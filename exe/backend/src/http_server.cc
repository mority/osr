#include "osr/backend/http_server.h"

#include <fstream>
#include <sstream>
#include <utility>

#include "boost/algorithm/string.hpp"
#include "boost/asio/post.hpp"
#include "boost/beast/core/string.hpp"
#include "boost/beast/version.hpp"
#include "boost/json.hpp"

#include "fmt/core.h"

#include "utl/enumerate.h"
#include "utl/pipes.h"
#include "utl/to_vec.h"

#include "net/web_server/responses.h"
#include "net/web_server/serve_static.h"
#include "net/web_server/web_server.h"

#include "osr/backend/area_routing.h"
#include "osr/geojson.h"
#include "osr/lookup.h"
#include "osr/routing/algorithms.h"
#include "osr/routing/parameters.h"
#include "osr/routing/profiles/bike.h"
#include "osr/routing/profiles/bike_sharing.h"
#include "osr/routing/profiles/car.h"
#include "osr/routing/profiles/car_parking.h"
#include "osr/routing/profiles/car_sharing.h"
#include "osr/routing/profiles/foot.h"
#include "osr/routing/route.h"
#include "osr/routing/with_profile.h"

using namespace net;
using net::web_server;

namespace http = boost::beast::http;
namespace fs = std::filesystem;
namespace json = boost::json;

namespace osr::backend {

template <typename Body>
void set_cors_headers(http::response<Body>& res) {
  using namespace boost::beast::http;
  res.set(field::access_control_allow_origin, "*");
  res.set(field::access_control_allow_headers,
          "X-Content-Type-Options, X-Requested-With, Content-Type, Accept, "
          "Authorization");
  res.set(field::access_control_allow_methods, "GET, POST, OPTIONS");
  res.set(field::access_control_max_age, "3600");
}

web_server::string_res_t json_response(
    web_server::http_req_t const& req,
    std::string const& content,
    http::status const status = http::status::ok) {
  auto res = net::string_response(req, content, status, "application/json");
  set_cors_headers(res);
  return res;
}

location parse_location(json::value const& v) {
  auto const& obj = v.as_object();
  return {obj.at("lat").as_double(), obj.at("lng").as_double(),
          obj.contains("level") ? level_t{obj.at("level").to_number<float>()}
                                : kNoLevel};
}

json::value to_json(std::vector<geo::latlng> const& polyline) {
  auto a = json::array{};
  for (auto const& p : polyline) {
    a.emplace_back(json::array{p.lng(), p.lat()});
  }
  return a;
}

struct http_server::impl {
  impl(boost::asio::io_context& ios,
       boost::asio::io_context& thread_pool,
       ways const& g,
       lookup const& l,
       platforms const* pl,
       elevation_storage const* elevations,
       std::string const& static_file_path,
       fs::path const& area_cells_file)
      : ioc_{ios},
        thread_pool_{thread_pool},
        w_{g},
        l_{l},
        pl_{pl},
        elevations_{elevations},
        server_{ioc_} {
    try {
      if (!static_file_path.empty() && fs::is_directory(static_file_path)) {
        static_file_path_ = fs::canonical(static_file_path).string();
        serve_static_files_ = true;
      }
    } catch (fs::filesystem_error const& e) {
      throw utl::fail("static file directory not found: {}", e.what());
    }
    load_areas(area_cells_file);
  }

  // Area cells as written by `osr-area-stats --cells-out`, grouped by area so
  // a request returns whole areas: every feature carries its area's "osm" id.
  void load_areas(fs::path const& p) {
    if (p.empty() || !fs::exists(p)) {
      return;
    }
    auto const content = [&]() {
      auto in = std::ifstream{p};
      auto ss = std::stringstream{};
      ss << in.rdbuf();
      return ss.str();
    }();

    auto by_osm = hash_map<std::string, std::size_t>{};
    auto const extend = [](auto&& self, geo::box& b,
                           json::value const& c) -> void {
      auto const& a = c.as_array();
      if (a.size() >= 2U && a[0].is_number()) {
        b.extend(geo::latlng{a[1].to_number<double>(),
                             a[0].to_number<double>()});
      } else {
        for (auto const& x : a) {
          self(self, b, x);
        }
      }
    };
    // Named, not iterated as a temporary: a range-for only extends the
    // lifetime of the last temporary in its range expression, so the parsed
    // document would be gone before the first iteration.
    auto doc = json::parse(content);
    area_routing_ = std::make_unique<area_routing>(
        w_, doc.as_object().at("features").as_array());
    std::cout << "area routing: " << area_routing_->n_areas() << " areas, "
              << area_routing_->n_hubs() << " hubs, "
              << area_routing_->n_connectors_ << " connectors, of which "
              << area_routing_->n_without_routing_node_
              << " are no routing node (left out)\n";
    for (auto& f : doc.as_object().at("features").as_array()) {
      auto const& props = f.as_object().at("properties").as_object();
      auto const osm = std::string{props.at("osm").as_string()};
      auto const [it, inserted] = by_osm.emplace(osm, areas_.size());
      if (inserted) {
        areas_.emplace_back();
      }
      auto& area = areas_[it->second];
      if (auto const* lv = props.if_contains("levels"); lv != nullptr) {
        area.levels_ = *lv;
      }
      extend(extend, area.bbox_,
             f.as_object().at("geometry").as_object().at("coordinates"));
      area.features_.emplace_back(std::move(f));
    }

    // The skeletons are held by area_routing, which groups the pieces of a
    // merged area separately; the features here are grouped by osm id.
    for (auto i = std::size_t{0U}; i != area_routing_->n_areas(); ++i) {
      if (auto const it = by_osm.find(area_routing_->osm_of(i));
          it != end(by_osm)) {
        areas_[it->second].skeletons_.push_back(i);
      }
    }
    std::cout << "loaded " << areas_.size() << " areas with cells from " << p
              << '\n';
  }

  void handle_areas(web_server::http_req_t const& req,
                    web_server::http_res_cb_t const& cb) {
    auto const query = boost::json::parse(req.body()).as_object();
    auto const waypoints = query.at("waypoints").as_array();
    auto const view = geo::box{
        geo::latlng{waypoints[1].as_double(), waypoints[0].as_double()},
        geo::latlng{waypoints[3].as_double(), waypoints[2].as_double()}};
    auto features = json::array{};
    for (auto const& a : areas_) {
      if (!a.bbox_.overlaps(view)) {
        continue;
      }
      for (auto const& f : a.features_) {
        features.push_back(f);
      }
      // The medial axis of the same area, for the Skeleton layer.
      for (auto const i : a.skeletons_) {
        auto const& osm = area_routing_->osm_of(i);
        for (auto const& [from, to] : area_routing_->skeleton(i).edges()) {
          auto props = json::object{{"kind", "skeleton"}, {"osm", osm}};
          if (!a.levels_.is_null()) {
            props["levels"] = a.levels_;
          }
          features.push_back(json::object{
              {"type", "Feature"},
              {"properties", std::move(props)},
              {"geometry",
               json::object{
                   {"type", "LineString"},
                   {"coordinates",
                    json::array{json::array{from.lng(), from.lat()},
                                json::array{to.lng(), to.lat()}}}}}});
        }
      }
    }
    cb(json_response(req, json::serialize(json::object{
                              {"type", "FeatureCollection"},
                              {"features", std::move(features)}})));
  }

  static search_profile get_search_profile_from_request(
      boost::json::object const& q) {
    auto const profile_it = q.find("profile");
    return profile_it == q.end() || !profile_it->value().is_string()
               ? search_profile::kFoot
               : to_profile(profile_it->value().as_string());
  }

  static routing_algorithm get_routing_algorithm_from_request(
      boost::json::object const& q) {
    auto const routing_it = q.find("routing");
    return routing_it == q.end() || !routing_it->value().is_string()
               ? routing_algorithm::kDijkstra
               : to_algorithm(routing_it->value().as_string());
  }

  void handle_route(web_server::http_req_t const& req,
                    web_server::http_res_cb_t const& cb) {
    auto const q = boost::json::parse(req.body()).as_object();
    auto const profile = get_search_profile_from_request(q);
    auto const direction_it = q.find("direction");
    auto const routing_algo = get_routing_algorithm_from_request(q);
    auto const dir = to_direction(direction_it == q.end() ||
                                          !direction_it->value().is_string()
                                      ? to_str(direction::kForward)
                                      : direction_it->value().as_string());
    auto const from = parse_location(q.at("start"));
    auto const to = parse_location(q.at("destination"));
    auto const max_it = q.find("max");
    auto const max = static_cast<cost_t>(
        max_it == q.end() ? 3600 : max_it->value().as_int64());
    auto const foot_speed_result =
        q.try_at("footSpeed")->try_to_number<float>();
    auto const params =
        profile == search_profile::kFoot && foot_speed_result.has_value()
            ? foot<false,
                   elevator_tracking>::parameters{.speed_meters_per_second_ =
                                                      foot_speed_result.value()}
            : get_parameters(profile);

    // Walking profiles cross areas through their cells.
    auto const areas =
        area_routing_ != nullptr && !area_routing_->empty() &&
                (profile == search_profile::kFoot ||
                 profile == search_profile::kWheelchair)
            ? std::optional{area_routing_->sharing()}
            : std::nullopt;
    auto const* sharing = areas.has_value() ? &*areas : nullptr;

    auto const p = route(params, w_, l_, profile, from, to, max, dir, 100,
                         nullptr, sharing, elevations_, routing_algo);

    auto const p1 = route(params, w_, l_, profile, from, std::vector{to}, max,
                          dir, 100, nullptr, sharing, elevations_);

    auto const print = [](char const* name, std::optional<path> const& p) {
      if (p.has_value()) {
        std::cout << name << " cost: " << p->cost_ << "\n";
      } else {
        std::cout << name << ": not found\n";
      }
    };
    print("p", p);
    print("p1", p1.at(0));

    if (!p.has_value()) {
      cb(json_response(req, "could not find a valid path",
                       http::status::not_found));
      return;
    }
    // Both the path as routed - an area crossing as straight lines through
    // the hubs of its cells, marked "area" - and, marked "geodesic", the
    // shortest walkable line each crossing stands for.
    auto fc = to_featurecollection_value(w_, p);
    if (sharing != nullptr) {
      auto& features = fc.at("features").as_array();
      for (auto i = std::size_t{0U}; i != p->segments_.size(); ++i) {
        auto const& s = p->segments_[i];
        if (sharing->is_additional_node(s.from_) ||
            sharing->is_additional_node(s.to_)) {
          auto& props = features[i].as_object()["properties"].as_object();
          props["area"] = true;
          // Path reconstruction puts every segment without a way on level
          // 0; an area's segments are on its floor.
          props["level"] =
              area_routing_
                  ->level_of(sharing->is_additional_node(s.to_) ? s.to_
                                                                : s.from_)
                  .to_float();
        }
      }
      for (auto const& c : area_routing_->crossings(*p)) {
        features.push_back(json::object{
            {"type", "Feature"},
            {"properties",
             {{"geodesic", true},
              {"level", c.level_.to_float()},
              {"area_osm", area_routing_->osm_of(c.area_)},
              {"model_distance", c.model_distance_},
              {"geodesic_distance", c.geodesic_distance_}}},
            {"geometry", to_line_string(c.geodesic_)}});
        // The same crossing over the medial axis: what the other approach
        // would have the walker do.
        if (auto const over_skeleton = area_routing_->skeleton_path(c);
            over_skeleton.size() >= 2U) {
          auto length = 0.0;
          for (auto i = std::size_t{1U}; i != over_skeleton.size(); ++i) {
            length += geo::distance(over_skeleton[i - 1U], over_skeleton[i]);
          }
          features.push_back(json::object{
              {"type", "Feature"},
              {"properties",
               {{"skeleton", true},
                {"level", c.level_.to_float()},
                {"area_osm", area_routing_->osm_of(c.area_)},
                {"skeleton_distance", length}}},
              {"geometry", to_line_string(over_skeleton)}});
        }
      }
    }
    cb(json_response(req, json::serialize(fc)));
  }

  void handle_levels(web_server::http_req_t const& req,
                     web_server::http_res_cb_t const& cb) {
    auto const query = boost::json::parse(req.body()).as_object();
    auto const waypoints = query.at("waypoints").as_array();
    auto const min = point::from_latlng(
        {waypoints[1].as_double(), waypoints[0].as_double()});
    auto const max = point::from_latlng(
        {waypoints[3].as_double(), waypoints[2].as_double()});
    auto levels = hash_set<level_t>{};
    l_.find({min, max}, [&](way_idx_t const x) {
      auto const p = w_.r_->way_properties_[x];
      levels.emplace(p.from_level());
      if (p.from_level() != p.to_level()) {
        levels.emplace(p.to_level());
      }
    });
    auto levels_sorted =
        utl::to_vec(levels, [](level_t const l) { return l.to_float(); });
    utl::sort(levels_sorted, [](auto&& a, auto&& b) { return a > b; });
    cb(json_response(req,
                     json::serialize(utl::all(levels_sorted)  //
                                     | utl::emplace_back_to<json::array>())));
  }

  void handle_graph(web_server::http_req_t const& req,
                    web_server::http_res_cb_t const& cb) {
    auto const query = boost::json::parse(req.body()).as_object();
    auto const waypoints = query.at("waypoints").as_array();
    auto const profile = get_search_profile_from_request(query);
    auto const min =
        geo::latlng{waypoints[1].as_double(), waypoints[0].as_double()};
    auto const max =
        geo::latlng{waypoints[3].as_double(), waypoints[2].as_double()};

    auto gj = geojson_writer{.w_ = w_};
    l_.find({min, max}, [&](way_idx_t const w) { gj.write_way(w); });

    with_profile(profile,
                 [&]<Profile P>(P&&) { send_graph_response<P>(req, cb, gj); });
  }

  template <Profile P>
  void send_graph_response(web_server::http_req_t const& req,
                           web_server::http_res_cb_t const& cb,
                           geojson_writer& gj) {
    gj.finish(&get_dijkstra<P>());
    cb(json_response(req, gj.string()));
  }

  void handle_static(web_server::http_req_t const& req,
                     web_server::http_res_cb_t const& cb) {
    if (auto res = net::serve_static_file(
            boost::beast::string_view{static_file_path_}, req);
        res.has_value()) {
      cb(std::move(*res));
    } else {
      namespace http = boost::beast::http;
      cb(net::web_server::string_res_t{http::status::not_found, req.version()});
    }
  }

  void handle_platforms(web_server::http_req_t const& req,
                        web_server::http_res_cb_t const& cb) {
    utl::verify(pl_ != nullptr, "no platforms");

    auto const query = boost::json::parse(req.body()).as_object();
    auto const level = query.contains("level")
                           ? level_t{query.at("level").to_number<float>()}
                           : kNoLevel;
    auto const waypoints = query.at("waypoints").as_array();
    auto const min = point::from_latlng(
        {waypoints[1].as_double(), waypoints[0].as_double()});
    auto const max = point::from_latlng(
        {waypoints[3].as_double(), waypoints[2].as_double()});

    auto gj = geojson_writer{.w_ = w_, .platforms_ = pl_};
    pl_->find(min, max, [&](platform_idx_t const i) {
      if (level == kNoLevel || pl_->get_level(w_, i) == level) {
        gj.write_platform(i);
      }
    });

    cb(json_response(req, gj.string()));
  }

  void handle_request(web_server::http_req_t const& req,
                      web_server::http_res_cb_t const& cb) {
    std::cout << "[" << req.method_string() << "] " << req.target() << '\n';
    switch (req.method()) {
      case http::verb::options: return cb(json_response(req, {}));
      case http::verb::post: {
        auto const& target = req.target();
        if (target.starts_with("/api/route")) {
          return run_parallel(
              [this](web_server::http_req_t const& req1,
                     web_server::http_res_cb_t const& cb1) {
                handle_route(req1, cb1);
              },
              req, cb);
        } else if (target.starts_with("/api/levels")) {
          return run_parallel(
              [this](web_server::http_req_t const& req1,
                     web_server::http_res_cb_t const& cb1) {
                handle_levels(req1, cb1);
              },
              req, cb);
        } else if (target.starts_with("/api/graph")) {
          return run_parallel(
              [this](web_server::http_req_t const& req1,
                     web_server::http_res_cb_t const& cb1) {
                handle_graph(req1, cb1);
              },
              req, cb);
        } else if (target.starts_with("/api/areas")) {
          return run_parallel(
              [this](web_server::http_req_t const& req1,
                     web_server::http_res_cb_t const& cb1) {
                handle_areas(req1, cb1);
              },
              req, cb);
        } else if (target.starts_with("/api/platforms")) {
          return run_parallel(
              [this](web_server::http_req_t const& req1,
                     web_server::http_res_cb_t const& cb1) {
                handle_platforms(req1, cb1);
              },
              req, cb);
        } else {
          return cb(json_response(req, R"({"error": "Not found"})",
                                  http::status::not_found));
        }
      }
      case http::verb::get:
      case http::verb::head: return handle_static(req, cb);
      default:
        return cb(json_response(req,
                                R"({"error": "HTTP method not supported"})",
                                http::status::bad_request));
    }
  }

  template <typename Fn>
  void run_parallel(
      Fn&& handler,  // NOLINT(cppcoreguidelines-missing-std-forward)
      web_server::http_req_t const& req,
      web_server::http_res_cb_t const& cb) {
    boost::asio::post(
        thread_pool_, [req, cb, h = std::forward<Fn>(handler), this]() {
          try {
            h(req, [req, cb, this](web_server::http_res_t&& res) {
              boost::asio::post(ioc_, [cb, req, res{std::move(res)}]() mutable {
                try {
                  cb(std::move(res));
                } catch (std::exception const& e) {
                  return cb(json_response(
                      req, fmt::format(R"({{"error": "{}"}})", e.what()),
                      http::status::internal_server_error));
                }
              });
            });
          } catch (std::exception const& e) {
            return cb(json_response(
                req, fmt::format(R"({{"error": "{}"}})", e.what()),
                http::status::internal_server_error));
          }
        });
  }

  void listen(std::string const& host, std::string const& port) {
    server_.on_http_request(
        [this](web_server::http_req_t const& req,
               web_server::http_res_cb_t const& cb,
               bool /*ssl*/) { return handle_request(req, cb); });

    boost::system::error_code ec;
    server_.init(host, port, ec);
    if (ec) {
      std::cerr << "server init error: " << ec.message() << "\n";
      return;
    }

    std::cout << "Listening on http://" << host << ":" << port
              << "/ and https://" << host << ":" << port << "/" << '\n';
    if (host == "0.0.0.0") {
      std::cout << "Local link: http://127.0.0.1:" << port << "/" << '\n';
    }
    server_.run();
  }

  void stop() { server_.stop(); }

private:
  boost::asio::io_context& ioc_;
  boost::asio::io_context& thread_pool_;
  ways const& w_;
  lookup const& l_;
  platforms const* pl_;
  elevation_storage const* elevations_;

  struct area_features {
    json::value levels_;  // the area's levels, for its skeleton features
    std::vector<std::size_t> skeletons_;  // its pieces in area_routing_
    geo::box bbox_;
    json::array features_;
  };
  std::vector<area_features> areas_;
  std::unique_ptr<area_routing> area_routing_;

  web_server server_;
  bool serve_static_files_{false};
  std::string static_file_path_;
};

http_server::http_server(boost::asio::io_context& ioc,
                         boost::asio::io_context& thread_pool,
                         ways const& w,
                         lookup const& l,
                         platforms const* pl,
                         elevation_storage const* elevation,
                         std::string const& static_file_path,
                         fs::path const& area_cells_file)
    : impl_{new impl(ioc,
                     thread_pool,
                     w,
                     l,
                     pl,
                     elevation,
                     static_file_path,
                     area_cells_file)} {}

http_server::~http_server() = default;

void http_server::listen(std::string const& host, std::string const& port) {
  impl_->listen(host, port);
}

void http_server::stop() { impl_->stop(); }

}  // namespace osr::backend
