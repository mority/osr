#include "osr/backend/http_server.h"

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
#include "osr/routing/reach_cells.h"
#include "osr/routing/route.h"

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
       std::string const& static_file_path)
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

    auto const p = route(params, w_, l_, profile, from, to, max, dir, 100,
                         nullptr, nullptr, elevations_, routing_algo);

    auto const p1 = route(params, w_, l_, profile, from, std::vector{to}, max,
                          dir, 100, nullptr, nullptr, elevations_);

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
    cb(json_response(req, to_featurecollection(w_, p)));
  }

  void handle_cells(web_server::http_req_t const& req,
                    web_server::http_res_cb_t const& cb) {
    auto const q = boost::json::parse(req.body()).as_object();
    auto const profile_it = q.find("profile");
    auto const profile =
        profile_it == q.end() || !profile_it->value().is_string()
            ? search_profile::kCar
            : to_profile(profile_it->value().as_string());
    auto const sources =
        utl::to_vec(q.at("sources").as_array(),
                    [](json::value const& v) { return parse_location(v); });

    auto p = reach_cells_params{};
    if (auto const it = q.find("max"); it != q.end()) {
      p.max_ = static_cast<cost_t>(it->value().to_number<std::int64_t>());
    }
    if (auto const it = q.find("cellSize"); it != q.end()) {
      p.cell_size_m_ = it->value().to_number<double>();
    }
    if (auto const it = q.find("fillSteps"); it != q.end()) {
      p.fill_steps_ = it->value().to_number<unsigned>();
    }
    auto const get_minutes = [&](char const* key) {
      auto const it = q.find(key);
      return it == q.end() || !it->value().is_number()
                 ? std::nullopt
                 : std::optional{duration_t{static_cast<duration_t::rep>(
                       std::clamp(it->value().to_number<double>() * 60.0, 0.0,
                                  18.0 * 3600.0))}};
    };
    p.meet_margin_ = get_minutes("meetMarginMinutes");
    p.total_margin_ = get_minutes("totalMarginMinutes");
    if (auto const it = q.find("bandMinutes"); it != q.end()) {
      p.band_size_ = duration_t{static_cast<duration_t::rep>(std::clamp(
          it->value().to_number<double>() * 60.0, 0.0, 18.0 * 3600.0))};
    }
    utl::verify(p.cell_size_m_ >= 10.0, "cellSize must be >= 10m");

    auto const cells = compute_reach_cells(get_parameters(profile), w_, l_,
                                           profile, sources, p, 100);

    auto features = json::array{};
    auto const add = [&](reach_polygon const& polygon, json::object props) {
      auto rings = json::array{};
      for (auto const& ring : polygon) {
        // 1e-5 degrees ~ 1m: plenty for cells, halves the response size
        auto coords = json::array{};
        coords.reserve(ring.size());
        for (auto const& x : ring) {
          coords.emplace_back(json::array{std::round(x.lng() * 1e5) / 1e5,
                                          std::round(x.lat() * 1e5) / 1e5});
        }
        rings.emplace_back(std::move(coords));
      }
      features.emplace_back(json::object{
          {"type", "Feature"},
          {"properties", std::move(props)},
          {"geometry",
           {{"type", "Polygon"}, {"coordinates", std::move(rings)}}}});
    };
    // Bands first, so the outlines are drawn on top within one layer order.
    for (auto const& b : cells.bands_) {
      add(b.polygon_, {{"kind", "band"},
                       {"source_idx", b.source_},
                       {"band", b.band_},
                       {"from_min", b.band_ * p.band_size_.count() / 60.0}});
    }
    for (auto const [source, polygons] : utl::enumerate(cells.sources_)) {
      for (auto const& polygon : polygons) {
        add(polygon, {{"kind", "outline"}, {"source_idx", source}});
      }
    }
    auto const add_meet = [&](std::optional<reach_cells::meet> const& m,
                              char const* zone_kind, char const* point_kind) {
      if (!m.has_value()) {
        return;
      }
      for (auto const& polygon : m->zone_) {
        add(polygon, {{"kind", zone_kind}});
      }
      auto per_source = json::array{};
      for (auto const d : m->per_source_) {
        per_source.emplace_back(d.count() / 60.0);
      }
      features.emplace_back(json::object{
          {"type", "Feature"},
          {"properties",
           {{"kind", point_kind},
            {"latest_min", m->latest_ / 60.0},
            {"total_min", m->total_ / 60.0},
            {"per_source_min", std::move(per_source)}}},
          {"geometry",
           {{"type", "Point"},
            {"coordinates", json::array{m->pos_.lng(), m->pos_.lat()}}}}});
    };
    add_meet(cells.meet_, "meet_zone", "meet_point");
    add_meet(cells.total_meet_, "total_zone", "total_point");
    cb(json_response(
        req, json::serialize(json::object{{"type", "FeatureCollection"},
                                          {"features", std::move(features)}})));
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
    auto const min =
        geo::latlng{waypoints[1].as_double(), waypoints[0].as_double()};
    auto const max =
        geo::latlng{waypoints[3].as_double(), waypoints[2].as_double()};

    auto gj = geojson_writer{.w_ = w_};
    l_.find({min, max}, [&](way_idx_t const w) { gj.write_way(w); });

    gj.finish();
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
        } else if (target.starts_with("/api/cells")) {
          return run_parallel(
              [this](web_server::http_req_t const& req1,
                     web_server::http_res_cb_t const& cb1) {
                handle_cells(req1, cb1);
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
                         std::string const& static_file_path)
    : impl_{new impl(ioc, thread_pool, w, l, pl, elevation, static_file_path)} {
}

http_server::~http_server() = default;

void http_server::listen(std::string const& host, std::string const& port) {
  impl_->listen(host, port);
}

void http_server::stop() { impl_->stop(); }

}  // namespace osr::backend
