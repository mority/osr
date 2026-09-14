#include "osr/backend/area_routing.h"

#include <algorithm>
#include <map>
#include <utility>

namespace json = boost::json;

namespace osr::backend {

namespace {

geo::latlng to_latlng(json::value const& c) {
  auto const& a = c.as_array();
  return {a[1].to_number<double>(), a[0].to_number<double>()};
}

std::vector<geo::latlng> to_line(json::value const& coords) {
  auto out = std::vector<geo::latlng>{};
  for (auto const& c : coords.as_array()) {
    out.push_back(to_latlng(c));
  }
  return out;
}

std::string_view str(json::object const& props, char const* key) {
  auto const* v = props.if_contains(key);
  return v != nullptr && v->is_string() ? std::string_view{v->get_string()}
                                        : std::string_view{};
}

}  // namespace

area_routing::area_routing(ways const& w, json::array const& features) {
  // Every feature carries the index of the area it belongs to.
  auto by_area = std::map<std::int64_t, std::vector<json::object const*>>{};
  for (auto const& f : features) {
    auto const& props = f.as_object().at("properties").as_object();
    if (auto const* idx = props.if_contains("area"); idx != nullptr) {
      by_area[idx->to_number<std::int64_t>()].push_back(&f.as_object());
    }
  }

  // Back from the file into what layer 3 produced, and the geometry layer 5
  // draws with.
  auto inputs = std::vector<area_crossing_input>{};
  auto geometries = std::vector<area_geometry>{};
  for (auto const& [_, fs] : by_area) {
    auto in = area_crossing_input{};
    auto geo = area_geometry{};
    auto meshed = false;
    auto levels = area_levels{.any_ = true};
    auto neighbours = std::vector<std::pair<std::size_t, std::size_t>>{};
    // Per-connector costs, if the file has them for every connector (older
    // files carry only the flat cost per cell).
    auto spokes = std::vector<float>{};
    auto link_costs = std::vector<float>{};
    auto every_spoke = true;
    for (auto const* f : fs) {
      auto const& props = f->at("properties").as_object();
      auto const& coords = f->at("geometry").as_object().at("coordinates");
      auto const kind = str(props, "kind");
      if (kind == "area") {
        meshed = str(props, "status") == "meshed";
        in.osm_ = geo.osm_ = str(props, "osm");
        for (auto const& r : coords.as_array()) {
          geo.rings_.push_back(to_line(r));
        }
        if (auto const* lv = props.if_contains("levels"); lv != nullptr) {
          levels.any_ = false;
          for (auto const& l : lv->as_array()) {
            levels.bits_ |= static_cast<level_bits_t>(1U)
                            << to_idx(level_t{l.to_number<float>()});
          }
        }
      } else if (kind == "hub") {
        auto const cell = props.at("cell").to_number<std::size_t>();
        in.cells_.cost_.resize(std::max(in.cells_.cost_.size(), cell + 1U));
        in.cells_.cost_[cell] = props.at("cost").to_number<float>();
      } else if (kind == "connector" && props.contains("cell")) {
        in.cells_.connector_cell_.push_back(
            static_cast<area_cells::cell_idx_t>(
                props.at("cell").to_number<std::size_t>()));
        in.connector_nodes_.push_back(props.at("node").to_number<std::int64_t>());
        in.connector_pos_.push_back(to_latlng(coords));
        if (auto const* s = props.if_contains("spoke"); s != nullptr) {
          spokes.push_back(s->to_number<float>());
        } else {
          every_spoke = false;
        }
      } else if (kind == "neighbour") {
        neighbours.emplace_back(props.at("from").to_number<std::size_t>(),
                                props.at("to").to_number<std::size_t>());
        auto const* c = props.if_contains("cost");
        link_costs.push_back(c != nullptr ? c->to_number<float>() : 0.F);
      } else if (kind == "direct") {
        in.cells_.direct_.push_back(
            {.a_ = props.at("from").to_number<std::uint32_t>(),
             .b_ = props.at("to").to_number<std::uint32_t>(),
             .cost_ = props.at("cost").to_number<float>()});
      } else if (kind == "barrier") {
        geo.barriers_.push_back(to_line(coords));
      }
    }
    if (!meshed || in.cells_.cost_.empty() || in.connector_nodes_.empty()) {
      continue;
    }
    for (auto const& [x, y] : neighbours) {
      in.cells_.set_neighbour(static_cast<area_cells::cell_idx_t>(x),
                              static_cast<area_cells::cell_idx_t>(y));
    }
    if (every_spoke && spokes.size() == in.connector_nodes_.size()) {
      in.cells_.spoke_ = std::move(spokes);
      for (auto i = std::size_t{0U}; i != neighbours.size(); ++i) {
        in.cells_.set_link(
            static_cast<area_cells::cell_idx_t>(neighbours[i].first),
            static_cast<area_cells::cell_idx_t>(neighbours[i].second),
            link_costs[i]);
      }
    }
    in.level_ = area_floor(levels);
    geo.connectors_ = in.connector_pos_;
    inputs.push_back(std::move(in));
    geometries.push_back(std::move(geo));
  }

  graph_ = std::make_unique<area_graph>(
      w.n_nodes(), inputs, [&](std::int64_t const osm) {
        return w.find_node_idx(to_osm_node_idx(osm));
      });
  drawer_ = std::make_unique<crossing_drawer>(*graph_, std::move(geometries));
  n_connectors_ = graph_->n_connectors_;
  n_without_routing_node_ = graph_->n_without_routing_node_;

  // The medial axis of every area, for the comparison the UI draws. Cheap
  // next to everything else here - all of Berlin's take well under a second.
  for (auto i = std::size_t{0U}; i != drawer_->n_areas(); ++i) {
    auto const& g = drawer_->geometry(i);
    skeletons_.push_back(std::make_unique<area_skeleton>(g.rings_, g.barriers_));
  }
}

}  // namespace osr::backend
