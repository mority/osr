#pragma once

#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "geo/latlng.h"

#include "osr/types.h"

namespace osr {

// Area routing, layer by layer:
//   0  what a walkable area is and what blocks it       this file
//   1  does it need area routing, or do mapped ways do?  osr/area/served.h
//   2  exact shortest paths, the ground truth            osr/area/geodesic.h
//   3  subdivision into cells                            osr/area/cells.h
//   4  routing through the cells                         osr/area/area_graph.h,
//                                                        profiles/foot.h
//   5  drawing a crossing as its shortest line           osr/area/crossing.h
//                                                        (reuses layer 2)
//
// Layer 0 decides from OSM data what a walkable area is, what stands in the
// way inside it, and cleans its geometry so that everything built on top can
// take "inside" literally.
//
// None of this is as simple as the tags make it look. Every rule here comes
// from a real case where the data meant something other than it seemed to,
// and area_walkable_test.cc keeps those cases as tests.

// Looks up the value of an OSM tag; nullopt if the object has no such key.
using tag_lookup =
    std::function<std::optional<std::string_view>(std::string_view key)>;

// Which levels something is on, and its layer.
//
// An object with no `level` tag matches everything, as foot routing already
// treats kNoLevel: most of OSM is untagged, and reading that as "ground floor
// only" would disconnect far more than it fixes.
struct area_levels {
  static area_levels of(tag_lookup const&);

  // `*this` is the area. Under the permissive rule an untagged object joins
  // anything, which is right outdoors and dangerous indoors: at Chatelet a
  // footway stub with no level tag would attach to an area three floors below
  // it. The strict rule makes an area that states its level demand the same
  // from whatever connects to it.
  bool matches(area_levels const& o, bool strict) const;

  // For things joined only by geometry - a way whose nodes fall inside an
  // area without sharing any of its nodes. A footbridge over a plaza lies
  // inside its polygon and connects to nothing, and `layer` is what says so.
  // Deliberately not applied where a node is genuinely shared: a mapper who
  // joined a staircase to a plaza at one node meant them to connect,
  // whatever layer the staircase carries.
  bool matches_spatially(area_levels const& o, bool strict) const;

  // Same levels and layer: areas that may be merged.
  bool same_as(area_levels const& o) const;

  void merge(area_levels const& o);

  bool any_{false};
  level_bits_t bits_{0U};
  int layer_{0};
};

struct walkable_options {
  // A place=square with no surface tag counts as walkable. Off: a bare
  // place=square names a square, it says nothing about where one can walk -
  // Hansaplatz has Altonaer Strasse running through it.
  bool place_square_{false};
};

// Is this area - a closed way (`from_way`) or a multipolygon - ground that
// pedestrians can cross in any direction?
bool is_walkable_area(tag_lookup const&,
                      bool from_way,
                      walkable_options const& = {});

// A way that carries routing in the ordinary way. Its nodes are what areas
// connect to.
bool is_linear_way(tag_lookup const&);

// Something inside an area that walking cannot pass: a building or a
// barrier.
struct blocker {
  enum class kind : std::uint8_t { kNone, kBuilding, kBarrier };

  // Does it block an area on `area`'s levels?
  bool stands_in(area_levels const& area, bool strict_levels) const;

  kind kind_{kind::kNone};
  area_levels levels_{};

  // Buildings: the lowest level the building occupies. Above 0 it stands on
  // columns - an upper floor over the plaza - and the ground underneath is
  // open.
  float min_level_{0.F};
};

blocker blocker_of(tag_lookup const&);

// The holes of an area - its inner rings and the buildings standing in it -
// as one clean set. Everything downstream decides inside and outside by the
// even-odd rule, and that rule turns overlap into walkable ground: a building
// that is both an inner ring and a building way cancels itself out, two
// overlapping buildings cancel where they overlap, and a building reaching
// past the outline makes the ground outside walkable. So the holes are
// unioned and clipped to the outline. A courtyard inside a building (a hole
// in a hole) cannot be reached from the area and is dropped, as is a hole
// that is not a valid polygon.
std::vector<std::vector<geo::latlng>> clean_holes(
    std::vector<geo::latlng> const& outer,
    std::vector<std::vector<geo::latlng>> const& holes);

struct area_ring {
  std::vector<std::int64_t> ids_{};  // OSM node ids
  std::vector<geo::latlng> points_{};
};

struct walkable_area {
  std::int64_t id_{0};
  bool from_way_{false};
  std::string name_{};
  area_levels levels_{};
  std::vector<area_ring> rings_{};  // rings_[0] is the outline

  // Set when overlapping areas were merged into this one: the areas it was
  // made of, as "way/1" / "relation/2", and the outline nodes of those areas
  // that now lie inside it.
  std::vector<std::string> members_{};
  std::vector<std::pair<std::int64_t, geo::latlng>> dissolved_{};
};

struct merge_stats {
  int n_members_{0};  // areas that took part in a merge
  int n_merged_{0};  // areas those became
  int n_failed_{0};  // groups left unmerged because the union failed
};

// The same ground is often mapped more than once: a place=square over the
// highway=pedestrian surface, a small area around a station entrance inside
// a plaza. Areas on the same levels and layer whose interiors overlap are
// unioned into one area. Areas that merely touch along an edge stay apart:
// they are joined through the connectors they share.
//
// A vertex that survives the union keeps its node id; points the union
// creates, where two outlines cross, get negative ids no way can have.
merge_stats merge_overlapping(std::vector<walkable_area>&);

}  // namespace osr
