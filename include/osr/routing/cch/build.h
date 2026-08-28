#pragma once

#include "cista/memory_holder.h"

#include "osr/routing/cch/cch.h"
#include "osr/routing/cch/order.h"
#include "osr/ways.h"

namespace osr {

// Ways that motorized vehicles (car / bus / hgv) can use. The metric
// independent hierarchy is built once for this sub graph and is shared by all
// customizations of these profiles.
bool is_cch_way(way_properties const&);

// Turn restrictions that hold for every motorized profile. Only those may be
// used to prune shortcuts during the metric independent contraction, profile
// specific restrictions are applied during customization.
bool is_restricted_for_all(ways::routing const&,
                           node_idx_t,
                           way_pos_t from,
                           way_pos_t to);

// `n_threads == 0` uses all cores. Only the contraction of dense nodes is
// spread over threads, everything else is sequential.
cista::wrapped<cch> build_cch(ways const&, unsigned n_threads = 0U);

}  // namespace osr
