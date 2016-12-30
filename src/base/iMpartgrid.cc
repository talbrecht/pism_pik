// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cmath>
#include <cstring>

#include "iceModel.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/pism_const.hh"
#include "coupler/PISMOcean.hh"
#include "earth/PISMBedDef.hh"
#include "base/util/pism_utilities.hh"
#include "base/part_grid_threshold_thickness.hh"

namespace pism {


//! @file iMpartgrid.cc Methods implementing geometry evolution with PIK option -part_grid [@ref
//! Albrechtetal2011].

//! Redistribute residual ice mass from subgrid-scale parameterization.
/*!
  See [\ref Albrechtetal2011].  Manages the loop.

  FIXME: Reporting!

  FIXME: repeatRedist should be config flag?

  FIXME: resolve fixed number (=3) of loops issue
*/
void IceModel::residual_redistribution(IceModelVec2S &H_residual) {
  const int max_loopcount = 3;

  bool done = false;
  for (int i = 0; not done and i < max_loopcount; ++i) {
    residual_redistribution_iteration(H_residual, done);
    m_log->message(4, "redistribution loopcount = %d\n", i);
  }
}


//! @brief Perform one iteration of the residual mass redistribution.
/**
 * @param[in,out] H_residual Residual Ice thickness. Updated in place.
 * @param[out] done set to 'true' if this was the last iteration we needed
 *
 * @return 0 on success
 */
void IceModel::residual_redistribution_iteration(IceModelVec2S &H_residual, bool &done) {

  bool reduce_frontal_thickness = m_config->get_boolean("geometry.part_grid.reduce_frontal_thickness");
  const double thickness_threshold = m_config->get_double("geometry.ice_free_thickness_standard");
  const double sea_level = m_ocean->sea_level_elevation();

  const IceModelVec2S &bed_topography = m_beddef->bed_elevation();

  GeometryCalculator gc(*m_config);
  gc.set_icefree_thickness(thickness_threshold);

  gc.compute_mask(sea_level, bed_topography, m_ice_thickness, m_cell_type);

  // First step: distribute residual ice thickness
  {
    // will be destroyed at the end of the block
    IceModelVec::AccessList list{&m_cell_type, &m_ice_thickness, &m_Href, &H_residual};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (H_residual(i,j) <= 0.0) {
        continue;
      }

      StarStencil<int> m = m_cell_type.int_star(i,j);
      int N = 0; // number of empty or partially filled neighbors
      StarStencil<bool> neighbors;
      neighbors.set(false);

      if (mask::ice_free_ocean(m.e)) {
        N++;
        neighbors.e = true;
      }
      if (mask::ice_free_ocean(m.w)) {
        N++;
        neighbors.w = true;
      }
      if (mask::ice_free_ocean(m.n)) {
        N++;
        neighbors.n = true;
      }
      if (mask::ice_free_ocean(m.s)) {
        N++;
        neighbors.s = true;
      }

      if (N > 0)  {
        // Remaining ice mass will be redistributed equally among all
        // adjacent partially-filled cells (is there a more physical
        // way?)
        if (neighbors.e) {
          m_Href(i + 1, j) += H_residual(i, j) / N;
        }
        if (neighbors.w) {
          m_Href(i - 1, j) += H_residual(i, j) / N;
        }
        if (neighbors.n) {
          m_Href(i, j + 1) += H_residual(i, j) / N;
        }
        if (neighbors.s) {
          m_Href(i, j - 1) += H_residual(i, j) / N;
        }

        H_residual(i, j) = 0.0;
      } else {
        // Conserve mass, but (possibly) create a "ridge" at the shelf
        // front
        m_ice_thickness(i, j) += H_residual(i, j);
        H_residual(i, j) = 0.0;
      }

    }
  }

  m_ice_thickness.update_ghosts();

  // The loop above updated ice_thickness, so we need to re-calculate the mask:
  gc.compute_mask(sea_level, bed_topography, m_ice_thickness, m_cell_type);
  // and the surface elevation:
  gc.compute_surface(sea_level, bed_topography, m_ice_thickness, m_ice_surface_elevation);

  double
    remaining_residual_thickness        = 0.0,
    remaining_residual_thickness_global = 0.0,
    dx = m_grid->dx();

  // Second step: we need to redistribute residual ice volume if
  // neighbors which gained redistributed ice also become full.
  {
    // will be destroyed at the end of the block
    IceModelVec::AccessList list{&m_ice_thickness, &m_ice_surface_elevation,
        &bed_topography, &m_cell_type};
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_Href(i,j) <= 0.0) {
        continue;
      }

      double H_threshold = part_grid_threshold_thickness(m_cell_type.int_star(i, j),
                                                         m_ice_thickness.star(i, j),
                                                         m_ice_surface_elevation.star(i, j),
                                                         bed_topography(i,j),
                                                         dx,
                                                         reduce_frontal_thickness);

      double coverage_ratio = 1.0;
      if (H_threshold > 0.0) {
        coverage_ratio = m_Href(i, j) / H_threshold;
      }
      if (coverage_ratio >= 1.0) {
        // The current partially filled grid cell is considered to be full
        H_residual(i, j) = m_Href(i, j) - H_threshold;
        remaining_residual_thickness += H_residual(i, j);
        m_ice_thickness(i, j) += H_threshold;
        m_Href(i, j) = 0.0;
      }
      if (m_ice_thickness(i, j) < 0) {
        m_log->message(1,
                   "PISM WARNING: at i=%d, j=%d, we just produced negative ice thickness.\n"
                   "  H_threshold: %f\n"
                   "  coverage_ratio: %f\n"
                   "  m_Href: %f\n"
                   "  H_residual: %f\n"
                   "  ice_thickness: %f\n", i, j, H_threshold, coverage_ratio,
                   m_Href(i, j), H_residual(i, j), m_ice_thickness(i, j));
      }

    }
  }

  // check if redistribution should be run once more
  remaining_residual_thickness_global = GlobalSum(m_grid->com, remaining_residual_thickness);

  if (remaining_residual_thickness_global > 0.0) {
    done = false;
  } else {
    done = true;
  }

  m_ice_thickness.update_ghosts();
}

} // end of namespace pism
