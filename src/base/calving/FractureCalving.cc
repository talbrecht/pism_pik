/* Copyright (C) 2016, 2017 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "FractureCalving.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/error_handling.hh"
#include "base/util/IceModelVec2CellType.hh"
#include "remove_narrow_tongues.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/stressbalance/ShallowStressBalance.hh"
#include "base/util/PISMVars.hh"
#include "base/rheology/FlowLaw.hh"
#include "base/util/pism_options.hh"

namespace pism {
namespace calving {

FractureCalving::FractureCalving(IceGrid::ConstPtr g,
                           stressbalance::StressBalance *stress_balance)
  : StressCalving(g, stress_balance, 2) {

  m_K = m_config->get_double("calving.eigen_calving.K");
}

FractureCalving::~FractureCalving() {
  // empty
}

void FractureCalving::init() {

  m_log->message(2,
                 "* Initializing the 'Fracture calving' mechanism...\n");

  if (fabs(m_grid->dx() - m_grid->dy()) / std::min(m_grid->dx(), m_grid->dy()) > 1e-2) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "-calving fracture_calving using a non-square grid cell is not implemented (yet);\n"
                                  "dx = %f, dy = %f, relative difference = %f",
                                  m_grid->dx(), m_grid->dy(),
                                  fabs(m_grid->dx() - m_grid->dy()) / std::max(m_grid->dx(), m_grid->dy()));
  }

  m_strain_rates.set(0.0);

}

//! \brief Uses principal strain rates to apply "eigencalving" with constant K.
/*!
  See equation (26) in [\ref Winkelmannetal2011].
*/
void FractureCalving::compute_calving_rate(const IceModelVec2CellType &mask,
                                          IceModelVec2S &result) const {

  using std::max;

  // Distance (grid cells) from calving front where strain rate is evaluated
  int offset = m_stencil_width;

  const double eigenCalvOffset = 0.0;

  double m_Knew = options::Real("-eigen_calving_K","Eigencalving constant used in Fracture Calving",m_K);

  double m_Fnew = options::Real("-fracture_calving_K","Fracturecalving constant",0.0);

  update_strain_rates();

  const IceModelVec2S &D = *m_grid->variables().get_2d_cell_type("fracture_density");

  IceModelVec::AccessList list{&mask, &result, &m_strain_rates, &D};


  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    // Find partially filled or empty grid boxes on the icefree ocean, which
    // have floating ice neighbors after the mass continuity step
    if (mask.ice_free_ocean(i, j) and mask.next_to_floating_ice(i, j)) {

      // Average of strain-rate eigenvalues in adjacent floating grid cells.
      double
        eigen1             = 0.0,
        eigen2             = 0.0,
        fdens              = 0.0;
      {
        int N = 0;
        for (int p = -1; p < 2; p += 2) {
          const int I = i + p * offset;

            if (mask.floating_ice(I, j) and not mask.ice_margin(I, j)) {
              eigen1 += m_strain_rates(I, j, 0);
              eigen2 += m_strain_rates(I, j, 1);
              fdens  += D(I, j);
              N += 1;
            }

        }

        for (int q = -1; q < 2; q += 2) {
          const int J = j + q * offset;

            if (mask.floating_ice(i, J) and not mask.ice_margin(i, J)) {
              eigen1 += m_strain_rates(i, J, 0);
              eigen2 += m_strain_rates(i, J, 1);
              fdens  += D(i, J);
              N += 1;
            }
        }

        if (N > 0) {
          eigen1             /= N;
          eigen2             /= N;
          fdens              /= N;
        }
      }

 
      // Eigen Calving law
      //
      // eigen1 * eigen2 has units [s^-2] and calving_rate_horizontal
      // [m*s^1] hence, eigen_calving_K has units [m*s]
      if (eigen2 > eigenCalvOffset and eigen1 > 0.0) {
        // spreading in all directions
        //result(i, j) = m_K * eigen1 * (eigen2 - eigenCalvOffset);
        result(i, j) = m_Knew * eigen1 * (eigen2 - eigenCalvOffset);
      } else {
        result(i, j) = 0.0;
      }

      if (fdens > 0.0){
        m_log->message(2,
                           "!!!! fracdens=%f at (%d, %d)\n",
                           fdens, i, j);
        result(i, j) = 100000.0;
      }


    } else { // end of "if (ice_free_ocean and next_to_floating)"
      result(i, j) = 0.0;
    }
  }   // end of loop over grid points
}

std::map<std::string, Diagnostic::Ptr> FractureCalving::diagnostics_impl() const {
  return {{"fracture_calving_rate",
        Diagnostic::Ptr(new CalvingRate(this, "fracture_calving_rate",
                                        "horizontal calving rate due to Fracture calving"))}};
}

} // end of namespace calving
} // end of namespace pism
