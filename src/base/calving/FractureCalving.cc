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
#include "base/util/PISMVars.hh"
#include "base/util/pism_options.hh"

namespace pism {
namespace calving {

FractureCalving::FractureCalving(IceGrid::ConstPtr g,
                           stressbalance::StressBalance *stress_balance)
  : StressCalving(g, stress_balance, 2) {

  m_K = m_config->get_double("calving.eigen_calving.K");


  //const unsigned int stencil_width = m_config->get_double("grid.max_stencil_width");
  calv_rate.create(m_grid, "fracture_calving_rate",
              WITHOUT_GHOSTS);
  //              WITH_GHOSTS, stencil_width);
  calv_rate.set_attrs("internal",
                 "potential fracture calving rate",
                 "m year-1", ""); 

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


PMC_potential_calving_rate::PMC_potential_calving_rate(const FractureCalving *m)
  : Diag<FractureCalving>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, "fracture_calving_rate_potential")};
  set_attrs("potential fracture calving rate",
            "potential_fracture_calving_rate", 
            "m year-1", "m year-1", 0);
}

IceModelVec::Ptr PMC_potential_calving_rate::compute_impl() {
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "fracture_calving_rate_potential", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  result->copy_from(model->calving_rate());

  return result;
}


std::map<std::string, Diagnostic::Ptr> FractureCalving::diagnostics_impl() const {

  return {
          {"fracture_calving_rate",Diagnostic::Ptr(new CalvingRate(this, "fracture_calving_rate",
                                             "horizontal calving rate due to fracture-calving"))},
          {"fracture_calving_rate_potential", Diagnostic::Ptr(new PMC_potential_calving_rate(this))},
         };

}

//const IceModelVec2S& FractureCalving::calving_rate(IceModelVec2S &calv_rate) const{
const IceModelVec2S& FractureCalving::calving_rate() const{

  //compute_calving_option(calvrate);

  return calv_rate;
}


void FractureCalving::compute_calving_option(IceModelVec2S &result) const{

  using std::max;

  // Distance (grid cells) from calving front where strain rate is evaluated
  int offset = m_stencil_width;

  const double eigenCalvOffset = 0.0;

  //double seconds_per_year = convert(m_unit_system, 1.0, "year", "seconds");
  double seconds_per_year = 3.15569259747e7;

  double Knew = options::Real("-eigen_calving_K","Eigencalving constant used in Fracture Calving",m_K);

  int Foption = options::Integer("-fracture_calving_opt","Fracture Calving Option",0);

  double F1 = options::Real("-fracture_calving_K1","Fracture Calving constant 1",0.0);
  double F2a = options::Real("-fracture_calving_K2a","Fracture Calving constant 2a",0.0);
  double F2b = options::Real("-fracture_calving_K2b","Fracture Calving constant 2b",0.0);
  double F3 = options::Real("-fracture_calving_K3","Fracture Calving constant 3",0.0);


  update_strain_rates();

  const IceModelVec2S &D = *m_grid->variables().get_2d_scalar("fracture_density");
  const IceModelVec2CellType &mask         = *m_grid->variables().get_2d_cell_type("mask");


  IceModelVec::AccessList list{&mask, &m_strain_rates, &D};
  list.add(result);

  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    if (mask.floating_ice(i, j) or mask.ice_free_ocean(i, j)) {
    //if (mask.ice_free_ocean(i, j) and mask.next_to_floating_ice(i, j)) {

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



      ///////////////////////////////////////////////////////////////
      // option 0: Eigen Calving law
      //
      // eigen1 * eigen2 has units [s^-2] and calving_rate_horizontal
      // [m*s^1] hence, eigen_calving_K has units [m*s]
      if (eigen2 > eigenCalvOffset and eigen1 > 0.0) {
        // spreading in all directions
        result(i, j) = Knew * eigen1 * (eigen2 - eigenCalvOffset);
      } else {
        result(i, j) = 0.0;
      }

      ///////////////////////////////////////////////////////////////
      // Fracture calving options
      if (Foption == 1) {
        if (fdens > 0.0){
          // option 1
          //m_log->message(2, "!!!! fracdens=%f at (%d, %d)\n", fdens, i, j);
          //m_log->message(2, "!!!! eigen calving=%f, additional calving=%f at (%d, %d)\n", result(i, j)*seconds_per_year,F1*fdens, i, j);
          result(i, j) += F1 * fdens / seconds_per_year;
        }

      /////////////////////////////////////////////////////////////////
      } else if (Foption == 2) {
          // option 2
          result(i, j) = F2a + (F2b - F2a) * fdens / seconds_per_year;

      /////////////////////////////////////////////////////////////////
      } else if (Foption == 3) {
          // option 3

          if (eigen2 > eigenCalvOffset and eigen1 > 0.0) {
            // spreading in all directions
            result(i, j) = (F3*fdens + Knew) * eigen1 * (eigen2 - eigenCalvOffset);
          } else {
            result(i, j) = 0.0;
          }
      } 


    } else { // end of "if (ice_free_ocean and next_to_floating)"
      result(i, j) = 0.0;
    }

  }   // end of loop over grid points

  //return result;
}

void FractureCalving::compute_calving_rate(const IceModelVec2CellType &mask,
                                          IceModelVec2S &result) const {

  compute_calving_option(result);

  /*
  IceModelVec::AccessList list{&mask, &result};

  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    bool calving_at_front = (mask.ice_free_ocean(i, j) and mask.next_to_floating_ice(i, j));

    if (!calving_at_front) {
      result(i, j) = 0.0;
    }
  }*/
}


} // end of namespace calving
} // end of namespace pism
