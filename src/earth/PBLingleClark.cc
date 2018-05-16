// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2017 Constantine Khroulev
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

#include "PBLingleClark.hh"

#include <gsl/gsl_math.h>       // GSL_NAN

#include "base/util/io/PIO.hh"
#include "base/util/PISMTime.hh"
#include "base/util/IceGrid.hh"
#include "base/util/pism_options.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/error_handling.hh"
#include "base/util/PISMVars.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/pism_utilities.hh"
#include "deformation.hh"

namespace pism {
namespace bed {

PBLingleClark::PBLingleClark(IceGrid::ConstPtr g)
  : BedDef(g) {

  // A work vector. This storage is used to put thickness change on rank 0 and to get the plate
  // displacement change back.
  m_bed_displacement.create(m_grid, "bed_displacement", WITHOUT_GHOSTS);
  m_bed_displacement.set_attrs("internal",
                               "total (viscous and elastic) displacement "
                               "in the Lingle-Clark bed deformation model",
                               "meters", "");

  m_work0 = m_bed_displacement.allocate_proc0_copy();

  m_relief.create(m_grid, "bed_relief", WITHOUT_GHOSTS);
  m_relief.set_attrs("internal",
                     "bed relief relative to the modeled bed displacement",
                     "meters", "");

  bool use_elastic_model = m_config->get_boolean("bed_deformation.lc.elastic_model");

  m_bdLC = NULL;

  const int
    Mx = m_grid->Mx(),
    My = m_grid->My(),
    Z  = 4,                     // use Z = 4 for now; to reduce global drift?
    Nx = Z*(Mx - 1) + 1,
    Ny = Z*(My - 1) + 1;

  const double
    Lx = Z * (m_grid->x0() - m_grid->x(0)),
    Ly = Z * (m_grid->y0() - m_grid->y(0));

  m_extended_grid = IceGrid::Shallow(m_grid->ctx(),
                                     Lx, Ly,
                                     m_grid->x0(), m_grid->y0(),
                                     Nx, Ny, NOT_PERIODIC);

  m_viscous_bed_displacement.create(m_extended_grid,
                                    "viscous_bed_displacement", WITHOUT_GHOSTS);
  m_viscous_bed_displacement.set_attrs("model state",
                                 "bed displacement in the viscous half-space bed deformation model; "
                                 "see BuelerLingleBrown", "meters", "");
  // coordinate variables of the extended grid should have different names
  m_viscous_bed_displacement.metadata().get_x().set_name("x_lc");
  m_viscous_bed_displacement.metadata().get_y().set_name("y_lc");
  m_viscous_bed_displacement0 = m_viscous_bed_displacement.allocate_proc0_copy();

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {
      m_bdLC = new BedDeformLC(*m_config, use_elastic_model,
                               Mx, My,
                               m_grid->dx(), m_grid->dy(),
                               Nx, Ny);
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();
}

PBLingleClark::~PBLingleClark() {
  if (m_bdLC != NULL) {
    delete m_bdLC;
    m_bdLC = NULL;
  }
}

/*!
 * Initialize the model by computing the viscous bed displacement using uplift and the elastic
 * response using ice thickness.
 *
 * Then compute the bed relief as the difference between bed elevation and total bed displacement.
 */
void PBLingleClark::bootstrap_impl(const IceModelVec2S &bed,
                                   const IceModelVec2S &bed_uplift,
                                   const IceModelVec2S &ice_thickness) {
  m_t_beddef_last = m_grid->ctx()->time()->start();

  m_t  = GSL_NAN;
  m_dt = GSL_NAN;

  m_topg_last.copy_from(bed);

  //petsc::Vec::Ptr thickness0 = ice_thickness.allocate_proc0_copy();

  const IceModelVec2S *load_thickness = m_grid->variables().get_2d_scalar("beddef_load");
  petsc::Vec::Ptr thickness0 = load_thickness->allocate_proc0_copy();


  // initialize the plate displacement
  {
    bed_uplift.put_on_proc0(*m_work0);
    //ice_thickness.put_on_proc0(*thickness0);
    load_thickness->put_on_proc0(*thickness0);


    ParallelSection rank0(m_grid->com);
    try {
      if (m_grid->rank() == 0) {
        PetscErrorCode ierr = 0;

        m_bdLC->bootstrap(*thickness0, *m_work0);

        ierr = VecCopy(m_bdLC->total_displacement(), *m_work0);
        PISM_CHK(ierr, "VecCopy");

        ierr = VecCopy(m_bdLC->viscous_displacement(), *m_viscous_bed_displacement0);
        PISM_CHK(ierr, "VecCopy");
      }
    } catch (...) {
      rank0.failed();
    }
    rank0.check();
  }

  m_viscous_bed_displacement.get_from_proc0(*m_viscous_bed_displacement0);

  m_bed_displacement.get_from_proc0(*m_work0);

  // compute bed relief
  m_topg.add(-1.0, m_bed_displacement, m_relief);
}

/*! Initialize the Lingle-Clark bed deformation model.
 *
 * Inputs:
 *
 * - bed topography,
 * - ice thickness,
 * - plate displacement (either read from a file or bootstrapped using uplift) and
 *   possibly re-gridded.
 */
void PBLingleClark::init_impl(const InputOptions &opts) {
  m_log->message(2, "* Initializing the Lingle-Clark bed deformation model...\n");

  // Initialize bed topography and uplift maps.
  BedDef::init_impl(opts);

  const IceModelVec2S *load_thickness = m_grid->variables().get_2d_scalar("beddef_load");
  //const IceModelVec2S *ice_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  m_topg_last.copy_from(m_topg);

  if (opts.type == INIT_RESTART) {
    // Set m_viscous_bed_displacement by reading from the input file.
    m_viscous_bed_displacement.read(opts.filename, opts.record);
  } else if (opts.type == INIT_BOOTSTRAP) {
    this->bootstrap(m_topg, m_uplift, *load_thickness);
  } else {
    // do nothing
  }

  // Try re-gridding plate_displacement.
  regrid("Lingle-Clark bed deformation model",
         m_viscous_bed_displacement, REGRID_WITHOUT_REGRID_VARS);

  // Now that m_viscous_bed_displacement is finally initialized, put it on rank 0 and initialize
  // m_bdLC itself.
  {
    load_thickness->put_on_proc0(*m_work0);
    m_viscous_bed_displacement.put_on_proc0(*m_viscous_bed_displacement0);

    ParallelSection rank0(m_grid->com);
    try {
      if (m_grid->rank() == 0) {  // only processor zero does the step
        PetscErrorCode ierr = 0;

        m_bdLC->init(*m_work0, *m_viscous_bed_displacement0);

        ierr = VecCopy(m_bdLC->total_displacement(), *m_work0);
        PISM_CHK(ierr, "VecCopy");
      }
    } catch (...) {
      rank0.failed();
    }
    rank0.check();
  }

  m_bed_displacement.get_from_proc0(*m_work0);

  // compute bed relief
  m_topg.add(-1.0, m_bed_displacement, m_relief);
}

MaxTimestep PBLingleClark::max_timestep_impl(double t) const {
  (void) t;
  // no time-step restriction
  return MaxTimestep("bed_def lc");
}

/*!
 * Get total bed displacement on the PISM grid.
 *
 * This method uses the fact that m_bed_displacement is used to store bed displacement
 */
const IceModelVec2S& PBLingleClark::total_displacement() const {
  return m_bed_displacement;
}

//! Update the Lingle-Clark bed deformation model.
void PBLingleClark::update_with_thickness_impl(const IceModelVec2S &load_thickness,
                                               double t, double dt) {

  if ((fabs(t - m_t)   < 1e-12) &&
      (fabs(dt - m_dt) < 1e-12)) {
    return;
  }

  m_t  = t;
  m_dt = dt;

  double t_final = m_t + m_dt;

  // Check if it's time to update:
  double dt_beddef = t_final - m_t_beddef_last; // in seconds
  if ((dt_beddef < m_config->get_double("bed_deformation.update_interval", "seconds") and
       t_final < m_grid->ctx()->time()->end()) or
      dt_beddef < 1e-12) {
    return;
  }

  m_t_beddef_last = t_final;

  load_thickness.put_on_proc0(*m_work0);

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {  // only processor zero does the step
      PetscErrorCode ierr = 0;

      m_bdLC->step(dt_beddef, *m_work0);

      ierr = VecCopy(m_bdLC->total_displacement(), *m_work0);
      PISM_CHK(ierr, "VecCopy");

      ierr = VecCopy(m_bdLC->viscous_displacement(), *m_viscous_bed_displacement0);
      PISM_CHK(ierr, "VecCopy");
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();

  m_viscous_bed_displacement.get_from_proc0(*m_viscous_bed_displacement0);

  m_bed_displacement.get_from_proc0(*m_work0);

  // Update bed elevation using bed displacement and relief.
  {
    m_bed_displacement.add(1.0, m_relief, m_topg);
    // Increment the topg state counter. SIAFD relies on this!
    m_topg.inc_state_counter();
  }

  //! Finally, we need to update bed uplift and topg_last.
  compute_uplift(dt_beddef);
  m_topg_last.copy_from(m_topg);
}

void PBLingleClark::define_model_state_impl(const PIO &output) const {
  BedDef::define_model_state_impl(output);
  m_viscous_bed_displacement.define(output);
}

void PBLingleClark::write_model_state_impl(const PIO &output) const {
  BedDef::write_model_state_impl(output);

  m_viscous_bed_displacement.write(output);
}

} // end of namespace bed
} // end of namespace pism
