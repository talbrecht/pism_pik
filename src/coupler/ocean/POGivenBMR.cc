// Copyright (C) 2011, 2012 Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#include "POGivenBMR.hh"
#include "IceGrid.hh"
#include "PISMVars.hh"
#include "pism_options.hh"



POGivenBMR::POGivenBMRConstants::POGivenBMRConstants(const PISMConfig &config) {

  T0_water                         = config.get("water_melting_point_temperature"); // K
  sea_water_density                = config.get("sea_water_density");
  ice_density                      = config.get("ice_density");
  standard_gravity                 = config.get("standard_gravity");
  secpera                          = config.get("seconds_per_year");
  beta_CC                          = config.get("beta_CC");
  beta_CC_grad                     = beta_CC *  ice_density * standard_gravity; // K m-1

  // coefficients of the exponential fit for the yearly shelfbmassflux 

  c[0] = -0.03337955 * ice_density;
  c[1] =  0.02736375 * ice_density;
  c[2] = -0.02269549 / ice_density;

}


POGivenBMR::POGivenBMR(IceGrid &g, const PISMConfig &conf)
  : PGivenClimate<POModifier,PISMOceanModel>(g, conf, NULL),
    shelfbtemp(g.get_unit_system())
    //shelfbmassflux(g.get_unit_system())
{
  PetscErrorCode ierr = allocate_POGivenBMR(); CHKERRCONTINUE(ierr);
  if (ierr != 0)
    PISMEnd();

}

POGivenBMR::~POGivenBMR() {
  // empty
}

PetscErrorCode POGivenBMR::allocate_POGivenBMR() {
  PetscErrorCode ierr;

  option_prefix   = "-ocean_bmr";

  // will be de-allocated by the parent's destructor
  shelfbmassflux = new IceModelVec2T;

  m_fields["shelfbmassflux"] = shelfbmassflux;

  ierr = process_options(); CHKERRQ(ierr);

  std::map<std::string, std::string> standard_names;
  ierr = set_vec_parameters(standard_names); CHKERRQ(ierr);

  ierr = shelfbmassflux->create(grid, "shelfbmassflux", false); CHKERRQ(ierr);
  ierr = shelfbmassflux->set_attrs("climate_forcing","ice mass flux from ice shelf base (positive flux is loss from ice shelf)","kg m-2 s-1", ""); CHKERRQ(ierr);
  ierr = shelfbmassflux->set_glaciological_units("kg m-2 year-1"); CHKERRQ(ierr);
  shelfbmassflux->write_in_glaciological_units = true;
  
  //ierr = shelfbtemp.create(grid, "shelfbtemp", WITHOUT_GHOSTS); CHKERRQ(ierr);
  //ierr = shelfbtemp.set_attrs("climate_forcing","absolute temperature at ice shelf base","Kelvin", ""); CHKERRQ(ierr);
  shelfbtemp.init_2d("shelfbtemp", grid);
  shelfbtemp.set_string("pism_intent", "climate_state");
  shelfbtemp.set_string("long_name","absolute temperature at ice shelf base");
  shelfbtemp.set_units("Kelvin");


  return 0;
}

PetscErrorCode POGivenBMR::init(PISMVars &vars) {
  PetscErrorCode ierr;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  bool regrid = true;
  int start = -1;

  ierr = verbPrintf(2, grid.com,"* Initializing the ocean model 'BMR' (which reads 'shelfbmassflux' from file) ...\n"); CHKERRQ(ierr);

  ierr = shelfbmassflux->init(filename, bc_period, bc_reference_time); CHKERRQ(ierr);

  //ierr = process_options(); CHKERRQ(ierr);

  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!ice_thickness) { SETERRQ(grid.com, 1, "ERROR: ice thickness is not available"); }

  // FIXME: update does not work, commented version runs
  // read time-independent data right away:
  if (shelfbmassflux->get_n_records() == 1) {
    ierr = update(grid.time->current(), 0); CHKERRQ(ierr); // dt is irrelevant
  }

  ierr = verbPrintf(2, grid.com,"* Sub-shelf mass flux will be adjusted according to reference ice shelf base elevation"); CHKERRQ(ierr); 

  std::string filenametemp=filename;
  ierr = find_pism_input(filename, regrid, start); CHKERRQ(ierr);
  ierr = melt_ref_thk.create(grid, "melt_ref_thk", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = melt_ref_thk.set_attrs("model_state","reference ice geometry","m",""); CHKERRQ(ierr); // no CF standard_name ??

  // read reference ice geometry from file
  ierr = verbPrintf(2, grid.com,"  - Reading reference ice geometry ('melt_ref_thk') from '%s' ... \n",filename.c_str()); CHKERRQ(ierr);
  if (regrid) {
    ierr = melt_ref_thk.regrid(filename.c_str(), CRITICAL); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = melt_ref_thk.read(filename.c_str(), start); CHKERRQ(ierr); // fails if not found!
  }
  std::string ref_shelfbaseelev_history = "read from " + filename + "\n";
  //ierr = melt_ref_thk.set_string("history", ref_shelfbaseelev_history); CHKERRQ(ierr); //FIXME: melt_ref_thk is not NCSpatialVariable
  //ierr = melt_ref_thk.set_attr("history", ref_shelfbaseelev_history); CHKERRQ(ierr);
  filename=filenametemp;
  

  return 0;
}

PetscErrorCode POGivenBMR::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr = update_internal(my_t, my_dt); CHKERRQ(ierr);

  ierr = shelfbmassflux->average(m_t, m_dt); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode POGivenBMR::shelf_base_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;

  POGivenBMRConstants c(config);

  ierr = ice_thickness->begin_access();   CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const double pressure = c.ice_density * c.standard_gravity * (*ice_thickness)(i,j); // FIXME task #7297

      // temp is set to melting point at depth
      result(i,j) = c.T0_water - c.beta_CC * pressure;
    }
  }
  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode POGivenBMR::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

  POGivenBMRConstants c(config);

  PetscScalar **H;
  ierr = shelfbmassflux->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = melt_ref_thk.begin_access(); CHKERRQ(ierr);
  ierr = ice_thickness->get_array(H);   CHKERRQ(ierr);


  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      PetscReal reference_thickness = melt_ref_thk(i,j),
      shelfbaseelev     = - ( c.ice_density / c.sea_water_density ) * H[i][j], // FIXME issue #15
      ref_shelfbaseelev = - ( c.ice_density / c.sea_water_density ) * reference_thickness;

      // first order correction to melt rate with
      // bmr(z) = bmr0(z) + dbmr/dz * (z-z0)
      // db/dz is a function of bmr predominantely, we use an exponential fit here
      // parameters for yearly melt rates
      PetscReal dbmrdz = c.c[0] + c.c[1] * exp( c.c[2] * (*shelfbmassflux)(i,j) * c.secpera );
      result(i,j) = (*shelfbmassflux)(i,j) + dbmrdz/c.secpera * (shelfbaseelev - ref_shelfbaseelev) ;

    }
  }

  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = melt_ref_thk.end_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux->end_access(); CHKERRQ(ierr);

  // PetscErrorCode ierr = shelfbmassflux->copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode POGivenBMR::sea_level_elevation(double &result) {
  result = sea_level;
  return 0;
}


PetscErrorCode POGivenBMR::melange_back_pressure_fraction(IceModelVec2S &result) {
  PetscErrorCode ierr = result.set(0.0); CHKERRQ(ierr);
  return 0;
}


void POGivenBMR::add_vars_to_output(std::string keyword, std::set<std::string> &result) {

  PGivenClimate<POModifier,PISMOceanModel>::add_vars_to_output(keyword, result);

  if (keyword != "none" && keyword != "small") {
    result.insert("shelfbtemp");
  }
}

PetscErrorCode POGivenBMR::define_variables(std::set<std::string> vars, const PIO &nc, PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  ierr = PGivenClimate<POModifier,PISMOceanModel>::define_variables(vars, nc, nctype); CHKERRQ(ierr);

  if (set_contains(vars, "shelfbtemp")) {
    ierr = shelfbtemp.define(nc, nctype, true); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode POGivenBMR::write_variables(std::set<std::string> vars, const PIO &nc) {
  PetscErrorCode ierr;
  IceModelVec2S tmp;

  ierr = PGivenClimate<POModifier,PISMOceanModel>::write_variables(vars, nc); CHKERRQ(ierr);


  if (set_contains(vars, "shelfbtemp")) {
    if (!tmp.was_created()) {
      ierr = tmp.create(grid, "tmp", WITHOUT_GHOSTS); CHKERRQ(ierr);
    }

    tmp.metadata() = shelfbtemp;
    ierr = shelf_base_temperature(tmp); CHKERRQ(ierr);
    ierr = tmp.write(nc); CHKERRQ(ierr);
  }
  

  ierr = melt_ref_thk.write(nc); CHKERRQ(ierr);

  return 0;
}
