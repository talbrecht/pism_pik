// Copyright (C) 2011 Constantine Khroulev
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

#ifndef _POGIVENBMR_H_
#define _POGIVENBMR_H_

//#include "PGivenBMR.hh"
#include "PGivenClimate.hh"
#include "POModifier.hh"
//#include "NCVariable.hh"

//class POGivenBMR : public PGivenBMR<POModifier,PISMOceanModel>
class POGivenBMR : public PGivenClimate<POModifier,PISMOceanModel>
{
public:

  POGivenBMR(IceGrid &g, const PISMConfig &conf);
  virtual ~POGivenBMR();


  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(double my_t, double my_dt);

  virtual PetscErrorCode sea_level_elevation(double &result) {
    result = sea_level;
    return 0;
  }

  virtual PetscErrorCode shelf_base_temperature(IceModelVec2S &result);

  virtual PetscErrorCode shelf_base_mass_flux(IceModelVec2S &result);

  virtual void add_vars_to_output(std::string keyword, std::set<std::string> &result);

  virtual PetscErrorCode define_variables(std::set<std::string> vars, const PIO &nc,PISM_IO_Type nctype);

  virtual PetscErrorCode write_variables(std::set<std::string> vars, const PIO &nc);


protected:
  IceModelVec2S *ice_thickness, melt_ref_thk;
  NCSpatialVariable shelfbtemp, shelfbmassflux;
  //IceModelVec2T *shelfbtemp, *shelfbmassflux;
  //IceModelVec2T *shelfbmassflux;


private:
  PetscErrorCode allocate_POGivenBMR();
};


#endif /* _POGIVENBMR_H_ */
