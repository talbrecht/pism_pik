/* Copyright (C) 2013, 2014, 2015, 2016, 2017 PISM Authors
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

#ifndef _PBLINGLECLARK_H_
#define _PBLINGLECLARK_H_

#include "PISMBedDef.hh"

namespace pism {
namespace bed {

class BedDeformLC;

//! A wrapper class around BedDeformLC.
class PBLingleClark : public BedDef {
public:
  PBLingleClark(IceGrid::ConstPtr g);
  virtual ~PBLingleClark();

  const IceModelVec2S& total_displacement() const;

protected:
  virtual void define_model_state_impl(const PIO &output) const;
  virtual void write_model_state_impl(const PIO &output) const;

  MaxTimestep max_timestep_impl(double t) const;
  void init_impl(const InputOptions &opts);
  void bootstrap_impl(const IceModelVec2S &bed_elevation,
                      const IceModelVec2S &bed_uplift,
                      const IceModelVec2S &load_thickness);
  void update_with_thickness_impl(const IceModelVec2S &load_thickness,
                                  double my_t, double my_dt);

  //! Total (viscous and elastic) bed displacement.
  IceModelVec2S m_bed_displacement;

  //! Storage on rank zero. Used to pass the load to the serial deformation model and get
  //! bed displacement back.
  petsc::Vec::Ptr m_work0;

  //! Bed relief relative to the bed displacement.
  IceModelVec2S m_relief;

  //! Serial viscoelastic bed deformation model.
  BedDeformLC *m_bdLC;

  //! extended grid for the viscous plate displacement
  IceGrid::Ptr m_extended_grid;

  //! Viscous displacement on the extended grid (part of the model state).
  IceModelVec2S m_viscous_bed_displacement;
  //! rank 0 storage using the extended grid
  petsc::Vec::Ptr m_viscous_bed_displacement0;
};

} // end of namespace bed
} // end of namespace pism

#endif /* _PBLINGLECLARK_H_ */
