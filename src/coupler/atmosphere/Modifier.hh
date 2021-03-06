// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017 PISM Authors
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

#ifndef _PAMODIFIER_H_
#define _PAMODIFIER_H_

#include "pism/coupler/AtmosphereModel.hh"
#include "pism/util/Modifier.hh"

namespace pism {
namespace atmosphere {

class PAModifier : public Modifier<AtmosphereModel>
{
public:
  PAModifier(IceGrid::ConstPtr g, AtmosphereModel* in)
    : Modifier<AtmosphereModel>(g, in) {}
  virtual ~PAModifier() {}

protected:
  virtual void mean_precipitation_impl(IceModelVec2S &result) const
  {
    if (m_input_model != NULL) {
      m_input_model->mean_precipitation(result);
    }
  }

  virtual void mean_annual_temp_impl(IceModelVec2S &result) const
  {
    if (m_input_model != NULL) {
      m_input_model->mean_annual_temp(result);
    }
  }

  virtual void begin_pointwise_access_impl() const
  {
    if (m_input_model != NULL) {
      m_input_model->begin_pointwise_access();
    }
  }

  virtual void end_pointwise_access_impl() const
  {
    if (m_input_model != NULL) {
      m_input_model->end_pointwise_access();
    }
  }

  virtual void temp_time_series_impl(int i, int j, std::vector<double> &result) const
  {
    if (m_input_model != NULL) {
      m_input_model->temp_time_series(i, j, result);
    }
  }

  virtual void precip_time_series_impl(int i, int j, std::vector<double> &result) const
  {
    if (m_input_model != NULL) {
      m_input_model->precip_time_series(i, j, result);
    }
  }

  virtual void init_timeseries_impl(const std::vector<double> &ts) const
  {
    if (m_input_model != NULL) {
      m_input_model->init_timeseries(ts);
    }

    m_ts_times = ts;
  }
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* _PAMODIFIER_H_ */
