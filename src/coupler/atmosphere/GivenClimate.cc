// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include "GivenClimate.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"

namespace pism {
namespace atmosphere {

Given::Given(IceGrid::ConstPtr g)
  : PGivenClimate<AtmosphereModel>(g, nullptr)
{
  m_fields["precipitation"] = nullptr;
  m_fields["air_temp"]      = nullptr;

  process_options("-atmosphere_given");

  set_vec_parameters(m_filename, {});

  m_precipitation = m_fields["precipitation"].get();
  m_air_temp      = m_fields["air_temp"].get();

  {
    m_air_temp->set_attrs("diagnostic", "mean annual near-surface air temperature",
                          "Kelvin", "", 0);
    m_air_temp->metadata(0).set_double("valid_min", 0.0);
    m_air_temp->metadata(0).set_double("valid_max", 323.15); // 50 C
  }
  {
    m_precipitation->set_attrs("model_state", "precipitation rate",
                               "kg m-2 second-1", "precipitation_flux", 0);
    m_precipitation->metadata(0).set_string("glaciological_units", "kg m-2 year-1");
  }
}

Given::~Given() {
  // empty
}

void Given::init_impl(const Geometry &geometry) {
  m_log->message(2,
             "* Initializing the atmosphere model reading near-surface air temperature\n"
             "  and ice-equivalent precipitation from a file...\n");

  m_air_temp->init(m_filename, m_bc_period, m_bc_reference_time);
  m_precipitation->init(m_filename, m_bc_period, m_bc_reference_time);

  // read time-independent data right away:
  if (m_air_temp->n_records() == 1 && m_precipitation->n_records() == 1) {
    update(geometry, m_grid->ctx()->time()->current(), 0); // dt is irrelevant
  }
}

void Given::update_impl(const Geometry &geometry, double t, double dt) {
  update_internal(geometry, t, dt);

  m_precipitation->average(t, dt);
  m_air_temp->average(t, dt);
}

const IceModelVec2S& Given::mean_precipitation_impl() const {
  return *m_precipitation;
}

const IceModelVec2S& Given::mean_annual_temp_impl() const {
  return *m_air_temp;
}

void Given::begin_pointwise_access_impl() const {

  m_air_temp->begin_access();
  m_precipitation->begin_access();
}

void Given::end_pointwise_access_impl() const {

  m_air_temp->end_access();
  m_precipitation->end_access();
}

void Given::temp_time_series_impl(int i, int j, std::vector<double> &result) const {

  m_air_temp->interp(i, j, result);
}

void Given::precip_time_series_impl(int i, int j, std::vector<double> &result) const {

  m_precipitation->interp(i, j, result);
}

void Given::init_timeseries_impl(const std::vector<double> &ts) const {

  m_air_temp->init_interpolation(ts);

  m_precipitation->init_interpolation(ts);

  m_ts_times = ts;
}


} // end of namespace atmosphere
} // end of namespace pism
