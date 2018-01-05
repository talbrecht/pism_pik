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

#include <algorithm>            // std::min
#include <gsl/gsl_math.h>

#include "TemperatureIndexITM.hh"
#include "localITM.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/Vars.hh"
#include "pism/util/Time.hh"
#include "pism/coupler/AtmosphereModel.hh"
#include "pism/util/Mask.hh"
#include "pism/util/io/PIO.hh"

#include "pism/util/error_handling.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {
namespace surface {

namespace diagnostics {

/*! @brief Snow cover depth. */
class ITM_snow_depth : public Diag<TemperatureIndexITM>
{
public:
  ITM_snow_depth(const TemperatureIndexITM *m);
protected:
  IceModelVec::Ptr compute_impl() const;
};

ITM_snow_depth::ITM_snow_depth(const TemperatureIndexITM *m)
  : Diag<TemperatureIndexITM>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "snow_depth")};

  set_attrs("snow cover depth (set to zero once a year)", "",
            "m", "m", 0);
}

IceModelVec::Ptr ITM_snow_depth::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "snow_depth", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  result->copy_from(model->snow_depth());

  return result;
}


 }// end of namespace diagnostics

///// PISM surface model implementing a ITM scheme.

TemperatureIndexITM::TemperatureIndexITM(IceGrid::ConstPtr g)
  : SurfaceModel(g) {

  m_mbscheme                   = NULL;
  m_melt_conversion_factor     = 1.; //m_config->get_double("surface.itm.melt_conversion_factor");
  m_refreeze_fraction          = 0.5; //m_config->get_double("surface.itm.refreeze");

  m_mbscheme = new ITMMassBalance(m_config, m_sys);
  
  m_climatic_mass_balance.create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
  m_climatic_mass_balance.set_attrs("diagnostic",
                                    "instantaneous surface mass balance (accumulation/ablation) rate",
                                    "kg m-2 s-1",
                                    "land_ice_surface_specific_mass_balance_flux");
  m_climatic_mass_balance.metadata().set_string("glaciological_units", "kg m-2 year-1");
  m_climatic_mass_balance.metadata().set_string("comment", "positive values correspond to ice gain");

  // diagnostic fields:

  {
    m_accumulation.create(m_grid, "saccum", WITHOUT_GHOSTS);
    m_accumulation.set_attrs("diagnostic", "surface accumulation (precipitation minus rain)",
                             "kg m-2", "");

    m_melt.create(m_grid, "smelt", WITHOUT_GHOSTS);
    m_melt.set_attrs("diagnostic", "surface melt", "kg m-2", "");

    m_runoff.create(m_grid, "srunoff", WITHOUT_GHOSTS);
    m_runoff.set_attrs("diagnostic", "surface meltwater runoff",
                       "kg m-2", "");
  }

  m_snow_depth.create(m_grid, "snow_depth", WITHOUT_GHOSTS);
  m_snow_depth.set_attrs("diagnostic",
                         "snow cover depth (set to zero once a year)",
                         "m", "");
  m_snow_depth.set(0.0);

  m_firn_depth.create(m_grid, "firn_depth", WITHOUT_GHOSTS);
  m_firn_depth.set_attrs("diagnostic",
                         "firn cover depth",
                         "m", "");
  m_firn_depth.set(0.0);
}

TemperatureIndexITM::~TemperatureIndexITM() {
  delete m_mbscheme;
}

void TemperatureIndexITM::init_impl() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  // call the default implementation (not the interface method init())
  SurfaceModel::init_impl();

  // report user's modeling choices
  {
    m_log->message(2,
                   "* Initializing the ITM-based surface processes scheme.\n"
                   "  Precipitation and 2m air temperature provided by atmosphere are inputs.\n"
                   "  Surface mass balance and ice upper surface temperature are outputs.\n");
  }

  // initializing the model state
  InputOptions input = process_input_options(m_grid->com);

  std::string firn_file = m_config->get_string("surface.pdd.firn_depth_file");

  if (input.type == INIT_RESTART) {
    if (not firn_file.empty()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "surface.pdd.firn_depth_file is not allowed when"
                                    " re-starting from a PISM output file.");
    }

    m_firn_depth.read(input.filename, input.record);
    m_snow_depth.read(input.filename, input.record);
  } else if (input.type == INIT_BOOTSTRAP) {

    m_snow_depth.regrid(input.filename, OPTIONAL, 0.0);

    if (firn_file.empty()) {
      m_firn_depth.regrid(input.filename, OPTIONAL, 0.0);
    } else {
      m_firn_depth.regrid(firn_file, CRITICAL);
    }
  } else {

    m_snow_depth.set(0.0);

    if (firn_file.empty()) {
      m_firn_depth.set(0.0);
    } else {
      m_firn_depth.regrid(firn_file, CRITICAL);
    }
  }

  {
    regrid("ITM surface model", m_snow_depth);
    regrid("ITM surface model", m_firn_depth);
  }

  // finish up
  {
    m_next_balance_year_start = compute_next_balance_year_start(m_grid->ctx()->time()->current());

    m_accumulation.set(0.0);
    m_melt.set(0.0);
    m_runoff.set(0.0);
  }
}

MaxTimestep TemperatureIndexITM::max_timestep_impl(double my_t) const {
  return m_atmosphere->max_timestep(my_t);
}

double TemperatureIndexITM::compute_next_balance_year_start(double time) {
  // compute the time corresponding to the beginning of the next balance year
  double
    balance_year_start_day = m_config->get_double("surface.pdd.balance_year_start_day"),
    one_day                = units::convert(m_sys, 1.0, "days", "seconds"),
    year_start             = m_grid->ctx()->time()->calendar_year_start(time),
    balance_year_start     = year_start + (balance_year_start_day - 1.0) * one_day;

  if (balance_year_start > time) {
    return balance_year_start;
  }
  return m_grid->ctx()->time()->increment_date(balance_year_start, 1);
}

void TemperatureIndexITM::update_impl(double t, double dt) {
  m_log->message(2, "here, right at begining of update_impl\n");//FIXME
  if ((fabs(t - m_t) < 1e-12) &&
      (fabs(dt - m_dt) < 1e-12)) {
    return;
  }

  m_t  = t;
  m_dt = dt;

  // update to ensure that temperature and precipitation time series are correct:
  m_atmosphere->update(t, dt);

  // set up air temperature and precipitation time series
  int N = m_mbscheme->get_timeseries_length(dt);
  m_log->message(2, "N = %d \n",N); //FIXME

  const double dtseries = dt / N;
  m_log->message(2, "dtseries = %d \n",dtseries); //FIXME

  std::vector<double> ts(N), T(N), P(N);
  double ITM_melt = 0.;
  for (int k = 0; k < N; ++k) {
    ts[k] = t + k * dtseries;
  }

  double insolation = 900.;
  double surface_elevation = 0.;

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");


  IceModelVec::AccessList list{&mask, &m_climatic_mass_balance, &m_firn_depth, &m_snow_depth, &m_accumulation, &m_melt, &m_runoff}; 
  //TODO : include surface altitude. 

  const double melt_conversion_factor   = m_melt_conversion_factor ; // aus config oben.

  m_atmosphere->init_timeseries(ts);

  m_atmosphere->begin_pointwise_access();

  const double ice_density = m_config->get_double("constants.ice.density");
  ParallelSection loop(m_grid->com);
  m_log->message(2, "here, right before loop\n"); //FIXME
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // m_log->message(2, "i = %d, j = %d, mask = %f\n",i,j,mask(i,j)); //FIXME
      // // if (i == 100 && j == 48 ){

      // m_log->message(2, "*****************************************\n\ni = %d, j = %d, mask = %f\n",i,j,mask(i,j)); //FIXME
      // m_log->message(2, //FIXME: test
      //        "* get old model parameters \n"
      //        "  m_accumulation = %f, m_melt = %f, m_runoff = %f, m_climatic_mass_balance = %f\n",m_accumulation(i,j), m_melt(i, j), m_runoff(i, j), m_climatic_mass_balance(i, j));
      
      // reset total accumulation, melt, and runoff, and SMB
      {
        m_accumulation(i, j)          = 0.0;
        m_melt(i, j)                  = 0.0;
        m_runoff(i, j)                = 0.0;
        m_climatic_mass_balance(i, j) = 0.0;
      }

      // the temperature time series from the AtmosphereModel and its modifiers
      m_atmosphere->temp_time_series(i, j, T);

      // the precipitation time series from AtmosphereModel and its modifiers
      m_atmosphere->precip_time_series(i, j, P);

      // convert precipitation from "kg m-2 second-1" to "m second-1" (ITMMassBalance expects
      // accumulation in m/second ice equivalent)
      for (int k = 0; k < N; ++k) {
        P[k] = P[k] / ice_density;
        // kg / (m^2 * second) / (kg / m^3) = m / second
      }


      // Use temperature time series to remove rainfall from precipitation
      m_mbscheme->get_snow_accumulation(T,  // air temperature (input)
                                        P); // precipitation rate (input-output)

      // Use ITM, albedo, and the snow precipitation to get surface mass
      // balance (and diagnostics: accumulation, melt, runoff)
      {
        double next_snow_depth_reset = m_next_balance_year_start;

        for (int k = 0; k < N; ++k) {
          if (ts[k] >= next_snow_depth_reset) {
            m_snow_depth(i,j)       = 0.0;
            while (next_snow_depth_reset <= ts[k]) {
              next_snow_depth_reset = m_grid->ctx()->time()->increment_date(next_snow_depth_reset, 1);
            }
          }

          const double accumulation = P[k] * dtseries;

          // get albedo
          
          // m_log->message(2, //FIXME: test
          //        "  accumulation = %f\n m_melt = %f, m_snow_depth = %f, m_firn_depth = %f, mask = %f\n",accumulation,m_melt(i, j), m_snow_depth(i, j), m_firn_depth(i, j), mask(i, j));
          
          double albedo = m_mbscheme->get_albedo(m_melt(i, j), m_snow_depth(i, j), m_firn_depth(i, j), mask(i, j));
          
          // m_log->message(2, //FIXME: test
          //        "* get albedo \n"
          //        "  albedo = %f\n", albedo);

          // compute melt
          ITM_melt = m_mbscheme->calculate_ITM_melt(dtseries, insolation, T[k], surface_elevation, albedo);

          // m_log->message(2, //FIXME: test
          //        "* calculate ITM_melt = %f \n",ITM_melt
          //        );
          // compute changes in snow, firn, etc. 
          LocalMassBalanceITM::Changes changes;
          changes = m_mbscheme->step(melt_conversion_factor, m_refreeze_fraction, ITM_melt, m_firn_depth(i, j), m_snow_depth(i, j), accumulation);

          // update firn depth
          m_firn_depth(i, j) += changes.firn_depth;
          // update snow depth
          m_snow_depth(i, j) += changes.snow_depth;

          // m_log->message(2, //FIXME: test
          //        "* firn_depth = %f, snow_depth = %f\n* changes.firn_depth = %f, changes.snow_depth = %f\n",m_firn_depth(i, j),m_snow_depth(i,j),changes.firn_depth, changes.snow_depth
          //        );
          // update total accumulation, melt, and runoff, converting from "meters, ice equivalent"
          // to "kg / meter^2"
          {
            m_accumulation(i, j) += accumulation * ice_density;
            m_melt(i, j)         += changes.melt * ice_density;
            m_runoff(i, j)       += changes.runoff * ice_density;
          }

          // m_climatic_mass_balance (unlike m_accumulation, m_melt, and m_runoff), is a rate.
          // m * (kg / m^3) / second = kg / m^2 / second
          m_climatic_mass_balance(i, j) += changes.smb * ice_density / m_dt;
        } // end of the time-stepping loop
      }

      if (mask.ocean(i,j)) {
        m_firn_depth(i,j) = 0.0;  // no firn over the ocean
        m_snow_depth(i,j) = 0.0;  // snow over the ocean does not stick
      }
           // m_log->message(2,"here\n");
    // }//FIXME remove this one. 
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
  m_log->message(2,"here, right after loop\n");
  m_atmosphere->end_pointwise_access();

  m_next_balance_year_start = compute_next_balance_year_start(m_grid->ctx()->time()->current());
  m_log->message(2, "end of update_impl\n"); //FIXME
}

void TemperatureIndexITM::mass_flux_impl(IceModelVec2S &result) const {
  m_log->message(2, "mass_flux_impl \n"); //FIXME
  result.copy_from(m_climatic_mass_balance);
}

void TemperatureIndexITM::temperature_impl(IceModelVec2S &result) const {
  m_log->message(2, "temperature_impl \n"); //FIXME
  m_atmosphere->mean_annual_temp(result);
}

const IceModelVec2S& TemperatureIndexITM::accumulation() const {
  m_log->message(2, "sccumulation \n"); //FIXME
  return m_accumulation;
}

const IceModelVec2S& TemperatureIndexITM::melt() const {
  m_log->message(2, "melt \n"); //FIXME
  return m_melt;
}

const IceModelVec2S& TemperatureIndexITM::runoff() const {
  m_log->message(2, "runoff \n"); //FIXME
  return m_runoff;
}

const IceModelVec2S& TemperatureIndexITM::firn_depth() const {
  m_log->message(2, "firn_depth \n"); //FIXME
  return m_firn_depth;
}

const IceModelVec2S& TemperatureIndexITM::snow_depth() const {
  m_log->message(2, "snow_depth \n"); //FIXME
  return m_snow_depth;
}


void TemperatureIndexITM::define_model_state_impl(const PIO &output) const {
  SurfaceModel::define_model_state_impl(output);
  m_firn_depth.define(output, PISM_DOUBLE);
  m_snow_depth.define(output, PISM_DOUBLE);
}

void TemperatureIndexITM::write_model_state_impl(const PIO &output) const {
  m_log->message(2, "here, in write model state \n"); //FIXME
  SurfaceModel::write_model_state_impl(output);
  m_firn_depth.write(output);
  m_snow_depth.write(output);
}

namespace diagnostics {

/*! @brief Report surface melt, averaged over the reporting interval */
class SurfaceMelt : public DiagAverageRate<TemperatureIndexITM>
{
public:
  SurfaceMelt(const TemperatureIndexITM *m)
    : DiagAverageRate<TemperatureIndexITM>(m, "smelt", TOTAL_CHANGE) {

    m_vars = {SpatialVariableMetadata(m_sys, "smelt")};
    m_accumulator.metadata().set_string("units", "kg m-2");

    set_attrs("surface melt, averaged over the reporting interval", "",
              "kg m-2 s-1", "kg m-2 year-1", 0);
    m_vars[0].set_string("cell_methods", "time: mean");

    double fill_value = units::convert(m_sys, m_fill_value,
                                       m_vars[0].get_string("glaciological_units"),
                                       m_vars[0].get_string("units"));
    m_vars[0].set_double("_FillValue", fill_value);
  }

protected:
  const IceModelVec2S& model_input() {
    return model->melt();
  }
};

/*! @brief Report surface runoff, averaged over the reporting interval */
class SurfaceRunoff : public DiagAverageRate<TemperatureIndexITM>
{
public:
  SurfaceRunoff(const TemperatureIndexITM *m)
    : DiagAverageRate<TemperatureIndexITM>(m, "srunoff", TOTAL_CHANGE) {

    m_vars = {SpatialVariableMetadata(m_sys, "srunoff")};
    m_accumulator.metadata().set_string("units", "kg m-2");

    set_attrs("surface runoff, averaged over the reporting interval", "",
              "kg m-2 s-1", "kg m-2 year-1", 0);
    m_vars[0].set_string("cell_methods", "time: mean");

    double fill_value = units::convert(m_sys, m_fill_value,
                                       m_vars[0].get_string("glaciological_units"),
                                       m_vars[0].get_string("units"));
    m_vars[0].set_double("_FillValue", fill_value);
  }

protected:
  const IceModelVec2S& model_input() {
    return model->runoff();
  }
};

/*! @brief Report accumulation (precipitation minus rain), averaged over the reporting interval */
class Accumulation : public DiagAverageRate<TemperatureIndexITM>
{
public:
  Accumulation(const TemperatureIndexITM *m)
    : DiagAverageRate<TemperatureIndexITM>(m, "saccum", TOTAL_CHANGE) {

    m_vars = {SpatialVariableMetadata(m_sys, "saccum")};
    m_accumulator.metadata().set_string("units", "kg m-2");

    set_attrs("accumulation (precipitation minus rain), averaged over the reporting interval", "",
              "kg m-2 s-1", "kg m-2 year-1", 0);
    m_vars[0].set_string("cell_methods", "time: mean");

    double fill_value = units::convert(m_sys, m_fill_value,
                                       m_vars[0].get_string("glaciological_units"),
                                       m_vars[0].get_string("units"));
    m_vars[0].set_double("_FillValue", fill_value);
  }

protected:
  const IceModelVec2S& model_input() {
    return model->accumulation();
  }
};

} // end of namespace diagnostics

std::map<std::string, Diagnostic::Ptr> TemperatureIndexITM::diagnostics_impl() const {
  using namespace diagnostics;

  std::map<std::string, Diagnostic::Ptr> result = {
    {"saccum",      Diagnostic::Ptr(new Accumulation(this))},
    {"smelt",       Diagnostic::Ptr(new SurfaceMelt(this))},
    {"srunoff",     Diagnostic::Ptr(new SurfaceRunoff(this))},
    {"snow_depth",  Diagnostic::Ptr(new ITM_snow_depth(this))},
    {"firn_depth",  Diagnostic::wrap(m_firn_depth)},
  };

  result = pism::combine(result, SurfaceModel::diagnostics_impl());

  return result;
}

} // end of namespace surface
} // end of namespace pism
