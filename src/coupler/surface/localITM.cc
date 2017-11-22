// Copyright (C) 2009, 2010, 2011, 2013, 2014, 2015, 2016, 2017 Ed Bueler and Constantine Khroulev and Andy Aschwanden
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

#include <cassert>
#include <ctime>  // for time(), used to initialize random number gen
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>       // M_PI
#include <cmath>                // for erfc() in CalovGreveIntegrand()
#include <algorithm>

#include "pism/util/pism_const.hh"
#include "pism/util/ConfigInterface.hh"
#include "localITM.hh"
#include "pism/util/IceGrid.hh"

namespace pism {
namespace surface {

LocalMassBalanceITM::Changes::Changes() {
  snow_depth = 0.0;
  melt       = 0.0;
  runoff     = 0.0;
  smb        = 0.0;
}

LocalMassBalanceITM::LocalMassBalanceITM(Config::ConstPtr myconfig, units::System::Ptr system)
  : m_config(myconfig), m_unit_system(system),
    m_seconds_per_day(86400) {
  // empty
}

LocalMassBalanceITM::~LocalMassBalanceITM() {
  // empty
}

std::string LocalMassBalanceITM::method() const {
  return m_method;
}

ITMMassBalance::ITMMassBalance(Config::ConstPtr config, units::System::Ptr system)
  : LocalMassBalanceITM(config, system) {
  precip_as_snow     = m_config->get_boolean("surface.pdd.interpret_precip_as_snow");
  Tmin               = m_config->get_double("surface.pdd.air_temp_all_precip_as_snow");
  Tmax               = m_config->get_double("surface.pdd.air_temp_all_precip_as_rain");
  refreeze_ice_melt  = m_config->get_boolean("surface.pdd.refreeze_ice_melt");

  m_method = "insolation temperature melt";
}


/*! \brief Compute the number of points for temperature and
    precipitation time-series.
 */
unsigned int ITMMassBalance::get_timeseries_length(double dt) {
  const unsigned int    NperYear = static_cast<unsigned int>(m_config->get_double("surface.pdd.max_evals_per_year"));
  const double dt_years = units::convert(m_unit_system, dt, "seconds", "years");

  return std::max(1U, static_cast<unsigned int>(ceil(NperYear * dt_years)));
}




/*  Determine albedo over ice cells (no ground or ocean cells respected here)
    if melt >= critical_melt, change albedo to melting snow
    TODO: put constants into config. 
*/
double ITMMassBalance::get_albedo(double melt, double snow_depth, double firn_depth){
  double albedo = 0. ;
  const double critical_melt = 100.;
  const double albedo_dry_snow = 0.8;
  const double albedo_melting_snow = 0.6;
  const double albedo_firn = 0.5;
  const double albedo_ice = 0.4; 

  if (snow_depth >= 0.) {
    if (melt < critical_melt) {
      albedo = albedo_dry_snow;
    }
    else {
      albedo = albedo_melting_snow;
    }
  }
  else if (firn_depth >= 0.) {
    albedo = albedo_firn;
  }
  else {
    albedo = albedo_ice;
  }
  return albedo;
}


//! 
/* compute melt by equation (16) from Robinson2010
 * @param dt_series length of the step for the time-series
 * @param T air temperature at time [k]
 * @param insolation at time [k]
 * @param surface_elevation
 * @param albedo which was should be figured by get_albedo (?)
 * @param[out] melt pointer to a pre-allocated array with N-1 elements
 * output in mm water equivalent
 */
void ITMMassBalance::calculate_ITM_melt(double dt_series,
                                         const double &insolation,
                                         const double &T,
                                         double &surface_elevation,
                                         double &albedo,  
                                         double &ITM_melt) {

  const double rho_w = 1.;    // mass density of water
  const double L_m = 1.;      // latent heat of ice melting
  double z = 1.;               // surface elevation
  double tau_a = 1. +  1. * z;  // transmissivity of the atmosphere, linear fit, plug in values
  const double itm_c = 1.;
  const double itm_lambda = 1.; 

  assert(dt_series > 0.0);

  const double h_days = dt_series / m_seconds_per_day;

  // TODO: check units!!!!
  ITM_melt = h_days / (rho_w * L_m) * (tau_a*(1. - albedo) * insolation + itm_c + itm_lambda * T); // hier muss irgendwo die Formel (16) from Robinson2010;

}


//! \brief Extract snow accumulation from mixed (snow and rain)
//! precipitation using the temperature time-series.
/** Uses the temperature time-series to determine whether the
 * precipitation is snow or rain. Rain is removed entirely from the
 * surface mass balance, and will not be included in the computed
 * runoff, which is meltwater runoff. There is an allowed linear
 * transition for Tmin below which all precipitation is interpreted as
 * snow, and Tmax above which all precipitation is rain (see, e.g.
 * [\ref Hock2005b]).
 *
 * Sets P[i] to the *solid* (snow) accumulation *rate*.
 *
 * @param[in,out] P precipitation rate (array of length N)
 * @param[in] T air temperature (array of length N)
 * @param[in] N array length
 */
void ITMMassBalance::get_snow_accumulation(const std::vector<double> &T,
                                           std::vector<double> &P) {

  assert(T.size() == P.size());
  const size_t N = T.size();

  // Following \ref Hock2005b we employ a linear transition from Tmin to Tmax
  for (unsigned int i = 0; i < N; i++) {
    // do not allow negative precipitation
    if (P[i] < 0.0) {
      P[i] = 0.0;
      continue;
    }

    if (precip_as_snow || T[i] <= Tmin) { // T <= Tmin, all precip is snow
      // no change
    } else if (T[i] < Tmax) { // linear transition from Tmin to Tmax
      P[i] *= (Tmax - T[i]) / (Tmax - Tmin);
    } else { // T >= Tmax, all precip is rain -- ignore it
      P[i] = 0.0;
    }
  }

}


//! \brief Compute the surface mass balance at a location from the amount of 
//! melted snow and the accumulation amount in a time interval.
/*!
 * This is a ITM scheme. The input parameter `melt_conversion_factor` is a conversion between
 * snow melting and ice melting per time period.
 *
 * `accumulation` has units "meter / second".
 *
 * - a fraction of the melted snow and ice refreezes, conceptualized
 *   as superimposed ice, and this is controlled by parameter \c
 *   refreeze_fraction
 *
 * - the excess of 'ITM_melt' is used to melt both the ice that came
 *   from refreeze and then any ice which is already present.
 *
 * The scheme here came from EISMINT-Greenland [\ref RitzEISMINT], but
 * is influenced by R. Hock (personal communication).
 */
ITMMassBalance::Changes ITMMassBalance::step(const double &melt_conversion_factor,
                                             const double &refreeze_fraction,
                                             double ITM_melt,
                                             double old_firn_depth,
                                             double old_snow_depth,
                                             double accumulation) {

  Changes result;

  double
    firn_depth      = old_firn_depth,
    snow_depth      = old_snow_depth + accumulation,
    max_snow_melted = ITM_melt,
    firn_melted     = 0.0,
    snow_melted     = 0.0,
    excess_melt     = 0.0;

  if (ITM_melt <= 0.0) {            // The "no melt" case.
    snow_melted = 0.0;
    firn_melted = 0.0,
    excess_melt = 0.0;
  } else if (max_snow_melted <= snow_depth) {
    // Some of the snow melted and some is left; in any case, all of
    // the energy available for melt was used up in melting snow.
    snow_melted = max_snow_melted;
    firn_melted = 0.0;
    excess_melt = 0.0;
  } else if (max_snow_melted <= firn_depth + snow_depth) {
    // All of the snow is melted but some firn is left; in any case, all of
    // the energy available for melt was used up in melting snow.
    snow_melted = snow_depth;
    firn_melted = max_snow_melted - snow_melted;
    excess_melt = 0.0;
  } else {
    // All (firn and snow_depth meters) of snow melted. Excess_pddsum is the
    // positive degree days available to melt ice.
    firn_melted = firn_depth;
    snow_melted = snow_depth;
    excess_melt = ITM_melt - (firn_melted + snow_melted) ; // 
  }

  double
    ice_melted              = excess_melt * melt_conversion_factor,
    melt                    = snow_melted + firn_melted + ice_melted,
    ice_created_by_refreeze = 0.0;

  if (refreeze_ice_melt) {
    ice_created_by_refreeze = melt * refreeze_fraction;
  } else {
    // Should this only be snow melted?
    ice_created_by_refreeze = (firn_melted + snow_melted) * refreeze_fraction;
  }

  firn_depth -= firn_melted;
  // FIXME: need to add snow that hasn't melted, is this correct?
  // firn_depth += (snow_depth - snow_melted);
  // Turn firn into ice at X times accumulation
  // firn_depth -= accumulation *  m_config->get_double("surface.pdd.firn_compaction_to_accumulation_ratio");

  if (firn_depth < 0.0) {
    firn_depth = 0.0;
  }
  snow_depth -= snow_melted;
  if (snow_depth < 0.0) {
    snow_depth = 0.0;
  }

  const double runoff = melt - ice_created_by_refreeze;

  result.firn_depth = firn_depth - old_firn_depth;
  result.snow_depth = snow_depth - old_snow_depth;
  result.melt       = melt;
  result.runoff     = runoff;
  result.smb        = accumulation - runoff;

  return result;
}




} // end of namespace surface
} // end of namespace pism
