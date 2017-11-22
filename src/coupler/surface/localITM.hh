// Copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2017 Ed Bueler and Constantine Khroulev
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

#ifndef __localITM_hh
#define __localITM_hh


namespace pism {
namespace surface {

//! \brief Base class for a model which computes surface mass flux rate (ice
//! thickness per time) from precipitation and temperature.
/*!
  This is a process model.  At each spatial location, it uses a 1D array, with a
  time dimension, for the temperature used in melting snow or ice.  At each spatial
  location it assumes the precipitation is time-independent.

  This process model does not know its location on the ice sheet, but
  simply computes the surface mass balance from three quantities:
  - the time interval \f$[t,t+\Delta t]\f$,
  - time series of values of surface temperature at \f$N\f$ equally-spaced
  times in the time interval
  - a scalar precipation rate which is taken to apply in the time interval.

  This model also uses degree day factors passed-in in DegreeDayFactors `ddf`,
  and the standard deviation `pddStdDev`.  The latter is the standard deviation of the
  modeled temperature away from the input temperature time series which contains
  the part of location-dependent temperature cycle on the time interval.

  @note 
  - Please avoid using `config.get...("...")` calls
  inside those methods of this class which are called inside loops over 
  spatial grids.  Doing otherwise increases computational costs.
  - This base class should be more general.  For instance, it could allow as
  input a time series for precipation rate.
*/
class LocalMassBalanceITM {
public:


  LocalMassBalanceITM(Config::ConstPtr config, units::System::Ptr system);
  virtual ~LocalMassBalanceITM();

  std::string method() const;

  virtual unsigned int get_timeseries_length(double dt) = 0;

  //! Count positive degree days (ITMs).  Returned value in units of K day.
  /*! Inputs T[0],...,T[N-1] are temperatures (K) at times t, t+dt_series, ..., t+(N-1)dt_series.
    Inputs `t`, `dt_series` are in seconds.  */
  /*virtual void get_ITMs(double dt_series,
                        //const std::vector<double> &albedo,
                        const std::vector<double> &T,
                        std::vector<double> &ITMs) = 0;*/

  /*! Remove rain from precipitation. */
  virtual void get_snow_accumulation(const std::vector<double> &T,
                                     std::vector<double> &precip_rate) = 0;

  class Changes {
  public:
    Changes();

    double firn_depth;
    double snow_depth;
    double melt;
    double runoff;
    double smb;
  };

  /** 
   * Take a step of the ITM model.
   *
   * @param[in] melt_conversion_factor , refreeze_fraction
   * @param[in] albedo
   * @param[in] melt
   * @param[in] old_firn_depth firn depth [ice equivalent meters]
   * @param[in] old_snow_depth snow depth [ice equivalent meters]
   * @param[in] accumulation total solid (snow) accumulation during the time-step [ice equivalent meters]
   */
  virtual Changes step(const double &melt_conversion_factor,
                       const double &refreeze_fraction,
                       double ITM_melt,
                       double old_firn_depth,
                       double old_snow_depth,
                       double accumulation) = 0;

protected:
  std::string m_method;

  const Config::ConstPtr m_config;
  const units::System::Ptr m_unit_system;
  const double m_seconds_per_day;
};


//! A ITM implementation which computes the local mass balance based on an expectation integral.
/*!
  The expected number of positive degree days is computed by an integral in \ref CalovGreve05.

  Uses degree day factors which are location-independent.
*/
class ITMMassBalance : public LocalMassBalanceITM {

public:
  ITMMassBalance(Config::ConstPtr config, units::System::Ptr system);
  virtual ~ITMMassBalance() {}

  virtual unsigned int get_timeseries_length(double dt);

  virtual double get_albedo();
  
  virtual void calculate_ITM_melt(double dt_series,
                        const std::vector<double> &insolation,
                        const std::vector<double> &T,
                        std::vector<double> &surface_elevation,
                        std::vector<double> &albedo,  
                        std::vector<double> &ITM_melt);

  virtual void get_snow_accumulation(const std::vector<double> &T,
                                     std::vector<double> &precip_rate);

  virtual Changes step(const double &melt_conversion_factor,
                       const double &refreeze_fraction,
                       double ITM_melt,
                       double firn_depth,
                       double snow_depth,
                       double accumulation);

protected:

  bool precip_as_snow,          //!< interpret all the precipitation as snow (no rain)
    refreeze_ice_melt;          //!< refreeze melted ice
  double Tmin,             //!< the temperature below which all precipitation is snow
    Tmax;             //!< the temperature above which all precipitation is rain
};


	} // end of namespace surface
} // end of namespace pism

#endif