Antarctica glacial-cycle spin-up example
=========

The scripts in this directory apply the PISM model to the Antarctic
Ice Sheet using climate forcing over the last two glacial cycles.

These scripts produce one selected simulation out of a parameter ensemble, which has been used as reference simulation in

* T. Albrecht, R. Winkelmann and A. Levermann, 
_Glacial-cycle simulations of the Antarctic Ice Sheet with the Parallel Ice Sheet Model (PISM) – Part 1: Boundary conditions and climatic forcing_, **The Cryosphere**, 14, 599–632, <https://doi.org/10.5194/tc-14-599-2020>, 2020


and 
* T. Albrecht, R. Winkelmann and A. Levermann, 
_Glacial-cycle simulations of the Antarctic Ice Sheet with the Parallel Ice Sheet Model (PISM) – Part 2: Parameter ensemble analysis_, **The Cryosphere**, 14, 633–656, <https://doi.org/10.5194/tc-14-633-2020>, 2020





Getting the input data
---------

The processed input data are provided on request (contact albrecht@pik-potsdam.de). They include

  * Topography: ice thickness, bed topography and ice surface elevation from [Bedmap2](https://www.bas.ac.uk/project/bedmap-2/).
  
  
  * Geothermal heatflux [Martos et al., 2017](https://doi.pangaea.de/10.1594/PANGAEA.882503).
  
  
  * Climate: precipitation and air temperature from [RACMOv.2.3p2](https://www.projects.science.uu.nl/iceclimate/models/racmo.php), forced with ERA-Interim, 1979-2016 (Wessem et al., 2017).
  
 
 * [Ocean](http://science.sciencemag.org/content/346/6214/1227): salinity and temperatures observations of Antarctic Continental Shelf Bottom Water (Schmidtko et al., 2014). 
 
 
 * Basins based on [Zwally drainage basins](http://imbie.org/imbie-2016/drainage-basins/) extended to Southern ocean
  
  
  * Temperature forcing: surface temperature reconstruction at [WDC](http://www.pnas.org/content/113/50/14249) (WAIS Divide Ice Core, contact C. Buizert), combined with [EDC](ftp://ftp.ncdc.noaa.gov/pub/data/paleo/icecore/antarctica/epica_domec/edc3deuttemp2007.txt) at EPICA Dome C before 67kyr
  
  
  * Sea-level forcing: eustatic sea-level reconstructions from [ICE-6G](http://www.atmosp.physics.utoronto.ca/~peltier/data.php) (contact D. Peltier), combined with [SPECMAP](https://doi.pangaea.de/10.1594/PANGAEA.734145) before 120kyr
  


Preprocessing of present-day input data for Antarctic simulations with PISM are described in https://github.com/pism/pism-ais 



Running the 16-km spin-up
---------

First set your environment paths to locate input data, output directory and pismr executable in:

    $ set_environment.sh

Next look at what the run script will do (i.e. a dry-run):

    $ PISM_DO=echo ./run-paleo.sh 16

Then actually do the run; this will take a number of processor-hours:

    $ ./run_paleo.sh 16

This is a 16km (initMIP) grid run with 16 processes which requires about one week to complete. The first stage essentially creates the initial guess of the enthalpy and bedrock displacement (< 1 yr), the second stage improves the enthalpy field (a "no mass continuity" run over 100kyr), and then the third stage optimizes the till friction angle for a best fit to present-day ice thickness (50kyr). Grounding line position and ice shelves are prescribed in this stage, while the 'full' physics are applied.

In a fourth stage the actual paleo-climate forced reference simulation is performed, with applied sea-level forcing, surface temperature and ocean temperature forcing (via PICO) over the last 205kyr. Grounding lines and calving fronts can freely evolve.
