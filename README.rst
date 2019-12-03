PISM, a Parallel Ice Sheet Model
================================

The Parallel Ice Sheet Model is an open source, parallel, high-resolution ice sheet model:

- hierarchy of available stress balances
- marine ice sheet physics, dynamic calving fronts
- polythermal, enthalpy-based conservation of energy scheme
- extensible coupling to atmospheric and ocean models
- verification and validation tools
- `documentation <pism-docs_>`_ for users and developers
- uses MPI_ and PETSc_ for parallel simulations
- reads and writes `CF-compliant <cf_>`_  NetCDF_ files



Latest PIK code improvements
----------------------------

The pism_pik_1.0 branch is based on the stable version 361b6d from November 8th, 2017 and adds methods that are used in paleo simulations of the Antarctic Ice Sheet

- sub-ice shelf melting calculations with PICO, with added diagnostics
- iterative optimization of till-friction angle to present-day grounded surface elevation
- air temperature parameterization based on a multiregression fit to ERA-Interim reanalysis data kill ocean
- non-linear lapse rate scaling for precipitation
- reads ocean_kill_mask from file to constrain maximum ice extent
- fix of the calculation of sea-level relevant volume in the output timeseries

- fix of linear interpolation scheme for changing sea-level height
- fix of the bed deformation model to exclude changes in ice shelf loads and include changes in the ocean load 
- improvement of damage formation based on Borstad et al., 2016

.. You find in the examples/paleo-antarctica folder a working example of a paleo spin-up using all added functionality.

Code release v1.0-paleo-ensemble was used for ``Glacial cycles simulation of the Antarctic Ice Sheet with PISM`` in Albrecht, Winkelmann and Levermann, The Cryosphere (2019). If you make use of this code, please cite the respective paper.

This code release builds on a `previous release <pism-rebound_>`_ for a publication in Nature by Kingslake, Scherer, Albrecht et al., "Extensive retreat and re-advance of the West Antarctic Ice Sheet during the Holocene." Nature 558, no. 7710 (2018): 430.


PISM is jointly developed at the `University of Alaska, Fairbanks (UAF) <uaf_>`_ and the
`Potsdam Institute for Climate Impact Research (PIK) <pik_>`_. UAF developers are based in
the `Glaciers Group <glaciers_>`_ at the `Geophysical Institute <gi_>`_.

Please see ``ACKNOWLEDGE.rst`` and ``doc/funding.csv`` for a list of grants supporting
PISM development.

Homepage
--------

    http://www.pism-docs.org/

Download and Install
--------------------

See the `Installing PISM <pism-installation_>`_ on ``pism-docs.org``.

Contributing
------------

Want to contribute? Great! See `Committing to PISM <pism-contribute_>`_.

.. URLs

.. _uaf: http://www.uaf.edu/
.. _pik: http://www.pik-potsdam.de/
.. _pism-docs: http://www.pism-docs.org/
.. _pism-stable: http://www.pism-docs.org/wiki/doku.php?id=stable_version
.. _pism-contribute: http://www.pism-docs.org/wiki/doku.php?id=committing
.. _pism-installation: http://pism-docs.org/sphinx/installation/
.. _mpi: http://www.mcs.anl.gov/research/projects/mpi/
.. _petsc: http://www.mcs.anl.gov/petsc/
.. _cf: http://cf-pcmdi.llnl.gov/
.. _netcdf: http://www.unidata.ucar.edu/software/netcdf/
.. _glaciers: http://www.gi.alaska.edu/snowice/glaciers/
.. _gi: http://www.gi.alaska.edu
.. _NASA-MAP: http://map.nasa.gov/
.. _NASA-Cryosphere: http://ice.nasa.gov/
.. _NSF-Polar: https://nsf.gov/geo/plr/about.jsp
.. _pism-rebound: https://github.com/pism/pism/releases/tag/pik-holocene-gl-rebound

..
   Local Variables:
   fill-column: 90
   End:
