#!/bin/bash

# Copyright (C) 2019 PISM authors
# created by matthias.mengel@pik-potsdam.de and torsten.albrecht@pik-potsdam.de

###############################################################################
# Paleo spinup of Antarctic ice sheet model using various input data, 
# that can be obtained from albrecht@pik-potsdam.de on request.
# See README.md for more information
###############################################################################

SCRIPTNAME="#(run_paleo.sh)"

set -e  # exit on error

echo "$SCRIPTNAME   Run 16km grid paleo-climate spinup script using -pik options"
echo "$SCRIPTNAME   Run as './run_paleo.sh NN' for NN procs and 16km grid"

# get user and platform-specific variables like working_dir, pismcodedir,
# pism_exec and mpi command
source set_environment.sh

runname=`echo $PWD | awk -F/ '{print $NF}'`
#echo $runname
#codever=pism0.7_pik
thisdir=`echo $PWD`
outdir=$working_dir/$runname
PISM_EXEC=$pism_exec
#echo $outdir
mkdir -p $outdir/log


###############################################################################
# input data can be optained from albrecht@pik-potsdam.de or https://github.com/pism/pism-ais

export bedmap2file=$input_data_dir/bedmap2/bedmap2_16km.nc
export bheatflxfile=$input_data_dir/bheatflx_martos/bheatflx_martos_initmip16km.nc
export climatefile=$input_data_dir/racmo_wessem/racmo_wessem_initmip16km_mean1986_2005.nc
export upliftfile=$input_data_dir/uplift/urate_B_16km.nc
export gridfile=$input_data_dir/cdo_remapgrids/initmip16km.nc

export origfile=$input_data_dir/bedmap2_bheatflx_racmo_wessem_16km.nc
export atmfile=$origfile

if [ ! -f $oceanfile ]; then

  ncks -O -4 -v thk,topg,usurf,lon,lat $bedmap2file $origfile

  # correct for Lake Vostok
  ncap2 -O -s 'thk(161:181,261:291)=usurf(161:181,261:291)-topg(161:181,261:291)' $origfile $origfile

  #python $pismcode_dir/util/fill_missing.py -v bheatflx -i 0.07 -e 1.0 -f $bheatflxfile $bheatflxfile
  ncks -A -4 -v bheatflx $bheatflxfile $origfile

  #cdo -O -f nc4c timmean -selyear,1986/2005 $input_data_dir/racmo_wessem/racmo_wessem_initmip16km.nc $climatefile
  ncks -A -4 -v precipitation,air_temp $climatefile $origfile

  ncks -A -4 -v uplift $upliftfile $origfile
  ncrename -O -v uplift,dbdt $origfile

  ncks -A -4 -v x,y $gridfile $origfile

  #cdo -O -f nc4c merge

  #python $pismcode_dir/util/fill_missing.py -v v_magnitude -i 0.0 --eps 1.0 -f rignotvel_16km.nc -o rignotvel_16km_filled.nc
  #ncks -A -4 -v v_magnitude /p/projects/tumble/albrecht/pism_input/data/rignotvel/rignotvel_16km_filled.nc

fi

export oceanfile=$input_data_dir/schmidtko/schmidtko_initmip16km_means_11+12.nc
#ncks -A -4 -v salinity_ocean,theata_ocean,basins $oceanfile $origfile
#ncap2 -O -s 'where(basins==11) basins=12' $oceanfile &oceanfile
#ncap2 -O -s 'where(salinity_ocean==0.0) salinity_ocean=34.5' $oceanfile &oceanfile

export tforcefile=$input_data_dir/timeseries_edc-wdc_temp.nc #temperature
export slforcefile=$input_data_dir/imbrie06peltier15_sl.nc #sea level
export toforcefile=$input_data_dir/schmidtko14_edc-wdc_oceantemp_basins_response_fit_16km.nc #ocean temperature

export paramfile=$input_data_dir/pism_config_override
ncgen3 ${paramfile}.cdl -o ${paramfile}.nc

#export FTTMASK=${input_data_dir}/ftt_mask.nc
export OKMASK=${input_data_dir}/okill_mask.nc

# parallelization #############################################################
NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "run_paleo.sh 8" then NN = 8
  NN="$1"
fi
echo "$SCRIPTNAME              NN = $NN"
set -e  # exit on error


###### use only MPI if job is submitted #######################################
if [ -n "${PISM_ON_CLUSTER:+1}" ]; then  # check if env var is set
  echo "This run was submitted, use MPI"
  PISM_MPIDO=$pism_mpi_do
else
  echo "This is interactive, skip use of MPI"
  PISM_MPIDO=""
  NN=""
fi

#check if env var PISM_DO was set (i.e. PISM_DO=echo for a 'dry' run)
if [ -n "${PISM_DO:+1}" ] ; then  # check if env var is already set
  echo "$SCRIPTNAME         PISM_DO = $PISM_DO  (already set)"
else
  PISM_DO=""
fi
DO=$PISM_DO

echo "PISM_MPIDO = $PISM_MPIDO"
PISM_DO="$PISM_MPIDO $NN $PISM_EXEC"


# get new pism code if fetch is argument
if [ "$1" = "fetch" ]; then
  mkdir -p $outdir/bin/
  mkdir -p $outdir/log/
  rsync -aCv $pismcode_dir/bin/pismr $outdir/bin/
  cd $pismcode_dir
  echo ------ `date` --- $RUNNAME ------                  >> $thisdir/log/versionInfo
  echo "commit $(git log --pretty=oneline --max-count=1)" >> $thisdir/log/versionInfo
  echo "branch $( git branch | grep \*)"                  >> $thisdir/log/versionInfo
  cd $thisdir
fi


###############################################################################
source set_physics.sh

GRID=$SIXTEENKMGRID
SKIP=$SKIPFIFTEENKM
GRIDNAME=16km


###### output settings ########################################################
bootlength=200
nomasslength=200000
fitlength=50000
paleolength=205000



# #############################################################################
# bootstrap and diagnostic run to initialize bed deformation model
# #############################################################################
stage=obs
INNAME=$origfile
RESNAMEOBS=${outdir}/${stage}_${GRIDNAME}.nc
echo
echo "$SCRIPTNAME  bootstrapping on $GRIDNAME"
cmd="$PISM_MPIDO $NN $PISM_EXEC -i ${INNAME} -bootstrap $GRID \
        $stress_sia_opts $pre_calv_opts $bed_opts \
        $pre_ocean_opts $pre_atm_opts \
        $tech_opts $config_opts \
        -y 0.0001 -o $RESNAMEOBS"
echo $DO $cmd
$DO $cmd >> ${outdir}/log/out_${stage}_${GRIDNAME}.out
#exit # <-- uncomment to stop here


# #############################################################################
# bootstrap and SHORT smoothing run to 200 years
# #############################################################################
stage=boot
INNAME=$origfile
RESNAMEBOOT=${outdir}/${stage}_${GRIDNAME}.nc
RUNTIME=$bootlength
echo
echo "$SCRIPTNAME  bootstrapping on $GRIDNAME grid plus SIA run for $RUNTIME years"
cmd="$PISM_MPIDO $NN $PISM_EXEC -i ${INNAME} -bootstrap $GRID \
        $stress_sia_opts $pre_calv_opts $pre_bed_opts \
	$pre_ocean_opts $pre_atm_opts \
	$tech_opts $config_opts \
        -y $RUNTIME -o $RESNAMEBOOT" #-skip -skip_max $SKIP
#echo $DO $cmd
#$DO $cmd >> ${outdir}/log/out_${stage}_${GRIDNAME}.out
#exit # <-- uncomment to stop here


# #############################################################################
# run with -no_mass (no surface change) on coarse grid for 200kyr
# #############################################################################
stage=nomass
INNAME=$RESNAMEOBS
RESNAMENOMASS=${outdir}/${stage}_${GRIDNAME}.nc
TSNAME=${outdir}/ts_${stage}_${GRIDNAME}.nc
RUNTIME=$nomasslength
EXTRANAME=${outdir}/extra_${stage}_${GRIDNAME}.nc
expnomass="-extra_times 0:2500:$RUNTIME -extra_vars bmelt,tillwat,velsurf_mag,temppabase,diffusivity,hardav"
echo
echo "$SCRIPTNAME  -no_mass (no surface change) SIA for $RUNTIME years"
cmd="$PISM_MPIDO $NN $PISM_EXEC -i $INNAME -no_mass  \
     $stress_sia_opts $pre_calv_opts $pre_bed_opts \
     $pre_ocean_opts $pre_atm_opts \
     $tech_opts $config_opts \
     -ys 0 -y $RUNTIME \
     -extra_file $EXTRANAME $expnomass \
     -o $RESNAMENOMASS"
echo $DO $cmd
$DO $cmd >> ${outdir}/log/out_${stage}_${GRIDNAME}.out
#exit # <-- uncomment to stop here

#module load nco
#if [ ! -f $FTTMASK ]; then
  # create ftt_mask for EAIS in input file, could be done in preprocessing already
  #ncks -A -v thk,usurf,bedunc $origfile $FTTMASK
  #ncrename -O -v bedunc,ftt_mask $FTTMASK
  #ncatted -O -a standard_name,ftt_mask,d,, $FTTMASK $FTTMASK
  #ncatted -O -a long_name,ftt_mask,o,c,"mask where -force_to_thickness is applied" $FTTMASK $FTTMASK
  #ncatted -O -a units,ftt_mask,d,, $FTTMASK $FTTMASK
  #ncap2 -O -s 'where(ftt_mask>0.0) ftt_mask=0.0' $FTTMASK $FTTMASK
  #ncap2 -O -s 'where(usurf>2500.0) ftt_mask=1.0' $FTTMASK $FTTMASK
#fi

# #############################################################################
# run into quasi steady state with constant climate forcing after 50kyr,
# optimize for till friction angle
# #############################################################################
stage=fit
INNAME=$RESNAMENOMASS
RESNAMEFIT=${outdir}/${stage}_${GRIDNAME}.nc
RUNTIME=$fitlength

TSNAME=${outdir}/ts_${stage}_${GRIDNAME}.nc
EXTRANAME=${outdir}/extra_${stage}_${GRIDNAME}.nc
exvars="thk,bmelt,climatic_mass_balance,\
ice_surface_temp,shelfbmassflux,shelfbtemp,dHdt,mask,ocean_kill_mask,\
diff_usurf,usurf,tillphi,tillwat,diff_mask,taub_mag,tauc,\
velbar_mag,velsurf_mag,cell_area"
expoutfit="-extra_times 0:1000:$RUNTIME -extra_vars $exvars"
#tsvars="dt,ivol,ivolf,ivolg,imass,slvol,iareag,iareaf,iarea,surface_ice_flux,\
#grounded_basal_ice_flux,sub_shelf_ice_flux,nonneg_rule_flux,discharge_flux,\
#max_hor_vel,max_diffusivity"
tsoutfit="-ts_times 0:yearly:$RUNTIME " #-ts_vars=$tsvars"

echo
echo "$SCRIPTNAME  run into quasi steady state with constant climate forcing for $RUNTIME years"
cmd="$PISM_MPIDO $NN $PISM_EXEC -i $INNAME \
     $stress_opts $tillphi_param \
     $pre_calv_opts $pre_bed_opts \
     $pre_ocean_opts $fit_atm_opts \
     $tech_opts $config_opts \
     -ys 0 -y $RUNTIME \
     -ts_file $TSNAME $tsoutfit \
     -extra_file $EXTRANAME $expoutfit \
     -o $RESNAMEFIT" #-skip -skip_max $SKIP
echo $DO $cmd
$DO $cmd >> ${outdir}/log/out_${stage}_${GRIDNAME}.out


# add modern uplift as dbdt
#ncks -A -v viscous_bed_displacement,dbdt $RESNAMEOBS $RESNAMEFIT


export RESNAMEFITMOD=${outdir}/${stage}_${GRIDNAME}_mod.nc
source set_physics.sh


if [ ! -f $RESNAMEFITMOD ]; then
   # initiate from nomass run with optimized till friction angle and initial bedrok uplift rates
   cp $RESNAMENOMASS $RESNAMEFITMOD
   ncks -A -v tillphi $RESNAMEFIT $RESNAMEFITMOD
   ncks -A -v viscous_bed_displacement,dbdt $RESNAMEOBS $RESNAMEFITMOD   
fi

if [ ! -f $OKMASK ]; then

  # create ocean_kill_mask to constrain maximum ice extent to the edge of the continental shelf
  ncks -A -v bedunc,topg,thk $bedmap2file $OKMASK
  ncrename -O -v bedunc,ocean_kill_mask $OKMASK
  ncatted -O -a standard_name,ocean_kill_mask,d,, $OKMASK $OKMASK
  ncatted -O -a long_name,ocean_kill_mask,o,c,"mask where -ocean_kill cuts off ice" $OKMASK $OKMASK
  ncatted -O -a units,ocean_kill_mask,d,, $OKMASK $OKMASK
  ncap2 -O -s "where(topg<-2000) ocean_kill_mask=1" $OKMASK $OKMASK
  ncap2 -O -s "where(topg>=-2000) ocean_kill_mask=0" $OKMASK $OKMASK
  ncap2 -O -s "where(thk>0) ocean_kill_mask=0" $OKMASK $OKMASK
fi


# #############################################################################
# run paleo simulation with modified climate forcing, 2 glacial cycles (205kyr)
# #############################################################################
stage=paleo
INNAME=$RESNAMEFITMOD

RUNTIME=$paleolength
RESNAME=${outdir}/${stage}_${GRIDNAME}.nc

TSNAME=${outdir}/ts_${stage}_${GRIDNAME}.nc
SNAPNAME=${outdir}/snap_${stage}_${GRIDNAME}
EXTRANAME=${outdir}/extra_${stage}_${GRIDNAME}.nc

extratm=$((-paleolength)):1000:0
timestm=$((-paleolength)):yearly:0
snapstm=$((-paleolength)):5000:0
exvars="thk,topg,usurf,mask,velsurf_mag,dHdt,dbdt,diffusivity,\
bmelt,ice_surface_temp,climatic_mass_balance,shelfbmassflux,shelfbtemp,\
velbar,velbase,velsurf,wvelsurf,wvelbase,tauc,taub_mag,\
tempbase,tempsurf,sftgif,sftgrf,sftflf,hfgeoubed,\
rank,cell_area,ocean_kill_mask"
##surface_mass_balance_average,basal_mass_balance_average,\ ,discharge_flux,
#tsvars="vol,ivolf,ivolg,imass,slvol,iareag,iareaf,iarea,surface_ice_flux,\
#grounded_basal_ice_flux,sub_shelf_ice_flux,nonneg_rule_flux,discharge_flux,\
#max_hor_vel,max_diffusivity"
expout="-extra_file $EXTRANAME -extra_times $extratm -extra_vars $exvars"
snapout="-save_file $SNAPNAME -save_times $snapstm -save_split -save_size medium"
tsout="-ts_file $TSNAME -ts_times $timestm " #-ts_vars $tsvars"

regrid_opts="-bootstrap $GRID -regrid_file $INNAME -regrid_vars topg,thk,usurf,\
Href,tillwat,bmelt,enthalpy,litho_temp,temp,tillphi"
regrid_opts=""

echo
echo "$SCRIPTNAME  run paleo simulation with modified climate forcing starting $RUNTIME BP"
cmd="$PISM_MPIDO $NN $PISM_EXEC -i $INNAME $regrid_opts \
     $stress_opts $calv_opts $bed_opts \
     $ocean_opts $atm_opts \
     $tech_opts $config_opts \
    -ye 50 -ys $((-RUNTIME)) \
     $snapout $tsout $expout \
    -o $RESNAME -o_size medium" #-skip -skip_max $SKIP
echo $DO $cmd
$DO $cmd >> ${outdir}/log/out_${stage}_${GRIDNAME}.out



