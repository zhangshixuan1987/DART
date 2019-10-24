#!/bin/csh 

set year       = 2009
set month      = 1
set day        = 1
set tot_days   = 4

set yr = `printf %04d $year`
set mn = `printf %02d $month`
set dy = `printf %02d $day`

set BUFR_dir  = /compyfs/zhan391/acme_init/Observations/PREBUFR
set BUFR_idir = ${BUFR_dir}/prepqm
set BUFR_xdir = ${BUFR_dir}/prepout
set BUFR_odir = ${BUFR_dir}/${yr}${mn}_6H_E3SM

# directory where DART prepbufr programs are located, relative to here.
set DART_exec_dir = /compyfs/zhan391/code/DART_E3SM/observations/obs_converters/NCEP/ascii_to_obs/work

cp $DART_exec_dir/input.nml $BUFR_odir/input.nml

sed -i "s/year       = 2007/year       = $year/g" $BUFR_odir/input.nml
sed -i "s/month      = 1/month      = $month/g" $BUFR_odir/input.nml
sed -i "s/day        = 1/day        = $day/g" $BUFR_odir/input.nml
sed -i "s/tot_days   = 31/tot_days   = $tot_days/g" $BUFR_odir/input.nml
sed -i "s/daily_file = .true./daily_file = .false./g" $BUFR_odir/input.nml
sed -i "s/ObsBase = '..\/..\/prep_bufr\/data\/temp_obs.'/ObsBase = '..\/prepout\/temp_obs.'/g" $BUFR_odir/input.nml
cd $BUFR_odir

$DART_exec_dir/create_real_obs

exit 0
