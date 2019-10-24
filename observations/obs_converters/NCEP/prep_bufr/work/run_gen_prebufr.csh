#!/bin/csh
set workdir = `pwd`
# directory where DART prepbufr programs are located, relative to here.
set DART_exec_dir = /compyfs/zhan391/code/DART_E3SM
set ncepdat_exec_dir = ${DART_exec_dir}/observations/obs_converters/NCEP
set prebufr_exec_dir = ${ncepdat_exec_dir}/prep_bufr/work
set dartobs_exec_dir = ${ncepdat_exec_dir}/ascii_to_obs/work

setenv My_BUFR_DATA /compyfs/zhan391/acme_init/Observations/PREBUFR

set year       = 2009
set month      = 1
set day        = 1
set tot_days   = 10

set iday = 11
while ($iday <= $tot_days )

  @ nday = $day + $iday - 1 
  ${prebufr_exec_dir}/prepbufr.csh $year $month $nday 

  @ iday++
end

set yr = `printf %04d $year`
set mn = `printf %02d $month`
set dy = `printf %02d $day`

set BUFR_idir = ${My_BUFR_DATA}/prepqm
set BUFR_xdir = ${My_BUFR_DATA}/prepout
set BUFR_odir = ${My_BUFR_DATA}/${yr}${mn}_6H_E3SM

@ tot_days = $tot_days - 1

cp $dartobs_exec_dir/input.nml $BUFR_odir/input.nml

sed -i "s/year       = 2007/year       = $year/g" $BUFR_odir/input.nml
sed -i "s/month      = 1/month      = $month/g" $BUFR_odir/input.nml
sed -i "s/day        = 1/day        = $day/g" $BUFR_odir/input.nml
sed -i "s/tot_days   = 31/tot_days   = $tot_days/g" $BUFR_odir/input.nml
sed -i "s/daily_file = .true./daily_file = .false./g" $BUFR_odir/input.nml
#sed -i "s/obs_time   = .true./obs_time   = .false./g" $BUFR_odir/input.nml
sed -i "s/ObsBase = '..\/..\/prep_bufr\/data\/temp_obs.'/ObsBase = '..\/prepout\/temp_obs.'/g" $BUFR_odir/input.nml
cd $BUFR_odir

$dartobs_exec_dir/create_real_obs

exit 0
~          

