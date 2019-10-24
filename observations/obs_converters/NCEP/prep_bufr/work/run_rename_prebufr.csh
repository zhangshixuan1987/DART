#!/bin/csh
set workdir = `pwd`
# directory where DART prepbufr programs are located, relative to here.
set DART_exec_dir = /compyfs/zhan391/code/DART_E3SM
set ncepdat_exec_dir = ${DART_exec_dir}/observations/obs_converters/NCEP
set prebufr_exec_dir = ${ncepdat_exec_dir}/prep_bufr/work
set dartobs_exec_dir = ${ncepdat_exec_dir}/ascii_to_obs/work
setenv My_BUFR_DATA /compyfs/zhan391/acme_init/Observations/PREBUFR

set yrst       = 2009
set monst      = 1
set yred       = 2009
set moned      = 1

set iyr = $yrst
while ($iyr <= $yred )

 set imon = $monst
 while ($imon <= $moned )

 set yr = `printf %04d $iyr`
 set mn = `printf %02d $imon`

 set BUFR_idir = ${My_BUFR_DATA}/prepqm
 set BUFR_xdir = ${My_BUFR_DATA}/prepout
 set BUFR_odir = ${My_BUFR_DATA}/${yr}${mn}_6H_E3SM

 cd $BUFR_odir
  
 foreach file (obs_seq*)
  
  set year  = `echo $file |cut -c8-11`
  set month = `echo $file |cut -c12-13`
  set day   = `echo $file |cut -c14-15`
  set hour  = `echo $file |cut -c16-17`

  if ( $hour == 24 ) then
   echo $year $month $day $hour
   echo $file
    @ dp1 = $day + 1 
    set hr1 = 0 
    set day1  = `printf %02d $dp1`
    set hour1 = `printf %02d $hr1`
    rm -rvf obs_seq$year$month$day1$hour1
    ln -sf $file obs_seq$year$month$day1$hour1
  endif
   
 end  
 

 @ imon++
 end 

@ iyr++
end



exit 0
~          

