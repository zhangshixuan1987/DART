#!/bin/csh
set workdir = `pwd`
# directory where DART prepbufr programs are located, relative to here.
set DART_exec_dir = /compyfs/zhan391/code/DART_E3SM
set ncepdat_exec_dir = ${DART_exec_dir}/observations/obs_converters/NCEP
set prebufr_exec_dir = ${ncepdat_exec_dir}/prep_bufr/work
set dartobs_exec_dir = ${ncepdat_exec_dir}/ascii_to_obs/work

setenv  My_BUFR_DATA  /compyfs/zhan391/acme_init/Observations/PREBUFR
set BUFR_idir    = ${My_BUFR_DATA}/prepqm
set BUFR_xdir    = ${My_BUFR_DATA}/prepout

if( ! -d $BUFR_idir ) then
  mkdir -rp $BUFR_idir
endif 

if( ! -d $BUFR_xdir ) then
  mkdir -rp $BUFR_xdir
endif

set YearList = (2009)
set MonList  = (1  2  3  4  5  6  7  8  9  10 11 12)
set DayList  = (31 28 31 30 31 30 31 31 30 31 30 31)

set nyr      = $#YearList

set nmon     = $#MonList

set iyS = 1
set iyE = 1 

set idS = 1
set idE = 28

set imS = 2
set imE = 2

set dp1 = 1
set mp1 = 1
set yp1 = 1

if ( $iyE > $nyr ) then
   set iyE = $nyr
endif

if ( $imE > $nmon ) then
   set imE = $nmon
endif

set iy = $iyS
while ( $iy <= $iyE )

 set year = ${YearList[$iy]}

 set im = $imS
 while ( $im <= $imE ) 

  set month      = ${MonList[$im]}
  set tot_days   = ${DayList[$im]}

  if ( $idE > $tot_days ) then
   set idE = $tot_days
  endif

  set ystr = `printf %04d $year`
  set mstr = `printf %02d $month`

  set BUFR_odir  = ${My_BUFR_DATA}/${ystr}${mstr}_12H_E3SM
  set BUFR_odir1 = ${My_BUFR_DATA}/${ystr}${mstr}_12H_E3SM-FV
  set BUFR_odir2 = ${My_BUFR_DATA}/${ystr}${mstr}_12H_E3SM-SE

  if ( ! -d $BUFR_odir ) then
   mkdir -rp $BUFR_odir
   ln -sf $BUFR_odir  $BUFR_odir1
   ln -sf $BUFR_odir  $BUFR_odir2
  endif


  set iday = $idS
  while ( $iday <= $idE )

   echo $year $month $iday
   #####prebufr to txt file
   ${prebufr_exec_dir}/prepbufr.csh $year $month $iday 

   @ iday++
  end
  
   #####txt file to dart file
   cp $dartobs_exec_dir/input.nml $BUFR_odir/input.nml
   sed -i "s/year       = 2007/year       =  $year/g"    $BUFR_odir/input.nml
   sed -i "s/month      = 1/month         = $month/g"    $BUFR_odir/input.nml
   sed -i "s/day        = 1/day           =   $idS/g"    $BUFR_odir/input.nml
   sed -i "s/tot_days   = 31/tot_days     =   $idE/g"    $BUFR_odir/input.nml
   sed -i "s/daily_file = .true./daily_file = .false./g" $BUFR_odir/input.nml
  #sed -i "s/obs_time   = .true./obs_time   = .false./g" $BUFR_odir/input.nml
   sed -i "s/ObsBase = '..\/..\/prep_bufr\/data\/temp_obs.'/ObsBase = '..\/prepout\/temp_obs.'/g" $BUFR_odir/input.nml

   cd $BUFR_odir
   $dartobs_exec_dir/create_real_obs
 
 #rename the files 
  set iday = $idS
  while ( $iday <= $idE )

   set dstr = `printf %02d $iday`

   @ dp1  = $iday + 1 
   
   if ($dp1 > $idE) then
     if ($month == 12 ) then
     @ yp1  = $year + 1
     set ypstr = `printf %02d $yp1`
     set mpstr = `printf %02d 1`
     set dpstr = `printf %02d 1`
     ln -sf obs_seq${ystr}${mstr}${dstr}24 obs_seq${ypstr}${mpstr}${dpstr}00
     else
     set ypstr = $ystr
     @ mp1  = $month + 1 
     set mpstr = `printf %02d $mp1`
     set dpstr = `printf %02d 1`
     ln -sf obs_seq${ystr}${mstr}${dstr}24 obs_seq${ypstr}${mpstr}${dpstr}00
     endif 
   else
     set dpstr = `printf %02d $dp1`
     ln -sf obs_seq${ystr}${mstr}${dstr}24 obs_seq${ystr}${mstr}${dpstr}00
   endif 
   @ iday++
  end
  
 @ im++
 end 

@ iy++
end 

exit 0
~          
