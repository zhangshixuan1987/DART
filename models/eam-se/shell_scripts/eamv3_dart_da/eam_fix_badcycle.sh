#!/bin/bash -el 
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=esmd
#SBATCH  --time=2:00:00
#SBATCH  --partition=short
#SBATCH  --job-name=e3sm_dart_ensda_cyc 
#SBATCH  --nodes=20
#SBATCH  --output=e3sm_dart_ensda_cyc.%j 
#SBATCH  --exclusive
#SBATCH  --no-kill
#SBATCH  --requeue

#For cshell:
#limit stacksize unlimited
#limit datasize unlimited

#For bash 
#ulimit -s unlimited
#ulimit -d unlimited

#export SLURM_NNODES=20
#export SLURM_NTASKS=800

echo == Start of eam_fix_badcycle.sh ==
date
echo ============================================

#source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh
source /qfs/people/zhan391/e3sm_dart_work/code/DART/models/eam-se/work/env_mach_specific.sh

my_wkdir=${PWD}

cd ${my_wkdir}
source ./create_and_setup_case.sh

DATA_ASSIMILATION_ATM=TRUE
DATA_ASSIMILATION_CYCLES=${my_dart_cycle}
DATA_ASSIMILATION_WINDOW=${my_dart_window}

CASE_ROOT=${my_modeldir}/EN01/case_scripts
RUN_ROOT=${my_modeldir}/EN01/run
ARCHIVE_DIR="${my_modeldir}/archive"

my_fixdate="2011-12-28"
my_fixtod="43200"
CUR_YMD=${my_fixdate}
CUR_TOD=${my_fixtod}

#determine the time for previous DA cycle 
CUR_DATE=( `echo ${CUR_YMD}-${CUR_TOD} | sed -e "s#-# #g"` )
CUR_YEAR=`echo "${CUR_DATE[0]}" | bc`
CUR_MONTH=`echo "${CUR_DATE[1]}" | bc`
CUR_DAY=`echo "${CUR_DATE[2]}" | bc`
CUR_HOUR=`echo "${CUR_DATE[3]}" / 3600 | bc`
CUR_SECONDS=`echo "${CUR_DATE[3]}" | bc`
echo "valid time for eam forecast cycle is $CUR_YEAR $CUR_MONTH $CUR_DAY $CUR_SECONDS (seconds)"

# First Step: Loop over members and run 6-hour ensembel forecast
for i in `seq 1 ${my_ensnum}`;do
  echo === Starting member ${i} ===
  ENSTR=EN`printf "%02d" ${i}`
  CASE_NAME=${my_casename}.${ENSTR}
  CASE_DIR=`echo ${CASE_ROOT} | sed "s/EN01/${ENSTR}/g"`
  RUN_DIR=`echo ${RUN_ROOT} | sed "s/EN01/${ENSTR}/g"`
  REF_DIR="${ARCHIVE_DIR}/rest/${CUR_YMD}-${CUR_TOD}"
  CASE_ARCHIVE_DIR=${ARCHIVE_DIR}/rest/${CUR_YMD}-${CUR_TOD}
  echo "Run Case: ${CASE_NAME}"
  echo "Run Directory: ${CASE_ARCHIVE_DIR}"
  if [ ! -d "${CASE_ARCHIVE_DIR}" ];then
    echo "archive data directory not fond ...."
    exit 99
  fi

  for scomp in "atm" "lnd" "rof" "ocn" "ice" "drv"; do
     echo === E3SM component ${scomp} ===
     cd ${CASE_ARCHIVE_DIR}
     if [[ ${scomp} == "atm" ]]; then
        smod="eam"
        CUR_DATE_EXT=${CUR_YMD}-${CUR_TOD}
        ATM_INITIAL_FILENAME=${CASE_NAME}.${smod}.i.${CUR_DATE_EXT}.nc
        echo "process ${smod} ic: ${ATM_INITIAL_FILENAME}"
        if [[ -f "${ATM_INITIAL_FILENAME}" && ! -f "${RUN_DIR}/${ATM_INITIAL_FILENAME}" ]];then 
          cp -rp ${ATM_INITIAL_FILENAME} ${RUN_DIR}/ || exit 01
        elif [[ ! -f "${ATM_INITIAL_FILENAME}" ]]; then 
          echo "ERROR: initial condition file not found: ${ATM_INITIAL_FILENAME}"
          exit
        fi
        ATM_REST_FILENAME="${CASE_NAME}.${smod}.r.${CUR_DATE_EXT}.nc"
        echo "${ATM_REST_FILENAME}"   >  ${RUN_DIR}/rpointer.atm
        if [[ -f "${ATM_REST_FILENAME}" && ! -f "${RUN_DIR}/${ATM_REST_FILENAME}" ]];then 
          cp -rp ${ATM_REST_FILENAME} ${RUN_DIR}/ || exit 02 
        elif [[ ! -f "${ATM_REST_FILENAME}" ]]; then
          echo "ERROR: initial condition file not found: ${ATM_REST_FILENAME}"
          exit
        fi
     elif [[ ${scomp} == "lnd" ]]; then
        smod="elm"
        CUR_DATE_EXT=${CUR_YMD}-${CUR_TOD}
        LND_INITIAL_FILENAME=${CASE_NAME}.${smod}.r.${CUR_DATE_EXT}.nc
        echo "process ${smod} ic: ${LND_INITIAL_FILENAME}"
        if [[ -f "${LND_INITIAL_FILENAME}" && ! -f "${RUN_DIR}/${LND_INITIAL_FILENAME}" ]];then
          cp -rp ${LND_INITIAL_FILENAME} ${RUN_DIR}/ || exit 03
        elif [[ ! -f "${LND_INITIAL_FILENAME}" ]]; then
          echo "ERROR: initial condition file not found: ${LND_INITIAL_FILENAME}"
          exit
        fi
        echo "./${LND_INITIAL_FILENAME}"   >  rpointer.lnd
     elif [[ ${scomp} == "rof" ]]; then
        smod="mosart"
        CUR_DATE_EXT=${CUR_YMD}-${CUR_TOD}
        ROF_INITIAL_FILENAME=${CASE_NAME}.${smod}.r.${CUR_DATE_EXT}.nc
        echo "process ${smod} ic: ${ROF_INITIAL_FILENAME}"
        if [[ -f "${ROF_INITIAL_FILENAME}" && ! -f "${RUN_DIR}/${ROF_INITIAL_FILENAME}" ]];then
          cp -rp ${ROF_INITIAL_FILENAME} ${RUN_DIR}/ || exit 04
        elif [[ ! -f "${ROF_INITIAL_FILENAME}" ]]; then
          echo "ERROR: initial condition file not found: ${ROF_INITIAL_FILENAME}"
          exit
        fi
        echo "./${ROF_INITIAL_FILENAME}"   >  rpointer.rof
     elif [[ ${scomp} == "ocn"  &&  ${my_runtype} == "Full-CPL" ]]; then
        smod="mpaso"
        CUR_DATE_EXT=${CUR_YMD}_${CUR_TOD}
        OCN_INITIAL_FILENAME=${CASE_NAME}.${smod}.rst.${CUR_DATE_EXT}.nc
        echo "process ${smod} ic: ${OCN_INITIAL_FILENAME}"
        if [[ -f "${OCN_INITIAL_FILENAME}" && ! -f "${RUN_DIR}/${OCN_INITIAL_FILENAME}" ]];then 
          cp -rp ${OCN_INITIAL_FILENAME} ${RUN_DIR}/ || exit 05
        elif [[ ! -f "${OCN_INITIAL_FILENAME}" ]]; then
          echo "ERROR: initial condition file not found: ${OCN_INITIAL_FILENAME}"
          exit
        fi
        echo "${CUR_DATE}_`printf "%02d" ${REF_HOUR}`:00:00"  > rpointer.ocn
     elif [[ ${scomp} == "ocn"  && ${my_runtype} == "AMIP" ]]; then
        smod="docn"
        CUR_DATE_EXT=${CUR_YMD}-${CUR_TOD}
        OCN_INITIAL_FILENAME1="${CASE_NAME}.${smod}.r.${CUR_DATE_EXT}.nc"
        OCN_INITIAL_FILENAME2="${CASE_NAME}.${smod}.rs1.${CUR_DATE_EXT}.bin"
        echo "process ${smod} ic: ${OCN_INITIAL_FILENAME1}"
        echo "process ${smod} ic: ${OCN_INITIAL_FILENAME2}"
        if [[ -f "${OCN_INITIAL_FILENAME2}" && ! -f "${RUN_DIR}/${OCN_INITIAL_FILENAME2}" ]];then
          cp -rp ${OCN_INITIAL_FILENAME2} ${RUN_DIR}/ || exit 07
        elif [[ ! -f "${OCN_INITIAL_FILENAME2}" ]]; then
          echo "ERROR: initial condition file not found: ${OCN_INITIAL_FILENAME2}"
          exit
        fi
        echo "${OCN_INITIAL_FILENAME1}"  >  rpointer.ocn
        echo "${OCN_INITIAL_FILENAME2}"  >> rpointer.ocn
     elif [[ ${scomp} == "ice" ]]; then
        smod="mpassi"
        CUR_DATE_EXT=${CUR_YMD}_${CUR_TOD}
        ICE_INITIAL_FILENAME=${CASE_NAME}.${smod}.rst.${CUR_DATE_EXT}.nc
        echo "process ${smod} ic: ${ICE_INITIAL_FILENAME}"
        if [[ -f "${ICE_INITIAL_FILENAME}" && ! -f "${RUN_DIR}/${ICE_INITIAL_FILENAME}" ]];then
          cp -rp ${ICE_INITIAL_FILENAME} ${RUN_DIR}/ || exit 08
        elif [[ ! -f "${ICE_INITIAL_FILENAME}" ]]; then
          echo "ERROR: initial condition file not found: ${ICE_INITIAL_FILENAME}"
          exit
        fi
        echo "${CUR_DATE}_`printf "%02d" ${REF_HOUR}`:00:00"  > rpointer.ice
     else
        smod="cpl"
        CUR_DATE_EXT=${CUR_YMD}-${CUR_TOD}
        CPL_INITIAL_FILENAME=${CASE_NAME}.${smod}.r.${CUR_DATE_EXT}.nc
        echo "process ${smod} ic: ${CPL_INITIAL_FILENAME}"
        if [[ -f "${CPL_INITIAL_FILENAME}" && ! -f "${RUN_DIR}/${CPL_INITIAL_FILENAME}" ]];then
          cp -rp ${CPL_INITIAL_FILENAME} ${RUN_DIR}/ || exit 09
        elif [[ ! -f "${CPL_INITIAL_FILENAME}" ]]; then
          echo "ERROR: initial condition file not found: ${CPL_INITIAL_FILENAME}"
          exit
        fi
        echo "${CPL_INITIAL_FILENAME}"  >  rpointer.drv
     fi
  done

  cd ${CASE_DIR}

  ./xmlchange run_exe="--kill-on-bad-exit=1 --job-name=${CASE_NAME} \${EXEROOT}/e3sm.exe "
  ./xmlchange RUN_TYPE="hybrid"
  ./xmlchange CONTINUE_RUN=FALSE
  ./xmlchange BFBFLAG=FALSE
  ./xmlchange RUN_STARTDATE="${CUR_YMD}"
  ./xmlchange START_TOD="${CUR_TOD}"
  ./xmlchange REST_OPTION="nhours"
  ./xmlchange REST_N="${DATA_ASSIMILATION_WINDOW}"
  ./xmlchange STOP_OPTION="nhours"
  ./xmlchange STOP_N="${DATA_ASSIMILATION_WINDOW}"
  ./xmlchange GET_REFCASE=FALSE
  ./xmlchange RUN_REFCASE="${CASE_NAME}"
  ./xmlchange RUN_REFDATE="${CUR_YMD}"
  ./xmlchange RUN_REFTOD="${CUR_TOD}"
  ./xmlchange RUN_REFDIR="${REF_DIR}"
  ./xmlchange DOUT_S=True  #FALSE
  ./xmlchange DOUT_S_ROOT="${ARCHIVE_DIR}"

  ./case.setup

  echo ============================

done

# That's all folks!
sleep 10

echo ===== End of eam_fix_badcycle.sh =====
date
echo =====================================
