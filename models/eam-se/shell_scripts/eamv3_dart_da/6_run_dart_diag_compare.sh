#!/bin/bash -el
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=esmd
#SBATCH  --time=2:00:00
#SBATCH  --partition=short
#SBATCH  --job-name=e3sm_dart_diag 
#SBATCH  --nodes=10
#SBATCH  --output=e3sm_dart_diag.%j 
#SBATCH  --exclusive
#SBATCH  --no-kill
#SBATCH  --requeue

echo == Start of e3sm dart diagnostic ==
date
echo ============================================

#source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh
source /qfs/people/zhan391/e3sm_dart_work/code/DART/models/eam-se/work/env_mach_specific.sh

#For cshell:
#limit stacksize unlimited
#limit datasize unlimited

#For bash 
#ulimit -s unlimited
#ulimit -d unlimited

#export SLURM_NNODES=20
#export SLURM_NTASKS=800

VERBOSE='-v'
MOVE='/usr/bin/mv'
COPY='/usr/bin/cp --preserve=timestamps'
LINK='/usr/bin/ln -fs'
LINKV=TRUE
LIST='/usr/bin/ls'
REMOVE='/usr/bin/rm'
LAUNCHCMD='srun -N 1 -n 1'

my_wkdir=${PWD}
scomp="eam"
cd ${my_wkdir}

source ./create_and_setup_case.sh

# Set paths
E3SM_ROOT=${my_e3sm_code}
DART_ROOT=${my_dart_code}
DART_MODEL=${my_dart_eam}
DART_SCPTDIR=${DART_ROOT}/models/${DART_MODEL}/shell_scripts
DART_WORKDIR=${DART_ROOT}/models/${DART_MODEL}/work
BASE_OBSDIR=${my_dart_obsdir}
BASE_PHIS=${my_e3sm_topo}
BASE_SEMAPS=${my_e3sm_semap}
BASE_CSGRID=${my_e3sm_csgrid}
CASE_ROOT=${my_modelcase}

# Run options
DART_ENSNUM=${my_ensnum}
DART_CASE=${my_casename}
DART_RUNDIR="${my_dart_runpath}/dart_en"`printf "%02d" ${DART_ENSNUM}`
DART_NTASKS=${my_dart_ntask}
DART_ON_PGRID=${my_dart_pgrid}

DART_RUNDIR1=${my_dart_rundir1}
DART_REFCASE1=${my_dart_refcase1}
DART_RUNDIR2=${my_dart_rundir2}
DART_REFCASE2=${my_dart_refcase2}

DATA_ASSIMILATION_CYCLES=${my_dart_cycle}
DATA_ASSIMILATION_WINDOW=${my_dart_window}
DATA_ASSIMILATION_ATM=TRUE

homme_map_file="SEMapping.nc"
cs_grid_file="SEMapping_cs_grid.nc"

CUR_YMD=${my_dartymds}
CUR_TOD=${my_darttods}
CUR_DATE=( `echo ${CUR_YMD}-${CUR_TOD} | sed -e "s#-# #g"` )
CUR_YEAR=`echo "${CUR_DATE[0]}" | bc`
CUR_MONTH=`echo "${CUR_DATE[1]}" | bc`
CUR_DAY=`echo "${CUR_DATE[2]}" | bc`
CUR_SECONDS=`echo "${CUR_DATE[3]}" | bc`
CUR_HOUR=`echo "${CUR_DATE[3]}" / 3600 | bc`
echo "valid time for eam forecast cycle is $CUR_YEAR $CUR_MONTH $CUR_DAY $CUR_SECONDS (seconds)"
echo "valid time for eam forecast cycle is $CUR_YEAR $CUR_MONTH $CUR_DAY $CUR_HOUR (hours)"

CURRENT_DADIR="${DART_RUNDIR}/dart_diagnostics"
if [ ! -d ${CURRENT_DADIR} ]; then
  mkdir -p ${CURRENT_DADIR}
fi
cd ${CURRENT_DADIR}

user_dart_nl() {
  if [ -e "${1}/diag_dart_input.nml" ]; then
    ${2} ${1}/diag_dart_input.nml input.nml  || exit 10
    sed -i "/#/d;/^\!/d;/^[ ]*\!/d" input.nml
    sed -i '1,1i\WARNING: Changes to this file will be ignored. \n Edit \$DART_WORKDIR/diag_dart_input.nml instead.\n\n\n'  input.nml
  else
    echo "ERROR ... DART required file ${1}/eam_dart_input.nml not found ... ERROR"
    exit 11
  fi

  xlist=`grep '^[ ]*vertical_localization_coord' input.nml`
  xlist=( `echo $xlist | sed -e "s#[=,']# #g"` )
  if [ "${xlist[1]}" == "SCALEHEIGHT" ]; then
     list1=`grep '^[ ]*vert_normalization_scale_height' input.nml `
     list1=( `echo $list1 | sed -e "s#[=,]##g"` )
     if [ "${list1[1]}" != "1.5" ]; then
        echo "WARNING!  input.nml is not using 1.5 for vert_normalization_scale_height."
        echo "          Use a different value only if you definitely want to. "
     fi
  else
     echo "WARNING!  input.nml is not using SCALEHEIGHT for vertical_localization_coord."
     echo "          SCALEHEIGHT is highly recommended for EAM"
  fi
  
  # If possible, use the round-robin approach to deal out the tasks.
  # This facilitates using multiple nodes for the simultaneous I/O operations.
  if [ -v ${3} ]; then
     if [ ${#3[@]} -gt 0 ]; then
        sed -i "s#layout.*#layout = 2#"  input.nml
        sed -i "s#tasks_per_node.*#tasks_per_node = ${3}#" input.nml
     fi
  fi
  
}

user_obs_sequence_nl() {
  local filename_seq_list="$1"
  local filename_out="$2"
  local qc_metadata="$3"
  local min_qc="$4"
  local max_qc="$5"
  local calendar="$6"
  local print_only="$7"
  local edit_copies="$8"
  local new_copy_index="$9"

  [[ -f input.nml ]] || { echo "input.nml not found"; return 1; }

  cat <<EOF | ex -s input.nml
g#^ *filename_seq *=#s#=.*#= '',#
g#^ *filename_seq_list *=#s#=.*#= ${filename_seq_list},#
g#^ *filename_out *=#s#=.*#= ${filename_out},#
g#^ *qc_metadata *=#s#=.*#= ${qc_metadata},#
g#^ *min_qc *=#s#=.*#= ${min_qc},#
g#^ *max_qc *=#s#=.*#= ${max_qc},#
g#^ *calendar *=#s#=.*#= '${calendar}',#
g#^ *print_only *=#s#=.*#= ${print_only},#
g#^ *edit_copies *=#s#=.*#= ${edit_copies},#
g#^ *new_copy_index *=#s#=.*#= ${new_copy_index},#
wq
EOF
}

user_obs_common_nl() {
  local num_to_compare_at_once="$1"
  local filename_seq_list="$2"
  local filename_out_suffix="$3"
  local print_every="$4"
  local dart_qc_threshold="$5"
  local calendar="$6"
  local print_only="$7"
  local eval_and_assim_can_match="$8"

  [[ -f input.nml ]] || { echo "input.nml not found"; return 1; }

  cat <<EOF | ex -s input.nml
g#^ *num_to_compare_at_once *=#s#=.*#= ${num_to_compare_at_once},#
g#^ *filename_seq *=#s#=.*#= '',#
g#^ *filename_seq_list *=#s#=.*#= ${filename_seq_list},#
g#^ *filename_out_suffix *=#s#=.*#= ${filename_out_suffix},#
g#^ *print_every *=#s#=.*#= ${print_every},#
g#^ *dart_qc_threshold *=#s#=.*#= ${dart_qc_threshold},#
g#^ *calendar *=#s#=.*#= '${calendar}',#
g#^ *print_only *=#s#=.*#= ${print_only},#
g#^ *eval_and_assim_can_match *=#s#=.*#= ${eval_and_assim_can_match},#
wq
EOF
}

cd ${CURRENT_DADIR}

############################################################
# run obs_diag (observation space diagnostics) 
############################################################
OBSCOMMON_DIR="${CURRENT_DADIR}/obs_common"
if [ ! -d "${OBSCOMMON_DIR}" ];then
   mkdir ${OBSCOMMON_DIR}
fi
cd ${OBSCOMMON_DIR}

user_dart_nl ${DART_WORKDIR} ${COPY} ${my_task_per_node}

#delete line that is deprecated by current version of DART
MYSTRING=`grep eam_template_filename input.nml`
MYSTRING=( `echo $MYSTRING | sed -e "s#[=,']# #g"` )
EAMINPUT=${MYSTRING[1]}

MYSTRING=`grep eam_phis_filename input.nml`
MYSTRING=( `echo $MYSTRING | sed -e "s#[=,']# #g"` )
EAM_PHIS=${MYSTRING[1]}
${LINK} ${BASE_PHIS} ${EAM_PHIS} || exit 100

#Now, Link the grid information files 
if [ ! -f ${BASE_SEMAPS} ] && [ ! -f ${BASE_CSGRID} ]; then
  echo "ERROR ... no mapping file ${homme_map_file}"
  echo "ERROR ... no gridinfo file ${cs_grid_file}"
  echo "ERROR ... must provide either of them"
  exit 91
else
  if [ -f ${BASE_SEMAPS} ]; then
    ${COPY} -r ${BASE_SEMAPS} ${homme_map_file} || exit 101
  fi
  if [ -f ${BASE_CSGRID} ]; then
    ${COPY} -r ${BASE_CSGRID} ${cs_grid_file} || exit 102
  fi
fi

#observation sequence file to process
obs_sequence_list1="${OBSCOMMON_DIR}/obs_seq1"
obs_sequence_list2="${OBSCOMMON_DIR}/obs_seq2"
obs_sequence_list3="${OBSCOMMON_DIR}/obs_seq3"

${LIST} ${DART_RUNDIR1}/*/${DART_REFCASE1}.dart.e.${scomp}_obs_seq_final.* > ${obs_sequence_list1}
${LIST} ${DART_RUNDIR2}/*/${DART_REFCASE2}.dart.e.${scomp}_obs_seq_final.* > ${obs_sequence_list2}
${LIST} ${DART_RUNDIR}/*/${DART_CASE}.dart.e.${scomp}_obs_seq_final.* > ${obs_sequence_list3}

[[ -f input.nml ]] || { echo "input.nml not found"; return 1; }

#first step: use obs_sequence_tool to convert the file 
${COPY} -f ${DART_WORKDIR}/obs_sequence_tool  ${OBSCOMMON_DIR} || exit 55

filename_seq_list="'obs_seq.tmp'"
filename_out="'obs_seq.processed'"
qc_metadata="''" #"'Data QC','DART quality control'"
print_every=10000
min_qc=-888888.0 #0
max_qc=-888888.0 #5
calendar='Gregorian'
print_only=".false."
edit_copies=".true."
new_copy_index="1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13"

user_obs_sequence_nl \
  "${filename_seq_list}" \
  "${filename_out}" \
  "${qc_metadata}" \
  "${min_qc}" \
  "${max_qc}" \
  "${calendar}" \
  "${print_only}" \
  "${edit_copies}" \
  "${new_copy_index}"

for obs_seq in ${obs_sequence_list1} ${obs_sequence_list2} ${obs_sequence_list3};do 
  cp ${obs_seq} obs_seq.tmp
  if [[ ${obs_seq} == ${obs_sequence_list1} ]];then 
     filename_out=${DART_REFCASE1}.obs_seq.processed   
  fi 
  if [[ ${obs_seq} == ${obs_sequence_list2} ]];then
     filename_out=${DART_REFCASE2}.obs_seq.processed
  fi
  if [[ ${obs_seq} == ${obs_sequence_list3} ]];then
     filename_out=${DART_CASE}.obs_seq.processed
  fi

  ${LAUNCHCMD} ${OBSCOMMON_DIR}/obs_sequence_tool || exit 149
done 

#second step: use obs_common_subset to generate diagnostics 
${COPY} -f ${DART_WORKDIR}/obs_common_subset  ${OBSCOMMON_DIR} || exit 56
num_to_compare_at_once=3
filename_seq_list="'obs_seq1','obs_seq2','obs_seq3'"
filename_out_suffix="'.common'"
print_every=10000
dart_qc_threshold=3
calendar='Gregorian'
print_only=".false."
eval_and_assim_can_match=".false."
user_obs_common_nl \
  "${num_to_compare_at_once}" \
  "${filename_seq_list}" \
  "${filename_out_suffix}" \
  "${print_every}" \
  "${dart_qc_threshold}" \
  "${calendar}" \
  "${print_only}" \
  "${eval_and_assim_can_match}"
${LAUNCHCMD} ${OBSCOMMON_DIR}/obs_common_subset || exit 150

echo ===== End of e3sm dart diagnostic =====
date
echo =====================================
