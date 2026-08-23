#!/bin/bash -el
#------------------------------------------------------------------------------
# SLURM Batch Directives
#------------------------------------------------------------------------------
#SBATCH --account=esmd
#SBATCH --time=2:00:00
#SBATCH --partition=short
#SBATCH --job-name=regrid_diag
#SBATCH --nodes=1
#SBATCH --output=runtmp/logs/regrid_diag.%j
#SBATCH --exclusive
#SBATCH --no-kill
#SBATCH --requeue

set -Eeo pipefail
shopt -s nullglob
if [[ "${POSTPROCESS_DRIVER_ACTIVE:-FALSE}" != "TRUE" ]]; then
  echo "ERROR: internal post-processing worker; run 7_run_post_hist.sh" >&2
  exit 1
fi
echo "== Start of DART diagnostic =="
date
echo "============================================"

# System utilities
MOVE='/usr/bin/mv'
COPY='/usr/bin/cp --preserve=timestamps'
LINK='/usr/bin/ln -fs'
REMOVE='/usr/bin/rm'
LIST='/usr/bin/ls'

# Environment setup (assumes these are exported externally or in create_and_setup_case.sh)
# E3SM_ROOT, DART_ROOT, my_modeldir, my_ensnum, my_casename, etc.

DART_ROOT="${my_eam_dart_code}"
DART_MODEL=${my_eam_dart_model}
DART_WORKDIR=${DART_ROOT}/models/${DART_MODEL}/work
ARCHIVE_DIR="${my_modeldir}/${POST_ENSTR}/archive"
MAP_FILE="${my_eam_post_map_file:?my_eam_post_map_file is not set}"

# Dates
ymds="${POST_START_DATE}"
ymde="${POST_END_DATE}"

read -r sy sm sd <<< "$(echo ${ymds} | tr '-' ' ')"
read -r ey em ed <<< "$(echo ${ymde} | tr '-' ' ')"
mday=(31 28 31 30 31 30 31 31 30 31 30 31)

hist="eam.h0"
freq="clim"

input="${ARCHIVE_DIR}/atm/hist"
outdir="${ARCHIVE_DIR}/post"

mkdir -p "${outdir}"

jobid=${SLURM_JOB_ID:-$$}

cd ${outdir}
workdir=$(mktemp -d tmp.${jobid}.XXXX)
cd ${workdir}

ENSTR="${POST_ENSTR}"
CASE_NAME="${POST_CASE_NAME}"

flist=$(mktemp ./input${ENSTR}.XXXXXX)

echo "=== Starting ensemble member ${CASE_NAME}.${hist}.${freq} ==="

for year in $(seq "${sy}" "${ey}"); do
  for month in $(seq 1 12); do
    # Skip months outside the desired range if first/last year
    if [ "$year" -eq "$sy" ] && [ "$month" -lt "$sm" ]; then continue; fi
    if [ "$year" -eq "$ey" ] && [ "$month" -gt "$em" ]; then continue; fi  # note: > not >=

    yymm=$(printf "%04d-%02d" "${year}" "${month}")

    # Link files like ${CASE_NAME}.${hist}.YYYY-MM*.nc
    for ff in ${input}/${CASE_NAME}.${hist}.${yymm}.nc; do
      ln -sf $ff .
    done
  done
done

ls ${CASE_NAME}.${hist}.????-??.nc > ${flist}

# === Step 1: Regrid monthly h0 climatology files ===
clim_dest="${outdir}/atm/180x360_aave/${freq}"
mkdir -p "${clim_dest}"
while IFS= read -r ff; do
  outfile=$(basename "${ff}")
  ncremap -m "${MAP_FILE}" -i "${ff}" -o "${clim_dest}/${outfile}.tmp.${jobid}"
  ncdump -h "${clim_dest}/${outfile}.tmp.${jobid}" >/dev/null 2>&1
  ${MOVE} "${clim_dest}/${outfile}.tmp.${jobid}" "${clim_dest}/${outfile}"
done < ${flist}

cd ..
rm -rf "${workdir}"

echo "===== End of DART diagnostic ====="
date
echo "==================================="
