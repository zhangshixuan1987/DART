#!/bin/bash -el
#------------------------------------------------------------------------------
# SLURM Batch Directives
#------------------------------------------------------------------------------
#SBATCH --account=esmd
#SBATCH --time=2:00:00
#SBATCH --partition=short
#SBATCH --job-name=regrid_diag
#SBATCH --nodes=1
#SBATCH --output=regrid_diag.%j
#SBATCH --exclusive
#SBATCH --no-kill
#SBATCH --requeue

echo "== Start of DART diagnostic =="
date
echo "============================================"

# Load conda environment
source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh

# System utilities
MOVE='/usr/bin/mv'
COPY='/usr/bin/cp --preserve=timestamps'
LINK='/usr/bin/ln -fs'
REMOVE='/usr/bin/rm'
LIST='/usr/bin/ls'

# Environment setup (assumes these are exported externally or in create_and_setup_case.sh)
# E3SM_ROOT, DART_ROOT, my_modeldir, my_ensnum, my_casename, etc.

DART_MODEL=${my_dart_eam}
DART_WORKDIR=${DART_ROOT}/models/${DART_MODEL}/work
ARCHIVE_DIR="${my_modeldir}/archive"
MAP_FILE="/compyfs/zhan391/acme_init/map_file/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"

# Dates
ymds="2012-01-01"
ymde="2012-03-01"

read -r sy sm sd <<< "$(echo ${ymds} | tr '-' ' ')"
read -r ey em ed <<< "$(echo ${ymde} | tr '-' ' ')"
mday=(31 28 31 30 31 30 31 31 30 31 30 31)

# Variables
var_list=(TCO SCO FLUT PRECT PS PSL TS TBOT TUQ TVQ UBOT VBOT ZBOT U10 TAUX TAUY TREFHT QREFHT TMQ CAPE CIN TTQ
           U1000 U975 U950 U925 U900 U850 U800 U700 U600 U500 U400 U300 U200 U100 U050 U010
           V1000 V975 V950 V925 V900 V850 V800 V700 V600 V500 V400 V300 V200 V100 V050 V010
           Z1000 Z975 Z950 Z925 Z900 Z850 Z800 Z700 Z600 Z500 Z400 Z300 Z200 Z100 Z050 Z010
           T1000 T975 T950 T925 T900 T850 T800 T700 T600 T500 T400 T300 T200 T100 T050 T010
           Q1000 Q975 Q950 Q925 Q900 Q850 Q800 Q700 Q600 Q500 Q400 Q300 Q200 Q100 Q050 Q010
           OMEGA1000 OMEGA975 OMEGA950 OMEGA925 OMEGA900 OMEGA850 OMEGA800 OMEGA700 OMEGA600 OMEGA500
           OMEGA400 OMEGA300 OMEGA200 OMEGA100 TREFHTMN TREFHTMX PRECTMX)

hist="eam.h1"
freq="daily"
input="${ARCHIVE_DIR}/atm/hist"
outdir="${ARCHIVE_DIR}/post"

mkdir -p "${outdir}"

jobid=${SLURM_JOBID}

cd ${outdir}
workdir=$(mktemp -d tmp.${jobid}.XXXX)
cd ${workdir}

ENSTR=$(printf "EN%02d" ${i})
CASE_NAME="${my_casename}.${ENSTR}"

flist=$(mktemp ./input${ENSTR}.XXXXXX)

echo "=== Starting ensemble member ${CASE_NAME}.${hist}.${freq} ==="

ts_dest="${outdir}/atm/180x360_aave/ts/${freq}"
mkdir -p "${ts_dest}"

for year in $(seq "${sy}" "${ey}"); do
  yyyy=$(printf "%04d" "${year}")

  # Per-year scratch
  SCRATCH=$(mktemp -d scrach.${jobid}.XXXX)
  echo "Building daily temps for ${ENSTR} ${yyyy} in ${SCRATCH}"

  for month in $(seq 1 12); do
    # honor start/end month bounds on first/last year
    if [ "${year}" -eq "${sy}" ] && [ "${month}" -lt "${sm}" ]; then continue; fi
    if [ "${year}" -eq "${ey}" ] && [ "${month}" -gt "${em}" ]; then continue; fi

    mm=$(printf "%02d" "${month}")
    yymm="${yyyy}-${mm}"

    # days in month (handle leap-year Feb only)
    nday=${mday[$((month - 1))]}
    if (( month == 2 )) && (( (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0) )); then
      nday=29
    fi

    # For each day, average all hist files into a single daily file
    for dd in $(seq 1 "${nday}"); do
      dd2=$(printf "%02d" "${dd}")

      # Gather that day's input files directly from $input; no symlinks
      day_glob="${input}/${CASE_NAME}.${hist}.${yymm}-${dd2}"*.nc
      day_files=( ${day_glob} )

      if ((${#day_files[@]} == 0)); then
        # No files for this day -> skip
        continue
      fi

      # One daily mean file per day
      out_day="${SCRATCH}/ftmp_${yyyy}-${mm}-${dd2}.nc"
      ncra -O "${day_files[@]}" "${out_day}"
    done
  done

  # Collect all daily files for the year and sort (null-safe, portable)
  ffiles1=()
  while IFS= read -r -d '' f; do
    ffiles1+=("$f")
  done < <(find "$SCRATCH" -maxdepth 1 -type f -name "ftmp_${yyyy}-*.nc" -print0 | sort -z)

  if [ ${#ffiles1[@]} -eq 0 ]; then
    echo "Warning: No daily files built for ${ENSTR} ${yyyy} — skipping yearly aggregate."
    rm -rf -- "$SCRATCH"
    continue
  fi

  # Aggregate per variable, then regrid
  for var in "${var_list[@]}"; do
    outfile="${var}.${ENSTR}.${yyyy}.nc"
    tmp_year="./${var}_${ENSTR}_${yyyy}_year.nc"

    # Find the first daily file that *actually contains* this variable
    probe_file=""
    for f in "${ffiles1[@]}"; do
      if ncks -C -H -v "$var" "$f" >/dev/null 2>&1; then
        probe_file="$f"
        break
      fi
    done

    if [[ -z "$probe_file" ]]; then
      echo "Skipping ${var} — not present in any daily file for ${ENSTR} ${yyyy}."
      continue
    fi

    echo "Variable ${var} exists — concatenating and regridding..."
    # Concatenate along the record (time) dimension
    if ! ncrcat -O -v "$var" "${ffiles1[@]}" "$tmp_year"; then
      echo "ERROR: ncrcat failed for ${var} — skipping."
      continue
    fi

    # Regrid
    if ! ncremap -m "$MAP_FILE" -i "$tmp_year" -o "${ts_dest}/${outfile}"; then
      echo "ERROR: ncremap failed for ${var} — skipping."
      continue
    fi

    # Optional lightweight compression of the final product
    # ncks -O -4 -L 1 "${ts_dest}/${outfile}" "${ts_dest}/${outfile}"

  done

  # Clean per-year scratch
  rm -rf "${SCRATCH}"
done

cd ..
rm -rf "${workdir}"

echo "===== End of DART diagnostic ====="
date
echo "==================================="
