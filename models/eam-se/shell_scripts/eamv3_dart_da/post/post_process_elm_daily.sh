#!/bin/bash -el
#------------------------------------------------------------------------------
# SLURM Batch Directives
#------------------------------------------------------------------------------
#SBATCH --account=esmd
#SBATCH --time=02:00:00
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
MAP_FILE="/compyfs/zhan391/acme_init/map_file/map_r05_to_cmip6_180x360_aave.20200901.nc"

# Dates
ymds="2012-01-01"
ymde="2012-03-01"

read -r sy sm sd <<< "$(echo ${ymds} | tr '-' ' ')"
read -r ey em ed <<< "$(echo ${ymde} | tr '-' ' ')"
mday=(31 28 31 30 31 30 31 31 30 31 30 31)

# Variables
var_list=(
  H2OSNO FSNO QRUNOFF QSNOMELT FSNO_EFF SNORDSL SNOW FSA FSDS FSR FLDS
  FIRE FIRA SOILWATER_10CM SOILLIQ SOILICE QSOIL U10 U10WITHGUSTS TWS
  TSOI_10CM TSA TLAI THBOT TSOI TSOI_ICE SOILLIQ_ICE SOILICE_ICE
  TAUX TAUY FSH HC HCSOI EFLX_LH_TOT SNOW_DEPTH RH2M RAIN QVEGE QVEGT
  QBOT Q2M H2OSOI H2OSFC ZWT ZBOT TBOT TG THBOT PBOT
)

hist="elm.h1"
freq="daily"
input="${ARCHIVE_DIR}/lnd/hist"
outdir="${ARCHIVE_DIR}/post"

mkdir -p "${outdir}"

jobid=${SLURM_JOBID}

cd ${outdir}
workdir=$(mktemp -d tmp.${jobid}.XXXX)
cd ${workdir}

ENSTR=$(printf "EN%02d" ${i})
CASE_NAME="${my_casename}.${ENSTR}"

echo "=== Starting ensemble member ${CASE_NAME}.${hist}.${freq} ==="

# === Step 2: Regrid var_list from hist (daily) ===
ts_dest="${outdir}/lnd/180x360_aave/ts/${freq}"
mkdir -p "${ts_dest}"

for year in $(seq "${sy}" "${ey}"); do
  yyyy=$(printf "%04d" "${year}")

  # Per-year scratch
  SCRATCH=$(mktemp -d "./tmp_hist_${ENSTR}_${yyyy}.XXXXXX")
  echo "Building daily temps for ${ENSTR} ${yyyy} in ${SCRATCH}"

  # Build daily means for months within bounds
  for month in $(seq 1 12); do
    # honor start/end bounds in first/last year
    if [ "${year}" -eq "${sy}" ] && [ "${month}" -lt "${sm}" ]; then continue; fi
    if [ "${year}" -eq "${ey}" ] && [ "${month}" -gt "${em}" ]; then continue; fi

    mm=$(printf "%02d" "${month}")
    yymm="${yyyy}-${mm}"

    # days in month (leap-year Feb only)
    nday=${mday[$((month - 1))]}
    if (( month == 2 )) && (( (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0) )); then
      nday=29
    fi

    for dd in $(seq 1 "${nday}"); do
      dd2=$(printf "%02d" "${dd}")
      day_glob="${input}/${CASE_NAME}.${hist}.${yymm}-${dd2}"*.nc
      day_files=( ${day_glob} )

      # no files for this day → skip
      ((${#day_files[@]}==0)) && continue

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

  # Probe once for landfrac presence (in the first daily file)
  have_landfrac=0
  if ncks -m "${ffiles1[0]}" | grep -qE "^landfrac\b"; then
    have_landfrac=1
  fi

  for var in "${var_list[@]}"; do
    outfile="${var}.${ENSTR}.${yyyy}.nc"
    tmp_year="./${var}_${ENSTR}_${yyyy}_year.nc"

    # Build list of variables we need in the concat
    need_vars=("$var")
    if (( have_landfrac )); then
      need_vars+=("landfrac")
    fi

    # Filter files that contain *all* needed variables
    files_ok=()
    for f in "${ffiles1[@]}"; do
      ok=1
      for v in "${need_vars[@]}"; do
        ncks -C -H -v "$v" "$f" >/dev/null 2>&1 || { ok=0; break; }
      done
      (( ok )) && files_ok+=("$f")
    done

    if ((${#files_ok[@]} == 0)); then
      echo "Skipping ${var} — none of the daily files contain all of: ${need_vars[*]}"
      continue
    fi

    # Compose comma-separated list for -v without spaces
    sel_vars=$(IFS=,; echo "${need_vars[*]}")

    # Yearly concat using only files that have all requested variables
    ncrcat -O -v "${sel_vars}" "${files_ok[@]}" "${tmp_year}"

    if (( have_landfrac )); then
      # ELM-aware regrid when landfrac exists
      ncremap -P elm -m "${MAP_FILE}" -i "${tmp_year}" -o "${ts_dest}/${outfile}"
    else
      echo "landfrac not found for ${ENSTR} ${yyyy}; using plain ncremap for ${var}."
      ncremap -m "${MAP_FILE}" -i "${tmp_year}" -o "${ts_dest}/${outfile}"
    fi

    # Optional compression:
    # ncks -O -4 -L 1 "${ts_dest}/${outfile}" "${ts_dest}/${outfile}"
  done

  rm -rf "${SCRATCH}"
done

wait 

cd ..
rm -rf "${workdir}"

echo "===== End of DART diagnostic ====="
date
echo "==================================="
