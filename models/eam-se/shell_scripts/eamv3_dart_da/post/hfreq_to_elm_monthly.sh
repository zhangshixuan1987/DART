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


# --- Config you likely already have ---
source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
: "${OMP_NUM_THREADS:=1}"; export OMP_NUM_THREADS

# Input window (inclusive)
ymds="2011-12-01"
ymde="2011-12-31"

ENSTR="EN01"
CASE_NAME="${my_casename}.${ENSTR}"

ARCHIVE_DIR="${my_modeldir}/archive"
ts_dest1="${ARCHIVE_DIR}/post/lnd/180x360_aave/ts/daily"

outdir="${ARCHIVE_DIR}/post/lnd/180x360_aave/monthly"
mkdir -p "${outdir}"

jobid=${SLURM_JOBID}
SCRATCH=$(mktemp -d "${outdir}/tmp.monthly.${jobid}.XXXX")
trap 'rm -rf "${SCRATCH}"' EXIT

echo "== Monthly mean builder =="
echo "Range: ${ymds} → ${ymde}   Ensemble: ${ENSTR}"
echo "Daily src:    ${ts_dest1}"
echo "Output:       ${outdir}"
date
echo "======================================"

monthly_mean_slice () {
  local f="$1" var="$2" start="$3" end="$4" out="$5"

  # 1) Ensure variable exists (ncks returns non-zero if var is missing)
  if ! ncks -m -v "$var" "$f" >/dev/null 2>&1; then
    echo "  - skip ${var}: not in $(basename "$f")"
    return 0
  fi

  # 2) Compute monthly mean over the requested window in one step.
  #    If the time range doesn't exist, ncra will fail (non-zero) and we skip.
  if ncra -O -d time,"$start","$end" -v "$var" "$f" "${out}.tmp.nc" >/dev/null 2>&1; then
    ncks -O -4 -L 1 "${out}.tmp.nc" "$out"
    rm -f "${out}.tmp.nc"
  else
    # No overlap or another ncra error → skip quietly
    rm -f "${out}.tmp.nc" 2>/dev/null || true
    return 0
  fi
}

# ---------- discover variables (any year) ----------
# Expect filenames like: VAR.ENxx.YYYY.nc
declare -A vars_seen
for f in "${ts_dest1}"/*.${ENSTR}.*.nc; do
  bn=$(basename "$f")               # VAR.ENxx.YYYY.nc
  var="${bn%%.*}"                   # up to first dot
  vars_seen["$var"]=1
done
vars=( "${!vars_seen[@]}" )
echo "Discovered variables: ${#vars[@]}"

# ---------- iterate month-by-month from ymds..ymde (GNU date) ----------
start_month=$(date -d "${ymds:0:7}-01" +%Y-%m-01)
end_month=$(date -d "${ymde:0:7}-01" +%Y-%m-01)

curr="$start_month"
while : ; do
  yyyy=$(date -d "$curr" +%Y)
  mm=$(date   -d "$curr" +%m)

  # month start/end (clamped to ymds/ymde if boundary months)
  mstart="${yyyy}-${mm}-01 00:00:0.0"
  mend_day=$(date -d "$curr +1 month -1 day" +%d)
  mend="${yyyy}-${mm}-${mend_day} 23:59:59.0"
  if [[ "${yyyy}-${mm}" == "${ymds:0:7}" ]]; then mstart="${ymds} 00:00:0.0"; fi
  if [[ "${yyyy}-${mm}" == "${ymde:0:7}" ]]; then mend="${ymde} 23:59:59.0"; fi

  echo ">> Month ${yyyy}-${mm}  window: [${mstart} , ${mend}]"

  for var in "${vars[@]}"; do
    src="${ts_dest1}/${var}.${ENSTR}.${yyyy}.nc"
    [[ -f "$src" ]] || continue

    out="${outdir}/${var}.${ENSTR}.${yyyy}-${mm}.nc"
    if [[ -s "$out" ]]; then
      echo "  ${var} ${yyyy}-${mm} exists → skip"
      continue
    fi

    monthly_mean_slice "$src" "$var" "$mstart" "$mend" "$out"
    [[ -s "$out" ]] && echo "  wrote ${var} ${yyyy}-${mm}" || echo "  ${var} ${yyyy}-${mm}: no samples"
  done

  [[ "$curr" == "$end_month" ]] && break
  curr=$(date -d "$curr +1 month" +%Y-%m-01)
done

echo "== Done =="
date

