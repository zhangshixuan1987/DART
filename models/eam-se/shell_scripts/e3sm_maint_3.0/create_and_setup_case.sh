# Shared configuration for all numbered workflow stages.
# This file is sourced; it should define settings without running workflow work.

################################################################################
# --- Workflow and runtime paths ---------------------------------------------
# The workflow root is resolved from this configuration file at source time.
_my_config_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export my_workflow_root="${_my_config_dir}"
unset _my_config_dir
export my_workflow_lib="${my_workflow_root}/workflow_lib"
export my_runtime_dir="${my_workflow_root}/runtmp"
export my_log_dir="${my_runtime_dir}/logs"
export my_status_dir="${my_runtime_dir}/status"
export my_lock_dir="${my_runtime_dir}/locks"
export my_handoff_dir="${my_runtime_dir}/handoff"
export my_run_script_dir="${my_runtime_dir}/run_scripts"

################################################################################
# Step 4 --nodes must match my_job_nnodes; other stages use their own allocations.
################################################################################
export my_machine=compy
export my_project="esmd"
export my_jobqueue="slurm"
export my_walltime="02:00:00"
export my_task_per_node=40
export my_job_nnodes=160
export my_layout="custom-4_1x6_nhours"

################################################################################
# --- Machine and Slurm defaults ---------------------------------------------
# Explicit analysis environment used by Steps 2 and 3 for NCO and related tools.
################################################################################
export my_conda_setup_file="/people/zhan391/miniforge3/etc/profile.d/conda.sh"
export my_analysis_conda_env="e3sm_analysis"
export my_dart_env_file="${my_workflow_lib}/env/env_${my_machine}_specific.sh"
export my_eam_filter_nml="${my_workflow_lib}/namelists/eam/filter.nml"
export my_eam_perturb_nml="${my_workflow_lib}/namelists/eam/perturb.nml"
export my_eam_diag_nml="${my_workflow_lib}/namelists/eam/diagnostics.nml"
export my_elm_filter_nml="${my_workflow_lib}/namelists/elm/filter.nml"

################################################################################
# --- Ensemble and cycle execution -------------------------------------------
# Initial model time, ensemble size, concurrency, retries, and cycle batching.
################################################################################
export my_ensnum=40
export my_casedate="2011-12-01"
export my_casetod="00000"
export my_forecast_ready_for_da="FALSE"
export my_nodes_per_member=4
export my_retry_forecast_timeout_sec=12600
export my_wait_poll_interval_sec=20
export my_skip_completed_members=FALSE
export my_max_parallel_setup=10
export my_max_parallel_handoff=10
# Cycles attempted in one Step 4 allocation; use 1 for one cycle per job.
export my_cycles_per_job=1
# Start another cycle only when this runtime plus the shutdown margin remains.
export my_min_cycle_time_sec=12600
export my_cycle_shutdown_margin_sec=900

################################################################################
# --- E3SM experiment identity and output paths -------------------------------
# my_runtype must be Full-CPL or AMIP; keep compset and reference files consistent.
################################################################################
export my_e3sm_code="/qfs/people/zhan391/e3sm_dart_work/code/E3SMv3"
export my_runtype="Full-CPL"
export my_compset="WCYCL20TR"
export my_resolution="ne30pg2_r05_IcoswISC30E3r5"
export my_runpath="/compyfs/zhan391/v3_dart_cda_scratch"
export my_casename="DARTEN${my_ensnum}_${my_compset}_${my_resolution}_${my_machine}"

################################################################################
# Model experiment directory and shared executable. Step 1 builds E3SM once
# and reuses that executable for the ensemble members.
################################################################################
export my_modeldir="${my_runpath}/${my_casename}"
export my_modelcase="${my_modeldir}/case_scripts"
export my_modelexe="${my_modeldir}/build/e3sm.exe"
export my_eam_topography_file="/compyfs/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4pg2_x6t-SGH.c20210614.nc"
export my_eam_se_mapping_file="/qfs/people/zhan391/e3sm_dart_work/code/HOMME/SEMapping.nc"
export my_eam_cs_grid_file="/qfs/people/zhan391/e3sm_dart_work/code/DART/models/eam-se/work/SEMapping_cs_grid_NE30.nc"
export my_eam_post_map_file="/compyfs/zhan391/acme_init/map_file/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"
export my_elm_post_map_file="/compyfs/zhan391/acme_init/map_file/map_r05_to_cmip6_180x360_aave.20200901.nc"
export my_amip_sst_data="/compyfs/zhan391/acme_init/SST_forcing/sst_ice_NOAA_AVHRR_E3SM_1x1_c20231225.nc"
export my_amip_sst_grid="/compyfs/zhan391/acme_init/SST_forcing/domain.ocn.1x1.111007.nc"

################################################################################
# --- Reference initial conditions -------------------------------------------
# Full-CPL uses the coupled component restart set below. AMIP does not require
# MPAS-O, but the active initialization checks still require MPAS-I and coupler.
################################################################################
export my_refcase="20241109.v3-LR.DATM.ne30pg2_r05_IcoswISC30E3r5.pm-cpu"
export my_refdate=${my_casedate}
export my_reftod=${my_casetod}
export my_refdir="/compyfs/zhan391/acme_init/E3SMv3_INT/${my_refdate}-${my_reftod}"
export my_refeam_in="/compyfs/zhan391/acme_init/E3SMv3_INT/0101-01-01-00000/v3.LR.piControl.eam.i.0101-01-01-00000.nc"
export my_refeam_ic="${my_refdir}/v3.LR.ne30pg2.ERA5.eam.i.${my_refdate}-${my_reftod}.nc"
export my_refelm_in="${my_refdir}/${my_refcase}.elm.r.${my_refdate}-${my_reftod}.nc"
export my_refrof_in="${my_refdir}/${my_refcase}.mosart.r.${my_refdate}-${my_reftod}.nc"
export my_refocn_in="${my_refdir}/${my_refcase}.mpaso.rst.${my_refdate}_${my_reftod}.nc"
export my_refice_in="${my_refdir}/${my_refcase}.mpassi.rst.${my_refdate}_${my_reftod}.nc"
export my_refcpl_in="${my_refdir}/${my_refcase}.cpl.r.${my_refdate}-${my_reftod}.nc"

################################################################################
# --- Shared E3SM cycling controls -------------------------------------------
# Live completed-cycle counter. Step 4 updates this value transactionally.
################################################################################
export my_e3sm_completed_cycles=0
export my_e3sm_cycle_hours=6
export my_raw_archive_layout="per_member"
export my_shared_archive_dir="${my_modeldir}/archive"
export my_dart_root="${my_modeldir}/dart_en$(printf '%02d' "${my_ensnum}")"
my_member_archive_dir() {
  local member="${1:-}"
  [[ "${member}" =~ ^EN[0-9][0-9]$ ]] || { echo "invalid ensemble member: ${member:-unset}" >&2; return 1; }
  case "${my_raw_archive_layout}" in
    shared) printf '%s\n' "${my_shared_archive_dir}" ;;
    per_member) printf '%s\n' "${my_modeldir}/${member}/archive" ;;
    *) echo "invalid my_raw_archive_layout: ${my_raw_archive_layout}" >&2; return 1 ;;
  esac
}
export -f my_member_archive_dir
export my_e3sm_start_date=${my_casedate}
export my_e3sm_start_tod=${my_casetod}
export my_e3sm_end_date="2012-01-04"
export my_e3sm_end_tod="00000"

################################################################################
# --- EAM DART assimilation --------------------------------------------------
# Step 4 derives component nodes from my_job_nnodes and both DA switches:
# both enabled = equal halves; one enabled = all nodes; both disabled = no DA launch.
################################################################################
export my_eam_dart_da="on"
export my_eam_dart_cycle_hours=6
export my_eam_dart_end_date="${my_e3sm_end_date}"
export my_eam_dart_end_tod="${my_e3sm_end_tod}"
export my_eam_dart_run_dir="${my_dart_root}/eam"
export my_eam_dart_model="eam-se"
export my_eam_dart_pgrid=".true."
export my_eam_dart_code="/qfs/people/zhan391/e3sm_dart_work/code/DART"
# Lowest EAM model level at which observations may be assimilated; levels
# above it (smaller indices) are excluded near the diffusive model top.
export my_eam_no_obs_assim_above_level=5
export my_eam_use_log_vertical_scale=".true."
export my_eam_vert_normalization_scale_height="1.5"
# Use a workflow-owned loader rather than a generated file in DARTs work tree.
export my_eam_dart_obsdir="/compyfs/zhan391/acme_init/Observations/NCEP+ACARS+GPS+AIRS"
# Optional per-cycle EAM DA settings. Keys combine the exact valid time and
# parameter name. Omitted parameters retain their normal defaults.
declare -Ag my_eam_cycle_overrides=(
)

validate_my_eam_cycle_overrides() {
  local key stamp parameter value ymd tod
  for key in "${!my_eam_cycle_overrides[@]}"; do
    if [[ ! "${key}" =~ ^([0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{5}):(localization_cutoff|inflation_damping|no_obs_assim_above_level)$ ]]; then
      echo "invalid my_eam_cycle_overrides key: ${key}" >&2
      return 1
    fi
    stamp=${BASH_REMATCH[1]}
    parameter=${BASH_REMATCH[2]}
    ymd=${stamp:0:10}
    tod=${stamp:11:5}
    date -d "${ymd}" +%F >/dev/null 2>&1 || { echo "invalid override date: ${key}" >&2; return 1; }
    (( 10#${tod} < 86400 )) || { echo "override time is outside 00000-86399: ${key}" >&2; return 1; }
    value=${my_eam_cycle_overrides[${key}]}
    case "${parameter}" in
      localization_cutoff)
        [[ "${value}" =~ ^[0-9]+([.][0-9]+)?$ ]] && awk -v v="${value}" 'BEGIN {exit !(v > 0)}' || { echo "invalid localization cutoff for ${stamp}: ${value}" >&2; return 1; }
        ;;
      inflation_damping)
        [[ "${value}" =~ ^[0-9]+([.][0-9]+)?$ ]] && awk -v v="${value}" 'BEGIN {exit !(v >= 0 && v <= 1)}' || { echo "invalid inflation damping for ${stamp}: ${value}" >&2; return 1; }
        ;;
      no_obs_assim_above_level)
        [[ "${value}" =~ ^[1-9][0-9]*$ ]] && (( value <= 72 )) || { echo "invalid model-top cutoff level for ${stamp}: ${value}" >&2; return 1; }
        ;;
    esac
  done
}



################################################################################
# --- ELM DART assimilation --------------------------------------------------
# ELM DA links single-record h1 history and h2 vector files from each member archive.
################################################################################
export my_elm_dart_da="on"
export my_elm_dart_cycle_hours=24
export my_elm_dart_end_date="${my_e3sm_end_date}"
export my_elm_dart_end_tod="${my_e3sm_end_tod}"
export my_elm_dart_run_dir="${my_dart_root}/elm"
export my_elm_dart_code="/qfs/people/zhan391/e3sm_dart_work/code/DART_SCP"
export my_elm_sourcemods_dir="${my_elm_dart_code}/models/elm/DART_SourceMods/e3sm_maint_3.0/src.elm"
export my_elm_dart_obsdir="/compyfs/zhan391/acme_init/Observations/SMAP"
export my_elm_dart_model="elm"
export my_elm_history_stream="h1"
export my_elm_vector_history_stream="h2"
declare -Ag my_elm_dart_cycle_overrides=(
)

################################################################################
# --- Strongly coupled DA setup ----------------------------------------------
################################################################################
export strongly_coupled_on="off"

export atm_da_compute_posterior=".false."
export atm_da_output_sequential_prior_post=".false."
export atm_da_use_sequential_prior_post=".false."
export atm_da_output_mean=".true."
export atm_da_output_sd=".true."
export atm_da_output_members=".true."
export atm_da_strongly_coupled=".false."
export atm_da_state_model="Atmosphere"
export atm_da_obs_model="Atmosphere"

export lnd_da_output_sequential_prior_post=".false."
export lnd_da_use_sequential_prior_post=".false."
export lnd_da_perturb_from_single_instance=".true."
export lnd_da_perturbation_amplitude="0.2"
export lnd_da_perturbation_method="uniform"
export lnd_da_obs_sequence_in_name="obs_seq.out"
export lnd_da_inf_flavor_prior="0"
export lnd_da_inf_flavor_posterior="0"
export lnd_da_cutoff="0.4"
export lnd_da_spread_restoration=".false."
export lnd_da_sampling_error_correction=".false."
export lnd_da_horiz_dist_only=".true."
export lnd_da_strongly_coupled=".false."
export lnd_da_state_model="Land"
export lnd_da_obs_model="Land"

################################################################################
# --- EAM DART diagnostic defaults ------------------------------------------
# Step 6 uses this range when DIAG_START and DIAG_END are empty.
################################################################################
export my_eam_dart_diag_start="2011-12-01-00000"
export my_eam_dart_diag_end="2011-12-28-00000"
# Diagnostic switches use Fortran logical strings because the workers pass them to DART tools.
export my_eam_dart_diag_use_custom_range=".true."
export my_eam_dart_diag_run_closest_member=".false."
export my_eam_dart_diag_run_obs2netcdf=".false."
export my_eam_dart_diag_run_obs_diag=".true."
