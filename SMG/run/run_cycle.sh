#!/usr/bin/env bash
#-----------------------------------------------------------------------------#
#BOP
# !SCRIPT: run_cycle.sh — Orchestrate OBS → GSI → BAM cycles (ProTex-documented)
# !DESCRIPTION:
#   Robust driver to run a cyclic NWP workflow composed of three main stages:
#     1) OBSMAKE  — Prepare observations for data assimilation
#     2) GSI      — Run the GSI analysis (optionally with radiance bias spin-up)
#     3) BAM      — Run the forecast model using the analysis produced by GSI
#
#   This refactor preserves the original logic while improving:
#     • Safety     → set -Eeuo pipefail, strict IFS, defensive sourcing
#     • Clarity    → cohesive help/usage, function-level logging, explicit checks
#     • Portability→ avoids non-portable tools; requires Bash 4+ and GNU `date`
#     • Extensibility→ bias-correction spin-up switch, feature toggles per stage
#
# !USAGE:
#   run_cycle.sh [options]
#
# !OPTIONS:
#   Model options (BAM background):
#     -t  <int>       Model background truncation (modelTrunc)              [REQ]
#     -l  <int>       Model background levels (modelNLevs)                  [REQ]
#     -p  <str>       Model background file prefix (modelPrefix)            [REQ]
#     -mnp <int>      MPI ranks for BAM (modelMPITasks)                     [64]
#     -fct <int>      Forecast length in hours (modelFCT)                   [09]
#
#   GSI options:
#     -gt  <int>      GSI analysis truncation (gsiTrunc) [default: modelTrunc]
#     -gnp <int>      MPI ranks for GSI (gsiMPITasks)                        [144]
#
#   Cycle window:
#     -I <yyyymmddhh> Initial analysis datetime (LABELI)                    [REQ]
#     -F <yyyymmddhh> Final   analysis datetime (LABELF)                    [REQ]
#
#   Feature toggles:
#     --no-obsmake    Skip observer pre-processing (do_obsmake=0)
#     --no-gsi        Skip GSI analysis (do_gsi=0)
#     --no-bam        Skip BAM forecast (do_bam=0)
#     -q|--quiet      Less logging (verbose=false)
#     -v|--verbose    More logging (verbose=true)
#
#   Bias-correction spin-up (radiances):
#     -bc   <N>          Repeat the *same* analysis N times to spin up satbias
#                         & angle-dependent coefficients at a given datetime.
#
# !NOTES:
#   • The spin-up follows common NCEP guidance for short experiments when a
#     long cycling period is not available. It reuses identical BG+OBS while
#     feeding back updated satbias/angle coefficients each iteration.
#   • GNU `date` is required for date arithmetic (e.g., +6 hours, +fct hours).
#   • Environment: this script honors SMG_ROOT when already defined (e.g., via
#     module/wrapper). If unset, it attempts a robust derivation from the script
#     location. Sourcing of config_smg.ksh is defensive to avoid killing the run
#     under `set -e` if the config returns non-zero.
#
# !EXAMPLES:
#   # Production cycle every +6h from I to F (OBS→GSI→BAM):
#   run_cycle.sh -t 299 -l 64 -p PREF -I 2025010100 -F 2025010300 -gnp 144 -mnp 64 -fct 9 -v
#
#   # Bias spin-up only (10 iterations at -I), skipping OBS and BAM:
#   run_cycle.sh -t 299 -l 64 -p PREF -I 2025010100 -F 2025010100 -gt 299 -bc 10 --no-obsmake --no-bam -v
#
#   # Bias spin-up with 8 iterations at first label:
#   run_cycle.sh -t 299 -l 64 -p PREF -I 2025010100 -F 2025010300 -gt 299 -bc 8
#
# !AUTHOR: GDAD/CPTEC/INPE
# !DATE:   2025-09-20
#EOP
#-----------------------------------------------------------------------------#
#BOC
# --- Library directory (folder where this file lives) ---
__SRC="${BASH_SOURCE[0]}"
while [[ -L "$__SRC" ]]; do
  __LINK="$(readlink -- "$__SRC")"
  if [[ "$__LINK" = /* ]]; then
    __SRC="$__LINK"
  else
    __SRC="$(cd -- "$(dirname -- "$__SRC")" && cd -- "$(dirname -- "$__LINK")" && pwd)/$(basename -- "$__LINK")"
  fi
done
export RUN_CYCLE_DIR="$(cd -- "$(dirname -- "$__SRC")" && pwd -P)"

# --- load init (which loads helpers) ---
__init_path="${RUN_CYCLE_DIR}/__init__.sh"
# shellcheck disable=SC1090
. "$__init_path" || :        # no extra prints; don't break under set -e

trap '_log_err "Interrupted (SIGINT)"; exit 130' INT
trap '_log_err "Terminated (SIGTERM)"; exit 143' TERM
trap '_log_err "Failure at line %d" "$LINENO"' ERR

#------------------------- Environment bootstrap -----------------------------#
# Prefer environment-provided SMG_ROOT. If missing, derive from script path.
# --- Load system variables ---
if [[ -z "${SMG_ROOT:-}" ]]; then
  _log_warn "SMG_ROOT not set; deriving from script location"

  # Start from the directory of this script
  cur_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

  while [[ "$cur_dir" != "/" ]]; do
    if [[ -f "${cur_dir}/config_smg.ksh" ]]; then
      SMG_ROOT="$cur_dir"
      break
    fi
    cur_dir="$(dirname -- "$cur_dir")"
  done

  # Fallback: if not found, use the script's own directory
  if [[ -z "${SMG_ROOT:-}" ]]; then
    _log_warn "config_smg.ksh not found while climbing; using script dir"
    SMG_ROOT="$cur_dir"
  fi
fi

export SMG_ROOT

: "${scripts_smg:=${SMG_ROOT}/run/scripts}"
[[ -d "${scripts_smg}" ]] || _log WARN "scripts_smg directory not found: %s" "${scripts_smg}"

# Defensive sourcing: do not kill the run if config returns non-zero
if [[ -f "${SMG_ROOT}/config_smg.ksh" ]]; then
  # shellcheck disable=SC1090
  if ! . "${SMG_ROOT}/config_smg.ksh" _env_export; then
    _log_warn "config_smg.ksh returned non-zero; continuing"
  fi
fi

command -v date >/dev/null 2>&1 || _die 3 "'date' not found in PATH"

#------------------------------ Usage helper ---------------------------------#
# !FUNCTION: usage
# !DESCRIPTION: Print concise usage extracted from the ProTex header.
usage() {
  sed -n '/^# !USAGE:/,/^# !AUTHOR:/p' "${BASH_SOURCE[0]}" | sed 's/^# //'
}

#------------------------------ Defaults -------------------------------------#
modelMPITasks=64
modelFCT=09
gsiMPITasks=144

do_obsmake=1
do_gsi=1
do_bam=1

modelTrunc=""
modelNLevs=""
modelPrefix=""
gsiTrunc=""

LABELI=""
LABELF=""

BcCycles=0

#------------------------------ Option parsing -------------------------------#
# !FUNCTION: _parse_args
# !DESCRIPTION: Parse CLI options in a portable way (Bash 4+).
_parse_args() {
  (($#)) || { usage; return 1; }
  while (($#)); do
    case "$1" in
      -t)   modelTrunc="${2:?missing -t}"; shift 2;;
      -l)   modelNLevs="${2:?missing -l}"; shift 2;;
      -p)   modelPrefix="${2:?missing -p}"; shift 2;;
      -mnp) modelMPITasks="${2:?missing -mnp}"; shift 2;;

      -I)   LABELI="${2:?missing -I}"; shift 2;;
      -F)   LABELF="${2:?missing -F}"; shift 2;;
      -fct) modelFCT="${2:?missing -fct}"; shift 2;;

      -gt)  gsiTrunc="${2:?missing -gt}"; shift 2;;
      -gnp) gsiMPITasks="${2:?missing -gnp}"; shift 2;;

      -bc)  BcCycles="${2:?missing -bc}"; shift 2;;

      --no-obsmake) do_obsmake=0; shift;;
      --no-gsi)     do_gsi=0; shift;;
      --no-bam)     do_bam=0; shift;;

      -q|--quiet)   verbose=false; shift;;
      -v|--verbose) verbose=true; shift;;
      -h|--help)    usage; exit 0;;
      *)            _log WARN "Unknown option ignored: %s" "$1"; shift;;
    esac
  done
}

#------------------------------- Helpers -------------------------------------#
# !FUNCTION: _validate_required
# !DESCRIPTION: Validate required options and derive defaults.
_validate_required() {
  [[ -n "${modelTrunc}"  ]] || _die 2 "Missing -t (modelTrunc)"
  [[ -n "${modelNLevs}"  ]] || _die 2 "Missing -l (modelNLevs)"
  [[ -n "${modelPrefix}" ]] || _die 2 "Missing -p (modelPrefix)"
  [[ -n "${LABELI}"      ]] || _die 2 "Missing -I (LABELI)"
  [[ -n "${LABELF}"      ]] || _die 2 "Missing -F (LABELF)"
  if [[ -z "${gsiTrunc}" ]]; then gsiTrunc="${modelTrunc}"; fi
}

# !FUNCTION: _run_obsmake
# !DESCRIPTION: Execute observer pre-processing for a given analysis label.
_run_obsmake() {
  local lbl="$1"
  _log_info "Executing OBSMAKE for %s" "${lbl}"
  SECONDS=0
  /bin/bash "${scripts_smg}/run_obsmake.sh" "${lbl}"
  local rc=$?
  ((rc==0)) || _die 10 "Observer failed for %s (rc=%d)" "${lbl}" "${rc}"
  _log_ok   "Observer finished in %dm%02ds" "$((SECONDS/60))" "$((SECONDS%60))"
}

#BOP
# !FUNCTION: _run_gsi
# !INTERFACE: _run_gsi LABEL [BC_CYCLE]
# !DESCRIPTION:
#   Execute a single GSI analysis for the given analysis label (LABEL).
#
#   • By default, runs a "normal" assimilation cycle.
#   • When BC_CYCLE > 0, the run is considered part of a bias-correction (BC)
#     cycling sequence, and the BC_CYCLE index is passed downstream so that
#     other functions (e.g. CopyOutputsForCycle) can manage satbias file
#     handoff correctly.
#
#   Arguments:
#     LABEL     → Analysis datetime label (YYYYMMDDHH)
#     BC_CYCLE  → [optional] integer cycle index; defaults to 0
#
#   Behavior:
#     - Calls the runGSI driver script with the required truncation,
#       levels, prefix, and MPI task count.
#     - On error, aborts with code 20.
#     - Logs elapsed runtime in minutes/seconds.
#
# !USAGE:
#   _run_gsi 2025010100
#   _run_gsi 2025010100 3
#   _run_gsi 2025010100 0
#
# !NOTES:
#   - BC_CYCLE is passed as options to runGSI so that bias
#     correction cycling can be handled consistently across scripts.
#EOP
#BOC
_run_gsi() {
  local lbl=${1:? "Missing LABEL"}
  local bc_cycle=${2:-0}

  _log_info "Executing GSI for %s" "${lbl}" 
  if (( bc_cycle > 0 ));then
     _log_action "running bias corretion, bc_cycle=%d" "${bc_cycle}"
  fi

  SECONDS=0
  /bin/bash "${scripts_smg}/runGSI" \
    -I "${lbl}" \
    -T "${gsiTrunc}" \
    -t "${modelTrunc}" \
    -l "${modelNLevs}" \
    -p "${modelPrefix}" \
    -np "${gsiMPITasks}" \
    -bc "${bc_cycle}"

  local rc=$?
  ((rc==0)) || _die 20 "GSI failed for %s (rc=%d)" "${lbl}" "${rc}"
  _log_info "GSI finished in %dm%02ds" "$((SECONDS/60))" "$((SECONDS%60))"
}
#EOC

# !FUNCTION: _run_bam
# !DESCRIPTION: Execute BAM forecast from a given analysis label to +modelFCT h.
_run_bam() {
  local lbl="$1"
  local fct_date
  fct_date="$(date -u +%Y%m%d%H -d "${lbl:0:8} ${lbl:8:2} +${modelFCT} hours")"
  _log_info "Executing BAM: start=%s  fct_end=%s" "${lbl}" "${fct_date}"
  SECONDS=0
  /bin/bash "${scripts_smg}/run_model.sh" "${lbl}" "${fct_date}" "${modelPrefix}" "${modelTrunc}" "${modelNLevs}" "${modelMPITasks}" "No"
  local rc=$?
  ((rc==0)) || _die 30 "BAM failed for %s (rc=%d). Check PRE/BAM/POS." "${lbl}" "${rc}"
  _log_ok   "BAM finished in %dm%02ds" "$((SECONDS/60))" "$((SECONDS%60))"
}

#BOP
# !FUNCTION: _spinup_bias_once
# !INTERFACE: _spinup_bias_once LABEL [N_ITERS]
# !DESCRIPTION:
#   Perform multiple consecutive GSI runs at a single analysis time, reusing the
#   same background (BG) and observations (OBS), in order to spin up radiance
#   bias and/or angle-dependent bias coefficients.
#
#   Arguments:
#     LABEL     → Analysis datetime label (YYYYMMDDHH)
#     N_ITERS   → [optional] Number of iterations; defaults to 10
#
#   Behavior:
#     - Executes N consecutive GSI runs at the same LABEL.
#     - Each run calls `_run_gsi` with the same inputs, optionally in BC mode.
#     - Useful for short spin-up experiments of bias correction when no long
#       cycling is available.
#
# !USAGE:
#   _spinup_bias_once 2025010100          # 10 iterations, normal mode
#   _spinup_bias_once 2025010100 5        # 5 iterations
#
#EOP
#BOC
_spinup_bias_once() {
  local lbl=${1:? "Missing LABEL"}
  local n=${2:-10}

  _log_info "Bias-correction spin-up: %d iterations @ %s" \
    "${n}" "${lbl}"

  for ((i=1; i<=n; i++)); do
    _log_info "Spin-up iteration %d/%d" "$i" "$n"
    _run_gsi "${lbl}" ${i}
  done

  _log_ok "Bias-correction spin-up completed"
}
#EOC

#------------------------------- Main ----------------------------------------#
_main() {
   
   local verbose=true

  _parse_args "$@" || return 1
  _validate_required

  _log_info "Cycle window: %s → %s (step +6h)" "${LABELI}" "${LABELF}"
  _log_info "Model: trunc=%s nlevs=%s prefix=%s mnp=%s fct=%sh" \
             "${modelTrunc}" "${modelNLevs}" "${modelPrefix}" "${modelMPITasks}" "${modelFCT}"
  _log_info "GSI:   trunc=%s gnp=%s" "${gsiTrunc}" "${gsiMPITasks}"

  # Optional bias-correction spin-up
  if (( BcCycles > 0 )); then
      _spinup_bias_once "${LABELI}" "${BcCycles}"
  fi

  # Production cycle
  local current="${LABELI}"
  while [[ "${current}" -le "${LABELF}" ]]; do
    _log_info "Submitting system for date %s" "${current}"

    (( do_obsmake )) && _run_obsmake "${current}"
    (( do_gsi     )) && _run_gsi     "${current}"
    (( do_bam     )) && _run_bam     "${current}"

    current="$(date -u +%Y%m%d%H -d "${current:0:8} ${current:8:2} +6 hours")"
  done

  _log_info "run_cycle completed successfully"
}

_main "$@"
#EOC
#-----------------------------------------------------------------------------#

