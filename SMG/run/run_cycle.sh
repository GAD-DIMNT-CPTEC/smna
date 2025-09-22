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
set -Eeuo pipefail
IFS=$' \t\n'

#-------------------------------- Logging ------------------------------------#
: "${verbose:=true}"     # default: print INFO
: "${COLOR:=1}"          # set COLOR=0 to disable ANSI coloring
if [[ -t 2 && "${COLOR}" = "1" ]]; then
  C_INFO=$'\033[1;34m'; C_OK=$'\033[1;32m'; C_WARN=$'\033[1;33m'; C_ERR=$'\033[1;31m'; C_RST=$'\033[0m'
else
  C_INFO=; C_OK=; C_WARN=; C_ERR=; C_RST=
fi

# !FUNCTION: _log
# !DESCRIPTION: Print standardized log lines. Shows INFO only when verbose.
_log() { # _log <LEVEL> <fmt> [args...]
  local lvl="$1"; shift || true
  $verbose || [[ "$lvl" =~ ^(OK|WARN|ERR)$ ]] || return 0
  local fmt="$1"; shift || true
  local tag="[${lvl}]"
  case "$lvl" in
    INFO) printf "%s%s%s %s\n" "$C_INFO" "$tag" "$C_RST" "$(printf "$fmt" "$@")" ;;
    OK)   printf "%s%s%s %s\n" "$C_OK"   "$tag" "$C_RST" "$(printf "$fmt" "$@")" ;;
    WARN) printf "%s%s%s %s\n" "$C_WARN" "$tag" "$C_RST" "$(printf "$fmt" "$@")" ;;
    ERR)  printf "%s%s%s %s\n" "$C_ERR"  "$tag" "$C_RST" "$(printf "$fmt" "$@")" ;;
    *)    printf "%s %s\n" "$tag" "$(printf "$fmt" "$@")" ;;
  esac
}
_die() { local code="${1:-1}"; shift || true; _log ERR "$@"; exit "$code"; }

trap '_log ERR "Interrupted (SIGINT)"; exit 130' INT
trap '_log ERR "Terminated (SIGTERM)"; exit 143' TERM
trap '_log ERR "Failure at line %d" "$LINENO"' ERR

#------------------------- Environment bootstrap -----------------------------#
# Prefer environment-provided SMG_ROOT. If missing, derive from script path.
if [[ -z "${SMG_ROOT:-}" ]]; then
  _log WARN "SMG_ROOT not set; deriving from script location"
  self_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
  if [[ -d "${self_dir}/../config" || -f "${self_dir}/../config_smg.ksh" ]]; then
    SMG_ROOT="$(cd -- "${self_dir}/.." && pwd)"
  else
    SMG_ROOT="${self_dir}"
  fi
fi
export SMG_ROOT

: "${scripts_smg:=${SMG_ROOT}/run/scripts}"

# Defensive sourcing: do not kill the run if config returns non-zero
if [[ -f "${SMG_ROOT}/config_smg.ksh" ]]; then
  # shellcheck disable=SC1090
  if ! . "${SMG_ROOT}/config_smg.ksh" vars_export; then
    _log WARN "config_smg.ksh returned non-zero; continuing"
  fi
fi

command -v date >/dev/null 2>&1 || _die 3 "'date' not found in PATH"
[[ -d "${scripts_smg}" ]] || _log WARN "scripts_smg directory not found: %s" "${scripts_smg}"

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
BcLABELI=""
BcLABELF=""

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
  _log INFO "Executing OBSMAKE for %s" "${lbl}"
  SECONDS=0
  /bin/bash "${scripts_smg}/run_obsmake.sh" "${lbl}"
  local rc=$?
  ((rc==0)) || _die 10 "Observer failed for %s (rc=%d)" "${lbl}" "${rc}"
  _log OK   "Observer finished in %dm%02ds" "$((SECONDS/60))" "$((SECONDS%60))"
}

#BOP
# !FUNCTION: _run_gsi
# !INTERFACE: _run_gsi LABEL [IS_BC] [BC_CYCLE]
# !DESCRIPTION:
#   Execute a single GSI analysis for the given analysis label (LABEL).
#
#   • By default, runs a "normal" assimilation cycle (IS_BC=false).
#   • When IS_BC=true, the run is considered part of a bias-correction (BC)
#     cycling sequence, and the BC_CYCLE index is passed downstream so that
#     other functions (e.g. CopyOutputsForCycle) can manage satbias file
#     handoff correctly.
#
#   Arguments:
#     LABEL     → Analysis datetime label (YYYYMMDDHH)
#     IS_BC     → [optional] "true" or "false"; defaults to "false"
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
#   _run_gsi 2025010100 true  3
#   _run_gsi 2025010100 false 0
#
# !NOTES:
#   - IS_BC and BC_CYCLE are passed as options to runGSI so that bias
#     correction cycling can be handled consistently across scripts.
#   - If IS_BC=false, BC_CYCLE is ignored.
#EOP
#BOC
_run_gsi() {
  local lbl=${1:? "Missing LABEL"}
  local bc_cycle=${2:-0}

  # derive is_bc from BcCycle
  local is_bc=false
  (( bc_cycle > 0 )) && is_bc=true

  _log INFO "Executing GSI for %s (is_bc=%s, bc_cycle=%d)" \
    "${lbl}" "${is_bc}" "${bc_cycle}"

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
  _log OK "GSI finished in %dm%02ds" "$((SECONDS/60))" "$((SECONDS%60))"
}
#EOC

# !FUNCTION: _run_bam
# !DESCRIPTION: Execute BAM forecast from a given analysis label to +modelFCT h.
_run_bam() {
  local lbl="$1"
  local fct_date
  fct_date="$(date -u +%Y%m%d%H -d "${lbl:0:8} ${lbl:8:2} +${modelFCT} hours")"
  _log INFO "Executing BAM: start=%s  fct_end=%s" "${lbl}" "${fct_date}"
  SECONDS=0
  /bin/bash "${scripts_smg}/run_model.sh" "${lbl}" "${fct_date}" "${modelPrefix}" "${modelTrunc}" "${modelNLevs}" "${modelMPITasks}" "No"
  local rc=$?
  ((rc==0)) || _die 30 "BAM failed for %s (rc=%d). Check PRE/BAM/POS." "${lbl}" "${rc}"
  _log OK   "BAM finished in %dm%02ds" "$((SECONDS/60))" "$((SECONDS%60))"
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
#   _spinup_bias_once 2025010100 5 true 0 # 5 iterations, BC mode (cycle 0)
#
# !NOTES:
#   - If IS_BC=true, BC_CYCLE is passed downstream but does not increment
#     automatically. For full multi-cycle BC handling, use outer cycling logic.
#EOP
#BOC
_spinup_bias_once() {
  local lbl=${1:? "Missing LABEL"}
  local n=${2:-10}

  _log INFO "Bias-correction spin-up: %d iterations @ %s (is_bc=%s, bc_cycle=%d)" \
    "${n}" "${lbl}" "${is_bc}" "${bc_cycle}"

  for ((i=1; i<=n; i++)); do
    _log INFO "Spin-up iteration %d/%d" "$i" "$n"
    _run_gsi "${lbl}" ${i}
  done

  _log OK "Bias-correction spin-up completed"
}
#EOC

#------------------------------- Main ----------------------------------------#
_main() {
  _parse_args "$@" || return 1
  _validate_required

  _log INFO "Cycle window: %s → %s (step +6h)" "${LABELI}" "${LABELF}"
  _log INFO "Model: trunc=%s nlevs=%s prefix=%s mnp=%s fct=%sh" \
             "${modelTrunc}" "${modelNLevs}" "${modelPrefix}" "${modelMPITasks}" "${modelFCT}"
  _log INFO "GSI:   trunc=%s gnp=%s" "${gsiTrunc}" "${gsiMPITasks}"

  # Optional bias-correction spin-up
  if (( BcCycles > 0 )); then
      _spinup_bias_once "${BcLABELI}" "${BcCycles}"
  fi

  # Production cycle
  local current="${LABELI}"
  while [[ "${current}" -le "${LABELF}" ]]; do
    _log INFO "Submitting system for date %s" "${current}"

    (( do_obsmake )) && _run_obsmake "${current}"
    (( do_gsi     )) && _run_gsi     "${current}"
    (( do_bam     )) && _run_bam     "${current}"

    current="$(date -u +%Y%m%d%H -d "${current:0:8} ${current:8:2} +6 hours")"
  done

  _log OK "run_cycle completed successfully"
}

_main "$@"
#EOC
#-----------------------------------------------------------------------------#

