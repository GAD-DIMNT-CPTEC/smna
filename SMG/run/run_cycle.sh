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
#BOP
# !FUNCTION: ensure_root
# !DESCRIPTION: Discover project root by locating ".smg_root", export SMG_ROOT,
#               and source "$SMG_ROOT/etc/__init__.sh" exactly once (idempotent).
# !USAGE: Place near the top of any script and call: ensure_root <ROOT_PATH_VAR> <INIT_LOADED_VAR>|| exit $?
# !NOTE: Requires bash; uses PWD/BASH_SOURCE and pwd -P (no readlink -f).
#EOP
#EOC
# Ensure SMG_ROOT exists and initialize once
ensure_root() {
  local root_var="${1:?root_var required}" loaded_var="${2:?load_var required}"
  local marker="${ROOT_MARKER:-.smg_root}" init_rel="${INIT_REL:-etc/__init__.sh}"
  local src dir found init

  # 1) Root discovery (only if not set or invalid)
  if [[ -z "${!root_var:-}" || ! -f "${!root_var}/$marker" ]]; then
    for src in "${BASH_SOURCE[@]:-"$0"}" "$PWD"; do
      [[ -n "$src" ]] || continue
      dir=$([[ -d "$src" ]] && cd -- "$src" && pwd -P || cd -- "$(dirname -- "$src")" && pwd -P) || return 1
      while [[ "$dir" != "/" && ! -f "$dir/$marker" ]]; do dir="${dir%/*}"; done
      if [[ -f "$dir/$marker" ]]; then found="$dir"; break; fi
    done
    [[ -n "${found:-}" ]] || { printf '[ERROR] %s not found\n' "$marker" >&2; return 1; }
    printf -v "$root_var" %s "$found"; export "$root_var"
  fi
  # 2) Sanity check for init
  init="${!root_var}/${init_rel}"
  [[ -r "$init" ]] || { printf '[ERROR] Missing %s\n' "$init" >&2; return 2; }

  # 3) Load init exactly once (idempotent)
  if [[ "${!load_var:-0}" != 1 ]]; then
    # shellcheck source=/dev/null
    . "$init" || { printf '[ERROR] Failed to load %s\n' "$init" >&2; return 3; }
    printf -v "$load_var" 1
    export "$load_var"
  fi
}
ensure_root SMG_ROOT SMG_INIT_LOADED || exit $?
#EOC

#------------------------------ Mini sanity ----------------------------------#
command -v date >/dev/null 2>&1 || _die 3 "'date' not found in PATH"

# Diretório de scripts
: "${scripts_smg:=${SMG_ROOT}/run/scripts}"
[[ -d "${scripts_smg}" ]] || _die 4 "scripts_smg directory not found: %s" "${scripts_smg}"


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

#------------------------------ Parse comum + local --------------------------#
# 1) Parser global (__helpers__.sh) — flags comuns e HPC
__parse_args__ "$@" || exit $?

# 2) Local parser — only cycle flags, using leftover_args from the global parser
# !FUNCTION: _parse_args
# !DESCRIPTION: Parse CLI options in a portable way (Bash 4+).
_parse_cycle_args() {
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

      -h|--help)    usage; exit 0;;
      --) shift; break;;
      -*) _log_warn "Unknown option ignored (cycle): %s" "$1"; shift;;
      *)  break;;
   esac
  done
# any positional arguments after "--" remain in leftover_args (from the global parser)
return 0

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

#BOP
# !FUNCTION: _run_obsmake
# !INTERFACE: _run_obsmake LABEL
# !DESCRIPTION:
#   Execute the OBSMAKE stage for a given analysis label (LABEL).
#
#   • Automatically translates global flags from __parse_args__ into
#     CLI options for run_obsmake.sh (--verbose/--quiet, --dry-run).
#   • Supports environment overrides for source/destination/mode/pattern
#     through OBSM_NCEP_ROOT, OBSM_DEST_DIR, OBSM_MODE, OBSM_PATTERN.
#   • Delegates dry-run and logging behavior to run_obsmake.sh itself,
#     instead of wrapping with eval/run() here.
#
#   Arguments:
#     LABEL   → Analysis datetime label (YYYYMMDDHH).
#
#   Behavior:
#     - Builds a CLI argument list (`obsmake_cli`) based on global variables.
#     - Calls run_obsmake.sh with LABEL and the constructed options.
#     - On error, aborts with code 10.
#     - Logs elapsed runtime in minutes/seconds.
#
# !USAGE:
#   _run_obsmake 2025010100
#
# !ENVIRONMENT:
#   OBSM_NCEP_ROOT  → Passed as --ncep-root if set.
#   OBSM_DEST_DIR   → Passed as --dest if set.
#   OBSM_MODE       → Passed as --mode (link|copy|hardlink) if set.
#   OBSM_PATTERN    → Passed as --pattern if set.
#   verbose         → From __parse_args__, mapped to --verbose/--quiet.
#   dry_run         → From __parse_args__, mapped to --dry-run.
#
# !RETURNS:
#   0 on success; aborts with code 10 on failure.
#EOP
#BOC
_run_obsmake() {
  local lbl="$1"

  # --- Build CLI argument list for run_obsmake.sh ---
  local obsmake_cli=()

  # Verbosity: prefer explicit CLI flag
  if [[ "${verbose:-false}" == true ]]; then
    obsmake_cli+=("--verbose")
  else
    obsmake_cli+=("--quiet")
  fi

  # Dry-run support
  [[ "${dry_run:-false}" == true ]] && obsmake_cli+=("--dry-run")

  # Optional environment overrides
  [[ -n "${OBSM_NCEP_ROOT:-}" ]] && obsmake_cli+=("--ncep-root" "${OBSM_NCEP_ROOT}")
  [[ -n "${OBSM_DEST_DIR:-}"  ]] && obsmake_cli+=("--dest"      "${OBSM_DEST_DIR}")
  [[ -n "${OBSM_MODE:-}"      ]] && obsmake_cli+=("--mode"      "${OBSM_MODE}")      # link|copy|hardlink
  [[ -n "${OBSM_PATTERN:-}"   ]] && obsmake_cli+=("--pattern"   "${OBSM_PATTERN}")

  _log_info "Executing OBSMAKE for %s" "${lbl}"
  SECONDS=0

  # --- Call run_obsmake.sh directly ---
  /bin/bash "${scripts_smg}/run_obsmake.sh" "${lbl}" "${obsmake_cli[@]}"
  local rc=$?

  (( rc == 0 )) || _die 10 "Observer failed for %s (rc=%d)" "${lbl}" "${rc}"
  _log_ok "Observer finished in %dm%02ds" "$((SECONDS/60))" "$((SECONDS%60))"
}
#EOC

#BOP
# !FUNCTION: _run_gsi
# !INTERFACE: _run_gsi LABEL [BC_CYCLE]
# !DESCRIPTION:
#   Execute a single GSI analysis for the given analysis label (LABEL).
#
#   • By default, runs a "normal" assimilation cycle.
#   • Consumes cycle-local settings (modelTrunc, modelNLevs, modelPrefix,
#     gsiTrunc, gsiMPITasks) and global HPC layout (mpi_tasks, omp_threads)
#     parsed previously by __parse_args__.
#   • Avoids passing unknown flags to runGSI. Only uses runGSI-supported
#     options (-I, -T, -t, -l, -p, -np, -bc). If OMP threads are provided by
#     the global parser, exports OMP_NUM_THREADS to benefit the GSI run.
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
#     - Resolves effective MPI tasks for GSI: prefer gsiMPITasks (-gnp),
#       otherwise fallback to global mpi_tasks (if provided).
#     - Exports OMP_NUM_THREADS when global --cpus-per-task was given.
#     - Invokes runGSI with a minimal, valid CLI (only supported flags).
#     - On failure, aborts with code 20.
#
# !USAGE:
#   _run_gsi 2025010100
#   _run_gsi 2025010100 3
#   _run_gsi 2025010100 0
#
# !ENVIRONMENT:
#   omp_threads     → If set (from __parse_args__), exported as OMP_NUM_THREADS.
#   verbose, dry_run→ Optionally exported for downstream scripts to consult;
#                     not passed as flags (runGSI does not declare them).
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

  # --- Resolve effective MPI task count for GSI (-np) ---
  # Prefer the cycle-local -gnp (gsiMPITasks); otherwise fallback to global --ntasks (mpi_tasks)
  local gnp=
  if [[ -n "${gsiMPITasks:-}" ]]; then
    gnp="${gsiMPITasks}"
  elif [[ -n "${mpi_tasks:-}" ]]; then
    gnp="${mpi_tasks}"
  fi
  
  SECONDS=0

  # --- Build runGSI CLI (only supported flags) ---
  # -I label, -T gsiTrunc (analysis trunc), -t modelTrunc (background trunc),
  # -l modelNLevs, -p prefix, -np MPI, -bc BC index (0 for normal mode).
  local argv=(
    -I "${lbl}"
    -T "${gsiTrunc}"
    -t "${modelTrunc}"
    -l "${modelNLevs}"
    -p "${modelPrefix}"
    -bc "${bc_cycle}"
  )
  [[ -n "${gnp:-}" ]] && argv+=(-np "${gnp}")

  # --- Invoke runGSI (no wrapper; let runGSI handle submission/logging) ---
  /bin/bash "${scripts_smg}/runGSI" "${argv[@]}"
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
 
  # 1) Parser global (__helpers__.sh) — flags comuns e HPC
  __parse_args__ "$@" || exit $?

  # 2) Local parser — only cycle flags, using leftover_args from the global parser
  _parse_cycle_args "${leftover_args[@]}" || exit $?

  # 3) Validate required options and derive defaults.
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

