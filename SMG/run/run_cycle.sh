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
# !DESCRIPTION:
#   Discover the project root by searching for a marker file (default: ".smg_root").
#   Defines/exports the environment variable ROOT_VAR pointing to the root path,
#   without sourcing any init scripts. Idempotent.
# !USAGE:
#   ensure_root <ROOT_VAR> [MARKER='.smg_root']
# !EXAMPLE:
#   ensure_root SMG_ROOT           || exit $?
#   ensure_root SMG_ROOT .smg_root || exit $?
# !NOTES:
#   - Requires bash; relies on BASH_SOURCE, PWD, and "pwd -P" (no readlink -f).
#   - Walks upward from ${BASH_SOURCE} and $PWD until the MARKER is found.
#EOP
#BOC
ensure_root() {
  local root_var="${1:?root_var required}" marker="${2:-${ROOT_MARKER:-.smg_root}}"
  local src dir found

  # variable names must be valid identifiers
  [[ "$root_var" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ ]] || {
    printf '[ERROR] invalid var name: %s\n' "$root_var" >&2; return 1; }

  # if already set and valid, keep it
  if [[ -n "${!root_var:-}" && -f "${!root_var}/$marker" ]]; then
    export "$root_var"
    return 0
  fi

  # search from BASH_SOURCE chain and PWD
  for src in "${BASH_SOURCE[@]:-"$0"}" "$PWD"; do
    [[ -n "$src" ]] || continue
    if [[ -d "$src" ]]; then
      dir="$(cd -- "$src" && pwd -P)" || return 1
    else
      dir="$(cd -- "$(dirname -- "$src")" && pwd -P)" || return 1
    fi
    while [[ "$dir" != "/" && ! -f "$dir/$marker" ]]; do dir="${dir%/*}"; done
    if [[ -f "$dir/$marker" ]]; then found="$dir"; break; fi
  done

  [[ -n "${found:-}" ]] || { printf '[ERROR] %s not found\n' "$marker" >&2; return 1; }
  printf -v "$root_var" %s "$found"
  export "$root_var"
}
#EOC
#-----------------------------------------------------------------------------#
#BOP
# !FUNCTION: smg_source_init_once
# !DESCRIPTION:
#   Sources "$<ROOT_VAR>/<INIT_REL>" exactly once, controlled by a flag LOADED_VAR.
#   Does not attempt to discover the project root: assumes <ROOT_VAR> was already
#   set/exported by smg_ensure_root (or another mechanism).
# !USAGE:
#   smg_source_init_once <ROOT_VAR> <LOADED_VAR> [INIT_REL='etc/__init__.sh']
# !EXAMPLE:
#   smg_ensure_root SMG_ROOT || exit $?
#   smg_source_init_once SMG_ROOT SMG_INIT_LOADED || exit $?
# !NOTES:
#   - Idempotent: if LOADED_VAR==1, it will not re-source the init.
#   - Return codes: 2 if init missing/unreadable, 3 if sourcing fails.
#EOP
#BOC
source_init_once() {
  local root_var="${1:?root_var required}" loaded_var="${2:?loaded_var required}"
  local init_rel="${3:-${INIT_REL:-etc/__init__.sh}}"
  local init

  # variable names must be valid identifiers
  [[ "$root_var"   =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ ]] || { printf '[ERROR] invalid var: %s\n' "$root_var" >&2; return 1; }
  [[ "$loaded_var" =~ ^[a-zA-Z_][a-zA-Z0-9_]*$ ]] || { printf '[ERROR] invalid var: %s\n' "$loaded_var" >&2; return 1; }

  # require root to be already set
  [[ -n "${!root_var:-}" ]] || { printf '[ERROR] %s is unset\n' "$root_var" >&2; return 1; }

  # simple idempotence check (1 == already loaded)
  if [[ "${!loaded_var:-0}" == 1 ]]; then
    return 0
  fi

  init="${!root_var}/${init_rel}"
  [[ -r "$init" ]] || { printf '[ERROR] Missing init: %s\n' "$init" >&2; return 2; }

  # shellcheck source=/dev/null
  . "$init" || { printf '[ERROR] Failed to load: %s\n' "$init" >&2; return 3; }

  printf -v "$loaded_var" '%d' 1
  export "$loaded_var"
}
#EOC

#-----------------------------------------------------------------------------#
#BOP
# !FUNCTION: _parse_cycle_args
# !DESCRIPTION:
#   Parse cycle-specific CLI options. Intended to be called after the global
#   parser (__parse_args__), consuming only the leftover arguments passed
#   down from it. Sets defaults, applies user overrides, and updates global
#   shell variables used by the run cycle driver.
#
# !USAGE:
#   _parse_cycle_args "${leftover_args[@]}" || exit $?
#
# !OPTIONS:
#   -t   <INT>    Model spectral truncation.
#   -l   <INT>    Number of vertical levels in the model.
#   -p   <STR>    Model run prefix (e.g., SMT, CPT).
#   -mnp <INT>    Number of MPI tasks for the model (default: 64).
#
#   -I   <YYYYMMDDHH>  Initial cycle datetime (required).
#   -F   <YYYYMMDDHH>  Final   cycle datetime (required).
#   -fct <INT>         Forecast length in hours (default: 09).
#
#   -gt  <INT>    GSI spectral truncation.
#   -gnp <INT>    Number of MPI tasks for GSI (default: 144).
#
#   -bc  <INT>    Number of bias-correction spin-up cycles (default: 0).
#
#   --no-obsmake  Disable OBSMAKE stage (default: enabled).
#   --no-gsi      Disable GSI stage (default: enabled).
#   --no-bam      Disable BAM stage (default: enabled).
#
#   -h, --help    Show usage help and exit 0.
#   --            End of options; remaining args are positional.
#
# !BEHAVIOR:
#   - Sets shell variables: modelTrunc, modelNLevs, modelPrefix,
#     modelMPITasks, modelFCT, gsiTrunc, gsiMPITasks,
#     LABELI, LABELF, BcCycles, do_obsmake, do_gsi, do_bam.
#   - Applies defaults if not specified on the command line.
#   - Emits a warning and ignores unrecognized options.
#
# !NOTES:
#   - Requires Bash 4+ for arrays and indirect expansion.
#   - All variables are global (intended for use by _main driver).
#   - Must be called after __parse_args__, consuming leftover_args.
#   - Returns non-zero if required args are missing.
#EOP
_parse_cycle_args() {

  #------------------------------ Defaults -----------------------------------#
  # These assignments only apply if variables are unset or empty.
  # This allows overriding defaults via environment variables.

  : ${modelMPITasks:=64}   # default MPI tasks for model stage
  : ${modelFCT:=09}        # default forecast length in hours
  : ${gsiMPITasks:=144}    # default MPI tasks for GSI stage

  : ${do_obsmake:=1}       # run OBSMAKE stage by default
  : ${do_gsi:=1}           # run GSI stage by default
  : ${do_bam:=1}           # run BAM stage by default

  : ${modelTrunc:=}        # model spectral truncation (user-specified)
  : ${modelNLevs:=}        # number of vertical levels (user-specified)
  : ${modelPrefix:=}       # model run prefix (e.g. SMT, CPT)
  : ${gsiTrunc:=}          # GSI truncation (user-specified)

  : ${LABELI:=}            # initial cycle date (YYYYMMDDHH, required)
  : ${LABELF:=}            # final   cycle date (YYYYMMDDHH, required)

  : ${BcCycles:=0}         # bias-correction spin-up cycles (default: none)

  #------------------------------ Parse loop ---------------------------------#
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
      --no-gsi)     do_gsi=0;     shift;;
      --no-bam)     do_bam=0;     shift;;

      -h|--help)    usage; exit 0;;
      --)           shift; break;;
      -*)           _log_warn "Unknown option ignored (cycle): %s" "$1"; shift;;
      *)            break;;
    esac
  done

  # Any leftover args after "--" remain for positional handling
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
  . "${scripts_smg}/run_obsmake.sh" "${lbl}" "${obsmake_cli[@]}"
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
#     options (-I, -T, -t, -l, -p, --gsi-ntasks, -bc). If OMP threads are provided by
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

  # --- Resolve effective MPI task count for GSI (--gsi-ntasks) ---
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
  # -l modelNLevs, -p prefix, --gsi-ntasks MPI, -bc BC index (0 for normal mode).
  local argv=(
    -I "${lbl}"
    -T "${gsiTrunc}"
    -t "${modelTrunc}"
    -l "${modelNLevs}"
    -p "${modelPrefix}"
    -bc "${bc_cycle}"
  )
  [[ -n "${gnp:-}" ]] && argv+=(--ntasks "${gnp}")

  # --- Invoke runGSI (no wrapper; let runGSI handle submission/logging) ---
  . "${scripts_smg}/runGSI" "${argv[@]}"
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
  . "${scripts_smg}/run_model.sh" "${lbl}" "${fct_date}" "${modelPrefix}" "${modelTrunc}" "${modelNLevs}" "${modelMPITasks}" "No"
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
  local lbl="${1:?Missing LABEL}"   # required label argument
  local -i n="${2:-10}"             # number of iterations, default = 10
  local -i i rc                     # loop index and return code (integers, local)

  # Log the start of the spin-up process
  _log_info "Bias-correction spin-up: %d iterations @ %s" "$n" "$lbl"

  # Loop over the specified number of iterations
  for (( i=1; i<=n; i++ )); do
    _log_info "Spin-up iteration %d/%d" "$i" "$n"

    # Run GSI for this iteration
    _run_gsi "$lbl" "$i"
    rc=$?   # capture the exit status immediately

    # Report the result of this iteration
    _log_info "Iteration %d/%d finished with rc=%d" "$i" "$n" "$rc"

    # Optional: abort if any iteration fails
    if (( rc != 0 )); then
      _log_err "GSI failed on iteration %d (rc=%d)" "$i" "$rc"
      return "$rc"
    fi
  done

  # Log the completion of the spin-up procedure
  _log_ok "Bias-correction spin-up completed"
}


#EOC

#-----------------------------------------------------------------------------#
#BOP
# !FUNCTION: _main
# !DESCRIPTION:
#   Orchestrates a full SMNA cycle driver:
#     1) Bootstraps root + environment init (idempotent),
#     2) Parses global + local CLI flags (with leftovers hand-off),
#     3) Validates required options and derives defaults,
#     4) Optionally performs a one-time bias-correction spin-up,
#     5) Iterates production cycles from LABELI to LABELF (step +6h),
#        invoking OBSMAKE → GSI → BAM sub-routines as requested.
# !USAGE:
#   _main "$@"
# !INPUT (via CLI / env; validated by parsers):
#   LABELI          Initial cycle label (YYYYMMDDHH).
#   LABELF          Final   cycle label (YYYYMMDDHH).
#   do_obsmake      1/0 toggle to run OBSMAKE stage.
#   do_gsi          1/0 toggle to run GSI stage.
#   do_bam          1/0 toggle to run BAM stage.
#   BcCycles        Non-negative integer; if >0, run bias spin-up once.
#   modelTrunc      Spectral truncation for model.
#   modelNLevs      Number of vertical levels for model.
#   modelPrefix     Run prefix (e.g., SMT/CPT/...).
#   modelMPITasks   MPI tasks for model stage.
#   modelFCT        Forecast length (hours).
#   gsiTrunc        Truncation for GSI grid (if applicable).
#   gsiMPITasks     MPI tasks for GSI stage.
# !DEPENDENCIES:
#   Functions: ensure_root, source_init_once, mini_sanity,
#              __parse_args__, _parse_cycle_args, _validate_required,
#              _spinup_bias_once, _run_obsmake, _run_gsi, _run_bam,
#              _log_info, _die (used indirectly by sanity/validators).
#   External:  date(1) with GNU -d support (used to add +6 hours in UTC).
# !BEHAVIOR:
#   - Global parser may set `leftover_args` (array) consumed by local parser.
#   - Cycle loop increments `current` in +6h steps, inclusive of LABELF if equal.
#   - Each stage is conditionally executed via arithmetic context (( do_* )).
# !RETURNS:
#   0 on success; non-zero on error (propagated from called helpers/parsers).
# !NOTES:
#   - Assumes environment was not previously initialized; if it was, init is
#     idempotent thanks to source_init_once’s LOADED flag.
#   - Ensure TZ-sensitive behavior is explicit: date is invoked with -u.
#EOP
_main() {

  # Local read of verbosity (do not force global default here)
  local verbose=${verbose:-true}

  # --- 0) Bootstrap: project root & one-time init ---------------------------#
  # Ensures SMG_ROOT points to the repository root (via marker) and exports it.
  ensure_root SMG_ROOT || exit $?
  # Sources $SMG_ROOT/etc/__init__.sh exactly once (guards with INIT_LOADED=1).
  source_init_once SMG_ROOT INIT_LOADED

  # --- 0.1) Minimal prerequisites ------------------------------------------#
  # Checks for required commands and existence of scripts root.
  command -v date >/dev/null 2>&1 || _die 3 "'date' not found in PATH"

  # Run scripts directory
  : "${scripts_smg:=${SMG_ROOT}/run/scripts}"
  [[ -d "${scripts_smg}" ]] || _die 4 "scripts_smg directory not found: %s" "${scripts_smg}"

  # --- 1) Global CLI parsing (common flags, HPC knobs) ----------------------#
  # Parses top-level options (e.g., --debug, --dry-run, queueing params, etc.).
  # On success, may populate `leftover_args` with unconsumed positionals/options.
  __parse_args__ "$@" || exit $?

  # --- 2) Local (cycle) CLI parsing ----------------------------------------#
  # Parse only the cycle-specific flags/positionals using leftovers from global.
  # This keeps parsers decoupled and composable.
  _parse_cycle_args "${leftover_args[@]}" || exit $?

  # --- 3) Validate required inputs and derive defaults ----------------------#
  # Ensures LABELI/LABELF and stage toggles are consistent; computes missing
  # defaults (e.g., truncation, levels, MPI tasks) if not provided.
  _validate_required

  # --- 3.1) Informational summary ------------------------------------------#
  _log_info "Cycle window: %s → %s (step +6h)" "${LABELI}" "${LABELF}"
  _log_info "Model: trunc=%s nlevs=%s prefix=%s mnp=%s fct=%sh" \
            "${modelTrunc}" "${modelNLevs}" "${modelPrefix}" \
            "${modelMPITasks}" "${modelFCT}"
  _log_info "GSI:   trunc=%s gnp=%s" "${gsiTrunc}" "${gsiMPITasks}"

  # --- 4) Optional bias-correction spin-up ----------------------------------#
  # If BcCycles>0, prepare VarBC coefficients once before the production loop.
  if (( BcCycles > 0 )); then
    _run_obsmake "${LABELI}"
    _spinup_bias_once "${LABELI}" "${BcCycles}"
  fi

  # --- 5) Production loop: run selected stages every +6h --------------------#
  local current="${LABELI}"
  while [[ "${current}" -le "${LABELF}" ]]; do
    _log_info "Submitting system for date %s" "${current}"

    # Stage toggles use arithmetic context: 0 (skip) / non-zero (run).
    (( do_obsmake )) && _run_obsmake "${current}"
    (( do_gsi     )) && _run_gsi     "${current}"
    (( do_bam     )) && _run_bam     "${current}"

    # Advance in UTC by +6 hours, preserving YYYYMMDDHH formatting.
    current="$(date -u +%Y%m%d%H -d "${current:0:8} ${current:8:2} +6 hours")"
  done

  # --- 6) Epilogue ----------------------------------------------------------#
  _log_info "run_cycle completed successfully"
}
watch_i() {
  echo "[DEBUG] i=$i (func=${FUNCNAME[1]} line=${BASH_LINENO[0]})"
}
trap watch_i DEBUG
_main "$@"
#EOC
#-----------------------------------------------------------------------------#

