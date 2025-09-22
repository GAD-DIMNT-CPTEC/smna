#!/usr/bin/env bash
#-----------------------------------------------------------------------------#
#           Group on Data Assimilation Development - GDAD/CPTEC/INPE          #
#-----------------------------------------------------------------------------#
#BOP
# !SCRIPT: runGSI_functions.sh
#
# !DESCRIPTION:
#   Helper library sourced by runGSI. Provides logging helpers, fixed-files
#   preparation, observation linking, bias files handling, submission wrappers
#   for GSI and angle-dependent bias correction, and various utility routines.
#
# !USAGE:
#   # Always source this file before calling its functions:
#   #   source /path/to/runGSI_functions.sh
#   #
#   # Minimal workflow (example):
#   #   constants                 # set HPC defaults, files, paths
#   #   ParseOpts "$@"            # parse CLI options (called by runGSI)
#   #   linkObs "$AnlDate" "$runDir"
#   #   FixedFiles "$AnlTrunc" "$AnlNLevs" "$runDir"
#   #   subGSI "$runDir" "$IMax" "$JMax" "$KMax" "$AnlTrunc" "$BkgTrunc" \
#   #          "$AnlDate" "$oneobtest" "$is_cold"
#
# !NOTES:
#   • This file is meant to be *sourced*; it must not change global shell flags.
#     Functions that require strict mode enable it locally and restore on return.
#   • All log helpers respect the global boolean 'verbose' (default: false).
#     Errors (_log_err/_log_fail) always print, regardless of verbosity.
#   • Expected boolean convention: "true" / "false" (lowercase).
#
# !REQUIREMENTS:
#   bash >= 4.0, coreutils, awk, sed, find, sort, grep, printf
#
# !REVISION HISTORY:
#   2018-08-20 - J. G. de Mattos - Initial version
#   2025-09-19 - Refactor: safer logging, local strict-mode, getopts parsing,
#                path discovery, CRTM/coeff search robust to nested trees.
#EOP
#-----------------------------------------------------------------------------#
#BOC

# --- Default logging verbosity (string boolean) ---
: "${verbose:=false}"

# --- per-rank post-merge action (PEs) (move, delete, keep)
: "${PES_ACTION:=move}"

# --- Staging mode (copy, symlink, hardlink)
: "${STAGE_LINK_MODE:=copy}"

# --- Library directory (folder where this file lives) ---
# Make the helper self-aware even when sourced from anywhere.
# Works for both symlinked and real paths.
__RGF_SRC="${BASH_SOURCE[0]}"
while [[ -L "$__RGF_SRC" ]]; do
  __RGF_LINK_TARGET="$(readlink -- "$__RGF_SRC")"
  if [[ "$__RGF_LINK_TARGET" = /* ]]; then
    __RGF_SRC="$__RGF_LINK_TARGET"
  else
    __RGF_SRC="$(cd -- "$(dirname -- "$__RGF_SRC")" && cd -- "$(dirname -- "$__RGF_LINK_TARGET")" && pwd)/$(basename -- "$__RGF_LINK_TARGET")"
  fi
done
export RUN_GSI_FUNCS_DIR="$(cd -- "$(dirname -- "$__RGF_SRC")" && pwd)"
#EOC

#-----------------------------------------------------------------------------#
#------------------------------ HELPERS LOADER --------------------------------#
#-----------------------------------------------------------------------------#
# -----------------------------------------------------------------------------#
#BOP
# !FUNCTION: source_nearby
# !INTERFACE: source_nearby <filename> [start_dir]
# !DESCRIPTION:
#   Search for a target file (e.g., libs.sh) starting from the given directory
#   (or from the caller script location if not provided), climbing up toward root,
#   and at each level performing a pruned recursive search. Honors LIBS_PATH env
#   as a fast-path hint. Once found, the file is sourced into the current shell.
#
# !USAGE:
#   source_nearby "libs.sh"
#   source_nearby "config.sh" /opt/myproj/run/scripts
#
# !ENVIRONMENT:
#   LIBS_PATH                Optional absolute path to the file or a directory
#                            containing it. If valid, takes precedence.
#   SOURCE_NEARBY_VERBOSE    true|false (default: true). Controls [INFO]/[WARN] output.
#   SOURCE_NEARBY_MAXDEPTH   Max depth per directory level during search (default: 4).
#   SOURCE_NEARBY_GLOBAL_ROOT  Path for fallback global search (default: "/").
#   SOURCE_NEARBY_PRUNE      Colon-separated list of dirs to skip in search (added
#                            to defaults like .git, build, node_modules, etc.).
#
# !BEHAVIOR:
#   • If LIBS_PATH is set and valid, sources directly from it.
#   • Otherwise, resolves start_dir (or script dir) and attempts:
#       1. Same directory as the script.
#       2. Iteratively climb parent dirs, in each doing a pruned recursive search.
#       3. As last fallback, global search under SOURCE_NEARBY_GLOBAL_ROOT.
#   • Always silences "Permission denied" errors during find.
#
# !RETURNS:
#   0 on success (file sourced),
#   1 if target not found,
#   >1 for unexpected internal errors.
#
# !NOTES:
#   • Safe under 'set -euo pipefail'.
#   • Uses 'find -xdev' to avoid crossing filesystem mounts.
#   • Recommended: export LIBS_PATH once found to speed up subsequent calls.
#EOP
#BOC
source_nearby() {
  local target="${1:?filename required}"
  local start="${2:-}"

  # --------------------- logging control ---------------------
  local _verbose="${SOURCE_NEARBY_VERBOSE:-true}"
  _log() { "$_verbose" && printf '%s\n' "$*" >&2 || true; }

  # ----------------- helper: resolve script dir ---------------
  _script_dir() {
    # Use BASH_SOURCE when available (robust for sourced scripts)
    local src="${BASH_SOURCE[0]:-}"
    if [[ -n "$src" && -e "$src" ]]; then
      cd -- "$(dirname -- "$src")" && pwd -P
    else
      pwd -P
    fi
  }

  # --------------------- resolve start dir -------------------
  local script_dir
  if [[ -n "$start" ]]; then
    script_dir="$(cd -- "$start" && pwd -P)"
  else
    script_dir="$(_script_dir)"
  fi

  # -------------------- explicit path case -------------------
  if [[ "$target" == */* && -f "$target" ]]; then
    # shellcheck disable=SC1090
    source "$target"
    _log "[INFO] Loaded (explicit path): $target"
    return 0
  fi

  # -------------------- 0) LIBS_PATH hint --------------------
  if [[ -n "${LIBS_PATH:-}" ]]; then
    if [[ -f "$LIBS_PATH" ]]; then
      # shellcheck disable=SC1090
      source "$LIBS_PATH"
      _log "[INFO] Loaded via LIBS_PATH (file): $LIBS_PATH"
      return 0
    elif [[ -d "$LIBS_PATH" && -f "$LIBS_PATH/$target" ]]; then
      # shellcheck disable=SC1090
      source "$LIBS_PATH/$target"
      _log "[INFO] Loaded via LIBS_PATH (dir): $LIBS_PATH/$target"
      return 0
    fi
    _log "[WARN] LIBS_PATH set but not valid for '$target': $LIBS_PATH"
  fi

  # ----------------- search parameters -----------------------
  local maxdepth="${SOURCE_NEARBY_MAXDEPTH:-4}"
  local global_root="${SOURCE_NEARBY_GLOBAL_ROOT:-/}"

  # Default prune dirs
  local -a _prune_default=(.git .svn .hg build dist node_modules __pycache__ .mypy_cache .venv venv .tox .cache .idea .vscode)
  IFS=':' read -r -a _prune_extra <<< "${SOURCE_NEARBY_PRUNE:-}"
  local -a PRUNE_DIRS=("${_prune_default[@]}" "${_prune_extra[@]}")

  # Helper: pruned find in base dir
  _find_here() {
    local base="$1"
    local -a prune_expr=()
    for d in "${PRUNE_DIRS[@]}"; do
      [[ -n "$d" ]] && prune_expr+=(-name "$d" -o)
    done
    ((${#prune_expr[@]})) && unset 'prune_expr[${#prune_expr[@]}-1]'

    # Use -xdev to avoid crossing mounts; stop at first match (-quit)
    find "$base" -xdev \
      \( -type d \( "${prune_expr[@]}" \) -prune \) -o \
      \( -type f -name "$target" -maxdepth "$maxdepth" -print -quit \) 2>/dev/null
  }

  # ----------------- 1) same directory -----------------------
  if [[ -f "$script_dir/$target" ]]; then
    # shellcheck disable=SC1091
    source "$script_dir/$target"
    _log "[INFO] Loaded: $script_dir/$target"
    return 0
  fi

  # ----------------- 2) climb parent dirs --------------------
  local cur="$script_dir"
  local hit=""
  while :; do
    hit="$(_find_here "$cur")" || true
    if [[ -n "$hit" && -f "$hit" ]]; then
      # shellcheck disable=SC1090
      source "$hit"
      _log "[INFO] Loaded: $hit"
      return 0
    fi
    [[ "$cur" == "/" ]] && break
    cur="$(dirname -- "$cur")"
  done

  # ----------------- 3) global fallback ----------------------
  hit="$(_find_here "$global_root")" || true
  if [[ -n "$hit" && -f "$hit" ]]; then
    # shellcheck disable=SC1090
    source "$hit"
    _log "[INFO] Loaded (global): $hit"
    return 0
  fi

  printf '[ERROR] Could not find "%s" from "%s" or under "%s".\n' \
    "$target" "$script_dir" "$global_root" >&2
  return 1
}
#EOP


source_nearby "__helpers__.sh"


#EOC

#-----------------------------------------------------------------------------#
#------------------------------- SMALL UTILS ---------------------------------#
#-----------------------------------------------------------------------------#


#BOP
# !FUNCTION: usage
# !DESCRIPTION:
#   Print usage banner extracted from the main driver script ($runGSI).
#EOP
#BOC
usage() {
  echo
  echo "Usage:"
  if [[ -n "${runGSI:-}" && -r "$runGSI" ]]; then
    sed -n '/^#BOP/,/^#EOP/{/^#BOP/d;/^#EOP/d;p}' -- "$runGSI"
  else
    echo "runGSI not set; cannot extract usage block."
  fi
}
#EOC
#BOP
# !FUNCTION: StageBackgrounds
# !INTERFACE: StageBackgrounds RUN_DIR BKG_PREFIX BKG_DATE BKG_MRES ANL_DATE SUBT_MODEL_BAM
# !DESCRIPTION:
#   Stage (copy or link) BAM background forecasts (+03,+06,+09) **and** their
#   descriptor files into RUN_DIR with canonical names `BAM.fct.{HH}` and
#   `BAM.dir.{HH}` for GSI. Fails fast on missing inputs.
#
#   The staging mode can be customized via the environment variable
#   STAGE_LINK_MODE:
#     - "copy"    (default) copy files (preserving perms/timestamps)
#     - "symlink" create symlinks
#     - "hardlink" create hardlinks (same filesystem required)
#
#   If a helper `_copy_with_progress SRC DST` is defined in the environment,
#   it will be used automatically for copy operations (otherwise `cp -pf`).
#
# !USAGE:
#   StageBackgrounds "${RunDir}" "${BkgPrefix}" "${BkgDate}" "${BkgMRES}" "${AnlDate}" "${subt_model_bam}"
#
# !ARGUMENTS:
#   RUN_DIR         Target directory where files will be staged
#   BKG_PREFIX      Background prefix used by BAM outputs (e.g., "PREF")
#   BKG_DATE        Cycle base date of the background (yyyymmddhh)
#   BKG_MRES        BAM resolution label (e.g., "TQ0299L064")
#   ANL_DATE        Analysis date (yyyymmddhh) that defines -3/0/+3 hours window
#   SUBT_MODEL_BAM  BAM root path that contains dataout/<mres>/DAS/<BKG_DATE>/*
#
# !RETURNS:
#   0 on success; non-zero on first missing input or staging failure.
#
# !NOTES:
#   - Requires GNU date.
#   - Creates RUN_DIR if it does not exist.
#   - Overwrites existing staged targets.
#EOP
#BOC
StageBackgrounds() {
  _with_strict_mode   # enable strict mode only for this function

  # ------------------------ Args & basic checks ------------------------------
  local run_dir=${1:? "Missing RUN_DIR"}
  local bkg_prefix=${2:? "Missing BKG_PREFIX"}
  local bkg_date=${3:? "Missing BKG_DATE (yyyymmddhh)"}
  local bkg_mres=${4:? "Missing BKG_MRES"}
  local anl_date=${5:? "Missing ANL_DATE (yyyymmddhh)"}
  local subt_bam=${6:? "Missing SUBT_MODEL_BAM"}

  # Validate date formats (basic 10-digit check)
  [[ "${bkg_date}" =~ ^[0-9]{10}$ ]] || { _log_err "BKG_DATE must be yyyymmddhh"; return 2; }
  [[ "${anl_date}" =~ ^[0-9]{10}$ ]] || { _log_err "ANL_DATE must be yyyymmddhh"; return 2; }

  # Staging mode
  local mode="${STAGE_LINK_MODE:-copy}"
  case "${mode}" in
    copy|symlink|hardlink) : ;;
    *) _log_err "STAGE_LINK_MODE must be one of: copy|symlink|hardlink"; return 2 ;;
  esac

  # Source directory that contains BAM outputs for this background date
  local src_base="${subt_bam}/dataout/${bkg_mres}/DAS/${bkg_date}"

  # Ensure target dir exists and is writable
  mkdir -p -- "${run_dir}"
  [[ -w "${run_dir}" ]] || { _log_err "RUN_DIR is not writable: ${run_dir}"; return 3; }

  _need() {
    # _need FILEPATH  → check existence & regular file
    local f="$1"
    [[ -f "${f}" ]] || { _log_err "Missing file: ${f}"; return 1; }
  }

  # --------------------------- Main staging loop -----------------------------
  local inc Date HH FFCT FDIR src_fct src_dir
  for inc in -3 0 3; do
    # Hour bucket names: -3→03, 0→06, +3→09 (relative to analysis)
    HH=$(printf "%02g" $((inc + 6)))
    Date=$(date -u +%Y%m%d%H -d "${anl_date:0:8} ${anl_date:8:2} ${inc} hours")

    FFCT="GFCT${bkg_prefix}${bkg_date}${Date}F.fct.${bkg_mres}"
    FDIR="GFCT${bkg_prefix}${bkg_date}${Date}F.dir.${bkg_mres}"

    src_fct="${src_base}/${FFCT}"
    src_dir="${src_base}/${FDIR}"

    _log_info -f " • Staging HH=${HH}  (Date=${Date})"
    _log_info -f "   - Forecast : ${src_fct}"
    _log_info -f "   - Descriptor: ${src_dir}"

    _need "${src_fct}" || return 4
    _need "${src_dir}" || return 4

    _copy_one_safe "${src_fct}" "${run_dir}/BAM.fct.${HH}" "${mode}"
    _copy_one_safe "${src_dir}" "${run_dir}/BAM.dir.${HH}" "${mode}"
  done

  _log_ok "Backgrounds staged into ${run_dir} (mode=${mode})."
  return 0
}
#EOC
#BOP
# !FUNCTION: CopyGSIExec
# !INTERFACE: CopyGSIExec EXEC_PATH RUN_DIR
# !DESCRIPTION:
#   Stage the GSI executable into RUN_DIR and ensure it is executable.
#   By default it copies the binary preserving timestamps/permissions.
#
#   Behavior can be customized via environment variables:
#     - STAGE_LINK_MODE
#         "copy"    (default) → copy file
#         "symlink"           → create/update a symlink
#         "hardlink"          → create/update a hardlink (same filesystem)
#     - GSI_TARGET_NAME
#         If set, use this filename for the staged binary (e.g., "gsi.x").
#         Otherwise, use the basename of EXEC_PATH.
#     - VERIFY_COPY
#         If set to "1" and STAGE_LINK_MODE=copy, verify the copy via sha256.
#     - If a helper `_copy_with_progress SRC DST` exists, it is used for copies.
#
# !USAGE:
#   CopyGSIExec "/path/to/gsi.x" "${RunDir}"
#
# !RETURNS:
#   0 on success;
#   1 if EXEC_PATH is missing or not a regular file;
#   2 if RUN_DIR is not writable or cannot be created;
#   3 on staging failure (copy/link) or verification mismatch.
#
# !NOTES:
#   - Overwrites existing target.
#   - Sets executable bit on the staged file (for copy/hardlink modes).
#EOP
#BOC
CopyGSIExec() {
  _with_strict_mode   # enable strict mode only for this function

  local exec_src=${1:? "Missing EXEC_PATH"}
  local run_dir=${2:? "Missing RUN_DIR"}

  # Validate source
  if [[ ! -f "${exec_src}" ]]; then
    _log_err "Executable not found: %s\n" "${exec_src}"
    return 1
  fi

  # Ensure target dir
  mkdir -p -- "${run_dir}" || {
    _log_err "Cannot create RUN_DIR: %s\n" "${run_dir}"
    return 2
  }
  [[ -w "${run_dir}" ]] || {
    _log_err "RUN_DIR not writable: %s\n" "${run_dir}"
    return 2
  }

  # Resolve target name/path
  local target_name="${GSI_TARGET_NAME:-$(basename -- "${exec_src}")}"
  local exec_dst="${run_dir%/}/${target_name}"

  # Select staging mode
  local mode="${STAGE_LINK_MODE:-copy}"

  _copy_one_safe "${exec_src}" "${exec_dst}" "${mode}"

  if [[ "${mode}" == "copy" ]]; then
    chmod u+x,go+rx -- "${exec_dst}" || true

    # Optional verification
    if [[ "${VERIFY_COPY:-0}" == "1" ]]; then
      local src_hash dst_hash
      src_hash=$(sha256sum "${exec_src}" | awk '{print $1}') || true
      dst_hash=$(sha256sum "${exec_dst}" | awk '{print $1}') || true
      if [[ -z "${src_hash}" || -z "${dst_hash}" || "${src_hash}" != "${dst_hash}" ]]; then
        _log_err "sha256 mismatch after copy\n  src=%s\n  dst=%s\n" "${src_hash:-<nil>}" "${dst_hash:-<nil>}"
        return 3
      fi
    fi
  fi
  _log_ok "GSI executable staged to %s (mode=%s)" "${exec_dst}" "${mode}"
  return 0
}
#EOC
#BOP
# !FUNCTION: CopyOutputsForCycle
# !INTERFACE: CopyOutputsForCycle RUN_DIR DATAOUT_DIR IS_BC BC_CYCLE
# !DESCRIPTION:
#   Handle the post-processing copy stage after a successful GSI run.
#
#   This function moves or copies relevant GSI outputs from the working
#   run directory (RUN_DIR) to the final data output location (DATAOUT_DIR).
#
#   Behavior depends on whether Bias Correction (BC) cycling mode is active:
#
#   • IS_BC=true
#       - Outputs are copied into a cycle-specific directory:
#           ${DATAOUT_DIR}/BC.C{cycle+1}
#       - If BC_CYCLE > 0, the function also retrieves the previous cycle’s
#         satbias files:
#           satbias*.out → renamed as satbias*.in
#         These are injected back into the current RUN_DIR so that GSI can
#         continue bias coefficient cycling across iterations.
#       - Supports optional angle-dependent bias correction files if
#         satbiasAngOu/satbiasAngIn are defined.
#
#   • IS_BC=false
#       - Outputs are copied directly into DATAOUT_DIR with no cycle nesting
#         or bias-correction file handoff.
#
#   Required environment variables:
#       satbiasIn     → filename for bias input coefficients
#       satbiasPCIn   → filename for preconditioning bias input
#       satbiasOu     → filename for bias output coefficients
#       satbiasPCOu   → filename for preconditioning bias output
#
#   Optional environment variables:
#       satbiasAngIn  → filename for angle-dependent bias input
#       satbiasAngOu  → filename for angle-dependent bias output
#
# !USAGE:
#   CopyOutputsForCycle /path/to/run_dir /path/to/dataout true  2
#   CopyOutputsForCycle /path/to/run_dir /path/to/dataout false 0
#
# !NOTES:
#   - Expects helper function `copyFiles SRC DST` to perform the actual copy.
#   - BC cycle numbering starts at 0 internally, but outputs are stored in
#     BC.C01, BC.C02, … directories (cycle+1).
#   - Use `bc_cycle=0` to indicate the first iteration (no back-propagation
#     of satbias*.out yet).
#EOP
#BOC
CopyOutputsForCycle() {
  _with_strict_mode   # enable strict mode only for this function

  local run_dir=${1:? "Missing RUN_DIR"}
  local dataout=${2:? "Missing DATAOUT_DIR"}
  local is_bc=${3:-false}
  local bc_cycle=${4:-0}

  if ${is_bc}; then
    local next=$((bc_cycle+1))
    local bc_dir_out
    printf -v bc_dir_out "%s/BC.C%02d" "${dataout}" "${next}"
    mkdir -p -- "${bc_dir_out}"
    copyFiles "${run_dir}" "${bc_dir_out}"

    # Bring satbias*.out from previous cycle as input for next one (if not first)
    if (( bc_cycle > 0 )); then
      local prev_dir
      printf -v prev_dir "%s/BC.C%02d" "${dataout}" "${bc_cycle}"
      [[ -f "${prev_dir}/${satbiasOu}"    ]] && cp -pfr "${prev_dir}/${satbiasOu}"    "${run_dir}/${satbiasIn}"
      [[ -f "${prev_dir}/${satbiasPCOu}"  ]] && cp -pfr "${prev_dir}/${satbiasPCOu}"  "${run_dir}/${satbiasPCIn}"
      # angle bias optional:
      if [[ -n "${satbiasAngOu:-}" && -n "${satbiasAngIn:-}" && -f "${prev_dir}/${satbiasAngOu}" ]]; then
        cp -pfr "${prev_dir}/${satbiasAngOu}" "${run_dir}/${satbiasAngIn}"
      fi
    fi
  else
    copyFiles "${run_dir}" "${dataout}"
  fi
}
#EOC

#-----------------------------------------------------------------------------#
#------------------------------- CONSTANTS -----------------------------------#
#-----------------------------------------------------------------------------#
#BOP
# !FUNCTION: constants
# !INTERFACE: constants
# !DESCRIPTION:
#   Initialize cluster defaults, byte order, and fixed-file paths used by runGSI.
#   When invoked from a call stack that includes a 'main' function, this helper
#   also resolves and exports the absolute path of the caller script (runGSI) and
#   its directory. Cluster-specific settings are applied based on $hpc_name.
#
# !USAGE:
#   constants
#
# !BEHAVIOR:
#   • Detects runGSI path:
#       - If 'main' is present in the call stack, exports:
#           runGSI   = real path to the caller script
#           runGSIDir= directory containing the caller script
#   • Sets default BYTE_ORDER if unset: Big_Endian
#   • Applies cluster presets according to $hpc_name:
#       egeon:
#         MaxCoresPerNode=60, MTasks (default 120), ThreadsPerMPITask=1,
#         TasksPerNode=MaxCoresPerNode/ThreadsPerMPITask,
#         PEs=MTasks/ThreadsPerMPITask,
#         Nodes=(MTasks+MaxCoresPerNode-1)/MaxCoresPerNode,
#         Queue=batch, WallTime=01:00:00, BcCycles=0,
#         MPICH_* tuning variables exported
#       XC50/xc50:
#         MaxCoresPerNode=40, MPITasks (default 120), ThreadsPerMPITask=1,
#         TasksPerNode=MaxCoresPerNode/ThreadsPerMPITask,
#         PEs=MPITasks/ThreadsPerMPITask,
#         Nodes=(MPITasks+MaxCoresPerNode-1)/MaxCoresPerNode,
#         Queue=pesq, WallTime=01:00:00, BcCycles=0,
#         MPICH_* tuning variables exported
#       other:
#         Emits a WARNING and leaves conservative defaults
#   • Exports executable and fixfile paths if not already set:
#       execGSI, parmGSI, obsGSI, cldRadInfo, SatBiasSample, SatBiasPCSample,
#       ScanInfo, SatBiasAngSample, execBCAng, parmBCAng
#   • Exports canonical satbias filenames used by other helpers:
#       satbiasIn, satbiasOu, satbiasPCIn, satbiasPCOu, satbiasAngIn, satbiasAngOu
#
# !ENVIRONMENT:
#   Inputs (optional; defaults applied if unset):
#     hpc_name            Cluster selector: egeon | XC50 (case-insensitive) | other
#     BYTE_ORDER          Big_Endian | Little_Endian (default: Big_Endian)
#     home_cptec          Base path for executables (used for execGSI, execBCAng)
#     home_gsi_fix        Base path for GSI fix files (parmGSI, obsGSI, etc.)
#     public_fix          Base path for public fix assets (SatBias* samples)
#     MTasks / MPITasks   Total MPI ranks (cluster-specific default = 120)
#     ThreadsPerMPITask   OpenMP threads per MPI rank (default = 1)
#     Queue, WallTime     Scheduler defaults (cluster-specific)
#   Outputs (exported):
#     runGSI, runGSIDir, BYTE_ORDER, cluster layout vars, MPICH_* tunables,
#     execGSI, parmGSI, obsGSI, cldRadInfo, SatBiasSample, SatBiasPCSample,
#     ScanInfo, SatBiasAngSample, execBCAng, parmBCAng,
#     satbiasIn, satbiasOu, satbiasPCIn, satbiasPCOu, satbiasAngIn, satbiasAngOu
#
# !RETURNS:
#   0 on completion. Emits a WARNING if hpc_name is unknown; no hard failure.
#
# !NOTES:
#   • Path resolution uses GNU readlink -e; on systems without -e, consider a
#     portability shim or 'realpath' fallback.
#   • MTasks vs MPITasks:
#       - egeon branch uses MTasks to compute PEs/Nodes
#       - XC50 branch uses MPITasks
#     Ensure the appropriate variable is set for your cluster.
#   • This helper is side-effectful by design (exports many variables).
#     Call it early in your runGSI workflow.
#EOP
#BOC
constants() {
  _with_strict_mode   # enable strict mode only for this function

  # detect runGSI path from the call stack if available
  for (( i=0; i<${#FUNCNAME[@]}; i++ )); do
    if [[ "${FUNCNAME[$i]}" == "main" ]]; then
      export runGSIDir
      runGSIDir="$(cd -- "$(dirname -- "$(readlink -e -- "${BASH_SOURCE[$i]}")")" && pwd)"
      export runGSI
      runGSI="$(readlink -e -- "${BASH_SOURCE[$i]}")"
      break
    fi
  done

  # Endianness (default)
  export BYTE_ORDER="${BYTE_ORDER:-Big_Endian}"

  case "${hpc_name:-unknown}" in
    egeon)
      export MaxCoresPerNode=60
      export MTasks="${MTasks:-120}"
      export ThreadsPerMPITask="${ThreadsPerMPITask:-1}"
      export TasksPerNode=$(( MaxCoresPerNode / ThreadsPerMPITask ))
      export PEs=$(( MTasks / ThreadsPerMPITask ))
      export Nodes=$(( (MTasks + MaxCoresPerNode - 1) / MaxCoresPerNode ))
      export Queue="${Queue:-batch}"
      export WallTime="${WallTime:-01:00:00}"
      export BcCycles="${BcCycles:-0}"
      export MPICH_UNEX_BUFFER_SIZE=100000000
      export MPICH_MAX_SHORT_MSG_SIZE=4096
      export MPICH_PTL_UNEX_EVENTS=50000
      export MPICH_PTL_OTHER_EVENTS=2496
      ;;
    XC50|xc50)
      export MaxCoresPerNode=40
      export MPITasks="${MPITasks:-120}"
      export ThreadsPerMPITask="${ThreadsPerMPITask:-1}"
      export TasksPerNode=$(( MaxCoresPerNode / ThreadsPerMPITask ))
      export PEs=$(( MPITasks / ThreadsPerMPITask ))
      export Nodes=$(( (MPITasks + MaxCoresPerNode - 1) / MaxCoresPerNode ))
      export Queue="${Queue:-pesq}"
      export WallTime="${WallTime:-01:00:00}"
      export BcCycles="${BcCycles:-0}"
      export MPICH_UNEX_BUFFER_SIZE=100000000
      export MPICH_MAX_SHORT_MSG_SIZE=4096
      export MPICH_PTL_UNEX_EVENTS=50000
      export MPICH_PTL_OTHER_EVENTS=2496
      ;;
    *)
      _log_warn -f "Unknown hpc_name=%s; using conservative defaults." "${hpc_name:-unset}"
      ;;
  esac

  # Executables & fix paths (caller must export these beforehand)
  export execGSI="${execGSI:-${home_cptec}/bin/gsi.x}"
  export parmGSI="${parmGSI:-${home_gsi_fix}/gsiparm.anl}"
  export obsGSI="${obsGSI:-${home_gsi_fix}/obsfiles.rc}"
  export cldRadInfo="${cldRadInfo:-${home_gsi_fix}/cloudy_radiance_info.txt}"
  export SatBiasSample="${SatBiasSample:-${public_fix}/comgsi_satbias_in}"
  export SatBiasPCSample="${SatBiasPCSample:-${public_fix}/comgsi_satbias_pc_in}"
  export ScanInfo="${ScanInfo:-${home_gsi_fix}/global_scaninfo.txt}"
  export SatBiasAngSample="${SatBiasAngSample:-${home_gsi_fix}/global_satangbias.txt}"
  export execBCAng="${execBCAng:-${home_cptec}/bin/global_angupdate}"
  export parmBCAng="${parmBCAng:-${home_gsi_fix}/global_angupdate.namelist}"

  # Satbias names
  export satbiasIn="satbias_in"      satbiasOu="satbias_out"
  export satbiasPCIn="satbias_pc"    satbiasPCOu="satbias_pc.out"
  export satbiasAngIn="satbias_ang.in" satbiasAngOu="satbias_ang.out"
}
#EOC

#-----------------------------------------------------------------------------#
#------------------------ GRID RESOLUTION (BAM) -------------------------------#
#-----------------------------------------------------------------------------#
BAM_CoordSize() {
  local trunc=${1:? "Missing trunc"}
  case "${trunc}" in
    21) TimeStep=3600; IMax=64;   JMax=32  ;;
    31) TimeStep=1800; IMax=96;   JMax=48  ;;
    42) TimeStep=1800; IMax=128;  JMax=64  ;;
    62) TimeStep=900;  IMax=192;  JMax=96  ;;
    106)TimeStep=900;  IMax=320;  JMax=160 ;;
    126)TimeStep=600;  IMax=384;  JMax=192 ;;
    133)TimeStep=600;  IMax=400;  JMax=200 ;;
    159)TimeStep=600;  IMax=480;  JMax=240 ;;
    170)TimeStep=450;  IMax=512;  JMax=256 ;;
    213)TimeStep=300;  IMax=640;  JMax=320 ;;
    254)TimeStep=255;  IMax=768;  JMax=384 ;;
    299)TimeStep=200;  IMax=900;  JMax=450 ;;
    319)TimeStep=225;  IMax=960;  JMax=480 ;;
    341)TimeStep=200;  IMax=1024; JMax=512 ;;
    382)TimeStep=180;  IMax=1152; JMax=576 ;;
    511)TimeStep=150;  IMax=1536; JMax=768 ;;
    533)TimeStep=150;  IMax=1600; JMax=800 ;;
    666)TimeStep=240;  IMax=2000; JMax=1000;;
    863)TimeStep=150;  IMax=2592; JMax=1296;;
    1279)TimeStep=20;  IMax=3840; JMax=1920;;
    1332)TimeStep=20;  IMax=4000; JMax=2000;;
    *) _log_fail "Unknown truncation: %s" "$trunc" ; return 2 ;;
  esac
}

#-----------------------------------------------------------------------------#
#------------------------------ CLI PARSER -----------------------------------#
#-----------------------------------------------------------------------------#
#BOP
# !FUNCTION: ParseOpts
# !INTERFACE: ParseOpts -- <args...>
# !DESCRIPTION:
#   Manual argv parser (no getopts) for runGSI options. Populates environment
#   variables consumed by downstream steps and computes derived fields.
#
# !USAGE:
#   ParseOpts -- -t 299 -l 64 -p PREF -I 2025010100 [-T 299] [-bc 0|1|2...] \
#                 [-np 120] [-N 60] [-d 1] [-om true|false] [-h]
#
# !OPTIONS:
#   -t <int>   Background truncation (BkgTrunc)                 [REQUIRED]
#   -l <int>   Background levels (BkgNLevs)                     [REQUIRED]
#   -p <str>   Background file prefix (BkgPrefix)               [REQUIRED]
#   -I <yyyymmddhh>  Analysis datetime label (AnlDate)          [REQUIRED]
#   -T <int>   Analysis truncation (AnlTrunc); defaults to BkgTrunc
#   -bc <int>  Bias-correction cycle index (BcCycle); default 0
#              • If BcCycle > 0 → IsBC=true ; else → IsBC=false
#   -np <int>  MPI tasks (MPITasks) / legacy alias MTasks
#   -N  <int>  TasksPerNode
#   -d  <int>  ThreadsPerMPITask
#   -om <bool> One-observation test flag (ObsMod), true|false; default false
#   -h         Print usage (extracted from $runGSI) and return 0
#
# !BEHAVIOR:
#   • Validates required options and basic formats
#   • Exports: BkgTrunc, BkgNLevs, BkgPrefix, AnlDate, AnlTrunc,
#              BcCycle, IsBC, ObsMod, (compat) BcCycles=BcCycle
#              and optional MPITasks, TasksPerNode, ThreadsPerMPITask
#   • Computes and exports:
#       - AnlNLevs = BkgNLevs
#       - BkgDate = inctime(AnlDate, -6h)
#       - FctDate = inctime(AnlDate, +6h)
#       - AnlMRES = TQ%04dL%03d (AnlTrunc, AnlNLevs)
#       - BkgMRES = TQ%04dL%03d (BkgTrunc, BkgNLevs)
#
# !RETURNS:
#   0 on success; non-zero on missing/invalid arguments or inctime failure.
#
# !NOTES:
#   Requires: inctime tool in PATH/env; helper 'usage' for -h.
#   Safety: runs under local strict mode (set -euo pipefail) and restores on return.
#   Compatibility: exports BcCycles as alias of BcCycle for older callers.
#EOP
#BOC
ParseOpts() {
  _with_strict_mode   # enable strict mode only for this function

  # defaults
  ObsMod="${ObsMod:-false}"
  BcCycle="${BcCycle:-0}"   # novo nome
  IsBC="${IsBC:-false}"

  # allow an explicit separator to avoid parsing caller flags
  if [[ "${1:-}" == "--" ]]; then shift; fi

  # manual argv loop
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -t) BkgTrunc="${2:?missing value for -t}"; shift 2 ;;
      -l) BkgNLevs="${2:?missing value for -l}"; shift 2 ;;
      -p) BkgPrefix="${2:?missing value for -p}"; shift 2 ;;
      -I) AnlDate="${2:?missing value for -I}"; shift 2 ;;
      -T) AnlTrunc="${2:?missing value for -T}"; shift 2 ;;
      -bc)
        BcCycle="${2:?missing value for -bc}"
        shift 2
        ;;
      -np) MPITasks="${2:?missing value for -np}"; shift 2 ;;
      -N)  TasksPerNode="${2:?missing value for -N}"; shift 2 ;;
      -d)  ThreadsPerMPITask="${2:?missing value for -d}"; shift 2 ;;
      -om) ObsMod="${2:?missing value for -om}"; shift 2 ;;
      -h)
        usage || true
        return 0
        ;;
      --) shift; break ;;          # end of options
      -*)
        _log_fail "Unknown option: %s" "$1"
        usage || true
        return 2
        ;;
      *)
        break ;;
    esac
  done

  # --- validations -----------------------------------------------------------
  : "${BkgTrunc:?BkgTrunc (-t) not set}"
  : "${BkgNLevs:?BkgNLevs (-l) not set}"
  : "${BkgPrefix:?BkgPrefix (-p) not set}"
  : "${AnlDate:?AnlDate (-I) not set}"
  AnlTrunc="${AnlTrunc:-$BkgTrunc}"
  export AnlNLevs="$BkgNLevs"

  [[ "$AnlDate" =~ ^[0-9]{10}$ ]] || { _log_fail "Invalid AnlDate format: %s (expected YYYYMMDDHH)" "$AnlDate"; return 2; }
  [[ "$BkgTrunc" =~ ^[0-9]+$ && "$BkgNLevs" =~ ^[0-9]+$ ]] || { _log_fail "Trunc/levels must be integers"; return 2; }
  [[ "${ObsMod}" =~ ^(true|false)$ ]] || { _log_fail "ObsMod must be 'true' or 'false'"; return 2; }
  [[ "${BcCycle}" =~ ^[0-9]+$ ]] || { _log_fail "-bc must be a non-negative integer (cycle index)"; return 2; }

  # derive IsBC from BcCycle
  if (( BcCycle > 0 )); then
    IsBC=true
  else
    IsBC=false
  fi

  # --- derived fields --------------------------------------------------------
  : "${inctime:?inctime tool not found in PATH/env}"

  export BkgDate
  BkgDate="$("$inctime" "$AnlDate" -6h %y4%m2%d2%h2)"
  export FctDate
  FctDate="$("$inctime" "$AnlDate" +6h %y4%m2%d2%h2)"

  export AnlMRES
  AnlMRES=$(printf 'TQ%4.4dL%3.3d' "$AnlTrunc" "$AnlNLevs")
  export BkgMRES
  BkgMRES=$(printf 'TQ%4.4dL%3.3d' "$BkgTrunc" "$BkgNLevs")

  # export parsed/derived values
  export BkgTrunc BkgNLevs BkgPrefix AnlDate AnlTrunc ObsMod
  export BcCycle IsBC

  [[ -n "${MPITasks:-}"           ]] && export MPITasks
  [[ -n "${TasksPerNode:-}"       ]] && export TasksPerNode
  [[ -n "${ThreadsPerMPITask:-}"  ]] && export ThreadsPerMPITask
}
#EOC

#BOP
# !FUNCTION: linkObs
# !INTERFACE: linkObs RUN_DATE RUN_DIR
# !DESCRIPTION:
#   Stage required observation BUFR files for the GSI run by resolving each
#   OBS_INPUT name declared in ${parmGSI} against search paths derived from
#   ${ncep_ext} and RUN_DATE. For every observation “name” present in both
#   ${parmGSI} and ${obsGSI}, the function computes the timestamped filename
#   via ${inctime} and copies the matching file into RUN_DIR under that name.
#
# !USAGE:
#   linkObs 2025010100 /path/to/run
#
# !ARGUMENTS:
#   RUN_DATE  Cycle datetime (YYYYMMDDHH) used to timestamp input filenames
#   RUN_DIR   Destination directory where selected observation files are copied
#
# !BEHAVIOR:
#   • Builds a colon-separated search list:
#       ${ncep_ext}/YYYY/MM/DD : ${ncep_ext}/YYYYMMDD00/dataout/NCEP
#   • Extracts unique observation names from ${parmGSI} between OBS_INPUT:: and ::
#     (skipping comments/empty lines)
#   • For each name present in ${obsGSI}, resolves its pattern and applies:
#       file="$(${inctime} ${RUN_DATE} +0h <pattern-from-obsGSI>)"
#   • Copies the first file found (in search-order) into RUN_DIR/<name>
#   • Prints a colorized one-line status per file and a summary count at the end
#
# !RETURNS:
#   0 on completion (copy failures are reported to stdout but do not abort)
#   Non-zero only on argument/env validation errors (due to strict mode)
#
# !NOTES:
#   Requires globals: ncep_ext, parmGSI, obsGSI, inctime
#   Copy semantics: cp -pf (preserve mode/mtime; overwrite if exists)
#   File matching relies on ${obsGSI} first column for each observation “name”
#   Safety: runs with localized strict mode (set -euo pipefail); does not modify IFS globally
#   Logging: colorized printf is used; set verbose_local=false to mute file paths
#EOP
#BOC
linkObs_() {
  _with_strict_mode   # enable strict mode only for this function

  local runDate=${1:?}
  local runDir=${2:?}
  local verbose_local=true

  local obsDir
  obsDir="${ncep_ext}/${runDate:0:4}/${runDate:4:2}/${runDate:6:2}"
  obsDir="${obsDir}:${ncep_ext}/${runDate:0:8}00/dataout/NCEP"

  IFS=":" read -r -a obsPath <<<"$obsDir"
  local names
  names="$(sed -n '/OBS_INPUT::/,/::/{/OBS_INPUT/d;/::/d;/^!/d;p}' "${parmGSI}" | awk '{print $1}' | sort -u | xargs)"

  local count=0
  for name in $names; do
    if ! grep -qi -w "$name" "$obsGSI"; then
      printf "\n%s not found in %s list\n\n" "$name" "$obsGSI"
      continue
    fi
    local i filemask file
    for (( i=0; i<${#obsPath[@]}; i++ )); do
      filemask="${obsPath[$i]}/$(grep -iw "$name" "$obsGSI" | awk '{print $1}')"
      file="$("$inctime" "$runDate" +0h "$filemask")" || true
      if [[ -f "$file" ]]; then
        cp -pf -- "$file" "$runDir/$name" 2>/dev/null && {
          printf "\033[34;1m[\033[0m\033[32;1m OK\033[0m"
          ((count++))
        } || {
          printf "\033[34;1m[\033[0m\033[31;1m FAIL\033[0m"
        }
        if $verbose_local; then
          printf "\033[34;1m ]\033[0m link \033[37;1m%s\033[0m @ [ %s ]\n" "$name" "$file"
        else
          printf "\033[34;1m ]\033[0m link \033[37;1m%s\033[0m\n" "$name"
        fi
        break
      elif (( i == ${#obsPath[@]}-1 )); then
        printf "\n\033[31;1m File not found \033[0m\033[34;1m %s\033[0m\n\n" "$(basename -- "$file")"
      fi
    done
  done
  printf "\033[34;1m Found\033[0m\033[32;1m %d\033[0m\033[34;1m observation files to use\033[0m\n" "$count"
}
#EOC
#BOP
# !FUNCTION: linkObs
# !INTERFACE: linkObs RUN_DATE RUN_DIR
# !DESCRIPTION:
#   Stage required observation BUFR files for the GSI run by resolving each
#   OBS_INPUT name declared in ${parmGSI} against search paths derived from
#   ${ncep_ext} and RUN_DATE. For every observation “name” present in both
#   ${parmGSI} and ${obsGSI}, the function computes the timestamped filename
#   via ${inctime} and copies the matching file into RUN_DIR under that name.
#
# !USAGE:
#   linkObs 2025010100 /path/to/run
#
# !ARGUMENTS:
#   RUN_DATE  Cycle datetime (YYYYMMDDHH) used to timestamp input filenames
#   RUN_DIR   Destination directory where selected observation files are copied
#
# !BEHAVIOR:
#   • Builds a colon-separated search list (in priority order):
#       ${ncep_ext}/YYYY/MM/DD : ${ncep_ext}/YYYYMMDD00/dataout/NCEP
#   • Extracts unique observation names from ${parmGSI} between OBS_INPUT:: and ::
#     (skipping comments/empty lines)
#   • Builds a map name→pattern from ${obsGSI} (first token per line; skips comments)
#   • For each <name> present em ambos:
#       candidate="$(${inctime} ${RUN_DATE} +0h <dir>/<pattern-from-obsGSI>)"
#       copia o primeiro arquivo existente para RUN_DIR/<name>
#   • Loga uma linha por arquivo e um sumário ao final
#
# !RETURNS:
#   0 on success (copy failures for individual files don’t abort the whole run)
#   Non-zero only on argument/env validation errors
#
# !REQUIRES:
#   Env vars: ncep_ext, parmGSI, obsGSI, inctime
#   Tools   : awk, sed, sort, xargs, cp
#
# !SAFETY:
#   Local strict mode (set -euo pipefail); IFS preservado; sem efeitos colaterais globais.
#EOP
#BOC
linkObs() {
  local RUN_DATE=${1:? "Missing RUN_DATE (YYYYMMDDHH)"}
  local RUN_DIR=${2:? "Missing RUN_DIR"}

  # Strict mode local
  local _old_opt
  _old_opt=$(set +o) ; set -euo pipefail

  # ---------- Validations ----------
  : "${ncep_ext:?ncep_ext is required}"
  : "${parmGSI:?parmGSI is required}"
  : "${obsGSI:?obsGSI is required}"
  : "${inctime:?inctime is required}"

  [[ "${#RUN_DATE}" -eq 10 ]] || { _log_err "RUN_DATE must be YYYYMMDDHH, got '%s'" "$RUN_DATE"; eval "$_old_opt"; return 2; }
  [[ -d "${RUN_DIR}" ]] || mkdir -p -- "${RUN_DIR}"

  # ---------- Build search paths (priority order) ----------
  local y=${RUN_DATE:0:4} m=${RUN_DATE:4:2} d=${RUN_DATE:6:2} ymd=${RUN_DATE:0:8}
  local -a SEARCH_PATHS=(
    "${ncep_ext}/${y}/${m}/${d}"
    "${ncep_ext}/${ymd}00/dataout/NCEP"
  )
  _log_info "Obs search paths: %s" "$(IFS=:; echo "${SEARCH_PATHS[*]}")"

  # ---------- Extract names from parmGSI (OBS_INPUT block) ----------
  local names
  names="$(
    sed -n '/^[[:space:]]*OBS_INPUT::/,$p' "${parmGSI}" \
    | sed -n '1,/^[[:space:]]*::/p' \
    | sed '1d;$d' \
    | sed -e 's/^\s*!.*$//' -e '/^\s*$/d' \
    | awk '{print $1}' \
    | sort -u \
    | xargs
  )"

  if [[ -z "${names}" ]]; then
    _log_warn -f "No OBS_INPUT names found in %s" "${parmGSI}"
    eval "$_old_opt"; return 0
  fi
  _log_info "Found %d OBS_INPUT names in parmGSI" "$(wc -w <<<"${names}")"

  # ---------- Build name→pattern map from obsGSI ----------
  # Lines like:  prepbufr.gdas.YYYYMMDD.tHHz.nr  PREPBUFR  ...
  # We take the *first* token as the pattern and the second as the name; tolerate either order by normalizing.
  declare -A PATTERN_OF
  # Accept “name pattern” or “pattern name”; ignore comments/blank
  while read -r a b _rest; do
    [[ -z "${a:-}" || "${a#\!}" != "${a}" ]] && continue
    [[ -z "${b:-}" ]] && continue
    if [[ "${a}" == *YYYY* || "${a}" == *%* || "${a}" == *MM* || "${a}" == *DD* || "${a}" == *HH* ]]; then
      # a=pattern, b=name
      PATTERN_OF["$b"]="$a"
    else
      # a=name, b=pattern
      PATTERN_OF["$a"]="$b"
    fi
  done < <(sed -e 's/#.*$//' -e 's/^\s*!.*$//' -e '/^\s*$/d' "${obsGSI}")

  # ---------- Resolve & copy ----------
  local copied=0 missing=0
  local name pattern dir mask cand found
  for name in ${names}; do
    pattern="${PATTERN_OF[$name]:-}"
    if [[ -z "${pattern}" ]]; then
      _log_warn "Name '%s' not present in obsGSI — skipping" "$name"
      ((missing++)) || true
      continue
    fi

    found=""
    for dir in "${SEARCH_PATHS[@]}"; do
      mask="${dir}/${pattern}"
      # inctime expands the timestamped filename from the pattern
      cand="$("${inctime}" "${RUN_DATE}" +0h "${mask}" 2>/dev/null || true)"
      if [[ -n "${cand}" && -f "${cand}" ]]; then
        if cp -pf -- "${cand}" "${RUN_DIR}/${name}" 2>/dev/null; then
          _log_ok "Linked %-18s ← %s" "${name}" "${cand}"
          found="yes"
          ((copied++)) || true
          break
        else
          _log_err "Copy failed for %s → %s" "${cand}" "${RUN_DIR}/${name}"
          found="err"
          break
        fi
      fi
    done

    if [[ -z "${found}" ]]; then
      _log_warn -f "File not found for %-18s (pattern=%s)" "${name}" "${pattern}"
      ((missing++)) || true
    fi
  done

  _log_info "Obs staged into %s" "${RUN_DIR}"
  _log_ok   "Found %d observation file(s); %d missing" "${copied}" "${missing}"

  # restore shell opts
  eval "$_old_opt"
  return 0
}
#EOC

#-----------------------------------------------------------------------------#
#BOP
# !FUNCTION: FixedFiles
# !INTERFACE: FixedFiles Trunc NLevs runDir
# !DESCRIPTION:
#   Copia (ou localiza e copia) arquivos fixos necessários para a execução do GSI,
#   incluindo berror_stats, tabelas e coeficientes do CRTM. Implementa busca
#   robusta por coeficientes sob ${public_crtm} respeitando ${BYTE_ORDER},
#   independente de níveis extras de diretórios.
#
# !USAGE:
#   FixedFiles 299 64 "/path/to/runDir"
#
# !BEHAVIOR:
#   • Valida variáveis de ambiente críticas: public_fix, public_crtm, home_gsi_fix, BYTE_ORDER
#   • Gera rótulo de resolução (ex.: TQ0299L064) e tenta copiar o berror_stats correspondente,
#     caindo para TQ0254L064 se não existir
#   • Copia arquivos fixos do GSI (ozinfo, pcpinfo, errtable, bufrtab, etc.)
#   • Copia arquivos do usuário (anavinfo, satinfo, convinfo, scaninfo)
#   • Localiza e copia coeficientes base do CRTM (EmisCoeff IR/VIS/MW, Aerosol, Cloud)
#     usando `find` sob .../${BYTE_ORDER}/..., primeiro match
#   • Prefere FASTEM6; se ausente, usa FASTEM5 (fallback)
#   • A partir dos sensores listados em `satinfo` (coluna 1), copia SpcCoeff e TauCoeff;
#     para TauCoeff prioriza ODPS (TauCoeff/ODPS) e, se não encontrado, tenta genérico
#   • Mantém o destino “achatado” em runDir (nomes exclusivos), preservando compatibilidade
#
# !EXAMPLE:
#   export public_fix=/pesq/share/das/public/GSI_NCEP/fix
#   export public_crtm=/pesq/share/das/public/CRTM_REL-2.4.0
#   export home_gsi_fix=$HOME/SMG/cptec/gsi/fix
#   export BYTE_ORDER=Big_Endian   # ou Little_Endian
#   FixedFiles 299 64 "/mnt/beegfs/$USER/cycle/2025010100/DAS"
#
# !NOTES:
#   • Requer Bash 4+ (uso de mapfile e <(process substitution))
#   • Integração com helpers de log é opcional; se _log_* não existir, usa printf
#   • Evita vazar flags de `set`: restaura estado do shell ao sair (trap RETURN)
#   • plus_crtm ≡ public_crtm por convenção
#EOP
#BOC
FixedFiles() {
  _with_strict_mode   # enable strict mode only for this function

  # ------------------------ Args & basic checks ------------------------------
  local Trunc=${1:? "Missing Trunc"}
  local NLevs=${2:? "Missing NLevs"}
  local runDir=${3:? "Missing runDir"}

  local mres
  mres=$(printf 'TQ%4.4dL%3.3d' "${Trunc}" "${NLevs}")

  : "${public_fix:?public_fix is unset}"
  : "${public_crtm:?public_crtm is unset}"
  : "${home_gsi_fix:?home_gsi_fix is unset}"
  : "${BYTE_ORDER:?BYTE_ORDER is unset}"   # expected: Big_Endian | Little_Endian

  # Convenção informada: plus_crtm == public_crtm
  local plus_crtm="${public_crtm}"

  # Garante diretório de execução
  mkdir -p -- "${runDir}"

  # Copia arquivo se existir; mensagem padronizada
  _cp_if_exists() { # _cp_if_exists SRC DST_NAME
    local src=$1 dst="$runDir/$2"
    if [[ -f "$src" ]]; then
      _copy_one_safe "$src" "$dst" copy
      return 0
    else
      _log_warn "Missing file: %s" "$src"
      return 1
    fi
  }

  # Primeiro match sob ${public_crtm} com filtro por BYTE_ORDER e, opcionalmente, um subcaminho (hint)
  _first_found_under_byteorder() { # _first_found_under_byteorder PATTERN [SUBPATH_HINT]
    local pattern=$1
    local hint="${2:-}"  # ex.: "EmisCoeff" | "SpcCoeff" | "TauCoeff"
    if [[ -n "$hint" ]]; then
      find "$public_crtm" -type f -path "*/${BYTE_ORDER}/*" -path "*/${hint}/*" -name "$pattern" -print -quit 2>/dev/null || true
    else
      find "$public_crtm" -type f -path "*/${BYTE_ORDER}/*" -name "$pattern" -print -quit 2>/dev/null || true
    fi
  }

  _log_info -f "Starting FixedFiles for %s (%s levels) into %s" "$Trunc" "$NLevs" "$runDir"

  # ------------------------- GSI fixed files ---------------------------------
  # berror_stats: tenta usar o que combina com mres; fallback para TQ0254L064
  if ! _cp_if_exists "${public_ber}/gsir4.berror_stats.gcv.BAM.${mres}" "berror_stats"; then
      _log_err "Could not copy berror_stats (%s missing)" "${mres}"
  fi

  _cp_if_exists "${public_fix}/atms_beamwidth.txt"          "atms_beamwidth.txt"
  _cp_if_exists "${public_fix}/global_ozinfo.txt"           "ozinfo"
  _cp_if_exists "${public_fix}/global_pcpinfo.txt"          "pcpinfo"
  _cp_if_exists "${public_fix}/prepobs_errtable.global"     "errtable"

  if [[ "${ObsMod:-false}" == "true" ]]; then
    _cp_if_exists "${public_fix}/prepobs_prep.bufrtable"    "prepobs_prep.bufrtable"
  fi
  _cp_if_exists "${public_fix}/bufrtab.012"                 "bftab_sstphr"

  # ------------------------- User fixed files --------------------------------
  _cp_if_exists "${home_gsi_fix}/global_anavinfo.l${NLevs}.txt" "anavinfo"
  _cp_if_exists "${home_gsi_fix}/global_satinfo.txt"            "satinfo"
  _cp_if_exists "${home_gsi_fix}/global_convinfo.txt"           "convinfo"
  _cp_if_exists "${home_gsi_fix}/global_scaninfo.txt"           "scaninfo"

  # -------------------------- CRTM base coeffs -------------------------------
  # Preferir FASTEM6 (mais atual); fallback: FASTEM5
  local f_IRwater f_IRice f_IRland f_IRsnow f_VISice f_VISland f_VISsnow f_VISwater f_FASTEM6 f_FASTEM5 f_Aerosol f_Cloud
  f_IRwater=$(_first_found_under_byteorder "Nalli.IRwater.EmisCoeff.bin" "EmisCoeff")
  f_IRice=$(_first_found_under_byteorder   "NPOESS.IRice.EmisCoeff.bin"  "EmisCoeff")
  f_IRland=$(_first_found_under_byteorder  "NPOESS.IRland.EmisCoeff.bin" "EmisCoeff")
  f_IRsnow=$(_first_found_under_byteorder  "NPOESS.IRsnow.EmisCoeff.bin" "EmisCoeff")
  f_VISice=$(_first_found_under_byteorder  "NPOESS.VISice.EmisCoeff.bin" "EmisCoeff")
  f_VISland=$(_first_found_under_byteorder "NPOESS.VISland.EmisCoeff.bin" "EmisCoeff")
  f_VISsnow=$(_first_found_under_byteorder "NPOESS.VISsnow.EmisCoeff.bin" "EmisCoeff")
  f_VISwater=$(_first_found_under_byteorder "NPOESS.VISwater.EmisCoeff.bin" "EmisCoeff")
  f_FASTEM6=$(_first_found_under_byteorder "FASTEM6.MWwater.EmisCoeff.bin" "EmisCoeff")
  f_FASTEM5=$(_first_found_under_byteorder "FASTEM5.MWwater.EmisCoeff.bin" "EmisCoeff")
  f_Aerosol=$(_first_found_under_byteorder "AerosolCoeff.bin" "AerosolCoeff")
  f_Cloud=$(_first_found_under_byteorder   "CloudCoeff.bin"   "CloudCoeff")

  [[ -n "$f_IRwater"  ]] && _cp_if_exists "$f_IRwater"  "Nalli.IRwater.EmisCoeff.bin"
  [[ -n "$f_IRice"    ]] && _cp_if_exists "$f_IRice"    "NPOESS.IRice.EmisCoeff.bin"
  [[ -n "$f_IRland"   ]] && _cp_if_exists "$f_IRland"   "NPOESS.IRland.EmisCoeff.bin"
  [[ -n "$f_IRsnow"   ]] && _cp_if_exists "$f_IRsnow"   "NPOESS.IRsnow.EmisCoeff.bin"
  [[ -n "$f_VISice"   ]] && _cp_if_exists "$f_VISice"   "NPOESS.VISice.EmisCoeff.bin"
  [[ -n "$f_VISland"  ]] && _cp_if_exists "$f_VISland"  "NPOESS.VISland.EmisCoeff.bin"
  [[ -n "$f_VISsnow"  ]] && _cp_if_exists "$f_VISsnow"  "NPOESS.VISsnow.EmisCoeff.bin"
  [[ -n "$f_VISwater" ]] && _cp_if_exists "$f_VISwater" "NPOESS.VISwater.EmisCoeff.bin"

  if [[ -n "$f_FASTEM6" ]]; then
    _cp_if_exists "$f_FASTEM6" "FASTEM6.MWwater.EmisCoeff.bin"
  elif [[ -n "$f_FASTEM5" ]]; then
    _cp_if_exists "$f_FASTEM5" "FASTEM5.MWwater.EmisCoeff.bin"
  else
    _log_warn -f "Neither FASTEM6 nor FASTEM5 EmisCoeff found under %s" "$public_crtm"
  fi

  [[ -n "$f_Aerosol" ]] && _cp_if_exists "$f_Aerosol" "AerosolCoeff.bin"
  [[ -n "$f_Cloud"   ]] && _cp_if_exists "$f_Cloud"   "CloudCoeff.bin"

  # ---------------- SpcCoeff / TauCoeff de acordo com satinfo ----------------
  # Extrai sensores da 1ª coluna (ignora linhas vazias e comentários iniciados por '!')
  if [[ -s "${runDir}/satinfo" ]]; then
    mapfile -t _sensors < <(awk 'NF>0 && $1 !~ /^!/ {print $1}' "${runDir}/satinfo" | sort -u)
    for sensor in "${_sensors[@]}"; do
      # SpcCoeff
      local spc
      spc=$(_first_found_under_byteorder "${sensor}.SpcCoeff.bin" "SpcCoeff")
      if [[ -n "$spc" ]]; then
        _cp_if_exists "$spc" "${sensor}.SpcCoeff.bin"
      else
        _log_warn -f "SpcCoeff not found for %s under %s/%s" "$sensor" "$public_crtm" "$BYTE_ORDER"
      fi

      # TauCoeff – prioriza ODPS; se não achar ODPS, tenta genérico
      # Observação (02/07/2025): ODPS é mais adequado a micro-ondas; p/ ATMS alguns
      # arquivos só existem sob ODPS (p. ex., crtm-2.4.0_emc.1)
      local tau_odps tau_any
      tau_odps=$(find "$public_crtm" -type f -path "*/${BYTE_ORDER}/*" -path "*/TauCoeff/ODPS/*" -name "${sensor}.TauCoeff.bin" -print -quit 2>/dev/null || true)
      tau_any=$(_first_found_under_byteorder "${sensor}.TauCoeff.bin" "TauCoeff")

      if [[ -n "$tau_odps" ]]; then
        _cp_if_exists "$tau_odps" "${sensor}.TauCoeff.bin"
      elif [[ -n "$tau_any" ]]; then
        _cp_if_exists "$tau_any" "${sensor}.TauCoeff.bin"
      else
        _log_warn -f "TauCoeff not found for %s (ODPS or generic) under %s/%s" "$sensor" "$public_crtm" "$BYTE_ORDER"
      fi
    done
  else
    _log_warn -f "satinfo not found or empty at %s; skipping Spc/TauCoeff copy." "${runDir}/satinfo"
  fi

  # -------------------------- Done -------------------------------------------
  _log_ok -f "FixedFiles completed for %s (L=%s) into %s" "$Trunc" "$NLevs" "$runDir"
}
#EOC

#BOP
# !FUNCTION: getInfoFiles
# !INTERFACE: getInfoFiles ANL_DATE
# !DESCRIPTION:
#   Select and stage the most recent metadata files (satinfo, convinfo, ozinfo)
#   whose filename suffix date is <= ANL_DATE from ${home_gsi_fix}/infoFiles.
#   The chosen files are copied into ${runDir} with canonical names:
#   ${runDir}/satinfo, ${runDir}/convinfo, ${runDir}/ozinfo.
#
# !USAGE:
#   getInfoFiles 2025010100
#
# !ARGUMENTS:
#   ANL_DATE   Analysis datetime label (YYYYMMDDHH). The function picks, for each
#              kind in {satinfo, convinfo, ozinfo}, the newest file whose suffix
#              “.YYYYMMDDHH” is <= ANL_DATE.
#
# !BEHAVIOR:
#   • Scans ${home_gsi_fix}/infoFiles for files matching: global_*<kind>*
#   • Sorts candidates by their suffix date (descending) and picks the first <= ANL_DATE
#   • Copies the selected file into ${runDir}/<kind> using cp -pf (preserving attrs)
#   • Exposes the chosen absolute path via exported variables: $satinfo, $convinfo, $ozinfo
#   • Logs one line per kind (OK when staged, ERROR when missing)
#
# !RETURNS:
#   0  if all three (satinfo, convinfo, ozinfo) were resolved and copied
#   N  number of kinds that could not be determined or copied (N > 0)
#
# !NOTES:
#   Requires globals: home_gsi_fix, runDir; helper: assign
#   Directory layout: expects ${home_gsi_fix}/infoFiles/global_*<kind>*.<YYYYMMDDHH>
#   Safety: runs under local strict mode (set -euo pipefail) and restores shell flags on return
#   Example filenames:
#     global_satinfo.something.2024123100
#     global_convinfo.other.2024123006
#     global_ozinfo.XYZ.2024010100
#EOP
#BOC
getInfoFiles() {
  _with_strict_mode   # enable strict mode only for this function

  local anlDate="${1:?missing ANL_DATE}"
  : "${home_gsi_fix:?home_gsi_fix is unset}"
  : "${runDir:?runDir is unset}"

  local infoRoot="${home_gsi_fix}/infoFiles"
  [[ -d "$infoRoot" ]] || { _log_fail "infoFiles dir not found: %s" "$infoRoot"; return 2; }

  local missing=0
  local infoFile kind file infoDate chosen

  for kind in satinfo convinfo ozinfo; do
    chosen=""
    # find candidate files, newest first by suffix date
    while IFS= read -r -d '' file; do
      infoDate="${file##*.}"
      # numeric compare: choose the newest <= anlDate
      if [[ "$anlDate" -ge "$infoDate" ]]; then
        chosen="$file"
        break
      fi
    done < <(find "$infoRoot" -maxdepth 1 -type f -name "global_*${kind}*" -print0 \
             | xargs -0 -I{} bash -c 'f="{}"; echo "${f##*.} ${f}"' \
             | sort -r -k1,1 \
             | awk '{print $2}' ORS='\0')

    if [[ -n "$chosen" ]]; then
      _assign "${kind}" "$chosen"
      cp -pf -- "$chosen" "${runDir}/${kind}"
      _log_ok "Using %s -> %s" "$chosen" "${runDir}/${kind}"
    else
      _log_err "No %s file <= %s under %s" "$kind" "$anlDate" "$infoRoot"
      ((missing++))
    fi
  done

  return $missing
}
#EOC

#BOP
# !FUNCTION: getSatBias
# !INTERFACE: getSatBias ANL_DATE BKG_DATE RUN_DIR DATAOUT_DIR
# !DESCRIPTION:
#   Prepare radiance bias files for the GSI run. The function looks for the
#   most recent satbias_out and satbias_pc.out either in the current RUN_DIR
#   or in DATAOUT_DIR/BKG_DATE. If none are found, seed the run from the
#   provided samples (SatBiasSample / SatBiasPCSample), creating satbias_in
#   and satbias_pc in RUN_DIR.
#
# !USAGE:
#   getSatBias 2025010100 2024123118 "/path/to/run" "/path/to/dataout"
#
# !ARGUMENTS:
#   ANL_DATE    Analysis datetime label (YYYYMMDDHH), used for logs only
#   BKG_DATE    Background datetime (YYYYMMDDHH), used to locate DATAOUT/BKG_DATE
#   RUN_DIR     Target directory where satbias_in and satbias_pc will be written
#   DATAOUT_DIR Root directory where previous cycles are stored
#
# !BEHAVIOR:
#   • Looks for satbias_out → copies to RUN_DIR/satbias_in (warm start)
#   • Looks for satbias_pc.out → copies to RUN_DIR/satbias_pc
#   • If not found, seeds from SatBiasSample / SatBiasPCSample (cold start)
#   • Chooses the “latest” file based on descending sort of matches
#   • Idempotent: safely overwrites targets (cp -pf)
#
# !RETURNS:
#   0  if warm start (previous satbias_out found)
#   1  if cold start (no satbias_out found, seeded from sample)
#   >1 on unexpected errors (e.g., missing environment variables; due to set -e)
#
# !NOTES:
#   Required global variables:
#     • satbiasIn, satbiasOu, satbiasPCIn, satbiasPCOu
#     • SatBiasSample, SatBiasPCSample
#   Safety/robustness:
#     • Runs under strict mode locally (set -euo pipefail), restoring shell state on exit
#     • Uses find + sort to resolve the most recent file
#   Integration:
#     • Caller can branch on return code to select cold/warm namelist, e.g.:
#       if getSatBias ...; then cold=".false."; else cold=".true."; fi
#EOP
#BOC
getSatBias() {
  _with_strict_mode   # enable strict mode only for this function

  local AnlDate="${1:?missing ANL_DATE}"
  local bkgDate="${2:?missing BKG_DATE}"
  local runDir_="${3:?missing RUN_DIR}"
  local dataOut="${4:?missing DATAOUT_DIR}"

  : "${SatBiasSample:?SatBiasSample unset}"
  : "${SatBiasPCSample:?SatBiasPCSample unset}"
  : "${satbiasIn:?satbiasIn unset}"
  : "${satbiasOu:?satbiasOu unset}"
  : "${satbiasPCIn:?satbiasPCIn unset}"
  : "${satbiasPCOu:?satbiasPCOu unset}"

  local cold=0
  local findDir1="${runDir_}"
  local findDir2="${dataOut}/${bkgDate}"

  # satbias_out -> satbias_in
  local fileSatbiasOu
  fileSatbiasOu="$(find "$findDir1" "$findDir2" -type f -iname "${satbiasOu}*" -print 2>/dev/null \
                  | sort -r | head -n1 || true)"

  if [[ -z "$fileSatbiasOu" ]]; then
    _log_warn -f "satbias_out not found; seeding satbias_in from sample"
    rm -f -- "${runDir_}/${satbiasIn}" || true
    cp -pf -- "${SatBiasSample}" "${runDir_}/${satbiasIn}"
    cold=1
  else
    cp -pf -- "${fileSatbiasOu}" "${runDir_}/${satbiasIn}"
    _log_ok "Copied %s -> %s" "$fileSatbiasOu" "${runDir_}/${satbiasIn}"
  fi

  # satbias_pc.out -> satbias_pc
  local fileSatbiasPCOu
  fileSatbiasPCOu="$(find "$findDir1" "$findDir2" -type f -iname "${satbiasPCOu}*" -print 2>/dev/null \
                    | sort -r | head -n1 || true)"
  if [[ -z "$fileSatbiasPCOu" ]]; then
    _log_warn -f "satbias_pc.out not found; seeding satbias_pc from sample"
    rm -f -- "${runDir_}/${satbiasPCIn}" || true
    cp -pf -- "${SatBiasPCSample}" "${runDir_}/${satbiasPCIn}"
  else
    cp -pf -- "${fileSatbiasPCOu}" "${runDir_}/${satbiasPCIn}"
    _log_ok "Copied %s -> %s" "$fileSatbiasPCOu" "${runDir_}/${satbiasPCIn}"
  fi

  return "$cold"
}
#EOC

#BOP
# !FUNCTION: subGSI
# !INTERFACE: subGSI <runDir> <IMAX> <JMAX> <KMAX> <AnlTrunc> <BkgTrunc> <AnlDate> <oneobtest> <cold>
# !DESCRIPTION:
#   Generate and submit a batch job script to run the GSI analysis under SLURM
#   (Egeon) or PBS (XC50). The function creates a job script with proper
#   scheduler directives and runtime environment setup, then submits it with
#   sbatch/qsub.
#
# !USAGE:
#   subGSI /path/to/run 960 480 64 299 299 2025010100 .false. .true.
#
# !ARGUMENTS:
#   runDir     Target directory where the job script will be created/executed
#   IMAX       Horizontal grid size (X-dimension)
#   JMAX       Horizontal grid size (Y-dimension)
#   KMAX       Vertical levels
#   AnlTrunc   Analysis truncation (spectral resolution)
#   BkgTrunc   Background truncation
#   AnlDate    Analysis datetime label (YYYYMMDDHH)
#   oneobtest  Flag: enable/disable one-observation test (.true./.false.)
#   cold       Flag: use cold start namelist (.true./.false.)
#
# !BEHAVIOR:
#   • Copies gsiparm template (cold/normal) into runDir and substitutes placeholders
#   • Generates job script (gsi.qsb) with scheduler settings (SLURM/PBS)
#   • Performs sanity check: TasksPerNode * ThreadsPerMPITask == MaxCoresPerNode
#   • Submits the job with sbatch -W (SLURM) or qsub -W block=true (PBS)
#   • Standard output/error logs are redirected to gsiStdout_* and .out files
#
# !NOTES:
#   Requires: env vars (execGSI, parmGSI, cldRadInfo, hpc_name, Nodes,
#             TasksPerNode, ThreadsPerMPITask, MaxCoresPerNode, MTasks, Queue)
#   Logs: All stdout/stderr routed to runDir/gsiStdout_${AnlDate}.<timestamp>.log
#   Safe under Bash strict mode (set -euo pipefail).
#EOP
#BOC
subGSI() {
  _with_strict_mode   # enable strict mode only for this function

  local runDir_="${1:?}"
  local imax="${2:?}" jmax="${3:?}" kmax="${4:?}"
  local atrc="${5:?}" btrc="${6:?}" andt="${7:?}" onet="${8:?}" cold="${9:?}"

  : "${parmGSI:?}"; : "${cldRadInfo:?}"; : "${execGSI:?}"; : "${hpc_name:?}"

  # prepara gsiparm e info de execução
  cp -pf -- "${cldRadInfo}" "${runDir_}/$(basename -- "${cldRadInfo}")"
  local gsiparm_dst="${runDir_}/$(basename -- "${parmGSI}")"
  if [[ "$cold" == ".true." || "$cold" == "true" || "$cold" == "1" ]]; then
    sed "s/#CENTER#/cptec/g" "${parmGSI}.cold" > "${gsiparm_dst}"
  else
    sed "s/#CENTER#/cptec/g" "${parmGSI}" > "${gsiparm_dst}"
  fi
  sed -i -e "s/#IMAX#/${imax}/g" \
         -e "s/#JMAX#/${jmax}/g" \
         -e "s/#KMAX#/${kmax}/g" \
         -e "s/#TRUNC#/${atrc}/g" \
         -e "s/#TRUNCB#/${btrc}/g" \
         -e "s/#NLVL#/$((${kmax}-1))/g" \
         -e "s/#NLV#/${kmax}/g" \
         -e "s/#LABELANL#/${andt}/g" \
         -e "s/#ONEOBTEST#/${onet}/" "${gsiparm_dst}"

  # sanity: TasksPerNode * ThreadsPerMPITask == MaxCoresPerNode
  local sanity=$(( ${TasksPerNode:?} * ${ThreadsPerMPITask:?} ))
  if (( sanity != ${MaxCoresPerNode:?} )); then
    printf '[ERROR] MPI/OpenMP layout mismatch: %d != %d\n' "$sanity" "$MaxCoresPerNode" >&2
    return 3
  fi

  local runTime; runTime="$(date '+runTime-%H:%M:%S')"
  local jobfile="${runDir_}/gsi.qsb"

  case "${hpc_name}" in
    egeon)
      # ------------------------ SLURM (egeon) ------------------------
      {
        printf '%s\n' '#!/bin/bash'
        printf '%s\n' "#SBATCH --output=${runDir_}/gsiStdout.${andt}.${runTime}.out"
        printf '%s\n' "#SBATCH --error=${runDir_}/gsiStdErr.${andt}.${runTime}.out"
        printf '%s\n' "#SBATCH --time=${WallTime}"
        printf '%s\n' "#SBATCH --nodes=${Nodes}"
        printf '%s\n' "#SBATCH --ntasks=${MTasks}"
        printf '%s\n' "#SBATCH --ntasks-per-node=${TasksPerNode}"
        printf '%s\n' "#SBATCH --cpus-per-task=${ThreadsPerMPITask}"
        printf '%s\n' "#SBATCH --job-name=GSI-%s\n" "${andt}"
        printf '%s\n' "#SBATCH --partition=${Queue}"
        printf '\n'
        printf '%s\n' "cd ${runDir_}"
        printf '%s\n' 'pwd'
        printf '\n'
        printf '%s\n' 'export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK'
        printf '%s\n' 'ulimit -c unlimited'
        printf '%s\n' 'ulimit -s unlimited'
        printf '%s\n' 'export I_MPI_DEBUG=15'
        printf '\n'
        printf '%s\n' "source ${home_gsi}/env.sh egeon ${compiler}"
        printf '%s\n' 'module list'
        printf '\n'
        printf '%s\n' "mpirun -np \$SLURM_NTASKS ./$(basename -- "${execGSI}")"
      } > "${jobfile}"
      ( cd "${runDir_}" && sbatch -W "$(basename -- "${jobfile}")" );;

    XC50|xc50)
      # ------------------------- PBS (XC50) --------------------------
      {
        printf '%s\n' '#!/bin/bash'
        printf '%s\n' "#PBS -o ${runDir_}/gsiAnl.${andt}.${runTime}.out"
        printf '%s\n' '#PBS -S /bin/bash'
        printf '%s\n' "#PBS -q ${Queue}"
        printf '%s\n' "#PBS -l nodes=${Nodes}:ppn=${MaxCoresPerNode}"
        printf '%s\n' '#PBS -N gsiAnl'
        printf '%s\n' '#PBS -A CPTEC'
        printf '\n'
        printf '%s\n' 'set -euo pipefail'
        printf '%s\n' 'ulimit -c unlimited'
        printf '%s\n' 'ulimit -s unlimited'
        printf '%s\n' 'export ATP_ENABLED=1'
        printf '%s\n' 'cd "${PBS_O_WORKDIR:-.}"'
        printf '\n'
        printf '%s\n' "aprun -n ${PEs} -N ${TasksPerNode} -d ${ThreadsPerMPITask} $(basename -- "${execGSI}") > gsiStdout_${andt}.${runTime}.log 2>&1"
      } > "${jobfile}"
      ( cd "${runDir_}" && qsub -W block=true "$(basename -- "${jobfile}")" );;

    *)
      printf '[ERROR] Unsupported hpc_name=%s\n' "${hpc_name}" >&2
      return 4;;
  esac
}
#EOC

#BOP
# !FUNCTION: subBCAng
# !INTERFACE: subBCAng <runDir> <AnlDate>
# !DESCRIPTION:
#   Generate and submit a batch job script to run the angle-dependent radiance
#   bias correction (`global_angupdate`) under SLURM (Egeon) or PBS (XC50).
#   The function prepares the namelist with date placeholders, creates a job
#   script with proper scheduler directives, and submits it.
#
# !USAGE:
#   subBCAng /path/to/run 2025010100
#
# !ARGUMENTS:
#   runDir     Target directory where the job script will be created/executed
#   AnlDate    Analysis datetime label (YYYYMMDDHH)
#
# !BEHAVIOR:
#   • Copies and customizes angle-update namelist (parmBCAng) with date fields
#   • Generates job script (angupdate.qsb) with scheduler settings
#   • Submits job with sbatch -W (SLURM) or qsub -W block=true (PBS)
#   • Redirects program output to gsiAngUpdateStdout_* and .out files
#
# !NOTES:
#   Requires: env vars (execBCAng, parmBCAng, hpc_name, Nodes, TasksPerNode,
#             ThreadsPerMPITask, MTasks, Queue)
#   Logs: Output captured in runDir/gsiAngUpdateStdout_${AnlDate}.<timestamp>.log
#   Cold/warm start handling is not relevant for angle bias correction.
#   Safe under Bash strict mode (set -euo pipefail).
#EOP
#BOC
subBCAng() {
  _with_strict_mode   # enable strict mode only for this function

  local runDir_="${1:?}"
  local andt="${2:?}"

  : "${parmBCAng:?}"; : "${execBCAng:?}"; : "${hpc_name:?}"

  # prepara namelist
  cp -pf -- "${parmBCAng}" "${runDir_}/$(basename -- "${parmBCAng}")"
  sed -i -e "s/#YEAR#/${andt:0:4}/g" \
         -e "s/#MONTH#/${andt:4:2}/g" \
         -e "s/#DAY#/${andt:6:2}/g" \
         -e "s/#HOUR#/${andt:8:2}/g" "${runDir_}/$(basename -- "${parmBCAng}")"

  local runTime; runTime="$(date '+runTime-%H:%M:%S')"
  local jobfile="${runDir_}/angupdate.qsb"

  case "${hpc_name}" in
    egeon)
      # ------------------------ SLURM (egeon) ------------------------
      {
        printf '%s\n' '#!/bin/bash'
        printf '%s\n' "#SBATCH --output=${runDir_}/gsiAngUpdateStdOut.${andt}.${runTime}.out"
        printf '%s\n' "#SBATCH --error=${runDir_}/gsiAngUpdateStdErr.${andt}.${runTime}.out"
        printf '%s\n' '#SBATCH --time=00:15:00'
        printf '%s\n' "#SBATCH --nodes=${Nodes}"
        printf '%s\n' "#SBATCH --ntasks=${MTasks}"
        printf '%s\n' "#SBATCH --ntasks-per-node=${TasksPerNode}"
        printf '%s\n' "#SBATCH --cpus-per-task=${ThreadsPerMPITask}"
        printf '%s\n' "#SBATCH --job-name=AngUp${andt:4:6}"
        printf '%s\n' "#SBATCH --partition=${Queue}"
        printf '\n'
        printf '%s\n' "cd ${runDir_}"
        printf '%s\n' 'pwd'
        printf '\n'
        printf '%s\n' 'export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK'
        printf '%s\n' 'ulimit -c unlimited'
        printf '%s\n' 'ulimit -s unlimited'
        printf '\n'
        printf '%s\n' "source ${home_gsi}/env.sh egeon ${compiler}"
        printf '%s\n' 'module list'
        printf '\n'
        printf '%s\n' "mpirun -np \$SLURM_NTASKS ./$(basename -- "${execBCAng}")"
      } > "${jobfile}"
      ( cd "${runDir_}" && sbatch -W "$(basename -- "${jobfile}")" );;

    XC50|xc50)
      # ------------------------- PBS (XC50) --------------------------
      {
        printf '%s\n' '#!/bin/bash'
        printf '%s\n' "#PBS -o ${runDir_}/gsiAngUpdate.${andt}.${runTime}.out"
        printf '%s\n' '#PBS -S /bin/bash'
        printf '%s\n' "#PBS -q ${Queue}"
        printf '%s\n' '#PBS -l nodes=1:ppn=40'
        printf '%s\n' '#PBS -l walltime=00:05:00'
        printf '%s\n' '#PBS -N gsiAngUpdate'
        printf '%s\n' '#PBS -A CPTEC'
        printf '\n'
        printf '%s\n' 'set -euo pipefail'
        printf '%s\n' 'export ATP_ENABLED=1'
        printf '%s\n' 'cd "${PBS_O_WORKDIR:-.}"'
        printf '\n'
        printf '%s\n' "aprun -n 40 $(basename -- "${execBCAng}") > gsiAngUpdateStdout_${andt}.${runTime}.log 2>&1"
      } > "${jobfile}"
      ( cd "${runDir_}" && qsub -W block=true "$(basename -- "${jobfile}")" );;

    *)
      printf '[ERROR] Unsupported hpc_name=%s\n' "${hpc_name}" >&2
      return 4;;
  esac
}
#EOC

#BOP
# !FUNCTION: mergeDiagFiles
# !INTERFACE: mergeDiagFiles RUN_DIR ANL_DATE
# !DESCRIPTION:
#   Merge per-rank GSI diagnostic files (peNNNN.*) into consolidated
#   files named "diag_<type>_<loop>.<ANL_DATE>" inside RUN_DIR.
#
#   After merging, per-rank parts (peNNNN.<type>_<loop>) are handled by
#   PES_ACTION (default: move → moved into RUN_DIR/diag). Alternatives:
#   delete (remove parts) or keep (leave as-is).
#
# !USAGE:
#   mergeDiagFiles "/path/to/rundir" "YYYYMMDDHH"
#
# !BEHAVIOR:
#   • Validates inputs and RUN_DIR existence
#   • Merges: pe????.<type>_<loop> → diag_<type>_<loop>.<ANL_DATE> (atomic)
#   • Handles per-rank parts according to PES_ACTION (move|delete|keep)
#
# !ENV:
#   PES_ACTION=move|delete|keep
#       Action for per-rank files after merge. Default: move.
#   COPY_PES=true|false  (legacy; mapped to PES_ACTION)
#       true  → move  (default behavior)
#       false → delete
#   parmGSI             (optional) path to GSI namelist; if set, miter is
#                       computed (informational only)
#
# !NOTES:
#   • Safe under `set -u` and handles filenames with spaces.
#   • Single-pass indexing; no repeated 'find' inside the loop.
#   • Atomic merge: write to *.tmp.$$ then mv → final file.
#EOP
#BOC
mergeDiagFiles() {
  _with_strict_mode   # enable strict mode only for this function

  local runDir=${1:? "mergeDiagFiles: missing RUN_DIR"}
  local AnDate=${2:? "mergeDiagFiles: missing ANL_DATE"}
  [[ -d "$runDir" ]] || { printf '[ERROR] RUN_DIR not found: %s\n' "$runDir" >&2; return 2; }

  # Optional: compute miter (+1) from parmGSI if available (informational)
  local miter=""
  if [[ -n "${parmGSI:-}" && -f "${parmGSI}" ]]; then
    miter=$(
      awk -F'[=, ]+' 'BEGIN{IGNORECASE=1}
        /miter[ \t]*=/{for(i=1;i<=NF;i++) if(tolower($i)=="miter"){print $(i+1)+0; exit}}
      ' "${parmGSI}" 2>/dev/null || true
    )
    [[ -n "$miter" ]] && printf '[INFO] Computed miter from parmGSI: %s (no relink)\n' "$((miter+1))"
  fi

  #------------------- per-rank post-merge action (PEs) -----------------------#
  # Back-compat mapping from COPY_PES → PES_ACTION (if not explicitly set)
  local PES_ACTION="${PES_ACTION:-}"
  if [[ -z "$PES_ACTION" && -n "${COPY_PES:-}" ]]; then
    if [[ "$COPY_PES" == "true" ]]; then
      PES_ACTION="move"
    else
      PES_ACTION="delete"
    fi
  fi
  # Default if still empty: move
  case "${PES_ACTION:-move}" in
    move|delete|keep) : ;;
    *) PES_ACTION="move" ;;
  esac

  #-------------------- index files & group by token --------------------------#
  # Expect pattern: peNNNN.<type>_<loop>
  declare -A tokens=()
  local f base token
  for f in "${runDir}"/pe*.*; do
    base=${f##*/}     # basename
    token=${base#*.}  # part after first dot
    [[ "$token" == "$base" ]] && continue   # skip if no dot
    tokens["$token"]=1
  done

  if (( ${#tokens[@]} == 0 )); then
    printf '[WARNING] No files matching pattern pe*.* in %s\n' "$runDir"
    return 0
  fi

  #------------------------------ merge loop ----------------------------------#
  local rc=0
  local token type loop outfile tmpfile diagdir
  local -a token_files abs_files

  for token in "${!tokens[@]}"; do
    # token is "<type>_<loop>"
    loop=${token##*_}
    type=${token%_*}

    if [[ -z "$type" || -z "$loop" ]]; then
      printf '[WARNING] Skipping malformed token: %s\n' "$token" >&2
      continue
    fi

    outfile="${runDir}/diag_${type}_${loop}.${AnDate}"
    tmpfile="${outfile}.tmp.$$"

    # Collect parts for this token, in natural order (pe0000, pe0001, …)
    mapfile -t token_files < <(
      find -P -O3 "$runDir" -maxdepth 1 -type f -name "pe*.${type}_${loop}" -printf '%f\n' \
        | LC_ALL=C sort -V
    )

    if (( ${#token_files[@]} == 0 )); then
      printf '[WARNING] No parts found for %s (pattern: pe*.%s_%s)\n' "$token" "$type" "$loop" >&2
      rc=1
      continue
    fi

    printf '\033[34;1mMerging diag files \033[m\033[32;1m%s\033[m\033[34;1m, outer loop \033[m\033[32;1m%s\033[m \033[34;1m[\033[m' "$type" "$loop"

    # Atomic merge: write to tmp then move to final
    : > "$tmpfile"
    {
      local part
      for part in "${token_files[@]}"; do
        cat -- "${runDir}/${part}" >> "$tmpfile"
      done
      mv -f -- "$tmpfile" "$outfile"
    } 2>/dev/null || true

    if [[ -s "$outfile" ]]; then
      printf '\033[32;1m OK \033[m\033[34;1m]\033[m\n'
    else
      printf '\033[31;1m FAIL \033[m\033[34;1m]\033[m\n' >&2
      rc=1
      [[ -e "$tmpfile" ]] && rm -f -- "$tmpfile"
      continue
    fi

    #--------------- handle per-rank parts as per PES_ACTION ------------------#
    # Build absolute paths from token_files (safe for spaces)
    abs_files=()
    if (( ${#token_files[@]} > 0 )); then
      local t
      for t in "${token_files[@]}"; do
        abs_files+=( "${runDir}/${t}" )
      done
    fi

    case "$PES_ACTION" in
      move)
        diagdir="${runDir}/diag"
        mkdir -p -- "$diagdir"
        # Try mv first (fast, atomic across same FS). Fallback to copy+remove.
        if ! mv -f -- "${abs_files[@]}" "$diagdir/" 2>/dev/null; then
          # Cross-FS or permission corner cases: copy then remove originals
          if command -v rsync >/dev/null 2>&1; then
            if [[ "${verbose:-false}" == "true" ]]; then
              rsync -a --human-readable --info=progress2,name0 -i -- "${abs_files[@]}" "$diagdir/" || true
            else
              rsync -a -q -- "${abs_files[@]}" "$diagdir/" || true
            fi
          else
            cp -f -- "${abs_files[@]}" "$diagdir/" 2>/dev/null || true
          fi
          # Remove originals only if copies exist at destination
          rm -f -- "${abs_files[@]}" 2>/dev/null || true
        fi
        ;;
      delete)
        rm -f -- "${abs_files[@]}" 2>/dev/null || true
        ;;
      keep)
        : # do nothing
        ;;
    esac

    # (Reference only; final-link disabled)
    # if [[ -n "$miter" && "$loop" == "$(printf "%02d" $((miter+1)))" ]]; then
    #   ln -sf -- "$outfile" "${runDir}/diag_${type}.${AnDate}"
    # fi
  done

  return "$rc"
}
#EOC

#BOP
# !FUNCTION: copyFiles
# !INTERFACE: copyFiles FROM_DIR TO_DIR
# !DESCRIPTION:
#   Archive outputs from a completed GSI run. Copies diagnostics, analysis
#   results, bias correction files, configuration files, and logs from the
#   run directory (FROM_DIR) into the archive directory (TO_DIR).
#
# !USAGE:
#   copyFiles /path/to/run /path/to/archive
#
# !ARGUMENTS:
#   FROM_DIR   Source run directory containing GSI outputs
#   TO_DIR     Destination archive directory (created if missing)
#
# !BEHAVIOR:
#   • Copies all diag_* and fort.* outputs
#   • Preserves legacy bundle of pe* diagnostics under diag/ if present
#   • Copies BAM.anl → renamed as GANL<BkgPrefix><AnlDate>S.unf.<BkgMRES>
#   • Copies satbias_in/out and satbias_pc/in/out if found
#   • Copies key config files (gsiparm.anl, anavinfo, satinfo, convinfo)
#   • Copies GSI stdout logs (gsiStdout*)
#   • All copies use cp -pf to preserve timestamps/permissions
#
# !RETURNS:
#   0 on success (all existing files copied)
#   Non-zero only if cp/find fail under strict mode
#
# !NOTES:
#   Requires globals: BkgPrefix, AnlDate, BkgMRES,
#                     satbiasIn, satbiasOu, satbiasPCIn, satbiasPCOu
#   Safe under Bash strict mode (set -euo pipefail) with shell state restored.
#   Diagnostic pe* shards are not merged here — see mergeDiagFiles.
#   By design, no files are removed from FROM_DIR; this is a copy-only archiver.
#EOP
#BOC
copyFiles() {
  _with_strict_mode   # enable strict mode only for this function

  local fromDir="${1:?}"; local toDir="${2:?}"
  : "${BkgPrefix:?BkgPrefix unset}"
  : "${AnlDate:?AnlDate unset}"
  : "${BkgMRES:?BkgMRES unset}"
  : "${satbiasIn:?}"; : "${satbiasOu:?}"; : "${satbiasPCIn:?}"; : "${satbiasPCOu:?}"

  mkdir -p -- "${toDir}"

  # GSI diagnostics & fortran unit files
  find -P -O3 "${fromDir}" -maxdepth 1 -type f -iname "diag_*"  -exec cp -pf -- {} "${toDir}" \;
  find -P -O3 "${fromDir}" -maxdepth 1 -type f -iname "fort.*"  -exec cp -pf -- {} "${toDir}" \;

  # legacy bundle of pe* diag files
  if [[ -d "${fromDir}/diag" ]]; then
    cp -prf -- "${fromDir}/diag" "${toDir}/"
  fi

  # analysis file
  if [[ -f "${fromDir}/BAM.anl" ]]; then
    cp -pf -- "${fromDir}/BAM.anl" "${toDir}/GANL${BkgPrefix}${AnlDate}S.unf.${BkgMRES}"
  fi

  # bias files
  for f in "${satbiasIn}" "${satbiasOu}" "${satbiasPCIn}" "${satbiasPCOu}"; do
    [[ -n "$f" && -f "${fromDir}/$f" ]] && cp -pf -- "${fromDir}/$f" "${toDir}/"
  done

  # configs & logs
  for f in gsiparm.anl anavinfo satinfo convinfo; do
    [[ -f "${fromDir}/${f}" ]] && cp -pf -- "${fromDir}/${f}" "${toDir}/"
  done
  find -P -O3 "${fromDir}" -maxdepth 1 -type f -name "gsiStdout*" -exec cp -pf -- {} "${toDir}" \;
}
#EOC


#-----------------------------------------------------------------------------#
# Simple colored logger (level-aware)
#-----------------------------------------------------------------------------#
#BOP
# !FUNCTION: .log
# !INTERFACE: .log LEVEL MESSAGE...
# !DESCRIPTION:
#   Minimal colorized logger; LEVEL is 0..7 (emerg..debug).
#   This helper ignores the global 'verbose'. Use _log_* for standardized logs.
#EOP
#BOC
function .log () {
  local LOG_LEVELS=([0]="emerg"  [1]="alert" [2]="crit"  [3]="err" \
                    [4]="warning"[5]="notice"[6]="info"  [7]="debug")
  local LOG_COLORS=([0]="\033[31;1m" [1]="\033[33;1m" [2]="\033[31;1m" [3]="\033[31;1m" \
                    [4]="\033[33;1m" [5]="\033[32;1m" [6]="\033[34;1m" [7]="\033[32;1m")
  local LEVEL="${1:?missing level}"; shift
  local msg="$*"
  if [[ "$LEVEL" =~ ^[0-7]$ ]]; then
    printf "%s[%s]\033[0m %s\n" "${LOG_COLORS[$LEVEL]}" "${LOG_LEVELS[$LEVEL]}" "$msg"
  else
    printf "[log] %s\n" "$msg"
  fi
}
#EOC

if [[ "${BASH_SOURCE[0]}" == "$0" && -z "${SMG_SETUP_AS_LIBRARY}" ]]; then
  
  # If a valid function is called, execute it; otherwise, show help
  if declare -F "$1" >/dev/null 2>&1; then
    cmd="$1"; shift
    "$cmd" "$@"
  else
    _log_err "Unknown command: $1"
    echo
    show_help
    exit 1
  fi
fi

