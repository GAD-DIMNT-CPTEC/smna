#!/usr/bin/env bash
#=============================================================================#
#                            GLOBAL INIT / HELPERS                            #
#=============================================================================#
# Modo estrito antecipado (antes de carregar helpers). Mantém robustez geral.
set -o errexit -o nounset -o pipefail

# Full silence if the first arg requests help (even forced messages)
if [[ "${1:-}" =~ ^(-h|h|--help)$ ]]; then
  export SUPPRESS_LOGS=true
fi

# --------------------- logging control ---------------------
# Minimal printf-based log until helpers are available
# (INFO/ERROR only)
__ml() {
  [[ "${SUPPRESS_LOGS:-false}" == true ]] && return 0
  local lvl="${1^^}"; shift
  local c_rst=$'\033[0m' c_info=$'\033[1;34m' c_err=$'\033[1;31m'
  local tag="[${lvl}]"
  [[ "$lvl" == "INFO"  ]] && printf '%s%s%s %s\n' "$c_info" "$tag" "$c_rst" "$*" && return
  [[ "$lvl" == "ERROR" ]] && printf '%s%s%s %s\n' "$c_err"  "$tag" "$c_rst" "$*" && return
  printf '%s\n' "$tag $*"
}

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
#   LIBS_PATH                  Optional absolute path to the file or a directory
#                              containing it. If valid, takes precedence.
#   SOURCE_NEARBY_VERBOSE      true|false (default: true). Controls [INFO]/[WARN].
#   SOURCE_NEARBY_MAXDEPTH     Max depth per directory level (default: 4).
#   SOURCE_NEARBY_GLOBAL_ROOT  Fallback root for global search (default: "/").
#   SOURCE_NEARBY_PRUNE        Colon-separated extra dirs to prune from search.
#
# !RETURNS:
#   0 on success; 1 if target not found; >1 for unexpected internal errors.
#
# !NOTES:
#   • Safe under 'set -euo pipefail'.
#   • Uses 'find -xdev -maxdepth N' to avoid crossing mounts and limit depth.
#   • Caches the discovered directory in LIBS_PATH, speeding future calls.
#EOP
#BOC
source_nearby() {
  local target="${1:?filename required}"
  local start="${2:-}"


  # ----------------- helper: resolve script dir ---------------
  # If this file is sourced, BASH_SOURCE[0] points to the current file;
  # if executed, it still works. Fall back to $PWD if unresolved.
  _script_dir() {
    local src="${BASH_SOURCE[0]:-}"
    if [[ -n "$src" && -e "$src" ]]; then
      ( cd -- "$(dirname -- "$src")" && pwd -P )
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
    export LIBS_PATH="$(cd -- "$(dirname -- "$target")" && pwd -P)"
    __ml info "Loaded (explicit path): $target"
    return 0
  fi

  # -------------------- 0) LIBS_PATH hint --------------------
  if [[ -n "${LIBS_PATH:-}" ]]; then
    if [[ -f "$LIBS_PATH" ]]; then
      # shellcheck disable=SC1090
      source "$LIBS_PATH"
      export LIBS_PATH="$(cd -- "$(dirname -- "$LIBS_PATH")" && pwd -P)"
      __ml info "Loaded via LIBS_PATH (file): $LIBS_PATH"
      return 0
    elif [[ -d "$LIBS_PATH" && -f "$LIBS_PATH/$target" ]]; then
      # shellcheck disable=SC1090
      source "$LIBS_PATH/$target"
      __ml info "Loaded via LIBS_PATH (dir): $LIBS_PATH/$target"
      return 0
    fi
    __ml error "LIBS_PATH set but not valid for '$target': $LIBS_PATH"
  fi

  # ----------------- search parameters -----------------------
  local maxdepth="${SOURCE_NEARBY_MAXDEPTH:-4}"
  local global_root="${SOURCE_NEARBY_GLOBAL_ROOT:-/}"

  # Default prune dirs + extras
  local -a PRUNE_DIRS=(.git .svn .hg build dist node_modules __pycache__ .mypy_cache .venv venv .tox .cache .idea .vscode)
  if [[ -n "${SOURCE_NEARBY_PRUNE:-}" ]]; then
    local IFS=':'; read -r -a _prune_extra <<< "${SOURCE_NEARBY_PRUNE}"
    PRUNE_DIRS+=("${_prune_extra[@]}")
  fi

  # Build prune expression array (no trailing -o)
  local -a PRUNE_EXPR=()
  if ((${#PRUNE_DIRS[@]})); then
    for d in "${PRUNE_DIRS[@]}"; do
      [[ -n "$d" ]] && PRUNE_EXPR+=(-name "$d" -o)
    done
    ((${#PRUNE_EXPR[@]})) && unset 'PRUNE_EXPR[${#PRUNE_EXPR[@]}-1]'
  fi

  # GNU/FreeBSD find with -maxdepth. (BSD/macOS find lacks -maxdepth.)
  _find_here() {
    local base="$1"
    find "$base" -xdev -maxdepth "$maxdepth" \
      \( -type d \( "${PRUNE_EXPR[@]}" \) -prune \) -o \
      \( -type f -name "$target" -print -quit \) 2>/dev/null
  }

  # ----------------- 1) same directory -----------------------
  if [[ -f "$script_dir/$target" ]]; then
    # shellcheck disable=SC1091
    source "$script_dir/$target"
    export LIBS_PATH="$script_dir"
    __ml info "Loaded: $script_dir/$target"
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
      export LIBS_PATH="$(cd -- "$(dirname -- "$hit")" && pwd -P)"
      __ml info "Loaded: $hit"
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
    export LIBS_PATH="$(cd -- "$(dirname -- "$hit")" && pwd -P)"
    __ml info "Loaded (global): $hit"
    return 0
  fi

  printf '[ERROR] Could not find "%s" from "%s" or under "%s".\n' \
    "$target" "$script_dir" "$global_root" >&2
  return 1
}
#EOC

#=============================================================================#
#BOP
# !FUNCTION: ensure_helpers_and_env
# !INTERFACE: ensure_helpers_and_env [<helpers_file>]
# !DESCRIPTION:
#   Ensure the base helpers are sourced via `source_nearby` and immediately
#   (re)export environment settings by calling `_env_export`. It also disables
#   any active Conda environment by invoking `disable_conda`.
#   This function NEVER falls back to POSIX dot sourcing; if `source_nearby`
#   is unavailable or fails, an error code is returned.
#
# !USAGE:
#   # In any script that depends on helpers/_env_export/_with_strict_mode:
#   ensure_helpers_and_env               # loads "__helpers__.sh" using source_nearby
#   ensure_helpers_and_env "../lib/__helpers__.sh"
#
# !BEHAVIOR:
#   • Requires: source_nearby <file>
#   • If __HELPERS_SH_LOADED=true, it skips re-sourcing helpers and just calls _env_export.
#   • Otherwise, sources the given helpers file via source_nearby, then calls _env_export.
#   • Always calls _env_export after ensuring helpers are present.
#   • Always calls disable_conda to ensure Conda is deactivated.
#
# !RETURNS:
#   0 on success; non-zero if:
#     - source_nearby is missing, or
#     - sourcing helpers fails, or
#     - _env_export is missing/failed, or
#     - disable_conda is missing/failed.
#
# !DEPENDENCIES:
#   • source_nearby: function that resolves and sources helper files relative
#     to the script location (not the current working directory).
#   • __helpers__.sh: must define _env_export, disable_conda, and (optionally)
#     _with_strict_mode and logging helpers (_log_info/_warn/_err).
#
# !EXAMPLE:
#   ensure_helpers_and_env || { echo "[ERROR] helpers/env setup failed"; exit 1; }
#
# !NOTES:
#   - This function does NOT use _with_strict_mode because it is provided by
#     the helpers being loaded here.
#   - Idempotent: safe to call multiple times; it will always re-run _env_export
#     and re-disable Conda if active.
#EOP
#BOC
ensure_helpers_and_env() {

  local helpers_file="${1:-__helpers__.sh}"
  local rc=0

  # Load helpers only if not loaded yet
  if [[ "${__HELPERS_SH_LOADED:-false}" != true ]]; then
    source_nearby "${helpers_file}"
    rc=$?
    if (( $rc != 0 )); then
      rc=$?
      __ml error "source_nearby '${helpers_file}' failed (rc=${rc})"
      return "${rc}"
    fi
    export __HELPERS_SH_LOADED=true
  fi

  # From this point, logging helpers may exist; fall back to _eh_log otherwise
  if declare -F _log_info >/dev/null; then
    _log_info "Helpers loaded: %s" "${helpers_file}"
  else
    __ml info "Helpers loaded: ${helpers_file}"
  fi

  # _env_export must exist and succeed; always run it (idempotent)
  if ! declare -F _env_export >/dev/null; then
    __ml error "_env_export not found after loading helpers"
    return 3
  fi
  _env_export
  rc=$?
  if (( $rc != 0 )); then
    if declare -F _log_err >/dev/null; then
      _log_err "Environment export failed (rc=%d)" "${rc}"
    else
      __ml error "Environment export failed (rc=${rc})"
    fi
    return "${rc}"
  fi

  # Checks if Conda is active and deactivates it if necessary
  if ! declare -F disable_conda >/dev/null; then
    _log_err "disable_conda not found after loading helpers"
    return 4
  fi
  disable_conda
  rc=$?
  if (( $rc != 0 )); then
     _log_err "disable_conda failed (rc=%d). Aborting." "${rc}"
    return "${rc}"
  fi

  return 0
}
#EOC

#=============================================================================#
# Initialize helpers + environment, then parse CLI arguments
#=============================================================================#

# Ensure that base helpers are loaded and environment variables are exported.
# This will source "__helpers__.sh" if not already loaded and run _env_export.
ensure_helpers_and_env __helpers__.sh

# Ensure leftover_args exists as a global array (avoids "unbound variable"
# errors under 'set -u' when resetting "$@" with its contents).
declare -ga leftover_args=()

# Call the argument parser provided by helpers. It will populate leftover_args[]
# with positional arguments not consumed by option parsing.
__parse_args__ "$@"

# Replace "$@" with the filtered leftover_args so downstream code only sees
# true positional arguments and not global/common flags.
set -- "${leftover_args[@]}"

# -- debug info --- #
_log_debug "queue=%s job=%s walltime=%s dry_run=%s verbose=%s debug=%s" \
              "${queue:-}" "${job_name:-}" "${walltime:-}" "${dry_run}" "${verbose}" "${debug}"
_log_debug "MPI=%s OMP=%s Nodes=%s C/Node=%s Procs=%s" \
              "${mpi_tasks:-}" "${omp_threads}" "${nodes:-}" "${cores_per_node:-}" "${total_procs:-}"
