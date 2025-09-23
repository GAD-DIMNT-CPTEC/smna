#!/usr/bin/env bash
#-----------------------------------------------------------------------------#
#------------------------------ HELPERS LOADER --------------------------------#
#-----------------------------------------------------------------------------#
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
#   • Caches the discovered path in LIBS_PATH (dir), speeding future calls.
#EOP
#BOC
source_nearby() {
  local target="${1:?filename required}"
  local start="${2:-}"

  # --------------------- logging control ---------------------
  local _verbose="${SOURCE_NEARBY_VERBOSE:-true}"
  _log() { [[ "${_verbose}" == true ]] && printf '%s\n' "$*" >&2 || true; }

  # ----------------- helper: resolve script dir ---------------
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
    _log "[INFO] Loaded (explicit path): $target"
    # cache directory
    export LIBS_PATH="$(cd -- "$(dirname -- "$target")" && pwd -P)"
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

  # Default prune dirs + extras
  local -a PRUNE_DIRS=(.git .svn .hg build dist node_modules __pycache__ .mypy_cache .venv venv .tox .cache .idea .vscode)
  if [[ -n "${SOURCE_NEARBY_PRUNE:-}" ]]; then
    local IFS=':'; read -r -a _prune_extra <<< "${SOURCE_NEARBY_PRUNE}"
    PRUNE_DIRS+=("${_prune_extra[@]}")
  fi
  
  # Constrói a expressão de prune como ARRAY, sem NUL e sem command substitution
  local -a PRUNE_EXPR=()
  if ((${#PRUNE_DIRS[@]})); then
    for d in "${PRUNE_DIRS[@]}"; do
      [[ -n "$d" ]] && PRUNE_EXPR+=(-name "$d" -o)
    done
    # remove o -o final
    ((${#PRUNE_EXPR[@]})) && unset 'PRUNE_EXPR[${#PRUNE_EXPR[@]}-1]'
  fi
  
  # Busca (GNU/FreeBSD find com -maxdepth)
  _find_here() {
    local base="$1"
    find "$base" -xdev -maxdepth "${SOURCE_NEARBY_MAXDEPTH:-4}" \
      \( -type d \( "${PRUNE_EXPR[@]}" \) -prune \) -o \
      \( -type f -name "$target" -print -quit \) 2>/dev/null
  }

  # ----------------- 1) same directory -----------------------
  if [[ -f "$script_dir/$target" ]]; then
    # shellcheck disable=SC1091
    source "$script_dir/$target"
    _log "[INFO] Loaded: $script_dir/$target"
    export LIBS_PATH="$script_dir"
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
      export LIBS_PATH="$(cd -- "$(dirname -- "$hit")" && pwd -P)"
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
    export LIBS_PATH="$(cd -- "$(dirname -- "$hit")" && pwd -P)"
    return 0
  fi

  printf '[ERROR] Could not find "%s" from "%s" or under "%s".\n' \
    "$target" "$script_dir" "$global_root" >&2
  return 1
}
#EOP

# --- colors for warnings ---
if [[ -t 2 ]]; then
  : "${C_WARN:=$'\033[1;33m'}"
  : "${C_RST:=$'\033[0m'}"
else
  : "${C_WARN:=}"
  : "${C_RST:=}"
fi
export C_WARN C_RST

# --- Load helpers only once ---
if ! ${__HELPERS_SH_LOADED:-false}; then
  if ! source_nearby "__helpers__.sh"; then   # seguro mesmo com 'set -e'
    rc=$?
    # source_nearby já emite [ERROR]/[WARN]; aqui só registramos o código
    printf "%s[WARN]%s __helpers__.sh load failed (rc=%d); continuing\n" \
           "$C_WARN" "$C_RST" "$rc" 1>&2
  fi
fi

