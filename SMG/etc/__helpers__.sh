#!/usr/bin/env bash
#===============================================================================
# Script: __helpers__.sh
#
# Purpose:
#   Common helper utilities (logging, argument parsing, progress reporting,
#   simple copying/linking with progress, environment assignment/expansion,
#   system detection, and cluster environment export). Designed to be *sourced*
#   by other Bash scripts.
#
# Usage:
#   # Recommended: source this file at the top of your script
#   source "/path/to/__helpers__.sh"
#
#   # Example pattern inside caller scripts:
#   my_command() {
#     _parse_args "$@"
#     $verbose && _log_info "Running my_command with args: %s" "${leftover_args[*]}"
#     # ...
#   }
#
# Notes:
#   - This library is intended to be sourced (not executed). The shebang is
#     present for consistency only; no code runs on direct execution.
#   - Global flag `verbose` defaults to `false`. Set `verbose=true` in the
#     caller to enable logging, or pass `-v/--verbose` to functions that call
#     `_parse_args`.
#   - Functions prefixed with `_` are internal helpers; they can be used by
#     callers, but their interface may evolve conservatively.
#   - Requires Bash 4+ (for namerefs with `local -n`).
#
# Environment:
#   verbose      : boolean string "true"/"false" (default: false)
#   auto_yes     : when true, skip confirmation prompts (default: false)
#   dry_run      : when true, only log actions (default: false)
#   do_restore   : when true, enable restore behavior if supported (default: false)
#   do_fix       : when true, allow creating dirs/symlinks/copies if needed (default: false)
#
# License:
#   LGPL-3.0-or-later (adjust as needed for your project)
#===============================================================================
# guard: avoid re-sourcing this file multiple times
__HELPERS_SH_LOADED=${__HELPERS_SH_LOADED:-false}

if $__HELPERS_SH_LOADED; then
  return 0
fi

__HELPERS_SH_LOADED=true

# --- Default logging verbosity (exported for downstream scripts) ---
export verbose=${verbose:-false}

# --- Default logging debugging (exported for downstream scripts) ---
export debug=${debug:-false}

# --- Absolute path to the directory where this helpers script resides ---
export helpers_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

#-----------------------------------------------------------------------------#
#----------------------- internal helpers (hidden from help) -----------------#
#-----------------------------------------------------------------------------#
#BOP
# !FUNCTION: _with_strict_mode
# !INTERFACE: _with_strict_mode
# !DESCRIPTION:
#   Enable "strict mode" (set -euo pipefail) inside the current function scope
#   only. It saves the current shell options with `set +o` and registers a trap
#   on RETURN to restore them when the function exits (either normally or due
#   to an error).
#
# !USAGE:
#   my_function() {
#     _with_strict_mode
#     # function logic here
#   }
#
# !BEHAVIOR:
#   • Activates safer execution by:
#       - `-e`  → exit immediately on error
#       - `-u`  → error on unset variables
#       - `-o pipefail` → pipelines fail if any command fails
#   • The original shell options are restored automatically after the function
#     finishes, ensuring no side effects on the caller environment.
#
# !EXAMPLE:
#   safe_copy() {
#     _with_strict_mode
#     cp -- "$1" "$2"
#   }
#
#   # In this example, any error in cp (missing file, permission issue, etc.)
#   # will stop execution, but once safe_copy returns the original shell
#   # options are restored.
#EOP
_with_strict_mode() {
  # Snapshot current options as shell code (e.g., 'set +o errexit; set +o nounset; ...')
  local __o; __o="$(set +o)"

  # Enable strict mode
  set -euo pipefail

  # Restore snapshot when this function returns (expand __o now, not later)
  trap -- "$(printf 'eval %q; trap - RETURN' "$__o")" RETURN
}


#BOP
# !FUNCTION: __parse_args__
# !INTERFACE: __parse_args__ "$@"
# !DESCRIPTION:
#   Parse common options and component selection flags for compilation.
#
#   Verbosity and utilities:
#     -v, --verbose      → verbose=true  (enable logging)
#     -q, --quiet        → verbose=false (disable logging)
#     -d, --debug        → debug=true    (enable logging)
#     -y, --yes          → auto_yes=true (auto-confirm prompts)
#     -f, --fix          → do_fix=true   (create dirs/symlinks/copies as needed)
#     --dry-run          → dry_run=true  (skip actual actions)
#     --restore          → do_restore=true (restore from *.bak when supported)
#
#   Notes:
#     1) Flags are evaluated in order: the last one wins (e.g. --verbose followed by --quiet disables verbosity).
#
#   Remaining positional arguments are preserved in leftover_args[].
#
# !EXAMPLE:
#   my_function() {
#     _parse_args "$@"
#     $verbose && _log_info "Components: GSI=$compgsi ANG=$compang BAM=$compbam INC=$compinctime"
#     if $dry_run; then
#       _log_info "DRY-RUN: would run with args: %s" "${leftover_args[*]}"
#     fi
#   }
#EOP
#BOC
__parse_args__() {
  _with_strict_mode   # enable strict mode only for this function

  # ---- Safe defaults (respect existing env vars, otherwise set to false) ----
  verbose=${verbose:-false}
  debug=${debug:-false}
  auto_yes=${auto_yes:-false}
  dry_run=${dry_run:-false}
  do_restore=${do_restore:-false}
  do_fix=${do_fix:-false}

  # Declare leftover_args explicitly (global) to avoid -u surprises
  # If you don't want to rely on 'declare -g', just omit it; the assignment below will create it.
  declare -g -a leftover_args 2>/dev/null || true
  leftover_args=()
  
  while [[ $# -gt 0 ]]; do
    case "$1" in
      # ---- Verbosity / utilities ----
      -v|--verbose) verbose=true ;;
      -q|--quiet)   verbose=false ;;
      -d|--debug)   debug=true ;;
      -y|--yes)     auto_yes=true ;;
      -f|--fix)     do_fix=true ;;
      --dry-run)    dry_run=true ;;
      --restore)    do_restore=true ;;

      # ---- Stop parsing options ----
      --) shift; break ;;   # explicit end of options
      -*) break ;;          # unknown short/long flag → stop
      *)  break ;;          # first positional arg → stop
    esac
    shift
  done

  # ---- Preserve leftover positional arguments ----
  leftover_args=("$@")
}
#EOC

#-----------------------------------------------------------------------------#
# Colored logging helpers (drop-in, keeps the same signatures)
#-----------------------------------------------------------------------------#

#BOP
# !FUNCTION: _log_colors_init
# !DESCRIPTION:
#   Initialize ANSI color tags when stdout is an interactive TTY and COLOR=1.
#   When not a TTY or COLOR != 1, colors are disabled (empty strings).
#   This keeps logs clean when redirected to files or piped.
#
# !NOTES:
#   • Uses stdout (fd 1) as reference for color enablement.
#   • Honors env var COLOR (default: 1). Set COLOR=0 to force no color.
#EOP
#BOC
_log_colors_init() {
  local _COLOR="${COLOR:-1}"
  if [[ -t 1 && "${_COLOR}" = "1" ]]; then
    C_INFO=$'\033[1;34m'   # bold blue
    C_DBG=$'\033[1;35m'    # bold purple
    C_OK=$'\033[1;32m'     # bold green
    C_WARN=$'\033[1;33m'   # bold yellow
    C_ERR=$'\033[1;31m'    # bold red
    C_ACT=$'\033[1;36m'    # bold cyan (ACTION)
    C_RST=$'\033[0m'
  else
    C_INFO=; C_OK=; C_WARN=; C_ERR=; C_ACT=; C_RST=
  fi
}
# Initialize at load time (no-op if not TTY)
_log_colors_init
#EOC


#BOP
# !FUNCTION: _log_msg
# !DESCRIPTION:
#   Print a standardized log message with a given level (INFO, OK, WARNING,
#   ACTION, ERROR, FAIL, etc.). Messages are printed only when the global
#   variable `verbose` is set to `true`, unless forced with `-f`.
#   Supports printf-style formatting. Safe with 'set -u'.
#
# !INTERFACE:
#   _log_msg <LEVEL> [-f] <format> [args...]
#
# !EXAMPLES:
#   _log_msg INFO "Starting step %s" "$step"
#   _log_msg WARNING -f "Low disk space: %s" "$mountpoint"
#   _log_msg OK "Built %d target(s) in %0.2f s" "$n" "$elapsed"
#
# !NOTES:
#   • Do not add a trailing newline to <format>; it’s appended automatically.
#   • `verbose` is the Bash boolean string `true` or `false` (default: false).
#   • Output goes to stdout by default; to route to stderr temporarily, wrap:
#       { _log_msg ERR "msg"; } 1>&2
#EOP
#BOC
_log_msg() {
  local level="$1"; shift
  local force=false
  if [[ "${1:-}" == "-f" ]]; then force=true; shift; fi

  # default: verbose=false
  local v="${verbose:-false}"
  # default: debug=false
  local d="${debug:-false}"

  # Unless explicitly forced:
  $force || { [[ "${level^^}" == "DEBUG" ]] && { $d || return 0; } \
              || { $v || [[ "$level" =~ ^(OK|WARNING|WARN|ERROR|ERR|FAIL|ACTION)$ ]] || return 0; }; }

  # map level → color tag + normalized label
  local tag color
  case "${level^^}" in
    DEBUG)    tag="[DEBUG]";   color="${C_DBG:-}";;
    INFO)     tag="[INFO]";    color="${C_INFO:-}";;
    OK)       tag="[OK]";      color="${C_OK:-}";;
    WARNING|WARN)
              tag="[WARNING]"; color="${C_WARN:-}";;
    ERROR|ERR)
              tag="[ERROR]";   color="${C_ERR:-}";;
    FAIL)     tag="[FAIL]";    color="${C_ERR:-}";;
    ACTION)   tag="[ACTION]";  color="${C_ACT:-}";;
    *)        tag="[${level}]"; color="";;
  esac

  # $1 is the format string; the rest are printf args
  local fmt; fmt="${1:-}"; shift || true

  if [[ -z "$fmt" ]]; then
    printf "%s%s%s\n" "${color}" "${tag}" "${C_RST}"
  else
    # colorize only the tag; message remains plain (better for grepping)
    # use -- to stop printf from interpretarando algo como opção
    printf "%s%s%s %s\n" "${color}" "${tag}" "${C_RST}" "$(printf -- "$fmt" "$@")"
  fi

}
#EOC

#BOP
# !FUNCTION: _log_debug, _log_info, _log_ok, _log_err, _log_warning, _log_action, _log_fail, _die
# !DESCRIPTION:
#   Convenience wrappers around `_log_msg` that set the appropriate level.
#   `_log_err` and `_log_fail` always force output (they behave as if `-f` was passed).
#   Other wrappers accept an optional `-f` to force output.
#
# !INTERFACE:
#   _log_debug   [-f] <format> [args...]
#   _log_info    [-f] <format> [args...]
#   _log_ok      [-f] <format> [args...]
#   _log_warning [-f] <format> [args...]
#   _log_action  [-f] <format> [args...]
#   _log_err          <format> [args...]   # forced output
#   _log_fail         <format> [args...]   # forced output
#   _die       [code] <format> [args...]
#
# !EXAMPLES:
#   _log_info "Preparing NCEP input copy (layout=%s): %s" "$layout" "$cycle"
#   _log_ok   "Artifacts available at %s" "$outdir"
#   _log_warning -f "Retrying download (%d/%d)..." "$i" "$max"
#   _log_err  "Failed to load cluster paths via vars_export"
#
# !NOTES:
#   • Formatting follows `printf` semantics (placeholders are expanded).
#   • Wrappers append a newline automatically; do not include one in <format>.
#EOP
#BOC


# Wrappers (keep behavior consistent)
_log_debug()   { _log_msg "DEBUG"   "$@"; }
_log_info()    { _log_msg "INFO"    "$@"; }
_log_ok()      { _log_msg "OK"      "$@"; }
_log_warn()    { _log_msg "WARNING" "$@"; }
_log_action()  { _log_msg "ACTION"  "$@"; }

# Errors should always print, regardless of verbose
_log_fail()    { _log_msg "FAIL"  -f "$@"; }
_log_err()     { _log_msg "ERROR" -f "$@"; }
_die()         { local code="${1:-1}"; shift || true; _log_err "$@"; exit "$code"; }

#EOC


#BOP
# !FUNCTION: _progress_bar
# !INTERFACE: _progress_bar <current> <total> <message>
# !DESCRIPTION:
#   Render a simple textual progress bar (width=40) with percentage and counters.
#
# !USAGE:
#   _progress_bar 3 10 "Copying files"
#
# !BEHAVIOR:
#   • Prints a single updating line with carriage return
#   • Shows percent complete and current/total count
#
# !EXAMPLE:
#   for i in {1..10}; do
#     _progress_bar "$i" 10 "Processing"
#     sleep 0.2
#   done
#EOP
#BOC
_progress_bar(){
  local current=${1:-0} total=${2:-0} msg="${3:-Processing}"
  local width=40
  (( total <= 0 )) && { printf "\r[INFO] %s [no files]        \n" "$msg"; return; }
  (( current < 0 )) && current=0
  (( current > total )) && current=$total
  local progress=$(( current * width / total ))
  local percent=$(( 100 * current / total ))
  local bar=""
  for ((i=0; i<width; i++)); do
    if (( i < progress )); then bar+="█"; else bar+="░"; fi
  done
  printf "\r$C_ACT[ACTION]$C_RST %s [%s] %3d%% (%d/%d) " "$msg" "$bar" "$percent" "$current" "$total"
}
#EOC
#BOP
# !FUNCTION: _copy_one_safe
# !INTERFACE: _copy_one_safe SRC DST [METHOD]
# !DESCRIPTION:
#   Safely copy or link a single file from SRC to DST.
#
#   METHOD:
#     - copy     : atomic copy (cp -p to temp, then mv), preserves mode/mtime
#     - hardlink : create a hard link (fails if across filesystems)
#     - symlink  : create/replace a symlink to SRC (absolute target)
#     Aliases: link → hardlink, simlink → symlink
#
#   Behavior:
#     - If DST is an existing directory, the file is placed as DST/basename(SRC).
#     - Parent directory of final destination is created as needed.
#     - Returns non-zero on errors; prints concise error messages to stderr.
#
# !USAGE:
#   _copy_one_safe /data/in/a.bin /data/out copy
#   _copy_one_safe /data/in/a.bin /data/out/another.bin symlink
#   _copy_one_safe /data/in/a.bin /data/out link
#
# !RETURNS:
#   0  ok
#   1  invalid arguments
#   2  source missing or not a regular file
#   3  hardlink across different filesystems
#   4  operation failed (cp/ln/mv errors)
#
# !NOTES:
#   - Portable stat device check (GNU/BSD) to detect cross-FS for hardlink.
#   - Uses absolute SRC path for symlink to avoid broken links when cwd changes.
#EOP
#BOC
_copy_one_safe() {
  _with_strict_mode   # enable strict mode only for this function

  local src="${1:-}"; local dst="${2:-}"; local method="${3:-copy}"
  [[ -n "${src}" && -n "${dst}" ]] || { _log_err "copy_one_safe: missing SRC or DST\n" >&2; return 1; }

  # normalize method aliases
  case "${method}" in
    copy|hardlink|symlink) ;;
    link)    method="hardlink" ;;   # alias
    simlink) method="symlink"  ;;   # alias (common typo)
    *) _log_err "copy_one_safe: unknown METHOD '%s' (expected copy|hardlink|symlink)\n" "${method}" >&2; return 1 ;;
  esac

  # source checks
  if [[ ! -e "${src}" ]]; then
    _log_err "copy_one_safe: source does not exist: %s\n" "${src}" >&2
    return 2
  fi
  if [[ ! -f "${src}" ]]; then
    _log_err "copy_one_safe: source is not a regular file: %s\n" "${src}" >&2
    return 2
  fi

  # resolve destination path: if DST is a directory, append basename(SRC)
  local dst_path="${dst}"
  if [[ -d "${dst_path}" ]]; then
    dst_path="${dst%/}/$(basename -- "${src}")"
  fi
  # ensure parent dir exists
  mkdir -p -- "$(dirname -- "${dst_path}")" || { _log_err "mkdir -p '%s' failed\n" "$(dirname -- "${dst_path}")" >&2; return 4; }

  # absolute source path (for symlink target)
  local src_abs
  src_abs="$(cd -- "$(dirname -- "${src}")" && pwd)/$(basename -- "${src}")"

  # device id helper (GNU/BSD)
  _dev_id() { stat -Lc %d "$1" 2>/dev/null || stat -f %d "$1" 2>/dev/null; }

  case "${method}" in
    copy)
      # atomic-ish copy: cp to temp in same dir, then mv
      local tmp="${dst_path}.tmp.$$"
      # -p preserve mode/mtime; -f overwrite temp if exists
      cp -pf -- "${src_abs}" "${tmp}" 2>/dev/null || { _log_err "cp '%s' -> '%s' failed\n" "${src_abs}" "${tmp}" >&2; rm -f -- "${tmp}" 2>/dev/null || true; return 4; }
      mv -f -- "${tmp}" "${dst_path}" 2>/dev/null || { _log_err "mv '%s' -> '%s' failed\n" "${tmp}" "${dst_path}" >&2; rm -f -- "${tmp}" 2>/dev/null || true; return 4; }
      ;;

    hardlink)
      # must be same filesystem
      local sdev ddev
      sdev="$(_dev_id "$(dirname -- "${src_abs}")")" || { _log_err "stat device(src) failed\n" >&2; return 4; }
      ddev="$(_dev_id "$(dirname -- "${dst_path}")")" || { _log_err  "stat device(dst) failed\n" >&2; return 4; }
      if [[ "${sdev}" != "${ddev}" ]]; then
        printf "[ERROR] hardlink across filesystems: %s -> %s\n" "${src_abs}" "${dst_path}" >&2
        return 3
      fi
      ln -f -- "${src_abs}" "${dst_path}" 2>/dev/null || { _log_err "ln (hardlink) failed: %s -> %s\n" "${src_abs}" "${dst_path}" >&2; return 4; }
      ;;

    symlink)
      # replace existing file/symlink
      ln -sfn -- "${src_abs}" "${dst_path}" 2>/dev/null || { _log_err "ln -s failed: %s -> %s\n" "${src_abs}" "${dst_path}" >&2; return 4; }
      ;;
  esac

  return 0
}
#EOC


#BOP
# !FUNCTION: _copy_with_progress
# !INTERFACE: _copy_with_progress <arrayname> <src_dir> <dest_dir> [copy|link] [msg]
# !DESCRIPTION:
#   Copy or symlink a list of files from a source directory to a destination,
#   showing a simple progress bar on TTY (fallback mode when rsync isn't used).
#   Skips safely when source and destination are the same.
#
# !USAGE:
#   files=("a.txt" "b.txt")
#   _copy_with_progress files /tmp/src /tmp/dest copy "Copying inputs"
#
# !OPTIONS:
#   copy   (default)  Use `cp -pf` to copy files
#   link              Use `ln -sf` to create symbolic links
#
# !BEHAVIOR:
#   • Ensures destination directory exists
#   • Iterates through file list and copies/links one by one
#   • On TTY, calls `_progress_bar` to render visual feedback; otherwise logs a single line
#
# !NOTES:
#   - Depends on `_progress_bar` when stdout is a TTY (`[[ -t 1 ]]`).
#   - The first argument is the *name* of the array (without [@]).
#   - action: copy|symlink|hardlink  (default: copy)
#   - Requires Bash 4+ (nameref) and optional _progress_bar/_log_* helpers.
#EOP
#BOC
_copy_with_progress() {
  _with_strict_mode 2>/dev/null || true  # enable strict mode only for this function if available

  local -n _files_ref="$1"; shift        # nameref to the input array
  local src_dir="$1"; shift
  local dest_dir="$1"; shift
  local action="${1:-copy}"; shift || true
  local progress_msg="${1:-Processing files}"

  mkdir -p -- "$dest_dir" || { _log_err "mkdir -p '%s' failed" "$dest_dir"; return 1; }

  local total=${#_files_ref[@]}
  if (( total == 0 )); then
    _log_info -f "%s [nothing to do]" "$progress_msg"
    return 0
  fi

  local count=0 f src dst

  for f in "${_files_ref[@]}"; do
    ((count++))
    src="${src_dir%/}/$f"
    dst="${dest_dir%/}/$f"

    if [[ ! -e "$src" ]]; then
      _log_warn "Missing: %s" "$src"
      continue
    fi

    # Skip if src and dst are the same file (same inode/device)
    if [[ -e "$dst" ]] && [[ "$src" -ef "$dst" ]]; then
      _log_debug "Same file (skipping): %s" "$dst"
      [[ -t 1 ]] && _progress_bar "$count" "$total" "$progress_msg"
      continue
    fi

    # Show progress/log before executing
    if [[ -t 1 ]]; then
      _progress_bar "$count" "$total" "$progress_msg"
    else
      _log_action -f "(%d/%d) %s -> %s" "$count" "$total" "$src" "$dst"
    fi

    case "$action" in
      copy)
        # If destination exists and contents are identical, skip quietly
        if [[ -e "$dst" ]] && cmp -s -- "$src" "$dst"; then
          _log_debug "Unchanged (skipping): %s" "$dst"
        else
          # Ensure parent dir exists (in case list has subpaths)
          mkdir -p -- "$(dirname -- "$dst")" || { _log_err "mkdir -p '%s' failed" "$(dirname -- "$dst")"; return 1; }
          cp -pf -- "$src" "$dst" || { _log_err "Failed to copy '%s' -> '%s'" "$src" "$dst"; return 1; }
        fi
        ;;
      symlink|link|ln|sym)
        mkdir -p -- "$(dirname -- "$dst")" || { _log_err "mkdir -p '%s' failed" "$(dirname -- "$dst")"; return 1; }
        ln -sfn -- "$src" "$dst" || { _log_err "Failed to symlink '%s' -> '%s'" "$src" "$dst"; return 1; }
        ;;
      hardlink|hlink|hln)
        mkdir -p -- "$(dirname -- "$dst")" || { _log_err "mkdir -p '%s' failed" "$(dirname -- "$dst")"; return 1; }
        # Hard links must be on the same filesystem; if fails, fall back to copy
        if ! ln -f -- "$src" "$dst" 2>/dev/null; then
          _log_warn "Hardlink failed (cross-device?), falling back to copy: %s" "$dst"
          cp -pf -- "$src" "$dst" || { _log_err "Failed to copy '%s' -> '%s' (fallback)" "$src" "$dst"; return 1; }
        fi
        ;;
      *)
        _log_err "Unknown action '%s' (expected: copy|symlink|hardlink)" "$action"
        return 2
        ;;
    esac
  done

  # In TTY, the bar draws inline; append a clean OK marker once.
  if [[ -t 1 ]]; then
    printf " - $C_OK[OK]$C_RST\n"
  else
    _log_ok "%s done" "$progress_msg"
  fi
  return 0
}
#EOC

#BOP
# !FUNCTION: _copy_with_progress_
# !INTERFACE: _copy_with_progress_ <array_name> <src_dir> <dest_dir> [copy|link] ["Message"]
# !DESCRIPTION:
#   Copy or link a list of files from a source directory to a destination
#   directory, with progress reporting. Supports two modes:
#     • If `rsync` is available → uses `rsync -a --human-readable --info=progress2`
#       for detailed progress (percentage, bytes, ETA).
#     • If `rsync` is not available → falls back to `_copy_with_progress`
#       which draws a simple progress bar (TTY only).
#   Also supports dry-run mode (global var `dry_run=true`) where operations are
#   only logged, not executed.
#
# !USAGE:
#   files=(file1 file2 file3)
#   _copy_with_progress_ files /path/src /path/dest copy "Copying input files"
#
# !OPTIONS:
#   array_name   Name of the bash array variable containing filenames (not expanded with [@]).
#   src_dir      Source directory containing the files.
#   dest_dir     Destination directory for copy/link.
#   action       Either "copy" (default) or "link".
#   message      Optional label for logging (default: "Copying files").
#
# !BEHAVIOR:
#   • Creates destination directory if missing.
#   • For each file in array:
#       - Verifies existence in src_dir.
#       - Executes copy/link or logs action if dry-run.
#       - Logs progress:
#           · With rsync: clean numeric progress (percentage, size, ETA).
#           · Without rsync: falls back to `_copy_with_progress` (TTY bar).
#   • Reports total processed vs expected at the end.
#
# !NOTES:
#   - Requires `_log_info`, `_log_warn`, `_log_action`, `_log_ok` helpers for logging.
#   - Respects global `dry_run` variable (true/false).
#   - Designed for internal use inside higher-level workflows.
#EOP
#BOC
_copy_with_progress_() {
  _with_strict_mode   # enable strict mode only for this function

  local files_name="$1"; shift                 # array name (e.g., files), without [@]
  local -n _files_ref="$files_name"            # nameref to iterate here
  local src_dir="$1"; shift
  local dest_dir="$1"; shift
  local action="${1:-copy}"; shift || true
  local progress_msg="${1:-Copying files}"
  local dry="${dry_run:-false}"

  mkdir -p -- "$dest_dir" || { _log_err "mkdir -p failed: %s" "$dest_dir"; return 1; }

  local total=${#_files_ref[@]}
  echo "total: ${total}"
  if (( total == 0 )); then
    _log_info "%s [nothing to do]" "$progress_msg"
    return 0
  fi

  # No rsync: delegate to simple bar (with dry-run handling)
  if ! command -v rsync >/dev/null 2>&1; then
    if $dry; then
      local f
      for f in "${_files_ref[@]}"; do
        _log_info "Dry-run: would %s %s -> %s/" "$action" "${src_dir%/}/$f" "$dest_dir"
      done
      _log_ok "%s: %d/%d processed (dry-run)" "$progress_msg" "$total" "$total"
      return 0
    fi
    _copy_with_progress "$files_name" "$src_dir" "$dest_dir" "$action" "$progress_msg"
    return $?
  fi

  # With rsync: show --info=progress2
  local n=0 ok=0 f src dst
  for f in "${_files_ref[@]}"; do
    src="${src_dir%/}/$f"
    dst="${dest_dir%/}/$f"

    if [[ ! -e "$src" ]]; then
      _log_warn "Missing: %s" "$src"
      continue
    fi

    ((n++))
    _log_action "(%d/%d) %s -> %s" "$n" "$total" "$src" "$dst"

    case "$action" in
      copy)
        if $dry; then
          _log_info "Dry-run: would copy %s -> %s/" "$src" "$dest_dir"
          ((ok++))
        else
          # preflight: will transfer?
          if rsync -ai --dry-run -- "$src" "$dest_dir"/ | grep -q '^[^\.]'; then
            rsync -a --human-readable --partial \
                  --info=progress2,name0 \
                  -i --out-format='%i %n' -- "$src" "$dest_dir"/ && ((ok++)) || _log_warn "rsync failed: %s" "$src"
          else
            printf "[INFO] Up-to-date: %s -> %s/\n" "$src" "$dest_dir"
            ((ok++))
          fi
        fi
        ;;
      link)
        if $dry; then
          _log_info "Dry-run: would link %s -> %s" "$src" "$dst"
          ((ok++))
        else
          if ln -sf -- "$src" "$dst"; then
            ((ok++))
          else
            _log_warn "link failed: %s" "$src"
          fi
        fi
        ;;
      *)
        _log_warn "Unknown action: %s (skipping)" "$action"
        ;;
    esac
  done

  _log_ok "%s: %d/%d processed" "$progress_msg" "$ok" "$total"
}
#EOC


#BOP
# !FUNCTION: _copy_dir_with_progress
# !INTERFACE: _copy_dir_with_progress <src_dir> <dest_dir> [msg] [action]
# !DESCRIPTION:
#   Enumerate all files in a directory and copy or link them to a destination,
#   showing progress. Uses `_copy_with_progress` internally.
#
# !USAGE:
#   _copy_dir_with_progress /src/path /dest/path ["message"] [copy|link]
#
# !OPTIONS:
#   copy   (default)  Copy files with `cp`
#   link              Symlink files with `ln`
#
# !BEHAVIOR:
#   • Collects all files in the source directory
#   • Builds a list of basenames
#   • Calls `_copy_with_progress` to handle the actual operation
#
# !EXAMPLE:
#   _copy_dir_with_progress ./PRE/datain/2019/sst ./datain/2019/sst "Copying SST files" copy
#
# !NOTES:
#   Depends on `_copy_with_progress`
#EOP
#BOC
_copy_dir_with_progress(){
  _with_strict_mode   # enable strict mode only for this function

  local src_dir="$1"
  local dest_dir="$2"
  local progress_msg="${3:-Copying files}"
  local action="${4:-copy}"

  mkdir -p -- "$dest_dir" || { _log_err "mkdir -p '$dest_dir' failed"; return 1; }

  local files=()
  while IFS= read -r -d '' f; do
    files+=( "$(basename "$f")" )
  done < <(find -L "$src_dir" -maxdepth 1 -type f -print0)

  _copy_with_progress files "$src_dir" "$dest_dir" "$action" "$progress_msg"
}
#EOC

#BOP
# !FUNCTION: _assign
# !INTERFACE: _assign <KEY> <VALUE...>
# !DESCRIPTION:
#   Define and export an environment variable from a KEY and VALUE pair,
#   expanding ${VAR} references using the current environment and variables
#   already set by previous _assign calls.
#
# !USAGE:
#   _assign HOME /home/${USER}
#   _assign subt_smg ${SUBMIT_HOME}/${nome_smg}
#
# !BEHAVIOR:
#   • Joins VALUE arguments into a single string
#   • Expands variable references like ${USER}, ${SUBMIT_HOME}
#   • If `envsubst` is available, uses it for expansion
#   • Otherwise, uses a safe eval/printf fallback (disallows command substitution)
#   • Exports KEY=expanded_value to the shell environment
#
# !NOTES:
#   • Order matters: referenced variables must be defined earlier
#   • Does not allow backticks, $(), or arithmetic $(( )) substitutions
#EOP
#BOC
_assign() {
  local key="$1"; shift
  # Join VALUE... preserving internal spaces
  local raw="$*"

  # Safe expansion using envsubst, if available
  if command -v envsubst >/dev/null 2>&1; then
    local expanded
    expanded="$(printf '%s' "$raw" | envsubst)" || return 1
    export "$key=$expanded"
    return 0
  fi

  # Fallback without envsubst: block command substitution
  if printf '%s' "$raw" | grep -Eq '(`|\$\(|\$\(\()'; then
    _log_err "_assign: forbidden command substitution in value for $key" >&2
    return 1
  fi

  # Perform only variable expansion ${VAR} / $VAR using eval+printf
  # shellcheck disable=SC2086
  local expanded
  expanded="$(eval "printf '%s' \"$raw\"")" || return 1
  export "$key=$expanded"
}
#EOC

#BOP
# !FUNCTION: _list_files_array
# !INTERFACE: _list_files_array VAR_NAME SRC_DIR [find_args...]
# !DESCRIPTION:
#   Populate the array named VAR_NAME with the basenames of files found in
#   SRC_DIR (optionally filtered with extra find(1) arguments).
#
# !USAGE:
#   _list_files_array filesArray /path/to/src -name '*.txt'
#
# !BEHAVIOR:
#   • Validates inputs and clears the target array
#   • Runs a pruned, single-level search: find -L SRC_DIR -maxdepth 1 -type f ...
#   • Stores only basenames into VAR_NAME
#   • Emits [DEBUG]/[INFO] logs if $verbose=true
#
# !NOTES:
#   • Requires Bash 4+ (nameref)
#   • Safe if no matches: returns 0 with an empty array
#EOP
#BOC
_list_files_array() {
  # optional strict mode if your helper exists
  if declare -F with_strict_mode >/dev/null 2>&1; then with_strict_mode; fi

  local __outvar="${1:?VAR_NAME required}"
  local src_dir="${2:?SRC_DIR required}"
  shift 2 || true

  # nameref to external array
  local -n __arr_ref="$__outvar" 2>/dev/null || {
    _log_err "Target array '%s' is not a valid name" "$__outvar"
    return 2
  }
  __arr_ref=()

  if [[ ! -d "$src_dir" ]]; then
    _log_warn "Source directory not found: %s" "$src_dir"
    return 0
  fi

  # Log filters (quoted to keep wildcards literal)
  _log_debug "Listing files in %s (filters: %q)" "$src_dir" "$*"

  # Collect matches (null-delimited)
  local -a __found=()
  if mapfile -d '' -t __found < <(find -L "$src_dir" -maxdepth 1 -type f "$@" -print0 2>/dev/null); then
    :
  else
    # find error (permissions, bad predicate, etc.)
    _log_warn "find failed in %s with args: %q" "$src_dir" "$*"
    return 0
  fi

  # Fill array with basenames
  local p
  for p in "${__found[@]}"; do
    __arr_ref+=( "$(basename -- "$p")" )
  done

  _log_info "Found %d file(s) in %s" "${#__arr_ref[@]}" "$src_dir"
  return 0
}
#EOC

#BOP
# !FUNCTION: _env_export
# !INTERFACE: _env_export
# !DESCRIPTION:
#   Load and export cluster-specific environment variables from a config file
#   in `etc/mach/${hpc_name}_paths.conf`. Each line has the format:
#     KEY   VALUE
#   Blank lines and lines starting with '#' are ignored. Leading/trailing
#   whitespace is trimmed. Uses _assign to perform variable expansion and export.
#
# !USAGE:
#   _env_export
#
# !BEHAVIOR:
#   • Builds config file path: $(dirname ${BASH_SOURCE})/mach/${hpc_name}_paths.conf
#   • Aborts if the file is missing
#   • Reads each non-empty, non-comment line
#   • Splits into KEY and VALUE
#   • Calls _assign KEY VALUE (which expands ${VAR} references)
#   • Exports all resulting variables
#
# !EXAMPLE:
#   # inside mach/egeon_paths.conf
#   HOME        /home/${USER}
#   subt_smg    ${SUBMIT_HOME}/${nome_smg}
#
#   # usage
#   _env_export   # sets HOME and subt_smg expanded
#
# !NOTES:
#   • Depends on: _assign
#   • Requires: hpc_name set to match a config file
#EOP
#BOC
_env_export(){
  _with_strict_mode   # enable strict mode only for this function
  # Avoid re-running if already done
  if [[ "${ENV_EXPORTED:-false}" == "true" ]]; then
    return 0
  fi

  # Ensure hpc_name is defined before using it
  if [[ -z "${hpc_name:-}" ]]; then
    _log_warn "hpc_name was not set before loading environment paths!"
    # Try to detect HPC system automatically
    if ! detect_hpc_system; then
      _log_err "detect_hpc_system failed. Aborting."
      return 1
    fi
    _log_info -f "hpc_name set to: %s" "$hpc_name"
  fi

  local confdir filepaths
  confdir="${helpers_dir}/mach"
  filepaths="${confdir}/${hpc_name}_paths.conf"

  [[ -f "$filepaths" ]] || { echo "[ERROR] Missing: $filepaths" >&2; return 1; }

  # Read lines "KEY  VALUE"; ignore comments/blank lines
  while IFS= read -r line; do
    [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
    # key = first token; value = the rest (left-trim)
    local key="${line%%[[:space:]]*}"
    local value="${line#"$key"}"
    value="${value#"${value%%[![:space:]]*}"}"   # trim leading spaces
    [[ -z "$key" || -z "$value" ]] && continue
    _assign "$key" "$value"
  done < "$filepaths"

  # Mark as already exported
  export ENV_EXPORTED=true
  
}
#EOC

#BOP
# !FUNCTION: _bootstrap_env_root
# !INTERFACE: _bootstrap_env_root <env_var> <anchor_file>
# !DESCRIPTION:
#   Ensure <env_var> (e.g., SMG_ROOT) is exported and points to the project root.
#   If <env_var> is already set, do nothing. Otherwise, climb parent directories
#   starting from the caller script location until <anchor_file> is found; export
#   its parent directory to <env_var>. If not found, fall back to the caller
#   script directory.
# !EXAMPLE:
#   _bootstrap_env_root SMG_ROOT "config_smg.sh"
#EOP
#BOC
_bootstrap_env_root() {
  _with_strict_mode 2>/dev/null || true

  local env_var="${1:?environment variable name required}"    # e.g., SMG_ROOT
  local anchor="${2:?anchor file required}"                   # e.g., config_smg.sh

  # If already defined, do nothing
  local current
  eval "current=\"\${$env_var:-}\""
  if [[ -n "$current" ]]; then
    _log_debug "%s already set → %s" "$env_var" "$current"
    return 0
  fi

  # Resolve the caller script absolute directory (follow symlinks)
  local this_file this_dir cur
  _resolve_this_file this_file
  this_dir="$(cd -- "$(dirname -- "$this_file")" && pwd -P)"

  cur="$this_dir"
  while [[ "$cur" != "/" ]]; do
    if [[ -f "$cur/$anchor" ]]; then
      eval "export $env_var=\"\$cur\""
      _log_info "Set %s → %s (found %s)" "$env_var" "$cur" "$anchor"
      return 0
    fi
    cur="$(dirname -- "$cur")"
  done

  # Fallback: not found; use the caller script directory
  _log_warn "%s not found while climbing; using script dir" "$anchor"
  eval "export $env_var=\"\$this_dir\""
  return 1
}
#EOC

#BOP
# !FUNCTION: _run_project_config
# !INTERFACE: _run_project_config <env_var> <config_filename> <target_fn> [args...]
# !DESCRIPTION:
#   Source the config file located at "<$env_var>/<config_filename>" and invoke
#   the given <target_fn> within that file, passing optional [args...].
#   Fail-soft: logs warnings and returns non-zero if anything goes wrong.
# !EXAMPLE:
#   _run_project_config SMG_ROOT "config_smg.sh" _env_export
#EOP
#BOC
_run_project_config() {
  _with_strict_mode 2>/dev/null || true

  local env_var="${1:?env var required}"                # e.g., SMG_ROOT
  local cfg_file="${2:?config filename required}"       # e.g., config_smg.sh
  local target="${3:?target function required}"         # e.g., _env_export
  shift 3 || true
  local -a target_args=( "$@" )

  # Resolve root from env_var
  local root
  eval "root=\"\${$env_var:-}\""
  if [[ -z "$root" ]]; then
    _log_err "%s is not set; call _bootstrap_env_root first" "$env_var"
    return 2
  fi

  local cfg_path="${root%/}/${cfg_file}"
  if [[ ! -f "$cfg_path" ]]; then
    _log_err "Config not found: %s" "$cfg_path"
    return 3
  fi

  # Source defensively and then call the target if defined
  _source_or_explain "$cfg_path" || true
  if ! declare -F "$target" >/dev/null 2>&1; then
    _log_err "Target not defined after sourcing config: %s" "$target"
    return 4
  fi

  # Call target
  if ! "$target" "${target_args[@]}"; then
    _log_warn "Target %s returned non-zero" "$target"
    return 5
  fi

  _log_ok "Config %s executed: %s" "$cfg_file" "$target"
  return 0
}
#EOC

#BOP
# !FUNCTION: _resolve_script_dir
# !INTERFACE: _resolve_script_dir <out_var> [bash_source]
# !DESCRIPTION:
#   Resolve the canonical absolute directory of a script, following symlinks.
#   Writes the directory path into <out_var>.
#
# !USAGE:
#   # from inside a script
#   _resolve_script_dir RUN_GSI_DIR "${BASH_SOURCE[0]}"
#   echo "Script dir is $RUN_GSI_DIR"
#
# !ARGUMENTS:
#   <out_var>     Name of variable to receive the resolved directory path.
#   [bash_source] Path to the script (default: "${BASH_SOURCE[0]}").
#
# !NOTES:
#   • Follows symlinks recursively until the real file is found.
#   • Uses pwd -P to resolve .. and symlinked dirs as well.
#   • Safe under set -euo pipefail.
#EOP
#BOC
_resolve_script_dir() {
  local __out="${1:?out var required}"
  local __src="${2:-${BASH_SOURCE[0]}}"

  # Follow symlinks
  while [[ -L "$__src" ]]; do
    local __link
    __link="$(readlink -- "$__src")"
    if [[ "$__link" = /* ]]; then
      __src="$__link"
    else
      __src="$(cd -- "$(dirname -- "$__src")" && cd -- "$(dirname -- "$__link")" && pwd)/$(basename -- "$__link")"
    fi
  done

  # Assign result = directory of resolved file
  local __dir
  __dir="$(cd -- "$(dirname -- "$__src")" && pwd -P)"
  printf -v "$__out" '%s' "$__dir"
}
#EOC

#BOP
# !FUNCTION: detect_hpc_system
# !INTERFACE: detect_hpc_system [-v|--verbose]
# !DESCRIPTION:
#   Identify the HPC system and set global environment flags accordingly.
#
# !USAGE:
#   detect_hpc_system [-v]
#
# !BEHAVIOR:
#   • Reads uname/hostname to infer platform
#   • Sets: hpc_system, hpc_name, is_egeon, is_cray, WRAPPER, LC_ALL
#
# !EXAMPLE:
#   detect_hpc_system -v
#
# !NOTES:
#   Relies on logging helpers: _log_info, _log_warn, _log_err, _log_action
#   Option parsing via _parse_args (optional)
#EOP
#BOC
detect_hpc_system(){
  _with_strict_mode   # enable strict mode only for this function

  # Parse common flags (-v/--verbose etc.)
  _parse_args "$@" 2>/dev/null || true
  set -- "${leftover_args[@]}"

  # Prevent verbosity changes from leaking
  local verbose=$verbose

  local sys_info short_hostname
  sys_info="$(uname -a)"
  short_hostname="$(hostname -s)"

  export is_egeon=false
  export is_cray=false

  if printf '%s' "$sys_info" | grep -qi 'cray_ari_s'; then
    export hpc_system="cray"
    export hpc_name="xc50"
    export is_cray=true
    export WRAPPER="ftn"
    _log_info 'Detected: Cray XC50'
    return 0

  elif printf '%s' "$sys_info" | grep -qi 'egeon'; then
    export hpc_system="linux"
    export hpc_name="egeon"
    export is_egeon=true
    export WRAPPER="mpif90"
    export LC_ALL="en_US.UTF-8"
    _log_info 'Detected: EGEON Cluster'
    return 0

  elif printf '%s' "$short_hostname" | grep -qi 'headnode'; then
    export hpc_system="linux"
    export hpc_name="egeon"
    export is_egeon=true
    export WRAPPER="mpif90"
    export LC_ALL="en_US.UTF-8"
    _log_warn 'Detected: HEADNODE of EGEON (build-only)'
    return 0

  else
    _log_err 'Unknown machine: %s' "$short_hostname"
    _log_action -f '1) Add the machine under etc/mach/'
    _log_action -f '2) Create an entry in copy_fixed_files (etc/smg_setup.sh)'
    return 1
  fi
}
#EOC

#BOP
# !FUNCTION: disable_conda
# !INTERFACE: disable_conda [-v|--verbose]
# !DESCRIPTION:
#   Deactivate a detected Conda environment (if any) to avoid toolchain conflicts,
#   and clean related environment variables.
#
# !USAGE:
#   disable_conda [-v]
#
# !BEHAVIOR:
#   • If CONDA_PREFIX is set:
#       - Try `conda deactivate` (or `source deactivate`)
#       - Unset common Conda-related variables
#   • Otherwise, no-op with OK message
#
# !EXAMPLE:
#   disable_conda -v
#
# !NOTES:
#   Safe to call even if Conda is not installed.
#EOP
#BOC
disable_conda(){
  # Parse common flags (-v/--verbose e etc.)
  _parse_args "$@" 2>/dev/null || true
  set -- "${leftover_args[@]}"

  # Prevent verbosity changes from leaking
  local verbose=$verbose

  if [[ -n "$CONDA_PREFIX" ]]; then
    _log_warn "Conda environment detected: %s" "$CONDA_PREFIX"
    _log_action  "Deactivating Conda..."

    # Try conda/mamba/micromamba deactivation; ignore errors
    if command -v conda >/dev/null 2>&1; then
      conda deactivate 2>/dev/null || source deactivate 2>/dev/null || true
    elif command -v mamba >/dev/null 2>&1; then
      mamba deactivate 2>/dev/null || true
    elif command -v micromamba >/dev/null 2>&1; then
      micromamba deactivate 2>/dev/null || true
    else
      _log_warn "Conda command not found, but CONDA_PREFIX is set."
    fi

    unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_PROMPT_MODIFIER
    _log_ok "Conda has been disabled."
  else
    _log_ok "No active Conda environment detected."
  fi
}
#EOC
