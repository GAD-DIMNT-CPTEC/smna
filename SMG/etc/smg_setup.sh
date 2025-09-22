#!/bin/bash -x
#-----------------------------------------------------------------------------#
#           Group on Data Assimilation Development - GDAD/CPTEC/INPE          #
#-----------------------------------------------------------------------------#
#BOP
# !SCRIPT: smg_setup.sh
#
# !DESCRIPTION:
#   This script sets up the Global Modeling System (SMG) on the HPC, including
#   exporting environment variables, copying fixed files, configuring the
#   environment, compiling, and running test cases.
#
# !USAGE:
#   ./smg_setup.sh <option>
#
# !OPTIONS:
#   configure        - Create directories, copy fixed files, adjust scripts.
#   compile          - Build BAM, GSI and utilities.
#   copy_ncep_inputs - Copy GDAS/NCEP analysis & SST inputs by cycle.
#   testcase         - Prepare and populate a test case environment.
#   help             - Show this documentation and a command list.
#
# !REVISION HISTORY:
#   20 Dec 2017 - J. G. de Mattos - Initial version
#   07 Feb 2025 - Updates for robustness and efficiency
#   15 Sep 2025 - Help system & descriptions improved; alias; small fixes
#EOP
#-----------------------------------------------------------------------------#

# --- Identify this script path even when sourced from another script ---
# Use a canonical absolute path if available; otherwise, keep the raw path.
if command -v readlink >/dev/null 2>&1; then
  SMG_DOC_FILE="$(readlink -f "${BASH_SOURCE[0]}")"
else
  SMG_DOC_FILE="${BASH_SOURCE[0]}"
fi
export SMG_DOC_FILE  # expose to wrappers or child scripts if needed

# --- Minimum required CMake version ---
export CMAKE_VERSION_MIN="4.1.1"

# --- Default logging verbosity ---
export verbose=false

#-----------------------------------------------------------------------------#
#----------------------- internal helpers (hidden from help) -----------------#
#-----------------------------------------------------------------------------#
#BOP
# !FUNCTION: _log_msg
# !DESCRIPTION:
#   Print a standardized log message with a given level (INFO, OK, WARNING, ACTION, etc.).
#   Messages are printed only when the global variable `verbose` is set to `true`,
#   unless output is explicitly forced with the `-f` flag.
#   Supports printf-style formatting (placeholders like %s, %d).
#   Safe under 'set -u' (defaults to `false` when `verbose` is unset).
#
# !INTERFACE:
#   _log_msg <LEVEL> [-f] <format> [args...]
#
# !EXAMPLES:
#   # Respect global verbosity (will print only if verbose=true)
#   _log_msg INFO "Starting step %s" "$step"
#
#   # Force a single message regardless of `verbose`
#   _log_msg WARNING -f "Low disk space: %s" "$mountpoint"
#
#   # Using numeric formatting
#   _log_msg OK "Built %d target(s) in %0.2f s" "$n_targets" "$elapsed"
#
# !NOTES:
#   • Do not add a trailing newline to <format>; the logger appends it automatically.
#   • `verbose` is expected to be the Bash boolean string `true` or `false`.
#   • Output goes to stdout; if you need stderr routing, adapt the implementation accordingly.
#EOP
#BOC
# Core: supports printf-style formatting and a -f flag to force output
_log_msg() {
  local level="$1"; shift
  local force=false
  if [[ "${1:-}" == "-f" ]]; then force=true; shift; fi

  # default is verbose=false (string "true"/"false")
  local v="${verbose:-false}"

  if $v || $force; then
    # $1 is the format string; the rest are printf args
    printf "[%s] $1\n" "$level" "${@:2}"
  fi
}
#EOC

#BOP
# !FUNCTION: _log_info, _log_ok, _log_err, _log_warning, _log_action, _log_fail
# !DESCRIPTION:
#   Convenience wrappers around `_log_msg` that set the appropriate level.
#   `_log_err` and `_log_fail` always force output (they behave as if `-f` was passed).
#   Other wrappers accept an optional `-f` to force output.
#
# !INTERFACE:
#   _log_info    [-f] <format> [args...]
#   _log_ok      [-f] <format> [args...]
#   _log_warning [-f] <format> [args...]
#   _log_action  [-f] <format> [args...]
#   _log_err          <format> [args...]   # forced output
#   _log_fail         <format> [args...]   # forced output
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
_log_info()    { _log_msg "INFO"    "$@"; }
_log_ok()      { _log_msg "OK"      "$@"; }
_log_warn()    { _log_msg "WARNING" "$@"; }
_log_action()  { _log_msg "ACTION"  "$@"; }

# Errors should always print, regardless of verbose
_log_fail()    { _log_msg "FAIL"  -f "$@"; }
_log_err()     { _log_msg "ERROR" -f "$@"; }
#EOC

#BOP
# !FUNCTION: _parse_args
# !INTERFACE: _parse_args "$@"
# !DESCRIPTION:
#   Parse common options and component selection flags for compilation.
#
#   Verbosity and utilities:
#     -v, --verbose      → verbose=true  (enable logging)
#     -q, --quiet        → verbose=false (disable logging)
#     -y, --yes          → auto_yes=true (auto-confirm prompts)
#     -f, --fix          → do_fix=true   (create dirs/symlinks/copies as needed)
#     --dry-run          → dry_run=true  (skip actual actions)
#     --restore          → do_restore=true (restore from *.bak when supported)
#
#   CMake:
#     -c, --cmake X      → set cmake_version=X
#     --cmake=X          → same as above
#
#   Component selection (compilation):
#     --all              → enable all components (GSI, ANGUPDATE, BAM, inctime)
#     --none             → disable all components
#     --gsi / --no-gsi
#     --ang / --no-ang           (ANGUPDATE)
#     --bam / --no-bam
#     --inctime / --no-inctime
#
#   Notes:
#     1) Default values may be inherited from the environment (compgsi/compang/compbam/compinctime).
#        If unset, defaults to "false".
#     2) Flags are evaluated in order: the last one wins (e.g. --all followed by --no-bam disables BAM only).
#
#   Remaining positional arguments are preserved in leftover_args[].
#
# !EXAMPLE:
#   my_function() {
#     _parse_args "$@"
#     $verbose && _log_info "Components: GSI=$compgsi ANG=$compang BAM=$compbam INC=$compinctime"
#     if $dry_run; then
#       _log_info "DRY-RUN: would run with args: ${leftover_args[*]}"
#     fi
#   }
#EOP
#BOC
_parse_args() {
  # ---- Safe defaults (respect existing env vars, otherwise set to false) ----
  verbose=${verbose:-false}
  auto_yes=${auto_yes:-false}
  dry_run=${dry_run:-false}
  do_restore=${do_restore:-false}
  do_fix=${do_fix:-false}
  cmake_version=${cmake_version:-""}

  # ---- Components: allow override from environment or default to false ----
  compgsi=${compgsi:-false}
  compang=${compang:-false}
  compbam=${compbam:-false}
  compinctime=${compinctime:-false}

  leftover_args=()

  while [[ $# -gt 0 ]]; do
    case "$1" in
      # ---- Verbosity / utilities ----
      -v|--verbose) verbose=true ;;
      -q|--quiet)   verbose=false ;;
      -y|--yes)     auto_yes=true ;;
      -f|--fix)     do_fix=true ;;
      --dry-run)    dry_run=true ;;
      --restore)    do_restore=true ;;

      # ---- CMake options ----
      -c|--cmake)   cmake_version="${2:?--cmake requires a version}"; shift ;;
      --cmake=*)    cmake_version="${1#*=}" ;;

      # ---- Component presets ----
      --all)        compgsi=true; compang=true; compbam=true; compinctime=true ;;
      --none)       compgsi=false; compang=false; compbam=false; compinctime=false ;;

      # ---- Individual component toggles ----
      --gsi)        compgsi=true ;;
      --no-gsi)     compgsi=false ;;
      --ang|--angupdate)        compang=true ;;
      --no-ang|--no-angupdate)  compang=false ;;
      --bam)        compbam=true ;;
      --no-bam)     compbam=false ;;
      --inctime)    compinctime=true ;;
      --no-inctime) compinctime=false ;;

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
    printf "\r[INFO] %s [%s] %3d%% (%d/%d) " "$msg" "$bar" "$percent" "$current" "$total"
}
#EOC

#BOP
# !FUNCTION: _copy_with_progress
# !INTERFACE: _copy_with_progress <arrayname> <src_dir> <dest_dir> [copy|link] [msg]
# !DESCRIPTION:
#   Copy or symlink a list of files from a source directory to a destination,
#   showing a simple progress bar on TTY (fallback mode when rsync isn't used).
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
#EOP
#BOC
_copy_with_progress(){
  local -n _files="$1"; shift
  local src_dir="$1"; shift
  local dest_dir="$1"; shift
  local action="${1:-copy}"; shift || true
  local progress_msg="${1:-Processing files}"

  mkdir -p -- "$dest_dir" || { echo "[ERROR] mkdir -p '$dest_dir' failed"; return 1; }

  local total=${#_files[@]}
  if (( total == 0 )); then
    printf "[INFO] %s [nothing to do]\n" "$progress_msg"
    return 0
  fi

  local count=0 f src dst
  for f in "${_files[@]}"; do
    ((count++))
    src="${src_dir%/}/$f"
    dst="${dest_dir%/}/$f"

    if [[ ! -e "$src" ]]; then
      printf "[WARNING] Missing: %s\n" "$src"
      continue
    fi

    # Log / barra antes de executar (opcional: mover para depois, se preferir)
    if [[ -t 1 ]]; then
      _progress_bar "$count" "$total" "$progress_msg"
    else
      printf "[ACTION] (%d/%d) %s -> %s\n" "$count" "$total" "$src" "$dst"
    fi

    if [[ "$action" == "copy" ]]; then
      cp -pf -- "$src" "$dest_dir"/ || { echo "[ERROR] Failed to copy '$f'"; return 1; }
    else
      ln -sf -- "$src" "$dest_dir"/ || { echo "[ERROR] Failed to link '$f'"; return 1; }
    fi
  done

  # Em TTY, a barra já atualizou; adiciona um OK ao final
  [[ -t 1 ]] && echo " - [OK]"
}
#EOC

#BOP
# !FUNCTION: _copy_with_progress_
# !INTERFACE: _copy_with_progress_ <array_name> <src_dir> <dest_dir> [copy|link] ["Mensagem"]
# !DESCRIPTION:
#   Copy or link a list of files from a source directory to a destination
#   directory, with progress reporting. Supports two modes:
#     • If `rsync` is available → uses `rsync -a --human-readable --info=progress2`
#       for detailed progress (percentage, bytes, ETA).
#     • If `rsync` is not available → falls back to `__copy_with_progress__`
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
#           · Without rsync: falls back to __copy_with_progress__, which draws a bar if TTY.
#   • Reports total processed vs expected at the end.
#
# !NOTES:
#   - Requires `_log_info`, `_log_warn`, `_log_action`, `_log_ok` helpers for logging.
#   - Respects global `dry_run` variable (true/false).
#   - Designed for internal use inside higher-level workflows (e.g., copy_ncep_input).
#   - Fallback function `__copy_with_progress__` must exist for the no-rsync case.
#EOP

_copy_with_progress_() {
  local files_name="$1"; shift                 # nome do array (ex.: files), sem [@]
  local -n _files_ref="$files_name"            # nameref para iterar aqui
  local src_dir="$1"; shift
  local dest_dir="$1"; shift
  local action="${1:-copy}"; shift || true
  local progress_msg="${1:-Copying files}"
  local dry="${dry_run:-false}"

  mkdir -p -- "$dest_dir" || { _log_err "mkdir -p failed: %s" "$dest_dir"; return 1; }

  local total=${#_files_ref[@]}
  if (( total == 0 )); then
    _log_info "%s [nothing to do]" "$progress_msg"
    return 0
  fi

  # Se não houver rsync, delega para a tua função (com manejo de dry-run)
  if ! command -v rsync >/dev/null 2>&1; then
    if $dry; then
      local f
      for f in "${_files_ref[@]}"; do
        _log_info "Dry-run: would %s %s -> %s/" "$action" "${src_dir%/}/$f" "$dest_dir"
      done
      _log_ok "%s: %d/%d processed (dry-run)" "$progress_msg" "$total" "$total"
      return 0
    fi
    # Usa tua barra antiga
    _copy_with_progress "$files_name" "$src_dir" "$dest_dir" "$action" "$progress_msg"
    return $?
  fi

  # Com rsync disponível: usa --info=progress2 (limpo em TTY e em log)
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
          # pré-voo: haverá transferência?
          if rsync -ai --dry-run -- "$src" "$dest_dir"/ | grep -q '^[^\.]'; then
            # sim: copia com nome + progresso + itemize
            rsync -a --human-readable --partial \
                  --info=progress2,name0 \
                  -i --out-format='%i %n' -- "$src" "$dest_dir"/ && ((ok++)) || _log_warn "rsync failed: %s" "$src"
          else
            # não: está igual
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
    local src_dir="$1"
    local dest_dir="$2"
    local progress_msg="${3:-Copying files}"
    local action="${4:-copy}"

    mkdir -p "$dest_dir" || { echo "[ERROR] mkdir -p '$dest_dir' failed"; return 1; }

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
    # Pass the entire current environment (including variables exported by _assign itself)
    local expanded
    expanded="$(printf '%s' "$raw" | envsubst)" || return 1
    export "$key=$expanded"
    return 0
  fi

  # Fallback without envsubst: block command substitution
  # (disallow `cmd`, $(cmd) or $((...)) in configuration values)
  if printf '%s' "$raw" | grep -Eq '(`|\$\(|\$\(\()'; then
    echo "[ERROR] _assign: forbidden command substitution in value for $key" >&2
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
# !FUNCTION: _use_local_cmake
# !INTERFACE: _use_local_cmake [-v] [-c VERSION]
# !DESCRIPTION:
#   Ensure a usable version of CMake is available under `utils/cmake`.
#
# !USAGE:
#   use_local_cmake [-v] [-c VERSION]
#
# !OPTIONS:
#   -v, --verbose       Enable verbose logging
#   -c, --cmake VERSION Specify the required CMake version
#
# !BEHAVIOR:
#   • If a suitable version (>= CMAKE_VERSION_MIN) is already installed → use it
#   • If missing or outdated → download and install into `utils/cmake`
#   • Falls back to system `cmake` if appropriate
#
# !EXAMPLE:
#   use_local_cmake -v -c 3.31.6
#
# !NOTES:
#   Relies on `_parse_args` and logging helpers (`_log_info`, `_log_warn`, `_log_ok`)
#   Downloads prebuilt binaries from Kitware GitHub releases
#EOP
#BOC
_use_local_cmake(){
    # Parse common flags (-v/--verbose e etc.)
    _parse_args "$@" 2>/dev/null || true
    set -- "${leftover_args[@]}"
  
    # Prevent verbosity changes from leaking
    local verbose=$verbose
  
    # Derive directories relative to the script
    local script_path="$(realpath "${BASH_SOURCE[0]}")"
    local script_dir="$(dirname "$script_path")"
    local project_dir="$(dirname "$script_dir")"

    # Target install path for local CMake
    local cmake_dir="$project_dir/utils/cmake"
    local cmake_bin="$cmake_dir/bin/cmake"

    # Minimum required version
    local version_min="${cmake_version:-$CMAKE_VERSION_MIN}"
    local url="https://github.com/Kitware/CMake/releases/download/v$version_min/cmake-$version_min-linux-x86_64.tar.gz"
    
    # Helper: check if version1 >= version2
    _version_ge() { printf '%s\n%s\n' "$2" "$1" | sort -V -C; }

    # --- Step 1: Check if a local CMake already exists ---
    if [[ -x "$cmake_bin" ]]; then
        local v
        v=$("$cmake_bin" --version | head -n1 | awk '{print $3}')
        if _version_ge "$v" "$version_min"; then
            _log_info "Using existing local CMake $v at: $cmake_bin"
        else
            _log_warn -f "Local CMake $v < $version_min, upgrading..."
            rm -rf "$cmake_dir"
        fi
    else
       _log_warn -f "Local CMake not found!"
    fi

    # --- Step 2: Download and install if missing or outdated ---
    if [[ ! -x "$cmake_bin" ]]; then
        _log_info -f "Downloading CMake %s ..." "$version_min"
        mkdir -p "$cmake_dir"
        wget -q --show-progress --progress=bar:force -O /tmp/cmake.tar.gz "$url"
    
        if [[ -t 1 ]]; then
            _log_info -f "Extracting CMake archive to: %s (this may take a while) ..." "$cmake_dir"
            tar --checkpoint=100 \
                --checkpoint-action=exec='echo -n "."' \
                -xzf /tmp/cmake.tar.gz -C "$cmake_dir" --strip-components=1
            echo " "    # garante quebra de linha após os pontos
        else
            _log_info -f "Extracting CMake archive to: %s (this may take a while, please wait) ..." "$cmake_dir"
            tar -xzf /tmp/cmake.tar.gz -C "$cmake_dir" --strip-components=1
            echo " "    # garante quebra de linha após os pontos
        fi
    
        rm -f /tmp/cmake.tar.gz
        _log_ok -f "CMake %s installed at %s" "$version_min" "$cmake_dir"
    fi

    # --- Step 3: Export PATH ---
    export PATH="$cmake_dir/bin:$PATH"
    _log_info -f "Using local CMake at: $cmake_bin"
}
#EOC

#BOP
# !FUNCTION: _list_files_array
# !INTERFACE: _list_files_array VAR_NAME SRC_DIR [find_args...]
# !DESCRIPTION:
#   Populate the array named VAR_NAME with the basenames of files
#   found in SRC_DIR (optionally filtered with extra find arguments).
#
# !USAGE:
#   list_files_array filesArray /path/to/src [-name '*.txt']
#
# !BEHAVIOR:
#   • Clears the target array
#   • Runs `find` with optional filters
#   • Stores only basenames into the array
#
# !EXAMPLE:
#   list_files_array myFiles ./data -name '*.form'
#   echo "${myFiles[@]}"
#
# !NOTES:
#   Requires Bash 4+ (for nameref `local -n`)
#EOP
#BOC
_list_files_array() {
  local __outvar=$1
  local src_dir=$2
  shift 2

  local -n arr="$__outvar"   # reference to the external array
  arr=()                     # reset the array

  while IFS= read -r -d '' f; do
    arr+=( "$(basename "$f")" )
  done < <(find -L "$src_dir" -maxdepth 1 -type f "$@" -print0)
}
#EOC

#-----------------------------------------------------------------------------#
#------------------------------ END internal helpers -------------------------#
#-----------------------------------------------------------------------------#
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
  # Parse common flags (-v/--verbose e etc.)
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
    # Shadow verbose locally to avoid leaking the flag
    _log_err 'Unknown machine: %s' "$short_hostname"
    _log_action -f '1) Add the machine under etc/mach/'
    _log_action -f '2) Create an entry in copy_fixed_files (etc/smg_setup.sh)'
    return 1
  fi
}
#EOC

#BOP
# !FUNCTION: vars_export
# !INTERFACE: vars_export
# !DESCRIPTION:
#   Load and export cluster-specific environment variables from a config file
#   in `etc/mach/${hpc_name}_paths.conf`. Each line has the format:
#     KEY   VALUE
#   Blank lines and lines starting with '#' are ignored. Leading/trailing
#   whitespace is trimmed. Uses _assign to perform variable expansion and export.
#
# !USAGE:
#   vars_export
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
#   vars_export   # sets HOME and subt_smg expanded
#
# !NOTES:
#   • Depends on: _assign
#   • Requires: hpc_name set to match a config file
#EOP
#BOC
vars_export(){
  local confdir filepaths
  confdir="$(dirname -- "${BASH_SOURCE[0]}")/mach"
  filepaths="${confdir}/${hpc_name}_paths.conf"

  [[ -f "$filepaths" ]] || { echo "[ERROR] Missing: $filepaths" >&2; return 1; }

  # Lê linhas "KEY  VALUE"; ignora comentários/linhas vazias
  while IFS= read -r line; do
    [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
    # chave = primeira palavra; valor = resto da linha (trim inicial)
    local key="${line%%[[:space:]]*}"
    local value="${line#"$key"}"
    value="${value#"${value%%[![:space:]]*}"}"   # trim leading spaces
    [[ -z "$key" || -z "$value" ]] && continue
    _assign "$key" "$value"
  done < "$filepaths"
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


#BOP
# !FUNCTION: copy_fixed_files
# !INTERFACE: copy_fixed_files [-v|--verbose]
# !DESCRIPTION:
#   Copy or symlink fixed BAM model files depending on the HPC environment.
#
# !USAGE:
#   copy_fixed_files [-v|--verbose]
#
# !BEHAVIOR:
#   • On Egeon-like hosts (hpc_name=egeon) → copy files with progress
#   • On Cray hosts (hostname starts with 'c') → create symlinks
#   • File lists are discovered automatically using `find`
#
# !EXAMPLE:
#   copy_fixed_files   # will populate datain/dataout/databcs automatically
#
# !NOTES:
#   Depends on helper functions: `_copy_with_progress`, `_copy_dir_with_progress`
#   Requires env vars: public_bam, subt_model_bam, subt_pre_bam
#EOP
#BOC
copy_fixed_files(){
  # Parse common flags (-v/--verbose e etc.)
  _parse_args "$@" 2>/dev/null || true
  set -- "${leftover_args[@]}"

  # Prevent verbosity changes from leaking
  local verbose=$verbose

  # 1) load env variables
  if ! vars_export; then
    _log_err "Failed to load cluster paths via vars_export"
    return 1
  fi

  # --- sanity checks ---
  if [[ -z "$public_bam" || -z "$subt_pre_bam" || -z "$subt_model_bam" ]]; then
    _log_err "Required environment vars missing. Need: public_bam, subt_pre_bam, subt_model_bam, home_pos_bam, subt_pos_bam"
    return 1
  fi

  local base_pre="${public_bam}/PRE/datain"
  if [[ ! -d "$base_pre" ]]; then
    _log_err "Directory not found: %s" "$base_pre"
    return 1
  fi

  # --- gather available years ---
  mapfile -t years < <(find -L "$base_pre" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' 2>/dev/null | sort)
  if (( ${#years[@]} == 0 )); then
    _log_err "No year directories found under %s" "$base_pre"
    return 1
  fi
  # --- end sanity checks ---

  # Lists of fixed files to be copied when on Egeon-like systems
  # MODEL/datain → subt_model_bam/datain
  _list_files_array filesModDataIn "${public_bam}/MODEL/datain"
  
  # PRE/datain (somente HybridLevels.*)
  _list_files_array filesPreDataIn "${public_bam}/PRE/datain" -name 'HybridLevels.*'
  
  # PRE/dataout
  _list_files_array filesPreDataOut "${public_bam}/PRE/dataout"


  # PRE/databcs (*.form & *.bin)
  _list_files_array filesPreDataBC "${public_bam}/PRE/databcs" \( -name '*.bin' -o -name '*.form' \)

  if $is_egeon; then
    # ------------------------------- COPY ----------------------------------
    _copy_with_progress filesModDataIn "${public_bam}/MODEL/datain" "${subt_model_bam}/datain" "copy" "Copying MODEL/datain"
    _copy_with_progress filesPreDataOut "${public_bam}/PRE/dataout"   "${subt_pre_bam}/dataout" "copy" "Copying PRE/dataout"
    _copy_with_progress filesPreDataBC  "${public_bam}/PRE/databcs"   "${subt_pre_bam}/databcs" "copy" "Copying PRE/databcs"
    _copy_with_progress filesPreDataIn  "${public_bam}/PRE/datain"    "${subt_pre_bam}/datain" "copy" "Copying PRE/datain (HybridLevels)"


  elif $is_cray; then
    # ------------------------------- LINK ----------------------------------
    _log_info "Linking fixed files..."

    mkdir -p \
      "${subt_model_bam}/datain" \
      "${subt_pre_bam}/datain" \
      "${subt_pre_bam}/dataout" \
      "${subt_pre_bam}/databcs" \
      "${subt_pre_bam}/dataco2" \
      "${subt_pre_bam}/datasst" \
      "${subt_pre_bam}/dataTop" || { _log_err "mkdir -p for Cray destinations failed"; return 1; }

    # Se desejar manter o comportamento antigo com globs:
    _copy_with_progress filesModDataIn "/cray_home/joao_gerd/BAMFIX/model/datain" "${subt_model_bam}/datain" "link" "Linking MODEL/datain"
    _copy_with_progress filesPreDataOut "/cray_home/joao_gerd/BAMFIX/pre/dataout"  "${subt_pre_bam}/dataout" "link" "Linking PRE/dataout"
    _copy_with_progress filesPreDataBC  "/cray_home/joao_gerd/BAMFIX/pre/databcs"  "${subt_pre_bam}/databcs" "link" "Linking PRE/databcs"
    _copy_with_progress filesPreDataIn  "/cray_home/joao_gerd/BAMFIX/pre/datain"   "${subt_pre_bam}/datain"  "link" "Linking PRE/datain"
    _log_info " Linking completed. [OK]"

  else
     _log_err "Unknown environment: HOSTNAME= %s, hpc_name= %s" "${HOSTNAME}" "${hpc_name}" >&2
     return 1
  fi
}

#EOC

#BOP
# !FUNCTION: copy_ncep_input
# !INTERFACE: copy_ncep_input [--dry-run] [-v|--verbose] [--src-root DIR] [--layout ncep-gfs|pre-year] YYYYMMDDHH
# !DESCRIPTION:
#   Copy GDAS/NCEP inputs (atmospheric analysis and SST) for a given cycle.
#   Layouts supported:
#     • ncep-gfs (default): SRC/YYYY/MM/DD/HH/{gdas…, rtgsst…}
#     • pre-year:          SRC/YYYY/ncep_anl/{gdas…} and SRC/YYYY/sst/{rtgsst…}
#
# !USAGE:
#   copy_ncep_input 2024021000
#   copy_ncep_input --src-root /data/ncep_gfs 2024021000
#   copy_ncep_input --layout pre-year --src-root "${public_bam}/PRE/datain" 2019010100
#
# !OPTIONS:
#   --dry-run              Validate and list expected files without copying
#   -v, --verbose          enable verbose logging
#   --src-root DIR         Override ncep_gfs as the root directory
#   --layout LAYOUT        One of: ncep-gfs (default) | pre-year
#
# !BEHAVIOR:
#   • Resolves source dir:
#       - ncep-gfs → $src_root/YYYY/MM/DD/HH
#       - pre-year → $src_root/YYYY/{ncep_anl|sst}
#   • For GDAS, tries the first existing filename among:
#       <prefix>.T<HH>Z.<suffix>.<YYYYMMDDHH>
#     where prefixes={gdas,gdas1,gblav}, suffixes={SAnl,atmanl.nemsio,atmanl.netcdf}
#   • For SST, tries (in order):
#       rtgssthr_grb_0.083.grib2.<YYYYMMDD>
#       rtgssthr_grb_0.5.grib2.<YYYYMMDD>
#       gdas1.T<HH>Z.sstgrb2.<YYYYMMDDHH>
#       gdas1.T<HH>Z.sstgrb.<YYYYMMDDHH>
#     (checks the SST source dir first, then the GDAS dir as fallback)
#   • If --dry-run → only print files to be copied
#   • Otherwise, copy files into $subt_pre_bam/datain with a progress bar
#
# !NOTES:
#   Requires env var: subt_pre_bam
#   Defaults: src_root=${ncep_gfs}, layout=ncep-gfs
#   Depends on: _parse_args, _log_*, _copy_with_progress
#EOP
#BOC
copy_ncep_input(){
  # Parse common flags (-v/--verbose e etc.)
  _parse_args "$@" 2>/dev/null || true
  set -- "${leftover_args[@]}"

  # Prevent verbosity changes from leaking
  local verbose=$verbose

  local src_root="${ncep_gfs}"
  local layout="ncep-gfs"

  # parse extras for this function
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --src-root) src_root="${2:?--src-root requires DIR}"; shift 2; continue ;;
      --src-root=*) src_root="${1#*=}"; shift; continue ;;
      --layout) layout="${2:?--layout requires value}"; shift 2; continue ;;
      --layout=*) layout="${1#*=}"; shift; continue ;;
      *) break ;;
    esac
  done

  # --- Args ---
  if [[ $# -ne 1 ]]; then
    _log_err "Usage: copy_ncep_input [--src-root DIR] [--layout ncep-gfs|pre-year] YYYYMMDDHH"
    return 1
  fi
  local yyyymmddhh="$1"
  if [[ ! "$yyyymmddhh" =~ ^[0-9]{10}$ ]]; then
    _log_err "Invalid cycle format: %s (expected YYYYMMDDHH)" "$yyyymmddhh"; return 1
  fi
  if [[ -z "$subt_pre_bam" ]]; then
    _log_err "Environment variable not set: subt_pre_bam"; return 1
  fi
  if [[ -z "$src_root" ]]; then
    _log_err "Missing source root (set ncep_gfs or pass --src-root)"; return 1
  fi

  # --- Derive names ---
  local yyyy=${yyyymmddhh:0:4} mm=${yyyymmddhh:4:2} dd=${yyyymmddhh:6:2} hh=${yyyymmddhh:8:2}
  local yyyymmdd=${yyyymmddhh:0:8}

  local src_dir_gdas src_dir_sst
  case "$layout" in
    ncep-gfs)
      src_dir_gdas="${src_root}/${yyyy}/${mm}/${dd}/${hh}"
      src_dir_sst="${src_dir_gdas}"
      ;;
    pre-year)
      src_dir_gdas="${src_root}/${yyyy}/ncep_anl"
      src_dir_sst="${src_root}/${yyyy}/sst"
      ;;
    *)
      _log_err "Unknown layout: %s (use ncep-gfs or pre-year)" "$layout"; return 1
      ;;
  esac

  # --- Check sources ---
  if [[ ! -d "$src_dir_gdas" ]]; then
    _log_err "Source directory not found (GDAS): %s" "$src_dir_gdas"; return 1
  fi
  if [[ ! -d "$src_dir_sst" ]]; then
    _log_err "Source directory not found (SST): %s" "$src_dir_sst"; return 1
  fi

  local dest_dir="${subt_pre_bam}/datain"
  $dry_run || mkdir -p -- "$dest_dir" || { _log_err "Failed to create destination: %s" "$dest_dir"; return 1; }

  # --- GDAS filename candidates (multiple naming variants) ---
  # Build candidate names like: <prefix>.T<HH>Z.<suffix>.<YYYYMMDDHH>
  local pfx=(gdas gdas1 gblav)
  local sfx=(atmanl.netcdf SAnl atmanl.nemsio SAnl)
  local gdas_candidates=()
  local p s
  for p in "${pfx[@]}"; do
    for s in "${sfx[@]}"; do
      gdas_candidates+=( "${p}.T${hh}Z.${s}.${yyyymmddhh}" )
    done
  done
  # Pick the first existing candidate in the GDAS source dir
  local gdas_pick=""
  for cand in "${gdas_candidates[@]}"; do
    if [[ -f "${src_dir_gdas}/${cand}" ]]; then
      gdas_pick="$cand"
      break
    fi
  done

  # --- SST filename candidates (multiple naming variants) ---
  # Prefer daily RTG SST (0.083 then 0.5); fallback to hourly gdas1 SST.
  local sst_candidates=(
    "rtgssthr_grb_0.083.grib2.${yyyymmdd}"
    "rtgssthr_grb_0.5.grib2.${yyyymmdd}"
    "gdas1.T${hh}Z.sstgrb2.${yyyymmddhh}"
    "gdas1.T${hh}Z.sstgrb.${yyyymmddhh}"
  )
  # Try SST dir first, then GDAS dir (handles layout differences)
  local sst_pick="" sst_pick_dir=""
  for cand in "${sst_candidates[@]}"; do
    if [[ -f "${src_dir_sst}/${cand}" ]]; then
      sst_pick="$cand"; sst_pick_dir="$src_dir_sst"; break
    elif [[ -f "${src_dir_gdas}/${cand}" ]]; then
      sst_pick="$cand"; sst_pick_dir="$src_dir_gdas"; break
    fi
  done

  _log_info -f "Preparing NCEP input copy (layout=%s): %s" "$layout" "$yyyymmddhh"
  _log_info "GDAS dir: %s" "$src_dir_gdas"
  _log_info "SST  dir: %s" "$src_dir_sst"
  _log_info "Trying GDAS candidates: %s" "${gdas_candidates[*]}"
  _log_info "Trying SST  candidates: %s" "${sst_candidates[*]}"

  # --- Build file list (warn if missing) ---
  local files=()
  if [[ -n "$gdas_pick" ]]; then
    files+=("$gdas_pick")
  else
    _log_warn "No GDAS file found in %s (tried: %s)" "$src_dir_gdas" "${gdas_candidates[*]}"
  fi
  if [[ -n "$sst_pick" ]]; then
    files+=("$sst_pick")
  else
    _log_warn "No SST file found in %s or %s (tried: %s)" "$src_dir_sst" "$src_dir_gdas" "${sst_candidates[*]}"
  fi
  if (( ${#files[@]} == 0 )); then
    _log_err "No expected NCEP files found for cycle %s (layout=%s)" "$yyyymmddhh" "$layout"
    return 1
  fi

  # --- Copy with progress (two sources possible) ---

  if [[ -n "$gdas_pick" ]]; then
    local g=( "$gdas_pick" )
    _copy_with_progress g "$src_dir_gdas" "$dest_dir" copy "Copying GDAS analysis"
  fi
  
  if [[ -n "$sst_pick" ]]; then
    local s=( "$sst_pick" )
    _copy_with_progress s "$sst_pick_dir" "$dest_dir" copy "Copying SST file"
  fi

  return 0
}
#EOC

#BOP
# !FUNCTION: modify_scripts
# !INTERFACE: modify_scripts [--dry-run] [--restore]
# !DESCRIPTION:
#   Update paths and environment variable references in SMG run scripts,
#   Makefiles, and namelists to point to the configured run-tree locations.
#   If --restore is provided, revert those files from their <file>.bak backups.
#
# !USAGE:
#   modify_scripts
#   modify_scripts --dry-run
#   modify_scripts --restore
#
# !OPTIONS:
#   --dry-run   Show what would be changed/restored, without modifying files
#   --restore   Restore all known targets from their <file>.bak backup
#
# !BEHAVIOR:
#   • Default mode: idempotent edits with one-time backup (<file>.bak)
#   • Restore mode: if <file>.bak exists, replace the edited file with its backup
#   • Skips missing files with a warning
#
# !NOTES:
#   Depends on: _parse_args, _log_info/_log_ok/_log_warn/_log_err/_log_action
#   Requires Bash 4+ (associative arrays). Uses POSIX sed/awk and .bak backups.
#EOP
#BOC
modify_scripts(){
  # Parse common flags (-v/--verbose e etc.)
  _parse_args "$@" 2>/dev/null || true
  set -- "${leftover_args[@]}"

  # Prevent verbosity changes from leaking
  local verbose=$verbose
  local rc=0

  # ------------------- target lists -------------------
  local smg_scripts=(
    "${run_smg}/run_cycle.sh"
    "${scripts_smg}/runGSI"
    "${scripts_smg}/run_model.sh"
    "${scripts_smg}/run_obsmake.sh"
    "${scripts_smg}/run_blsdas.sh"
  )
  local makefiles=(
    "${home_model_bam}/source/Makefile"
    "${home_pos_bam}/source/Makefile"
  )
  local nml="${home_run_bam}/PostGridHistory.nml"
  local envfile="${home_run_bam}/EnvironmentalVariables"

  # Collate all targets to simplify restore loop
  local all_targets=("${smg_scripts[@]}" "${makefiles[@]}" "$nml" "$envfile")

  # ------------------- helpers -------------------
  _ensure_backup_once() {
    local f="$1"
    [[ -f "$f" ]] || return 0
    [[ -f "${f}.bak" ]] && return 0
    if $dry_run; then
      _log_info "Dry-run: would create backup %s" "${f}.bak"
      return 0
    fi
    cp -p -- "$f" "${f}.bak" || { _log_err "Failed to create backup for %s" "$f"; return 1; }
    _log_ok "Backup created: %s.bak" "$f"
  }

  _restore_file() {
    local f="$1" b="${f}.bak"
    if [[ ! -f "$b" ]]; then
      _log_warn -f "No backup found for restore: %s" "$f"
      return 0
    fi
    if $dry_run; then
      _log_info "Dry-run: would restore %s from %s" "$f" "$b"
      return 0
    fi
    # Optionally keep a timestamped safety backup of current file before restore
    if [[ -f "$f" ]]; then
      cp -p -- "$f" "${f}.pre-restore.$(date +%Y%m%d%H%M%S)" >/dev/null 2>&1 || true
    fi
    mv -f -- "$b" "$f" || { _log_err "Failed to restore %s" "$f"; return 1; }
    _log_ok "Restored %s" "$f"
    return 0
  }

  _sed_idempotent_insert_after() {
    local f="$1" pat="$2" ins="$3"
    [[ -f "$f" ]] || { _log_warn "Missing file: %s" "$f"; return 0; }
    grep -Fqx "$ins" "$f" && { _log_info "Already configured: %s" "$f"; return 0; }
    $dry_run && { _log_info "Dry-run: would insert after '%s' in %s → %s" "$pat" "$f" "$ins"; return 0; }
    sed -i -e "/$pat/a\\$ins" "$f" || { _log_err "sed insert failed for %s" "$f"; return 1; }
    _log_ok "Inserted config line in %s" "$f"
  }

  _sed_delete_next_line_after_match() {
    local f="$1" pat="$2"
    [[ -f "$f" ]] || { _log_warn "Missing file: %s" "$f"; return 0; }
    $dry_run && { _log_info "Dry-run: would delete line after '%s' in %s" "$pat" "$f"; return 0; }
    awk -v pat="$pat" '
      BEGIN{del=0}
      {
        if(!del && $0 ~ pat){ print; getline; del=1; next }
        print
      }' "$f" > "${f}.tmp" && mv "${f}.tmp" "$f" || { _log_err "Failed deleting after match in %s" "$f"; return 1; }
    _log_ok "Deleted line after '%s' in %s" "$pat" "$f"
  }

  _sed_replace_or_insert() {
    local f="$1" key="$2" val="$3" anchor="$4"
    [[ -f "$f" ]] || { _log_warn "Missing file: %s" "$f"; return 0; }
    local new="export ${key}=\"${val}\""
    if grep -E -q "^export[[:space:]]+${key}=" "$f"; then
      if $dry_run; then
        _log_info "Dry-run: would replace %s in %s → %s" "$key" "$f" "$new"
        return 0
      fi
      sed -i -E "s|^export[[:space:]]+${key}=.*|${new}|" "$f" || { _log_err "sed replace failed for %s (%s)" "$f" "$key"; return 1; }
      _log_ok "Updated %s in %s" "$key" "$f"
    else
      if $dry_run; then
        _log_info "Dry-run: would insert %s after anchor '%s' in %s" "$new" "$anchor" "$f"
        return 0
      fi
      if grep -Fq "$anchor" "$f"; then
        sed -i "/$anchor/a\\$new" "$f" || { _log_err "sed insert failed for %s (%s)" "$f" "$key"; return 1; }
      else
        printf '\n%s\n' "$new" >> "$f" || { _log_err "append failed for %s" "$f"; return 1; }
      fi
      _log_ok "Inserted %s into %s" "$key" "$f"
    fi
  }

  # ------------------- restore mode -------------------
  if $do_restore; then
    _log_action "Restoring modified files from backups..."
    local f
    for f in "${all_targets[@]}"; do
      [[ -n "$f" ]] || continue
      [[ -e "$f" || -e "${f}.bak" ]] || { _log_warn "Target path not found (skipping): %s" "$f"; continue; }
      _restore_file "$f" || rc=1
    done
    (( rc == 0 )) && _log_ok "Restore completed." || _log_warn "Restore completed with warnings/errors."
    return $rc
  fi

  # ------------------- default modify mode -------------------
  local marker_scripts="# Loading system variables"
  local source_line="source \"${home_smg}/config_smg.ksh\" vars_export"

  # 1) RUN SCRIPTS
  for f in "${smg_scripts[@]}"; do
    [[ -f "$f" ]] || { _log_warn "Missing script: %s" "$f"; continue; }
    _ensure_backup_once "$f" || rc=1
    _sed_delete_next_line_after_match "$f" "$marker_scripts" || rc=1
    _sed_idempotent_insert_after "$f" "$marker_scripts" "$source_line" || rc=1
  done

  # 2) MAKEFILES
  local make_anchor="# Path where the Model executable should be located"
  for f in "${makefiles[@]}"; do
    [[ -f "$f" ]] || { _log_warn "Missing Makefile: %s" "$f"; continue; }
    _ensure_backup_once "$f" || rc=1
    _sed_delete_next_line_after_match "$f" "$make_anchor" || rc=1
    local new_path2="PATH2=\"${f%/*}/exec\""
    _sed_idempotent_insert_after "$f" "$make_anchor" "$new_path2" || rc=1
  done

  # 3) NAMELIST
  if [[ -f "$nml" ]]; then
    _ensure_backup_once "$nml" || rc=1
    if ! $dry_run; then
      awk '!/^DirInPut='\''\// && !/^DirOutPut='\''\// && !/^DirMain='\''\//' "$nml" > "${nml}.tmp" \
        && mv "${nml}.tmp" "$nml" || { _log_err "Failed to sanitize namelist %s" "$nml"; rc=1; }
    else
      _log_info "Dry-run: would sanitize absolute Dir* lines in %s" "$nml"
    fi
    declare -A nml_replacements=(
      ["DirInPut"]="\"${subt_model_bam}/dataout\"/TQ0042L028"
      ["DirOutPut"]="\"${subt_grh_bam}/dataout\"/TQ0042L028"
      ["DirMain"]="\"${subt_bam}\""
    )
    for key in "${!nml_replacements[@]}"; do
      local anchor="!${key}=subt_bam"
      local line="${key}='${nml_replacements[$key]}'"
      if grep -Fqx "$line" "$nml"; then
        _log_info "Already configured: %s in %s" "$key" "$nml"
      else
        if $dry_run; then
          _log_info "Dry-run: would insert %s after '%s' in %s" "$key" "$anchor" "$nml"
        else
          if grep -Fq "$anchor" "$nml"; then
            sed -i "/$anchor/a\\$line" "$nml" || { _log_err "Failed to insert %s in %s" "$key" "$nml"; rc=1; }
          else
            printf '\n%s\n' "$line" >> "$nml" || { _log_err "Append failed for %s" "$nml"; rc=1; }
          fi
          _log_ok "Set %s in %s" "$key" "$nml"
        fi
      fi
    done
  else
    _log_warn "Missing namelist: %s" "$nml"
  fi

  # 4) ENVIRONMENTAL VARIABLES
  if [[ -f "$envfile" ]]; then
    _ensure_backup_once "$envfile" || rc=1
    declare -A env_replacements=(
      ["HOMEBASE"]="${home_bam}"
      ["SUBTBASE"]="${subt_bam}"
      ["WORKBASE"]="${work_bam}"
    )
    local anchor_env="# BAM path in HOME"
    for key in "${!env_replacements[@]}"; do
      _sed_replace_or_insert "$envfile" "$key" "${env_replacements[$key]}" "$anchor_env" || rc=1
    done
  else
    _log_warn "Missing env file: %s" "$envfile"
  fi

  (( rc == 0 )) && _log_ok "Script modifications completed." \
               || _log_warn "Script modifications completed with warnings/errors."
  return $rc
}
#EOC

#BOP
# !FUNCTION: configure
# !INTERFACE: configure [--dry-run] [-y|--yes]
# !DESCRIPTION:
#   Create required directories, ensure local CMake is available, copy fixed
#   files, and adjust SMG scripts/Makefiles with correct paths and variables.
#   Prompts for confirmation unless AUTO_ACCEPT=yes or -y/--yes is provided.
#
# !USAGE:
#   configure
#   configure -y
#   configure --dry-run
#
# !OPTIONS:
#   -y, --yes     Auto-confirm (skip prompt)
#   --dry-run     Validate and list actions, but do not apply changes
#
# !BEHAVIOR:
#   • Load cluster paths with vars_export
#   • Prompt for confirmation (unless auto-accepted)
#   • Create all needed directories
#   • Ensure CMake availability via use_local_cmake
#   • Copy fixed files (copy_fixed_files)
#   • Patch scripts/Makefiles (modify_scripts), if available
#
# !NOTES:
#   Depends on: _parse_args, _log_info/_log_ok/_log_warn/_log_err/_log_action,
#               vars_export, use_local_cmake, copy_fixed_files, modify_scripts
#   Uses env var AUTO_ACCEPT=yes as auto-confirmation.
#EOP
configure(){
  # Parse common flags (-v/--verbose e etc.)
  _parse_args "$@" 2>/dev/null || true
  set -- "${leftover_args[@]}"

  # Prevent verbosity changes from leaking
  local verbose=$verbose

  # 1) load env variables
  if ! vars_export; then
    _log_err "Failed to load cluster paths via vars_export"
    return 1
  fi

  _log_info -f "Configuring SMG..."

  # 2) Confirmação (pula se auto_yes, AUTO_ACCEPT=yes, ou --dry-run)
  if ! $auto_yes && [[ "${AUTO_ACCEPT}" != "yes" ]] && ! $dry_run; then
    if [[ -t 0 ]]; then
      read -r -p "Do you want to continue? (Y/N) " -n 1 REPLY; echo
      [[ "$REPLY" =~ ^[Nn]$ ]] && { _log_info "Aborted by user."; return 0; }
    else
      _log_warn -f "Non-interactive session and no auto-confirm flag; proceeding by default."
    fi
  fi

  # 3) Diretórios necessários
  local dirs=(
    "${subt_dataout}" "${subt_bam}" "${subt_gsi}" "${subt_obs_dataout}" "${subt_run}"
    "${subt_run_gsi}" "${subt_obs_run}" "${subt_pre_bam}/datain" "${subt_pre_bam}/dataout"
    "${subt_pre_bam}/databcs" "${subt_pre_bam}/datasst" "${subt_pre_bam}/dataco2"
    "${subt_pre_bam}/dataTop" "${subt_pre_bam}/exec" "${subt_model_bam}/datain"
    "${subt_model_bam}/dataout" "${subt_model_bam}/exec" "${subt_pos_bam}/datain"
    "${subt_pos_bam}/dataout" "${subt_pos_bam}/exec" "${subt_grh_bam}/datain"
    "${subt_grh_bam}/dataout" "${subt_gsi}/datain" "${subt_gsi}/dataout"
  )

  # 4) Dry-run: apenas lista as ações
  if $dry_run; then
    _log_info "Dry-run: would create the following directories if missing:"
    printf '  - %s\n' "${dirs[@]}"
    _log_info "Dry-run: would run: use_local_cmake ; copy_fixed_files ; modify_scripts (if defined)"
    return 0
  fi

  # 5) Criar diretórios
  _log_action "Creating necessary directories..."
  local d rc=0
  for d in "${dirs[@]}"; do
    if [[ -d "$d" ]]; then
      _log_info "Exists: %s" "$d"
    else
      if mkdir -p -- "$d"; then
        _log_ok "Created: %s" "$d"
      else
        _log_err "Failed to create: %s" "$d"
        rc=1
      fi
    fi
  done
  (( rc != 0 )) && return $rc

  # 6) Garantir CMake local
  if declare -F _use_local_cmake >/dev/null; then
    if ! _use_local_cmake "$@"; then
      _log_err "use_local_cmake failed"
      return 1
    fi
  else
    _log_warn "use_local_cmake is not defined; skipping CMake setup."
  fi

  # 7) Copiar arquivos fixos
  if declare -F copy_fixed_files >/dev/null; then
    if ! copy_fixed_files -v "$@"; then
      _log_err "copy_fixed_files failed"
      return 1
    fi
  else
    _log_warn "copy_fixed_files is not defined; skipping fixed files step."
  fi

  # 8) Ajustar scripts/Makefiles
  if declare -F modify_scripts >/dev/null; then
    if ! modify_scripts "$@"; then
      _log_err "modify_scripts failed"
      return 1
    fi
  else
    _log_warn "modify_scripts is not defined; skipping modifications."
  fi

  _log_ok "Configuration completed successfully."
  return 0
}
#EOC

#BOP
# !FUNCTION: verify_executables
# !INTERFACE: verify_executables [--fix]
# !DESCRIPTION:
#   Verifies whether expected executables were built under $HOME and are present
#   under /mnt/beegfs for runtime. Honors the following booleans (exported):
#     compgsi, compang, compbam, compinctime  (default: false)
#   With --fix, copies missing runtime executables from $HOME → /mnt/beegfs.
#
# !USAGE:
#   verify_executables
#   verify_executables --fix
#
# !OPTIONS:
#   -v, --verbose       Enable verbose logging
#   -f, --fix           creates target dirs and copies HOME → BEEGFS for BUILT_ONLY items
#                       (rsync -a --info=progress2,name0 -i if available; else cp -pf).
#
# !BEHAVIOR:
#   • For each enabled component, checks:
#       - HOME path (build artifact)
#       - BEEGFS path (runtime location)
#   • Status per item:
#       READY        → present on /mnt/beegfs (OK to run)
#       BUILT_ONLY   → present only in HOME (compiled but not copied)
#       RUNTIME_ONLY → present only in /mnt/beegfs (unusual/manual copy)
#       MISSING      → missing in both locations
#
# !RETURNS:
#   0 if everything is READY or RUNTIME_ONLY
#   1 if there are BUILT_ONLY items and --fix was not used
#   2 if there are MISSING items
#
# !NOTES:
#   • Relies on: _log_info, _log_ok, _log_warn, _log_err, _log_action
#   • Executable check uses -x (file exists and is executable).
#EOP
verify_executables() {
  # Parse common flags (-v/--verbose e etc.)
  _parse_args "$@" 2>/dev/null || true
  set -- "${leftover_args[@]}"

  # Prevent verbosity changes from leaking
  local verbose=$verbose

  # Prevent variable changes from leaking
  local do_fix=$do_fix

  # 1) load env variables
  if ! vars_export; then
    _log_err "Failed to load cluster paths via vars_export"
    return 1
  fi

  local ver="${nome_smg}"

  # Base paths (version may be overridden via SMNA_VER)
  local base_home="${home_smg}"
  local base_beeg="${subt_smg}"

  # Flags (booleans as strings)
  local compgsi="${compgsi:-false}"
  local compang="${compang:-false}"
  local compbam="${compbam:-false}"
  local compinctime="${compinctime:-false}"

  # rsync options
  if $verbose; then
    # Verbose mode: show detailed progress, human-readable sizes, and itemized changes
    rsync_opts=(-a --human-readable --info=progress2,name0 -i)
  else
    # Quiet mode: keep it minimal, suppress most of the output
    rsync_opts=(-a -q)
  fi

  # Build the checklist: name|home_path|beegfs_path
  local items=()
  [[ "$compgsi"     == true ]] && items+=("gsi.x|${base_home}/cptec/gsi/build/src/gsi/gsi.x|${base_beeg}/cptec/bin/gsi.x")
  [[ "$compang"     == true ]] && items+=("gsi_angupdate.exe|${base_home}/cptec/gsi/build/src/angle/gsi_angupdate.exe|${base_beeg}/cptec/bin/gsi_angupdate.exe")
  [[ "$compinctime" == true ]] && items+=("inctime|${base_home}/cptec/bin/inctime|${base_beeg}/cptec/bin/inctime")
  [[ "$compbam"     == true ]] && items+=("bam_pre_ParPre_MPI|${base_home}/cptec/bam/pre/build/ParPre_MPI|${base_beeg}/bam/pre/exec/ParPre_MPI")
  [[ "$compbam"     == true ]] && items+=("bam_model_ParModel_MPI|${base_home}/cptec/bam/model/build/ParModel_MPI|${base_beeg}/bam/model/exec/ParModel_MPI")
  [[ "$compbam"     == true ]] && items+=("bam_pos_PostGrib|${base_home}/cptec/bam/pos/build/PostGrib|${base_beeg}/bam/pos/exec/PostGrib")

  if ((${#items[@]}==0)); then
    _log_warn "No items enabled by flags (compgsi/compang/compbam/compinctime)."
    return 0
  fi

  _log_info -f "Checking executables for version %s" "$ver"

  local ready=0 built_only=0 runtime_only=0 missing=0 fixed=0
  local IFS='|' name hpath bpath

  for row in "${items[@]}"; do
    read -r name hpath bpath <<<"$row"
    local have_home=false have_beeg=false
    [[ -x "$hpath" ]] && have_home=true
    [[ -x "$bpath" ]] && have_beeg=true

    if $have_beeg && $have_home; then
      _log_ok "READY        | %s | home=%s | beegfs=%s" "$name" "$hpath" "$bpath"
      ((ready++))
    elif $have_beeg && ! $have_home; then
      _log_info "RUNTIME_ONLY | %s | home=%s | beegfs=%s" "$name" "$hpath" "$bpath"
      ((runtime_only++))
    elif ! $have_beeg && $have_home; then
      _log_warn "BUILT_ONLY   | %s | home=%s | beegfs=%s" "$name" "$hpath" "$bpath"
      ((built_only++))

      if $do_fix; then
        _log_action "Copying to runtime (beegfs) → %s" "$bpath"
        mkdir -p -- "$(dirname "$bpath")" || { _log_err "mkdir failed: %s" "$(dirname "$bpath")"; continue; }
        if command -v rsync >/dev/null 2>&1; then
          if rsync "${rsync_opts[@]}" -- "$hpath" "$(dirname "$bpath")/"; then
            _log_ok "Copied via rsync: %s → %s" "$hpath" "$bpath"
            ((fixed++))
          else
            _log_err "rsync failed: %s → %s" "$hpath" "$bpath"
          fi
        else
          if cp -pf -- "$hpath" "$bpath"; then
            _log_ok "Copied via cp: %s → %s" "$hpath" "$bpath"
            ((fixed++))
          else
            _log_err "cp failed: %s → %s" "$hpath" "$bpath"
          fi
        fi
      fi
    else
      _log_err "MISSING      | %s | home=%s | beegfs=%s" "$name" "$hpath" "$bpath"
      ((missing++))
    fi
  done

  _log_info -f "Summary: READY=%d, BUILT_ONLY=%d, RUNTIME_ONLY=%d, MISSING=%d, FIXED=%d" \
            "$ready" "$built_only" "$runtime_only" "$missing" "$fixed"

  # Return codes as specified
  if ((missing>0)); then
    _log_warn "There are MISSING items: build/copy required before running on beegfs."
    return 2
  fi
  if ((built_only>0)) && ! $do_fix; then
    _log_warn "There are BUILT_ONLY items. Run: verify_executables --fix"
    return 1
  fi
  return 0
}
#EOC

#BOP
# !FUNCTION: compile
# !INTERFACE: compile [--dry-run] [-v|--verbose]
# !DESCRIPTION:
#   Build the SMG stack: GSI (+angupdate), BAM and utilities. Verifies outputs
#   and copies executables to ${home_cptec}/bin when successful.
#
# !USAGE:
#   compile
#   compile --dry-run
#
# !OPTIONS:
#   --dry-run   Show planned actions and exit without building
#
# !BEHAVIOR:
#   • Load cluster paths (vars_export)
#   • Ensure ${home_cptec}/bin exists
#   • Conditionally build components based on boolean flags:
#       compgsi/compang/compbam/compinctime (true|false)
#   • Log to component-specific compile.log when applicable
#   • Verify expected artifacts and copy them to ${home_cptec}/bin
#
# !NOTES:
#   Depends on: _parse_args, _log_info/_log_ok/_log_warn/_log_err/_log_action, vars_export, verify_executables
#   Expects env flags (boolean): compgsi, compang, compbam, compinctime (default: false)
#   Expects env: compiler, hpc_name, home_* paths
#EOP
#BOC
compile(){
  # Parse common flags (-v/--verbose e etc.)
  _parse_args "$@" 2>/dev/null || true
  set -- "${leftover_args[@]}"

  # Prevent verbosity changes from leaking
  local verbose=$verbose

  local rc=0

  # --- Defaults for boolean flags if unset ---
  : "${compgsi:=false}"
  : "${compang:=false}"
  : "${compbam:=false}"
  : "${compinctime:=false}"

  # --- Load paths ---
  if ! vars_export; then
    _log_err "Failed to load cluster paths via vars_export"
    return 1
  fi

  # --- Ensure bin dir ---
  if [[ -z "$home_cptec" ]]; then
    _log_err "home_cptec is not set"
    return 1
  fi
  if $dry_run; then
    _log_info "Dry-run: would ensure directory: %s/bin" "$home_cptec"
  else
    mkdir -p -- "${home_cptec}/bin" || { local verbose=true; _log_err "Cannot create %s/bin" "$home_cptec"; return 1; }
  fi

  _log_info "Starting compilation..."

  # --------------- GSI ----------------
  if $compgsi; then
    if [[ -z "$home_gsi" || -z "$home_gsi_src" || -z "$home_gsi_bin" ]]; then
      _log_err "GSI paths not set (home_gsi/home_gsi_src/home_gsi_bin)"; rc=1
    else
      _log_action -f "Compiling GSI …"
      _log_info   "PATH %s" "$home_gsi"
      if $dry_run; then
        _log_info "Dry-run: would run: (cd %s && . ./compile.sh -C %s | tee compile.log)" "$home_gsi" "$compiler"
      else
        pushd "$home_gsi" >/dev/null || { local verbose=true; _log_err "cd failed: %s" "$home_gsi"; return 1; }
        . ./compile.sh -C "${compiler}" 2>&1 | tee "${home_gsi}/compile.log"
        popd >/dev/null || true
        if [[ ! -e "${home_gsi_src}/gsi.x" ]]; then
          _log_err "GSI compilation failed. See %s" "${home_gsi}/compile.log"; rc=1
        else
          $dry_run && _log_info "Dry-run: would copy gsi.x to %s and %s" "$home_gsi_bin" "${home_cptec}/bin" \
                  || { cp -pvf -- "${home_gsi_src}/gsi.x" "${home_gsi_bin}/" && cp -pvf -- "${home_gsi_src}/gsi.x" "${home_cptec}/bin/"; } || rc=1
        fi
      fi
    fi
  else
    _log_info -f "Skipping GSI (compgsi=false)"
  fi

  # -------------- ANGUPDATE --------------
  if $compang; then
    if [[ -z "$home_gsi" || -z "$hpc_name" || -z "$compiler" ]]; then
      _log_err "ANGUPDATE needs home_gsi/hpc_name/compiler"; rc=1
    else
      _log_action -f "Compiling GSI bias correction (global_angupdate) …"
      if $dry_run; then
        _log_info "Dry-run: would source %s/env.sh %s %s" "$home_gsi" "$hpc_name" "$compiler"
        _log_info "Dry-run: would (cd %s/util/global_angupdate && ln -sf Makefile.conf.%s-%s Makefile.conf && make clean && make)" "$home_gsi" "$hpc_name" "$compiler"
        _log_info "Dry-run: would copy global_angupdate to %s/bin" "$home_cptec"
      else
        if [[ -f "${home_gsi}/env.sh" ]]; then
          # shellcheck disable=SC1090
          source "${home_gsi}/env.sh" "${hpc_name}" "${compiler}"
        else
          _log_warn "env.sh not found under %s; proceeding without." "$home_gsi"
        fi
        pushd "${home_gsi}/util/global_angupdate" >/dev/null || { _log_err "cd failed: %s" "${home_gsi}/util/global_angupdate"; return 1; }
        ln -sf "Makefile.conf.${hpc_name}-${compiler}" Makefile.conf
        make -f Makefile clean && make -f Makefile
        popd >/dev/null || true

        if [[ ! -e "${home_gsi}/util/global_angupdate/global_angupdate" ]]; then
          _log_err "ANGUPDATE compilation failed."
          rc=1
        else
          cp -pvf -- "${home_gsi}/util/global_angupdate/global_angupdate" "${home_cptec}/bin/" || rc=1
        fi
      fi
    fi
  else
    _log_info -f "Skipping ANGUPDATE (compang=false)"
  fi

  # --------------- BAM ----------------
  if $compbam; then
    if [[ -z "$home_bam" ]]; then
      _log_err "home_bam is not set"; rc=1
    else
      _log_action -f "Compiling BAM …"
      _log_info   "PATH %s" "$home_bam"
      if $dry_run; then
        _log_info "Dry-run: would run: (cd %s && . ./compile.sh | tee compile.log)" "$home_bam"
      else
        pushd "$home_bam" >/dev/null || { _log_err "cd failed: %s" "$home_bam"; return 1; }
        . ./compile.sh 2>&1 | tee "${home_bam}/compile.log"
        popd >/dev/null || true
        # (Adicione validação/cópias de artefatos do BAM se necessário)
      fi
    fi
  else
    _log_info -f "Skipping BAM (compbam=false)"
  fi

  # ----------- INCTIME utility ----------
  if $compinctime; then
    if [[ -z "$util_inctime" ]]; then
      _log_err "util_inctime is not set"; rc=1
    else
      _log_action -f "Compiling inctime utility …"
      _log_info   "PATH %s" "$util_inctime"
      if $dry_run; then
        _log_info "Dry-run: would module swap gnu9/9.4.0 intel/2021.4.0"
        _log_info "Dry-run: would (cd %s/src && ARCH=Darwin_intel make)" "$util_inctime"
        _log_info "Dry-run: would copy inctime to %s/bin" "$home_cptec"
      else
        command -v module >/dev/null 2>&1 && module swap gnu9/9.4.0 intel/2021.4.0 || _log_warn "'module' command not available; skipping swap"
        pushd "${util_inctime}/src" >/dev/null || { _log_err "cd failed: %s" "${util_inctime}/src"; return 1; }
        export ARCH=Darwin_intel
        make
        popd >/dev/null || true
        if [[ ! -e "${util_inctime}/src/inctime" ]]; then
          _log_err "inctime compilation failed."
          rc=1
        else
          cp -pvf -- "${util_inctime}/src/inctime" "${home_cptec}/bin/" || rc=1
        fi
      fi
    fi
  else
    _log_info -f "Skipping inctime (compinctime=false)"
  fi

  # ----------- VERIFY IF EXECUTABLES ----------
  if $dry_run; then
     _log_info -f "Dry-run: Checking executables for version"
  else
     verify_executables --fix
  fi

  (( rc == 0 )) && _log_ok "Compilation completed successfully." \
               || _log_warn "Compilation finished with errors."
  return $rc
}
#EOC

#BOP
# !FUNCTION: testcase
# !INTERFACE: testcase [--dry-run] [--year YYYY]
# !DESCRIPTION:
#   Prepare a sample testcase by selecting a year under $public_bam/PRE/datain,
#   then copy only the NCEP input files for one cycle using copy_ncep_input
#   with the 'pre-year' layout. Fixed files remain the responsibility of
#   copy_fixed_files.
#
# !USAGE:
#   testcase
#   testcase --year 2019
#
# !NOTES:
#   Requires: public_bam, subt_pre_bam
#   Depends on: _parse_args, _log_*, copy_ncep_input
#EOP
#BOC
testcase(){
  # Parse common flags (-v/--verbose e etc.)
  _parse_args "$@" 2>/dev/null || true
  set -- "${leftover_args[@]}"

  # Prevent verbosity changes from leaking
  local verbose=$verbose

  local chosen_year=""
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --year) chosen_year="${2:?--year requires YYYY}"; shift 2; continue ;;
      --year=*) chosen_year="${1#*=}"; shift; continue ;;
      *) break ;;
    esac
  done
  # Monte a flag condicionalmente
  local dr=()
  $dry_run && dr+=(--dry-run)

  if [[ -z "$public_bam" || -z "$subt_pre_bam" ]]; then
    _log_err "Missing env: public_bam, subt_pre_bam"; return 1
  fi

  local base_pre="${public_bam}/PRE/datain"
  if [[ ! -d "$base_pre" ]]; then
    _log_err "Directory not found: %s" "$base_pre"; return 1
  fi

  # listar anos
  mapfile -t years < <(find -L "$base_pre" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' 2>/dev/null | sort)
  (( ${#years[@]} > 0 )) || { _log_err "No year directories under %s" "$base_pre"; return 1; }

  if [[ -z "$chosen_year" ]]; then
    if [[ -t 0 ]]; then
      echo -e "\e[34;1m Select a testcase year:\e[0m"
      local i=1; for y in "${years[@]}"; do printf "\e[31;1m[%2d]\e[0m \e[32;1m%s\e[0m\n" "$i" "$y"; ((i++)); done
      local answer; read -r answer
      [[ "$answer" =~ ^[0-9]+$ ]] && (( answer>=1 && answer<=${#years[@]} )) || { _log_err "Invalid selection"; return 1; }
      chosen_year="${years[$((answer-1))]}"
    else
      chosen_year="${years[-1]}"; _log_info "No TTY; selecting latest year: %s" "$chosen_year"
    fi
  else
    [[ -d "${base_pre}/${chosen_year}" ]] || { _log_err "Year not found: %s" "${base_pre}/${chosen_year}"; return 1; }
  fi

  # escolher um ciclo a partir de um arquivo GDAS de análise (suporta várias variantes)
  local src_anl_dir="${base_pre}/${chosen_year}/ncep_anl"
  if [[ ! -d "$src_anl_dir" ]]; then
    _log_err "Missing directory: %s" "$src_anl_dir"; return 1
  fi

  local cycle="" picked_file=""
  # Liste nomes e escolha o primeiro que combine com os padrões de análise
  while IFS= read -r fname; do
    # precisa terminar com .YYYYMMDDHH
    if [[ "$fname" =~ \.([0-9]{10})$ ]]; then
      local cand_cycle="${BASH_REMATCH[1]}"
      # precisa conter T??Z e ser um arquivo de análise (não SST)
      if [[ "$fname" =~ T[0-9]{2}Z ]] && [[ "$fname" =~ (SAnl|atmanl\.nemsio|atmanl\.netcdf) ]]; then
        cycle="$cand_cycle"
        picked_file="$fname"
        break
      fi
    fi
  done < <(find -L "$src_anl_dir" -maxdepth 1 -type f -printf '%f\n' | sort)

  if [[ -z "$cycle" ]]; then
    _log_err "No GDAS analysis file found under %s (looked for *T??Z.*{SAnl,atmanl.nemsio,atmanl.netcdf}.*.<YYYYMMDDHH>)" "$src_anl_dir"
    return 1
  fi

  _log_action -f "Preparing testcase using year=%s cycle=%s" "$chosen_year" "$cycle"

  # chama a função genérica usando o layout pre-year
  copy_ncep_input "${dr[@]}" -- --layout pre-year --src-root "$base_pre" "$cycle"
  
  local rc=$?
  (( rc==0 )) && _log_ok "Testcase ready: year=%s cycle=%s" "$chosen_year" "$cycle"
  return $rc
}
#EOC


# ================= FRIENDLY HELP + SUGGESTIONS + DESCRIPTIONS ===================
#BOP
# !FUNCTION: list_funcs
# !INTERFACE: list_funcs [FILE]
# !DESCRIPTION:
#   List function names defined in a Bash file (defaults to this script).
#
# !USAGE:
#   list_funcs
#   list_funcs path/to/other.sh
#
# !BEHAVIOR:
#   • Uses awk to match Bash function headers like:  name() { … }
#   • Prints one function name per line (no leading underscores filtered here)
#
# !EXAMPLE:
#   list_funcs | sort
#
# !NOTES:
#   Defaults to $SMG_DOC_FILE when FILE is not provided.
#EOP
#BOC
list_funcs(){
  local file="${1:-${SMG_DOC_FILE}}"
  awk '
  /^[ \t]*[a-zA-Z_][a-zA-Z0-9_]*[ \t]*\(\)[ \t]*\{/ {
      name=$1
      sub(/[ \t]*\(\)[ \t]*\{[ \t]*$/, "", name)
      print name
  }' "$file"
}
#EOC

#BOP
# !FUNCTION: get_doc
# !INTERFACE: get_doc <FUNCTION_NAME> [FILE]
# !DESCRIPTION:
#   Extract the DESCRIPTION block of a function from ProTex headers (#BOP/#EOP).
#
# !USAGE:
#   get_doc configure
#   get_doc copy_fixed_files ./smg_setup.sh
#
# !BEHAVIOR:
#   • Scans for a #BOP…#EOP block whose !FUNCTION or !INTERFACE matches the name
#   • Captures lines under !DESCRIPTION: (until the next !-tag or #EOP)
#   • Strips the leading ‘#’ and returns the plain text
#
# !EXAMPLE:
#   desc=$(get_doc compile) && printf "%s\n" "$desc"
#
# !NOTES:
#   Expected ProTex format:
#     #BOP
#     # !FUNCTION: name
#     # !INTERFACE: name
#     # !DESCRIPTION:
#     #   lines…
#     #EOP
#   Defaults to $SMG_DOC_FILE when FILE is not provided.
#EOP
#BOC
get_doc(){
  local funcname="$1"
  local file="${2:-$SMG_DOC_FILE}"
  awk -v f="$funcname" '
  /^#BOP/ { inblk=1; matchok=0; cap=0; desc=""; next }
  /^#EOP/ {
      if (inblk && matchok) {
          sub(/[ \t\r\n]+$/, "", desc)
          print desc
      }
      inblk=0; matchok=0; cap=0; desc=""
      next
  }
  inblk {
      line=$0
      if (line ~ /^#[ \t]*!FUNCTION:[ \t]*/) {
          tmp=line; sub(/^#[ \t]*!FUNCTION:[ \t]*/, "", tmp)
          if (tmp==f) matchok=1
          next
      }
      if (line ~ /^#[ \t]*!INTERFACE:[ \t]*/) {
          tmp=line; sub(/^#[ \t]*!INTERFACE:[ \t]*/, "", tmp)
          if (tmp==f) matchok=1
          next
      }
      if (matchok && line ~ /^#[ \t]*!DESCRIPTION:/) {
          cap=1
          tmp=line; sub(/^#[ \t]*!DESCRIPTION:[ \t]*/, "", tmp)
          if (length(tmp)) desc = desc tmp "\n"
          next
      }
      if (matchok && cap) {
          if (line ~ /^#[ \t]*!/) { cap=0; next }
          if (line ~ /^#/) {
              tmp=line
              sub(/^#[ \t]?/, "", tmp)
              desc = desc tmp "\n"
          }
      }
  }
  ' "$file"
}
#EOC

#BOP
# !FUNCTION: get_help_block
# !INTERFACE: get_help_block <FUNCTION_NAME> [FILE]
# !DESCRIPTION:
#   Extract the full ProTex help block (#BOP … #EOP) for a given function,
#   stripping leading '#' so the block prints cleanly.
#
# !USAGE:
#   get_help_block configure
#   get_help_block copy_fixed_files ./smg_setup.sh
#
# !BEHAVIOR:
#   • Scans the file for a #BOP…#EOP section
#   • Matches the function name by !FUNCTION or !INTERFACE
#   • Collects all commented lines inside the block and removes leading “#”
#   • Prints the entire block to stdout
#
# !EXAMPLE:
#   ./smg_setup.sh help configure   # internally calls get_help_block
#
# !NOTES:
#   Defaults to $SMG_DOC_FILE if FILE is not given.
#   Used by show_help_func to render detailed help for one function.
#EOP
#BOC
get_help_block(){
  local funcname="$1"
  local file="${2:-$SMG_DOC_FILE}"
  awk -v f="$funcname" '
  /^#BOP/ { inblk=1; matchok=0; out=""; next }
  /^#EOP/ {
      if (inblk && matchok) {
          sub(/[[:space:]\r\n]+$/, "", out)
          print out
      }
      inblk=0; matchok=0; out=""
      next
  }
  inblk {
      line=$0
      # detect function match via !FUNCTION or !INTERFACE
      if (line ~ /^#[[:space:]]*!FUNCTION:[[:space:]]*/) {
          tmp=line; sub(/^#[[:space:]]*!FUNCTION:[[:space:]]*/, "", tmp)
          if (tmp==f) matchok=1
      }
      if (line ~ /^#[[:space:]]*!INTERFACE:[[:space:]]*/) {
          tmp=line; sub(/^#[[:space:]]*!INTERFACE:[[:space:]]*/, "", tmp)
          if (tmp==f) matchok=1
      }
      if (matchok && line ~ /^#/) {
          sub(/^#[[:space:]]?/, "", line)
          out = out line "\n"
      }
  }
  ' "$file"
}
#EOC

#BOP
# !FUNCTION: show_help_func
# !INTERFACE: show_help_func <FUNCTION_NAME> [FILE]
# !DESCRIPTION:
#   Print a detailed help message for a single function, showing the full
#   ProTex block (#BOP…#EOP). Falls back to DESCRIPTION only if no block found.
#
# !USAGE:
#   show_help_func configure
#   show_help_func compile ./smg_setup.sh
#
# !BEHAVIOR:
#   • Verifies that the function exists in the file (via list_funcs)
#   • Calls get_help_block to fetch the full #BOP/#EOP block
#   • If block is missing, prints only the DESCRIPTION (from get_doc)
#   • Wraps output with separators for readability
#
# !EXAMPLE:
#   ./smg_setup.sh help copy_fixed_files
#
# !NOTES:
#   Intended for CLI dispatcher: “help <function>”
#   Complements show_help (list of all functions with short description).
#EOP

show_help_func(){
  local name="$1"
  local file="${2:-$SMG_DOC_FILE}"
  if [[ -z "$name" ]]; then
    echo "[ERROR] Usage: help <function>"
    return 1
  fi
  # verifica se a função existe no arquivo
  if ! list_funcs "$file" | grep -qx "$name"; then
    echo "[ERROR] Unknown function: $name"
    return 1
  fi
  local block
  block="$(get_help_block "$name" "$file")"
  if [[ -z "$block" ]]; then
    echo "[WARN] No #BOP/#EOP help block found for: $name"
    # fallback: mostra só a DESCRIPTION se existir
    local desc; desc="$(get_doc "$name" "$file")"
    [[ -n "$desc" ]] && { echo "$name -"; printf "%s\n" "$desc"; }
    return 0
  fi
  # impressão “bonita”
  echo "-------------------------------------------------------------------------------"
  echo "$block"
  echo "-------------------------------------------------------------------------------"
}

#BOP
# !FUNCTION: show_help
# !INTERFACE: show_help [FILE]
# !DESCRIPTION:
#   Print a CLI help message listing available commands and their descriptions.
#
# !USAGE:
#   show_help
#   show_help ./smg_setup.sh
#
# !BEHAVIOR:
#   • Prints Usage and Available commands sections
#   • Uses list_funcs to enumerate functions in FILE (default: $SMG_DOC_FILE)
#   • Skips internal helpers (show_help, get_doc, list_funcs) and names starting with “_”
#   • For each function, fetches its DESCRIPTION with get_doc; falls back to “(no description)”
#
# !EXAMPLE:
#   ./smg_setup.sh help
#
# !NOTES:
#   Intended as the user-facing help for this script.
#EOP
#BOC
show_help(){
  local file="${1:-$SMG_DOC_FILE}"
  echo "Usage:"
  echo "./$(basename "$file") <option>"
  echo
  echo "Available commands:"
  echo

  while read -r fname; do
    case "$fname" in
      help|show_help|get_doc|list_funcs|get_help_block|show_help_func) continue ;;  # skip internal helpers
      _*) continue ;;                            # skip functions starting with "_"
    esac
    local desc
    desc=$(get_doc "$fname" "$file")
    if [ -z "$desc" ]; then
      desc="(no description)"
    fi
    echo "$fname -"
    printf "%s\n" "$desc" | sed 's/^/\t/'
    echo
  done < <(list_funcs "$file")
}
#EOC

# --------------------------- CLI dispatcher -------------------------------------
#BOP
# !SECTION: CLI dispatcher
# !INTERFACE: ./$(basename "$SMG_DOC_FILE") <command> [args...]
# !DESCRIPTION:
#   Entry point that routes CLI calls to functions defined in this script.
#   Provides a ‘help’ command, validates environment (hpc_name), and invokes
#   detect_hpc_system on demand.
#
# !USAGE:
#   ./smg_setup.sh help
#   ./smg_setup.sh configure -y
#   ./smg_setup.sh copy_ncep_input --dry-run 2025010100
#
# !BEHAVIOR:
#   • If no command or “help” → prints help and exits 0
#   • Ensures hpc_name is set; if not, calls detect_hpc_system and aborts on failure
#   • If <command> maps to a declared function → executes it with remaining args
#   • Otherwise, prints an error, shows help, and exits 1
#
# !EXAMPLE:
#   ./smg_setup.sh modify_scripts --dry-run
#
# !NOTES:
#   Exits the shell (this section is the executable CLI, not a function).
#EOP
#BOC

help(){
  local target="$1"

  if [[ -n "$target" ]]; then
    # era "$2"; correto é "$1"
    show_help_func "$target" "$SMG_DOC_FILE"
  else
    show_help "$SMG_DOC_FILE"
    echo "Tip: use 'help <function>' for detailed help of a single command."
  fi

  # Se executado como script, encerra; se sourceado, apenas retorna
  if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    exit 0
  else
    return 0
  fi
}

# Ensure hpc_name is defined before using it
if [[ -z "${hpc_name:-}" ]]; then
    _log_warn -f "hpc_name was not set before loading smg_setup.sh!"
    
    # Try to detect HPC system automatically
    detect_hpc_system
    exit_code=$?
    
    if [[ $exit_code -ne 0 ]]; then
        _log_err -f "detect_hpc_system failed (exit code %d). Aborting." "$exit_code"
        exit "$exit_code"
    fi

    _log_info -f "hpc_name set to: %s" "$hpc_name"
fi

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
#EOC
# ============================================================================


