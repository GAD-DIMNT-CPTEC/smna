#!/usr/bin/env bash
#-----------------------------------------------------------------------------#
#           Group on Data Assimilation Development - GDAD/CPTEC/INPE          #
#-----------------------------------------------------------------------------#
#BOP
# !SCRIPT: run_obsmake.sh — Stage observation files for a given START_DATE
#
# !DESCRIPTION:
#   Stages the BUFR observation files required by GSI for a specific cycle
#   datetime (START_DATE, format YYYYMMDDHH). The script locates files under
#   ${ncep_ext}/YYYY/MM/DD (default layout used in SMNA) and places them in
#   ${subt_gsi_datain_obs}, using symbolic links by default (copy/hardlink
#   available via --mode).
#
#   Robust, portable and idempotent:
#     • Strict mode and safe IFS.
#     • Validates inputs, directories and environment.
#     • Handles empty sources gracefully and returns non-zero on no files.
#     • Works with large directories (find + while read -d '').
#
# !USAGE:
#   run_obsmake.sh START_DATE [--ncep-root DIR] [--dest DIR]
#                              [--mode link|copy|hardlink]
#                              [--pattern 'gdas.*'] [--dry-run] [-v|-q] [-h]
#
# !ARGUMENTS:
#   START_DATE         Cycle datetime label, e.g., 2015043006 (YYYYMMDDHH).
#
# !OPTIONS:
#   --ncep-root DIR    Root of NCEP external tree (defaults to $ncep_ext).
#   --dest DIR         Destination directory (defaults to $subt_gsi_datain_obs).
#   --mode MODE        How to stage files: link (default), copy, hardlink.
#   --pattern GLOB     Source filename glob (default: 'gdas.*').
#   --dry-run          Print actions without performing them.
#   -v, --verbose      Verbose messages.
#   -q, --quiet        Only errors.
#   -h, --help         Show this help extracted from #BOP/#EOP.
#
# !EXIT CODES:
#   0  Success and >=1 files staged
#   1  Usage / invalid arguments / unexpected error
#   2  No files found to stage
#
# !ENVIRONMENT:
#   ncep_ext               Default source root if --ncep-root not provided.
#   subt_gsi_datain_obs    Default destination if --dest not provided.
#
# !EXAMPLES:
#   ./run_obsmake.sh 2025010106 --mode link
#   ./run_obsmake.sh 2025010106 --ncep-root /data/NCEP --dest /run/gsi/obs --pattern 'gdas.*'
#
# !REVISION HISTORY:
#   2025-09-21  Refactor for robustness, options, and ProTex documentation.
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
export RUN_OBSMAKE_DIR="$(cd -- "$(dirname -- "$__SRC")" && pwd -P)"

# --- load init (which loads helpers) ---
__init_path="${RUN_OBSMAKE_DIR}/__init__.sh"
# shellcheck disable=SC1090
. "$__init_path" || :        # no extra prints; don't break under set -e

usage() {
  sed -n '/^#BOP/,/^#EOP/{/^#BOP/d;/^#EOP/d;p}' "${BASH_SOURCE[0]}"
}

#------------------------------- Defaults ------------------------------------#
verbose=false
quiet=false
dry_run=false
pattern='gdas.*'
mode='link'        # link|copy|hardlink
ncep_root="${ncep_ext:-}"
dest_dir="${subt_gsi_datain_obs:-}"

#----------------------------- Parse options ---------------------------------#
_with_strict_mode

[[ $# -ge 1 ]] || { _log_err "Missing START_DATE."; usage; exit 1; }
START_DATE="$1"; shift || true

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ncep-root) ncep_root="${2:-}"; shift 2 ;;
    --dest)      dest_dir="${2:-}"; shift 2 ;;
    --mode)      mode="${2:-}"; shift 2 ;;
    --pattern)   pattern="${2:-}"; shift 2 ;;
    --dry-run)   dry_run=true; shift ;;
    -v|--verbose) verbose=true; quiet=false; shift ;;
    -q|--quiet)  quiet=true; verbose=false; shift ;;
    -h|--help)   usage; exit 0 ;;
    *) err "Unknown option: $1"; usage; exit 1 ;;
  esac
done

# Try to auto-load SMG config if variables are missing and config is reachable
if [[ -z "${ncep_root}" || -z "${dest_dir}" ]]; then
  if command -v realpath >/dev/null 2>&1; then
    SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"
    RootDir="$(dirname "$SCRIPT_PATH")"
    SMG_ROOT="${SMG_ROOT:-$RootDir}"
    if [[ -f "${SMG_ROOT}/../../config_smg.ksh" ]]; then
      # shellcheck disable=SC1090
      . "${SMG_ROOT}/../../config_smg.ksh" vars_export || true
      # Refill defaults after sourcing
      : "${ncep_root:=${ncep_ext:-}}"
      : "${dest_dir:=${subt_gsi_datain_obs:-}}"
    fi
  fi
fi

#------------------------------ Validation -----------------------------------#
[[ "${#START_DATE}" -eq 10 && "${START_DATE}" =~ ^[0-9]{10}$ ]] \
  || die "START_DATE must be YYYYMMDDHH (10 digits): got '${START_DATE}'"

[[ -n "${ncep_root}" ]] || _die "--ncep-root not set and \$ncep_ext is empty."
[[ -n "${dest_dir}"  ]] || _die "--dest not set and \$subt_gsi_datain_obs is empty."

Y=${START_DATE:0:4}
M=${START_DATE:4:2}
D=${START_DATE:6:2}
#HH=${START_DATE:8:2}  # not used by current layout
src_dir="${ncep_root}/${Y}/${M}/${D}"

[[ -d "${src_dir}" ]] || _die "Source directory not found: ${src_dir}"

# Normalize destination and create it
if [[ ! -d "${dest_dir}" ]]; then
  $dry_run || mkdir -p -- "${dest_dir}"
  _log_info "Created destination: ${dest_dir}"
fi

#----------------------------- Staging files ---------------------------------#
_log_info "Scanning: ${src_dir} (pattern: ${pattern})"
count=0

# Use find -print0 and while-read loop to handle arbitrary filenames safely
while IFS= read -r -d '' path; do
  base=$(basename -- "$path")
  target="${dest_dir}/${base}"

  case "${mode}" in
    link)
      $dry_run || ln -sfn -- "$path" "$target"
      _log_info "ln -sfn -- '$path' '$target'"
      ;;
    hardlink)
      $dry_run || ln -fn -- "$path" "$target"
      _log_info "ln -fn -- '$path' '$target'"
      ;;
    copy)
      $dry_run || cp -pf -- "$path" "$target"
      _log_info "cp -pf -- '$path' '$target'"
      ;;
    *)
      _die "Invalid --mode '${mode}'. Use link|copy|hardlink."
      ;;
  esac
  ((count++))
done < <(find -P "${src_dir}" -type f -size +0c -name "${pattern}" -print0)

if (( count > 0 )); then
  _log_ok "Staged ${count} file(s) into ${dest_dir}."
  exit 0
else
  _log_err "No files found in ${src_dir} matching '${pattern}'."
  exit 2
fi

#EOP

