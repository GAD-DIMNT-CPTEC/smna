#!/usr/bin/env bash
#------------------------------------------------------------------------------#
#BOP
# !SCRIPT: check_paths.sh
# !DESCRIPTION:
#   Validate and report all directory paths defined in a SMG configuration file
#   (e.g., egeon_paths.conf). Expands variables, checks existence, and displays
#   colorized results for easy inspection. Optionally creates missing paths and
#   logs results for audit purposes.
#
# !USAGE:
#   ./check_paths.sh [pathfile] [--create-missing] [--yes] [--log <file>]
#
# !EXAMPLES:
#   ./check_paths.sh etc/mach/egeon_paths.conf
#   ./check_paths.sh etc/mach/egeon_paths.conf --create-missing
#   ./check_paths.sh etc/mach/egeon_paths.conf --create-missing --yes
#   ./check_paths.sh etc/mach/egeon_paths.conf --log ~/checkpaths.log
#
# !RETURNS:
#   Exit code 0 if all paths exist, >0 otherwise.
#
# !NOTES:
#   • Requires Bash 4+
#   • Only lines with pattern "<key> <value>" are parsed
#   • Comments (#...) and blank lines are ignored
#EOP
#------------------------------------------------------------------------------#

set -euo pipefail

#=== Colors ===================================================================#
if [[ -t 1 ]]; then
  C_OK=$'\033[1;32m'   # bold green
  C_WARN=$'\033[1;33m' # bold yellow
  C_ERR=$'\033[1;31m'  # bold red
  C_H=$'\033[1;34m'    # bold blue
  C_RST=$'\033[0m'
else
  C_OK= C_WARN= C_ERR= C_H= C_RST=
fi

#=== Defaults =================================================================#
CREATE_MISSING=false
AUTO_YES=false
LOG_FILE=
PATH_FILE=

#=== Parse arguments ==========================================================#
while [[ $# -gt 0 ]]; do
  case "$1" in
    --create-missing) CREATE_MISSING=true ;;
    --yes) AUTO_YES=true ;;
    --log)
      shift || { echo "${C_ERR}[ERROR]${C_RST} Missing argument for --log"; exit 1; }
      LOG_FILE="$1"
      ;;
    -*)
      echo "${C_ERR}[ERROR]${C_RST} Unknown option: $1" >&2
      exit 1
      ;;
    *) PATH_FILE="$1" ;;
  esac
  shift
done

PATH_FILE="${PATH_FILE:-etc/mach/egeon_paths.conf}"

if [[ ! -f "$PATH_FILE" ]]; then
  echo "${C_ERR}[ERROR]${C_RST} Path file not found: $PATH_FILE" >&2
  exit 1
fi

#=== Initialize log ===========================================================#
if [[ -n "${LOG_FILE}" ]]; then
  mkdir -p "$(dirname "$LOG_FILE")"
  {
    echo "------------------------------------------------------------------"
    echo "# SMG Path Check Log — $(date '+%Y-%m-%d %H:%M:%S')"
    echo "# File: ${PATH_FILE}"
    echo "# Options: create-missing=${CREATE_MISSING}, yes=${AUTO_YES}"
    echo "------------------------------------------------------------------"
  } > "$LOG_FILE"
fi

log() { [[ -n "$LOG_FILE" ]] && echo "$*" >> "$LOG_FILE"; }

#=== Header ===================================================================#
echo "${C_H}#------------------------------------------------------------------------------#${C_RST}"
echo "${C_H}# Validating paths from:${C_RST} $PATH_FILE"
echo "${C_H}#------------------------------------------------------------------------------#${C_RST}"
log "Validating paths from: $PATH_FILE"

#=== Counters =================================================================#
total=0 ok=0 warn=0 fail=0 created=0

#=== Load config ==============================================================#
# shellcheck disable=SC1090
source "$PATH_FILE"

#=== Extract and validate variables ===========================================#
vars=$(grep -E '^[[:space:]]*[a-zA-Z_]+\s+[^\#]+' "$PATH_FILE" \
  | awk '{print $1}' \
  | grep -E '_(smg|bam|gsi|blsdas|obs|sfc|bkg|cptec|util|public|data|home|subt|work)' \
  | sort -u)

printf "%-30s %-70s %s\n" "Variable" "Expanded Path" "Status"
printf "%-30s %-70s %s\n" "--------" "--------------" "------"
log "Variable,ExpandedPath,Status"

for var in $vars; do
  ((total++))
  val=$(eval echo "\${$var:-}")
  if [[ -z "$val" ]]; then
    printf "%-30s %-70s %s\n" "$var" "(empty)" "${C_WARN}[WARN]${C_RST}"
    log "$var,(empty),WARN"
    ((warn++))
    continue
  fi

  if [[ -d "$val" ]]; then
    printf "%-30s %-70s %s\n" "$var" "$val" "${C_OK}[OK]${C_RST}"
    log "$var,$val,OK"
    ((ok++))
  else
    printf "%-30s %-70s %s" "$var" "$val" "${C_ERR}[MISSING]${C_RST}"
    log "$var,$val,MISSING"

    if $CREATE_MISSING; then
      if $AUTO_YES; then
        mkdir -p "$val" && ((created++))
        echo " ${C_OK}[CREATED]${C_RST}"
        log "$var,$val,CREATED"
      else
        echo
        read -rp " → Create this directory? [y/N]: " ans
        if [[ "${ans,,}" == y ]]; then
          mkdir -p "$val" && ((created++))
          echo " ${C_OK}[CREATED]${C_RST}"
          log "$var,$val,CREATED"
        fi
      fi
    else
      echo
    fi
    ((fail++))
  fi
done

#=== Summary ==================================================================#
echo
echo "${C_H}#------------------------------------------------------------------------------#${C_RST}"
printf "Total: %d | ${C_OK}OK:%d${C_RST} | ${C_WARN}WARN:%d${C_RST} | ${C_ERR}MISSING:%d${C_RST}" \
  "$total" "$ok" "$warn" "$fail"
if $CREATE_MISSING; then
  printf " | ${C_OK}CREATED:%d${C_RST}\n" "$created"
else
  echo
fi

log "------------------------------------------------------------------"
log "Summary: Total=$total OK=$ok WARN=$warn MISSING=$fail CREATED=$created"
log "------------------------------------------------------------------"

[[ $fail -eq 0 ]] && exit 0 || exit 2

