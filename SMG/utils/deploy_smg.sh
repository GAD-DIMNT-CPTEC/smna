#!/usr/bin/env bash
#-----------------------------------------------------------------------------#
#           Group on Data Assimilation Development - GDAD/CPTEC/INPE          #
#-----------------------------------------------------------------------------#
#BOP
# !SCRIPT: deploy_smg.sh
# !DESCRIPTION:
#   Deploy modified SMG files to the remote EGEON server using rsync,
#   preserving directory structure. By default, only files flagged as
#   "modified" by SVN are transferred.
#
# !CALLING SEQUENCE:
#   ./deploy_smg.sh [--dry-run] [--all]
#
# !OPTIONS:
#   --dry-run   Show what would be transferred, without sending files.
#   --all       Send full predefined file list instead of SVN changes only.
#
# !REVISION HISTORY:
#   03 Oct 2025 - J. G. de Mattos - Initial Version
#   03 Oct 2025 - GPT-5 Review    - Error handling, options, logging improved
# !REMARKS:
#   Requires rsync and svn in PATH.
#EOP
#-----------------------------------------------------------------------------#
set -euo pipefail

#------------------------------- config --------------------------------------#
dest="joao.gerd@egeon:/home/joao.gerd/SMNA_TESTAR/SMG"

# Optional static file list (fallback for --all)
static_files=(
  config_smg.ksh
  etc/__helpers__.sh
  etc/__init__.sh
  etc/smg_setup.sh
  run/run_cycle.sh
  run/scripts/gsi_scripts/runGSI_Functions.sh
  run/scripts/runGSI
  run/scripts/run_obsmake.sh
  cptec/bam/run/MODELIN.das
  cptec/bam/run/MODELIN.default
  cptec/bam/run/runModel
  cptec/bam/run/runPre
)

#------------------------------- functions ----------------------------------#
_log() { printf '%s\n' "$*" >&2; }
_log_info() { _log "[INFO] $*"; }
_log_ok()   { _log "[OK]   $*"; }
_log_err()  { _log "[ERROR] $*"; }

usage() {
  cat <<EOF
Usage: $0 [--dry-run] [--all]

Deploy SMG files to remote server.

Options:
  --dry-run   Only show what would be transferred
  --all       Deploy static file list (ignores svn status)
EOF
}

#------------------------------- parse args ---------------------------------#
dry_run=""
use_all=false
for arg in "$@"; do
  case "$arg" in
    --dry-run) dry_run="--dry-run" ;;
    --all)     use_all=true ;;
    -h|--help) usage; exit 0 ;;
    *) _log_err "Unknown option: $arg"; usage; exit 1 ;;
  esac
done

#------------------------------- file list ----------------------------------#
if $use_all; then
  files=("${static_files[@]}")
else
  # take only modified files from SVN
  if ! mapfile -t files < <(svn status | awk '/^M/{print $2}'); then
    _log_err "Failed to run svn status"
    exit 2
  fi
fi

if [[ ${#files[@]} -eq 0 ]]; then
  _log_info "No files to transfer"
  exit 0
fi

#------------------------------- execution ----------------------------------#
_log_info "Starting rsync to: $dest"
rsync -avR $dry_run "${files[@]}" "$dest"
_log_ok "Transfer completed successfully"
#EOC
#-----------------------------------------------------------------------------#

