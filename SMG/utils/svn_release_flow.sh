#!/usr/bin/env bash
#-----------------------------------------------------------------------------#
#           Group on Data Assimilation Development - GDAD/CPTEC/INPE          #
#-----------------------------------------------------------------------------#
#BOP
# !SCRIPT: svn_release_flow.sh
# !DESCRIPTION:
#   Given an existing SVN branch, create:
#     1) a tag from that branch,
#     2) the trunk (if missing),
#     3) optionally a new improvements branch from trunk,
#     4) optionally close (delete) the original branch.
#
#   Layout assumed: <REPO_ROOT>/{trunk,branches,tags}
#
# !CALLING SEQUENCE:
#   ./svn_release_flow.sh \
#     --repo https://svn.server/path/project \
#     --source-branch SMNA_v3.0.x \
#     --tag SMNA_v3.0.0 \
#     [--impr-branch SMNA_v3.0.x-improvements] \
#     [--trunk-subpath SMNA/SMG] \
#     [--message "Release SMNA_v3.0.0"] \
#     [--overwrite-trunk] [--dry-run] [--yes]
#
# !BEHAVIOR:
#   - --impr-branch is OPTIONAL.
#   - If --impr-branch is NOT provided:
#       * improvements branch is NOT created
#       * source branch is NOT closed
#       * warnings are emitted
#
# !EXIT CODES:
#   0 OK, 1 usage, 2 config/missing, 3 svn error
#EOP
#-----------------------------------------------------------------------------#
set -euo pipefail

#------------------------------- logging -------------------------------------#
_log()  { printf '%s\n' "$*" >&2; }
_info() { _log "[INFO] $*"; }
_ok()   { _log "[OK]   $*"; }
_warn() { _log "[WARN] $*"; }
_err()  { _log "[ERROR] $*"; }
_die()  { _err "$*"; exit 3; }

#------------------------------- usage ---------------------------------------#
usage() {
  sed -n '1,/^#EOP/p' "$0" | sed -n 's/^# \{0,1\}//p'
  cat <<'EOF'

Usage:
  svn_release_flow.sh --repo URL --source-branch NAME --tag NAME
                      [--impr-branch NAME]
                      [--trunk-subpath SUBDIR]
                      [--message MSG]
                      [--overwrite-trunk]
                      [--dry-run] [--yes]

EOF
}

#------------------------------- args ----------------------------------------#
repo=""
src_branch=""
tag_name=""
impr_branch=""
trunk_subpath=""
message="Automated release flow"
dry_run=false
yes_mode=false
overwrite_trunk=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo)            repo="${2:-}"; shift 2 ;;
    --source-branch)   src_branch="${2:-}"; shift 2 ;;
    --tag)             tag_name="${2:-}"; shift 2 ;;
    --impr-branch)     impr_branch="${2:-}"; shift 2 ;;
    --trunk-subpath)   trunk_subpath="${2:-}"; shift 2 ;;
    --message)         message="${2:-}"; shift 2 ;;
    --dry-run)         dry_run=true; shift ;;
    --yes)             yes_mode=true; shift ;;
    --overwrite-trunk) overwrite_trunk=true; shift ;;
    -h|--help)         usage; exit 0 ;;
    *) _err "Unknown option: $1"; usage; exit 1 ;;
  esac
done

[[ -n "$repo" && -n "$src_branch" && -n "$tag_name" ]] || {
  _err "Missing required options."
  usage
  exit 1
}

repo="${repo%/}"

#------------------------------- helpers SVN ---------------------------------#
svn_exists() { svn info "$1" >/dev/null 2>&1; }

run() {
  if $dry_run; then
    printf '[DRY-RUN] %s\n' "$*" >&2
  else
    eval "$@"
  fi
}

detect_layout() {
  local brc="${repo}/branches" brs="${repo}/branch"
  local tgc="${repo}/tags"     tgs="${repo}/tag"

  if svn_exists "$brc"; then BR="branches"; else BR="branch"; fi
  if svn_exists "$tgc"; then TG="tags";     else TG="tag";     fi
}

#------------------------------- layout --------------------------------------#
detect_layout

branches_url="${repo}/${BR}"
tags_url="${repo}/${TG}"
trunk_url="${repo}/trunk"

if [[ -n "$trunk_subpath" ]]; then
  trunk_full="${trunk_url}/${trunk_subpath}"
else
  trunk_full="${trunk_url}"
fi

src_url="${branches_url}/${src_branch}"
tag_url="${tags_url}/${tag_name}"

if [[ -n "$impr_branch" ]]; then
  impr_url="${branches_url}/${impr_branch}"
else
  impr_url=""
fi

#------------------------------- info ----------------------------------------#
_info "Repository root:     ${repo}"
_info "Layout:              ${BR}/${TG}"
_info "Source branch:       ${src_url}"
_info "Tag:                 ${tag_url}"
_info "Trunk path:          ${trunk_full}"

if [[ -n "$impr_branch" ]]; then
  _info "Improvements branch:${impr_url}"
else
  _warn "No --impr-branch provided."
  _warn "Improvements branch will NOT be created."
  _warn "Source branch will NOT be closed."
fi

[[ -n "$message" ]] && _info "Message:             ${message}"
$overwrite_trunk && _warn "Will OVERWRITE trunk path if it exists."

#------------------------------- checks --------------------------------------#
svn_exists "$src_url" || _die "Source branch not found: ${src_url}"
svn_exists "$tag_url" && _die "Tag already exists: ${tag_url}"

if [[ -n "$impr_url" ]] && svn_exists "$impr_url"; then
  _die "Improvements branch already exists: ${impr_url}"
fi

if svn_exists "$trunk_full" && ! $overwrite_trunk; then
  _warn "Trunk path already exists: ${trunk_full}"
  _warn "Use --overwrite-trunk to overwrite it. Aborting."
  exit 2
fi

#------------------------------- confirm -------------------------------------#
if ! $yes_mode; then
  echo
  _info "Planned operations:"
  echo "  1) svn copy ${src_url} -> ${tag_url}"

  if $overwrite_trunk && svn_exists "$trunk_full"; then
    echo "  2) svn rm   ${trunk_full}"
  fi

  echo "  3) ensure trunk hierarchy"
  echo "  4) svn copy ${src_url} -> ${trunk_full}"

  if [[ -n "$impr_branch" ]]; then
    echo "  5) svn copy ${trunk_full} -> ${impr_url}"
    echo "  6) svn rm   ${src_url}"
  else
    echo "  5) (skipped) improvements branch"
    echo "  6) (skipped) source branch closure"
  fi

  read -r -p "Proceed? [y/N]: " ans
  [[ "${ans:-}" =~ ^[Yy]$ ]] || { _warn "Cancelled by user."; exit 0; }
fi

#------------------------------- 1) tag --------------------------------------#
_info "Creating tag..."
run "svn copy '${src_url}' '${tag_url}' -m $(printf %q "${message}: create tag ${tag_name}")"
_ok "Tag created."

#------------------------------- 2) trunk ------------------------------------#
if $overwrite_trunk && svn_exists "$trunk_full"; then
  _warn "Removing existing trunk path..."
  run "svn rm '${trunk_full}' -m $(printf %q "${message}: remove existing trunk path")"
fi

if [[ -n "$trunk_subpath" ]]; then
  IFS='/' read -r -a parts <<< "$trunk_subpath"
  partial="${trunk_url}"
  for p in "${parts[@]}"; do
    partial="${partial}/${p}"
    if ! svn_exists "$partial"; then
      run "svn mkdir '${partial}' -m $(printf %q "${message}: create trunk path ${partial#${repo}/}")"
    fi
  done
fi

_info "Initializing trunk from source branch..."
run "svn copy '${src_url}' '${trunk_full}' -m $(printf %q "${message}: initialize trunk from ${src_branch}")"
_ok "Trunk ready."

#------------------------------- 3) improvements -----------------------------#
if [[ -n "$impr_branch" ]]; then
  _info "Creating improvements branch..."
  run "svn copy '${trunk_full}' '${impr_url}' -m $(printf %q "${message}: create improvements branch ${impr_branch}")"
  _ok "Improvements branch created."
else
  _warn "Skipping improvements branch creation."
fi

#------------------------------- 4) close source -----------------------------#
if [[ -n "$impr_branch" ]]; then
  _info "Closing source branch..."
  run "svn rm '${src_url}' -m $(printf %q "${message}: close source branch ${src_branch}")"
  _ok "Source branch closed."
else
  _warn "Source branch preserved."
fi

_ok "Release flow completed successfully."
#EOC
#-----------------------------------------------------------------------------#

