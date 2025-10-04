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
#     3) a new improvements branch from trunk,
#     4) and then close (delete) the original branch.
#
#   Layout assumed: <REPO_ROOT>/{trunk,branches,tags}
#
# !CALLING SEQUENCE:
#   ./svn_release_flow.sh \
#     --repo https://svn.server/path/project \
#     --source-branch SMNA_v3.0.x \
#     --tag SMNA_v3.0.0 \
#     --impr-branch SMNA_v3.0.x-improvements \
#     [--trunk-subpath SMNA/SMG] [--message "Release SMNA_v3.0.0"] \
#     [--dry-run] [--yes]
#
# !OPTIONS:
#   --repo URL               Repository root (parent of trunk/branches/tags)
#   --source-branch NAME     Existing branch name under branches/
#   --tag NAME               New tag name under tags/
#   --impr-branch NAME       New branch name under branches/ for improvements
#   --trunk-subpath SUBDIR   Optional subpath inside trunk/ (e.g., "SMNA/SMG")
#   --message MSG            Commit message base
#   --dry-run                Show commands only (no changes to server)
#   --yes                    Non-interactive; skip confirmation prompt
#
# !BEHAVIOR & SAFETY:
#   - Aborts if tag or improvement branch already exist (unless you change logic).
#   - Creates trunk if missing; aborts if trunk exists (to evitar overwrite).
#   - "Close branch" = 'svn rm' on the source branch URL with a close message.
#   - Requires 'svn' CLI configured (auth via cache/agent).
#
# !EXAMPLES:
#   ./svn_release_flow.sh --repo https://svn.cptec.inpe.br/smna \
#     --source-branch SMNA_v3.0.x \
#     --tag SMNA_v3.0.0 \
#     --impr-branch SMNA_v3.0.x-improvements \
#     --trunk-subpath SMNA/SMG \
#     --message "Release SMNA_v3.0.0 (tag), create trunk and improvements branch" \
#     --dry-run
#
# !REQUIRES:
#   Bash 4+, svn in PATH
#
# !EXIT CODES:
#   0 OK, 1 usage, 2 config/missing, 3 svn error
#EOP
#-----------------------------------------------------------------------------#
set -euo pipefail

#------------------------------- logging -------------------------------------#
_log()      { printf '%s\n' "$*" >&2; }
_info()     { _log "[INFO] $*"; }
_ok()       { _log "[OK]   $*"; }
_warn()     { _log "[WARN] $*"; }
_err()      { _log "[ERROR] $*"; }
_die()      { _err "$*"; exit 3; }

#------------------------------- usage ---------------------------------------#
usage() {
  sed -n '1,/^#EOP/p' "$0" | sed -n 's/^# \{0,1\}//p'
  cat <<'EOF'

Usage:
  $0 --repo URL --source-branch NAME --tag NAME --impr-branch NAME
     [--trunk-subpath SUBDIR] [--message MSG] [--dry-run] [--yes] [--overwrite-trunk]

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
    --repo)           repo="${2:-}"; shift 2 ;;
    --source-branch)  src_branch="${2:-}"; shift 2 ;;
    --tag)            tag_name="${2:-}"; shift 2 ;;
    --impr-branch)    impr_branch="${2:-}"; shift 2 ;;
    --trunk-subpath)  trunk_subpath="${2:-}"; shift 2 ;;
    --message)        message="${2:-}"; shift 2 ;;
    --dry-run)        dry_run=true; shift ;;
    --yes)            yes_mode=true; shift ;;
    --overwrite-trunk) overwrite_trunk=true; shift ;;
    -h|--help)        usage; exit 0 ;;
    *) _err "Unknown option: $1"; usage; exit 1 ;;
  esac
done

[[ -n "$repo" && -n "$src_branch" && -n "$tag_name" && -n "$impr_branch" ]] || { _err "Missing required options."; usage; exit 1; }
repo="${repo%/}"

# Normaliza repo sem barra final
repo="${repo%/}"

#------------------------------- helpers SVN ---------------------------------#
svn_exists() { svn info "$1" >/dev/null 2>&1; }           # URL exists?
svn_ls()     { svn ls "$1" >/dev/null 2>&1; }

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

#------------------------------- sanity checks -------------------------------#
detect_layout

branches_url="${repo}/${BR}"
tags_url="${repo}/${TG}"
trunk_url="${repo}/trunk"
[[ -n "$trunk_subpath" ]] && trunk_full="${trunk_url}/${trunk_subpath}" || trunk_full="${trunk_url}"

src_url="${branches_url}/${src_branch}"
tag_url="${tags_url}/${tag_name}"
impr_url="${branches_url}/${impr_branch}"

_info "Repository root:     ${repo}"
_info "Layout:              ${BR}/${TG}"
_info "Source branch (URL): ${src_url}"
_info "Tag (URL):           ${tag_url}"
_info "Trunk (URL):         ${trunk_full}"
_info "Improvements branch: ${impr_url}"
[[ -n "$message" ]] && _info "Message:             ${message}"
$overwrite_trunk && _warn "Will OVERWRITE trunk path: ${trunk_full}"

svn_exists "$src_url" || _die "Source branch not found: ${src_url}"
if svn_exists "$tag_url";  then _die "Tag already exists: ${tag_url}"; fi
if svn_exists "$impr_url"; then _die "Improvements branch already exists: ${impr_url}"; fi

# Trunk existence rule
if svn_exists "$trunk_full"; then
  if ! $overwrite_trunk; then
    _warn "Trunk path already exists: ${trunk_full}"
    _warn "Use --overwrite-trunk to delete and recreate it safely. Aborting."
    exit 2
  fi
fi

# Confirm
if ! $yes_mode; then
  echo
  _info "About to perform the following operations:"
  echo "  1) svn copy ${src_url} -> ${tag_url}"
  if $overwrite_trunk && svn_exists "$trunk_full"; then
    echo "  2) svn rm   ${trunk_full} (overwrite trunk path)"
  fi
  echo "  3) ensure trunk path hierarchy exists"
  echo "  4) svn copy ${src_url} -> ${trunk_full}"
  echo "  5) svn copy ${trunk_full} -> ${impr_url}"
  echo "  6) svn rm   ${src_url} (close source branch)"
  read -r -p "Proceed? [y/N]: " ans
  [[ "${ans:-}" =~ ^[Yy]$ ]] || { _warn "Cancelled by user."; exit 0; }
fi

#------------------------------- 1) tag --------------------------------------#
_info "Creating tag from source branch..."
run "svn copy '${src_url}' '${tag_url}' -m $(printf %q "${message}: create tag ${tag_name} from ${src_branch}")"
_ok "Tag created: ${tag_url}"

#------------------------------- 2) trunk ------------------------------------#
# overwrite trunk path if requested
if $overwrite_trunk && svn_exists "$trunk_full"; then
  run "svn rm '${trunk_full}' -m $(printf %q "${message}: remove existing trunk subpath to overwrite")"
fi

_info "Creating trunk (subpath) from source branch..."
# Garante que diretórios intermediários no trunk existam (quando trunk_subpath tem múltiplos níveis)
if [[ -n "$trunk_subpath" ]]; then
  # Cria skeleton do trunk se necessário (svn requires explicit dir creation)
  # Ex.: trunk/SMNA e trunk/SMNA/SMG
  IFS='/' read -r -a parts <<< "$trunk_subpath"
  partial="${trunk_url}"
  for p in "${parts[@]}"; do
    partial="${partial}/${p}"
    if ! svn_exists "$partial"; then
      run "svn mkdir '${partial}' -m $(printf %q "${message}: create trunk path ${partial#${repo}/}")"
    fi
  done
fi

# Copia conteúdo do branch para o trunk_full (último nível)
run "svn copy '${src_url}' '${trunk_full}' -m $(printf %q "${message}: initialize trunk from ${src_branch}")"
_ok "Trunk ready at: ${trunk_full}"

#------------------------------- 3) improvements branch ----------------------#
_info "Creating improvements branch from trunk..."
run "svn copy '${trunk_full}' '${impr_url}' -m $(printf %q "${message}: create improvements branch ${impr_branch} from trunk")"
_ok "Improvements branch created: ${impr_url}"

#------------------------------- 4) close original branch --------------------#
_info "Closing original source branch..."
run "svn rm '${src_url}' -m $(printf %q "${message}: close source branch ${src_branch} after creating tag/trunk/impr-branch")"
_ok "Source branch closed: ${src_url}"

_ok "All done."
#EOC
#-----------------------------------------------------------------------------#

