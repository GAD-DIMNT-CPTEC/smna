#!/usr/bin/env bash
#------------------------------------------------------------------------------#
#BOP
# !SCRIPT: clean_dist.sh
# !DESCRIPTION:
#   Safely clean the SMG distribution by removing build artifacts, temporary
#   files, and (optionally) binary data files (IEEE, GRIB, NetCDF, HDF, etc.)
#   without touching source code or text-based scripts. Optionally removes such
#   files from SVN.
#
# !USAGE:
#   ./utils/clean_dist.sh [--dry-run] [--svn] [--yes] [--bin]
#
# !OPTIONS:
#   --dry-run   Show what would be removed without deleting anything
#   --svn       Also remove matching files from SVN repository
#   --yes       Non-interactive mode (no confirmation)
#   --bin       Include scientific and binary files (IEEE, GRIB, NetCDF, HDF, etc.)
#
# !AUTHOR:  J. G. de Mattos (GDAD/CPTEC/INPE)
# !REVISION: 05 Oct 2025
#EOP
#------------------------------------------------------------------------------#

set -uo pipefail
C_BLUE=$'\033[1;34m'; C_GREEN=$'\033[1;32m'; C_WARN=$'\033[1;33m'; C_RST=$'\033[0m'

echo "${C_BLUE}#------------------------------------------------------------------------------#${C_RST}"
echo "${C_BLUE}# Cleaning SMG distribution...${C_RST}"
echo "${C_BLUE}#------------------------------------------------------------------------------#${C_RST}"

DRY_RUN=false
SVN_MODE=false
AUTO_YES=false
CLEAN_BINARIES=false

for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY_RUN=true ;;
    --svn) SVN_MODE=true ;;
    --yes) AUTO_YES=true ;;
    --bin) CLEAN_BINARIES=true ;;
    *)
      echo "${C_WARN}[WARN] Ignoring unknown argument: $arg${C_RST}"
      ;;
  esac
done

# Patterns to clean (build + temporários)
patterns=(
  "*.o" "*.mod" "*.a" "*.so" "*.out" "*.err" "*.tmp" "*.bak"
  "*~" ".DS_Store" "nohup.out" "*.swp" "*.log" "__pycache__"
  "core.*" "ecbuild.log" "CMakeCache.txt" "cmake_install.cmake"
)

removed_local=0
removed_svn=0
removed_bin=0

#------------------------------------------------------------------------------#
# STEP 1: Remove arquivos conhecidos pelos padrões
#------------------------------------------------------------------------------#
for pat in "${patterns[@]}"; do
  echo "[INFO] Searching for: $pat"
  mapfile -d '' files < <(find . -name "$pat" -print0 2>/dev/null)
  (( ${#files[@]} == 0 )) && continue

  svn_files=()
  local_files=()

  for f in "${files[@]}"; do
    if svn info "$f" &>/dev/null; then
      svn_files+=("$f")
    else
      local_files+=("$f")
    fi
  done

  # --- Local files ---
  if ((${#local_files[@]} > 0)); then
    if $DRY_RUN; then
      printf "  [DRY-RUN] Would remove %d local files:\n" "${#local_files[@]}"
      printf "    %s\n" "${local_files[@]}"
    else
      printf "  [INFO] Removing %d local files:\n" "${#local_files[@]}"
      printf "    %s\n" "${local_files[@]}"
      rm -rf "${local_files[@]}" 2>/dev/null || true
      ((removed_local+=${#local_files[@]}))
    fi
  fi

  # --- SVN files ---
  if $SVN_MODE && ((${#svn_files[@]} > 0)); then
    printf "  [INFO] Found %d versioned files matching '%s'\n" "${#svn_files[@]}" "$pat"
    if $DRY_RUN; then
      printf "    [DRY-RUN] Would svn remove: %s\n" "${svn_files[@]}"
    else
      if $AUTO_YES; then
        svn rm --force "${svn_files[@]}" || true
        ((removed_svn+=${#svn_files[@]}))
      else
        read -rp "  → Remove these from SVN? [y/N]: " ans
        if [[ "${ans,,}" == y ]]; then
          svn rm --force "${svn_files[@]}" || true
          ((removed_svn+=${#svn_files[@]}))
        fi
      fi
    fi
  fi
done

#------------------------------------------------------------------------------#
# STEP 2: Detecta binários científicos (IEEE, GRIB, NetCDF, HDF, etc.)
#------------------------------------------------------------------------------#
if $CLEAN_BINARIES; then
  echo
  echo "[INFO] Scanning for binary data files (IEEE, GRIB, NetCDF, HDF, etc.)..."

  mapfile -d '' bin_candidates < <(find . -type f -size +1k ! -path "./.svn/*" -print0 2>/dev/null)

  for f in "${bin_candidates[@]}"; do
    [[ ! -s "$f" ]] && continue
    desc=$(file "$f" 2>/dev/null || true)

    # ignora textos e scripts
    if grep -qiE "text|script|source|ASCII|UTF-8" <<<"$desc"; then
      continue
    fi

    # detecta binários científicos
    if grep -qiE "IEEE|GRIB|NetCDF|HDF|Hierarchical|CDF|binary data|compiled object" <<<"$desc"; then
      if $DRY_RUN; then
        echo "  [DRY-RUN] Binary candidate: $f"
        echo "             → $desc"
      else
        if svn info "$f" &>/dev/null; then
          if $SVN_MODE; then
            if $AUTO_YES; then
              svn rm --force "$f" && ((removed_bin++))
            else
              read -rp "  → Remove binary file from SVN ($f)? [y/N]: " ans
              [[ "${ans,,}" == y ]] && svn rm --force "$f" && ((removed_bin++))
            fi
          fi
        else
          echo "  [INFO] Removing binary data file: $f"
          echo "         → $desc"
          rm -f "$f" && ((removed_bin++))
        fi
      fi
    fi
  done
fi

#------------------------------------------------------------------------------#
# STEP 3: Summary
#------------------------------------------------------------------------------#
echo
echo "${C_GREEN}[OK]${C_RST} Cleanup completed."
echo "   Local files removed: $removed_local"
echo "   SVN files removed:   $removed_svn"
if $CLEAN_BINARIES; then
  echo "   Binary data removed:  $removed_bin"
fi

if $DRY_RUN; then
  echo "${C_WARN}[NOTE]${C_RST} Dry-run mode — no files actually deleted."
fi

