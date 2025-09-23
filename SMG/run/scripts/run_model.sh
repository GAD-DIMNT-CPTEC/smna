#!/usr/bin/env bash
#BOP
# !SCRIPT: run_bam_cycle.sh
# !DESCRIPTION:
#   Orquestra o pré-processamento (chopping de ozônio + pré completo), execução do
#   modelo BAM e pós-processamento opcional, consumindo a análise do GSI.
#   Mantém a convenção de argumentos posicionais:
#     1) LABELANL (YYYYMMDDHH)
#     2) LABELFCT (YYYYMMDDHH)
#     3) PREFIX   (ex.: SMT/GA... conforme teu fluxo)
#     4) TRC      (truncamento espectral, ex.: 299)
#     5) NLV      (n níveis, ex.: 64)
#     6) NPROC    (ntasks; default: egeon=64, XC50=480)
#     7) RUNPOS   (yes|y para ativar pós; default yes)
#
# !USAGE:
#   ./run_bam_cycle.sh 2025010100 2025010200 SMT 299 64 64 yes
#
# !NOTES:
#   - Este script resolve o diretório do runPre.func e carrega o arquivo
#     EnvironmentalVariables “irmão” do runPre.func.
#   - Requer que getBAMSize e getMPIinfo estejam disponíveis após o source.
#   - Usa GNU date para aritmética de datas.
#EOP


# -------------------- Localização do root do projeto ----------------- #
resolve_root() {
__SRC="${BASH_SOURCE[0]}"
while [[ -L "$__SRC" ]]; do
  __LINK="$(readlink -- "$__SRC")"
  if [[ "$__LINK" = /* ]]; then
    __SRC="$__LINK"
  else
    __SRC="$(cd -- "$(dirname -- "$__SRC")" && cd -- "$(dirname -- "$__LINK")" && pwd)/$(basename -- "$__LINK")"
  fi
done
export RUN_MODEL_DIR="$(cd -- "$(dirname -- "$__SRC")" && pwd -P)"

# --- load init (which loads helpers) ---
__init_path="${RUN_MODEL_DIR}/__init__.sh"
# shellcheck disable=SC1090
. "$__init_path" || :        # no extra prints; don't break under set -e

}

# --------------- Carregar runPre.func e EnvironmentalVariables ------- #
load_runpre_and_env() {
  # runPre.func pode estar no ${home_run_bam} (definido por env) ou ao lado
  if [[ -n "${home_run_bam:-}" && -f "${home_run_bam}/runPre.func" ]]; then
    RUNPRE_FUNC="${home_run_bam}/runPre.func"
  elif [[ -f "./runPre.func" ]]; then
    RUNPRE_FUNC="./runPre.func"
  else
    die "runPre.func not found (searched in \${home_run_bam} and ./)"
  fi
  # shellcheck source=/dev/null
  source "${RUNPRE_FUNC}"

  # EnvironmentalVariables “irmão” do runPre.func
  RUNPRE_DIR="$(cd -- "$(dirname -- "${RUNPRE_FUNC}")" && pwd -P)"
  ENV_FILE="${RUNPRE_DIR}/EnvironmentalVariables"
  [[ -f "${ENV_FILE}" ]] || die "EnvironmentalVariables not found at ${ENV_FILE}"
  # shellcheck source=/dev/null
  source "${ENV_FILE}"

  # Sanidade: funções requeridas
  for fn in getBAMSize getMPIinfo; do
    type -t "${fn}" >/dev/null 2>&1 || die "Required function '${fn}' not found after sourcing env."
  done
}

# ---------------------- Parse dos argumentos ------------------------ #
parse_args() {
  if [[ $# -lt 5 ]]; then
    cat >&2 <<'USAGE'
[ERROR] Missing required arguments.

USAGE:
  run_bam_cycle.sh LABELANL LABELFCT PREFIX TRC NLV [NPROC] [RUNPOS]

ARGS:
  LABELANL  Analysis datetime (YYYYMMDDHH)
  LABELFCT  Forecast end datetime (YYYYMMDDHH)
  PREFIX    Background prefix (e.g., SMT)
  TRC       Spectral truncation (e.g., 299)
  NLV       Number of levels (e.g., 64)
  NPROC     (optional) ntasks; default depends on HPC (egeon=64, XC50=480)
  RUNPOS    (optional) yes|y to enable post; default yes
USAGE
    exit 3
  fi

  LABELANL="$1"; export LABELANL
  LABELFCT="$2"; export LABELFCT
  PREFIX="$3";   export PREFIX
  TRC="$4";      export TRC
  NLV="$5";      export NLV

  # NPROC: default por plataforma
  if [[ $# -ge 6 && -n "${6:-}" ]]; then
    NPROC="$6"
  else
    case "${hpc_name}" in
      egeon) NPROC=64  ;;
      XC50)  NPROC=480 ;;
      *)     NPROC=64  ;; # fallback
    esac
    _warn "NPROC not provided; defaulting to ${NPROC} for ${hpc_name}"
  fi
  export NPROC

  # RUNPOS normalizado (default yes)
  if [[ $# -ge 7 && -n "${7:-}" ]]; then
    RUNPOS="$(printf '%s' "${7}" | tr '[:upper:]' '[:lower:]')"
  else
    RUNPOS="yes"
  fi
  export RUNPOS
}

# --------------------- Resolver layout MPI/OMP ---------------------- #
resolve_mpi_layout() {
  # Defaults locais (respeita overrides do ambiente se existirem)
  : "${tasks_per_node:=16}"
  : "${cpus_per_task:=8}"

  # Pede ao helper calcular e validar
  getMPIinfo -np "${NPROC}" -N "${tasks_per_node}" -d "${cpus_per_task}" || \
    die "getMPIinfo failed"

  # Aderir aos nomes resolvidos pelo helper
  tasks_per_node="${TasksPerNode}"
  cpus_per_task="${ThreadsPerMPITask}"
}

# -------------------------- Execução principal ---------------------- #
main() {
  resolve_root
  load_runpre_and_env
  parse_args "$@"
  resolve_mpi_layout

  # getBAMSize define JM (usado para postfix)
  getBAMSize "${TRC}"

  postfix="$(printf 'G%5.5dL%3.3d' "${JM}" "${NLV}")"
  export postfix

  MRES="$(printf 'TQ%4.4dL%3.3d' "${TRC}" "${NLV}")"
  export MRES

  ANL="GANL${PREFIX}${LABELANL}S.unf.${MRES}"
  export ANL

  # LABELFGS = LABELANL + 6h; YESTERDAY = YYYYMMDD00 do dia anterior
  LABELFGS="$(date -u +%Y%m%d%H -d "${LABELANL:0:8} ${LABELANL:8:2} +6 hours")"
  YESTERDAY="$(date -u +%Y%m%d00 -d "${LABELANL:0:8} ${LABELANL:8:2} -1 days")"
  export LABELFGS YESTERDAY

  # Pastas de I/O (vindas do ENV_FILE)
  : "${subt_model_bam:?subt_model_bam is required}"
  : "${subt_gsi_dataout:?subt_gsi_dataout is required}"

  modelDataIn="${subt_model_bam}/datain"
  gsiDataOut="${subt_gsi_dataout}/${LABELANL}"
  export modelDataIn gsiDataOut

  printf '\n\033[34;1m >> Submetendo o MCGA:\033[m \033[31;1m%s\033[m\n\n' "${LABELANL}"
  printf '\033[34;1m > Resolucao (espectral) :\033[m \033[31;1m%s\033[m\n' "${MRES}"
  printf '\033[34;1m > Resolucao (grade)     :\033[m \033[31;1m%s\033[m\n' "${postfix}"
  printf '\033[34;1m > Condicao Inicial      :\033[m \033[31;1m%s\033[m\n' "${LABELANL}"
  printf '\033[34;1m > Previsao ate          :\033[m \033[31;1m%s\033[m\n' "${LABELFCT}"
  printf '\033[34;1m > Pos-proc. Ativado     :\033[m \033[31;1m%s\033[m\n' "${RUNPOS}"

  # Vai para o diretório de execução do BAM
  cd -- "${home_run_bam:?home_run_bam is required}"

  # 1) Pré (apenas chopping de ozônio); saída em bam/model/datain/
  /bin/bash runPre -v -t "${TRC}" -l "${NLV}" -I "${LABELANL}" -s -n chp -O
  _log_info "1st call to runPre completed."

  # Copia OZONSMT -> OZON${PREFIX}
  cp -f -- "${modelDataIn}/OZONSMT${LABELANL}S.grd.${postfix}" \
            "${modelDataIn}/OZON${PREFIX}${LABELANL}S.grd.${postfix}"

  # Limpa arquivos desnecessários
  rm -f -- "${modelDataIn}/GANLSMT${LABELANL}S.unf."* || true
  rm -f -- "${modelDataIn}/OZONSMT${LABELANL}S.unf."* || true

  # 2) Copia análise do GSI para o Model/DataIn
  cp -fp -- "${gsiDataOut}/GANL${PREFIX}${LABELANL}S.unf.${MRES}" "${modelDataIn}/"
  _log_info "GSI analysis staged into modelDataIn."

  # 3) Pré completo (usando a análise do GSI)
  /bin/bash runPre -v -t "${TRC}" -l "${NLV}" -I "${LABELANL}" -p "${PREFIX}" -n das
  _log_info "2nd call to runPre (full pre) completed."

  # 4) Rodar o Modelo
  /bin/bash runModel -das -v \
    -np "${NPROC}" -N "${tasks_per_node}" -d "${cpus_per_task}" \
    -t "${TRC}" -l "${NLV}" -I "${LABELANL}" -F "${LABELFCT}" -W "${LABELFCT}" \
    -p "${PREFIX}" -s sstwkl -ts 3 -r -tr 6 -i 2

  # 5) Pós-processamento opcional
  case "${RUNPOS}" in
    yes|y)
      # Garante executável de pós em ${subt_pos_bam}/exec/PostGrib
      if [[ ! -x "${subt_pos_bam}/exec/PostGrib" ]]; then
        cp -v -- "${home_pos_bam}/exec/PostGrib" "${subt_pos_bam}/exec/"
      fi

      cd -- "${home_run_bam}"
      _log_info "Starting runPos..."
      /bin/bash runPos -np 15 -N 3 -d 8 \
        -t "${TRC}" -l "${NLV}" -I "${LABELANL}" -F "${LABELFCT}" -p "${PREFIX}" \
        > /dev/null 2>&1
      _log_info "runPos completed."
      ;;
    *)
      _log_info "Post-processing disabled by RUNPOS='${RUNPOS}'."
      ;;
  esac

  _log "Workflow finished successfully."
}

main "$@"

