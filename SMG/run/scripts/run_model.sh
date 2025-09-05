#! /bin/bash -x

# Descomentar para debugar
#set -o xtrace

### Define hpc_name below is needed because in some cases (ex.: First Run),
### this script may be called from the command line.
### 
lognode=`cat /proc/sys/kernel/hostname | cut  -b 1-6`

case $lognode in

  clogin)
    STR=`uname -a`
    SUB='cray'
    if [[ "$STR" == *"$SUB"* ]]; then
      echo -n "This will run on cray XC50 ..."
      export hpc_name="XC50"
    fi
    ;;

  headno)
    STR=`uname -a`
    SUB='egeon'
    if [[ "$STR" == *"$SUB"* ]]; then
      echo -n "This will run on EGEON Cluster ..."
      export hpc_name="egeon"
    fi
    ;;

  egeon-)
    STR=`uname -a`
    SUB='egeon'
    if [[ "$STR" == *"$SUB"* ]]; then
      echo -n "This will run on EGEON Cluster ..."
      export hpc_name="egeon"
    fi
    ;;

  *)
    mach=`cat /proc/sys/kernel/hostname`
    echo -n "The configurations for "$mach" is not defined yet !"
    echo -n "1) Add the machine to the defined systems in etc/mach ; and"
    echo -n "2) add an option for it in the function copy_fixed_files in etc/functions"
    exit
    ;;
esac

# Carregando as variaveis do sistema
SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"
RootDir="$(dirname "$SCRIPT_PATH")"
export SMG_ROOT=${RootDir}
# Localiza e carrega o runPre.func (pode estar no home_run_bam ou ao lado deste script)
if [ -f "${home_run_bam}/runPre.func" ]; then
  RUNPRE_FUNC="${home_run_bam}/runPre.func"
elif [ -f "./runPre.func" ]; then
  RUNPRE_FUNC="./runPre.func"
else
  echo "[ERROR] runPre.func not found (looked in \${home_run_bam} and ./)"; exit 1
fi
# shellcheck source=/dev/null
source "${RUNPRE_FUNC}"

# === NOVO: Carregar EnvironmentalVariables no MESMO diretório do runPre.func ===
# Resolve o diretório do runPre.func para achar o EnvironmentalVariables “irmão”
RUNPRE_DIR="$(cd -- "$(dirname -- "${RUNPRE_FUNC}")" && pwd)"
ENV_FILE="${RUNPRE_DIR}/EnvironmentalVariables"
if [ -f "${ENV_FILE}" ]; then
  # shellcheck source=/dev/null
  source "${ENV_FILE}"
else
  echo "[ERROR] EnvironmentalVariables not found at ${ENV_FILE}"; exit 1
fi

# Sanidade: garantir que funções que usaremos existam
for fn in getBAMSize getMPIinfo; do
  if ! type -t "$fn" >/dev/null 2>&1; then
    echo "[ERROR] function '$fn' not found after sourcing EnvironmentalVariables"; exit 1
  fi
done


# carregando funcoes do pre-processamento
source ${RUNPRE_FUNC}

# Verificando argumentos de entrada
if [ -z "${1}" ]
then
  echo "LABELANL is not set" 
  exit 3
else
  export LABELANL=${1}
fi
if [ -z "${2}" ]
then
  echo "LABELFCT is not set" 
  exit 3
else
  export LABELFCT=${2} 
fi
if [ -z "${3}" ]
then
  echo "PREFIX is not set" 
  exit 3
else
  export PREFIX=${3}
fi
if [ -z "${4}" ]
then
  echo "TRC is not set" 
  exit 3
else
  export TRC=${4}
fi
if [ -z "${5}" ]
then
  echo "NLV is not set" 
  exit 3
else
  export NLV=${5}
fi
if [ -z "${6}" ]
then
  echo "NPROC is not set"
  case ${hpc_name} in
     	egeon)  echo "setting to 64"
             	export NPROC=64  # ntasks
	;;
  	XC50)  	echo "setting to 480"
  		export NPROC=480             	
	;;
  esac
else
  export NPROC=${6}
fi
if [ "$#" == 7 ]
then 
  export RUNPOS=$(echo ${7} | tr [:upper:] [:lower:])   
else 
  export RUNPOS="yes"
fi

## Resolve MPI/OpenMP layout from current HPC and user overrides
# Defaults (mantém o comportamento atual se não vier do caller)
: "${tasks_per_node:=16}"
: "${cpus_per_task:=8}"

# Faz a resolução completa (detecta máquina, checa limites, etc.)
# Aceita também -np/-N/-d do próprio script, se você repassar.
getMPIinfo -np "${NPROC}" -N "${tasks_per_node}" -d "${cpus_per_task}" || {
  echo "[ERROR] getMPIinfo failed"; exit 1;
}
# Opcional: sincronizar com nomes locais caso você use essas variáveis adiante
tasks_per_node="${TasksPerNode}"
cpus_per_task="${ThreadsPerMPITask}"


getBAMSize ${TRC}
export postfix=$(printf "G%5.5dL%3.3d \n" $JM $NLV)
export MRES=`echo ${TRC} ${NLV} | awk '{printf("TQ%4.4dL%3.3d\n",$1,$2)}'`
export ANL="GANL${PREFIX}${LABELANL}S.unf.${MRES}"
### export LABELFGS=`${inctime} ${LABELANL} +6h %y4%m2%d2%h2`
export LABELFGS=`date -u +%Y%m%d%H -d "${LABELANL:0:8} ${LABELANL:8:2} +6 hours" `
### export YESTERDAY=$(${inctime} ${LABELANL} -1d %y4%m2%d200)
export YESTERDAY=`date -u +%Y%m%d00 -d "${LABELANL:0:8} ${LABELANL:8:2} -1 days" `
export modelDataIn=${subt_model_bam}/datain
export gsiDataOut=${subt_gsi_dataout}/${LABELANL}

echo -e ""
echo -e "\033[34;1m >> Submetendo o MCGA:\033[m \033[31;1m${LABELANL}\033[m"
echo -e " "
echo -e "\033[34;1m > Resolucao (espectral) : \033[m \033[31;1m${MRES}\033[m"
echo -e "\033[34;1m > Resolucao (grade)     : \033[m \033[31;1m${postfix}\033[m"
echo -e "\033[34;1m > Condicao Inicial      : \033[m \033[31;1m${LABELANL}\033[m"
echo -e "\033[34;1m > Previsao ate          : \033[m \033[31;1m${LABELFCT}\033[m"
echo -e "\033[34;1m > Pos-proc. Ativado     : \033[m \033[31;1m${RUNPOS}\033[m"

#
# mudando para diretorio dos scripts do BAM
#

cd ${home_run_bam}

#
# rodando o somente o Chopping do pré para pegar o arquivo de ozônio
# saida gerada no bam/model/datain/ 

/bin/bash runPre -v -t ${TRC} -l ${NLV} -I ${LABELANL} -s -n chp -O

STATUS=$?
echo "1st call to runPre. Status: "${STATUS}
if [ ${STATUS} -ne 0 ];then
   exit ${STATUS}
fi

cp -f ${modelDataIn}/OZONSMT${LABELANL}S.grd.${postfix} ${modelDataIn}/OZON${PREFIX}${LABELANL}S.grd.${postfix}

#
# remove arquivos desnecessarios
#

rm -fr ${modelDataIn}/GANLSMT${LABELANL}S.unf.*
rm -fr ${modelDataIn}/OZONSMT${LABELANL}S.unf.*

#
# Copiando Analise do GSI para o Model/DataIn
#

cp -pfr ${gsiDataOut}/GANL${PREFIX}${LABELANL}S.unf.${MRES} ${modelDataIn}

#
# Rodando os demais processos do pré e usando a análise do GSI
#
# /bin/bash runPre -v -t 299 -l 64 -I ${LABELANL}  -n 0 -p SMT -s -O -T -G -Gp gblav -Gt Grid

/bin/bash runPre -v -t ${TRC} -l ${NLV} -I ${LABELANL} -p ${PREFIX} -n das
STATUS=$?
echo "2nd call to runPre. Status: "${STATUS}
if [ ${STATUS} -ne 0 ];then
   exit ${STATUS}
fi

# Rodando o Modelo

/bin/bash runModel -das -v -np ${NPROC} -N ${tasks_per_node} -d ${cpus_per_task} \
                   -t ${TRC} -l ${NLV} -I ${LABELANL} -F ${LABELFCT} -W ${LABELFCT} \
                   -p ${PREFIX} -s sstwkl -ts 3 -r -tr 6 -i 2


# Pos-processa as previsoes caso a variavel RUNPOS possua o valor Yes ou Y
if [ ${RUNPOS} == "yes" -o ${RUNPOS} == "y" ]
then

  # Verifica se o executavel se encontra em ${subt_pos_bam_run}
  if [ ! -e ${subt_pos_bam}/exec/PostGrib ]
  then

    cp -v ${home_pos_bam}/exec/PostGrib ${subt_pos_bam}/exec

  fi

  cd ${home_run_bam}
  echo   "./runPos -np 15 -N 3 -d 8 -t ${TRC} -l ${NLV} -I ${LABELANL} -F ${LABELFCT} -p ${PREFIX} > /dev/null 2>&1"
  /bin/bash runPos -np 15 -N 3 -d 8 -t ${TRC} -l ${NLV} -I ${LABELANL} -F ${LABELFCT} -p ${PREFIX} > /dev/null 2>&1
  STATUS=$?
   if [ ${STATUS} -ne 0 ];then
      exit ${STATUS}
   fi

fi

exit 0
