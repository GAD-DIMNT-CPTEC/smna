#! /bin/bash

# Descomentar para debugar
#set -o xtrace


WhereIam=$(dirname ${BASH_SOURCE})

# Carregando as variaveis do sistema
source ${WhereIam}/../../config_smg.ksh vars_export


# carregando funcoes do pre-processamento

source ${home_run_bam}/runPre.func




args=${@} #save arguments temporarily
while (( $# )); do
   opt=$1
   case ${opt} in
       -t) TRC=$2; shift 2;; # shift twice [p.e. -t 299 ]
       -l) NLV=$2; shift 2;;
       -p) PREFIX=$2; shift 2;;
       -I) LABELANL=$2; shift 2;;
       -F) LABELFCT=$2; shift 2;;
       -np) MPITasks=$2; shift 2 ;;
       -N) TasksPerNode=$2; shift 2;;
       -d) ThreadsPerMPITask=$2; shift 2;;
       -pw) walltime=$2;shift 2;;
       -pq) queue=$2; shift 2;;
       -pn) queue_name=$2; shift 2;;
       -pos) RUNPOS='.TRUE.';shift 1;;
       -h) cat < ${0} | sed -n '/^#BOP/,/^#EOP/p' ; exit 0;;
       *) echo -e "\033[31;1mWarning:\033[m Unknown argument:\033[33;1m $opt\033[m"; shift 1;;
   esac
done
set -- $args #restore arguments


# Verificando argumentos de entrada
# Data da Condição Inicial (cold start)
if [ -z "${LABELANL}" ];then
  echo "LABELANL is not set" 
  exit 1
fi

# Data final das previsões
if [ -z "${LABELFCT}" ];then
  echo "LABELFCT is not set" 
  exit 2
fi

# prefixo dos arquivos
if [ -z "${PREFIX}" ];then
  export PREFIX=CPT
fi

# truncamento
if [ -z "${TRC}" ];then
  echo "TRC is not set" 
  exit 3
fi

# numero de niveis verticais
if [ -z "${NLV}" ];then
  echo "NLV is not set" 
  exit 4
fi

# Numero de processadores que serao utilizados no Job
if [ -z ${MPITasks} ];then
   MPITasks=80
fi

# Numero de processadores utilizados por tarefas MPI
if [ -z ${TasksPerNode} ];then
   TasksPerNode=40
fi

# Number of cores hosting OpenMP threads
if [ -z ${ThreadsPerMPITask} ]; then
   ThreadsPerMPITask=1
fi

# Number of cores de cada nó do sistam
# Tupa = 24
# XC50 = 40
if [ -z ${CoresPerNode} ]; then
   CoresPerNode=40
fi


# define PBS walltime 
if [ -z ${walltime} ];then
   walltime=00:45:00
fi

# define PBS queue
if [ -z ${queue} ];then
   queue=pesq
fi

# define PBS queue
if [ -z ${queue_name} ];then
   queue_name="BAM${TRC}"
fi


if [ -z ${RUNPOS} ];then
   RUNPOS='.FALSE.'
fi



getBAMSize ${TRC}
export postfix=$(printf "G%5.5dL%3.3d \n" $JM $NLV)
export MRES=`echo ${TRC} ${NLV} | awk '{printf("TQ%4.4dL%3.3d\n",$1,$2)}'`
export ANL="GANL${PREFIX}${LABELANL}S.unf.${MRES}"
export LABELFGS=`${inctime} ${LABELANL} +6h %y4%m2%d2%h2`
export YESTERDAY=$(${inctime} ${LABELANL} -1d %y4%m2%d200)
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


# Rodando o Modelo

/bin/bash runModel -das -f -v -pw ${walltime} -pn ${queue_name} -np ${MPITasks} -N ${TasksPerNode} -d ${ThreadsPerMPITask} -t ${TRC} -l ${NLV} -I ${LABELANL} -F ${LABELFCT} -W  ${LABELFCT} -p ${PREFIX} -u -s sstwkl -ts 6 -i 2
STATUS=$?
if [ ${STATUS} -ne 0 ];then
   exit ${STATUS}
fi

# Pos-processa as previsoes caso a variavel RUNPOS possua o valor Yes ou Y
if [ ${RUNPOS} == ".TRUE." ];then

  # Verifica se o executavel se encontra em ${subt_pos_bam_run}
  if [ ! -e ${subt_pos_bam}/exec/PostGrib ]
  then 
  
    cp -v ${home_pos_bam}/source/PostGrib ${subt_pos_bam}/exec
  
  fi

  cd ${home_run_bam}
  echo   "./runPos -np 120 -N 12 -d 1 -t ${TRC} -l ${NLV} -I ${LABELANL} -F ${LABELFCT} -p ${PREFIX} > /dev/null 2>&1"
  /bin/bash runPos -np 120 -N 12 -d 1 -t ${TRC} -l ${NLV} -I ${LABELANL} -F ${LABELFCT} -p ${PREFIX} > /dev/null 2>&1
  STATUS=$?
   if [ ${STATUS} -ne 0 ];then
      exit ${STATUS}
   fi

fi

exit 0
