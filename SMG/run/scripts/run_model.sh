#! /bin/bash

# Descomentar para debugar
#set -o xtrace


# Carregando as variaveis do sistema
source /home/jose.aravequia/SMNA_v3.0.0.t11889/SMG/config_smg.ksh vars_export

# carregando funcoes do pre-processamento

source ${home_run_bam}/runPre.func

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
  echo "setting to 480"
  export NPROC=480
else
  export NPROC=${6}
fi
if [ "$#" == 7 ]
then
  export RUNPOS=$(echo ${7} | tr [:upper:] [:lower:])
else
  export RUNPOS="yes"
fi

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
#


/bin/bash runPre -v -t ${TRC} -l ${NLV} -I ${LABELANL} -s -n 2 -O
STATUS=$?
if [ ${STATUS} -ne 0 ];then
   exit ${STATUS}
fi


mv -f ${modelDataIn}/OZONSMT${LABELANL}S.grd.${postfix} ${modelDataIn}/OZON${PREFIX}${LABELANL}S.grd.${postfix}

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

/bin/bash runPre -v -t ${TRC} -l ${NLV} -I ${LABELANL} -p CPT -n 3
STATUS=$?
if [ ${STATUS} -ne 0 ];then
   exit ${STATUS}
fi

# Rodando o Modelo

/bin/bash runModel -das -v -np ${NPROC} -N 10 -d 4 -t ${TRC} -l ${NLV} -I ${LABELANL} -F ${LABELFCT} -W  ${LABELFCT} -p ${PREFIX} -s sstwkl -ts 3 -r -tr 6 -i -3
STATUS=$?
if [ ${STATUS} -ne 0 ];then
   exit ${STATUS}
fi

# Pos-processa as previsoes caso a variavel RUNPOS possua o valor Yes ou Y
if [ ${RUNPOS} == "yes" -o ${RUNPOS} == "y" ]
then

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
