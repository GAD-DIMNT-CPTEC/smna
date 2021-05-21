#! /bin/bash

# Descomentar para debugar
#set -o xtrace

# Carregando as variaveis do sistema
source /lustre_xc50/joao_gerd/SMG/config_smg.ksh vars_export

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

export MRES=`echo ${TRC} ${NLV} | awk '{printf("TQ%4.4dL%3.3d\n",$1,$2)}'`
export ANL="GANL${PREFIX}${LABELANL}S.unf.${MRES}"
export LABELFGS=`${inctime} ${LABELANL} +6h %y4%m2%d2%h2`
export YESTERDAY=$(${inctime} ${LABELANL} -1d %y4%m2%d200)
export runmodel=${subt_model_bam_exec}_${PREFIX}
export modelpos=${subt_pos_bam_exec}_${PREFIX}

echo -e ""
echo -e "\033[34;1m >> Submetendo o MCGA:\033[m \033[31;1m${LABELANL}\033[m"
echo -e " "
echo -e "\033[34;1m > Resolucao         : \033[m \033[31;1m${MRES}\033[m"
echo -e "\033[34;1m > Condicao Inicial  : \033[m \033[31;1m${LABELANL}\033[m"
echo -e "\033[34;1m > Previsao ate      : \033[m \033[31;1m${LABELFCT}\033[m"
echo -e "\033[34;1m > Pos-proc. Ativado : \033[m \033[31;1m${RUNPOS}\033[m"

#
# Copiando Analise do GSI para o Model/DataIn
#

cp -pfvr ${work_gsi_dataout}/${LABELANL}/GANL${PREFIX}${LABELANL}S.unf.${MRES} ${subt_model_bam_datain}

## Copiando/linkando demais arquivos do pré
#
## SST está disponível com 1 dia de atraso
#ln -s ${ncep_ext}/${YESTERDAY}/dataout/NCEP/rtgssthr_grb_0.083.grib2.${YESTERDAY:0:8} ${subt_pre_bam_datain}/
#
## Umidade do solo tem as 00UTC e as 12UTC
#DATA=${LABELANL}
#file="${ncep_ext}/${DATA}/dataout/Umid_Solo/GL_SM.GPNR.${DATA}.vfm"
#c=0
#until [ -e ${file} -o $c -gt 3 ]; do
#   DATA=$(${inctime} ${DATA} -6h %y4%m2%d2%h2)
#   file="${ncep_ext}/${DATA}/dataout/Umid_Solo/GL_SM.GPNR.${DATA}.vfm"
#   c=$((c+1))
#done
#if [ ! -e ${file} ];then
#   echo -ne "\033[32;1m Soil moisture file for\033[m: \033[34;1m${LABELANL}\033[m \033[32;1m Not Found\033[m"
#   exit -1
#else
#   ln -s ${file} ${subt_pre_bam_datain}/GL_SM.GPNR.${LABELANL}.vfm
#fi

# Rodando o Pre 
echo -ne "\033[34;1m > Pre-Processamento : \033[m"

cd ${home_smg}/run/scripts/bam_scripts

#echo "./runPre ${TRC} ${NLV} ${LABELANL} ${PREFIX} 0 F F 1534 64"
#/bin/bash runPre ${TRC} ${NLV} ${LABELANL} ${PREFIX} 0 F F ${TRC} ${NLV} #> /dev/null 2>&1
#/bin/bash runPre ${TRC} ${NLV} ${LABELANL} ${PREFIX} 0 F F 1534 64
echo "/bin/bash runPre -t ${TRC} -l ${NLV} -I ${LABELANL} -p CPT -n 3 -ti 299 -li 64"
/bin/bash runPre -t ${TRC} -l ${NLV} -I ${LABELANL} -p CPT -n 3 -ti 299 -li 64

if [ $? -eq 0 ]
then
  echo -e "\033[34;1m [\033[m\033[32;2m OK \033[m\033[34;1m]\033[m"
else
  echo -e "\033[34;1m[\033[m\033[31;1m Falhou \033[m\033[34;1m]\033[m"
  exit 1
fi

# Rodando o Modelo
cd ${home_smg}/run/scripts/bam_scripts

date

echo ". runModel -np ${NPROC} -N 10 -d 4 -t ${TRC} -l ${NLV} -I ${LABELANL} -F ${LABELFCT} -W  ${LABELFCT} -p ${PREFIX} -s sstwkl -ts 3 -r -tr 6 -i -3"
/bin/bash runModel -np ${NPROC} -N 10 -d 4 -t ${TRC} -l ${NLV} -I ${LABELANL} -F ${LABELFCT} -W  ${LABELFCT} -p ${PREFIX} -s sstwkl -ts 3 -r -tr 6 -i -3

if [ $? -eq 0 ]
then
 echo -e "\033[34;1m [\033[m\033[32;2m OK \033[m\033[34;1m]\033[m"
else
  echo -e "\033[34;1m[\033[m\033[31;1m Falhou \033[m\033[34;1m]\033[m"
  exit 1
fi

# Pos-processa as previsoes caso a variavel RUNPOS possua o valor Yes ou Y
if [ ${RUNPOS} == "yes" -o ${RUNPOS} == "y" ]
then

  # Verifica se o executavel se encontra em ${subt_pos_bam_run}
  if [ ! -e ${subt_pos_bam_run}/PostGrib ]
  then 
  
    cp -v ${home_pos_bam_source}/PostGrib ${subt_pos_bam_run}
  
  fi

  cd ${home_smg}/run/scripts/bam_scripts
  echo   "./runPos -np 120 -N 12 -d 1 -t ${TRC} -l ${NLV} -I ${LABELANL} -F ${LABELFCT} -p ${PREFIX} > /dev/null 2>&1"
  /bin/bash runPos -np 120 -N 12 -d 1 -t ${TRC} -l ${NLV} -I ${LABELANL} -F ${LABELFCT} -p ${PREFIX} > /dev/null 2>&1
fi

exit 0
