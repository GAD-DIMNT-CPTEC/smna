#! /bin/bash
#-----------------------------------------------------------------------------#
#           Group on Data Assimilation Development - GDAD/CPTEC/INPE          #
#-----------------------------------------------------------------------------#
#BOP
#
# !SCRIPT: script utilizado para extrair as observações para os ciclos de 
#          assimilação de dados
#
# !DESCRIPTION:
#
# !CALLING SEQUENCE:
#
#   ./run_obsmake.sh <opções>
#
#      As <opções> válidas são
#          * START_DATE : Data da condição inicial 
#
#          exemplo:
#          ./run_obsmake.sh 2015043006
# 
# !REVISION HISTORY: 
#
# !REMARKS:
#
#EOP
#-----------------------------------------------------------------------------#
#BOC

WhereIam=$(dirname ${BASH_SOURCE})

# Carregando as variaveis do sistema
source ${WhereIam}/../../config_smg.ksh vars_export

subwrd ( ) {
   str=$(echo "${@}" | awk '{ for (i=1; i<=NF-1; i++) printf("%s ",$i)}')
   n=$(echo "${@}" | awk '{ print $NF }')
   echo "${str}" | awk -v var=${n} '{print $var}'
}

# parse options

args=${@} #save arguments temporarily
while (( $# )); do
   opt=$1
   case ${opt} in
       -I) runDate=$2; shift 2;;
       -u) userSSH=$2; shift 2;;
       -h) cat < ${0} | sed -n '/^#BOP/,/^#EOP/p' ; exit 0;;
       *) echo -e "\033[31;1mWarning:\033[m Unknown argument:\033[33;1m $opt\033[m"; shift 1;;
   esac
done
set -- $args #restore arguments


if [ -z ${runDate} ];then
   echo -e "\033[31;1m runDate not set \033[m"
   exit 1
fi


 
# define parameters
parmGSI=${home_gsi_fix}/gsiparm.anl
obsGSI=${home_gsi_fix}/obsfiles.rc

# define dirs to find observations
obsDir=${ncep_ext}/${runDate}/dataout/NCEP
obsDir=${obsDir}:${subt_obs_run}/${runDate}


# define scp command
SCP="scp -pC"
if [ ! -z ${userSSH} ];then
   SCP+=" ${userSSH}"
fi
SCP+="@egeon:/oper/dados/preproc/brutos/model/gdas/"



IFS=":"; read -a obsPath < <(echo "${obsDir}")
nPaths=${#obsPath[@]}
IFS=" "

names=$(sed -n '/OBS_INPUT::/,/::/{/OBS_INPUT/d;/::/d;/^!/d;p}' ${parmGSI} | awk '{print $1}' | sort -u | xargs)
for name in ${names};do
   i=0
   while [ $i -le $((nPaths-1)) ];do
      filemask=$(grep -iw ${name} ${obsGSI} | awk '{print $1}'; exit ${PIPESTATUS[0]} )
      if [ $? -eq 0 ];then
         file=$(${inctime} ${runDate} +0h ${filemask})
         if [ -e ${obsPath[$i]}/${file} ];then break ;fi

         if [ $i -eq $((nPaths-1)) ];then
            dirOut=${subt_obs_run}/${runDate}
            dirIn=${runDate:0:4}/${runDate:4:2}/${runDate:6:2}
            if [ ! -e ${dirOut} ];then
               mkdir -p ${dirOut}
            fi
            ${SCP}${dirIn}/${file%%.${runDate:0:8}}* ${dirOut}/${file}
         fi
      fi
      i=$((i+1))
   done
done
