#!/bin/bash
#-----------------------------------------------------------------------------#
#           Group on Data Assimilation Development - GDAD/CPTEC/INPE          #
#-----------------------------------------------------------------------------#
#BOP
#
# !SCRIPT: config_smg.sh
#
# !DESCRIPTION: Script utilizado para a configuracao e instacao do SMG
#
# !CALLING SEQUENCE: ./config_smg.sh <opcoes>
#                    nota: Para saber mais sobre as opcoes disponiveis
#                          execute ./config_smg.sh ajuda
#
# !REVISION HISTORY:
#
#    ?? ??? ???? - ??????????????? - Initial Version
#    20 Dec 2017 - J. G. de Mattos - separa em diferentes arquivos
#                                    * etc/paths.sh -> caminhos
#                                    * etc/functions.sh -> funcoes
#
#    Oct 2023 - J. A. Aravequia - acrescenta identificação do sistema de HPC
#                                 acrescenta a variável global hpc_name
#                                permite definições para diferentes computadores:
#                                    * etc/mach/XC50_paths.conf
#                                    * etc/mach/egeon_paths.conf
#                                    * etc/mach/new_system_avail.conf
#
# !REMARKS:
#
#    * Os caminhos do SMG estao contidos no arquivo etc/paths.sh
#    * As funcoes que podem ser executadas estao no arquivo etc/functions.sh
#
#EOP
#-----------------------------------------------------------------------------#
#BOC

RootDir=$(dirname ${BASH_SOURCE})

export SMG_ROOT=$RootDir

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
      echo "It's there."
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

. ${RootDir}/etc/functions.sh

#
# Verifica os argumentos passados junto com o script
#

echo -e ""
echo -e "\e[36;1m >>> ${BASH_SOURCE##*/} executado a partir de\e[m \e[32;1m${0##*/}\e[m"

if [ $# = 0 ];then
  echo -e ""
  echo -e "\e[31;1m > Nao foi passado nenhum argumento! \e[m"
  ajuda
  banner
  exit -1
fi

echo -en "\e[34;1m Opcao escolhida: \e[m \e[37;1m ${1} \e[m"

f=0
for function in $(grep -i '(){$' ${RootDir}/etc/functions.sh | sed 's/(){//g');do
   if [ ${1} == ${function} ];then
     f=1
     echo -e "\e[37;1m[\e[m\e[32;1m OK \e[m\e[37;1m]\e[m"
     vars_export
     ${1}
   fi
done

#------------------------------------------------------------#
# Pipes create SubShells, so the while read is running on a
# different shell than your script, that makes my f variable
# never changes (only the one inside the pipe subshell).
# Passing the data in a sub-shell instead, like it's a file
# before the while loop. This assumes that I don't want some
# intermittent file:

#while read function; do
#   if [ ${1} == ${function} ];then
#     f=1
#     echo -e "\e[37;1m[\e[m\e[32;1m OK \e[m\e[37;1m]\e[m"
#     vars_export
#     ${1}
#     banner
#   fi
#done < <(grep -i '(){$' ${RootDir}/etc/functions.sh | sed 's/(){//g')

#------------------------------------------------------------#

if [ ${f} -eq 0 ];then

  echo -e "\e[37;1m[\e[m\e[31;1m FAIL \e[m\e[37;1m]\e[m"
  echo -e ""
  echo -e "\e[37;1m Opcao desconhecida, <ajuda>: \e[m"
  vars_export
  ajuda
  banner
fi

#EOC
#-----------------------------------------------------------------------------#
