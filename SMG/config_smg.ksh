#!/bin/bash
#-----------------------------------------------------------------------------#
#           Group on Data Assimilation Development - GDAD/CPTEC/INPE          #
#-----------------------------------------------------------------------------#
#BOP
#
# !SCRIPT: config_smg.sh
#
# !DESCRIPTION: Script used for SMG configuration and installation
#
# !CALLING SEQUENCE: ./config_smg.sh <options>
#                    note: To learn more about the available options
#                          execute ./config_smg.sh help
#
# !REVISION HISTORY:
#
#    ?? ??? ???? - ??????????????? - Initial Version
#    20 Dec 2017 - J. G. de Mattos - split into different files
#                                    * etc/paths.sh -> paths
#                                    * etc/functions.sh -> functions
#
#    Oct 2023 - J. A. Aravequia - added HPC system identification
#                                 added the global variable hpc_name
#                                 allows definitions for different computers:
#                                    * etc/mach/XC50_paths.conf
#                                    * etc/mach/egeon_paths.conf
#                                    * etc/mach/new_system_avail.conf
#
# !REMARKS:
#
#    * SMG paths are contained in the etc/paths.sh file
#    * Functions that can be executed are in the etc/functions.sh file
#
#EOP
#-----------------------------------------------------------------------------#
#BOC

RootDir=$(dirname ${BASH_SOURCE})

export SMG_ROOT=${HOME}/SMNA_v3.0.0.t11889/SMG
echo "Installation path: SMG_ROOT=$SMG_ROOT"

lognode=$(cat /proc/sys/kernel/hostname | cut -b 1-6)

case $lognode in

  clogin)
    STR=$(uname -a)
    SUB='cray'
    if [[ "$STR" == *"$SUB"* ]]; then
      echo "This will run on Cray XC50 ..."
      export hpc_name="XC50"
    fi
    ;;

  headno)
    STR=$(uname -a)
    SUB='egeon'
    if [[ "$STR" == *"$SUB"* ]]; then
      echo "This will run on EGEON Cluster ..."
      export hpc_name="egeon"
    fi
    ;;

  *)
    mach=$(cat /proc/sys/kernel/hostname)
    echo "The configurations for $mach are not defined yet!"
    echo "1) Add the machine to the defined systems in etc/mach; and"
    echo "2) Add an option for it in the function copy_fixed_files in etc/functions"
    exit
    ;;
esac

. ${SMG_ROOT}/etc/functions.sh

# Check the arguments passed with the script
echo -e ""
echo -e "\e[36;1m >>> ${BASH_SOURCE##*/} executed from \e[m \e[32;1m${0##*/}\e[m"

if [ $# -eq 0 ]; then
  echo -e ""
  echo -e "\e[31;1m > No arguments were passed! \e[m"
  help
  banner
  exit -1
fi

echo -en "\e[34;1m Selected option: \e[m \e[37;1m ${1} \e[m"

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
  echo -e "\e[37;1m Unknown option, <help>: \e[m"
  vars_export
  help
  banner
fi

#EOC
#-----------------------------------------------------------------------------#

