#!/bin/bash -x
#-----------------------------------------------------------------------------#
#           Group on Data Assimilation Development - GDAD/CPTEC/INPE          #
#-----------------------------------------------------------------------------#
#BOP
#
# !SCRIPT:
#
# !DESCRIPTION:
#
# !CALLING SEQUENCE:
#
# !REVISION HISTORY:
# 20 Dec 2017 - J. G. de Mattos - Initial Version
#
# !REMARKS:
#
#
#EOP
#-----------------------------------------------------------------------------#
#BOC

assign(){
#DESCRIPTION: Defines and exports an environment variable from the provided arguments
  eval export $1=$2
}

vars_export(){
#DESCRIPTION: Exports necessary environment variables from a configuration file specific to the HPC system in use
  FilePaths=$(dirname ${BASH_SOURCE})/mach/${hpc_name}_paths.conf
  while read line; do
    assign $line
  done < <(sed -r 's/^\s*(.*\S)*\s*$/\1/;/^$/d;/^#/d' ${FilePaths})
}

copy_fixed_files(){
#DESCRIPTION: Copies required fixed files for any run, adapting to the system in use (EGEON or Cray XC50)
  vars_export

  if [ ${HOSTNAME:0:1} = 'e' ] || [ ${hpc_name} = "egeon" ]; then
    cp -pf ${home_model_bam}/datain/* ${subt_model_bam}/datain/

    cp -pf ${public_bam}/MODEL/datain/AeroVar.Tab ${subt_model_bam}/datain/
    cp -pf ${public_bam}/MODEL/datain/ETAMPNEW_DATA ${subt_model_bam}/datain/
    cp -pf ${public_bam}/MODEL/datain/F_nwvl200_mu20_lam50_res64_t298_c080428.bin ${subt_model_bam}/datain/
    cp -pf ${public_bam}/MODEL/datain/iceoptics_c080917.bin ${subt_model_bam}/datain/
    cp -pf ${public_bam}/MODEL/datain/ocnalbtab24bnd.bin ${subt_model_bam}/datain/

    # Pre-processing
    echo "home_pre_bam: ${home_pre_bam} --->>> subt_pre_bam: ${subt_pre_bam}"

    cp -pf ${home_pre_bam}/datain/* ${subt_pre_bam}/datain/
    cp -pf ${home_pre_bam}/dataout/* ${subt_pre_bam}/dataout/

    cp -pf ${public_bam}/PRE/dataout/WaterNavy.dat ${subt_pre_bam}/dataout/
    cp -pf ${public_bam}/PRE/dataout/TopoNavy.dat ${subt_pre_bam}/dataout/
    cp -pf ${public_bam}/PRE/dataout/HPRIME.dat ${subt_pre_bam}/dataout/

    cp -pf ${home_pre_bam}/databcs/* ${subt_pre_bam}/databcs/

    cp -pf ${public_bam}/PRE/databcs/sib2soilms.form ${subt_pre_bam}/databcs/
    cp -pf ${public_bam}/PRE/databcs/FluxCO2.bin ${subt_pre_bam}/databcs/
    cp -pf ${public_bam}/PRE/databcs/FluxCO2.ctl ${subt_pre_bam}/databcs/
    cp -pf ${public_bam}/PRE/databcs/claymsk.form ${subt_pre_bam}/databcs/
    cp -pf ${public_bam}/PRE/databcs/clmt.form ${subt_pre_bam}/databcs/
    cp -pf ${public_bam}/PRE/databcs/deltat.form ${subt_pre_bam}/databcs/
    cp -pf ${public_bam}/PRE/databcs/ersst.bin ${subt_pre_bam}/databcs/
    cp -pf ${public_bam}/PRE/databcs/ibismsk.form ${subt_pre_bam}/databcs/
    cp -pf ${public_bam}/PRE/databcs/ndviclm.form ${subt_pre_bam}/databcs/
    cp -pf ${public_bam}/PRE/databcs/sandmsk.form ${subt_pre_bam}/databcs/
    cp -pf ${public_bam}/PRE/databcs/sib2msk.form ${subt_pre_bam}/databcs/
    cp -pf ${public_bam}/PRE/databcs/soiltext.form ${subt_pre_bam}/databcs/
  elif [ ${HOSTNAME:0:1} = 'c' ]; then
    ln -s /cray_home/joao_gerd/BAMFIX/model/datain/* ${subt_model_bam}/datain/
    ln -s /cray_home/joao_gerd/BAMFIX/pre/datain/* ${subt_pre_bam}/datain/
    ln -s /cray_home/joao_gerd/BAMFIX/pre/dataout/* ${subt_pre_bam}/dataout/
    ln -s /cray_home/joao_gerd/BAMFIX/pre/databcs/* ${subt_pre_bam}/databcs/
    ln -s cray_home/joao_gerd/BAMFIX/pre/dataco2/* ${subt_pre_bam}/dataco2/
    ln -s cray_home/joao_gerd/BAMFIX/pre/datasst/* ${subt_pre_bam}/datasst/
    ln -s cray_home/joao_gerd/BAMFIX/pre/dataTop/* ${subt_pre_bam}/dataTop/
  fi
}

configurar(){
#DESCRIPTION: Configures the modeling system, creating folders and symbolic links, as well as copying fixed files and modifying scripts and Makefiles as needed
  vars_export

  echo ""
  echo -e "\033[34;1m > The variable nome_smg has the value \033[36;1m${nome_smg}\033[m \033[m"
  echo -e "\033[34;1m > The SMG configurator will create folders with the name \033[36;1m${nome_smg}\033[m \033[m"
  echo -e "\033[34;1m > In the scratch[in,out] disks below. \033[m"
  echo ""
  echo -e "\033[34;1m > SMG HOME: \033[36;1m${home_smg}\033[m \033[m"
  echo -e "\033[34;1m > SMG SUBM: \033[36;1m${subt_smg}\033[m \033[m"
  echo -e "\033[34;1m > SMG WORK: \033[36;1m${work_smg}\033[m \033[m"

  if [ "/"${home_smg}"/" == "/"${HOME}/${nome_smg}"/" ]; then
     echo ""
     echo -e "\033[31;1m > The system detected that you are trying to install the \033[m"
     echo -e "\033[31;1m > SMG in home, which is not recommended, as there may be insufficient space. \033[m"
  fi
  echo ""

  read -p "Do you want to continue? (Y/N) " -n 1 -r
  echo    # (optional) move to a new line

  if [[ $REPLY = [Yy] ]]; then

    echo ""
    echo -e "\033[34;1m > Configuring the SMG... \033[m"
    echo ""

    # Creating directory structure
    if test ! -s ${subt_dataout}; then mkdir -p ${subt_dataout}; fi
    if test ! -e ${subt_bam}; then mkdir -p ${subt_bam}; fi
    if test ! -e ${subt_gsi}; then mkdir -p ${subt_gsi}; fi

    if test ! -e ${subt_obs_dataout}; then mkdir -p ${subt_obs_dataout}; fi

    if test ! -s ${subt_run}; then mkdir -p ${subt_run}; fi
    if test ! -s ${subt_run_gsi}; then mkdir -p ${subt_run_gsi}; fi
    if test ! -s ${subt_obs_run}; then mkdir -p ${subt_obs_run}; fi

    for dir in exec datain dataout; do
      if [ ! -e ${subt_model_bam}/${dir} ]; then
        mkdir -p ${subt_model_bam}/${dir}
      fi
    done

    for dir in exec datain dataout datasst databcs dataco2 dataTop; do
      if [ ! -e ${subt_pre_bam}/${dir} ]; then
        mkdir -p ${subt_pre_bam}/${dir}
      fi
    done

    for dir in exec datain dataout; do
      if [ ! -e ${subt_pos_bam}/${dir} ]; then
        mkdir -p ${subt_pos_bam}/${dir}
      fi
    done

    for dir in exec datain dataout; do
      if [ ! -e ${subt_grh_bam}/${dir} ]; then
        mkdir -p ${subt_grh_bam}/${dir}
      fi
    done

    if test ! -e ${subt_gsi}; then mkdir -p ${subt_gsi}; fi

    # Modify SMG scripts and Makefiles
    sed -i "/# Loading system variables/{n;d}" ${run_smg}/run_cycle.sh ${scripts_smg}/runGSI ${scripts_smg}/run_model.sh ${scripts_smg}/run_obsmake.sh ${scripts_smg}/run_blsdas.sh
    sed -i "/# Loading system variables/a\source "${home_smg}"\/config_smg.ksh vars_export" ${run_smg}/run_cycle.sh ${scripts_smg}/runGSI ${scripts_smg}/run_model.sh ${scripts_smg}/run_obsmake.sh ${scripts_smg}/run_blsdas.sh

    sed -i "/# Path where the Model executable should be located/{n;d}" ${home_model_bam}/source/Makefile
    sed -i "/# Path where the Model executable should be located/a\PATH2="${home_model_bam}/exec"" ${home_model_bam}/source/Makefile

    sed -i "/# Path where the Pos executables should be located/{n;d}" ${home_pos_bam}/source/Makefile
    sed -i "/# Path where the Pos executables should be located/a\PATH2="${home_pos_bam}/exec"" ${home_pos_bam}/source/Makefile

    sed -i "/DirInPut='\//,1d" ${home_run_bam}/PostGridHistory.nml
    sed -i "/!DirInPut=subt_model_bam_dataout/a\DirInPut=\'"${subt_model_bam}/dataout"\/TQ0042L028\'" ${home_run_bam}/PostGridHistory.nml

    sed -i "/DirOutPut='\//,1d" ${home_run_bam}/PostGridHistory.nml
    sed -i "/!DirOutPut=subt_grh_bam_dataout/a\DirOutPut=\'"${subt_grh_bam}/dataout"\/TQ0042L028\'" ${home_run_bam}/PostGridHistory.nml

    sed -i "/DirMain='\//,1d" ${home_run_bam}/PostGridHistory.nml
    sed -i "/!DirMain=subt_bam/a\DirMain=\'"${subt_bam}"\'" ${home_run_bam}/PostGridHistory.nml

    sed -i "/export HOMEBASE=\//,1d" ${home_run_bam}/EnvironmentalVariables
    sed -i "/# BAM path in HOME/a\export HOMEBASE="${home_bam}"" ${home_run_bam}/EnvironmentalVariables

    sed -i "/export SUBTBASE=\//,1d" ${home_run_bam}/EnvironmentalVariables
    sed -i "/# BAM path in scratchin (SUBMIT_HOME)/a\export SUBTBASE="${subt_bam}"" ${home_run_bam}/EnvironmentalVariables

    # Copy necessary files
    copy_fixed_files

    echo ""
    echo -e "\033[34;1m > SMG configuration complete \033[m"

    banner

  elif [[ ${REPLY} = [Nn] ]]; then
    echo -e "\033[34;1m > Exiting the configurator. \033[m"
    exit 0
  else
    exit 1
  fi
}

compilar(){
#DESCRIPTION: Compiles the complete system, including utilities, GSI, BAM, and other necessary components
  vars_export
  echo "Compiler : " $compiler

  if [ ! -e ${home_cptec}/bin ]; then
    mkdir -p ${home_cptec}/bin
  fi

  if [ ${HOSTNAME:0:1} = 'e' ]; then
    if [ ${HOSTNAME} != "eslogin01" -a ${HOSTNAME} != "eslogin02" ]; then
      echo "#####################################################################"
      echo "#                                                                   #"
      echo "#               You are logged in to ${HOSTNAME}                    #"
      echo "#                                                                   #"
      echo "# Before proceeding with Installation, log in to one of these servers: #"
      echo "#                                                                   #"
      echo "# $ ssh eslogin01 -XC                                               #"
      echo "#                                                                   #"
      echo "#  or                                                               #"
      echo "#                                                                   #"
      echo "# $ ssh eslogin02 -XC                                               #"
      echo "#                                                                   #"
      echo "#####################################################################"

      exit
    fi
  fi

  compgsi=1
  compang=1
  compbam=0

  echo ""
  echo "%%%"
  echo " Compiling SMG utilities:  "
  echo "%%%"

  cd ${home_gsi}

  echo ""
  echo "%%%%%%%%%%%%%%%%%%%%"
  echo " Compiling GSI:  "
  echo "%%%%%%%%%%%%%%%%%%%%"
  echo ""

  if [ ${compgsi} -eq 1 ]; then

    echo ""
    echo "+++++++++++++++++++++++++"
    echo "   Compiling GSI:     "
    echo "+++++++++++++++++++++++++"
    echo ""
    echo " Time: "  `date`
    echo " !!! Estimated compilation time: 30 minutes !!!"
    sleep 1
    echo ""

    cd ${home_gsi}
    pwd
    ./compile.sh -G -C ${compiler} 2>&1 | tee ${home_gsi}/compile.log

    if [ -e ${home_gsi_src}/gsi.x ]; then
      echo ""
      echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      echo "!                                                !"
      echo "!      GSI compilation successful                !"
      echo "!                                                !"
      echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      echo ""
    else
      echo ""
      echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      echo "!                                          !"
      echo "! GSI compilation failed !!!!               !"
      echo "! Check configurations or disk space!       !"
      echo "! More information in the compile.log file  !"
      echo "!                                          !"
      echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      echo ""
      exit -9
    fi
  fi

  if [ ${compang} -eq 1 ]; then

    if [ ${hpc_name} = 'XC50' ]; then
      source ./env.sh xc50 ${compiler}
    elif [ ${hpc_name} = 'egeon' ]; then
      source ./env.sh egeon ${compiler}
    fi

    echo ""
    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "   Compiling GSI bias correction utility:     "
    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

    cd ${home_gsi}/util/global_angupdate
    echo ln -sf Makefile.conf.${hpc_name}-${compiler} Makefile.conf
    ln -sf Makefile.conf.${hpc_name}-${compiler} Makefile.conf
    make -f Makefile clean
    make -f Makefile

    if [ -e ${home_gsi}/util/global_angupdate/global_angupdate ]; then
      echo ""
      echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      echo "!                                                                                           !"
      echo "!      GSI angle bias correction utility compilation successful                             !"
      echo "!                                                                                           !"
      echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      echo ""

      cp -pfvr ${home_gsi}/util/global_angupdate/global_angupdate ${home_cptec}/bin/global_angupdate

    else
      echo ""
      echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      echo "!                                                                      !"
      echo "! GSI angle bias correction utility compilation failed !!!!            !"
      echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      echo ""
      exit -9
    fi
  fi

  if [ ${compbam} -eq 1 ]; then
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%"
    echo " Compiling BAM:  "
    echo "%%%%%%%%%%%%%%%%%%%%"
    echo ""

    export mkname=${compiler}_${SUB}
    if [ ${hpc_name} = "egeon" ]; then
      module purge
      module load intel/2021.4.0
      module load mpi/2021.4.0 impi/2021.4.0
      module load netcdf/4.7.4 pnetcdf/1.12.2 netcdf-fortran/4.5.3
      module list
    fi

    echo "+++++++++++++++++++++++++"
    echo "   Compiling pre:        "
    echo "+++++++++++++++++++++++++"

    if [ ${SUB} = "cray" ]; then 
      export NETCDF_FORTRAN_DIR=${NETCDF_DIR}
    fi

    cd ${home_pre_bam}/build
    make clean ${mkname}
    make ${mkname}
    make install

    echo "+++++++++++++++++++++++++"
    echo "   Compiling pos:        "
    echo "+++++++++++++++++++++++++"

    cd ${home_pos_bam}/source
    make clean 
    make clean ${mkname}

    echo "+++++++++++++++++++++++++"
    echo "   Compiling model:      "
    echo "+++++++++++++++++++++++++"

    cd ${home_model_bam}/source
    make clean $mkname

    echo -e "\033[34;1m ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \033[m"
    echo -e "\033[34;1m > SMG compilation complete. Check for possible errors in the log file, if created. \033[m"
    echo -e "\033[34;1m ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \033[m"

  fi   

  cp -pfr ${home_gsi_src}/gsi.x ${home_cptec}/bin

  cp -pfr ${home_pre_bam}/exec/* ${subt_pre_bam}/exec
  cp -pfr ${home_model_bam}/exec/* ${subt_model_bam}/exec
  cp -pfr ${home_pos_bam}/exec/* ${subt_pos_bam}/exec

  echo -e "\033[34;1m > SMG compilation complete. Check for possible errors in the log file, if created. \033[m"

  banner
}

testcase(){
#DESCRIPTION: Unpacks testcase files and prepares the environment for test execution
  vars_export

  copy_fixed_files

  echo -e "\e[34;1m Choose one of the available options for the testcase:\e[m"
  i=1
  for line in $(ls -1 ${public_bam}/PRE/datain); do
    year[$i]=${line}
    opts[$i]=${i}
    echo -e "\e[31;1m [$i]\e[m\e[37;1m - Testcase for \e[m\e[32;1m${year[$i]}\e[m"
    i=$((i+1))
  done
  read answer

  anlfile=$(ls -1 ${public_bam}/PRE/datain/${year[$answer]}/ncep_anl/gdas1.*.SAnl.*|head -n 1)
  cp -pvfrL ${anlfile} ${subt_pre_bam}/datain/
  cp -pvfrL ${public_bam}/PRE/datain/${year[${answer}]}/sst/* ${subt_pre_bam}/datain/
  cp -pvfrL ${public_bam}/PRE/datain/${year[${answer}]}/sno/* ${subt_pre_bam}/datain/
  cp -pvfrL ${public_bam}/PRE/datain/${year[${answer}]}/smc/*.vfm ${subt_pre_bam}/datain/

  cp -pvfr ${public_bam}/PRE/dataout/* ${subt_pre_bam}/dataout/
  cp -pvfr ${public_bam}/PRE/databcs/* ${subt_pre_bam}/databcs/
  cp -pvfr ${public_bam}/PRE/datasst/* ${subt_pre_bam}/datasst/

  ln -s ${home_pos_bam}/datain/* ${subt_pos_bam}/datain/
}

banner(){
#DESCRIPTION: Prints information about the script and where to find more details
  echo -e ""
  echo -e "\e[34;1m > For more information about this SMG distribution, read the file: \e[m"
  echo -e "\e[32;1m > ${home_smg}/README \e[m"
  echo -e ""
}

ajuda(){
#DESCRIPTION: Displays help with available options
  echo ""
  echo " Usage.....: ${0##*/} <option>"
  echo ""
  First=1
  grep -i '(){$' ${BASH_SOURCE} | sed 's/(){//g' | while read function; do
    dsc=$(sed -n "/${function}(){$/ {n;p}" ${BASH_SOURCE} | sed "s/#*[DdEeSsCcRrIiPpTtIiOoNn]*://g")
    if [ ${First} -eq 1 ]; then
      First=0
      printf " Options..:%2s* \e[1;31m%s\e[m \e[1;37;1m-->\e[m\e[1;34m%s\e[m\n" " " "$function" "${dsc}"
    else
      printf "%12s* \e[1;31m%s\e[m \e[1;37;1m-->\e[m\e[1;34m%s\e[m\n" " " "$function" "${dsc}"
    fi
  done
  echo ""
  First=1
  grep -i '(){$' ${BASH_SOURCE} | sed 's/(){//g' | while read function; do
    dsc=$(sed -n "/${function}(){$/ {n;p}" ${BASH_SOURCE} | sed "s/#*[DdEeSsCcRrIiPpTtIiOoNn]*://g")
    if [ ${First} -eq 1 ]; then
      First=0
      printf " Examples:%2s\e[1;37m${0##*/}\e[m \e[1;31m%s\e[m\n" " " "$function"
    else
      printf "%12s\e[1;37m${0##*/}\e[m \e[1;31m%s\e[m\n" " " "$function"
    fi
  done
}

teste(){
#DESCRIPTION: Simple test function to verify the operation of functions
   echo "HELLO WORLD"
}

#EOC
#-----------------------------------------------------------------------------#

