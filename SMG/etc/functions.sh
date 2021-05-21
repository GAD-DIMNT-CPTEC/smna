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
Assing(){
eval export $1=$2
}

# Variaveis principais
vars_export(){
#DESCRIPTION: Exporta as variáveis de ambiente do SMG
FilePaths=$(dirname ${BASH_SOURCE})/paths.conf
while read line;do

   Assing $line

done < <(sed -r 's/^\s*(.*\S)*\s*$/\1/;/^$/d;/^#/d' ${FilePaths})
}

copy_fixed_files(){
#DESCRIPTION: Copia os arquivos fixos necessarios para qualquer rodada

  vars_export

  #
  # Copiando arquivos fixos necessários para qualquer rodada, independente do período
  # estes arquivos não estão sob controle de versão pois são binários

  #
  # BAM
  # 

  #-----------------------------------------------------------------------------------------
  # model/datain
  if [ ${HOSTNAME:0:1} = 'e' ];then
     cp -pf ${home_model_bam_datain}/* ${subt_model_bam_datain}/
   
     cp -pf ${public_bam}/MODEL/datain/AeroVar.Tab ${subt_model_bam_datain}/
     cp -pf ${public_bam}/MODEL/datain/ETAMPNEW_DATA ${subt_model_bam_datain}/
     cp -pf ${public_bam}/MODEL/datain/F_nwvl200_mu20_lam50_res64_t298_c080428.bin ${subt_model_bam_datain}/
     cp -pf ${public_bam}/MODEL/datain/iceoptics_c080917.bin ${subt_model_bam_datain}/
     cp -pf ${public_bam}/MODEL/datain/ocnalbtab24bnd.bin ${subt_model_bam_datain}/
   
     #-----------------------------------------------------------------------------------------
     # pre/datain
   
     cp -pf ${home_pre_bam_datain}/* ${subt_pre_bam_datain}/
   
     #-----------------------------------------------------------------------------------------
     # pre/dataout
   
     cp -pf ${home_pre_bam_dataout}/* ${subt_pre_bam_dataout}/
   
     cp -pf ${public_bam}/PRE/dataout/WaterNavy.dat ${subt_pre_bam_dataout}/
     cp -pf ${public_bam}/PRE/dataout/TopoNavy.dat ${subt_pre_bam_dataout}/
     cp -pf ${public_bam}/PRE/dataout/HPRIME.dat ${subt_pre_bam_dataout}/
   
     #-----------------------------------------------------------------------------------------
     # pre/databcs
   
     cp -pf ${home_pre_bam_databcs}/* ${subt_pre_bam_databcs}/
   
     cp -pf ${public_bam}/PRE/databcs/sib2soilms.form ${subt_pre_bam_databcs}/
     cp -pf ${public_bam}/PRE/databcs/FluxCO2.bin ${subt_pre_bam_databcs}/
     cp -pf ${public_bam}/PRE/databcs/FluxCO2.ctl ${subt_pre_bam_databcs}/
     cp -pf ${public_bam}/PRE/databcs/claymsk.form ${subt_pre_bam_databcs}/
     cp -pf ${public_bam}/PRE/databcs/clmt.form ${subt_pre_bam_databcs}/
     cp -pf ${public_bam}/PRE/databcs/deltat.form ${subt_pre_bam_databcs}/
     cp -pf ${public_bam}/PRE/databcs/ersst.bin ${subt_pre_bam_databcs}/
     cp -pf ${public_bam}/PRE/databcs/ibismsk.form ${subt_pre_bam_databcs}/
     cp -pf ${public_bam}/PRE/databcs/ndviclm.form ${subt_pre_bam_databcs}/
     cp -pf ${public_bam}/PRE/databcs/sandmsk.form ${subt_pre_bam_databcs}/
     cp -pf ${public_bam}/PRE/databcs/sib2msk.form ${subt_pre_bam_databcs}/
     cp -pf ${public_bam}/PRE/databcs/soiltext.form ${subt_pre_bam_databcs}/
  elif [ ${HOSTNAME:0:1} = 'c' ];then

     cp -pf ${home_model_bam_datain}/* ${subt_model_bam_datain}/
   
     scp ${USER//_/.}@tupa:${public_bam}/MODEL/datain/AeroVar.Tab ${subt_model_bam_datain}/
     scp ${USER//_/.}@tupa:${public_bam}/MODEL/datain/ETAMPNEW_DATA ${subt_model_bam_datain}/
     scp ${USER//_/.}@tupa:${public_bam}/MODEL/datain/F_nwvl200_mu20_lam50_res64_t298_c080428.bin ${subt_model_bam_datain}/
     scp ${USER//_/.}@tupa:${public_bam}/MODEL/datain/iceoptics_c080917.bin ${subt_model_bam_datain}/
     scp ${USER//_/.}@tupa:${public_bam}/MODEL/datain/ocnalbtab24bnd.bin ${subt_model_bam_datain}/
   
     #-----------------------------------------------------------------------------------------
     # pre/datain
   
     cp -pf ${home_pre_bam_datain}/* ${subt_pre_bam_datain}/
   
     #-----------------------------------------------------------------------------------------
     # pre/dataout
   
     cp -pf ${home_pre_bam_dataout}/* ${subt_pre_bam_dataout}/
   
     scp ${USER//_/.}@tupa:${public_bam}/PRE/dataout/WaterNavy.dat ${subt_pre_bam_dataout}/
     scp ${USER//_/.}@tupa:${public_bam}/PRE/dataout/TopoNavy.dat ${subt_pre_bam_dataout}/
     scp ${USER//_/.}@tupa:${public_bam}/PRE/dataout/HPRIME.dat ${subt_pre_bam_dataout}/
   
     #-----------------------------------------------------------------------------------------
     # pre/databcs
   
     cp -pf ${home_pre_bam_databcs}/* ${subt_pre_bam_databcs}/
   
     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/sib2soilms.form ${subt_pre_bam_databcs}/
     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/FluxCO2.bin ${subt_pre_bam_databcs}/
     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/FluxCO2.ctl ${subt_pre_bam_databcs}/
     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/claymsk.form ${subt_pre_bam_databcs}/
     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/clmt.form ${subt_pre_bam_databcs}/
     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/deltat.form ${subt_pre_bam_databcs}/
     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/ersst.bin ${subt_pre_bam_databcs}/
     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/ibismsk.form ${subt_pre_bam_databcs}/
     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/ndviclm.form ${subt_pre_bam_databcs}/
     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/sandmsk.form ${subt_pre_bam_databcs}/
     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/sib2msk.form ${subt_pre_bam_databcs}/
     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/soiltext.form ${subt_pre_bam_databcs}/
  fi

}

# Funcao configurar
configurar(){
#DESCRIPTION: Configura o Sistema de modelagem (cria links e pastas)

# Exportando variaveis
  vars_export

  echo ""
  echo -e "\033[34;1m > A variavel nome_smg possui o valor \033[36;1m${nome_smg}\033[m \033[m"
  echo -e "\033[34;1m > o configurador do SMG ira criar pastas com o nome \033[36;1m${nome_smg}\033[m \033[m"
  echo -e "\033[34;1m > nos discos scratch[in,out] abaixo. \033[m"
  echo ""
  echo -e "\033[34;1m > SMG HOME: \033[36;1m${home_smg}\033[m \033[m"
  echo -e "\033[34;1m > SMG SUBM: \033[36;1m${subt_smg}\033[m \033[m"
  echo -e "\033[34;1m > SMG WORK: \033[36;1m${work_smg}\033[m \033[m" 
  if [ "/"${home_smg}"/" == "/"${HOME}/${nome_smg}"/" ] 
  then 
     echo ""
     echo -e "\033[31;1m > O sistema detectou que voce esta tentando instalar o \033[m"
     echo -e "\033[31;1m > SMG no home, o que nao e recomendado, pois podera haver falta de espaco. \033[m"
  fi
  echo ""

  read -p "Deseja continuar? (S/N) " -n 1 -r
  echo    # (optional) move to a new line

  if [[ $REPLY = [Ss] ]];then

# Mensagem ao usuario
  echo ""
  echo -e "\033[34;1m > Configurando o SMG ... \033[m"
  echo ""

# Criando estrutura de diretorios

# Cria datainout em diante (bam, gsi, obs, interface)
  if test ! -s ${subt_dataout}; then mkdir -p ${subt_dataout}; fi
  if test ! -s ${work_dataout}; then mkdir -p ${work_dataout}; fi

#  if test ! -s ${work_bam}; then mkdir -p ${work_bam}; fi
#  if test ! -e ${subt_bam}; then ln -sf ${work_bam} ${subt_bam}; fi
  if test ! -e ${subt_bam}; then mkdir -p ${subt_bam}; fi
  if test ! -s ${work_bam}; then mkdir -p ${work_bam}; fi

#  if test ! -s ${work_gsi}; then mkdir -p ${work_gsi}; fi
#  if test ! -e ${subt_gsi}; then ln -sf ${work_gsi} ${subt_gsi}; fi
  if test ! -e ${subt_gsi}; then mkdir -p ${subt_gsi}; fi
  if test ! -s ${work_gsi}; then mkdir -p ${work_gsi}; fi

  if test ! -s ${work_gsi_dataout}; then mkdir -p ${work_gsi_dataout}; fi
  if test ! -e ${subt_gsi_dataout}; then ln -sf ${work_gsi_dataout} ${subt_gsi_dataout}; fi

#  if test ! -s ${work_gsi_datain}; then mkdir -p ${work_gsi_datain}; fi
#  if test ! -e ${subt_gsi_datain}; then ln -sf ${work_gsi_datain} ${subt_gsi_datain}; fi
#  if test ! -e ${subt_gsi_datain}; then mkdir -p ${subt_gsi_datain}; fi
#  if test ! -s ${subt_gsi_datain_obs}; then mkdir -p ${subt_gsi_datain_obs}; fi

#  if test ! -s ${work_gsi_datain}; then ln -sf ${subt_gsi_datain} ${work_gsi_datain}; fi

  if test ! -s ${work_gsi_datain_obs}; then mkdir -p ${work_gsi_datain_obs}; fi
  if test ! -s ${subt_gsi_datain}; then ln -s ${work_gsi_datain} ${subt_gsi_datain}; fi
#  if test ! -s ${work_gsi_datain_bkg}; then mkdir -p ${work_gsi_datain_bkg}; fi

#  if test ! -s ${work_obs_dataout}; then mkdir -p ${work_obs_dataout}; fi
#  if test ! -e ${subt_obs_dataout}; then ln -sf ${work_obs_dataout} ${subt_obs_dataout}; fi
  if test ! -e ${subt_obs_dataout}; then mkdir -p ${subt_obs_dataout}; fi
  if test ! -s ${work_obs_dataout}; then ln -sf ${subt_obs_dataout} ${work_obs_dataout}; fi

# Cria datarun em diante (bam, gsi, obs, interface)
  if test ! -s ${subt_run}; then mkdir -p ${subt_run}; fi
  if test ! -s ${work_run}; then mkdir -p ${work_run}; fi

  if test ! -s ${subt_run_gsi}; then mkdir -p ${subt_run_gsi}; fi
  if test ! -e ${work_run_gsi}; then ln -sf ${subt_run_gsi} ${work_run_gsi}; fi

#  if test ! -s ${subt_sfc_run}; then mkdir -p ${subt_sfc_run}; fi
#  if test ! -e ${work_sfc_run}; then ln -sf ${subt_sfc_run} ${work_sfc_run}; fi

  if test ! -s ${subt_obs_run}; then mkdir -p ${subt_obs_run}; fi
  if test ! -e ${work_obs_run}; then ln -sf ${subt_obs_run} ${work_obs_run}; fi

# Cria pastas do MCGA
  if test ! -s ${work_pre_bam}; then mkdir -p ${work_pre_bam}; fi
  if test ! -s ${subt_pre_bam_datain}; then mkdir -p ${subt_pre_bam_datain}; fi
  if test ! -e ${work_pre_bam_datain}; then ln -sf ${subt_pre_bam_datain} ${work_pre_bam_datain}; fi
  if test ! -s ${work_pre_bam_dataout}; then mkdir -p ${work_pre_bam_dataout}; fi
  if test ! -e ${subt_pre_bam_dataout}; then ln -sf ${work_pre_bam_dataout} ${subt_pre_bam_dataout}; fi

  if test ! -s ${subt_pre_bam_datasst}; then mkdir -p ${subt_pre_bam_datasst}; fi
#  if test ! -s ${oiv2monthly}; then mkdir -p ${oiv2monthly}; fi

  if test ! -s ${subt_pre_bam_run}; then mkdir -p ${subt_pre_bam_run}; fi
  if test ! -e ${work_pre_bam_run}; then ln -sf ${subt_pre_bam_run} ${work_pre_bam_run}; fi
  if test ! -s ${subt_pre_bam_databcs}; then mkdir -p ${subt_pre_bam_databcs}; fi
  if test ! -e ${work_pre_bam_databcs}; then ln -sf ${subt_pre_bam_databcs} ${work_pre_bam_databcs}; fi

  if test ! -s ${work_model_bam}; then mkdir -p ${work_model_bam}; fi
  if test ! -s ${subt_model_bam_datain}; then mkdir -p ${subt_model_bam_datain}; fi
  if test ! -e ${work_model_bam_datain}; then ln -sf ${subt_model_bam_datain} ${work_model_bam_datain}; fi
  if test ! -s ${work_model_bam_dataout}; then mkdir -p ${work_model_bam_dataout}; fi
  if test ! -e ${subt_model_bam_dataout}; then ln -sf ${work_model_bam_dataout} ${subt_model_bam_dataout}; fi
  if test ! -s ${subt_model_bam_run}; then mkdir -p ${subt_model_bam_run}; fi
  if test ! -e ${work_model_bam_run}; then ln -sf ${subt_model_bam_run} ${work_model_bam_run}; fi

  if test ! -s ${work_pos_bam}; then mkdir -p ${work_pos_bam}; fi
  if test ! -s ${subt_pos_bam_datain}; then mkdir -p ${subt_pos_bam_datain}; fi
  if test ! -e ${work_pos_bam_datain}; then ln -sf ${subt_pos_bam_datain} ${work_pos_bam_datain}; fi
  if test ! -s ${work_pos_bam_dataout}; then mkdir -p ${work_pos_bam_dataout}; fi
  if test ! -e ${subt_pos_bam_dataout}; then ln -sf ${work_pos_bam_dataout} ${subt_pos_bam_dataout}; fi
  if test ! -s ${subt_pos_bam_run}; then mkdir -p ${subt_pos_bam_run}; fi
  if test ! -e ${work_pos_bam_run}; then ln -sf ${subt_pos_bam_run} ${work_pos_bam_run}; fi

  if test ! -s ${subt_grh_bam}; then mkdir -p ${subt_grh_bam}; fi
  if test ! -s ${work_grh_bam_dataout}; then mkdir -p ${work_grh_bam_dataout}; fi
  if test ! -e ${subt_grh_bam_dataout}; then ln -sf ${work_grh_bam_dataout} ${subt_grh_bam_dataout}; fi
  if test ! -s ${subt_grh_bam_run}; then mkdir -p ${subt_grh_bam_run}; fi
  if test ! -e ${work_grh_bam_run}; then ln -sf ${subt_grh_bam_run} ${work_grh_bam_run}; fi

  if test ! -e ${subt_gsi}; then mkdir -p ${subt_gsi}; fi

# Pasta de saida do GSI para os arquivos do testcase
#  if test ! -e ${subt_gsi_dataout}/2013123118; then mkdir -p ${subt_gsi_dataout}/2013123118; fi

# Altera scripts e Makefiles do SMG
  sed -i "/# Carregando as variaveis do sistema/{n;d}" ${run_smg}/run_cycle.sh ${scripts_smg}/runGSI ${scripts_smg}/run_model.sh ${scripts_smg}/run_obsmake.sh ${scripts_smg}/run_blsdas.sh
  sed -i "/# Carregando as variaveis do sistema/a\source "${home_smg}"\/config_smg.ksh vars_export"  ${run_smg}/run_cycle.sh ${scripts_smg}/runGSI ${scripts_smg}/run_model.sh ${scripts_smg}/run_obsmake.sh ${scripts_smg}/run_blsdas.sh

  sed -i "/# Caminho onde deverao ficar os arquivos executaveis do Pre/{n;d}" ${home_pre_bam_source}/Makefile
  sed -i "/# Caminho onde deverao ficar os arquivos executaveis do Pre/a\PATH2="${home_pre_bam_run}"" ${home_pre_bam_source}/Makefile

  sed -i "/# Caminho onde devera fica o arquivo executavel do Model/{n;d}" ${home_model_bam_source}/Makefile
  sed -i "/# Caminho onde devera fica o arquivo executavel do Model/a\PATH2="${home_model_bam_run}"" ${home_model_bam_source}/Makefile

  sed -i "/# Caminho onde deverao ficar os arquivos executaveis do Pos/{n;d}" ${home_pos_bam_source}/Makefile
  sed -i "/# Caminho onde deverao ficar os arquivos executaveis do Pos/a\PATH2="${home_pos_bam_run}"" ${home_pos_bam_source}/Makefile

  sed -i "/DirInPut='\//,1d" ${home_run_bam}/PostGridHistory.nml
  sed -i "/!DirInPut=subt_model_bam_dataout/a\DirInPut=\'"${subt_model_bam_dataout}"\/TQ0042L028\'" ${home_run_bam}/PostGridHistory.nml

  sed -i "/DirOutPut='\//,1d" ${home_run_bam}/PostGridHistory.nml
  sed -i "/!DirOutPut=subt_grh_bam_dataout/a\DirOutPut=\'"${subt_grh_bam_dataout}"\/TQ0042L028\'" ${home_run_bam}/PostGridHistory.nml

  sed -i "/DirMain='\//,1d" ${home_run_bam}/PostGridHistory.nml
  sed -i "/!DirMain=subt_bam/a\DirMain=\'"${subt_bam}"\'" ${home_run_bam}/PostGridHistory.nml

  sed -i "/export HOMEBASE=\//,1d" ${home_run_bam}/EnvironmentalVariables
  sed -i "/# Caminho do BAM no HOME/a\export HOMEBASE="${home_bam}"" ${home_run_bam}/EnvironmentalVariables

  sed -i "/export SUBTBASE=\//,1d" ${home_run_bam}/EnvironmentalVariables
  sed -i "/# Caminho do BAM no scratchin (SUBMIT_HOME)/a\export SUBTBASE="${subt_bam}"" ${home_run_bam}/EnvironmentalVariables

  sed -i "/export WORKBASE=\//,1d" ${home_run_bam}/EnvironmentalVariables
  sed -i "/# Caminho do BAM no scratchout (SUBMIT_WORK)/a\export WORKBASE="${work_bam}"" ${home_run_bam}/EnvironmentalVariables

  sed -i "/export PATHBASE=\//,1d" ${home_run_bam}/EnvironmentalVariablesMCGA
  sed -i "/#export PATHBASE=home_bam/a\export PATHBASE="${home_bam}"" ${home_run_bam}/EnvironmentalVariablesMCGA

  sed -i "/echo \${DK}; export DK=\//,1d" ${home_run_bam}/EnvironmentalVariablesMCGA
  sed -i "/#echo \${DK}; export DK=subt_bam/a\echo \${DK}; export DK="${subt_bam}"" ${home_run_bam}/EnvironmentalVariablesMCGA

  sed -i "/echo \${DK2}; export DK2=\//,1d" ${home_run_bam}/EnvironmentalVariablesMCGA
  sed -i "/#echo \${DK2}; export DK2=subt_bam/a\echo \${DK2}; export DK2="${subt_bam}"" ${home_run_bam}/EnvironmentalVariablesMCGA

  # Copiando arquivos necessários
  copy_fixed_files

  echo ""
  echo -e "\033[34;1m > Configuracao SMG completa \033[m"

  banner

elif [[ ${resposta} = N || ${resposta} = "n" ]]
then
  echo -e "\033[34;1m > Saindo do configurador. \033[m"
  exit 0
else
  exit 1
fi

}

# Funcao compilar
compilar(){
#DESCRIPTION: Compila o sistema completo (cria executaveis)

# Exportando variaveis principais
   vars_export

# Criando diretório onde estarão os executáveis
   if [ ! -e ${home_cptec_bin} ];then
       mkdir -p ${home_cptec_bin}
   fi
# Verificando se esta logado no eslogin01
   if [ ${HOSTNAME:0:1} = 'e' ];then
      if [ ${HOSTNAME} != "eslogin01" -a ${HOSTNAME} != "eslogin02" ];then
           echo "#####################################################################"
           echo "#                                                                   #"
           echo "#               Voce esta logado no ${HOSTNAME}                       #"
           echo "#                                                                   #"
           echo "# Antes de proceder com a Instalacao logar em um destes servidores: #"
           echo "#                                                                   #"
           echo "# $ ssh eslogin01 -XC                                               #"
           echo "#                                                                   #"
           echo "#  ou                                                               #"
           echo "#                                                                   #"
           echo "# $ ssh eslogin02 -XC                                               #"
           echo "#                                                                   #"
           echo "#####################################################################"
      
           exit
      fi
   fi
 
   echo "" 
   echo "%%%"
   echo " Compilando utilitarios do SMG:  "
   echo "%%%"


## Compilando utilitario inctime
   echo ""
   echo "%%% Inctime:"
   cd ${util_inctime} 
   make

## Compilando utilitário spectral to grid
   echo ""
   echo "%%% spc2grd:"
   cd ${util_spc2grd}
   ./compile.sh

#
#############################################################################
# Compilacao gsi
##############################################################################
   
   cd ${home_gsi}

   echo ""
   echo "%%%%%%%%%%%%%%%%%%%%"
   echo " Compilando o GSI:  "
   echo "%%%%%%%%%%%%%%%%%%%%"
   echo ""
# Pacote NETCDF
   if [ ${HOSTNAME:0:1} = 'e' ];then
      . /opt/modules/default/etc/modules.sh
      module load netcdf/4.2.0
   fi

# NETCDF_DIR=/opt/cray/netcdf/4.2.0/pgi/119
   export NETCDF=${NETCDF_DIR}
   echo "Diretorio NETCDF:    " $NETCDF

# Bilioteca LAPACK
   if [ ${HOSTNAME:0:1} = 'e' ];then
      export LAPACK_PATH=/opt/lapack/3.1.1
      echo "Bilioteca LAPACK:    " $LAPACK_PATH
   fi

# Configurando o GSI
   echo ""
   echo "+++++++++++++++++++++++++"
   echo "Limpando o diretorio GSI:"
   echo "+++++++++++++++++++++++++"
   echo ""
   echo " com ./clean -a"
   echo "" 
   ./clean -a

   if [ ${HOSTNAME:0:1} = 's' ];then
      compiler=PGI
   elif [ ${HOSTNAME:0:1} = 'c' ];then
      compiler=gfortran
   fi

   echo "1" | perl arch/Config.pl -corepath=$(pwd) -os=Cray -mach=x86_64 -comp=${compiler} -usewrf=0 -netcdf=${NETCDF}
 
# Compilando o GSI
   echo ""
   echo "+++++++++++++++++++++++++"
   echo "   Compilando o GSI:     "
   echo "+++++++++++++++++++++++++" 
   echo ""
   echo " Time: "  `date`
   echo " !!! Tempo previsto para compilacao 30 minutos !!!"
   sleep 1  
   echo ""
   
   cd ${home_gsi}
   ./compile 2>&1 | tee ${home_gsi}/compile.log
 
   if [ -e ${home_gsi_src}/gsi.exe  ]; then
     echo ""
     echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     echo "!                                                !"
     echo "!      Compilacao do GSI com Sucesso             !"
     echo "!                                                !"
     echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     echo ""
   else 
     echo ""
     echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     echo "!                                          !"
     echo "! Compilacao do GSI falhou !!!!            !"
     echo "! Rever configuracoes ou espaço em disco!  !"
     echo "! Mais informacoes no arquivo compile.log  !"
     echo "!                                          !"
     echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     echo ""
     exit -9
   fi
 
# Compilando do utilitario para geracão de bias relacionado com o ângulo
 
   echo ""
   echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
   echo "   Compilando utilitário de gereção de correção de bias do GSI:     "
   echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" 
 
   cd ${home_gsi}/util/gsi_angupdate
   make
 
 
   if [ -e ${home_gsi}/util/gsi_angupdate/gsi_angupdate.exe  ]; then
     echo ""
     echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     echo "!                                                                                           !"
     echo "!      Compilacao do utilitário de gereção de correção de bias do ângulo do GSI com Sucesso !"
     echo "!                                                                                           !"
     echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     echo ""
 
     cp -pfvr ${home_gsi}/util/gsi_angupdate/gsi_angupdate.exe ${home_cptec_bin}/gsi_angupdate.exe
 
   else 
     echo ""
     echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     echo "!                                                                      !"
     echo "! Compilacao do utilitário de gereção de correção de bias  falhou !!!! !"
     echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     echo ""
     exit -9
   fi


#############################################################################
# Compilacao bam
#############################################################################

   echo ""
   echo "%%%%%%%%%%%%%%%%%%%%"
   echo " Compilando o BAM:  "
   echo "%%%%%%%%%%%%%%%%%%%%"
   echo ""


#########################################
# Compilacao pre/pos
#########################################

#  . /opt/modules/default/etc/modules.sh
#  module swap pgi pgi/11.10.0

   echo "+++++++++++++++++++++++++"
   echo "   Compilando o pre:     "
   echo "+++++++++++++++++++++++++" 
 
   cd ${home_pre_bam_source}
   make clean gnu_cray
   
   echo "+++++++++++++++++++++++++"
   echo "   Compilando o pos:     "
   echo "+++++++++++++++++++++++++"   
   
   cd ${home_pos_bam_source}
   make clean gnu_cray
  
  
#########################################
# Compilacao model
#########################################

   echo "+++++++++++++++++++++++++"
   echo "   Compilando o model:   "
   echo "+++++++++++++++++++++++++"   

#  . /opt/modules/default/etc/modules.sh
#  module swap PrgEnv-pgi PrgEnv-cray

   cd ${home_model_bam_source}
#   make clean cray_cray32
   make clean gnu_cray


   echo -e "\033[34;1m ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \033[m"
   echo -e "\033[34;1m > Compilacao SMG completa.# Verifique possiveis erros no arquivo de log, caso o tenha criado. \033[m"
   echo -e "\033[34;1m ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \033[m"
#########################################
# Compilacao do blsdas
#########################################

   echo "+++++++++++++++++++++++++"
   echo "   Compilando o BLSDAS:   "
   echo "+++++++++++++++++++++++++"   
 
   . /opt/modules/default/etc/modules.sh
   module swap PrgEnv-pgi PrgEnv-gnu

   export ARCH=linux_gnu
   
   cd ${home_blsdas}
   make clean
   make

#   # compilando bibliotecas necessárias
#   for lib in mpeu bufr w3lib sigioBAM BilinInterp;do
#      cd ${home_blsdas_lib}/${lib}
#      make clean
#      make
#   done
# 
#   # compilando codigo fonte
#   cd ${home_blsdas_src}
#   make clean
#   make

#########################################
# Copiando executaveis
#########################################

# Copia todos os executaveis para a pasta cptec/bin


   cp -pfr ${home_gsi_src}/gsi.exe ${home_cptec_bin}
   cp -pfr ${util_inctime}/inctime ${home_cptec_bin}
   cp -pfr ${util_spc2grd}/spc2grd.x ${home_cptec_bin}

#
# exec BAM
#
   cp -pfr ${home_pre_bam_run}/* ${subt_pre_bam_run}
   cp -pfr ${home_model_bam_run}/* ${subt_model_bam_run}
   cp -pfr ${home_pos_bam_run}/* ${subt_pos_bam_run}



   echo -e "\033[34;1m > Compilacao SMG completa. Verifique possiveis erros no arquivo de log, caso o tenha criado. \033[m"

   banner

}

# Funcao aloca tescase
testcase(){
#DESCRIPTION: Descompacta os arquivos do TestCase

  vars_export

  #
  # copiando arquivos fixos (Já devem estar pois são copiados durante a configuracão)
  #

  copy_fixed_files

  #
  # Pre-processamento BAM
  #
  echo -e "\e[34;1m Escolha uma das opções disponíveis para o testcase:\e[m"
  i=1
  for line in $(ls -1 ${public_bam}/PRE/datain);do
     year[$i]=${line}
     opts[$i]=${i}
     echo -e "\e[31;1m [$i]\e[m\e[37;1m - Testcase para \e[m\e[32;1m${year[$i]}\e[m"
     i=$((i+1))
  done
  read answer

  anlfile=$(ls -1 ${public_bam}/PRE/datain/${year[$answer]}/ncep_anl/gdas1.*.SAnl.*|head -n 1)
  cp -pvfrL ${anlfile}  ${subt_pre_bam_datain}/
  cp -pvfrL ${public_bam}/PRE/datain/${year[${answer}]}/sst/* ${subt_pre_bam_datain}/
  cp -pvfrL ${public_bam}/PRE/datain/${year[${answer}]}/sno/* ${subt_pre_bam_datain}/
  cp -pvfrL ${public_bam}/PRE/datain/${year[${answer}]}/smc/*.vfm ${subt_pre_bam_datain}/

  #
  # independentes
  #

  cp -pvfr ${public_bam}/PRE/dataout/* ${subt_pre_bam_dataout}/
  cp -pvfr ${public_bam}/PRE/databcs/* ${subt_pre_bam_databcs}/
  cp -pvfr ${public_bam}/PRE/datasst/* ${subt_pre_bam_datasst}/

  ln -s ${home_pos_bam_datain}/* ${subt_pos_bam_datain}/

  #
  # GSI fix
  #

#  cp -v ${public_fix}/xrt004.ana_satbang_rst.20120417_09z.txt ${subt_gsi_dataout}/2012123118/satbias_ang.out
#  cp -v ${public_fix}/xrt004.ana_satbias_rst.20120417_09z.txt ${subt_gsi_dataout}/2012123118/satbias_out

}

# Final
banner(){
#DESCRIPTION: Imprime algumas informacoes sobre o script

#  vars_export

  echo -e ""
  echo -e "\e[34;1m > Para mais informacoes sobre esta distribuicao do SMG, leia o arquivo: \e[m"
  echo -e "\e[32;1m > ${home_smg}/README \e[m"
  echo -e ""

}
ajuda(){
#DESCRIPTION: Mostra a ajuda
  echo ""
  echo " Uso.....: ${0##*/} <opcao>"
  echo ""
  First=1
  grep -i '(){$' ${BASH_SOURCE} | sed 's/(){//g'|while read function; do
     dsc=$(sed -n "/${function}(){$/ {n;p}" ${BASH_SOURCE} | sed "s/#*[DdEeSsCcRrIiPpTtIiOoNn]*://g")
     if [ ${First} -eq 1 ];then
        First=0
        printf " Opcoes..:%2s* \e[1;31m%s\e[m \e[1;37;1m-->\e[m\e[1;34m%s\e[m\n" " " "$function" "${dsc}"
     else
        printf "%12s* \e[1;31m%s\e[m \e[1;37;1m-->\e[m\e[1;34m%s\e[m\n" " " "$function" "${dsc}"
     fi
  done
  echo ""
  First=1
  grep -i '(){$' ${BASH_SOURCE} | sed 's/(){//g'|while read function; do
     dsc=$(sed -n "/${function}(){$/ {n;p}" ${BASH_SOURCE} | sed "s/#*[DdEeSsCcRrIiPpTtIiOoNn]*://g")
     if [ ${First} -eq 1 ];then
        First=0
        printf " Exemplos:%2s\e[1;37m${0##*/}\e[m \e[1;31m%s\e[m\n" " " "$function"
     else
        printf "%12s\e[1;37m${0##*/}\e[m \e[1;31m%s\e[m\n" " " "$function"
     fi
  done

}

teste(){
#DESCRIPTION: Funcao para testar uso de funcoes
   echo "OLA MUNDO"
}


#EOC
#-----------------------------------------------------------------------------#

