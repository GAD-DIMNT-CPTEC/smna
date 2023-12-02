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
eval export $1=$2
}

# Variaveis principais
vars_export(){
#DESCRIPTION: Exporta as variáveis de ambiente do SMG
FilePaths=$(dirname ${BASH_SOURCE})/mach/${hpc_name}_paths.conf
while read line;do

   assign $line

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
  if [ ${HOSTNAME:0:1} = 'e' ] || [ ${hpc_name} = "egeon" ];then
     cp -pf ${home_model_bam}/datain/* ${subt_model_bam}/datain/

     cp -pf ${public_bam}/MODEL/datain/AeroVar.Tab ${subt_model_bam}/datain/
     cp -pf ${public_bam}/MODEL/datain/ETAMPNEW_DATA ${subt_model_bam}/datain/
     cp -pf ${public_bam}/MODEL/datain/F_nwvl200_mu20_lam50_res64_t298_c080428.bin ${subt_model_bam}/datain/
     cp -pf ${public_bam}/MODEL/datain/iceoptics_c080917.bin ${subt_model_bam}/datain/
     cp -pf ${public_bam}/MODEL/datain/ocnalbtab24bnd.bin ${subt_model_bam}/datain/

     #-----------------------------------------------------------------------------------------
     # pre/datain

     echo " home_pre_bam :" ${home_pre_bam} " --->>>  subt_pre_bam :"${subt_pre_bam}

     cp -pf ${home_pre_bam}/datain/* ${subt_pre_bam}/datain/

     #-----------------------------------------------------------------------------------------
     # pre/dataout

     cp -pf ${home_pre_bam}/dataout/* ${subt_pre_bam}/dataout/

     cp -pf ${public_bam}/PRE/dataout/WaterNavy.dat ${subt_pre_bam}/dataout/
     cp -pf ${public_bam}/PRE/dataout/TopoNavy.dat ${subt_pre_bam}/dataout/
     cp -pf ${public_bam}/PRE/dataout/HPRIME.dat ${subt_pre_bam}/dataout/

     #-----------------------------------------------------------------------------------------
     # pre/databcs

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
  elif [ ${HOSTNAME:0:1} = 'c' ];then

     ln -s /cray_home/joao_gerd/BAMFIX/model/datain/* ${subt_model_bam}/datain/
     ln -s /cray_home/joao_gerd/BAMFIX/pre/datain/* ${subt_pre_bam}/datain/
     ln -s /cray_home/joao_gerd/BAMFIX/pre/dataout/* ${subt_pre_bam}/dataout/
     ln -s /cray_home/joao_gerd/BAMFIX/pre/databcs/* ${subt_pre_bam}/databcs/
     ln -s cray_home/joao_gerd/BAMFIX/pre/dataco2/* ${subt_pre_bam}/dataco2/
     ln -s cray_home/joao_gerd/BAMFIX/pre/datasst/* ${subt_pre_bam}/datasst/
     ln -s cray_home/joao_gerd/BAMFIX/pre/dataTop/* ${subt_pre_bam}/dataTop/


#      echo -e '\033[31;1m-------------------------------------------------------------------\033[m'
#      echo -e '\033[31;1mATENÇÃO:\033[m\033[33;1m É necessário verificar com o administrador\033[m'
#      echo -e '\033[33;1m         onde estão os arquivos fixos para o modelo BAM\033[m'
#      echo -e ''
#      echo -e '\033[33;1m Entre em contato com \033[m\033[34;1meduardo.khamis@inpe.br\033[m ou \033[34;1mjoao.gerd@inpe.br\033[m'
#      echo -e '\033[31;1m-------------------------------------------------------------------\033[m'
#
#     cp -pf ${home_model_bam_datain}/* ${subt_model_bam}/datain/
#
#     scp ${USER//_/.}@tupa:${public_bam}/MODEL/datain/AeroVar.Tab ${subt_model_bam}/datain/
#     scp ${USER//_/.}@tupa:${public_bam}/MODEL/datain/ETAMPNEW_DATA ${subt_model_bam}/datain/
#     scp ${USER//_/.}@tupa:${public_bam}/MODEL/datain/F_nwvl200_mu20_lam50_res64_t298_c080428.bin ${subt_model_bam}/datain/
#     scp ${USER//_/.}@tupa:${public_bam}/MODEL/datain/iceoptics_c080917.bin ${subt_model_bam}/datain/
#     scp ${USER//_/.}@tupa:${public_bam}/MODEL/datain/ocnalbtab24bnd.bin ${subt_model_bam}/datain/
#
#     #-----------------------------------------------------------------------------------------
#     # pre/datain
#
#     cp -pf ${home_pre_bam_datain}/* ${subt_pre_bam}/datain/
#
#     #-----------------------------------------------------------------------------------------
#     # pre/dataout
#
#     cp -pf ${home_pre_bam_dataout}/* ${subt_pre_bam}/dataout/
#
#     scp ${USER//_/.}@tupa:${public_bam}/PRE/dataout/WaterNavy.dat ${subt_pre_bam}/dataout/
#     scp ${USER//_/.}@tupa:${public_bam}/PRE/dataout/TopoNavy.dat ${subt_pre_bam}/dataout/
#     scp ${USER//_/.}@tupa:${public_bam}/PRE/dataout/HPRIME.dat ${subt_pre_bam}/dataout/
#
#     #-----------------------------------------------------------------------------------------
#     # pre/databcs
#
#     cp -pf ${home_pre_bam_databcs}/* ${subt_pre_bam}/databcs/
#
#     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/sib2soilms.form ${subt_pre_bam}/databcs/
#     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/FluxCO2.bin ${subt_pre_bam}/databcs/
#     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/FluxCO2.ctl ${subt_pre_bam}/databcs/
#     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/claymsk.form ${subt_pre_bam}/databcs/
#     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/clmt.form ${subt_pre_bam}/databcs/
#     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/deltat.form ${subt_pre_bam}/databcs/
#     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/ersst.bin ${subt_pre_bam}/databcs/
#     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/ibismsk.form ${subt_pre_bam}/databcs/
#     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/ndviclm.form ${subt_pre_bam}/databcs/
#     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/sandmsk.form ${subt_pre_bam}/databcs/
#     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/sib2msk.form ${subt_pre_bam}/databcs/
#     scp ${USER//_/.}@tupa:${public_bam}/PRE/databcs/soiltext.form ${subt_pre_bam}/databcs/
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
  if test ! -e ${subt_bam}; then mkdir -p ${subt_bam}; fi
  if test ! -e ${subt_gsi}; then mkdir -p ${subt_gsi}; fi

  if test ! -e ${subt_obs_dataout}; then mkdir -p ${subt_obs_dataout}; fi

# Cria datarun em diante (bam, gsi, obs, interface)
  if test ! -s ${subt_run}; then mkdir -p ${subt_run}; fi
  if test ! -s ${subt_run_gsi}; then mkdir -p ${subt_run_gsi}; fi
  if test ! -s ${subt_obs_run}; then mkdir -p ${subt_obs_run}; fi

# Cria pastas do MCGA

  # model
  for dir in exec datain dataout; do
    if [ ! -e ${subt_model_bam}/${dir} ];then
       mkdir -p ${subt_model_bam}/${dir}
    fi
  done

  # pre
  for dir in exec datain dataout datasst databcs dataco2 dataTop; do
    if [ ! -e ${subt_pre_bam}/${dir} ];then
       mkdir -p ${subt_pre_bam}/${dir}
    fi
  done

  # pos
  for dir in exec datain dataout; do
    if [ ! -e ${subt_pos_bam}/${dir} ];then
       mkdir -p ${subt_pos_bam}/${dir}
    fi
  done

  # grh
  for dir in exec datain dataout; do
    if [ ! -e ${subt_grh_bam}/${dir} ];then
       mkdir -p ${subt_grh_bam}/${dir}
    fi
  done


  if test ! -e ${subt_gsi}; then mkdir -p ${subt_gsi}; fi


# Altera scripts e Makefiles do SMG
  sed -i "/# Carregando as variaveis do sistema/{n;d}" ${run_smg}/run_cycle.sh ${scripts_smg}/runGSI ${scripts_smg}/run_model.sh ${scripts_smg}/run_obsmake.sh ${scripts_smg}/run_blsdas.sh
  sed -i "/# Carregando as variaveis do sistema/a\source "${home_smg}"\/config_smg.ksh vars_export"  ${run_smg}/run_cycle.sh ${scripts_smg}/runGSI ${scripts_smg}/run_model.sh ${scripts_smg}/run_obsmake.sh ${scripts_smg}/run_blsdas.sh

  sed -i "/# Caminho onde devera fica o arquivo executavel do Model/{n;d}" ${home_model_bam}/source/Makefile
  sed -i "/# Caminho onde devera fica o arquivo executavel do Model/a\PATH2="${home_model_bam}/exec"" ${home_model_bam}/source/Makefile

  sed -i "/# Caminho onde deverao ficar os arquivos executaveis do Pos/{n;d}" ${home_pos_bam}/source/Makefile
  sed -i "/# Caminho onde deverao ficar os arquivos executaveis do Pos/a\PATH2="${home_pos_bam}/exec"" ${home_pos_bam}/source/Makefile

  sed -i "/DirInPut='\//,1d" ${home_run_bam}/PostGridHistory.nml
  sed -i "/!DirInPut=subt_model_bam_dataout/a\DirInPut=\'"${subt_model_bam}/dataout"\/TQ0042L028\'" ${home_run_bam}/PostGridHistory.nml

  sed -i "/DirOutPut='\//,1d" ${home_run_bam}/PostGridHistory.nml
  sed -i "/!DirOutPut=subt_grh_bam_dataout/a\DirOutPut=\'"${subt_grh_bam}/dataout"\/TQ0042L028\'" ${home_run_bam}/PostGridHistory.nml

  sed -i "/DirMain='\//,1d" ${home_run_bam}/PostGridHistory.nml
  sed -i "/!DirMain=subt_bam/a\DirMain=\'"${subt_bam}"\'" ${home_run_bam}/PostGridHistory.nml

  sed -i "/export HOMEBASE=\//,1d" ${home_run_bam}/EnvironmentalVariables
  sed -i "/# Caminho do BAM no HOME/a\export HOMEBASE="${home_bam}"" ${home_run_bam}/EnvironmentalVariables

  sed -i "/export SUBTBASE=\//,1d" ${home_run_bam}/EnvironmentalVariables
  sed -i "/# Caminho do BAM no scratchin (SUBMIT_HOME)/a\export SUBTBASE="${subt_bam}"" ${home_run_bam}/EnvironmentalVariables


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
   echo "Compiler : " $compiler
# Criando diretório onde estarão os executáveis
   if [ ! -e ${home_cptec}/bin ];then
       mkdir -p ${home_cptec}/bin
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

## Liga ou desliga compilação de componentes
# 

# 0: Não compila;
# 1: Compila.

compgsi=0
compang=1
compbam=0

#################




   echo ""
   echo "%%%"
   echo " Compilando utilitarios do SMG:  "
   echo "%%%"

#### Comando `date ...` do Linux substitui o inctime 
## Compilando utilitario inctime
#   echo ""
#   echo "%%% Inctime:"
#   cd ${util_inctime}
#   make clean
#   ./autogen.sh
#   ./configure --prefix=${home_cptec}
#   if [ $? -ne 0 ];then
#      echo ''
#      echo -e '\033[31;1mErro ao compilar inctime[m'
#      echo -e '\033[32;1mProvavelmente não esta definido FC (fortran compiler)\033[m'
#      echo ''
#   fi
#   if [ ${hpc_name} = "egeon" ];then
#      # remove extra "-module ../include"  that messed with linking
#      sed -i 's/\-module\ \.\.\/include/ /g' src/Makefile
#   fi

#   make
#   make install

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

if [ ${compgsi} -eq 1 ]; then

# Configurando o GSI

   if [ ${hpc_name} = 'XC50' ];then
        source ./env.sh xc50 ${compiler}
   elif [ ${hpc_name} = 'egeon' ];then
        source ./env.sh egeon ${compiler}
   fi

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
   pwd
   ./compile.sh 2>&1 | tee ${home_gsi}/compile.log

   if [ -e ${home_gsi_src}/gsi.x ]; then
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

fi

if [ ${compang} -eq 1 ]; then

# Compilando do utilitario para geracão de bias relacionado com o ângulo

   if [ ${hpc_name} = 'XC50' ];then
        source ./env.sh xc50 ${compiler}
   elif [ ${hpc_name} = 'egeon' ];then
        source ./env.sh egeon ${compiler}
   fi
   echo ""
   echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
   echo "   Compilando utilitário de gereção de correção de bias do GSI:     "
   echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

   ## OLD GSI version : cd ${home_gsi}/util/gsi_angupdate
   cd ${home_gsi}/util/global_angupdate
   echo ln -sf Makefile.conf.${hpc_name}-${compiler} Makefile.conf
   ln -sf Makefile.conf.${hpc_name}-${compiler} Makefile.conf
   make -f Makefile clean
   make -f Makefile
   
   ## OLD GSI version :if [ -e ${home_gsi}/util/gsi_angupdate/gsi_angupdate.exe  ]; then
   if [ -e ${home_gsi}/util/global_angupdate/global_angupdate ]; then
     echo ""
     echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     echo "!                                                                                           !"
     echo "!      Compilacao do utilitário de gereção de correção de bias do ângulo do GSI com Sucesso !"
     echo "!                                                                                           !"
     echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     echo ""

     cp -pfvr ${home_gsi}/util/global_angupdate/global_angupdate ${home_cptec}/bin/global_angupdate

   else
     echo ""
     echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     echo "!                                                                      !"
     echo "! Compilacao do utilitário de gereção de correção de bias  falhou !!!! !"
     echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     echo ""
     exit -9
   fi

fi
#############################################################################
# Compilacao bam
#############################################################################

if [ ${compbam} -eq 1 ]; then
   echo ""
   echo "%%%%%%%%%%%%%%%%%%%%"
   echo " Compilando o BAM:  "
   echo "%%%%%%%%%%%%%%%%%%%%"
   echo ""
 
   export mkname=${compiler}_${SUB}
   if [ ${hpc_name} = "egeon" ];then
      . /mnt/beegfs/jose.aravequia/libs_impi/env_libs_impi.sh
   fi
#########################################
# Compilacao pre/pos
#########################################

#  . /opt/modules/default/etc/modules.sh
#  module swap pgi pgi/11.10.0

   echo "+++++++++++++++++++++++++"
   echo "   Compilando o pre:     "
   echo "+++++++++++++++++++++++++"

  ## As instalacoes NETCDF separam C da Fortran , mas no XC50 estão no mesmo diretório
   if [ ${SUB} = "cray" ];then 
      export NETCDF_FORTRAN_DIR = ${NETCDF_DIR}
   fi

   cd ${home_pre_bam}/build
   make clean ${mkname}
   make ${mkname}
   make install

   echo "+++++++++++++++++++++++++"
   echo "   Compilando o pos:     "
   echo "+++++++++++++++++++++++++"

   cd ${home_pos_bam}/source
   make clean 

   make clean ${mkname}


#########################################
# Compilacao model
#########################################

   echo "+++++++++++++++++++++++++"
   echo "   Compilando o model:   "
   echo "+++++++++++++++++++++++++"

#  . /opt/modules/default/etc/modules.sh
#  module swap PrgEnv-pgi PrgEnv-cray

   cd ${home_model_bam}/source
#   make clean cray_cray32
   make clean $mkname


   echo -e "\033[34;1m ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \033[m"
   echo -e "\033[34;1m > Compilacao SMG completa.# Verifique possiveis erros no arquivo de log, caso o tenha criado. \033[m"
   echo -e "\033[34;1m ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \033[m"

fi   
#########################################
# Compilacao do blsdas
#########################################
#
#   echo "+++++++++++++++++++++++++"
#   echo "   Compilando o BLSDAS:   "
#   echo "+++++++++++++++++++++++++"
#
#   . /opt/modules/default/etc/modules.sh
#   module swap PrgEnv-pgi PrgEnv-gnu
#
#   export ARCH=linux_gnu
#
#   cd ${home_blsdas}
#   make clean
#   make
#
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


   cp -pfr ${home_gsi_src}/gsi.x ${home_cptec}/bin

#
# exec BAM
#
   cp -pfr ${home_pre_bam}/exec/* ${subt_pre_bam}/exec
   cp -pfr ${home_model_bam}/exec/* ${subt_model_bam}/exec
   cp -pfr ${home_pos_bam}/exec/* ${subt_pos_bam}/exec



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
  cp -pvfrL ${anlfile}  ${subt_pre_bam}/datain/
  cp -pvfrL ${public_bam}/PRE/datain/${year[${answer}]}/sst/* ${subt_pre_bam}/datain/
  cp -pvfrL ${public_bam}/PRE/datain/${year[${answer}]}/sno/* ${subt_pre_bam}/datain/
  cp -pvfrL ${public_bam}/PRE/datain/${year[${answer}]}/smc/*.vfm ${subt_pre_bam}/datain/

  #
  # independentes
  #

  cp -pvfr ${public_bam}/PRE/dataout/* ${subt_pre_bam}/dataout/
  cp -pvfr ${public_bam}/PRE/databcs/* ${subt_pre_bam}/databcs/
  cp -pvfr ${public_bam}/PRE/datasst/* ${subt_pre_bam}/datasst/

  ln -s ${home_pos_bam}/datain/* ${subt_pos_bam}/datain/

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
