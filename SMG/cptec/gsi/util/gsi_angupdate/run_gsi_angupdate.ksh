#!/bin/ksh
#####################################################
# machine set up (users should change this part)
#####################################################
#
# GSIPROC = processor number used for GSI analysis
#------------------------------------------------
  GSIPROC=1
  ARCH='LINUX_PBS'
# Supported configurations:
            # IBM_LSF,
            # LINUX, LINUX_LSF, LINUX_PBS,
            # DARWIN_PGI
#
#####################################################
# case set up (users should change this part)
#####################################################
#
# ANAL_TIME= analysis time  (YYYYMMDDHH)
# WORK_ROOT= working directory, where angupdate executable runs
# GSI_WORK_ROOT= GSI working directory, where GSI runs
# GSI_ANGUPDATE_EXE  = path and name of the gsi angupdate executable 
  ANAL_TIME=2011032212
  WORK_ROOT=/comGSI_v3.3/run/angupdate_${ANAL_TIME}
  GSI_WORK_ROOT=/comGSI_v3.3/run/arw_2011032212
  GSI_ANGUPDATE_EXE=/comGSI_v3.3/util/gsi_angupdate/gsi_angupdate.exe
#
#
#####################################################
# Users should NOT change script after this point
#####################################################
#
case $ARCH in
   'IBM_LSF')
      ###### IBM LSF (Load Sharing Facility)
      BYTE_ORDER=Big_Endian
      RUN_COMMAND="mpirun.lsf " ;;

   'IBM_LoadLevel')
      ###### IBM LoadLeve 
      BYTE_ORDER=Big_Endian
      RUN_COMMAND="poe " ;;

   'LINUX')
      BYTE_ORDER=Little_Endian
      if [ $GSIPROC = 1 ]; then
         #### Linux workstation - single processor
         RUN_COMMAND=""
      else
         ###### Linux workstation -  mpi run
        RUN_COMMAND="mpirun -np ${GSIPROC} -machinefile ~/mach "
      fi ;;

   'LINUX_LSF')
      ###### LINUX LSF (Load Sharing Facility)
      BYTE_ORDER=Little_Endian
      RUN_COMMAND="mpirun.lsf " ;;

   'LINUX_PBS')
      BYTE_ORDER=Little_Endian
      #### Linux cluster PBS (Portable Batch System)
      RUN_COMMAND="mpiexec_mpt -np ${GSIPROC} " ;;
##      RUN_COMMAND="mpirun -np ${GSIPROC} " ;;

   'DARWIN_PGI')
      ### Mac - mpi run
      BYTE_ORDER=Little_Endian
      if [ $GSIPROC = 1 ]; then
         #### Mac workstation - single processor
         RUN_COMMAND=""
      else
         ###### Mac workstation -  mpi run
         RUN_COMMAND="mpirun -np ${GSIPROC} -machinefile ~/mach "
      fi ;;

   * )
     print "error: $ARCH is not a supported platform configuration."
     exit 1 ;;
esac


##################################################################################
# Check GSI needed environment variables are defined and exist
#
 
# Make sure ANAL_TIME is defined and in the correct format
if [ ! "${ANAL_TIME}" ]; then
  echo "ERROR: \$ANAL_TIME is not defined!"
  exit 1
fi

# Make sure WORK_ROOT is defined and exists
if [ ! "${WORK_ROOT}" ]; then
  echo "ERROR: \$WORK_ROOT is not defined!"
  exit 1
fi

# Make sure the GSI executable exists
if [ ! -x "${GSI_ANGUPDATE_EXE}" ]; then
  echo "ERROR: ${GSI_ANGUPDATE_EXE} does not exist!"
  exit 1
fi

# Check to make sure the number of processors for running GSI was specified
if [ -z "${GSIPROC}" ]; then
  echo "ERROR: The variable $GSIPROC must be set to contain the number of processors to run GSI"
  exit 1
fi

#
##################################################################################
# Create the ram work directory and cd into it

workdir=${WORK_ROOT}
echo " Create working directory:" ${workdir}

if [ -d "${workdir}" ]; then
  rm -rf ${workdir}
fi
mkdir -p ${workdir}
cd ${workdir}

#
##################################################################################

# Loop over first and last outer loops to generate innovation
# diagnostic files for indicated observation types (groups)
#
ls -l ${GSI_WORK_ROOT}/diag_* > listpe
#  Collect diagnostic files for obs types (groups) below
   listall="hirs2_n14 msu_n14 sndr_g08 sndr_g11 sndr_g12 sndr_g13 sndr_g08_prep sndr_g11_prep sndr_g12_prep sndr_g13_prep sndrd1_g11 sndrd2_g11 sndrd3_g11 sndrd4_g11 sndrd1_g12 sndrd2_g12 sndrd3_g12 sndrd4_g12 sndrd1_g13 sndrd2_g13 sndrd3_g13 sndrd4_g13 sndrd1_g14 sndrd2_g14 sndrd3_g14 sndrd4_g14 sndrd1_g15 sndrd2_g15 sndrd3_g15 sndrd4_g15 hirs3_n15 hirs3_n16 hirs3_n17 amsua_n15 amsua_n16 amsua_n17 amsub_n15 amsub_n16 amsub_n17 hsb_aqua airs_aqua amsua_aqua imgr_g08 imgr_g11 imgr_g12 ssmi_f13 ssmi_f14 imgr_g14 imgr_g15 ssmi_f15 hirs4_n18 hirs4_metop-a amsua_n18 amsua_metop-a mhs_n18 mhs_metop-a amsre_low_aqua amsre_mid_aqua amsre_hig_aqua ssmis_las_f16 ssmis_uas_f16 ssmis_img_f16 ssmis_env_f16 ssmis_las_f17 ssmis_uas_f17 ssmis_img_f17 ssmis_env_f17 ssmis_las_f18 ssmis_uas_f18 ssmis_img_f18 ssmis_env_f18 ssmis_las_f19 ssmis_uas_f19 ssmis_img_f19 ssmis_env_f19 ssmis_las_f20 ssmis_uas_f20 ssmis_img_f20 ssmis_env_f20 iasi_metop-a hirs4_n19 amsua_n19 mhs_n19 seviri_m08 seviri_m09 seviri_m10 cris_npp atms_npp hirs4_metop-b amsua_metop-b mhs_metop-b iasi_metop-b"
   for type in $listall; do
      count=`grep diag_${type}_ges.${ANAL_TIME} listpe | wc -l`
      if [[ $count -gt 0 ]]; then
         ln -fs ${GSI_WORK_ROOT}/diag_${type}_ges.${ANAL_TIME} ./diag_${type}.${ANAL_TIME}
      fi
   done


##################################################################################

echo " Copy GSI executable, background file, and link observation bufr to working directory"

# Save a copy of the GSI executable in the workdir
cp ${GSI_ANGUPDATE_EXE} gsi_angupdate.exe

## if [[ -s satbias_ang.in ]]; then
   cp ${GSI_WORK_ROOT}/satbias_angle ./satbias_ang.in
## fi

iy=$(echo $ANAL_TIME |cut -c1-4)
im=$(echo $ANAL_TIME |cut -c5-6)
id=$(echo $ANAL_TIME |cut -c7-8)
ih=$(echo $ANAL_TIME |cut -c9-10)

#
##################################################################################
# Set some parameters for use by the GSI executable and to build the namelist
echo " Build the namelist "

# Build the GSI namelist on-the-fly
cat << EOF > gsi_angupdate.namelist
 &setup
  jpch=2680,nstep=90,nsize=20,wgtang=0.008333333,wgtlap=0.0,
  iuseqc=1,dtmax=1.0,
  iyy1=${iy},imm1=${im},idd1=${id},ihh1=${ih},
  iyy2=${iy},imm2=${im},idd2=${id},ihh2=${ih},
  dth=01,ndat=50
 /
 &obs_input
  dtype(01)='hirs3',     dplat(01)='n17',       dsis(01)='hirs3_n17',
  dtype(02)='hirs4',     dplat(02)='metop-a',   dsis(02)='hirs4_metop-a',
  dtype(03)='goes_img',  dplat(03)='g11',       dsis(03)='imgr_g11',
  dtype(04)='goes_img',  dplat(04)='g12',       dsis(04)='imgr_g12',
  dtype(05)='airs',      dplat(05)='aqua',      dsis(05)='airs281SUBSET_aqua',
  dtype(06)='amsua',     dplat(06)='n15',       dsis(06)='amsua_n15',
  dtype(07)='amsua',     dplat(07)='n18',       dsis(07)='amsua_n18',
  dtype(08)='amsua',     dplat(08)='metop-a',   dsis(08)='amsua_metop-a',
  dtype(09)='amsua',     dplat(09)='aqua',      dsis(09)='amsua_aqua',
  dtype(10)='mhs',       dplat(10)='n18',       dsis(10)='mhs_n18',
  dtype(11)='mhs',       dplat(11)='metop-a',   dsis(11)='mhs_metop-a',
  dtype(12)='ssmi',      dplat(12)='f15',       dsis(12)='ssmi_f15',
  dtype(13)='amsre_low', dplat(13)='aqua',      dsis(13)='amsre_aqua',
  dtype(14)='amsre_mid', dplat(14)='aqua',      dsis(14)='amsre_aqua',
  dtype(15)='amsre_hig', dplat(15)='aqua',      dsis(15)='amsre_aqua',
  dtype(16)='ssmis_las', dplat(16)='f16',       dsis(16)='ssmis_f16',
  dtype(17)='ssmis_uas', dplat(17)='f16',       dsis(17)='ssmis_f16',
  dtype(18)='ssmis_img', dplat(18)='f16',       dsis(18)='ssmis_f16',
  dtype(19)='ssmis_env', dplat(19)='f16',       dsis(19)='ssmis_f16',
  dtype(20)='sndrd1',    dplat(20)='g12',       dsis(20)='sndrD1_g12',
  dtype(21)='sndrd2',    dplat(21)='g12',       dsis(21)='sndrD2_g12',
  dtype(22)='sndrd3',    dplat(22)='g12',       dsis(22)='sndrD3_g12',
  dtype(23)='sndrd4',    dplat(23)='g12',       dsis(23)='sndrD4_g12',
  dtype(24)='sndrd1',    dplat(24)='g11',       dsis(24)='sndrD1_g11',
  dtype(25)='sndrd2',    dplat(25)='g11',       dsis(25)='sndrD2_g11',
  dtype(26)='sndrd3',    dplat(26)='g11',       dsis(26)='sndrD3_g11',
  dtype(27)='sndrd4',    dplat(27)='g11',       dsis(27)='sndrD4_g11',
  dtype(28)='sndrd1',    dplat(28)='g13',       dsis(28)='sndrD1_g13',
  dtype(29)='sndrd2',    dplat(29)='g13',       dsis(29)='sndrD2_g13',
  dtype(30)='sndrd3',    dplat(30)='g13',       dsis(30)='sndrD3_g13',
  dtype(31)='sndrd4',    dplat(31)='g13',       dsis(31)='sndrD4_g13',
  dtype(32)='iasi',      dplat(32)='metop-a',   dsis(32)='iasi616_metop-a',
  dtype(33)='hirs4',     dplat(33)='n19',       dsis(33)='hirs4_n19',
  dtype(34)='amsua',     dplat(34)='n19',       dsis(34)='amsua_n19',
  dtype(35)='mhs',       dplat(35)='n19',       dsis(35)='mhs_n19',
  dtype(36)='amsub',     dplat(36)='n17',       dsis(36)='amsub_n17',
  dtype(37)='hirs4',     dplat(37)='metop-b',   dsis(37)='hirs4_metop-b',
  dtype(38)='amsua',     dplat(38)='metop-b',   dsis(38)='amsua_metop-b',
  dtype(39)='mhs',       dplat(39)='metop-b',   dsis(39)='mhs_metop-b',
  dtype(40)='iasi',      dplat(40)='metop-b',   dsis(40)='iasi616_metop-b',
  dtype(41)='atms',      dplat(41)='npp',       dsis(41)='atms_npp',
  dtype(42)='cris',      dplat(42)='npp',       dsis(42)='cris_npp',
  dtype(43)='sndrd1',    dplat(43)='g14',       dsis(43)='sndrD1_g14',
  dtype(44)='sndrd2',    dplat(44)='g14',       dsis(44)='sndrD2_g14',
  dtype(45)='sndrd3',    dplat(45)='g14',       dsis(45)='sndrD3_g14',
  dtype(46)='sndrd4',    dplat(46)='g14',       dsis(46)='sndrD4_g14',
  dtype(47)='sndrd1',    dplat(47)='g15',       dsis(47)='sndrD1_g15',
  dtype(48)='sndrd2',    dplat(48)='g15',       dsis(48)='sndrD2_g15',
  dtype(49)='sndrd3',    dplat(49)='g15',       dsis(49)='sndrD3_g15',
  dtype(50)='sndrd4',    dplat(50)='g15',       dsis(50)='sndrD4_g15',
 /
EOF

#
###################################################
#  run  GSI
###################################################

case $ARCH in
   'IBM_LSF'|'IBM_LoadLevel')
      ${RUN_COMMAND} ./gsi_angupdate.exe < gsi_angupdate.namelist > stdout 2>&1  ;;

   * )
      ${RUN_COMMAND} ./gsi_angupdate.exe < gsi_angupdate.namelist > stdout 2>&1  ;;
esac

#
##################################################################
#

exit 0
