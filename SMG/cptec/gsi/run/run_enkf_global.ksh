#!/bin/ksh
#####################################################
# machine set up (users should change this part)
#####################################################
set -x

module load intel netcdf mpt

#####################################################
##case set up (users should change this part)
#####################################################
#
# GFSCASE = cases used for DTC test
#           T574, T254, T126, T62, enkf_glb_t254 
# ANAL_TIME= analysis time  (YYYYMMDDHH)
# WORK_ROOT= working directory, where GSI runs
# PREPBURF = path of PreBUFR conventional obs
# BK_ROOT  = path of background files
# OBS_ROOT = path of observations files
# FIX_ROOT = path of fix files
# GSI_EXE  = path and name of the gsi executable 
  ANAL_TIME=2014040506
  GFSCASE=enkf_glb_t254
  WORK_ROOT=comGSIv3.4-EnKFv1.0/run/${GFSCASE}
  DIAG_ROOT=comGSIv3.4-EnKFv1.0/run/${GFSCASE}-observer
  BK_ROOT=data/global_case/enkf_glb_t254/bk
  GSI_ROOT=comGSIv3.4-EnKFv1.0
  OBS_ROOT=data/global_case/enkf_glb_t254/obs
  FIX_ROOT=${GSI_ROOT}/fix
  ENKF_EXE=${GSI_ROOT}/src/main/enkf/global_enkf
  CRTM_ROOT=CRTM_REL-2.1.3
  ENKF_NAMELIST=${GSI_ROOT}/run/enkf_wrf_namelist.sh

# Note:  number of pe >= NMEM_ENKF
NMEM_ENKF=10
LEVS=64
NVARS=5

# Set the JCAP resolution which you want.
# All resolutions use LEVS=64
if [[ "$GFSCASE" = "T62" ]]; then
  JCAP=62
  JCAP_B=62
elif [[ "$GFSCASE" = "T126" ]]; then
  JCAP=126
  JCAP_B=126
elif [[ "$GFSCASE" = "enkf_glb_t254" ]]; then
  JCAP=254
  JCAP_B=254
elif [[ "$GFSCASE" = "T254" ]]; then
  JCAP=254
  JCAP_B=574
elif [[ "$GFSCASE" = "T574" ]]; then
  JCAP=574
  JCAP_B=1534
else
   echo "INVALID case = $GFSCASE"
   exit
fi

# Given the requested resolution, set dependent resolution parameters
if [[ "$JCAP" = "382" ]]; then
   LONA=768
   LATA=384
   DELTIM=180
   resol=1
elif [[ "$JCAP" = "574" ]]; then
   LONA=1152
   LATA=576
   DELTIM=1200
   resol=2
elif [[ "$JCAP" = "254" ]]; then
   LONA=512
   LATA=256
   DELTIM=1200
   resol=2
elif [[ "$JCAP" = "126" ]]; then
   LONA=256
   LATA=128
   DELTIM=1200
   resol=2
elif [[ "$JCAP" = "62" ]]; then
   LONA=192
   LATA=94
   DELTIM=1200
   resol=2
else
   echo "INVALID JCAP = $JCAP"
   exit
fi
NLAT=` expr $LATA + 2 `

ncp=/bin/cp


# Given the analysis date, compute the date from which the
# first guess comes.  Extract cycle and set prefix and suffix
# for guess and observation data files
# gdate=`$ndate -06 $ANAL_TIME`
gdate=$ANAL_TIME
PDYa=`echo $ANAL_TIME | cut -c1-8`
cyca=`echo $ANAL_TIME | cut -c9-10`
PDYg=`echo $gdate | cut -c1-8`
cycg=`echo $gdate | cut -c9-10`
prefix_tbc=gdas1.t${cycg}z


# Directories for test case
dirtbc=$BK_ROOT

# Fixed files
CONVINFO=${FIX_ROOT}/global_convinfo.txt
SATINFO=${FIX_ROOT}/global_satinfo.txt
SCANINFO=${FIX_ROOT}/global_scaninfo.txt
OZINFO=${FIX_ROOT}/global_ozinfo.txt
LOCINFO=${FIX_ROOT}/global_hybens_locinfo.l64.txt

# temporary directory
tmpdir=${WORK_ROOT}/${GFSCASE}_tmp

# Set up $tmpdir
rm -rf $tmpdir
mkdir -p $tmpdir
cd $tmpdir

NAM_ENKF=""

cat << EOF > enkf.nml
 &nam_enkf
  datestring="$ANAL_TIME",datapath="$tmpdir/",
  analpertwtnh=0.85,analpertwtsh=0.85,analpertwttr=0.85,
  covinflatemax=1.e2,covinflatemin=1,pseudo_rh=.true.,iassim_order=0,
  corrlengthnh=2000,corrlengthsh=2000,corrlengthtr=2000,
  lnsigcutoffnh=2.0,lnsigcutoffsh=2.0,lnsigcutofftr=2.0,
  lnsigcutoffpsnh=2.0,lnsigcutoffpssh=2.0,lnsigcutoffpstr=2.0,
  lnsigcutoffsatnh=2.0,lnsigcutoffsatsh=2.0,lnsigcutoffsattr=2.0,
  obtimelnh=1.e30,obtimelsh=1.e30,obtimeltr=1.e30,
  saterrfact=1.0,numiter=3,
  sprd_tol=1.e30,paoverpb_thresh=0.98,
  nlons=$LONA,nlats=$LATA,nlevs=$LEVS,nanals=$NMEM_ENKF,nvars=$NVARS,
  deterministic=.true.,sortinc=.true.,lupd_satbiasc=.false.,
  reducedgrid=.true.,readin_localization=.true.,
  $NAM_ENKF
 /
 &END
 &satobs_enkf
  sattypes_rad(1) = 'amsua_n15',     dsis(1) = 'amsua_n15',
  sattypes_rad(2) = 'amsua_n18',     dsis(2) = 'amsua_n18',
  sattypes_rad(3) = 'amsua_n19',     dsis(3) = 'amsua_n19',
  sattypes_rad(4) = 'amsub_n16',     dsis(4) = 'amsub_n16',
  sattypes_rad(5) = 'amsub_n17',     dsis(5) = 'amsub_n17',
  sattypes_rad(6) = 'amsua_aqua',    dsis(6) = 'amsua_aqua',
  sattypes_rad(7) = 'amsua_metop-a', dsis(7) = 'amsua_metop-a',
  sattypes_rad(8) = 'airs_aqua',     dsis(8) = 'airs281SUBSET_aqua',
  sattypes_rad(9) = 'hirs3_n17',     dsis(9) = 'hirs3_n17',
  sattypes_rad(10)= 'hirs4_n19',     dsis(10)= 'hirs4_n19',
  sattypes_rad(11)= 'hirs4_metop-a', dsis(11)= 'hirs4_metop-a',
  sattypes_rad(12)= 'mhs_n18',       dsis(12)= 'mhs_n18',
  sattypes_rad(13)= 'mhs_n19',       dsis(13)= 'mhs_n19',
  sattypes_rad(14)= 'mhs_metop-a',   dsis(14)= 'mhs_metop-a',
  sattypes_rad(15)= 'goes_img_g11',  dsis(15)= 'imgr_g11',
  sattypes_rad(16)= 'goes_img_g12',  dsis(16)= 'imgr_g12',
  sattypes_rad(17)= 'goes_img_g13',  dsis(17)= 'imgr_g13',
  sattypes_rad(18)= 'goes_img_g14',  dsis(18)= 'imgr_g14',
  sattypes_rad(19)= 'goes_img_g15',  dsis(19)= 'imgr_g15',
  sattypes_rad(20)= 'avhrr3_n16',    dsis(20)= 'avhrr3_n16',
  sattypes_rad(21)= 'avhrr3_n17',    dsis(21)= 'avhrr3_n17',
  sattypes_rad(22)= 'avhrr3_n18',    dsis(22)= 'avhrr3_n18',
  sattypes_rad(23)= 'amsre_aqua',    dsis(23)= 'amsre_aqua',
  sattypes_rad(24)= 'ssmis_f16',     dsis(24)= 'ssmis_f16',
  sattypes_rad(25)= 'ssmis_f17',     dsis(25)= 'ssmis_f17',
  sattypes_rad(26)= 'ssmis_f18',     dsis(26)= 'ssmis_f18',
  sattypes_rad(27)= 'ssmis_f19',     dsis(27)= 'ssmis_f19',
  sattypes_rad(28)= 'ssmis_f20',     dsis(28)= 'ssmis_f20',
  sattypes_rad(29)= 'sndrd1_g11',    dsis(29)= 'sndrD1_g11',
  sattypes_rad(30)= 'sndrd2_g11',    dsis(30)= 'sndrD2_g11',
  sattypes_rad(31)= 'sndrd3_g11',    dsis(31)= 'sndrD3_g11',
  sattypes_rad(32)= 'sndrd4_g11',    dsis(32)= 'sndrD4_g11',
  sattypes_rad(33)= 'sndrd1_g12',    dsis(33)= 'sndrD1_g12',
  sattypes_rad(34)= 'sndrd2_g12',    dsis(34)= 'sndrD2_g12',
  sattypes_rad(35)= 'sndrd3_g12',    dsis(35)= 'sndrD3_g12',
  sattypes_rad(36)= 'sndrd4_g12',    dsis(36)= 'sndrD4_g12',
  sattypes_rad(37)= 'sndrd1_g13',    dsis(37)= 'sndrD1_g13',
  sattypes_rad(38)= 'sndrd2_g13',    dsis(38)= 'sndrD2_g13',
  sattypes_rad(39)= 'sndrd3_g13',    dsis(39)= 'sndrD3_g13',
  sattypes_rad(40)= 'sndrd4_g13',    dsis(40)= 'sndrD4_g13',
  sattypes_rad(41)= 'sndrd1_g14',    dsis(41)= 'sndrD1_g14',
  sattypes_rad(42)= 'sndrd2_g14',    dsis(42)= 'sndrD2_g14',
  sattypes_rad(43)= 'sndrd3_g14',    dsis(43)= 'sndrD3_g14',
  sattypes_rad(44)= 'sndrd4_g14',    dsis(44)= 'sndrD4_g14',
  sattypes_rad(45)= 'sndrd1_g15',    dsis(45)= 'sndrD1_g15',
  sattypes_rad(46)= 'sndrd2_g15',    dsis(46)= 'sndrD2_g15',
  sattypes_rad(47)= 'sndrd3_g15',    dsis(47)= 'sndrD3_g15',
  sattypes_rad(48)= 'sndrd4_g15',    dsis(48)= 'sndrD4_g15',
  sattypes_rad(49)= 'iasi_metop-a',  dsis(49)= 'iasi616_metop-a',
  sattypes_rad(50)= 'seviri_m08',    dsis(50)= 'seviri_m08',
  sattypes_rad(51)= 'seviri_m09',    dsis(51)= 'seviri_m09',
  sattypes_rad(52)= 'seviri_m10',    dsis(52)= 'seviri_m10',
  sattypes_rad(53)= 'amsua_metop-b', dsis(53)= 'amsua_metop-b',
  sattypes_rad(54)= 'hirs4_metop-b', dsis(54)= 'hirs4_metop-b',
  sattypes_rad(55)= 'mhs_metop-b',   dsis(15)= 'mhs_metop-b',
  sattypes_rad(56)= 'iasi_metop-b',  dsis(56)= 'iasi616_metop-b',
  sattypes_rad(57)= 'avhrr3_metop-b',dsis(56)= 'avhrr3_metop-b',
  sattypes_rad(58)= 'atms_npp',      dsis(58)= 'atms_npp',
  sattypes_rad(59)= 'cris_npp',      dsis(59)= 'cris_npp',
  $SATOBS_ENKF
 /
 &END
 &ozobs_enkf
  sattypes_oz(1) = 'sbuv2_n16',
  sattypes_oz(2) = 'sbuv2_n17',
  sattypes_oz(3) = 'sbuv2_n18',
  sattypes_oz(4) = 'sbuv2_n19',
  sattypes_oz(5) = 'omi_aura',
  sattypes_oz(6) = 'gome_metop-a',
  sattypes_oz(7) = 'gome_metop-b',
  $OZOBS_ENKF
 /
 &END
EOF

cat enkf.nml

$ncp $ENKF_EXE        ./enkf.x

$ncp $CONVINFO        ./convinfo
$ncp $SATINFO         ./satinfo
$ncp $SCANINFO        ./scaninfo
$ncp $OZINFO          ./ozinfo
$ncp $LOCINFO         ./hybens_locinfo

$ncp $DIAG_ROOT/satbias_in ./satbias_in
$ncp $DIAG_ROOT/satbias_angle ./satbias_angle

# get mean
$ncp $BK_ROOT/sfg_${gdate}_fhr06_ensmean ./sfg_${ANAL_TIME}_fhr06_ensmean
list="conv amsua_metop-a amsua_n18 amsua_n15"
for type in $list; do
   $ncp $DIAG_ROOT/diag_${type}_ges.ensmean .
done

# get each member
imem=1
while [[ $imem -le $NMEM_ENKF ]]; do
   member="mem"`printf %03i $imem`
   $ncp $BK_ROOT/sfg_${gdate}_fhr06_${member} ./sfg_${ANAL_TIME}_fhr06_${member}
   list="conv amsua_metop-a amsua_n18 amsua_n15"
   for type in $list; do
      $ncp $DIAG_ROOT/diag_${type}_ges.${member} .
   done
   (( imem = $imem + 1 ))
done

eval "mpirun.lsf $tmpdir/enkf.x < enkf.nml > stdout"
rc=$?

exit
