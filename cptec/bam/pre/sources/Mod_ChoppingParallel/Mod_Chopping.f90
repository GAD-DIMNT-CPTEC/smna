!# @info
!# ---
!# INPE/CPTEC, DIDMD, Modelling and Development Division
!# ---
!# </br>
!#
!# **Module**: Mod_Chopping </br></br>
!#
!# **Brief**: Module used for Running Chopping in Parallel. Uses stand alone
!# execution mode in an new MPI comminicator (check README), through use of
!# Mod_Parallelism_Group_Chopping module. </br>
!#
!# The Chopping_parallel.f90 algorithm reads a data set (topography, surface 
!# pressure Ln, virtual temperature, divergence, vorticity, specific humidity,
!# ozone and tracers) from the gdas1.ThhZ.Sanl.YYYYMMDDHH file. This file is 
!# called atmosphere analysis, and is generated in NCEP through the data 
!# assimilation algorithm. For vertical interpolation the algorithm reads the
!# DeltaSigma.LZZZ file which contains the variation of the sigma level which
!# one wishes to interpolate. Another file that is used in the process is the 
!# Topography.TQXXXX spectral topography, which is intended to replace the low
!# resolution topography that is present in the analysis file. </br>
!#
!# Task 6061 - using subroutines to return erros to Mod_Pre; precision 
!# adjustments: fixes TopoWaterGT30 reprodutibility. </br></br>
!#
!# **Files in:**
!#
!# &bull; pre/datain/gdas1.ThhZ.SAnl.YYYYMMDDHH (Ex.: pre/datain/gdas1.T00Z.SAnl.2015043000) </br>
!# &bull; pre/datain/DeltaSigma.LZZZ (Ex.: pre/datain/DeltaSigma.L064) </br>
!# &bull; pre/dataout/Topography.TQXXXX (Ex.: pre/dataout/Topography.TQ0299) </br></br>
!#
!# **Files out:**
!#
!# &bull; model/datain/GANLNMCYYYYMMDDHHS.unf.TQXXXXLZZZ (Ex.: model/datain/GANLNMC2015043000S.unf.TQ0299L064) </br>
!! &bull; model/datain/OZONNMCYYYYMMDDHHS.unf.GZZZZZLZZZ (Ex.: model/datain/OZONNMC2015043000S.unf.G00450L064) </br>
!! &bull; model/datain/TRACNMCYYYYMMDDHHS.unf.GZZZZZLZZZ </br>
!! &bull; pre/dataout/GANLNMCYYYYMMDDHHS.unf.TQXXXXLZZZ.GrADS </br>
!! &bull; pre/dataout/GANLNMCYYYYMMDDHHS.unf.TQXXXXLZZZ.GrADS.ctl </br></br>
!# 
!# &bull; model/datain/GANLSMTYYYYMMDDHHS.unf.TQXXXXLZZZ </br>
!# &bull; model/datain/OZONSMTYYYYMMDDHHS.unf.GZZZZZLZZZ </br>
!# &bull; model/datain/TRACSMTYYYYMMDDHHS.unf.GZZZZZLZZZ </br>
!# &bull; pre/dataout/GANLSMTYYYYMMDDHHS.unf.TQYYYYLZZZ.GrADS </br>
!# &bull; pre/dataout/GANLSMTYYYYMMDDHHS.unf.TQYYYYLZZZ.GrADS.ctl </br></br>
!# 
!# When in the names of these files there is the prefix SMT instead of NMC, is because they suffered dampening (smoothing).
!# 
!# **Author**: Paulo Bonatti </br>
!#
!# **Version**: 2.2.0 </br></br>
!# @endinfo
!#
!# @changes
!# <ul type="disc">
!#  <li>13-11-2004 - Paulo Bonatti  - version: 1.0.0 </li>
!#  <li>01-08-2017 - Simone Tomita  - version: 1.1.1 </li>
!#  <li>26-04-2019 - Denis Eiras    - version: 2.0.0 </li>
!#  <li>12-10-2019 - Eduardo Khamis - version: 2.1.0 - changing for operational Chopping </li>
!#  <li>12-11-2019 - Denis Eiras    - version: 2.2.0 </li>
!# </ul>
!# @endchanges
!#
!# @bug
!# <ul type="disc">
!#  <li>Bug 6241 - incorrect pressure levels and variable number in GrADSOut.ctl e OZONEGrADSOut.ctl </li>
!#  <li>Bug 6239 - Grads ozon file generated with errors </li>
!# </ul>
!# @endbug
!#
!# @todo
!# <ul type="disc">
!#  <li>None items at this time </li>
!# </ul>
!# @endtodo
!#
!# @documentation
!#
!# For theoretical information, please visit the following link: </br> <http://urlib.net/8JMKD3MGP3W34R/3SME6J2> </br>
!# **&#9993;**<mailto:atende.cptec@inpe.br> </br></br>
!# @enddocumentation
!#
!# @warning
!# Copyright Under GLP-3.0
!# &copy; https://opensource.org/licenses/GPL-3.0
!# @endwarning
!#
!#---


module Mod_Chopping

  use nemsio_gfs
  use nemsio_module &
    , c_nemsio_init => nemsio_init &
    , c_nemsio_open => nemsio_open &
    , c_nemsio_getfilehead => nemsio_getfilehead &
    , c_nemsio_close => nemsio_close &
    , c_nemsio_finalize => nemsio_finalize &
    , c_nemsio_gfile => nemsio_gfile

  use Mod_Parallelism_Group_Chopping, only : &
    initializeParallelismVariables &
    , mpiCommGroup &
    , myId &                  
    ! MPI process rank
    , mpiMasterProc &         
    ! MPI master process rank of this group
    , isMasterProc &          
    !  MPI process rank is master in this group
    , getMyIdString &
    , fatalError &
    , msgOutMaster &
    , msgDump

  use Mod_Parallelism_Fourier, only : &
    myId_four

  !  USE Dumpgraph, ONLY:   &
  !       dumpgra, Writectl

  use Mod_InitChoppingParallel, only : &
    Initall

  use Mod_InputParameters, only : InitParameters, &
    KmaxInpp, KmaxOutp, &
    ImaxOut, JmaxOut, TruncInp, TruncOut, &
    Mnwv2Inp, Mnwv2Out, Mnwv3Inp, Mnwv3Out, &
    NTracers, Kdim, ICaseDec, IcaseRec, &
    IDVCInp, IDSLInp, ForecastDay, TimeOfDay, cTv, MonChar, &
    DataCPT, DataInp, DataOut, DataOup, DataTop, DataTopG, DataSig, DataSigInp, &
    OzonInp, TracInp, OzonOut, TracOut, &
    RoCp, RoCp1, RoCpR, p0mb, p0Pa, iMaxInp, jmaxInp &
    , ChoppingNameListData, fillDataOutFileName

  use Mod_InputArrays, only : GetArrays, ClsArrays, GetSpHuTracers, ClsSpHuTracers, &
    DateInitial, DateCurrent, DelSInp, SigIInp, SigLInp, SigIOut, &
    SigLOut, qWorkInp, qWorkInp3D, qWorkOut, DelSigmaInp, SigInterInp, &
    SigLayerInp, DelSigmaOut, SigInterOut, SigLayerOut, &
    qTopoInp, qLnPsInp, qTopoOut, qLnPsOut, &
    qDivgInp, qVortInp, qTvirInp, qDivgOut, qVortOut, &
    qTvirOut, qUvelInp, qVvelInp, qSpHuInp, qSpHuOut, &
    qUvelOut, qVvelOut, qWorkOut1, qworkprout, &
    gworkprout, qtorto, gpresaux, &
    gWorkOut, gTopoInp, gTopoOut, qTopoOutSpec, gTopoOutGaus, gTopoOutGaus8, gTopoDel, gLnPsInp, &
    gPsfcInp, gLnPsOut, gPsfcOut, gUvelInp, gVvelInp, &
    gTvirInp, gDivgInp, gVortInp, gPresInp, gUvelOut, &
    gPresInpp, gWorkprInp, &
    gVvelOut, gVvelInp, gTvirOut, gTvirInp, gPresOut, gSpHuInp, gSpHuOut, &
    gSpHuInp, qWorkInOut, qWorkInOut1, qSpHuInp, qDivgInp, qVortInp, qTvirInp, &
    qVvelInp

  use Mod_VerticalInterpolation, only : VertSigmaInter

  use Mod_Communications, only : Collect_Grid_Full, Collect_Spec, Clear_Communications

  use Mod_Transform, only : DepositSpecToGrid, CreateSpecToGrid, DoSpecToGrid, &
    DestroySpecToGrid, CreateGridToSpec, DoGridToSpec, &
    DepositGridToSpec, DestroyGridToSpec, DepositGridToSpec_PK, Clear_Transform

  use Mod_SpecDynamics, only : dztouv, uvtodz, Clear_SpecDynamics

  use Mod_Utils, only : NewSigma, SigmaInp, NewPs, coslat, rcl, lati, IJtoIBJB, CyclicLinear_inter, Clear_Utils

  use Mod_Sizes, only : jbmax, Ibmaxperjb, mymmax, msinproc, mnmap, mymnmap, &
    mnmap_out, myfirstlev, mylastlev, ibmax, mynMap, &
    ThreadDecomp, ReshapeVerticalGroups, havesurf, kmaxloc, &
    mnmax, mymnextmax, mymnmax, kmaxloc_out, kMaxloc_In, kmax, lm2m, &
    gridmap, ibperij, jbperij, mmax, &
    iperijb, jperijb, myfirstlon, mylastlon, mnmax_out, &
    imax, jmax, Ibmaxperjb, myfirstlat, myfirstlat_diag, mylastlat_diag, mylastlat, Clear_Sizes

  use Mod_FileManager, only: &
    getFileUnit,  &
    openFile, &
    fileExists
  
  use Mod_Messages, only : &
    msgInLineFormatOut &
    , msgNewLine
  
  implicit none

  public :: generateChopping
  public :: getNameChopping
  public :: initChopping
  public :: shouldRunChopping

  private
  include 'precision.h'
  include 'pre.h'
  include 'messages.h'
  include 'mpif.h'
  include 'files.h'

  integer :: ios, nRec, IOL

  integer :: i, j, k, nt, ierror, ib, jb, i1, i2, m, mm, nn, mw

  logical :: GetNewTop = .false., GetNewSig = .false., ExistGDAS = .false., &
    ExistGANLCPT = .false., ExistGANLSMT = .false., ExistGANL = .false., &
    VerticalInterp = .true.

  character (len = 12) :: Tdef = '  z         '
  real (kind = p_r8) :: MagPs

  type(ChoppingNameListData) :: var
  namelist /ChoppingNameList/ var

  character(len = *), parameter :: headerMsg = 'Chopping            | '


contains


  function getNameChopping() result(returnModuleName)
    !# Returns Chopping Module Name
    !# ---
    !# @info
    !# **Brief:** Returns Chopping Module Name. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 <br>
    !# @endin
    implicit none
    character (len = maxModuleNameLength) :: returnModuleName
    
    returnModuleName = "Chopping"
  end function getNameChopping


  subroutine initChopping(nameListFileUnit)
    !# Initializes Chopping Module
    !# ---
    !# @info
    !# **Brief:** Initialization of Chopping module, defined in PRE_run.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 <br>
    !# @endin
    implicit none
    integer, intent(in) :: nameListFileUnit

    read(unit = nameListFileUnit, nml = ChoppingNameList)

  end subroutine initChopping


  function shouldRunChopping() result(shouldRun)
    !# Returns true if Module Should Run as a dependency
    !# ---
    !# @info
    !# **Brief:** Returns true if Module Should Run as a dependency, when it
    !# does not generated its out files and was not marked to run. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: jul/2019 <br>
    !# @endin
    implicit none
    logical :: shouldRun

    shouldRun = .not. fileExists(getOutFileName())
  end function shouldRunChopping


  function getOutFileName() result(outFilename)
    !# Gets Out Filename
    !# ---
    !# @info
    !# **Brief:** Gets Out Filename. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: may/2019 <br>
    !# @endin
    implicit none
    character(len = maxPathLength) :: outFilename

    call fillDataOutFileName(var)
    outFilename = trim(var%dirout) // trim(DataOut)
  end function getOutFileName


  subroutine printNameList()
    !# Prints namelist of Chopping Module
    !# ---
    !# @info
    !# **Brief:** Prints namelist of Chopping Module, defined in PRE.nml. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: mar/2019 <br>
    !# @endin  
    implicit none

    if (isMasterProc()) then
      write (unit = p_nfprt, fmt = '(A)')      ' '
      write (unit = p_nfprt, fmt = '(A)')      ' &ChoppingNameList'
      write (unit = p_nfprt, fmt = '(A,I5)')   '  mEndInp            = ', var%mEndInp
      write (unit = p_nfprt, fmt = '(A,I5)')   '  mEndOut            = ', var%mEndOut
      write (unit = p_nfprt, fmt = '(A,I5)')   '  kMaxOut            = ', var%kMaxOut
      write (unit = p_nfprt, fmt = '(A,I5)')   '  mEndMin            = ', var%mEndMin
      write (unit = p_nfprt, fmt = '(A,I5)')   '  mEndCut            = ', var%mEndCut
      write (unit = p_nfprt, fmt = '(A,I5)')   '  iter               = ', var%iter
      write (unit = p_nfprt, fmt = '(A,F7.3)') '  smthPerCut         = ', var%smthPerCut
      write (unit = p_nfprt, fmt = '(A,L6)')   '  getOzone           = ', var%getOzone
      write (unit = p_nfprt, fmt = '(A,L6)')   '  getTracers         = ', var%getTracers
      write (unit = p_nfprt, fmt = '(A,L6)')   '  grADS              = ', var%grADS
      write (unit = p_nfprt, fmt = '(A,L6)')   '  grADSOnly          = ', var%grADSOnly
      write (unit = p_nfprt, fmt = '(A,L6)')   '  gdasOnly           = ', var%gdasOnly
      write (unit = p_nfprt, fmt = '(A,A )')   '  dataGDAS           = ', trim(var%dataGDAS)
      write (unit = p_nfprt, fmt = '(A,L6)')   '  smoothTopo         = ', var%smoothTopo
      write (unit = p_nfprt, fmt = '(A,L6)')   '  rmGANL             = ', var%rmGANL
      write (unit = p_nfprt, fmt = '(A,L6)')   '  linearGrid         = ', var%linearGrid
      write (unit = p_nfprt, fmt = '(2A)')     '  dateLabel          = ', var%dateLabel
      write (unit = p_nfprt, fmt = '(2A)')     '  utc                = ', var%utc
      write (unit = p_nfprt, fmt = '(2A)')     '  NCEPName           = ', trim(var%NCEPName)
      write (unit = p_nfprt, fmt = '(2A)')     '  StrFormat          = ', trim(var%StrFormat)     
      write (unit = p_nfprt, fmt = '(A,L6)')   '  givenfouriergroups = ', var%givenfouriergroups
      write (unit = p_nfprt, fmt = '(A,I5)')   '  nProc_vert         = ', var%nProc_vert
      write (unit = p_nfprt, fmt = '(A,I5)')   '  ibdim_size         = ', var%ibdim_size
      write (unit = p_nfprt, fmt = '(A,I5)')   '  tamBlock           = ', var%tamBlock
      write (unit = p_nfprt, fmt = '(2A)')     '  dGDInp             = ', trim(var%dGDInp)
      write (unit = p_nfprt, fmt = '(2A)')     '  dirInp             = ', trim(var%dirInp)
      write (unit = p_nfprt, fmt = '(2A)')     '  dirOut             = ', trim(var%dirOut)
      write (unit = p_nfprt, fmt = '(2A)')     '  dirTop             = ', trim(var%dirTop)
      write (unit = p_nfprt, fmt = '(2A)')     '  dirSig             = ', trim(var%dirSig)
      write (unit = p_nfprt, fmt = '(2A)')     '  dirGrd             = ', trim(var%dirGrd)


      !      write (unit = p_nfprt, fmt = '(2A)')     '  prefix             = ', trim(var%prefix)
      !      write (unit = p_nfprt, fmt = '(2A)')     '  DirMain    = ', trim(DirMain)
      !      write (unit = p_nfprt, fmt = '(2A)')     '  DirHome    = ', trim(DirHome)
      write (unit = p_nfprt, fmt = '(A)')      ' /'
    endif
    !    WRITE (UNIT = p_nfprt, FMT = '(/,A)')  ' &SoilMoistureWeeklyNameList'
    !    WRITE (UNIT = p_nfprt, FMT = '(A,I6)') '      xDim = ', xDim
    !    WRITE (UNIT = p_nfprt, FMT = '(A,I6)') '      yDim = ', yDim
    !    WRITE (UNIT = p_nfprt, FMT = '(A,I6)') '      zDim = ', zDim
    !    WRITE (UNIT = p_nfprt, FMT = '(A,I6)') '      xMax = ', xMax
    !    WRITE (UNIT = p_nfprt, FMT = '(A,I6)') '      yMax = ', yMax
    !    WRITE (UNIT = p_nfprt, FMT = '(A,L6)') '     grADS = ', grADS
    !    WRITE (UNIT = p_nfprt, FMT = '(A,L6)') '    linear = ', linear
    !
    !    WRITE (UNIT = p_nfprt, FMT = '(A)') '  dirPreIn = ' // TRIM(dirPreIn)
    !    WRITE (UNIT = p_nfprt, FMT = '(A)') 'dirPreTemp = ' // TRIM(dirPreTemp)
    !    WRITE (UNIT = p_nfprt, FMT = '(A)') ' dirPreOut = ' // TRIM(dirPreOut)
    !    WRITE (UNIT = p_nfprt, FMT = '(A)') '   varName = ' // TRIM(varName)
    !    WRITE (UNIT = p_nfprt, FMT = '(A)') '      date = ' // TRIM(date)
  end subroutine printNameList


  function generateChopping() result(isExecOk)
    !# Generates Chopping
    !# ---
    !# @info
    !# **Brief:** Generates Chopping output. This subroutine is the main method
    !# for use this module. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 <br>
    !# @endin
    implicit none

    integer :: nftop, nfozr, nfsig, nftrr, lrec
    logical :: isExecOk
    integer :: mEndInpAux, mEndOutAux, kMaxInpAux, kMaxOutAux, mEndCutAux

    isExecOk = .false.

    call initializeParallelismVariables()

    mEndInpAux = var%mEndInp
    kMaxInpAux = var%kMaxInp
    mEndOutAux = var%MendOut
    kMaxOutAux = var%kMaxOut
    mEndCutAux = var%mEndCut
    !  if (isMasterProc()) call printNameList(var)

    if(.not. InitParameters(1, var)) return
    call InitAll (ImaxOut, var)
    call GetArrays(var%kMaxInp, var%kMaxOut)

  
    if(trim(var%dataGDAS) == 'Grid')then
      inquire (file = trim(var%dGDInp) // trim(var%gdasinp), EXIST = ExistGDAS)
      inquire (file = trim(var%dirInp) // trim(DataCPT), EXIST = ExistGANLCPT)
      inquire (file = trim(var%dirInp) // trim(DataInp), EXIST = ExistGANLSMT)

      if (isMasterProc()) then
        write (unit = p_nfprt, FMT = '(A,L6)') trim(var%dGDInp) // trim(var%gdasinp), ExistGDAS
        write (unit = p_nfprt, FMT = '(A,L6)') trim(var%dirInp) // trim(DataCPT), ExistGANLCPT
        write (unit = p_nfprt, FMT = '(A,L6)') trim(var%dirInp) // trim(DataInp), ExistGANLSMT
      endif
      if (.not.ExistGDAS .and. .not.ExistGANLCPT .and. .not.ExistGANLSMT) then
        call fatalError(headerMsg, &
          ' The NCEP Input file Does not Exist and' // &
          ' $The CPTEC Input file Does not Exist Also and' // &
          ' $The CPTEC Topo-Smoothed Input file Does not Exist Also' // &
          ' $No NCEP or GANL File' )
        return
      end if

      ExistGANL = ExistGANLCPT .or. ExistGANLSMT
      if ((ExistGANLSMT).and. (.not.ExistGDAS)) then
        if (var%smoothTopo) var%smoothTopo = .false.
      else
        if (ExistGANLCPT) then
          DataInp = DataCPT
          if (.not.var%smoothTopo .and. var%mEndOut > var%mEndMin) var%smoothTopo = .true.
        end if
      end if
      call MPI_BARRIER(mpiCommGroup, ierror)

      if (ExistGDAS .and. .not.ExistGANL) then
        call GetSpHuTracers(var%kMaxInp, var%kMaxOut)

        if ( .not. GDAStoGANL2() ) return

      end if
    end if

    call ClsArrays
    call ClsSpHuTracers
    call Clear_Sizes()
    call Clear_Utils()
    call Clear_Communications()
    call Clear_SpecDynamics()
    call Clear_Transform()
    call MPI_BARRIER(mpiCommGroup, ierror)

    var%mEndInp = mEndInpAux
    var%kMaxInp = kMaxInpAux
    var%MendOut = mEndOutAux
    var%kMaxOut = kMaxOutAux
    var%mEndCut = mEndCutAux
    !  if (isMasterProc()) call printNameList(var)

    if(.not. InitParameters(0, var)) return
    call InitAll (ImaxOut, var)
    call GetArrays(var%kMaxInp, var%kMaxOut)

    !  if (isMasterProc()) call writectl(1)
    inquire (file = trim(var%dGDInp) // trim(var%gdasinp), EXIST = ExistGDAS)
    inquire (file = trim(var%dirInp) // trim(DataCPT), EXIST = ExistGANLCPT)
    inquire (file = trim(var%dirInp) // trim(DataInp), EXIST = ExistGANLSMT)

    if (isMasterProc()) then
      write (unit = p_nfprt, FMT = '(A,L6)') trim(var%dGDInp) // trim(var%gdasinp), ExistGDAS
      write (unit = p_nfprt, FMT = '(A,L6)') trim(var%dirInp) // trim(DataCPT), ExistGANLCPT
      write (unit = p_nfprt, FMT = '(A,L6)') trim(var%dirInp) // trim(DataInp), ExistGANLSMT
    endif

    if (.not.ExistGDAS .and. .not.ExistGANLCPT .and. .not.ExistGANLSMT) then
      call fatalError(headerMsg, &
          ' The NCEP Input file Does not Exist and' // &
          ' $The CPTEC Input file Does not Exist Also and' // &
          ' $The CPTEC Topo-Smoothed Input file Does not Exist Also' // &
          ' $No NCEP or GANL File' )
        return
    end if

    ExistGANL = ExistGANLCPT .or. ExistGANLSMT
    if ((ExistGANLSMT).and. (.not.ExistGDAS)) then
      if (var%smoothTopo) var%smoothTopo = .false.
    else
      if (ExistGANLCPT) then
        DataInp = DataCPT
        if (.not.var%smoothTopo .and. var%mEndOut > var%mEndMin) var%smoothTopo = .true.
      end if
    end if

    call MPI_BARRIER(mpiCommGroup, ierror)
    if (ExistGDAS .and. .not.ExistGANL) then
      if(trim(var%dataGDAS) == 'Spec')then
        if(isMasterProc()) then
          if ( .not. GDAStoGANL() ) return
        else
          if (var%getOzone) then
            NTracers = NTracers + 1
            if (var%getTracers) NTracers = NTracers + 1
          endif
        endif
        call MPI_BARRIER(mpiCommGroup, ierror)
        call GetSpHuTracers(var%kMaxInp, var%kMaxOut)
      end if
    else

      if (isMasterProc()) write (unit = p_nfprt, FMT = '(/,A)') ' GANL file Already Exists'

      inquire (file = trim(var%dirInp) // trim(OzonInp), EXIST = var%getOzone)
      inquire (file = trim(var%dirInp) // trim(TracInp), EXIST = var%getTracers)
      if (var%getOzone) then
        NTracers = NTracers + 1
        if (var%getTracers) then
          if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A)') &
            ' Considering Just One Other Tracer Than Ozone'
          NTracers = NTracers + 1
        end if
      end if
      call GetSpHuTracers(var%kMaxInp, var%kMaxOut)

      mw = min(var%mEndOut+1,var%mEndInp+1)

      if (var%getOzone) then

        if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A)') ' Get Ozone from ' // trim(OzonInp)
        nfozr = openFile(trim(var%dirInp) // trim(OzonInp), 'unformatted', 'sequential', -1, 'read', 'old')
        if(nfozr < 0) return
        IF(TRIM(var%dataGDAS) == 'Grid')THEN
           nt=2
           DO k=1,var%kMaxInp
              READ (UNIT=nfozr) qWorkInp
              IF (k.ge.myfirstlev.and.k.le.mylastlev) THEN
                 DO mm=1,mymmax
                    m = msinproc(mm,myid_four)
                    i2 = 2*mymnmap(mm,m)-1
                    IF (m.gt.var%mEndInp+1) THEN 
                       DO nn=0,2*(var%mEndOut+1-m)+1
                        qSpHuInp(i2+nn,k+1-myfirstlev,nt) = 0.0_p_r8
                       ENDDO
                    ELSE
                       i1 = 2*mnmap(m,m)-1
                       DO nn=0,2*(mw-m)+1
                          qSpHuInp(i2+nn,k+1-myfirstlev,nt) = qWorkInp(i1+nn)
                       ENDDO
                       DO nn=2*(mw-m)+2,2*(var%mEndOut+1-m)+1
                          qSpHuInp(i2+nn,k+1-myfirstlev,nt) = 0.0_p_r8
                       ENDDO
                    ENDIF
                 ENDDO
              END IF
           END DO
        ELSE
           nt=2
           DO k=1,var%kMaxInp
              READ (UNIT=nfozr) qWorkInp
              IF (k.ge.myfirstlev.and.k.le.mylastlev) THEN
                 DO mm=1,mymmax
                    m = msinproc(mm,myid_four)
                    i2 = 2*mymnmap(mm,m)-1
                    IF (m.gt.var%mEndInp+1) THEN 
                       DO nn=0,2*(var%mEndOut+1-m)+1
                        qSpHuInp(i2+nn,k+1-myfirstlev,nt) = 0.0_p_r8
                       ENDDO
                    ELSE
                       i1 = 2*mnmap(m,m)-1
                       DO nn=0,2*(mw-m)+1
                          qSpHuInp(i2+nn,k+1-myfirstlev,nt) = qWorkInp(i1+nn)
                       ENDDO
                       DO nn=2*(mw-m)+2,2*(var%mEndOut+1-m)+1
                          qSpHuInp(i2+nn,k+1-myfirstlev,nt) = 0.0_p_r8
                       ENDDO
                    ENDIF
                 ENDDO
              END IF
           END DO
        END IF
        CLOSE(UNIT=nfozr)
   
   
        mw = min(var%mEndOut+1,var%mEndInp+1)

!        nt = 2
!        do k = 1, var%kMaxInp
!          read (unit = nfozr) qWorkprOut
!          !qSpHuInp(:,k,nt)=qWorkprOut
!          if (k.ge.myfirstlev.and.k.le.mylastlev) then
!            do mm = 1, mymmax
!              m = msinproc(mm, myid_four)
!              i1 = 2 * mnMap_out(m, m) - 1
!              i2 = 2 * mymnmap(mm, m) - 1
!              do nn = 0, 2 * (mMax - m) + 1
!                qSpHuInp(i2 + nn, k + 1 - myfirstlev, nt) = qWorkprOut(i1 + nn)
!              enddo
!            enddo
!          end if
!
!        end do
!        close(unit = nfozr)

        if (var%getTracers) then
          if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A)') ' Get Tracers from ' // trim(TracInp)
          nftrr = openFile(trim(var%dirInp) // trim(TracInp), 'unformatted', 'sequential', -1, 'read', 'old')
          if(nftrr < 0) return

          IF(TRIM(var%dataGDAS) == 'Grid')THEN
             Do nt=3,NTracers
                DO k=1,var%kMaxInp
                   READ (UNIT=nftrr) qWorkInp
                   IF (k.ge.myfirstlev.and.k.le.mylastlev) THEN
                   DO mm=1,mymmax
                      m = msinproc(mm,myid_four)
                      i2 = 2*mymnmap(mm,m)-1
                      IF (m.gt.var%mEndInp+1) THEN 
                         DO nn=0,2*(var%mEndOut+1-m)+1
                            qSpHuInp(i2+nn,k+1-myfirstlev,nt) = 0.0_p_r8
                         ENDDO
                      ELSE
                         i1 = 2*mnmap(m,m)-1
                         DO nn=0,2*(mw-m)+1
                           qSpHuInp(i2+nn,k+1-myfirstlev,nt) = qWorkInp(i1+nn)
                         ENDDO
                         DO nn=2*(mw-m)+2,2*(var%mEndOut+1-m)+1
                            qSpHuInp(i2+nn,k+1-myfirstlev,nt) = 0.0_p_r8
                         ENDDO
                      ENDIF
                   ENDDO
                   END IF
                END DO
             END DO
          ELSE
             Do nt=3,NTracers
                DO k=1,var%kMaxInp
                   READ (UNIT=nftrr) qWorkInp
                   IF (k.ge.myfirstlev.and.k.le.mylastlev) THEN
                   DO mm=1,mymmax
                      m = msinproc(mm,myid_four)
                      i2 = 2*mymnmap(mm,m)-1
                      IF (m.gt.var%mEndInp+1) THEN 
                         DO nn=0,2*(var%mEndOut+1-m)+1
                            qSpHuInp(i2+nn,k+1-myfirstlev,nt) = 0.0_p_r8
                         ENDDO
                      ELSE
                         i1 = 2*mnmap(m,m)-1
                         DO nn=0,2*(mw-m)+1
                           qSpHuInp(i2+nn,k+1-myfirstlev,nt) = qWorkInp(i1+nn)
                         ENDDO
                         DO nn=2*(mw-m)+2,2*(var%mEndOut+1-m)+1
                            qSpHuInp(i2+nn,k+1-myfirstlev,nt) = 0.0_p_r8
                         ENDDO
                      ENDIF
                   ENDDO
                   END IF
                END DO
             END DO
          END IF
          CLOSE(UNIT=nftrr)

!          nt = 3
!          do k = 1, var%kMaxInp
!            read (unit = nftrr) qWorkprOut
!            !  qSpHuInp(:,k,nt)=qWorkInp
!            if (k.ge.myfirstlev.and.k.le.mylastlev) then
!              do mm = 1, mymmax
!                m = msinproc(mm, myid_four)
!                i1 = 2 * mnMap_out(m, m) - 1
!                i2 = 2 * mymnmap(mm, m) - 1
!                do nn = 0, 2 * (mMax - m) + 1
!                  qSpHuInp(i2 + nn, k + 1 - myfirstlev, nt) = qWorkprOut(i1 + nn)
!                enddo
!              enddo
!            end if
!          end do
!          close(unit = nftrr)
        else
          if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A)') ' Other Tracers file Does not Exist'
        end if

      else

        if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A)') &
          ' Ozone file Does not Exist. Ignore Other Tracers'
        var%getTracers = .false.
        NTracers = 1

      end if

    end if
    if (isMasterProc()) write (unit = p_nfprt, FMT = '(/,A,I5)') ' NTracers = ', NTracers

    call MPI_BARRIER(mpiCommGroup, ierror)

    if (.not. ICRead_and_Chop(GetNewTop) ) return


    !TopoOut for New Grid Topography:

    inquire(IOLENGTH = lrec)gWorkprOut
    inquire (file = trim(var%dirtop) // trim(DataTopG), EXIST = GetNewTop)
    if(GetNewTop)then

      nftop = openFile(trim(var%dirtop) // trim(DataTopG), 'unformatted', 'direct', lrec, 'read', 'old')
      if(nftop < 0) return

      read  (unit = nftop, rec = 1) gWorkprOut
      do j = 1, jmaxout
        do i = 1, imaxout
          gTopoOutGaus8(i, j) = real(gWorkprOut(i, jmaxout + 1 - j), kind = p_r8)
        end do
      end do

      if(isMasterProc()) then
        write (unit = p_nfprt, FMT = '(/,A)') ' TopoOut for New Grid Topography:'
        write (unit = p_nfprt, FMT = '(1P3G12.5)') gWorkprOut(1, 1), &
          minval(gWorkprOut), maxval(gWorkprOut)
      endif
      close (unit = nftop)

      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb)
          gTopoOutGaus(ib, jb) = gTopoOutGaus8(i, j)
        end do
      end do
      do jb = 1, jbMax
        do ib = 1, Ibmaxperjb(jb)
          gTopoOut(ib, jb) = gTopoOutGaus(ib, jb)
        end do
      end do
      if(isMasterProc()) then
        write (unit = p_nfprt, FMT = '(/,A)') ' TopoOut for New Grid Block Topography:'
        write (unit = p_nfprt, FMT = '(1P3G12.5)') gTopoOut(1, 1), &
          minval(gTopoOut), maxval(gTopoOut)
      endif

    end if

    inquire (file = trim(var%dirSig)//'/'// trim(DataSig), EXIST = GetNewSig)
    if (GetNewSig) then

      if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A)') ' Getting New Delta Sigma'
        nfsig = openFile(trim(var%dirSig) // trim(DataSig), 'formatted', 'sequential', -1, 'read', 'old')
        if(nfsig < 0) return
            if(isMasterProc()) write(unit = p_nfprt, FMT = '(A)') ' SigInterOut(k), SigLayerOut(k)'
            do k = 1, var%kMaxOut + 1
              read(unit = nfsig, FMT = *) SigInterOut(k), SigLayerOut(k)
              SigInterOut(k) = SigInterOut(k) / 100._p_r4  ! transform from Pa to mbar
              if(isMasterProc()) write(unit = p_nfprt, FMT = *) SigInterOut(k), SigLayerOut(k)
            enddo
      !PK_sig read  (unit = nfsig, fmt = '(5F9.6)') DelSigmaOut

      close (unit = nfsig)

      if (var%kMaxOut == var%kMaxInp) then
        !IF (MAXVAL(ABS(DelSigmaOut-DelSigmaInp)) < 1.0E-04_p_r8) THEN
        if (isMasterProc()) then
          write (unit = p_nfprt, fmt = '(/,A)') ' MaxOut = kMaxInp And DelSima Is Quite The Same'
          !WRITE (UNIT=p_nfprt, FMT='(A,1PG12.5,/)') ' MAXVAL(ABS(DelSigmaOut-DelSigmaInp)) = ', &
          !                                 MAXVAL(ABS(DelSigmaOut-DelSigmaInp))
        end if
        !PK_sig call NewSigma(var%kMaxOut)
        SigInterInp=SigInterOut
        SigLayerInp=SigLayerOut
        DelSigmaInp = DelSigmaOut
        !ELSE
        !   CALL NewSigma
        !ENDIF
      else
        !call NewSigma(var%kMaxOut)
      end if

    else

      if (var%kMaxOut /= var%kMaxInp) then
        call fatalError(headerMsg,  'Error in Getting New Sigma: kMaxInp /= kMaxOut')
        call msgInLineFormatOut(headerMsg // ' kMaxInp = ',  '(A)')
        call msgInLineFormatOut(var%kMaxInp, '(I5)')
        call msgInLineFormatOut(', kMaxOut = ',  '(A)')
        call msgInLineFormatOut(var%kMaxOut, '(I5)')
        call msgNewLine()
        !        stop
        return
      end if
      if (IDVCInp == 0_p_i8 .or. IDVCInp == 1_p_i8) then
         ! Sigma Interfaces   (kmax+1) - SigIInp
         ! Sigma Layers       (kmax  ) - SigLInp
         !PK_sig call NewSigma(var%kMaxOut)
         SigInterOut = SigInterInp
         SigLayerOut = SigLayerInp
         DelSigmaOut = DelSigmaInp

      else if (IDVCInp == 2_p_i8) then
         ! Hybrid Interface A (kmax+1) - SigIInp
         ! Hybrid Interface B (kmax+1) - SigLInp
         SigInterOut = SigInterInp
         SigLayerOut = SigLayerInp
         DelSigmaOut = DelSigmaInp
      else
      !      stop ' ** (Error) **'
      call fatalError(headerMsg, 'Invalid IDVCInp. Aborting Chopping!')
      return
    end if

    end if

    if (var%kMaxOut == var%kMaxInp .and. &
      (.not.GetNewTop .and. .not.var%smoothTopo)) VerticalInterp = .false.
    if (GetNewTop)VerticalInterp = .true.

    !WRITE (UNIT=p_nfprt, FMT=*) VerticalInterp .OR. GrADS,VerticalInterp , GrADS
    if (VerticalInterp .or. var%grADS) then

      call ICRecomposition

      if (var%grADS) then
        write (Tdef(1:2), '(I2.2)') DateCurrent(1)
        write (Tdef(4:5), '(I2.2)') DateCurrent(3)
        write (Tdef(6:8), '(A3)')   MonChar(DateCurrent(2))
        write (Tdef(9:12), '(I4.4)') DateCurrent(4)
        if (.not. GetGrADSInp() ) return
        if (var%grADSOnly) then
          isExecOk = .true.
          return
        endif
      end if

    end if

    if (VerticalInterp) then

      if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A)') 'Doing Vertical Interpolation'

      if (var%gdasonly) then
       ! Sigma case, calculate pressure for each layer directly
        do j = 1, Jbmax
          do k = 1, var%kMaxInp
            do I = 1, Ibmaxperjb(j)
              gPresInp(i, k, j) = gPsfcInp(i, j) * SigLayerInp(k)
            end do
          end do
        end do
	
      else	
        ! Calculate pressure for each layer according to vertical coordinate
        if (IDVCInp == 0_p_i8 .or. IDVCInp == 1_p_i8) then
          ! Sigma case, calculate pressure for each layer directly
          do j = 1, Jbmax
            do k = 1, var%kMaxInp
              do I = 1, Ibmaxperjb(j)
                gPresInp(i, k, j) = gPsfcInp(i, j) * SigLayerInp(k)
              end do
            end do
          end do
        else if (IDVCInp == 2_p_i8) then
          ! Sigma-p case, first calculate pressure for each interface and then
          ! calculate pressure on each layer according to vertical structure
          do j = 1, Jbmax
            do k = 1, KmaxInpp
              do I = 1, Ibmaxperjb(j)
                gPresInpp(i, k, j) = SigInterInp(k) + gPsfcInp(i, j) * SigLayerInp(k)
              end do
            end do
          end do
          if (IDSLInp == 2_p_i8) then
            ! Mean over two interfaces
            do j = 1, Jbmax
              do k = 1, var%kMaxInp
                do I = 1, Ibmaxperjb(j)
                  gPresInp(i, k, j) = (gPresInpp(i, k, j) + gPresInpp(i, k + 1, j)) / 2.0_p_r8
                end do
              end do
            end do
          else
            ! Phillips interpolation over two interfaces
            do j = 1, Jbmax
              do k = 1, var%kMaxInp
                do I = 1, Ibmaxperjb(j)
                  gPresInp(i, k, j) = (gPresInpp(i, k, j) + gPresInpp(i, k + 1, j)) / 2.0_p_r8

                 !                      gPresInp(i,k,j) = ((gPresInpp(i,k,j)**RoCp1-gPresInpp(i,k+1,j)**RoCp1)/ &
                 !                           (RoCp1*(gPresInpp(i,k,j)-gPresInpp(i,k+1,j))))**RoCpR
                end do
              end do
            end do
          end if
        end if
      end if
      
      call NewPs(var%kMaxInp)

      IF ( var%gdasonly) THEN

          ! Sigma-p case, first calculate pressure for each interface and then
          ! calculate pressure on each layer according to vertical structure
          !
          DO j=1,Jbmax
             DO k=1,KmaxOutp
                DO I=1,Ibmaxperjb(j)
                   gPresaux(i,k)=SigInterOut(k)+gPsfcOut(i,j)*SigLayerOut(k)
                END DO
             END DO
             IF ( IDSLInp == 2_p_i8 ) THEN
               ! Mean over two interfaces
                DO k=1,var%kMaxOut
                   DO I=1,Ibmaxperjb(j)
                      gPresOut(i,k,j)= 0.5_p_r8 * (gPresaux(i,k)+gPresaux(i,K+1))
                   END DO
                END DO
             ELSE
               ! Phillips interpolation over two interfaces
                DO k=1,var%kMaxOut
                   DO I=1,Ibmaxperjb(j)
                      gPresOut(i,k,j)= ((gPresaux(i,k)**RoCp1-gPresaux(i,k+1)**RoCp1)/ &
                                       (RoCp1*(gPresaux(i,k)-gPresaux(i,k+1))))**RoCpR
                   END DO
                END DO
             ENDIF
          END DO
      ELSE
        ! Calculate pressure for each layer according to vertical coordinate
        IF ( IDVCInp == 0_p_i8 .OR. IDVCInp == 1_p_i8 ) THEN
          ! Sigma case, calculate pressure for each layer directly
          do j = 1, Jbmax
            do k = 1, var%kMaxOut
              do I = 1, Ibmaxperjb(j)
                gPresOut(i,k,j)=gPsfcOut(i,j)*SigLayerOut(k)
                !MagPs = maxval(gPsfcOut)
                !if(MagPs > p0mb) then
                !  gPresOut(i, k, j) = SigInterOut(k) * gPsfcOut(i, j) + SigLayerOut(k) * gPsfcOut(i, j)
                !else
                !  gPresOut(i, k, j) = SigInterOut(k) * gPsfcOut(i, j) + SigLayerOut(k) * gPsfcOut(i, j)
                !end if
                !if(isMasterProc()) print*, k, gPsfcOut(i, j), SigInterOut(k) * gPsfcOut(i, j), SigLayerOut(k) * gPsfcOut(i, j), gPresOut(i, k, j)
                !                gPresOutInter(nlevp-k+1)  = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*ps(i,j)
              end do
            end do
          end do
        ELSE IF ( IDVCInp == 2_p_i8 ) THEN
             ! Sigma-p case, first calculate pressure for each interface and then
             ! calculate pressure on each layer according to vertical structure
      
             DO j=1,Jbmax
                DO k=1,KmaxOutp
                   DO I=1,Ibmaxperjb(j)
                      gPresaux(i,k)=SigInterOut(k)+gPsfcOut(i,j)*SigLayerOut(k)
                   END DO
                END DO
                IF ( IDSLInp == 2_p_i8 ) THEN
                   DO k=1,var%kMaxOut
                      DO I=1,Ibmaxperjb(j)
                         gPresOut(i,k,j)= 0.5_p_r8 * (gPresaux(i,k)+gPresaux(i,K+1))
                      END DO
                   END DO
                ELSE
                   DO k=1,var%kMaxOut
                      DO I=1,Ibmaxperjb(j)
                        gPresOut(i,k,j)= 0.5_p_r8 * (gPresaux(i,k)+gPresaux(i,K+1))
!                         gPresOut(i,k,j)= ((gPresaux(i,k)**RoCp1-gPresaux(i,k+1)**RoCp1)/ &
!                                         (RoCp1*(gPresaux(i,k)-gPresaux(i,k+1))))**RoCpR
                      END DO
                   END DO
                ENDIF
             END DO
         END IF
      END IF

      !   DO j=1,Jbmax
      !       DO i=1,Ibmaxperjb(j)
      !          IF (iPerIJB(i,j).EQ.123.AND.jPerIJB(i,j).EQ.56) THEN
      !             DO k = var%kMaxInp,1,-1
      !                WRITE(70,*) gPresInp(i,k,j),gUvelInp(i,k,j)
      !                WRITE(72,*) gPresInp(i,k,j),gVvelInp(i,k,j)
      !                WRITE(74,*) gPresInp(i,k,j),gTvirInp(i,k,j)
      !             ENDDO
      !             DO k = var%kMaxOut,1,-1
      !                WRITE(71,*) gPresOut(i,k,j),gUvelOut(i,k,j)
      !                WRITE(73,*) gPresOut(i,k,j),gVvelOut(i,k,j)
      !                WRITE(75,*) gPresOut(i,k,j),gTvirOut(i,k,j)
      !             ENDDO
      !          ENDIF
      !       END DO
      !   END DO

      IF(TRIM(var%dataGDAS) == 'Grid')THEN
        DO j=1,Jbmax
           DO k=1,var%kMaxInp
              DO i=1,Ibmaxperjb(j)
                 gTvirInp(i,k,j)=gTvirInp(i,k,j)!/(1.0_p_r8+cTv*gSpHuInp(i,k,j,1))
              END DO
           END DO
        END DO
      ELSE
        DO j=1,Jbmax
           DO k=1,var%kMaxInp
              DO i=1,Ibmaxperjb(j)
                 gTvirInp(i,k,j)=gTvirInp(i,k,j)/(1.0_p_r8+cTv*gSpHuInp(i,k,j,1))
              END DO
           END DO
        END DO
      END IF

      do j = 1, Jbmax
        call VertSigmaInter (ibmax, Ibmaxperjb(j), &
          var%kMaxInp, var%kMaxOut, NTracers, &
          gPresInp(:, :, j), gUvelInp(:, :, j), gVvelInp(:, :, j), &
          gTvirInp(:, :, j), gSpHuInp(:, :, j, :), &
          gPresOut(:, :, j), gUvelOut(:, :, j), gVvelOut(:, :, j), &
          gTvirOut(:, :, j), gSpHuOut(:, :, j, :))
      end do

      do j = 1, Jbmax
        do k = 1, var%kMaxInp
          do i = 1, Ibmaxperjb(j)
            gTvirInp(i, k, j) = gTvirInp(i, k, j) * (1.0_p_r8 + cTv * gSpHuInp(i, k, j, 1))
          end do
        end do
      end do
      do j = 1, Jbmax
        do k = 1, var%kMaxOut
          do i = 1, Ibmaxperjb(j)
            gTvirOut(i, k, j) = gTvirOut(i, k, j) * (1.0_p_r8 + cTv * gSpHuOut(i, k, j, 1))
          end do
        end do
      end do

      if (var%grADS) then
        write (Tdef(1:2), '(I2.2)') DateCurrent(1)
        write (Tdef(4:5), '(I2.2)') DateCurrent(3)
        write (Tdef(6:8), '(A3)')   MonChar(DateCurrent(2))
        write (Tdef(9:12), '(I4.4)') DateCurrent(4)
        if(.not. GetGrADSOut() ) return
        if (var%grADSOnly) then
          isExecOk = .true.
          return
        endif
      end if

      call ICDecomposition

    end if

    if (.not. ICWrite() ) return

    call ClsArrays

    call ClsSpHuTracers

    isExecOk = .true.
  end function generateChopping


  function GDAStoGANL2() result(isExecOk)
    !# Converts from GDAS to GANL2
    !# ---
    !# @info
    !# **Brief:** Converts from GDAS to GANL2. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 <br>
    !# @endin
    implicit none
    logical :: isExecOk

    real (kind = p_r4), dimension (2 * 100 + 1) :: SiSl
    integer, parameter :: nemsio_realkind = 4
    integer :: iret, nelements
    type(c_nemsio_gfile) :: gfile
    integer :: ios, ios1, ij
    integer :: im, jm, nsoil, fieldsize, ntrac, ierr
    type(nemsio_head) :: gfshead
    type(nemsio_headv) :: gfsheadv
    real (kind = p_r4), allocatable :: buff(:)

    isExecOk = .false.

    call c_nemsio_init(iret = iret)
    ! Inicializa a lib nemsio
    
    ! call nemsio_open(gfile,TRIM(var%dGDInp)//TRIM(var%gdasinp),'READ',MPI_COMM_WORLD,iret=iret)
    call c_nemsio_open(gfile, trim(var%dGDInp) // trim(var%gdasinp), 'READ', iret = iret)
    ! Abertura do arquivo nemsio

    if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A)') ' Getting GANL from GDAS NCEP File'

    !  open (read) nemsio grid file headers
    call c_nemsio_getfilehead(gfile, &
      idate = gfshead%idate, &
      nfhour = gfshead%nfhour, &
      nfminute = gfshead%nfminute, &
      nfsecondn = gfshead%nfsecondn, &
      nfsecondd = gfshead%nfsecondd, &
      version = gfshead%version, &
      nrec = gfshead%nrec, &
      dimx = gfshead%dimx, &
      dimy = gfshead%dimy, &
      dimz = gfshead%dimz, &
      jcap = gfshead%jcap, &
      ntrac = gfshead%ntrac, &
      ncldt = gfshead%ncldt, &
      nsoil = gfshead%nsoil, &
      idsl = gfshead%idsl, &
      idvc = gfshead%idvc, &
      idvm = gfshead%idvm, &
      idrt = gfshead%idrt, &
      extrameta = gfshead%extrameta, &
      nmetavari = gfshead%nmetavari, &
      nmetavarr = gfshead%nmetavarr, &
      nmetavarl = gfshead%nmetavarl, &
      nmetavarr8 = gfshead%nmetavarr8, &
      nmetaaryi = gfshead%nmetaaryi, &
      nmetaaryr = gfshead%nmetaaryr, &
      iret = ios)

    call nemsio_getheadvar(gfile, 'fhour', gfshead%fhour, iret = ios)
    if(ios/=0) gfshead%fhour = gfshead%nfhour + gfshead%nfminute / 60.          &
      & + gfshead%nfsecondn / (3600. * gfshead%nfsecondd)

    !         call nemsio_getheadvar(gfile,'dimx',    gfshead%latb,iret=ios)
    !         call nemsio_getheadvar(gfile,'dimy',    gfshead%LONB,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'LEVS',    gfshead%LEVS,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'ITRUN',   gfshead%ITRUN,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'IORDER',  gfshead%IORDER,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'IREALF',  gfshead%IREALF,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'IGEN',    gfshead%IGEN,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'LATF',    gfshead%LATF,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'LONF',    gfshead%LONF,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'LATR',    gfshead%LATR,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'LONR',    gfshead%LONR,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'ICEN2',   gfshead%ICEN2,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'IENS',    gfshead%IENS,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'IDPP',    gfshead%IDPP,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'IDVT',    gfshead%IDVT,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'IDRUN',   gfshead%IDRUN,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'IDUSR',   gfshead%IDUSR,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'PDRYINI', gfshead%PDRYINI,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'IXGR',    gfshead%IXGR,IRET=ios)
    !         CALL NEMSIO_GETHEADVAR(GFILE,'NVCOORD', gfshead%NVCOORD,IRET=ios)

    call nemsio_gfs_alheadv(gfshead, gfsheadv)

    call c_nemsio_getfilehead(GFILE                                    &
      &, RECNAME=gfsheadv%RECNAME              &
      &, RECLEVTYP=gfsheadv%RECLEVTYP          &
      !     &,                         RECLEV=gfsheadv%RECLEV                &
      &, VCOORD = gfsheadv%VCOORD                &
      !     &,                         LAT=gfsheadv%LAT                      &
      !     &,                         LON=gfsheadv%LON                      &
      !     &,                         CPI=gfsheadv%CPI                      &
      !     &,                         RI=gfsheadv%RI                        &
      !     &,                         variname=gfsheadv%variname            &
      !     &,                         varrname=gfsheadv%varrname            &
      !     &,                         varlname=gfsheadv%varlname            &
      !     &,                         varival=gfsheadv%varival              &
      !     &,                         varrval=gfsheadv%varrval              &
      !     &,                         varlval=gfsheadv%varlval              &
      !     &,                         aryiname=gfsheadv%aryiname            &
      !     &,                         aryrname=gfsheadv%aryrname            &
      !     &,                         aryilen=gfsheadv%aryilen              &
      !     &,                         aryrlen=gfsheadv%aryrlen              &
      &, IRET = ios1)

    !GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-GFS-

    ! Descriptor Label (See NMC Office Note 85) (Not Used at CPTEC)

    ! Forecast Time               - TimeOfDay
    ! Initial Date                - DateInitial(1:4)
    !  o Initial Hour             - DateInitial(1)
    !  o Initial Month            - DateInitial(2)
    !  o Initial Day              - DateInitial(3)
    !  o Initial Year             - DateInitial(4)
    ! Sigma Interfaces and Layers - SiSl(1:201)
    ! Extra Information           - Extra(1:44)
    !  o ID Sigma Structure       - Extra(18)
    !    = 1 Phillips
    !    = 2 Mean
    !  o ID Vertical Coordinate   - Extra(19)
    !    = 1 Sigma (0 for old files)
    !    = 2 Sigma-p  READ (UNIT=nfnmc) TimeOfDay, DateInitial, SigIInp, SigLInp

    ! PK READ (UNIT=nfnmc) TimeOfDay, DateInitial, SiSl, Extra
    TimeOfDay = gfshead%nfsecondn
    DateInitial   (1) = gfshead%idate(4)
    DateInitial   (2) = gfshead%idate(2)
    DateInitial   (3) = gfshead%idate(3)
    DateInitial   (4) = gfshead%idate(1)

    SiSl(1:var%kMaxInp + 1) = gfsheadv%VCOORD(1:var%kMaxInp + 1, 1, 1)
    SiSl(var%kMaxInp + 2:2 * var%kMaxInp + 2) = gfsheadv%VCOORD(1:var%kMaxInp + 1, 2, 1)

    if(isMasterProc())then
      write (unit = p_nfprt, FMT = '(A,1X,F15.4)') 'TimeOfDay:', TimeOfDay
      write (unit = p_nfprt, FMT = '(A,1X,4I5)') 'DateInitial:', DateInitial
      write (unit = p_nfprt, FMT = '(A)') 'SiSl-1:'
      write (unit = p_nfprt, FMT = '(7F13.6)')  SiSl(1:var%kMaxInp + 1)
      write (unit = p_nfprt, FMT = '(A)') 'SiSl-2:'
      write (unit = p_nfprt, FMT = '(7F13.6)')  SiSl(var%kMaxInp + 2:2 * var%kMaxInp + 2)

      print*, ' VCOORD(:,1,1)=', gfsheadv%VCOORD(:, 1, 1)
      print *, "-------------------"

      print*, ' VCOORD(:,2,1)=', gfsheadv%VCOORD(:, 2, 1)
      print *, "-------------------"

      print*, ' VCOORD(:,3,1)=', gfsheadv%VCOORD(:, 3, 1)
      print *, "-------------------"

      print*, ' VCOORD(:,1,2)=', gfsheadv%VCOORD(:, 1, 2)
      print *, "-------------------"

      print*, ' VCOORD(:,2,2)=', gfsheadv%VCOORD(:, 2, 2)
      print *, "-------------------"

      print*, ' VCOORD(:,3,2)=', gfsheadv%VCOORD(:, 3, 2)
      print *, "-------------------"

      !  print *,' RECNAME  =',gfsheadv%RECNAME  
      !  print *, "-------------------"

      !  print *,' RECLEVTYP=',gfsheadv%RECLEVTYP
      !  print *, "-------------------"
    end if

    IDSLInp = gfshead%idsl
    IDVCInp = gfshead%idvc

    if (IDVCInp == 0_p_i8 .or. IDVCInp == 1_p_i8) then
      ! Sigma Interfaces   (kmax+1) - SigIInp
      ! Sigma Layers       (kmax  ) - SigLInp
      SigIInp(1:var%kMaxInp + 1) = SiSl(1:var%kMaxInp + 1)
      SigLInp(1:var%kMaxInp) = SiSl(var%kMaxInp + 2:2 * var%kMaxInp + 1)
    else if (IDVCInp == 2_p_i8) then
      ! Hybrid Interface A (kmax+1) - SigIInp
      ! Hybrid Interface B (kmax+1) - SigLInp
      SigIInp(1:var%kMaxInp + 1) = SiSl(1:var%kMaxInp + 1) / 100.0_p_r4  ! conversion from Pa to mbar
      SigLInp(1:var%kMaxInp + 1) = SiSl(var%kMaxInp + 2:2 * var%kMaxInp + 2)
    else
      !      stop ' ** (Error) **'
      call fatalError(headerMsg, 'Invalid IDVCInp. Aborting Chopping!')
      return
    end if
    if(isMasterProc())then
      write (unit = p_nfprt, FMT = *) 'IDSLInp ', IDSLInp
      write (unit = p_nfprt, FMT = *) 'IDVCINP ', IDVCINP
      if(IDVCInp == 2_p_i8)then
        write (unit = p_nfprt, FMT = '(/,A)')   'a_hybr  (in Pa):'
        write (unit = p_nfprt, FMT = *) SigIInp * 100.0_p_r4
        write (unit = p_nfprt, FMT = '(/,A)')   'b_hybr  '
        write (unit = p_nfprt, FMT = *) SigLInp
      else
        write (unit = p_nfprt, FMT = '(/,A)')   'SigIInp:'
        write (unit = p_nfprt, FMT = *) SigIInp
        write (unit = p_nfprt, FMT = '(/,A)')   'SigLInp:'
        write (unit = p_nfprt, FMT = *) SigLInp
      end if
    end if

    DateCurrent = DateInitial
    ! Forecast Day      - ForecastDay
    ! Time of Day       - TimeOfDay
    ! Initial Date      - DateInitial
    ! Current Date      - DateCurrent

    ! 
    !---read out data from nemsio file
    ! 
    call c_nemsio_getfilehead(gfile, dimx = im, dimy = jm, nsoil = nsoil, ntrac = ntrac, iret = ierr)
    if((isMasterProc()) .and. (ierr /= 0)) then
      call fatalError(headerMsg, 'cannot get dimension from gfile')
      call msgInLineFormatOut(headerMsg // 'Error code:', '(A)')
      call msgInLineFormatOut(ierr, '(I3)')
      call msgNewLine()
    endif
    if (ierr /= 0) return

    fieldsize = im * jm
    if(im * jm/=iMaxInp * jmaxInp) then
      print*, iMaxInp, jmaxInp, im, jm
      call fatalError(headerMsg, 'ERROR: dimension not match')
      call msgInLineFormatOut(headerMsg // 'Dimensions iMaxInp, jmaxInp, im, jm', '(A)')
      call msgInLineFormatOut(iMaxInp, '(I8)')
      call msgInLineFormatOut(jmaxInp, '(I8)')
      call msgInLineFormatOut(im, '(I8)')
      call msgInLineFormatOut(jm, '(I8)')
      call msgNewLine()
      return
    endif

    allocate(buff(fieldsize))
    call MPI_BARRIER(mpiCommGroup, ierr)
    nelements = size(buff)


    ! Spectral Coefficients of Orography (m)
    if(isMasterProc())then
      call nemsio_readrecv(gfile, 'hgt', 'sfc', 1, buff, iret = ierr)
    end if
    call MPI_BCAST(buff, nelements, MPI_REAL, mpiMasterProc, mpiCommGroup, ierr)
    call MPI_BARRIER(mpiCommGroup, ierr)
    if(isMasterProc()) then
      write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, buff(1), &
        minval(buff), &
        maxval(buff)
    endif
    if(minval(buff) < -1000.0 .and.maxval(buff) > 10000.0)then
      call fatalError(headerMsg, "Values out of range, aborting Chopping!")
      !      stop
      return
    end if
    ij = 0
    do j = 1, jmaxInp
      do i = 1, iMaxInp
        ij = ij + 1
        if(buff(ij)<0.5)buff(ij) = 0.0
        gWorkprInp(i, j) = buff(ij)
      end do
    end do
    do jb = 1, jbMax
      do ib = 1, ibMaxPerJB(jb)
        i = iPerIJB(ib, jb)
        j = jPerIJB(ib, jb)
        gTopoInp(ib, jb) = gWorkprInp(i, j)
      end do
    end do

    ! Spectral coefficients of ln(Ps) (ln(hPa)/10)

    if(isMasterProc())then
      call nemsio_readrecv(gfile, 'pres', 'sfc', 1, buff(:), iret = iret)
    end if
    call MPI_BCAST(buff, nelements, MPI_REAL, mpiMasterProc, mpiCommGroup, ierr)
    call MPI_BARRIER(mpiCommGroup, ierr)
    if(isMasterProc()) then
      write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, buff(1), &
        minval(buff), &
        maxval(buff)
    endif
    if(minval(buff) < 0.0 .and.maxval(buff) > 200000.0)then
      call fatalError(headerMsg, "Values out of range, aborting Chopping!")
      !      stop
      return
    end if
    ij = 0
    do j = 1, jmaxInp
      do i = 1, iMaxInp
        ij = ij + 1
        gWorkprInp(i, j) = buff(ij)
      end do
    end do
    do jb = 1, jbMax
      do ib = 1, ibMaxPerJB(jb)
        i = iPerIJB(ib, jb)
        j = jPerIJB(ib, jb)
        gPsfcInp(ib, jb) = gWorkprInp(i, j) / 100.0_p_r8!convert Pa to hPa
        gLnpsInp(ib, jb) = log(gPsfcInp(ib, jb) / 10.0_p_r8)!convert hPa to log(cPa)
      end do
    end do

    ! Spectral Coefficients of Virtual Temp (K)

    do k = 1, var%kMaxInp
      if(isMasterProc())then
        call nemsio_readrecv(gfile, 'tmp', 'mid layer', k, buff(:), iret = iret)
      end if
      call MPI_BCAST(buff, nelements, MPI_REAL, mpiMasterProc, mpiCommGroup, ierr)
      call MPI_BARRIER(mpiCommGroup, ierr)
      if(isMasterProc()) then
        write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, buff(1), &
          minval(buff), &
          maxval(buff)
      endif
      if(minval(buff) < 0.0 .and.maxval(buff) > 400.0)then
        call fatalError(headerMsg, "Values out of range, aborting Chopping!")
        !        stop
        return
      end if
      ij = 0
      do j = 1, jmaxInp
        do i = 1, iMaxInp
          ij = ij + 1
          gWorkprInp(i, j) = buff(ij)
        end do
      end do
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb)
          gTvirInp(ib, k, jb) = gWorkprInp(i, j)
        end do
      end do
    end do


    ! Spectral Coefficients of zonal wind (m/seg)

    do k = 1, var%kMaxInp
      if(isMasterProc())then
        call nemsio_readrecv(gfile, 'ugrd', 'mid layer', k, buff(:), iret = iret)
      endif
      call MPI_BCAST(buff, nelements, MPI_REAL, mpiMasterProc, mpiCommGroup, ierr)
      call MPI_BARRIER(mpiCommGroup, ierr)
      if(isMasterProc()) then
        write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, buff(1), &
          minval(buff), &
          maxval(buff)
      endif
      if(minval(buff) < -200.0 .and.maxval(buff) > 200.0)then
        call fatalError(headerMsg, "Values out of range, aborting Chopping!")
        !        stop
        return
      end if
      ij = 0
      do j = 1, jmaxInp
        do i = 1, iMaxInp
          ij = ij + 1
          gWorkprInp(i, j) = buff(ij)
        end do
      end do
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb)
          gUvelInp(ib, k, jb) = gWorkprInp(i, j)
        end do
      end do
    end do


    ! Spectral Coefficients of meridional wind (m/seg)

    do k = 1, var%kMaxInp
      if(isMasterProc())then
        call nemsio_readrecv(gfile, 'vgrd', 'mid layer', k, buff(:), iret = iret)
      endif
      call MPI_BCAST(buff, nelements, MPI_REAL, mpiMasterProc, mpiCommGroup, ierr)
      call MPI_BARRIER(mpiCommGroup, ierr)
      if(isMasterProc()) then
        write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, buff(1), &
          minval(buff), &
          maxval(buff)
      endif
      if(minval(buff) < -200.0 .and.maxval(buff) > 200.0)then
        call fatalError(headerMsg, "Values out of range, aborting Chopping!")
        !        stop
        return
      end if
      ij = 0
      do j = 1, jmaxInp
        do i = 1, iMaxInp
          ij = ij + 1
          gWorkprInp(i, j) = buff(ij)
        end do
      end do
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb)
          gVvelInp(ib, k, jb) = gWorkprInp(i, j)
        end do
      end do
    end do


    ! Spectral Coefficients of Specific Humidity (g/g)

    do k = 1, var%kMaxInp
      if(isMasterProc())then
        call nemsio_readrecv(gfile, 'spfh', 'mid layer', k, buff(:), iret = iret)
      endif

      call MPI_BCAST(buff, nelements, MPI_REAL, mpiMasterProc, mpiCommGroup, ierr)
      call MPI_BARRIER(mpiCommGroup, ierr)
      if(isMasterProc()) then
        write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, buff(1), &
          minval(buff), &
          maxval(buff)
      endif
      if(minval(buff) < -1.0e-2 .and.maxval(buff) > 0.4)then
        call fatalError(headerMsg, "Values out of range, aborting Chopping!")
        !        stop
        return
      end if
      ij = 0
      do j = 1, jmaxInp
        do i = 1, iMaxInp
          ij = ij + 1
          gWorkprInp(i, j) = max(buff(ij), 1e-21)
        end do
      end do
      do jb = 1, jbMax
        do ib = 1, ibMaxPerJB(jb)
          i = iPerIJB(ib, jb)
          j = jPerIJB(ib, jb)
          gSpHuInp(ib, k, jb, 1) = gWorkprInp(i, j)
        end do
      end do
    end do

    ! Spectral Coefficients of ozonio (g/g)

    if (var%getOzone) then
      NTracers = NTracers + 1
      do k = 1, var%kMaxInp
        if(isMasterProc())then
          call nemsio_readrecv(gfile, 'o3mr', 'mid layer', k, buff(:), iret = iret)
        endif
        call MPI_BCAST(buff, nelements, MPI_REAL, mpiMasterProc, mpiCommGroup, ierr)
        call MPI_BARRIER(mpiCommGroup, ierr)
        if(isMasterProc()) then
          write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, buff(1), &
            minval(buff), &
            maxval(buff)
        endif
        if(minval(buff) < -1.0e-2 .and.maxval(buff) > 0.4)then
          call fatalError(headerMsg, "Values out of range, aborting Chopping!")
          !          stop
          return
        end if
        ij = 0
        do j = 1, jmaxInp
          do i = 1, iMaxInp
            ij = ij + 1
            gWorkprInp(i, j) = max(buff(ij), 1e-21)
          end do
        end do
        do jb = 1, jbMax
          do ib = 1, ibMaxPerJB(jb)
            i = iPerIJB(ib, jb)
            j = jPerIJB(ib, jb)
            gSpHuInp(ib, k, jb, 2) = gWorkprInp(i, j)
          end do
        end do
      end do
      if (var%getTracers) then
        ! Spectral Coefficients of liquid cloud water (g/g)
        NTracers = NTracers + 1
        do k = 1, var%kMaxInp
          if(isMasterProc())then
            call nemsio_readrecv(gfile, 'clwmr', 'mid layer', k, buff(:), iret = iret)
          endif
          call MPI_BCAST(buff, nelements, MPI_REAL, mpiMasterProc, mpiCommGroup, ierr)
          call MPI_BARRIER(mpiCommGroup, ierr)
          if(isMasterProc()) then
            write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, buff(1), &
              minval(buff), &
              maxval(buff)
          endif
          if(minval(buff) < -1.0e-2 .and.maxval(buff) > 0.4)then
            call fatalError(headerMsg, "Values out of range, aborting Chopping!")
            !            stop
            return
          end if
          ij = 0
          do j = 1, jmaxInp
            do i = 1, iMaxInp
              ij = ij + 1
              gWorkprInp(i, j) = max(buff(ij), 1e-21)
            end do
          end do
          do jb = 1, jbMax
            do ib = 1, ibMaxPerJB(jb)
              i = iPerIJB(ib, jb)
              j = jPerIJB(ib, jb)
              gSpHuInp(ib, k, jb, 3) = gWorkprInp(i, j)
            end do
          end do
        end do
      end if
    end if
    if(isMasterProc())then
      print*, 'getTracers,NTracers', var%getTracers, NTracers
    end if

    if (var%grADS) then
      write (Tdef(1:2), '(I2.2)') DateCurrent(1)
      write (Tdef(4:5), '(I2.2)') DateCurrent(3)
      write (Tdef(6:8), '(A3)')   MonChar(DateCurrent(2))
      write (Tdef(9:12), '(I4.4)') DateCurrent(4)
      if(.not. GetGrADSInp_GDAS() ) return
      if (var%grADSOnly) then
        isExecOk = .true.
        return
      endif
    end if

    call ICDecompositionInput()

    if (.not. ICWriteGDAS() ) return

    !Fecha o arquivo nemsio
    call c_nemsio_close(gfile, iret = iret)

    !Finaliza
    call c_nemsio_finalize()

    deallocate(buff)

    !close(unit = nfnmc)

    if (var%gdasonly) then
      isExecOk = .true.
      return
    endif

    isExecOk = .true.

  end function GDAStoGANL2


  function GDAStoGANL() result(isExecOk)
    !# Converts from GDAS to GANL
    !# ---
    !# @info
    !# **Brief:** Converts from GDAS to GANL. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 <br>
    !# @endin
    implicit none
    logical :: isExecOk

    character (len = 1), dimension (32) :: Descriptor
    real (kind = p_r4), dimension (2 * 100 + 1) :: SiSl
    real (kind = p_r4), dimension (44) :: Extra
    integer :: nfnmc, nfozw, nftrw, nfcpt

    isExecOk = .false.

    if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A)') ' Getting GANL from GDAS NCEP File'
    nfnmc = openFile(trim(var%dGDInp) // trim(var%gdasInp), 'unformatted', 'sequential', -1, 'read', 'old')
    if(nfnmc < 0) return
    nfcpt = openFile(trim(var%dirInp) // trim(DataInp), 'unformatted', 'sequential', -1, 'write', 'replace')
    if(nfcpt < 0) return

    ! Descriptor Label (See NMC Office Note 85) (Not Used at CPTEC)
    read (unit = nfnmc) Descriptor

    ! Forecast Time               - TimeOfDay
    ! Initial Date                - DateInitial(1:4)
    !  o Initial Hour             - DateInitial(1)
    !  o Initial Month            - DateInitial(2)
    !  o Initial Day              - DateInitial(3)
    !  o Initial Year             - DateInitial(4)
    ! Sigma Interfaces and Layers - SiSl(1:201)
    ! Extra Information           - Extra(1:44)
    !  o ID Sigma Structure       - Extra(18)
    !    = 1 Phillips
    !    = 2 Mean
    !  o ID Vertical Coordinate   - Extra(19)
    !    = 1 Sigma (0 for old files)
    !    = 2 Sigma-p
    read (unit = nfnmc) TimeOfDay, DateInitial, SiSl, Extra

    if(isMasterProc())then
      write (unit = p_nfprt, FMT = '(A,1X,F15.4)') 'TimeOfDay:', TimeOfDay
      write (unit = p_nfprt, FMT = '(A,1X,4I5)') 'DateInitial:', DateInitial
      write (unit = p_nfprt, FMT = '(A)') 'SiSl:'
      write (unit = p_nfprt, FMT = '(7F13.6)')  SiSl
      write (unit = p_nfprt, FMT = '(A)') 'Extra:'
      write (unit = p_nfprt, FMT = '(7F13.6)')  Extra
    end if

    IDSLInp = int(Extra(18), p_i8)
    IDVCInp = int(Extra(19), p_i8)

    if (IDVCInp == 0_p_i8 .or. IDVCInp == 1_p_i8) then
      ! Sigma Interfaces   (kmax+1) - SigIInp
      ! Sigma Layers       (kmax  ) - SigLInp
      SigIInp(1:var%kMaxInp + 1) = SiSl(1:var%kMaxInp + 1)
      SigLInp(1:var%kMaxInp) = SiSl(var%kMaxInp + 2:2 * var%kMaxInp + 1)
    else if (IDVCInp == 2_p_i8) then
      ! Hybrid Interface A (kmax+1) - SigIInp
      ! Hybrid Interface B (kmax+1) - SigLInp
      SigIInp(1:var%kMaxInp + 1) = SiSl(1:var%kMaxInp + 1) / 100.0_p_r4  ! conversion from Pa to mbar
      SigLInp(1:var%kMaxInp + 1) = SiSl(var%kMaxInp + 2:2 * var%kMaxInp + 2)
    else
      call fatalError(headerMsg, "Invalid IDVCInp ! Aborting Chopping")
      return
    end if

    write (unit = p_nfprt, FMT = *) 'IDSLInp ', IDSLInp
    write (unit = p_nfprt, FMT = *) 'IDVCINP ', IDVCINP
    if(isMasterProc().and.IDVCInp == 2_p_i8)then
      write (unit = p_nfprt, FMT = '(/,A)')   'a_hybr  (in Pa):'
      write (unit = p_nfprt, FMT = *) SigIInp * 100.0_p_r4
      write (unit = p_nfprt, FMT = '(/,A)')   'b_hybr  '
      write (unit = p_nfprt, FMT = *) SigLInp
    else
      write (unit = p_nfprt, FMT = '(/,A)')   'SigIInp:'
      write (unit = p_nfprt, FMT = *) SigIInp
      write (unit = p_nfprt, FMT = '(/,A)')   'SigLInp:'
      write (unit = p_nfprt, FMT = *) SigLInp
    end if


    !  IF (ANY(SigIInp < 0.0_r8 .OR. SigIInp > 1.0_r8)) THEN
    !    WRITE (UNIT=p_nferr, FMT='(/,A)') ' SigI and SIgLi will be recalculated based on DelSInp'
    !    INQUIRE (FILE=TRIM(var%dirSig)//TRIM(DataSigInp), EXIST=GetNewSig)
    !    IF (GetNewSig) THEN
    !      WRITE (UNIT=p_nfprt, FMT='(/,A)') ' Getting New Delta Sigma'
    !      OPEN  (UNIT=nfsig, FILE=TRIM(var%dirSig)//TRIM(DataSigInp), FORM='FORMATTED', &
    !            ACCESS='SEQUENTIAL', ACTION='READ', STATUS='OLD', IOSTAT=ios)
    !      IF (ios /= 0) THEN
    !        WRITE (UNIT=p_nferr, FMT='(3A,I4)') ' ** (Error) ** Open file ', &
    !                                          TRIM(TRIM(var%dirSig)//TRIM(DataSigInp)), &
    !                                          ' returned IOStat = ', ios
    !        STOP ' ** (Error) **'
    !      END IF
    !      READ  (UNIT=nfsig, FMT='(5F9.6)') DelSInp
    !      CLOSE (UNIT=nfsig)
    !      CALL SigmaInp
    !    ELSE
    !      WRITE (UNIT=p_nferr, FMT='(A)') ' There is no file : '//TRIM(var%dirSig)//TRIM(DataSigInp)
    !    END IF
    !  END IF

    DateCurrent = DateInitial
    ! Forecast Day      - ForecastDay
    ! Time of Day       - TimeOfDay
    ! Initial Date      - DateInitial
    ! Current Date      - DateCurrent
    write (unit = nfcpt) ForecastDay, TimeOfDay, DateInitial, &
      DateCurrent, SigIInp, SigLInp, IDVCInp, IDSLInp


    ! Spectral Coefficients of Orography (m)
    read  (unit = nfnmc) qWorkInp
    write (unit = nfcpt) qWorkInp

    ! Spectral coefficients of ln(Ps) (ln(hPa)/10)
    read  (unit = nfnmc) qWorkInp
    write (unit = nfcpt) qWorkInp

    ! Spectral Coefficients of Virtual Temp (K)
    do k = 1, var%kMaxInp
      read  (unit = nfnmc) qWorkInp
      write (unit = nfcpt) qWorkInp
    end do

    ! Spectral Coefficients of Divergence and Vorticity (1/seg)

    if(TRIM(var%StrFormat) == 'old') then
      !
      !Spectral Coefficients of Divergence and Vorticity
      !
      do k = 1, var%kMaxInp
         !Divergence
         read  (unit = nfnmc) qWorkInp
         write (unit = nfcpt) qWorkInp
         !Vorticity
         read  (unit = nfnmc) qWorkInp
         write (unit = nfcpt) qWorkInp
      end do
    else if(TRIM(var%StrFormat) == 'new') then
    !
    !Spectral Coefficients of Divergence
    !
      do k = 1, var%kMaxInp
         !Divergence
         read  (unit = nfnmc) qWorkInp
         qWorkInp3D(:,k,1)=qWorkInp
         !Vorticity
         read  (unit = nfnmc) qWorkInp
         qWorkInp3D(:,k,2)=qWorkInp
      end do

      do k = 1, var%kMaxInp
         !Divergence
         qWorkInp=qWorkInp3D(:,k,1)
         write (unit = nfcpt) qWorkInp
      end do

      do k = 1, var%kMaxInp
         !Vorticity
         qWorkInp = qWorkInp3D(:,k,2)
         write (unit = nfcpt) qWorkInp
      end do
    else
      call fatalError(headerMsg, 'Invalid StrFormat. see at the namelist, Aborting Chopping!')
    endif
 
    ! Spectral Coefficients of Specific Humidity (g/g)
    do k = 1, var%kMaxInp
      read  (unit = nfnmc) qWorkInp
      write (unit = nfcpt) qWorkInp
    end do

    close(unit = nfcpt)

    if (var%getOzone) then

      ! Spectral Coefficients of Ozone (?)
      nfozw = openFile(trim(var%dirInp) // trim(OzonInp), 'unformatted', 'sequential', -1, 'write', 'replace')
      if(nfozw < 0) return

      do k = 1, var%kMaxInp
        read  (unit = nfnmc) qWorkInp
        write (unit = nfozw) qWorkInp
      end do
      close(unit = nfozw)
      NTracers = NTracers + 1

      if (var%getTracers) then
        ! Spectral Coefficients of Tracers (?)
        nftrw = openFile(trim(var%dirInp) // trim(TracInp), 'unformatted', 'sequential', -1, 'write', 'replace')
        if(nftrw < 0) return

        ios = 0
        Tracer : do
          do k = 1, var%kMaxInp
            read (unit = nfnmc, iostat = ios) qWorkInp
            if (ios /= 0) then
              if (ios == -1) then
                write (unit = p_nfprt, FMT = '(/,A,I5,A)') ' end of file Found - NTracers = ', &
                  NTracers, '  in:'
              else
                write (unit = p_nfprt, FMT = '(/,A,I5,A)') ' Reading Error - ios = ', ios, '  in:'
              end if
              write (unit = p_nfprt, FMT = '(1X,A,/)') trim(var%dGDInp) // trim(var%gdasInp)
              EXIT Tracer
            end if
            write (unit = nftrw) qWorkInp
          end do
          NTracers = NTracers + 1
        end do Tracer
        close(unit = nftrw)
      end if

    end if

    close(unit = nfnmc)

    if (var%gdasonly) then
      isExecOk = .true.
      print *, ' gdasOnly = .true. '
      return
    end if
    !      stop ' gdasOnly = .true. '
    isExecOk = .true.    

  end function GDAStoGANL


  function ICRead_and_Chop(GetNewTop) result(isExecOk)
    !# Chopps Input CPTEC No Topo-Smoothed Analysis file
    !# ---
    !# @info
    !# **Brief:** Chopps Input CPTEC No Topo-Smoothed Analysis file. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 <br>
    !# @endin
    implicit none
    logical, intent(inout) :: GetNewTop
    logical :: isExecOk
    
    integer :: mm, nn, i1, i2, k, m, mw
    integer :: nficr, nfozr, nftop, nftrr

    isExecOk = .false.
    
    ! REAL(KIND=r8) ::gTopoOutGaus8(iMaxout,jMaxout)
    nficr = openFile(trim(var%dirInp) // trim(DataInp), 'unformatted', 'sequential', -1, 'read', 'old') 
    if(nficr < 0) return

    mw = min(var%mEndOut + 1, var%mEndInp + 1)
    if (var%gdasonly) then
      ! Old file with sigma level
      read (unit = nficr) ForecastDay, TimeOfDay, DateInitial(1:4), &
        DateCurrent(1:4), SigIInp(1:KmaxInpp), SigLInp(1:var%kMaxInp)
    else

      read (unit = nficr) ForecastDay, TimeOfDay, DateInitial(1:4), &
        DateCurrent(1:4), SigIInp(1:KmaxInpp), SigLInp(1:KmaxInpp), IDVCInp, IDSLInp
    end if
    SigInterInp = SigIInp
    SigLayerInp = SigLInp

    if(isMasterProc())  write (unit = p_nfprt, FMT = '(/,A,I5,A,F15.4)') ' ForecastDay = ', ForecastDay, &
      ' TimeOfDay = ', TimeOfDay
    if(isMasterProc())write (unit = p_nfprt, FMT = '(/,A,4I5)') ' DateInitial = ', DateInitial
    if(isMasterProc())write (unit = p_nfprt, FMT = '(/,A,4I5)') ' DateCurrent = ', DateCurrent

    if (var%gdasonly) then
      do k = 1, var%kMaxInp
        DelSigmaInp(k) = SigInterInp(k) - SigInterInp(k + 1)
      end do
      if(isMasterProc())  write (unit = p_nfprt, FMT = '(/,A)')  ' DelSigmaInp:'
      if(isMasterProc())  write (unit = p_nfprt, FMT = '(7F10.6)') DelSigmaInp
      if(isMasterProc())  write (unit = p_nfprt, FMT = '(/,A)')  ' SigInterInp:'
      if(isMasterProc())  write (unit = p_nfprt, FMT = '(7F10.6)') SigInterInp
      if(isMasterProc())  write (unit = p_nfprt, FMT = '(/,A)')  ' SigLayerInp:'
      if(isMasterProc())  write (unit = p_nfprt, FMT = '(7F10.6)') SigLayerInp
    else
      if(isMasterProc())  write (unit = p_nfprt, FMT = '(/,A)')  ' a_hybr:'
      if(isMasterProc())  write (unit = p_nfprt, FMT = '(7F12.6)') SigInterInp
      if(isMasterProc())  write (unit = p_nfprt, FMT = '(/,A)')  ' b_hybr:'
      if(isMasterProc())  write (unit = p_nfprt, FMT = '(7F12.6)') SigLayerInp
    endif

    read (unit = nficr) qWorkInp
    do mm = 1, mymmax
      m = msinproc(mm, myid_four)
      i2 = 2 * mymnmap(mm, m) - 1
      if (m.gt.var%mEndInp + 1) then
        do nn = 0, 2 * (var%mEndOut + 1 - m) + 1
          qTopoInp(i2 + nn) = 0.0_p_r8
        enddo
      else
        i1 = 2 * mnmap(m, m) - 1
        do nn = 0, 2 * (mw - m) + 1
          qTopoInp(i2 + nn) = qWorkInp(i1 + nn)
        enddo
        do nn = 2 * (mw - m) + 2, 2 * (var%mEndOut + 1 - m) + 1
          qTopoInp(i2 + nn) = 0.0_p_r8
        enddo
      endif
    enddo
    if(isMasterProc()) then
      write (unit = p_nfprt, FMT = '(/,A)') ' TopoInp:'
      write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') 0, qWorkInp(1), &
        minval(qWorkInp(2:)), &
        maxval(qWorkInp(2:))
    endif
    !-----------------------------------------------------------------------------
    !    spectral topography (m)
    !-----------------------------------------------------------------------------
    qTopoOut = qTopoInp
    qTopoOutSpec = qTopoInp
    gTopoDel = 0.0_p_r8
    inquire (file = trim(var%dirtop) // trim(DataTop), EXIST = GetNewTop)
    if (GetNewTop.and.havesurf) then
      write (unit = p_nfprt, FMT = '(/,A)')' Getting New Topography'
      nftop = openFile(trim(var%dirtop) // trim(DataTop), 'unformatted', 'sequential', -1, 'read', 'old') 
      if(nftop < 0) return

      read  (unit = nftop) qWorkprOut
      close (unit = nftop)
      !
      !i1 = 1
      !DO m=1,mEndOut+1
      !   i2 = m
      !   DO nn=mEndOut+1,m,-1
      !      qtorto(2*i2-1) = qWorkprOut(2*i1-1)
      !      qtorto(2*i2  ) = qWorkprOut(2*i1  )
      !      i1 = i1+1
      !      i2 = i2 + nn
      !   ENDDO
      !ENDDO
      !qWorkprOut = qtorto
      !
      do mm = 1, mymmax
        m = msinproc(mm, myid_four)
        i2 = 2 * mymnmap(mm, m) - 1
        if (m.GT.var%mEndInp + 1) then
          do nn = 0, 2 * (var%mEndOut + 1 - m) + 1
            qTopoOut(i2 + nn) = 0.0_p_r8
            qTopoOutSpec(i2 + nn) = 0.0_p_r8
          enddo
        else
          i1 = 2 * mnmap_out(m, m) - 1
          if (isMasterProc()) write(*, *) ' reading topo i1 i2 ', i1, i2
          do nn = 0, 2 * (mw - m) + 1
            qTopoOut(i2 + nn) = qWorkprOut(i1 + nn)
            qTopoOutSpec(i2 + nn) = qWorkprOut(i1 + nn)
          enddo
          if (isMasterProc()) then
            write(98, *) ' qtopoinp for m', m
            write(98, *) (qtopoinp(i2 + nn), nn = 0, 2 * (mw - m) - 1)
            write(99, *) ' qtopoout for m', m
            write(99, *) (qtopoout(i2 + nn), nn = 0, 2 * (mw - m) - 1)
          endif
          do nn = 2 * (mw - m) + 2, 2 * (var%mEndOut + 1 - m) + 1
            qTopoOut(i2 + nn) = 0.0_p_r8
            qTopoOutSpec(i2 + nn) = 0.0_p_r8
          enddo
        endif
      enddo

    else
      if (var%smoothTopo.and.havesurf) then
        if(isMasterProc()) then
          write (unit = p_nfprt, FMT = '(/,A)') ' Chopping Old Topography for Smoothing'
        endif
        qTopoOut = qTopoInp
      end if
    end if
    if(isMasterProc()) then
      write (unit = p_nfprt, FMT = *) ' Chopping Old Topography for Smoothing', 'smoothTopo=', var%smoothTopo
    endif
    if (var%smoothTopo.and.havesurf) then
      call SmoothCoef()
    end if
    !-----------------------------------------------------------------------------
    !    LN Surface Pressure (cb)
    !-----------------------------------------------------------------------------
    read (unit = nficr) qWorkInp
    do mm = 1, mymmax
      m = msinproc(mm, myid_four)
      i2 = 2 * mymnmap(mm, m) - 1
      if (m.gt.var%mEndInp + 1) then
        do nn = 0, 2 * (var%mEndOut + 1 - m) + 1
          qLnPsInp(i2 + nn) = 0.0_p_r8
        enddo
      else
        i1 = 2 * mnmap(m, m) - 1
        do nn = 0, 2 * (mw - m) + 1
          qLnPsInp(i2 + nn) = qWorkInp(i1 + nn)
        enddo
        do nn = 2 * (mw - m) + 2, 2 * (var%mEndOut + 1 - m) + 1
          qLnPsInp(i2 + nn) = 0.0_p_r8
        enddo
      endif
    enddo
    !  IF (var%smoothTopo.and.havesurf) THEN
    !    CALL SmoothCoefAtm(qLnPsInp (1:Mnwv2Inp))
    !  END IF

    if(isMasterProc()) then
      write (unit = p_nfprt, FMT = '(/,A)') ' LnPsInp:'
      write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') 0, qWorkInp(1), &
        minval(qWorkInp(2:)), &
        maxval(qWorkInp(2:))
    endif
    !-----------------------------------------------------------------------------
    !    VIRTUAL TEMPERATURE (k)
    !-----------------------------------------------------------------------------
    if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A)') ' TvirInp:'
    do k = 1, var%kMaxInp
      read (unit = nficr) qWorkInp
      if (k.ge.myfirstlev.and.k.le.mylastlev) then
        do mm = 1, mymmax
          m = msinproc(mm, myid_four)
          i2 = 2 * mymnmap(mm, m) - 1
          if (m.gt.var%mEndInp + 1) then
            do nn = 0, 2 * (var%mEndOut + 1 - m) + 1
              qTvirInp(i2 + nn, k + 1 - myfirstlev) = 0.0_p_r8
            enddo
          else
            i1 = 2 * mnmap(m, m) - 1
            do nn = 0, 2 * (mw - m) + 1
              qTvirInp(i2 + nn, k + 1 - myfirstlev) = qWorkInp(i1 + nn)
            enddo
            do nn = 2 * (mw - m) + 2, 2 * (var%mEndOut + 1 - m) + 1
              qTvirInp(i2 + nn, k + 1 - myfirstlev) = 0.0_p_r8
            enddo
          endif
        enddo
        !       IF (var%smoothTopo.and.havesurf) THEN
        !          CALL SmoothCoefAtm(qTvirInp (1:Mnwv2Inp,k+1-myfirstlev))
        !       END IF
      end if

      if(isMasterProc()) then
        write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, qWorkInp(1), &
          minval(qWorkInp(2:)), &
          maxval(qWorkInp(2:))
      endif
    end do
    !-----------------------------------------------------------------------------
    !    DIVERGENCY AND VORTICITY (1/s)
    !-----------------------------------------------------------------------------
    if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A)') ' DivgInp - VortInp:'
    if(TRIM(var%StrFormat) == 'old') then
      !
      !Spectral Coefficients divergency and vorticity (1/s)
      !
      do k = 1, var%kMaxInp
         !
         ! Spectral Coefficients Divergency
         !
        read (unit = nficr) qWorkInp
        if (k.ge.myfirstlev.and.k.le.mylastlev) then
          do mm = 1, mymmax
            m = msinproc(mm, myid_four)
            i2 = 2 * mymnmap(mm, m) - 1
            if (m.gt.var%mEndInp + 1) then
              do nn = 0, 2 * (var%mEndOut + 1 - m) + 1
                qDivgInp(i2 + nn, k + 1 - myfirstlev) = 0.0_p_r8
              enddo
            else
              i1 = 2 * mnmap(m, m) - 1
              do nn = 0, 2 * (mw - m) + 1
                qDivgInp(i2 + nn, k + 1 - myfirstlev) = qWorkInp(i1 + nn)
              enddo
              do nn = 2 * (mw - m) + 2, 2 * (var%mEndOut + 1 - m) + 1
                qDivgInp(i2 + nn, k + 1 - myfirstlev) = 0.0_p_r8
              enddo
            endif
          enddo
          !       IF (var%smoothTopo.and.havesurf) THEN
          !          CALL SmoothCoefAtm(qDivgInp (1:Mnwv2Inp,k+1-myfirstlev))
          !       END IF
        end if

        if(isMasterProc()) then
          write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, qWorkInp(1), &
            minval(qWorkInp(2:)), &
            maxval(qWorkInp(2:))
        endif
         !
         ! Spectral Coefficients Vorticity
         !
        read (unit = nficr) qWorkInp
        if (k.ge.myfirstlev.and.k.le.mylastlev) then
          do mm = 1, mymmax
            m = msinproc(mm, myid_four)
            i2 = 2 * mymnmap(mm, m) - 1
            if (m.gt.var%mEndInp + 1) then
              do nn = 0, 2 * (var%mEndOut + 1 - m) + 1
                qVortInp(i2 + nn, k + 1 - myfirstlev) = 0.0_p_r8
              enddo
            else
              i1 = 2 * mnmap(m, m) - 1
              do nn = 0, 2 * (mw - m) + 1
                qVortInp(i2 + nn, k + 1 - myfirstlev) = qWorkInp(i1 + nn)
              enddo
              do nn = 2 * (mw - m) + 2, 2 * (var%mEndOut + 1 - m) + 1
                qVortInp(i2 + nn, k + 1 - myfirstlev) = 0.0_p_r8
              enddo
            endif
          enddo
          !       IF (var%smoothTopo.and.havesurf) THEN
          !          CALL SmoothCoefAtm(qVortInp (1:Mnwv2Inp,k+1-myfirstlev))
          !       END IF
        end if

        if(isMasterProc()) then
          write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, qWorkInp(1), &
            minval(qWorkInp(2:)), &
            maxval(qWorkInp(2:))
        endif
      end do
    else if(TRIM(var%StrFormat) == 'new') then
      !
      !Spectral Coefficients of Divergence
      !
      do k = 1, var%kMaxInp
         !
         ! Spectral Coefficients Divergence
         !
        read (unit = nficr) qWorkInp
        if (k.ge.myfirstlev.and.k.le.mylastlev) then
          do mm = 1, mymmax
            m = msinproc(mm, myid_four)
            i2 = 2 * mymnmap(mm, m) - 1
            if (m.gt.var%mEndInp + 1) then
              do nn = 0, 2 * (var%mEndOut + 1 - m) + 1
                qDivgInp(i2 + nn, k + 1 - myfirstlev) = 0.0_p_r8
              enddo
            else
              i1 = 2 * mnmap(m, m) - 1
              do nn = 0, 2 * (mw - m) + 1
                qDivgInp(i2 + nn, k + 1 - myfirstlev) = qWorkInp(i1 + nn)
              enddo
              do nn = 2 * (mw - m) + 2, 2 * (var%mEndOut + 1 - m) + 1
                qDivgInp(i2 + nn, k + 1 - myfirstlev) = 0.0_p_r8
              enddo
            endif
          enddo
          !       IF (var%smoothTopo.and.havesurf) THEN
          !          CALL SmoothCoefAtm(qDivgInp (1:Mnwv2Inp,k+1-myfirstlev))
          !       END IF
        end if

        if(isMasterProc()) then
          write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, qWorkInp(1), &
            minval(qWorkInp(2:)), &
            maxval(qWorkInp(2:))
        endif
      end do
      !
      !Spectral Coefficients of Vorticity
      !
      do k = 1, var%kMaxInp
         !
         ! Spectral Coefficients Vorticity
         !
        read (unit = nficr) qWorkInp
        if (k.ge.myfirstlev.and.k.le.mylastlev) then
          do mm = 1, mymmax
            m = msinproc(mm, myid_four)
            i2 = 2 * mymnmap(mm, m) - 1
            if (m.gt.var%mEndInp + 1) then
              do nn = 0, 2 * (var%mEndOut + 1 - m) + 1
                qVortInp(i2 + nn, k + 1 - myfirstlev) = 0.0_p_r8
              enddo
            else
              i1 = 2 * mnmap(m, m) - 1
              do nn = 0, 2 * (mw - m) + 1
                qVortInp(i2 + nn, k + 1 - myfirstlev) = qWorkInp(i1 + nn)
              enddo
              do nn = 2 * (mw - m) + 2, 2 * (var%mEndOut + 1 - m) + 1
                qVortInp(i2 + nn, k + 1 - myfirstlev) = 0.0_p_r8
              enddo
            endif
          enddo
          !       IF (var%smoothTopo.and.havesurf) THEN
          !          CALL SmoothCoefAtm(qVortInp (1:Mnwv2Inp,k+1-myfirstlev))
          !       END IF
        end if

        if(isMasterProc()) then
          write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, qWorkInp(1), &
            minval(qWorkInp(2:)), &
            maxval(qWorkInp(2:))
        endif
      end do
    else
      call fatalError(headerMsg, 'Invalid StrFormat. see at the namelist, Aborting Chopping!')
    endif
    !-----------------------------------------------------------------------------
    !    Specific Humidity (kg/kg)
    !-----------------------------------------------------------------------------
    if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A)') ' SpHuInp:'
    do k = 1, var%kMaxInp
      read (unit = nficr) qWorkInp
      if (k.ge.myfirstlev.and.k.le.mylastlev) then
        do mm = 1, mymmax
          m = msinproc(mm, myid_four)
          i2 = 2 * mymnmap(mm, m) - 1
          if (m.gt.var%mEndInp + 1) then
            do nn = 0, 2 * (var%mEndOut + 1 - m) + 1
              qSpHuInp(i2 + nn, k + 1 - myfirstlev, 1) = 0.0_p_r8
            enddo
          else
            i1 = 2 * mnmap(m, m) - 1
            do nn = 0, 2 * (mw - m) + 1
              qSpHuInp(i2 + nn, k + 1 - myfirstlev, 1) = qWorkInp(i1 + nn)
            enddo
            do nn = 2 * (mw - m) + 2, 2 * (var%mEndOut + 1 - m) + 1
              qSpHuInp(i2 + nn, k + 1 - myfirstlev, 1) = 0.0_p_r8
            enddo
          endif
        enddo
        !       IF (var%smoothTopo.and.havesurf) THEN
        !          CALL SmoothCoefAtm(qSpHuInp (1:Mnwv2Inp,k+1-myfirstlev,1))
        !       END IF
      end if

      if(isMasterProc()) then
        write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, qWorkInp(1), &
          minval(qWorkInp(2:)), &
          maxval(qWorkInp(2:))
      endif
    end do
    if(var%gdasonly)then
      close(unit = nficr, status = 'KEEP')
    else
      close(unit = nficr, status = 'KEEP')
    end if
    !-----------------------------------------------------------------------------
    !    TRacers (kg/kg)
    !-----------------------------------------------------------------------------
    if (var%getOzone) then

      nfozr = openFile(trim(var%dirInp) // trim(OzonInp), 'unformatted', 'sequential', -1, 'read', 'old') 
      if(nfozr < 0) return

      nt = 2
      do k = 1, var%kMaxInp
        read (unit = nfozr) qWorkInp
        if (k.ge.myfirstlev.and.k.le.mylastlev) then
          do mm = 1, mymmax
            m = msinproc(mm, myid_four)
            i2 = 2 * mymnmap(mm, m) - 1
            if (m.gt.var%mEndInp + 1) then
              do nn = 0, 2 * (var%mEndOut + 1 - m) + 1
                qSpHuInp(i2 + nn, k + 1 - myfirstlev, nt) = 0.0_p_r8
              enddo
            else
              i1 = 2 * mnmap(m, m) - 1
              do nn = 0, 2 * (mw - m) + 1
                qSpHuInp(i2 + nn, k + 1 - myfirstlev, nt) = qWorkInp(i1 + nn)
              enddo
              do nn = 2 * (mw - m) + 2, 2 * (var%mEndOut + 1 - m) + 1
                qSpHuInp(i2 + nn, k + 1 - myfirstlev, nt) = 0.0_p_r8
              enddo
            endif
          enddo
          !         IF (var%smoothTopo.and.havesurf) THEN
          !            CALL SmoothCoefAtm(qSpHuInp (1:Mnwv2Inp,k+1-myfirstlev,nt))
          !         END IF
        end if

      end do
      close(unit = nfozr)

      if (var%getTracers) then
        nftrr = openFile(trim(var%dirInp) // trim(TracInp), 'unformatted', 'sequential', -1, 'read', 'old') 
        if(nftrr < 0) return

        do nt = 3, NTracers
          do k = 1, var%kMaxInp
            read (unit = nftrr) qWorkInp
            if (k.ge.myfirstlev.and.k.le.mylastlev) then
              do mm = 1, mymmax
                m = msinproc(mm, myid_four)
                i2 = 2 * mymnmap(mm, m) - 1
                if (m.gt.var%mEndInp + 1) then
                  do nn = 0, 2 * (var%mEndOut + 1 - m) + 1
                    qSpHuInp(i2 + nn, k + 1 - myfirstlev, nt) = 0.0_p_r8
                  enddo
                else
                  i1 = 2 * mnmap(m, m) - 1
                  do nn = 0, 2 * (mw - m) + 1
                    qSpHuInp(i2 + nn, k + 1 - myfirstlev, nt) = qWorkInp(i1 + nn)
                  enddo
                  do nn = 2 * (mw - m) + 2, 2 * (var%mEndOut + 1 - m) + 1
                    qSpHuInp(i2 + nn, k + 1 - myfirstlev, nt) = 0.0_p_r8
                  enddo
                endif
              enddo
              !            IF (var%smoothTopo.and.havesurf) THEN
              !               CALL SmoothCoefAtm(qSpHuInp (1:Mnwv2Inp,k+1-myfirstlev,nt))
              !            END IF
            end if
          end do
        end do
        close(unit = nftrr)
      end if

    end if
    
    isExecOk = .true.

  end function ICRead_and_Chop


  function ICWrite() result(isExecOk)
    !# Writes IC ...?
    !# ---
    !# @info
    !# **Brief:** Writes IC ...?. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 <br>
    !# @endin
    implicit none
    logical :: isExecOk
    
    integer :: nficw, nfozg, nftrg

    isExecOk = .false.
    IF (  var%gdasOnly) THEN
       ! Sigma-p case, first calculate pressure for each interface and then
         SigIOut=SigInterOut*100.0_p_r8  ! convert back from mb to Pa
         SigLOut=SigLayerOut
    ELSE
      IF ( IDVCInp == 0_p_i4 .OR. IDVCInp == 1_p_i4 ) THEN
         SigIOut(1:KmaxOutp)    = SigInterOut(1:KmaxOutp)
         SigLOut(1:var%kMaxOut) = SigLayerOut(1:var%kMaxOut)
      ELSE IF ( IDVCInp == 2_p_i4 ) THEN
       ! Sigma-p case, first calculate pressure for each interface and then
         SigIOut=SigInterOut*100.0_p_r8  ! convert back from mb to Pa
         SigLOut=SigLayerOut
      END IF
    END IF
    if (isMasterProc()) then
      nficw = openFile(getOutFileName(), 'unformatted', 'sequential', -1, 'write', 'replace')
      if(nficw < 0) return
      IF (  var%gdasOnly) THEN
         write (unit = nficw) ForecastDay, TimeOfDay, DateInitial, &
           DateCurrent, SigIOut, SigLOut
      else
         if ( IDVCInp == 0_p_i4 .OR. IDVCInp == 1_p_i4 ) then
            write (unit = nficw) ForecastDay, TimeOfDay, DateInitial, &
               DateCurrent, SigIOut(1:KmaxOutp), SigLOut(1:var%kMaxOut)
         else if ( IDVCInp == 2_p_i4 ) then
            write (unit = nficw) ForecastDay, TimeOfDay, DateInitial, &
              DateCurrent, SigIOut, SigLOut
         endif
      endif 
      write (unit = p_nfprt, FMT = '(/,A,I5,A,F15.4)') ' ForecastDay = ', ForecastDay, &
        ' TimeOfDay = ', TimeOfDay
      write (unit = p_nfprt, FMT = '(/,A,4I5)') ' DateInitial = ', DateInitial
      write (unit = p_nfprt, FMT = '(/,A,4I5)') ' DateCurrent = ', DateCurrent
      write (unit = p_nfprt, FMT = '(/,A)')  ' DelSigmaOut:'
      write (unit = p_nfprt, FMT = '(7F10.6)') DelSigmaOut(1:var%kMaxOut)
      write (unit = p_nfprt, FMT = '(/,A)')  ' SigInterOut:'
      write (unit = p_nfprt, FMT = '(7F10.6)') SigInterOut(1:KmaxOutp)
      write (unit = p_nfprt, FMT = '(/,A)')  ' SigLayerOut:'
      write (unit = p_nfprt, FMT = '(7F10.6)') SigLayerOut(1:var%kMaxOut)
    end if
    !
    !Spectral Coefficients of Topography (m)
    !
    if (VerticalInterp) then
      if (havesurf) call Collect_Spec(qTopoOut, qWorkOut(:, 1), 1, 1, 0)
    else
      if (havesurf) call Collect_Spec(qTopoInp, qWorkOut(:, 1), 1, 1, 0)
    endif

    if (isMasterProc()) then
      qWorkprOut = qWorkOut(:, 1)
      write (unit = nficw) qWorkprOut
      write (unit = p_nfprt, FMT = '(/,A)') ' TopoOut:'
      write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') 0, qWorkOut(1, 1), &
        minval(qWorkOut(2:, 1)), maxval(qWorkOut(2:, 1))
    end if
    !
    !Spectral Coefficients of ln Surface Pressure (cb)
    !
    if (VerticalInterp) then
      if (havesurf) call Collect_Spec(qLnpsOut, qWorkOut(:, 1), 1, 1, 0)
    else
      if (havesurf) call Collect_Spec(qLnpsInp, qWorkOut(:, 1), 1, 1, 0)
    endif
    if (isMasterProc()) then
      qWorkprOut = qWorkOut(:, 1)
      write(p_nfprt, *) 'qLnpsOut  ', (qWorkprOut(i), i = 1, 20)
      write (unit = nficw) qWorkprOut
      write (unit = p_nfprt, FMT = '(/,A)') ' LnPsOut:'
      write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') 0, qWorkOut(1, 1), &
        minval(qWorkOut(2:, 1)), maxval(qWorkOut(2:, 1))

      write (unit = p_nfprt, FMT = '(/,A)') ' TvirOut:'
    end if
    !
    !Spectral Coefficients of Virtual Temperature (K)
    !
    if (VerticalInterp) then
      call Collect_Spec(qTvirOut, qWorkOut, kmaxloc, var%kMaxOut, 0)
    else
      call Collect_Spec(qTvirInp, qWorkOut, kmaxloc, var%kMaxOut, 0)
    endif

    if (isMasterProc()) then
      do k = 1, var%kMaxOut
        qWorkprOut = qWorkOut(:, k)
        write (unit = nficw) qWorkprOut
        write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, qWorkOut(1, k), &
          minval(qWorkOut(2:, k)), &
          maxval(qWorkOut(2:, k))
      end do
    end if
    !
    !Spectral Coefficients of Divergence and Vorticity (1/s)
    !
    if (VerticalInterp) then
      call Collect_Spec(qDivgOut, qWorkOut, kmaxloc, var%kMaxOut, 0)
      call Collect_Spec(qVortOut, qWorkOut1, kmaxloc, var%kMaxOut, 0)
    else
      call Collect_Spec(qDivgInp, qWorkOut, kmaxloc, var%kMaxOut, 0)
      call Collect_Spec(qVortInp, qWorkOut1, kmaxloc, var%kMaxOut, 0)
    endif

    if (isMasterProc()) then
      write (unit = p_nfprt, FMT = '(/,A)') ' DivgOut - VortOut:'
      ! Get NCEP Analysis and dumping to 'old' format or 'new' format
      if(TRIM(var%StrFormat) == 'old') then
        !
        !Spectral Coefficients of Divergence and Vorticity
        !
        do k = 1, var%kMaxOut
	  !Divergence
          qWorkprOut = qWorkOut(:, k)
          write (unit = nficw) qWorkprOut
          write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, qWorkOut(1, k), &
          minval(qWorkOut(2:, k)), &
          maxval(qWorkOut(2:, k))
          !Vorticity
          qWorkprOut = qWorkOut1(:, k)  
          write (unit = nficw) qWorkprOut
          write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, qWorkOut1(1, k), &
          minval(qWorkOut1(2:, k)), &
          maxval(qWorkOut1(2:, k))
        end do
      else if(TRIM(var%StrFormat) == 'new') then
        !
        !Spectral Coefficients of Divergence
	!
        do k = 1, var%kMaxOut
	  !Divergence
          qWorkprOut = qWorkOut(:, k)
          write (unit = nficw) qWorkprOut
          write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, qWorkOut(1, k), &
          minval(qWorkOut(2:, k)), &
          maxval(qWorkOut(2:, k))
        end do
        do k = 1, var%kMaxOut
          !Vorticity
          qWorkprOut = qWorkOut1(:, k)  
          write (unit = nficw) qWorkprOut
          write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, qWorkOut1(1, k), &
          minval(qWorkOut1(2:, k)), &
          maxval(qWorkOut1(2:, k))
        end do
        !Spectral Coefficients of Vorticity
      else
        call fatalError(headerMsg, 'Invalid StrFormat. see at the namelist, Aborting Chopping!')
      endif
    end if
    !
    !Spectral Coefficients of Specific Humidity (kg/kg)
    !
    if (VerticalInterp) then
      call Collect_Spec(qSpHuOut(:, :, 1), qWorkOut, kmaxloc, var%kMaxOut, 0)
    else
      call Collect_Spec(qSpHuInp(:, :, 1), qWorkOut, kmaxloc, var%kMaxOut, 0)
    endif
    if (isMasterProc()) then
      write (unit = p_nfprt, FMT = '(/,A)') ' SpHuOut:'
      do k = 1, var%kMaxOut
        qWorkprOut = qWorkOut(:, k)
        write (unit = nficw) qWorkprOut
        write (unit = p_nfprt, FMT = '(I5,1P3G12.5)') k, qWorkOut(1, k), &
          minval(qWorkOut(2:, k)), &
          maxval(qWorkOut(2:, k))
      end do
      write (unit = p_nfprt, FMT = '(A)') ' '

      close(unit = nficw)
    end if
    !
    !Tracers (kg/kg)
    !
    if (var%getOzone) then

      nt = 2
      if (VerticalInterp) then
        call Collect_Grid_Full(gSpHuOut(:, :, :, nt), gWorkOut, var%kMaxOut, 0)
      else
        call Collect_Grid_Full(gSpHuInp(:, :, :, nt), gWorkOut, var%kMaxOut, 0)
      endif

      if (isMasterProc()) then
        inquire (IOLENGTH = IOL) gWorkprOut
        nfozg = openFile(trim(var%dirout) // trim(OzonOut), 'unformatted', 'direct', IOL, 'write', 'replace')
        if(nfozg < 0) return

        do k = 1, var%kMaxOut
          gWorkprOut = gWorkOut(:, :, k)
          write (unit = nfozg, REC = k) gWorkprOut
        end do
        close(unit = nfozg)
      end if

      if (var%getTracers) then

        if (isMasterProc()) then
          nftrg = openFile(trim(var%dirout) // trim(TracOut), 'unformatted', 'sequential', -1, 'write', 'replace')
          if(nftrg < 0) return
        end if

        do nt = 3, Ntracers
          if (VerticalInterp) then
            call Collect_Grid_Full(gSpHuOut(:, :, :, nt), gWorkOut, var%kMaxOut, 0)
          else
            call Collect_Grid_Full(gSpHuInp(:, :, :, nt), gWorkOut, var%kMaxOut, 0)
          endif
          if (isMasterProc()) then
            do k = 1, var%kMaxOut
              gWorkprOut = gWorkOut(:, :, k)
              write (unit = nftrg)gWorkprOut
            end do
          end if
        end do
        if  (isMasterProc()) close(unit = nftrg)
      end if

    end if
    
    isExecOk = .true.    

  end function ICWrite


  function ICWriteGDAS() result(isExecOk)
    !# Writes GDAS file
    !# ---
    !# @info
    !# **Brief:** Writes GDAS file. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 <br>
    !# @endin
    implicit none
    logical :: isExecOk
    
    integer :: nfozw, nftrw
    integer :: nfcpt
    
    isExecOk = .false.

    if (isMasterProc()) then
      nfcpt = openFile(trim(var%dirInp) // trim(DataInp), 'unformatted', 'sequential', -1, 'write', 'replace')
      if(nfcpt < 0) return
      write (unit = nfcpt) ForecastDay, TimeOfDay, DateInitial, &
        DateCurrent, SigIInp, SigLInp, IDVCInp, IDSLInp

      write (unit = p_nfprt, FMT = '(/,A,I5,A,F15.4)') ' ForecastDay = ', ForecastDay, &
        ' TimeOfDay = ', TimeOfDay
      write (unit = p_nfprt, FMT = '(/,A,4I5)') ' DateInitial = ', DateInitial
      write (unit = p_nfprt, FMT = '(/,A,4I5)') ' DateCurrent = ', DateCurrent
    end if

    ! Spectral Coefficients of Orography (m)

    if (havesurf) call Collect_Spec(qTopoInp, qWorkOut(:, 1), 1, 1, 0)
    if (isMasterProc()) then
      qWorkprOut = qWorkOut(:, 1)
      write (unit = nfcpt) qWorkprOut
    end if

    ! Spectral coefficients of ln(Ps) (ln(hPa)/10)

    if (havesurf) call Collect_Spec(qLnpsInp, qWorkOut(:, 1), 1, 1, 0)
    if (isMasterProc()) then
      qWorkprOut = qWorkOut(:, 1)
      write (unit = nfcpt) qWorkprOut
    end if

    ! Spectral Coefficients of Virtual Temp (K)

    call Collect_Spec(qTvirInp, qWorkInOut, kmaxloc_in, var%kMaxInp, 0)
    if (isMasterProc()) then
      do k = 1, var%kMaxInp
        qWorkprOut = qWorkInOut(:, k)
        write (unit = nfcpt) qWorkprOut
      end do
    end if

    ! Spectral Coefficients of Divergence and Vorticity (1/seg)

    call Collect_Spec(qDivgInp, qWorkInOut, kmaxloc_in, var%kMaxInp, 0)
    call Collect_Spec(qVortInp, qWorkInOut1, kmaxloc_in, var%kMaxInp, 0)
 
   ! Get NCEP Analysis and dumping to 'old' format or 'new' format
    if(TRIM(var%StrFormat) == 'old') then
      if (isMasterProc()) then
        do k = 1, var%kMaxInp
          qWorkprOut = qWorkInOut(:, k)
          write (unit = nfcpt) qWorkprOut
          qWorkprOut = qWorkInOut1(:, k)
          write (unit = nfcpt) qWorkprOut
        end do
      end if
    else if(TRIM(var%StrFormat) == 'new') then
      !Spectral Coefficients of Divergence
      if (isMasterProc()) then
        do k = 1, var%kMaxInp
          qWorkprOut = qWorkInOut(:, k)
          write (unit = nfcpt) qWorkprOut
        end do
      end if
      !Spectral Coefficients of Vorticity
      if (isMasterProc()) then
        do k = 1, var%kMaxInp
          qWorkprOut = qWorkInOut1(:, k)
          write (unit = nfcpt) qWorkprOut
        end do
      end if
    else
      call fatalError(headerMsg, 'Invalid StrFormat. see at the namelist, Aborting Chopping!')
    endif

    ! Spectral Coefficients of Specific Humidity (g/g)

    call Collect_Spec(qSpHuInp(:, :, 1), qWorkInOut, kmaxloc_in, var%kMaxInp, 0)
    if (isMasterProc()) then
      do k = 1, var%kMaxInp
        qWorkprOut = qWorkInOut(:, k)
        write (unit = nfcpt) qWorkprOut
      end do
    end if

    if (isMasterProc()) then
      close(unit = nfcpt)
    endif

    if (var%getOzone) then

      if (isMasterProc()) then
        ! Spectral Coefficients of Ozone (?)
        nfozw = openFile(trim(var%dirInp) // trim(OzonInp), 'unformatted', 'sequential', -1, 'write', 'replace')
        if(nfozw < 0) return
      end if
      nt = 2
      call Collect_Spec(qSpHuInp(:, :, nt), qWorkInOut, kmaxloc_in, var%kMaxInp, 0)
      if (isMasterProc()) then
        do k = 1, var%kMaxInp
          qWorkprOut = qWorkInOut(:, k)
          write (unit = nfozw) qWorkprOut
        end do
        close(unit = nfozw)
      end if

      if (var%getTracers) then
        if (isMasterProc()) then
          ! Spectral Coefficients of Tracers (?)
          nftrw = openFile(trim(var%dirInp) // trim(TracInp), 'unformatted', 'sequential', -1, 'write', 'replace')
          if(nftrw < 0) return
        end if

        do nt = 3, Ntracers
          call Collect_Spec(qSpHuInp(:, :, nt), qWorkInOut, kmaxloc_in, var%kMaxInp, 0)
          if (isMasterProc()) then
            do k = 1, var%kMaxInp
              qWorkprOut = qWorkInOut(:, k)
              write (unit = nftrw)qWorkprOut
            end do
          end if
        end do
        if  (isMasterProc()) close(unit = nftrw)
      end if
    end if

    isExecOk = .true.
    
  end function ICWriteGDAS


  subroutine ICRecomposition
    !# ICRecomposition ... ?
    !# ---
    !# @info
    !# **Brief:** ICRecomposition ... ? </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 <br>
    !# @endin
    implicit none
    integer :: k, ib, jb, ns, nf
    integer :: jbFirst
    integer :: jbLast
    integer :: ibFirst
    integer :: ibLast
    integer :: kFirst
    integer :: kLast
    integer :: mnRIExtFirst
    integer :: mnRIExtLast

    call ThreadDecomp(1, jbMax, jbFirst, jbLast, "ICRecomp ")
    call ThreadDecomp(1, ibMax, ibFirst, ibLast, "ICRecomp ")
    call ThreadDecomp(1, kMaxloc, kFirst, kLast, "ICRecomp ")
    call ThreadDecomp(1, 2 * mymnExtMax, mnRIExtFirst, mnRIExtLast, "ICRecomp ")
    !
    call dztouv(qDivgInp, qVortInp, qUvelInp, qVvelInp, mnRIExtFirst, mnRIExtLast)

    ns = 2
    nf = 3 + Ntracers
    if (GetNewTop .or. var%smoothTopo)  ns = ns + 1
    if (var%grADS) nf = nf + 2
    if (isMasterProc()) write(p_nfprt, *) 'a nf, ns ', nf, ns
    !$OMP BARRIER
    !$OMP SINGLE
    call CreateSpecToGrid(nf, ns, nf, ns)

    call DepositSpecToGrid(qTopoInp, gTopoInp)
    call DepositSpecToGrid(qLnpsInp, gLnpsInp)
    if (GetNewTop .or. var%smoothTopo) then
      call DepositSpecToGrid(qTopoOut, gTopoOut)
    endif
    if (var%grADS) then
      call DepositSpecToGrid(qDivgInp, gDivgInp)
      call DepositSpecToGrid(qVortInp, gVortInp)
    endif
    call DepositSpecToGrid(qTvirInp, gTvirInp)
    call DepositSpecToGrid(qUvelInp, gUvelInp)
    call DepositSpecToGrid(qVvelInp, gVvelInp)
    do nt = 1, Ntracers
      call DepositSpecToGrid(qSpHuInp(:, :, nt), gSpHuInp(:, :, :, nt))
    enddo
    !$OMP END SINGLE
    call DoSpecToGrid(var%tamBlock)
    !$OMP BARRIER
    !$OMP SINGLE
    call DestroySpecToGrid()
    !$OMP END SINGLE
    call msgOutMaster(headerMsg // 'ICRdecomp ', ' after spectogrid ')
    !   m1 =0.
    !   m2 = 1.e7
    do j = 1, Jbmax
      do I = 1, Ibmaxperjb(j)
        gPsfcInp(i, j) = 10.0_p_r8 * exp(gLnPsInp(i, j))
      end do
    end do

    if (GetNewTop .or. var%smoothTopo) then
      if(GetNewTop)then
        gTopoOut = 0.0
        do jb = 1, jbMax
          do ib = 1, Ibmaxperjb(jb)
            gTopoOut(ib, jb) = gTopoOutGaus(ib, jb)
          end do
        end do
      end if
      do j = 1, Jbmax
        do I = 1, Ibmaxperjb(j)
          gTopoDel(i, j) = gTopoOut(i, j) - gTopoInp(i, j)
        end do
      end do

    else

      do j = 1, Jbmax
        do I = 1, Ibmaxperjb(j)
          if(gTopoInp(i, j)<0.5)gTopoInp(i, j) = 0.0
          gTopoOut(i, j) = gTopoInp(i, j)
        end do
      end do

    end if

    do j = 1, Jbmax
      do k = 1, var%kMaxInp
        do i = 1, Ibmaxperjb(j)
          gUvelInp(i, k, j) = gUvelInp(i, k, j) / coslat(i, j)
          gVvelInp(i, k, j) = gVvelInp(i, k, j) / coslat(i, j)
        end do
      end do
    end do
    !   m1 =0.
    !   m2 = 1.e7
    !   m3 =0.
    !   m4 = 1.e7
    !   DO j=1,Jbmax
    !      DO I=1,Ibmaxperjb(j)
    !       m1 = max(m1,gtvirInp(i,1,j))
    !       m2 = min(m2,gtvirInp(i,1,j))
    !       m3 = max(m3,gtvirInp(i,var%kMaxInp,j))
    !       m4 = min(m4,gtvirInp(i,var%kMaxInp,j))
    !     END DO
    !   END DO
    !    write(p_nfprt,*) myId , 'max tvir 1',m1, ' min ', m2
    !    write(p_nfprt,*) myId , 'max tvir kmax',m3, ' min ', m4

  end subroutine ICRecomposition


  subroutine ICDecomposition()
    !# ICDecomposition ... ?
    !# ---
    !# @info
    !# **Brief:** ICDecomposition ... ? </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 <br>
    !# @endin
    implicit none

    integer :: ib, jb, k

    call ReshapeVerticalGroups(var%kMaxOut, kMaxloc_out)
    kmaxloc = kMaxloc_out
    kmax = var%kMaxOut

    do jb = 1, jbmax
      do k = 1, var%kMaxOut
        do ib = 1, ibMaxperjb(jb)
          gUvelOut(ib, k, jb) = gUvelOut(ib, k, jb) * coslat(ib, jb) * rcl(ib, jb)
          gVvelOut(ib, k, jb) = gVvelOut(ib, k, jb) * coslat(ib, jb) * rcl(ib, jb)
        end do
      end do
    end do
    !$OMP BARRIER
    !$OMP SINGLE
    call CreateGridToSpec(3 + Ntracers, 2)
    call DepositGridToSpec(qTopoOut, gTopoOut)
    call DepositGridToSpec(qLnpsOut, gLnpsOut)
    call DepositGridToSpec(qUvelOut, gUvelOut)
    call DepositGridToSpec(qVvelOut, gVvelOut)
    call DepositGridToSpec(qTvirOut, gTvirOut)
    do nt = 1, Ntracers
      call DepositGridToSpec(qSpHuOut(:, :, nt), gSpHuOut(:, :, :, nt))
    enddo
    !$OMP END SINGLE
    call DoGridToSpec(var%tamBlock)
    !$OMP BARRIER
    !$OMP SINGLE
    call DestroyGridToSpec()
    !$OMP END SINGLE
    
    call Uvtodz(qUvelOut, qVvelOut, qDivgOut, qVortOut, 1, 2 * mymnmax)
    ! obtain div and vort tendencies
    
    if (GetNewTop) then
      qTopoOut = qTopoOutSpec
    end if

  end subroutine ICDecomposition


  subroutine ICDecompositionInput()
    !# ICDecompositionInput ... ?
    !# ---
    !# @info
    !# **Brief:** ICDecompositionInput ... ? </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 <br>
    !# @endin
    implicit none

    integer :: ib, jb, k

    call ReshapeVerticalGroups(var%kMaxInp, kMaxloc_In)
    kmaxloc = kMaxloc_In
    kmax = var%kMaxInp

    do jb = 1, jbmax
      do k = 1, var%kMaxInp
        do ib = 1, ibMaxperjb(jb)
          gUvelInp(ib, k, jb) = gUvelInp(ib, k, jb) * coslat(ib, jb) * rcl(ib, jb)
          gVvelInp(ib, k, jb) = gVvelInp(ib, k, jb) * coslat(ib, jb) * rcl(ib, jb)
        end do
      end do
    end do
    !$OMP BARRIER
    !$OMP SINGLE
    !    CALL CreateGridToSpec(3+Ntracers, 2)
    call CreateGridToSpec(6, 2)
    call DepositGridToSpec_PK(qTopoInp, gTopoInp)
    call DepositGridToSpec_PK(qLnpsInp, gLnpsInp)
    call DepositGridToSpec_PK(qUvelInp, gUvelInp)
    call DepositGridToSpec_PK(qVvelInp, gVvelInp)
    call DepositGridToSpec_PK(qTvirInp, gTvirInp)
    call DepositGridToSpec_PK(qSpHuInp(:, :, 1), gSpHuInp(:, :, :, 1))
    call DepositGridToSpec_PK(qSpHuInp(:, :, 2), gSpHuInp(:, :, :, 2))
    call DepositGridToSpec_PK(qSpHuInp(:, :, 3), gSpHuInp(:, :, :, 3))

    !$OMP END SINGLE
    call DoGridToSpec(var%tamBlock)
    !$OMP BARRIER
    !$OMP SINGLE
    call DestroyGridToSpec()
    !$OMP END SINGLE
    ! 
    !   obtain div and vort tendencies
    ! 
    call Uvtodz(qUvelInp, qVvelInp, qDivgInp, qVortInp, 1, 2 * mymnmax)

    if (GetNewTop) then
      qTopoInp = qTopoOutSpec
    end if

  end subroutine ICDecompositionInput


  subroutine SmoothCoef()
    !# Smoothes Spherical Harmonics Coefficients
    !# ---
    !# @info
    !# **Brief:** Smoothes Spherical Harmonics Coefficients. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 <br>
    !# @endin

    !          Using Hoskin's(?) Filter

    integer :: mn, n

    real (kind = p_r8) :: rm, rn, cx, red, rmn, ck

    rm = real(2 * var%mendcut, p_r8) * real(2 * var%mendcut + 1, p_r8)
    rn = real(var%mendcut - 1, p_r8) * real(var%mendcut, p_r8)
    red = (var%smthpercut)**(-(rm * rm) / (rn * rn) / var%iter)
    cx = -LOG(red) / (rm * rm)
    do mn = 1, mymnmax
      n = mynMap(mn)
      rmn = real(n - 1, p_r8) * real(n, p_r8)
      ck = (exp(cx * rmn * rmn))**var%iter
      qTopoOut(2 * mn - 1) = qTopoOut(2 * mn - 1) * ck
      qTopoOut(2 * mn) = qTopoOut(2 * mn) * ck
    enddo

  end subroutine SmoothCoef


  subroutine SmoothCoef2()
    !# Smoothes Spherical Harmonics Coefficients 2
    !# ---
    !# @info
    !# **Brief:** Smoothes Spherical Harmonics Coefficients 2. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 <br>
    !# @endin  
  
    !          Using Hoskin's(?) Filter
  
    integer :: mn, n
  
    real (kind = p_r8) :: rm, rn, cx, red, rmn, ck
 
    rm = real(2 * var%mendcut, p_r8) * real(2 * var%mendcut + 1, p_r8)
    rn = real(var%mendcut - 1, p_r8) * real(var%mendcut, p_r8)
    red = (var%smthpercut)**(-(rm * rm) / (rn * rn) / var%iter)
    cx = -log(red) / (rm * rm)
    do mn = 1, mymnmax
      n = mynMap(mn)
      rmn = real(n-1, p_r8) * real(n, p_r8)
      ck = (exp(cx * rmn * rmn))**var%iter
      qTopoInp(2 * mn - 1) = qTopoInp(2 * mn - 1) * ck
      qTopoInp(2 * mn) = qTopoInp(2 * mn) * ck
    enddo
  
  end subroutine SmoothCoef2


  function GetGrADSInp() result(isExecOk)
    !# Writes GradsInp file
    !# ---
    !# @info
    !# **Brief:** Writes GradsInp file. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 <br>
    !# @endin 
    implicit none
    logical :: isExecOk
    integer :: nfctl, nfgrd

    isExecOk = .false.

    if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A,L6,/)') ' GrADS     = ', var%grADS
    if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A,L6,/)') ' grADSOnly = ', var%grADSOnly

    if (isMasterProc()) then
      inquire (IOLENGTH = IOL) gWorkprOut
      nfgrd = openFile(trim(var%dirGrd) // trim(DataOut) // '.GrADS', 'unformatted', 'direct', IOL, 'write', 'replace')
      if(nfgrd < 0) return
    endif

    call Collect_Grid_Full(gTopoInp, gWorkOut(:, :, 1), 1, 0)
    if(isMasterProc()) then
      gWorkprOut = gWorkOut(:, :, 1)
      nRec = 1
      write (unit = nfgrd, REC = nRec) gWorkprOut
    endif
    call Collect_Grid_Full(gPsfcInp, gWorkOut(:, :, 1), 1, 0)
    if(isMasterProc()) then
      gWorkprOut = gWorkOut(:, :, 1)
      nRec = 2
      write (unit = nfgrd, REC = nRec) gWorkprOut
    endif
    call Collect_Grid_Full(gTvirInp, gWorkOut, var%kMaxInp, 0)
    if(isMasterProc()) then
      do k = 1, var%kMaxInp
        gWorkprOut = gWorkOut(:, :, k)
        nRec = 2 + int(k)
        write (unit = nfgrd, REC = nRec) gWorkprOut
      end do
    endif
    call Collect_Grid_Full(gDivgInp, gWorkOut, var%kMaxInp, 0)
    if(isMasterProc()) then
      do k = 1, var%kMaxInp
        gWorkprOut = gWorkOut(:, :, k)
        nRec = 2 + int(k + var%kMaxInp)
        write (unit = nfgrd, REC = nRec) gWorkprOut
      end do
    endif
    call Collect_Grid_Full(gVortInp, gWorkOut, var%kMaxInp, 0)
    if(isMasterProc()) then
      do k = 1, var%kMaxInp
        gWorkprOut = gWorkOut(:, :, k)
        nRec = 2 + int(k + 2 * var%kMaxInp)
        write (unit = nfgrd, REC = nRec) gWorkprOut
      end do
    endif
    call Collect_Grid_Full(gSpHuInp(:, :, :, 1), gWorkOut, var%kMaxInp, 0)
    if(isMasterProc()) then
      do k = 1, var%kMaxInp
        gWorkprOut = gWorkOut(:, :, k)
        nRec = 2 + int(k + 3 * var%kMaxInp)
        write (unit = nfgrd, REC = nRec) gWorkprOut
      end do
    endif
    call Collect_Grid_Full(gUvelInp, gWorkOut, var%kMaxInp, 0)
    if(isMasterProc()) then
      do k = 1, var%kMaxInp
        gWorkprOut = gWorkOut(:, :, k)
        nRec = 2 + int(k + 4 * var%kMaxInp)
        write (unit = nfgrd, REC = nRec) gWorkprOut
      end do
    endif
    call Collect_Grid_Full(gVvelInp, gWorkOut, var%kMaxInp, 0)
    if(isMasterProc()) then
      do k = 1, var%kMaxInp
        gWorkprOut = gWorkOut(:, :, k)
        nRec = 2 + int(k + 5 * var%kMaxInp)
        write (unit = nfgrd, REC = nRec) gWorkprOut
      end do
    endif
    if (var%getOzone) then
      nt = 2
      call Collect_Grid_Full(gSpHuInp(:, :, :, nt), gWorkOut, var%kMaxInp, 0)
      if(isMasterProc()) then
        do k = 1, var%kMaxInp
          gWorkprOut = gWorkOut(:, :, k)
          nRec = 2 + int(k + 6 * var%kMaxInp)
          write (unit = nfgrd, REC = nRec) gWorkprOut
        end do
      endif
    end if
    if (var%getTracers) then
      do nt = 3, NTracers
        call Collect_Grid_Full(gSpHuInp(:, :, :, nt), gWorkOut, var%kMaxInp, 0)
        if(isMasterProc()) then
          do k = 1, var%kMaxInp
            gWorkprOut = gWorkOut(:, :, k)
            nRec = 2 + int(k + (4 + nt) * var%kMaxInp)
            write (unit = nfgrd, REC = nRec) gWorkprOut
          end do
        endif
      end do
    end if

    if (isMasterProc()) then
      close (unit = nfgrd)
    endif

    if (isMasterProc()) then
      nfctl = openFile(trim(var%dirGrd) // trim(DataOut) // '.GrADS.ctl', 'formatted', 'sequential', -1, 'write', 'replace')
      if(nfctl < 0) return
      write (unit = nfctl, FMT = '(A)') 'DSET ^' // trim(DataOut) // '.GrADS'
      write (unit = nfctl, FMT = '(A)') 'options yrev big_endian'
      write (unit = nfctl, FMT = '(A)') 'undef 1e20'
      write (unit = nfctl, FMT = '(A,I5,A,2F12.6)') 'xdef ', ImaxOut, ' linear ', &
        0.0_p_r8, 360.0_p_r8 / real(ImaxOut, p_r8)
      write (unit = nfctl, FMT = '(A,I5,A)') 'ydef ', JmaxOut, ' levels '
      write (unit = nfctl, FMT = '(6F10.3)') lati(JmaxOut:1:-1)
      write (unit = nfctl, FMT = '(A,I5,A)') 'zdef ', var%kMaxInp, ' levels '

      if (var%gdasonly) then 
        ! Sigma-p case, first calculate pressure for each interface and then
        ! calculate pressure on each layer according to vertical structure
        write (unit = nfctl, FMT='(6F10.3)') SigInterInp(1:var%kMaxInp)+1000.0_p_r8*SigLayerInp(1:var%kMaxInp)
      else
        !      write (unit = nfctl, FMT = '(6F10.3)') 1000.0_p_r8 * SigLayerInp
        if ( IDVCInp == 0_p_i4 .or. IDVCInp == 1_p_i4 ) then
           ! Sigma case, calculate pressure for each layer directly
            write (unit = nfctl, FMT='(6F10.3)') 1000.0_p_r8*SigLayerInp(1:var%kMaxInp)
        else if ( IDVCInp == 2_p_i4 ) then
            ! Sigma-p case, first calculate pressure for each interface and then
            ! calculate pressure on each layer according to vertical structure
            write (unit = nfctl, FMT='(6F10.3)') SigInterInp(1:var%kMaxInp)+1000.0_p_r8*SigLayerInp(1:var%kMaxInp)
        end if
      endif
      write (unit = nfctl, FMT = '(3A)') 'tdef 1 linear ', Tdef, ' 1dy'
      if (NTracers == 1) then
        write (unit = nfctl, FMT = '(A)') 'vars 8'
      else
        write (unit = nfctl, FMT = '(A,I5)') 'vars ', 7 + NTracers
      end if
      write (unit = nfctl, FMT = '(A)') 'topo   0 99 Topography        ' // TruncInp // ' (m)'
      write (unit = nfctl, FMT = '(A)') 'pslc   0 99 Surface Pressure  ' // TruncInp // ' (hPa)'
      write (unit = nfctl, FMT = '(A,I3,A)') 'tvir ', var%kMaxInp, ' 99 Virt Temperature  ' // &
        TruncInp // ' (K)'
      write (unit = nfctl, FMT = '(A,I3,A)') 'divg ', var%kMaxInp, ' 99 Divergence        ' // &
        TruncInp // ' (1/s)'
      write (unit = nfctl, FMT = '(A,I3,A)') 'vort ', var%kMaxInp, ' 99 Vorticity         ' // &
        TruncInp // ' (1/s)'
      write (unit = nfctl, FMT = '(A,I3,A)') 'umes ', var%kMaxInp, ' 99 Specific Humidity ' // &
        TruncInp // ' (kg/kg)'
      write (unit = nfctl, FMT = '(A,I3,A)') 'uvel ', var%kMaxInp, ' 99 Zonal Wind        ' // &
        TruncInp // ' (m/s)'
      write (unit = nfctl, FMT = '(A,I3,A)') 'vvel ', var%kMaxInp, ' 99 Meridional Wind   ' // &
        TruncInp // ' (m/s)'
      if (var%getOzone) then
        write (unit = nfctl, FMT = '(A,I3,A)') 'ozon ', var%kMaxInp, ' 99 Ozone             ' // &
          TruncInp // ' (?)'
      end if
      if (var%getTracers) then
        do nt = 3, NTracers
          write (unit = nfctl, FMT = '(A,I1,A,I3,A)') 'trc', nt - 2, ' ', var%kMaxInp, &
            ' 99 Tracer            ' // TruncInp // ' (?)'
        end do
      end if
      write (unit = nfctl, FMT = '(A)') 'endvars'
      close (unit = nfctl)
    endif

    isExecOk = .true.

  end function GetGrADSInp


  function GetGrADSOut() result(isExecOk)
    !# Writes GradsOut file
    !# ---
    !# @info
    !# **Brief:** Writes GradsOut file. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 <br>
    !# @endin 
    implicit none
    logical :: isExecOk
    integer :: nfctl, nfgrd

    isExecOk = .false.

    if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A,L6,/)') ' GrADS     = ', var%grADS
    if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A,L6,/)') ' GrADSOnly = ', var%grADSOnly

    if (isMasterProc()) then
      inquire (IOLENGTH = IOL) gWorkprOut
      nfgrd = openFile(trim(var%dirGrd) // trim(DataOut) // '.GrADSOut', 'unformatted', 'direct', IOL, 'write', 'replace')
      if(nfgrd < 0) return
    endif
    call Collect_Grid_Full(gTopoOut, gWorkOut(:, :, 1), 1, 0)
    if(isMasterProc()) then
      gWorkprOut = gWorkOut(:, :, 1)
      nRec = 1
      write (unit = nfgrd, REC = nRec)gWorkprOut
    endif
    call Collect_Grid_Full(gPsfcOut, gWorkOut(:, :, 1), 1, 0)
    if(isMasterProc()) then
      gWorkprOut = gWorkOut(:, :, 1)
      nRec = 2
      write (unit = nfgrd, REC = nRec) gWorkprOut
    endif
    call Collect_Grid_Full(gTvirOut, gWorkOut, var%kMaxOut, 0)
    if(isMasterProc()) then
      do k = 1, var%kMaxOut
        gWorkprOut = gWorkOut(:, :, k)
        nRec = 2 + int(k)
        write (unit = nfgrd, REC = nRec) gWorkprOut
      end do
    endif
    call Collect_Grid_Full(gpresOut, gWorkOut, var%kMaxOut, 0)
    if(isMasterProc()) then
      do k = 1, var%kMaxOut
        gWorkprOut = gWorkOut(:, :, k)
        nRec = 2 + int(k + var%kMaxOut)
        write (unit = nfgrd, REC = nRec) gWorkprOut
      end do
    endif
    call Collect_Grid_Full(gUvelOut, gWorkOut, var%kMaxOut, 0)
    if(isMasterProc()) then
      do k = 1, var%kMaxOut
        gWorkprOut = gWorkOut(:, :, k)
        nRec = 2 + int(k + 2 * var%kMaxOut)
        write (unit = nfgrd, REC = nRec) gWorkprOut
      end do
    endif
    call Collect_Grid_Full(gVvelOut, gWorkOut, var%kMaxOut, 0)
    if(isMasterProc()) then
      do k = 1, var%kMaxOut
        gWorkprOut = gWorkOut(:, :, k)
        nRec = 2 + int(k + 3 * var%kMaxOut)
        write (unit = nfgrd, REC = nRec) gWorkprOut
      end do
    endif
    call Collect_Grid_Full(gSpHuOut(:, :, :, 1), gWorkOut, var%kMaxOut, 0)
    if(isMasterProc()) then
      do k = 1, var%kMaxOut
        gWorkprOut = gWorkOut(:, :, k)
        nRec = 2 + int(k + 4 * var%kMaxOut)
        write (unit = nfgrd, REC = nRec) gWorkprOut
      end do
    endif
    if (isMasterProc()) then
      close (unit = nfgrd)
    endif

    if (isMasterProc()) then
      nfctl = openFile(trim(var%dirGrd) // trim(DataOut) // '.GrADSOut.ctl', 'formatted', 'sequential', -1, 'write', 'replace')
      if(nfctl < 0) return
      write (unit = nfctl, FMT = '(A)') 'DSET ^' // trim(DataOut) // '.GrADSOut'
      write (unit = nfctl, FMT = '(A)') 'options yrev big_endian'
      write (unit = nfctl, FMT = '(A)') 'undef 1e20'
      write (unit = nfctl, FMT = '(A,I5,A,2F12.6)') 'xdef ', ImaxOut, ' linear ', &
        0.0_p_r8, 360.0_p_r8 / real(ImaxOut, p_r8)
      write (unit = nfctl, FMT = '(A,I5,A)') 'ydef ', JmaxOut, ' levels '
      write (unit = nfctl, FMT = '(6F10.3)') lati(JmaxOut:1:-1)
      write (unit = nfctl, FMT = '(A,I5,A)') 'zdef ', var%kMaxOut, ' levels '
      
      if (var%gdasonly) then
         ! Sigma-p case, first calculate pressure for each interface and then
         ! calculate pressure on each layer according to vertical structure
         !write (unit = nfctl, FMT='(6F10.3)') 1000.0_p_r8 * SigLayerOut(1:var%kMaxOut)
         write (unit = nfctl, FMT='(6F10.3)') SigInterOut(1:var%kMaxOut)+1000.0_p_r8*SigLayerOut(1:var%kMaxOut)
      else
         ! bug: level > 1000 and 1 level more
         ! WRITE (UNIT=nfctl, FMT='(6F10.3)') SigInterOut+1000.0_r8*SigInterOut
      
         ! using solution from file Grads.ctl ... but using "Out" variables !
         if ( IDVCInp == 0_p_i4 .or. IDVCInp == 1_p_i4 ) then
             ! Sigma case, calculate pressure for each layer directly
             write (unit = nfctl, FMT='(6F10.3)') 1000.0_p_r8 * SigLayerOut(1:var%kMaxOut)
         else if ( IDVCInp == 2_p_i4 ) then
             ! Sigma-p case, first calculate pressure for each interface and then
             ! calculate pressure on each layer according to vertical structure
             !write (unit = nfctl, FMT='(6F10.3)') 1000.0_p_r8 * SigLayerOut(1:var%kMaxOut)

             write (unit = nfctl, FMT='(6F10.3)') SigInterOut(1:var%kMaxOut)+1000.0_p_r8*SigLayerOut(1:var%kMaxOut)
         end if
      endif
      write (unit = nfctl, FMT = '(3A)') 'tdef 1 linear ', Tdef, ' 1dy'
      if (NTracers == 1) then
        write (unit = nfctl, FMT = '(A)') 'vars 7'
      else
        write (unit = nfctl, FMT = '(A,I5)') 'vars ', 4 + 3
      end if
      write (unit = nfctl, FMT = '(A)') 'topo   0 99 Topography        ' // TruncInp // ' (m)'
      write (unit = nfctl, FMT = '(A)') 'pslc   0 99 Surface Pressure  ' // TruncInp // ' (hPa)'
      write (unit = nfctl, FMT = '(A,I3,A)') 'tvir ', var%kMaxOut, ' 99 Virt Temperature  ' // &
        TruncInp // ' (K)'
      write (unit = nfctl, FMT = '(A,I3,A)') 'pres ', var%kMaxOut, ' 99 Pressure ' // &
        TruncInp // ' (kg/kg)'
      write (unit = nfctl, FMT = '(A,I3,A)') 'uvel ', var%kMaxOut, ' 99 Zonal Wind        ' // &
        TruncInp // ' (m/s)'
      write (unit = nfctl, FMT = '(A,I3,A)') 'vvel ', var%kMaxOut, ' 99 Meridional Wind   ' // &
        TruncInp // ' (m/s)'
      write (unit = nfctl, FMT = '(A,I3,A)') 'umes ', var%kMaxOut, ' 99 Specific Humidity Wind   ' // &
        TruncInp // ' (m/s)'
      write (unit = nfctl, FMT = '(A)') 'endvars'
      close (unit = nfctl)
    endif

    if (var%getOzone) then
      if (isMasterProc()) then
        inquire (iolength = IOL) gWorkprOut
        nfgrd = openFile(trim(var%dirGrd)//trim(DataOut)//'.OZONEGrADSOut', 'unformatted', 'direct', IOL , 'write', 'replace')
        if (nfgrd < 0) return
      endif
      call Collect_Grid_Full(gSpHuOut(:,:,:,2), gWorkOut, var%kMaxOut, 0)
      if (isMasterProc()) then
        do k = 1, var%kMaxOut
            gWorkprOut = gWorkOut(:,:,k)
            nRec = 0 + int(k)
            write (unit = nfgrd, rec = nRec) gWorkprOut
        end do
        close (unit = nfgrd)
      endif

      if (isMasterProc()) then
        nfctl = openFile(trim(var%dirGrd)//trim(DataOut)//'.OZONEGrADSOut.ctl', 'formatted', 'sequential', -1, 'write', 'replace')
        if (nfctl < 0) return
        write (unit = nfctl, fmt = '(A)') 'dset ^'//trim(DataOut)//'.OZONEGrADSOut'
        write (unit = nfctl, fmt = '(A)') 'options yrev big_endian'
        write (unit = nfctl, fmt = '(A)') 'undef 1e20'
        write (unit = nfctl, fmt = '(A,I5,A,2F12.6)') 'xdef ',ImaxOut,' linear ', &
              0.0_p_r8, 360.0_p_r8/real(ImaxOut,p_r8)
        write (unit = nfctl, fmt = '(A,I5,A)') 'ydef ',JmaxOut,' levels '
        write (unit = nfctl, fmt = '(6F10.3)') lati(JmaxOut:1:-1)
        write (unit = nfctl, fmt = '(A,I5,A)') 'zdef ',var%kMaxOut,' levels '

        if (var%gdasonly) then        
              ! Sigma-p case, first calculate pressure for each interface and then
              ! calculate pressure on each layer according to vertical structure
              !write (unit = nfctl, FMT='(6F10.3)') SigInterOut(1:var%kMaxOut)+1000.0_p_r8*SigLayerOut(1:var%kMaxOut)
              ! write (unit = nfctl, FMT='(6F10.3)') 1000.0_p_r8 * SigLayerOut(1:var%kMaxOut)
                write (unit = nfctl, FMT='(6F10.3)') SigInterOut(1:var%kMaxOut)+1000.0_p_r8*SigLayerOut(1:var%kMaxOut)
        else
           ! bug: level > 1000 and 1 level more
           ! WRITE (UNIT=nfctl, FMT='(6F10.3)') SigInterOut+1000.0_r8*SigInterOut
        
           ! using solution from file Grads.ctl ... but using "Out" variables !
           if ( IDVCInp == 0_p_i4 .or. IDVCInp == 1_p_i4 ) then
              ! Sigma case, calculate pressure for each layer directly
              write (unit = nfctl, FMT='(6F10.3)') 1000.0_p_r8 * SigLayerOut(1:var%kMaxOut)
           else if ( IDVCInp == 2_p_i4 ) then
              ! Sigma-p case, first calculate pressure for each interface and then
              ! calculate pressure on each layer according to vertical structure
                write (unit = nfctl, FMT='(6F10.3)') SigInterOut(1:var%kMaxOut)+1000.0_p_r8*SigLayerOut(1:var%kMaxOut)
              ! write (unit = nfctl, FMT='(6F10.3)') 1000.0_p_r8 * SigLayerOut(1:var%kMaxOut)
           end if
        endif
        write (unit = nfctl, fmt = '(3A)') 'tdef 1 linear ', Tdef,' 1dy'
        if (NTracers == 1) then
          write (unit = nfctl, fmt = '(A)') 'vars 1'
        else
          write (unit = nfctl, fmt = '(A,I5)') 'vars ',0+1
        end if
        write (unit = nfctl, fmt = '(A,I3,A)') 'ozon ',var%kMaxOut,' 99 OZONE   '// &
            TruncInp//' (kg/kg)'
        write (unit = nfctl, fmt = '(A)') 'endvars'
        close (unit = nfctl)
      endif
    endif
    
    isExecOk = .true.
    
  end function GetGrADSOut


  function GetGrADSInp_GDAS() result(isExecOk)
    !# Writes GrADS_PK file
    !# ---
    !# @info
    !# **Brief:** Writes GrADS_PK file. </br>
    !# **Authors**: </br>
    !# &bull; Denis Eiras </br>
    !# **Date**: apr/2019 <br>
    !# @endin
    implicit none
    logical :: isExecOk

    integer :: nfctl, nfgrd

    isExecOk = .false.
        
    if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A,L6,/)') ' GrADS     = ', var%grADS
    if(isMasterProc()) write (unit = p_nfprt, FMT = '(/,A,L6,/)') ' GrADSOnly = ', var%grADSOnly

    if (isMasterProc()) then
      inquire (IOLENGTH = IOL) gWorkprOut
      nfgrd = openFile(trim(var%dirGrd) // trim(DataOut) // '.GrADS_PK', 'unformatted', 'direct', IOL, 'write', 'replace')
      if(nfgrd < 0) return
    endif

    call Collect_Grid_Full(gTopoInp, gWorkOut(:, :, 1), 1, 0)
    if(isMasterProc()) then
      gWorkprOut = gWorkOut(:, :, 1)
      nRec = 1
      write (unit = nfgrd, REC = nRec) gWorkprOut
    endif
    call Collect_Grid_Full(gPsfcInp, gWorkOut(:, :, 1), 1, 0)
    if(isMasterProc()) then
      gWorkprOut = gWorkOut(:, :, 1)
      nRec = 2
      write (unit = nfgrd, REC = nRec) gWorkprOut
    endif
    call Collect_Grid_Full(gTvirInp, gWorkOut, var%kMaxInp, 0)
    if(isMasterProc()) then
      do k = 1, var%kMaxInp
        gWorkprOut = gWorkOut(:, :, k)
        nRec = 2 + int(k)
        write (unit = nfgrd, REC = nRec) gWorkprOut
      end do
    endif
    call Collect_Grid_Full(gUvelInp, gWorkOut, var%kMaxInp, 0)
    if(isMasterProc()) then
      do k = 1, var%kMaxInp
        gWorkprOut = gWorkOut(:, :, k)
        nRec = 2 + int(k + var%kMaxInp)
        write (unit = nfgrd, REC = nRec) gWorkprOut
      end do
    endif
    call Collect_Grid_Full(gVvelInp, gWorkOut, var%kMaxInp, 0)
    if(isMasterProc()) then
      do k = 1, var%kMaxInp
        gWorkprOut = gWorkOut(:, :, k)
        nRec = 2 + int(k + 2 * var%kMaxInp)
        write (unit = nfgrd, REC = nRec) gWorkprOut
      end do
    endif
    call Collect_Grid_Full(gSpHuInp(:, :, :, 1), gWorkOut, var%kMaxInp, 0)
    if(isMasterProc()) then
      do k = 1, var%kMaxInp
        gWorkprOut = gWorkOut(:, :, k)
        nRec = 2 + int(k + 3 * var%kMaxInp)
        write (unit = nfgrd, REC = nRec) gWorkprOut
      end do
    endif
    call Collect_Grid_Full(gUvelInp, gWorkOut, var%kMaxInp, 0)
    if(isMasterProc()) then
      do k = 1, var%kMaxInp
        gWorkprOut = gWorkOut(:, :, k)
        nRec = 2 + int(k + 4 * var%kMaxInp)
        write (unit = nfgrd, REC = nRec) gWorkprOut
      end do
    endif
    call Collect_Grid_Full(gVvelInp, gWorkOut, var%kMaxInp, 0)
    if(isMasterProc()) then
      do k = 1, var%kMaxInp
        gWorkprOut = gWorkOut(:, :, k)
        nRec = 2 + int(k + 5 * var%kMaxInp)
        write (unit = nfgrd, REC = nRec) gWorkprOut
      end do
    endif
    if (var%getOzone) then
      nt = 2
      call Collect_Grid_Full(gSpHuInp(:, :, :, nt), gWorkOut, var%kMaxInp, 0)
      if(isMasterProc()) then
        do k = 1, var%kMaxInp
          gWorkprOut = gWorkOut(:, :, k)
          nRec = 2 + int(k + 6 * var%kMaxInp)
          write (unit = nfgrd, REC = nRec) gWorkprOut
        end do
      endif
    end if
    if (var%getTracers) then
      do nt = 3, NTracers
        call Collect_Grid_Full(gSpHuInp(:, :, :, nt), gWorkOut, var%kMaxInp, 0)
        if(isMasterProc()) then
          do k = 1, var%kMaxInp
            gWorkprOut = gWorkOut(:, :, k)
            nRec = 2 + int(k + (4 + nt) * var%kMaxInp)
            write (unit = nfgrd, REC = nRec) gWorkprOut
          end do
        endif
      end do
    end if

    if (isMasterProc()) then
      close (unit = nfgrd)
    endif

    if (isMasterProc()) then
      nfctl = openFile(trim(var%dirGrd) // trim(DataOut) // '.GrADS_PK.ctl', 'formatted', 'sequential', -1, 'write', 'replace')
      if(nfctl < 0) return
      write (unit = nfctl, FMT = '(A)') 'DSET ^' // trim(DataOut) // '.GrADS_PK'
      write (unit = nfctl, FMT = '(A)') 'options yrev big_endian'
      write (unit = nfctl, FMT = '(A)') 'undef 1e20'
      write (unit = nfctl, FMT = '(A,I5,A,2F12.6)') 'xdef ', ImaxOut, ' linear ', &
        0.0_p_r8, 360.0_p_r8 / real(ImaxOut, p_r8)
      write (unit = nfctl, FMT = '(A,I5,A)') 'ydef ', JmaxOut, ' levels '
      write (unit = nfctl, FMT = '(6F10.3)') lati(JmaxOut:1:-1)
      write (unit = nfctl, FMT = '(A,I5,A)') 'zdef ', var%kMaxInp, ' levels '
      if (var%gdasonly) then        
         ! Sigma-p case, first calculate pressure for each interface and then
         ! calculate pressure on each layer according to vertical structure
               write (unit = nfctl, FMT='(6F10.3)') SigIInp(1:var%kMaxInp)+1000.0_p_r8*SigLInp(1:var%kMaxInp)
      else
           ! bug: level > 1000 and 1 level more
           ! WRITE (UNIT=nfctl, FMT='(6F10.3)') SigInterOut+1000.0_r8*SigInterOut
        
           ! using solution from file Grads.ctl ... but using "Out" variables !
           if ( IDVCInp == 0_p_i4 .or. IDVCInp == 1_p_i4 ) then
              ! Sigma case, calculate pressure for each layer directly
              write (unit = nfctl, FMT='(6F10.3)') 1000.0_p_r8 * SigLInp(1:var%kMaxInp)
           else if ( IDVCInp == 2_p_i4 ) then
              ! Sigma-p case, first calculate pressure for each interface and then
              ! calculate pressure on each layer according to vertical structure
              !write (unit = nfctl, FMT='(6F10.3)') SigInterOut(1:var%kMaxOut)+1000.0_p_r8*SigLayerOut(1:var%kMaxOut)
               write (unit = nfctl, FMT='(6F10.3)') SigIInp(1:var%kMaxInp)+1000.0_p_r8*SigLInp(1:var%kMaxInp)
           end if
      endif
      write (unit = nfctl, FMT = '(3A)') 'tdef 1 linear ', Tdef, ' 1dy'
      if (NTracers == 1) then
        write (unit = nfctl, FMT = '(A)') 'vars 8'
      else
        write (unit = nfctl, FMT = '(A,I5)') 'vars ', 7 + NTracers
      end if
      write (unit = nfctl, FMT = '(A)') 'topo   0 99 Topography        ' // TruncInp // ' (m)'
      write (unit = nfctl, FMT = '(A)') 'pslc   0 99 Surface Pressure  ' // TruncInp // ' (hPa)'
      write (unit = nfctl, FMT = '(A,I3,A)') 'tvir ', var%kMaxInp, ' 99 Virt Temperature  ' // &
        TruncInp // ' (K)'
      write (unit = nfctl, FMT = '(A,I3,A)') 'divg ', var%kMaxInp, ' 99 Divergence        ' // &
        TruncInp // ' (1/s)'
      write (unit = nfctl, FMT = '(A,I3,A)') 'vort ', var%kMaxInp, ' 99 Vorticity         ' // &
        TruncInp // ' (1/s)'
      write (unit = nfctl, FMT = '(A,I3,A)') 'umes ', var%kMaxInp, ' 99 Specific Humidity ' // &
        TruncInp // ' (kg/kg)'
      write (unit = nfctl, FMT = '(A,I3,A)') 'uvel ', var%kMaxInp, ' 99 Zonal Wind        ' // &
        TruncInp // ' (m/s)'
      write (unit = nfctl, FMT = '(A,I3,A)') 'vvel ', var%kMaxInp, ' 99 Meridional Wind   ' // &
        TruncInp // ' (m/s)'
      if (var%getOzone) then
        write (unit = nfctl, FMT = '(A,I3,A)') 'ozon ', var%kMaxInp, ' 99 Ozone             ' // &
          TruncInp // ' (?)'
      end if
      if (var%getTracers) then
        do nt = 3, NTracers
          write (unit = nfctl, FMT = '(A,I1,A,I3,A)') 'trc', nt - 2, ' ', var%kMaxInp, &
            ' 99 Tracer            ' // TruncInp // ' (?)'
        end do
      end if
      write (unit = nfctl, FMT = '(A)') 'endvars'
      close (unit = nfctl)
    endif

    isExecOk = .true.
    
  end function GetGrADSInp_GDAS


end module Mod_Chopping
