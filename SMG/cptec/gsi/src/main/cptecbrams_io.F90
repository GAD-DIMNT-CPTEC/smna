!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !MODULE: cptecBRAMS_io -- Implements a CPTEC BRAMS interface to Read/Write
!                           atmospheric and surface at GSI.
!
! !DESCRIPTON:This version reads a vfm file created by BRAMS output.
!             A later version will read directly from the BRAMS restart file.
!             The guess is read in by complete horizontal fields, one field
!             per processor, in parallel.  Each horizontal input field is 
!             converted from the staggered c-grid to an unstaggered a-grid.
!             On the c-grid, u is shifted 1/2 point in the negative x direction
!             and v 1/2 point in the negative y direction, but otherwise the
!             three grids are regular.  When the fields are read in, nothing
!             is done to mass variables, but wind variables are interpolated to
!             mass points.
!
!             
!                 
!\\
!\\
! !INTERFACE:
!
module cptecBRAMS_io

! GSI kinds
   use kinds, only: r_kind

   implicit none
   private
!
! !PUBLIC MEMBER FUNCTIONS:
!

   public :: readBrams    ! read BRAMS atmopheric guess fields
   public :: readBramsSfc ! read BRAMS surface guess fields
   public :: writeBrams   ! write BRAMS atmospheric analysis fields
!   public :: write_bramsfc! write BRAMS surface analysis fields

!
! !PRIVATE MEMBER VARIABLES:
!
   !
   ! Transfers Variables
   !

   real(r_kind), allocatable :: g_ps      ( :,: ) ! Surface Pressure (kPa)
   real(r_kind), allocatable :: g_z       ( :,: ) ! Orography (m)
   real(r_kind), allocatable :: g_tv      (:,:,:) ! Virtural Temperature (K)
   real(r_kind), allocatable :: g_u       (:,:,:) ! Zonal Wind (m/s)
   real(r_kind), allocatable :: g_v       (:,:,:) ! Meridional Wind (m/s)
   real(r_kind), allocatable :: g_q       (:,:,:) ! Specific Humidy (kg/kg)
   real(r_kind), allocatable :: g_pi      (:,:,:) ! pressure (exner function) (kPa)
   real(r_kind), allocatable :: g_tk      (:,:,:) ! Absolute Temperature (K)
   real(r_kind), allocatable :: g_gh      (:,:,:) ! Geopotential Height (m)
   real(r_kind), allocatable :: g_rh      (:,:,:) ! Relative Humidity ()

   real(r_kind), allocatable :: g_prsi      (:,:,:) ! 
   real(r_kind), allocatable :: g_prsl      (:,:,:) !

   real(r_kind), allocatable :: g_ql      (:,:,:) ! Mixing ratio for cloud (kg/kg)
   real(r_kind), allocatable :: g_qi      (:,:,:) ! Mixing ratio for ice (kg/kg)
   real(r_kind), allocatable :: g_qr      (:,:,:) ! Mixing ratio for rain (kg/kg)
   real(r_kind), allocatable :: g_qs      (:,:,:) ! Mixing ratio for snow (kg/kg)
   real(r_kind), allocatable :: g_qg      (:,:,:) ! Mixing ratio for graupel (kg/kg)
   real(r_kind), allocatable :: g_qnr     (:,:,:) ! Rain number concentration ()

   real(r_kind), allocatable :: gs_th2    ( :,: ) ! 
   real(r_kind), allocatable :: gs_q2     ( :,: ) !
   real(r_kind), allocatable :: gs_tsk    ( :,: ) !
   real(r_kind), allocatable :: gs_soilt1 ( :,: ) !
   real(r_kind), allocatable :: gs_tslb   (:,:,:) !
   real(r_kind), allocatable :: gs_smois  (:,:,:) !

!
! write out units
!

   integer, parameter :: stdinp = 5 ! standard input
   integer, parameter :: stdout = 6 ! standard output
   integer, parameter :: stderr = 0 ! standard error
!
! Some parameters
!

      integer, parameter :: iUdef = -9999
      real,    parameter :: Undef = -9.99999996E+11 !-9.99E9

!
! !REVISION HISTORY:
!
!   26 Mar 2019 - J. G. de Mattos -  Initial code.
!
! !REMARKS:
!          need to pay special attention to various surface fields to make sure
!          they are correct (check units).  fact10 needs to be computed from 
!          10m wind, which is included in wrf mass interface file.
!
!          Cloud water currently set to zero.  Need to fix before doing precip
!          assimilation--wait until global adopts Brad Ferrier's multi-species 
!          scheme??
!
!          Ozone currently set to zero.  No ozone variable in wrf mass core--climatology
!          instead.  Do we use climatology?
!
!          No background bias yet. (biascor ignored)
!
!EOP
!-----------------------------------------------------------------------------!

   character(len=64),parameter :: myname='cptecBRAMS_io'

   contains

!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: GenReadBRAMS
!              
! !DESCRIPTON: rotine to get atmospheric guess fields from BRAMS and return to
!              GSI
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine readBrams( mype )
!
! !USES:
!
      ! GSI kinds
      use kinds, only: i_kind, r_kind

      !GSI bundle
      use gsi_metguess_mod, only: gsi_metguess_bundle, gsi_metguess_get
      use gsi_bundlemod, only: gsi_bundlegetpointer, gsi_bundleprint

      use guess_grids, only: nfldsig, ntguessig, ifilesig
      
      use rapidrefresh_cldsurf_mod, only: l_cloud_analysis

      implicit none
!
! !INPUT PARAMETERS:
!
      integer, intent(in   ) :: mype ! mpi task id

! !REVISION HISTORY:
!
!   26 Mar 2019 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      character(len=64),parameter :: myname_=trim(myname)//' :: read_BRAMS( )'

      !
      ! GSI bungle guess fields
      !

      real(r_kind), pointer, dimension(  :,:) :: ges_ps_it     => NULL()

      real(r_kind), pointer, dimension(  :,:) :: ges_th2_it    => NULL()
      real(r_kind), pointer, dimension(  :,:) :: ges_q2_it     => NULL()
      real(r_kind), pointer, dimension(  :,:) :: ges_tsk_it    => NULL()
      real(r_kind), pointer, dimension(  :,:) :: ges_soilt1_it => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_tslb_it   => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_smois_it  => NULL()
      
      real(r_kind), pointer, dimension(:,:  ) :: ges_z_it      => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_u_it      => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_v_it      => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_tv_it     => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_q_it      => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_pi_it     => NULL()

      real(r_kind), pointer, dimension(:,:,:) :: ges_ql_it     => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_qi_it     => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_qr_it     => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_qs_it     => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_qg_it     => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: ges_qnr_it    => NULL()


      !
      ! Auxiliary vars
      !

      integer         :: istatus
      integer         :: it, i, ier
      integer(i_kind) :: n_actual_clouds
      character(len=20) :: FileAnlFG

#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      ! Inquire about cloud guess fields
      call gsi_metguess_get('clouds::3d',n_actual_clouds,istatus)
      if (n_actual_clouds>0) then
         ! Get pointer for each of the hydrometeors from guess at time index "it"
         ier=0
         it=ntguessig
         call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'ql', ges_ql_it, istatus );ier=ier+istatus
         call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qi', ges_qi_it, istatus );ier=ier+istatus
         call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qr', ges_qr_it, istatus );ier=ier+istatus
         call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qs', ges_qs_it, istatus );ier=ier+istatus
         call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qg', ges_qg_it, istatus );ier=ier+istatus
!         call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qnr',ges_qnr_it,istatus );ier=ier+istatus
         if (ier/=0) n_actual_clouds=0
      end if


      !
      ! Loop over all time gues fields
      !

      do it = 1, nfldsig

         !
         ! Allocating Atmospheric transfer fields
         !

         call AllocateAtmFields( )
 
         !
         ! Read BRAMS fields
         !

         call ReadAtmBRAMS ( mype, it )

         !
         ! Put fields at GSI Bundle
         !

         ! Orography

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'z', ges_z_it, istatus) 
         if(istatus.eq.0) ges_z_it = g_z
         
         ! Surface Pressure

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'ps', ges_ps_it, istatus) 
         if(istatus.eq.0) ges_ps_it = g_ps 

         ! pressure (exner function)

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'pi', ges_pi_it, istatus) 
         if(istatus.eq.0) ges_pi_it = g_pi

         ! Virtual Temperature

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'tv', ges_tv_it, istatus) 
         if(istatus.eq.0) ges_tv_it = g_tv

         ! Specific Humidity

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'q', ges_q_it, istatus) 
         if(istatus.eq.0) ges_q_it = g_q

         ! Zonal wind component

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'u', ges_u_it, istatus) 
         if(istatus.eq.0) ges_u_it = g_u

         ! Meridional wind component

         call gsi_bundlegetpointer (gsi_metguess_bundle(it), 'v', ges_v_it, istatus) 
         if(istatus.eq.0) ges_v_it = g_v

!         write(FileAnlFG,'("FileAnlFG.",I2.2)') ifilesig(it)
!         call GenWriteBRAMS_FG (FileAnlFG, mype, 0, it )

         call DeallocateAtmFields ( )

!         ! surface moisture observation
!         
!         if(i_use_2mq4b>0) then
!            call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it),'q2m',ges_q2_it,istatus )
!            if (ier/=0) call die(trim(myname),'cannot get pointers for q2m, ier =',istatus)
!         endif

         ! hydrometeors

         if(l_cloud_analysis .or. n_actual_clouds>0) then

            !
            ! Allocating Atmospheric transfer fields
            !
   
            call AllocateCloudFields( )
    
            !
            ! Read BRAMS fields
            !
   
            call ReadCloudBRAMS ( mype, it )
   
            !
            ! Put fields at GSI Bundle
            !

            ! Mixing ratio for cloud
            call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'ql', ges_ql_it, istatus )
            if(istatus.eq.0) ges_ql_it = g_ql

            ! Mixing ratio for ice
            call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qi', ges_qi_it, istatus )
            if(istatus.eq.0) ges_qi_it = g_qi

            ! Mixing ratio for rain
            call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qr', ges_qr_it, istatus )
            if(istatus.eq.0) ges_qr_it = g_qr

            ! Mixing ratio for snow
            call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qs', ges_qs_it, istatus )
            if(istatus.eq.0) ges_qs_it = g_qs

            ! Mixing ratio for graupel
            call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qg', ges_qg_it, istatus )
            if(istatus.eq.0) ges_qg_it = g_qg

!            ! Rain number concentration
!            call GSI_BundleGetPointer ( GSI_MetGuess_Bundle(it), 'qnr',ges_qnr_it,istatus )
!            if(istatus.eq.0) ges_qnr_it = g_qnr

            !
            ! Clear ges transfer fields
            !
   
            call DeallocateCloudFields ( )

         endif


     enddo

  end subroutine
!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: ReadAtmBRAMS
!
! !DESCRIPTON: rotine to read atmospheric guess fields
!
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine ReadAtmBRAMS( mype, it )
!
! !USES:
!
      ! Model read 
      use BramsMod, only : BramsFile
      use ModConstants, only: r4, r8

      ! GSI grid informations
      use gridmod, only : grd => grd_a

      ! GSI file ges information
      use guess_grids, only: ifilesig 

      !GSI kinds
      use kinds, only : r_kind

      !MPI
      use mpimod, only : NPe


      implicit none
!
! !INPUT PARAMETERS:
!
      integer, intent(in   ) :: mype ! mpi task id
      integer, intent(in   ) :: it   ! gues field number
!
! !PARAMTERS:
!

      integer, parameter                     :: atmNVars = 7
      character(len=40), dimension(atmNVars) :: atmVarName = [                             &
                                                              'TOPT                      ',& !  1 - tropography (m)
                                                              'SFC_PRESS                 ',& !  2 - surface pressure (mb)
                                                              'VTMP                      ',& !  3 - pegar THETA + PI + RV(mixing Ratio)
                                                              'UMES                      ',& !  4 - pegar RV(mixing ratio)
                                                              'UP                        ',& !  5 - u-wind (m/s)
                                                              'VP                        ',& !  6 - v-wind (m/s)
                                                              'PI                        ' & !  7 - pressure (exner function)(mb)
                                                              ]
  

! !REVISION HISTORY:
!
!   26 Mar 2019 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC      
      character(len=64), parameter :: myname_=trim(myname)//' :: ReadAtmBRAMS( )'

      type(BramsFile) :: Brams  ! BRAMS files data type

      real(r4), allocatable :: grid(:,:)
      real(r_kind)          :: work(grd%itotsub)
!      real(r_kind), allocatable  :: work(:)

      !
      ! Auxiliary Variables
      !

      integer :: istat

      integer :: WVar(NPe)
      integer :: WLev(NPe)


      integer :: Nfields
      integer :: nlon
      integer :: nlat
      integer :: nlevs
      integer :: VarNlevs

      integer :: i, j, k, ij
      integer :: ii,jj,kk
      integer :: icount
      integer :: ivar
      integer :: iret
      integer :: ilev
      integer :: iPe
      character(len=20) :: FileGuess
      character(len=20) :: FileHeader

      integer, parameter :: MyPe_Out = 0

      
      !
      ! Define files names
      ! vfm => sequential binary ieee file (grid points field)
      ! head => ascii file with field information
      !

      !
      ! Files to read and get first-guess
      ! !! hardwire files! not a best way !!
      !

#ifdef DEBUG
      WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif


      write(FileGuess,'("BRAMS.vfm.",I2.2)') ifilesig(it)
      write(FileHeader,'("BRAMS.head.",I2.2)') ifilesig(it)

      !
      ! All tasks open and read header
      !

      iret = brams%BrOpen(FileHeader, FileGuess)

      if (iret.ne.0)then

         write(stdout,'(2(A,1x),I4,1x,A)')trim(myname_),'Problem to open/read BRAMS files at ',MyPe, 'rank!'

         stop
      endif

      !
      ! consistency check
      !

      nlon  = brams%GetDim('nnxp')
      nlat  = brams%GetDim('nnyp')
      nlevs = brams%GetDim('nnzp')

      if (                            &
           (nlon .ne. grd%nlon .or.   &
            nlat .ne. grd%nlat .or.   &
            nlevs.ne. grd%nsig) .and. &
            MyPe .eq. 0               &
         )then

         write (6,'(3(1x,A))') trim(myname_),':: *** ERROR *** reading ', trim(FileGuess)
         write (6,'(2(1x,A,1x,I4))')'<nlon>', nlon,'<->.',grd%nlon
         write (6,'(2(1x,A,1x,I4))')'<nlat>', nlat,'<->.',grd%nlat
         write (6,'(2(1x,A,1x,I4))')'<nlev>',nlevs,'<->.',grd%nsig

         stop
      endif


      !
      ! Process guess atmospheric fields according to type of input file. BRAMS files are
      ! are in a rotated stereographical polar projection and need to be transformed 
      ! to a regular la/lon grid.
      ! Once on the grid, fields need to be scattered from the full domain to 
      ! sub-domains.
      !
      NFields = (atmNVars-2)*nlevs + 2
      iCount  = 1
      ivar    = 1
      do ivar = 1, atmNVars

         VarNlevs = brams%getvarnlevels(trim(atmvarname(ivar)),iret)

         if(iret.ne.0)then
            write(*,*)trim(myname_),':: *** ERROR *** '
            write(*,*)'Requested field not found', trim(atmVarName(ivar))
            stop
         endif

         do ilev = 1, VarNLevs
         
!            i  = ceiling(float(icount)/float(NPe))
!            iPe = ( ( icount + NPe ) - NPe * i ) - 1
            iPe = iCount - 1

            !----------------------------------------
            ! This is used to scatter between all Pe's
            ! Como cada Pe ira pegar um campo
            ! em um nivel, estas duas variaveis dirao
            ! qual variavel é lida em cada Pe
            !
            wvar(iCount) = ivar ! What Variable
            wlev(iCount) = ilev ! What Level
            !----------------------------------------

            if ( iPe .eq. MyPe )then

               allocate(grid(nlon, nlat) )

               call brams%GetField(grid, trim(atmVarName(ivar)), ilev)

               !
               ! some necessary adjusts to Specific Humidity field
               !

               if (trim(atmVarName(ivar)) .eq. 'UMES' .or. trim(atmVarName(ivar)) .eq. 'Q2MT' )then
                  do j=1,nlat
                     do i=1,nlon
!                        grid(i,j) = grid(i,j) * 0.001_r4 ! convert from g/kg to kg/kg
                        grid(i,j) = MAX(1.0e-12_r8,grid(i,j))
                     enddo
                  enddo
               endif

               ! Convert Surface Pressure from hPa to kPa (mb to cb). 
               ! NOTA: GSI use ps in millibar but this conversion
               !       is made by internals subroutines. So at this
               !       point ps need be in centibar.

               if(trim(atmVarName(iVar)) .eq. 'SFC_PRESS' .or. trim(atmVarName(iVar)) .eq. 'PI' )then
                  grid = grid * 0.1_r4
               endif

               !
               ! reorganize grid to be used by GSI
               !
               !  - reorder the output array so that it is a one-dimensional 
               !    array read in an order consistently with that assumed for total
               !    domain gsi grids.
               !

!               allocate (work(grd%itotsub))
!               allocate(work(nlon*nlat))

               call ReordField(grid,work)

               deallocate(grid)

            endif

            if ( iCount .eq. NPe .or. ( ivar .eq. atmNVars .and. ilev .eq. NLevs) )then

               !
               ! Transfer contents of 2-d array global to 3-d subdomain array
               ! Fields are scattered from the full domain to sub-domains.
               !
               
               call GenReload( work,  & ! Input Field read by Pe
                               WVar,  & !       Position in Var list
                               WLev,  & !       Level
                               icount,& !       MaxCount until here
                               g_z,   & ! OutPut Topography
                               g_ps,  & !        Surface Pressure
                               g_tv,  & !        Virtural Temperature
                               g_q,   & !        Specific Humidity
                               g_u,   & !        Zonal Wind
                               g_v,   & !        Meridional Wind
                               g_pi   & !        Pressure (exner function)
                               )

               !
               ! after transfer and scatter reset all counters
               !
               
               WVar   = 0
               WLev   = 0
               icount = 0

!               deallocate(work)

            endif

            icount = icount + 1
         enddo

      enddo

      iret = brams%BrClose( )

      if(iret.ne.0)then
         write(*,*)trim(myname_),':: *** WARNING *** '
         write(*,*)'Some problem to Release bramsFile memory'
         stop
      endif

      return
   end subroutine
!EOC
!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: ReadAtmBRAMS
!
! !DESCRIPTON: rotine to read atmospheric guess fields
!
!             
!                 
!\\
!\\
! !INTERFACE:
!
   subroutine ReadCloudBRAMS( mype, it )
!
! !USES:
!
      ! Model read 
      use BramsMod, only : BramsFile
      use ModConstants, only: r4, r8

      ! GSI grid informations
      use gridmod, only: grd => grd_a

      ! GSI file ges information
      use guess_grids, only: ifilesig 

      !GSI kinds
      use kinds, only : r_kind

      !MPI
      use mpimod, only : NPe


      implicit none
!
! !INPUT PARAMETERS:
!
      integer, intent(in   ) :: mype ! mpi task id
      integer, intent(in   ) :: it   ! gues field number
!
! !PARAMTERS:
!

      integer, parameter                       :: cloudNVars = 5
      character(len=40), dimension(cloudNVars) :: cloudVarName = [                             &
                                                                  'CLOUD                     ',& !  1 - Cloud Mixing Ratio (g/kg) -> ql - cloud_liquid
                                                                  'ICE                       ',& !  2 - ice Mixing Ratio (g/kg) -> qi - cloud_ice
                                                                  'RAIN                      ',& !  3 - Rain Mixing Ratio (g/kg) -> qr - rain
                                                                  'SNOW                      ',& !  4 - Snow Mixing Ratio (g/kg) -> qs - snow
                                                                  'GRAUPEL                   ' & !  5 - Graupel Mixing Ratio (g/kg) -> qg - graupel
                                                                  ] 

  

! !REVISION HISTORY:
!
!   26 Mar 2019 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC      
      character(len=64), parameter :: myname_=trim(myname)//' :: ReadCloudBRAMS( )'

      type(BramsFile) :: Brams  ! BRAMS files data type

      real(r4), pointer :: grid(:,:)
      real(r_kind) :: work(grd%itotsub)

      !
      ! Auxiliary Variables
      !

      integer :: istat

      integer :: WVar(NPe)
      integer :: WLev(NPe)


      integer :: Nfields
      integer :: nlon
      integer :: nlat
      integer :: nlevs
      integer :: VarNlevs

      integer :: i, j, k, ij
      integer :: ii,jj,kk
      integer :: icount
      integer :: ivar
      integer :: iret
      integer :: ilev
      integer :: iPe
      character(len=20) :: FileGuess
      character(len=20) :: FileHeader

      integer, parameter :: MyPe_Out = 0

      
      !
      ! Define files names
      ! fct => sequential binary ieee file (spectral and grid points field)
      ! dir => ascii file with field information
      !

      !
      ! Files to read and get first-guess
      ! !! hardwire files! not a best way !!
      !

#ifdef DEBUG
      WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      write(FileGuess,'("BRAMS.vfm.",I2.2)') ifilesig(it)
      write(FileHeader,'("BRAMS.head.",I2.2)') ifilesig(it)

      !
      ! All tasks open and read header
      !

      iret = brams%BrOpen(FileHeader, FileGuess)

      if (iret.ne.0)then

         write(stdout,'(2(A,1x),I4,1x,A)')trim(myname_),'Problem to open/read BRAMS files at ',MyPe, 'rank!'

         stop
      endif

      !
      ! consistency check
      !

      nlon  = brams%GetDim('nnxp')
      nlat  = brams%GetDim('nnyp')
      nlevs = brams%GetDim('nnzp')

      if (                             & 
           (nlon  .ne. grd%nlon .or.   &
            nlat  .ne. grd%nlat .or.   &
            nlevs .ne. grd%nsig) .and. &
            MyPe  .eq. 0               &
         )then

         write (6,'(3(1x,A))') trim(myname_),':: *** ERROR *** reading ', trim(FileGuess)
         write (6,'(2(1x,A,1x,I4))')'<nlon>', nlon,'<->.',grd%nlon
         write (6,'(2(1x,A,1x,I4))')'<nlat>', nlat,'<->.',grd%nlat
         write (6,'(2(1x,A,1x,I4))')'<nlev>',nlevs,'<->.',grd%nsig

         stop
      endif


      !
      ! Process guess atmospheric fields according to type of input file. BRAMS files are
      ! are in a rotated stereographical polar projection and need to be transformed 
      ! to a regular la/lon grid.
      ! Once on the grid, fields need to be scattered from the full domain to 
      ! sub-domains.
      !
      NFields = CloudNVars*nlevs
      iCount  = 1
      ivar    = 1

      do ivar = 1, CloudNVars

         VarNLevs = brams%GetVarNlevels(trim(CloudVarName(ivar)),iret)

         if(iret.ne.0)then
            write(*,*)trim(myname_),':: *** ERROR *** '
            write(*,*)'Requested field not found', trim(CloudVarName(ivar))
            stop
         endif

         do ilev = 1, VarNLevs
         
!            i  = ceiling(float(icount)/float(NPe))
!            iPe = ( ( icount + NPe ) - NPe * i ) - 1
            iPe = iCount - 1

            !----------------------------------------
            ! This is used to scatter between all Pe's
            ! Como cada Pe ira pegar um campo
            ! em um nivel, estas duas variaveis dirao
            ! qual variavel é lida em cada Pe
            !
            wvar(iCount) = ivar ! What Variable
            wlev(iCount) = ilev ! What Level
            !----------------------------------------


            if ( iPe .eq. MyPe )then
               allocate(grid(nlon, nlat) )

               call brams%GetField(grid, trim(cloudVarName(ivar)), ilev)

               !
               ! reorganize grid to be used by GSI
               !
               !  - reorder the output array so that it is a one-dimensional 
               !    array read in an order consistently with that assumed for total
               !    domain gsi grids.
               !

!               allocate (work(grd%itotsub))

               call ReordField(grid,work)

               deallocate(grid)

            endif

            if ( iCount .eq. NPe .or. ( ivar .eq. CloudNVars .and. ilev .eq. NLevs) )then

               !
               ! Transfer contents of 2-d array global to 3-d subdomain array
               ! Fields are scattered from the full domain to sub-domains.
               !
               
               call CloudReload( work,  & ! Input Field read by Pe
                                 WVar,  & !       Position in Var list
                                 WLev,  & !       Level
                                 icount,& !       MaxCount until here
                                 g_ql,  & ! OutPut Cloud Mixing Ratio
                                 g_qi,  & !        Ice Mixing Ratio
                                 g_qr,  & !        Rain Mixing Ratio
                                 g_qs,  & !        Snow Mixing Ratio
                                 g_qg   & !        Graupel Mixing Ratio
                               )

               !
               ! after transfer and scatter reset all counters
               !
               
               WVar   = 0
               WLev   = 0
               icount = 0

!               deallocate(work)

            endif

            icount = icount + 1
         enddo

      enddo

      iret = brams%BrClose( )
      if(iret.ne.0)then
         write(*,*)trim(myname_),':: *** WARNING *** '
         write(*,*)'Some problem to Release bramsFile memory'
         stop
      endif

      return
   end subroutine

!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: AllocateAtmFields
!
! !DESCRIPTON: rotine to allocate ges transfer fields
!
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine AllocateAtmFields( )

      use gridmod, only: lat2, lon2, nsig
      implicit none

! !REVISION HISTORY:
!
!   26 Mar 2019 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC  
      character(len=64), parameter :: myname_=trim(myname)//' :: AllocateAtmFields( )'


#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      allocate( g_z  (lat2, lon2) )
      allocate( g_ps (lat2, lon2) )
      allocate( g_tv (lat2, lon2, nsig) )
      allocate( g_u  (lat2, lon2, nsig) )
      allocate( g_v  (lat2, lon2, nsig) )
      allocate( g_q  (lat2, lon2, nsig) )
      allocate( g_pi (lat2, lon2, nsig) )

   end subroutine

!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: DeallocateAtmFields
!
! !DESCRIPTON: rotine to Deallocate ges transfer fields
!
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine DeallocateAtmFields( )

      implicit none

! !REVISION HISTORY:
!
!   26 Mar 2019 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC      
      character(len=64), parameter :: myname_=trim(myname)//' :: DeAllocateAtmFields( )'


#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      deallocate( g_z  )
      deallocate( g_ps )
      deallocate( g_tv )
      deallocate( g_u  )
      deallocate( g_v  )
      deallocate( g_q  )
      deallocate( g_pi )


   end subroutine

!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: AllocateCloudFields
!
! !DESCRIPTON: rotine to allocate ges Cloud transfer fields
!
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine AllocateCloudFields( )

      use gridmod, only: lat2, lon2, nsig
      implicit none

! !REVISION HISTORY:
!
!   16 Dec 2019 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC  
      character(len=64), parameter :: myname_=trim(myname)//' :: AllocateCloudFields( )'


#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      allocate( g_ql  (lat2, lon2, nsig) )
      allocate( g_qi  (lat2, lon2, nsig) )
      allocate( g_qr  (lat2, lon2, nsig) )
      allocate( g_qs  (lat2, lon2, nsig) )
      allocate( g_qg  (lat2, lon2, nsig) )
!      allocate( g_qnr (lat2, lon2, nsig) )

   end subroutine

!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: DeallocateCloudFields
!
! !DESCRIPTON: rotine to Deallocate ges cloud transfer fields
!
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine DeallocateCloudFields( )

      implicit none

! !REVISION HISTORY:
!
!   16 Dec 2019 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC      
      character(len=64), parameter :: myname_=trim(myname)//' :: DeAllocateCloudFields( )'


#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      deallocate( g_ql  )
      deallocate( g_qi  )
      deallocate( g_qr  )
      deallocate( g_qs  )
      deallocate( g_qg  )
!      deallocate( g_qnr )

   end subroutine

!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: AllocateAtmFields
!
! !DESCRIPTON: rotine to allocate ges transfer fields
!
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine AllocateAnlFields( )

      use gridmod, only: lat2, lon2, nsig
      implicit none

! !REVISION HISTORY:
!
!   26 Mar 2019 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC  
      character(len=64), parameter :: myname_=trim(myname)//' :: AllocateAtmFields( )'


#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      allocate( g_z  (lat2, lon2) )

      allocate( g_u  (lat2, lon2, nsig) )
      allocate( g_v  (lat2, lon2, nsig) )
      allocate( g_tk (lat2, lon2, nsig) )
      allocate( g_rH (lat2, lon2, nsig) )
      allocate( g_gh (lat2, lon2, nsig) )


      allocate( g_ps   (lat2, lon2) )
      allocate( g_tv   (lat2, lon2, nsig) )
      allocate( g_q    (lat2, lon2, nsig) )
      allocate( g_pi   (lat2, lon2, nsig) )
      allocate( g_prsi (lat2, lon2, nsig) )
      allocate( g_prsl (lat2, lon2, nsig) )

   end subroutine

!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: DeallocateAtmFields
!
! !DESCRIPTON: rotine to Deallocate ges transfer fields
!
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine DeallocateAnlFields( )

      implicit none

! !REVISION HISTORY:
!
!   26 Mar 2019 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC      
      character(len=64), parameter :: myname_=trim(myname)//' :: DeAllocateAtmFields( )'


#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      deallocate( g_z )

      deallocate( g_u  )
      deallocate( g_v  )
      deallocate( g_tk )
      deallocate( g_rH )
      deallocate( g_gh )

      deallocate( g_ps  )
      deallocate( g_tv  )
      deallocate( g_q   )
      deallocate( g_pi  )
      deallocate( g_prsi)
      deallocate( g_prsl)

   end subroutine

!EOC

!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: GenReload
!
! !DESCRIPTON: rotine to transfer contents of 2-d array global to 3-d subdomain array
!              and scatter from the full domain to sub-domains.
!             
!                 
!\\
!\\
! !INTERFACE:
!

subroutine GenReload( work, wvar, wlev, icount, &
                      g_z, g_ps, g_tv, g_q, g_u, g_v, g_pi)

! !USES:

  use kinds, only: r_kind,i_kind
  use mpimod, only: npe,mpi_comm_world,ierror,mpi_rtype
  use gridmod, only: grd => grd_a
  implicit none

! !INPUT PARAMETERS:

!  real(r_kind),dimension(grd%itotsub)   ,intent(in   ) :: work
  real(r_kind),dimension(:)             ,intent(in   ) :: work
  integer(i_kind),dimension(npe)        ,intent(in   ) :: wvar
  integer(i_kind),dimension(npe)        ,intent(in   ) :: wlev
  integer(i_kind)                       ,intent(in   ) :: icount

! !OUTPUT PARAMETERS:

  real(r_kind),dimension(:,:),  intent(  out) :: g_z
  real(r_kind),dimension(:,:),  intent(  out) :: g_ps
  real(r_kind),dimension(:,:,:),intent(  out) :: g_tv
  real(r_kind),dimension(:,:,:),intent(  out) :: g_q
  real(r_kind),dimension(:,:,:),intent(  out) :: g_u
  real(r_kind),dimension(:,:,:),intent(  out) :: g_v
  real(r_kind),dimension(:,:,:),intent(  out) :: g_pi


! !REVISION HISTORY:
!   2004-05-14  treadon
!   2019-12-16  de Mattos  - adapt to BRAMS model
!
! !REMARKS:
!
!   language: f90
!
!
!EOP
!-------------------------------------------------------------------------
!BOC
      character(len=64), parameter :: myname_=trim(myname)//' :: GenReload( )'

  integer(i_kind) :: i,j,k,ij
  integer(i_kind) :: Pe
  integer(i_kind) :: var
  integer(i_kind) :: lev
  real(r_kind),dimension(grd%lat2*grd%lon2,npe):: sub


#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

  call mpi_alltoallv(work,          & 
                     grd%ijn_s,     &
                     grd%displs_s,  &
                     mpi_rtype,     &
                     sub,           &
                     grd%irc_s,     &
                     grd%ird_s,     &
                     mpi_rtype,     &
                     mpi_comm_world,&
                     ierror)

!$omp parallel do  schedule(dynamic,1) private(Pe,i,j,ij,var,lev)
  do Pe = 1, icount

     var = wvar(Pe)
     lev = wlev(Pe)

     select case ( var )

        case (1) ! Topography

           ij=0
           do j=1,grd%lon2
              do i=1,grd%lat2
                 ij=ij+1
                 g_z(i,j) = sub(ij,Pe)
              end do
           end do

        case (2) ! Surface Pressure

           ij=0
           do j=1,grd%lon2
              do i=1,grd%lat2
                 ij=ij+1
                 g_ps(i,j) = sub(ij,Pe)
              end do
           end do

        case (3) ! Virtual Temperature

           ij=0
           do j=1,grd%lon2
              do i=1,grd%lat2
                 ij=ij+1
                 g_tv(i,j,lev) = sub(ij,Pe)
              end do
           end do

        case (4) ! Specific Humidity

           ij=0
           do j=1,grd%lon2
              do i=1,grd%lat2
                 ij=ij+1
                 g_q(i,j,lev) = sub(ij,Pe)
              end do
           end do

        case (5) ! zonal wind (UP)

          ij=0
          do j=1,grd%lon2
             do i=1,grd%lat2
                ij=ij+1
                g_u(i,j,lev) = sub(ij,Pe)
             end do
          end do

        case (6) ! meridional wind (VP)

          ij=0
          do j=1,grd%lon2
             do i=1,grd%lat2
                ij=ij+1
                g_v(i,j,lev) = sub(ij,Pe)
             end do
          end do

        case (7) ! Pressure [exner function] (PI)

          ij=0
          do j=1,grd%lon2
             do i=1,grd%lat2
                ij=ij+1
                g_pi(i,j,lev) = sub(ij,Pe)
             end do
          end do


     end select
  enddo


  return
end subroutine GenReload
!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: CloudReload
!
! !DESCRIPTON: rotine to transfer contents of 2-d array global to 3-d subdomain array
!              and scatter from the full domain to sub-domains.
!             
!                 
!\\
!\\
! !INTERFACE:
!

subroutine CloudReload( work, wvar, wlev, icount, &
                      g_ql, g_qi, g_qr, g_qs, g_qg)

! !USES:

  use kinds, only: r_kind,i_kind
  use mpimod, only: npe,mpi_comm_world,ierror,mpi_rtype
  use gridmod, only: grd => grd_a
  implicit none

! !INPUT PARAMETERS:

  real(r_kind),dimension(grd%itotsub)   ,intent(in   ) :: work
  integer(i_kind),dimension(npe)        ,intent(in   ) :: wvar
  integer(i_kind),dimension(npe)        ,intent(in   ) :: wlev
  integer(i_kind)                       ,intent(in   ) :: icount

! !OUTPUT PARAMETERS:

  real(r_kind),dimension(:,:,:),intent(  out) :: g_ql
  real(r_kind),dimension(:,:,:),intent(  out) :: g_qi
  real(r_kind),dimension(:,:,:),intent(  out) :: g_qr
  real(r_kind),dimension(:,:,:),intent(  out) :: g_qs
  real(r_kind),dimension(:,:,:),intent(  out) :: g_qg


! !REVISION HISTORY:
!   2004-05-14  treadon
!   2019-12-16  de Mattos  - adapt to BRAMS model
!
! !REMARKS:
!
!   language: f90
!
!
!EOP
!-------------------------------------------------------------------------
!BOC
      character(len=64), parameter :: myname_=trim(myname)//' :: GenReload( )'

  integer(i_kind) :: i,j,k,ij
  integer(i_kind) :: Pe
  integer(i_kind) :: var
  integer(i_kind) :: lev
  real(r_kind),dimension(grd%lat2*grd%lon2,npe):: sub


#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

  call mpi_alltoallv(work,          & 
                     grd%ijn_s,     &
                     grd%displs_s,  &
                     mpi_rtype,     &
                     sub,           &
                     grd%irc_s,     &
                     grd%ird_s,     &
                     mpi_rtype,     &
                     mpi_comm_world,&
                     ierror)

!$omp parallel do  schedule(dynamic,1) private(Pe,i,j,ij,var,lev)
  do Pe = 1, icount

     var = wvar(Pe)
     lev = wlev(Pe)

     select case ( var )

        case (1) ! Cloud Mixing Ratio

           ij=0
           do j=1,grd%lon2
              do i=1,grd%lat2
                 ij=ij+1
                 g_ql(i,j,lev) = sub(ij,Pe)
              end do
           end do

        case (2) ! Ice Mixing Ratio

           ij=0
           do j=1,grd%lon2
              do i=1,grd%lat2
                 ij=ij+1
                 g_qi(i,j,lev) = sub(ij,Pe)
              end do
           end do

        case (3) ! Rain Mixing Ratio

          ij=0
          do j=1,grd%lon2
             do i=1,grd%lat2
                ij=ij+1
                g_qr(i,j,lev) = sub(ij,Pe)
             end do
          end do

        case (4) ! Snow Mixing Ratio

          ij=0
          do j=1,grd%lon2
             do i=1,grd%lat2
                ij=ij+1
                g_qs(i,j,lev) = sub(ij,Pe)
             end do
          end do

        case (5) ! Graupel Mixing Ratio

          ij=0
          do j=1,grd%lon2
             do i=1,grd%lat2
                ij=ij+1
                g_qg(i,j,lev) = sub(ij,Pe)
             end do
          end do


     end select
  enddo


  return
end subroutine CloudReload
!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: ReordField
!
! !DESCRIPTION: This routine reorder the output array so that it is a
!               one-dimensional array read in an order consisten with 
!               that assumed for total domain gsi grids.
!
!               The assumed order for the input grid is longitude as
!               the first dimension with array index increasing from
!               east to west.  The second dimension is latitude with
!               the index increasing from north to south.  This ordering
!               differs from that used in the GSI.
!
!               The GSI ordering is latitude first with the index
!               increasing from south to north.  The second dimension is
!               longitude with the index increasing from east to west.
!
!               Thus, the code below also rearranges the indexing and
!               order of the dimensions to make the output grid
!               consistent with that which is expected in the rest of
!               gsi.
!\\
!\\
! !INTERFACE:
!

 subroutine ReordField(GridIn,GridOut)

! !USES:

   use ModConstants, only : r4
   use kinds, only: r_kind,i_kind
   use gridmod, only: grd => grd_a

   implicit none

! !INPUT PARAMETERS:

   real(r4),            intent(in   ) :: GridIn (:,:) ! input grid <imax,jmax> (nlon,nlat)
   real(r_kind),        intent(  out) :: GridOut( : ) ! output grid <itotsub>

!
!
! !REVISION HISTORY:
!   2019-12-16  J. G. Z. de Mattos
!
! !REMARKS:
!   language: f90
!
!
!EOP
!-------------------------------------------------------------------------
!BOC
      character(len=64), parameter :: myname_=trim(myname)//' :: ReordField( )'


!  Declare local variables
   real(r_kind), allocatable :: swap(:,:)
   integer(i_kind) :: i, j, k
   integer(i_kind) :: nlat, nlon

#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

   nlat = grd%nlat
   nlon = grd%nlon

!  Transfer local work array to output grid
!   do k = 1, grd%itotsub
!
!      j          = grd%nlat-grd%ltosi_s(k)
!      i          = grd%ltosj_s(k)
!      if(i.ge.1 .and. i .le. ipt .and. j .ge. 1 .and. j .le. jpt)then
!         GridOut(k) = real(GridIn(i,j),r_kind)
!      else
!         GridOut(k) = 0.0
!      endif
!
!   end do

   allocate(swap(nlat,nlon))
   call SFCTrans(GridIn, swap)

!   do j = 1, nlat
!      do i = 1, nlon
!         swap(j,i) = GridIn(i,j)
!      enddo
!   enddo

   GridOut = 0.0
   do k=1,grd%itotsub
      i = grd%ltosi_s(k) ! latitude
      j = grd%ltosj_s(k) ! longitude
      if(i.ge.1 .and. i .le. nlat .and. j .ge. 1 .and. j .le. nlon)then
!         GridOut(k) = GridIn(j,i)
         GridOut(k) = swap(i,j)
      endif
   enddo

   deallocate(swap)

   return
 end subroutine ReordField

!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !FUNCTION: GenReadBRAMSSFC
!
! !DESCRIPTON: subroutine to get BRAMS sfc fields and tranfer to GSI
!
!             
!                 
!\\
!\\
! !INTERFACE:
!
   subroutine readBramsSfc(  iope,  mype,  &
                             fact10,       &! 10-meter wind factor
                             td0,          &! surface temperature [k]
                             sheleg,       &! snow depth [m]
                             vtype,        &! vegetation type  [1-13]
                             vcover,       &! vegetation cover [-]
                             stype,        &! soil type [1-9]
                             tg0,          &! surface soil temperature [k]
                             w0,           &! surface soil moisture [fraction]
                             lsimsk,       &! land sea mask
                             Z0,           &! surface roughness length [cm]
                             use_sfc_any   &
                             )  
!
!  !USES:
!
      ! GSI guess number actual count of sfc in-cache time slots 
      use guess_grids, only: nfldsfc

      ! GSI grid - SFC grid sizes
      use gridmod, only: nlat_sfc, nlon_sfc

      ! GSI kinds
      use kinds, only: r_kind,i_kind

      ! MPI GSI 
      use mpimod, only: mpi_itype,mpi_rtype,mpi_comm_world

      implicit none

!
! !INPUT PARAMETERS:
!
    integer(i_kind), intent(in   ) :: iope        ! mpi task handling i/o
    integer(i_kind), intent(in   ) :: mype        ! mpi task id
    logical,         intent(in   ) :: use_sfc_any !

!
! !OUTPUT PARAMETERS:
!

    integer(i_kind), dimension(:,:)  , intent(inout) :: lsimsk ! land sea mask 
    real(r_kind)   , dimension(:,:)  , intent(inout) :: vtype  ! vegetation type  [1-13]
    real(r_kind)   , dimension(:,:)  , intent(inout) :: stype  ! soil type [1-9]
    real(r_kind)   , dimension(:,:,:), intent(inout) :: fact10 ! 10-meter wind factor
    real(r_kind)   , dimension(:,:,:), intent(inout) :: td0    ! surface temperature [k]
    real(r_kind)   , dimension(:,:,:), intent(inout) :: sheleg ! snow depth [m]
    real(r_kind)   , dimension(:,:,:), intent(inout) :: vcover ! vegetation cover [-]
    real(r_kind)   , dimension(:,:,:), intent(inout) :: tg0    ! surface soil temperature [k]
    real(r_kind)   , dimension(:,:,:), intent(inout) :: w0     ! surface soil moisture [fraction]
    real(r_kind)   , dimension(:,:,:), intent(inout) :: Z0     ! surface roughness length [cm]   


! !REVISION HISTORY:
!
!   06 Jul 2016 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

      character(len=100), parameter :: myname_=trim(myname)//':: Read_BRAMSSFC( )'

      !
      ! Auxiliary vars
      !

      integer :: istat
      integer :: it
      integer :: npts
      integer :: nptsall


#ifdef DEBUG
      WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      npts    = nlat_sfc * nlon_sfc
      nptsall = npts * nfldsfc


      !
      ! Loop over all time sfc gues fields
      !
      do it = 1, nfldsfc

         !
         ! Read BRAMS fields
         !

         call GenReadBRAMSSFC  (   iope,  mype,     &
                                 lsimsk (:,:),    &
                                 vtype  (:,:),    &
                                 stype  (:,:),    &
                                 fact10 (:,:,it), &
                                 td0    (:,:,it), &
                                 sheleg (:,:,it), &
                                 vcover (:,:,it), &
                                 tg0    (:,:,it), &
                                 w0     (:,:,it), &
                                 Z0     (:,:,it), &
                                 use_sfc_any,     &
                                 it               &
                             )

      enddo

      !
      ! Load onto all processors
      !


      call mpi_bcast( lsimsk,    npts, mpi_itype, iope, mpi_comm_world, istat )
      call mpi_bcast(    td0, nptsall, mpi_rtype, iope, mpi_comm_world, istat )
      call mpi_bcast( fact10, nptsall, mpi_rtype, iope, mpi_comm_world, istat )
      call mpi_bcast( sheleg, nptsall, mpi_rtype, iope, mpi_comm_world, istat )
      call mpi_bcast(     Z0, nptsall, mpi_rtype, iope, mpi_comm_world, istat )

      if(use_sfc_any)then

         call mpi_bcast(  vtype,    npts, mpi_rtype, iope, mpi_comm_world, istat )
         call mpi_bcast(  stype,    npts, mpi_rtype, iope, mpi_comm_world, istat )
         call mpi_bcast( vcover, nptsall, mpi_rtype, iope, mpi_comm_world, istat )
         call mpi_bcast(    tg0, nptsall, mpi_rtype, iope, mpi_comm_world, istat )
         call mpi_bcast(     w0, nptsall, mpi_rtype, iope, mpi_comm_world, istat )

      end if

   end subroutine

!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !FUNCTION: GenReadBRAMSSFC
!
! !DESCRIPTON: subroutine to read sfc files from BRAMS fct files
!
!             
!                 
!\\
!\\
! !INTERFACE:
!
   subroutine GenReadBRAMSSFC( iope,  mype,               &
                               lsimsk, vtype, stype,      & 
                               fact10,td0, sheleg,vcover, &
                               tg0, w0, Z0,               &
                               use_sfc_any, it )
       ! Model read 
       use BramsMod, only: BramsFile
       use ModConstants, only:  r4, r8
!       use sigio_BRAMSMod, only : BRAMS_GetAvailUnit
       ! GSI kinds
       use kinds, only: r_kind,i_kind

       ! GSI guess files names
       use guess_grids, only: ifilesfc

      ! GSI grid informations
      use gridmod, only : grd=>grd_a

       implicit none

!
! !INPUT PARAMETERS:
!
    integer(i_kind), intent(in   ) :: iope        ! mpi task handling i/o
    integer(i_kind), intent(in   ) :: mype        ! mpi task id
    logical,         intent(in   ) :: use_sfc_any !
    integer(i_kind), intent(in   ) :: it
!
! !OUTPUT PARAMETERS:
!

    integer(i_kind), dimension(:,:), intent(inout) :: lsimsk ! land sea mask 
    real(r_kind)   , dimension(:,:), intent(inout) :: vtype  ! vegetation type  [1-13]
    real(r_kind)   , dimension(:,:), intent(inout) :: stype  ! soil type [1-9]
    real(r_kind)   , dimension(:,:), intent(inout) :: fact10 ! 10-meter wind factor
    real(r_kind)   , dimension(:,:), intent(inout) :: td0    ! surface temperature [k]
    real(r_kind)   , dimension(:,:), intent(inout) :: sheleg ! snow depth [m]
    real(r_kind)   , dimension(:,:), intent(inout) :: vcover ! vegetation cover [-]
    real(r_kind)   , dimension(:,:), intent(inout) :: tg0    ! surface soil temperature [k]
    real(r_kind)   , dimension(:,:), intent(inout) :: w0     ! surface soil moisture [fraction]
    real(r_kind)   , dimension(:,:), intent(inout) :: Z0     ! surface roughness length [cm]

!
! !PARAMTERS:
!

       integer, parameter                  :: NVars = 10
       character(len=40), dimension(NVars) :: VName = [                             &
                                                       'PATCH_AREA                ',& !  1- land sea ice mask [0,1,2]
                                                       'LEAF_CLASS                ',& !  2- Vegetation Type [ type numbers ]
                                                       'SOIL_TEXT                 ',& !  3- soil type [1-9]
                                                       'FACT10M                   ',& !  4- 10-meter wind factor <- definir
                                                       'VEG_TEMP                  ',& !  5- Surface Temperature [k]
                                                       'snowdepthj                ',& !  6- snow depth [m]
                                                       'VEG_FRACAREA              ',& !  7- vegetation cover [-]
                                                       'SOIL_TEMPJ                ',& !  8- surface soil temperature [k]
                                                       'SOIL_WATER                ',& !  9- surface soil moisture [fraction]
                                                       'PATCH_ROUGH               ' & ! 10- surface roughness length [cm]
                                                      ]


! !REVISION HISTORY:
!
!   13 Jan 2020 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      character(len=100), parameter :: myname_=trim(myname)//' :: GenReadSFCBRAMS( )'

      type(BRAMSFile) :: BRAMS  ! BRAMS files data type

      real(r_kind), allocatable :: TmpGrid(:,:)

      integer :: nlon
      integer :: nlat
      integer :: npatch
      integer :: nzg

      integer :: ivar
      integer :: istat
      integer :: iret
      integer :: i, j

      character(len=20) :: FileGuess
      character(len=20) :: FileHeader



      !
      ! Auxiliary Variables
      !

      real(r4), dimension( :,: ), allocatable :: u
      real(r4), dimension( :,: ), allocatable :: v
      real(r4), dimension( :,: ), allocatable :: speed
      real(r4), dimension( :,: ), allocatable :: u10m
      real(r4), dimension( :,: ), allocatable :: v10m
      real(r4), dimension( :,: ), allocatable :: speed10m
      real(r4), dimension( :,: ), allocatable :: field2d
      real(r4), dimension(:,:,:), allocatable :: field3d
      real(r4), dimension(:,:,:), allocatable :: patch_area

      !
      ! Define files names
      ! vfm => sequential binary ieee file (grid points field)
      ! head => ascii file with field information
      !

      !
      ! Files to read and get first-guess
      ! !! hardwire files! not a best way !!
      !

#ifdef DEBUG
      WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      write(FileGuess,'("BRAMS.vfm.",I2.2)') ifilesfc(it)
      write(FileHeader,'("BRAMS.head.",I2.2)') ifilesfc(it)

      !
      ! All tasks open and read header
      !

      iret = brams%BrOpen(FileHeader, FileGuess)

      if (iret.ne.0)then

         write(stdout,'(2(A,1x),I4,1x,A)')trim(myname_),'Problem to open/read BRAMS files at ',MyPe, 'rank!'

         stop
      endif

      !
      ! consistency check
      !

      nLon   = brams%GetDim('nnxp')
      nLat   = brams%GetDim('nnyp')
      npatch = brams%GetDim('npatch')
      nzg    = brams%GetDim('nzg')

      if (                            &
           (nLon .ne. grd%nlon .or.   &
            nLat .ne. grd%nlat).and.  &
            MyPe .eq. 0               &
         )then

         write (6,'(3(1x,A))') trim(myname_),':: *** ERROR *** reading ', trim(FileGuess)
         write (6,'(2(1x,A,1x,I4))')'<nlon>', nlon,'<->.',grd%nlon
         write (6,'(2(1x,A,1x,I4))')'<nlat>', nlat,'<->.',grd%nlat

         stop
      endif

      !-----------------------------------------------------------------!
      ! get patch_area
      !

      allocate(patch_area(nLon, nLat, npatch))
      do i=2,npatch
         call brams%GetField(patch_area(:,:,i),'PATCH_AREA', patch = i)
      enddo
      patch_area(:,:,1) = sum (patch_area(:,:,2:npatch), dim=3)
      where(patch_area(:,:,2) .lt. 0) patch_area(:,:,1) = 0
      where(patch_area(:,:,2) .eq. undef) patch_area(:,:,1) = undef

      !
      !-----------------------------------------------------------------!
      


      do ivar = 1, NVars

         if (                   &
             ( ivar .eq. 2 .or. & ! VEGETATION TYPE [LEAF_CLASS]
               ivar .eq. 3 .or. & ! SOIL TYPE [SOIL_TEXT]
               ivar .eq. 7 .or. & ! VEGETATION COVER [VEG_FRACAREA]
               ivar .eq. 8 .or. & ! SURFACE SOIL TEMPERATURE [SOIL_TEMPJ]
               ivar .eq. 9      & ! SURFACE SOIL MOISTURE [SOIL_WATER]
             ) .and. .not. use_sfc_any &
            ) cycle

         select case ( ivar )

            case ( 1) ! Land Sea Ice Mask [1-Land,0-Sea]

               allocate(field2d(nLon, nLat))

               where(patch_area(:,:,1) .gt. 0.0) field2d = 0.0
               where(patch_area(:,:,1) .le. 0.0) field2d = 1.0
               where(patch_area(:,:,2) .le. undef) field2d = undef

               allocate(TmpGrid(nLat,nLon))
               call SFCTrans(Field2d,TmpGrid)

               lsimsk = nint(abs(TmpGrid))

               deallocate(TmpGrid)
               deallocate(field2d)

            case ( 2) ! Vegetation Type [ type numbers ]

               allocate(field2d(nLon, nLat))
               call brams%GetField(field2d,trim(VName(ivar)), patch = 2) ! patch 2 have higest area

               where(patch_area(:,:,2) .le. 0.0   ) field2d = 0.0
               where(patch_area(:,:,2) .eq. undef ) field2d = undef

               call SFCTrans(field2d,vtype)

               deallocate(field2d)

            case ( 3) ! Soil Type [type numbers]

               allocate(field2d(nLon, nLat))
               call brams%GetField(field2d,trim(VName(ivar)), patch = 2) ! patch 2 have higest area
               
               where(patch_area(:,:,2) .le. 0.0   ) field2d = 0.0
               where(patch_area(:,:,2) .eq. undef ) field2d = undef

               call SFCTrans(field2d,stype)

               deallocate(field2d)

            case ( 4) ! 10-meter Wind Factor [ - ]

               allocate(u10m(nLon,nLat))
               call brams%GetField(u10m,'U10MJ', level = 1)
               allocate(v10m(nLon,nLat))
               call brams%GetField(v10m,'V10MJ', level = 1)
               allocate(speed10m(nLon,nLat))
               speed10m = sqrt(u10m*u10m + v10m*v10m)
               where(u10m .eq. undef .or. v10m .eq. undef) speed10m = undef
   
               deallocate(u10m)
               deallocate(v10m)
   
               allocate(u(nLon,nLat))
               call brams%GetField(u,'UP', level = 1)
               allocate(v(nLon,nLat))
               call brams%GetField(v,'VP', level = 1)
               allocate(speed(nLon,nLat))
               speed = sqrt(u*u + v*v)
   
               deallocate(u)
               deallocate(v)
   
               allocate(field2d(nLon,nLat))
   
               field2d = speed10m/speed
   
               where(patch_area(:,:,2) .le. 0.0   ) field2d = 0.0
               where(patch_area(:,:,2) .eq. undef ) field2d = undef

               call SFCTrans(field2d,fact10)

               deallocate(field2d)

            case ( 5) ! Surface Temperature [ k ]

               allocate(field2d(nLon, nLat))
               call brams%GetField(field2d,trim(VName(ivar)), patch = 2) ! patch 2 have higest area
               
               where(patch_area(:,:,2) .eq. undef ) field2d = undef

               call SFCTrans(field2d,td0)

               deallocate(field2d)

            case ( 6) ! Snow Depth [ m ]

               allocate(field2d(nLon, nLat))
               call brams%GetField(field2d,trim(VName(ivar)), patch = 2) ! patch 2 have higest area
               
               where(patch_area(:,:,2) .le. 0.0   ) field2d = 0.0
               where(patch_area(:,:,2) .eq. undef ) field2d = undef

               call SFCTrans(field2d,sheleg)

               deallocate(field2d)

            case ( 7) ! Vegetation Cover [ - ]

               allocate(field2d(nLon, nLat))
               call brams%GetField(field2d,trim(VName(ivar)), patch = 2) ! patch 2 have higest area
               
               where(patch_area(:,:,2) .le. 0.0   ) field2d = 0.0
               where(patch_area(:,:,2) .eq. undef ) field2d = undef

               call SFCTrans(field2d,vcover)

               deallocate(field2d)

            case ( 8) ! Surface Soil Temperature [ k ]

               allocate(field2d(nLon, nLat))
               call brams%GetField(field2d,trim(VName(ivar)), level = nzg, patch = 2) ! patch 2 have higest area
               
               where(patch_area(:,:,2) .eq. undef ) field2d = undef

               call SFCTrans(field2d,tg0)

               deallocate(field2d)

            case ( 9) ! Surface soil moisture [ fraction ]

               allocate(field2d(nLon, nLat))
               call brams%GetField(field2d,trim(VName(ivar)), level = nzg, patch = 2) ! patch 2 have higest area
               
               where(patch_area(:,:,2) .eq. undef ) field2d = undef

               call SFCTrans(field2d,w0)

               deallocate(field2d)

            case (10) ! Roughness Length [ cm ]
               allocate(field2d(nLon, nLat))
               call brams%GetField(field2d,trim(VName(ivar)), patch = 2) ! patch 2 have higest area
               
               where(patch_area(:,:,2) .le. 0.0   ) field2d = 0.0
               where(patch_area(:,:,2) .eq. undef ) field2d = undef

               call SFCTrans(field2d,Z0)

               deallocate(field2d)

         end select

      enddo

      !
      ! Close BRAMS Files
      !

      iret = brams%BrClose( )
      if(iret.ne.0)then
         write(*,*)trim(myname_),':: *** WARNING *** '
         write(*,*)'Some problem to Release bramsFile memory'
         stop
      endif

      return
   end subroutine

!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !FUNCTION: LUAvail
!
! !DESCRIPTON: function to return next available logical unit
!
!
!             
!                 
!\\
!\\
! !INTERFACE:
!
   subroutine SFCTrans( GridIn, GridOut )
! !USES:
      use ModConstants, only : r4
      use kinds, only: r_kind, i_kind
      use constants, only: zero,one

      implicit none
! !IMPUT PARAMETERS:

      real(r4),            intent(in   ) :: GridIn (:,:) ! input grid  ->  IMax,JMax
      real(r_kind),        intent(inout) :: GridOut(:,:) ! output grid ->  JMax,IMax

! !DESCRIPTION: This routine reorder the output array so that it is a 
!               two-dimensional array read in an order consisten with 
!               that assumed for total domain gsi grids.
!
!               The assumed order for the input grid is longitude as
!               the first dimension with array index increasing from
!               east to west.  The second dimension is latitude with
!               the index increasing from north to south.  This ordering
!               differs from that used in the GSI.
!
!               The GSI ordering is latitude first with the index
!               increasing from south to north.  The second dimension is
!               longitude with the index increasing from east to west.
!
!               Thus, the code below also rearranges the indexing and
!               order of the dimensions to make the output grid
!               consistent with that which is expected in the rest of
!               gsi.
!
!
! !REVISION HISTORY:
!   06 Jul 2016 - J. G. de Mattos -  Initial code.
!
!
!EOP
!-------------------------------------------------------------------------
!BOC
      integer         :: iMax, jMax
      integer(i_kind) :: i, j

      character(len=64), parameter :: myname_=trim(myname)//' :: SFCTrans( )'


#ifdef DEBUG
    WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif
!  Compute mean along southern and northern latitudes

      iMax = size(GridIn,1) ! lon
      jMax = size(GridIn,2) ! lat

!  Transfer local work array to output grid

      do j = 1, iMax
         do i = 1, jMax
            GridOut(i,j) = GridIn( j, i )
         end do
      end do

   return


   end subroutine
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: Write_BRAMS
!              
! !DESCRIPTON: rotine to write atmospheric analysis fields from GSI to BRAMS 
!
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine writeBRAMS(increment, mype, mype_atm, mype_sfc )
!
! !USES:
!
      ! Model read 
      use BramsMod, only : bramsFile

      ! GSI kinds
      use kinds, only: r_kind, i_kind

      !GSI bundle
      use gsi_metguess_mod, only: gsi_metguess_bundle
      use gsi_bundlemod, only: gsi_bundlegetpointer
      
      !GSI Guess informations
      use guess_grids, only: ges_prsl, ges_prsi
      use guess_grids, only: ntguessig,ntguessfc,ifilesig,nfldsig

      use gsi_4dvar, only: lwrite4danl

      implicit none
!
! !INPUT PARAMETERS:
!
      integer(i_kind), intent(in   ) :: increment
      integer(i_kind), intent(in   ) :: mype
      integer(i_kind), intent(in   ) :: mype_atm
      integer(i_kind), intent(in   ) :: mype_sfc

! !REVISION HISTORY:
!
!   06 Jul 2016 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      character(len=64), parameter :: myname_=trim(myname)//' :: writeBRAMS( )'

      type(BRAMSFile) :: brams  ! BRAMS files data type

      !
      ! GSI bungle guess fields
      !

      real(r_kind), pointer, dimension(:,:  ) :: anl_z_it    => NULL()
      real(r_kind), pointer, dimension(:,:  ) :: anl_ps_it   => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_pi_it   => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_u_it    => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_v_it    => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_tv_it   => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_q_it    => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_oz_it   => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_cw_it   => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_ql_it   => NULL()
      real(r_kind), pointer, dimension(:,:,:) :: anl_qi_it   => NULL()
  
      !
      ! Auxiliary vars
      !

      integer :: istatus
      integer :: ntlevs
      integer :: it
      integer :: itout
      integer :: i

      character(len=20) :: FileAnl
      character(len=20) :: FileAnlFG

#ifdef DEBUG
      WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
      return
#endif

      !
      ! Write atmospheric analysis file
      ! Loop over all time analysis fields
      !

      if (.not.lwrite4danl) then
         ntlevs = 1
      else
         ntlevs = nfldsig
      end if
      
      do it = 1, ntlevs

         !
         ! define anl file name
         !

         if (increment>0) then
            FileAnl   = 'BRAMS.inc'
            FileAnlFG = 'BRAMS.inc.rg'
            itout     = increment
            if(mype.eq.0) write(6,*) 'WRITE_BRAMS: writing time slot ', itout
         else if (.not.lwrite4danl) then
            FileAnl   = 'BRAMS.anl'
            FileAnlFG = 'BRAMS.anl.rg'
            itout     = ntguessig
            if(mype.eq.0) write(6,*) 'WRITE_BRAMS: writing single analysis state for F ', itout
         else
            write(FileAnl,'("BRAMS.anl.",I2.2)') ifilesig(it)
            write(FileAnlFG,'("BRAMS.anl.rg.",I2.2)') ifilesig(it)
            itout = it
            if(mype.eq.0) write(6,*) 'WRITE_BRAMS: writing full analysis state for F ', itout
         endif


         !
         ! Allocating anl transfer fields
         !

         call AllocateAnlFields( )

         !
         ! Get fields from GSI Bundle
         !

         ! Orography
         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'z', anl_z_it, istatus)
         if(istatus.eq.0) g_z = anl_z_it

         ! Surface Pressure
         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'ps', anl_ps_it, istatus)
         if(istatus.eq.0) g_ps = anl_ps_it

         ! Pressure (Exner Function)
         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'pi', anl_pi_it, istatus)
         if(istatus.eq.0) g_pi = anl_pi_it

         ! Zonal wind component
         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'u', anl_u_it, istatus)
         if(istatus.eq.0) g_u = anl_u_it

         ! Meridional wind component
         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'v', anl_v_it, istatus)
         if(istatus.eq.0) g_v = anl_v_it
      
         ! Virtual Temperature
         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'tv', anl_tv_it, istatus)
         if(istatus.eq.0) g_tv = anl_tv_it

         ! Specific Humidity
         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'q', anl_q_it, istatus)
         if(istatus.eq.0) g_q = anl_q_it

!         ! Ozone
!         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'oz', anl_oz_it, istatus)
!         if(istatus.eq.0) g_oz = anl_oz_it
!
!         ! Total Cloud Water Content
!         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'cw', anl_cw_it, istatus)
!         if(istatus.eq.0) g_cw = anl_cw_it
!
!         ! Cloud Liquid Water Content
!
!         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'ql',anl_ql_it, istatus)
!         if(istatus.eq.0) g_ql = anl_ql_it
!
!         ! Cloud Ice Water Content
!
!         call gsi_bundlegetpointer (gsi_metguess_bundle(itout), 'qi', anl_qi_it, istatus)
!         if(istatus.eq.0) g_qi = anl_qi_it

         g_prsi = ges_prsi(:,:,:,itout)
         g_prsl = ges_prsl(:,:,:,itout)

         ! Compute geopotential Heights
         
         call comp_geop_hgt(g_z, g_q, g_tv, g_prsl, g_prsi, g_gh)

         ! Compute Relative Humidity
      
         call comp_rh(g_q, g_tv, g_prsl, g_prsi, g_rH)
         
         ! Compute Absolute Temperature (K)

         call comp_tk(g_q, g_tv, g_tk)

         !
         ! Write BRAMS ANL fields
         !
         if(mype.eq.0) write(6,*) 'write BRAMS ANL fields' 
         call GenWriteBRAMS (FileAnl, mype, mype_atm, itout )

         !
         ! Write BRAMS ANL fields (same fields as First Guess at regular grid)
         !

         if(mype.eq.0) write(6,*) 'write BRAMS ANL same as first guess fields at regular grid' 
         call GenWriteBRAMS_FG (FileAnlFG, mype, mype_atm, itout )

         !
         ! Clear anl transfer fields
         !

         call DeallocateAnlFields ( )

     enddo

  end subroutine
!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: GenWriteBRAMS
!
! !DESCRIPTON: rotine to read atmospheric guess fields
!
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine GenWriteBRAMS( FileAnl, MyPe, MyPe_Out, it )
!
! !USES:
!
      ! Model read 
      use BRAMSMod, only : bramsFile
      use ModConstants, only: i4, r4
      
      ! GSI grid informations
      use gridmod, only : grd=>grd_a

      ! GSI file guess information
      use guess_grids, only: ntguessig, ifilesig

      ! GSI analysis time
      use obsmod, only: iadate

      ! GSI 4dvar case
      use gsi_4dvar, only: ibdate, nhr_obsbin, lwrite4danl

      !GSI kinds
      use kinds, only : r_kind, i_kind

      !GSI constants
      use constants, only:zero 

      !MPI
      use mpimod, only : NPe

!
! !INPUT PARAMETERS:
!
      character(len=*), intent(in   ) :: FileAnl
      integer,          intent(in   ) :: MyPe     ! current mpi task id
      integer,          intent(in   ) :: MyPe_Out ! mpi task id to write
      integer,          intent(in   ) :: it

!
! !PARAMTERS:
!

      integer, parameter                  :: NVars = 5
!      character(len=40), dimension(NVars) :: VName = [                            &
!                                                      'UP                        ',& !  1 - u-wind (m/s)
!                                                      'VP                        ',& !  2 - v-wind (m/s) 
!                                                      'TEMPK                     ',& !  3 - temperature (k)
!                                                      'GEOPOTENTIAL              ',& !  4 - Geopotential H
!                                                      'RELATIVE HUMIDITY         ' & !  5 - Relative Humidity
!                                                     ] 

! !REVISION HISTORY:
!
!   26 Mar 2019 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC      
      character(len=100), parameter :: myname_=':: GenWriteBRAMS()'

      type(BRAMSFile) :: brams  ! BRAMS files data type

      !
      ! Variables to write GFCT spectral file
      !

      integer (i4)                          :: ifday, rc, nymd, nhms
      real    (r4)                          :: tod
      integer (r4), dimension(4)            :: idate, idatec

      real(r_kind), dimension(:),   allocatable :: work1 ! Physical Space
      real(r4),     dimension(:,:), allocatable :: work2 ! Physical Space
      real(r_kind), dimension(:,:), allocatable :: tmp   ! Physical Space
  
      !
      ! Auxiliary Variables
      !

      integer(i_kind),dimension(5) :: mydate
      integer(i_kind),dimension(8) :: ida,jda
      real(r_kind),   dimension(5) :: fha


      integer :: istat

      integer :: WVar(NPe)
      integer :: WLev(NPe)


      integer :: Nflds
      integer :: nLon
      integer :: nLat
      integer :: nlevs

      integer :: i, j, k, ij
      integer :: icount
      integer :: ivar
      integer :: iret
      integer :: ilev
      integer :: iPe

      integer :: iret_write

      character(len=20) :: FileGuess
      character(len=20) :: FileHeader

      integer :: OutUnit, iunit, ios
      logical :: isopen

      

#ifdef DEBUG
      WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      !
      ! Open File to write
      !

      if (myPe .eq. MyPe_Out)then

        ! Find a unit available to write anl
         find_unit:do iunit = 10, 254
   
            inquire (Unit = iunit, opened = isopen, iostat = ios )
            
            if (.not.isopen.and.ios.eq.0)then
               outUnit = iunit
               exit
            endif
   
            if(iunit .eq. 254)then
               write(*,'(A,1x,A)')trim(myname_),'Units from 10 to 254 are already in use!'
               stop
            endif
   
         end do find_unit

         !
         ! Open File to write
         !

         Open(Unit = OutUnit,       &
              File = trim(FileAnl),&
              Form = 'unformatted')

      endif

      !
      ! Set dir guess file name
      !

      write(FileGuess,'("BRAMS.vfm.",I2.2)') ifilesig(ntguessig)
      write(FileHeader,'("BRAMS.head.",I2.2)') ifilesig(ntguessig)

      !
      ! All tasks open and read header
      !

      Nflds = NVars * grd%nsig

      ! all tasks get ges header information

      iret = brams%BrOpen(FileHeader, FileGuess)

      if (iret.ne.0)then

         write(stdout,'(2(A,1x),I4,1x,A)')trim(myname_),'Problem to open/read BRAMS files at ',MyPe, 'rank!'

         stop
      endif


      !
      ! consistency check
      !

      nLon  = brams%GetDim('nnxp')
      nLat  = brams%GetDim('nnyp')
      nLevs = brams%GetDim('nnzp')

      if (                            &
           (nLon .ne. grd%nlon .or.   &
            nLat .ne. grd%nlat .or.   &
            nLevs.ne. grd%nsig) .and. &
            MyPe .eq. 0               &
         )then

         write (6,'(3(1x,A))') trim(myname_),':: *** ERROR *** reading ', trim(FileGuess)
         write (6,'(2(1x,A,1x,I4))')'<nlon>', nLon,'<->.',grd%nlon
         write (6,'(2(1x,A,1x,I4))')'<nlat>', nLat,'<->.',grd%nlat
         write (6,'(2(1x,A,1x,I4))')'<nlev>',nLevs,'<->.',grd%nsig

         stop
      endif


! Load date
      if (.not.lwrite4danl) then
        mydate = iadate
      else
!     increment mydate ...

         mydate = ibdate
         fha(:) = zero ; ida=0; jda=0
         fha(2) = real(nhr_obsbin*(it-1))  ! relative time interval in hours
         ida(1) = mydate(1) ! year
         ida(2) = mydate(2) ! month
         ida(3) = mydate(3) ! day
         ida(4) = 0         ! time zone
         ida(5) = mydate(4) ! hour

! Move date-time forward by nhr_assimilation hours
         call w3movdat(fha,ida,jda)
         mydate(1) = jda(1)
         mydate(2) = jda(2)
         mydate(3) = jda(3)
         mydate(4) = jda(5)
      end if

!      if (MyPe .eq. MyPe_Out)then
!         call BRAMS_WriteAnlHeader(BRAMS, mydate, iret)
!         iret_write = iret_write + iret
!      endif

      !
      !------------------------------------------!
      !
      wvar   = 0
      ilev   = 0
      icount = 1
      ivar   = 1
      do while (ivar .le. nvars)

         do ilev = 1, grd%nsig

            !----------------------------------------
            ! This is used to scatter between all Pe's
            ! Como cada Pe ira pegar um campo
            ! em um nivel, estas duas variaveis dirao
            ! qual variavel é lida em cada Pe
            !
            wvar(icount) = ivar ! What Variable
            wlev(icount) = ilev ! What Level
            !----------------------------------------

            if ( icount .eq. NPe .or. ( ivar .eq. NVars .and. ilev .eq. NLevs) )then

               !
               ! Transfer contents of 3-d subdomain array to 2-d array global
               !

               allocate(work1(nLon*nLat))

               call GenGatherBRAMS(&
                                    g_u,   & ! Input  Zonal wind component
                                    g_v,   & !        Meridional wind component
                                    g_tk,  & !        Absolute Temperature
                                    g_gh,  & !        Geopotential Height
                                    g_rH,  & !        Relative Humidity
                                    WVar,  & !        Position in Var list
                                    WLev,  & !        Level
                                    icount,& !        MaxCount until here
                                    MyPe,  & !
                                    work1  & ! OutPut Field read by Pe
                                  )

               if (MyPe .lt. iCount)then

                  allocate(work2(nLon,nLat))

                  call load_grid(work1, work2)

                  deallocate(work1)


                  if (MyPe .eq. MyPe_Out)then

                     call WriteField(OutUnit, work2, MyPe_Out, icount)

                  else   ! send to MyPe_out

                     call SendField(MyPe, MyPe_Out, work2)

                  endif

                  if(allocated(work2)) deallocate(work2)

               endif
               !
               ! reset all counters
               !
               
               WVar   = 0
               WLev   = 0
               icount = 0

            endif

            !
            ! Transform from Grid to Spec
            !

            icount = icount + 1
         enddo

         ivar = ivar + 1
      enddo

      !
      ! Close BRAMS Files
      !

      return
   end subroutine
!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: GenWriteBRAMS
!
! !DESCRIPTON: rotine to read atmospheric guess fields
!
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine GenWriteBRAMS_FG( FileAnl, MyPe, MyPe_Out, it )
!
! !USES:
!
      ! Model read 
      use BRAMSMod, only : bramsFile
      use ModConstants, only: i4, r4
      
      ! GSI grid informations
      use gridmod, only : grd=>grd_a

      ! GSI file guess information
      use guess_grids, only: ntguessig, ifilesig, ges_prsi, ges_prsl

      ! GSI analysis time
      use obsmod, only: iadate

      ! GSI 4dvar case
      use gsi_4dvar, only: ibdate, nhr_obsbin, lwrite4danl

      !GSI kinds
      use kinds, only : r_kind, i_kind

      !GSI constants
      use constants, only:zero 

      !MPI
      use mpimod, only : NPe

!
! !INPUT PARAMETERS:
!
      character(len=*), intent(in   ) :: FileAnl
      integer,          intent(in   ) :: MyPe     ! current mpi task id
      integer,          intent(in   ) :: MyPe_Out ! mpi task id to write
      integer,          intent(in   ) :: it

!
! !PARAMTERS:
!

      integer, parameter                  :: NVars = 8
!      character(len=40), dimension(NVars) :: VName = [                            &
!                                                      'UP                        ',& !  1 - u-wind (m/s)
!                                                      'VP                        ',& !  2 - v-wind (m/s) 
!                                                      'TEMPK                     ',& !  3 - temperature (k)
!                                                      'GEOPOTENTIAL              ',& !  4 - Geopotential H
!                                                      'RELATIVE HUMIDITY         ' & !  5 - Relative Humidity
!                                                     ] 

! !REVISION HISTORY:
!
!   26 Mar 2019 - J. G. de Mattos -  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC      
      character(len=100), parameter :: myname_=':: GenWriteBRAMS_FG()'

      type(BRAMSFile) :: brams  ! BRAMS files data type

      !
      ! Variables to write GFCT spectral file
      !

      real(r_kind), dimension(:),   allocatable :: work1
      real(r4),     dimension(:,:), allocatable :: work2
      real(r4),     dimension(:,:), allocatable :: work3
      real(r_kind), dimension(:,:), allocatable :: tmp
  
      !
      ! Auxiliary Variables
      !

      integer(i_kind),dimension(5) :: mydate
      integer(i_kind),dimension(8) :: ida,jda
      real(r_kind),   dimension(5) :: fha


      integer :: istat

      integer :: WVar(NPe)
      integer :: WLev(NPe)


      integer :: Nflds
      integer :: nnxp
      integer :: nnyp
      integer :: nLon
      integer :: nLat
      integer :: nlevs
      integer :: varNLevs

      integer :: i, j, k, ij
      integer :: x, y
      integer :: icount
      integer :: ivar
      integer :: iret
      integer :: ilev
      integer :: iPe

      integer :: iret_write

      character(len=20) :: FileGuess
      character(len=20) :: FileHeader

      integer :: OutUnit, iunit, ios
      logical :: isopen

      real(r4), allocatable :: xmap(:,:), ymap(:,:)
      real(r4), allocatable :: weight(:,:,:)
      real(r4) :: sumWeight
      

#ifdef DEBUG
      WRITE(stdout,'(     2A)')'Hello from ', trim(myname_)
#endif

      !
      ! Open File to write
      !

      if (myPe .eq. MyPe_Out)then

        ! Find a unit available to write anl
         find_unit:do iunit = 10, 254
   
            inquire (Unit = iunit, opened = isopen, iostat = ios )
            
            if (.not.isopen.and.ios.eq.0)then
               outUnit = iunit
               exit
            endif
   
            if(iunit .eq. 254)then
               write(*,'(A,1x,A)')trim(myname_),'Units from 10 to 254 are already in use!'
               stop
            endif
   
         end do find_unit

         !
         ! Open File to write
         !

         Open(Unit = OutUnit,       &
              File = trim(FileAnl),&
              Form = 'unformatted')

      endif

      !
      ! Set dir guess file name
      !

      write(FileGuess,'("BRAMS.vfm.",I2.2)') ifilesig(ntguessig)
      write(FileHeader,'("BRAMS.head.",I2.2)') ifilesig(ntguessig)

      !
      ! All tasks open and read header
      !

      Nflds = (NVars-1) * grd%nsig + 1

      ! all tasks get ges header information

      iret = brams%BrOpen(FileHeader, FileGuess)

      if (iret.ne.0)then

         write(stdout,'(2(A,1x),I4,1x,A)')trim(myname_),'Problem to open/read BRAMS files at ',MyPe, 'rank!'

         stop
      endif


      !
      ! consistency check
      !

      nnxp  = brams%GetDim('nnxp')
      nnyp  = brams%GetDim('nnyp')
      nLevs = brams%GetDim('nnzp')

      nLon  = brams%GetDim('nlon')
      nLat  = brams%GetDim('nlat')

      allocate(xMap(nLon,nLat))
      allocate(yMap(nLon,nLat))
      allocate(weight(nLon,nLat,4))
   
      call brams%GetMapCoord('xmap', xMap)
      call brams%GetMapCoord('ymap', yMap)
      call brams%GetInterpWeight(weight)

      if (                            &
           (nnxp .ne. grd%nlon .or.   &
            nnyp .ne. grd%nlat .or.   &
            nLevs.ne. grd%nsig) .and. &
            MyPe .eq. 0               &
         )then

         write (6,'(3(1x,A))') trim(myname_),':: *** ERROR *** reading ', trim(FileGuess)
         write (6,'(2(1x,A,1x,I4))')'<nlon>', nnxp,'<->.',grd%nlon
         write (6,'(2(1x,A,1x,I4))')'<nlat>', nnyp,'<->.',grd%nlat
         write (6,'(2(1x,A,1x,I4))')'<nlev>',nLevs,'<->.',grd%nsig

         stop
      endif


! Load date
      if (.not.lwrite4danl) then
        mydate = iadate
      else
!     increment mydate ...

         mydate = ibdate
         fha(:) = zero ; ida=0; jda=0
         fha(2) = real(nhr_obsbin*(it-1))  ! relative time interval in hours
         ida(1) = mydate(1) ! year
         ida(2) = mydate(2) ! month
         ida(3) = mydate(3) ! day
         ida(4) = 0         ! time zone
         ida(5) = mydate(4) ! hour

! Move date-time forward by nhr_assimilation hours
         call w3movdat(fha,ida,jda)
         mydate(1) = jda(1)
         mydate(2) = jda(2)
         mydate(3) = jda(3)
         mydate(4) = jda(5)
      end if

!      if (MyPe .eq. MyPe_Out)then
!         call BRAMS_WriteAnlHeader(BRAMS, mydate, iret)
!         iret_write = iret_write + iret
!      endif

      ! some adjustments
      ! from centibar to milibar (kPa to hPa)
      ! 
      g_ps = g_ps*10.0_r_kind
      g_pi = g_pi*10.0_r_kind
      g_prsi = g_prsi*10.0_r_kind
      g_prsl = g_prsl*10.0_r_kind

      !
      !------------------------------------------!
      !
      wvar   = 0
      ilev   = 0
      icount = 1
      ivar   = 1
      do while (ivar .le. nvars)
         if (ivar .eq. 1)then
            varNLevs = 1
         else
            varNLevs = grd%nsig
         endif
         do ilev = 1, varNLevs!grd%nsig

            !----------------------------------------
            ! This is used to scatter between all Pe's
            ! Como cada Pe ira pegar um campo
            ! em um nivel, estas duas variaveis dirao
            ! qual variavel é lida em cada Pe
            !
            wvar(icount) = ivar ! What Variable
            wlev(icount) = ilev ! What Level
            !----------------------------------------

            if ( icount .eq. NPe .or. ( ivar .eq. NVars .and. ilev .eq. NLevs) )then

               !
               ! Transfer contents of 3-d subdomain array to 2-d array global
               !

               allocate(work1(nnxp*nnyp))

               call GenGatherBRAMS_FG (       & ! Input:
                                       g_ps,  & !  * Surface Pressure
                                       g_tv,  & !  * Virtual temperature
                                       g_q,   & !  * Specific Humidity
                                       g_u,   & !  * Zonal wind component
                                       g_v,   & !  * Meridional wind component
                                       g_pi,  & !  * Pressure (Exner Function)
                                       g_prsi,& !  * Pressure at interface
                                       g_prsl,& !  * Pressure level
                                       WVar,  & !  * Position in Var list
                                       WLev,  & !  * Level
                                       icount,& !  * MaxCount until here
                                       MyPe,  & !
                                       work1  & ! OutPut Field read by Pe
                                     )

               if (MyPe .lt. iCount)then

                  allocate(work2(nnxp,nnyp))

                  call load_grid(work1, work2)

                  deallocate(work1)

                  ! interpolate to regular lat/lon grid
                  allocate(work3(nlon,nlat))

                  do j=1,nlat
                     do i=1,nlon
                        x = xMap(i,j)
                        y = yMap(i,j)

                        sumWeight = sum(weight(i,j,:))
                        if (sumWeight .eq. 0)then
                           work3(i,j) = undef
                        else
                           work3(i,j) = work2(x  ,y  ) * weight(i,j,1) + &
                                        work2(x+1,y  ) * weight(i,j,2) + &
                                        work2(x  ,y+1) * weight(i,j,3) + &
                                        work2(x+1,y+1) * weight(i,j,4)
                        endif
                     enddo
                  enddo

                  if(allocated(work2)) deallocate(work2)

                  if (MyPe .eq. MyPe_Out)then

                     call WriteField(OutUnit, work3, MyPe_Out, icount)

                  else   ! send to MyPe_out

                     call SendField(MyPe, MyPe_Out, work3)

                  endif

                  if(allocated(work3)) deallocate(work3)

               endif
               !
               ! reset all counters
               !
               
               WVar   = 0
               WLev   = 0
               icount = 0

            endif

            !
            ! Transform from Grid to Spec
            !

            icount = icount + 1
         enddo

         ivar = ivar + 1
      enddo

      !
      ! Close BRAMS Files
      !

      return
   end subroutine
!EOC
!-------------------------------------------------------------------------
!    NOAA/NCEP, National Centers for Environmental Prediction GSI        !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  load_grid --- strip off south/north latitude rows
!
! !INTERFACE:
!
 subroutine load_grid(grid_in,grid_out)

! !USES:
   use gridmod, only: grd => grd_a
   use kinds, only : r_kind, i_kind
   use ModConstants, only : r4
   implicit none

! !INPUT PARAMETERS:

   real(r_kind),dimension(max(grd%iglobal,grd%itotsub)), intent(in   ) :: grid_in  ! input grid
   real(r4),    dimension(grd%nlon,grd%nlat),            intent(  out) :: grid_out ! output grid

! !DESCRIPTION: This routine prepares grids for use in splib
!               grid to spectral tranforms.  This preparation
!               entails to two steps
!                  1) reorder indexing of the latitude direction.
!                     The GSI ordering is south to north.
!                     
!                  2) put on (nLon,nLat) intead (nLat,nLon) of GSI
!
! !REVISION HISTORY:
!   2004-08-27  treadon
!   2013-10-25  todling - move from gridmod to this module
!   2020-01-24  J.G.Z. de Mattos - adapt to brams model
!
! !REMARKS:
!   language: f90
!   machine:  ibm rs/6000
!
! !AUTHOR:
!   treadon          org: np23                date: 2004-08-27
!
!EOP
!-------------------------------------------------------------------------
   integer(i_kind) i,j,k,nlatm1,jj
   real(r_kind),dimension(grd%nlon,grd%nlat) :: grid
   real(r_kind),dimension(grd%nlat,grd%nlon) :: swap

!  Transfer input grid from 1d to 2d local array.  

   do k=1,grd%iglobal
      i = grd%ltosi(k) ! latitude
      j = grd%ltosj(k) ! longitude

      if(i .ge. 1 .and. i .le. grd%nLat .and. j .ge. 1 .and. j .le. grd%nLon)then
         swap(i,j) = grid_in(k)
      endif

   end do
   
   do j=1,grd%nlat
      do i=1,grd%nlon
         grid_out(i,j) = swap(j,i)
      enddo
   enddo
   
!  Transfer contents of local array to output array.
!   nlatm1=grd%nlat-1
!   do j=1,grd%nlat
!      jj=j-1
!      do i=1,grd%nlon
!         grid_out(i,jj)=real(grid(i,j),r4)
!      end do
!   end do
   
   return
 end subroutine load_grid
!
!EOC
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: GenGatherBRAMS
!              
! !DESCRIPTION: Transfer contents of 3d subdomains to 2d work arrays over pes
!
!             
!                 
!\\
!\\
! !INTERFACE:
!

   subroutine GenGatherBRAMS( g_u, g_v, g_tk, g_gh, g_rH, &
                            wvar, wlev, icount, MyPe, work  )
!
! !USES:
!

      use kinds, only: r_kind,i_kind
      use mpimod, only: npe,mpi_comm_world,ierror,mpi_rtype
      use gridmod, only: strip
      use gridmod, only: grd => grd_a

      implicit none
!
! !INPUT PARAMETERS:
!
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_u
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_v
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_tk
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_gh
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_rH

  integer(i_kind),dimension(npe), intent(in   ) :: wvar
  integer(i_kind),dimension(npe), intent(in   ) :: wlev
  integer(i_kind),                intent(in   ) :: icount
  integer(i_kind),                intent(in   ) :: MyPe


!
! !OUTPUT PARAMETERS:
!

  real(r_kind), dimension(:),     intent(inout) :: work

!
! !REVISION HISTORY:
!   2013-06-19  treadon
!   2013-10-24  todling   - update interface to strip
!   2016-08-04  de Mattos - Adpat to BRAMS model
!
!EOP
!-------------------------------------------------------------------------
!BOC

  real(r_kind),dimension(grd%lat1*grd%lon1,npe) :: sub

  integer :: ivar
  integer :: ilev
  integer :: Pe
  integer :: var



!$omp parallel do  schedule(dynamic,1) private(Pe,ivar,ilev)

      do Pe = 1, icount

         ivar = wvar(Pe)
         ilev = wlev(Pe)

         select case ( ivar )
            case (1) ! Zonal wind component

               call strip ( g_u(:,:,ilev), sub(:,Pe) )

            case (2) ! Meridional wind component

               call strip ( g_v(:,:,ilev), sub(:,Pe) )

            case (3) ! absotute temperature

               call strip ( g_tk(:,:,ilev), sub(:,Pe) )
           
            case (4) ! Geopotential Heigth

               call strip ( g_gh(:,:,ilev), sub(:,Pe) )

            case (5) ! Relative Humidity

               call strip ( g_rH(:,:,ilev), sub(:,Pe) )

         end select
      enddo

      call mpi_alltoallv(sub,           &
                         grd%isc_g,     &
                         grd%isd_g,     &
                         mpi_rtype,     &
                         work,          &
                         grd%ijn,       &
                         grd%displs_g,  &
                         mpi_rtype,     &
                         mpi_comm_world,&
                         ierror         &
                         )
   end subroutine
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOP
!
! !ROUTINE: GenGatherBRAMS
!              
! !DESCRIPTION: Transfer contents of 3d subdomains to 2d work arrays over pes
!
!             
!                 
!\\
!\\
! !INTERFACE:
!
!   subroutine GenGatherBRAMS_FG(g_ps, g_tv, g_q, g_u, g_v, g_pi,&
!                            wvar, wlev, icount, MyPe, work  )

   subroutine GenGatherBRAMS_FG(g_ps, g_tv, g_q, g_u, g_v, &
                                g_pi, g_prsi, g_prsl,&
                            wvar, wlev, icount, MyPe, work  )
!
! !USES:
!

      use kinds, only: r_kind,i_kind
      use mpimod, only: npe,mpi_comm_world,ierror,mpi_rtype
      use gridmod, only: strip
      use gridmod, only: grd => grd_a

      implicit none
!
! !INPUT PARAMETERS:
!

  real(r_kind),dimension(:,:  ),  intent(in   ) :: g_ps
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_tv
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_q
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_u
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_v
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_pi
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_prsi
  real(r_kind),dimension(:,:,:),  intent(in   ) :: g_prsl

  integer(i_kind),dimension(npe), intent(in   ) :: wvar
  integer(i_kind),dimension(npe), intent(in   ) :: wlev
  integer(i_kind),                intent(in   ) :: icount
  integer(i_kind),                intent(in   ) :: MyPe


!
! !OUTPUT PARAMETERS:
!

  real(r_kind), dimension(:),     intent(inout) :: work

!
! !REVISION HISTORY:
!   2013-06-19  treadon
!   2013-10-24  todling   - update interface to strip
!   2016-08-04  de Mattos - Adpat to BRAMS model
!
!EOP
!-------------------------------------------------------------------------
!BOC

  real(r_kind),dimension(grd%lat1*grd%lon1,npe) :: sub

  integer :: ivar
  integer :: ilev
  integer :: Pe
  integer :: var



!$omp parallel do  schedule(dynamic,1) private(Pe,ivar,ilev)

      do Pe = 1, icount

         ivar = wvar(Pe)
         ilev = wlev(Pe)

         select case ( ivar )
            case (1) ! Surface Pressure

               call strip ( g_ps(:,:), sub(:,Pe) )

            case (2) ! Virtual Temperature

               call strip ( g_tv(:,:,ilev), sub(:,Pe) )

            case (3) ! Specific Humidity

               call strip ( g_q(:,:,ilev), sub(:,Pe) )

            case (4) ! Zonal wind component

               call strip ( g_u(:,:,ilev), sub(:,Pe) )

            case (5) ! Meridional wind component

               call strip ( g_v(:,:,ilev), sub(:,Pe) )

            case (6) ! Pressure (Exner Function)

               call strip ( g_pi(:,:,ilev), sub(:,Pe) )

            case (7) ! Pressure at interface

               call strip ( g_prsi(:,:,ilev), sub(:,Pe) )

            case (8) ! Pressure at boundary

               call strip ( g_prsl(:,:,ilev), sub(:,Pe) )

         end select
      enddo

      call mpi_alltoallv(sub,           &
                         grd%isc_g,     &
                         grd%isd_g,     &
                         mpi_rtype,     &
                         work,          &
                         grd%ijn,       &
                         grd%displs_g,  &
                         mpi_rtype,     &
                         mpi_comm_world,&
                         ierror         &
                         )
   end subroutine
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: BAM_SendField_ -
!
! 
! !DESCRIPTION:
!
!
! !INTERFACE:
!   
  subroutine SendField(MyPe, toPe, Field, istat)
     use ModConstants, only: i4, r4
     use MPI

     implicit none
!
! !INPUT PARAMETERS:
!
     ! Source Pe
     integer(i4),           intent(in   ) :: MyPe

     ! Target Pe
     integer(i4),           intent(in   ) :: toPe

     ! Field to be send
     real(r4),              intent(in   ) :: field(:,:)
!
! !OUTPUT PARAMETERS:
! 
     integer(i4), optional, intent(  out) :: istat
!
! !REVISION HISTORY: 
!
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
     character(len=100), parameter :: myname_=':: SendField( ... )'

     integer :: iret
     integer :: sizefield

     if(present(istat)) istat = 0

     sizefield = size(field)
      
     call mpi_send(field, sizefield, MPI_FLOAT, toPe, 100, MPI_COMM_WORLD, iret)

     if(iret .ne. 0)then
        write(*,'(2A,I5,1x,A,1x,I5)')trim(myname_),': ERROR to send field from',MyPe,'to',toPe
        if(present(istat)) istat = iret
        return
     endif

  end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: WriteField_MPIr4 - write fields of BAM files from a MPI Pe.
!
! 
! !DESCRIPTION: Esta rotina escreve um campo do modelo BAM 
!               
!
! !INTERFACE:
! 

  subroutine WriteField(OutUnit, field, OutPe, iCount, istat)
     use ModConstants, only: i4, r4
     use MPI

     implicit none
!
! !INPUT PARAMETERS:
! 
     ! Output logical unit
     integer(i4),           intent(in   ) :: OutUnit

     ! Field to be writed
     real(r4),              intent(in   ) :: field(:,:)

     ! How many Pe's are working
     integer(i4),           intent(in   ) :: iCount

     ! What Pe will write
     integer(i4),           intent(in   ) :: OutPe
!
! !OUTPUT PARAMETERS:
! 

     integer(i4), optional, intent(  out) :: istat

!
! !REVISION HISTORY: 
!
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
     character(len=100), parameter :: myname_=':: WriteFields( ... )'

     real(r4), allocatable :: buff(:,:)
     integer :: sizebuff
     integer :: Pe
     integer :: iret
     integer :: isize, jsize
     integer :: status(MPI_Status_size)
   
     if(present(istat)) istat = 0

     isize    = size(field,1)
     jsize    = size(field,2)
     sizebuff = isize*jsize

     allocate(buff(isize,jsize))

      

     do Pe = 0, iCount-1

        if(Pe .eq. OutPe)then

           call WriteField_Serial(OutUnit, field, iret)

           if(iret .ne. 0)then
              write(stdout,*)trim(myname_),': error to write field, ',iret
              if(present(istat)) istat = iret
              stop
           endif

        else

           call mpi_recv(buff, sizebuff, MPI_FLOAT, Pe, 100, MPI_COMM_WORLD, status, iret)

           call WriteField_Serial(OutUnit, buff, iret)
           if(iret .ne. 0)then
              if(present(istat)) istat = iret
              stop
           endif
        endif

     end do
      
     deallocate(buff)

  end subroutine
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: WriteField_Serialr4 - 
!
! 
! !DESCRIPTION: 
!              
!
! !INTERFACE:
!   

  subroutine WriteField_Serial(OutUnit, Field, istat)
     use ModConstants, only: i4, r4

     implicit none
!
! !INPUT PARAMETERS:
!
     ! Output logical unit
     integer(i4),           intent(in   ) :: OutUnit

     ! Field to be writed
     real(r4),              intent(in   ) :: Field(:,:)
!
! !OUTPUT PARAMETERS:
! 

     integer(i4), optional, intent(  out) :: istat

!
! !REVISION HISTORY: 
!
!  11 Oct 2016 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
     character(len=100), parameter :: myname_=':: WriteFields_Serial( ... )'

     integer :: iret
     integer :: siz

     siz=size(Field)

     if(present(istat)) istat = 0

     write( OutUnit, iostat = iret ) Field

     if(iret .ne. 0)then
        write(*,'(1A)')trim(myname_),': ERROR to write BRAMS file'
        if(present(istat)) istat = iret
        return
     endif

     return

  end subroutine
!
!EOC
!-------------------------------------------------------------------------
!    NOAA/NCEP, National Centers for Environmental Prediction GSI        !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: load_geop_hgt --- Populate guess geopotential height
!
! !INTERFACE:
!
  subroutine comp_geop_hgt(z, q, tv, prsl, prsi, geop_hgt)
     use ModConstants, only: i4, r4
     use guess_grids, only: use_compress

! !USES:

    use constants, only: one,eps, rd, grav, half, t0c, fv
    use constants, only: cpf_a0, cpf_a1, cpf_a2, cpf_b0, cpf_b1, cpf_c0, cpf_c1, cpf_d, cpf_e
    use constants, only: psv_a, psv_b, psv_c, psv_d
    use constants, only: ef_alpha, ef_beta, ef_gamma
    use gridmod, only: lat2, lon2, nsig

    implicit none

! !INPUT PARAMETERS:
    real(r_kind),dimension(:,:,:) :: tv
    real(r_kind),dimension(:,:,:) :: q
    real(r_kind),dimension(:,:,:) :: prsl
    real(r_kind),dimension(:,:,:) :: prsi
    real(r_kind),dimension(:,:  ) :: z

! !OUTPUT PARAMETER:
    real(r_kind),dimension(:,:,:) ::geop_hgt

! !DESCRIPTION: Obtain geopotential height fields
!
! !REVISION HISTORY:
!   2020-09-28  de Mattos - Version adapted from load_geop_ght
!                           to be used by BRAMS anl
!
! !REMARKS:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origi11 n 2000; Compaq/HP
!
! !AUTHOR:
!   treadon          org: w/nmc20      date: 2003-10-15
!
!EOP
!-------------------------------------------------------------------------

    character(len=*),parameter::myname_=myname//'comp_geop_hgt'
    real(r_kind),parameter:: thousand = 1000.0_r_kind

    integer(i4) i,j,k,jj,ier,istatus
    real(r_kind) h,dz,rdog
    real(r_kind),dimension(nsig+1):: height
    real(r_kind) cmpr, x_v, rl_hm, fact, pw, tmp_K, tmp_C, prs_sv, prs_a, ehn_fct, prs_v


    rdog = rd/grav

    if (use_compress) then

!      Compute compressibility factor (Picard et al 2008) and geopotential heights at midpoint 
!      of each layer

       do j=1,lon2
          do i=1,lat2
             k  = 1
             fact    = one + fv * q(i,j,k)
             pw      = eps + q(i,j,k)*( one - eps )
             tmp_K   = tv(i,j,k) / fact
             tmp_C   = tmp_K - t0c
             prs_sv  = exp(psv_a*tmp_K**2 + psv_b*tmp_K + psv_c + psv_d/tmp_K)        ! Pvap sat, eq A1.1 (Pa)
             prs_a   = thousand * exp(half*(log(prsi(i,j,k)) + log(prsl(i,j,k))))     ! (Pa) 
             ehn_fct = ef_alpha + ef_beta*prs_a + ef_gamma*tmp_C**2 ! enhancement factor (eq. A1.2)
             prs_v   = q(i,j,k) * prs_a / pw   ! vapor pressure (Pa)
             rl_hm   = prs_v / prs_sv    ! relative humidity
             x_v     = rl_hm * ehn_fct * prs_sv / prs_a     ! molar fraction of water vapor (eq. A1.3)
 
             ! Compressibility factor (eq A1.4 from Picard et al 2008)
             cmpr = one - (prs_a/tmp_K) * (cpf_a0 + cpf_a1*tmp_C + cpf_a2*tmp_C**2 &
                        + (cpf_b0 + cpf_b1*tmp_C)*x_v + (cpf_c0 + cpf_c1*tmp_C)*x_v**2 ) &
                        + (prs_a**2/tmp_K**2) * (cpf_d + cpf_e*x_v**2)

             h  = rdog * tv(i,j,k)
             dz = h * cmpr * log(prsi(i,j,k)/prsl(i,j,k))
             height(k) = z(i,j) + dz   

             do k=2,nsig
                fact    = one + fv * half * (q(i,j,k-1)+q(i,j,k))
                pw      = eps + half * (q(i,j,k-1)+q(i,j,k)) * (one - eps)
                tmp_K   = half * (tv(i,j,k-1)+tv(i,j,k)) / fact
                tmp_C   = tmp_K - t0c
                prs_sv  = exp(psv_a*tmp_K**2 + psv_b*tmp_K + psv_c + psv_d/tmp_K)  ! eq A1.1 (Pa)
                prs_a   = thousand * exp(half*(log(prsl(i,j,k-1))+log(prsl(i,j,k))))   ! (Pa)
                ehn_fct = ef_alpha + ef_beta*prs_a + ef_gamma*tmp_C**2 ! enhancement factor (eq. A1.2)
                prs_v   = half*(q(i,j,k-1)+q(i,j,k) ) * prs_a / pw   ! (Pa)
                rl_hm   = prs_v / prs_sv    ! relative humidity
                x_v     = rl_hm * ehn_fct * prs_sv / prs_a     ! molar fraction of water vapor (eq. A1.3)
                cmpr    = one - (prs_a/tmp_K) * ( cpf_a0 + cpf_a1*tmp_C + cpf_a2*tmp_C**2 &
                          + (cpf_b0 + cpf_b1*tmp_C)*x_v + (cpf_c0 + cpf_c1*tmp_C)*x_v**2 ) &
                          + (prs_a**2/tmp_K**2) * (cpf_d + cpf_e*x_v**2)
                h       = rdog * half * (tv(i,j,k-1)+tv(i,j,k))
                dz      = h * cmpr * log(prsl(i,j,k-1)/prsl(i,j,k))

                height(k) = height(k-1) + dz
             end do

             do k=1,nsig
                geop_hgt(i,j,k)=height(k) - z(i,j)
             end do
          enddo
       enddo

    else

!      Compute geopotential height at midpoint of each layer
       do j=1,lon2
          do i=1,lat2
             k  = 1
             h  = rdog * tv(i,j,k)
             dz = h * log(prsi(i,j,k)/prsl(i,j,k))
             height(k) = z(i,j) + dz
 
             do k=2,nsig
                h  = rdog * half * (tv(i,j,k-1)+tv(i,j,k))
                dz = h * log(prsl(i,j,k-1)/prsl(i,j,k))
                height(k) = height(k-1) + dz
             end do

             do k=1,nsig
                geop_hgt(i,j,k)=height(k) - z(i,j)
             end do
          end do
       end do

    endif

    return
  end subroutine comp_geop_hgt
!
!EOC
!-------------------------------------------------------------------------
!    NOAA/NCEP, National Centers for Environmental Prediction GSI        !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: load_rh --- calculte relative Humidty
!
! !INTERFACE:
!
  subroutine comp_rH(q, tv, prsl, prsi, rH)

! !USES:
    use ModConstants, only: i4, r4
    use constants, only: half, one, eps, t0c, fv
    use constants, only: psv_a, psv_b, psv_c, psv_d
    use gridmod, only: lat2, lon2, nsig

    implicit none

! !INPUT PARAMETERS:
    real(r_kind),dimension(:,:,:) :: tv
    real(r_kind),dimension(:,:,:) :: q
    real(r_kind),dimension(:,:,:) :: prsl
    real(r_kind),dimension(:,:,:) :: prsi

! !OUTPUT PARAMETER:
    real(r_kind),dimension(:,:,:) :: rH

! !DESCRIPTION: Obtain geopotential height fields
!
! !REVISION HISTORY:
!   2020-09-28  de Mattos - Version adapted from load_geop_ght
!                           to be used by BRAMS anl
!
! !REMARKS:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origi11 n 2000; Compaq/HP
!
! !AUTHOR:
!   treadon          org: w/nmc20      date: 2003-10-15
!
!EOP
!-------------------------------------------------------------------------

    character(len=*),parameter::myname_=myname//'comp_rh'
    real(r_kind),parameter:: thousand = 1000.0_r_kind

    integer(i4) :: i,j,k
    real(r_kind)    :: fact, pw
    real(r_kind)    :: tmp_K, tmp_C
    real(r_kind)    :: prs_sv, prs_a, prs_v

    do j=1,lon2
       do i=1,lat2

          k         = 1
          fact      = one + fv * q(i,j,k)
          pw        = eps + q(i,j,k)*( one - eps )
          tmp_K     = tv(i,j,k) / fact
          tmp_C     = tmp_K - t0c
          prs_sv    = exp(psv_a*tmp_K**2 + psv_b*tmp_K + psv_c + psv_d/tmp_K)        ! Pvap sat, eq A1.1 (Pa)
          prs_a     = thousand * exp(half*(log(prsi(i,j,k)) + log(prsl(i,j,k))))     ! (Pa) 
          prs_v     = q(i,j,k) * prs_a / pw   ! vapor pressure (Pa)
          rH(i,j,k) = prs_v / prs_sv    ! relative humidity
 
          do k=2,nsig

             fact      = one + fv * half * (q(i,j,k-1)+q(i,j,k))
             pw        = eps + half * (q(i,j,k-1)+q(i,j,k)) * (one - eps)
             tmp_K     = half * (tv(i,j,k-1)+tv(i,j,k)) / fact
             tmp_C     = tmp_K - t0c
             prs_sv    = exp(psv_a*tmp_K**2 + psv_b*tmp_K + psv_c + psv_d/tmp_K)  ! eq A1.1 (Pa)
             prs_a     = thousand * exp(half*(log(prsl(i,j,k-1))+log(prsl(i,j,k))))   ! (Pa)
             prs_v     = half*(q(i,j,k-1)+q(i,j,k) ) * prs_a / pw   ! (Pa)
             rH(i,j,k) = prs_v / prs_sv    ! relative humidity

          end do

       enddo
    enddo

    return
  end subroutine comp_rH
  !
  !EOC
!-------------------------------------------------------------------------
!    NOAA/NCEP, National Centers for Environmental Prediction GSI        !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: load_rh --- calculte relative Humidty
!
! !INTERFACE:
!
  subroutine comp_tk(q, tv, tk)

! !USES:
    use ModConstants, only: i4, r4
    use constants, only: one, fv
    use gridmod, only: lat2, lon2, nsig

    implicit none

! !INPUT PARAMETERS:
    real(r_kind),dimension(:,:,:) :: tv
    real(r_kind),dimension(:,:,:) :: q

! !OUTPUT PARAMETER:
    real(r_kind),dimension(:,:,:) :: tk

! !DESCRIPTION: Obtain absolute temperature
!
! !REVISION HISTORY:
!   2020-09-28  de Mattos - Version adapted from load_geop_ght
!                           to be used by BRAMS anl
!
! !REMARKS:
!   language: f90
!   machine:  ibm rs/6000 sp; SGI Origi11 n 2000; Compaq/HP
!
! !AUTHOR:
!   treadon          org: w/nmc20      date: 2003-10-15
!
!EOP
!-------------------------------------------------------------------------

    character(len=*),parameter::myname_=myname//'comp_tk'

    integer(i4) :: i,j,k


    do k=1,nsig
       do j=1,lon2
          do i=1,lat2
             tk(i,j,k) = tv(i,j,k)/(one+fv*q(i,j,k))
          enddo
      enddo
    enddo

    return
  end subroutine comp_tk
  !
  !EOC
  !-------------------------------------------------------------------------


end module
