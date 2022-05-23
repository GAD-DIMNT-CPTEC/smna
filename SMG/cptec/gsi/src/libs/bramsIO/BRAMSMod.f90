!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!BOI
!
! !TITLE: Input/Output BRAMS interface Documentation \\ Version 1.0.0
!
! !AUTHORS: João Gerd Zell de Mattos
!
! !AFFILIATION: Modeling and Development Division, CPTEC/INPE
!
! !DATE: October 11, 2019
!
! !INTRODUCTION:
!      Input/Output BRAMS interface (BRAMMod) is a Fortran 90 collection of 
!      routines/functions for accessing the forecasting data files of the
!      Brazilian developments on the Regional Atmospheric Modelling System.
!
! \subsection{BRAMS data files}
!
! Os arquivos de previsão do modelo BRAMS são constituídos por dois arquivos em
! formatos distintos: 
!   \begin{enumerate}
!      \item Um arquivo no formato ASCII, denominado arquivo {\bf head} contém um 
!            cabeçalho descrevendo algumas informações da simulação, tais
!            como: data inicial, tempo de simulação, dimensão horizontal e vertical
!            e etc ... Este arquivo também possui uma tabela contendo as variáveis 
!            disponíveis no arquivo de previsão, identificando o tipo, dimensão e 
!            tamanho em bytes de cada variável;
!      \item Um arquivo no formato VFM ...

!   \end{enumerate}
!   
! \newpage
! \subsection{Principais Rotinas/Funções}
!
! \begin{verbatim}
!  ------------------------+-----------------------------------------
!     Routine/Function     |             Description
!  ------------------------+-----------------------------------------
!                          |
!   BrOpen                 | Open a BRAMS file
!   BrClose                | Close a BRAMS file
!   GetField               | Return a BRAMS field from a file
!   InqVar                 |
!   GetTimeInfo            |
!   GetDim                 |
!   GetMapCoord            |
!   GetVerticalInfo        !
!   GetVarNLevels          |
!                          |
!  ------------------------+-----------------------------------------
! \end{verbatim}
!
! \subsection{Exemplo de uso}
!
! O primeiro passo é carregar este módulo no programa fortran e definir uma estrutura 
! de dados contendo as informações do BRAMS.
! 
! \begin{enumerate}
!
!    \item Defina no início do programa fortran o uso do módulo BRAMSMod:
!
!    \begin{verbatim}
!       use BRAMSMod    
!    \end{verbatim}
!
!    \item Defina uma variável que conterá a estrutura de dados e uma variável
!          para receber o status de saída da função:
!
!    \begin{verbatim}
!       type(BramsFile) :: bramsArq
!       integer         :: ierr
!    \end{verbatim}
!
!    \item Abra o arquivo para leitura. Os arquivos do Brams são compostos por
!          um arquivo de cabecalho (head) e um arquivo ASCII (vfm), será necessário
!          passar somente o arquivo de cabecalho:
! 
!    \begin{verbatim}
!        ierr = bramsArq%Open('OPQUE-A-2019-06-05-000000-head.txt')
!    \end{verbatim}
! 
! 
!    \item Faça a leitura do campo disponível no BRAMS, defina antes o nome da 
!          variável, o nível vertical e aloque um vetor do tamanho necessário para
!          retornar o campo solicitado:
! 
!    \begin{verbatim}
!    
!       Allocate(grid(192*96))
!       
!       VName = 'UP'
!       ilev  = 1
!       
!       grid = bramsArq%GetField(trim(VName(ivar)), ilev, iret)
!   
!    \end{verbatim}
! 
!    \item Depois que forem lidos todos os campos necessários do modelo BRAMS feche o arquivo:
!    
!    \begin{verbatim}
!        ierr = bramsArq%Close( )
!    \end{verbatim}
! 
! \end{enumerate}
! 
!Veja os prólogos na próxima seção para obter detalhes adicionais.
!
!EOI
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!
Module BRAMSMod
  use m_inpak90, only : I90_LoadF, I90_fullRelease, I90_perr, I90_lcase
  use m_inpak90, only : I90_Label, I90_Gint, I90_GLine, I90_GFloat
  use ModConstants, only : deg2rad, earthRadius, earthDiameter
  use ModConstants, only : i4, i8, r4, r8, strlen, stdout, stdinp, stderr, undef
  use ModConstants, only : cp, cpor, p00
  use ModConstants, only : r60inv
  use m_time, only: cal2jul, jul2cal
  implicit none
  !  private

  integer, parameter :: VarNameMaxLength = 256
  integer, parameter :: iosize = 4 ! single precision file

  type BramsHeader
     character(len=16) :: VarName    ! field name
     integer           :: idim_type  ! field dimensionality (coded)
     integer           :: ngrid      ! grid number
     integer(kind=i8)  :: npointer   ! on ASCII coded files, field starting position
     integer(kind=i8)  :: nvalues    ! field size
     integer(kind=i8)  :: ipos       ! position inside binary file
     !----------------------------------------------------------------------------------------!
     ! this is for to use procedure pointer to
     ! get extra field derived from native field
     logical           :: extra
     procedure(FieldMethod), pointer, nopass :: func => null() ! function to obtain variable
     ! procedure( ),pointer, nopass :: f => null() ! subroutine to obtain variable
     !
     !----------------------------------------------------------------------------------------!

     type(BramsHeader), pointer :: next => null()
  end type BramsHeader

  type BramsFile
     private
     logical                   :: isOpen = .false.
     character(len=512)        :: Fhead  ! Header File Name
     character(len=512)        :: Fvfm   ! VFm file Name
     integer                   :: Funit
     integer                   :: nvars  ! total number of variables at vfm brams file
     integer                   :: nvars_extra ! total number of derived variables
     integer                   :: nnxp   ! Number of grid cells in the x-direction
     integer                   :: nnyp   ! Number of grid cells in the y-direction
     integer                   :: nnzp   ! Number of grid cells in the z-direction
     integer                   :: nzg    ! Number of soil layers
     integer                   :: nzs    ! Maximum number of snowpack layers
     integer                   :: npatch ! Number of subgrid patches
     integer                   :: nwave
     integer                   :: iyr    ! initial forecast year
     integer                   :: imo    ! initial forecast month
     integer                   :: idy    ! initial forecast day
     integer                   :: ihr    ! initial forecast hour
     integer                   :: imn    ! initial forecast minute
     integer                   :: isc    ! initial forecast second
     integer                   :: fyr    ! end forecast year
     integer                   :: fmo    ! end forecast month
     integer                   :: fdy    ! end forecast day
     integer                   :: fhr    ! end forecast hour
     integer                   :: fmn    ! end forecast minute
     integer                   :: fsc    ! end forecast second
     integer                   :: nLon        ! number of longitudes of regular lat/lon grid;
     integer                   :: nLat        ! number of latitudes of regular lat/lon grid;
     integer                   :: mcphys_type
     integer                   :: idimSize(7) ! array size of each idim_type
     integer                   :: jdim   ! all horizontal grids are 1D (jdim=0) or 2D (jdim=1);
     integer, pointer          :: xMap(:,:) => null( )  ! BRAMS x index, indexed by post grid indices
     integer, pointer          :: yMap(:,:) => null( )  ! BRAMS y index, indexed by post grid indices
     real,    pointer          :: x(:) => null( )       ! BRAMS x index, at regular lat/lon grid
     real,    pointer          :: y(:) => null( )       ! BRAMS y index, at regular lat/lon grid

     real                      :: fct_time
     real                      :: polelon
     real                      :: polelat
     real                      :: centlon     ! grid center longitude (degrees); from RAMSIN
     real                      :: centlat     ! grid center latitude (degrees); from RAMSIN
     real                      :: lonMin      ! start point of longitude (degree)
     real                      :: lonMax      ! end point of longitude (degree)
     real                      :: latMin      ! start point of latitude (degree)
     real                      :: latMax      ! end point of latitude (degree)
     real                      :: deltaxn     ! delta x on polar stereographic projection;
     real                      :: deltayn     ! delta y on polar stereographic projection;
     real                      :: delLon      ! delta longitude of regular lat/lon grid;
     real                      :: delLat      ! delta latitude of regular lat/lon grid;
     real                      :: ztop        ! altitude of top, km above sea level;
     real, pointer             :: xtn(:) => null( )     ! x coordinate of cell center on polar stereographic projection;
     real, pointer             :: ytn(:) => null( )     ! y coordinate of cell center on polar stereographic projection;
     real, pointer             :: xmn(:) => null( )     ! x coordinate of higher cell boundary on polar stereographic projection;
     real, pointer             :: ymn(:) => null( )     ! y coordinate of higher cell boundary on polar stereographic projection;
     real, pointer             :: dzt(:) => null( )     ! inverse of the atmosphere thickness (DZT = 1 / (z2- z1);
     real, pointer             :: glon(:,:) => null( )  ! longitudes on BRAMS rotated stereographic polar projection;
     real, pointer             :: glat(:,:) => null( )  ! latitudes on BRAMS rotated stereographic polar projection;
     real, pointer             :: fmapt(:,:) => null( ) ! Map scale Factor at t points
     real, pointer             :: dxt(:,:) => null( )     ! 
     real, pointer             :: dyt(:,:) => null( )     ! 

     real, pointer             :: weight(:,:,:) => null( )! weight of four neighbour BRAMS points;
     real, pointer             :: lat(:) => null( )     ! latitudes of regular lat/lon grid;
     real, pointer             :: lon(:) => null( )     ! longitudes of regular lat/lon grid.

     type(BramsHeader), pointer :: Header
     
     contains

        procedure, public  :: BrOpen => OpenFile_
        procedure, public  :: BrClose => CloseBHead_
        procedure, public  :: InqVar => InqVar_
        procedure, public  :: GetField => GetField_
        procedure, public  :: GetTimeInfo => GetTimeInfo_
        procedure, public  :: GetDim => GetDim_
        procedure, public  :: GetMapCoord => GetMapCoord_
        procedure, public  :: GetVerticalInfo => GetVerticalInfo_
        procedure, public  :: GetMapInfo => GetMapInfo_
        procedure, public  :: GetVarNLevels => GetNLevels_
        procedure, private :: OpenBramsHeader => OpenBHead_
        procedure, private :: RegisterField => RegisterField_
        procedure, public  :: printVars => printVars_
        procedure, public  :: GetInterpWeight => GetInterpWeight_
        
  end type BramsFile

  public :: bramsFile
  public :: xy_ll
  public :: ll_xy
  public :: distang
  public :: haversine

  character(len=*),parameter :: myname = 'BramsMod'

  !---------------------------------------------------------!
  ! - Procedure pointer a Fortran 2003 feature.
  !   This feature is used to create
  !   registerField subroutine
  !
!  type :: RealFunc
!     character(len=40) :: VarName   ! Name of variable
!     integer           :: idim_type ! field dimensionality (coded)
!     procedure(func),pointer, nopass :: f => null() ! function to obtain variable
!     type(RealFunc), pointer :: root => null()
!     type(RealFunc), pointer :: next => null()
!  end type ExtraHeader
!  interface
!     function func(x,y) result(z)
!        import BramsFile
!        type(BramsFile), intent(in) :: x
!        integer,         intent(in) :: y
!        real, pointer               :: z(:,:)
!     end function func
!  end interface
  interface
     subroutine FieldMethod(x,y,z)
        import bramsFile
        type(bramsFile), intent(in   ) :: x
        integer,         intent(in   ) :: y
        real,            intent(inout) :: z(:,:)
     end subroutine FieldMethod
  end interface
  !
  !---------------------------------------------------------!

  interface GetVal
     module procedure GetFloat_, GetInt_, &
                      GetArrayFloat_, GetArrayInt_
  end interface
  
  interface GetField
     module procedure GetField_
  end interface GetField

contains
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: OpenFile_ - routine to open Brams files 
!
! 
! !DESCRIPTION: Esta rotina serve de interface para o usuário acessar o arquivo
!               vfm do modelo BRAMS.
!               É necessário informar somente o nome do arquivo header que descreve
!               o conteúdo do arquivo vfm. É importante informar que o arquivo header
!               e o vfm devem possuir nomes compatíveis como no exemplo:
!              
!                  header: OPQUE-A-2019-06-05-000000-head.txt
!                  vfm   : OPQUE-A-2019-06-05-000000-g1.vfm
!
!               caso contrário é necessário informar o nome dos dois arquivos.
!
! !INTERFACE:
!

  function OpenFile_(self, HeaderFile, vfmFile) result(istatus)
  
    class(BramsFile),           intent(inout) :: self
!
! !INPUT PARAMETERS:
!
    ! Header file name
    character(len=*),           intent(in   ) :: HeaderFile
    ! vfm file name
    character(len=*), optional, intent(in   ) :: vfmFile
!
! !OUTPUT PARAMETERS:
!
    integer                                   :: istatus 
!    
! !SEE ALSO:
!   function/routines from inpak90 module
!
!    - I90_perr( )  - send a simple error message to _stderr_ 
!
!   function/routines from self module
!
!    - OpenBramsHeader( ) - read header auxiliary file of vfm file
!
! !REVISION HISTORY: 
!  11 Nov 2019 - J. G. de Mattos - Initial Version
!
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!

    character(len=*),parameter :: myname_=trim(myname)//' :: OpenFile_(...)'

    integer :: iret
    integer :: idx
    logical :: exists

    istatus = 0
    if(.not.present(vfmFile))then
       ! get vfm name file
       idx  = index(HeaderFile,'.',back=.true.) - 5
       self%Fvfm = HeaderFile(1:idx)//'g1.vfm'
    else
       self%Fvfm = vfmFile
    endif

    inquire(file=trim(self%Fvfm), exist=exists)
    if(.not.exists)then
       call i90_perr(myname_,'not found: '//trim(self%Fvfm),-1)
       stop
    endif
    
    ! Read Head file
    call self%OpenBramsHeader( HeaderFile, iret )
    if(iret.ne.0)then
       call i90_perr(myname_,'OpenBHead_("'//trim(HeaderFile)//'")',iret)
    else
       self%isOpen = .true.
    endif
    
    istatus = iret

  end function OpenFile_
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: OpenBHead_ - routine to read BRAMS run info contained in the header
!                         BRAMS file.
!
! 
! !DESCRIPTION: Esta rotina é reponsável pela leitura das informações contidas
!               no arquivo Head das simulações do BRAMS. As informações são
!               incluídas na estrutura bramsFile e posteriormente utilizadas pelas
!               demais rotinas deste modulo. Esta é uma rotina interna privada do
!               próprio módulo. 
!
! !INTERFACE:
!
  subroutine OpenBHead_(self, HeaderFile, istatus)
!
! !INPUT PARAMETERS:
!      
    ! BRAMS file structure
    class(BramsFile), intent(inout) :: self

    ! Header file name
    character(len=*), intent(in   ) :: HeaderFile
!
! !OUTPUT PARAMETERS:
!

    integer, optional               :: istatus

!    
! !SEE ALSO:
!   function/routines from inpak90 module
!
!    - I90_perr( )  - send a simple error message to _stderr_ 
!
!   function/routines from self module
!
!    - GetVal ( ) - 
!    - xy_ll ( ) -
!    - ll_xy ( ) -
!    - BuildAxis ( ) -
!    - RegisterField ( ) -
!
! !REVISION HISTORY: 
!  11 Nov 2019 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!

    character(len=*),parameter :: myname_=trim(myname)//' :: OpenBHead_(...)'

    type(bramsHeader), pointer :: Head => null()

    integer(kind=i4) :: i, j, iret
    integer(kind=i4) :: ymd, hms
    real(kind=r4)    :: incr
    integer(kind=i8) :: idx
    logical :: exists
    character(len=strlen) :: VarName

    real(kind=r4) :: deltaSum
    real(kind=r4) :: distXLow, distYLow
    real(kind=r4) :: xt2, yt2
    real(kind=r4) :: c1
!    real(kind=r4) :: x, y
    integer(kind=i4) :: xPost, yPost
    integer(kind=i4) :: xBrams, yBrams

    if(present(istatus)) istatus = 0


    self%Fhead = trim(HeaderFile)
    inquire(file=trim(HeaderFile), exist=exists)
    if(.not.exists)then
       call i90_perr(myname_,'not found: '//trim(HeaderFile),-1)
       stop
    endif

    call I90_LoadF(trim(self%Fhead), iret)
    if (iret .ne. 0 )then
       call i90_perr(myname_,'i90_LoadF("'//trim(self%Fhead)//'")',iret)
       stop
    endif

    ! BRAMS options
    call GetVal('__mcphys_type', self%mcphys_type)

    ! BRAMS grid
    call GetVal(  '__nnxp',   self%nnxp) ! Number of grid cells in the x-direction
    call GetVal(  '__nnyp',   self%nnyp) ! Number of grid cells in the y-direction
    call GetVal(  '__nnzp',   self%nnzp) ! Number of grid cells in the z-direction
    call GetVal(   '__nzg',    self%nzg) ! Number of soil layers
    call GetVal(   '__nzs',    self%nzs) ! Maximum number of snowpack layers
    call GetVal('__npatch', self%npatch) ! Number of subgrid patches
!    call GetVal( '__nwave',  self%nwave)

    !BRAMS time
    call GetVal( '__iyear1',      self%iyr)
    call GetVal('__imonth1',      self%imo)
    call GetVal( '__idate1',      self%idy)
    call GetVal( '__itime1',           hms)
    call GetVal(   '__time', self%fct_time)

    !
    ! BRAMS itime is like 0000 0600 1200 1800
    !
    self%ihr = INT ( hms / 100 )
    self%imn = MOD ( hms,  100 )
    self%isc = 0

    ! calcular aqui as datas da previsão utilizando o fct_time

    ymd  = self%iyr*10000 + self%imo*100 + self%idy + hms
    incr = self%fct_time * (r60inv*r60inv) / 24.0

    call jul2cal(cal2jul(ymd,hms) + incr, &
                 self%fyr, & ! Forecast year
                 self%fmo, & ! Forecast month
                 self%fdy, & ! Forecast day
                 self%fhr, & ! Forecast hour
                 self%fmn, & ! Forecast minute
                 self%fsc  & ! Forecast second
                 )

    !
    ! geographic info
    !

    call GetVal('__polelat', self%polelat)
    call GetVal('__polelon', self%polelon)

    call GetVal('__centlat', self%centlat)
    call GetVal('__centlon', self%centlon)

    call GetVal('__deltaxn', self%deltaxn)
    call GetVal('__deltayn', self%deltayn)

    call GetVal('__ztop', self%ztop)

    call GetVal('__jdim', self%jdim)

    allocate(self%xtn(self%nnxp))
    call GetVal(  '__xtn',   self%xtn)

    allocate(self%ytn(self%nnyp))
    call GetVal(  '__ytn',   self%ytn)

    allocate(self%xmn(self%nnxp))
    call GetVal(  '__xmn',   self%xmn)

    allocate(self%ymn(self%nnyp))
    call GetVal(  '__ymn',   self%ymn)

    allocate(self%dzt(self%nnzp))
    call GetVal(  '__dztn01',   self%dzt)

    call i90_FullRelease( iret )
    if (iret .ne. 0 )then
       call i90_perr(myname_,'i90_FullRelease( )',iret)
       stop
    endif    

    allocate(self%glon(self%nnxp,self%nnyp))
    allocate(self%glat(self%nnxp,self%nnyp))
    allocate(self%fmapt(self%nnxp,self%nnyp))

    !  Calculates map factor and geographical lat/lon at t-points for a given polar stereographic grid
    c1 = (2.0 * EarthRadius) ** 2.0
    do j=1,self%nnyp
       do i=1,self%nnxp

          xt2 = self%xtn(i)*self%xtn(i)
          yt2 = self%ytn(j)*self%ytn(j)

          self%fmapt(i,j) = 1.0 + (xt2+yt2) / c1

          call xy_ll(self%glat(i,j),self%glon(i,j),self%polelat,self%polelon,self%xtn(i),self%ytn(j))

       enddo
    enddo
 
    allocate(self%dxt(self%nnxp,self%nnyp))
    allocate(self%dyt(self%nnxp,self%nnyp))

     do j = 1, self%nnyp
        do i = 2,self%nnxp
           self%dxt(i,j)=self%fmapt(i,j)/(self%xmn(i)-self%xmn(i-1))
        enddo
        self%dxt(1,j)=self%dxt(2,j)*self%fmapt(1,j)/self%fmapt(2,j)
     enddo

     if (self%jdim == 1) then
        do i = 1,self%nnxp
           do j = 2,self%nnyp
              self%dyt(i,j)=self%fmapt(i,j)/(self%ymn(j)-self%ymn(j-1))
           enddo
           self%dyt(i,1)=self%dyt(i,2)*self%fmapt(i,1)/self%fmapt(i,2)
        enddo
     else
        do j=1,self%nnyp
           do i=1,self%nnxp
              self%dyt(i,j)=1./self%deltayn
           enddo
        enddo
     endif

     self%dxt = 1.0/self%dxt
     self%dyt = 1.0/self%dyt
    


    ! compute average delta longitude over all longitudes:
    ! sum of average delta longitude over all latitudes
    ! divided by number of latitudes
    deltaSum = sum( &
         (self%glon(self%nnxp,:) - self%glon(1,:))/&
         (self%nnxp-1)&
         )
    self%delLon = deltaSum / self%nnyp

    ! same procedure for latitudes:
    ! sum of average delta latitude over all longitudes
    ! divided by number of longitudes
    deltaSum = sum( &
         (self%glat(:,self%nnyp) - self%glat(:,1))/&
         (self%nnyp-1)&
         )
    self%delLat = deltaSum / self%nnxp

    !   find BRAMS grid extremes

    self%lonMin = minval(self%glon)
    self%lonMax = maxval(self%glon)
    self%latMin = minval(self%glat)
    self%latMax = maxval(self%glat)

    ! post grid number of points (#intervals + 1)

    self%nLon = 1 + ceiling((self%lonMax-self%lonMin)/self%delLon)
    self%nLat = 1 + ceiling((self%latMax-self%latMin)/self%delLat)

    allocate(self%Lon(self%nLon))
    call BuildAxis(self%nLon,self%lonMin, self%delLon, self%Lon)
    allocate(self%Lat(self%nLat))
    call BuildAxis(self%nLat,self%latMin, self%delLat, self%Lat)

    allocate(self%x(self%nLon))
    allocate(self%y(self%nLat))
    allocate(self%xMap(self%nLon, self%nLat))
    allocate(self%yMap(self%nLon, self%nLat))
    allocate(self%weight(self%nLon, self%nLat,4))
    self%xMap   = -1
    self%yMap   = -1
    self%weight = 0.0
    do yPost = 1, self%nLat
       do xPost = 1,self%nLon

          ! project lat-lon post grid point onto BRAMS grid

          call ll_xy(self%lat(yPost), self%lon(xPost), &
               self%polelat, self%polelon, &
               self%x(xPost), self%y(yPost))

          xBrams = floor((self%x(xPost)-self%xtn(1))/self%deltaxn) + 1
          yBrams = floor((self%y(yPost)-self%ytn(1))/self%deltayn) + 1
          if ( 1 <= xBrams .and. xBrams < self%nnxp .and. &
               1 <= yBrams .and. yBrams < self%nnyp ) then
             self%xMap(xPost,yPost) = xBrams
             self%yMap(xPost,yPost) = yBrams

             distXLow = (self%x(xPost) - self%xtn(xBrams))/self%deltaxn
             distYLow = (self%y(yPost) - self%ytn(yBrams))/self%deltayn
             self%weight(xPost,yPost,1) = (1.0-distXLow)*(1.0-distYLow)  !(i  , j  )
             self%weight(xPost,yPost,2) = (    distXLow)*(1.0-distYLow)  !(i+1, j  )
             self%weight(xPost,yPost,3) = (1.0-distXLow)*(    distYLow)  !(i  , j+1)
             self%weight(xPost,yPost,4) = (    distXLow)*(    distYLow)  !(i+1, j+1)

             ! minimize rounding errors: 
             ! the sum of all weights should be 1.0, but it is usualy not,
             ! due to rounding. Try to minimize rounding by replacing one
             ! of the weights by 1-sum(remaining weights)

             self%weight(xPost,yPost,4) = 1.0 - &
                  sum(self%weight(xPost,yPost,1:3))

          end if
       end do
    end do


    ! define idimSize for each idim_type
    self%idimSize(:) = 0
    self%idimSize(2) = self%nnxp*self%nnyp
    self%idimSize(3) = self%nnxp*self%nnyp*self%nnzp
    self%idimSize(4) = self%nnxp*self%nnyp*self%npatch*self%nzg
    self%idimSize(5) = self%nnxp*self%nnyp*self%npatch*self%nzs
    self%idimSize(6) = self%nnxp*self%nnyp*self%npatch
    self%idimSize(7) = self%nnxp*self%nnyp*self%nwave

    !
    ! Open header to read variable name and positions
    !

    Open(unit=98, file=trim(self%Fhead), form='formatted')
    read(98,*)self%nvars

    allocate(Head)
    self%Header => Head

    idx = 1
    do i=1,self%nvars

       read(98,*) VarName, Head%npointer, Head%idim_type, Head%ngrid, Head%nvalues
       Head%VarName = trim(adjustl(I90_lcase(VarName)))
       Head%extra   = .false.

       Head%ipos = idx + iosize
       idx       = Head%ipos + Head%nvalues*iosize + iosize

       if (i.lt.self%nvars)then
          allocate(Head%next)
          Head => Head%next
       endif

    enddo
    close(98)


    !registra os métodos para gerar variaveis dependentes
    self%nvars_extra = 0 ! inicializa com zero
    call self%registerField(      'VTMP', 3,  virtualTemperature)
    call self%registerField(      'UMES', 3,    SpecificHumidity)
    call self%registerField(      'Q2MT', 2, SpecificHumidity_2m)
    call self%registerField(     'CLOUD', 3,       CloudMixRatio)
    call self%registerField(   'GRAUPEL', 3,     GraupelMixRatio)
    call self%registerField(      'SNOW', 3,        SnowMixRatio)
    call self%registerField(      'RAIN', 3,        RainMixRatio)
    call self%registerField(       'ICE', 3,         IceMixRatio)
    call self%registerField( 'SFC_PRESS', 2,        CompSfcPress)


  end subroutine OpenBHead_
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: CloseBHead_ - routine to release bramsFile structure.
!
! 
! !DESCRIPTION: Esta rotina é reponsável pela leitura das informações contidas
 
!
! !INTERFACE:
!

  function CloseBHead_(self) result(istatus)
!
! !INPUT PARAMETERS:
!      
    ! BRAMS file structure
    class(BramsFile), intent(inout) :: self
!
! !OUTPUT PARAMETERS:
!
    integer                         :: istatus
! !REVISION HISTORY: 
!  11 Nov 2019 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!

    character(len=*),parameter :: myname_=trim(myname)//' ::CloseBHead(...)'

    integer :: ierr

    type(bramsHeader), pointer :: Head => null()

    istatus = 0

    self%nvars  = -1

    self%nnxp   = -1
    self%nnyp   = -1
    self%nnzp   = -1
    self%nzg    = -1
    self%nzs    = -1
    self%npatch = -1
    self%nwave  = -1

    self%iyr    = -1
    self%imo    = -1
    self%idy    = -1
    self%ihr    = -1
    self%imn    = -1
    self%isc    = -1

    self%fyr    = -1
    self%fmo    = -1
    self%fdy    = -1
    self%fhr    = -1
    self%fmn    = -1
    self%fsc    = -1

    self%polelat= -1
    self%polelon= -1
    self%deltaxn= -1
    self%deltayn= -1

    self%isOpen = .false.

    deallocate(self%xtn, stat=ierr)
    if (ierr .ne. 0 ) istatus = istatus + ierr
    deallocate(self%ytn, stat=ierr)
    if (ierr .ne. 0 ) istatus = istatus + ierr
    deallocate(self%xmn, stat=ierr)
    if (ierr .ne. 0 ) istatus = istatus + ierr
    deallocate(self%ymn, stat=ierr)
    if (ierr .ne. 0 ) istatus = istatus + ierr
    deallocate(self%dzt, stat=ierr)
    if (ierr .ne. 0 ) istatus = istatus + ierr
    deallocate(self%glon, stat=ierr)
    if (ierr .ne. 0 ) istatus = istatus + ierr
    deallocate(self%glat, stat=ierr)
    if (ierr .ne. 0 ) istatus = istatus + ierr
    deallocate(self%Lon, stat=ierr)
    if (ierr .ne. 0 ) istatus = istatus + ierr
    deallocate(self%Lat, stat=ierr)
    if (ierr .ne. 0 ) istatus = istatus + ierr
    deallocate(self%x, stat=ierr)
    if (ierr .ne. 0 ) istatus = istatus + ierr
    deallocate(self%y, stat=ierr)
    if (ierr .ne. 0 ) istatus = istatus + ierr
    deallocate(self%xMap, stat=ierr)
    if (ierr .ne. 0 ) istatus = istatus + ierr
    deallocate(self%yMap, stat=ierr)
    if (ierr .ne. 0 ) istatus = istatus + ierr
    deallocate(self%weight, stat=ierr)
    if (ierr .ne. 0 ) istatus = istatus + ierr
    deallocate(self%fmapt, stat=ierr)
    if (ierr .ne. 0 ) istatus = istatus + ierr
    deallocate(self%dxt, stat=ierr)
    if (ierr .ne. 0 ) istatus = istatus + ierr
    deallocate(self%dyt, stat=ierr)
    if (ierr .ne. 0 ) istatus = istatus + ierr

    Head => self%Header%next
    do
       deallocate(self%Header, stat=ierr)
       if (ierr .ne. 0 ) istatus = istatus + ierr
       if(.not.associated(Head))exit
       self%Header => Head
       Head => Head%next
    enddo

     if(istatus.ne.0)then
       call i90_perr(myname_,'some problem to release bramsFile ...',istatus)
    endif

  end function CloseBHead_
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP

  function InqVar_(self, vname) result(VarInfo)
    class(BramsFile), intent(in   ) :: self
    character(len=*), intent(in   ) :: vname

    type(bramsHeader), pointer :: VarInfo

    integer :: i
    character(len=strlen) :: VarName

    if (.not.self%isOpen) then
       call i90_perr(myname,':: BRAMS fille not opened !',-99)
    endif

    VarName = trim(adjustl(I90_lcase(vname)))
    
    VarInfo => self%Header

    do while(associated(VarInfo))
       if(trim(VarName).eq.trim(VarInfo%VarName))return
       VarInfo=>VarInfo%next
    enddo
    
    return

  end function InqVar_


  recursive subroutine GetField_(self, fld, vname, level, patch, nwave)
    class(BramsFile),  intent(in   ) :: self
    real,              intent(inout) :: fld(:,:)
    character(len=*),  intent(in   ) :: vname
    integer, optional, intent(in   ) :: level
    integer, optional, intent(in   ) :: patch
    integer, optional, intent(in   ) :: nwave

    character(len=*),parameter :: myname_= myname//' :: GetField_( ... )'

    
    character(len=256) :: msg

    real, pointer     :: a( : ) => null()
    real, allocatable :: b( : )
    real, allocatable :: c(:,:)

    real, pointer :: tmp(:,:,:,:) => null()

    type(bramsHeader), pointer :: Find => null()
    integer :: i, j
    integer :: x, y
    integer :: idx, pt
    integer :: iret
    integer :: nx, ny, npts

    real :: sumWeight
    character(len=strlen) :: VarName
    character(len=4) :: i_char, j_char

    if (.not.self%isOpen) then
       call i90_perr(myname,':: BRAMS file not opened yet!',-99)
       return
    endif

    VarName = trim(adjustl(I90_lcase(vname)))

!---------------------------------------------------------------------!
! Some sanity check

    Find => self%Inqvar(trim(VarName))
    if(.not.associated(Find))then
       call I90_perr(myname_,'Error: Field not found:'//trim(VName))
       stop
    endif

    if(present(level))then
       select case (Find%idim_type)
          case (3)
             if (level < 1 .or. level > self%nnzp)then
                write(i_char,'(I4)')level
                write(j_char,'(I4)')self%nnzp
                write(msg,'(3(1x,A),A)')'level',trim(adjustl(i_char)), &
                       'is out of range, bounds are 1:',trim(adjustl(j_char))
                call i90_perr(myname_,trim(msg))
                stop
             endif

          case (4)
             if (level < 1 .or. level > self%nzg)then
                write(i_char,'(I4)')level
                write(j_char,'(I4)')self%nzg
                write(msg,'(3(1x,A),A)')'level',trim(adjustl(i_char)), &
                       'is out of range, bounds are 1:',trim(adjustl(j_char))
                call i90_perr(myname_,trim(msg))
                stop
             endif

          case (5)
             if (level < 1 .or. level > self%nzs)then
                write(i_char,'(I4)')level
                write(j_char,'(I4)')self%nzs
                write(msg,'(3(1x,A),A)')'level',trim(adjustl(i_char)), &
                       'is out of range, bounds are 1:',trim(adjustl(j_char))
                call i90_perr(myname_,trim(msg))
                stop
             endif
       end select
    endif
    if(present(patch))then
       if (patch < 1 .or. patch > self%npatch)then
                write(i_char,'(I4)')patch
                write(j_char,'(I4)')self%npatch
                write(msg,'(3(1x,A),A)')'patch',trim(adjustl(i_char)), &
                       'is out of range, bounds are 1:',trim(adjustl(j_char))
          call i90_perr(myname_,trim(msg))
          stop
       endif
    endif

    if(present(nwave))then
       if (nwave < 1 .or. nwave > self%nwave)then
                write(i_char,'(I4)')nwave
                write(j_char,'(I4)')self%nwave
                write(msg,'(A,3(1x,A))')'nwave',trim(adjustl(i_char)), &
                       'is out of range, bounds are 1:',trim(adjustl(j_char))
          call i90_perr(myname_,trim(msg))
          
          stop
       endif
    endif
! end sanity check
!---------------------------------------------------------------------!
    if(Find%extra)then
       call Find%func(self, level, fld)
       return
    endif

    allocate(a(Find%nValues))
    allocate(b(Find%nValues))
    call RAMS_c_open(trim(self%Fvfm)//char(0),'r'//char(0))
    call vfirecr(10,a,Find%nvalues,'LIN', b, Find%npointer)

    deallocate(b)
    call RAMS_c_close()

    !
    ! A organização da matriz dentro do arquivo vfm é
    ! diferente da forma utilizada no arquivo history
    allocate(c(self%nnxp,self%nnyp))

    select case (Find%idim_type)

       case(2)
          ! fld(1,BRAMS%nnxp,BRAMS%nnyp,1) ! history and vfm
          !         fld => b
          tmp(1:1,1:self%nnxp,1:self%nnyp,1:1) => a
          do j=1,self%nnyp
             do i=1,self%nnxp
                c(i,j) = tmp(1,i,j,1)
             enddo
          enddo
   
       case(3)

          !fld(BRAMS%nnzp,BRAMS%nnxp,BRAMS%nnyp,1)) ! history
          !fld(BRAMS%nnxp,BRAMS%nnyp,BRAMS%nnzp,1)) ! vfm
          if(.not.present(level))then
             call i90_perr(myname_,'level no set!',-1)
             c = undef
             return
          endif

          tmp(1:self%nnxp,1:self%nnyp,1:self%nnzp,1:1) => a
          do j=1,self%nnyp
             do i=1,self%nnxp
                c(i,j) = tmp(i,j,level,1)
             enddo
          enddo

       case(4)
          ! fld(BRAMS%nzg,BRAMS%nnxp,BRAMS%nnyp,BRAMS%npatch) hist
          ! fld(BRAMS%nnxp,BRAMS%nnyp,BRAMS%nzg,BRAMS%npatch) vfm
!          tmp(1:self%nzg,1:self%nnxp,1:self%nnyp,1:self%npatch) => a
          tmp(1:self%nnxp,1:self%nnyp,1:self%nzg, 1:self%npatch) => a


          if(present(level).and.present(patch))then

             do j=1,self%nnyp
                do i=1,self%nnxp
                   c(i,j) = tmp(i,j,level,patch)
                enddo
             enddo

          else if(present(level) .and. .not. present(patch))then
             ! Resultado depende da variável:
             !  * incluir um flag especificando método a ser utilizado.
             !    Neste momento estou colocando o somatório do valor, mas
             !    pode-se definir o valor baseado na classe dominante,
             !    por exemplo, tipo de vegetação, tipo de solo ....
             !    Outro ponto é que o valor pode ser a média dos valores
             !    de cada patch
             !    
             do j=1,self%nnyp
                do i=1,self%nnxp
                   c(i,j) = sum(tmp(i, j, level, :))
                enddo
             enddo

          else if(present(patch) .and. .not. present(level))then
             ! Resultado depende da variável:
             !  * incluir um flag especificando método a ser utilizado.
             !    Neste momento estou colocando o valor médio, mas
             !    pode-se definir que o valor pode ser o somatório
             !    dos valores de cada nivel, como no caso dos patchs acima.
             !    
             do j=1,self%nnyp
                do i=1,self%nnxp
                   c(i,j) = sum(tmp(i,j,:,patch))/self%nzg
                enddo
             enddo

          else
             write(msg,'(A)') 'neither `level` nor `patch` are present!'
             call i90_perr(myname_,msg,-99)
             c = undef
             return
          endif
   
       case(5)
          ! fld(BRAMS%nzs,BRAMS%nnxp,BRAMS%nnyp,BRAMS%npatch)
          tmp(1:self%nzs,1:self%nnxp,1:self%nnyp,1:self%npatch) => a
          if(present(level).or.present(patch))then

             if(present(level).and.present(patch))then

                do j=1,self%nnyp
                   do i=1,self%nnxp
                      c(i,j) = tmp(level,i,j,patch)
                   enddo
                enddo

             else if(present(level) .and. .not. present(patch))then
                ! Resultado depende da variável:
                !  * incluir um flag especificando método a ser utilizado
                !    neste momento estou colocando o somatório do valor, mas
                !    pode-se definir o valor baseado na classe dominante,
                !    por exemplo, tipo de vegetação, tipo de solo ....
                !    Outro ponto é que o valor pode ser a média dos valores
                !    de cada patch
                !    
                do j=1,self%nnyp
                   do i=1,self%nnxp
                      c(i,j) = sum(tmp(level,i,j,:))
                   enddo
                enddo

             else if(present(patch) .and. .not. present(level))then
                ! Resultado depende da variável:
                !  * incluir um flag especificando método a ser utilizado
                !    neste momento estou colocando o valor médio, mas
                !    pode-se definir que o valor pode ser o somatório
                !    dos valores de cada nivel, como no caso dos patchs acima.
                !    
                do j=1,self%nnyp
                   do i=1,self%nnxp
                      c(i,j) = sum(tmp(:,i,j,patch))/self%nzs
                   enddo
                enddo

             endif
          else
             write(msg,'(A)') 'neither `level` nor `patch` are present!'
             call i90_perr(myname_,msg,-99)
             c = undef
             return
          endif
   
       case(6)
          ! fld(1,BRAMS%nnxp,BRAMS%nnyp,BRAMS%npatch)
          if(.not.present(patch))then
             call i90_perr(myname_,'patch no set!',-1)
             c = undef
             return
          endif

          tmp(1:1,1:self%nnxp,1:self%nnyp,1:self%npatch) => a

          do j=1,self%nnyp
             do i=1,self%nnxp
                c(i,j) = tmp(1, i, j, patch)
             enddo
          enddo

       case(7)
          ! fld(1,BRAMS%nnxp,BRAMS%nnyp,BRAMS%nwave)
          if(.not.present(nwave))then
             call i90_perr(myname_,'patch no set!',-1)
             c = undef
             return
          endif

          tmp(1:1, 1:self%nnxp,1:self%nnyp,1:self%nwave) => a
          do j=1,self%nnyp
             do i=1,self%nnxp
                c(i,j) = tmp(1, i, j, nwave)
             enddo
          enddo

    end select

    deallocate(a)

    nx   = size(fld,1)
    ny   = size(fld,2)
    npts = nx*ny

    if(npts .eq. self%nnxp*self%nnyp)then
       ! brams native field size/projection
       fld = c
    else if(npts .eq. self%nLon*self%nLat)then
       ! Do interpolation to a regular lat/lon grid
       fld = undef
       do j=1,self%nLat
          do i=1,self%nLon
             x = self%xMap(i,j)
             y = self%yMap(i,j)
             sumWeight = sum(self%weight(i,j,:))
             if (sumWeight .eq. 0)then
                fld(i,j) = undef
             else
                fld(i,j) = c(x  ,y  ) * self%weight(i,j,1) + &
                           c(x+1,y  ) * self%weight(i,j,2) + &
                           c(x  ,y+1) * self%weight(i,j,3) + &
                           c(x+1,y+1) * self%weight(i,j,4)
             endif
          enddo
       enddo
    else

       write(i_char,'(I4)')nwave
       write(j_char,'(I4)')self%nwave
       write(msg,'(A,1x,2I4)')'unknown requested field size',nx, ny
       call i90_perr(myname_,trim(msg))
       stop
    endif
    deallocate(c)

    return
  end subroutine GetField_
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: GetTimeInfo_ - routine to return some information about the time
!                           of the BRAMS file 
!
! 
! !DESCRIPTION: Esta rotina retorna informacões sobre a data dos arquivos do
!               modelo BRAMS. Podem ser retornadas as seguintes informacões:
!               
!               1. ihr: Hora da condicão inicial utilizada para a simulacão
!               2. iyr: Ano da condicão inicial utilizada para a simulacão
!               3. idy: dia da condicão inicial utilizada para a simulacão
!               4. imo: mês da condicão inicial utilizada para a simulacão
!               5. fhr: hora da previsão da simulacão
!               6. fyr: ano da previsão da simulacão
!               7. fdy: dia da previsão da simulacão
!               8. fmo: mês da previsão da simulacão
!
! !INTERFACE:
!

   function GetTimeInfo_(self, DName) result(dt)
      
      implicit none
!
! !INPUT PARAMETERS:
!      
      ! BRAMS file structure
      class(BramsFile),     intent(in   ) :: self

      ! BRAMS time request
      character(len=*),  intent(in   ) :: DName
      !
      ! Can be:
      !   1. ihr: request hour of initial condition
      !   2. iyr: request year of initial condition
      !   3. idy: request day of initial condition
      !   4. imo: request month of initial condition
      !   5. fhr: request hour of forecast
      !   6. fyr: request year of forecast
      !   7. fdy: request day of forecast
      !   8. fmo: request month of forecast
!
! !OUTPUT PARAMETERS:
!
      ! time of simulation
      integer                          :: dt
!
! !REVISION HISTORY: 
!  11 Nov 2019 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
      character(len=100), parameter :: myname_=':: GetTimeInfo( ... )'

      character(len=10) :: DimName

      if (.not.self%isOpen) then
         call i90_perr(myname,':: BRAMS file not opened yet!',-99)
         return
      endif

      DimName = Trim(Adjustl(I90_lcase(DName)))

      select case (DimName)
         case ('iyr' , 'iyear')
            dt = self%iyr
         case ('imo' , 'imonth')
            dt = self%imo
         case ('idy' , 'iday')
            dt = self%idy
         case ('ihr' , 'ihour')
            dt = self%ihr
         case ('imn' , 'iminute')
            dt = self%imn
         case ('isc' , 'isecond')
            dt = self%isc
         case ('fyr' , 'fyear')
            dt = self%fyr
         case ('fmo' , 'fmonth')
            dt = self%fmo
         case ('fdy' , 'fday')
            dt = self%fdy
         case ('fhr' , 'fhour')
            dt = self%fhr
         case ('fmn' , 'fminute')
            dt = self%fmn
         case ('fsc' , 'fsecond')
            dt = self%fsc
         case ('fct_time')
            dt = self%fct_time
         case default
            write(*,'(4A)')trim(myname_),': wrong dimension <',trim(DimName),'>'
            write(*,'( A)')'try:'
            write(*,'( A)')' * ihr, idy, imo, iyr: Initial Time'
            write(*,'( A)')' * fhr, fdy, fmo, fyr: Forecast Time'
            write(*,'( A)')''
            dt = -1
      end select

   end function
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: GetDim_ - return information about one dimension of BAM file.
!
!
! 
! !DESCRIPTION: Rotina que retorma informacão sobre um das dimensões do arquivo
!               do modelo BAM. Podem ser retornadas informacões sobre:
!               
!               1. nvars    : total number of variables at vfm brams file
!               2. nnxp/nlon: Number of grid cells in the x-direction (longitude)
!               3. nnyp/nlat: Number of grid cells in the y-direction (latitude)
!               4. nnzp/nlev: Number of grid cells in the z-direction (altitude)
!               5. nzg/nsoil: Number of soil layers (deep)
!               6. nzs/nsnow: Maximum number of snowpack layers
!               7. npatch   : Number of subgrid patches
!               8. nwave
!
! !INTERFACE:
!

   function GetDim_(self, DName) result(dim)
      
      implicit none
!
! !INPUT PARAMETERS:
!      
      ! BRAMS file structure
      class(bramsFile),     intent(in   ) :: self

      !BAM dimension name
      character(len=*),  intent(in   ) :: DName
      !
      ! Can be:
      !
      !     1. nvars      - total number of variables at vfm brams file
      !     2. nnxp/nlon  - Number of grid cells in the x-direction
      !     3. nnyp/nlat  - Number of grid cells in the y-direction
      !     4. nnzp/nlev  - Number of grid cells in the z-direction
      !     5. nzg/nsoil  - Number of soil layers
      !     6. nzs/nsnow  - Maximum number of snowpack layers
      !     7. npatch     - Number of subgrid patches
      !     8. nwave
      !

!
! !OUTPUT PARAMETERS:
!
      integer                          :: dim
!
! !REVISION HISTORY: 
!  12 Nov 2019 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!

      character(len=100), parameter :: myname_=':: BAM_GetOneDim( ... )'

      character(len=6) :: DimName

      if (.not.self%isOpen) then
         call i90_perr(myname,':: BRAMS file not opened yet!',-99)
         return
      endif

      DimName = Trim(Adjustl(I90_lcase(DName)))

      select case (DimName)
         case ('nvars')
            dim = self%nvars
         case ('nnxp')
            dim = self%nnxp
         case ('nnyp')
            dim = self%nnyp
         case ('nlon')
            dim = self%nlon
         case ('nlat')
            dim = self%nlat
         case ('nnzp')
            dim = self%nnzp
         case ('nzg'  , 'nsoil')
            dim = self%nzg
         case ('nzs'  , 'nsnow')
            dim = self%nzs
         case ('npatch')
            dim = self%npatch
         case ('nwave')
            dim = self%nwave
         case default
            write(*,'(4A)')trim(myname_),': wrong dimension <',trim(DimName),'>'
            write(*,'( A)')'try: nnxp, nnyp, nnzp, nzg, nzs, npatch, nwave'
            write(*,'( A)')'     nlon, nlat, nlev, nsoil, nsnow'
            dim = -1
      end select

   end function
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: GetMapCoord_ - routine to return informations world coordinates
!                        of BRAMS model file
!
! 
! !DESCRIPTION: Esta rotina retorna um vetor contendo o valor de uma determinada
!               informacão geográfica do modelo BAM. Podem ser retornadas as 
!               seguintes informacões:
!
!               1. glon : real longitude at polar stereographic projection
!               2. glat : real latitude at polar stereographic projection  
!
! !INTERFACE:
!

   subroutine  GetMapCoord_(self, what, coord)
      implicit none
!
! !INPUT PARAMETERS:
!      
      ! BRAMS file structure
      class(BramsFile),  intent(in   ) :: self

      ! Name of inquired Word Coordinate
      character(len=*),  intent(in   ) :: what
      ! 
      ! Can be:
      !        1. glon : real longitude at polar stereographic projection
      !        2. glat : real latitude at polar stereographic projection
      !        3. xtn  : delta X at polar stereographic projection
      !        4  ytn  : delta Y at polar stereographic projection
      !        5. rlon : longitude at regular lat/lon grid'
      !        6. rlat : latitude at regular lat/lon grid'
      !        7. dtx  : delta X at regular lat/lon grid'
      !        8. dty  : delta Y at regular lat/lon grid'

!
! !OUTPUT PARAMETERS:
!
      real(r4),          intent(  out) :: coord(:,:)
!
! !REVISION HISTORY: 
!  13 Oct 2019 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!

      character(len=100), parameter :: myname_=':: GetMapCoord_( ... )'

      integer :: lenIn
      integer :: lenOu
      integer :: i
      integer :: j

      character(len=10) :: DimName

      if (.not.self%isOpen) then
         call i90_perr(myname,':: BRAMS file not opened yet!',-99)
         return
      endif

      DimName = Trim(Adjustl(I90_lcase(what)))

      select case (trim(DimName))
         case('glon') ! real longitude

            !check consistensy
            lenIn = size(coord)
            lenOu = self%nnxp*self%nnyp
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            coord = self%glon


         case('glat') ! real latitude

            !check consistensy
            lenIn = size(coord)
!            lenOu = BFile%head%JMax
            lenOu = self%nnxp*self%nnyp
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            coord = self%glat

         case('xtn') ! real longitude

            !check consistensy
            lenIn = size(coord)
            lenOu = self%nnxp*self%nnyp
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            do j=1,self%nnyp
               do i=1,self%nnxp
                  coord(i,j) = self%xtn(i)
               enddo
            enddo

         case('ytn') ! real latitude

            !check consistensy
            lenIn = size(coord)
!            lenOu = BFile%head%JMax
            lenOu = self%nnxp*self%nnyp
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            do j=1,self%nnyp
               do i=1,self%nnxp
                  coord(i,j) = self%ytn(i)
               enddo
            enddo

         case('dxt')

            !check consistensy
            lenIn = size(coord)
            lenOu = self%nnxp*self%nnyp
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            coord = self%dxt

         case('dyt')

            !check consistensy
            lenIn = size(coord)
!            lenOu = BFile%head%JMax
            lenOu = self%nnxp*self%nnyp
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            coord = self%dyt

         case('fmapt')

            !check consistensy
            lenIn = size(coord)
!            lenOu = BFile%head%JMax
            lenOu = self%nnxp*self%nnyp
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            coord = self%fmapt


         case('rlon') ! longitude at regula lat/lon grid

            !check consistensy
            lenIn = size(coord)
            lenOu = self%nLon*self%nLat
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            do j=1,self%nLat
               do i=1,self%nLon
                  coord(i,j) = self%lon(i)
               enddo
            enddo

         case('rlat') ! latitude at regular lat/lon grid

            !check consistensy
            lenIn = size(coord)
!            lenOu = BFile%head%JMax
            lenOu = self%nLon*self%nLat
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            do j=1,self%nLat
               do i=1,self%nLon
                  coord(i,j) = self%lat(j)
               enddo
            enddo

         case('xmap') 

            !check consistensy
            lenIn = size(coord)
!            lenOu = BFile%head%JMax
            lenOu = self%nLon*self%nLat
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            coord = self%xMap

         case('ymap') 

            !check consistensy
            lenIn = size(coord)
!            lenOu = BFile%head%JMax
            lenOu = self%nLon*self%nLat
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            coord = self%yMap


         case default
            write(stdout,'(2(1x,A))')trim(myname_),': ERROR'
            write(stdout,'(2(1x,A))')'Wrong coordinate inquired :',trim(what)
            write(stdout,'(2(1x,A))')'validy coordinates:'
            write(stdout,'(2(1x,A))')'glon : real longitude at polar stereographic projetion'
            write(stdout,'(2(1x,A))')'glat : real latitude at polar stereographic proction'
            write(stdout,'(2(1x,A))')'xtn  : delta X at polar stereographic projection'
            write(stdout,'(2(1x,A))')'ytn  : delta Y at polar stereographic projection'
            write(stdout,'(2(1x,A))')'dtx  : delta X at t points'
            write(stdout,'(2(1x,A))')'dty  : delta Y at t points'
            write(stdout,'(2(1x,A))')'rlon : longitude at regular lat/lon grid'
            write(stdout,'(2(1x,A))')'rlat : latitude at regular lat/lon grid'

            coord = undef
            return 

      end select

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
! !IROUTINE: GetVerticalInfo_ - routine to return informations about vertical
!                               coordinates of BRAMS model file
!
! 
! !DESCRIPTION: Esta rotina retorna um vetor contendo o valor de uma determinada
!               informacão sobre a coordenada vertical modelo BRAMS.
!               Podem ser retornadas as seguintes informacões:
!
!               1. dzt  : inverso da espessura da camada
!
! !INTERFACE:
!

   subroutine  GetVerticalInfo_(self, what, coord)
      implicit none
!
! !INPUT PARAMETERS:
!      
      ! BRAMS file structure
      class(BramsFile),  intent(in   ) :: self

      ! Name of inquired Word Coordinate
      character(len=*),  intent(in   ) :: what
      ! 
      ! Can be:
      !        1. dzt  : inverse of the atmosphere thickness

!
! !OUTPUT PARAMETERS:
!
      real(r4),          intent(  out) :: coord(:)
!
! !REVISION HISTORY: 
!  18 Nov 2020 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
      character(len=100), parameter :: myname_=':: GetVertivalInfo_( ... )'

      integer :: lenIn
      integer :: lenOu
      integer :: i
      integer :: j

      character(len=10) :: DimName

      if (.not.self%isOpen) then
         call i90_perr(myname,':: BRAMS file not opened yet!',-99)
         return
      endif

      DimName = Trim(Adjustl(I90_lcase(what)))

      select case (trim(DimName))
         case('dzt') ! inverse of the atmosphere thickness

            !check consistensy
            lenIn = size(coord)
            lenOu = self%nnzp
            if (lenIn.ne.lenOu) then
               write (stdout,'(2(1x,A))')trim(myname_),':ERROR, wrong size of coordinate array:'
               write (stdout,'(I4,1x,A,1x,I4)')lenIn,' .ne. ',lenOu
               return
            endif

            coord = self%dzt

         case default
            write(stdout,'(2(1x,A))')trim(myname_),': ERROR'
            write(stdout,'(2(1x,A))')'Wrong vertical info inquired :',trim(what)
            write(stdout,'(2(1x,A))')'validy vertical info:'
            write(stdout,'(2(1x,A))')'dzt : inverse of atmosphere thickness'

            coord = undef
            return 

      end select

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
! !IROUTINE: GetMapInfo_ - routine to return informations world coordinates
!                        of BRAMS model file
!
! 
! !DESCRIPTION: Esta rotina retorna o valor de uma determinada informação geográfica
!               do modelo BRAMS. Podem ser retornadas as seguintes informações:
!               
!
!               1. lonMin : start point of longitude at polar stereographic projection
!               2. lonMax : end point of longitude at polar stereographic projection  
!               3. latMin : start point of latitude at polar stereographic projection
!               4. latMax : end point of latitude at polar stereographic projection
!               5. deltaxn: delta longitude of polar stereographic projection;
!               6. deltayn: delta latitude ofpolar stereographic projection;
!               7. delLon : delta longitude of regular lat/lon grid;
!               8. delLat : delta latitude of regular lat/lon grid;
!               9. ztop   : altitude of top, km above sea level 
!

!
! !INTERFACE:
!

   function  GetMapInfo_(self, what) result(coord)
      implicit none
!
! !INPUT PARAMETERS:
!      
      ! BRAMS file structure
      class(BramsFile),     intent(in   ) :: self

      ! Name of inquired Word Coordinate
      character(len=*),  intent(in   ) :: what
      ! 
      ! Can be:
      !        1. lonMin : start point of longitude at polar stereographic projection
      !        2. lonMax : end point of longitude at polar stereographic projection  
      !        3. latMin : start point of latitude at polar stereographic projection
      !        4. latMax : end point of latitude at polar stereographic projection 
      !        5. deltaxn: delta longitude of polar stereographic projection;
      !        6. deltayn: delta latitude ofpolar stereographic projection;
      !        7. delLon : delta longitude of regular lat/lon grid;
      !        8. delLat : delta latitude of regular lat/lon grid;
      !        9. ztop   : altitude of top, km above sea level;
!
! !OUTPUT PARAMETERS:
!
      real(r4)                         :: coord
!
! !REVISION HISTORY: 
!  13 Oct 2019 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!

      character(len=100), parameter :: myname_=':: GetMapInfo__( ... )'

      character(len=10) :: DimName

      if (.not.self%isOpen) then
         call i90_perr(myname,':: BRAMS file not opened yet!',-99)
         return
      endif

      DimName = Trim(Adjustl(I90_lcase(what)))

      select case (trim(DimName))
         case('lonmin' , 'lon1')
            coord = self%lonMin
         case('lonmax')
            coord = self%lonMax
         case('latmin' , 'lat1')
            coord = self%latMin
         case('latmax')
            coord = self%latMax
         case('deltaxn')
            coord = self%deltaxn
         case('deltayn')
            coord = self%deltayn
         case('dellon')
            coord = self%delLon
         case('dellat')
            coord = self%delLat
         case('ztop')
            coord = self%ztop
         case default
            write(stdout,'(2(1x,A))')trim(myname_),': ERROR'
            write(stdout,'(2(1x,A))')'Wrong coordinate inquired :',trim(what)
            write(stdout,'(2(1x,A))')'validy coordinates:'
            write(stdout,'(2(1x,A))')' '
            write(stdout,'(2(1x,A))')'lonmin : start point of longitude at polar stereographic projetion'
            write(stdout,'(2(1x,A))')'lonmax : end point of longitude at polar stereographic projetion'
            write(stdout,'(2(1x,A))')'latmin : start point of latitude at polar stereographic proction'
            write(stdout,'(2(1x,A))')'latmax : end point of latitude at polar stereographic proction'
            write(stdout,'(2(1x,A))')'deltaxn: delta longitude of polar stereographic projection'
            write(stdout,'(2(1x,A))')'deltayn: delta latitude ofpolar stereographic projection'
            write(stdout,'(2(1x,A))')'delLon : delta longitude of regular lat/lon grid'
            write(stdout,'(2(1x,A))')'delLat : delta latitude of regular lat/lon grid'
            write(stdout,'(2(1x,A))')'ztop   : altitude of top, km above sea level'

            coord = undef
            return 

      end select

   end function
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IFUNCTION: GetNLevels_ - return how many vertical levels has one
!                           givem variable
!
! 
! !DESCRIPTION: Esta funcão retorna o número de níveis verticais que uma dada
!               variável possui.
!
!
!
! !INTERFACE:
!

   function GetNlevels_(self,VName,istat) result(nlevs)

!
! !INPUT PARAMETERS:
!      
      ! BRAMS file structure
      class(bramsFile),   intent(in   ) :: self

      ! Name of a BAM variable
      character(len=*),   intent(in   ) :: VName
      
!
! !OUTPUT PARAMETERS:
!
      integer, optional,  intent(  out) :: istat
!
! !RETURN VALUE:
!
      integer                           :: NLevs
!
! !REVISION HISTORY: 
!  21 Nov 2019 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
      character(len=100), parameter :: myname_=':: GetNLevels_( ... )'

      character(len=strlen)  :: VarName
      type(bramsHeader), pointer :: Head => null()

      if (.not.self%isOpen) then
         call i90_perr(myname,':: BRAMS file not opened yet!',-99)
         return
      endif
      

      if(present(istat)) istat = 0

      VarName = trim(adjustl(I90_lcase(VName)))
      
      Head => self%Header
      do while(associated(Head))
         if(trim(VarName).eq.trim(Head%VarName))exit
         Head=>Head%next
      enddo
  
      if(.not.associated(Head))then
         write(*,'(4A)')trim(myname_),': variable not found ! <',trim(VarName),'>'
         if(present(istat)) istat = -1
         return
      endif

      select case(Head%idim_type)
         case(2)
           NLevs = 1 
         case(3)
           NLevs = self%nnzp
         case(4)
           NLevs = self%nzg   ! o que fazer com o npatch nesse caso?
         case(5)
           NLevs = self%nzs   ! o que fazer com o npatch nesse caso?
         case(6)
           NLevs = self%npatch ! sera que considera como nivel?
         case(7)
           NLevs = self%nwave  ! sera que considera como nivel?
      end select

      return

   end function

   subroutine GetInterpWeight_( self, weight )
      class(bramsfile)        :: self
      real(r4), intent(inout) :: weight(:,:,:)

      weight = self%weight

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
! !IROUTINE: GetFloat_ - return information about one dimension of BRAMS file.
!
!
! 
! !DESCRIPTION: Esta rotina busca pelo label no arquivo Header do BRAMS
!               e retorna o valor real atribuído ao label buscado.
!               
!
! !INTERFACE:
!
  subroutine GetFloat_(label, val)

!
! !INPUT PARAMETERS:
!      

    character(len=*), intent(in   ) :: label
!
! !OUTPUT PARAMETERS:
!     
    real,             intent(  out) :: val
!
!
! !SEE ALSO:
!   function/routines from inpak90 module
!
!    - I90_Label( )  - selects a label (key)
!    - I90_Gint( )   - return next integer number (function)
!    - I90_GFloat( ) - return next float number (function)
!    - I90_GLine( )  - selects next line
!
! !REVISION HISTORY: 
!  12 Nov 2019 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
    character(len=*),parameter :: myname_= myname//' :: GetFloat_( ... )'
    character(len=256) :: msg

    integer :: iret
    integer :: nlines

    call I90_Label( trim(label), iret )
    if(iret .ne. 0)then
       write(msg,'(3A)')'i90_label("',trim(label),'")'
       call i90_perr(myname_,trim(msg),iret)
    endif

    call I90_GLine( iret )
    if(iret .ne. 0)then
       write(msg,'(3A)')'i90_GLine("',trim(label),'")'
       call i90_perr(myname_,trim(msg),iret)
    endif

    nlines = I90_GInt( iret )
    if(iret .ne. 0)then
       write(msg,'(3A)')'i90_GLine("',trim(label),'")'
       call i90_perr(myname_,trim(msg),iret)
    endif

    call I90_GLine( iret )
    if(iret .ne. 0)then
       write(msg,'(3A)')'i90_GLine("',trim(label),'")'
       call i90_perr(myname_,trim(msg),iret)
    endif

    val=I90_GFloat( iret )

    if(iret .ne. 0)then
       write(msg,'(3A)')'i90_GFloat("',trim(label),'")'
       call i90_perr(myname_,trim(msg),iret)
    endif

    return
  end subroutine GetFloat_
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: GetInt_ - return information about one dimension of BRAMS file.
!
!
! 
! !DESCRIPTION: Esta rotina busca pelo label no arquivo Header do BRAMS
!               e retorna o valor inteiro atribuído ao label buscado.
!               
!
! !INTERFACE:
!
  subroutine GetInt_(label, val)
!
! !INPUT PARAMETERS:
! 
    character(len=*), intent(in   ) :: label
!
! !OUTPUT PARAMETERS:
!  
    integer,          intent(  out) :: val
!
!
! !SEE ALSO:
!   function/routines from inpak90 module
!
!    - I90_Label( )  - selects a label (key)
!    - I90_Gint( )   - return next integer number (function)
!    - I90_GFloat( ) - return next float number (function)
!    - I90_GLine( )  - selects next line
!
! !REVISION HISTORY: 
!  12 Nov 2019 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
    character(len=*),parameter :: myname_= myname//' :: GetInt_( ... )'
    character(len=256) :: msg

    integer :: iret
    integer :: nlines

    call I90_Label( trim(label), iret )
    if(iret .ne. 0)then
       write(msg,'(3A)')'i90_label("',trim(label),'")'
       call i90_perr(myname_,trim(msg),iret)
    endif

    call I90_GLine( iret )
    if(iret .ne. 0)then
       write(msg,'(3A)')'i90_GLine("',trim(label),'")'
       call i90_perr(myname_,trim(msg),iret)
    endif

    nlines = I90_GInt( iret )
    if(iret .ne. 0)then
       write(msg,'(3A)')'i90_GLine("',trim(label),'")'
       call i90_perr(myname_,trim(msg),iret)
    endif

    call I90_GLine( iret )
    if(iret .ne. 0)then
       write(msg,'(3A)')'i90_GLine("',trim(label),'")'
       call i90_perr(myname_,trim(msg),iret)
    endif

    val=I90_Gint( iret )
    if(iret .ne. 0)then
       write(msg,'(3A)')'i90_Gint("',trim(label),'")'
       call i90_perr(myname_,trim(msg),iret)
    endif

    return
  end subroutine GetInt_
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: GetArrayFloat_ - return information about one dimension of BRAMS file.
!
!
! 
! !DESCRIPTION: Esta rotina busca pelo label no arquivo Header do BRAMS
!               e retorna um array com os valores reais atribuídos ao label buscado.
!               
!
! !INTERFACE:
!
  subroutine GetArrayFloat_(label, val)
!
! !INPUT PARAMETERS:
!  
    character(len=*), intent(in   ) :: label
!
! !OUTPUT PARAMETERS:
!   
    real,             intent(  out) :: val(:)
!
!
! !SEE ALSO:
!   function/routines from inpak90 module
!
!    - I90_Label( )  - selects a label (key)
!    - I90_Gint( )   - return next integer number (function)
!    - I90_GFloat( ) - return next float number (function)
!    - I90_GLine( )  - selects next line
!
! !REVISION HISTORY: 
!  12 Nov 2019 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
    character(len=*),parameter :: myname_= myname//' :: GetFloat_( ... )'
    character(len=256) :: msg

    integer :: iret
    integer :: nlines
    integer :: vsize
    integer :: i

    call I90_Label( trim(label), iret )
    if(iret .ne. 0)then
       write(msg,'(3A)')'i90_label("',trim(label),'")'
       call i90_perr(myname_,trim(msg),iret)
    endif

    call I90_GLine( iret )
    if(iret .ne. 0)then
       write(msg,'(3A)')'i90_GLine("',trim(label),'")'
       call i90_perr(myname_,trim(msg),iret)
    endif

    nlines = I90_GInt( iret )
    if(iret .ne. 0)then
       write(msg,'(3A)')'i90_GLine("',trim(label),'")'
       call i90_perr(myname_,trim(msg),iret)
    endif

    vsize  = size(val)
    if(nlines.ne.vsize)then
       print*, 'Array with different size:'
       print*, 'Requested ', vsize
       print*, trim(label)//' has', nlines
       val = -999.9
       return
    endif

    do i=1,nlines
       call I90_GLine( iret )
       if(iret .ne. 0)then
          write(msg,'(3A)')'i90_GLine("',trim(label),'")'
          call i90_perr(myname_,trim(msg),iret)
       endif
   
       val(i)=I90_GFloat( iret )
       if(iret .ne. 0)then
          write(msg,'(3A)')'i90_GFloat("',trim(label),'")'
          call i90_perr(myname_,trim(msg),iret)
       endif
    enddo

    return
  end subroutine GetArrayFloat_
!
!EOC
!
!-----------------------------------------------------------------------------!
!             Modeling and Development Division - DMD/CPTEC/INPE              !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: GetArrayInt_ - return information about one dimension of BRAMS file.
!
!
! 
! !DESCRIPTION: Esta rotina busca pelo label no arquivo Header do BRAMS
!               e retorna um array com os valores inteiros atribuídos ao label buscado.
!               
!
! !INTERFACE:
!
  subroutine GetArrayInt_(label, val)
!
! !INPUT PARAMETERS:
! 
    character(len=*),     intent(in   ) :: label
!
! !OUTPUT PARAMETERS:
! 
    integer, Allocatable, intent(  out) :: val(:)
!
!
! !SEE ALSO:
!   function/routines from inpak90 module
!
!    - I90_Label( )  - selects a label (key)
!    - I90_Gint( )   - return next integer number (function)
!    - I90_GFloat( ) - return next float number (function)
!    - I90_GLine( )  - selects next line
!
! !REVISION HISTORY: 
!  12 Nov 2019 - J. G. de Mattos - Initial Version
!
!
!EOP
!-----------------------------------------------------------------------------!
!BOC
!
    character(len=*),parameter :: myname_= myname//' :: GetInt_( ... )'
    character(len=256) :: msg

    integer :: iret
    integer :: nlines
    integer :: vsize
    integer :: i

    call I90_Label( trim(label), iret )
    if(iret .ne. 0)then
       write(msg,'(3A)')'i90_label("',trim(label),'")'
       call i90_perr(myname_,trim(msg),iret)
    endif

    call I90_GLine( iret )
    if(iret .ne. 0)then
       write(msg,'(3A)')'i90_GLine("',trim(label),'")'
       call i90_perr(myname_,trim(msg),iret)
    endif

    nlines = I90_GInt( iret )
    if(iret .ne. 0)then
       write(msg,'(3A)')'i90_GLine("',trim(label),'")'
       call i90_perr(myname_,trim(msg),iret)
    endif
    
    vsize  = size(val)

    if(nlines.ne.vsize)then
       print*, 'Array with different size:'
       print*, 'Requested ', vsize
       print*, trim(label)//' has', nlines
       val = -999
       return
    endif

    do i=1,nlines
       call I90_GLine( iret )
       if(iret .ne. 0)then
          write(msg,'(3A)')'i90_GLine("',trim(label),'")'
          call i90_perr(myname_,trim(msg),iret)
       endif
   
       val(i)=I90_Gint( iret )
       if(iret .ne. 0)then
          write(msg,'(3A)')'i90_Gint("',trim(label),'")'
          call i90_perr(myname_,trim(msg),iret)
       endif
    enddo

    return
  end subroutine GetArrayInt_
!
!EOC
!
!-----------------------------------------------------------------------------!
! ATENTION: Function to obtain derived variables
!-----------------------------------------------------------------------------!
  subroutine virtualTemperature(self, ilev, vtmp)
     type(bramsfile), intent(in   ) :: self
     integer,         intent(in   ) :: ilev
     real,            intent(inout) :: vtmp(:,:)

     real, allocatable :: theta(:,:)
     real, allocatable :: tempk(:,:)
     real, allocatable :: pi(:,:)
     real, allocatable :: rv(:,:)
     real, allocatable :: q(:,:)
     integer :: nx, ny

     character(len=512) :: myname_='virtualTemperature()'

     nx = size(vtmp,1)
     ny = size(vtmp,2)

     allocate(rv(nx,ny))
     call self%GetField(rv, 'RV',level=ilev)

     allocate(q(nx,ny))

     q = rv / (1.0 + rv)
     where(rv.eq.undef) q = undef

     deallocate(rv)

     allocate(theta(nx,ny))
     call self%GetField(theta, 'THETA',level=ilev)

     allocate(pi(nx,ny))
     call self%GetField(pi, 'PI',level=ilev)

     allocate(tempk(nx,ny))
     tempk = theta * pi / cp
     where(theta.eq.undef .or. pi.eq.undef) tempk = undef

     deallocate(theta)
     deallocate(pi)

     vtmp = tempk * ( 1 + 0.61 * q) 
     where(tempk.eq.undef)vtmp = undef

     deallocate(tempk)

  end subroutine

  subroutine SpecificHumidity(self, ilev, umes)
     type(bramsfile), intent(in   ) :: self
     integer,         intent(in   ) :: ilev
     real,            intent(inout) :: umes(:,:)

     real, allocatable :: rv(:,:)
     integer :: nx, ny
     character(len=512) :: myname_='SpecificHumidity()'

     nx = size(umes,1)
     ny = size(umes,2)

     allocate(rv(nx,ny))
     call self%GetField(rv, 'RV',level=ilev)

     umes = rv / (1.0 + rv)
     where(rv.eq.undef) umes = undef

     deallocate(rv)

  end subroutine


  subroutine SpecificHumidity_2m(self, ilev, umes)
     type(bramsfile), intent(in   ) :: self
     integer,         intent(in   ) :: ilev
     real,            intent(inout) :: umes(:,:)

     real, allocatable :: rv(:,:)
     integer :: nx, ny
     character(len=512) :: myname_='SpecificHumidity_2m()'

     nx = size(umes,1)
     ny = size(umes,2)

     allocate(rv(nx,ny))
     call self%GetField(rv, 'RV2MJ',level=ilev)


     umes = rv / (1.0 + rv)
     where(rv.eq.undef) umes = undef

     deallocate(rv)

  end subroutine


  subroutine CloudMixRatio(self, ilev, cloud)
     type(bramsfile), intent(in   ) :: self
     integer,         intent(in   ) :: ilev
     real,            intent(inout) :: cloud(:,:)

     logical, allocatable :: bitMap(:,:)
     integer :: nx, ny
     character(len=512) :: myname_='CloudMixRatio()'

     nx = size(cloud,1)
     ny = size(cloud,2)

     call self%GetField(cloud, 'RCP', level=ilev)
     
     allocate(bitMap(nx,ny))
     bitMap = .false.
     where(cloud .eq. undef) bitMap = .true.

     cloud = cloud * 1.e3
     cloud = max(cloud, 0.0)

     where(bitMap) cloud = undef
     deallocate(bitMap)

  end subroutine

  subroutine GraupelMixRatio(self, ilev, graupel)
     type(bramsfile), intent(in   ) :: self
     integer,         intent(in   ) :: ilev
     real,            intent(inout) :: graupel(:,:)

     logical, allocatable :: bitMap(:,:)
     integer :: nx,ny
     character(len=512) :: myname_='GraupeldMixRatio()'

     nx = size(graupel,1)
     ny = size(graupel,2)

     call self%GetField(graupel, 'RGP', level=ilev)
     
     allocate(bitMap(nx,ny))
     bitMap = .false.
     where(graupel .eq. undef) bitMap = .true.

     graupel = graupel * 1.e3
     graupel = max(graupel, 0.0)

     where(bitMap) graupel = undef
     deallocate(bitMap)

  end subroutine

  subroutine SnowMixRatio(self, ilev, snow)
     type(bramsfile), intent(in   ) :: self
     integer,         intent(in   ) :: ilev
     real,            intent(inout) :: snow(:,:)

     logical, allocatable :: bitMap(:,:)
     integer :: nx, ny
     character(len=512) :: myname_='SnowdMixRatio()'

     nx = size(snow,1)
     ny = size(snow,2)

     call self%GetField(snow, 'RSP', level=ilev)
     
     allocate(bitMap(nx,ny))
     bitMap = .false.
     where(snow .eq. undef) bitMap = .true.

     snow = snow * 1.e3
     snow = max(snow, 0.0)

     where(bitMap) snow = undef
     deallocate(bitMap)

  end subroutine

  subroutine RainMixRatio(self, ilev, rain)
     type(bramsfile), intent(in   ) :: self
     integer,         intent(in   ) :: ilev
     real,            intent(inout) :: rain(:,:)

     logical,allocatable :: bitMap(:,:)
     integer :: nx, ny
     character(len=512) :: myname_='RaindMixRatio()'

     nx = size(rain,1)
     ny = size(rain,2)

     call self%GetField(rain, 'RRP', level=ilev)
     
     allocate(bitmap(nx,ny))
     bitMap = .false.
     where(rain .eq. undef) bitMap = .true.

     rain = rain * 1.e3
     rain = max(rain, 0.0)

     where(bitMap) rain = undef
     deallocate(bitMap)

  end subroutine


  subroutine IceMixRatio(self, ilev, ice)
     type(bramsfile), intent(in   ) :: self
     integer,         intent(in   ) :: ilev
     real,            intent(inout) :: ice(:,:)

     real, allocatable :: RPP(:,:)
     real, allocatable :: RSP(:,:)
     real, allocatable :: RAP(:,:)
     real, allocatable :: RGP(:,:)
     real, allocatable :: RHP(:,:)
     real, allocatable :: Q6(:,:)
     real, allocatable :: Q7(:,:)
     integer :: npts
     integer :: nx, ny
     character(len=512) :: myname_='IcedMixRatio()'

     nx = size(ice,1)
     ny = size(ice,2)

     ice = 0.0

     allocate(RPP(nx,ny))
     call self%GetField(RPP, 'RPP',level=ilev)
     ice = ice + RPP
     where(RPP.eq.undef) ice = undef
     deallocate(RPP)

     allocate(RSP(nx,ny))
     call self%GetField(RSP, 'RSP',level=ilev)
     ice = ice + RSP
     where(RSP.eq.undef) ice = undef
     deallocate(RSP)

     if(self%mcphys_type .le. 1)then
        allocate(RAP(nx,ny))
        call self%GetField(RAP, 'RAP',level=ilev)
        ice = ice + RAP
        where(RAP.eq.undef) ice = undef
        deallocate(RAP)
     endif

     allocate(RGP(nx,ny))
     call self%GetField(RGP, 'RGP', level=ilev)

     if(self%mcphys_type .le. 1)then
        allocate(Q6(nx,ny))
        call self%GetField(Q6, 'Q6',level=ilev)
        Q6  = max(Q6,      0.0)
        Q6  = min(Q6, 334000.0)
        Q6  = 1. - ( Q6 / 334000.0 )
        where(RGP .eq. undef) Q6 = undef
        RGP = RGP * Q6
        where(Q6.eq.undef) RGP = undef
        deallocate(Q6)
     endif
     ice = ice + RGP
     where(RGP .eq. undef) ice = undef
     deallocate(RGP)
    
     if(self%mcphys_type .le. 1) then
        allocate(Q7(nx,ny))
        call self%GetField(Q7, 'Q7',level=ilev)
        Q7  = max(Q7,      0.0)
        Q7  = min(Q7, 334000.0)
        Q7  = 1. - ( Q7 / 334000.0 )

        allocate(RHP(nx,ny))
        call self%GetField(RHP, 'RHP', level=ilev)

        where(RHP .eq. undef) Q7 = undef
        RHP = RHP * Q7
        where(Q7.eq.undef) RHP = undef
        deallocate(Q7)

        ice = ice + RHP
        where(RHP .eq. undef) ice = undef
        deallocate(RHP)

     endif

     ice = ice * 1.e3
     ice = max(ice, 0.0)

  end subroutine

  subroutine CompSfcPress(self, ilev, sfcPress)

     type(bramsfile), intent(in   ) :: self
     integer,         intent(in   ) :: ilev
     real,            intent(inout) :: sfcPress(:,:) ! mb

     real, allocatable :: pi01(:,:)
     real, allocatable :: pi02(:,:)
     integer :: nx, ny
     character(len=512) :: myname_='SfcPress()'
     
     nx = size(sfcPress,1)
     ny = size(sfcPress,2)

     allocate(pi01(nx,ny))
     call self%GetField(PI01, 'PI',level=1)

     allocate(pi02(nx,ny))
     call self%GetField(PI02, 'PI',level=2)

     sfcPress = (0.5 * (PI01 + PI02) / cp)**cpor * p00 * .01

     where(PI01.eq.undef) sfcPress = undef
     where(PI02.eq.undef) sfcPress = undef

  end subroutine

   subroutine registerField_(self, FName, idim_type, fnct)
      class(bramsfile), intent(inout) :: self
      integer,          intent(in   ) :: idim_type
      character(len=*), intent(in   ) :: FName
      procedure(FieldMethod)          :: fnct

      type(bramsHeader), pointer :: Header => null()
      type(bramsHeader), pointer :: Find => null()
      character(len=strlen) :: VarName

      VarName = trim(adjustl(I90_lcase(FName)))      

      ! procura se o FName já está associado
      Find => self%Header
      do while(associated(Find))
         if(trim(Find%VarName) .eq. trim(VarName))then
            print*,'erro: Function Name already exist:',trim(VarName)
            stop
         endif
         Header => Find
         Find => Find%next
      enddo

      !se não encontrou continua

      self%nvars_extra = self%nvars_extra + 1

      allocate(Header%next)
      Header => Header%next

      Header%VarName   = trim(VarName)
      Header%idim_type = idim_type
      Header%extra     = .true.
      Header%func      => Fnct

      return
   end subroutine

  subroutine xy_ll(qlat,qlon,polelat,polelon,x,y)
    real, intent(OUT) :: qlat
    real, intent(OUT) :: qlon
    real, intent(IN)  :: polelat
    real, intent(IN)  :: polelon
    real, intent(IN)  :: x
    real, intent(IN)  :: y
    real              :: sinplat,sinplon
    real              :: cosplat,cosplon
    real              :: x3p,y3p,z3p
    real              :: x3q,y3q,z3q,r3q
    real              :: xq,yq,zq
    real              :: d,t
    real              :: alpha

    !    include "post_rconfig.h"

    ! Evaluate sine and cosine of latitude and longitude of pole point p.
    sinplat = sin (polelat * deg2rad)
    cosplat = cos (polelat * deg2rad) 
    sinplon = sin (polelon * deg2rad)
    cosplon = cos (polelon * deg2rad)

    ! Compute (x3,y3,z3) coordinates of the pole point where the origin is the
    ! center of the earth, the z axis is the north pole, the x axis is the
    ! equator and prime meridian, and the y axis is the equator and 90 E.
    x3p     = earthRadius * cosplat * cosplon
    y3p     = earthRadius * cosplat * sinplon
    z3p     = earthRadius * sinplat

    ! Compute distance d from given point R on the polar stereographic plane
    ! to the pole point P:
    d       = sqrt (x ** 2 + y ** 2)

    ! Compute angle QCP where C is the center of the Earth.  This is twice
    ! angle QAP where A is the antipodal point.  Angle QAP is the same as
    ! angle RAP:

    alpha   = 2. * atan2 (d,real(earthDiameter,r4))

    ! Compute zq, the height of Q relative to the polar stereographic plane:
    zq      = earthRadius * (cos (alpha) - 1.)

    ! Compute the parameter t which is the the distance ratio AQ:AR
    t       = (earthDiameter + zq) / earthDiameter

    ! Compute xq and yq, the x and y coordinates of Q in polar stereographic space:
    xq      = t * x
    yq      = t * y

    ! Transform location of Q from (x,y,z) coordinates to (x3,y3,z3):
    x3q     = x3p - xq * sinplon - yq * cosplon * sinplat  &
         + zq * cosplat * cosplon
    y3q     = y3p + xq * cosplon - yq * sinplon * sinplat  &
         + zq * cosplat * sinplon
    z3q     = z3p + yq * cosplat + zq * sinplat

    ! Compute the latitude and longitude of Q:
    qlon    = atan2 (y3q,x3q) / deg2rad
    r3q     = sqrt (x3q ** 2 + y3q ** 2)

    qlat    = atan2 (z3q,r3q) / deg2rad

  end subroutine xy_ll

  subroutine ll_xy(qlat, qlon, polelat, polelon, x, y)
    real, intent(in) :: qlat
    real, intent(in) :: qlon
    real, intent(in) :: polelat
    real, intent(in) :: polelon
    real, intent(out) :: x
    real, intent(out) :: y
    real :: sinplat, sinplon, sinqlat, sinqlon
    real :: cosplat, cosplon, cosqlat, cosqlon
    real :: x3p, y3p, z3p
    real :: x3q, y3q, z3q
    real :: xq, yq, zq
    real :: t

    !      include "post_rconfig.h"
    !-----------------------------------------------------------

    ! Evaluate sine and cosine of latitude and longitude of pole point p and
    ! input point q.
    sinplat = sin (polelat * deg2rad)
    cosplat = cos (polelat * deg2rad)
    sinplon = sin (polelon * deg2rad)
    cosplon = cos (polelon * deg2rad)

    sinqlat = sin (qlat * deg2rad)
    cosqlat = cos (qlat * deg2rad)
    sinqlon = sin (qlon * deg2rad)
    cosqlon = cos (qlon * deg2rad)

    ! Compute (x3,y3,z3) coordinates where the origin is the center of the earth,
    ! the z axis is the north pole, the x axis is the equator and prime
    ! meridian, and the y axis is the equator and 90 E.

    ! For the pole point, these are:
    x3p = earthRadius * cosplat * cosplon
    y3p = earthRadius * cosplat * sinplon
    z3p = earthRadius * sinplat

    ! For the given lat,lon point, these are:
    z3q = earthRadius * sinqlat
    x3q = earthRadius * cosqlat * cosqlon
    y3q = earthRadius * cosqlat * sinqlon

    ! Transform q point from (x3,y3,z3) coordinates in the above system to
    ! polar stereographic coordinates (x,y,z):
    xq = - sinplon * (x3q - x3p) + cosplon * (y3q - y3p)
    yq = cosplat * (z3q - z3p)                                      &
         - sinplat * (cosplon * (x3q - x3p) + sinplon * (y3q - y3p))
    zq = sinplat * (z3q - z3p)                                      &
         + cosplat * (cosplon * (x3q - x3p) + sinplon * (y3q - y3p))

    ! Parametric equation for line from antipodal point at (0,0,-2 earthRadius) to
    ! point q has the following parameter (t) value on the polar stereographic
    ! plane:
    t = earthDiameter / (earthDiameter + zq)

    ! This gives the following x and y coordinates for the projection of point q
    ! onto the polar stereographic plane:
    x = xq * t
    y = yq * t

  end subroutine ll_xy

  real function distang(dlat1,dlon1,dlat2,dlon2)

    real, intent(in) :: dlat1,dlon1,dlat2,dlon2  ! input coordinates in degrees

    real :: lat1,lon1,lat2,lon2 ! coordinates in radians
    real :: xxd, xxm
    if (dlon1.eq.dlon2.and.dlat1.eq.dlat2) then
       distang=0.
       return
    endif
    ! convert to radians
    lon1=dlon1*deg2rad
    lat1=dlat1*deg2rad
    lon2=dlon2*deg2rad
    lat2=dlat2*deg2rad

    xxd=cos(lon1-lon2)*cos(lat1)*cos(lat2) + sin(lat1)*sin(lat2)
    xxm=amin1(1.0,amax1(-1.0,xxd))
    distang=acos(xxm)*earthRadius

    return
  end function distang

  function haversine(ilon1, ilat1, ilon2, ilat2) result(km)

    real,intent(in) :: ilon1, ilat1
    real,intent(in) :: ilon2, ilat2

    real :: lon1, lat1
    real :: lon2, lat2
    real :: dlon, dlat
    real :: a, c
    real :: km
    !    """
    !    Calculate the great circle distance between two points 
    !    on the earth (specified in decimal degrees)
    !    """
    ! convert decimal degrees to radians 

    lon1=ilon1*deg2rad
    lat1=ilat1*deg2rad
    lon2=ilon2*deg2rad
    lat2=ilat2*deg2rad

    ! haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = (sin(dlat/2)*sin(dlat/2)) + cos(lat1) * cos(lat2) * (sin(dlon/2)*sin(dlon/2))
    c = 2 * asin(sqrt(a)) 
    km = earthRadius * c
    return
  end function haversine

  !
  ! Variables
  !
  !    subroutine InsertRealScalar_(var, vname, data)
  !       type(variable),   intent(inout) :: var
  !       character(len=*), intent(in   ) :: vname
  !       real,           intent(in   )   :: data
  !
  !       allocate(var)
  !
  !       
  !
  !    end subroutine

  subroutine toll()
  end subroutine toll

  subroutine BuildAxis(nVal, firstVal, delVal, val)
    integer, intent(in) :: nVal
    real, intent(in) :: firstVal
    real, intent(in) :: delVal
    real, intent(out) :: val(:)

    integer :: i

    do i = 1, nVal
       val(i) = firstVal + real(i-1)*delVal
    end do
  end subroutine BuildAxis

  function printVars_(self) result(iret)
     class(bramsFile), intent(in) :: self
     type(bramsHeader), pointer :: Head => null()
     integer :: iret
     iret = 0
     Head => self%Header
     do while(associated(Head))
        print*,trim(Head%VarName)
        iret = iret + 1
        Head=>Head%next
    enddo
     
  end function
end module BRAMSMod
