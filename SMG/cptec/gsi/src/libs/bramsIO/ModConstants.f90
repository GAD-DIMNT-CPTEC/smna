module ModConstants
   implicit none
   private

   !
   !precisao das variaveis
   integer, public, parameter :: I4 = SELECTED_INT_KIND( 9)  ! Kind for 32-bits Integer Numbers
   integer, public, parameter :: I8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers

   integer, public, parameter :: R4 = SELECTED_REAL_KIND( 6) ! Kind for 32-bits Real Numbers
   integer, public, parameter :: R8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers

   !Logical Units 
   integer,  public, parameter :: stderr = 0 ! Error Unit
   integer,  public, parameter :: stdinp = 5 ! Input Unit
   integer,  public, parameter :: stdout = 6 ! Output Unit

   !Undefine value
   real(kind=r4), public, parameter :: Undef = -1e12

   ! Constants
   real(kind=r8), public, parameter :: cp      = 1004.0_r8     ! Specific heat of dry air at constant pressure (Joules kg^-1 K^-1)
   real(kind=r8), public, parameter :: Rd      = 287.04_r8     ! Constante dos Gases para o ar seco [ Joutes kg^-1 K^-1]
   real(kind=r8), public, parameter :: cpor    = cp/Rd
   real(kind=r8), public, parameter :: p00     = 1.e5_r8       ! Sea level standard atmospheric pressure [Pa] ??

   real(kind=r8), public, parameter :: pi      = 4.0_r8*ATAN(1.0_r8)
   real(kind=r8), public, parameter :: deg2rad = pi / 180.0_r8 ! convert from degree to radians
   real(kind=r8), public, parameter :: rad2deg = 180.0_r8 / pi ! convert from radians to degree

   real(kind=r4), public, parameter :: r60inv  = 1.0/60.0

   real(kind=r8), public, parameter :: earthRadius   = 6.37E6_r8
   real(kind=r8), public, parameter :: earthDiameter = 2*earthRadius
   real(kind=r8), public, parameter :: earthRadius2  = earthRadius*earthRadius
   real(kind=r8), public, parameter :: earthRadius1  = 1.0_r8/earthRadius
   real(kind=r8), public, parameter :: earthRadius12 = earthRadius1*earthRadius1

   !------------------------------------------------!
   ! Word Geodetic System 1984 (WGS84)
   real(kind=r8), public, parameter :: EarthRadiusMajorAxis = 6378137.0000_r8
   real(kind=r8), public, parameter :: EarthRadiusMinorAxis = 6356752.3142_r8
   real(kind=r8), public, parameter :: flattening           = 1.0/298.257223563

   ! Mean Earth Radius in m.  The value below is consistent
   ! with NCEP's routines and grids.
   real (kind=r8), public, parameter :: A_WGS84  = 6378137.0_r8
   real (kind=r8), public, parameter :: B_WGS84  = 6356752.314_r8
   real (kind=r8), public, parameter :: RE_WGS84 = A_WGS84
   real (kind=r8), public, parameter :: E_WGS84  = 0.081819192_r8

   !
   !------------------------------------------------!
   ! other parameters
   integer, public, parameter :: strlen = 512!1024
   
end module ModConstants
