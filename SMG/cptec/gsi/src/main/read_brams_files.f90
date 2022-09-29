!---------------------------------------------------------------------!
!BOP  
! !IROUTINE: read_brams_files 
!
! !DESCRIPTION: same as read_files.f90, but for files used in brams
!  
!               Get info about atm & sfc guess files.
!               figure out available time levels of background fields 
!               for later input to GSI.
!
! !INTERFACE:

subroutine read_brams_files(mype)
   use kinds, only: r_kind, i_kind, i_llong
   use constants, only: r60inv
   use gsi_4dvar, only: l4dvar, iwinbgn, winlen, nhr_assimilation
   use obsmod, only: iadate
   use read_obsmod, only: gsi_inquire
   use mpimod, only: mpi_comm_world, ierror,mpi_rtype, mpi_itype, npe
   use BRAMSMod, only: bramsFile

   use guess_grids, only: nfldsig,      & ! actual count of in-cache time slots for sigma
                          nfldsfc,      & ! actual count of in-cache time slots for sfc
                          nfldnst,      & ! actual count of in-cache time slots for nst FCST fields ???
                          ntguessig,    & ! location of actual guess time for sigma fields
                          ntguessfc,    & ! location of actual guess time for sfc fields
                          ifilesig,     & ! array used to open the correct sigma guess files
                          ifilesfc,     & ! array used to open the correct surface guess files
                          hrdifsig,     & ! times for cached sigma guess_grid
                          hrdifsfc,     & ! times for cached surface guess_grid
                          hrdifsig_all, & ! a list of all times for sigma files
                          hrdifsfc_all, & ! a list of all times for sfc files
                          create_gesfinfo ! Allocate guess-files information arrays

   implicit none

!
! !INPUT PARAMETERS:
!
   integer(i_kind), intent(in) :: mype
!
! !LOCAL VARIABLES:
!
   character(len=64), parameter :: myname_= 'read_brams_files( )'
   integer, parameter           :: MaxFiles = 100
   real(r_kind)       :: t4dv
   integer(i_kind)    :: workPe
   integer(i_kind)    :: nhrHalf
   integer(i_kind)    :: nfld
   integer(i_kind)    :: ihour
   integer(i_kind)    :: inmings ! minutes of initial time (from 1978 see below)
   integer(i_kind)    :: fnmings ! minutes of forecast time (from 1978 see below)
   integer(i_kind)    :: nminanl ! minutes of analysis time (from 1978 see below)
   integer(i_kind)    :: iamana
   integer(i_kind)    :: ndiff
   integer(i_kind)    :: i
   integer(i_kind)    :: iret
   integer(i_kind)    :: npem1
   integer(i_llong)   :: lenbytes
   logical            :: FExist
!   type(BramsFile), pointer    :: BRAMS =>null()
   type(BramsFile)    :: BRAMS
   character(len=256) :: WhatPrint
   character(len= 20) :: FileGuess
   character(len= 20) :: FileHeader

   real(r_kind),    allocatable :: tFiles(:) ! time of files
   integer(i_kind), allocatable :: iFiles(:) ! i order of files
   integer(i_kind),dimension(5) :: itime
   integer(i_kind),dimension(5) :: ftime
   integer(i_kind)              :: fct_time

   type adtime
      real    :: tval
      integer :: ival
      type(adtime), pointer :: next => null()
   end type
   type(adtime), pointer :: timeRoot => null()
   type(adtime), pointer :: time => null()
!
! !REVISION HISTORY:
!   2019-03-19 - J.G.Z. de Mattos - Initial version
!
! !REMARKS:
!   Adapted from read_files.f90, read_wrf_nmm_files.f90
!EOP
!---------------------------------------------------------------------!
!BOC




   nfld    = 0
   workPe  = npe - 1
   nhrHalf = nhr_assimilation / 2
   if (nhrHalf*2 < nhr_assimilation) nhrHalf = nhrHalf + 1

   ! Convert analysis time to minutes relative to fixed date
   call w3fs21(iadate,nminanl)
   write(6,*)trim(myname_),' :  analysis date,minutes ',iadate,nminanl

   
   if (mype == workPe)then

      allocate(timeRoot)
      time => timeRoot
      do ihour = 0, MaxFiles-1
      
         write(FileGuess,'(''BRAMS.vfm.'',i2.2)')ihour
         call gsi_inquire(lenbytes, fexist, FileGuess, mype)

         if(FExist .and. lenBytes > 0)then
            write(FileHeader,'("BRAMS.head.",I2.2)') ihour
            iret = brams%BrOpen(FileHeader, FileGuess)
            if (iret .ne. 0)then
               WhatPrint='ERROR: Some error to open/read BRAMS Head file:'//trim(FileHeader)
               write(6,'(2(1x,A))')trim(myname_),trim(WhatPrint)
               call stop2(99)
            endif

            !
            ! check for consistency of times from brams guess files
            !

            ! Initial contition time

            itime(:) = 0
            itime(1) = brams%GetTimeInfo('iyr')
            itime(2) = brams%GetTimeInfo('imo')
            itime(3) = brams%GetTimeInfo('idy')
            itime(4) = brams%GetTimeInfo('ihr')
            itime(5) = brams%GetTimeInfo('imn')
   
            !
            ! Convert time to minutes relative to fixed date
            ! itime(1) -> Year
            ! itime(2) -> Month
            ! itime(3) -> Day
            ! itime(4) -> Hour
            ! itime(5) -> Minute
            ! inmings   -> Integer number of minutes since 1 jan 1978
   
            call w3fs21(itime,inmings)
            ! Forecast time

            ftime(:) = 0
            ftime(1) = brams%GetTimeInfo('fyr')
            ftime(2) = brams%GetTimeInfo('fmo')
            ftime(3) = brams%GetTimeInfo('fdy')
            ftime(4) = brams%GetTimeInfo('fhr')
            ftime(5) = brams%GetTimeInfo('fmn')
   
            !
            ! Convert time to minutes relative to fixed date
            ! ftime(1) -> Year
            ! ftime(2) -> Month
            ! ftime(3) -> Day
            ! ftime(4) -> Hour
            ! ftime(5) -> Minute
            ! fnmings   -> Integer number of minutes since 1 jan 1978
   
            call w3fs21(ftime,fnmings) 

            ! NOTA:
            ! iwinbgn is the time since ref at start of 4dvar window (hours)
            ! this parameter was defined in gesinfo
            !

            t4dv = real((fnmings - iwinbgn),r_kind)*r60inv ! position on time window [ -3, 0, 3 ]
            if (l4dvar)then
               if (t4dv < 0 .or. t4dv > winlen ) cycle
            else
               ndiff = (fnmings - nminanl)*r60inv
               if(abs(ndiff) > nhrHalf) cycle
            endif          

            time%tval = t4dv
            time%ival = ihour
            nfld      = nfld + 1
            if (fnmings == nminanl) iamana = nfld

            allocate(time%next)
            time => time%next
            
            iret = brams%BrClose( )


         endif
      enddo

      if(nfld > 0) then
         time => timeRoot
         allocate(tFiles(nfld), iFiles(nfld))
!         write(6,'(A5,1x,I3.1,1x,A10)')'Found',nfld,'BRAMS files:'
         do i = 1, nfld
!            write(*,'(A4,1x,I3.1,1x,A10,I2.2)')'File',i,'BRAMS.fct.',time%ival
            tFiles(i) = time%tval
            iFiles(i) = time%ival
            time => time%next
         enddo
      else
         write(6,*)trim(myname_),' ***ERROR*** NO valid BRAMS fields; aborting...'
         call stop2(169)
      end if

      !deallocate time array
      time=>timeRoot%next
      do
         deallocate(timeRoot)
         if(.not.associated(time)) exit
         timeRoot => time
         time => time%next
      enddo

   endif

   ! Broadcast guess file information to all tasks
   call mpi_bcast(  nfld, 1,mpi_itype,workPe,mpi_comm_world,ierror) ! number of files
   call mpi_bcast(iamana, 1,mpi_itype,workPe,mpi_comm_world,ierror)

   if (.not. allocated(tFiles)) allocate(tFiles(nfld))
   if (.not. allocated(iFiles)) allocate(iFiles(nfld))

   call mpi_bcast(tFiles,nfld,mpi_rtype,workPe,mpi_comm_world,ierror) ! time of each file
   call mpi_bcast(iFiles,nfld,mpi_itype,workPe,mpi_comm_world,ierror) ! number of each file

   ! Allocate space for guess information files
   ntguessig = iamana
   ntguessfc = iamana
   nfldsig   = nfld
   nfldsfc   = nfld
   nfldnst   = nfld
   
   call create_gesfinfo


   ! transfer guess information to GSI system
   do i=1,nfld
      hrdifsig(i) = tFiles(i)
      ifilesig(i) = iFiles(i)
      
      hrdifsfc(i) = tFiles(i)
      ifilesfc(i) = iFiles(i)

      hrdifsig_all(i) = tFiles(i)
      hrdifsfc_all(i) = iFiles(i)
   end do

   if(mype == 0)then
      write(6,*)trim(myname_),':  atm/sfc fcst files used in analysis  :  '
      do i=1,nfld
       if(i==iamana)then
          write(*,'(A4,1x,I3.1,1x,A10,I2.2,1x,A)')'File',i,'BRAMS.fct.',ifilesig(i), '<-- Analysis time '
       else
          write(*,'(A4,1x,I3.1,1x,A10,I2.2)')'File',i,'BRAMS.fct.',ifilesig(i)
       endif
      enddo
   endif

end subroutine
