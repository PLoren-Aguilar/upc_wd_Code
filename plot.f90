      subroutine plot
!========================================================================
!  This subroutine dumps data in ascii format in order to plot it
!
!  Last revision: 14/March/2015
!========================================================================
      use mod_parameters
      use mod_commons
! 
!--Force to declare EVERYTHING
!
      implicit none
!
!--Local variables
!
      integer       :: AllocateStatus, p, k, readnstep
      character(30) :: numstep, filename1, filename2
!
!--Declare structures in order to read data
!
      type :: outdata
        real(4), dimension(ndim+2) :: ioxyzhm
        real(4), dimension(ndim+2) :: iovxyzut
        real(4) :: iorho
        real(4) :: ioka1
        integer :: iostep
        integer(1) :: iopartype
        integer(1) :: iostar
      end type
      type(outdata), allocatable, dimension(:) :: data_array
!
      type :: coutdata
        real(4), dimension(nel) :: comp
      end type
      type(coutdata), allocatable, dimension(:) :: cdata_array
!
!--Number of datafile to read?
!
      print*, 'Number of datafile to read?'
      read*, readnstep
!
!--Open data unit
!
      write(numstep,'(i4.4)') readnstep
      filename1 = 'bodi'//trim(numstep)//'.out'
      filename2 = 'comp'//trim(numstep)//'.out'
      open (1, FILE=filename1, form='unformatted')
      open (2, FILE=filename2, form='unformatted')
!
!--Read time and npart data
!
      read(1) tnow, npart, nbody, ndead, nbody1, nbody2, Omega0
!
!--Allocate space for the bodi and comp data
!
      allocate (data_array(npart), STAT = AllocateStatus)
      if (AllocateStatus /= 0) then
         print*, "Not enough memory for data_array!!!"
         stop
      endif
!
      allocate (cdata_array(npart), STAT = AllocateStatus)
      if (AllocateStatus /= 0) then
         print*, "Not enough memory for data_array!!!"
         stop
      endif
!
!--Read bodi and comp data
!
      read(1) data_array
      read(2) cdata_array
!
!--Close files
!
      close(unit=1)
!
!--Open data unit to write
!
      write(numstep,'(i4.4)') readnstep
      filename1 = 'asciibodi'//trim(numstep)//'.out'
      filename2 = 'asciicomp'//trim(numstep)//'.out'
      open (1, FILE=filename1)
      open (2, FILE=filename2)
!
!--Write time and npart data
!
      write(1,*) 'ASPLASH_HEADERLINE_TIME'
      write(1,*) tnow
      write(1,*) 'ASPLASH_HEADERLINE_GAMMA'
      write(1,*) '1.667'

! Find out number of particles of each type
      nHe = 0
      nCO = 0
      ndead = 0
      do p = 1,npart
         select case (data_array(p)%iopartype)
           case(0)
             nCO = nCO + 1
           case(1)
             nHe = nHe + 1
           case(2)
             ndead = ndead + 1
         end select
      enddo 
      write(1,*) 'nCO, nHe, ndead'
      write(1,*) nCO, nHe, ndead

! Change reference frame from co-rotation to intertial
      !do p=1,npart
      !   data_array(p)%iovxyzut(1) = data_array(p)%iovxyzut(1) - data_array(p)%ioxyzhm(2)*Omega0
      !   data_array(p)%iovxyzut(2) = data_array(p)%iovxyzut(2) + data_array(p)%ioxyzhm(1)*Omega0
      !enddo

! Print CO particles first 
      do p=1,npart
         if (data_array(p)%iopartype == 0) then
           write(1,'(12(1ES12.4),i13.10,i3.1,i3.1)') (data_array(p)%ioxyzhm(k),k=1,5), &
           data_array(p)%iorho, (data_array(p)%iovxyzut(k),k=1,5),data_array(p)%ioka1,&
           data_array(p)%iostep,data_array(p)%iopartype+1,data_array(p)%iostar
         endif
      enddo
! Then He particles
      do p=1,npart
         if (data_array(p)%iopartype == 1) then
           write(1,'(12(1ES12.4),i13.10,i3.1,i3.1)') (data_array(p)%ioxyzhm(k),k=1,5), &
           data_array(p)%iorho, (data_array(p)%iovxyzut(k),k=1,5),data_array(p)%ioka1,&
           data_array(p)%iostep,data_array(p)%iopartype+1,data_array(p)%iostar
         endif
      enddo
! And finally Dead particles
      do p=1,npart
         if (data_array(p)%iopartype == 2) then
           write(1,'(12(1ES12.4),i13.10,i3.1,i3.1)') (data_array(p)%ioxyzhm(k),k=1,5), &
           data_array(p)%iorho, (data_array(p)%iovxyzut(k),k=1,5),data_array(p)%ioka1,&
           data_array(p)%iostep,data_array(p)%iopartype+1,data_array(p)%iostar
         endif
      enddo
!
!--Write composition as well
!
      do p=1,npart
         write(2,'(16(1ES12.4))') (cdata_array(p)%comp(k),k=1,nel)
      enddo
!
!--Close data units
!
      close (1)
      close (2)
!      
      end subroutine plot
