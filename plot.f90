      SUBROUTINE plot
!========================================================================
!  This subroutine dumps data in ascii format in order to plot it
!
!  Last revision: 14/March/2015
!========================================================================
      USE mod_parameters
      USE mod_commons
! 
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local variables
!
      INTEGER :: AllocateStatus, p, k, readnstep
      CHARACTER(30) :: numstep, filename1, filename2
!
!--Declare structures in order to read data
!
      TYPE :: outdata
        REAL(4), DIMENSION(ndim+2) :: ioxyzhm
        REAL(4), DIMENSION(ndim+2) :: iovxyzut
        REAL(4) :: iorho
        REAL(4) :: ioka1
        INTEGER(1) :: iopartype
        INTEGER(1) :: iostar
      END TYPE
      TYPE(outdata), ALLOCATABLE, DIMENSION(:) :: data_array
!
      TYPE :: coutdata
        REAL(4), DIMENSION(nel) :: comp
      END TYPE
      TYPE(coutdata), ALLOCATABLE, DIMENSION(:) :: cdata_array
!
!--Number of datafile to read?
!
      PRINT*, 'Number of datafile to read?'
      READ*, readnstep
!
!--Open data unit
!
      WRITE(numstep,'(i4.4)') readnstep
      filename1 = 'bodi'//trim(numstep)//'.out'
      filename2 = 'comp'//trim(numstep)//'.out'
      !filename2 = 'comp0001.out'
      OPEN (1, FILE=filename1, FORM='unformatted')
      OPEN (2, FILE=filename2, FORM='unformatted')
!
!--Read time and npart data
!
      READ(1) tnow, nstep, npart, nbody, ndead, nbody1, nbody2, Omega0
!
!--Allocate space for the bodi and comp data
!
      ALLOCATE (data_array(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for data_array!!!"
         STOP
      ENDIF
!
      ALLOCATE (cdata_array(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for data_array!!!"
         STOP
      ENDIF
!
!--Read bodi and comp data
!
      READ(1) data_array
      READ(2) cdata_array
!
!--Close files
!
      CLOSE(UNIT=1)
!
!--Open data unit to write
!
      WRITE(numstep,'(i4.4)') readnstep
      filename1 = 'asciibodi'//trim(numstep)//'.out'
      filename2 = 'asciicomp'//trim(numstep)//'.out'
      !filename2 = 'asciicomp0001.out'
      OPEN (1, FILE=filename1)
      OPEN (2, FILE=filename2)
!
!--Write time and npart data
!
      WRITE(1,*) 'ASPLASH_HEADERLINE_TIME'
      WRITE(1,*) tnow
      WRITE(1,*) 'ASPLASH_HEADERLINE_GAMMA'
      WRITE(1,*) '1.667'
      !WRITE(1,*) 'Omega0=', Omega0
      PRINT*, 'Omega0=', Omega0

! Find out number of particles of each type
      nHe = 0
      nCO = 0
      ndead = 0
      DO p = 1,npart
         SELECT CASE (data_array(p)%iopartype)
           CASE(0)
             nCO = nCO + 1
           CASE(1)
             nHe = nHe + 1
           CASE(2)
             ndead = ndead + 1
         END SELECT
      ENDDO 
      WRITE(1,*) 'nCO, nHe, ndead'
      WRITE(1,*) nCO, nHe, ndead

! Change reference frame from co-rotation to intertial
      DO p=1,npart
         data_array(p)%iovxyzut(1) = data_array(p)%iovxyzut(1) - data_array(p)%ioxyzhm(2)*Omega0
         data_array(p)%iovxyzut(2) = data_array(p)%iovxyzut(2) + data_array(p)%ioxyzhm(1)*Omega0
      ENDDO

! Print CO particles first 
      DO p=1,npart
         IF (data_array(p)%iopartype == 0) THEN
           WRITE(1,'(12(1ES12.4),i3.1,i3.1)') (data_array(p)%ioxyzhm(k),k=1,5), &
           data_array(p)%iorho, (data_array(p)%iovxyzut(k),k=1,5),data_array(p)%ioka1,&
           data_array(p)%iopartype,data_array(p)%iostar
         ENDIF
      ENDDO
! Then He particles
      DO p=1,npart
         IF (data_array(p)%iopartype == 1) THEN
           WRITE(1,'(12(1ES12.4),i3.1,i3.1)') (data_array(p)%ioxyzhm(k),k=1,5), &
           data_array(p)%iorho, (data_array(p)%iovxyzut(k),k=1,5),data_array(p)%ioka1,&
           data_array(p)%iopartype,data_array(p)%iostar
         ENDIF
      ENDDO
! And finally Dead particles
      DO p=1,npart
         IF (data_array(p)%iopartype == 2) THEN
           WRITE(1,'(12(1ES12.4),i3.1,i3.1)') (data_array(p)%ioxyzhm(k),k=1,5), &
           data_array(p)%iorho, (data_array(p)%iovxyzut(k),k=1,5),data_array(p)%ioka1,&
           data_array(p)%iopartype,data_array(p)%iostar
         ENDIF
      ENDDO
!
!--Write composition as well
!
      DO p=1,npart
         !IF (data_array(p)%iorho <= 1.) PRINT*, p,data_array(p)%iorho,cdata_array(p)%comp(2) 
         WRITE(2,'(16(1ES12.4))') (cdata_array(p)%comp(k),k=1,nel)
      ENDDO
!
!--Close data units
!
      CLOSE (1)
      CLOSE (2)
!      
      END SUBROUTINE plot
