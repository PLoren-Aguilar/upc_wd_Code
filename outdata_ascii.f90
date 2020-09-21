      SUBROUTINE outdata_ascii
!===================================================================
!  This subroutine writes data into disk.
!
!  Last revision: 15/March/2015
!===================================================================
!
!--Load modules
!
      USE mod_parameters
      USE mod_commons
!
!--Force do declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local variables
!
      REAL, DIMENSION(npart) :: temporal
      INTEGER, DIMENSION(npart) :: itemporal
      INTEGER :: p, k, num
      CHARACTER(30) :: filename1, filename2, numstep
!
!--Open data unit
!
      num = nstep/nout + 1
      WRITE(numstep,'(i4.4)') num
      filename1 = 'asciibodi'//trim(numstep)//'.out'
      filename2 = 'asciicomp'//trim(numstep)//'.out'
      OPEN (1, FILE=filename1, FORM='unformatted')
      OPEN (2, FILE=filename2, FORM='unformatted')
!
      WRITE(1,*) 'ASPLASH_HEADERLINE_TIME'
      WRITE(1,*) tnow
      WRITE(1,*) 'ASPLASH_HEADERLINE_GAMMA'
      WRITE(1,*) '1.667'
      WRITE(1,*) 'NPART'
      WRITE(1,*) npart
      WRITE(2,*) tnow
      DO p=1,npart
         WRITE(1,'(16(1ES12.4),i3.1,i3.1)') (xyzhm(k,p),k=1,5),rho(p),  &
        (vxyzut(k,p),k=1,5),ka1(p),fh(p),tscdyn(p),tscnuc(p),enuc(p),   &
         partype(p),star(p)
         tscdyn(p) = 0.0
         tscnuc(p) = 0.0
         enuc(p)   = 0.0
!
         WRITE(2,'(16(1ES12.4))') (xss(k,p),k=1,nel)
      ENDDO
!
      CLOSE (UNIT=1)
      CLOSE (UNIT=2)
!
      END SUBROUTINE outdata_ascii
