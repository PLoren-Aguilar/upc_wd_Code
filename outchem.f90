      SUBROUTINE outchem
!===================================================================
!
!  This subroutine writes data into disk.
!
!  Last revision: 15/March/2015
!
!===================================================================
!
!--Load modules
!
      USE mod_commons, ONLY : xyzhm,rho,vxyzut,ka1,fh,tscdyn, tscnuc,   &
      enuc, xss, nb, partype, gxyzu, star, nstep, nout, tnow, npart
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
      CHARACTER(30) :: filename1,filename2,filename3
!
!--Open data unit
!
      num=nstep/nout + 1
      WRITE(filename1,'(a,i4.4,a)') 'He',num,'.out'
      WRITE(filename2,'(a,i4.4,a)') 'C',num,'.out'
      WRITE(filename3,'(a,i4.4,a)') 'O',num,'.out'
      OPEN (UNIT=1, FILE=filename1, STATUS='new')
      OPEN (UNIT=2, FILE=filename2, STATUS='new')
      OPEN (UNIT=3, FILE=filename3, STATUS='new')
!
      WRITE(1,*) 'ASPLASH_HEADERLINE_TIME'
      WRITE(1,*) tnow
      WRITE(1,*) 'ASPLASH_HEADERLINE_GAMMA'
      WRITE(1,*) '1.667'
      WRITE(1,*) 'NPART'
      WRITE(1,*) npart
      DO p=1,npart
         WRITE(1,'(16(1ES12.4),i3.1,i3.1)') (xyzhm(k,p),k=1,5),xss(2,p),  &
        (vxyzut(k,p),k=1,5),ka1(p),fh(p),tscdyn(p),tscnuc(p),enuc(p),   &
         partype(p),star(p)
      ENDDO
!
      WRITE(2,*) 'ASPLASH_HEADERLINE_TIME'
      WRITE(2,*) tnow
      WRITE(2,*) 'ASPLASH_HEADERLINE_GAMMA'
      WRITE(2,*) '1.667'
      WRITE(2,*) 'NPART'
      WRITE(2,*) npart
      DO p=1,npart
         WRITE(2,'(16(1ES12.4),i3.1,i3.1)') (xyzhm(k,p),k=1,5),xss(3,p),  &
        (vxyzut(k,p),k=1,5),ka1(p),fh(p),tscdyn(p),tscnuc(p),enuc(p),   &
         partype(p),star(p)
      ENDDO
!
      WRITE(3,*) 'ASPLASH_HEADERLINE_TIME'
      WRITE(3,*) tnow
      WRITE(3,*) 'ASPLASH_HEADERLINE_GAMMA'
      WRITE(3,*) '1.667'
      WRITE(3,*) 'NPART'
      WRITE(3,*) npart
      DO p=1,npart
         WRITE(3,'(16(1ES12.4),i3.1,i3.1)') (xyzhm(k,p),k=1,5),xss(4,p),  &
        (vxyzut(k,p),k=1,5),ka1(p),fh(p),tscdyn(p),tscnuc(p),enuc(p),   &
         partype(p),star(p)
      ENDDO
!
      CLOSE (UNIT=1)
      CLOSE (UNIT=2)
      CLOSE (UNIT=3)
!
      END SUBROUTINE outchem
