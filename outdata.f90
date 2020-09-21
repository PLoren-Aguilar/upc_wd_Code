      SUBROUTINE outdata
!===================================================================
!  This subroutine writes data into disk.
!
!  Last revision: 15/March/2015
!===================================================================
!
!--Load modules
!
      USE mod_parameters, ONLY : ndim, nel
      USE mod_commons, ONLY : xyzhm, rho, vxyzut, ka1, fh,xss, nb,      &
      partype, gxyzu, star, nstep, nout, npart, tnow, RELFLAG, nbody,   &
      nbody1, nbody2, ndead, tend, SIMTYPE, trelax, Omega0, mode
!
!--Force do declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local variables
!
      INTEGER :: p, k, num
      CHARACTER(30) :: filename
!
!--Declare structures in order to dump data
!
      TYPE :: iodata
        REAL(4), DIMENSION(ndim+2) :: ioxyzhm
        REAL(4), DIMENSION(ndim+2) :: iovxyzut
        REAL(4) :: iorho
        REAL(4) :: ioka1
        INTEGER(1) :: iopartype
        INTEGER(1) :: iostar
      END TYPE iodata
      TYPE(iodata), DIMENSION(npart) :: data_array
!
      TYPE :: ciodata
        REAL(4), DIMENSION(nel) :: comp
      END TYPE ciodata
      TYPE(ciodata), DIMENSION(npart) :: cdata_array
!
!--Load data into data arrays
!
      DO p = 1,npart
         data_array(p)%ioxyzhm(1:ndim+2)  = xyzhm(1:ndim+2,p)
         data_array(p)%iovxyzut(1:ndim+2) = vxyzut(1:ndim+2,p)
         data_array(p)%iorho              = rho(p)
         data_array(p)%ioka1              = ka1(p)
         data_array(p)%iopartype          = partype(p)
         data_array(p)%iostar             = star(p)
         cdata_array(p)%comp(1:nel)       = xss(1:nel,p)
      ENDDO
!
!--Open bodi unit and write data
!
      num = nstep/nout + 1
      WRITE(filename,'(a,i4.4,a)') 'bodi',num,'.out'
      OPEN (1, FILE=filename, FORM='unformatted', STATUS='new')
      WRITE(1) tnow, nstep, npart, nbody, ndead, nbody1, nbody2, Omega0
      WRITE(1) data_array
      CLOSE(1)
!
!--Open comp unit and write data
!
      IF (mode(1) == 'binary') THEN
         WRITE(filename,'(a,i4.4,a)') 'comp',num,'.out'
         OPEN (2, FILE=filename, FORM='unformatted', STATUS='new')
         WRITE(2) cdata_array
         CLOSE(2)
      ELSEIF (mode(1) == 'run') THEN
         IF ((RELFLAG.EQV..false.) .AND. (MOD(nstep,10*nout) == 0)) THEN
           WRITE(filename,'(a,i4.4,a)') 'comp',num,'.out'
           OPEN (2, FILE=filename, FORM='unformatted', STATUS='new')
           WRITE(2) cdata_array
           CLOSE(2)
        ENDIF
      ENDIF
!
!--Open treepars file and write data
!
      OPEN (3, FILE='treepars', STATUS='unknown')
      WRITE(3,'(i4.4,a)') num,' :: nstep'
      WRITE(3,'(1(1ES12.4),a)') tend,' :: final time in first run'
      WRITE(3,'(i4.4,a)') nout,' :: steps between outputs'
      WRITE(3,'(i4.4,a)') SIMTYPE,' :: SIMTYPE (1:Single star, 2: Binary system)'
      WRITE(3,'(L1,a)') RELFLAG,' :: RELFLAG'
      WRITE(3,'(f8.4,a)') trelax,' :: TRELAX'
      CLOSE (3)
! 
      END SUBROUTINE outdata
