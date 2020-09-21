      SUBROUTINE outdata
!===================================================================
!  This subroutine writes data into disk.
!
!  Last revision: 15/March/2015
!===================================================================
!
!--Load modules
!
      USE mod_essentials
      USE mod_commons, ONLY : xyzhm,rho,vxyzut,ka1,fh,tscdyn, tscnuc,   &
                              enuc, xss, nb, partype, gxyzu, star
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
      CHARACTER(30) :: filename
!
!--Declare structures in order to dump data
!
      TYPE :: outdata
        REAL, DIMENSION(ndim+2) :: ioxyzhm
        REAL, DIMENSION(ndim+2) :: iovxyzut
        REAL :: iorho
        REAL :: ioka1
        INTEGER(1) :: iopartype
        INTEGER(1) :: iostar
      END TYPE
      TYPE(outdata), DIMENSION(npart) :: data_array
!
      TYPE :: coutdata
        REAL, DIMENSION(nel) :: comp
      END TYPE
      TYPE(coutdata), DIMENSION(npart) :: cdata_array
!
!--Load data into data arrays
!
      DO p = 1,npart
         data_array(p)%ioxyzhm(1:ndim+2)  = xyzhm(1:ndim+2,p)
         data_array(p)%iovxyzut(1:ndim+2) = vxyzut(1:ndim+2,p)
         data_array(p)%iorho              = rho(p)
         data_array(p)%ioka1              = ka1(p)
         data_array(p)%iopartype          = partype(p)
         data_array(p)%iopartype          = star(p)
         cdata_array(p)%comp(1:nel)       = xss(1:nel,p)
      ENDDO
!
!--Open bodi unit and write data
!
      num = nstep/nout + 1
      WRITE(filename,'(a,i4.4,a)') 'bodi',num,'.out'
      OPEN (1, FILE=filename, FORM='unformatted', STATUS='new')
      WRITE(1) tnow, npart, data_array
      CLOSE(1)
!
!--Open comp unit and write data
!
      IF ((RELFLAG.EQV..false.) .AND. (MOD(10*nstep,nout) == 0)) THEN
         WRITE(filename,'(a,i4.4,a)') 'comp',num,'.out'
         OPEN (2, FILE=filename, FORM='unformatted', STATUS='new')
         WRITE(2) cdata_array
         CLOSE(2)
      ENDIF
!
      tscdyn(1:nbody) = 0.0
      tscnuc(1:nbody) = 0.0
      enuc(1:nbody)   = 0.0
!
      END SUBROUTINE outdata
