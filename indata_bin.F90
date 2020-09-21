      SUBROUTINE indata
!===================================================================
!  This subroutine reads the initial data.
!
!  Last revision: 15/March/2015
!===================================================================
!
!--Load modules
!
      USE mod_essentials
      USE mod_commons
!     
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local variables
!
      REAL, DIMENSION(npart) :: temporal 
      REAL :: real_factor,cmx,cmy,cmz
      INTEGER, DIMENSION(npart) :: itemporal 
      INTEGER ::  p, k, nums, m
      CHARACTER(30) :: filename
!
!--Declare structures in order to read data
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
!--Open data unit
!
      WRITE(filename,'(a,i4.4,a)') 'bodi',nstep_ini,'.out'
      OPEN (1, FILE=filename, FORM='unformatted')
      WRITE(filename,'(a,i4.4,a)') 'comp',nstep_ini,'.out'
      OPEN (2, FILE=filename, FORM='unformatted')
!
!--Read bodi data
!
      READ(1) tnow, nums, data_array
!
!--Check if the number of particles is consistent
!
      IF (nums /= npart) THEN
         IF (rank == MASTER) PRINT*, 'nparts NE nbody!!!'
         STOP
      ENDIF
!
!--Composition data
!
      READ(2) cdata_array
!
!--Close files
!
      CLOSE(UNIT=1)
      CLOSE(UNIT=2)
!
!--Transfer data into the appropriate locations
!
      DO p = 1,npart
         xyzhm(1:ndim+2,p)  = data_array(p)%ioxyzhm(1:ndim+2)
         vxyzut(1:ndim+2,p) = data_array(p)%iovxyzut(1:ndim+2)
         rho(p)             = data_array(p)%iorho
         ka1(p)             = data_array(p)%ioka1
         partype(p)         = data_array(p)%iopartype
         star(p)            = data_array(p)%iopartype
         xss(1:nel,p)       = cdata_array(p)%comp(1:nel)
      ENDDO
!
!--Start up variables
!
      vsigmax = 1.0d0
      cur     = 0.0d0
      div     = 0.0d0
      luminuc = 0.0d0
      enuc    = 0.0d0
      tscdyn  = 0.0d0
      tscnuc  = 0.0d0
!
      eps = 1d30
      DO p = 1, npart
         eps = MIN(eps,xyzhm(4,p))
         ilist(p)  = p
         plist(p)  = p
      ENDDO
      eps  = 1.4d0*2.d0*eps
      eps3 = eps*eps*eps
!
!--Sort particle list to get rid of dead particles
!
      CALL sort
#ifdef debug
      IF (rank == MASTER) PRINT*, 'sort called'
#endif
!
      END SUBROUTINE indata
