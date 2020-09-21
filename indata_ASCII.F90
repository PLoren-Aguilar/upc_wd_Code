      SUBROUTINE indata
!===================================================================
!
!  This subroutine reads the initial data.
!
!  Last revision: 15/March/2015
!
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
!--Open data unit
!
      WRITE(filename,'(a,i4.4,a)') 'bodi',nstep_ini,'.out'
      OPEN (UNIT=1, FILE=filename)
      WRITE(filename,'(a,i4.4,a)') 'comp',nstep_ini,'.out'
      OPEN (UNIT=2, FILE=filename)
!
!--Read time header
!
      READ(1,*)
      READ(1,*) tnow
      READ(2,*) 
!
!--Read Gamma 
!
      READ(1,*)
      READ(1,*) gamma
!
!--Read number of dust and gas particles
!
      READ(1,*)
      READ(1,*) nums
!
!--Check if the number of particles is consistent
!
      IF (nums /= npart) THEN
         IF (rank == MASTER) PRINT*, 'nparts NE nbody!!!'
         STOP
      ENDIF
!
      eps = 1d30
      partloop : DO p = 1,npart
!
!--SPH particle data
!
         READ(1,*) (xyzhm(k,p),k=1,5),rho(p),(vxyzut(k,p),k=1,5),ka1(p),&
         fh(p),tscdyn(p),tscnuc(p),enuc(p),partype(p),star(p)
!
!--Read old code format
!
!         READ(1,*) (xyzhm(k,p),k=1,5),rho(p),(vxyzut(k,p),k=1,5),ka1(p),&
!         fh(p),tscdyn(p),tscnuc(p),enuc(p)
!         IF (p.LE.nbody1) star(p)=1
!         IF (p.GT.nbody1) star(p)=2
!         partype(p) = 0

!
!         IF (partype(p).EQ.-1) partype(p)=2
!         IF (partype(p).EQ.1) vxyzut(5,p)=1e5
!
         ilist(p)  = p
         plist(p)  = p
         tscdyn(p) = 0.0
         tscnuc(p) = 0.0
         enuc(p)   = 0.0
         eps = MIN(eps,xyzhm(4,p))
!
!--Composition data
!
         READ(2,*) (xss(k,p),k=1,nel)
      END DO partloop
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
      eps  = 1.4d0*2.d0*eps
      eps3 = eps*eps*eps
!
      CLOSE(UNIT=1)
      CLOSE(UNIT=2)
!
!--Sort particle list to get rid of dead particles
!
      CALL sort
#ifdef debug
      IF (rank == MASTER) PRINT*, 'sort called'
#endif
!
      END SUBROUTINE indata
