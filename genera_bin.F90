      SUBROUTINE genera_bin
!=====================================================================
!     PROGRAM TO GENERATE THE INITIAL STATE OF THE PARTICLES THAT WILL
!     BE USED IN THE SPH PROGRAM
!
!     Last revision: 15/March/2015 by PLA
!=====================================================================
!
!--Load modules
!
      USE mod_parameters
      USE mod_commons
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local variables
!
      REAL :: m1, m2, x1, x2, v1, v2, rad1, rad2, omega, vrel
      INTEGER ::  i, p, k, nb1, nb2, AllocateStatus
      CHARACTER(30) :: file1, file2, filename
!
!--Read simulation parameters
!
      OPEN (1,FILE='treepars',STATUS='old')
      READ(1,*) nstep_infile
      READ(1,*) tend
      READ(1,*) nout
      READ(1,*) SIMTYPE
      READ(1,*) RELFLAG
      READ(1,*) trelax
      CLOSE(1) 
!
!--Read total number of particles
!
      PRINT*, 'Total number of particles in the binary simulation?'
      READ*, npart
!
!--Allocate particle data arrays and vectors
!
      ALLOCATE (xyzhm(ndim+2,npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for xyzhm!!!"
         STOP
      ENDIF
!
      ALLOCATE (vxyzut(ndim+2,npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for vxyzut!!!"
         STOP
      ENDIF
!
      ALLOCATE (rho(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for rho!!!"
         STOP
      ENDIF
!
      ALLOCATE (ka1(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ka1!!!"
         STOP
      ENDIF
!
      ALLOCATE (fh(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for fh!!!"
         STOP
      ENDIF
!
      ALLOCATE (tscdyn(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for tscdyn!!!"
         STOP
      ENDIF
!
      ALLOCATE (tscnuc(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ka1!!!"
         STOP
      ENDIF
!
      ALLOCATE (enuc(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for enuc!!!"
         STOP
      ENDIF
!
      ALLOCATE (partype(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for partype!!!"
         STOP
      ENDIF
!
      ALLOCATE (star(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for star!!!"
         STOP
      ENDIF
!
      ALLOCATE (xss(nel,npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for xss!!!"
         STOP
      ENDIF
!
!--Read first star data file
!
      PRINT*,'Name of the file for the first star?'
      READ*,file1
!
!--Read old ASCII format
!
      OPEN (UNIT=1,FILE=file1,STATUS='old')
      READ(1,*)
      READ(1,*) tnow
      READ(1,*)
      READ(1,*) 
      READ(1,*) 
      READ(1,*) nb1
!
      m1=0.0d0
      DO 10 p=1,nb1
         READ(1,*) (xyzhm(k,p),k=1,5),rho(p),(vxyzut(k,p),k=1,5),       &
         ka1(p),fh(p),tscdyn(p),tscnuc(p),enuc(p),partype(p),star(p)
!
!--Read even older code format
!
!         READ(1,*) (xyzhm(k,p),k=1,5),rho(p),(vxyzut(k,p),k=1,5),       &
!         ka1(p),fh(p),tscdyn(p),tscnuc(p),enuc(p)

         fh(p) = 1.0d0
         star(p) = 1
         m1 = m1 + xyzhm(5,p)
10    ENDDO
      CLOSE(1)
!
!--Read body 1 composition
!
!      OPEN(UNIT=2,FILE='cdata1',STATUS='old')
!      DO 20 k=1,nel
!         READ(2,*) xss(k,1)
!         DO 30 p=1,nb1
!            xss(k,p) = xss(k,1)
!30       ENDDO
!20    ENDDO
!      CLOSE(2)
!
!--When setting up a Helium layer calculation, setup
!  abundances by hand
!
      DO 20 p=1,nb1
         DO 30 k=1,nel
            xss(k,p) = 0.0d0
30       ENDDO
         IF (partype(p).EQ.0) THEN
            xss(3,p) = 0.4
            xss(4,p) = 0.6
         ELSEIF (partype(p).EQ.1) THEN
            xss(2,p) = 1.0
         ELSE
            xss(3,p) = 0.4
            xss(4,p) = 0.6
         ENDIF
20    ENDDO
!
!--Read second star data file
!
      PRINT*,'Name of the file for the second star?'
      READ*,file2
!
      OPEN (UNIT=3,FILE=file2,STATUS='old')
      READ(3,*)
      READ(3,*) tnow
      READ(3,*)
      READ(3,*)
      READ(3,*)
      READ(3,*) nb2
!
      m2=0.0d0
      DO 40 p=1+nb1,nb1+nb2
         READ(3,*) (xyzhm(k,p),k=1,5),rho(p),(vxyzut(k,p),k=1,5),       &
         ka1(p),fh(p),tscdyn(p),tscnuc(p),enuc(p),partype(p),star(p)
!
!--Read old code format
!
!         READ(3,*) (xyzhm(k,p),k=1,5),rho(p),(vxyzut(k,p),k=1,5),       &
!         ka1(p),fh(p),tscdyn(p),tscnuc(p),enuc(p)

         fh(p) = 1.0d0
         star(p) = 2
         m2 = m2 + xyzhm(5,p)
40    ENDDO
      CLOSE(3)
!
!--Read body 2 composition
!
!      OPEN(UNIT=4,FILE='cdata2',STATUS='old')
!      DO 50 k=1,nel
!         READ(4,*) xss(k,nb1+1)
!         DO 60 p=nb1+1,nb1+nb2
!            xss(k,p) = xss(k,nb1+1)
!60       ENDDO
!50    ENDDO
!      CLOSE(4)
!
!--When setting up a Helium layer calculation, setup
!  abundances by hand
!
      DO 50 p=nb1+1,nb1+nb2
         DO 60 k=1,nel
            xss(k,p) = 0.0d0
60       ENDDO
         IF (partype(p).EQ.0) THEN
            xss(3,p) = 0.4
            xss(4,p) = 0.6
         ELSEIF (partype(p).EQ.1) THEN
            xss(2,p) = 1.0
         ELSE
            xss(3,p) = 0.4
            xss(4,p) = 0.6
         ENDIF
50    ENDDO
!
!--Setup binary system initial positions and velocities
!
      PRINT*,'Distance between objects?'
      READ*,rad2
!
      x1 = m2*rad2/(m2+m1)
      x2 = rad2-x1
      WRITE(*,'(4(1pe12.4))') x1,x2,m1,m2
!
      DO 70 p=1,nb1
         xyzhm(1,p) = xyzhm(1,p) - x1
70    ENDDO
      DO 80 p=nb1+1,nb1+nb2
         xyzhm(1,p) = xyzhm(1,p) + x2
80    ENDDO
!
!--Non-corotating initial conditions
!
!      omega = DSQRT((m1+m2)/rad2**3)
!      DO 90 p=1,nb1+nb2
!         vxyzut(1,p) = -xyzhm(2,p)*omega
!         vxyzut(2,p) =  xyzhm(1,p)*omega
!         vxyzut(3,p) =  0.0d0
!90    ENDDO
!
!--Initial velocity conditions for a collision 
!
!      PRINT*, 'Relative velocity between the stars? (Km/s)'
!      READ*, vrel
!      vrel = vrel*1.0E5*(unt/unl)
!
!      v1 = vrel
!      v2 = 0.0
!      DO 90 p=1,nb1+nb2
!         vxyzut(2,p) =  0.0d0
!         vxyzut(3,p) =  0.0d0
!         IF (p <= nb1) THEN
!            vxyzut(1,p) = vrel
!         ELSE
!            vxyzut(1,p) = 0.0
!         ENDIF
!90    ENDDO
!
!--Write to disk initial data
!
      nstep  = 0
      nout   = 1
      nbody1 = nb1
      nbody2 = nb2
      Omega0 = 0.0
      CALL outdata
!
      END SUBROUTINE genera_bin
