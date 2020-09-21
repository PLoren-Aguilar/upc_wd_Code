      SUBROUTINE genera_single
!====================================================================
!  PROGRAM TO GENERATE THE INITIAL STATE OF THE PARTICLES THAT WILL
!  BE USED IN THE SPH PROGRAM
!
!  Last revision 15/March/2015 by PLA
!====================================================================
!
!--Load modules
!
      USE mod_parameters
      USE mod_commons
      USE IFPORT
      USE mod_EOS
! 
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
#ifdef MPI
      INCLUDE 'mpif.h'
#endif
!
!--Include EOS definitions
!
!      INCLUDE 'vector_eos.dek'
!
!--Local definitions
!
      DOUBLE PRECISION, DIMENSION(maxit) :: radle, densle 
      REAL :: xc, yc, zc, dh, dh2, rx, ry, rz, rnk, radmax,          &
                 dr2, rad2, masstot, radmin, temp, massp, sum1, sum2,   &
                 rhoc_min, rhoc_max, rhoc, mass_min, mass_max, mass,    &
                 REL, abar, zbar, t1, t2
      REAL, PARAMETER :: tolerance = 0.01
      INTEGER :: p, k, i, j, l, nbod, NNNZ, endit, iseed, AllocateStatus
      CHARACTER(13) :: ONC
      LOGICAL :: done
!
!--Startout MPI, if necessary
!
#ifdef MPI
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
#else
      size = 1
      rank = MASTER
#endif
      nprocs = size
      !IF (size.NE.nprocs) THEN
      !   PRINT*, 'fatal error: nprocs',nprocs,' != ',size,' np'
      !   STOP
      !ENDIF
!
!--Read initial model parameters
!
      OPEN (UNIT=1,FILE='genparams',STATUS='old')
      READ(1,*) nbody
      READ(1,*) masstot
      READ(1,*) rnk
      READ(1,*) temp
      READ(1,*) tnow
      CLOSE(1)
!
!--Check if the number of particles is consistent
!
      npart  = nbody
      ndead  = 0
      nbody1 = npart
      nbody2 = 0

      nbmax = 200
      nbmin = 5
      mmax  = 2*npart + 2
      !IF (nbod.NE.nbody) THEN
      !   PRINT*, 'nbod =! nbody'
      !   STOP
      !ENDIF
!
!--Allocate particle data arrays and vectors
!
      ALLOCATE (xyzhm(ndim+2,2*npart+2), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for xyzhm!!!"
         STOP
      ENDIF
!
      ALLOCATE (xyzhmp(ndim+2,2*npart+2), STAT = AllocateStatus)
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
      ALLOCATE (vxyzutp(ndim+2,npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for vxyzut!!!"
         STOP
      ENDIF
!
     ALLOCATE (axyzut(ndim+2,npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for vxyzut!!!"
         STOP
      ENDIF
!
     ALLOCATE (axyzutp(ndim+2,npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for vxyzut!!!"
         STOP
      ENDIF
!
     ALLOCATE (gxyzu(ndim+1,npart), STAT = AllocateStatus)
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
      ALLOCATE (div(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ka1!!!"
         STOP
      ENDIF
!
      ALLOCATE (divt(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ka1!!!"
         STOP
      ENDIF
!
      ALLOCATE (cur(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ka1!!!"
         STOP
      ENDIF
!
      ALLOCATE (css(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ka1!!!"
         STOP
      ENDIF
!
      ALLOCATE (cvs(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ka1!!!"
         STOP
      ENDIF
!
      ALLOCATE (cps(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ka1!!!"
         STOP
      ENDIF
!
      ALLOCATE (uintprev(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ka1!!!"
         STOP
      ENDIF
!
      ALLOCATE (luminuc(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ka1!!!"
         STOP
      ENDIF
!
      ALLOCATE (luminucp(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ka1!!!"
         STOP
      ENDIF
!
      ALLOCATE (dhdt(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ka1!!!"
         STOP
      ENDIF
!
      ALLOCATE (dhdtp(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ka1!!!"
         STOP
      ENDIF
!
      ALLOCATE (vsigmax(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ka1!!!"
         STOP
      ENDIF
!
      ALLOCATE (enucp(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ka1!!!"
         STOP
      ENDIF
!
      ALLOCATE (dPdT(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ka1!!!"
         STOP
      ENDIF
!
      ALLOCATE (press(npart), STAT = AllocateStatus)
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
     ALLOCATE (ilist(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ilist!!!"
         STOP
      ENDIF
!
      ALLOCATE (plist(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for plist!!!"
         STOP
      ENDIF
!
      ALLOCATE (xss(nel,npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for xss!!!"
         STOP
      ENDIF
!
      ALLOCATE (nb(nbmax,npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for nb!!!"
         STOP
      ENDIF
!
!--Read and setup WD composition
!
      OPEN(UNIT=2,FILE='cdata1',STATUS='old')
      DO 10 k=1,nel
          READ(2,*) xss(k,1)
          DO p=2,nbody
             xss(k,p) = xss(k,1)
20        ENDDO
10    ENDDO
      CLOSE(2)
!
!--Use a Bisection strategy in order to solve Lane-Embden 
!  equation for the given mass
!
      WRITE(*,*) ' '
      WRITE(*,*) '************************************************** '
      WRITE(*,*) 'Using a bisection scheme to solve the lane-embden'
      WRITE(*,*) 'equation and generate a first density vs radius'
      WRITE(*,*) 'profile of the WD...'
      WRITE(*,*) '************************************************** '
      WRITE(*,*) ' '
      done = .FALSE.
      rhoc_min = 1.0d4
      rhoc_max = 1.0d8
      DO 30 WHILE (done.EQV..FALSE.)
!
         CALL lanembden(maxit,rhoc_min,endit,mass,radle,densle)
         mass_min = mass
!
         CALL lanembden(maxit,rhoc_max,endit,mass,radle,densle)
         mass_max = mass
!
         rhoc = 0.5*(rhoc_min + rhoc_max)
         CALL lanembden(maxit,rhoc,endit,mass,radle,densle)
!
         IF (masstot.GT.mass) THEN
            rhoc_min = rhoc 
         ELSEIF (masstot.LT.mass) THEN
            rhoc_max = rhoc
         ENDIF 
!               
         REL = ABS(masstot-mass)/masstot
         IF (REL.LE.tolerance) THEN
            done = .TRUE.
         ENDIF
         WRITE(*,35) 'M_WD=',masstot,' M_Bisec=',mass
30    ENDDO
35    FORMAT (A5,1(1F7.4),A8,1F7.4)
!
!--Call lanembden subroutine a last time, saving this time the
!  values for density, radius, pressure, etc
      CALL lanembden(maxit,rhoc,endit,mass,radle,densle)
!
!--Generate random particle distribution
!
      radmax = 10.*radle(endit)
      rad2   = radmax*radmax
      massp  = mass/FLOAT(npart)
!$OMP PARALLEL DEFAULT(none) shared(xyzhm,vxyzut,radle,densle,rho,npart)   &
!$OMP shared(massp,rad2,radmax,temp,endit,ka1,fh) private(p,dr2,xc,yc,zc,i) &
!$OMP private(done,iseed)
!$OMP DO SCHEDULE(runtime)
      DO 40 p=1,npart
         dr2 = 1.0d30
         CALL SEED(p*time())         
         DO 50 WHILE (dr2.GT.rad2)
            xc  = 2.*radmax*RAND()-radmax
            yc  = 2.*radmax*RAND()-radmax
            zc  = 2.*radmax*RAND()-radmax
            dr2 = xc*xc+yc*yc+zc*zc
50       ENDDO
         xyzhm(1,p)  = xc
         xyzhm(2,p)  = yc
         xyzhm(3,p)  = zc
         vxyzut(1,p) = 0.
         vxyzut(2,p) = 0.
         vxyzut(3,p) = 0.
         ka1(p) = ka_min
         fh(p) = 1.
!
!--Compute density using Lane-Embden equation solution
!
         i = 0
         done = .FALSE.
         DO 60 WHILE (done.EQV..FALSE.)
            i = i + 1
!
            IF ((10.*radle(i).LE.SQRT(dr2)).AND.                       &
               (SQRT(dr2).LE.10.*radle(i+1))) THEN
               xyzhm(5,p)  = massp
               rho(p)      = densle(i)/uden
               xyzhm(4,p)  = hfact*((xyzhm(5,p)/rho(p))**(1./3.))
               vxyzut(5,p) = temp
               done = .TRUE.
            ENDIF
!
            IF (i.EQ.endit) THEN
               PRINT*, 'Particle',p,' not found'
               STOP
            ENDIF
60       ENDDO
40    ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
!--Sort particles. Since now factor is a dynamical quantity we need to
!  call sort in order to define it !!
!
      CALL sort
!
!--Now, use the tree to self-consistently calculate the rho-h relation
!
      PRINT*,'Trying to compute h. May take a long time'
      nstep = 1
      nout  = 1
      CALL iter_rhoh
      PRINT*, 'h calculation finished'
!
!--Read EOS tables
!
      CALL read_helm_table
!
      OPEN(UNIT=4,FILE='ZHELI.10b',status='old')
      aion(1) = 1.0d0
      zion(1) = 1.0d0
      DO 90 k=1,nel-2
         READ(4,8050) NNNZ,ONC,aion(k+1),zion(k+1)
90    ENDDO
      CLOSE(UNIT=4)
8050  FORMAT(I4,1X,A5,1X,F4.0,1X,F4.0)
!
!--Calculate thermal energy
!
      CALL EOS0
      PRINT*,'EOS finished'
!
!--Write data to disk
!
      nstep  = 0
      nout   = 1
      Omega0 = 0.0
      CALL outdata
!
      END SUBROUTINE genera_single
!
      SUBROUTINE lanembden(maxit,rhoc,endit,masa,r,rho)
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--I/O variables
!
      DOUBLE PRECISION, DIMENSION(maxit), INTENT(OUT) :: rho, r
      DOUBLE PRECISION, INTENT(IN) :: rhoc
      INTEGER, INTENT(IN)  :: maxit
      INTEGER, INTENT(OUT) :: endit
!
!--Local variables
!
      DOUBLE PRECISION, DIMENSION(2) :: y, yout
      DOUBLE PRECISION :: x, h, t, alpha, gama, K, pres, masa, pi,      &
                          alpha0
      INTEGER :: oup, flag, i, j
!
!--Physical parameters in cgs
!
      REAL, PARAMETER :: plank=6.6261d-27, me=9.1094d-28,            &
                            mp=1.6726d-24, G=6.6726d-8, msol=1.9891d33, &
                            rsol=6.955d10, mue=2.0d0, npol=1.5d0
!
!--External subroutines
!
      EXTERNAL derivs
!
!--Initial variables
!
      gama = (npol+1.)/npol
      pi   = 4.0*ATAN(1.0)
      K    = ((3.)**(2./3.)*(plank**2))/(20.*(pi**(2./3.))*             &
      me*((mp*mue)**(gama)))
      alpha0 = (((npol+1.)*K)/(4.*pi*G))**(1./2.)
      alpha  = alpha0*rhoc**((1.-npol)/(2.*npol))
      y(1)   = 1.0
      y(2)   = 0.0
!
!--Start up the intergration
!
      t = 0.0
      h = 1.0e-6 
      r(1)   = alpha*t/rsol
      rho(1) = rhoc*(y(1)**npol)
      pres   = K*(rho(1)**gama) 
      masa   = -4*pi*(alpha0**3.)*(rhoc**((3.-npol)/(2*npol)))*         &
      (t**2.)*(y(2))
      masa   = masa/msol
!
!--Integrate lane-embden equation
!
      j = 0
      i = 0
      DO WHILE (y(1).GT.0)
         j = j + 1
! 
         CALL rk4(npol,y,2,t,h,yout,derivs)
!
         t = t + h
         y = yout
         IF (mod(j,10000).EQ.0) THEN
           i = i + 1
           r(i)   = alpha*t/rsol
           rho(i) = rhoc*(y(1)**npol)
           pres   = K*(rho(i)**gama)               
           masa   = -4*pi*(alpha0**3.)*(rhoc**((3.-npol)/(2*npol)))*    &
      (t**2.)*(y(2))
           masa   = masa/msol
         ENDIF
      ENDDO
      endit = i
!
      END SUBROUTINE lanembden
!
      SUBROUTINE derivs(npol,x,y,dydx)
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--I/O variables
!
      DOUBLE PRECISION, DIMENSION(*) :: y, dydx
      DOUBLE PRECISION :: npol, x
!
      dydx(1) = y(2)
      IF (y(1) >= 0) THEN
        IF (dabs(x)<1.0e-9) THEN
          dydx(2) = -(y(1)**npol)/3.0
        ELSE
          dydx(2) = -y(1)**npol-2.0*y(2)/x
        ENDIF
      ELSE
        dydx(2) = 0.0
      ENDIF
!
      END SUBROUTINE derivs
!
      SUBROUTINE rk4(npol,y,n,x,h,yout,derivs)
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--I/O variables
!
      REAL, DIMENSION(n), INTENT(IN) :: y
      REAL, INTENT(IN) :: x, h, npol
      REAL, DIMENSION(n), INTENT(OUT) :: yout
!
!--Local variables
!
      REAL, DIMENSION(size(y)) :: k1, k2, k3, k4, yt
      REAL :: h6, hh, xh
      INTEGER :: n
!
      hh=h*0.5
      h6=h/6.0
!
      call derivs(npol,x,y,k1)
      xh=x+hh
      yt=y+hh*k1
!
      call derivs(npol,xh,yt,k2)
      yt=y+hh*k2
!
      call derivs(npol,xh,yt,k3)
      yt=y+h*k3
!
      call derivs(npol,x+h,yt,k4)
      yout=y+h6*(k1+k4+2.0*(k2+k3))
!
      END SUBROUTINE rk4
