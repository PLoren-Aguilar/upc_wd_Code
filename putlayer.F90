      SUBROUTINE putlayer
!========================================================================
!  This subroutine is the main driver of the code. It evolves the system
!  in time using a given integrator method, and prints diagnostics,
!  whenever necessary
!
!  Last revision: 15/March/2015
!========================================================================
!
!--Load modules
!
      USE mod_parameters
      USE mod_commons
      USE IFPORT
! 
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local variables
!
      REAL, DIMENSION(nel,npart) :: xss_He
      REAL, DIMENSION(5,npart)   :: xyzhm_He, vxyzut_He
      REAL, DIMENSION(npart)     :: rho_He, fh_He, ka1_He, star_He
      REAL :: rmax, dr2, mHe, xc, yc, zc
      INTEGER, DIMENSION(npart) :: partype_He
      INTEGER :: m, p, k, num,  number
      CHARACTER(30) :: filename
!
!--Startout the code reading the necessary data files and parameters
!
      CALL startout
!
!--How many He particles do you want in the spherical layer? 
!
      WRITE(*,*) 'How many He particles do you want?'
      READ(*,*) nHe 
      WRITE(*,*) 'Mass of the He layer (in solar masses)?'
      READ(*,*) mHe 
!
!--Find out the maximum radius of the WD
!
      rmax = 0.0
      DO p=1,nbody
         rmax = MAX(rmax,xyzhm(1,p)**2 + xyzhm(2,p)**2 + xyzhm(3,p)**2)
      ENDDO
!
!--Now put the spherical layer of He particles around the WD
!
      DO p=1,nHe
         dr2 = 1.0d30
         CALL SEED(p*time())
         DO WHILE (dr2.LT.rmax.OR.dr2.GT.1.5*rmax)
!            xc  = 2.*DSQRT(rmax)*RAND()-DSQRT(rmax)
!            yc  = 2.*DSQRT(rmax)*RAND()-DSQRT(rmax)
!            zc  = 2.*DSQRT(rmax)*RAND()-DSQRT(rmax)
            xc = RAND()*SIN(2.0*pi*RAND())*COS(2.0*pi*RAND())
            yc = RAND()*SIN(2.0*pi*RAND())*SIN(2.0*pi*RAND())
            zc = RAND()*COS(2.0*pi*RAND())
            dr2 = xc*xc+yc*yc+zc*zc
         ENDDO
         xyzhm_He(1,p) = xc
         xyzhm_He(2,p) = yc
         xyzhm_He(3,p) = zc
!
         vxyzut_He(1,p) = 0.0
         vxyzut_He(2,p) = 0.0
         vxyzut_He(3,p) = 0.0
         vxyzut_He(4,p) = 1.0d6
!
         ka1_He(p)     = ka_min
         fh_He(p)      = 1.
         partype_He(p) = 1
         star_He(p)    = 1
!
         xyzhm_He(5,p) = mHe/nHe
         rho_He(p)     = 1.0
         xyzhm_He(4,p) = hfact*((xyzhm_He(5,p)/rho_He(p))**(1./3.))
!
         DO k=1,nel
            xss_He(k,p) = 0.0
         ENDDO
         xss_He(2,p) = 1.0
      ENDDO
!
!--Now, write data into file
!
      num = 2
      WRITE(filename,'(a,i4.4,a)') 'bodi',num,'.out'
      OPEN (UNIT=1, FILE=filename, STATUS='new')
      WRITE(filename,'(a,i4.4,a)') 'comp',num,'.out'
      OPEN (UNIT=2, FILE=filename, STATUS='new')
!
!--Write CO data
!
      number = nbody + nHe
      WRITE(1,*) 'ASPLASH_HEADERLINE_TIME'
      WRITE(1,*) tnow
      WRITE(1,*) 'ASPLASH_HEADERLINE_GAMMA'
      WRITE(1,*) '1.667'
      WRITE(1,*) 'NPART'
      WRITE(1,*) number
      WRITE(2,*) tnow
      DO p=1,nbody
         WRITE(1,'(16(1ES12.4),i3.1,i3.1)') (xyzhm(k,p),k=1,5),rho(p),  &
        (vxyzut(k,p),k=1,5),ka1(p),fh(p),tscdyn(p),tscnuc(p),enuc(p),   &
         partype(p), star(p)
!
         WRITE(2,'(16(1ES12.4))') (xss(k,p),k=1,nel)
      ENDDO
!
!--Write He data
!
      DO p=1,nHe
         WRITE(1,'(16(1ES12.4),i3.1,i3.1)') (xyzhm_He(k,p),k=1,5),rho_He(p), &
        (vxyzut_He(k,p),k=1,5),ka1_He(p),fh_He(p),0.0,0.0,0.0,          &
        partype_He(p), star(p)
!
         WRITE(2,'(16(1ES12.4))') (xss_He(k,p),k=1,nel)
      ENDDO
!
      CLOSE (UNIT=1)
      CLOSE (UNIT=2)
!
      END SUBROUTINE putlayer
