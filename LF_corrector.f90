      SUBROUTINE corrector
!============================================================
!
!  This subroutine corrects the the integrated quantities 
!
!  Last revision: 15/March/2015
!============================================================
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
      REAL :: fa, dt1, dt2, dt3, denom, dtmpmin_old
      REAL, PARAMETER :: dtfact_b=0.3d0, dtfact_h=0.3d0, sigma=1.0d0,&
                            f1 = 1./256., f2 = 255./256.
      INTEGER, DIMENSION(npart) :: itemporal
      INTEGER :: p, k
!
!--Save dtmp_min for later use
!
      dtmpmin_old = dtmp_min  
!
!--Correct positions
!
      DO 10 k=1,ndim
!$OMP PARALLEL DEFAULT(none) shared(nbody,xyzhm,xyzhmp,vxyzut,vxyzutp)  &
!$OMP shared(axyzut,axyzutp,dtmp_min,k) private(p)
!$OMP DO SCHEDULE(runtime)
         DO 20 p=1,nbody
            xyzhm(k,p) = xyzhm(k,p) + 0.1667d0*(axyzut(k,p)-            &
                         axyzutp(k,p))*dtmp_min*dtmp_min
20       ENDDO
!$OMP END DO
!$OMP END PARALLEL
 10   ENDDO
!
!--Correct velocities
!
      DO 30 k=1,ndim
!$OMP PARALLEL DEFAULT(none) shared(vxyzut,vxyzutp,axyzut,axyzutp,dtmp_min,k)   &
!$OMP shared(nbody) private(p)
!$OMP DO SCHEDULE(runtime)
         DO 40 p=1,nbody
            vxyzut(k,p) = vxyzut(k,p) + 0.5*(axyzut(k,p)-               &
                          axyzutp(k,p))*dtmp_min
40       ENDDO
!$OMP END DO
!$OMP END PARALLEL
30    ENDDO
!
!--Correct thermal energies and temperatures if necessary
!
!$OMP PARALLEL DEFAULT(none) shared(vxyzut,vxyzutp,axyzut,axyzutp,dtmp_min)     &
!$OMP shared(nbody,enuc,enucp,luminuc,luminucp,RELFLAG) private(p)
!$OMP DO SCHEDULE(runtime)
      DO 50 p=1,nbody
         vxyzut(4,p) = vxyzut(4,p) + 0.5*(axyzut(4,p)-                  &
                       axyzutp(4,p))*dtmp_min + enuc(p) - enucp(p)
50    ENDDO
!$OMP END DO
      IF (RELFLAG.EQV..false.) THEN
!$OMP DO SCHEDULE(runtime)
         DO 60 p=1,nbody
            vxyzut(5,p) = vxyzut(5,p) + 0.5*(axyzut(5,p)-axyzutp(5,p))* &
                          dtmp_min + luminuc(p) - luminucp(p)
60       ENDDO
!$OMP END DO
      ENDIF
!$OMP END PARALLEL
!!
!!--Correct smoothing lenghts
!!
!!$OMP PARALLEL DEFAULT(none) shared(nbody,xyzhm,xyzhmp,dhdt,dhdtp,dtmp_min)  &
!!$OMP private(p)
!!$OMP DO SCHEDULE(runtime)
!      DO 70 p=1,nbody
!         xyzhm(4,p) = xyzhm(4,p) + 0.5*(dhdt(p) - dhdtp(p))*dtmp_min
!70    ENDDO
!!$OMP END DO
!!$OMP END PARALLEL
!
!--Save quantities for later use
!
      DO 80 k=1,5 
!$OMP PARALLEL DEFAULT(none) shared(axyzut,axyzutp,vxyzut,vxyzutp,xyzhm,xyzhmp) &
!$OMP shared(k,nbody) private(p)
!$OMP DO SCHEDULE(runtime)
         DO 90 p=1,nbody
            xyzhmp(k,p)  = xyzhm(k,p)
            vxyzutp(k,p) = vxyzut(k,p)
            axyzutp(k,p) = axyzut(k,p)
90       ENDDO
!$OMP END DO
!$OMP END PARALLEL
80    ENDDO
!
!$OMP PARALLEL DEFAULT(none) shared(enucp,enuc,luminucp,luminuc,dhdtp)  &
!$OMP shared(dhdt,nbody) private(p)
!$OMP DO SCHEDULE(runtime)
      DO 100 p=1,nbody
         enucp(p)    = enuc(p)
         luminucp(p) = luminuc(p)
         dhdtp(p)    = dhdt(p)
100   ENDDO
!$OMP END DO
!$OMP END PARALLEL
      RETURN
!
!--Calculate next time-step
!
      dt1 = 1.0d30
      dt2 = 1.0d30
!$OMP PARALLEL DEFAULT(none) shared(nbody,axyzut,xyzhm,vsigmax) &
!$OMP private(p,fa) reduction(MIN:dt1) reduction(MIN:dt2)
!$OMP DO SCHEDULE(runtime)
      DO 110 p=1,nbody
         fa = DSQRT(axyzut(1,p)*axyzut(1,p)+axyzut(2,p)*axyzut(2,p) +   &
              axyzut(3,p)*axyzut(3,p))
         dt1 = MIN(dt1,DSQRT(xyzhm(4,p)/fa))
         IF (vsigmax(p).NE.0.0d0) THEN
            dt2 = MIN(dt2,sigma*xyzhm(4,p)/vsigmax(p))
         ENDIF
110   ENDDO
!$OMP END DO
!$OMP END PARALLEL
      dtmp_min = dtfact_h*MIN(dt1,dt2)
!
      dt3 = 1.0d30
      IF (nstep.GT.nstep_ini+1) THEN
!$OMP PARALLEL DEFAULT(none) shared(nbody,vxyzut,uintprev,rho)  &
!$OMP private(p) reduction(MIN:dt3)
!$OMP DO SCHEDULE(runtime)
         DO 120 p=1,nbody
            IF (DABS(vxyzut(4,p)-uintprev(p)).GT.0.0d0.AND.rho(p).GT.   &
                5.0d0) THEN
                dt3 = MIN(DABS(uintprev(p)/(vxyzut(4,p)-uintprev(p))),  &
                      dt3)
            ENDIF
120      ENDDO
!$OMP END DO
!$OMP END PARALLEL
         dt3 = dtfact_b*dtmpmin_old*dt3
      ENDIF
      dtmp_min = MIN(dtmp_min,dt3)
!
!--Typical time-scales calculation
!
!$OMP PARALLEL DEFAULT(none) shared(nbody,tscnuc,cps,vxyzut,enuc,dtnuc) &
!$OMP shared(tscdyn,press,gxyzu,css,rho) private(p,denom)
!$OMP DO SCHEDULE(runtime)
      DO 130 p=1,nbody
         IF ((dtnuc(p).GT.0.0d0).AND.(enuc(p).GT.0.0d0)) THEN
            tscnuc(p) = cps(p)*vxyzut(5,p)/enuc(p)*dtnuc(p)
         ELSE
            tscnuc(p) = 0.0d0
         ENDIF
!
         denom = DSQRT(gxyzu(1,p)*gxyzu(1,p) + gxyzu(2,p)*gxyzu(2,p)+   &
                 gxyzu(3,p)*gxyzu(3,p))*rho(p)*css(p)
         IF (denom.NE.0.0d0) THEN
             tscdyn(p) = press(p)/denom
         ELSE
             tscdyn(p) = 0.0d0
         ENDIF
130   ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
      END SUBROUTINE corrector
