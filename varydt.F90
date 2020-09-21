       SUBROUTINE varydt
!===================================================================
!  This subroutine calculates the new time step
!
!  Last revision: 15/March/2015
!===================================================================
!
!--Load modules
!
      USE mod_commons, ONLY : axyzut, xyzhm, vsigmax, vxyzut, uintprev, &
      rho, partype, nbody, dtmp_min, nstep, nstep_ini, gxyzu, enuc,     &
      tscdyn, tscnuc, cps, dtnuc, press, css, tnow, eosflag, rank
      USE mod_parameters, ONLY : MASTER
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local definitions
!
      REAL :: fa, ft, dt1, dt2, dt3, dt4, dt5, dt6, dtmpmin_old, denom
      REAL, PARAMETER :: dtfact=0.3, sigma=1.0, dtlow = 1.e-5
      INTEGER :: p
!
!--Save old time-step for later use and reset time-step
!
      dtmpmin_old = dtmp_min
!
!--New time-step calculation
!
      dtmp_min = 1.0e30
      dt1 = 1.0e30
      dt2 = 1.0e30
      dt3 = 1.0e30
      dt4 = 1.0e30
      dt5 = 1.0e30
      dt6 = 1.0e30
!$OMP PARALLEL DEFAULT(none) shared(axyzut,xyzhm,vsigmax,vxyzut)     &
!$OMP shared(nbody,uintprev,rho,nstep,partype) private(p,fa) &
!$OMP reduction(MIN:dt1,dt2,dt3)
!$OMP DO SCHEDULE(runtime)
      partloop : DO p=1,nbody
#ifdef Helium
         IF (partype(p) /= 0) CYCLE
#else
         IF (partype(p) == 2) CYCLE
#endif

! Dynamical Courant condition
         fa = SQRT(axyzut(1,p)*axyzut(1,p) + axyzut(2,p)*axyzut(2,p) +  &
            axyzut(3,p)*axyzut(3,p))
         IF (fa /= 0.0) THEN
            dt1 = MIN(SQRT(xyzhm(4,p)/fa),dt1)
         ENDIF

! Viscous signal condition
         IF (vsigmax(p) /= 0.0) THEN
            dt2 = MIN(sigma*xyzhm(4,p)/vsigmax(p),dt2)
         ENDIF

      ENDDO partloop
!$OMP END DO
!$OMP END PARALLEL
      dtmp_min = dtfact*MIN(dt1,dt2)
!
!--Return if first simulation step
!
      IF (nstep == nstep_ini) RETURN
!
!--Thermal Courant condition
!
!$OMP PARALLEL DEFAULT(none) shared(vxyzut,uintprev,rho,nbody,partype)  &
!$OMP shared(axyzut,eosflag) private(p,ft) reduction(MIN:dt3,dt4)
!$OMP DO SCHEDULE(runtime)
      partloop2 : DO p = 1,nbody
#ifdef Helium
         IF (partype(p) /= 0) CYCLE
#else
         IF (partype(p) == 2) CYCLE
#endif

! Thermal energy condition
         ft  = ABS(axyzut(4,p))
         IF (ft /= 0.0) THEN

! Limit the minimum time step given by the energy condition. If the predicted
! time-step becomes too small, the particle is killed
            IF (vxyzut(4,p)/ft >= dtlow) THEN
               dt3 = MIN(dt3,vxyzut(4,p)/ft)
            ELSE
               partype(p) = 2
            ENDIF
         ENDIF

! Temperature condition
         ft  = ABS(axyzut(5,p))
         IF (ft /= 0.0) THEN

! Same thing with temperature
            IF (vxyzut(5,p)/ft >= dtlow) THEN
               dt4 = MIN(dt4,vxyzut(5,p)/ft)
            ELSE
               partype(p) = 2
            ENDIF
         ENDIF
      ENDDO partloop2
!$OMP END DO
!$OMP END PARALLEL
      dt3 = dtfact*dt3
      dt4 = dtfact*dt4
      dtmp_min = MIN(dtmp_min,dt3)
      dtmp_min = MIN(dtmp_min,dt4)
!
!--Typical time-scales calculation
!
!$OMP PARALLEL DEFAULT(none) shared(nbody,tscnuc,cps,vxyzut,enuc,dtnuc) &
!$OMP shared(tscdyn,press,gxyzu,css,rho) private(p,denom) reduction(MIN:dt5,dt6)
!$OMP DO SCHEDULE(runtime)
      partloop3 : DO p=1,nbody

! Nuclear time-scale
         IF ((dtnuc(p) > 0.0) .AND. (enuc(p) > 0.0)) THEN
            tscnuc(p) = cps(p)*vxyzut(5,p)/enuc(p)*dtnuc(p)
         ELSE
            tscnuc(p) = 0.0
         ENDIF
         dt5 = MIN(dt5, tscnuc(p))

! Dynamical time-scale
         denom = DSQRT(gxyzu(1,p)*gxyzu(1,p) + gxyzu(2,p)*gxyzu(2,p)+   &
                 gxyzu(3,p)*gxyzu(3,p))*rho(p)*css(p)
         IF (denom /= 0.0) THEN
             tscdyn(p) = press(p)/denom
         ELSE
             tscdyn(p) = 0.0
         ENDIF
         dt6 = MIN(dt6, tscdyn(p))

      ENDDO partloop3
!$OMP END DO
!$OMP END PARALLEL
!
!--Print time-scales
!
      IF (rank == MASTER) THEN 
        OPEN (1,FILE='timescales.out',STATUS='unknown',access='append')
        WRITE(1,'(i7,7(1pe13.5))') nstep, tnow, dt1, dt2, dt3, dt4, dt5, dt6
        CLOSE (1)
      ENDIF
! 
      END SUBROUTINE varydt
