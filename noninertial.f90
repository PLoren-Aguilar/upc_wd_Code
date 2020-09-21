       SUBROUTINE noninertial
!===================================================================
!  This subroutine calculates the effect of non-inertial
!  forces 
!
!  Last revision: 9/October/2015
!===================================================================
!
!--Load modules
!
      USE mod_parameters, ONLY : ndim
      USE mod_commons, ONLY : xyzhm, vxyzut, rotforc, star, Omega0,     &
      nbody, star, nstep, nstep_ini, rotforc, dtmp_min, nout, nstep,    &
      RELFLAG
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local variables
!
      REAL, DIMENSION(ndim) :: cmp1, cmp2
      REAL :: mass, mtot1, mtot2, cmd, Omega, Omegadot
      INTEGER :: k, p
!
!--Calculate angular velocity of the system
!
       mtot1 = 0.0
!$OMP PARALLEL DEFAULT(none) shared(xyzhm,star,nbody) private(p) &
!$OMP reduction(+:mtot1)
!$OMP DO SCHEDULE(runtime)
       DO p = 1,nbody
          IF (star(p) == 1) mtot1 = mtot1 + xyzhm(5,p)
       ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
       cmp1 = 0.0
       DO k=1,ndim
!$OMP PARALLEL DEFAULT(none) shared(k,xyzhm,nbody,star) private(p,mass) &
!$OMP reduction(+:cmp1)
!$OMP DO SCHEDULE(runtime)
          DO p=1,nbody
             IF (star(p) == 1) THEN
                mass = xyzhm(5,p)
                cmp1(k) = cmp1(k) + mass*xyzhm(k,p)
             ENDIF
          ENDDO
!$OMP END DO
!$OMP END PARALLEL
          cmp1(k) = cmp1(k)/mtot1
       ENDDO
!
       mtot2 = 0.0
!$OMP PARALLEL DEFAULT(none) shared(nbody,star,xyzhm) private(p) &
!$OMP reduction(+:mtot2) 
!$OMP DO SCHEDULE(runtime)
       DO p=1,nbody
          IF (star(p) == 2) mtot2 = mtot2 + xyzhm(5,p)
       ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
       cmp2 = 0.0
       DO k=1,ndim
!$OMP PARALLEL DEFAULT(none) shared(nbody,k,xyzhm,star) private(p,mass) &
!$OMP reduction(+:cmp2)
!$OMP DO SCHEDULE(runtime)
          DO p=1,nbody
             IF (star(p) == 2) THEN
                mass    = xyzhm(5,p)
                cmp2(k) = cmp2(k) + mass*xyzhm(k,p)
             ENDIF
          ENDDO
!$OMP END DO
!$OMP END PARALLEL
          cmp2(k) = cmp2(k)/mtot2
       ENDDO
!
!--Centre of mass distance and rotational angular velocity
!
       cmd   = SQRT((cmp1(1)-cmp2(1))**2 + (cmp1(2)-cmp2(2))**2 +       &
             (cmp1(3)-cmp2(3))**2)
       Omega = SQRT((mtot1+mtot2)/cmd**3)
!
!--Change Omega only at the start of the simulation or during a
!  relaxation
!
       IF (RELFLAG.EQV..true.) THEN
         IF (MOD(nstep,nout).EQ.0) THEN
           Omega0 = Omega
         ELSEIF (Omega0 == 0.0) THEN
           Omega0 = Omega
         ENDIF
       ELSE
         IF (Omega0 == 0.0) THEN
           Omega0 = Omega
         ENDIF
       ENDIF
!
!--This should be added in case of having a dynamical rotating frame. I
!  haven't been able to make it work. Needs more testing
!
       !Omegadot = (Omega - Omega0)/dtmp_min
       Omegadot = 0.0 
!
!--Calculate Coriolis and centrifugal forces
!
!$OMP PARALLEL DEFAULT(none) shared(nbody,rotforc,xyzhm,vxyzut,Omega0) &
!$OMP shared(Omegadot) private(p)
!$OMP DO SCHEDULE(runtime)
      DO p=1,nbody
         rotforc(1,p) = xyzhm(1,p)*Omega0**2 + 2.*vxyzut(2,p)*Omega0 +  &
                        xyzhm(2,p)*Omegadot
         rotforc(2,p) = xyzhm(2,p)*Omega0**2 - 2.*vxyzut(1,p)*Omega0 -  &
                        xyzhm(1,p)*Omegadot
         rotforc(3,p) = 0.0
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
      END SUBROUTINE noninertial
