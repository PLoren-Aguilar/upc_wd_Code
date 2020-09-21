      SUBROUTINE predictor
!============================================================
!
!  This subroutine predicts in time the integrated quantities 
!
!  Last revision: 15/March/2015
!
!============================================================
!
!--Load modules
!
      USE mod_parameters, ONLY : ndim
      USE mod_parameters, ONLY : tmin, tmax
      USE mod_commons, ONLY : uintprev, vxyzut, axyzutp, dtmp_min,      &
                              enucp, luminucp, xyzhm, dhdtp, eps, eps3, &
                              xyzhmp, vxyzutp, nbody, RELFLAG
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local variables
!
      INTEGER :: p, k
!
!--Save thermal energy for later use
!
!$OMP PARALLEL DEFAULT(none) shared(nbody,uintprev,vxyzutp) private(p)
!$OMP DO SCHEDULE(runtime)
      DO p=1,nbody
         uintprev(p) = vxyzutp(4,p)
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
!--Evolve positions
!
      DO k=1,ndim
!$OMP PARALLEL DEFAULT(none) shared(nbody,xyzhm,vxyzut,vxyzutp,dtmp_min) &
!$OMP shared(xyzhmp,axyzutp,k) private(p)
!$OMP DO SCHEDULE(runtime)
         DO p=1,nbody
#ifdef Helium
            xyzhm(k,p) = xyzhmp(k,p) + 0.5*vxyzutp(k,p)*dtmp_min
#else
            xyzhm(k,p) = xyzhm(k,p) + vxyzut(k,p)*dtmp_min +            &
                         0.5d0*axyzutp(k,p)*dtmp_min*dtmp_min
#endif
         ENDDO
!$OMP END DO
!$OMP END PARALLEL
      ENDDO 
!
!--Evolve velocities
!
      DO k=1,ndim
!$OMP PARALLEL DEFAULT(none) shared(nbody,vxyzut,axyzutp,dtmp_min,k) &
!$OMP shared(vxyzutp) private(p)
!$OMP DO SCHEDULE(runtime)
         DO p=1,nbody
#ifdef Helium
            vxyzut(k,p) = vxyzutp(k,p) + 0.5*axyzutp(k,p)*dtmp_min
#else
            vxyzut(k,p) = vxyzut(k,p) + axyzutp(k,p)*dtmp_min
#endif
         ENDDO
!$OMP END DO
!$OMP END PARALLEL
      ENDDO
!
!--Evolve thermal energies and temperatures, if necessary
!
!$OMP PARALLEL DEFAULT(none) shared(nbody,vxyzut,axyzutp,dtmp_min) &
!$OMP shared(enucp,luminucp,RELFLAG,vxyzutp) private(p)
!$OMP DO SCHEDULE(runtime)
      DO p=1,nbody
#ifdef Helium
         vxyzut(4,p) = vxyzutp(4,p) + 0.5*axyzutp(4,p)*dtmp_min +       &
                       enucp(p)
#else
         IF ((vxyzut(5,p) > tmin).AND.(vxyzut(5,p) < tmax)) THEN
            vxyzut(4,p) = vxyzut(4,p) + axyzutp(4,p)*dtmp_min +         &
                       enucp(p)
         ENDIF
#endif
      ENDDO
!$OMP END DO

      IF (RELFLAG.EQV..false.) THEN
!$OMP DO SCHEDULE(runtime)
         DO p=1,nbody
#ifdef Helium
            vxyzut(5,p) = vxyzutp(5,p) + 0.5*axyzutp(5,p)*dtmp_min +    &
                          luminucp(p)
#else
            vxyzut(5,p) = vxyzut(5,p) + axyzutp(5,p)*dtmp_min +         &
                          luminucp(p)
#endif
         ENDDO
!$OMP END DO
      ENDIF
!$OMP END PARALLEL
!
!--Save eps for late use
!
      eps = 1.0d30
!$OMP PARALLEL DEFAULT(none) shared(nbody,xyzhm) private(p)  &
!$OMP reduction(MIN:eps)
!$OMP DO SCHEDULE(runtime)
      DO p=1,nbody
         eps  = MIN(eps,xyzhm(4,p))
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      eps  = 1.4d0*2.0d0*eps
      eps3 = eps*eps*eps
!
!--Correct He layer velocity
!
#ifdef Helium
         CALL norm_layer
         CALL layer
#endif
!
      END SUBROUTINE predictor
