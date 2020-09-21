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
      USE mod_essentials
      USE mod_commons, ONLY : uintprev, vxyzut, axyzutp, dtmp_min,      &
                              enucp, luminucp, xyzhm, dhdtp, eps, eps3, &
                              vxyzutp
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
      DO 10 p=1,nbody
         uintprev(p) = vxyzutp(4,p)
10    ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
!--Evolve positions
!
      DO 20 k=1,ndim
!$OMP PARALLEL DEFAULT(none) shared(nbody,xyzhm,vxyzut,vxyzutp,dtmp_min) &
!$OMP shared(axyzutp,k) private(p)
!$OMP DO SCHEDULE(runtime)
         DO 30 p=1,nbody
            xyzhm(k,p) = xyzhm(k,p) + vxyzut(k,p)*dtmp_min +            &
                         0.5d0*axyzutp(k,p)*dtmp_min*dtmp_min
30       ENDDO
!$OMP END DO
!$OMP END PARALLEL
20    ENDDO 
!
!--Evolve velocities
!
      DO 40 k=1,ndim
!$OMP PARALLEL DEFAULT(none) shared(nbody,vxyzut,axyzutp,dtmp_min,k) &
!$OMP private(p)
!$OMP DO SCHEDULE(runtime)
         DO 50 p=1,nbody
            vxyzut(k,p) = vxyzut(k,p) + axyzutp(k,p)*dtmp_min
50       ENDDO
!$OMP END DO
!$OMP END PARALLEL
40    ENDDO
!
!--Evolve thermal energies and temperatures, if necessary
!
!$OMP PARALLEL DEFAULT(none) shared(nbody,vxyzut,axyzutp,dtmp_min) &
!$OMP shared(enucp,luminucp,RELFLAG) private(p)
!$OMP DO SCHEDULE(runtime)
      DO 60 p=1,nbody
         vxyzut(4,p) = vxyzut(4,p) + axyzutp(4,p)*dtmp_min +            &
                       enucp(p)
60    ENDDO
!$OMP END DO
      IF (RELFLAG.EQV..false.) THEN
!$OMP DO SCHEDULE(runtime)
         DO 70 p=1,nbody
            vxyzut(5,p) = vxyzut(5,p) + axyzutp(5,p)*dtmp_min +         &
                          luminucp(p)
70       ENDDO
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
      DO 90 p=1,nbody
         eps  = MIN(eps,xyzhm(4,p))
90    ENDDO
!$OMP END DO
!$OMP END PARALLEL
      eps  = 1.4d0*2.0d0*eps
      eps3 = eps*eps*eps
!
      END SUBROUTINE predictor
