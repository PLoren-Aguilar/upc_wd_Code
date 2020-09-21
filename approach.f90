       SUBROUTINE approach(frac,FINISH)
!===========================================================================
!
!      THIS SUBROUTINE SUBSTRACTS CM VELOCITY IN ORDER TO CORRECTLY
!      RELAX
!
!      Last revision: 15/March/2015
!
!===========================================================================
!
!--Load modules
!
       USE mod_parameters, ONLY : ndim
       USE mod_commons, ONLY : xyzhm, nbody, nbody1, nbody2
!
!--Force to declare EVERYTHING
!
       IMPLICIT NONE
!
!--I/O variables
!
       REAL, INTENT(IN)  :: frac
       LOGICAL, INTENT(OUT) :: FINISH
!
!--Local variables
!
       REAL, DIMENSION(ndim) :: cmp1, cmp2
       REAL :: mtot1, mtot2, m, d1, d2, gravpot1, gravpot2 
       INTEGER :: p, k
!
!--Loop over particles to calculate center of mass velocity
!  and position of the first star
!
       cmp1=0.0d0
       mtot1=0.0d0
       DO k=1,ndim
!$OMP PARALLEL DEFAULT(none) shared(k,xyzhm,nbody1) private(p,m) &
!$OMP reduction (+:cmp1) reduction(+:mtot1)
!$OMP DO SCHEDULE(runtime)
          DO p=1,nbody1
             m = xyzhm(5,p)
             mtot1 = mtot1 + m
             cmp1(k) = cmp1(k) + m*xyzhm(k,p)
          ENDDO
!$OMP END DO
!$OMP END PARALLEL
       ENDDO
       cmp1 = cmp1/mtot1
       d1   = SQRT(cmp1(1)**2+cmp1(2)**2+cmp1(3)**2)
!
       cmp2=0.0d0
       mtot2=0.0d0
       DO k=1,ndim
!$OMP PARALLEL DEFAULT(none) shared(nbody,nbody1,k,xyzhm) private(p,m) &
!$OMP reduction (+:cmp2) reduction(+:mtot2)
!$OMP DO SCHEDULE(runtime)
          DO p=nbody1+1,nbody
             m = xyzhm(5,p)
             mtot2 = mtot2 + m
             cmp2(k) = cmp2(k) + m*xyzhm(k,p)
          ENDDO
!$OMP END DO
!$OMP END PARALLEL
       ENDDO
       cmp2 = cmp2/mtot2
!
!--Loop over particles to reduce stellar distance
!
!$OMP PARALLEL DEFAULT(none) shared(nbody,nbody1,xyzhm,cmp1,cmp2,d1,d2,frac) &
!$OMP private(p)
!$OMP DO SCHEDULE(runtime)
       DO p=1,nbody1
          xyzhm(1,p) = xyzhm(1,p) - frac*ABS(cmp1(1)-cmp2(2))*cmp1(1)/d1
          xyzhm(2,p) = xyzhm(2,p) - frac*ABS(cmp1(2)-cmp2(2))*cmp1(2)/d1
          xyzhm(3,p) = xyzhm(3,p) - frac*ABS(cmp1(3)-cmp2(3))*cmp1(3)/d1
       ENDDO
!$OMP END DO
!
!$OMP DO SCHEDULE(runtime)
       DO p=nbody1+1,nbody
          xyzhm(1,p) = xyzhm(1,p) - frac*ABS(cmp1(1)-cmp2(2))*cmp2(1)/d2
          xyzhm(2,p) = xyzhm(2,p) - frac*ABS(cmp1(2)-cmp2(2))*cmp2(2)/d2
          xyzhm(3,p) = xyzhm(3,p) - frac*ABS(cmp1(3)-cmp2(3))*cmp2(3)/d2
       ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
!--Calculate gravitational potentials and check if there is RL overflow
!
       FINISH = .false.
       !DO p=1,nbody1
       !   gravpot1 = mtot1/SQRT((xyzhm(1,p)-cmp1(1))**2 +              &
       !                (xyzhm(2,p)-cmp1(2))**2 +                        &
       !                (xyzhm(3,p)-cmp1(3))**2)
       !   gravpot2 = mtot2/SQRT((xyzhm(1,p)-cmp2(1))**2 +              &
       !                (xyzhm(2,p)-cmp2(2))**2 +                        &
       !                (xyzhm(3,p)-cmp2(3))**2)
!!
!          IF ((gravpot1-gravpot2).LE.0.0.AND.ABS(gravpot1-gravpot2)  &
       !       .GT.5d-3) THEN
       !      FINISH=.true.
       !   ENDIF
       !ENDDO
!
       END SUBROUTINE
