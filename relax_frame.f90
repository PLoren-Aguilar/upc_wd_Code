       SUBROUTINE relax_frame
!===========================================================================
!      THIS SUBROUTINE SUBSTRACTS CM VELOCITY IN ORDER TO CORRECTLY RELAX
!
!      Last revision: 15/March/2015
!===========================================================================
!
!--Load modules
!
       USE mod_parameters, ONLY : ndim
       USE mod_commons, ONLY : xyzhm, vxyzut, nbody
!
!--Force to declare EVERYTHING
!
       IMPLICIT NONE
!
!--Local variables
!
       REAL, DIMENSION(3) :: cmp(3),cmv(3)
       REAL :: mtot, mass
       INTEGER :: k, m, p
!
!--Loop over particles to calculate center of mass velocity and position of the system
!
       mtot = 0.0
       cmp  = 0.0
       cmv  = 0.0
!$OMP PARALLEL DEFAULT(none) shared(xyzhm,nbody) private(m,p) &
!$OMP reduction(+:mtot)
!$OMP DO SCHEDULE(runtime)
       DO p=1,nbody
          mtot = mtot + xyzhm(5,p)
       ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
       DO k=1,ndim
!$OMP PARALLEL DEFAULT(none) shared(k,xyzhm,vxyzut,nbody) private(p,mass) &
!$OMP reduction(+:cmp) reduction(+:cmv)
!$OMP DO SCHEDULE(runtime)
          DO p=1,nbody
             mass = xyzhm(5,p)
             cmp(k) = cmp(k) + mass*xyzhm(k,p)
             cmv(k) = cmv(k) + mass*vxyzut(k,p)
          ENDDO
!$OMP END DO
!$OMP END PARALLEL
       ENDDO
       cmp = cmp/mtot
       cmv = cmv/mtot
!
!--Remove any possible center-of-mass offset
!
       DO k=1,ndim
!$OMP PARALLEL DEFAULT(none) shared(nbody,k,xyzhm,vxyzut,cmp,cmv) private(m,p)
!$OMP DO SCHEDULE(runtime)
          DO p=1,nbody
             xyzhm(k,p)  = xyzhm(k,p)  - cmp(k)
             vxyzut(k,p) = vxyzut(k,p) - cmv(k)
          ENDDO
!$OMP END DO
!$OMP END PARALLEL
       ENDDO
!
       END SUBROUTINE relax_frame
