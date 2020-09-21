      SUBROUTINE forces
!===================================================================
!
!  This subroutine adds up the hydro and gravitational forces.
!  It also adds coriolis and dumping forces in a relaxation.
!
!  Last revision: 9/October/2015
!
!===================================================================
!
!--Locad modules
!
      USE mod_parameters, ONLY : ndim
      USE mod_commons, ONLY : axyzut, gxyzu, xyzhm, vxyzut, rotforc,    &
      nbody, rotforc, SIMTYPE, trelax
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local variables
!
      INTEGER :: k, p, m
!
!--Calculate non-inertial forces
!
      SELECT CASE (SIMTYPE) 
         CASE (1)
            rotforc = 0.0
         CASE (2)
            CALL noninertial
         CASE DEFAULT
            PRINT*, 'Bad SIMTYPE !!!', SIMTYPE
            STOP
      END SELECT
!
!--Add gravitational, hydrodynamical, and inertial forces
!
      DO k=1,ndim
!$OMP PARALLEL DEFAULT(none) shared(nbody,axyzut,gxyzu,k,trelax,vxyzut) &
!$OMP shared(rotforc) private(p)
!$OMP DO SCHEDULE(runtime)
         DO p=1,nbody
            axyzut(k,p) = axyzut(k,p) + gxyzu(k,p) -                    &
                          trelax*vxyzut(k,p) + rotforc(k,p)
         ENDDO
!$OMP END DO
!$OMP END PARALLEL
      ENDDO
!
      END SUBROUTINE 
