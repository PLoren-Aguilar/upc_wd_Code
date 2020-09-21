      SUBROUTINE detonation(outdet1,outdet2)
!===========================================================
!
!  This subroutine checks the detonation conditions
!
!  Last revision: 15/March/2015
!===========================================================
!
!--Load modules
!
      USE mod_essentials
      USE mod_commons, ONLY : vxyzut, rho
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--I/O variables
!
      LOGICAL, INTENT(OUT) :: outdet1, outdet2
!
!--Check detonation condition 1
!
      outdet1=.FALSE.
      outdet2=.FALSE.
!$OMP PARALLEL DEFAULT(none) shared(nbody,rho,vxyzut,outdet1,outdet2) &
!$OMP private(p)
!$OMP DO SCHEDULE(runtime)
      DO p=1,nbody
         IF (rho(p)*uden.GT.2.0d6) THEN
            IF (vxyzut(5,p).GT.2.5d9)  THEN
!$OMP ATOMIC
               outdet1 = .TRUE.
            ENDIF
         ENDIF
         IF ((tscnuc(p)-tscdyn(p)).LE.0.0d0) THEN
!$OMP ATOMIC
            outdet2 = .TRUE.
         ENDIF
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
! 
      END SUBROUTINE detonation
