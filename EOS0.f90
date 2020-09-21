      SUBROUTINE EOS0
!========================================================================
!     THIS   SUBROUTINE   ALLOWS   TO  SELECT   BETWEEN  THE DIFFERENT  
!     EQUATIONS  OF  STATE  OF  THE  CODE. IT'S IMPORTANT  TO  REALIZE  
!     THAT  EVERY  EOS  SUBROUTINE MUST PROVIDE (AT LEAST):
!
!     * TEMPERATURE, PRESSURE, Cv AND SPEED OF SOUND
!
!     Last revision: 15/March/2015!
!=========================================================================
!
!--Load modules
!
      USE mod_commons, ONLY : nbody, partype, rho
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local definitions
!
      INTEGER :: p
!
!--Call the EOS for each active particle
!
      DO p = 1, nbody
         IF (rho(p) == 0) PRINT*, p,rho(p),partype(p)
         IF (partype(p) /= 2) CALL degenerate(p) 
      ENDDO
!
      END SUBROUTINE EOS0
