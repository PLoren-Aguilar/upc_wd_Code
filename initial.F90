      SUBROUTINE initial(cmp1,masa1,cmp2,masa2,omega,cm_d)
!===================================================================
!
!  This subroutine sets up initial conditions for the binary system.
!
!  Last revision: 15/March/2015
!
!===================================================================
!
!--Load modules
!
      USE mod_essentials
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--I/O definitions
!
      REAL, INTENT(OUT) :: cmp1, masa1, cmp2, masa2, omega, cm_d

!
!--A small kick is included to force the merger. Should only be used to
!  force a merger if initial conditions are set too far
!
      CALL approach(0.5)
      CALL relax_frame
!
#ifdef debug
      IF (rank.EQ.MASTER) PRINT*, 'initial called'
#endif
!
      END SUBROUTINE initial
