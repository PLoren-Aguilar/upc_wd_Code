subroutine initial(cmp1,masa1,cmp2,masa2,omega,cm_d)
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
  use mod_essentials
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--I/O definitions
!
  real, intent(out) :: cmp1, masa1, cmp2, masa2, omega, cm_d

!
!--A small kick is included to force the merger. Should only be used to
!  force a merger if initial conditions are set too far
!
  call approach(0.5)
  call relax_frame
!
#ifdef debug
  if (rank == MASTER) print*, 'initial called'
#endif
!
end subroutine initial
