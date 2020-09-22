subroutine eos0
!========================================================================
!     This subroutine allows to select between the different equations  
!     of state in the code. It's important to realize tha every eos must
!     provide at least:
!
!     * Temperature, pressure, cv and speed of sound
!
!     Last revision: 06/April/2019 
!=========================================================================
!
!--Load modules
!
  use mod_commons, only : nbody, partype, rho, istep0, istep, step

  implicit none
!
!--Local definitions
!
  integer :: p
!
!--Call the eos for each active particle
!
  do p = 1, nbody
     if (partype(p) /= 2) call degenerate(p)
  enddo
!
end subroutine eos0
