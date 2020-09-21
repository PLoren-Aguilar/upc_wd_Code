subroutine setup(omega)
!===========================================================================
!
! This subroutine setup up initial conditions
!
! Last revision: 9/April/2019 
!
!===========================================================================
!
!--Load modules
!
  use mod_essentials
  use mod_commons, only : xyzhm, vxyzut
!
!--Force to declare EVERYTHING
!
  implicit NONE
!
!--I/O variables
!
  real, intent(in) :: omega
!
!--Local variables
!
  integer :: p
!
!--Adjust body data to the initial conditions of the simulation
!
  do p=1,nbody
    vxyzut(1,p) = -xyzhm(2,p)*omega
    vxyzut(2,p) =  xyzhm(1,p)*omega
    vxyzut(3,p) =  0.0d0
  enddo

end subroutine setup
