       SUBROUTINE setup(omega)
!===========================================================================
!
!      THIS SUBROUTINE SETS UP INITIAL CONDITIONS
!
!      Last revision: 15/March/2015
!
!===========================================================================
!
!--Load modules
!
       USE mod_essentials
       USE mod_commons, ONLY : xyzhm, vxyzut
!
!--Force to declare EVERYTHING
!
       IMPLICIT NONE
!
!--I/O variables
!
       REAL, INTENT(IN) :: omega
!
!--Local variables
!
       INTEGER :: p
!
!--Adjust body data to the initial conditions of the simulation
!
       DO p=1,nbody
          vxyzut(1,p) = -xyzhm(2,p)*omega
          vxyzut(2,p) =  xyzhm(1,p)*omega
          vxyzut(3,p) =  0.0d0
       ENDDO
!
       END SUBROUTINE setup
