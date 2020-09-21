subroutine separate
!========================================================================
!
!  This subroutine is a call to various analysis subroutines
!
!  Last revision: 14/March/2015
!========================================================================
!
!--Load modules
!
  use mod_commons, only : nstep, nout
! 
!--Force to declare EVERYTHING
!
  implicit none
!
!--Startout the code reading the necessary data files and parameters
!
  call startout
!
!--Analize the results
!
  call separation
!
!--Write results
!
  nstep = 1
  nout  = 1
  call outdata
!
end subroutine separate
!
subroutine separation
!
!--Load modules
!
  use mod_commons,    only : xyzhm, vxyzut, rho, nbody1, nbody
  use mod_parameters, only : unm, unl, unt, uden
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--Local definitions
!
  real    :: masa ,mass, dens, temp, rin, rout, rmax, dr, r,        &
          xcm, ycm ,zcm, angx, angy, angz, Imom, RWD, omega, v
  integer :: p, k, n
!
!--Position all data around the centre of the primary star
!
  do p = 1, nbody1
    xyzhm(1,p) = xyzhm(1,p) - 0.03
  enddo
!
  do p = nbody1+1,nbody
    xyzhm(1,p) = xyzhm(1,p) + 0.03
  enddo

end subroutine separation
