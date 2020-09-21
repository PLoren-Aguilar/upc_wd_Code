      SUBROUTINE separate
!========================================================================
!
!  This subroutine is a call to various analysis subroutines
!
!  Last revision: 14/March/2015
!========================================================================
!
!--Load modules
!
      USE mod_commons, ONLY : nstep, nout
! 
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Startout the code reading the necessary data files and parameters
!
      CALL startout
!
!--Analize the results
!
      CALL separation
!
!--Write results
!
      nstep = 1
      nout  = 1
      CALL outdata
!
      END SUBROUTINE separate
!
      SUBROUTINE separation
!
!--Load modules
!
      USE mod_commons, ONLY : xyzhm, vxyzut, rho, nbody1, nbody
      USE mod_parameters, ONLY : unm, unl, unt, uden
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local definitions
!
      REAL :: masa ,mass, dens, temp, rin, rout, rmax, dr, r,        &
                 xcm, ycm ,zcm, angx, angy, angz, Imom, RWD, omega, v
      INTEGER :: p, k, n
!
!--Position all data around the centre of the primary star
!
      DO p = 1, nbody1
         xyzhm(1,p) = xyzhm(1,p) - 0.03
      ENDDO
!
      DO p = nbody1+1,nbody
         xyzhm(1,p) = xyzhm(1,p) + 0.03
      ENDDO

      END SUBROUTINE separation
