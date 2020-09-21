module mod_essentials
!===============================================================
!  This module contains all the "essential" variables
!  and parameters of the code
!
! * tnow     = Global time-step of the simulation
! * tend     = Final time of the simulation
! * dtmp_min = Integration time increment
! * trelax   = Relaxation coeficient
! * nbody1   = Number of particles of the first star
! * nbody2   = Number of particles of the second star
! * npart    = Total number of particles in the simulation
! * nprocs   = Number of nodes in an MPI calculation
! * ndim     = Number of spatial dimensions 
! * nel      = Number of elements in the Nuclear network
! * nstep    = Number of performed integration steps
! * nout     = Number of steps necessary to write an output
! * SIMTYPE  = Type of simulation
! * SETTYPE  = Setup flag
! * RELFLAG  = Relaxation flag
! * balsara  = Balsara switch flag
! * visvar   = Artificial viscosity flag
!===============================================================
!
!--Reals
!
 real :: tnow, tend, trelax, gamma
!
!--Integer parameters
!
 integer(1), parameter :: ndim=3
 integer(1), parameter :: nel=16
 integer(1), parameter :: MASTER=0
!
!--Integer lists
!
 integer, dimension(npart) :: ilist, plist
!
!--Integers
!
 integer    :: nstep, nstep_ini, nout, nbody, nCO, nHe, ndead, factor
 integer(1) :: SIMTYPE, SETTYPE, ierr, size, rank
!
!--Logical
!
 logical :: RELFLAG, balsara, visvar
!
end module mod_essentials
