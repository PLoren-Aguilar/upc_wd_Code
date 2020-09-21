        MODULE mod_essentials
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
        REAL :: tnow, tend, dtmp_min, trelax, gamma
!
!--Integer parameters
!
        !INTEGER, PARAMETER :: nbody1 = 350000, nbody2 = 350000,         &
        !npart  = nbody1 + nbody2
        !INTEGER, PARAMETER :: maxnode = npart + npart + 2
        !INTEGER(1), PARAMETER :: nprocs=1
        INTEGER(1), PARAMETER :: ndim=3
        INTEGER(1), PARAMETER :: nel=16
        INTEGER(1), PARAMETER :: MASTER=0
!
!--Integer lists
!
        INTEGER, DIMENSION(npart) :: ilist, plist
!
!--Integers
!
        INTEGER    :: nstep, nstep_ini, nout, nbody, nCO, nHe, ndead,   &
                      factor
        INTEGER(1) :: SIMTYPE, SETTYPE, ierr, size, rank
!
!--Logical
!
        LOGICAL :: RELFLAG, balsara, visvar
!
        END MODULE mod_essentials
