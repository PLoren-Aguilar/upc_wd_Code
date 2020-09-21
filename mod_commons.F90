module mod_commons
!===================================================================
!  This module contains all the common variables of the code
!
! * xyzhm   = Contains x, y, z, mass and h of the particles
! * xss     = Contains the chemical mass abundances of the particles
! * vxyzu   = Contains vx, vy, vz and uint of the particles
! * axyzu   = Contains ax, ay, az and dudt
! * axyzup  = Contains axp, ayp, azp, and dudtp
! * gxyzu   = Contains gx, gy, gz and phi
! * cur     = Contains the norm of curl of the velocity
! * div     = Contains the divergence of the velocity
! * dtem    = Contains the time-variation of the temperature
! * dtemp   = Contains the time-variation of the temperature at
!             the previous time-step
! * vsigmax = Contains the maximum signal velocity. Necessary
!             for the use of artificial viscosity
! * ka1     = Artificial viscosity coefficient
! * fh      = Contains the grah-h correction term
! * cvs     = Contains the specific heat per unit mass, at 
!             constant volume
! * css     = Contains the sound speed
! * press   = Contains the gas pressure
! * tem     = Contains the gas temperature
! * dPdT    = Contains the variation of gas pressure with
!             respect to temperature
! * enuc    = Contains the nuclear energy per unit mass
! * eneutr  = Contains the neutrino energy per unit mass?
!===================================================================
  use mod_parameters, only : nel, ndim
!
!--Real matrices
!
  real, allocatable, dimension(:,:) :: xyzhm, xyzhmp, xss, vxyzut,&
        vxyzutp, axyzut, axyzutp,  gxyzu, gxyzup, rotforc, fix
!
!--Real vectors
!
  real, allocatable, dimension(:) :: cur, div, divt, dtem, dtemp,       &
        vsigmax, ka1, fh, cvs, css, press, tem, dPdT, enuc, eneutr, rho,&
        tscdyn, tscnuc, luminuc, cps, uintprev, enucp, luminucp, dhdtp, &
        dhdt, dtnuc, norm, rhoG
  real, dimension(nel-1) :: aion, zion
!
!--Real numbers
!
  real :: eps, eps3, globmax, Omega01, Omega02, dtmax, dtmaxin, tnow, tend, &
          trelax, gamma, dt0, dt0in, fixx1, fixy1, fixz1, fixx2,            &
          fixy2, fixz2, Omega0, ttot, tpred, tirho, tEOS,                   &
          thydro, tforc, tcorr, tdt, tsort, ttree1, ttree2, tnuc
!
!--Integer matrices
!
  integer, allocatable, dimension(:,:) :: nb 
!
!--Integer vectors
!
  integer, allocatable, dimension(:) :: ilist, plist, step, istep0
  integer, allocatable, dimension(:) :: partype, eosflag, star
!
!--Default integer numbers
!
  integer :: globnmin, globnmax, globnvec, globdone, nstep,                  &
             nstep_ini, nstep_infile, nout, nbody, nCO, nHe, ndead, factor,  &
             nbody1, nbody2, nbmax, nbmin, mmax, npart, istep, nactive,      &
             npart0, npart1, nextstep
!
!--Integer numbers
!
  integer :: SIMTYPE, SETTYPE, ierr, size, rank, nprocs
!
!--Logical
!
  logical, allocatable, dimension(:) :: active, pseudoactive
  logical :: RELFLAG, balsara, visvar, reset, SYNC
!
!--Character
!
  character(30) :: mode(1)
end module mod_commons
