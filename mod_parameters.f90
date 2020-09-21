module mod_parameters
!==========================================================================
!--This module contains all the parameters of the code
!
! * unl   = Code unit of lenght  = 0.1*R_Sun
! * unm   = Code unit of mass    = 1*M_sun
! * unt   = Code unit of time    = 50 s
! * uen   = Code unit of energy 
! * up    = Code unit of power
! * uden  = Code unit of density
! * uv    = Code unit of velocity
! * mmax  = Maximum number of nodes in the tree
! * nbmax = Maximum number of neighbours per particle
!==========================================================================
!
!--Mathematical constants
!
  real, parameter :: pi=3.141592654
!
!--Physical constants (in cgs!!!)
!
  real, parameter :: unl=6.96e9, unm=1.989e33, g=6.6726e-8,                   &
                     unt2=unl**3/unm/g, unt=SQRT(unt2), uen=unm*unl*unl/unt2, &
                     up=uen/(unl**3), uden=unm/(unl**3), uv=unl/unt
!
!--Mathematical constants
!
  real,    parameter :: log2=0.30103
!
!--Tree calculation
!
  real,    parameter :: tol=0.7, bigno = 1.0e30, hashfac=1.0
  integer, parameter :: maxpas=10000, maxit = 1000000
!
!--Limits
!
  real,    parameter :: hfact=1.2,  hmin=1.0e-10, hmax=1.0e1
  real,    parameter :: rhomin = 1e-5, rhomax = 1e10
  real,    parameter :: tmax=1.0e10, tmin=2.0e3
  real,    parameter :: xbox = 10.0, ybox=10. , zbox=10.
  integer, parameter :: maxstep=1073741824
!
!--Viscosity
!
  real,    parameter :: ka_min=0.1, ka_max=1.5
!
!--Code dimensions
!
  integer, parameter :: ndim=3
  integer, parameter :: nel=16
!
!--MPI stuff
!
  integer, parameter :: MASTER=0

end module mod_parameters
