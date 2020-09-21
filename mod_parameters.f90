        MODULE mod_parameters
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
!--Load modules
!
        !USE mod_essentials
!
!--Mathematical constants
!
        REAL, PARAMETER :: pi=3.141592654
!
!--Physical constants (in cgs!!!)
!
        REAL, PARAMETER :: unl=6.96e9, unm=1.989e33, g=6.6726e-8,       &
        unt2=unl**3/unm/g, unt=SQRT(unt2), uen=unm*unl*unl/unt2,       &
        up=uen/(unl**3), uden=unm/(unl**3), uv=unl/unt
!
!--Tree calculation
!
        REAL, PARAMETER :: tol=0.7, bigno = 1.0e30, hashfac=1.0
        !INTEGER, PARAMETER :: mmax = 2*npart + 2, nbmax = 1000,        &
        !                      maxpas=10000, maxit = 1000000
        INTEGER, PARAMETER :: maxpas=10000, maxit = 1000000
!
!--Limits
!
        REAL, PARAMETER :: hfact=1.2,  hmin=1.0e-6, hmax=1.0e-1
        REAL, PARAMETER :: rhomin = 1.0, rhomax = 1e8
        REAL, PARAMETER :: tmax=1.0e10, tmin=2.0e3 
!
!--Viscosity
!
        REAL, PARAMETER :: ka_min=0.1, ka_max=1.5
!
!--Code dimensions
!
        INTEGER, PARAMETER :: ndim=3
        INTEGER, PARAMETER :: nel=16
!
!--MPI stuff
!
        INTEGER, PARAMETER :: MASTER=0
!
        END MODULE mod_parameters
