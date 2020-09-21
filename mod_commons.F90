        MODULE mod_commons
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
        USE mod_parameters, ONLY : nel
!
!--Real matrices
!
        REAL, ALLOCATABLE, DIMENSION(:,:) :: xyzhm, xyzhmp, xss, vxyzut,&
        vxyzutp, axyzut, axyzutp,  gxyzu, gxyzup, rotforc

        !REAL, DIMENSION(5,maxnode) :: xyzhm, xyzhmp
        !REAL, DIMENSION(nel,npart) :: xss
        !REAL, DIMENSION(5,npart)   :: vxyzut, vxyzutp, axyzut,          &
        !                              axyzutp
        !REAL, DIMENSION(4,npart)   :: gxyzu, gxyzup
        !REAL, DIMENSION(3,npart)   :: rotforc
!
!--Real vectors
!
        REAL, ALLOCATABLE, DIMENSION(:) :: cur, div, divt, dtem, dtemp, &
        vsigmax, ka1, fh, cvs, css, press, tem, dPdT, enuc, eneutr, rho,&
        tscdyn, tscnuc, luminuc, cps, uintprev, enucp, luminucp, dhdtp, &
        dhdt, dtnuc, norm, rhoG
!
        REAL, DIMENSION(nel-1) :: aion, zion    

        !REAL, DIMENSION(npart) :: cur, div, divt, dtem, dtemp,          &
        !vsigmax, ka1, fh, cvs, css, press, tem, dPdT, enuc, eneutr, rho,&
        !tscdyn, tscnuc, luminuc, cps, uintprev, enucp, luminucp, dhdtp, &
        !dhdt, dtnuc, error, norm, rhoG
        !REAL, DIMENSION(nel-1) :: aion, zion
!
!--Real numbers
!
        REAL :: eps, eps3, globmax, Omega0, tnow, tend, dtmp_min,       &
        trelax, gamma
!
!--Integer matrices
!
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nb 

        !INTEGER, DIMENSION(nbmax,npart)  :: nb 
!
!--Integer vectors
!
        INTEGER, ALLOCATABLE, DIMENSION(:) :: ilist, plist
        INTEGER, ALLOCATABLE, DIMENSION(:) :: partype, eosflag, star
!        
        !INTEGER(1), DIMENSION(npart) :: partype, star, ilist, plist
!
!--Default integer numbers
!
        INTEGER :: globnmin, globnmax, globnvec, globdone, nstep,       &
        nstep_ini, nstep_infile, nout, nbody, nCO, nHe, ndead, factor,  &
        nbody1, nbody2, nbmax, nbmin, mmax, npart
!
!--Integer numbers
!
        INTEGER :: SIMTYPE, SETTYPE, ierr, size, rank, nprocs
!
!--Logical
!
        LOGICAL :: RELFLAG, balsara, visvar
!
!--Character
!
        CHARACTER(30) :: mode(1)

        END MODULE mod_commons
