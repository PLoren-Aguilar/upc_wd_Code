subroutine sph
!========================================================================
!  This subroutine is the main driver of the code. It evolves the system
!  in time using a given integrator method, and prints diagnostics,
!  whenever necessary
!
!  Last revision: 06/April/2019
!========================================================================
!
!--Load modules
!
  use mod_parameters, only : MASTER, ndim, maxstep
  use mod_commons, only : axyzut, axyzutp, enuc, enucp,             &
      luminuc, luminucp, dhdt, dhdtp, dtnuc, vxyzut, xyzhm, xss,    &
      globnmin, globnmax, globnvec, globmax, globdone, rank, size,  &
      nprocs, nstep, nout, tend, tnow, nbody, npart, ndead, partype,&
      dtmax, RELFLAG, istep, nactive
! 
!--Force to declare EVERYTHING
!
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!--Local definitions
!
  real, dimension(ndim)  :: cmp1, cmp2
  real    :: omega, cm_d, deldis, masa1, masa2, t1, t2, avgt
  integer :: p, k, ierr, nCO, nHe
  logical :: finish
!
!--Startout MPI, if necessary
!
#ifdef MPI
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
#else
  size = 1
  rank = MASTER
#endif
  nprocs = size
!
!--Startout the code reading the necessary data files and parameters
!
  call startout
  nCO   = 0
  nHe   = 0
  ndead = 0
!$omp parallel default(none) shared(npart,partype) private(p) &
!$omp reduction(+:nCO) reduction(+:nHe) reduction(+:ndead) 
!$omp do schedule(runtime)
  do p=1,npart
    if (partype(p) == 0) nCO   = nCO   + 1
    if (partype(p) == 1) nHe   = nHe   + 1
    if (partype(p) == 2) ndead = ndead + 1
  end do
!$omp end do
!$omp end parallel
  nbody = npart - ndead
  nstep = 0
!
!--Kickstart the integrator
!
  call kickstart
!
!--Relax the system a first time
!
  call relax
!
!--Main time evolution loop
!
  nstep  = nstep + 1
  finish = .false.
  timeloop : do while (finish.eqv..false.)
!
!--Advance time
!
    tnow  = tnow  + dtmax
!
!--Evolve system a step
!
    call predcorr
!
!--Sort particles by particle type to avoid killed ones
!
    call sort
!
!--Calculate diagnostics and write outputs
!
    call diagnostics
!
!--Print time-step information and advance time-step
!
    if (rank == MASTER) then
      write(*,'(a,i4,a,1pe12.4,a,1pe12.4)') 'step = ', nstep,' tnow =',tnow,' dtmax =',dtmax
    endif
    nstep = nstep + 1
!
!--Check if simulation has ended
!
    if (tnow > tend) finish = .true.
  end do timeloop 
!
!--Ouput data for the last time
!
  call diagnostics
!
!--Stop MPI
!
#ifdef MPI
  call MPI_FINALIZE(ierr)
#endif
end subroutine sph
