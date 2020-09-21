subroutine kickstart
!==============================================================================================
!  This subroutine call the necessary subroutines to start the
!  integration process
!
!  Last revision: 9/April/2019
!==============================================================================================
!
!--Load modules
!
  use mod_parameters, only : MASTER, ndim, maxstep, log2
  use mod_commons   , only : axyzut, axyzutp, enuc, enucp, luminuc, luminucp, dhdt, dhdtp, &
                            xyzhm, xyzhmp, vxyzut, vxyzutp, nbody, rank, nstep, istep,     &
                            istep0, step, dtmax, partype, npart, nstep_infile, dt0,        &
                            reset, nstep_ini, active, nextstep
! 
!--Force to declare EVERYTHING
!
  implicit none
!
!--Local definitions
!
  integer, dimension(npart) :: step_ini
  integer :: k, m, n, p, minstep, newstep
!
!--Calculate the initial time-step for all particles to kickstart the
!  ITS scheme if we are on the first step
!
  if ((nstep_infile == 0) .or. (reset.eqv..true.)) then
    if (dt0 > dtmax) dt0 = dtmax
    n = int(log10(dtmax/dt0)/log10(2.)) + 1
    if (n > 31) then
      write(*,*) 'dtini too small!'
      stop
    endif
    step(1:npart) = maxstep/2**(n-1)
  endif
  active(1:npart) = .true.
  istep0(1:npart) = 0
  istep = 0
  reset = .false.
!
!--Now search for the next time-step. Careful, you need to find the
!  smallest istep !!!
!
  nextstep = maxstep
  minstep  = maxstep
!$omp parallel default(none) shared(nbody,partype,istep0,istep,step) private(p,newstep) &
!$omp reduction(MIN:nextstep) reduction(MIN:minstep) 
!$omp do schedule(runtime)
  do p = 1,nbody
    if (partype(p) == 2) cycle
    if (istep < maxstep) then
      newstep  = istep0(p)+step(p)
    elseif (istep == maxstep) then
      newstep  = step(p)
    endif
    nextstep = min(nextstep,newstep)
    minstep  = min(minstep,step(p))
  enddo
!$omp end do
!$omp end parallel
  istep = nextstep
  if (nextstep > maxstep) then
    write(*,*) 'Something wrong with nextstep. nextstep > maxstep'
    stop
  endif
!
!--If using global time-steps, set all particles to the minimum
!  time-step
#ifdef global
  do p=1,nbody
    if (partype(p) == 2) cycle
    step(p) = minstep
  enddo
#endif
!
!--Calculation of density, smoothing lenght and tree
!
  call iter_rhoh
#ifdef debug
   if (rank == MASTER) print*, 'iter_rhoh  called'
#endif
!
!--Calculate EOS
!
  if (nstep_infile == 0) then
    call eos0
  else
    call eos
  endif
#ifdef debug
  if (rank == MASTER) print*, 'EOS called'
#endif
!
!--Calculation of hydrodynamical quantities
!
  call hydro_rs
#ifdef debug
  if (rank == MASTER) print*, 'hydro_rs called'
#endif
!
!--Sum forces. Add imaginary forces if necessary
!
  call forces
#ifdef debug
  if (rank == MASTER) print*, 'forces called'
#endif      
!
!--Variables update
!
!$omp parallel default(none) shared(axyzutp,axyzut,enucp,enuc)  &
!$omp shared(luminucp,luminuc,xyzhmp,xyzhm,vxyzutp,vxyzut)      &
!$omp shared(nbody,istep0,step,istep,nstep) private(m,p,k)
!$omp do schedule(runtime)
  do p=1,nbody
    if ((nstep == 0) .or. (istep0(p) + step(p) == istep)) then
      do k=1,ndim+2
        axyzutp(k,p) = axyzut(k,p)
      enddo
      enucp(p)    = enuc(p)
      luminucp(p) = luminuc(p)
    endif
  enddo
!$omp end do
!$omp end parallel
!
!--First time-steps calculation
!
   !call varydt
#ifdef debug
  if (rank == MASTER) print*, 'varydt called'
#endif
!
end subroutine kickstart
