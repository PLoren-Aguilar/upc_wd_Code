subroutine corrector
!============================================================
!  This subroutine corrects the the integrated quantities 
!
!  Last revision: 06/April/2019
!============================================================
!
!--Load modules
!
  use mod_parameters, only : tmin, tmax, ndim, maxstep, xbox, ybox, zbox, MASTER
  use mod_commons,    only : xyzhm, vxyzut, axyzut, xyzhmp, vxyzutp, axyzutp,             &
                             nbody, enuc, enucp, luminuc, luminucp, RELFLAG, dhdt,        &
                             dhdtp, npart, istep0, istep, nextstep, dtmax, step, partype, &
                             SIMTYPE, active, nactive, rank, nstep
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--Local variables
!
  real, dimension(npart) :: temporal
  real :: fa, dt, dt1, dt2, dt3, denom
  real, PARAMETER :: dtfact_b=0.3d0, dtfact_h=0.3d0, sigma=1.0d0,   &
                     f1 = 1./256., f2 = 255./256.
  integer, dimension(npart) :: itemporal
  integer :: p, k
!
!--Correct positions
!
  do k=1,ndim
!$omp parallel default(none) shared(nbody,xyzhm,xyzhmp,vxyzut,vxyzutp)  &
!$omp shared(axyzut,axyzutp,k,dtmax,istep0,step,istep) private(p,dt)
!$omp do schedule(runtime)
     do p=1,nbody
        if (istep0(p) + step(p) == istep) then
           dt = dtmax*float(step(p))/float(maxstep)
#ifdef Helium
           xyzhm(k,p)=xyzhmp(k,p)+(f1*vxyzutp(k,p)+f2*vxyzut(k,p))*dt
#else
           xyzhm(k,p)=xyzhm(k,p)+0.1667*(axyzut(k,p)-axyzutp(k,p))*dt*dt
#endif
        endif
     enddo
!$omp end do
!$omp end parallel
  enddo
!
!--Correct velocities
!
  do k=1,ndim
!$omp parallel default(none) shared(vxyzut,vxyzutp,axyzut,axyzutp,k)            &
!$omp shared(nbody,istep,istep0,step,dtmax) private(p,dt)
!$omp do schedule(runtime)
     do p=1,nbody
        if (istep0(p) + step(p) == istep) then
           dt = dtmax*float(step(p))/float(maxstep)
#ifdef Helium
           vxyzut(k,p)=vxyzutp(k,p)+(f1*axyzutp(k,p)+f2*axyzut(k,p))*dt
#else
           vxyzut(k,p)=vxyzut(k,p)+0.5*(axyzut(k,p)-axyzutp(k,p))*dt
#endif
        endif
     enddo
!$omp end do
!$omp end parallel
  enddo
!
!--Correct thermal energies and temperatures if necessary
!
!$omp parallel default(none) shared(vxyzut,vxyzutp,axyzut,axyzutp,partype)      &
!$omp shared(nbody,enuc,enucp,luminuc,luminucp,RELFLAG,istep0,istep,step,dtmax) &
!$omp shared(SIMTYPE) private(p,dt)
!$omp do schedule(runtime)
  do p=1,nbody
     if (istep0(p) + step(p) == istep) then
        dt = dtmax*float(step(p))/float(maxstep)
#ifdef Helium
        vxyzut(4,p)=vxyzutp(4,p)+(f1*axyzutp(4,p)+f2*axyzut(4,p))*dt + enuc(p)-enucp(p)
#else
        if ((vxyzut(5,p) > tmin).and.(vxyzut(5,p) < tmax)) then
           vxyzut(4,p) = vxyzut(4,p) + 0.5*(axyzut(4,p)-axyzutp(4,p))*dt + enuc(p)-enucp(p)
        endif
#endif
     endif
  enddo
!$omp end do
  if (RELFLAG.eqv..false.) then
!$omp do schedule(runtime)
     do p=1,nbody
        if (istep0(p) + step(p) == istep) then
           dt = dtmax*float(step(p))/float(maxstep)
#ifdef Helium
           vxyzut(5,p) = vxyzutp(5,p) + (f1*axyzutp(5,p)+f2*axyzut(5,p))*dt + luminuc(p)-luminucp(p)
#else   
           vxyzut(5,p) = vxyzut(5,p) + 0.5*(axyzut(5,p)-axyzutp(5,p))*dt + luminuc(p)-luminucp(p)
#endif
        endif
     enddo
!$omp end do
  endif
!$omp end parallel
!
!--Correct He layer velocity
!
!#ifdef Helium
  call norm_layer
  call layer
!#endif
!
!--Save quantities for later use
!
!$omp parallel default(none) shared(axyzut,axyzutp,vxyzut,vxyzutp,xyzhm,xyzhmp) &
!$omp shared(nbody,istep0,step,istep,partype) private(p)
!$omp do schedule(runtime)
  do p=1,nbody
    if (istep0(p) + step(p) == istep) then
      xyzhmp(1:ndim+2,p)  = xyzhm(1:ndim+2,p)
      vxyzutp(1:ndim+2,p) = vxyzut(1:ndim+2,p)
      axyzutp(1:ndim+2,p) = axyzut(1:ndim+2,p)
    endif
  enddo
!$omp end do
!$omp end parallel
!
!$omp parallel default(none) shared(enucp,enuc,luminucp,luminuc,dhdtp)  &
!$omp shared(dhdt,nbody,istep0,step,istep) private(p)
!$omp do schedule(runtime)
  do p=1,nbody
     if (istep0(p) + step(p) == istep) then
        enucp(p)    = enuc(p)
        luminucp(p) = luminuc(p)
        dhdtp(p)    = dhdt(p)
     endif
  enddo
!$omp end do
!$omp end parallel
!
!--Count the number of active particles in the step
!
  active  = .false.
  nactive = 0
!$omp parallel default(none) shared(nbody,partype,istep0,step,istep,nstep,active) &
!$omp private(p) reduction(+:nactive) 
!$omp do schedule(runtime)
  do p=1,nbody
     if (partype(p) == 2) cycle
     if (istep0(p) + step(p) == istep) then
        nactive   = nactive + 1
        active(p) = .true.
     endif
  enddo
!$omp end do
!$omp end parallel
!
  if (rank == MASTER) then
     print*, 'istep=',istep,' out of', maxstep
     print*, nactive,' particles moved'
  endif
!
  if (rank == MASTER) then
     if (nactive == 0) then
        do p=1,nbody
           if (partype(p) == 2) cycle
           print*, nstep,p,istep0(p),step(p),istep,maxstep
        enddo
        stop
     endif
  endif

end subroutine corrector
