subroutine varydt
!===================================================================
!  This subroutine calculates the new time step
!
!  Last revision: 15/March/2015
!===================================================================
!
!--Load modules
!
  use mod_commons,    only : axyzutp, xyzhmp, vsigmax, vxyzutp, uintprev,                         &
                             rho, partype, nbody, nstep, nstep_ini, gxyzu, enuc, tscdyn,          &
                             tscnuc, cps, dtnuc, press, css, tnow, eosflag, rank, RELFLAG,        &
                             step, dtmax, istep, istep0, nactive, rank, SIMTYPE, active, nextstep,&
                             SYNC
  use mod_parameters, only : MASTER, maxstep, log2
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--Local definitions
!
  real :: fa, ft, dt, dt1, dt2, dt3, dt4, dt5, dt6, dtmpmin_old, denom, dtmin, frac, dtnew
  real, parameter :: dtfact=0.3, sigma=1.0, dtlow = 1.e-8
  integer :: p, n, newstep, syncstep, minstep
!
!--If reached the end of a integration cycle, activate the flag to stop
!
  if (istep == maxstep) SYNC = .true.
!
!--New time-step calculation
!
!$omp parallel default(none) shared(axyzutp,xyzhmp,vsigmax,vxyzutp)      &
!$omp shared(nbody,uintprev,rho,partype,istep0,step,istep,nstep)         &
!$omp shared(RELFLAG,dtmax,rank,SIMTYPE) private(p,fa,dtmin,newstep,syncstep) &
!$omp private(dt1,dt2,dt3,dt4,ft,frac,n,dtnew)
!$omp do schedule(runtime)
  partloop : do p=1,nbody
     if ((nstep /=0) .and. (istep0(p) + step(p) /= istep)) cycle
     if (partype(p) == 2) cycle
!#ifdef Helium
!    if (partype(p) /= 0) cycle
!#else
!    if (partype(p) == 2) cycle
!#endif

! Reset all time-scales
     dt1 = 1.0e30
     dt2 = 1.0e30
     dt3 = 1.0e30
     dt4 = 1.0e30

! Dynamical Courant condition
     fa = sqrt(axyzutp(1,p)*axyzutp(1,p) + axyzutp(2,p)*axyzutp(2,p) + axyzutp(3,p)*axyzutp(3,p))
     if (fa /= 0.0) then
        dt1 = min(sqrt(xyzhmp(4,p)/fa),dt1)
     endif

! Viscous signal condition
     if (vsigmax(p) /= 0.0) then
        dt2 = min(sigma*xyzhmp(4,p)/vsigmax(p),dt2)
     endif

! Thermal energy condition
     if (RELFLAG.eqv..false.) then
        ft = abs(axyzutp(4,p))
        if ((ft /= 0.0) .and. (vxyzutp(4,p) > 0.0)) then
           dt3 = min(dt3,vxyzutp(4,p)/ft)
        endif

! Temperature condition

        ft = abs(axyzutp(5,p))
        if (ft /= 0.0) then
          dt4 = min(dt4,vxyzutp(5,p)/ft)
        endif
     endif

! Calculate the shortest time-scale
     if ((dt1 == 0) .and. (dt2 /= 0)) then
        dtmin = dt2
     elseif ((dt1 /= 0) .and. (dt2 == 0)) then
        dtmin = dt1
     else
        dtmin = min(dt1,dt2)
     endif
     dtmin = min(dtmin,dt3)
     dtmin = min(dtmin,dt4)
     dtmin = dtfact*dtmin
     if (dtmin > dtmax) then
        dtmin = dtmax
        !if (rank == MASTER) then
        !  write(*,*) 'ERROR in time-steps. dtmin > dtmax !!!',dt1,dt2,dt3,dt4,dtmin,dtmax
        !  write(*,*) 'This is **not** the optimal way to run'
        !  stop
        !endif
     elseif (dtmin < dtlow) then
        partype(p) = 2
     endif

! Update the individual time-step
     frac = dtmin/(dtmax*step(p)/maxstep)
     if (frac >= 2.) then
        syncstep = maxstep - istep
        if (syncstep == 0) syncstep = maxstep
        newstep  = 2.*step(p)
        if (mod(syncstep,newstep) == 0) then
              step(p) = min(newstep,maxstep)
        endif
     elseif (frac <= 0.5) then
        if (dtmin <= 0) then
           PRint*, 'dtmin zeroooo',dt1,dt2,dt3,dt4
           stop
        endif
        dtnew = dtmax*step(p)/maxstep*frac
        n = int(log10(dtmax/dtnew)/log10(2.0)) + 1
        if (n > 30) then  ! Kill the particle if the time-step is too small
           partype(p) = 2
        else
           step(p) = maxstep/2**(n-1)
        endif
     endif

! Update the present time-step for the particle
     if (istep < maxstep) then
        istep0(p) = istep
     elseif (istep == maxstep) then
        istep0(p) = 0
     else
        write(*,*) 'Varydt: ERROR with istep !!!'
        stop
     endif
  enddo partloop
!$omp end do
!$omp end parallel
!
!--Now search for the next time-step. Careful, you need to find the
!  smallest istep !!!
!
  nextstep = maxstep
  minstep  = maxstep
!$omp parallel default(none) shared(nbody,partype,istep0,istep,step) private(p,newstep) &
!$omp reduction(min:nextstep) reduction(min:minstep) 
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
end subroutine varydt
