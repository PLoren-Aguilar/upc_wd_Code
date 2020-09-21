subroutine predictor
!============================================================
!  This subroutine predicts in time the integrated quantities 
!
!  Last revision: 6/April/2019
!============================================================
!
!--Load modules
!
  use mod_parameters, only : ndim
  use mod_parameters, only : tmin, tmax, maxstep
  use mod_commons,    only : uintprev, vxyzut, axyzutp, enucp, luminucp, xyzhm, &
                             dhdtp, eps, eps3, xyzhmp, vxyzutp, nbody, RELFLAG, &
                             step, istep, istep0, dtmax, SIMTYPE
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--Local variables
!
  real    :: dt
  integer :: p, k
!
!--Save thermal energy for later use
!
!$omp parallel default(none) shared(nbody,uintprev,vxyzutp) private(p)
!$omp do schedule(runtime)
  do p=1,nbody
     uintprev(p) = vxyzutp(4,p)
  enddo
!$omp end do
!$omp end parallel
!
!--Evolve positions
!
  do k=1,ndim
!$omp parallel default(none) shared(nbody,xyzhm,vxyzut,vxyzutp) &
!$omp shared(xyzhmp,axyzutp,k,dtmax,istep,istep0,step) private(p,dt)
!$omp do schedule(runtime)
    do p=1,nbody
       dt = dtmax*float(istep - istep0(p))/float(maxstep)
#ifdef Helium
       xyzhm(k,p) = xyzhmp(k,p) + 0.5*vxyzutp(k,p)*dt
#else
       xyzhm(k,p) = xyzhmp(k,p) + vxyzutp(k,p)*dt + 0.5*axyzutp(k,p)*dt*dt
#endif
    enddo
!$omp end do
!$omp end parallel
  enddo 
!
!--Evolve velocities
!
  do k=1,ndim
!$omp parallel default(none) shared(nbody,vxyzut,axyzutp,k) &
!$omp shared(vxyzutp,dtmax,istep,istep0,step) private(p,dt)
!$omp do schedule(runtime)
     do p=1,nbody
        dt = dtmax*float(istep - istep0(p))/float(maxstep)
#ifdef Helium
        vxyzut(k,p) = vxyzutp(k,p) + 0.5*axyzutp(k,p)*dt
#else
        vxyzut(k,p) = vxyzutp(k,p) + axyzutp(k,p)*dt
#endif
     enddo
!$omp end do
!$omp end parallel
  enddo
!
!--Evolve thermal energies and temperatures, if necessary
!
!$omp parallel default(none) shared(nbody,vxyzut,axyzutp,SIMTYPE) &
!$omp shared(enucp,luminucp,RELFLAG,vxyzutp,dtmax,istep,istep0,step) private(p,dt)
!$omp do schedule(runtime)
  do p=1,nbody
     dt = dtmax*float(istep - istep0(p))/float(maxstep)
#ifdef Helium
     vxyzut(4,p) = vxyzutp(4,p) + 0.5*axyzutp(4,p)*dt + enucp(p)
#else
     if ((vxyzut(5,p) > tmin).and.(vxyzut(5,p) < tmax)) then
        vxyzut(4,p) = vxyzutp(4,p) + axyzutp(4,p)*dt + enucp(p)
     endif
#endif
  enddo
!$omp end do

  if (RELFLAG.eqv..false.) then
!$omp do schedule(runtime)
     do p=1,nbody
        dt = dtmax*float(istep - istep0(p))/float(maxstep)
#ifdef Helium
        vxyzut(5,p) = vxyzutp(5,p) + 0.5*axyzutp(5,p)*dt + luminucp(p)
#else
        vxyzut(5,p) = vxyzutp(5,p)  + axyzutp(5,p)*dt + luminucp(p)
#endif
     enddo
!$omp end do
  endif
!$omp end parallel
!
!--Save eps for late use
!
  eps = 1.0d30
!$omp parallel default(none) shared(nbody,xyzhm) private(p)  &
!$omp reduction(MIN:eps)
!$omp do schedule(runtime)
  do p=1,nbody
     eps  = MIN(eps,xyzhm(4,p))
  enddo
!$omp end do
!$omp end parallel
  eps  = 1.4d0*2.0d0*eps
  eps3 = eps*eps*eps
!
!--Correct He layer velocity
!
!#ifdef Helium
  call norm_layer
  call layer
!#endif
!
end subroutine predictor
