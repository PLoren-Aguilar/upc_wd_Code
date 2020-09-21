subroutine noninertial
!===================================================================
!  This subroutine calculates the effect of non-inertial
!  forces 
!
!  Last revision: 9/April/2019
!===================================================================
!
!--Load modules
!
  use mod_parameters, only : ndim
  use mod_commons,    only : xyzhm, vxyzut, rotforc, star, Omega01, Omega02, nbody, star, nstep, &
                             nstep_ini, rotforc, dtmax, dtmaxin, nout, nstep, RELFLAG, axyzut,   &
                             gxyzu, Omega0, fixx1, fixx2
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--Local variables
!
  real, dimension(ndim) :: cmp1, fcm1, cmp2, fcm2
  real    :: mass, mtot1, mtot2, cmd, Omega, Omegadot, cmp1mod, cmp2mod, fcm1mod, fcm2mod
  integer :: k, p
!
!--Calculate angular velocity of the system
!
  mtot1 = 0.0
!$omp parallel default(none) shared(xyzhm,star,nbody) private(p) &
!$omp reduction(+:mtot1)
!$omp do schedule(runtime)
  do p = 1,nbody
    if (star(p) == 1) mtot1 = mtot1 + xyzhm(5,p)
  enddo
!$omp end do
!$omp end parallel

  cmp1 = 0.0
  fcm1 = 0.0
  do k=1,ndim
!$omp parallel default(none) shared(k,xyzhm,axyzut,gxyzu,nbody,star) private(p,mass) &
!$omp reduction(+:cmp1) reduction(+:fcm1)
!$omp do schedule(runtime)
    do p=1,nbody
      if (star(p) == 1) then
        mass = xyzhm(5,p)
        cmp1(k) = cmp1(k) + mass*xyzhm(k,p)
        fcm1(k) = fcm1(k) + mass*(axyzut(k,p)+gxyzu(k,p))
      endif
    enddo
!$omp end do
!$omp end parallel
     cmp1(k) = cmp1(k)/mtot1
     fcm1(k) = fcm1(k)/mtot1
  enddo

  mtot2 = 0.0
!$omp parallel default(none) shared(nbody,star,xyzhm) private(p) &
!$omp reduction(+:mtot2) 
!$omp do schedule(runtime)
  do p=1,nbody
    if (star(p) == 2) mtot2 = mtot2 + xyzhm(5,p)
  enddo
!$omp end do
!$omp end parallel

  cmp2 = 0.0
  fcm2 = 0.0
  do k=1,ndim
!$omp parallel default(none) shared(nbody,k,xyzhm,axyzut,gxyzu,star) private(p,mass) &
!$omp reduction(+:cmp2) reduction(+:fcm2)
!$omp do schedule(runtime)
    do p=1,nbody
      if (star(p) == 2) then
        mass    = xyzhm(5,p)
        cmp2(k) = cmp2(k) + mass*xyzhm(k,p)
        fcm2(k) = fcm2(k) + mass*(axyzut(k,p)+gxyzu(k,p))
      endif
    enddo
!$omp end do
!$omp end parallel
    cmp2(k) = cmp2(k)/mtot2
    fcm2(k) = fcm2(k)/mtot2
  enddo

!
!--Centre of mass distance and rotational angular velocity
!
! if (RELFLAG.eqv..true.) then
!   cmp1mod = sqrt(cmp1(1)**2 + cmp1(2)**2 + cmp1(3)**2)
!   fcm1mod = sqrt(fcm1(1)**2 + fcm1(2)**2 + fcm1(3)**2)
!   Omega01 = sqrt(fcm1mod/(mtot1*abs(fixx1)))
!
!   cmp2mod = sqrt(cmp2(1)**2 + cmp2(2)**2 + cmp2(3)**2)
!   fcm2mod = sqrt(fcm2(1)**2 + fcm2(2)**2 + fcm2(3)**2)
!   Omega02 = sqrt(fcm2mod/(mtot2*abs(fixx2)))
!
!   WRITE(*,*) 'Angular velocities are', Omega01, Omega02
!   Omega0  = (Omega01 + Omega02)/2.
! endif
!
!--Calculate Omega0
!
  Omega0 = sqrt((mtot1+mtot2)/(abs(fixx1) + abs(fixx2))**3)
!
!--Setup dtmax
!
  if (dtmax == 0.0) dtmax = dtmaxin*Omega0
!
!--This should be added in case of having a dynamical rotating frame. I
!  haven't been able to make it work. Needs more testing
!
  !Omegadot = (Omega - Omega0)/dt
  Omegadot = 0.0 
!
!--Calculate Coriolis and centrifugal forces
!
!$omp parallel default(none) shared(nbody,rotforc,xyzhm,vxyzut,Omega01,Omega02) &
!$omp shared(Omega0,Omegadot,RELFLAG,star) private(p,Omega)
!$omp do schedule(runtime)
  do p=1,nbody
    if (RELFLAG.eqv..false.) then
      rotforc(1,p) = xyzhm(1,p)*Omega0**2 + 2.*vxyzut(2,p)*Omega0 + xyzhm(2,p)*Omegadot
      rotforc(2,p) = xyzhm(2,p)*Omega0**2 - 2.*vxyzut(1,p)*Omega0 - xyzhm(1,p)*Omegadot
      rotforc(3,p) = 0.0
    elseif (RELFLAG.eqv..true.) then
      rotforc(1,p) = xyzhm(1,p)*Omega0**2 !+ 2.*vxyzut(2,p)*Omega0
      rotforc(2,p) = xyzhm(2,p)*Omega0**2 !- 2.*vxyzut(1,p)*Omega0
      rotforc(3,p) = 0.0
    endif
  enddo
!$omp end do
!$omp end parallel
!
end subroutine noninertial
