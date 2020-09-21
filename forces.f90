subroutine forces
!===================================================================
!
!  This subroutine adds up the hydro and gravitational forces.
!  It also adds coriolis and dumping forces in a relaxation.
!
!  Last revision: 6/April/2019  
!
!===================================================================
!
!--Locad modules
!
  use mod_parameters, only : ndim
  use mod_commons,    only : axyzut, gxyzu, xyzhm, vxyzut, rotforc, &
                             nbody, rotforc, SIMTYPE, trelax,       &
                             istep0, istep, step, nstep
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--Local variables
!
  integer :: k, p, m
!
!--Calculate non-inertial forces
!
  select case (SIMTYPE) 
    case (1)
      rotforc = 0.0
    case (2)
      call noninertial
    case default
      print*, 'Bad SIMTYPE !!!', SIMTYPE
      stop
  end select
!
!--Add gravitational, hydrodynamical, and inertial forces
!
  do k=1,ndim
!$omp parallel default(none) shared(nbody,axyzut,gxyzu,k,trelax,vxyzut) &
!$omp shared(rotforc,istep0,istep,step,nstep) private(p)
!$omp do schedule(runtime)
    do p=1,nbody
      if ((nstep == 0) .or. (istep0(p) + step(p) == istep)) then
        axyzut(k,p) = axyzut(k,p) + gxyzu(k,p) - trelax*vxyzut(k,p) + rotforc(k,p)
      endif
    enddo
!$omp end do
!$omp end parallel
  enddo
!
!--EXPERIMENTAL: I will try to add a dissipation into the temperature "acceleration"
!
!$omp parallel default(none) shared(nbody,axyzut,gxyzu,k,trelax,vxyzut) &
!$omp shared(rotforc,istep0,istep,step,nstep) private(p)
!$omp do schedule(runtime)
  do p=1,nbody
    if ((nstep == 0) .or. (istep0(p) + step(p) == istep)) then
      axyzut(5,p) = axyzut(5,p) - trelax*(vxyzut(5,p)-1.0e7)
    endif
  enddo
!$omp end do
!$omp end parallel
!
end subroutine forces
