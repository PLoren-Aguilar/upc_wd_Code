subroutine relax_frame
!===========================================================================
! This subroutine substracts cm velocity in order to correctly relax
!
! Last revision: 95/April/2019
!===========================================================================
!
!--Load modules
!
  use mod_parameters, only : ndim
  use mod_commons,    only : xyzhm, vxyzut, nbody
!
!--Force to declare EVERYTHING
!
  IMPLICIT NONE
!
!--Local variables
!
  real, dimension(3) :: cmp(3),cmv(3)
  real    :: mtot, mass
  integer :: k, m, p
!
!--Loop over particles to calculate center of mass velocity and position of the system
!
  mtot = 0.0
  cmp  = 0.0
  cmv  = 0.0
!$omp parallel default(none) shared(xyzhm,nbody) private(m,p) &
!$omp reduction(+:mtot)
!$omp do schedule(runtime)
  do p=1,nbody
     mtot = mtot + xyzhm(5,p)
  enddo
!$omp end do
!$omp end parallel
!
  do k=1,ndim
!$omp parallel default(none) shared(k,xyzhm,vxyzut,nbody) private(p,mass) &
!$omp reduction(+:cmp) reduction(+:cmv)
!$omp do schedule(runtime)
     do p=1,nbody
        mass = xyzhm(5,p)
        cmp(k) = cmp(k) + mass*xyzhm(k,p)
        cmv(k) = cmv(k) + mass*vxyzut(k,p)
     enddo
!$omp end do
!$omp end parallel
  enddo
  cmp = cmp/mtot
  cmv = cmv/mtot
!
!--Remove any possible center-of-mass offset
!
  do k=1,ndim
!$omp parallel default(none) shared(nbody,k,xyzhm,vxyzut,cmp,cmv) private(m,p)
!$omp do schedule(runtime)
     do p=1,nbody
        xyzhm(k,p)  = xyzhm(k,p)  - cmp(k)
        vxyzut(k,p) = vxyzut(k,p) - cmv(k)
     enddo
!$omp end do
!$omp end parallel
  enddo

end subroutine relax_frame
