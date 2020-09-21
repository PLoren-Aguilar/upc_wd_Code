subroutine approach(frac,finish)
!===========================================================================
!
!      This subroutine approaches the stars while keeping them correctly
!      positioned around the com
!
!      Last revision: 02/March/2017
!
!===========================================================================
!
!--Load modules
!
  use mod_parameters, only : ndim, MASTER
  use mod_commons,    only : xyzhmp, nbody, star, rank, partype, Omega0, &
                             fix, fixx1, fixx2, nstep_infile
!
!--Force to declare everything
!
  implicit none
!
!--I/O variables
!
  real,    intent(in)  :: frac
  logical, intent(out) :: finish
!
!--Local variables
!
  real, dimension(ndim) :: cmp1, cmp2
  real    :: mtot1, mtot2, m, d1, d2, gravpot1, gravpot2, dist,    &
             x1, x2, rad2, rx, ry, rz, cmd
  integer :: p, k
!
!--Loop over particles to calculate the center of mass positions for the stars
!
  cmp1  = 0.0d0
  cmp2  = 0.0d0
  mtot1 = 0.0d0
  mtot2 = 0.0d0
!$omp parallel default(none) shared(xyzhmp,nbody,star) private(p,k) &
!$omp reduction(+:cmp1) reduction(+:cmp2) reduction(+:mtot1) reduction(+:mtot2)
!$omp do schedule(runtime)
  do p=1,nbody
    if (star(p) == 1) then
      mtot1 = mtot1 + xyzhmp(5,p)
      do k=1,ndim
        cmp1(k) = cmp1(k) + xyzhmp(5,p)*xyzhmp(k,p)
      enddo
    elseif (star(p) == 2) then
      mtot2 = mtot2 + xyzhmp(5,p)
      do k=1,ndim
        cmp2(k) = cmp2(k) + xyzhmp(5,p)*xyzhmp(k,p)
      enddo
    endif
  enddo
!$omp end do
!$omp end parallel
  cmp1 = cmp1/mtot1
  cmp2 = cmp2/mtot2
!
!--Reduce/reposition the stellar distance by the desired amount
!
!  dist  = (1.0-frac)*(ABS(fixx1)+ABS(fixx2))
!  fixx1 = -mtot2*dist/(mtot2+mtot1)
!  fixx2 = dist-abs(fixx1)
!  if (rank == MASTER) then
!    PRINT*, 'Old positions = ',cmp1(1),cmp2(1)
!    PRINT*, 'New positions = ',fixx1,fixx2
!  endif
!
!--Reposition both stars around the CM given the new distance
!
!$omp parallel default(none) shared(nbody,star,xyzhmp,cmp1,cmp2,x1,x2) &
!$omp shared(fixx1,fixx2) private(p)
!$omp do schedule(runtime)
  do p=1,nbody
    if (star(p) == 1) then
      xyzhmp(1,p) = xyzhmp(1,p) - cmp1(1) + fixx1
      xyzhmp(2,p) = xyzhmp(2,p) - cmp1(2)
      xyzhmp(3,p) = xyzhmp(3,p) - cmp1(3)
    elseif (star(p) == 2) then
      xyzhmp(1,p) = xyzhmp(1,p) - cmp2(1) + fixx2
      xyzhmp(2,p) = xyzhmp(2,p) - cmp2(2)
      xyzhmp(3,p) = xyzhmp(3,p) - cmp2(3)
    endif
  enddo
!$omp end do
!$omp end parallel
!
!--Calculate cm position for both stars after repositioning
!
  cmp1  = 0.0d0
  cmp2  = 0.0d0
  mtot1 = 0.0d0
  mtot2 = 0.0d0
!$omp parallel default(none) shared(xyzhmp,nbody,star) private(p,k) &
!$omp reduction(+:cmp1) reduction(+:cmp2) reduction(+:mtot1) reduction(+:mtot2)
!$omp do schedule(runtime)
  do p=1,nbody
    if (star(p) == 1) then
      mtot1 = mtot1 + xyzhmp(5,p)
      do k=1,ndim
        cmp1(k) = cmp1(k) + xyzhmp(5,p)*xyzhmp(k,p)
      enddo
    elseif (star(p) == 2) then
      mtot2 = mtot2 + xyzhmp(5,p)
      do k=1,ndim
        cmp2(k) = cmp2(k) + xyzhmp(5,p)*xyzhmp(k,p)
      enddo
    endif
  enddo
!$omp end do
!$omp end parallel
  cmp1 = cmp1/mtot1
  cmp2 = cmp2/mtot2
!
!--Now remove the CO particles from the He layer of the stars
!
  do p=1,nbody
    if ((star(p) == 1) .and. (partype(p) == 0)) then
      rad2 = (xyzhmp(1,p)-cmp1(1))**2 + (xyzhmp(2,p)-cmp1(2))**2 +   &
             (xyzhmp(3,p)-cmp1(3))**2
      if (rad2 > 0.09*0.09) then
        rx = (xyzhmp(1,p)-cmp1(1))/sqrt(rad2)
        ry = (xyzhmp(2,p)-cmp1(2))/sqrt(rad2)
        rz = (xyzhmp(3,p)-cmp1(3))/sqrt(rad2)
!
        xyzhmp(1,p) = xyzhmp(1,p) - (0.09/3.)*rx
        xyzhmp(2,p) = xyzhmp(2,p) - (0.09/3.)*ry
        xyzhmp(3,p) = xyzhmp(3,p) - (0.09/3.)*rz
      endif
    elseif ((star(p) == 2) .and. (partype(p) == 0)) then
      rad2 = (xyzhmp(1,p)-cmp2(1))**2 + (xyzhmp(2,p)-cmp2(2))**2 +   &
             (xyzhmp(3,p)-cmp2(3))**2
      if (rad2 > 0.06*0.06) then
        rx = (xyzhmp(1,p)-cmp2(1))/sqrt(rad2)
        ry = (xyzhmp(2,p)-cmp2(2))/sqrt(rad2)
        rz = (xyzhmp(3,p)-cmp2(3))/sqrt(rad2)
!
        xyzhmp(1,p) = xyzhmp(1,p) - (0.06/3.)*rx
        xyzhmp(2,p) = xyzhmp(2,p) - (0.06/3.)*ry
        xyzhmp(3,p) = xyzhmp(3,p) - (0.06/3.)*rz
      endif
    endif
  enddo
!
!--Calculate gravitational potentials and check if there is RL overflow
!
  finish = .false.
  !do p=1,nbody1
  !   gravpot1 = mtot1/SQRT((xyzhm(1,p)-cmp1(1))**2 +              &
  !                (xyzhm(2,p)-cmp1(2))**2 +                        &
  !                (xyzhm(3,p)-cmp1(3))**2)
  !   gravpot2 = mtot2/SQRT((xyzhm(1,p)-cmp2(1))**2 +              &
  !                (xyzhm(2,p)-cmp2(2))**2 +                        &
  !                (xyzhm(3,p)-cmp2(3))**2)

  !  if ((gravpot1-gravpot2).LE.0.0.and.ABS(gravpot1-gravpot2).GT.5d-3) then
  !      finish=.true.
  !   endif
  !enddo
end subroutine
