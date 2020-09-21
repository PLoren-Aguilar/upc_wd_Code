subroutine norm_layer
!============================================================
!  This subroutine calculates the hydrodynamical 
!  accelerations
!
!  Last revision: 9/April/2019
!============================================================
!
!--Load modules
!
  use mod_parameters, only : ndim
  use mod_commons,    only : xyzhm, vxyzut, rho, partype, nb, norm, &
                             rhoG, rank, factor, nbody, ierr, npart,&
                             istep, istep0, step, nstep
  use mod_functions,  only : wk
!
!--Force to declare EVERYTHING
!
  implicit none 
#ifdef MPI
  INCLUDE 'mpif.h'
#endif 
!
!--Local variables
!
#ifdef MPI
  real, dimension(2*npart) :: sentarray, receivearray
#endif 
  real, dimension(ndim) :: vxyz
  real :: massp, r2, u2p, hp, scalar, dr, drx, dry, drz, dvx,       &
          dvy, dvz, rhop, rhoq, chi, wkp, DD_wkp
  real, parameter :: norma=10./9.
  integer ::  i, q, qlocal, p, k, m
  integer, parameter :: nu=3
!
!--Calculate normalization
!
!$omp parallel default(none) shared(xyzhm,rho,rank,partype,nb,norm)         & 
!$omp shared(rhoG,factor,nbody,nstep,istep,istep0,step) private(i,p,q,rhoq) &
!$omp private(rhop,hp,massp,u2p,wkp,DD_wkp,qlocal,dr,drx,dry,drz)
!$omp do schedule(runtime)
  partloop : do qlocal = 1, factor
    q = qlocal + rank*factor
!
    if (q > nbody) cycle
    if ((nstep /= 0) .and. (istep0(q) + step(q) /= istep)) cycle
!
!--Go only through active He particles
!
    if (partype(q) /= 1)  cycle
    if (istep0(q) + step(q) /= istep) cycle
!
!--Store for later use q-particle related variables 
!
    rhoq   = rho(q)
!
!--Loop over active CO neighbours
!
    rhoG(q) = 0.0
    norm(q) = 0.0
    neiloop : do i = 1,nb(1,qlocal)
      p = nb(i+1,qlocal)
!
      if (partype(p) /= 0) cycle
!
!--p-particle variables
!
      rhop   = rho(p)
      hp     = xyzhm(4,p)
      massp  = xyzhm(5,p)
!
!--Calculate kernel related variables
!
      drx = xyzhm(1,q) - xyzhm(1,p)
      dry = xyzhm(2,q) - xyzhm(2,p)
      drz = xyzhm(3,q) - xyzhm(3,p)
      dr  = SQRT(drx*drx + dry*dry + drz*drz)
!
      u2p    = dr*dr/(hp*hp)
      wkp = norma*wk(u2p)*u2p/(hp**ndim)
      !wkp    = wk(u2p)/(hp**ndim)
!
!--Calculate normalization factor
!
      norm(qlocal) = norm(qlocal) + (massp/rhop)*wkp
      rhoG(qlocal) = rhoG(qlocal) + massp*wkp
    enddo  neiloop
  enddo  partloop
!$omp end do
!$omp end parallel
!
!--Force synchronization
!
#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
#endif
!
!--Tranfer of variables across MPI processes, if necessary
! 
#ifdef MPI
!$omp parallel default(none) shared(sentarray,norm,rhoG,factor)         &
!$omp private(p)
!$omp do schedule(runtime)
  do p = 1, factor
    sentarray((p-1)*2 + 1) = norm(p) 
    sentarray((p-1)*2 + 2) = rhoG(p) 
  enddo
!$omp end do
!$omp end parallel
!
  call MPI_ALLGATHER(sentarray,2*factor,MPI_doUBLE_PRECISION,       &
                     receivearray,2*factor,MPI_doUBLE_PRECISION,    &
                     MPI_COMM_WORLD,ierr)
!
!$omp parallel default(none) shared(nbody,receivearray,norm,rhoG,nstep,istep0,step,istep) &
!$omp private(p)
!$omp do schedule(runtime)
  do p = 1, nbody
    if ((nstep /= 0) .and. (istep0(p) + step(p) /= istep)) cycle
    norm(p) = receivearray((p-1)*2 + 1)
    rhoG(p) = receivearray((p-1)*2 + 2)
  enddo
!$omp end do
!$omp end parallel
#endif MPI
!
end subroutine norm_layer
