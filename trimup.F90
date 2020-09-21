subroutine trimup
!===========================================================================
!
! This subroutine eliminates the CO particles in the He layer
!
! Last revision: 09/April/2019
!
!===========================================================================
!
!--Load modules
!
  use mod_commons, only : xyzhm, xyzhmp, nbody, star, partype, nb,     &
                          rank, factor
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
  integer :: m, p, plocal, q, nHe, nCO, ierr
  logical, dimension(nbody)  :: change
  logical, dimension(factor) :: tobechanged
#ifdef MPI
  logicaL, dimension(nbody) :: sentarray, receivearray
#endif
!
!--Now remove the CO particles from the He layer of the stars
!
  tobechanged(1:factor) = .false.
!$omp parallel default(none) shared(nbody,nb,partype,xyzhmp,rank,factor) &
!$omp shared(tobechanged) private(p,plocal,q,m,nCO,nHe)
!$omp do schedule(runtime)
  do plocal=1,factor
    p = plocal + rank*factor
    if (partype(p) /= 0) cycle
    nCO = 0
    nHe = 0
    do m=1,nb(1,plocal)
      q = nb(m+1,plocal)
      if (partype(q) == 0) nCO = nCO + 1
        if (partype(q) == 1) nHe = nHe + 1
    enddo

    if (nHe > 10*nCO) then
      tobechanged(plocal) = .true.
    endif 
  enddo
!$omp end do
!$omp end parallel
!
!--Force synchronization
!
#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
!--Transfer tobetranfered
! 
#ifdef MPI
!$omp parallel default(none) shared(factor,sentarray,tobechanged) &
!$omp private(p)
!$omp do schedule(runtime)
  do p = 1, factor
    sentarray(p) = tobechanged(p)
  enddo
!$omp end do
!$omp end parallel
!
  call MPI_ALLGATHER(sentarray,factor,MPI_LOGICAL,receivearray,     &
                     factor,MPI_LOGICAL,MPI_COMM_WORLD,ierr)
!
!$omp parallel default(none) shared(nbody,receivearray,change) &
!$omp private(p)
!$omp do schedule(runtime)
  do p = 1, nbody
    change(p) = receivearray(p)
  enddo
!$omp end do
!$omp end parallel
#else
!$omp parallel default(none) shared(nbody,tobechanged,change) &
!$omp private(p)
!$omp do schedule(runtime)
  do p = 1, nbody
    change(p) = tobechanged(p)
  enddo
!$omp end do
!$omp end parallel
#endif MPI
!
!--Now change the particles that need to be changed
!
  do p=1,nbody
    if (change(p) .eqv. .true.) THEN
      partype(p) = 2
    endif
  enddo
!
end subroutine
