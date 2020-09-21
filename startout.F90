subroutine startout
!===========================================================================
!
!     This subroutine calls all the necessary subroutines in order to
!     start the code
!
!     Last revision: 15/March/2015
!
!===========================================================================
!
!--Load modules
!
  use mod_parameters, only : MASTER
  use mod_commons,    only : rank, nstep_infile, SIMTYPE, relflag,  &
                             nprocs, nbody1, nbody2, tnow
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--Local variables
!
#ifdef openmp
  integer :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
#endif
!
!--Read code parameters needed for startup
!
  call inparams
#ifdef debug
  if (rank == MASTER) print*, 'imparams called'
#endif
!
!--Read SPH initial data
!
  call indata(0)
#ifdef debug
  if (rank == MASTER) print*, 'indata called'
#endif
!
!--Print simulation summary
!
  if (rank == MASTER) then
     print*, '====================================================='
     print*, '            SIMULATION INFORMATION'
     print*, ''
     write(*,'(a,i4.4,a)')'Starting from bodi',nstep_infile,'.out'
     write(*,'(a,1ES12.4)')'Starting from time=',tnow
     write(*,'(a,i1,a,L1)') 'Number of WDs = ',SIMTYPE
     write(*,'(a,L1)') 'Relaxation = ', RELFLAG
#ifdef MPI
     write(*,'(a,i3)') 'Number of MPI processes = ',nprocs
#endif
#ifdef openmp
!$OMP PARALLEL
     if (OMP_GET_THREAD_NUM() == 0) write(*,'(a,i3)')               &
        'Number of OpenMP processes = ', OMP_GET_NUM_THREADS()
!$OMP end PARALLEL
#endif
#ifndef global
     write(*,'(a)') 'With individual time-steps'
#endif
     write(*,'(a,i6)') 'Number of particles wd #1 = ', nbody1
     write(*,'(a,i6)') 'Number of particles wd #2 = ', nbody2
     print*,'======================================================'
  endif

end subroutine startout
