       subroutine gettime(time)
!===================================================================
!
!  This subroutine gives the present time
!
!  Last revision: 6/April/2019
!
!===================================================================
!
!--Force to declare everything
!
       implicit none
!
!--I/O definitions
!
       real :: time
!
!--Local definitions
!
#ifdef openmp 
      integer :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
      real :: OMP_GET_WTIME
#endif
!
#ifdef openmp
!$omp parallel
       if (OMP_GET_THREAD_NUM() == 0) time = OMP_GET_WTIME()
!$omp end parallel
#else
       call cpu_time(time)
#endif
!
       end subroutine gettime    
