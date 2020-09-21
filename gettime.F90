       SUBROUTINE gettime(time)
!===================================================================
!
!  This subroutine gives the present time
!
!  Last revision: 15/March/2015
!
!===================================================================
!
!--Force to declare EVERYTHING
!
       IMPLICIT NONE
!
!--I/O definitions
!
       REAL :: time
!
!--Local definitions
!
#ifdef openmp 
      INTEGER :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
      REAL :: OMP_GET_WTIME
#endif
!
#ifdef openmp
!$OMP PARALLEL
       IF (OMP_GET_THREAD_NUM().EQ.0) time = OMP_GET_WTIME()
!$OMP END PARALLEL
#else
       CALL cpu_time(time)
#endif
!
       END SUBROUTINE gettime    
