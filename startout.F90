      SUBROUTINE startout
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
      USE mod_parameters, ONLY : MASTER
      USE mod_commons, ONLY : rank, nstep_infile, SIMTYPE, relflag, nprocs,&
                              nbody1, nbody2
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local variables
!
#ifdef openmp
      INTEGER :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
#endif
!
!--Read code parameters needed for startup
!
      CALL inparams
#ifdef debug
      IF (rank == MASTER) PRINT*, 'imparams called'
#endif
!
!--Read SPH initial data
!
      CALL indata
#ifdef debug
      IF (rank == MASTER) PRINT*, 'indata called'
#endif
!
!--Print simulation summary
!
      IF (rank == MASTER) THEN
         PRINT*, '====================================================='
         PRINT*, '            SIMULATION INFORMATION'
         PRINT*, ''
         WRITE(*,'(a,i4.4,a)')'Starting from bodi',nstep_infile,'.out'
         WRITE(*,'(a,i1,a,L1)') 'Number of WDs = ',SIMTYPE
         WRITE(*,'(a,L1)') 'Relaxation = ', RELFLAG
#ifdef MPI
         WRITE(*,'(a,i3)') 'Number of MPI processes = ',nprocs
#endif
#ifdef openmp
!$OMP PARALLEL
         IF (OMP_GET_THREAD_NUM() == 0) WRITE(*,'(a,i3)')               &
            'Number of OpenMP processes = ', OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
#endif
         WRITE(*,'(a,i6)') 'Number of particles wd #1 = ', nbody1
         WRITE(*,'(a,i6)') 'Number of particles wd #2 = ', nbody2
         PRINT*,'======================================================'
      ENDIF
!
      END SUBROUTINE startout
