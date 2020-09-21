      SUBROUTINE sph
!========================================================================
!  This subroutine is the main driver of the code. It evolves the system
!  in time using a given integrator method, and prints diagnostics,
!  whenever necessary
!
!  Last revision: 15/March/2015
!========================================================================
!
!--Load modules
!
      USE mod_parameters, ONLY : MASTER, ndim
      USE mod_commons, ONLY : axyzut, axyzutp, enuc, enucp,             &
      luminuc, luminucp, dhdt, dhdtp, dtnuc, vxyzut, xyzhm, xss,        &
      globnmin, globnmax, globnvec, globmax, globdone, rank, size,&
      nprocs, nstep, nout, dtmp_min, tend, tnow, nbody, npart, ndead,   &
      partype
! 
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
#ifdef MPI
      INCLUDE 'mpif.h'
#endif
!
!--Local definitions
!
      REAL, DIMENSION(ndim)  :: cmp1, cmp2
      REAL    :: omega, cm_d, deldis, masa1, masa2, t1, t2, avgt
      INTEGER :: p, k, nwrite, ierr, nCO, nHe
      LOGICAL :: FINISH
!
!--Startout MPI, if necessary
!
#ifdef MPI
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
#else
      size = 1
      rank = MASTER
#endif
      nprocs = size
      !IF (size.NE.nprocs) THEN
      !   PRINT*, 'fatal error: nprocs',nprocs,' != ',size,' np'
      !   STOP
      !ENDIF
!
!--Startout the code reading the necessary data files and parameters
!
      CALL startout
      nCO   = 0
      nHe   = 0
      ndead = 0
      DO p=1,npart
         IF (partype(p) == 0) nCO   = nCO + 1
         IF (partype(p) == 1) nHe   = nHe + 1
         IF (partype(p) == 2) ndead = ndead + 1
      END DO
      nbody = npart - ndead
      PRINT*, 'nCO=',nCO,' nHe=',nHe,' ndead=',ndead

!
!--Kickstart the integrator
!
      CALL kickstart
!
!--Get time
!
      CALL gettime(t1)
!
!--Main time evolution loop
!
      nwrite = nstep
      FINISH = .false.
      timeloop : DO WHILE (FINISH.EQV..false.)
!
         tnow  = tnow + dtmp_min
         nstep = nstep + 1
!
!--Print timestep information
!
         IF (rank == MASTER) THEN
            WRITE(*,'(a,i4,a,i4,a,1pe12.4,a,1pe12.4)') 'step = ',       &
            nstep-nwrite,' out of ',nout,' tnow =',tnow,' dt =',dtmp_min
         ENDIF
!
!--Evolve system a step
!
         CALL predcorr
!
!--Relax system if necessary
!
         CALL relax
!
!--Calculate diagnostics and write outputs
!
         CALL gettime(t2)
         CALL diagnostics(t1,t2)
         CALL gettime(t1)
!
!--Check if simulation has ended
!
         IF (tnow > tend) FINISH = .true.
!
!--Update nstep
!
         IF (nstep-nwrite == nout) nwrite = nstep
      END DO timeloop 
!
!--Ouput data for the last time
!
      IF (nstep-nwrite /= nout) nstep = nwrite + nout
      CALL gettime(t2)
      CALL diagnostics(t1,t2)
!
!--Stop MPI
!
#ifdef MPI
      CALL MPI_FINALIZE(ierr)
#endif
      END SUBROUTINE sph
