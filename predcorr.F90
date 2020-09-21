subroutine predcorr
!===================================================================
!  This subroutine evolves the system a single step using an 
!  Adams-Bashford predictor-corrector integrator:
!
!  Serna A., Alimi J.-M, Chieze J.-P, 1995, Astrophysical Journal,
!  461, 884
!
!  The present version can run using individual time-steps
!
!  Last revision: 21/March/2017
!===================================================================
!
!--Load modules
!
  use mod_parameters, only : MASTER, maxstep, log2
  use mod_commons,    only : rank, RELFLAG, nbody, nactive, istep,              &
      istep, istep0, step, npart, partype, nextstep, SYNC, dtmax, ttot,         &
      tpred, ttree1, ttree2, tirho, tEOS, thydro, tforc, tcorr, tdt, tsort, tnuc
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--Local definitions
!
  real :: dt, t1, t2
  integer, dimension(30) :: nbin
  integer :: i, p, n
!
!--Calculation of predictor-corrector cycle. It will loop until all particles
!  arrive to a synchronised step.
!
  SYNC   = .false.
  ttot   = 0.
  tpred  = 0.
  tirho  = 0.
  tEOS   = 0.
  thydro = 0.
  tforc  = 0.
  tnuc   = 0.
  tcorr  = 0.
  tdt    = 0.
  do while (SYNC .eqv. .false.) 
!
!--Evolve positions, velocities, etc. Predictor phase
!
     call predictor

#ifdef debug
     if (rank == MASTER) print*, 'predictor called'
#endif
!
!--Calculation of density, smoothing length and tree+gravitational forces
!
     call iter_rhoh

#ifdef debug
     if (rank == MASTER) print*, 'iter_rhoh called'
#endif
!
!--Calculate EOS
!
     call EOS

#ifdef debug
     if (rank == MASTER) print*, 'EOS called'
#endif
!
!--Calculation of hydrodynamical quantities
!
     call hydro_rs

#ifdef debug
     if (rank == MASTER) print*, 'hydro_rs called'
#endif
!
!--Sum forces. Add imaginary forces if necessary
!
     call forces

#ifdef debug
     if (rank == MASTER) print*, 'forces called'
#endif
!
!--Calculation of nuclear burning
!
     if (RELFLAG.eqv..false.) then
!        call burn

#ifdef debug
        if (rank == MASTER) print*, 'burn called'
#endif
     endif
!
!--Correction calculation
!
     call corrector

#ifdef debug
     if (rank == MASTER) print*, 'corrector called'
#endif
!
!--Next time-step calculation
!
     call varydt

#ifdef debug
     if (rank == MASTER) print*, 'varydt called'
#endif
!
!--Relax the system if necessary (only if all particles are synchronised)
!
     if ((RELFLAG.eqv..true.) .and. (SYNC.eqv..true.)) then
        call relax
#ifdef debug
        if (rank == MASTER) print*, 'relax called'
#endif
     endif
  enddo
end subroutine predcorr
