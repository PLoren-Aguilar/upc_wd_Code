         SUBROUTINE predcorr
!===================================================================
!  This subroutine evolves the system a single step using an 
!  Adams-Bashford predictor-corrector integrator:
!
!  Serna A., Alimi J.-M, Chieze J.-P, 1995, Astrophysical Journal,
!  461, 884
!
!  Last revision: 15/March/2015
!===================================================================
!
!--Load modules
!
         USE mod_parameters, ONLY : MASTER
         USE mod_commons, ONLY : rank, RELFLAG
!
!--Force to declare EVERYTHING
!
         IMPLICIT NONE
!
!--Calculation of predictor step
!
         CALL predictor
#ifdef debug
         IF (rank.EQ.MASTER) PRINT*, 'predictor called'
#endif
!
!--Calculation of density, smoothing lenght and tree+gravitational forces
!
         CALL iter_rhoh
#ifdef debug
         IF (rank.EQ.MASTER) PRINT*, 'iter_rhoh called'
#endif
!
!--Calculate EOS
!
         CALL EOS
#ifdef debug
         IF (rank.EQ.MASTER) PRINT*, 'EOS called'
#endif
!
!--Calculation of hydrodynamical quantities
!
         CALL hydro_rs
#ifdef debug
         IF (rank.EQ.MASTER) PRINT*, 'hydro_rs called'
#endif
!
!--Sum forces. Add imaginary forces if necessary
!
         CALL forces
#ifdef debug
         IF (rank.EQ.MASTER) PRINT*, 'forces called'
#endif
!
!--Calculation of nuclear burning

         IF (RELFLAG.EQV..false.) CALL burn
#ifdef debug
         IF (rank.EQ.MASTER) PRINT*, 'burn called'
#endif
!
!--Correction calculation
!
         CALL corrector
!
!--Sort particles to avoid killed ones
!
         CALL sort
#ifdef debug
      IF (rank.EQ.MASTER) PRINT*, 'sort called'
#endif
!
!--Next time-step calculation
!
         CALL varydt
#ifdef debug
         IF (rank.EQ.MASTER) PRINT*, 'varydt called'
#endif
!
         END SUBROUTINE predcorr
