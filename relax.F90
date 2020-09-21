      SUBROUTINE relax
!===================================================================
!
!  This subroutine relaxes the system whenever necessary.
!
!  Last revision: 15/March/2015
!
!===================================================================
!
!--Load modules
!
      USE mod_parameters, ONLY : MASTER
      USE mod_commons, ONLY : RELFLAG, rank, SIMTYPE, nstep, nout
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local variables
!
      LOGICAL FINISH
!
!--Choose between single and binary system relaxation
!
      IF (RELFLAG.EQV..true.) THEN
!
!--In case of binary relaxation, slowly approach stars 
!  until Roche-Lobe overflow is achieved
!
         IF (SIMTYPE.EQ.2) THEN
!
!--Note that after removing the distance, the systems needs to
!  relax, so we will only reduce the distance every nout time
!  steps
!
            IF (MOD(nstep,nout).EQ.0) THEN
               CALL approach(0.01,FINISH)
            ENDIF
!
            IF (FINISH.EQV..TRUE.) THEN
               PRINT*, 'Binary relaxation finished'
               STOP
            ENDIF
         ENDIF           ! End of SIMTYPE IF
!
!--Remove any center of mass displacement
!
         CALL relax_frame
      ENDIF           ! End of RELFLAG IF 
#ifdef debug
      IF (rank.EQ.MASTER) PRINT*, 'relaxations called'
#endif
!
      END SUBROUTINE relax
