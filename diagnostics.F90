      SUBROUTINE diagnostics(t1,t2)
!===========================================================
!  This subroutine outputs all the necessary information
!
!  Last revision: 15/March/2015
!===========================================================
!
!--Load modules
!
      USE mod_parameters, ONLY : MASTER
      USE mod_commons, ONLY : globnmin, globnmax, globnvec, globmax,    &
      globdone, nstep, nout, rank, nCO, nHe, ndead
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--I/O variables
!
      REAL, INTENT(IN) :: t1, t2
!
      IF (rank.EQ.MASTER) THEN
         CALL energy
!
         IF (MOD(nstep,nout).EQ.0) THEN
            CALL outdata
!
            PRINT*, ''
            PRINT*, '============================================'
            WRITE(*,'(a,i4.4,a)') 'bodi',nstep/nout+1,'.out written'
            PRINT*, ''
            WRITE(*,'(a,i6,a,i6,a,i6)') 'NMIN:',globnmin, ' NMAX:',     &
                  globnmax,' AVG:',globnvec
            WRITE(*,'(a,i7)') 'Notdone=',globdone
            WRITE(*,'(a,1pe12.4,a)') 'Maxdiff',globmax*100.,' %'
            PRINT*, ''
            WRITE(*,'(a,1pe12.4)') 'Average iteration time=',           &
                  (t2-t1)/nout
            WRITE(*,'(a,i6,a,i6,a,i6)') 'nCO=',nCO,' nHe=',nHe,         &
                  ' ndead=',ndead
            PRINT*, '============================================'
            PRINT*, ''
         ENDIF
      ENDIF
!
      END SUBROUTINE diagnostics
