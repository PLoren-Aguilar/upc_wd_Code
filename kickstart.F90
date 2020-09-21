      SUBROUTINE kickstart
!===================================================================
!  This subroutine call the necessary subroutines to start the
!  integration process
!
!  Last revision: 15/March/2015
!
!===================================================================
!
!--Load modules
!
      USE mod_parameters, ONLY : MASTER, ndim
      USE mod_commons, ONLY : axyzut, axyzutp, enuc, enucp, luminuc,    &
      luminucp, dhdt, dhdtp, xyzhm, xyzhmp, vxyzut, vxyzutp, nbody,     &
      rank, nstep
! 
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local definitions
!
      INTEGER :: k, m, p
!
!--Calculation of density, smoothing lenght and tree
!
      CALL iter_rhoh
#ifdef debug
      IF (rank.EQ.MASTER) PRINT*, 'iter_rhoh  called'
#endif
!
!--Calculate EOS
!
      IF (nstep == 1) THEN
         CALL EOS0
      ELSE
         CALL EOS
      ENDIF
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
!--Variables initialization
!
!$OMP PARALLEL DEFAULT(none) shared(axyzutp,axyzut,enucp,enuc)        &
!$OMP shared(luminucp,luminuc,dhdtp,dhdt,xyzhmp,xyzhm,vxyzutp,vxyzut) &
!$OMP shared(nbody) private(m,p,k)
!$OMP DO SCHEDULE(runtime)
      DO p=1,nbody
         DO k=1,ndim+2
            xyzhmp(k,p)  = xyzhm(k,p)
            vxyzutp(k,p) = vxyzut(k,p)
            axyzutp(k,p) = axyzut(k,p)
         ENDDO
         enucp(p)    = enuc(p)
         luminucp(p) = luminuc(p)
         dhdtp(p)    = dhdt(p)
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
!--First time-steps calculation
!
      CALL varydt
#ifdef debug
      IF (rank.EQ.MASTER) PRINT*, 'varydt called'
#endif
!
      END SUBROUTINE kickstart
