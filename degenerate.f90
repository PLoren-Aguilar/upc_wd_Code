      SUBROUTINE degenerate(p)
!========================================================================
!     THIS   SUBROUTINE   ALLOWS   TO  SELECT   BETWEEN  THE DIFFERENT  
!     EQUATIONS  OF  STATE  OF  THE  CODE. IT'S IMPORTANT  TO  REALIZE  
!     THAT  EVERY  EOS  SUBROUTINE MUST PROVIDE (AT LEAST):
!
!     * TEMPERATURE, PRESSURE, Cv AND SPEED OF SOUND
!
!     Last revision: 15/March/2015!
!=========================================================================
!
!--Load modules
!
      USE mod_parameters, ONLY : uden, unm, uen, up, uv, unl, tmin, tmax,&
      ndim, nel, MASTER
      USE mod_commons, ONLY : rho, xss, aion, zion, vxyzut, press, css, &
      cvs, dPdT, cps, rank
      USE mod_EOS
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Helmholtz EOS definitions
!
!      INCLUDE 'vector_eos.dek'
!
!--I/O variables
!
      INTEGER, INTENT(IN) :: p
!
!--Local variables
!
!      REAL(8), DIMENSION(1) :: ewant_row
      REAL(8)  :: abar, zbar
      INTEGER :: k
!
!--Set initial data for the EOS
!
      IF (vxyzut(5,p).GT.tmax) vxyzut(5,p) = tmax
      IF (vxyzut(5,p).LT.tmin) vxyzut(5,p) = tmin
!
      temp_row(1)  = vxyzut(5,p)
      den_row(1)   = rho(p)*uden    ! EOS works in cgs units!!!
!
      abar = 0.0
      zbar = 0.0
      DO k=1,nel-1
         abar = abar + xss(k,p)/aion(k)
         zbar = zbar + xss(k,p)*zion(k)/aion(k)
      ENDDO
      abar_row(1) = 1.0/abar
      zbar_row(1) = zbar/abar
!
!--Set the number of elements sent into the EOS
!
      jlo_eos=1
      jhi_eos=1
!
!--Call EOS to obtain a first estimate for energy, pressure and derivatives
!
      CALL helmeos
      IF (rank.EQ.MASTER.AND.eosfail.EQV..true.) PRINT*,'EOScorr fail'
!
      vxyzut(4,p) = etot_row(1)*(unm/uen)
      press(p) = ptot_row(1)/up
      css(p)   = cs_row(1)/uv
      cvs(p)   = cv_row(1)*(unm/uen)
      dPdT(p)  = dpt_row(1)*(unl**3/uen)
      cps(p)   = cp_row(1)/(uen/unm)
!
      END SUBROUTINE degenerate
