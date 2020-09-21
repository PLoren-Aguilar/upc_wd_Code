      SUBROUTINE EOS0
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
      USE mod_parameters, ONLY : uden, unm, uen, up, uv, unl, tmin,     &
      tmax, nel, ndim, MASTER
      USE mod_commons, ONLY : rho, xss, aion, zion, vxyzut, press, css, &
      cvs, dPdT, cps, error, rank, nbody
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Helmholtz EOS definitions
!
      INCLUDE 'vector_eos.dek'
!
!--Local variables
!
      REAL, DIMENSION(nrowmax) :: ewant_row, temp_row, den_row, abar_row, &
      zbar_row, det_row, ptot_row, cs_row, cv_row, dpt_row, cp_row, etot_row
      REAL    :: abar, zbar
      INTEGER :: p, k, jlo_eos, jhi_eos
      LOGICAL :: eosfail = .false.
!
!--Set initial data for the EOS
!
!$OMP PARALLEL default(none) shared(nbody,vxyzut,rho,xss,aion,zion)   &
!$OMP shared(rank,press,css,cvs,dPdT,cps) private(p,temp_row,den_row) &
!$OMP private(abar,zbar,abar_row,zbar_row,jlo_eos,jhi_eos,eosfail)    &
!$OMP private(etot_row,ptot_row,cs_row,cv_row,dpt_row,cp_row,det_row)
!$OMP DO SCHEDULE(runtime)
      partloop : DO p = 1,nbody
        IF (vxyzut(5,p) > tmax) vxyzut(5,p) = tmax
        IF (vxyzut(5,p) < tmin) vxyzut(5,p) = tmin
!
        temp_row(1)  = vxyzut(5,p)
        den_row(1)   = rho(p)*uden    ! EOS works in cgs units!!!
!
        abar = 0.0
        zbar = 0.0
        DO k=1,nel-1
           abar = abar + xss(k,1)/aion(k)
           zbar = zbar + xss(k,1)*zion(k)/aion(k)
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
!$OMP CRITICAL
        CALL helmeos(jlo_eos,jhi_eos,temp_row,den_row,abar_row,zbar_row, &
            det_row,ptot_row,cs_row,cv_row,dpt_row,cp_row,etot_row,eosfail)
!$OMP END CRITICAL
        IF ((rank == MASTER)  .AND. (eosfail.EQV..true.)) THEN
           PRINT*,'EOScorr fail'
           STOP
        ENDIF
!
        vxyzut(4,p) = etot_row(1)*(unm/uen)
        press(p) = ptot_row(1)/up
        css(p)   = cs_row(1)/uv
        cvs(p)   = cv_row(1)*(unm/uen)
        dPdT(p)  = dpt_row(1)*(unl**3/uen)
        cps(p)   = cp_row(1)/(uen/unm)
      ENDDO partloop
!$OMP END DO
!$OMP END PARALLEL
!
      END SUBROUTINE EOS0
