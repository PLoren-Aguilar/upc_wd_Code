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
      cvs, dPdT, cps, rank, nbody, partype
      USE mod_EOS
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Helmholtz EOS definitions
!
      !INCLUDE 'vector_eos.dek'
!
!--Local variables
!
      !REAL, DIMENSION(nbody) :: ewant_row, temp_row, den_row, abar_row, &
      !zbar_row, det_row, ptot_row, cs_row, cv_row, dpt_row, cp_row, etot_row
      !INTEGER :: jlo_eos, jhi_eos
      !LOGICAL :: eosfail = .false.
      REAL    :: abar, zbar
      INTEGER :: p, k, nactive
!
!--Set initial data for the EOS
!
!$OMP PARALLEL default(none) shared(nbody,vxyzut,rho,xss,aion,zion)   &
!$OMP shared(rank,press,css,cvs,dPdT,cps,temp_row,den_row,partype)    & 
!$OMP shared(jlo_eos,jhi_eos,eosfail,etot_row,ptot_row,cs_row)        &
!$OMP shared(cv_row,dpt_row,cp_row,det_row,abar_row,zbar_row,eoslist) &
!$OMP private(p,k,abar,zbar)
!$OMP DO SCHEDULE(runtime)
      partloop : DO p = 1,nbody
        IF (partype(p) == 2) CYCLE
        IF (vxyzut(5,p) > tmax) vxyzut(5,p) = tmax
        IF (vxyzut(5,p) < tmin) vxyzut(5,p) = tmin
!
        temp_row(p)  = vxyzut(5,p)
        den_row(p)   = rho(p)*uden    ! EOS works in cgs units!!!
!
        abar = 0.0
        zbar = 0.0
        DO k=1,nel-1
           abar = abar + xss(k,p)/aion(k)
           zbar = zbar + xss(k,p)*zion(k)/aion(k)
        ENDDO
        abar_row(p) = 1.0/abar
        zbar_row(p) = zbar/abar
      ENDDO partloop
!$OMP END DO
!$OMP END PARALLEL
!
!--Call EOS to obtain a first estimate for energy, pressure and
!derivatives
!

!--Construct the new list of active particles
       nactive = 0
       DO p = 1,nbody
         IF (partype(p) /= 2) THEN
            nactive = nactive + 1
            eoslist(nactive) = p
         ENDIF
      ENDDO
      jlo_eos = 1
      jhi_eos = nactive
      CALL helmeos
      PRINT*, 'Helmeos called'
      IF ((rank == MASTER)  .AND. (eosfail.EQV..true.)) THEN
         PRINT*,'EOScorr fail'
         STOP
      ENDIF
!
!$OMP PARALLEL default(none) shared(nbody,vxyzut,press,css,cvs,dPdT,cps)  &
!$OMP shared(etot_row,ptot_row,cs_row,cv_row,dpt_row,cp_row,partype) private(p)
!$OMP DO schedule(runtime)
      partloop2 : DO p = 1, nbody
        IF (partype(p) == 2) CYCLE
        vxyzut(4,p) = etot_row(p)*(unm/uen)
        press(p) = ptot_row(p)/up
        css(p)   = cs_row(p)/uv
        cvs(p)   = cv_row(p)*(unm/uen)
        dPdT(p)  = dpt_row(p)*(unl**3/uen)
        cps(p)   = cp_row(p)/(uen/unm)
      ENDDO partloop2
!$OMP END DO
!$OMP END PARALLEL
!
      END SUBROUTINE EOS0
