      SUBROUTINE EOS
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
      USE mod_commons,    ONLY : rho, xss, aion, zion, vxyzut, press,   &
      css, cvs, dPdT, cps, partype, rank, RELFLAG, nbody
      USE mod_EOS
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local variables
!
      REAL, PARAMETER :: eos_tol = 1.0e-8
      REAL, DIMENSION(nbody) :: tnew, error
      REAL :: rel, abar, zbar
      INTEGER  :: p, k, m, newton, eosflag, alldone, nactive, new_nactive 
      INTEGER, PARAMETER :: max_newton = 50
      LOGICAL, DIMENSION(nbody) :: done
!
!--Set initial data for the EOS
!
      nactive     = 0
!$OMP PARALLEL DEFAULT(none) shared(nbody,vxyzut,xss,aion,zion,rho) & 
!$OMP shared(partype,ewant_row,temp_row,den_row,abar_row,zbar_row,done)  &
!$OMP shared(eoslist) private(p,k,abar,zbar) reduction(+:nactive)
!$OMP DO SCHEDULE(runtime)
      partloop : DO p = 1, nbody
        IF (partype(p) == 2) CYCLE
        IF (vxyzut(5,p) > tmax) vxyzut(5,p) = tmax
        IF (vxyzut(5,p) < tmin) vxyzut(5,p) = tmin
!
        ewant_row(p) = vxyzut(4,p)*(uen/unm) ! EOS works in cgs units!!!
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
!
        done(p) = .false.
        eoslist(p) = p
        nactive = nactive + 1
      ENDDO partloop
!$OMP END DO
!$OMP END PARALLEL
!
!--Call EOS to obtain a first estimate for energy, pressure and derivatives
!
      jlo_eos = 1
      jhi_eos = nactive
      eosfail = .false.
      CALL helmeos
      IF ((rank == MASTER) .AND. (eosfail.EQV..true.)) THEN
         PRINT*,'EOScorr fail'
         STOP
      ENDIF
      IF (RELFLAG.EQV..true.) RETURN
!
!--Now, loop over NR iterations untill all particles are done or the maximum
!  number of iterations is reached
!
      newton  = 0
      alldone = 0
      NRloop : DO
        newton = newton + 1
        IF ((newton > max_newton) .OR. (alldone >= nbody)) EXIT NRloop
!
!--Create the initial condition. Calculate new temperature
!
!$OMP PARALLEL default(none) shared(nbody,tnew,temp_row,etot_row,ewant_row) &
!$OMP shared(det_row,error,done,nactive,new_nactive,eoslist,partype) private(m,p)
!$OMP DO schedule(runtime)
        partloop2 : DO m = 1,nactive
          p = eoslist(m)
          IF ((done(p).EQV..true.) .OR. (partype(p) == 2)) CYCLE

! Update temperature
          tnew(p) = temp_row(p)-(etot_row(p)-ewant_row(p))/det_row(p)

! Calculate the error
          error(p) = ABS(tnew(p)-temp_row(p))/temp_row(p)

! Store new temperature
          temp_row(p) = tnew(p)

! If freezing, keep temperature inside EOS tables
          IF (temp_row(p) < tmin) THEN
             temp_row(p) = tmin
             error(p)    = 0.1*eos_tol
          ENDIF
        ENDDO partloop2
!$OMP END DO
!$OMP END PARALLEL
!
!--Call EOS
!
        jlo_eos = 1
        jhi_eos = nactive
        CALL helmeos
!
!--New temperature
!
!$OMP PARALLEL default(none) shared(nactive,tnew,temp_row,etot_row,ewant_row) &
!$OMP shared(det_row,error,done,eoslist,new_nactive,partype) private(m,p) reduction(+:alldone)
!$OMP DO schedule(runtime)
        partloop3 : DO m = 1,nactive
          p = eoslist(m)
          IF ((done(p).EQV..true.) .OR. (partype(p) == 2)) CYCLE

!--Update temperature
          tnew(p) = temp_row(p) - (etot_row(p)-ewant_row(p))/det_row(p)

!--Do no allow temperature to change more than an order of magnitude per
!  iteration
          IF (tnew(p) > 10.*temp_row(p)) tnew(p) = 10.0*temp_row(p)
          IF (tnew(p) < 0.1*temp_row(p)) tnew(p) = 0.1*temp_row(p)

!--Calculate the error
          error(p) = ABS(tnew(p)-temp_row(p))/temp_row(p)

!--Store new temperature
          temp_row(p) = tnew(p)

!--If freezing, keep temperature inside EOS tables
          IF (temp_row(p) < tmin) THEN
             temp_row(p) = tmin
             error(p)    = 0.1*eos_tol
          ENDIF

!--If too hot, keep temperature inside EOS tables
          IF (temp_row(p) > tmax) THEN
             temp_row(p) = tmax
             error(p) = 0.1*eos_tol
          ENDIF

!--Decide whether particle is done or not
          IF ((error(p) < eos_tol) .AND. (done(p).EQV..false.)) THEN
             done(p) = .true.
             alldone = alldone + 1
          ENDIF
        ENDDO partloop3
!$OMP END DO
!$OMP END PARALLEL

!--Construct the new list of active particles
        nactive = 0
        DO p = 1,nbody
           IF (done(p).EQV..false.) THEN
              nactive = nactive + 1
              eoslist(nactive) = p
           ENDIF
        ENDDO     
      ENDDO NRloop
!
!--Do a last loop to store data. If the Newton-Rapshon fails to find 
!  a valid temperature, keep it  constant. Check also if temperature 
!  and the temperature predicted from the internal energy differe more 
!  than a 5%
!$OMP PARALLEL default(none) shared(nbody,vxyzut,done,press) &
!$OMP shared(css,cvs,dPdT,cps,temp_row,etot_row,ptot_row,cs_row) &
!$OMP shared(cv_row,dpt_row,cp_row,partype) private(p,eosflag,rel)
!$OMP DO schedule(runtime)
      partloop4 : DO p = 1,nbody
        IF (partype(p) == 2) CYCLE
!
        rel = ABS(temp_row(p)-vxyzut(5,p))/vxyzut(5,p)
        IF ((done(p).EQV..false.) .OR. (rel > 0.05)) THEN
           eosflag = 2
        ELSE
           eosflag = 1
        END IF
!
!--If eosflag=1 we store temp_row, whereas for eosflag=2 we store 
!  etot_row
!
        IF (eosflag == 2) THEN  ! Temperature as input
           vxyzut(4,p) = etot_row(p)*(unm/uen)
        ELSE                    ! Internal energy as input
           vxyzut(5,p) = temp_row(p)
        END IF
        press(p) = ptot_row(p)/up
        css(p)   = cs_row(p)/uv
        cvs(p)   = cv_row(p)*(unm/uen)
        dPdT(p)  = dpt_row(p)*(unl**3/uen)
        cps(p)   = cp_row(p)/(uen/unm)
      ENDDO partloop4
!$OMP END DO
!$OMP END PARALLEL
!
      END SUBROUTINE EOS
