      subroutine degenerate2(p)
!========================================================================
!     THIS   subroutine   ALLOWS   TO  SELECT   BETWEEN  THE DifFERENT  
!     EQUATIONS  OF  STATE  OF  THE  CODE. IT'S IMPorTANT  TO  realIZE  
!     THAT  EVERY  EOS  subroutine MUST PROVIDE (AT LEAST):
!
!     * TEMPERATURE, PRESSURE, Cv and SPEED OF SOUND
!
!     Last revision: 21/March/2017
!=========================================================================
!
!--Load modules
!
      use mod_parameters, only : uden, unm, uen, up, uv, unl, tmin, tmax,&
      ndim, nel, MASTER
      use mod_commons,    only : rho, xss, aion, zion, vxyzut, press,   &
      css, cvs, dPdT, cps, partype, rank, RELFLAG, eosflag, SIMTYPE
      use mod_EOS
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--I/O variables
!
      integer, intent(in) :: p
!
!--Local variables
!
      real, parameter :: eos_tol=1.0d-8
      real  :: abar, zbar, tnew, densp, temp, ewantp, errorp, etotp, &
                  detp, cvp, rel
      integer  :: k, newton
      integer, parameter :: max_newton=50
!
!--Set initial data for the EOS
!
      if (vxyzut(5,p) > tmax) vxyzut(5,p) = tmax
      if (vxyzut(5,p) < tmin) vxyzut(5,p) = tmin
!
      ewant_row(1) = vxyzut(4,p)*(uen/unm) ! EOS works in cgs units!!!
      temp_row(1)  = vxyzut(5,p)
      den_row(1)   = rho(p)*uden    ! EOS works in cgs units!!!
!
!--Composition
!
      abar = 0.0
      zbar = 0.0
      do k=1,nel-1
         abar = abar + xss(k,p)/aion(k)
         zbar = zbar + xss(k,p)*zion(k)/aion(k)
      enddo
      abar_row(1) = 1.0/abar
      zbar_row(1) = zbar/abar
!
!--Set the number of elements sent into the EOS
!
      jlo_eos = 1
      jhi_eos = 1
!
!--Call EOS to obtain a first estimate for energy, pressure and derivatives
!
      call helmeos
      if ((rank == MASTER) .and. (eosfail.eqv..true.)) then
         write(*,*) 'EOScorr fail'
         stop
      endif
!
!--Finish EOS in case of a relaxation
!
      if (RELFLAG.eqv..true.) then 
        vxyzut(4,p) = etot_row(1)*(unm/uen)
        press(p) = ptot_row(1)/up
        css(p)   = cs_row(1)/uv
        cvs(p)   = cv_row(1)*(unm/uen)
        dPdT(p)  = dpt_row(1)*(unl**3/uen)
        cps(p)   = cp_row(1)/(uen/unm)
        return
      endif
!
!--Create the initial condition. Calculate new temperature
!
      tnew = temp_row(1) - (etot_row(1) - ewant_row(1))/det_row(1)
!
!--Calculate the error
!
      errorp = abs(tnew-temp_row(1))/temp_row(1)
!
!--Store new temperature
!
      temp_row(1) = tnew
!
!--If freezing, keep temperature inside EOS tables
!
      if (temp_row(1) < tmin) then
          temp_row(1) = tmin
          errorp = 0.1*eos_tol
      endif
!
!--Loop over particles doing the Newton-Raphson iteration
!
      newton = 0
      whileloop : do while ((errorp > eos_tol).and.(newton < max_newton))
         newton = newton + 1
!
!--Call EOS
!
         jlo_eos = 1
         jhi_eos = 1
         call helmeos
         if ((rank == MASTER) .and. (eosfail.eqv..true.)) print*,'EOS in&
             particle: ',p,' failed'
!
!--New temperature
!
         tnew = temp_row(1) - (etot_row(1)-ewant_row(1))/det_row(1)
!
!--Do no allow temperature to change more than an order of magnitude per
!  iteration
!
         if (tnew > 10.*temp_row(1)) tnew = 10.0*temp_row(1)
         if (tnew < 0.1*temp_row(1)) tnew = 0.1*temp_row(1)
!
!--Calculate the error
!
         errorp = abs(tnew-temp_row(1))/temp_row(1)
! 
!--Store new temperature
!
         temp_row(1) = tnew
!
!--If freezing, keep temperature inside EOS tables
!
         if (temp_row(1) < tmin) then
             temp_row(1) = tmin
             errorp = 0.1*eos_tol
         endif
!
!--If too hot, keep temperature inside EOS tables
!
         if (temp_row(1) > tmax) then
             temp_row(1) = tmax
             errorp = 0.1*eos_tol
         endif
      enddo whileloop   ! End Newton-Raphson loop
!
!--If the Newton-Rapshon fails to find a valid temperature, keep it 
!  constant. Check also if temperature and the temperature predicted
!  from the internal energy differe more than a 5%
      rel = dabs(temp_row(1)-vxyzut(5,p))/vxyzut(5,p)
      if ((newton >= max_newton).or.(rel > 0.05)) then
        eosflag(p) = 2
      else
        eosflag(p) = 1
      end if
!
!--If eosflag=1 we store temp_row, whereas for eosflag=2 we store 
!  etot_row
!
      if (eosflag(p) == 2) then  ! Temperature as input
        vxyzut(4,p) = etot_row(1)*(unm/uen)
        !error(p) = 0.0d0
      else                    ! Internal energy as input
        vxyzut(5,p) = temp_row(1)
        !error(p) = errorp
      end if
      press(p) = ptot_row(1)/up
      css(p)   = cs_row(1)/uv
      cvs(p)   = cv_row(1)*(unm/uen)
      dPdT(p)  = dpt_row(1)*(unl**3/uen)
      cps(p)   = cp_row(1)/(uen/unm)
!
      end subroutine degenerate2
