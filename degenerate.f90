      subroutine degenerate(p)
!========================================================================
!     THIS   subroutine CALLS THE EQUATION OF  STATE  OF  THE  CODE. 
!     IT'S IMPORTANT  TO  realIZE  THAT  EVERY  EOS  subroutine MUST 
!     PROVIDE (AT LEAST):
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
      use mod_commons, only : rho, xss, aion, zion, vxyzut, press, css, &
      cvs, dPdT, cps, rank
      use mod_EOS
!
!--Force to declare EVERYTHING
!
      implicit none
!
!--I/O variables
!
      integer, intent(in) :: p
!
!--Local variables
!
      real  :: abar, zbar
      integer  :: k
!
!--Set temperature and density
!
      if (vxyzut(5,p) > tmax) vxyzut(5,p) = tmax
      if (vxyzut(5,p) < tmin) vxyzut(5,p) = tmin
      temp_row(1)  = vxyzut(5,p)
      den_row(1)   = rho(p)*uden    ! EOS works in cgs units!!!
!
!--Set composition
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
      jlo_eos=1
      jhi_eos=1
!
!--Call EOS to obtain a first estimate for energy, pressure and derivatives
!
      call helmeos
      if ((rank == MASTER) .and. (eosfail.eqv..true.)) then
        write(*,*) 'EOS fail'
        stop
      endif
!
      vxyzut(4,p) = etot_row(1)*(unm/uen)
      press(p) = ptot_row(1)/up
      css(p)   = cs_row(1)/uv
      cvs(p)   = cv_row(1)*(unm/uen)
      dPdT(p)  = dpt_row(1)*(unl**3/uen)
      cps(p)   = cp_row(1)/(uen/unm)
!
      end subroutine degenerate
