subroutine analysis
!========================================================================
!
!  This subroutine is a call to various analysis subroutines
!
!  Last revision: 14/March/2015
!========================================================================
! 
!--Force to declare EVERYTHinG
!
  implicit none
!
!--Startout the code reading the necessary data files and parameters
!
  call startout
!
!--Analize the results
!
  call analyze
!
!--Do some chemical element data files in order to render
!
  call outchem
!
end subroutine analysis
!
subroutine analyze
!
!--Load modules
!
  use mod_commons
  use mod_parameters
!
!--Force to declare EVERYTHinG
!
  implicit none
!
!--Local definitions
!
  real, dimension(nel) :: xssp
  real :: masa ,mass, dens, temp, rin, rout, rmax, dr, r,             &
          xcm, ycm ,zcm, angx, angy, angz, Imom, RWD, omega, v,       &
          rest, cvsp, cpsp, dPdTp, gg, presion, gama, vxcm, vycm, vzcm
  real, PARAMETER :: alpha = 1.0, LL = 0.1*unl, GNew = 6.67e-8
  integer :: p, k, n, l, num, inum
  
  type :: outdata
     real, dimension(ndim+2) :: ioxyzhm
     real, dimension(ndim+2) :: iovxyzut
     real :: iorho
     real :: ioka1
     integer :: iopartype
     integer :: iostar
  end type
  type(outdata), dimension(npart) :: data_array
!
  type :: coutdata
     real, dimension(nel) :: comp
  end type
  type(coutdata), dimension(npart) :: cdata_array
!
!--Do a test to find the most efficient way to write bodi binary files
!
  do p = 1,npart
    data_array(p)%ioxyzhm(1:ndim+2)  = xyzhm(1:ndim+2,p)
    data_array(p)%iovxyzut(1:ndim+2) = vxyzut(1:ndim+2,p)
    data_array(p)%iorho = rho(p)
    data_array(p)%ioka1 = ka1(p)
    data_array(p)%iopartype = partype(p)
    data_array(p)%iopartype = star(p)
    cdata_array(p)%comp(1:nel) = xss(1:nel,p)
  enddo
! 
  open (1, file='test_bin', form="unformatted")
  write(1) tnow,npart,data_array
  close(1)
!
  do p=1, npart
    data_array(p)%ioxyzhm(1:ndim+2) = 0.0
  enddo
!
  open (1, file='test_bin', form="unformatted")
  read(1) tnow,num,data_array
  close(1)
  print*, tnow, num, data_array(100)%ioxyzhm(1:ndim+2)
!
!--Now let's go for comp binary files
!
  open (2, file='ctest_bin', form="unformatted")
  write(2) tnow, cdata_array
  close(2)
  stop      
!
!--Build a bodi file with showing the density of a given
!  chemical element
!
!      do p = 1, npart
!         rest = 0.0
!         rest = rest + xss(2,p)
!         do k=5,nel
!            rest = rest + xss(k,p)
!         enddo
!         !if (xss(2,p).GT.1e-10) then
!         !   PRinT*, p,rest,(xss(k,p),k=1,nel)
!         !   STOP
!         !endif
!         if (partype(p).EQ.0) then
!            vxyzut(4,p) = rest
!         ELSE
!!            vxyzut(4,p) = 0.0
!         endif
!      enddo
!
!      call noninertial
!
!--Write results
!
!      nstep = nstep + 1
!      call outdata
!      RETURN
!
!--Open I/O files
!
  open (unit=1, file='profiles.dat', status='new')
  open (unit=2, file='totals.dat',   status='new')
!
!--This sould be the headers for the profiles.dat file
!
!# rad(0.1Rsun)  Temp(K)  Dens(g/cm^3)   Omega     Mass(M_sun) cp(erg/gK) dPdT(dyn/K) qc(erg/cm^2s) 1-1/gamma  Press (dyn)
!#===================================================================================================================================================#
!# rad(0.1Rsun)  Temp(K)  Dens(g/cm^3) Mass(M_sun) xcm(0.1Rsun)  ycm()      zcm()     vxcm (Km/s)    vycm ()     vzcm()  v()  He(fraction)  C()  O()
!#===================================================================================================================================================#
!
!--Position all data around the centre of the primary star
!
  mass  = 0.0
  xcm   = 0.0
  ycm   = 0.0
  zcm   = 0.0
  vxcm  = 0.0
  vycm  = 0.0
  vzcm  = 0.0
  do p = 1, nbody
    !if (star(p) == 2) then 
    if (rho(p) > 1.0) then 
      mass = mass + xyzhm(5,p)
      xcm  = xcm  + xyzhm(5,p)*xyzhm(1,p)
      ycm  = ycm  + xyzhm(5,p)*xyzhm(2,p)
      zcm  = zcm  + xyzhm(5,p)*xyzhm(3,p)
!
      vxcm  = vxcm  + xyzhm(5,p)*vxyzut(1,p)
      vycm  = vycm  + xyzhm(5,p)*vxyzut(2,p)
      vzcm  = vzcm  + xyzhm(5,p)*vxyzut(3,p)
    endif
  enddo
  xcm = xcm/mass
  ycm = ycm/mass
  zcm = zcm/mass
!
  vxcm = vxcm/mass
  vycm = vycm/mass
  vzcm = vzcm/mass

  do p = 1, nbody
    xyzhm(1,p) = xyzhm(1,p) - xcm
    xyzhm(2,p) = xyzhm(2,p) - ycm
    xyzhm(3,p) = xyzhm(3,p) - zcm
!
    vxyzut(1,p) = vxyzut(1,p) - vxcm
    vxyzut(2,p) = vxyzut(2,p) - vycm
    vxyzut(3,p) = vxyzut(3,p) - vzcm
  enddo

  do p = 1, nbody
    r = sqrt(xyzhm(1,p)**2  + xyzhm(2,p)**2  + xyzhm(3,p)**2)
    v = sqrt(vxyzut(1,p)**2 + vxyzut(2,p)**2 + vxyzut(3,p)**2)
    !print*, r,v*(1e-5*unl/unt)
  enddo
!
!--Calculate radial profiles of density, temperature... etc in a
!  cross-section around the mid-plane
!
  n    = 1000
  rmax = log10(4.0)
  rin  = log10(1e-6)
  dr   = (rmax-rin)/n
  rout = rin + dr
  masa = 0.0
     
  do k = 1,n
    num   = 0.0 
    mass  = 0.0
    dens  = 0.0
    temp  = 0.0
    omega = 0.0
    xssp  = 0.0
    do p=1, nbody
      r = sqrt(xyzhm(1,p)**2  + xyzhm(2,p)**2 + xyzhm(3,p)**2)
      v = sqrt(vxyzut(1,p)**2 + vxyzut(2,p)**2 + vxyzut(3,p)**2)
      if ((log10(r) > rin) .and. (log10(r) <= rout)) then 
        masa = masa + xyzhm(5,p)
        num= num + 1
        if (abs(xyzhm(3,p)) < 0.005) then
          mass  = mass  + xyzhm(5,p)
          temp  = temp  + xyzhm(5,p)*vxyzut(5,p)
          dens  = dens  + xyzhm(5,p)*rho(p)
          omega = omega + xyzhm(5,p)*v
          do l=1,nel
            xssp(l) = xssp(l) + xyzhm(5,p)*xss(l,p)
          enddo
        endif
      endif
    enddo
!
    if (mass > 0.0) then
      temp  = temp/mass
      dens  = dens/mass
      omega = omega/mass
      do l=1,nel
        xssp(l) = xssp(l)/mass
      enddo
!
!--Call EOS to calculate the value of the specific heat
!
      call deg(dens,temp,presion,xssp,cvsp,cpsp,dPdTp,gama)
!
!--Write results
!
      gg = (GNew*mass*unm/(unl*10**((rin+rout)/2.))**2)**0.5
      !write(1,'(11(1pe12.4))') 10**((rin+rout)/2.),temp,dens*uden,&
      !                        (omega/unt)/(10**((rin+rout)/2.)),  &
      !                        masa,cvsp,dPdTp,                    &
      !                        alpha**2*(dens*uden)*cvsp*temp*     &
      !                        (gg*LL)**0.5*(2**1.5),1 - cvsp/cpsp,&
      !                        presion,gama
      !if (dens > 1.0) then
      write(1,'(16(1pe12.4))') 10**((rin+rout)/2.),temp,dens*uden,     &
                               masa,xcm,ycm,zcm,vxcm*(1.e-5*unl/unt),   &
                               vycm*(1e-5*unl/unt),vzcm*(1e-5*unl/unt), &
                               omega*(1e-5*unl/unt),xssp(4)/xssp(2),    &
                               xssp(4)/xssp(3),1.,xssp(4)/xssp(5),      &
                               xssp(4)/xssp(6)
      !endif
    endif
    rin  = rout
    rout = rout + dr
  enddo
!
!--Calculate TOTAL values for the angular momentum and moment of inertia
!  for the WD and disk
!
  RWD  = 0.1
  angx = 0.0
  angy = 0.0
  angz = 0.0
  Imom = 0.0
  mass = 0.0
  do p = 1, nbody
    r = sqrt(xyzhm(1,p)**2 + xyzhm(2,p)**2)
    if (r <= RWD) then
      mass = mass + xyzhm(5,p)
      angx = angx + xyzhm(5,p)*(xyzhm(2,p)*vxyzut(3,p) -         &
             xyzhm(3,p)*vxyzut(2,p))
      angy = angy + xyzhm(5,p)*(xyzhm(1,p)*vxyzut(3,p) -         &
             xyzhm(3,p)*vxyzut(1,p))
      angz = angz + xyzhm(5,p)*(xyzhm(1,p)*vxyzut(2,p) -         &
             xyzhm(2,p)*vxyzut(1,p))
      Imom = Imom + xyzhm(5,p)*r*r
    endif
  enddo
  write(2,*) 'Values for the WD'
  write(2,'(5(1pe12.4))') angx*unm*unl**2/unt,angy*unm*unl**2/unt,  &
                              angz*unm*unl**2/unt,Imom*unm*unl**2, mass
!
  angx = 0.0
  angy = 0.0
  angz = 0.0
  Imom = 0.0
  mass = 0.0
  do p = 1, nbody
    r = sqrt(xyzhm(1,p)**2 + xyzhm(2,p)**2)
    if (r > RWD) then
      mass = mass + xyzhm(5,p)
      angx = angx + xyzhm(5,p)*(xyzhm(2,p)*vxyzut(3,p) -         &
            xyzhm(3,p)*vxyzut(2,p))
      angy = angy + xyzhm(5,p)*(xyzhm(1,p)*vxyzut(3,p) -         &
            xyzhm(3,p)*vxyzut(1,p))
      angz = angz + xyzhm(5,p)*(xyzhm(1,p)*vxyzut(2,p) -         &
             xyzhm(2,p)*vxyzut(1,p))
      Imom = Imom + xyzhm(5,p)*r*r
    endif
  enddo
  write(2,*) 'Values for the Disk'
  write(2,'(5(1pe12.4))') angx*unm*unl**2/unt,angy*unm*unl**2/unt,  &
                            angz*unm*unl**2/unt,Imom*unm*unl**2, mass
!
!--Close I/O files
!
  close(1)
  close(2)
!
end subroutine analyze

subroutine deg(rhop,temp,presion,xssp,cvsp,cpsp,dPdTp,gama)
!========================================================================
!     Call to the equation of state
!
!     Last revision: 9/April/2019
!=========================================================================
!
!--Load modules
!
  use mod_parameters
  use mod_commons
  use mod_EOS
!
!--Force to declare EVERYTHinG
!
  implicit none
!
!--Helmholtz EOS definitions
!
  !inCLUDE 'vector_eos.dek'

!
!--I/O variables
!
  real, dimension(nel), intent(in) :: xssp
  real, intent(in) :: rhop, temp
  real, intent(out) :: gama, presion, cvsp, cpsp, dPdTp
!
!--Local variables
!
  !integer :: jlo_eos, jhi_eos
  !real, dimension(nrowmax) :: ewant_row, temp_row, den_row, abar_row, &
  !zbar_row, det_row, ptot_row, cs_row, cv_row, dpt_row, cp_row, etot_row
  !LOGICAL :: eosfail
  real    :: abar, zbar
  integer :: p, k
!
!--Set initial data for the EOS
!
  temp_row(1) = temp
  if (temp > tmax) temp_row(1) = tmax
  if (temp < tmin) temp_row(1) = tmin
  den_row(1) = rhop*uden    ! EOS works in cgs units!!!
!
  abar = 0.0
  zbar = 0.0
  do k=1,nel-1
    abar = abar + xssp(k)/aion(k)
    zbar = zbar + xssp(k)*zion(k)/aion(k)
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
  !call helmeos(jlo_eos,jhi_eos,temp_row,den_row,abar_row,zbar_row,  &
  !      det_row,ptot_row,cs_row,cv_row,dpt_row,cp_row,eosfail)
  call helmeos
  if (rank == MASTER .and. eosfail.eqv..true.) then
    print*,'EOScorr fail'
    stop
  endif

  presion  = ptot_row(1)
  cvsp   = cv_row(1)
  cpsp   = cp_row(1)
  dPdTp  = dpt_row(1)
  gama   = nabad_gas_row(1)

end subroutine deg
