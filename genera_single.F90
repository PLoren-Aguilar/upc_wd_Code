subroutine genera_single
!====================================================================
!  Program to generate the initial state of the partciles that will
!  be used in the sph program
!
!  Last revision 9/April/2019  by PLA
!====================================================================
!
!--Load modules
!
  use mod_parameters
  use mod_commons
  !use ifPORT
  use mod_EOS
! 
!--Force to declare EVERYTHING
!
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif
!
!--Local definitions
!
  real, dimension(maxit) :: radle, densle, normdensle 
  real :: xc, yc, zc, dh, dh2, rx, ry, rz, radmax,                  &
             dr2, rad, masstot, radmin, temp, massp, sum1, sum2,      &
             rhoc_min, rhoc_max, rhoc, mass_min, mass_max, mass,       &
             REL, abar, zbar, t1, t2, rand1, rand2, theta, phi, rnddens, rnddens2, densmax
  real, parameter :: tolerance = 0.01
  integer :: p, k, i, j, l, nbod, NNNZ, endit, iseed, AllocateStatus, ndone, i0, i1
  character(13) :: ONC
  logical :: done
!
!--Startout MPI, if necessary
!
#ifdef MPI
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
#else
  size = 1
  rank = MASTER
#endif
  nprocs = size
!
!--Read initial model parameters
!
  open (unit=1,file='genparams',status='old')
  read(1,*) nbody
  read(1,*) masstot
  read(1,*) temp
  read(1,*) tnow
  close(1)
!
!--Setup the number of particles
!
  npart  = nbody
  ndead  = 0
  nbody1 = npart
  nbody2 = 0
!
!--Setup the neighbour list matrix size
!
  nbmax = 2000
  nbmin = 1
  mmax  = 2*npart + 2
!
!--Allocate particle data arrays and vectors
!
  allocate (xyzhm(ndim+2,2*npart+2), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for xyzhm!!!"
    stop
  endif
!
  allocate (xyzhmp(ndim+2,2*npart+2), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for xyzhm!!!"
    stop
  endif
!
  allocate (vxyzut(ndim+2,npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for vxyzut!!!"
    stop
  endif
!
  allocate (vxyzutp(ndim+2,npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for vxyzut!!!"
    stop
  endif
!
  allocate (axyzut(ndim+2,npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for vxyzut!!!"
    stop
  endif
!
  allocate (axyzutp(ndim+2,npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for vxyzut!!!"
    stop
  endif
!
  allocate (gxyzu(ndim+1,npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for vxyzut!!!"
    stop
  endif
!
  allocate (rho(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for rho!!!"
    stop
  endif
!
  allocate (ka1(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for ka1!!!"
    stop
  endif
! 
  allocate (div(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for div!!!"
    stop
  endif
!
  allocate (divt(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for divt!!!"
    stop
  endif
!
  allocate (cur(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for cur!!!"
    stop
  endif
!
  allocate (css(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for css!!!"
    stop
  endif
!
  allocate (cvs(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for cvs!!!"
    stop
  endif
!
  allocate (cps(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for cps!!!"
    stop
  endif
!
  allocate (uintprev(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for uintprev!!!"
    stop
  endif
!
  allocate (luminuc(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for luminuc!!!"
    stop
  endif
!
  allocate (luminucp(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for luminucp!!!"
    stop
  endif
!
  allocate (dhdt(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for dhdt!!!"
    stop
  endif
!
  allocate (dhdtp(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
  print*, "Not enough memory for dhdtp!!!"
  stop
      endif
!
  allocate (vsigmax(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for vsigmax!!!"
    stop
  endif
!
  allocate (enucp(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for enucp!!!"
    stop
  endif
! 
  allocate (dPdT(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for dPdT!!!"
    stop
  endif
!
  allocate (press(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for press!!!"
    stop
  endif
!
  allocate (fh(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for fh!!!"
    stop
  endif
!
  allocate (tscdyn(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for tscdyn!!!"
    stop
  endif
!
  allocate (tscnuc(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for tscnuc!!!"
    stop
  endif
!
  allocate (enuc(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for enuc!!!"
    stop
  endif
!
  allocate (partype(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for partype!!!"
    stop
  endif
!
  allocate (star(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for star!!!"
    stop
  endif
!
  allocate (ilist(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for ilist!!!"
    stop
  endif
!
  allocate (plist(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for plist!!!"
    stop
  endif
!
  allocate (xss(nel,npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for xss!!!"
    stop
  endif
!
  allocate (nb(nbmax,npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for nb!!!"
    stop
  endif
!
  allocate (step(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for step!!!"
    stop
  endif
!
  allocate (istep0(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for istep0!!!"
    stop
  endif
!
!--Setup WD composition
!
  xss = 0.0
  if (masstot < 0.4) then
    do p=1,nbody
      xss(2,p) = 1.0
    enddo
  elseif ((masstot >= 0.4) .and. (masstot < 1.15)) then
    do p=1,nbody
      xss(3,p) = 0.5
      xss(4,p) = 0.5
    enddo
  elseif (masstot >= 1.15) then
    do p=1,nbody
      xss(4,p) = 0.8
      xss(5,p) = 0.2
    enddo
  endif
!
!--Use a Bisection strategy in order to solve Lane-Embden 
!  equation for the given mass
!
  write(*,*) ' '
  write(*,*) '************************************************** '
  write(*,*) 'Using a bisection scheme to solve the lane-embden'
  write(*,*) 'equation and generate a first density vs radius'
  write(*,*) 'profile of the WD...'
  write(*,*) '************************************************** '
  write(*,*) ' '
  done     = .false.
  rhoc_min = 1.0d4
  rhoc_max = 1.0d8
  do while (done.eqv..false.)

    call lanembden(maxit,rhoc_min,endit,mass,radle,densle)
    mass_min = mass

    call lanembden(maxit,rhoc_max,endit,mass,radle,densle)
    mass_max = mass

    rhoc = 0.5*(rhoc_min + rhoc_max)
    call lanembden(maxit,rhoc,endit,mass,radle,densle)

    if (masstot > mass) then
      rhoc_min = rhoc 
    elseif (masstot < mass) then
      rhoc_max = rhoc
    endif 
!               
    REL = ABS(masstot-mass)/masstot
    if (REL <= tolerance) then
      done = .true.
    endif
    write(*,'(A5,1(1F7.4),A8,1F7.4)') 'M_WD=',masstot,' M_Bisec=',mass
  enddo
!
!--Call lanembden subroutine a last time, saving this time the
!  values for density, radius, pressure, etc
  call lanembden(maxit,rhoc,endit,mass,radle,densle)
  massp  = mass/FLOAT(npart)
!
!--Randomly distribute the particles according to the found density profile
!
  !call random_seed()
  call srand(time())
  p = 0
  partloop : do while (p < nbody)
    rnddens2 =  0.0
    rnddens  = -1.0
    do while (rnddens2 > rnddens)
      !    
      !Generate a uniformely distributed random number to select a model density
      !call random_number(rand1)
      rand1 = rand()
      i = int(endit*rand1) + 1
      if (i == endit + 1) i = endit
      rnddens = densle(i)
      !
      !Now, generate a second random number to select a radius for the particle
      !call random_number(rand2)
      rand2 = rand()
      rnddens2 = rand2*densle(1)
    enddo
    !
    !And finally, put the particle at the found radius
    p = p + 1
  
    i0   = 1
    i1   = endit
    done = .false.
    do while (densle(i0) > rnddens2) 
      i0 = i0 + 1
    enddo
    i1  = i0
    i0  = i0-1
    rad = radle(i0) + (rnddens2-densle(i0))*(radle(i1)-radle(i0))/(densle(i1)-densle(i0))
    rad = 10.*rad
    
    !call random_number(rand1)
    !call random_number(rand2)
    rand1 = rand()
    rand2 = rand()
    theta       = 2.0*pi*rand1
    phi         = acos(1.0-2.0*rand2)
    xyzhm(1,p)  = rad*sin(phi)*cos(theta)
    xyzhm(2,p)  = rad*sin(phi)*sin(theta)
    xyzhm(3,p)  = rad*cos(phi)
    xyzhm(4,p)  = hfact*(massp*uden/rnddens2)**0.333
    xyzhm(5,p)  = massp
    rho(p)      = rnddens2/uden
    !
    !And complete the rest of the particle's properties
    vxyzut(1,p) = 0.
    vxyzut(2,p) = 0.
    vxyzut(3,p) = 0.
    vxyzut(5,p) = temp
    ka1(p)      = ka_min
    fh(p)       = 1.
    step(p)     = 0
    istep0(p)   = 0
    partype(p)  = 0
  enddo partloop
!
!--Save eps for late use. Store the p's variable for fiagnostic purposes
!
  eps = 1.0d30
!$omp parallel default(none) shared(nbody,xyzhm,xyzhmp,vxyzut,vxyzutp) private(p)  &
!$omp reduction(min:eps)
!$omp do schedule(runtime)
  do p=1,nbody
     eps  = min(eps,xyzhm(4,p))
     xyzhmp(:,p)  = xyzhm(:,p)
     vxyzutp(:,p) = vxyzut(:,p)
  enddo
!$omp end do
!$omp end parallel
  eps  = 1.4d0*2.0d0*eps
  eps3 = eps*eps*eps
!
!--Sort particles. Since now factor is a dynamical quantity we need to
!  call sort in order to define it !!
!
  call sort
!
!--Now, use the tree to self-consistently calculate the rho-h relation
!
  print*,'Trying to compute h. May take a long time'
  nstep = 1
  step  = 0
  call iter_rhoh
  print*, 'h calculation finished'
!
!--Read EOS tables
!
  call read_helm_table
!
  open(unit=4,file='ZHELI.10b',status='old')
  aion(1) = 1.0d0
  zion(1) = 1.0d0
  do k=1,nel-2
    read(4,'(I4,1X,A5,1X,F4.0,1X,F4.0)') NNNZ,ONC,aion(k+1),zion(k+1)
  enddo
  close(unit=4)
!
!--Calculate thermal energy
!
  call eos0
  print*,'eos finished'
!
!--Write data to disk
!
  nstep  = 0
  nout   = 1
  Omega0 = 0.0
  call diagnostics(0.,0.)

end subroutine genera_single
!
subroutine lanembden(maxit,rhoc,endit,masa,r,rho)
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--I/O variables
!
  real, dimension(maxit), intent(out) :: rho, r
  real, intent(in) :: rhoc
  integer, intent(in)  :: maxit
  integer, intent(out) :: endit
!
!--Local variables
!
  real, dimension(2) :: y, yout
  real    :: x, h, t, alpha, gama, K, pres, masa, pi, alpha0
  integer :: oup, flag, i, j
!
!--Physical parameters in cgs
!
   real, parameter :: plank=6.6261d-27, me=9.1094d-28, mp=1.6726d-24, &
                      G=6.6726d-8, msol=1.9891d33, rsol=6.955d10, mue=2.0d0, npol=1.5d0, tol=1e-4
!
!--External subroutines
!
  external derivs
!
!--Initial variables
!
  gama   = (npol+1.)/npol
  pi     = 4.0*atan(1.0)
  K      = ((3.)**(2./3.)*(plank**2))/(20.*(pi**(2./3.))*me*((mp*mue)**(gama)))
  alpha0 = (((npol+1.)*K)/(4.*pi*G))**(1./2.)
  alpha  = alpha0*rhoc**((1.-npol)/(2.*npol))
  y(1)   = 1.0
  y(2)   = 0.0
!
!--Start up the intergration
!
  t      = 0.0
  h      = 1.0e-6 
  r(1)   = alpha*t/rsol
  rho(1) = rhoc*(y(1)**npol)
  pres   = K*(rho(1)**gama) 
  masa   = -4*pi*(alpha0**3.)*(rhoc**((3.-npol)/(2*npol)))*(t**2.)*(y(2))
  masa   = masa/msol
!
!--Integrate lane-embden equation
!
  j = 0
  i = 0
  do while (y(1) > 0)

    j = j + 1

    call rk4(npol,y,2,t,h,yout,derivs)

    t = t + h
    y = yout
    if (mod(j,100000) == 0) then
      i      = i + 1
      r(i)   = alpha*t/rsol
      rho(i) = rhoc*(y(1)**npol)
      pres   = K*(rho(i)**gama)               
      masa   = -4*pi*(alpha0**3.)*(rhoc**((3.-npol)/(2*npol)))*(t**2.)*(y(2))
      masa   = masa/msol
      !print*, i,r(i)
    endif
  enddo
  endit = i

end subroutine lanembden
!
subroutine derivs(npol,x,y,dydx)
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--I/O variables
!
  real, dimension(*) :: y, dydx
  real :: npol, x
!
  dydx(1) = y(2)
  if (y(1) >= 0) then
    if (abs(x)<1.0e-9) then
      dydx(2) = -(y(1)**npol)/3.0
    else
      dydx(2) = -y(1)**npol-2.0*y(2)/x
    endif
  else
    dydx(2) = 0.0
  endif

end subroutine derivs
!
subroutine rk4(npol,y,n,x,h,yout,derivs)
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--I/O variables
!
  real, dimension(n), intent(in) :: y
  real, intent(in) :: x, h, npol
  real, dimension(n), intent(out) :: yout
!
!--Local variables
!
  real, dimension(size(y)) :: k1, k2, k3, k4, yt
  real    :: h6, hh, xh
  integer :: n
!
  hh=h*0.5
  h6=h/6.0
!
  call derivs(npol,x,y,k1)
  xh=x+hh
  yt=y+hh*k1
!
  call derivs(npol,xh,yt,k2)
  yt=y+hh*k2
!
  call derivs(npol,xh,yt,k3)
  yt=y+h*k3
!
  call derivs(npol,x+h,yt,k4)
  yout=y+h6*(k1+k4+2.0*(k2+k3))
!
end subroutine rk4
