subroutine genera_bin
!=====================================================================
!  Subroutine that generates the initial state of the particles that
!  will be used in the binary simulation
!
!  Last revision: 9/April/2017 by PLA
!=====================================================================
!
!--Load modules
!
  use mod_parameters
  use mod_commons
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--Local variables
!
  real, dimension(3) :: cmp1, cmp2
  real :: m1, m2, x1, x2, v1, v2, rad1, rad2, omega, vrel, cmd
  integer ::  i, p, k, nb1, nb2, AllocateStatus
  character(30) :: file1, file2, filename
!
!--Read simulation parameters
!
  open (1,file='treepars',status='old')
  read(1,*) nstep_infile
  read(1,*) tend
  read(1,*) dtmaxin
  read(1,*) dt0in
  read(1,*) nout
  read(1,*) SIMTYPE
  read(1,*) RELFLAG
  read(1,*) reset
  read(1,*) trelax
  read(1,*)
  read(1,*)
  close(1) 
!
!--Read total number of particles
!
  print*, 'Total number of particles in the binary simulation?'
  read*, npart
!
!--Read data for the first star
!
  write(*,*) 'Reading data for the first star'
  call indata(0)
  nb1 = npart1

! Identify the first star and calculate its mass
  m1 = 0.0
  do p=1,nb1
    fh(p)   = 1.0
    star(p) = 1
    m1      = m1 + xyzhm(5,p)
  enddo 
!
!--Setup body 1 composition
!
  xss = 0.0
  if (m1 < 0.4) then
    do p=1,nb1
      xss(2,p) = 1.0
    enddo
  elseif ((m1 >= 0.4) .AND. (m1 < 1.15)) then
    do p=1,nb1
      if (partype(p) == 0) then
        xss(3,p) = 0.5
        xss(4,p) = 0.5
      elseif (partype(p) == 1) then
        xss(2,p) = 1.0
      endif
    enddo
  elseif (m1 >= 1.15) then
    do p=1,nb1
      if (partype(p) == 0) then
        xss(4,p) = 0.8
        xss(5,p) = 0.2
      elseif (partype(p) == 1) then
        xss(2,p) = 1.0
      endif
    enddo
  endif
!
!--Read data for the second star
!
  write(*,*) 'Reading data for the second star'
  call indata(1)
  nb2 = npart1

! Identify the second star and calculate its mass
  m2 = 0.0
  do p=nb1,nb1+nb2
    fh(p)   = 1.0
    star(p) = 2
    m2      = m2 + xyzhm(5,p)
  enddo
!
!--Setup body 2 composition
!
  print*, nb1,nb2,m1,m2
  if (m2 < 0.4) then
    do p=nb1,nb1+nb2
      xss(2,p) = 1.0
    enddo
  elseif ((m2 >= 0.4) .AND. (m2 < 1.15)) then
    do p=nb1,nb1+nb2
      if (partype(p) == 0) then
        xss(3,p) = 0.5
        xss(4,p) = 0.5
      elseif (partype(p) == 1) then
        xss(2,p) = 1.0
      endif
    enddo
  elseif (m2 >= 1.15) then
    do p=nb1,nb1+nb2
      if (partype(p) == 0) then
        xss(4,p) = 0.8
        xss(5,p) = 0.2
      elseif (partype(p) == 1) then
        xss(2,p) = 1.0
      endif
    enddo
  endif
!
!--Setup binary system initial positions and velocities
!
  print*,'Distance between objects?'
  read*,rad2
!
  x1 = m2*rad2/(m2+m1)
  x2 = rad2-x1
  write(*,'(4(1pe12.4))') x1,x2,m1,m2
!
  do p=1,nb1
    xyzhm(1,p) = xyzhm(1,p) - x1
  enddo
  do p=nb1+1,nb1+nb2
    xyzhm(1,p) = xyzhm(1,p) + x2
  enddo
!
  cmp1  = 0.0
  cmp2  = 0.0
  m1    = 0.0
  m2    = 0.0
!$omp parallel default(none) shared(xyzhm,nb1,nb2,star) private(p,k)      &
!$omp reduction(+:cmp1) reduction(+:cmp2) reduction(+:m1) reduction(+:m2)
!$omp do schedule(runtime)
  do p=1,nb1+nb2
    if (star(p) == 1) then
      m1 = m1 + xyzhm(5,p)
      do k=1,ndim
        cmp1(k) = cmp1(k) + xyzhm(5,p)*xyzhm(k,p)
      enddo
    elseif (star(p) == 2) then
      m2 = m2 + xyzhm(5,p)
      do k=1,ndim
        cmp2(k) = cmp2(k) + xyzhm(5,p)*xyzhm(k,p)
      enddo
    endif
  enddo
!$omp end do
!$omp end parallel
  cmp1 = cmp1/m1
  cmp2 = cmp2/m2
  if (rank == MASTER) then
    print*, 'm1=',m1, 'm2=',m2
    print*, 'Where they really are',cmp1(1),cmp2(1)
  endif
!
!--Centre of mass distance and rotational angular velocity
!
  cmd     = sqrt((cmp1(1)-cmp2(1))**2 + (cmp1(2)-cmp2(2))**2 + (cmp1(3)-cmp2(3))**2)
  Omega01 = sqrt((m1+m2)/cmd**3)
  Omega02 = Omega01
  Omega0  = Omega01
!
!--Assign fixed position for each star (necessary for binary relaxation)
!
  fixx1 = -x1
  fixy1 = 0.
  fixz1 = 0.
  fixx2 = x2
  fixy2 = 0.
  fixz2 = 0.
!
!--Corotating initial conditions
!
  omega = sqrt((m1+m2)/rad2**3)
  do p=1,nb1+nb2
    vxyzut(1,p) = 0.0
    vxyzut(2,p) = 0.0
    vxyzut(3,p) = 0.0
    vxyzut(5,p) = 1.0e7
  enddo
!
!--Non-corotating initial conditions
!
!  omega = Dsqrt((m1+m2)/rad2**3)
!  do p=1,nb1+nb2
!    vxyzut(1,p) = -xyzhm(2,p)*omega
!    vxyzut(2,p) =  xyzhm(1,p)*omega
!    vxyzut(3,p) =  0.0d0
!  enddo
!
!--Initial velocity conditions for a collision 
!
!  print*, 'Relative velocity between the stars? (Km/s)'
!  read*, vrel
!  vrel = vrel*1.0E5*(unt/unl)
!
!  v1 = vrel
!  v2 = 0.0
!  do p=1,nb1+nb2
!    vxyzut(2,p) =  0.0d0
!    vxyzut(3,p) =  0.0d0
!    if (p <= nb1) then
!      vxyzut(1,p) = vrel
!    else
!      vxyzut(1,p) = 0.0
!    endif
!  enddo
!
!--Write to disk initial data
!
  nstep  = 0
  nout   = 1
  nbody1 = nb1
  nbody2 = nb2
  call outdata
end subroutine genera_bin
