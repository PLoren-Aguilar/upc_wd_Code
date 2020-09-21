subroutine indata(alloc)
!===================================================================
!  This subroutine reads the initial data.
!
!  Last revision: 6/April/2019
!===================================================================
!
!--Load modules
!
  use mod_parameters, only : ndim, nel, MASTER
  use mod_commons
  !use ifPorT
!     
!--Force to declare EVERYTHING
!
  implicit none
!
!--I/O variables
!
  integer, intent(in) :: alloc
!
!--Local variables
!
  real          :: r2
  integer       :: i, p, k, m, nums
  integer       :: AllocateStatus=0, DeallocateStatus=0
  character(30) :: filename1, filename2, numstep
!
!--Declare structures in order to read data
!
  type :: outdata
    real(4), dimension(ndim+2) :: ioxyzhm
    real(4), dimension(ndim+2) :: iovxyzut
    real(4) :: iorho
    real(4) :: ioka1
    integer :: iostep
    integer(1) :: iopartype
    integer(1) :: iostar
  end type
  type(outdata), allocatable, dimension(:) :: data_array
!
  type :: coutdata
     real(4), dimension(nel) :: comp
  end type
  type(coutdata), allocatable, dimension(:) :: cdata_array
!
!--Open data unit
!
  if ((mode(1) == 'run') .or. (mode(1) == 'putlayer')) then
    write(numstep,'(i4.4)') nstep_infile
    filename1 = 'bodi'//trim(numstep)//'.out'
    filename2 = 'comp'//trim(numstep)//'.out'
    open (1, FILE=filename1, form='unformatted')
    open (2, FILE=filename2, form='unformatted')
  endif
!
  if (mode(1) == 'binary') then
    write(*,*) 'Name of the file?'
    read*, filename1
    open (1, FILE=filename1, form='unformatted')
  endif
!
!--Read time and npart data
!
  read(1) tnow, npart1, nbody, ndead, nbody1, nbody2, Omega0
!
  if ((mode(1) == 'binary') .and. (alloc == 0)) npart0 = npart1
  if ((mode(1) == 'run') .or. (mode(1) == 'putlayer')) npart = npart1
!
!--Allocate space for bodi and comp data
!
  if (mode(1) == 'binary') then
    if (allocated(data_array)) deallocate(data_array, stat=DeallocateStatus)
    if (DeallocateStatus /= 0 ) then
      print*, "Error deallocating data_array"
      stop
    endif
  endif
!  
  allocate (data_array(npart1), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for data_array!!!"
    stop
  endif
!
  if ((mode(1) == 'run') .or. (mode(1) == 'putlayer')) then
    allocate (cdata_array(npart1), stat = AllocateStatus)
    if (AllocateStatus /= 0) then
      print*, "Not enough memory for cdata_array!!!"
      stop
    endif
  endif
!
!--Read bodi and comp data and close files
!
  read(1) data_array
  close(unit=1)

  if ((mode(1) == 'run') .or. (mode(1) == 'binary') .or. (mode(1) == 'putlayer')) then
    read(2) cdata_array
    close(unit=2)
  endif
!
!--Allocate particle data arrays and vectors
!
  if (alloc == 0) then
    allocate (xyzhm(ndim+2,2*npart+2), xyzhmp(ndim+2,2*npart+2), stat = AllocateStatus)
    if (AllocateStatus /= 0) then
      print*, "Not enough memory for xyzhm!!!"
      stop
    endif
!
    allocate (vxyzut(ndim+2,npart), vxyzutp(ndim+2,npart), stat = AllocateStatus)
    if (AllocateStatus /= 0) then
      print*, "Not enough memory for vxyzut!!!"
      stop
    endif
!
    allocate (axyzut(ndim+2,npart), axyzutp(ndim+2,npart), stat =  AllocateStatus)
    if (AllocateStatus /= 0) then
      print*, "Not enough memory for axyzut!!!"
      stop
    endif
!
    allocate (gxyzu(ndim+1,npart), gxyzup(ndim+1,npart), stat = AllocateStatus)
    if (AllocateStatus /= 0) then
      print*, "Not enough memory for axyzut!!!"
      stop
    endif
!
    allocate (rho(npart), rhoG(npart), norm(npart), stat = AllocateStatus)
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
    allocate (partype(npart), stat = AllocateStatus)
    if (AllocateStatus /= 0) then
      print*, "Not enough memory for partype!!!"
      stop
    endif
!
    allocate (eosflag(npart), stat = AllocateStatus)
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
    allocate (rotforc(ndim,npart), stat = AllocateStatus)
    if (AllocateStatus /= 0) then
      print*, "Not enough memory for rotforc!!!"
      stop
    endif
!
    allocate (fh(npart), stat = AllocateStatus)
    if (AllocateStatus /= 0) then
      print*, "Not enough memory for fh!!!"
      stop
    endif
!
    allocate (dhdt(npart), dhdtp(npart), stat = AllocateStatus)
    if (AllocateStatus /= 0) then
      print*, "Not enough memory for fh!!!"
      stop
    endif
!
    allocate (xss(nel,npart), stat = AllocateStatus)
    if (AllocateStatus /= 0) then
      print*, "Not enough memory for xss!!!"
      stop
    endif
!
    allocate (step(npart), stat = AllocateStatus)
    if (AllocateStatus /= 0) then
      print*, "Not enough memory for step!!!"
      stop
    endif
!
    allocate (luminuc(npart), luminucp(npart), stat = AllocateStatus)
    if (AllocateStatus /= 0) then
      print*, "Not enough memory for luminuc!!!"
      stop
    endif
!
    allocate (enuc(npart), enucp(npart), stat = AllocateStatus)
    if (AllocateStatus /= 0) then
      print*, "Not enough memory for enuc!!!"
      stop
    endif
  endif
!
!--Transfer data into the appropriate locations
!
  if (alloc == 1) then 
    p = npart0
  elseif (alloc == 0) then
    p = 0
  else
    write(*,*) 'Wrong value for alloc !'
  endif
!
  do i = 1,npart1
    p = p + 1
    xyzhmp(1:ndim+2,p)  = data_array(i)%ioxyzhm(1:ndim+2)
    vxyzutp(1:ndim+2,p) = data_array(i)%iovxyzut(1:ndim+2)
    xyzhm(1:ndim+2,p)   = xyzhmp(1:ndim+2,p)
    vxyzut(1:ndim+2,p)  = vxyzutp(1:ndim+2,p)
    rho(p)              = data_array(i)%iorho
    ka1(p)              = data_array(i)%ioka1
    step(p)             = data_array(i)%iostep
    partype(p)          = data_array(i)%iopartype
    star(p)             = data_array(i)%iostar
    if (mode(1) /= 'binary') then
      xss(1:nel,p)      = cdata_array(i)%comp(1:nel)
    endif
!
!    r2 = xyzhm(1,p)**2 + xyzhm(2,p)**2 + xyzhm(3,p)**2
!    if (sqrt(r2) > 0.12) then
!      xyzhmp(1,p) = 0.06 + rand()*0.06
!      xyzhmp(2,p) = 0.06 + rand()*0.06
!      xyzhmp(3,p) = 0.06 + rand()*0.06
!      vxyzutp(1,p) = 0.0
!      vxyzutp(2,p) = 0.0
!      vxyzutp(3,p) = 0.0
!      xyzhm(1,p) = xyzhmp(1,p)
!      xyzhm(2,p) = xyzhmp(2,p)
!      xyzhm(3,p) = xyzhmp(3,p)
!      vxyzut(1,p) = 0.0
!      vxyzut(2,p) = 0.0
!      vxyzut(3,p) = 0.0
!    endif
  enddo
!
!--Deallocate I/O arrays
!
  if (allocated(data_array))  deallocate(data_array, stat=DeallocateStatus)
  if (allocated(cdata_array)) deallocate(cdata_array, stat=DeallocateStatus)
!
!--Exit if setting up a binary system
!
  if (mode(1) == 'binary') return
!
!--Start up the rest of variables
!
  allocate (vsigmax(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for vsigmax!!!"
    stop 
  endif
!     
  allocate (cur(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for cur!!!"
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
  allocate (css(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for css!!!"
    stop
  endif
!
  allocate (press(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for press!!!"
    stop
  endif
!
  allocate (cvs(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for cvs!!!"
    stop
  endif
!
  allocate (dPdT(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for dPdT!!!"
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
  allocate (dtnuc(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for dtnuc!!!"
    stop
  endif
!
  allocate (istep0(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for istep0!!!"
    stop
  endif
!
  allocate (active(npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for active!!!"
    stop
  endif
!
  allocate (fix(ndim,npart), stat = AllocateStatus)
  if (AllocateStatus /= 0) then
    print*, "Not enough memory for fix!!!"
    stop
  endif
!
  nstep_ini = 0
  nstep     = 0
  vsigmax   = 1.0
  cur       = 0.0
  div       = 0.0
  luminuc   = 0.0
  enuc      = 0.0
  tscdyn    = 0.0
  tscnuc    = 0.0
  dtnuc     = 0.0
  if (SIMTYPE == 1)  then
    dtmax = dtmaxin
  else
    dtmax = dtmaxin*Omega0
  endif
  dt0 = dt0in*dtmax
!
  eps = 1d30
  do p = 1, npart
    eps       = min(eps,xyzhm(4,p))
    ilist(p)  = p
    plist(p)  = p
    vxyzut(5,p) = 1.0e7
  enddo
  eps  = 1.4d0*2.d0*eps
  eps3 = eps*eps*eps
!
  do p = 1, npart
    if (star(p) == 1) then
      fix(1,p) = fixx1 
      fix(2,p) = fixy1 
      fix(3,p) = fixz1
    elseif (star(p) == 2) then
      fix(1,p) = fixx2
      fix(2,p) = fixy2
      fix(3,p) = fixz2
    endif
  enddo 
!
!--Sort particle list to get rid of dead particles
!
  call sort
#ifdef debug
  if (rank == MASTER) print*, 'sort called'
#endif
!
!--Set the maximum number of neighbours and nodes in the tree and
!  allocate space for the neighbour list
!
  nbmax = 2000
  nbmin = 1
  mmax   = 2*npart + 2
end subroutine indata
