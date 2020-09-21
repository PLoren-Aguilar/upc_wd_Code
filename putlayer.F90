subroutine putlayer
!========================================================================
!  This subroutine puts a He buffer on the top of an ideally relaxed WD
!
!  Last revision: 05/June/2017
!========================================================================
!
!--Load modules
!
  use mod_parameters
  use mod_commons
  !use ifPORT
! 
!--Force to declare EVERYTHING
!
  IMPLICIT NONE
!
!--Local variables
!
  real(4), dimension(nel,1000000) :: xss_He
  real(4), dimension(5,100000)   :: xyzhm_He, vxyzut_He
  real(4), dimension(1000000)     :: rho_He, fh_He, ka1_He
  real(4) :: rmax, rdown, rup, rad, theta, phi, mHe, xc, yc, zc
  integer(1), dimension(1000000) :: partype_He, star_He
  integer, dimension(1000000) :: step_He
  integer :: m, p, k, num, ntot, AllocateStatus=0, DeAllocateStatus=0
  character(30) :: filename
!
!--Declare structures in order to dump data
!
  type :: iodata
     real(4), dimension(ndim+2) :: ioxyzhm
     real(4), dimension(ndim+2) :: iovxyzut
     real(4) :: iorho
     real(4) :: ioka1
     integer :: iostep
     integer(1) :: iopartype
     integer(1) :: iostar
  end type iodata
  type(iodata), dimension(:), allocatable :: data_array
!
  type :: ciodata
     real(4), dimension(nel) :: comp
  end type ciodata
  type(ciodata), dimension(:), allocatable :: cdata_array
!
!--Startout the code reading the necessary data files and parameters
!
  nprocs = 1
  CALL startout
!
!--How many He particles do you want in the spherical layer? 
!
  write(*,*) 'How many He particles do you want?'
  READ(*,*) nHe 
  write(*,*) 'Mass of the He layer (in solar masses)?'
  READ(*,*) mHe 
!
!--Find out the maximum radius of the WD
!
  rmax = 0.15*0.15
  !do p=1,nbody
  !  rmax = MAX(rmax,xyzhm(1,p)**2 + xyzhm(2,p)**2 + xyzhm(3,p)**2)
  !enddo
!
!--Now put the spherical layer of He particles around the WD
!
  rdown = sqrt(rmax) + 0.005
  rup   = 0.05
  call srand(time())
  do p=1,nHe
    rad   = rdown + rand()*rup
    theta = 2.0*pi*rand()
    phi   = acos(1.0-2.0*rand())
    xc = rad*sin(phi)*cos(theta)
    yc = rad*sin(phi)*sin(theta)
    zc = rad*cos(phi)
    xyzhm_He(1,p) = xc
    xyzhm_He(2,p) = yc
    xyzhm_He(3,p) = zc
!
    vxyzut_He(1,p) = 0.0
    vxyzut_He(2,p) = 0.0
    vxyzut_He(3,p) = 0.0
    vxyzut_He(4,p) = 1.0
    vxyzut_He(5,p) = 1.0d7
!
    ka1_He(p)     = ka_min
    fh_He(p)      = 1.
    partype_He(p) = 1
    star_He(p)    = 1
    step_He(p)    = maxstep
!
    xyzhm_He(5,p) = mHe/nHe
    rho_He(p)     = 0.2
    xyzhm_He(4,p) = hfact*(xyzhm_He(5,p)/rho_He(p))**0.333
    do k=1,nel
      xss_He(k,p) = 0.0
    enddo
    xss_He(2,p) = 1.0
  enddo
!
!--Allocate the dimension of the data structures
!
   ntot = npart + nHe
   ALLOCATE (data_array(ntot), STAT = AllocateStatus)
   if (AllocateStatus /= 0) THEN
     write(*,*) "Not enough memory for data_array!!!"
     STOP
   endif
!
   ALLOCATE (cdata_array(ntot), STAT = AllocateStatus)
   if (AllocateStatus /= 0) THEN
     write(*,*) "Not enough memory for data_array!!!"
     STOP
   endif
!
!--Now, load CO data into I/O arrays
!
   do p = 1,npart
     data_array(p)%ioxyzhm(1:ndim+2)  = xyzhm(1:ndim+2,p)
     data_array(p)%iovxyzut(1:ndim+2) = vxyzut(1:ndim+2,p)
     data_array(p)%iorho              = rho(p)
     data_array(p)%ioka1              = ka1(p)
     data_array(p)%iostep             = step(p)
     data_array(p)%iopartype          = partype(p)
     data_array(p)%iostar             = star(p)
     cdata_array(p)%comp(1:nel)       = xss(1:nel,p)
   enddo
!
!--And now, load He layer data into I/O arrays
!
   do p = 1,nHe
     data_array(p+npart)%ioxyzhm(1:ndim+2)  = xyzhm_He(1:ndim+2,p)
     data_array(p+npart)%iovxyzut(1:ndim+2) = vxyzut_He(1:ndim+2,p)
     data_array(p+npart)%iorho              = rho_He(p)
     data_array(p+npart)%ioka1              = ka1_He(p)
     data_array(p+npart)%iostep             = step_He(p)
     data_array(p+npart)%iopartype          = partype_He(p)
     data_array(p+npart)%iostar             = star_He(p)
     cdata_array(p+npart)%comp(1:nel)       = xss_He(1:nel,p)
   enddo
!
!--Now, write data into files
!
   num    = 2
   Omega0 = 0.0
   nbody  = ntot
   ndead  = 0
   nbody1 = ntot
   nbody2 = 0
   write(filename,'(a,i4.4,a)') 'bodi',num,'.out'
   open (1, file=filename, form='unformatted', status='new')
   write(1) tnow, ntot, nbody, ndead, nbody1, nbody2, Omega0
   write(1) data_array
   close(1)
!
   write(filename,'(a,i4.4,a)') 'comp',num,'.out'
   open (2, file=filename, form='unformatted', status='new')
   write(2) cdata_array
   close(2)
!
end subroutine putlayer
