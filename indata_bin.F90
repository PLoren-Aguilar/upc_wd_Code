      subroutine indata
!===================================================================
!  This subroutine reads the initial data.
!
!  Last revision: 15/March/2015
!===================================================================
!
!--Load modules
!
      use mod_essentials
      use mod_commons
!     
!--Force to declare EVERYTHING
!
      implicit none
!
!--Local variables
!
      real, dimension(npart) :: temporal 
      real :: real_factor,cmx,cmy,cmz
      integer, dimension(npart) :: itemporal 
      integer ::  p, k, nums, m
      character(30) :: filename
!
!--Declare structures in order to read data
!
      type :: outdata
        real, dimension(ndim+2) :: ioxyzhm
        real, dimension(ndim+2) :: iovxyzut
        real :: iorho
        real :: ioka1
        integer(1) :: iopartype
        integer(1) :: iostar
      end type
      type(outdata), dimension(npart) :: data_array
!
      type :: coutdata
        real, dimension(nel) :: comp
      end type
      type(coutdata), dimension(npart) :: cdata_array
!
!--Open data unit
!
      write(filename,'(a,i4.4,a)') 'bodi',nstep_ini,'.out'
      open (1, FILE=filename, form='unformatted')
      write(filename,'(a,i4.4,a)') 'comp',nstep_ini,'.out'
      open (2, FILE=filename, form='unformatted')
!
!--Read bodi data
!
      read(1) tnow, nums, data_array
!
!--Check if the number of particles is consistent
!
      IF (nums /= npart) then
         IF (rank == MASTER) print*, 'nparts NE nbody!!!'
         stop
       endif
!
!--Composition data
!
      read(2) cdata_array
!
!--Close files
!
      close(unit=1)
      close(unit=2)
!
!--Transfer data into the appropriate locations
!
      do p = 1,npart
         xyzhm(1:ndim+2,p)  = data_array(p)%ioxyzhm(1:ndim+2)
         vxyzut(1:ndim+2,p) = data_array(p)%iovxyzut(1:ndim+2)
         rho(p)             = data_array(p)%iorho
         ka1(p)             = data_array(p)%ioka1
         partype(p)         = data_array(p)%iopartype
         star(p)            = data_array(p)%iopartype
         xss(1:nel,p)       = cdata_array(p)%comp(1:nel)
      enddo
!
!--Start up variables
!
      vsigmax = 1.0d0
      cur     = 0.0d0
      div     = 0.0d0
      luminuc = 0.0d0
      enuc    = 0.0d0
      tscdyn  = 0.0d0
      tscnuc  = 0.0d0
!
      eps = 1d30
      do p = 1, npart
         eps = MIN(eps,xyzhm(4,p))
         ilist(p)  = p
         plist(p)  = p
      enddo
      eps  = 1.4d0*2.d0*eps
      eps3 = eps*eps*eps
!
!--Sort particle list to get rid of dead particles
!
      call sort
#ifdef debug
      IF (rank == MASTER) print*, 'sort called'
#endif
!
      end subroutine indata
