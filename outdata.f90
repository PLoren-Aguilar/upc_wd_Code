subroutine outdata
!===================================================================
!  This subroutine writes data into disk.
!
!  Last revision: 15/March/2015
!===================================================================
!
!--Load modules
!
  use mod_parameters, only : ndim, nel
  use mod_commons,    only : xyzhm, rho, vxyzut, ka1, fh,xss, nb,       &
      partype, gxyzu, star, nstep, nout, npart, tnow, RELFLAG, nbody,   &
      nbody1, nbody2, ndead, tend, SIMtype, trelax, Omega0, mode, step, &
      dtmaxin, dt0in, reset, nstep_infile, fixx1, fixy1, fixz1, fixx2,  &
      fixy2, fixz2
!
!--Force do declare EVERYTHING
!
  implicit none
!
!--Local variables
!
  integer :: p, k, num
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
  type(iodata), dimension(npart) :: data_array
!
  type :: ciodata
     real(4), dimension(nel) :: comp
  end type ciodata
  type(ciodata), dimension(npart) :: cdata_array
!
!--Load data into data arrays
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
!--Open bodi unit and write data
!
  !num = nstep/nout
  num = nstep_infile+nstep/nout
  write(filename,'(a,i4.4,a)') 'bodi',num,'.out'
  open (1, file=filename, form='unformatted', status='new')
  write(1) tnow, npart, nbody, ndead, nbody1, nbody2, Omega0
  write(1) data_array
  close(1)
!
!--Open comp unit and write data
!
  write(filename,'(a,i4.4,a)') 'comp',num,'.out'
  open (2, file=filename, form='unformatted', status='new')
  write(2) cdata_array
  close(2)
!
!--Open treepars file and write data
!
  open (3, file='treepars', status='unknown')
  write(3,'(i4.4,a)') num,' :: nstep'
  write(3,'(1(1ES12.4),a)') tend,' :: final time in first run'
  write(3,'(1(1ES12.4),a)') dtmaxin,' :: maximum time step'
  write(3,'(1(1ES12.4),a)') dt0in,' :: initial time-step (in units of dtmax)'
  write(3,'(i4.4,a)') nout, ' :: Number of steps before output'
  write(3,'(i4.4,a)') SIMtype,' :: SIMtype (1:Single star, 2: Binary system)'
  write(3,'(L1,a)') RELFLAG,' :: RELFLAG'
  write(3,'(L1,a)') reset,' :: resetflag'
  write(3,'(f8.4,a)') trelax,' :: TRELAX'
  write(3,'(3(1ES12.4),a)') fixx1,fixy1,fixz1,' :: fixed position for the first star'
  write(3,'(3(1ES12.4),a)') fixx2,fixy2,fixz2,' :: fixed position for the second star'
  close (3)
! 
end subroutine outdata
