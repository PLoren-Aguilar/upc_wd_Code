subroutine inparams           
!=====================================================
!
!  This subroutine reads the initial code paameters
!
!  Last revision: 9/April/2019 
!
!=====================================================
!
!--Load modules
!
  use mod_parameters, only : nel
  use mod_commons   , only : nstep_infile, tend, nout, SIMTYPE, trelax, &
                             RELFLAG, nstep, aion, zion, dtmaxin, dtmax,&
                             dt0, dt0in, reset, fixx1, fixy1, fixz1,    &
                             fixx2, fixy2, fixz2
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--Local variables
!
  integer      :: k, NNNZ
  character(5) :: ONC
!
!--Read parameters
!
  open (1,FILE='treepars',status='old')
  read(1,*) nstep_infile
  read(1,*) tend
  read(1,*) dtmaxin
  read(1,*) dt0in
  read(1,*) nout
  read(1,*) SIMTYPE
  read(1,*) RELFLAG
  read(1,*) reset
  read(1,*) trelax
  read(1,*) fixx1,fixy1,fixz1
  read(1,*) fixx2,fixy2,fixz2
  close (1)
!
!--Integration parameters
!
  dtmax = 0.0
  dt0   = 0.0
!
!--Read isotopes (neglecting gamma photons)
!
  open(4,file='ZHELI.10b',status='old')
  aion(1) = 1.0
  zion(1) = 1.0
  isoloop : do k=1,nel-2
    read(4,8050) NNNZ,ONC,aion(k+1),zion(k+1)
  enddo isoloop
8050 format(I4,1X,A5,1X,F4.0,1X,F4.0)
  close(unit=4)
!
!--Read EOS data
!
  call read_helm_table
!
!--Read nuclear tables
!
!  call RPARAM
!  call RNETWORK
end subroutine inparams
