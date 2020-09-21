subroutine diagnostics
!===========================================================
!  This subroutine outputs all the necessary information
!
!  Last revision: 15/March/2015
!===========================================================
!
!--Load modules
!
  use mod_parameters, only : MASTER, maxstep, log2, maxstep
  use mod_commons,    only : globnmin, globnmax, globnvec, globmax,   &
  globdone, nstep, nout, rank, nCO, nHe, ndead, step, partype, nbody, &
  dtmax, nactive, nstep_infile
!
!--Force to declare EVERYTHING
!
  IMPLICIT NONE
!
!--Local variables
!
  integer, dimension(30) :: nbin
  integer                :: n, p, i
  real :: dt
!
  if (rank == MASTER) then
    if (mod(nstep,nout) == 0) then
      CALL energy
      CALL outdata

      print*, ''
      print*, '============================================'
      write(*,'(a,i4.4,a)') 'bodi', nstep/nout,'.out written'
      print*, ''
      write(*,'(a,i6,a,i6,a,i6)') 'NMIN:',globnmin, ' NMAX:',       &
            globnmax,' AVG:',globnvec
      write(*,'(a,i7)') 'Notdone=',globdone
      write(*,'(a,1pe12.4,a)') 'Maxdiff',globmax*100.,' %'
      print*, ''
      write(*,'(a,i6,a,i6,a,i6)') 'nCO=',nCO,' nHe=',nHe,           &
            ' ndead=',ndead
      print*, '============================================'
      print*, ''
    endif
  endif
!
!--Print timestepping information
!
  if ((rank == MASTER) .and. dtmax /= 0) then
    nbin = 0
    do p=1,nbody
      if (partype(p) == 2) cycle
      dt      = dtmax*float(step(p))/float(maxstep)
      n       = int(log10(dtmax/dt)/log10(2.)) + 1          
      nbin(n) = nbin(n) + 1
    enddo
    write(*,'(a)') '         Timestep distribution'
    write(*,'(a)') '========================================'
    do i=1,30
      write(*,'(i2,i11,1pe12.4,i7)') i, maxstep/2**(i-1),dtmax/2**(i-1), nbin(i)
    enddo
  endif
!
end subroutine diagnostics
