subroutine relax
!===================================================================
!
!  This subroutine relaxes the system whenever necessary.
!
!  Last revision: 15/March/2015
!
!===================================================================
!
!--Load modules
!
  use mod_parameters, only : MASTER
  use mod_commons,    only : RELFLAG, rank, SIMTYPE, nstep, nout
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--Local variables
!
  logical :: FINISH
!
!--In case of binary relaxation, slowly approach stars 
!  until Roche-Lobe overflow is achieved
!
  if (SIMTYPE == 2) then
! Note that after removing the distance, the systems needs to
! relax, so we will only reduce the distance every nout time steps
    call approach(0.0,FINISH)
!
    if (FINISH.eqv..true.) then
      print*, 'Binary relaxation finished'
      stop
    endif
  endif           ! End of SIMTYPE if
!
!--Remove any center of mass displacement
!
  call relax_frame
#ifdef debug
  if (rank == MASTER) print*, 'relax_frame called'
#endif
!
end subroutine relax
