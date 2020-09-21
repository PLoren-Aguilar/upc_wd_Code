program main
!----------------------------------------------------------------
!
!     Program for hydrodynamicsl calculations using SPH and a 
!     binary tree for gravity and neighbour searching. The present
!     version supports OpenMP+MPI parallelization, the presence
!     of a He layer and invividual time-steps
!
!     AUTHORS: * Jordi Isern Vilaboy
!              * Enrique Garcia-Berro
!              * Josep Guerrero
!              * Pablo Loren-Aguilar
!              * Alba Gutierrez-Pedemonte
!              * Gabriela Aznar-Siguan
!
!     The code has several working modes:
!
!     * single: Generates a initial WD model by solving the 
!               Lane-Embden equation. A parameter file called
!               genparams must be provided.
!
!     * binary: Constructs a binary initial model from two WD
!               models. The initial parameters for the binary
!               system will be entered by hand
!
!     * run:    Starts the main SPH code. The code will be run in
!               different modes according to the parameters provided
!               through the parameters file called treepars 
!
!     * analysis: Calls the subroutine "analysis", useful to perform
!                 any analysis to a particular time-step
!
!     * putlayer: Puts a He layer on the top of a (ideally relaxed) WD
!
!     * separate: Used to separate a pair of WDs on a binary system
!
!     * plot:     Transforms the binary output into an ASCII file
!                 suitable for plotting
!
!     Last revision: 6/April/2019
!
!----------------------------------------------------------------
!
!--Load modules
!
  use mod_commons, only : mode
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--Setup running mode
!
  call getarg(1,mode(1))
  if ( (mode(1) /= 'single')   .and. (mode(1) /= 'binary')   .and. &
       (mode(1) /= 'run')      .and. (mode(1) /= 'analysis') .and. &
       (mode(1) /= 'putlayer') .and. (mode(1) /= 'separate') .and. &
       (mode(1) /= 'plot') ) then
    print*, 'BAD RUNNING MODE!'
    stop
  endif

  if (mode(1) == 'single')   call genera_single
  if (mode(1) == 'binary')   call genera_bin
  if (mode(1) == 'run')      call sph
  if (mode(1) == 'analysis') call analysis
  if (mode(1) == 'putlayer') call putlayer
  if (mode(1) == 'separate') call separate
  if (mode(1) == 'plot')     call plot

end program main
