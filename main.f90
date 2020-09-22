        PROGRAM main
!----------------------------------------------------------------
!
!     PROGRAM FOR HYDRODYNAMICAL CALCULATIONS USING SPH AND A 
!     BINARY TREE FOR GRAVITY AND NEIGHBOR SEARCHING. THE PRESENT
!     VERSION SUPPORTS OpenMP+MPI PARALLELIZATION

!     The code has tree working modes:
!
!     * single: Generates a initial WD model by solving the 
!               Lane-Embden equation. A parameter file called
!               genparams must be provided.
!
!     * binary: Constructs a binary initial model from two WD
!               models. The initial parameters for the binary
!               system will be entered by hand
!
!     * run: Starts the main SPH code. The code will be run in
!            different modes according to the parameters provided
!            through the parameters file called treepars 
!
!     Last revision: 15/March/2015
!
!----------------------------------------------------------------
!
!--Load modules
!
        USE mod_commons, ONLY : mode
!
!--Force to declare EVERYTHING
!
        IMPLICIT NONE
!
!--Setup running mode
!
        CALL getarg(1,mode(1))
        IF (mode(1).NE.'single'.AND.mode(1).NE.'binary'.AND.mode(1).NE.'ejecta'.AND.  &
            mode(1).NE.'run'.AND.mode(1).NE.'analysis'.AND.mode(1).NE.  &
            'putlayer'.AND.mode(1).NE.'separate'.AND.mode(1).NE.'plot') THEN
            PRINT*, 'BAD RUNNING MODE!'
            STOP
        ENDIF
!
        IF (mode(1).EQ.'single')   CALL genera_single
        IF (mode(1).EQ.'binary')   CALL genera_bin
        IF (mode(1).EQ.'run')      CALL sph
        IF (mode(1).EQ.'analysis') CALL analysis
        IF (mode(1).EQ.'putlayer') CALL putlayer
        IF (mode(1).EQ.'separate') CALL separate
        IF (mode(1).EQ.'plot')     CALL plot
!
        END
