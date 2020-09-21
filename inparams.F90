        SUBROUTINE inparams           
!=====================================================
!
!  This subroutine reads the initial code paameters
!
!  Last revision: 15/March/2015
!
!=====================================================
!
!--Load modules
!
        USE mod_parameters, ONLY : nel
        USE mod_commons, ONLY : nstep_infile, tend, nout, SIMTYPE, trelax, &
                                RELFLAG, nstep, aion, zion
!
!--Force to declare EVERYTHING
!
        IMPLICIT NONE
!
!--Local variables
!
        INTEGER      :: k, NNNZ
        CHARACTER(5) :: ONC
!
!--Open parameter file
!
        OPEN (1,FILE='treepars',STATUS='old')
!
!--Read parameters
!
        READ(1,*) nstep_infile
        READ(1,*) tend
        READ(1,*) nout
        READ(1,*) SIMTYPE
        READ(1,*) RELFLAG
        READ(1,*) trelax   
!
!--Initializations        
!
!        nstep = nout*(nstep_ini-1)
!
!--Close parameter file
!
        CLOSE (1)
!
!--Integration parameters
!
        !dtfact   = 0.3d0
        !dtmp_max = 0.05d0
!
!--Read isotopes (neglecting gamma photons)
!
        OPEN(4,FILE='ZHELI.10b',STATUS='old')
!
        aion(1) = 1.0
        zion(1) = 1.0
        isoloop : DO k=1,nel-2
           READ(4,8050) NNNZ,ONC,aion(k+1),zion(k+1)
        ENDDO isoloop
8050    FORMAT(I4,1X,A5,1X,F4.0,1X,F4.0)
!
        CLOSE(UNIT=4)
!
!--Read EOS data
!
        CALL read_helm_table
!
!--Read nuclear tables
!
        CALL RPARAM
        CALL RNETWORK
!
        END SUBROUTINE inparams
