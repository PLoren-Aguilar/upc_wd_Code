      SUBROUTINE dissolve
!============================================================
!  This subroutine disolves the CO particles inside the He
!  layer
!
!  Last revision: 04/March/2017
!============================================================
!
!--Load modules
!
      USE mod_parameters, ONLY : ndim
      USE mod_commons,    ONLY : xyzhm, rho, partype, nb, nbody, rank,  &
                                 factor, npart, ierr, vxyzut
      USE mod_functions,  ONLY : wk
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
#ifdef MPI
      INCLUDE 'mpif.h'
#endif 
!
!--Local variables
!
      INTEGER, DIMENSION(npart) :: isentarray, ireceivearray
      REAL, DIMENSION(2*npart)  :: sentarray, receivearray
      INTEGER, DIMENSION(npart) :: localpartype
      REAL, DIMENSION(npart)    :: localtemp, localmass
      REAL    :: massHe, tempHe, r2, u2p, hp, dr, drx, dry, drz, wkp
      INTEGER :: i, q, qlocal, p, k, m, nCO, nHe
!
!--Figure out the mass of the He(layer) particles
!
      massHe = 0.0
      massloop : DO p=1,nbody
         IF (partype(p) == 1) THEN
           massHe = xyzhm(5,p)
           EXIT massloop
         ENDIF
      ENDDO massloop
!
!--If the CO particle has more He(layer) neighbours than CO neighbours,
!  dissolve the CO particle. 
!
!$OMP PARALLEL DEFAULT(none) shared(nbody,partype,localpartype,localtemp,localmass,nb,rank,factor)  &
!$OMP shared(xyzhm,vxyzut,massHe) private(p,q,qlocal,i,nCO,nHe,tempHe)
!$OMP DO SCHEDULE(runtime)
      partloop : DO qlocal = 1, factor
        q = qlocal + rank*factor
!
        IF ((q >= nbody) .OR. (partype(q) /= 0)) CYCLE
!
        nCO = 0
        nHe = 0
        tempHe = 0.0
        neiloop : DO i = 1,nb(1,qlocal)
          p = nb(i+1,qlocal)
!
          IF (partype(p) == 0) nCO = nCO + 1
          IF (partype(p) == 1) THEN
             nHe = nHe + 1
             tempHe = tempHe + xyzhm(5,p)*vxyzut(5,p)
          ENDIF
!
        ENDDO neiloop
        IF (nHe > nCO) THEN
          localpartype(qlocal) = 1
          localtemp(qlocal)    = tempHe/(nHe*massHe)
          localmass(qlocal)    = massHe
        ELSE
          localpartype(qlocal) = partype(q)
          localtemp(qlocal)    = vxyzut(5,q)
          localmass(qlocal)    = xyzhm(5,q)
        ENDIF
      ENDDO partloop
!$OMP END DO
!$OMP END PARALLEL
!
!--Transfer of particle types MPI processes, if necessary
! 
#ifdef MPI
!$OMP PARALLEL DEFAULT(none) shared(sentarray,localpartype,localtemp,localmass,factor)        &
!$OMP shared(isentarray) private(p)
!$OMP DO SCHEDULE(runtime)
      DO p = 1, factor
        isentarray(p) = localpartype(p)
        sentarray(2*(p-1) + 1) = localtemp(p)
        sentarray(2*(p-1) + 2) = localmass(p)
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
      CALL MPI_ALLGATHER(isentarray,factor,MPI_INTEGER,ireceivearray,   &
                         factor,MPI_INTEGER,MPI_COMM_WORLD,ierr)
      CALL MPI_ALLGATHER(sentarray,2*factor,MPI_DOUBLE_PRECISION,       &
                         receivearray,2*factor,MPI_DOUBLE_PRECISION,    &
                         MPI_COMM_WORLD,ierr)
!
!$OMP PARALLEL DEFAULT(none) shared(nbody,ireceivearray,receivearray,partype,vxyzut,xyzhm)  &
!$OMP private(p)
!$OMP DO SCHEDULE(runtime)
      DO p = 1, nbody
        partype(p)  = ireceivearray(p)
        vxyzut(5,p) = receivearray(2*(p-1)+1)
        xyzhm(5,p)  = receivearray(2*(p-1)+2)
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
#else
!$OMP PARALLEL DEFAULT(none) shared(nbody,localpartype,partype,localpartype,massHe,xyzhm)  &
!$OMP shared(vxyzut) private(p)
!$OMP DO SCHEDULE(runtime)
      DO p = 1, nbody
        xyzhm(5,p)  = massHe
        vxyzut(5,p) = localtemp(p)
        partype(p)  = localpartype(p)
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
#endif MPI

      END SUBROUTINE dissolve
