      SUBROUTINE norm_layer
!============================================================
!  This subroutine calculates the hydrodynamical 
!  accelerations
!
!  Last revision: 15/March/2015
!============================================================
!
!--Load modules
!
      USE mod_parameters, ONLY : ndim
      USE mod_commons, ONLY    : xyzhm, vxyzut, rho, partype, nb, norm, &
                                 rhoG, rank, factor, nbody, ierr, npart
      USE mod_functions, ONLY  : wk
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
#ifdef MPI
      REAL, DIMENSION(2*npart) :: sentarray, receivearray
#endif 
      REAL, DIMENSION(ndim) :: vxyz
      REAL :: massp, r2, u2p, hp, scalar, dr, drx, dry, drz, dvx,    &
                 dvy, dvz, rhop, rhoq, chi, wkp, DD_wkp
      REAL, PARAMETER :: norma=10./9.
      INTEGER ::  i, q, qlocal, p, k, m
      INTEGER, PARAMETER :: nu=3
!
!--Calculate normalization
!
!$OMP PARALLEL DEFAULT(none) shared(xyzhm,rho,rank,partype,nb,norm)     & 
!$OMP shared(rhoG,factor,nbody) private(i,p,q,rhoq,rhop,hp,massp,u2p,wkp) &
!$OMP private(DD_wkp,qlocal,dr,drx,dry,drz)
!$OMP DO SCHEDULE(runtime)
      DO 10 qlocal = 1, factor
         q = qlocal + rank*factor
!
         IF (q.GT.nbody) GOTO 15
!
!--Go only through active He particles
!
         IF (partype(q).NE.1)  GOTO 15
!
!--Store for later use q-particle related variables 
!
         rhoq   = rho(q)
!
!--Loop over active CO neighbours
!
         rhoG(q) = 0.0
         norm(q) = 0.0
         DO 20 i = 1,nb(1,qlocal)
            p = nb(i+1,qlocal)
!
            IF (partype(p).NE.0)  GOTO 25
!
!--p-particle variables
!
            rhop   = rho(p)
            hp     = xyzhm(4,p)
            massp  = xyzhm(5,p)
!
!--Calculate kernel related variables
!
            drx = xyzhm(1,q) - xyzhm(1,p)
            dry = xyzhm(2,q) - xyzhm(2,p)
            drz = xyzhm(3,q) - xyzhm(3,p)
            dr = SQRT(drx*drx + dry*dry + drz*drz)
!
            u2p    = dr*dr/(hp*hp)
            DD_wkp = norma*wk(u2p)*u2p/(hp**ndim)
            wkp    = wk(u2p)/(hp**ndim)
!
!--Calculate normalization factor
!
            norm(qlocal) = norm(qlocal) + (massp/rhop)*DD_wkp
            rhoG(qlocal) = rhoG(qlocal) + massp*wkp
25          CONTINUE
20       ENDDO      ! End of neighbours loop
15       CONTINUE
10    ENDDO      !End of particles loop
!$OMP END DO
!$OMP END PARALLEL
!
!--Force synchronization
!
#ifdef MPI
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
#endif
!
!--Tranfer of variables across MPI processes, if necessary
! 
#ifdef MPI
!$OMP PARALLEL DEFAULT(none) shared(sentarray,norm,rhoG,factor)         &
!$OMP private(p)
!$OMP DO SCHEDULE(runtime)
      DO p = 1, factor
         sentarray((p-1)*2 + 1) = norm(p) 
         sentarray((p-1)*2 + 2) = rhoG(p) 
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
      CALL MPI_ALLGATHER(sentarray,2*factor,MPI_DOUBLE_PRECISION,       &
                         receivearray,2*factor,MPI_DOUBLE_PRECISION,    &
                         MPI_COMM_WORLD,ierr)
!
!$OMP PARALLEL DEFAULT(none) shared(nbody,receivearray,norm,rhoG) &
!$OMP private(p)
!$OMP DO SCHEDULE(runtime)
      DO p = 1, nbody
         norm(p) = receivearray((p-1)*2 + 1)
         rhoG(p) = receivearray((p-1)*2 + 2)
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
#endif MPI
!
      END SUBROUTINE norm_layer
