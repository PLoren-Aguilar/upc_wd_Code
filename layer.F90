      SUBROUTINE layer
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
                                 rhoG, nbody, rank, factor, npart, ierr
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
      REAL, DIMENSION(ndim+2, npart)  :: newvxyzut
      REAL, DIMENSION((ndim+2)*npart) :: sentarray, receivearray
      REAL, DIMENSION(ndim) :: vxyz
      REAL :: massp, r2, u2p, hp, scalar, dr, drx, dry, drz, dvx,       &
                 dvy, dvz, rhop, rhoq, chi, wkp, chichi, tCO, uCO, temp,&
                 normq
      REAL, PARAMETER :: norma = 10./9.
      INTEGER ::  i, q, qlocal, p, k, m
      INTEGER, PARAMETER :: nu=1
!
!--Mean velocity flow
!
      newvxyzut = 0.0
!$OMP PARALLEL DEFAULT(none) shared(xyzhm,vxyzut,rho,rank,partype,nb,norm,factor) & 
!$OMP shared(rhoG,nbody,newvxyzut) private(i,p,q,rhoq,rhop,hp,scalar,chi,chichi)  &
!$OMP private(dr,drx,dry,drz,dvx,dvy,dvz,massp,u2p,wkp,qlocal,vxyz,tCO,uCO,temp,normq)
!$OMP DO SCHEDULE(runtime)
      DO 10 qlocal = 1, factor
         q = qlocal + rank*factor
         IF (q.GT.nbody) GOTO 15
!
!--Go only through active He particles
!
         IF (partype(q).NE.1) GOTO 15
!
!--Store for later use q-particle related variables 
!
         normq   = 1.0
         rhoq    = rho(q)
!
!--Loop over active CO neighbours
!
         vxyz= 0.0d0
         tCO = 0.0d0
         uCO = 0.0d0
!
         DO 20 i = 1,nb(1,qlocal)
            p = nb(i+1,qlocal)
            IF (partype(p).NE.0) GOTO 25
!
!--Calculate unit vectors and velocity scalar product
!
            drx = xyzhm(1,q) - xyzhm(1,p)
            dry = xyzhm(2,q) - xyzhm(2,p)
            drz = xyzhm(3,q) - xyzhm(3,p)
            dr  = SQRT(drx*drx + dry*dry + drz*drz)
            IF (dr.GT.0.0) THEN
               drx = drx/dr
               dry = dry/dr
               drz = drz/dr
            ELSE
               drx = 0.0
               dry = 0.0
               drz = 0.0
            ENDIF
!
            dvx = vxyzut(1,q) - vxyzut(1,p)
            dvy = vxyzut(2,q) - vxyzut(2,p)
            dvz = vxyzut(3,q) - vxyzut(3,p)
            scalar = dvx*drx + dvy*dry + dvz*drz
!
!--p-particle variables
!
            rhop   = rho(p)
            temp   = vxyzut(5,p)
            hp     = xyzhm(4,p)
            massp  = xyzhm(5,p)
!
!--Calculate kernel related variables
!
            u2p = dr*dr/(hp*hp)
            wkp = wk(u2p)/(hp**ndim)
            !wkp = norma*wk(u2p)*u2p/(hp**ndim)
!
!--Density coefficient
!
            chi = 1.0/(rhoq + rhop)
!
!--Make He layer follow mean CO flow
!
            vxyz(1) = vxyz(1) - nu*massp*chi*scalar*drx*wkp
            vxyz(2) = vxyz(2) - nu*massp*chi*scalar*dry*wkp
            vxyz(3) = vxyz(3) - nu*massp*chi*scalar*drz*wkp
!
!--Make He layer to adopt mean thermal equilibrium with CO
!
            uCO = uCO + 0.5*nu*massp*chi*scalar*scalar*wkp
            tCO = tCO + (massp/rhop)*temp*wkp
25          CONTINUE
20       ENDDO      ! End of neighbours loop
!
!--Define "switch" function
!
         IF (rhoG(q).EQ.0.0) THEN
            chichi = 0.0
         ELSE
            chichi = EXP(-rhoq/rhoG(q))
            !chichi = DEXP(-(xyzhm(5,p)/xyzhm(5,q))*rhoq/rhoG(q))
         ENDIF
!
!--Update velocity and thermal energy
!
         DO k=1,3
            newvxyzut(k,qlocal) = (1.-chichi)*vxyzut(k,q) +             &
                               chichi*vxyz(k)/normq
         ENDDO
         newvxyzut(4,qlocal) =  (1.-chichi)*vxyzut(4,q) +               &
                               chichi*uCO/normq
         newvxyzut(5,qlocal) =  (1.-chichi)*vxyzut(5,q) +               &
                               chichi*tCO/normq
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
!--Tranfer of accelerations across MPI processes, if necessary
! 
#ifdef MPI
!$OMP PARALLEL DEFAULT(none) shared(sentarray,newvxyzut,factor)        &
!$OMP private(p,k)
!$OMP DO SCHEDULE(runtime)
      DO p = 1, factor
         DO k = 1, ndim+2
            sentarray((p-1)*(ndim+2)+k) = newvxyzut(k,p)
         ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
      CALL MPI_ALLGATHER(sentarray,(ndim+2)*factor,MPI_DOUBLE_PRECISION,&
                         receivearray,(ndim+2)*factor,                  &
                         MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
!
!$OMP PARALLEL DEFAULT(none) shared(nbody,receivearray,vxyzut,partype)  &
!$OMP private(p,k)
!$OMP DO SCHEDULE(runtime)
      DO p = 1, nbody
         DO k = 1, ndim+2
            IF (partype(p).EQ.1) THEN
               vxyzut(k,p)=receivearray((p-1)*(ndim+2)+k)
            ENDIF
         ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
#else
!$OMP PARALLEL DEFAULT(none) shared(nbody,newvxyzut,vxyzut,partype)  &
!$OMP private(p,k)
!$OMP DO SCHEDULE(runtime)
      DO p = 1, nbody
         DO k = 1, ndim+2
            IF(partype(p).EQ.1) vxyzut(k,p) = newvxyzut(k,p)
         ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
#endif MPI
!
      END SUBROUTINE layer
