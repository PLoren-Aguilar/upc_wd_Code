       SUBROUTINE iter_rhoh
!=============================================================================
!      THIS SUBROUTINE CALCULATES NEW DENSITY AND SMOOTHING LENGTH ITERATIVELY 
!
!      Last revision: 15/March/2015
!=============================================================================
!
!--Load modules
!
       USE mod_commons,    ONLY : xyzhm, rho, fh, div, divt, cur, dhdt, &
                                  globmax, globdone, rho, partype, nb,  &
                                  npart, factor, rank, nbody, ierr, size 
       USE mod_functions,  ONLY : wk, dk
       USE mod_parameters, ONLY : rhomin, rhomax, hmin, hmax, hfact,    &
                                  nel, ndim, MASTER
!
!--Force to declare everything 
!
       IMPLICIT NONE
#ifdef MPI
       INCLUDE 'mpif.h'
#endif
!
!--Local variables
!
       REAL, DIMENSION(7*npart) :: sentarray, receivearray
       REAL, DIMENSION(nbody) :: haux
       REAL  :: mp, hp, hnew, rhoh, rhop, drhopdh, dhdrhoh, omega,      &
                fp, funch, dfunchdh, diff, max_diff, drhodtp, divp,     &
                divtp, curp
       REAL, PARAMETER :: tolerance=1.0e-3
       INTEGER  :: i, k, p, m, plocal, notdone
       INTEGER, PARAMETER :: itermax=250 
       LOGICAL  :: converged
!
!--Calculate neighbour lists and gravitational forces
!
       CALL treeconst
       !DO plocal = 1, factor
       !   DO i = 1, nb(1,plocal)
       !      !PRINT*, plocal,i,nb(1,plocal),nb(i+1,plocal)
       !      IF (nb(i+1,plocal).EQ.0) THEN
       !         !PRINT*, 'Aqui tambien!!!',factor,rank,plocal,i,partype(plocal+rank*factor)
       !         PRINT*, 'Aqui tambien!!!'
       !         STOP
       !      ENDIF
       !   ENDDO
       !ENDDO
#ifdef debug
      IF (rank == MASTER) THEN
         PRINT*, 'treeconst called',rank
      ENDIF

#endif
!
!--Loop through particles
!
       notdone  = 0
       max_diff = 0
!$OMP PARALLEL DEFAULT(none) shared(xyzhm,rho,fh,dhdt,div,divt,cur)     &
!$OMP shared(rank,haux,partype,factor,nbody,nb)                 &
!$OMP private(m,p,plocal,i,mp,hp,rhop,dhdrhoh,omega,fp,funch,dfunchdh)  &
!$OMP private(hnew,divp,divtp,curp,drhopdh,drhodtp,rhoh,diff,converged) & 
!$OMP reduction(MAX:max_diff) reduction(+:notdone)
!$OMP DO SCHEDULE(runtime)
       partloop : DO plocal = 1, factor
          p = plocal + rank*factor
          IF ((p > nbody) .OR. (partype(p) == 2)) CYCLE
!
!--Old h and mass
!
          hp = xyzhm(4,p)
          mp = xyzhm(5,p)
!
!--Main NR loop
!
          converged = .false.
          iterloop : DO i=1,itermax
!
!--Calculate density using h, for t
!
             CALL new_rho(p,plocal,mp,hp,rhop,drhopdh,drhodtp,divp,     &
                          divtp,curp)
!
!--rho compatible with the old h
!
             rhoh    = mp*(hfact/hp)**ndim
             dhdrhoh = - hp/(ndim*rhop)
             omega   = 1.0 - dhdrhoh*drhopdh
!
!--If Newton-Raphson does not converge, we keep initial data
!
             IF (omega < 1.0e-16) THEN     
               fp = 1.0
               EXIT iterloop
             END IF
!
!--function we want to find the zeros from
!
             funch    = rhoh - rhop
             dfunchdh = omega/dhdrhoh
             fp       = 1.0/omega
!
!--Calculate new h
!
             hnew = hp - funch/dfunchdh
!
!--Overwrite if iteration is going wrong
!
             IF ((hnew < 0.0) .OR. (fp < 1.0e-6)) THEN
                hnew = hfact*(mp/rhop)**(1.0/REAL(ndim))
             END IF
!
!--Limit the change in h... just in case
!
             IF (hnew > 1.2*hp) hnew = 1.2*hp
             IF (hnew < 0.8*hp) hnew = 0.8*hp
             IF (hnew > hmax)   hnew = hmax
             IF (hnew < hmin)   hnew = hmin
!               
!--Look at convergence
!
             diff = ABS(rhop-rhoh)/rhop
             !diff = ABS(hnew-hp)/hp
             IF ((diff < tolerance) .AND. (omega > 0.)) converged = .true.
!
!--Update h
!
             hp = hnew
! 
!--Abandon the loop if finished
!
             IF (converged.EQV..true.) EXIT iterloop 
          ENDDO  iterloop 
!
!--Check change in h and update particle values
!
          IF (diff > tolerance) THEN
             notdone = notdone + 1
          ELSE
             max_diff = MAX(max_diff,diff)
          ENDIF
!
          rho(plocal)  = rhop
          fh(plocal)   = fp
          dhdt(plocal) = dhdrhoh*drhodtp*fp
          haux(plocal) = hp
          div(plocal)  = divp
          divt(plocal) = divtp
          cur(plocal)  = curp
       ENDDO partloop
!$OMP END DO
!$OMP END PARALLEL
!
!--Force synchronization
!
#ifdef MPI     
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
!--Transfer variables across MPI processess if necessary
!
#ifdef MPI
!$OMP PARALLEL DEFAULT(none) shared(sentarray,rho,fh,dhdt,haux,div,divt)&
!$OMP shared(cur,factor) private(p)
!$OMP DO SCHEDULE(runtime)
      DO p = 1, factor
         sentarray((p-1)*7 + 1) = rho(p)
         sentarray((p-1)*7 + 2) = fh(p)
         sentarray((p-1)*7 + 3) = dhdt(p)
         sentarray((p-1)*7 + 4) = haux(p)
         sentarray((p-1)*7 + 5) = div(p)
         sentarray((p-1)*7 + 6) = divt(p)
         sentarray((p-1)*7 + 7) = cur(p)
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!      
      CALL MPI_ALLGATHER(sentarray,7*factor,MPI_DOUBLE_PRECISION,       &
                         receivearray,7*factor,MPI_DOUBLE_PRECISION,    &
                         MPI_COMM_WORLD,ierr)
!
!$OMP PARALLEL DEFAULT(none) shared(receivearray,rho,fh,dhdt,haux,div)  &
!$OMP shared(divt,cur,xyzhm,nbody) private(p)
!$OMP DO SCHEDULE(runtime)
      DO p = 1, nbody
         rho(p)     = receivearray((p-1)*7 + 1)
         fh(p)      = receivearray((p-1)*7 + 2)
         dhdt(p)    = receivearray((p-1)*7 + 3)
         xyzhm(4,p) = receivearray((p-1)*7 + 4)
         div(p)     = receivearray((p-1)*7 + 5)
         divt(p)    = receivearray((p-1)*7 + 6)
         cur(p)     = receivearray((p-1)*7 + 7)
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
#else
!$OMP PARALLEL DEFAULT(none) shared(nbody,xyzhm,haux,partype) private(p)
!$OMP DO SCHEDULE(runtime)
      DO p = 1, nbody
         IF (partype(p) /= 2) xyzhm(4,p) = haux(p)
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
#endif MPI
!
!--Print status
!
#ifdef MPI
      CALL MPI_REDUCE(max_diff, globmax, 1, MPI_DOUBLE_PRECISION,       &
                      MPI_MAX, MASTER, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(notdone, globdone, 1, MPI_INTEGER, MPI_SUM,       &
                      MASTER, MPI_COMM_WORLD, ierr)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) ! Force synchronization
#else
      globmax  = max_diff
      globdone = notdone
#endif
!
      !DO p = 1,nbody
      !   IF (rho(p) == 0.0) PRINT*, p,rho(p),partype(p)
      !ENDDO

       END SUBROUTINE iter_rhoh
!
       SUBROUTINE new_rho(q,qlocal,mq,hq,rhofunc,gradh,drhodt,divvfunc, &
                          divtfunc,rotn)
!=====================================================================
!
!  This subroutine calculates the new density value for
!  particle p
!
!=====================================================================
!
!--Load modules
!
      USE mod_commons, ONLY    : xyzhm, vxyzut, nb, partype
      USE mod_functions, ONLY  : wk, dk, dk_h
      USE mod_parameters, ONLY : rhomin, rhomax, ndim
!
!--Force do declare EVERYTHING
!
      IMPLICIT NONE
!
!--I/O variables
!
      REAL, INTENT(IN)  :: mq, hq
      REAL, INTENT(OUT) :: rhofunc, gradh, divvfunc, divtfunc, rotn,    &
                              drhodt
      INTEGER, INTENT(IN)  :: q, qlocal
!
!--Local variables
!
      REAL, DIMENSION(ndim) :: vel0, pos0, dv, dp, rot
      REAL :: massp, dterm, r2, u2q, u2p, hp, hq3, hp3, hq5, hq2,       &
      hp2, hp5, vijrij, vijj, dkhp, dkq, dkp, wkq, wkp, dkhq
      INTEGER :: p, k, m, i
!
!--Variables start up
!
      pos0(1:ndim) = xyzhm(1:ndim,q)
      vel0(1:ndim) = vxyzut(1:ndim,q)
      rhofunc   = 0.0
      divvfunc  = 0.0
      divtfunc  = 0.0
      rot       = 0.0
      drhodt    = 0.0
      gradh     = 0.0
!
      hq2 = hq*hq
      hq3 = hq2*hq
      hq5 = hq3*hq2
!
      neiloop : DO i = 1,nb(1,qlocal)
         p    = nb(i+1,qlocal)
         IF (partype(q) /= partype(p)) CYCLE neiloop
         hp   = xyzhm(4,p)
         hp2  = hp*hp
         hp3  = hp2*hp
         hp5  = hp3*hp2
!
         massp = xyzhm(5,p)
!
         vijrij = 0.0
         dimloop : DO k = 1,ndim
            dv(k)  = vel0(k) - vxyzut(k,p)
            dp(k)  = pos0(k) - xyzhm(k,p)
            vijrij = vijrij  + dv(k)*dp(k)
         ENDDO dimloop
         r2 = dp(1)*dp(1)+dp(2)*dp(2)+dp(3)*dp(3)
!
         u2q   = r2/hq2
         u2p   = r2/hp2
         wkq   = wk(u2q)/hq3            
         wkp   = wk(u2p)/hp3            
         dkq   = dk(u2q)/hq5
         dkp   = dk(u2p)/hp5
         dkhq  = dk_h(u2q)/(hq2*hq2)
         dkhp  = dk_h(u2p)/(hp2*hp2)
         dterm = massp*0.5*(dkq+dkp)
!
         rhofunc  = rhofunc  + massp*wkq
         gradh    = gradh    + massp*dkhq
!
         drhodt   = drhodt + massp*vijrij*dkq
         divvfunc = divvfunc + vijrij*dterm
         divtfunc = divtfunc + (vxyzut(5,q)-vxyzut(5,p))*dterm
         rot(1)   = rot(1)+(dv(2)*dp(3)-dv(3)*dp(2))*dterm
         rot(2)   = rot(2)+(dv(3)*dp(1)-dv(1)*dp(3))*dterm
         rot(3)   = rot(3)+(dv(1)*dp(2)-dv(2)*dp(1))*dterm
      ENDDO neiloop
!
!--Limit density to keep it inside the tables
!
      IF (rhofunc > rhomax) rhofunc = rhomax
      IF (rhofunc < rhomin) rhofunc = rhomin
!
      divvfunc = divvfunc/rhofunc
      rotn     = SQRT(rot(1)**2+rot(2)**2+rot(3)**2)/rhofunc
!
      END SUBROUTINE new_rho
