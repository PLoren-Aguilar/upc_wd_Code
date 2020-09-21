subroutine iter_rhoh
!=============================================================================
!      This subroutine calculates new density and smoothing length iteratively 
!      It also calculates the tree structure, neighbour lists and gravitational
!      forces by calling treeconst
!
!      Last revision: 9/April/2019
!=============================================================================
!
!--Load modules
!
  use mod_commons,    only : xyzhm, rho, fh, div, divt, cur, dhdt, globmax, globdone, rho, partype, nb,       &
                             npart, factor, rank, nbody, ierr, size, step, istep0, istep, nstep, ttot, tirho, &
                             eps, eps3
  use mod_functions,  only : wk_Q5, dk_Q5, wk, dk
  use mod_parameters, only : rhomin, rhomax, hmin, hmax, hfact, nel, ndim, MASTER
!
!--Force to declare everything 
!
  implicit none 
#ifdef MPI
  include 'mpif.h'
#endif
!
!--Local variables
!
  real, dimension(7*npart) :: sentarray, receivearray
  real, dimension(nbody)   :: haux
  real  :: mp, hp, hnew, rhoh, rhop, drhopdh, dhdrhoh, omega, fp, funch, dfunchdh, diff, max_diff, drhodtp, &
           divp, divtp, curp, t1, t2
  real, parameter :: tolerance=1.0e-3
  integer  :: i, k, p, m, plocal, notdone
  integer, parameter :: itermax=250
  logical  :: converged
!
!--Calculate neighbour lists and gravitational forces
!
  call treeconst
#ifdef debug
  if (rank == MASTER) then
    print*, 'treeconst called',rank
  endif
#endif
!
!--Loop through particles
!
  notdone  = 0
  max_diff = 0
!$omp parallel default(none) shared(xyzhm,rho,fh,dhdt,div,divt,cur)     &
!$omp shared(rank,haux,partype,factor,nbody,nb,istep,istep0,step,nstep) &
!$omp private(m,p,plocal,i,mp,hp,rhop,dhdrhoh,omega,fp,funch,dfunchdh)  &
!$omp private(hnew,divp,divtp,curp,drhopdh,drhodtp,rhoh,diff,converged) & 
!$omp reduction(MAX:max_diff) reduction(+:notdone)
!$omp do schedule(runtime)
  partloop : do plocal = 1, factor
    p = plocal + rank*factor
    if ((p > nbody) .or. (partype(p) == 2)) cycle
    if ((nstep /= 0) .and. (istep0(p) + step(p) /= istep)) cycle
!
!--Old h and mass
!
    hp = xyzhm(4,p)
    mp = xyzhm(5,p)
!
!--Main NR loop
!
    diff = 1e30
    converged = .false.
    iterloop : do i=1,itermax
!
!--Calculate density using h, for t
!
      call new_rho(p,plocal,mp,hp,rhop,drhopdh,drhodtp,divp,divtp,curp)
!
!--rho compatible with the old h
!
      rhoh    = mp*(hfact/hp)**ndim
      dhdrhoh = - hp/(ndim*rhop)
      omega   = 1.0 - dhdrhoh*drhopdh
!
!--If Newton-Raphson does not converge, we keep initial data
!
      if (omega < 1.0e-16) then     
        fp = 1.0
        exit iterloop
      end if
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
      if ((hnew < 0.0) .OR. (fp < 1.0e-6)) then
        hnew = hfact*(mp/rhop)**(1.0/real(ndim))
      end if
!
!--Limit the change in h... just in case
!
      if (hnew > 1.2*hp) hnew = 1.2*hp
      if (hnew < 0.8*hp) hnew = 0.8*hp
      if (hnew > hmax)   hnew = hmax
      if (hnew < hmin)   hnew = hmin
!               
!--Look at convergence
!
      !diff = abs(rhop-rhoh)/rhop
      diff = ABS(hnew-hp)/hp
      if ((diff < tolerance) .and. (omega > 0.)) converged = .true.
!
!--Update h
!
      hp = hnew
! 
!--Abandon the loop if finished
!
      if (converged.eqv..true.) exit iterloop 
    enddo  iterloop 
!
!--Check change in h and update particle values
!
    if (diff > tolerance) then
      notdone = notdone + 1
    else
      max_diff = max(max_diff,diff)
    endif
!
    rho(plocal)  = rhop
    fh(plocal)   = fp
    dhdt(plocal) = dhdrhoh*drhodtp*fp
    haux(plocal) = hp
    div(plocal)  = divp
    divt(plocal) = divtp
    cur(plocal)  = curp
  enddo partloop
!$omp end do
!$omp end parallel
!
!--Force synchronization
!
#ifdef MPI     
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
!--Transfer variables across MPI processess if necessary
!
#ifdef MPI
!$omp parallel default(none) shared(sentarray,rho,fh,dhdt,haux,div,divt)&
!$omp shared(cur,factor) private(p)
!$omp do schedule(runtime)
  do p = 1, factor
    sentarray((p-1)*7 + 1) = rho(p)
    sentarray((p-1)*7 + 2) = fh(p)
    sentarray((p-1)*7 + 3) = dhdt(p)
    sentarray((p-1)*7 + 4) = haux(p)
    sentarray((p-1)*7 + 5) = div(p)
    sentarray((p-1)*7 + 6) = divt(p)
    sentarray((p-1)*7 + 7) = cur(p)
  enddo
!$omp end do
!$omp end parallel
!     
  call MPI_ALLGATHER(sentarray,7*factor,MPI_doUBLE_PRECISION,receivearray,7*factor, &
                     MPI_doUBLE_PRECISION,MPI_COMM_WORLD,ierr)
!
!$omp parallel default(none) shared(receivearray,rho,fh,dhdt,haux,div)  &
!$omp shared(divt,cur,xyzhm,nbody,istep,istep0,step,nstep) private(p)
!$omp do schedule(runtime)
  do p = 1, nbody
    if ((nstep == 0) .or. (istep0(p) + step(p) == istep)) then
      rho(p)     = receivearray((p-1)*7 + 1)
      fh(p)      = receivearray((p-1)*7 + 2)
      dhdt(p)    = receivearray((p-1)*7 + 3)
      xyzhm(4,p) = receivearray((p-1)*7 + 4)
      div(p)     = receivearray((p-1)*7 + 5)
      divt(p)    = receivearray((p-1)*7 + 6)
      cur(p)     = receivearray((p-1)*7 + 7)
    endif
  enddo
!$omp end do
!$omp end parallel
#else
!$omp parallel default(none) shared(nbody,xyzhm,haux,partype,istep,istep0,step,nstep) &
!$omp private(p)
!$omp do schedule(runtime)
  do p = 1, nbody
    if ((nstep == 0) .or. (istep0(p) + step(p) == istep)) then
      if (partype(p) /= 2) xyzhm(4,p) = haux(p)
    endif
  enddo
!$omp end do
!$omp end parallel
#endif MPI
!
!--Print status
!
#ifdef MPI
  call MPI_REDUCE(max_diff, globmax, 1, MPI_doUBLE_PRECISION, MPI_MAX, MASTER, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(notdone, globdone, 1, MPI_integer, MPI_SUM, MASTER, MPI_COMM_WORLD, ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! Force synchronization
#else
  globmax  = max_diff
  globdone = notdone
#endif
!
end subroutine iter_rhoh
!
subroutine new_rho(q,qlocal,mq,hq,rhofunc,gradh,drhodt,divvfunc, divtfunc,rotn)
!=====================================================================
!
!  This subroutine calculates the new density value for
!  particle p
!
!=====================================================================
!
!--Load modules
!
  use mod_commons,    only : xyzhm, vxyzut, nb, partype, istep0, istep, step
  use mod_functions,  only : wk_Q5, dk_Q5, wk, dk
  use mod_parameters, only : rhomin, rhomax, ndim
!
!--Force do declare EVERYTHING
!
  implicit none
!
!--I/O variables
!
  real, intent(in)     :: mq, hq
  real, intent(out)    :: rhofunc, gradh, divvfunc, divtfunc, rotn, drhodt
  integer, intent(in)  :: q, qlocal
!
!--Local variables
!
  real, dimension(ndim) :: vel0, pos0, dv, dp, rot
  real :: massp, dterm, r2, u2q, u2p, hp, hq3, hp3, hq5, hq2, hp2, hp5, vijrij, &
          vijj, dkhp, dkq, dkp, wkq, wkp, dkhq
  integer :: p, k, m, i
!
!--Variables start up and definitions
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
!--Loop over neighbouring particles
!
  neiloop : do i = 1,nb(1,qlocal)
    p    = nb(i+1,qlocal)
    !if (partype(q) /= partype(p)) cycle neiloop

! Kernel-related definitions
    hp   = xyzhm(4,p)
    hp2  = hp*hp
    hp3  = hp2*hp
    hp5  = hp3*hp2
    massp = xyzhm(5,p)

! Distance between particles
    vijrij = 0.0
    dimloop : do k = 1,ndim
      dv(k)  = vel0(k) - vxyzut(k,p)
      dp(k)  = pos0(k) - xyzhm(k,p)
      vijrij = vijrij  + dv(k)*dp(k)
    enddo dimloop
    r2 = dp(1)*dp(1)+dp(2)*dp(2)+dp(3)*dp(3)

! Kernel calculation
    u2q   = r2/hq2
    u2p   = r2/hp2

! Cubic Kernel
    wkq   = wk(u2q)/hq3            
    wkp   = wk(u2p)/hp3            
    dkq   = dk(u2q)/hq5
    dkp   = dk(u2p)/hp5
    dkhq  = -(3.*wk(u2q) + sqrt(u2q)*dk(u2q))/(hq2*hq2)
    dkhp  = -(3.*wk(u2p) + sqrt(u2p)*dk(u2p))/(hp2*hp2)

! Quintic Kernel
    !wkq   = wk_Q5(u2q)/hq3
    !wkp   = wk_Q5(u2p)/hp3
    !dkq   = dk_Q5(u2q)/hq5
    !dkp   = dk_Q5(u2p)/hp5
    !dkhq  = -(3.*wk_Q5(u2q) + sqrt(u2q)*dk_Q5(u2q))/(hq2*hq2)
    !dkhp  = -(3.*wk_Q5(u2p) + sqrt(u2p)*dk_Q5(u2p))/(hp2*hp2)
    dterm = massp*0.5*(dkq+dkp)

! Summation over neighbours to calculate density and h-gradient
    rhofunc  = rhofunc  + massp*wkq
    gradh    = gradh    + massp*dkhq

! Rest of the summations to calculate change in density, velocity divergence
! and rotational
    drhodt   = drhodt + massp*vijrij*dkq
    divvfunc = divvfunc + vijrij*dterm
    divtfunc = divtfunc + (vxyzut(5,q)-vxyzut(5,p))*dterm
    rot(1)   = rot(1)+(dv(2)*dp(3)-dv(3)*dp(2))*dterm
    rot(2)   = rot(2)+(dv(3)*dp(1)-dv(1)*dp(3))*dterm
    rot(3)   = rot(3)+(dv(1)*dp(2)-dv(2)*dp(1))*dterm
  enddo neiloop
!
!--Limit density to keep it inside the tables
!
  if (rhofunc > rhomax) rhofunc = rhomax
  if (rhofunc < rhomin) rhofunc = rhomin
!
  divvfunc = divvfunc/rhofunc
  rotn     = sqrt(rot(1)**2+rot(2)**2+rot(3)**2)/rhofunc
!
end subroutine new_rho
