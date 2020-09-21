subroutine layer
!===================================================================
!  This subroutine keep the helium layer stick to the top of the wd
!
!  Last revision: 9/April/2019
!===================================================================
!
!--Load modules
!
  use mod_parameters, only : ndim
  use mod_commons,    only : xyzhm, vxyzut, rho, partype, nb, norm, &
                             rhoG, nbody, rank, factor, npart, ierr,&
                             istep, istep0, step, nstep, RELFLAG
  use mod_functions,  only : wk
!
!--Force to declare EVERYTHING
!
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif 
!
!--Local variables
!
  real, dimension(ndim+2, npart)  :: newvxyzut
  real, dimension((ndim+2)*npart) :: sentarray, receivearray
  real, dimension(ndim) :: vxyz
  real :: massp, r2, u2p, hp, scalar, dr, drx, dry, drz, dvx,       &
          dvy, dvz, rhop, rhoq, chi, wkp, chichi, tCO, uCO, temp,   &
          normq
  real, parameter :: norma = 10./9.
  integer ::  i, q, qlocal, p, k, m
  integer, parameter :: nu=1
!
!--Mean velocity flow
!
  newvxyzut = 0.0
!$omp parallel default(none) shared(xyzhm,vxyzut,rho,rank,partype,nb,norm,factor,RELFLAG) & 
!$omp shared(rhoG,nbody,newvxyzut,nstep,istep,istep0,step) private(i,p,q,rhoq,rhop,hp)  &
!$omp private(scalar,chi,chichi,dr,drx,dry,drz,dvx,dvy,dvz,massp,u2p,wkp,qlocal)  &
!$omp private(vxyz,tCO,uCO,temp,normq)
!$omp do schedule(runtime)
  partloop : do qlocal = 1, factor
     q = qlocal + rank*factor
     if (q > nbody) cycle
!
!--Go only through active He particles
!
     if (partype(q) /= 1) cycle
     !if (istep0(q) + step(q) /= istep) cycle
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
     neiloop : do i = 1,nb(1,qlocal)
        p = nb(i+1,qlocal)
        if (partype(p) /= 0) cycle
!
!--Calculate unit vectors and velocity scalar product
!
        drx = xyzhm(1,q) - xyzhm(1,p)
        dry = xyzhm(2,q) - xyzhm(2,p)
        drz = xyzhm(3,q) - xyzhm(3,p)
        dr  = SQRT(drx*drx + dry*dry + drz*drz)
        if (dr > 0.0) then
           drx = drx/dr
           dry = dry/dr
           drz = drz/dr
        else
           drx = 0.0
           dry = 0.0
           drz = 0.0
        endif
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
        !wkp = wk(u2p)/(hp**ndim)
        wkp = norma*wk(u2p)*u2p/(hp**ndim)
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
        !uCO = uCO + 0.5*nu*massp*chi*scalar*scalar*wkp
        tCO = tCO + (massp/rhop)*temp*wkp
     enddo neiloop
!
!--Define "switch" function
!
     if (rhoG(q) == 0.0) then
       chichi = 0.0
     else
       chichi = EXP(-0.1*rhoq/rhoG(q))
       !chichi = DEXP(-(xyzhm(5,p)/xyzhm(5,q))*rhoq/rhoG(q))
     endif
!
!--Update velocity and thermal energy
!
     do k=1,3
       newvxyzut(k,qlocal) = (1.-chichi)*vxyzut(k,q) + chichi*vxyz(k)/normq
     enddo

     if (RELFLAG .eqv. .true.) then
       newvxyzut(4,qlocal) = vxyzut(4,q)
       newvxyzut(5,qlocal) = vxyzut(5,q)
     else  
       newvxyzut(4,qlocal) = (1.-chichi)*vxyzut(4,q) +  chichi*uCO/normq
       newvxyzut(5,qlocal) = (1.-chichi)*vxyzut(5,q) +  chichi*tCO/normq
     endif
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
!--Tranfer of accelerations across MPI processes, if necessary
! 
#ifdef MPI
!$omp parallel default(none) shared(sentarray,newvxyzut,factor)        &
!$omp private(p,k)
!$omp do schedule(runtime)
  do p = 1, factor
    do k = 1, ndim+2
      sentarray((p-1)*(ndim+2)+k) = newvxyzut(k,p)
    enddo
  enddo
!$omp end do
!$omp end parallel
!
  call MPI_ALLGATHER(sentarray,(ndim+2)*factor,MPI_doUBLE_PRECISION,&
                     receivearray,(ndim+2)*factor,                  &
                     MPI_doUBLE_PRECISION,MPI_COMM_WORLD,ierr)
!
!$omp parallel default(none) shared(nbody,receivearray,vxyzut,partype)  &
!$omp private(p,k)
!$omp do schedule(runtime)
  do p = 1, nbody
    do k = 1, ndim+2
      if (partype(p) == 1) then
        vxyzut(k,p)=receivearray((p-1)*(ndim+2)+k)
      endif
    enddo
  enddo
!$omp end do
!$omp end parallel
#else
!$omp parallel default(none) shared(nbody,newvxyzut,vxyzut,partype)  &
!$omp private(p,k)
!$omp do schedule(runtime)
  do p = 1, nbody
    do k = 1, ndim+2
      if(partype(p) == 1) vxyzut(k,p) = newvxyzut(k,p)
    enddo
  enddo
!omp end do
!$omp end parallel
#endif MPI
!
end subroutine layer
