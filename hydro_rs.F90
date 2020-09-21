subroutine hydro_rs
!===============================================================
!  This subroutine calculates the hydrodynamical accelerations 
!  within the ITS scheme
!
!  Last revision: 21/March/2017
!===============================================================
!
!--Load modules
!
  use mod_parameters, only : ndim, maxstep
  use mod_commons,    only : xyzhm, vxyzut, axyzut, ka1, css, cvs,  &
                             rho, press, fh, dPdT, cur, vsigmax,    &
                             div, nb, partype, rank, factor,        &
                             nbody, balsara, ierr, npart, istep0,   &
                             istep, step, dtmax, nstep, active,     &
                             pseudoactive, axyzutp, vxyzutp
  use mod_functions,  only : dk_Q5, wk_Q5, wk, dk
  use mod_parameters, only : ka_min, ka_max
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
#ifdef MPI
  real, dimension((ndim+3)*npart)::sentarray,receivearray
#endif 
  real, dimension(ndim)  :: vel0, pos0, posp
  real :: massp, r2, u2q, u2p, hq, hp, hqp, vqprqp,                 &
          vijj, rhop, rhoq, fq, fp, dterm, dtermp, dtermq, cqp,     &
          rhoqp, rij2, cq, cp, pterm, vsig1, vsig2, qterm1, qterm2, &
          qterm3, kq, kp, kqp, sqp, uterm, vaux, presq, presp,      &
          uintq, uintp, vsig, fqp, tempq, tempp, cvq, cvp, cvqp,    &
          tterm, ptermq, ptermp, fhq, fhp, fhqp, vsigu, S, depdtq,  &
          depdtp, dt, wkp
  real, parameter :: length=0.2d0
  integer ::  i, q, qlocal, p, k, m
!
!--Viscosity Switch
!
!$omp parallel default(none) shared(div,ka1,fh,css,xyzhm)    &
!$omp shared(nbody,partype,istep0,istep,step,dtmax,nstep)    &  
!$omp private(p,S,dt)
!$omp do schedule(runtime)
  do p=1,nbody
     if (partype(p) == 2) cycle
     if ((nstep /=0) .and. (istep0(p) + step(p) /= istep)) cycle
     dt = dtmax*float(step(p))/float(maxstep)

     if (div(p)*(ka_max-ka1(p)) < 0.0) then
        S = -div(p)*(ka_max-ka1(p))
     else
        S = 0.0
     endif

     ka1(p) = ka1(p)  + (-(ka1(p)-ka_min)*(2.0*length*css(p))/ &
             xyzhm(4,p) + S)*dt

     if (ka1(p) > ka_max) ka1(p) = ka_max
     if (ka1(p) < ka_min) ka1(p) = ka_min
  enddo
!$omp end do
!$omp end parallel
!
!--Hydro acceleration
!
  axyzut           = 0.0
  vsigmax(1:npart) = 0.0
  !pseudoactive(1:npart) = .false.
!$omp parallel default(none) shared(xyzhm,vxyzut,vxyzutp,axyzut,vsigmax,partype) &
!$omp shared(div,cur,fh,nb,press,dPdT,css,cvs,ka1,balsara,rho,rank)              &
!$omp shared(nbody,factor,istep,istep0,step,nstep,active,pseudoactive,dtmax)     &
!$omp private(p,q,pos0,vel0,vaux,rhoq,rhop,hq,hp,presq,presp,uintq)              &
!$omp private(uintp,tempq,depdtq,cq,cp,cvq,cvp,fq,fp,kq,kp,fhq)                  &
!$omp private(fhp,massp,tempp,posp,depdtp,u2q,u2p,dtermq,dtermp,dterm)           &
!$omp private(ptermq,ptermp,tterm,uterm,hqp,rhoqp,cqp,cvqp,kqp,vqprqp)           &
!$omp private(r2,vijj,fqp,vsig,qterm1,sqp,vsigu,qterm2,qlocal,m,dt)
!$omp do schedule(runtime)
  partloop : do qlocal = 1, factor
     q = qlocal + rank*factor
!
     if ((q > nbody) .or. (partype(q) == 2)) cycle
     if ((nstep /= 0) .and. (istep0(q) + step(q) /= istep)) cycle
!
!--Variables start up
!
     do k=1,ndim
        pos0(k) = xyzhm(k,q)
        vel0(k) = vxyzut(k,q)
     enddo
     vaux = 0.0
!
!--Store for later use q-particle related variables 
!
     rhoq   = rho(q)
     hq     = xyzhm(4,q)
     presq  = press(q)
     uintq  = vxyzut(4,q)
     tempq  = vxyzut(5,q)
     depdtq = dPdT(q)
     cq     = css(q)
     cvq    = cvs(q)
     if (div(q) /= 0.0) then 
        fq   = abs(div(q))/(abs(div(q))+cur(q)+1.0e-4*cq/hq)
     else
        fq   = 0.0
     endif
     kq     = ka1(q)
     fhq    = fh(q)
!
!--Loop over neighbours
!
     neiloop : do i = 1,nb(1,qlocal)
        p = nb(i+1,qlocal)

        if (partype(p) == 2) cycle
!
!--Save p-particle variables for later use
!
        vqprqp = 0.0
        do k=1,ndim
           posp(k) = pos0(k) - xyzhm(k,p)
           vqprqp  = vqprqp + (vel0(k) - vxyzut(k,p))*posp(k)
        enddo
        r2=posp(1)**2+posp(2)**2+posp(3)**2

        rhop   = rho(p)
        hp     = xyzhm(4,p)
        massp  = xyzhm(5,p)
        presp  = press(p)
        uintp  = vxyzut(4,p)
        tempp  = vxyzut(5,p)
        depdtp = dPdT(p)
        cp     = css(p)
        cvp    = cvs(p)
        if (div(p) /= 0.0) then
           fp   = abs(div(p))/(abs(div(p))+cur(p)+1.0e-4*cp/hp)
        else
           fp   = 0.0
        endif
        kp     = ka1(p)
        fhp    = fh(p)
!
!--Calculate kernel related variables
!
        u2q    = r2/(hq*hq)
        u2p    = r2/(hp*hp)

! Cubic Kernel
        dtermq = fhq*dk(u2q)/hq**(ndim+2)
        dtermp = fhq*dk(u2p)/hp**(ndim+2)

! Quintic Kernel
        !dtermq = dk_Q5(u2q)/(hq**(ndim+2))*fhq
        !dtermp = dk_Q5(u2p)/(hp**(ndim+2))*fhp
        dterm  = 0.5*(dtermp+dtermq)
!
!--Pressure and thermal energy 
!
        ptermq = presq/rhoq**2
        ptermp = presp/rhop**2
        tterm  = tempq*depdtq/rhoq**2
        uterm  = ptermq
!
!--Symmetrizations
!
        hqp    = 0.5*(hq + hp)
        rhoqp  = 0.5*(rhoq + rhop)
        cqp    = 0.5*(cq + cp)
        cvqp   = 0.5*(cvq + cvp)
        kqp    = 0.5*(kq + kp)
        fqp    = 0.5*(fq + fp)
!
!--Calculate viscosity parameters
!
        if (vqprqp < 0.0) then
           vijj   =  vqprqp/(sqrt(r2) + 0.01*hq*hq)
           vsig   =  cqp - vijj
           qterm1 = -fqp*kqp*vsig*vijj/rhoqp
        else
           vsig   =  cqp
           qterm1 =  0.0
        endif
!
!--Thermal conductivity
!
        sqp    =  0.0
        vsigu  =  sqrt(abs(presq-presp)/rhoqp)
        qterm2 = -fqp*sqp*vsigu*(uintq-uintp)/rhoqp*sqrt(r2)
!
!--Velocity signal
!
        vsigmax(qlocal) = max(vsig,vaux)
        vaux = vsigmax(qlocal)      
!
!--Viscous momentum variation
!
        do k=1,ndim
           axyzut(k,qlocal) = axyzut(k,qlocal) - massp*(ptermq*dtermq + ptermp*dtermp + qterm1*dterm)*posp(k)
        enddo
!
!--Viscous thermal energy/temperature variation
!
        axyzut(4,qlocal) = axyzut(4,qlocal) + massp*(uterm*dtermq + 0.5*qterm1*dterm)*vqprqp + massp*qterm2*dterm
        axyzut(5,qlocal) = axyzut(5,qlocal) + massp*(tterm*dtermq + 0.5*qterm1*dterm)*vqprqp/cvq
    enddo neiloop
  enddo partloop
!$omp end do
!$omp end parallel
!
!--Force synchronization
!
#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WorLD,ierr) 
#endif
!
!--Tranfer of accelerations across MPI processes, if necessary
! 
#ifdef MPI
!$omp parallel default(none) shared(sentarray,axyzut,vsigmax,factor)    &
!$omp private(p,k)
!$omp do schedule(runtime)
  do p = 1, factor
     do k = 1, ndim+2
        sentarray((p-1)*(ndim+3)+k) = axyzut(k,p)
     enddo
     sentarray(p*(ndim+3)) = vsigmax(p) 
  enddo
!$omp end do
!$omp end parallel

  call MPI_ALLGATHER(sentarray,(ndim+3)*factor,MPI_doUBLE_PRECISION,&
                     receivearray,(ndim+3)*factor,                  &
                     MPI_doUBLE_PRECISION,MPI_COMM_WorLD,ierr)

!$omp parallel default(none) shared(receivearray,axyzut,vsigmax,nbody)  &
!$omp shared(istep,istep0,step,nstep) private(p,k)
!$omp do schedule(runtime)
  do p = 1, nbody
     if ((nstep == 0) .or. (istep0(p) + step(p) == istep)) then
        do k = 1, ndim+2
           axyzut(k,p) = receivearray((p-1)*(ndim+3)+k)
        enddo
        vsigmax(p) = receivearray(p*(ndim+3))
     endif
  enddo
!$omp end do
!$omp end parallel
#endif MPI

end subroutine hydro_rs
