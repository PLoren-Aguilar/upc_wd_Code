      SUBROUTINE hydro_rs
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
      USE mod_commons, ONLY : xyzhm, vxyzut, axyzut, ka1, css, cvs, rho,&
      press, fh, dPdT, cur, vsigmax, div, nb, partype, rank, factor,    &
      dtmp_min, nbody, balsara, ierr, npart
      USE mod_functions, ONLY  : dk
      USE mod_parameters, ONLY : ka_min, ka_max
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
      REAL, DIMENSION((ndim+3)*npart)::sentarray,receivearray
#endif 
      REAL, DIMENSION(ndim) :: vel0, pos0, posp
      REAL :: massp, r2, u2q, u2p, hq, hp, hqp, vqprqp,                 &
                 vijj, rhop, rhoq, fq, fp, dterm, dtermp, dtermq, cqp,  &
                 rhoqp, rij2, cq, cp, pterm, vsig1, vsig2, qterm1,      &
                 qterm2, qterm3, kq, kp, kqp, sqp, uterm, vaux, presq,  &
                 presp, uintq, uintp, vsig, fqp, tempq, tempp, cvq, cvp,&
                 cvqp, tterm, ptermq, ptermp, fhq, fhp, fhqp, vsigu, S, &
                 depdtq, depdtp
      REAL, PARAMETER :: length=0.2d0
      INTEGER ::  i, q, qlocal, p, k, m
!
!--Viscosity Switch
!
!$OMP PARALLEL DEFAULT(none) shared(div,ka1,fh,css,xyzhm,dtmp_min) &
!$OMP shared(nbody,partype) private(p,S)
!$OMP DO SCHEDULE(runtime)
      DO p=1,nbody
         IF (partype(p) == 2) CYCLE
!
         IF (div(p)*(ka_max-ka1(p)) < 0.0) THEN
            S = -div(p)*(ka_max-ka1(p))
         ELSE
            S = 0.0
         ENDIF
!
         ka1(p) = ka1(p)  + (-(ka1(p)-ka_min)*(2.0*length*css(p))/ &
                  xyzhm(4,p) + S)*dtmp_min
!
         IF (ka1(p) > ka_max) ka1(p) = ka_max
         IF (ka1(p) < ka_min) ka1(p) = ka_min
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
!--Hydro acceleration
!
      axyzut  = 0.0
      vsigmax = 0.0
!$OMP PARALLEL DEFAULT(none) shared(xyzhm,vxyzut,axyzut,vsigmax,partype) &
!$OMP shared(div,cur,fh,nb,press,dPdT,css,cvs,ka1,balsara,rho,rank)      &
!$OMP shared(nbody,factor)                                               &
!$OMP private(p,q,pos0,vel0,vaux,rhoq,rhop,hq,hp,presq,presp,uintq)      &
!$OMP private(uintp,tempq,depdtq,cq,cp,cvq,cvp,fq,fp,kq,kp,fhq)          &
!$OMP private(fhp,massp,tempp,posp,depdtp,u2q,u2p,dtermq,dtermp,dterm)   &
!$OMP private(ptermq,ptermp,tterm,uterm,hqp,rhoqp,cqp,cvqp,kqp,vqprqp)   &
!$OMP private(r2,vijj,fqp,vsig,qterm1,sqp,vsigu,qterm2,qlocal,m)
!$OMP DO SCHEDULE(runtime)
      partloop : DO qlocal = 1, factor
         q = qlocal + rank*factor
!
         IF ((q > nbody) .OR. (partype(q) == 2)) CYCLE
!
!--Variables start up
!
         DO k=1,ndim
            pos0(k) = xyzhm(k,q)
            vel0(k) = vxyzut(k,q)
         ENDDO
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
         fq     = ABS(div(q))/(ABS(div(q))+cur(q)+1.0e-4*cq/hq)
         kq     = ka1(q)
         fhq    = fh(q)
!
!--Loop over neighbours
!
         neiloop : DO i = 1,nb(1,qlocal)
            p = nb(i+1,qlocal)
!
            IF (partype(p) == 2) CYCLE
!
!--Save p-particle variables for later use
!
            vqprqp = 0.0
            DO k=1,ndim
               posp(k) = pos0(k) - xyzhm(k,p)
               vqprqp  = vqprqp + (vel0(k) - vxyzut(k,p))*posp(k)
            ENDDO
            r2=posp(1)**2+posp(2)**2+posp(3)**2
!
            rhop   = rho(p)
            hp     = xyzhm(4,p)
            massp  = xyzhm(5,p)
            presp  = press(p)
            uintp  = vxyzut(4,p)
            tempp  = vxyzut(5,p)
            depdtp = dPdT(p)
            cp     = css(p)
            cvp    = cvs(p)
            fp     = ABS(div(p))/(ABS(div(p))+cur(p)+1.0e-4*cp/hp)
            kp     = ka1(p)
            fhp    = fh(p)
!
!--Calculate kernel related variables
!
            u2q    = r2/(hq*hq)
            u2p    = r2/(hp*hp)
            dtermq = dk(u2q)/(hq**(ndim+2))*fhq
            dtermp = dk(u2p)/(hp**(ndim+2))*fhp
            dterm  = 0.5*(dtermp+dtermq)
!
!--Pressure and thermal energy 
!
            ptermq = presq/rhoq**2
            ptermp = presp/rhop**2
            tterm  = tempq*depdtq/(rhoq**2) 
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
            IF (vqprqp < 0.0) THEN
               !IF (r2 == 0.0d0) THEN
               !   vijj = 0.0d0
               !ELSE
                  vijj = vqprqp/(SQRT(r2) + 0.1*hq)
               !ENDIF 
               vsig    = cqp - vijj
               qterm1  = -fqp*kqp*vsig*vijj/rhoqp
            ELSE
               vsig    = cqp
               qterm1  = 0.0
            ENDIF
!
!--Thermal conductivity
!
            sqp    = 0.0
            vsigu  = SQRT(ABS(presq-presp)/rhoqp)
            qterm2 = -fqp*sqp*vsigu*(uintq-uintp)/rhoqp*SQRT(r2)
!
!--Velocity signal
!
            vsigmax(qlocal) = MAX(vsig,vaux)
            vaux = vsigmax(qlocal)      
!
!--Viscous momentum variation
!
            DO k=1,ndim
                  axyzut(k,qlocal) = axyzut(k,qlocal) - massp*(ptermq*  &
      dtermq + ptermp*dtermp + qterm1*dterm)*posp(k)
            ENDDO
!
!--Viscous thermal energy/temperature variation
!
            axyzut(4,qlocal) = axyzut(4,qlocal) + massp*(uterm*dtermq + &
                        0.5*qterm1*dtermq)*vqprqp + massp*qterm2*dterm
            axyzut(5,qlocal) = axyzut(5,qlocal) + massp*(tterm*dtermq + &
                        0.5*qterm1*dtermq)*vqprqp/cvq
         END DO neiloop
      END DO partloop
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
!$OMP PARALLEL DEFAULT(none) shared(sentarray,axyzut,vsigmax,factor)    &
!$OMP private(p,k)
!$OMP DO SCHEDULE(runtime)
      DO p = 1, factor
         DO k = 1, ndim+2
            sentarray((p-1)*(ndim+3)+k) = axyzut(k,p)
         ENDDO
         sentarray(p*(ndim+3)) = vsigmax(p) 
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
      CALL MPI_ALLGATHER(sentarray,(ndim+3)*factor,MPI_DOUBLE_PRECISION,&
                         receivearray,(ndim+3)*factor,                  &
                         MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
!
!$OMP PARALLEL DEFAULT(none) shared(receivearray,axyzut,vsigmax,nbody)  &
!$OMP private(p,k)
!$OMP DO SCHEDULE(runtime)
      DO p = 1, nbody
         DO k = 1, ndim+2
            axyzut(k,p) = receivearray((p-1)*(ndim+3)+k)
         ENDDO
         vsigmax(p) = receivearray(p*(ndim+3))
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
#endif MPI
!
      END SUBROUTINE hydro_rs
