subroutine burn
!=========================================================================
! This subroutine calculates nuclear burning
!
! Last revision: 6/April/2019
!=========================================================================
!
!--Load modules
!
  use mod_EOS
  use mod_parameters, only : unm, uden, unt, uen, nel, ndim, MASTER,  &
                             maxstep
  use mod_commons,    only : rho, aion, zion, vxyzut, xss, cvs, enuc, &
                             luminuc, dtnuc, partype, nstep, tnow,    &
                             rank, nbody, istep0, istep, step, dtmax
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--Helmholtz EOS definitions
!
  !INCLUDE 'vector_eos.dek'
  !real, dimension(nrowmax) :: ewant_row, temp_row, den_row, abar_row, &
  !zbar_row, det_row, ptot_row, cs_row, cv_row, dpt_row, cp_row, etot_row
  !integer :: jlo_eos, jhi_eos
  !LOGICAL :: eosfail
!
!--Local variables
!
  real, dimension(nel-1) :: xss2                 
  real :: rhop, temp, cvp, enucp, luminucp, sumAE2, t1, t2, abar,   &
          zbar, sumdt, dt
  real, PARAMETER :: rhotiny=5.0
  integer :: i, p, JK, k, m, iread, iteration, itermax, nstep0
  integer, PARAMETER :: iter_max=1000, NIS=nel-2, NSP=NIS+1, NRE=29,&
                            ngrid=60
!
!--Nuclear network variables
!
  real, dimension(ngrid,NRE,11) :: vgrid 
  real, dimension(ngrid,11) :: aNegrid
  real, dimension(ngrid) :: tgrid
  real, dimension(NRE)   :: V2, flag, ndata, n1, n2, n3, n4, a1,    &
                            a2, a3, a4, Qrad, Qnu, QVAL, AE2
  real, dimension(NSP,1) :: XXTOT
  real, dimension(NSP)   :: AN, ZN, BE, YYB, XXXT
  real :: DTMOLDP, SUMYY, DTMNEWP, BE0
  integer, dimension(NRE)   :: K1, K2, K3, K4, K5, K6, K7, K8
  integer :: ireac, kgrid, NSNUC, ICO
  character(5),  dimension(NSP) :: ONC
  character(6),  dimension(NRE) :: z1, z2, z3, z4
  character(37), dimension(NRE) :: reaction
!
  common /cvit/V2,tgrid,aNegrid
  common /cvgrid/vgrid,flag,ndata
  common /network/ireac,kgrid,n1,n2,n3,n4,z1,z2,z3,z4,a1,a2,a3,a4,  &
                  reaction
  common /Q/Qrad,Qnu
  common /CM/K1,K2,K3,K4,K5,K6,K7,K8
  common /CNAME/ONC
  common /CNETW/AN,ZN,QVAL,BE,BE0
  common /ABUND/XXTOT
  common /CSTEP/NSNUC,ICO
  common /CSNC2/YYB
  common /NUCLEOS/XXXT
  common /CODER/AE2
!
!--Initializations
!
  if (NSP /= nel-1) then
    print*,'Cambiar parametros burn!!'
    return
  endif
!
  iread   = 0
  ICO     = nstep
  itermax = 0 
!
!!$OMP PARALLEL DEFAULT(none) shared(vxyzut,cvs,rho,xss,iread)         &       
!!$OMP shared(dtmp_min,enuc,luminuc,aion,zion,tnow) private(m,rhop,temp,cvp) &
!!$OMP private(enucp,luminucp,iteration,sumdt,DTMNEWP,DTMOLDP,sumAE2)        &
!!$OMP private(abar,zbar,temp_row,den_row,cv_row,abar_row,zbar_row,jlo_eos)  &
!!$OMP private(jhi_eos,xss2,XXTOT,XXXT,SUMYY,i,p) reduction(MAX:itermax)
!!$OMP do SCHEDULE(runtime)
  partloop : do p = 1,nbody
     if (partype(p) == 2) cycle
     if (istep0(p) + istep /= istep) cycle
!
!--Density, temperature and composition initialisation (to cgs units)
!  Skip hydrogen: xss(1) == He ... xss(nel-1) == photons
!
     rhop = rho(p)*uden
     temp = vxyzut(5,p)
     cvp  = cvs(p)
     do i = 1,nel-1
       xss2(i) = xss(i+1,p)
     enddo 
!
!--Iteration over the nuclear time-steps
!
     iteration = 0
     sumdt     = 0.0
     enucp     = 0.0
     luminucp  = 0.0
     dt = dtmax*float(step(p))/float(maxstep) 
     burnloop : do
       if ((sumdt >= dt*unt) .or. (iteration > iter_max))  exit burnloop
!
!--Activate nuclear burning only for T>1.0d7
!
       if (temp < 1.0e7) exit burnloop
!
!--Low density particle burning
!
       if (rhop < rhotiny) then 
         print*,'Careful, incinerated particle !!!',p,temp
         exit burnloop
       endif
!
!--Store chemical abundances composition
!
       do k = 1,NSP
         XXTOT(k,1) = xss2(k)
       enddo
!
!--Avoid excessively small time steps
!
       if (sumdt == 0.0) then
         DTMNEWP = dt*unt
       else
         DTMNEWP = min(DTMNEWP,dt*unt-sumdt)
       endif
!
!--Abundances calculation
!
       DTMOLDP = 0.0
       SUMYY   = 0.0
       call SNUC(tnow,DTMNEWP,rhop,temp,SUMYY,DTMOLDP,1,iread)
       if (iread == 0) iread = 1
!
!--If nuclear is bigger than SPH time-step, adopt SPH time-step 
!  (shouldnt happen since dtmold <= dtmnew)
!
       DTMOLDP = min(DTMOLDP,dt*unt-sumdt)
       sumdt   = sumdt + DTMOLDP
!
!--Abundances change
!
       do k = 1,NSP
         xss2(k) = XXXT(k)
       enddo
!
!--Total nuclear luminosity released. AE2 is in erg g-1 s-1, so 
!  change to code units
!
       sumAE2 = sum(AE2)*DTMOLDP*(unm/uen)  
       enucp  = enucp + sumAE2            
!
!--Luminosity/cvs*dtmold calculation
!
       luminucp = luminucp + sumAE2/cvp
!
!--Change temperature only if photodesintegration is present
!  Recalculate specific heats if temperature has changed
!
       if (sumAE2 < 0.0) then
         temp = temp + sumAE2/cvp
         temp_row(1) = temp
         den_row(1)  = rhop
         abar = 0.0
         zbar = 0.0
         do k = 1,nel-2
           abar = abar + xss2(k)/aion(k+1)
           zbar = zbar + xss2(k)*zion(k+1)/aion(k+1)
         enddo
         abar_row(1) = 1.0/abar
         zbar_row(1) = zbar/abar
!
         eoslist(1) = 1
         jlo_eos = 1
         jhi_eos = 1
         call helmeos
         cvp = cv_row(1)*(unm/uen)
       endif
!
       iteration = iteration + 1 
     enddo burnloop 
!
!--Statistics
!
     itermax = max(iteration,itermax)
!         
!--Store changes
!
     vxyzut(5,p) = temp
     cvs(p)      = cvp
     enuc(p)     = enucp
     luminuc(p)  = luminucp
     dtnuc(p)    = sumdt/unt
     do k = 1,nel-1
       xss(k+1,p) = xss2(k)
     end do
  enddo partloop
!!$OMP end do 
!!$OMP end PARALLEL

  if (rank == MASTER) print*,'burn: max. num. iteraciones:',itermax

end subroutine burn
