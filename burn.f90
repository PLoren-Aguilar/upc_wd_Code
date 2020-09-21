      SUBROUTINE burn
!=========================================================================
! This subroutine calculates nuclear burning
!
! Last revision: 15/March/2015
!=========================================================================
!
!--Load modules
!
      USE mod_parameters, ONLY : unm, uden, unt, uen, nel, ndim, MASTER
      USE mod_commons, ONLY : rho, aion, zion, vxyzut, xss, cvs, enuc,  &
      luminuc, dtnuc, partype, nstep, dtmp_min, tnow, rank, nbody
      USE mod_EOS
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Helmholtz EOS definitions
!
      !INCLUDE 'vector_eos.dek'

!
      !REAL, DIMENSION(nrowmax) :: ewant_row, temp_row, den_row, abar_row, &
      !zbar_row, det_row, ptot_row, cs_row, cv_row, dpt_row, cp_row, etot_row
      !INTEGER :: jlo_eos, jhi_eos
      !LOGICAL :: eosfail
!
!--Local variables
!
      REAL, DIMENSION(nel-1) :: xss2                 
      REAL :: rhop, temp, cvp, enucp, luminucp, sumAE2, t1, t2, abar,&
                 zbar, sumdt
      REAL, PARAMETER :: rhotiny=5.0
      INTEGER :: i, p, JK, k, m, iread, iteration, itermax, nstep0
      INTEGER, PARAMETER :: iter_max=1000, NIS=nel-2, NSP=NIS+1, NRE=29,&
                            ngrid=60
!
!--Nuclear network variables
!
      REAL, DIMENSION(ngrid,NRE,11) :: vgrid 
      REAL, DIMENSION(ngrid,11) :: aNegrid
      REAL, DIMENSION(ngrid) :: tgrid
      REAL, DIMENSION(NRE)   :: V2, flag, ndata, n1, n2, n3, n4, a1, &
                                   a2, a3, a4, Qrad, Qnu, QVAL, AE2
      REAL, DIMENSION(NSP,1) :: XXTOT
      REAL, DIMENSION(NSP)   :: AN, ZN, BE, YYB, XXXT
      REAL :: DTMOLDP, SUMYY, DTMNEWP, BE0
      INTEGER, DIMENSION(NRE)   :: K1, K2, K3, K4, K5, K6, K7, K8
      INTEGER :: ireac, kgrid, NSNUC, ICO
      CHARACTER(5),  DIMENSION(NSP) :: ONC
      CHARACTER(6),  DIMENSION(NRE) :: z1, z2, z3, z4
      CHARACTER(37), DIMENSION(NRE) :: reaction
!
      COMMON /cvit/V2,tgrid,aNegrid
      COMMON /cvgrid/vgrid,flag,ndata
      COMMON /network/ireac,kgrid,n1,n2,n3,n4,z1,z2,z3,z4,a1,a2,a3,a4,  &
                      reaction
      COMMON /Q/Qrad,Qnu
      COMMON /CM/K1,K2,K3,K4,K5,K6,K7,K8
      COMMON /CNAME/ONC
      COMMON /CNETW/AN,ZN,QVAL,BE,BE0
      COMMON /ABUND/XXTOT
      COMMON /CSTEP/NSNUC,ICO
      COMMON /CSNC2/YYB
      COMMON /NUCLEOS/XXXT
      COMMON /CODER/AE2
!
!--Initializations
!
      IF (NSP /= nel-1) THEN
        PRINT*,'Cambiar parametros burn!!'
        RETURN
      END IF
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
!!$OMP DO SCHEDULE(runtime)
      partloop : DO p = 1,nbody
         IF (partype(p) == 2) CYCLE
!
!--Density, temperature and composition initialisation (to cgs units)
!  Skip hydrogen: xss(1) == He ... xss(nel-1) == photons
!
         rhop = rho(p)*uden
         temp = vxyzut(5,p)
         cvp  = cvs(p)
         DO i = 1,nel-1
            xss2(i) = xss(i+1,p)
         ENDDO 
!
!--Iteration over the nuclear time-steps
!
         iteration = 0
         sumdt     = 0.0
         enucp     = 0.0
         luminucp  = 0.0
         burnloop : DO
            IF ((sumdt >= dtmp_min*unt) .OR. (iteration > iter_max))  EXIT burnloop
!
!--Activate nuclear burning only for T>1.0d7
!
            IF (temp < 1.0e7) EXIT burnloop
!
!--Low density particle burning
!
            IF (rhop < rhotiny) THEN 
               PRINT*,'Careful, incinerated particle !!!',p,temp
               EXIT burnloop
            END IF
!
!--Store chemical abundances composition
!
            DO k = 1,NSP
               XXTOT(k,1) = xss2(k)
            ENDDO
!
!--Avoid excessively small time steps
!
            IF (sumdt == 0.0) THEN
               DTMNEWP = dtmp_min*unt
            ELSE
               DTMNEWP = MIN(DTMNEWP,dtmp_min*unt-sumdt)
            ENDIF
!
!--Abundances calculation
!
            DTMOLDP = 0.0
            SUMYY   = 0.0
            CALL SNUC(tnow,DTMNEWP,rhop,temp,SUMYY,DTMOLDP,1,iread)
            IF (iread == 0) iread = 1
!
!--If nuclear is bigger than SPH time-step, adopt SPH time-step 
!  (shouldnt happen since dtmold <= dtmnew)
!
            DTMOLDP = MIN(DTMOLDP,dtmp_min*unt-sumdt)
            sumdt   = sumdt + DTMOLDP
!
!--Abundances change
!
            DO k = 1,NSP
               xss2(k) = XXXT(k)
            ENDDO
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
            IF (sumAE2 < 0.0) THEN
               temp = temp + sumAE2/cvp
               temp_row(1) = temp
               den_row(1)  = rhop
               abar = 0.0
               zbar = 0.0
               DO k = 1,nel-2
                  abar = abar + xss2(k)/aion(k+1)
                  zbar = zbar + xss2(k)*zion(k+1)/aion(k+1)
               ENDDO
               abar_row(1) = 1.0/abar
               zbar_row(1) = zbar/abar
!
               eoslist(1) = 1
               jlo_eos = 1
               jhi_eos = 1
               CALL helmeos
               cvp = cv_row(1)*(unm/uen)
            ENDIF
!
            iteration = iteration + 1 
         ENDDO burnloop 
!
!--Statistics
!
         itermax = MAX(iteration,itermax)
!         
!--Store changes
!
         vxyzut(5,p) = temp
         cvs(p)      = cvp
         enuc(p)     = enucp
         luminuc(p)  = luminucp
         dtnuc(p)    = sumdt/unt
         DO k = 1,nel-1
            xss(k+1,p) = xss2(k)
         END DO
       ENDDO partloop
!!$OMP END DO 
!!$OMP END PARALLEL
!
      IF (rank == MASTER) PRINT*,'burn: max. num. iteraciones:',itermax
!
      END SUBROUTINE burn
