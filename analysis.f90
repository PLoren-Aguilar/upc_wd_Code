      SUBROUTINE analysis
!========================================================================
!
!  This subroutine is a call to various analysis subroutines
!
!  Last revision: 14/March/2015
!========================================================================
! 
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Startout the code reading the necessary data files and parameters
!
      CALL startout
!
!--Analize the results
!
      CALL analyze
!
!--Do some chemical element data files in order to render
!
      CALL outchem
!
      END SUBROUTINE analysis
!
      SUBROUTINE analyze
!
!--Load modules
!
      USE mod_commons
      USE mod_parameters
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local definitions
!
      REAL, DIMENSION(nel) :: xssp
      REAL :: masa ,mass, dens, temp, rin, rout, rmax, dr, r,        &
                 xcm, ycm ,zcm, angx, angy, angz, Imom, RWD, omega, v,  &
                 rest, cvsp, cpsp, dPdTp, gg, presion, gama, vxcm, vycm,&
                 vzcm
      REAL, PARAMETER :: alpha = 1.0, LL = 0.1*unl, GNew = 6.67e-8
      INTEGER :: p, k, n, l, num, inum
!      TYPE :: outdata
!        REAL, DIMENSION(ndim+2,npart) :: ioxyzhm
!        REAL, DIMENSION(ndim+2,npart) :: iovxyzut
!        REAL, DIMENSION(npart) :: iorho
!        REAL, DIMENSION(npart) :: ioka1
!        INTEGER, DIMENSION(npart) :: iopartype
!        INTEGER, DIMENSION(npart) :: iostar
!      END TYPE
!      TYPE(outdata) :: data_array
!
      TYPE :: outdata
        REAL, DIMENSION(ndim+2) :: ioxyzhm
        REAL, DIMENSION(ndim+2) :: iovxyzut
        REAL :: iorho
        REAL :: ioka1
        INTEGER :: iopartype
        INTEGER :: iostar
      END TYPE
      TYPE(outdata), DIMENSION(npart) :: data_array
!
      TYPE :: coutdata
        REAL, DIMENSION(nel) :: comp
      END TYPE
      TYPE(coutdata), DIMENSION(npart) :: cdata_array
!
!--Do a test to find the most efficient way to write bodi binary files
!
      DO p = 1,npart
         data_array(p)%ioxyzhm(1:ndim+2)  = xyzhm(1:ndim+2,p)
         data_array(p)%iovxyzut(1:ndim+2) = vxyzut(1:ndim+2,p)
         data_array(p)%iorho = rho(p)
         data_array(p)%ioka1 = ka1(p)
         data_array(p)%iopartype = partype(p)
         data_array(p)%iopartype = star(p)
         cdata_array(p)%comp(1:nel) = xss(1:nel,p)
      ENDDO
! 
      OPEN (1, FILE='test_bin', FORM="unformatted")
      WRITE(1) tnow,npart,data_array
      CLOSE(1)
!
      DO p=1, npart
         data_array(p)%ioxyzhm(1:ndim+2) = 0.0
      ENDDO
!
      OPEN (1, FILE='test_bin', FORM="unformatted")
      READ(1) tnow,num,data_array
      CLOSE(1)
      PRINT*, tnow, num, data_array(100)%ioxyzhm(1:ndim+2)
!
!--Now let's go for comp binary files
!
      OPEN (2, FILE='ctest_bin', FORM="unformatted")
      WRITE(2) tnow, cdata_array
      CLOSE(2)
      STOP      
!
!--Build a bodi file with showing the density of a given
!  chemical element
!
!      DO p = 1, npart
!         rest = 0.0
!         rest = rest + xss(2,p)
!         DO k=5,nel
!            rest = rest + xss(k,p)
!         ENDDO
!         !IF (xss(2,p).GT.1e-10) THEN
!         !   PRINT*, p,rest,(xss(k,p),k=1,nel)
!         !   STOP
!         !ENDIF
!         IF (partype(p).EQ.0) THEN
!            vxyzut(4,p) = rest
!         ELSE
!!            vxyzut(4,p) = 0.0
!         ENDIF
!      ENDDO
!
!      CALL noninertial
!
!--Write results
!
!      nstep = nstep + 1
!      CALL outdata
!      RETURN
!
!--Open I/O files
!
      OPEN (unit=1, file='profiles.dat', status='new')
      OPEN (unit=2, file='totals.dat', status='new')
!
!--This sould be the headers for the profiles.dat file
!
!# rad(0.1Rsun)  Temp(K)  Dens(g/cm^3)   Omega     Mass(M_sun) cp(erg/gK) dPdT(dyn/K) qc(erg/cm^2s) 1-1/gamma  Press (dyn)
!#===================================================================================================================================================#
!# rad(0.1Rsun)  Temp(K)  Dens(g/cm^3) Mass(M_sun) xcm(0.1Rsun)  ycm()      zcm()     vxcm (Km/s)    vycm ()     vzcm()  v()  He(fraction)  C()  O()
!#===================================================================================================================================================#
!
!--Position all data around the centre of the primary star
!
      mass  = 0.0
      xcm   = 0.0
      ycm   = 0.0
      zcm   = 0.0
      vxcm  = 0.0
      vycm  = 0.0
      vzcm  = 0.0
      DO p = 1, nbody
         !IF (star(p) == 2) THEN 
         IF (rho(p).GT.1.0) THEN 
            mass = mass + xyzhm(5,p)
            xcm  = xcm  + xyzhm(5,p)*xyzhm(1,p)
            ycm  = ycm  + xyzhm(5,p)*xyzhm(2,p)
            zcm  = zcm  + xyzhm(5,p)*xyzhm(3,p)
!
            vxcm  = vxcm  + xyzhm(5,p)*vxyzut(1,p)
            vycm  = vycm  + xyzhm(5,p)*vxyzut(2,p)
            vzcm  = vzcm  + xyzhm(5,p)*vxyzut(3,p)
         ENDIF
      ENDDO
      xcm = xcm/mass
      ycm = ycm/mass
      zcm = zcm/mass
!
      vxcm = vxcm/mass
      vycm = vycm/mass
      vzcm = vzcm/mass

      DO p = 1, nbody
         xyzhm(1,p) = xyzhm(1,p) - xcm
         xyzhm(2,p) = xyzhm(2,p) - ycm
         xyzhm(3,p) = xyzhm(3,p) - zcm
!
         vxyzut(1,p) = vxyzut(1,p) - vxcm
         vxyzut(2,p) = vxyzut(2,p) - vycm
         vxyzut(3,p) = vxyzut(3,p) - vzcm
      ENDDO

      DO p = 1, nbody
         r = SQRT(xyzhm(1,p)**2  + xyzhm(2,p)**2 + xyzhm(3,p)**2)
         v = SQRT(vxyzut(1,p)**2 + vxyzut(2,p)**2 + vxyzut(3,p)**2)
         !PRINT*, r,v*(1e-5*unl/unt)
      ENDDO
!
!--Calculate radial profiles of density, temperature... etc in a
!  cross-section around the mid-plane
!
      n    = 1000
      rmax = log10(4.0)
      rin  = log10(1e-6)
      dr   = (rmax-rin)/n
      rout = rin + dr
      masa = 0.0
     
      DO k = 1,n
         num   = 0.0 
         mass  = 0.0
         dens  = 0.0
         temp  = 0.0
         omega = 0.0
         xssp  = 0.0
         DO p=1, nbody
            r = SQRT(xyzhm(1,p)**2  + xyzhm(2,p)**2 + xyzhm(3,p)**2)
            v = SQRT(vxyzut(1,p)**2 + vxyzut(2,p)**2 + vxyzut(3,p)**2)
            IF ((log10(r) > rin) .AND. (log10(r) <= rout)) THEN 
               masa = masa + xyzhm(5,p)
               num= num + 1
               IF (ABS(xyzhm(3,p)).LT.0.005) THEN
                  mass  = mass  + xyzhm(5,p)
                  temp  = temp  + xyzhm(5,p)*vxyzut(5,p)
                  dens  = dens  + xyzhm(5,p)*rho(p)
                  omega = omega + xyzhm(5,p)*v
                  DO l=1,nel
                     xssp(l) = xssp(l) + xyzhm(5,p)*xss(l,p)
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
!
         IF (mass > 0.0) THEN
            temp  = temp/mass
            dens  = dens/mass
            omega = omega/mass
            DO l=1,nel
               xssp(l) = xssp(l)/mass
            ENDDO
!
!--Call EOS to calculate the value of the specific heat
!
            CALL deg(dens,temp,presion,xssp,cvsp,cpsp,dPdTp,gama)
!
!--Write results
!
            gg = (GNew*mass*unm/(unl*10**((rin+rout)/2.))**2)**0.5
            !WRITE(1,'(11(1pe12.4))') 10**((rin+rout)/2.),temp,dens*uden,&
            !                        (omega/unt)/(10**((rin+rout)/2.)),  &
            !                        masa,cvsp,dPdTp,                    &
            !                        alpha**2*(dens*uden)*cvsp*temp*     &
            !                        (gg*LL)**0.5*(2**1.5),1 - cvsp/cpsp,&
            !                        presion,gama
            !IF (dens > 1.0) THEN
            WRITE(1,'(16(1pe12.4))') 10**((rin+rout)/2.),temp,dens*uden,     &
                                    masa,xcm,ycm,zcm,vxcm*(1.e-5*unl/unt),   &
                                    vycm*(1e-5*unl/unt),vzcm*(1e-5*unl/unt), &
                                    omega*(1e-5*unl/unt),xssp(4)/xssp(2),    &
                                    xssp(4)/xssp(3),1.,xssp(4)/xssp(5),      &
                                    xssp(4)/xssp(6)
            !ENDIF
         ENDIF
         rin  = rout
         rout = rout + dr
      ENDDO
!
!--Calculate TOTAL values for the angular momentum and moment of inertia
!  for the WD and disk
!
      RWD  = 0.1
      angx = 0.0
      angy = 0.0
      angz = 0.0
      Imom = 0.0
      mass = 0.0
      DO p = 1, nbody
         r = SQRT(xyzhm(1,p)**2 + xyzhm(2,p)**2)
         IF (r <= RWD) THEN
             mass = mass + xyzhm(5,p)
             angx = angx + xyzhm(5,p)*(xyzhm(2,p)*vxyzut(3,p) -         &
                    xyzhm(3,p)*vxyzut(2,p))
             angy = angy + xyzhm(5,p)*(xyzhm(1,p)*vxyzut(3,p) -         &
                    xyzhm(3,p)*vxyzut(1,p))
             angz = angz + xyzhm(5,p)*(xyzhm(1,p)*vxyzut(2,p) -         &
                    xyzhm(2,p)*vxyzut(1,p))
             Imom = Imom + xyzhm(5,p)*r*r
         ENDIF
      ENDDO
      WRITE(2,*) 'Values for the WD'
      WRITE(2,'(5(1pe12.4))') angx*unm*unl**2/unt,angy*unm*unl**2/unt,  &
                              angz*unm*unl**2/unt,Imom*unm*unl**2, mass
!
      angx = 0.0
      angy = 0.0
      angz = 0.0
      Imom = 0.0
      mass = 0.0
      DO p = 1, nbody
         r = SQRT(xyzhm(1,p)**2 + xyzhm(2,p)**2)
         IF (r > RWD) THEN
             mass = mass + xyzhm(5,p)
             angx = angx + xyzhm(5,p)*(xyzhm(2,p)*vxyzut(3,p) -         &
                    xyzhm(3,p)*vxyzut(2,p))
             angy = angy + xyzhm(5,p)*(xyzhm(1,p)*vxyzut(3,p) -         &
                    xyzhm(3,p)*vxyzut(1,p))
             angz = angz + xyzhm(5,p)*(xyzhm(1,p)*vxyzut(2,p) -         &
                    xyzhm(2,p)*vxyzut(1,p))
             Imom = Imom + xyzhm(5,p)*r*r
         ENDIF
      ENDDO
      WRITE(2,*) 'Values for the Disk'
      WRITE(2,'(5(1pe12.4))') angx*unm*unl**2/unt,angy*unm*unl**2/unt,  &
                              angz*unm*unl**2/unt,Imom*unm*unl**2, mass
!
!--Close I/O files
!
      CLOSE(1)
      CLOSE(2)
!
      END SUBROUTINE analyze

      SUBROUTINE deg(rhop,temp,presion,xssp,cvsp,cpsp,dPdTp,gama)
!========================================================================
!     THIS   SUBROUTINE   ALLOWS   TO  SELECT   BETWEEN  THE DIFFERENT  
!     EQUATIONS  OF  STATE  OF  THE  CODE. IT'S IMPORTANT  TO  REALIZE  
!     THAT  EVERY  EOS  SUBROUTINE MUST PROVIDE (AT LEAST):
!
!     * TEMPERATURE, PRESSURE, Cv AND SPEED OF SOUND
!
!     Last revision: 15/March/2015!
!=========================================================================
!
!--Load modules
!
      USE mod_parameters
      USE mod_commons
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
!--I/O variables
!
      REAL, DIMENSION(nel), INTENT(IN) :: xssp
      REAL, INTENT(IN) :: rhop, temp
      REAL, INTENT(OUT) :: gama, presion, cvsp, cpsp, dPdTp
!
!--Local variables
!
      !INTEGER :: jlo_eos, jhi_eos
      !REAL, DIMENSION(nrowmax) :: ewant_row, temp_row, den_row, abar_row, &
      !zbar_row, det_row, ptot_row, cs_row, cv_row, dpt_row, cp_row, etot_row
      !LOGICAL :: eosfail
      REAL    :: abar, zbar
      INTEGER :: p, k
!
!--Set initial data for the EOS
!
      temp_row(1)  = temp
      IF (temp.GT.tmax) temp_row(1) = tmax
      IF (temp.LT.tmin) temp_row(1) = tmin
      den_row(1)   = rhop*uden    ! EOS works in cgs units!!!
!
      abar = 0.0
      zbar = 0.0
      DO k=1,nel-1
         abar = abar + xssp(k)/aion(k)
         zbar = zbar + xssp(k)*zion(k)/aion(k)
      ENDDO
      abar_row(1) = 1.0/abar
      zbar_row(1) = zbar/abar
!
!--Set the number of elements sent into the EOS
!
      jlo_eos=1
      jhi_eos=1
!
!--Call EOS to obtain a first estimate for energy, pressure and derivatives
!
      !CALL helmeos(jlo_eos,jhi_eos,temp_row,den_row,abar_row,zbar_row,  &
      !      det_row,ptot_row,cs_row,cv_row,dpt_row,cp_row,eosfail)
      CALL helmeos
      IF (rank.EQ.MASTER.AND.eosfail.EQV..true.) PRINT*,'EOScorr fail'
!
      presion  = ptot_row(1)
      cvsp   = cv_row(1)
      cpsp   = cp_row(1)
      dPdTp  = dpt_row(1)
      gama   = nabad_gas_row(1)
!
      END SUBROUTINE deg
