         SUBROUTINE gw_radiation(ipas)
!
!--Load modules
!
         USE mod_essentials
         USE mod_parameters, ONLY : maxpas
         USE mod_commons, ONLY : xyzhm, vxyzut, axyzut
!
!--Force to declare EVERYTHING
!
         IMPLICIT  NONE
!
!--I/O variables
!
         INTEGER, INTENT(IN) :: ipas
!
!--Local variables
!
         REAL(8), DIMENSION(ndim,ndim,ndim,ndim) :: pijkl
         REAL(8), DIMENSION(ndim,ndim,maxpas) :: ddi
         REAL(8), DIMENSION(ndim,ndim) :: deltaij, dddi
         REAL(8), DIMENSION(ndim) :: ni
         REAL(8) :: tpas,tpasant,tbefore,trdddi,hp00,hx00,hpp0,hxp0,    &
                    degdt,dddi2,time1
         INTEGER :: i,j,k,l,k2,p
!
!--Physical parameters in code units 
!  
!  c = 3x10^10 cm/s = 215.49
!  G = 6.67x10^{-8} cm^3 g^-1 s^-2 = 1.0
!  distobs = 10kpc = 3.08x10E+22 cm
!
         REAL(8), PARAMETER :: Gdivc4=4.64d-10, distobs=4.43d+12,       &
                               Gdiv5c5=4.30d-13
!
!--cgs-code conversion factors
!
         REAL(8), PARAMETER :: potfact=7.7d+47, velfact=1.4d+8,         &
                               lumisol=3.9d+33
!
!--Open I/O files
!
         OPEN(UNIT=1,FILE='degdt',STATUS='unknown',ACCESS='append')
         OPEN(UNIT=2,FILE='h+00', STATUS='unknown',ACCESS='append')
         OPEN(UNIT=3,FILE='hx00', STATUS='unknown',ACCESS='append')
         OPEN(UNIT=4,FILE='h+p0', STATUS='unknown',ACCESS='append')
         OPEN(UNIT=5,FILE='hxp0', STATUS='unknown',ACCESS='append')
!
         IF (ipas.EQ.1) THEN
            WRITE(1,*) 'seconds     degdt(Lsol)'
            WRITE(2,*) 'seconds     h+00'
            WRITE(3,*) 'seconds     hx00'
            WRITE(4,*) 'seconds     h+p0'
            WRITE(5,*) 'seconds     hxp0'
         ENDIF
!
!--Start the GW calculation.Construct Kronecker's delta
!
         DO i=1,3
            DO j=1,3
               deltaij(i,j) = 0.0d0
            ENDDO
         ENDDO
         DO i=1,3
            deltaij(i,i) = 1.0d0
         ENDDO
!
!--Construct PIJKL and NI tensors
!
         DO i=1,3
            ni(i) = deltaij(3,i)
         ENDDO
         DO i=1,3
            DO j=1,3
               DO k=1,3
                  DO l=1,3
                     pijkl(i,j,k,l) = (deltaij(i,k)-ni(i)*ni(k))*       &
                       (deltaij(j,l)-ni(j)*ni(l))-                      &
                       0.5d0*(deltaij(i,j)-ni(i)*ni(j))*                &
                       (deltaij(k,l)-ni(k)*ni(l))         
                  ENDDO
               ENDDO
            ENDDO
        ENDDO
!
        tbefore = tnow-dtmp_min
        tpasant = tbefore
        tpas    = tnow
        time1   = tpas*50.4d0
!
!--Variables startup
!    
        DO k=1,3
           DO k2=1,3
              ddi(k,k2,ipas) = 0.0d0
              dddi(k,k2)     = 0.0d0
           ENDDO
        ENDDO
        trdddi = 0.0d0
!
!--DDI calculation
!
        DO k=1,3
           DO l=1,3
              DO p=1,nbody
                 ddi(k,l,ipas) = ddi(k,l,ipas) + xyzhm(5,p)*(2.0d0*     &
                 vxyzut(k,p)*vxyzut(l,p) + xyzhm(k,p)*axyzut(l,p) +     &
                 xyzhm(l,p)*axyzut(k,p))
              ENDDO
           ENDDO
        ENDDO
!
!--H+ & Hx calculation
!
        hp00 = (Gdivc4)*(ddi(1,1,ipas)-ddi(2,2,ipas))/distobs
        hx00 = (Gdivc4)*(2.0d0*ddi(1,2,ipas))/distobs
        hpp0 = (Gdivc4)*(ddi(3,3,ipas)-ddi(2,2,ipas))/distobs
        hxp0 = (Gdivc4)*(-2.0d0*ddi(2,3,ipas))/distobs
!
!--Construct DI derivatives and the radiated energy
!
        IF (ipas.GT.1) THEN
           DO i=1,3
              DO j=1,3
                 dddi(i,j) = (ddi(i,j,ipas)-ddi(i,j,ipas-1))/           &
                           (tpas-tpasant)
              ENDDO
           ENDDO
!
           DO i=1,3
              trdddi = trdddi + dddi(i,i)
           ENDDO
!
           DO i=1,3
              DO j=1,3
                 dddi(i,j) = dddi(i,j)-(0.3333d0)*deltaij(i,j)*trdddi
              ENDDO
           ENDDO
!
           dddi2=0.0d0
           DO i=1,3
              DO j=1,3 
                 dddi2 = dddi2+dddi(i,j)*dddi(i,j)
              ENDDO
           ENDDO
           degdt = Gdiv5c5*dddi2
!
!--Code units-cgs conversion
!
           degdt = potfact*degdt
           DO i=1,3
              DO j=1,3
                 dddi(i,j) = velfact*dddi(i,j)
              ENDDO
           ENDDO    
!
!--End of GW calculation. Write results and end program.
!
           WRITE(1,'(2(1pe12.4))') time1,degdt/lumisol
!
        ENDIF
!
        WRITE(2,'(2(1pe12.4))') time1, hp00
        WRITE(3,'(2(1pe12.4))') time1, hx00
        WRITE(4,'(2(1pe12.4))') time1, hpp0
        WRITE(5,'(2(1pe12.4))') time1, hxp0
!
        CLOSE(1)
        CLOSE(2)
        CLOSE(3)
        CLOSE(4)
        CLOSE(5)
!
        CALL flush (1)
        CALL flush (2)
        CALL flush (3)
        CALL flush (4)
        CALL flush (5)
!
        END SUBROUTINE gw_radiation
