      subroutine neutrins

      include 'b3defs.f'
      real*8 t9,a(14),zn(14),enrg2,ui,drho(1),du(1),dt(1),
     $      dycc(1,17),ya(14),dgam,dp(1),ds(1),dep(1),ye,duu
      integer i,j,k,iqn,input,iout,ipart
      parameter (iqn=17)
      data a/4.0d0,12.0d0,16.0d0,20.0d0,24.0d0,28.0d0,32.0d0,
     &     36.0d0,40.0d0,44.0d0,48.0d0,52.0d0,56.0d0,60.0d0/
      data zn/2.0d0,6.0d0,8.0d0,10.0d0,12.0d0,14.0d0,16.0d0,18.0d0,
     &     20.0d0,22.0d0,24.0d0,26.0d0,28.0d0,30.0d0/

      do i=1,nbody

       input=1
       iout=1
       ipart=1
       drho(1)=rho(i)*unm/unl/unl/unl
       if(temp(i).lt.0.0d0) temp(i)=1.0d4 
       dt(1)=temp(i)
       duu=uint(i) 
       du(1)=duu*drho(1)*uen/unm   
       do j=1,nel
           ya(j)=x(i,j)/a(j)
       enddo
         ye=0.0d0
         ui=0.0d0
         duu=0.0d0
       do 20 j=1,nel
           ye=ye+zn(j)*ya(j)
 20    continue
       dycc(1,1)=0.5d0
       dycc(1,2)=0.0d0
       dycc(1,3)=0.0d0
       do j=4,iqn
         dycc(1,j)=ya(j-3)
       enddo

       call eosda(drho,dt,du,dudt,dp,ds,dycc,dgam,1,1,
     1     input,iout,ipart,dep)
       duu=du(1)/drho(1)
       ui=duu
       cv(i)=dudt(1)/uen*unm
       pres(i)=dp(1)/uen*unl*unl*unl
       enrg2=0.0d0
       t9=temp(i)  
       call ntrno(drho,t9,ye,a,zn,enrg2)            

         if(enrg2.ne.0.0d0)then
            input=2
            iout=1
            ipart=1
            ui = ui - enrg2*(dtmp_min*50.0d0)
c           ui = ui*uen
            drho(1)=rho(i)*unm/unl**3
            du(1)=drho(1)*ui
            dt(1)=t9
            dycc(1,1)=0.5d0
            dycc(1,2)=0.0d0
            dycc(1,3)=0.0d0
            do j=4,iqn
               dycc(1,j)=ya(j-3)
            enddo
            call eosda(drho,dt,du,dudt,dp,ds,dycc,dgam,1,1,
     1           input,iout,ipart,dep)
            pres(i)=dp(1)/uen*unl**3
            cv(i)=dudt(1)/uen*unm
            enrg2=enrg2*(dsqrt(unt2)*unm/uen)
            eneutr(i)=enrg2*dtmp_min
         endif
      enddo

11    continue

      return
      end

C----------------BEGINING OF NEUTRINOS PART-----------------------------

      subroutine ntrno(rhok,t9,ye,a,z,energ)
C=======================================================================
C
C       Esta subrutina calcula el ritmo de perdida de energia por 
C       neutrinos en erg/g/s.
C
C       Refs.: Itoh, N., Adachi, A., Nakagawa, M., Kohyama, Y., y 
C              Munakata, H., Ap. J., 339, 354 (1989). 
C              Itoh, N., Kohyama, Y., Ap. J., 275, 858, (1983). 
C              Munakata, K., Kohyama, Y., Itoh, N., Ap. J., 316, 708 
C              (1987). 
C
C       Revisada el 21.09.94 por Enrique Garcia-Berro
C       
C-----------------------------------------------------------------------
C
C       Parametros de entrada:
C
C       rho:    densidad 
C       ye:     numero electronico molar
C       t:      temperatura
C       z:      carga atomica
C       z2oa:   <Z^2/A>
C       hntrno: metodo de calculo, japoneses es 'j'
C
C-----------------------------------------------------------------------
C
C       Parametros de salida:
C
C       entrno: energia emitida por neutrinos (erg/g/s)
C
C=======================================================================
      implicit double precision (a-h,o-z)

C     ---   Declaracion de variables   ---

      dimension z(14),z2oa(14),a(14)
      real*8 rhok,t9,energ

       rho=rhok 
       do i=1,14
         z2oa(i)=z(i)**2/a(i)
       enddo
       t=t9
       energ=0.0d0 
       entrno1=0.0d0
       entrno2=0.0d0
       entrno3=0.0d0 

C     ---   Controles de temperatura   ---

       do i=1,14

         if (t .le. 1.0d7 .or. t .gt. 1.0d11) then
         energ=0.0d0
         return
         else
         continue
         end if

C     ---   Controles de densidad   ---

         if (rho .le. 1.0d0 .or. rho .gt. 1.0d14) then
         energ=0.0d0
         return
         else
         continue
         end if

C     ---   Ritmos de los japoneses   ---
         za=z(i)
         call nppps(rho,t,ye,enppps)
         call nbrms(rho,t,ye,z,enbrms)
         entrno1=enppps+enbrms
        
C     ---   Ritmos de Beaudet, Petrosian y Salpeter   ---

        call neutri(rho,t,ye,z,entrno2)

C     ---   Ritmos de I2Jr   ---

        call eweak(rho,t,ye,z2oa,entrno3,denudt,denudr)

      energ=energ+entrno1+entrno2+entrno3
      enddo
c
      return 
c
      end


C***********************************************************************
C***********************************************************************


      subroutine nppps(rho,t,ye,enppps)
C=======================================================================
C
C       Esta subrutina calcula el ritmo de perdida de energia por 
C       neutrinos en regimen de plasma neutrino, pair neutrino y 
C       photo neutrino en erg/g/s.
C
C       Refs.: Itoh, N., Adachi, A., Nakagawa, M., Kohyama, Y., y 
C              Munakata, H., Ap. J., 339, 354 (1989). 
C              Corrections appeared in Ap.J. 360:741 (1990) 
C              have been included.
C
C-----------------------------------------------------------------------
C
C       Nota:  se han tomado los siguientes parametros, 
C 
C              sin^2 zw= 0.23
C              ca      = 0.50
C              cv      = 0.50 + 2 sin 2 zw = 0.96
C              cpv     = 1 - cv = 0.04
C              cpa     = 1 - ca = 0.50
C              n           =  2
C              cv2 + ca2   =  1.1716
C              cv2 - ca2   =  0.6716
C              cpv2 - cpa2 = -0.2484
C              cpv2 + cpa2 =  0.2516
C
C-----------------------------------------------------------------------
C
C       Revisada el 21.10.94 por Enrique Garcia-Berro.
C
C-----------------------------------------------------------------------
C
C       Parametros de entrada:
C
C       rho:    densidad 
C       ye:     numero electronico molar
C       t:      temperatura
C
C-----------------------------------------------------------------------
C
C       Parametros de salida:
C
C       enppps: energia emitida por neutrinos (erg/g/s)
C
C=======================================================================
      implicit double precision (a-h,o-z)

C     ---   Auxiliary variables (photo neutrino process)   ---

      dimension ci7(3,7),ci8(3,7),ci9(3,7),
     &          di7(3,5),di8(3,5),di9(3,5)

C     ---   Auxiliary variables (pair neutrino process)   ---

      dimension bii1(2),bii2(2),bii3(2),cii(2)

C     ---   Auxiliary variables (plasma neutrino process)   ---

      dimension alfiii(4),betiii(4),gamiii(4),biiib(3)

C     ---   Parameters   ---

      parameter (pi=3.141592654d0)

C     ---   Data for the photoneutrino process   ---

      data bi1 /6.290d-3/,
     &     bi2 /7.483d-3/,
     &     bi3 /3.061d-4/

      data ci7 / 1.008d11, 8.156d10, 1.067d11,
     &           0.000d00, 9.728d08,-9.782d09,
     &           0.000d00,-3.806d09,-7.193d09,
     &           0.000d00,-4.384d09,-6.936d09,
     &           0.000d00,-5.774d09,-6.893d09,
     &           0.000d00,-5.249d09,-7.041d09,
     &           0.000d00,-5.153d09,-7.193d09/

      data ci8 / 9.889d10, 1.813d11, 9.750d10,
     &          -4.524d08,-7.556d09, 3.484d10,
     &          -6.088d06,-3.304d09, 5.199d09,
     &           4.269d07,-1.031d09,-1.695d09,
     &           5.172d07,-1.764d09,-2.865d09,
     &           4.910d07,-1.851d09,-3.395d09,
     &           4.388d07,-1.928d09,-3.418d09/

      data ci9 / 9.581d10, 1.459d12, 2.424d11,
     &           4.107d08, 1.314d11,-3.669d09,
     &           2.305d08,-1.169d11,-8.691d09,
     &           2.236d08,-1.765d11,-7.967d09,
     &           1.580d08,-1.867d11,-7.932d09,
     &           2.165d08,-1.983d11,-7.987d09,
     &           1.721d08,-1.896d11,-8.333d09/

      data di7 / 0.0d0,-1.879d10,-2.919d10,
     &           0.0d0,-9.667d09,-1.185d10,
     &           0.0d0,-5.602d09,-7.270d09,
     &           0.0d0,-3.370d09,-4.222d09,
     &           0.0d0,-1.825d09,-1.560d09/

      data di8 /-1.135d08, 1.652d09,-1.549d10,
     &           1.256d08,-3.119d09,-9.338d09,
     &           5.149d07,-1.839d09,-5.899d09,
     &           3.436d07,-1.458d09,-3.035d09,
     &           1.005d07,-8.956d08,-1.598d09/

      data di9 / 4.724d08,-7.094d11,-2.254d10,
     &           2.976d08,-3.697d11,-1.551d10,
     &           2.242d08,-2.189d11,-7.793d09,
     &           7.937d07,-1.273d11,-4.489d09,
     &           4.859d07,-5.705d10,-2.185d09/

C     ---   Data for the pair neutrino process   ---

      data aii0 /6.002d19/,
     &     aii1 /2.084d20/,
     &     aii2 /1.872d21/

      data bii1 / 9.383d-1, 1.2383d0/,
     &     bii2 /-4.141d-1,-0.8141d0/,
     &     bii3 / 5.829d-2, 0.0000d0/

      data cii / 5.5924d0, 4.9924d0/

C     ---   Data for the plasma neutrino process   ---

      data aiii0a /2.320d-7/,aiii1a /8.449d-8/,
     &     aiii2a /1.787d-8/,aiii3a /0.000d+0/,
     &     biii1a /2.581d-2/,biii2a /1.734d-2/,
     &     biii3a /6.990d-4/,ciiia  /0.56457d0/

      data alfiii / 2.787d-7, 5.408d-7,-1.812d-7, 1.639d-8/,
     &     betiii /-3.936d-6,-7.715d-6, 2.731d-6,-2.445d-7/,
     &     gamiii / 1.408d-5, 2.751d-5,-1.024d-5, 9.137d-7/

      data ciiib1 /-7.683d-1/,
     &     ciiib2 / 1.798d-1/
      data biiib / 1.148d-2,1.209d-1,2.432d-4/

C     ---   General definitions   ---

      al =1.6862838d-10*t
      ali=1.0d0/al
      tl =dlog10(t)
      csi=ali*(rho*ye*1.0d-9)**(1.0d0/3.0d0)

C     ---   (I) Photoneutrino process   ---

         if (tl .lt. 8.0d0) then
         c=0.5654d0+tl-7.0d0
         else
         c=1.5654d0
         end if

         if (tl .ge. 9.0d0) then
         go to 4
         else
         continue
         end if

         if (tl .ge. 8.0d0) then
         go to 2
         else
         continue
         end if

      tau=tl-7.0d0
      aux=cos(1.0d1*pi*tau)
      ai0=0.5d0*(ci7(1,1)+ci7(1,7)*aux)
      ai1=0.5d0*(ci7(2,1)+ci7(2,7)*aux)
      ai2=0.5d0*(ci7(3,1)+ci7(3,7)*aux)

         do 1 j=2,6
         aux1=cos(5.0d0*pi*(j-1)*tau/3.0d0)
         aux2=sin(5.0d0*pi*(j-1)*tau/3.0d0)
         ai0=ai0+ci7(1,j)*aux1+di7(1,j-1)*aux2
         ai1=ai1+ci7(2,j)*aux1+di7(2,j-1)*aux2
         ai2=ai2+ci7(3,j)*aux1+di7(3,j-1)*aux2
    1    continue

      go to 6

    2 continue
      tau=tl-8.0d0
      aux=cos(1.0d1*pi*tau)
      ai0=0.5d0*(ci8(1,1)+ci8(1,7)*aux)
      ai1=0.5d0*(ci8(2,1)+ci8(2,7)*aux)
      ai2=0.5d0*(ci8(3,1)+ci8(3,7)*aux)

         do 3 j=2,6
         aux1=cos(5.0d0*pi*(j-1)*tau/3.0d0)
         aux2=sin(5.0d0*pi*(j-1)*tau/3.0d0)
         ai0=ai0+ci8(1,j)*aux1+di8(1,j-1)*aux2
         ai1=ai1+ci8(2,j)*aux1+di8(2,j-1)*aux2
         ai2=ai2+ci8(3,j)*aux1+di8(3,j-1)*aux2
    3    continue

      go to 6

    4 continue

      tau=tl-9.0d0
      aux=cos(1.0d1*pi*tau)
      ai0=0.5d0*(ci9(1,1)+ci9(1,7)*aux)
      ai1=0.5d0*(ci9(2,1)+ci9(2,7)*aux)
      ai2=0.5d0*(ci9(3,1)+ci9(3,7)*aux)
         do 5 j=2,6
         aux1=cos(5.0d0*pi*(j-1)*tau/3.0d0)
         aux2=sin(5.0d0*pi*(j-1)*tau/3.0d0)
         ai0=ai0+ci9(1,j)*aux1+di9(1,j-1)*aux2
         ai1=ai1+ci9(2,j)*aux1+di9(2,j-1)*aux2
         ai2=ai2+ci9(3,j)*aux1+di9(3,j-1)*aux2
    5    continue

    6 continue

      f=((ai0+csi*(ai1+csi*ai2))*dexp(-c*csi))/
     &  (csi*csi*csi+((bi3*ali+bi2)*ali+bi1)*ali)

      q=0.666d0*(1.0d0+2.045d0*al)**(-2.066d0)/(1.0d0+rho*ye/
     &  ((((-1.604d8*al+8.499d8)*al+1.653d8)*al+1.875d8)*al))

      enui=0.8374d0*(1.0d0-0.1044d0*q)*(rho*ye)*f*al**5.0d0

C     --- (II) Pair neutrino process   ---

      nn=1+amin0(1,idint(t*1.0d-10))

      g=(((918.6d0*al*al+1534.0d0)*al*al+133.5d0)*
     &  al*al-13.04d0)*al*al+1.0d0

      f=(aii0+csi*(aii1+aii2*csi))*dexp(-csi*cii(nn))/
     &  (csi*csi*csi+((bii3(nn)*ali+bii2(nn))*ali+bii1(nn))*ali)

      q=1.0d0/(1.0050d0+0.3967d0*al**0.5d0+10.7480d0*al*al)*
     & (1.0d0+rho*ye/(7.692d7*al*al*al+9.715d6*al**0.5d0))**(-0.3d0)

      enuii=0.8374d0*(1.0d0+0.1044d0*q)*g*f*dexp(-2.0d0*ali)

C     ---   (III) Plasma neutrino process   ---

         if (tl .ge. 7.8d0) then
        
         aiii0=aiii0a
         aiii1=aiii1a
         aiii2=aiii2a
         aiii3=aiii3a

         biii1=biii1a
         biii2=biii2a
         biii3=biii3a

         ciii=ciiia

         else

         aiii0=(alfiii(1)*tl+betiii(1))*tl+gamiii(1)
         aiii1=(alfiii(2)*tl+betiii(2))*tl+gamiii(2)
         aiii2=(alfiii(3)*tl+betiii(3))*tl+gamiii(3)
         aiii3=(alfiii(4)*tl+betiii(4))*tl+gamiii(4)

         biii1=biiib(1)
         biii2=biiib(2)
         biii3=biiib(3)

         ciii=ciiib1+ciiib2*tl

         end if

      f=((aiii0+csi*(aiii1+csi*(aiii2+csi*aiii3)))*dexp(-ciii*csi))/
     &  (csi*csi*csi+((biii3*ali+biii2)*ali+biii1)*ali)

      enuiii=0.92480d0*f*(rho*ye)**3.0d0

C     ---   Total rate   ---
      enppps=(enui+enuii+enuiii)/rho

      return
      end

C***********************************************************************
C***********************************************************************


      subroutine nbrms(rho,t,ye,z,enbrms)
C=======================================================================
C
C       Esta subrutina calcula el ritmo de perdida de energia por 
C       neutrinos en regimen de neutrino bremsstrahlung en erg/g/s.
C
C       Refs.: Itoh, N., Adachi, A., Nakagawa, M., Kohyama, Y., y 
C              Munakata, H., Ap. J., 339, 354 (1989). 
C
C-----------------------------------------------------------------------
C
C       Nota:  se han tomado los siguientes parametros, 
C 
C              sin^2 zw= 0.23
C              ca      = 0.50
C              cv      = 0.50 + 2 sin 2 zw = 0.96
C              cpv     = 1 - cv = 0.04
C              cpa     = 1 - ca = 0.50
C              n           =  2
C              cv2 + ca2   =  1.1716
C              cv2 - ca2   =  0.6716
C              cpv2 - cpa2 = -0.2484
C              cpv2 + cpa2 =  0.2516
C
C-----------------------------------------------------------------------
C
C       Revisada el 15.10.91 por Abulafia
C
C-----------------------------------------------------------------------
C
C       Parametros de entrada:
C
C       rho:    densidad 
C       ye:     numero electronico molar
C       t:      temperatura
C       z:      carga atomica
C
C-----------------------------------------------------------------------
C
C       Parametros de salida:
C
C       enbrms: energia emitida por neutrinos bremsstrahlung (erg/g/s)
C
C=======================================================================
      implicit double precision (a-h,o-z)

C     ---   Parametros   ---

      parameter (pi=3.141592654d0)

C     ---   Datos para la emision de neutrinos por bremsstrahlung   ---

      dimension za(9)

      data za/2.0d0,6.0d0,8.0d0,
     &        1.0d1,1.2d1,1.4d1,
     &        1.6d1,2.0d1,2.8d1/

C     ---   Definiciones generales   ---

      ul=2.0d0*pi*((dlog10(rho)-3.0d0)/1.0d1)
      us=2.0d0*pi*((dlog10(rho)-3.0d0)/9.0d0)
      gamma=(2.275d5/t)*(z**(5.0d0/3.0d0))*((rho*ye)**(1.0d0/3.0d0))
      tbdpn=1.779d9*(dsqrt(1.0d0+1.018d-4*((rho*ye)**(2.0d0/3.0d0)))-
     &      1.0d0)

C     ---   Bremsstrahlung parcialmente degenerado   ---

         if (t .gt. tbdpn) then
         call nbrpd(t,rho,f,g)
         go to 3
         else
         continue
         end if

C     ---   Fase solida   ---

         if (gamma .gt. 1.78d2) then

C     ---   Helio y Z < 2   ---

            if (z .le. za(1)) then
            call nsolid(us,1,gamma,f,g)
            else
            continue
            end if

C     ---   Hierro y Z > 28   ---

            if (z .ge. za(9)) then
            call nsolid(us,9,gamma,f,g)
            else
            continue
            end if

C     ---   Elementos entre el He y el Fe   ---

            if (z .gt. za(1) .and. z .lt. za(9)) then

C     ---   Interpolacion   ---

               do 1 i=1,8

C     ---   Busqueda del punto de interpolacion   ---

                  if (z .gt. za(i) .and. z .le. za(i+1)) then

C     ---   Calculo de los coeficientes de interpolacion   ---

                  call nsolid(us,i,gamma,fd,gd)
                  call nsolid(us,i+1,gamma,fu,gu)

C     ---   Calculo de la contribucion al ritmo   ---

                  dfdz=(fu-fd)/(za(i+1)-za(i))
                  dgdz=(gu-gd)/(za(i+1)-za(i))
                  f=fd+dfdz*(z-za(i))
                  g=gd+dgdz*(z-za(i))

                  else
                  continue
                  end if

    1          continue

            else
            continue
            end if

         else
         continue
         end if

C     ---   Fase liquida   ---

         if (gamma .le. 1.78d2) then

C     ---   Helio y Z < 2   ---

            if (z .le. za(1)) then
            call nliquid(ul,1,gamma,f,g)
            else
            continue
            end if

C     ---   Hierro y Z > 28   ---

            if (z .ge. za(9)) then
            call nliquid(ul,9,gamma,f,g)
            else
            continue
            end if
C     ---   Elementos entre el He y el Fe   ---

            if (z .gt. za(1) .and. z .lt. za(9)) then

C     ---   Interpolacion   ---

               do 2 i=1,8

C     ---   Busqueda del punto de interpolacion   ---

                  if (z .gt. za(i) .and. z .le. za(i+1)) then

C     ---   Calculo de los coeficientes de interpolacion   ---

                  call nliquid(ul,i,gamma,fd,gd)
                  call nliquid(ul,i+1,gamma,fu,gu)

C     ---   Calculo de la contribucion al ritmo   ---

                  dfdz=(fu-fd)/(za(i+1)-za(i))
                  dgdz=(gu-gd)/(za(i+1)-za(i))
                  f=fd+dfdz*(z-za(i))
                  g=gd+dgdz*(z-za(i))
                  else
                  continue
                  end if

    2          continue

            else
            continue
            end if

         else
         continue
         end if

    3 continue

C     ---   Energia emitida por neutrino bremsstrahlung   ---

      enbrms=0.5738d0*(0.8374d0*f+0.0874d0*g)*z*ye*(t*1.0d-8)**6.0d0

      return
      end

C***********************************************************************
C***********************************************************************

      subroutine nliquid(u,j,gamma,fliqd,gliqd)
C=======================================================================
C
C       Esta subrutina calcula la contribucion al ritmo de perdida de 
C       energia por neutrinos en la zona liquida.
C
C       Refs.: Itoh, N., Kohyama, Y., Ap. J., 275, 858, (1983). 
C
C       Revisada el 23.09.91 por Abulafia
C
C-----------------------------------------------------------------------
C
C       Parametros de entrada:
C
C       u:     parametro adimensional.
C       j:     especie atomica
C
C
C              1 ........   4He
C              2 ........   12C
C              3 ........   16O
C              4 ........   20Ne
C              5 ........   24Mg
C              6 ........   28Si
C              7 ........   32S
C              8 ........   40Ca
C              9 ........   56Fe
C
C       gamma: constante de acoplamiento electronica. 
C
C-----------------------------------------------------------------------
C
C       Parametros de salida:
C
C       fliqd: coeficiente f en fase solida.
C       gliqd: coeficiente g en fase solida.
C
C=======================================================================
      implicit double precision (a-h,o-z)
C     ---   Declaracion de variables   ---

      double precision i0,i1,i2,i3,i4,i5,j1,j2,j3,j4

C     ---   Dimensionado de variables   ---

      dimension a0(9),a1(9),a2(9),a3(9),a4(9),a5(9),
     &          e0(9),e1(9),e2(9),e3(9),e4(9),e5(9),
     &          i0(9),i1(9),i2(9),i3(9),i4(9),i5(9),
     &          p0(9),p1(9),p2(9),p3(9),p4(9),p5(9),
     &                      f1(9),f2(9),f3(9),f4(9),
     &                      b1(9),b2(9),b3(9),b4(9),
     &                      j1(9),j2(9),j3(9),j4(9),
     &                      q1(9),q2(9),q3(9),q4(9),
     &                                    c(9),d(9),
     &                                    g(9),h(9),
     &                                    k(9),l(9),
     &                                    r(9),s(9),
     &              alf0(9),alf1(9),alf2(9),alf3(9),
     &              bet0(9),bet1(9),bet2(9),bet3(9)

C     ---   Datos generales   ---
      data a0/ 0.09037d0, 0.17946d0, 0.20933d0, 0.23425d0, 0.25567d0,
     &         0.27445d0, 0.29120d0, 0.32001d0, 0.34888d0/
      data a1/-0.03009d0,-0.05821d0,-0.06740d0,-0.07503d0,-0.08158d0,
     &        -0.08734d0,-0.09250d0,-0.10142d0,-0.11076d0/
      data a2/-0.00564d0,-0.01089d0,-0.01293d0,-0.01472d0,-0.01632d0,
     &        -0.01777d0,-0.01910d0,-0.02145d0,-0.02349d0/
      data a3/-0.00544d0,-0.01147d0,-0.01352d0,-0.01522d0,-0.01667d0,
     &        -0.01793d0,-0.01905d0,-0.02094d0,-0.02283d0/
      data a4/-0.00290d0,-0.00656d0,-0.00776d0,-0.00872d0,-0.00952d0,
     &        -0.01019d0,-0.01076d0,-0.01168d0,-0.01250d0/
      data a5/-0.00224d0,-0.00519d0,-0.00613d0,-0.00688d0,-0.00748d0,
     &        -0.00798d0,-0.00839d0,-0.00902d0,-0.00971d0/

      data b1/-0.02148d0,-0.04969d0,-0.05950d0,-0.06776d0,-0.07491d0,
     &        -0.08120d0,-0.08682d0,-0.09651d0,-0.10661d0/
      data b2/-0.00817d0,-0.01584d0,-0.01837d0,-0.02045d0,-0.02220d0,
     &        -0.02370d0,-0.02500d0,-0.02716d0,-0.02860d0/
      data b3/-0.00300d0,-0.00504d0,-0.00567d0,-0.00616d0,-0.0653d0,
     &        -0.00683d0,-0.00706d0,-0.00738d0,-0.00785d0/
      data b4/-0.00170d0,-0.00281d0,-0.00310d0,-0.00331d0,-0.00345d0,
     &        -0.00356d0,-0.00363d0,-0.00370d0,-0.00385d0/

      data  c/ 0.00671d0, 0.00945d0, 0.00952d0, 0.00932d0, 0.00899d0,
     &         0.00858d0, 0.00814d0, 0.00721d0, 0.00766d0/

      data  d/ 0.28130d0, 0.34529d0, 0.36029d0, 0.37137d0, 0.38006d0,
     &         0.38714d0, 0.39309d0, 0.40262d0, 0.40991d0/

      data e0/-0.02006d0, 0.06781d0, 0.09304d0, 0.11465d0, 0.13455d0,
     &         0.15315d0, 0.17049d0, 0.20051d0, 0.23159d0/
      data e1/ 0.01790d0,-0.00944d0,-0.01656d0,-0.02253d0,-0.02828d0,
     &        -0.03391d0,-0.03930d0,-0.04877d0,-0.05891d0/
      data e2/-0.00783d0,-0.01289d0,-0.01489d0,-0.01680d0,-0.01846d0,
     &        -0.01988d0,-0.02113d0,-0.02331d0,-0.02531d0/
      data e3/-0.00021d0,-0.00589d0,-0.00778d0,-0.00942d0,-0.01087d0,
     &        -0.01218d0,-0.01338d0,-0.01541d0,-0.01747d0/
      data e4/ 0.00024d0,-0.00404d0,-0.00520d0,-0.00613d0,-0.00693d0,
     &        -0.00763d0,-0.00825d0,-0.00925d0,-0.01021d0/
      data e5/-0.00014d0,-0.00330d0,-0.00418d0,-0.00488d0,-0.00547d0,
     &        -0.00596d0,-0.00639d0,-0.00705d0,-0.00778d0/

      data f1/ 0.00538d0,-0.02213d0,-0.03076d0,-0.03824d0,-0.04508d0,
     &        -0.05142d0,-0.05729d0,-0.06744d0,-0.07767d0/
      data f2/-0.00175d0,-0.01136d0,-0.01390d0,-0.01601d0,-0.01782d0,
     &        -0.01941d0,-0.02082d0,-0.02311d0,-0.02516d0/
      data f3/-0.00346d0,-0.00467d0,-0.00522d0,-0.00571d0,-0.00607d0,
     &        -0.00631d0,-0.00648d0,-0.00670d0,-0.00711d0/
      data f4/-0.00031d0,-0.00131d0,-0.00161d0,-0.00183d0,-0.00200d0,
     &        -0.00212d0,-0.00221d0,-0.00233d0,-0.00260d0/

      data  g/-0.02199d0,-0.02342d0,-0.02513d0,-0.02697d0,-0.02832d0,
     &        -0.02919d0,-0.02978d0,-0.03070d0,-0.03076d0/

      data  h/ 0.17300d0, 0.24819d0, 0.27480d0, 0.29806d0, 0.31541d0,
     &         0.32790d0, 0.33756d0, 0.35242d0, 0.36908d0/

      data i0/ 0.00192d0, 0.00766d0, 0.00951d0, 0.01103d0, 0.01231d0,
     &         0.01342d0, 0.01440d0, 0.01606d0, 0.01892d0/
      data i1/-0.00301d0,-0.00710d0,-0.00838d0,-0.00942d0,-0.01030d0,
     &        -0.01106d0,-0.01172d0,-0.01285d0,-0.01493d0/
      data i2/-0.00073d0,-0.00028d0,-0.00011d0, 0.00004d0, 0.00016d0,
     &         0.00027d0, 0.00037d0, 0.00055d0, 0.00125d0/
      data i3/ 0.00182d0, 0.00232d0, 0.00244d0, 0.00252d0, 0.00259d0,
     &         0.00264d0, 0.00269d0, 0.00276d0, 0.00262d0/
      data i4/ 0.00037d0, 0.00044d0, 0.00046d0, 0.00047d0, 0.00048d0,
     &         0.00049d0, 0.00050d0, 0.00051d0, 0.00055d0/
      data i5/ 0.00116d0, 0.00158d0, 0.00168d0, 0.00176d0, 0.00183d0,
     &         0.00188d0, 0.00193d0, 0.00201d0, 0.00209d0/

      data j1/ 0.01706d0, 0.02300d0, 0.02455d0, 0.02573d0, 0.02669d0,
     &         0.02748d0, 0.02816d0, 0.02927d0, 0.03034d0/
      data j2/-0.00753d0,-0.01078d0,-0.01167d0,-0.01236d0,-0.01291d0,
     &        -0.01338d0,-0.01379d0,-0.01445d0,-0.01519d0/
      data j3/ 0.00066d0, 0.00118d0, 0.00132d0, 0.00144d0, 0.00154d0,
     &         0.00162d0, 0.00169d0, 0.00181d0, 0.00204d0/
      data j4/-0.00060d0,-0.00089d0,-0.00097d0,-0.00103d0,-0.00108d0,
     &        -0.00112d0,-0.00116d0,-0.00122d0,-0.00135d0/

      data k/-0.01021d0,-0.01259d0,-0.01314d0,-0.01354d0,-0.01386d0,
     &       -0.01411d0,-0.01432d0,-0.01465d0,-0.01494d0/

      data l/ 0.06417d0, 0.07917d0, 0.08263d0, 0.08515d0, 0.08711d0,
     &        0.08869d0, 0.09001d0, 0.09209d0, 0.09395d0/

      data p0/-0.01112d0,-0.00769d0,-0.00700d0,-0.00649d0,-0.00583d0,
     &        -0.00502d0,-0.00415d0,-0.00255d0,-0.00011d0/
      data p1/ 0.00603d0, 0.00356d0, 0.00295d0, 0.00246d0, 0.00192d0,
     &         0.00132d0, 0.00070d0,-0.00042d0,-0.00222d0/
      data p2/-0.00149d0,-0.00184d0,-0.00184d0,-0.00183d0,-0.00179d0,
     &        -0.00173d0,-0.00166d0,-0.00151d0,-0.00104d0/
      data p3/ 0.00047d0, 0.00146d0, 0.00166d0, 0.00181d0, 0.00193d0,
     &         0.00203d0, 0.00211d0, 0.00224d0, 0.00225d0/
      data p4/ 0.00040d0, 0.00031d0, 0.00032d0, 0.00033d0, 0.00034d0,
     &         0.00034d0, 0.00034d0, 0.00034d0, 0.00037d0/
      data p5/ 0.00028d0, 0.00069d0, 0.00082d0, 0.00093d0, 0.00102d0,
     &         0.00110d0, 0.00116d0, 0.00126d0, 0.00140d0/

      data q1/ 0.00422d0, 0.01052d0, 0.01231d0, 0.01379d0, 0.01501d0,
     &         0.01604d0, 0.01693d0, 0.01841d0, 0.02013d0/
      data q2/-0.00009d0,-0.00354d0,-0.00445d0,-0.00518d0,-0.00582d0,
     &        -0.00639d0,-0.00691d0,-0.00778d0,-0.00883d0/
      data q3/-0.00066d0,-0.00014d0, 0.00002d0, 0.00013d0, 0.00024d0,
     &         0.00035d0, 0.00045d0, 0.00062d0, 0.00090d0/
      data q4/-0.00003d0,-0.00018d0,-0.00026d0,-0.00033d0,-0.00038d0,
     &        -0.00044d0,-0.00048d0,-0.00056d0,-0.00071d0/

      data  r/-0.00561d0,-0.00829d0,-0.00921d0,-0.01000d0,-0.01059d0,
     &        -0.01103d0,-0.01136d0,-0.01188d0,-0.01249d0/

      data  s/ 0.03522d0, 0.05211d0, 0.05786d0, 0.06284d0, 0.06657d0,
     &         0.06928d0, 0.07140d0, 0.07468d0, 0.07850d0/

      data alf0/-0.07980d0,-0.05483d0,-0.06597d0,-0.06910d0,-0.07003d0,
     &          -0.07023d0,-0.06814d0,-0.06934d0,-0.07608d0/
      data alf1/ 0.17057d0,-0.01946d0, 0.06048d0, 0.07685d0, 0.07808d0,
     &           0.07660d0, 0.05799d0, 0.06443d0, 0.11559d0/
      data alf2/ 1.51980d0, 1.86310d0, 1.74860d0, 1.74280d0, 1.75870d0,
     &           1.77560d0, 1.81880d0, 1.82500d0, 1.75730d0/
      data alf3/-0.61058d0,-0.78873d0,-0.74308d0,-0.75047d0,-0.76675d0,
     &          -0.78191d0,-0.80866d0,-0.82010d0,-0.79677d0/

      data bet0/-0.05881d0,-0.06711d0,-0.07356d0,-0.07123d0,-0.06960d0,
     &          -0.06983d0,-0.06880d0,-0.07255d0,-0.08034d0/
      data bet1/ 0.00165d0, 0.06859d0, 0.10865d0, 0.08264d0, 0.06577d0,
     &           0.06649d0, 0.05775d0, 0.08529d0, 0.14368d0/
      data bet2/ 1.82700d0, 1.74360d0, 1.70150d0, 1.76760d0, 1.81180d0,
     &           1.82170d0, 1.84540d0, 1.81230d0, 1.73140d0/
      data bet3/-0.76993d0,-0.74498d0,-0.73653d0,-0.77896d0,-0.80797d0,
     &          -0.81839d0,-0.83439d0,-0.82499d0,-0.79467d0/

C     ---   Comienzo del calculo   ---

      gam3=1.0d0/(gamma**(1.0d0/3.0d0))

C     ---   Variables de interpolacion   ---

      v=alf0(j)+gam3*(alf1(j)+gam3*(alf2(j)+gam3*alf3(j)))
      w=bet0(j)+gam3*(bet1(j)+gam3*(bet2(j)+gam3*bet3(j)))

C     ---   Parametros de la interpolacion   ---

      u1=u
      u2=2.0d0*u
      u3=3.0d0*u
      u4=4.0d0*u
      u5=5.0d0*u

      fu0001=(a0(j)/2.0d0)+a1(j)*dcos(u1)+a2(j)*dcos(u2)+
     &                     a3(j)*dcos(u3)+a4(j)*dcos(u4)+
     &                     a5(j)*dcos(u5)+
     &                     b1(j)*dsin(u1)+b2(j)*dsin(u2)+
     &                     b3(j)*dsin(u3)+b4(j)*dsin(u4)+
     &                     c(j)*u+d(j)

      fu0160=(e0(j)/2.0d0)+e1(j)*dcos(u1)+e2(j)*dcos(u2)+
     &                     e3(j)*dcos(u3)+e4(j)*dcos(u4)+
     &                     e5(j)*dcos(u5)+
     &                     f1(j)*dsin(u1)+f2(j)*dsin(u2)+
     &                     f3(j)*dsin(u3)+f4(j)*dsin(u4)+
     &                     g(j)*u+h(j)

      gu0001=(i0(j)/2.0d0)+i1(j)*dcos(u1)+i2(j)*dcos(u2)+
     &                     i3(j)*dcos(u3)+i4(j)*dcos(u4)+
     &                     i5(j)*dcos(u5)+
     &                     j1(j)*dsin(u1)+j2(j)*dsin(u2)+
     &                     j3(j)*dsin(u3)+j4(j)*dsin(u4)+
     &                     k(j)*u+l(j)

      gu0160=(p0(j)/2.0d0)+p1(j)*dcos(u1)+p2(j)*dcos(u2)+
     &                     p3(j)*dcos(u3)+p4(j)*dcos(u4)+
     &                     p5(j)*dcos(u5)+
     &                     q1(j)*dsin(u1)+q2(j)*dsin(u2)+
     &                     q3(j)*dsin(u3)+q4(j)*dsin(u4)+
     &                     r(j)*u+s(j)

C     ---   Calculo de los coeficientes f y g en fase liquida   ---

      fliqd=v*fu0001+(1.0d0-v)*fu0160
      gliqd=w*gu0001+(1.0d0-w)*gu0160

      return
      end


C***********************************************************************
C***********************************************************************


     subroutine nsolid(u,j,gamma,fsold,gsold)
C=======================================================================
C
C       Esta subrutina calcula la contribucion al ritmo de perdida de 
C       energia por neutrinos en la zona solida (lattice + phonon).
C
C       Refs.: Itoh, N., Kohyama, Y., Matsumoto, N., Seki, M., Ap. J., 
C              285, 304, (1984). Las correcciones aparecidas en Ap. J.,
C              322, 584, (1987) estan incluidas.
C
C       Revisada el 23.09.91 por Abulafia
C
C-----------------------------------------------------------------------
C
C       Parametros de entrada:
C
C       u:     parametro adimensional.
C       j:     especie atomica
C
C
C              1 ........   4He
C              2 ........   12C
C              3 ........   16O
C              4 ........   20Ne
C              5 ........   24Mg
C              6 ........   28Si
C              7 ........   32S
C              8 ........   40Ca
C              9 ........   56Fe
C
C       gamma: constante de acoplamiento electronica. 
C
C-----------------------------------------------------------------------
C
C       Parametros de salida:
C
C       fsold: coeficiente f en fase solida.
C       gsold: coeficiente g en fase solida.
C
C=======================================================================
      implicit double precision (a-h,o-z)
C     ---   Declaracion de variables   ---

      double precision i0l,i1l,i2l,i3l,i4l,j1l,j2l,j3l,kcl,lcl
      double precision i0p,i1p,i2p,i3p,i4p,j1p,j2p,j3p,kcp,lcp

C     ---   Dimensionado de variables   ---

      dimension i0l(9),i1l(9),i2l(9),i3l(9),i4l(9),
     &          i0p(9),i1p(9),i2p(9),i3p(9),i4p(9),
     &          a0l(9),a1l(9),a2l(9),a3l(9),a4l(9),
     &          a0p(9),a1p(9),a2p(9),a3p(9),a4p(9),
     &               e0(9),e1(9),e2(9),e3(9),e4(9),
     &               p0(9),p1(9),p2(9),p3(9),p4(9),
     &                        b1l(9),b2l(9),b3l(9),
     &                        b1p(9),b2p(9),b3p(9),
     &                        j1l(9),j2l(9),j3l(9),
     &                        j1p(9),j2p(9),j3p(9),
     &                           f1(9),f2(9),f3(9),
     &                           q1(9),q2(9),q3(9),
     &                               kcl(9),lcl(9),
     &                               kcp(9),lcp(9),
     &                                 cl(9),dl(9),
     &                                 cp(9),dp(9),
     &                                   g(9),h(9),
     &                                   r(9),s(9),
     &         alf0l(9),alf1l(9),alf2l(9),alf3l(9),
     &         bet0l(9),bet1l(9),bet2l(9),bet3l(9),
     &         alf0p(9),alf1p(9),alf2p(9),alf3p(9),
     &         bet0p(9),bet1p(9),bet2p(9),bet3p(9)

C     ---   Datos generales   ---

      data a0l/-0.02296d0, 0.03677d0, 0.03232d0, 0.03224d0, 0.03191d0,
     &          0.03318d0, 0.03471d0, 0.03754d0, 0.04192d0/
      data a1l/ 0.01601d0,-0.01066d0,-0.00874d0,-0.01010d0,-0.01083d0,
     &         -0.01218d0,-0.01341d0,-0.01521d0,-0.01768d0/
      data a2l/-0.00433d0,-0.00458d0,-0.00413d0,-0.00285d0,-0.00211d0,
     &         -0.00145d0,-0.00098d0,-0.00052d0,-0.00007d0/
      data a3l/ 0.00015d0,-0.00177d0,-0.00190d0,-0.00193d0,-0.00192d0,
     &         -0.00197d0,-0.00204d0,-0.00218d0,-0.00241d0/
      data a4l/-0.00034d0,-0.00138d0,-0.00139d0,-0.00123d0,-0.00109d0,
     &         -0.00099d0,-0.00093d0,-0.00086d0,-0.00080d0/

      data b1l/ 0.01558d0,-0.00244d0,-0.00344d0,-0.00607d0,-0.00753d0,
     &         -0.00940d0,-0.01107d0,-0.01349d0,-0.01705d0/
      data b2l/ 0.00191d0,-0.00206d0,-0.00261d0,-0.00279d0,-0.00281d0,
     &         -0.00281d0,-0.00280d0,-0.00280d0,-0.00268d0/
      data b3l/-0.00055d0,-0.00037d0,-0.00070d0,-0.00078d0,-0.00087d0,
     &         -0.00093d0,-0.00100d0,-0.00115d0,-0.00141d0/

      data  cl/-0.01694d0,-0.01093d0,-0.00791d0,-0.00365d0,-0.00110d0,
     &          0.00123d0, 0.00304d0, 0.00531d0, 0.00818d0/

      data  dl/ 0.10649d0, 0.12431d0, 0.13980d0, 0.13861d0, 0.14075d0,
     &          0.13959d0, 0.13847d0, 0.13927d0, 0.13629d0/

      data e0/-0.03654d0, 0.04719d0, 0.04421d0, 0.05766d0, 0.06145d0,
     &         0.06964d0, 0.07593d0, 0.08106d0, 0.09256d0/
      data e1/ 0.02395d0,-0.01353d0,-0.00883d0,-0.01613d0,-0.01751d0,
     &        -0.02211d0,-0.02568d0,-0.02720d0,-0.03290d0/
      data e2/-0.00448d0,-0.00619d0,-0.00857d0,-0.00739d0,-0.00750d0,
     &        -0.00656d0,-0.00575d0,-0.00613d0,-0.00523d0/
      data e3/-0.00033d0,-0.00211d0,-0.00257d0,-0.00309d0,-0.00340d0,
     &        -0.00379d0,-0.00413d0,-0.00460d0,-0.00539d0/
      data e4/-0.00088d0,-0.00176d0,-0.00214d0,-0.00222d0,-0.00231d0,
     &        -0.00236d0,-0.00240d0,-0.00260d0,-0.00276d0/

      data f1/ 0.01730d0, 0.00456d0, 0.00629d0,-0.00079d0,-0.00322d0,
     &        -0.00863d0,-0.01332d0,-0.01665d0,-0.02574d0/
      data f2/ 0.00402d0,-0.00174d0,-0.00210d0,-0.00328d0,-0.00386d0,
     &        -0.00450d0,-0.00495d0,-0.00559d0,-0.00630d0/
      data f3/-0.00005d0,-0.00031d0,-0.00099d0,-0.00115d0,-0.00148d0,
     &        -0.00162d0,-0.00177d0,-0.00229d0,-0.00285d0/

      data  g/-0.02222d0,-0.02259d0,-0.02610d0,-0.01999d0,-0.01798d0,
     &        -0.01325d0,-0.00913d0,-0.00698d0, 0.00022d0/

      data  h/ 0.13969d0, 0.20343d0, 0.26993d0, 0.27099d0, 0.29179d0,
     &         0.29081d0, 0.28946d0, 0.31671d0, 0.31871d0/

      data i0l/-0.00647d0, 0.00106d0, 0.00199d0, 0.00362d0, 0.00468d0,
     &          0.00573d0, 0.00662d0, 0.00789d0, 0.00977d0/
      data i1l/ 0.00440d0,-0.00048d0,-0.00112d0,-0.00224d0,-0.00297d0,
     &         -0.00370d0,-0.00432d0,-0.00520d0,-0.00653d0/
      data i2l/-0.00110d0,-0.00022d0,-0.00003d0, 0.00031d0, 0.00053d0,
     &          0.00076d0, 0.00096d0, 0.00123d0, 0.00171d0/
      data i3l/ 0.00001d0, 0.00019d0, 0.00014d0, 0.00009d0, 0.00004d0,
     &          0.00000d0,-0.00004d0,-0.00010d0,-0.00024d0/
      data i4l/-0.00007d0,-0.00001d0, 0.00001d0, 0.00003d0, 0.00005d0,
     &          0.00007d0, 0.00009d0, 0.00012d0, 0.00017d0/

      data j1l/ 0.00294d0, 0.00658d0, 0.00745d0, 0.00778d0, 0.00810d0,
     &          0.00826d0, 0.00839d0, 0.00865d0, 0.00869d0/
      data j2l/ 0.00059d0,-0.00180d0,-0.00209d0,-0.00241d0,-0.00261d0,
     &         -0.00279d0,-0.00293d0,-0.00312d0,-0.00323d0/
      data j3l/-0.00018d0, 0.00036d0, 0.00044d0, 0.00054d0, 0.00059d0,
     &          0.00065d0, 0.00068d0, 0.00073d0, 0.00075d0/

      data kcl/-0.00337d0,-0.00398d0,-0.00447d0,-0.00444d0,-0.00451d0,
     &         -0.00448d0,-0.00445d0,-0.00448d0,-0.00439d0/

      data lcl/ 0.02116d0, 0.02499d0, 0.02811d0, 0.02793d0, 0.02840d0,
     &          0.02820d0, 0.02800d0, 0.02820d0, 0.02766d0/

      data p0/-0.00938d0,-0.00047d0,-0.00111d0, 0.00240d0, 0.00384d0,
     &         0.00636d0, 0.00848d0, 0.01025d0, 0.01464d0/
      data p1/ 0.00610d0, 0.00063d0, 0.00110d0,-0.00124d0,-0.00219d0,
     &        -0.00389d0,-0.00534d0,-0.00652d0,-0.00957d0/
      data p2/-0.00114d0,-0.00064d0,-0.00074d0,-0.00019d0, 0.00005d0,
     &         0.00050d0, 0.00091d0, 0.00123d0, 0.00222d0/
      data p3/-0.00010d0, 0.00030d0, 0.00024d0, 0.00022d0, 0.00017d0,
     &         0.00013d0, 0.00008d0, 0.00001d0,-0.00022d0/
      data p4/-0.00018d0,-0.00006d0,-0.00004d0, 0.00001d0, 0.00004d0,
     &         0.00008d0, 0.00011d0, 0.00016d0, 0.00025d0/

      data q1/ 0.00320d0, 0.01013d0, 0.01286d0, 0.01396d0, 0.01526d0,
     &         0.01587d0, 0.01632d0, 0.01790d0, 0.01867d0/
      data q2/ 0.00107d0,-0.00247d0,-0.00281d0,-0.00366d0,-0.00412d0,
     &        -0.00465d0,-0.00506d0,-0.00556d0,-0.00615d0/
      data q3/-0.00008d0, 0.00052d0, 0.00057d0, 0.00077d0, 0.00086d0,
     &         0.00100d0, 0.00111d0, 0.00119d0, 0.00133d0/

      data  r/-0.00442d0,-0.00650d0,-0.00861d0,-0.00866d0,-0.00933d0,
     &        -0.00931d0,-0.00928d0,-0.01015d0,-0.01023d0/

      data  s/ 0.02775d0, 0.04087d0, 0.05414d0, 0.05448d0, 0.05868d0,
     &         0.05857d0, 0.05837d0, 0.06388d0, 0.06442d0/


      data a0p/-0.01373d0, 0.02231d0, 0.01599d0, 0.01672d0, 0.01608d0,
     &          0.01767d0, 0.01927d0, 0.02137d0, 0.02584d0/
      data a1p/ 0.00957d0,-0.00589d0,-0.00191d0,-0.00325d0,-0.00329d0,
     &         -0.00455d0,-0.00567d0,-0.00661d0,-0.00894d0/
      data a2p/-0.00204d0,-0.00279d0,-0.00330d0,-0.00248d0,-0.00222d0,
     &         -0.00177d0,-0.00143d0,-0.00134d0,-0.00097d0/
      data a3p/-0.00005d0,-0.00073d0,-0.00075d0,-0.00078d0,-0.00078d0,
     &         -0.00084d0,-0.00091d0,-0.00103d0,-0.00125d0/
      data a4p/-0.00003d0,-0.00043d0,-0.00047d0,-0.00045d0,-0.00042d0,
     &         -0.00041d0,-0.00041d0,-0.00043d0,-0.00045d0/

      data b1p/ 0.00661d0,-0.00095d0, 0.00088d0,-0.00080d0,-0.00108d0,
     &         -0.00244d0,-0.00365d0,-0.00477d0,-0.00739d0/
      data b2p/ 0.00135d0,-0.00059d0,-0.00098d0,-0.00128d0,-0.00138d0,
     &         -0.00152d0,-0.00162d0,-0.00175d0,-0.00190d0/
      data b3p/-0.00035d0, 0.00002d0,-0.00036d0,-0.00037d0,-0.00044d0,
     &         -0.00044d0,-0.00045d0,-0.00054d0,-0.00062d0/

      data  cp/-0.00811d0,-0.00729d0,-0.00776d0,-0.00500d0,-0.00399d0,
     &         -0.00239d0,-0.00113d0,-0.00029d0, 0.00167d0/

      data  dp/ 0.05098d0, 0.06630d0, 0.08995d0, 0.08939d0, 0.09525d0,
     &          0.09466d0, 0.09407d0, 0.10004d0, 0.09950d0/

      data i0p/-0.00338d0, 0.00024d0,-0.00017d0, 0.00092d0, 0.00137d0,
     &          0.00213d0, 0.00277d0, 0.00343d0, 0.00480d0/
      data i1p/ 0.00231d0, 0.00018d0, 0.00055d0,-0.00017d0,-0.00044d0,
     &         -0.00094d0,-0.00137d0,-0.00178d0,-0.00271d0/
      data i2p/-0.00047d0,-0.00028d0,-0.00038d0,-0.00022d0,-0.00015d0,
     &         -0.00003d0, 0.00008d0, 0.00018d0, 0.00047d0/
      data i3p/-0.00003d0, 0.00012d0, 0.00011d0, 0.00011d0, 0.00011d0,
     &          0.00010d0, 0.00009d0, 0.00007d0, 0.00001d0/
      data i4p/ 0.00000d0,-0.00004d0,-0.00003d0,-0.00003d0,-0.00003d0,
     &         -0.00003d0,-0.00003d0,-0.00003d0,-0.00001d0/

      data j1p/ 0.00111d0, 0.00339d0, 0.00429d0, 0.00461d0, 0.00500d0,
     &          0.00520d0, 0.00535d0, 0.00579d0, 0.00604d0/
      data j2p/ 0.00042d0,-0.00082d0,-0.00088d0,-0.00115d0,-0.00129d0,
     &         -0.00146d0,-0.00159d0,-0.00176d0,-0.00197d0/
      data j3p/-0.00010d0, 0.00015d0, 0.00014d0, 0.00022d0, 0.00025d0,
     &          0.00030d0, 0.00034d0, 0.00038d0, 0.00044d0/

      data kcp/-0.00161d0,-0.00212d0,-0.00287d0,-0.00286d0,-0.00304d0,
     &         -0.00303d0,-0.00301d0,-0.00321d0,-0.00320d0/

      data lcp/ 0.01013d0, 0.01332d0, 0.01803d0, 0.01796d0, 0.01915d0,
     &          0.01906d0, 0.01896d0, 0.02018d0, 0.02012d0/

      data  alf0l/    1.6449d0,   0.6252d0,   0.4889d0,   0.4993d0,
     &                0.5104d0,   0.5502d0,   0.5962d0,   0.5793d0,
     &                0.6798d0/
      data  alf1l/  -23.2588d0,  10.6819d0,  16.1962d0,  15.8082d0,
     &               16.1851d0,  15.4934d0,  14.4670d0,  15.1152d0,
     &               12.7527d0/
      data  alf2l/  272.1670d0, -70.6879d0,-138.4860d0,-134.0100d0,
     &             -147.1310d0,-148.0760d0,-144.0690d0,-151.6510d0,
     &             -140.1800d0/
      data  alf3l/-1074.7000d0, -44.3349d0, 185.7060d0, 171.0490d0,
     &              230.5050d0, 250.3260d0, 251.8580d0, 276.8910d0,
     &              268.8290d0/

      data  bet0l/    1.6443d0,   0.6307d0,   0.5111d0,   0.5366d0,
     &                0.5855d0,   0.6380d0,   0.6814d0,   0.7331d0,
     &                0.7783d0/
      data  bet1l/  -23.2414d0,  10.4966d0,  15.4195d0,  15.4573d0,
     &               14.5626d0,  13.3510d0,  12.3005d0,  11.2258d0,
     &               10.2315d0/
      data  bet2l/  272.0080d0, -68.7973d0,-130.1540d0,-141.5810d0,
     &             -141.7690d0,-136.2650d0,-130.5920d0,-127.8900d0,
     &             -124.2640d0/
      data  bet3l/-1074.2500d0, -50.0581d0, 159.6050d0, 217.6680d0,
     &              237.9640d0, 235.7860d0, 229.2470d0, 238.5240d0,
     &              241.3060d0/

      data alf0p/   -0.1394d0,   0.5481d0,   0.3173d0,   0.3167d0,
     &               0.2524d0,   0.2448d0,   0.2617d0,   0.2427d0,
     &               0.2847d0/
      data alf1p/    7.0680d0, -20.4731d0, -14.4048d0, -14.2426d0,
     &             -12.4235d0, -12.1275d0, -12.4639d0, -12.2010d0,
     &             -13.0828d0/
      data alf2p/ -115.5940d0, 223.9220d0, 186.9100d0, 183.2260d0,
     &             172.5400d0, 173.0720d0, 177.6570d0, 180.4920d0,
     &             192.5030d0/
      data alf3p/  619.9170d0,-534.9400d0,-476.8100d0,-461.2490d0,
     &            -446.9220d0,-457.7020d0,-475.7000d0,-496.3580d0,
     &            -543.1480d0/

      data bet0p/   -0.1394d0,   0.5413d0,   0.3073d0,   0.2475d0,
     &               0.2521d0,   0.2809d0,   0.3088d0,   0.3205d0,
     &               0.3221d0/
      data bet1p/    7.0664d0, -20.2069d0, -13.7973d0, -12.0132d0,
     &             -12.0769d0, -12.7225d0, -13.3769d0, -13.7451d0,
     &             -13.8640d0/
      data bet2p/ -115.5800d0, 220.7060d0, 176.9940d0, 165.9970d0,
     &             170.6080d0, 177.8420d0, 184.6140d0, 194.3840d0,
     &             201.0700d0/
      data bet3p/  619.8790d0,-524.1240d0,-438.7520d0,-422.3750d0,
     &            -446.8310d0,-472.0400d0,-494.2780d0,-539.2810d0,
     &            -573.1160d0/

C     ---   Comienzo del calculo   ---

      gam3=1.0d0/(gamma**(1.0d0/3.0d0))

C     ---   Calculo de las variables de interpolacion   ---

      vl=alf0l(j)+gam3*(alf1l(j)+gam3*(alf2l(j)+gam3*alf3l(j)))
      wl=bet0l(j)+gam3*(bet1l(j)+gam3*(bet2l(j)+gam3*bet3l(j)))
      vp=alf0p(j)+gam3*(alf1p(j)+gam3*(alf2p(j)+gam3*alf3p(j)))
      wp=bet0p(j)+gam3*(bet1p(j)+gam3*(bet2p(j)+gam3*bet3p(j)))

C     ---   Parametros de la interpolacion   ---

      u1=u
      u2=2.0d0*u
      u3=3.0d0*u
      u4=4.0d0*u

      fu0171=(a0l(j)/2.0d0)+a1l(j)*dcos(u1)+a2l(j)*dcos(u2)+
     &                      a3l(j)*dcos(u3)+a4l(j)*dcos(u4)+
     &                      b1l(j)*dsin(u1)+b2l(j)*dsin(u2)+
     &                      b3l(j)*dsin(u3)+cl(j)*u+dl(j)

      fup171=(a0p(j)/2.0d0)+a1p(j)*dcos(u1)+a2p(j)*dcos(u2)+
     &                      a3p(j)*dcos(u3)+a4p(j)*dcos(u4)+
     &                      b1p(j)*dsin(u1)+b2p(j)*dsin(u2)+
     &                      b3p(j)*dsin(u3)+cp(j)*u+dp(j)

      gu0171=(i0l(j)/2.0d0)+i1l(j)*dcos(u1)+i2l(j)*dcos(u2)+
     &                      i3l(j)*dcos(u3)+i4l(j)*dcos(u4)+
     &                      j1l(j)*dsin(u1)+j2l(j)*dsin(u2)+
     &                      j3l(j)*dsin(u3)+kcl(j)*u+lcl(j)

      gup171=(i0p(j)/2.0d0)+i1p(j)*dcos(u1)+i2p(j)*dcos(u2)+
     &                      i3p(j)*dcos(u3)+i4p(j)*dcos(u4)+
     &                      j1p(j)*dsin(u1)+j2p(j)*dsin(u2)+
     &                      j3p(j)*dsin(u3)+kcp(j)*u+lcp(j)

      fu5000=(e0(j)/2.0d0)+e1(j)*dcos(u1)+e2(j)*dcos(u2)+
     &                     e3(j)*dcos(u3)+e4(j)*dcos(u4)+
     &                     f1(j)*dsin(u1)+f2(j)*dsin(u2)+
     &                     f3(j)*dsin(u3)+g(j)*u+h(j)

      gu5000=(p0(j)/2.0d0)+p1(j)*dcos(u1)+p2(j)*dcos(u2)+
     &                     p3(j)*dcos(u3)+p4(j)*dcos(u4)+
     &                     q1(j)*dsin(u1)+q2(j)*dsin(u2)+
     &                     q3(j)*dsin(u3)+r(j)*u+s(j)

C     ---   Contribucion de la red (lattice)   ---

      fl=(1.0d0-vl)*fu0171+vl*fu5000
      gl=(1.0d0-wl)*gu0171+wl*gu5000

C     ---   Contribucion de los fonones (phonon)   ---

      fp=vp*fup171
      gp=wp*gup171

C     ---   Contribucion total   ---

      fsold=fl+fp
      gsold=gl+gp

      return
      end


C***********************************************************************
C***********************************************************************


      subroutine nbrpd(rho,t,fbrpd,gbrpd)
C=======================================================================
C
C       Esta subrutina calcula la contribucion al ritmo de perdida de 
C       energia por neutrinos de un gas de electrones parcialmente de-
C       generado en regimen de neutrino-pair bremsstrahlung. 
C
C       Refs.: Munakata, K., Kohyama, Y., Itoh, N., Ap. J., 316, 708 
C              (1987). 
C
C       Revisada el 23.09.91 por Abulafia
C
C-----------------------------------------------------------------------
C
C       Parametros de entrada:
C
C       rho: densidad.
C       t:   temperatura
C
C-----------------------------------------------------------------------
C
C       Parametros de salida:
C       fbrpd: coeficiente f para el bremsstrahlung.
C       gbrpd: coeficiente g para el bremsstrahlung.
C
C=======================================================================
      implicit double precision (a-h,o-z)

      data a0/ 2.35d+1/, a1/6.83d+4/, a2/7.81d+8/
      data a3/ 2.30d+2/, a4/6.70d+5/, a5/7.66d+9/
      data b1/ 1.47d+0/, b2/3.29d-2/

C     ---   Calculo de los parametros necesarios   ---

      t8=1.0d-8*t

      b3=7.75d5*(t8**1.50d0)+2.47d2*(t8**3.85d0)
      b4=4.07d0+2.40d-2*(t8**1.40d0)
      b5=4.59d-5/(t8**0.11d0)

      exim1=((7.05d6+5.12d4*(t8**1.5d0))*(t8**1.5d0))/rho

C     ---   Calculo de los coeficientes del bremsstrahlung   ---

      fbrpd=1.0d0/(a0+a1*(t8**(-2.0d0))+a2*(t8**(-5.0d0)))+
     &      1.26d0*(1.0d0+exim1)/(1.0d0+exim1*(b1+b2*exim1))

      gbrpd=1.0d0/((1.0d0+1.0d-9*rho)*(a3+a4*(t8**(-2.0d0))+
     &      a5*(t8**(-5.0d0))))+1.0d0/((b3/rho)+b4+b5*
     &      (rho**0.656d0))

      return
      end


C***********************************************************************
C***********************************************************************


      subroutine neutri(rho,t,ye,z,enu)
C=======================================================================
C
C     Par, plasma y fotoneutrinos a la Beaudet, Petrosian y Salpeter, 
C     brehmsstrahlung a la De Zotti.
C
C-----------------------------------------------------------------------
C
C     Revisada el 24.01.95 por E. Garcia-Berro.
C
C-----------------------------------------------------------------------
C
C       Parametros de entrada:
C
C       rho:    densidad 
C       ye:     numero electronico molar
C       t:      temperatura
C       z:      carga atomica
C
C-----------------------------------------------------------------------
C
C       Parametros de salida:
C
C       enu:    energia emitida por neutrinos (erg/g/s)
C
C=======================================================================
      implicit double precision (a-h,o-z)

C     ---   Dimensionado de variables   ---

      dimension a0(3),a1(3),a2(3),b1(3),b2(3),b3(3),c(3),f(3)

C     ---   Datos   ---

      data a0 / 6.002d+19, 2.320d-07, 4.886d+10/
      data a1 / 2.084d+20, 8.449d-08, 7.580d+10/
      data a2 / 1.872d+21, 1.789d-08, 6.023d+10/
      data b1 / 9.383d-01, 2.581d-02, 6.290d-03/
      data b2 /-4.141d-01, 1.734d-02, 7.483d-03/
      data b3 / 5.829d-02, 6.990d-04, 3.061d-04/

      data c /5.5924d0, 0.5646d0, 1.5654d0/

      data abr0 /0.5167d0/, bbr0 /7.577d0/, bbr1 /-2.801d0/, 
     &     bbr2 /0.267d0/, cbr0 /0.4775d0/, cbr1 /1.707d0/, 
     &     cbr2 /-9.090d0/, cbr3 /5.250d0/

C     ---   Parametros del calculo   ---

      rhoam=rho*ye
      bl=t/5.9302d09
      xi=1.0d-03*rhoam**(1.0d0/3.0d0)/bl
      gl=1.0d0-13.04d0*bl**2+133.5d0*bl**4+1534.0d0*bl**6+918.60d0*bl**8

C     ---   Llamada a la funcion f(x,l)   ---

         do 1 i=1,3
         call funf(a0(i),a1(i),a2(i),b1(i),b2(i),b3(i),c(i),bl,xi,f(i))
    1    continue

C     ---   Par, plasma y fotoneutrinos   ---

      fpar=f(1)
      fpls=f(2)
      ffot=f(3)

      vblamb=2.0d0/bl

C     ---   Control de underflows   ---

         if (fpar .lt. 1.0d-60) then
         fpar=0.0d0
         else
         continue
         end if

         if (fpls .lt. 1.0d-60) then
         fpls=0.0d0
         else
         continue
         end if

         if (ffot .lt. 1.0d-60) then
         ffot=0.0d0
         else
         continue
         end if

C     ---   Par neutrinos solo a temperaturas altas   ---

      qq=(rhoam**3.0d0)*fpls+rhoam*(bl**5.0d0)*ffot

         if (vblamb .gt. 7.0d2) then
         continue
         else
         qq=qq+gl*dexp(-vblamb)*fpar
         end if

      enu1=qq/rho

C     ---   Neutrino brehmsstrahlung   ---

      xbr=dlog10(rhoam)
      ybr=dlog10(t)-8.0d0

      gbr=cbr0*datan(cbr1*(xbr-ybr)+cbr2)+cbr3
      fbr=abr0*xbr-dsqrt(bbr0+bbr1*xbr+bbr2*(xbr**2.0d0))+ybr*gbr
      enu2=7.87d-3*(2.0d0*ye)*z*(1.0d1**(fbr))

C     ---   Totales   ---

      enu=enu1+enu2

      return
      end


C***********************************************************************
C***********************************************************************


      subroutine funf(a0,a1,a2,b1,b2,b3,c,bl,xi,flxi)
C=======================================================================
C
C     Funcion auxiliar para el calculo de los neutrinos.
C
C-----------------------------------------------------------------------
C
C     Revisada el 24.01.95 por E. Garcia-Berro.
C
C-----------------------------------------------------------------------
C
C       Parametros de entrada:
C
C       a's:    Coeficientes a_i de la funcion
C       b's:    Coeficientes b_i de la funcion
C       c:      Coeficiente c de la funcion
C       bl:     Lambda
C       xi:     Exi
C
C-----------------------------------------------------------------------
C
C       Parametros de salida:
C
C       flxi:   Valor de la funcion.
C
C=======================================================================
      implicit double precision (a-h,o-z)

      cxi=c*xi

         if (cxi .gt. 7.0d2) then
         flxi=0.0d0
         else
         flxi=(a0+a1*xi+a2*(xi**2.0d0))*dexp(-cxi)/(xi**3.0d0+(b1/bl)+
     &        (b2/(bl**2.0d0))+(b3/(bl**3.0d0)))
         end if

      return
      end 


C***********************************************************************
C***********************************************************************


      subroutine eweak(rho,t,ye,z2oa,eneut,denudt,denudr)
C=======================================================================
C
C     Neutrino losses a la Beaudet, Petrosian and Salpeter. Ap. J., 150,
C     979 (1967), brehmsstrahlung a la Festa & Ruderman.
C
C-----------------------------------------------------------------------
C
C     Input parameters:
C
C     rho:     Density.
C     t:       Temperature.
C     ye:      Electron molecular weight.
C     z2oa:    <Z^2/A>
C
C-----------------------------------------------------------------------
C
C     Output parameters:
C
C     eneut:   Neutrino losses.
C     denudt:  Temperature derivative.
C     denudr:  Density derivative.
C
C=======================================================================
      implicit double precision (a-h,o-z)

C     ---   Necessary parameters   ---

      t6=t/1.0d6
      rhmu=rho*ye
      elsal=t6/5930.2d0
      xsal=(5.9302d0/t6)*((rhmu)**(1.0d0/3.0d0))

C     ---   ...and their powers   ---

      xsal2=xsal*xsal
      xsal3=xsal2*xsal
      elsal2=elsal*elsal
      elsal3=elsal2*elsal

C     ---   Compute g and its derivative   ---

      gsal  =elsal2*(-13.04d0+elsal2*(133.5d0+elsal2*(1534.0d0+918.6d0*
     &       elsal2)))
      gsal=gsal+1.0d0
      dlgdlt=elsal2*(-26.08d0+elsal2*(534.0d0+elsal2*(9204.0d0+7348.8d0
     &       *elsal2)))
      dlgdlt=dlgdlt/gsal

C     ---  Compute the denominator of interpolation formula   ---

      b1sal=2.581d-2/elsal
      b2sal=b1sal+1.734d-2/elsal2
      btsal=xsal3+b2sal+6.99d-4/elsal3

C     ---   Compute the rest of interpolation formula   ---

      a1sal=8.449d-8*xsal
      a2sal=1.787d-8*xsal2
      atsal=2.32d-7+a1sal+a2sal
      ahsal=(a1sal+a2sal+a2sal)/atsal

C     ---   Compute plasma neutrinos   ---

      csal=0.56457d0*xsal
      qplsal=rhmu*rhmu*rhmu*atsal*dexp(-csal)/btsal

C     ---    ...and its derivatives   ---

      dfdrpl=(ahsal-csal)/3.0d0-xsal3/btsal
      dfdtpl=csal+3.0d0-ahsal-(b1sal+b2sal)/btsal

C     ---  Compute the denominator of interpolation formula   ---

      b1sal=6.29d-3/elsal
      b2sal=b1sal+7.483d-3/elsal2
      btsal=xsal3+b2sal+3.061d-4/elsal3

C     ---   Compute the rest of interpolation formula   ---

      a1sal=7.580d+10*xsal
      a2sal=6.023d+10*xsal2
      atsal=4.886d+10+a1sal+a2sal
      ahsal=(a1sal+a2sal+a2sal)/atsal

C     ---   Compute photo neutrinos   ---

      csal=1.5654d0*xsal
      qphsal=rhmu*elsal2*elsal3*atsal*dexp(-csal)/btsal

C     ---    ...and its derivatives   ---

      dfdrph=(ahsal-csal)/3.0d0-xsal3/btsal
      dfdtph=csal+3.0d0-ahsal-(b1sal+b2sal)/btsal

C     ---  Compute the denominator of interpolation formula   ---

      b1sal=0.9383d0/elsal
      b2sal=b1sal-0.4141d0/elsal2
      btsal=xsal3+b2sal+5.829d-2/elsal3

C     ---   Compute the rest of interpolation formula   ---

      a1sal=2.084d+20*xsal
      a2sal=1.872d+21*xsal2
      atsal=6.002d+19+a1sal+a2sal
      ahsal=(a1sal+a2sal+a2sal)/atsal

C     ---   Compute pair neutrinos   ---

      csal=5.5924d0*xsal
      qpasal=gsal*dexp(-csal-2.0d0/elsal)*atsal/btsal

C     ---    ...and its derivatives   ---
      dfdrpa=(ahsal-csal)/3.0d0-xsal3/btsal
      dfdtpa=csal+3.0d0-ahsal-(b1sal+b2sal)/btsal

C     ---   Combine the three contributions   ---

      eneut=(qplsal+qphsal+qpasal)/rho
      denudr=(qplsal*(2.0d0+dfdrpl)+qphsal*dfdrph+qpasal*
     1       (dfdrpa-1.0d0))/(rho)
      denudt=(qplsal*dfdtpl+qphsal*(5.0d0+dfdtph)+qpasal*
     1       (dlgdlt+2.0d0/elsal+dfdtpa))/(rho)

C***********************************************************************
C       ***   Festa and Ruderman, Phys. Rev., 180, 1227, 1969   ***
C    ***   Neutrino-pair bremstrahlung (relativistic degeneracy)   ***
C***********************************************************************

C     ---   Necessary parameters   ---

      tt2=t6*t6
      tt6=tt2*tt2*tt2
      deblr=25.3286d0-dlog(rho)

C     ---   Neutrino bremsstrahlung ala Festa & Ruderman   ---

      eneufr=7.6d-13*z2oa*tt6

C     ---   Correction factor   ---

         if (deblr .gt. 0.0d0) then
         afac=dexp(-deblr*deblr/224.0d0)
         else
         afac=1.0d0
         end if

      eneufr=eneufr*afac

C      ---   Derivatives   ---

         if (deblr .gt. 0.0d0) then
         afac=deblr/112.0d0
         else
         afac=0.0d0
         end if

      defrt=6.0d0*eneufr
      defrr=afac*eneufr

C***********************************************************************
C       ***   Remember Gandelman and Pinaev, JETP, 10, 764, 1960   ***
C ***   Non relativistic degeneracy  eneut=0.85e-5*rho*z2oa*(t8**4.5)***
C***********************************************************************

C     ---   Neutrino and nuclear contributions   ---

      eneut=eneut+eneufr
      denudt=defrt+denudt
      denudr=defrr+denudr

      return
      end


C***********************************************************************
C***********************************************************************



