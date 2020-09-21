         subroutine gw_radiation(ipas)

!--Load modules

         use mod_essentials
         use mod_parameters, ONLY : maxpas
         use mod_commons,    ONLY : xyzhm, vxyzut, axyzut

!--Force to declare EVERYTHING

         implicit none 

!--I/O variables

         integer, intent(in) :: ipas

!--Local variables

         real(8), dimension(ndim,ndim,ndim,ndim) :: pijkl
         real(8), dimension(ndim,ndim,maxpas)    :: ddi
         real(8), dimension(ndim,ndim)           :: deltaij, dddi
         real(8), dimension(ndim)                :: ni
         real(8)  :: tpas,trdddi,hp00,hx00,hpp0,hxp0, degdt,dddi2,time1
         integer  :: i,j,k,l,k2,p

!--Physical parameters in code units 
!  
!  c = 3x10^10 cm/s = 215.49
!  G = 6.67x10^{-8} cm^3 g^-1 s^-2 = 1.0
!  distobs = 10kpc = 3.08x10E+22 cm

         real(8), parameter :: Gdivc4=4.64d-10, distobs=4.43d+12,       &
                               Gdiv5c5=4.30d-13

!--cgs-code conversion factors

         real(8), parameter :: potfact=7.7d+47, velfact=1.4d+8,         &
                               lumisol=3.9d+33

!--Open I/O files

         open(unit=1,file='degdt',status'unknown', access='append')
         open(unit=2,file='h+00', status='unknown',access='append')
         open(UNIT=3,file='hx00', status='unknown',access='append')
         open(UNIT=4,file='h+p0', status='unknown',access='append')
         open(UNIT=5,file='hxp0', status='unknown',access='append')

         if (ipas == 1) then
           write(1,*) 'seconds     degdt(Lsol)'
           write(2,*) 'seconds     h+00'
           write(3,*) 'seconds     hx00'
           write(4,*) 'seconds     h+p0'
           write(5,*) 'seconds     hxp0'
         endif

!--Start the GW calculation.Construct Kronecker's delta

         do i=1,3
           do j=1,3
             deltaij(i,j) = 0.0d0
           enddo
         enddo

         do i=1,3
            deltaij(i,i) = 1.0d0
         enddo

!--Construct PIJKL and NI tensors

         do i=1,3
            ni(i) = deltaij(3,i)
         enddo
         do i=1,3
            do j=1,3
               do k=1,3
                  do l=1,3
                     pijkl(i,j,k,l) = (deltaij(i,k)-ni(i)*ni(k))*       &
                       (deltaij(j,l)-ni(j)*ni(l))-                      &
                       0.5d0*(deltaij(i,j)-ni(i)*ni(j))*                &
                       (deltaij(k,l)-ni(k)*ni(l))         
                  enddo
               enddo
            enddo
        enddo

        time1   = tnow*50.4d0 !=> time1 is the physical time. In my old code each code unit of time was 50.4 seconds.

!--Variables startup
    
        do k=1,3
           do k2=1,3
              ddi(k,k2,ipas) = 0.0d0
              dddi(k,k2)     = 0.0d0
           enddo
        enddo
        trdddi = 0.0d0

!--DDI calculation

        do k=1,3
           do l=1,3
              do p=1,nbody
                 ddi(k,l,ipas) = ddi(k,l,ipas) + xyzhm(5,p)*(2.0d0*     &
                 vxyzut(k,p)*vxyzut(l,p) + xyzhm(k,p)*axyzut(l,p) +     &
                 xyzhm(l,p)*axyzut(k,p))
              enddo
           enddo
        enddo

!--H+ & Hx calculation

        hp00 = (Gdivc4)*(ddi(1,1,ipas)-ddi(2,2,ipas))/distobs
        hx00 = (Gdivc4)*(2.0d0*ddi(1,2,ipas))/distobs
        hpp0 = (Gdivc4)*(ddi(3,3,ipas)-ddi(2,2,ipas))/distobs
        hxp0 = (Gdivc4)*(-2.0d0*ddi(2,3,ipas))/distobs

!--Construct DI derivatives and the radiated energy

        if (ipas > 1) then
           do i=1,3
              do j=1,3
                 dddi(i,j) = (ddi(i,j,ipas)-ddi(i,j,ipas-1))/dt !=> dt is the integration time-step
              enddo
           enddo

           do i=1,3
              trdddi = trdddi + dddi(i,i)
           enddo

           do i=1,3
              do j=1,3
                 dddi(i,j) = dddi(i,j)-(0.3333d0)*deltaij(i,j)*trdddi
              enddo
           enddo

           dddi2=0.0d0
           do i=1,3
              do j=1,3 
                 dddi2 = dddi2+dddi(i,j)*dddi(i,j)
              enddo
           enddo
           degdt = Gdiv5c5*dddi2

!--Code units-cgs conversion

           degdt = potfact*degdt
           do i=1,3
              do j=1,3
                 dddi(i,j) = velfact*dddi(i,j)
              enddo
           enddo    

!--end of GW calculation. write results and end program.

           write(1,'(2(1pe12.4))') time1,degdt/lumisol

        endif

        write(2,'(2(1pe12.4))') time1, hp00
        write(3,'(2(1pe12.4))') time1, hx00
        write(4,'(2(1pe12.4))') time1, hpp0
        write(5,'(2(1pe12.4))') time1, hxp0

        close(1)
        close(2)
        close(3)
        close(4)
        close(5)

        call flush (1)
        call flush (2)
        call flush (3)
        call flush (4)
        call flush (5)
!
        end subroutine gw_radiation
