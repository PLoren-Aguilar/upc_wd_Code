      SUBROUTINE helio(relflag,simtype,nbody1,nbody,nel,posx,posy,posz,mass,xss,rank)

      implicit none

      logical, intent(in) :: relflag
      integer, intent(in) :: nbody1,nbody,nel,simtype,rank
      real(8), intent(in) :: mass(nbody),posx(nbody),posy(nbody),posz(nbody)
      real(8), intent(in out) :: xss(nel,nbody)

      integer :: p,i,j
      integer :: ntheta,nphi,nher,cont,nhe
      real(8), parameter :: masshe1=0.01d0,masshe2=0.01d0                 ! masa de helio de la wd
      real(8) :: massparhe,massparco,mass1,mass2,mass_comp_he

! ajustar esta variable, segun valor deseado
      nhe=3000

      IF((simtype.eq.1).and.(relflag.eqv..true.)) THEN

        ! numero de particulas de helio
        mass1=sum(mass(1:nbody1))
        massparhe=masshe1/nhe
        massparco=(mass1-masshe1)/(nbody1-nhe)

        ntheta=12
        nphi=ntheta*2
        nher=nhe/nphi/ntheta
        if(ntheta*nphi*nher.lt.nhe) then
          nher=nher+1                                 
          nhe=nher*nphi*ntheta
          massparhe=masshe1/nhe
          massparco=(mass1-masshe1)/(nbody1-nhe)
!        if(ntheta*nphi*nher.ne.nhe) then
!          nphi=nphi+(nhe-ntheta*nphi*nher)/ntheta+1
        end if

        if(rank.eq.0) print*,'num particulas de He:',ntheta*nphi*nher,ntheta,nphi,nher
        if(rank.eq.0) print*,'masa particulas He:',massparhe,' C=:',massparco

        call pongo_he(nbody,nel,ntheta,nphi,nher,1,nbody1,mass,posx,posy,posz,xss,massparhe,massparco)

        ! cuento particulas He puestas
        mass_comp_he=0.0d0
        cont=0
        do p=1,nbody
           if(xss(2,p).gt.0.5d0) then
             mass_comp_he=mass_comp_he+mass(p)
             cont=cont+1
             if(mass(p).ne.massparhe) then
               if(rank.eq.0) print*,'MAL HE!!'
             end if
           else
             if(mass(p).ne.massparco) then
               if (rank.eq.0) print*,'MAL CO!!'
             end if
           end if
        end do
        if(rank.eq.0) print*,'num particulas de He puestas:',ntheta*nphi*nher
        if(rank.eq.0) print*,'comprobacion massa total: M = ',sum(mass(:)),'MHe = ',mass_comp_he

      ELSE IF ((simtype.eq.2.).and.(relflag.eqv..true.)) then

        ! numero de particulas de helio
        mass1=sum(mass(1:nbody1))
        massparhe=masshe1/nhe
        massparco=(mass1-masshe1)/(nbody1-nhe)

        ntheta=12
        nphi=ntheta*2
        nher=nhe/nphi/ntheta
        if(ntheta*nphi*nher.ne.nhe) then
          nphi=nphi+(nhe-ntheta*nphi*nher)/ntheta+1
        end if
        call pongo_he(nbody,nel,ntheta,nphi,nher,1,nbody1,mass,posx,posy,posz,xss,massparhe,massparco)

        ! numero de particulas de helio
        mass2=sum(mass(nbody1+1:nbody))
        massparhe=masshe2/nhe
        massparco=(mass2-masshe2)/((nbody-nbody1)-nhe)

        ntheta=12
        nphi=ntheta*2
        nher=nhe/nphi/ntheta
        if(ntheta*nphi*nher.ne.nhe) then
          nphi=nphi+(nhe-ntheta*nphi*nher)/ntheta+1
        end if
        call pongo_he(nbody,nel,ntheta,nphi,nher,nbody1+1,nbody,mass,posx,posy,posz,xss,massparhe,massparco)

      END IF

      END SUBROUTINE


      SUBROUTINE ordeno_rad(nbody,length,particulas,posx,posy,posz,cmp,orden)

      implicit none
      
      integer, intent(in) :: length,particulas(length+1),nbody
      real(8), intent(in) :: posx(nbody),posy(nbody),posz(nbody),cmp(3)
      integer, intent(out) :: orden(length)

      real(8) :: radio(length),radioord(length)
      integer :: p,j,i

      orden=0
      radioord=0.0d0

      ! compruebo datos estan ordenados
      if(particulas(length+1).ne.0) print*, 'ordeno_rad: malos datos'

      ! ORDENO PARTICULAS POR RADIO
      do i=1,length
         p=particulas(i)
         radio(i)=dsqrt((posx(p)-cmp(1))**2+(posy(p)-cmp(2))**2+(posz(p)-cmp(3))**2)
         if(i.eq.1) then
           orden(i)=p
           radioord(i)=radio(i)
         else if (i.eq.2) then
           if(radio(i).gt.radio(i-1)) then
              radioord(i)=radio(i)
              orden(i)=p
           else
              radioord(2)=radioord(1)
              orden(2)=orden(1)
              radioord(1)=radio(i)
              orden(1)=p
           end if
         else ! las p-1 posiciones anteriores estan ordenadas
           if(radioord(i-1).le.radio(i)) then
             radioord(i)=radio(i)
             orden(i)=p
           else 
             do j=i-2,1,-1
                if(radioord(j).lt.radio(i)) then
                  ! muevo array radioord desde posicion j+1 hasta al final (posicion p) una posicion a la derecha y coloco radio(p) en posicion j
                  call muevo(nbody,radioord,orden,j+1,i)
                  orden(j+1)=p
                  radioord(j+1)=radio(i)
                  go to 654
                end if
             end do
             ! radio(p) es el mas pequeño de todos, lo pongo al principio
             call muevo(nbody,radioord,orden,1,i)
             orden(1)=p
             radioord(1)=radio(i)
           end if
         end if
 654 continue
      end do
     
      END SUBROUTINE


      SUBROUTINE muevo(nbody,radioord,ordeno,j,p)

      integer, intent(in) :: nbody,j,p
      integer, intent(in out) :: ordeno(nbody)
      real(8), intent(in out) :: radioord(nbody)

      do i=p,j+1,-1
         radioord(i)=radioord(i-1)
         ordeno(i)=ordeno(i-1)
      end do

      END SUBROUTINE


      SUBROUTINE pongo_he(nbody,nel,ntheta,nphi,nher,nbody0,nbody1,mass,posx,posy,posz,xss,massparhe,massparco)
        
      implicit none

      real(8), parameter :: pi=3.141592654d0
      integer, intent(in) :: nbody0,nbody1,nbody,nel,ntheta,nphi,nher
      real(8), intent(in) :: posx(nbody),posy(nbody),posz(nbody),massparhe,massparco
      real(8), intent(in out) :: xss(nel,nbody),mass(nbody)

      integer :: contpar(nphi*ntheta),cuadrante(nphi*ntheta,nbody/2)
      integer :: inphi,intheta,inmatrix
      integer :: p,i,j,cont,orden(nbody/2)
      real(8) :: masstot,cmp(3),masshe,radio,theta,phi
 

      masstot=sum(mass(nbody0:nbody1))
      cmp(1)=sum(posx(nbody0:nbody1)*mass(nbody0:nbody1))/masstot
      cmp(2)=sum(posy(nbody0:nbody1)*mass(nbody0:nbody1))/masstot
      cmp(3)=sum(posz(nbody0:nbody1)*mass(nbody0:nbody1))/masstot

      ! asigno queso a cada particula -> construyo martiz cuadrante
      ! la matriz cuadrante contiene en cada fila las particulas que pertenecen a ese cuadrante
      ! cuadrante esta ordenada por theta  --> phi

      contpar=0
      cuadrante=0
      do p=nbody0,nbody1
             
         radio=dsqrt((posx(p)-cmp(1))**2+(posy(p)-cmp(2))**2+(posz(p)-cmp(3))**2)
         phi=datan((posy(p)-cmp(2))/(posx(p)-cmp(1)))
         if((posx(p)-cmp(1)).lt.0.0d0) then ! 2º o 3º cuadrante
           phi=phi+pi
         else if ((posy(p)-cmp(2)).lt.0.0d0) then ! 4º cuadrante
           phi=2.0d0*pi+phi
         end if
         theta=datan(dsqrt((posx(p)-cmp(1))**2+(posy(p)-cmp(2))**2)/(posz(p)-cmp(3)))
         if((posz(p)-cmp(3)).lt.0.0d0) then
           theta=pi+theta
         end if

         inphi=phi*nphi/2.0d0/pi+1
         intheta=theta*ntheta/pi+1

         if((inphi.gt.nphi).or.(inphi.lt.1)) print*,'MAL phi',inphi
         if((intheta.gt.ntheta).or.(intheta.lt.1)) print*,'MAL theta',intheta

         inmatrix=(intheta-1)*nphi+inphi ! numero de queso al que pertenece
         if((inmatrix.lt.1).or.(inmatrix.gt.ntheta*nphi)) print*,'MAL inmatrix'
         contpar(inmatrix)=contpar(inmatrix)+1
         cuadrante(inmatrix,contpar(inmatrix))=p 

      end do   
 
      ! recorro cuadrantes: en cada cuadrante ordeno particulas por radio
      do i=1,ntheta*nphi

         ! ordeno particulas por radio de este cuadrante (de menor a mayor)
         call ordeno_rad(nbody,contpar(i),cuadrante(i,1:contpar(i)+1),posx,posy,posz,cmp,orden)

         ! pongo helio en este cuadrante
         cont=0
         do j=contpar(i),1,-1  ! recorro de radio mayor a menor en este cuadrante
            cont=cont+1
            p=orden(j)
            if((p.lt.nbody0).or.(p.gt.nbody1)) print*,'particula mal contada'
            if(cont.le.nher) then  ! pongo helio
              xss(:,p)=0.0d0
              xss(2,p)=1.0d0
              mass(p)=massparhe
            else ! no toca helio
              xss(:,p)=0.0d0
              xss(3,p)=0.5d0
              xss(4,p)=0.5d0
              mass(p)=massparco
            end if
         end do 

         orden=0

      end do

      END SUBROUTINE


      ! calculo masa acumulada
      
!      massacum=0.0d0
!      do i=1,nbody
!         p=orden(i)
!         pp=p+nbody0-1
!         if(p.gt.nbody) then
!           print*,'problemon',p
!           stop
!         end if
!         massacum=massacum+mass(pp)
!         if(massacum.ge.(masstot-masshe)) then
!            xss(pp,:)=0.0d0
!            xss(pp,2)=1.0d0
!         else
!            xss(pp,:)=0.0d0
!            xss(pp,3)=0.5d0
!            xss(pp,4)=0.5d0
!         end if
!      end do




!          intheta=i/nphi+1
!          inphi=mod(i,nphi)

!          if(inphi.eq.0) inphi=nphi
  
!          phi1=2.0d0*pi/nphi*(inphi-1)
!          phi2=2.0d0*pi/nphi*inphi
!          theta1=pi/ntheta*(intheta-1)
!          theta2=pi/ntheta*intheta
 


