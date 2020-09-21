      subroutine dissolve
!============================================================
!  This subroutine disolves the CO particles inside the He
!  layer
!
!  Last revision: 06/April/2019
!============================================================
!
!--Load modules
!
      use mod_parameters, only : ndim
      use mod_commons,    only : xyzhm, rho, partype, nb, nbody, rank,  &
                                 factor, npart, ierr, vxyzut
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
      integer, dimension(npart) :: isentarray, ireceivearray
      real, dimension(2*npart)  :: sentarray, receivearray
      integer, dimension(npart) :: localpartype
      real, dimension(npart)    :: localtemp, localmass
      real    :: massHe, tempHe, r2, u2p, hp, dr, drx, dry, drz, wkp
      integer :: i, q, qlocal, p, k, m, nCO, nHe
!
!--Figure out the mass of the He(layer) particles
!
      massHe = 0.0
      massloop : do p=1,nbody
         if (partype(p) == 1) then
           massHe = xyzhm(5,p)
           exit massloop
         endif
      enddo massloop
!
!--If the CO particle has more He(layer) neighbours than CO neighbours,
!  dissolve the CO particle. 
!
!$OMP parallel default(none) shared(nbody,partype,localpartype,localtemp,localmass,nb,rank,factor)  &
!$OMP shared(xyzhm,vxyzut,massHe) private(p,q,qlocal,i,nCO,nHe,tempHe)
!$OMP do schedule(runtime)
      partloop : do qlocal = 1, factor
        q = qlocal + rank*factor
!
        if ((q >= nbody) .or. (partype(q) /= 0)) cycle
!
        nCO = 0
        nHe = 0
        tempHe = 0.0
        neiloop : do i = 1,nb(1,qlocal)
          p = nb(i+1,qlocal)
!
          if (partype(p) == 0) nCO = nCO + 1
          if (partype(p) == 1) then
             nHe = nHe + 1
             tempHe = tempHe + xyzhm(5,p)*vxyzut(5,p)
          endif
!
        enddo neiloop
        if (nHe > nCO) then
          localpartype(qlocal) = 1
          localtemp(qlocal)    = tempHe/(nHe*massHe)
          localmass(qlocal)    = massHe
        ELSE
          localpartype(qlocal) = partype(q)
          localtemp(qlocal)    = vxyzut(5,q)
          localmass(qlocal)    = xyzhm(5,q)
        endif
      enddo partloop
!$OMP end do
!$OMP end parallel
!
!--Transfer of particle types MPI processes, if necessary
! 
#ifdef MPI
!$OMP parallel default(none) shared(sentarray,localpartype,localtemp,localmass,factor)        &
!$OMP shared(isentarray) private(p)
!$OMP do schedule(runtime)
      do p = 1, factor
        isentarray(p) = localpartype(p)
        sentarray(2*(p-1) + 1) = localtemp(p)
        sentarray(2*(p-1) + 2) = localmass(p)
      enddo
!$OMP end do
!$OMP end parallel
!
      call MPI_ALLGATHER(isentarray,factor,MPI_INTEGER,ireceivearray,   &
                         factor,MPI_INTEGER,MPI_COMM_WorLD,ierr)
      call MPI_ALLGATHER(sentarray,2*factor,MPI_DOUBLE_PRECISION,       &
                         receivearray,2*factor,MPI_DOUBLE_PRECISION,    &
                         MPI_COMM_WORLD,ierr)
!
!$OMP parallel default(none) shared(nbody,ireceivearray,receivearray,partype,vxyzut,xyzhm)  &
!$OMP private(p)
!$OMP do schedule(runtime)
      do p = 1, nbody
        if (partype(p) == 0) then
          partype(p)  = ireceivearray(p)
          vxyzut(5,p) = receivearray(2*(p-1)+1)
          xyzhm(5,p)  = receivearray(2*(p-1)+2)
        endif
      enddo
!$OMP end do
!$OMP end parallel
#else
!$OMP parallel default(none) shared(nbody,localpartype,partype,localpartype,massHe,xyzhm)  &
!$OMP shared(vxyzut) private(p)
!$OMP do schedule(runtime)
      do p = 1, nbody
        if (partype(p) == 0) then
          xyzhm(5,p)  = massHe
          vxyzut(5,p) = localtemp(p)
          partype(p)  = localpartype(p)
        endif
      enddo
!$OMP end do
!$OMP end parallel
#endif MPI

      end subroutine dissolve
