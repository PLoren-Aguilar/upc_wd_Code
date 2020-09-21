      SUBROUTINE treeconst
!===========================================================================
!     This subroutine builds the tree and calculates neighbour lists
!
!     Last revision: 15/March/2015
!===========================================================================
!
!--Load modules
!
      USE mod_parameters, ONLY : ndim, MASTER, hmax, hmin
      USE mod_commons,    ONLY : xyzhm, gxyzu, eps, eps3, nb, globnmin, &
      globnmax, globnvec, partype, nbody, npart, factor, rank, nstep,   &
      nout, nbmax, nbmin, mmax, ierr, size, istep, step, istep0, nstep
                                  
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
#ifdef MPI
      INCLUDE 'mpif.h'
#endif MPI
!
!--Pointers
!
      REAl, DIMENSION(:), POINTER :: rx, ry, rz, qxx, qyy, qzz, qxy,    &
                                        qyz, qzx
      INTEGER, DIMENSION(:), POINTER :: list, nay
!
!--Local variables
!
      REAL, DIMENSION(6*mmax), TARGET :: quad
      REAL, DIMENSION(3*mmax), TARGET :: r
      REAL, DIMENSION(mmax)  :: em, qrad, hnode
      REAL, DIMENSION(nbody) :: xx, yy, zz, emm
#ifdef MPI
      REAL, DIMENSION((ndim+1)*npart) :: sentarray, receivearray
#endif MPI
      REAL, DIMENSION(100) :: xmap, ymap, zmap
      REAL :: t1, t2, t3, rr, rr05, rr2, rrx, rry, rrz, hm, dxi, dyi,   &
              dzi, rcomp2, qradn2, hsum2, hn, eps2, frx, fry, frz,      &
              pot, frxl, fryl, frzl, potl, difx, dify, difz,            &
              accpar, fl, fll, emred, ddd, xnp, ynp, znp,d, dmin, hnew, &
              hup, hdown
      REAL, PARAMETER :: tol=0.7, bigno = 1.0e30, hashfac=1.0
      INTEGER, DIMENSION(npart) :: sentarray2, receivearray2
      INTEGER, DIMENSION(mmax) :: isubnode, node, isib, ipar, idau, key,&
                                  next, ihash
      INTEGER, DIMENSION(mmax), TARGET :: listga, listgn
      INTEGER, DIMENSION(nbody) :: nearl, newnearl
      INTEGER :: i, j, k, p, q, s, m, n, nlist, nlistga, nlistgn,       &
                 nnodes, knodes, nroot, iprint, iterm, idisk1, idisk2,  &
                 idisk3, nactive, new, newold, l, nn, ll,               &
                 message, icbrt, icbrt1, nhash, jj, ibin, itemp, ibin1, &
                 nj, np, mm, iz, iy, ix, nc, lllx, llux, llly, lluy,    &
                 lllz, lluz, llx, lux, lly, luy, llz, luz, jx, icjx, jy,&
                 icjy, jz, icjz, lk, nwal, mlocal, newlist, niter
      INTEGER :: nvec, nsig, nmin, nmax, notdone, AllocateStatus=0,     &
                 DeallocateStatus=0
      INTEGER, PARAMETER :: itermax=1000
#ifdef openmp 
      INTEGER :: OMP_GET_THREAD_NUM
#endif
      LOGICAL :: done
!
!--Initialize pointers
!
      list => listga
      nay  => listgn
      rx   => r(1:mmax)
      ry   => r(mmax+1:2*mmax)
      rz   => r(2*mmax+1:3*mmax)
      qxx  => quad(1:mmax)
      qyy  => quad(mmax+1:2*mmax)
      qzz  => quad(2*mmax+1:3*mmax)
      qxy  => quad(3*mmax+1:4*mmax)
      qyz  => quad(4*mmax+1:5*mmax)
      qzx  => quad(5*mmax+1:6*mmax)
!
!--Initialize the matrix containing the list of neighbours for each particle. In MPI calculations
!  there will be a copy of the matrix for every MPI process, containing each copy the list of neighbours
!  for the particles allocated in each MPI process. Since the number of particles in each process
!  (factor) is a dynamical quantity, we have to redefine the size of the matrix everytime we create
!  and store the tree structure.
      IF (ALLOCATED(nb)) DEALLOCATE(nb, STAT=DeallocateStatus)
      IF (DeallocateStatus /= 0 ) THEN
        PRINT*, "Error deallocating nb"
        STOP
      ENDIF
!
      ALLOCATE (nb(nbmax,factor), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
        PRINT*, "Treeconst: Not enough memory for nb!!!"
        STOP
      ENDIF
!
!--Here starts the tree building. Initialize arrays
!
      accpar = tol**2
!$OMP PARALLEL DEFAULT(none) shared(list,xyzhm,rx,ry,rz,em,isib) &
!$OMP shared(idau,qrad,hnode,qxx,qyy,qzz,qxy,qyz,qzx,nbody) &
!$OMP shared(partype) private(j)
!$OMP DO SCHEDULE(runtime)
      DO 10 j=1,nbody
         list(j)  = j
         rx(j)    = xyzhm(1,j)
         ry(j)    = xyzhm(2,j)
         rz(j)    = xyzhm(3,j)
         em(j)    = xyzhm(5,j)
         isib(j)  = 0
         idau(j)  = 0
         qrad(j)  = 0.
         hnode(j) = 4.*xyzhm(4,j)*xyzhm(4,j)
         qxx(j)   = 0.
         qyy(j)   = 0.
         qzz(j)   = 0.
         qxy(j)   = 0.
         qyz(j)   = 0.
         qzx(j)   = 0.
10    ENDDO
!$OMP END DO
!$OMP END PARALLEL
      nactive = nbody
      new     = nbody + 1
!
!--Build binary tree including quadrupole moments
!  return here to fill in each new level of hierarchy:
!
1     CONTINUE
!
!--Find nearest neighbor
!
      icbrt=int(hashfac*float(nactive)**(1./3.))
      icbrt1=icbrt-1
      nhash=icbrt**3.0
      IF (nhash.GT.nactive) THEN
         WRITE(*,*)'nhash greater than n'
         STOP
      ENDIF
!
      IF (icbrt.GT.100) THEN
         WRITE(*,*)'icbrt greater than ncbrt'
         STOP
      ENDIF
!
!--Bin the points by X and fill XMAP
!
      CALL indexx(nactive,list,rx,next)
!
      jj=1
      DO 20 ibin=1,icbrt
        i = (ibin*nactive)/icbrt
        IF (i.NE.nactive) xmap(ibin) = 0.5d0*(rx(list(next(i))) +       &
                    rx(list(next(i+1))))
        itemp = icbrt*(ibin-1)
        DO 30 j=jj,i
           key(next(j))=itemp
30      ENDDO
        jj=i+1
20    ENDDO
!
!--Bin the points by Y and fill YMAP
!
      CALL indexx(nactive,list,ry,next)
!
      jj=1
      DO 40 ibin=1,icbrt
         i = (ibin*nactive)/icbrt
         IF (i.NE.nactive) ymap(ibin) = 0.5d0*(ry(list(next(i))) +      &
                     ry(list(next(i+1))))
         ibin1 = ibin-1
         DO 50 j=jj,i
            nj      = next(j)
            key(nj) = icbrt*(key(nj)+ibin1)
50       ENDDO
         jj=i+1
40    ENDDO
!
!--Bin the points by Z and fill ZMAP
!
      CALL indexx(nactive,list,rz,next)
!
      jj=1
      DO 60 ibin=1,icbrt
         i = (ibin*nactive)/icbrt
         IF (i.NE.nactive) zmap(ibin) = 0.5d0*(rz(list(next(i))) +      &
                     rz(list(next(i+1))))
         DO 70 j=jj,i
            nj      = next(j)
            key(nj) = key(nj) + ibin
70       ENDDO
         jj=i+1
60    ENDDO
!
!--Now fill the head-of-list table and the linked list
!
      DO 80 m=1,nhash
         ihash(m) = 0
80    ENDDO
!
      DO 90 j=1,nactive
         m        = key(j)
         next(j)  = ihash(m)
         ihash(m) = j
90    ENDDO
!
!--Now loop over the particles to find nearest neighbors
!
!$OMP PARALLEL DEFAULT(none) shared(nactive,key,icbrt,icbrt1,list,nay,rx,ry,rz)  &
!$OMP shared(ihash,next,xmap,ymap,zmap,rank,size) private(np,mm,iz,iy,ix,nc)     &
!$OMP private(xnp,ynp,znp,lllx,llux,llly,lluy,lllz,lluz,llx,lux,lly,luy,llz,luz) &
!$OMP private(jx,jy,jz,icjx,icjy,icjz,m,k,lk,d,ddd,dmin,nwal)
!$OMP DO SCHEDULE(runtime)
      DO 100 np=1,nactive
         !PRINT*, 'np=', np, 'rank=', rank
!
!--Find its cell
!
         mm = key(np)-1
         iz = mod(mm,icbrt)
         mm = mm/icbrt
         iy = mod(mm,icbrt)
         ix = mm/icbrt
!
!--Set closest distance so far and coordinates
!
         nc  = 0
         ddd = bigno
         xnp = rx(list(np))
         ynp = ry(list(np))
         znp = rz(list(np))
!
!--Initialize the limits of the search loop
!
         lllx = ix
         llux = ix
         llly = iy
         lluy = iy
         lllz = iz
         lluz = iz
         llx  = ix
         lux  = ix
         lly  = iy
         luy  = iy
         llz  = iz
         luz  = iz
!
!--Loop for finding closest particle in the volume of new cells
!
2        CONTINUE
!
         DO 110 jx=llx,lux
            icjx = icbrt*jx
            DO 120 jy=lly,luy
               icjy = icbrt*(jy+icjx)+1
               DO 130 jz=llz,luz
                  m = icjy+jz
                  k = ihash(m)
                  IF (k.NE.0) THEN
135                  CONTINUE
                     IF (k.NE.np) THEN
                        lk = list(k)
                        d  = (rx(lk)-xnp)*(rx(lk)-xnp) +                &
                             (ry(lk)-ynp)*(ry(lk)-ynp) +                &
                             (rz(lk)-znp)*(rz(lk)-znp)
                        IF (d.LT.ddd) THEN
                           ddd = d
                           nc  = k
                        ENDIF
                     ENDIF
                     k = next(k)
                     IF (k.NE.0) GOTO 135
                  ENDIF
130             ENDDO
120         ENDDO
110      ENDDO
!
!--We must now find the closest wall, if any
!
         dmin = bigno
         nwal = 0
         IF (lllx.GT.0) THEN
            d = xnp-xmap(lllx)
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 1
               ENDIF
            ENDIF
         ENDIF
!
         IF (llux.LT.icbrt1) THEN
            d = xmap(llux+1)-xnp
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 2
               ENDIF
            ENDIF
         ENDIF
!
         IF (llly.GT.0) THEN
            d = ynp-ymap(llly)
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 3
               ENDIF
            ENDIF
         ENDIF
!
         IF (lluy.LT.icbrt1) THEN
            d=ymap(lluy+1)-ynp
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 4
               ENDIF
            ENDIF
         ENDIF
!
         IF (lllz.GT.0) THEN
            d = znp-zmap(lllz)
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 5
               ENDIF
            ENDIF
         ENDIF
!
         IF (lluz.LT.icbrt1) THEN
            d = zmap(lluz+1)-znp
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 6
               ENDIF
            ENDIF
         ENDIF
!
!--Reset the search volume to its augmented value
!
         llx = lllx
         lux = llux
         lly = llly
         luy = lluy
         llz = lllz
         luz = lluz
!
!--Augment it and set the next search according to which wall is closest
!
         IF (nwal.NE.0) THEN
            IF (nwal.EQ.1) THEN
               lllx = lllx-1
               llx  = lllx
               lux  = lllx
            ELSE IF (nwal.EQ.2) THEN
               llux = llux+1
               llx  = llux
               lux  = llux
            ELSE IF (nwal.EQ.3) THEN
               llly = llly-1
               lly  = llly
               luy  = llly
            ELSE IF (nwal.EQ.4) THEN
               lluy = lluy+1
               lly  = lluy
               luy  = lluy
            ELSE IF (nwal.EQ.5) THEN
               lllz = lllz-1
               llz  = lllz
               luz  = lllz
            ELSE IF (nwal.EQ.6) THEN
               lluz = lluz+1
               llz  = lluz
               luz  = lluz
            ENDIF
            GOTO 2
         ENDIF
         nay(np) = nc
100   ENDDO
!$OMP END DO
!$OMP END PARALLEL
      newold = new
!
!--Find new nodes (It would be nice to do this loop in parallel, but I
!  don't see how to do it wihout introducing too many changes)
!
      DO 140 j=1,nactive
         l=list(j)
         IF (isib(l).EQ.0) THEN
            n = nay(j)
            IF (nay(n).EQ.j) THEN
               ll       = list(n)
               isib(ll) = l
               isib(l)  = ll
               ipar(ll) = new
               ipar(l)  = new
               em(new)  = em(l)+em(ll)
               fl       = em(l)/em(new)
               fll      = em(ll)/em(new)
               emred    = fl*fll*em(new)
               difx     = rx(ll)-rx(l)
               dify     = ry(ll)-ry(l)
               difz     = rz(ll)-rz(l)
               rx(new)  = rx(l)+fll*difx
               ry(new)  = ry(l)+fll*dify
               rz(new)  = rz(l)+fll*difz
!
!--find radius and maximum h
!
               rr           = sqrt(difx*difx+dify*dify+difz*difz)
               qrad(new)    = max(fll*rr+qrad(l),fl*rr+qrad(ll))
               xyzhm(4,new) = max(xyzhm(4,l),xyzhm(4,ll))
               hnode(new)   = max(hnode(l),hnode(ll))
!
!--compute quadrupole moments
!
               qxx(new) = (emred*difx)*difx
               qxy(new) = (emred*difx)*dify
               qzx(new) = (emred*difx)*difz
               qyy(new) = (emred*dify)*dify
               qyz(new) = (emred*dify)*difz
               qzz(new) = emred*(difz*difz)
               IF (l.GT.nbody) THEN
                  qxx(new) = qxx(new)+qxx(l)
                  qyy(new) = qyy(new)+qyy(l)
                  qzz(new) = qzz(new)+qzz(l)
                  qxy(new) = qxy(new)+qxy(l)
                  qyz(new) = qyz(new)+qyz(l)
                  qzx(new) = qzx(new)+qzx(l)
               ENDIF
               IF (ll.GT.nbody) THEN
                  qxx(new) = qxx(new)+qxx(ll)
                  qyy(new) = qyy(new)+qyy(ll)
                  qzz(new) = qzz(new)+qzz(ll)
                  qxy(new) = qxy(new)+qxy(ll)
                  qyz(new) = qyz(new)+qyz(ll)
                  qzx(new) = qzx(new)+qzx(ll)
               ENDIF
               isib(new) = 0
               idau(new) = l
               new       = new + 1
               IF (new.GT.mmax) STOP 'new.gt.mmax'
            ENDIF
         ENDIF
140   ENDDO
!
!--Compactify list:
!
      n = 1
      DO 150 j=1,nactive
         IF (isib(list(j)).EQ.0) THEN
            list(n)=list(j)
            n=n+1
         ENDIF
150   ENDDO
      DO 160 j=newold,new-1
         list(n)=j
         n=n+1
160   ENDDO
      nactive=n-1
      IF (nactive.GT.1) THEN
         GOTO 1
      ENDIF
      nroot=new-1
      message=6

      NULLIFY(nay)
      NULLIFY(list)
!
!--End of tree building. Calculate neighbour lists and gravitational forces
!
!$OMP PARALLEL DEFAULT(none) shared(rx,ry,rz,nb,eps,eps3)                   &
!$OMP shared(nroot,qrad,xyzhm,idau,em,hnode,accpar,isib,qxx,qyy,qzz)        &
!$OMP shared(qxy,qzx,qyz,gxyzu,partype,factor,nbody,nbmax,nbmin,rank)       &
!$OMP shared(istep0,istep,step,nstep)                                       &
!$OMP private(m,s,nlist,p,q,rrx,rry,rrz,hm,mlocal)                          &
!$OMP private(pot,eps2,nnodes,node,knodes,i,n,hn,dxi,dyi,dzi,rr2)           &
!$OMP private(rcomp2,qradn2,j,isubnode,nearl,k,rr,rr05,newnearl,newlist)    &
!$OMP private(frx,fry,frz,frxl,fryl,frzl,potl,hup,hdown,hnew,done,niter)
!$OMP DO SCHEDULE(runtime)
      partloop : DO mlocal = 1, factor
         m = mlocal + rank*factor
         IF ((m > nbody) .OR. (partype(m) == 2)) CYCLE
         !IF ((nstep /=0) .AND. (istep0(m) + step(m) /= istep)) CYCLE
!
!--Initialize quantities
!
         rrx     = rx(m)
         rry     = ry(m)
         rrz     = rz(m)
         hm      = hnode(m)
         nlist   = 0
         frx     = 0.
         fry     = 0.
         frz     = 0.
         pot     = 0.
         eps2    = 2.*eps
!
!--Walk the tree level by level to find neighbors and
!  construct list of distant nodes and particles
!
         nnodes  = 1
         node(1) = nroot
3        IF (nnodes.GT.0) THEN
            knodes=0
            DO 180 i=1,nnodes
               n   = node(i)
               hn  = hnode(n)
               dxi = rx(n)-rrx
               dyi = ry(n)-rry
               dzi = rz(n)-rrz
               rr2 = dxi*dxi + dyi*dyi + dzi*dzi
               IF (n.GT.nbody) THEN
                  rcomp2 = (qrad(n)+2.*MAX(xyzhm(4,n),xyzhm(4,m)))**2
                  qradn2 = qrad(n)+eps2
                  qradn2 = qradn2*qradn2
                  IF (rr2.LT.rcomp2) THEN
                     j                = idau(n)
                     knodes           = knodes + 1
                     isubnode(knodes) = j
                  ELSEIF (accpar*rr2.LT.qradn2) THEN
                     j                = idau(n)
                     knodes           = knodes + 1
                     isubnode(knodes) = j
                  ELSE
                     CALL grav_node(rx(n)-rrx,ry(n)-rry,rz(n)-rrz,em(n),&
                                    qxx(n),qyy(n),qzz(n),qxy(n),qzx(n), &
                                    qyz(n),frxl,fryl,frzl,potl)
                     frx = frx + frxl
                     fry = fry + fryl
                     frz = frz + frzl
                     pot = pot + potl
                  ENDIF
               ELSE
                  IF (n.NE.m) THEN
                     CALL grav_dist(rx(n)-rrx,ry(n)-rry,rz(n)-rrz,em(n),&
                                    eps,eps3,frxl,fryl,frzl,potl)
                     frx = frx + frxl
                     fry = fry + fryl
                     frz = frz + frzl
                     pot = pot + potl
                  ENDIF
                  IF (sqrt(rr2) < 0.5*(sqrt(hm)+sqrt(hn))) THEN
                     nlist        = nlist + 1
                     nearl(nlist) = n
                  ENDIF
               ENDIF
180         ENDDO
!
!--Initialise for next walk
!
            k=0
            DO 190 i=1,knodes
               k       = k + 1
               j       = isubnode(i)
               node(k) = j
               k       = k + 1
               node(k) = isib(j)
190         ENDDO
            nnodes=k
            GOTO 3
         ENDIF
!
!--Store gravitational forces and potential
!
         gxyzu(1,mlocal) = frx
         gxyzu(2,mlocal) = fry
         gxyzu(3,mlocal) = frz
         gxyzu(4,mlocal) = -pot
!
!--Store neighbour list
!
         IF ((nlist >= nbmax) .OR. (nlist <= nbmin)) THEN
           ! If the particle has an excess or too few neighbours we will use a
           ! bisection method to recalculate its smoothing length,
           ! so it has a propper number of neighbours
           PRINT*, 'Switching to bisection for particle',m,xyzhm(1,m),  &
           xyzhm(2,m),xyzhm(3,m),xyzhm(4,m),nlist
           if (nlist > nbmax) then
             hup   = xyzhm(4,m)
             hdown = hmin
           elseif (nlist < nbmin) then
             hup   = hmax
             hdown = xyzhm(4,m)
           endif
           done  = .false.
           niter = 0
           iterloop : DO WHILE (done.EQV..false.)
             niter = niter+1
!
             hnew = 0.5*(hup+hdown)
             newlist = 0
             DO p=1,nbody
               hn  = xyzhm(4,p)
               dxi = xyzhm(1,p)-xyzhm(1,m)
               dyi = xyzhm(2,p)-xyzhm(2,m)
               dzi = xyzhm(3,p)-xyzhm(3,m)
               rr2 = dxi**2 + dyi**2 + dzi**2
               IF (sqrt(rr2) < 0.5*(hnew+hn)) THEN
                  newlist = newlist + 1
                  newnearl(newlist) = p
               ENDIF
             ENDDO
!
             IF ((newlist <= nbmax) .AND. (newlist >= nbmin)) THEN
               done = .true.
               PRINT*, 'Bisection ended', m,' hnew=,',hnew, ' nlist=',newlist
             ELSEIF (newlist >= nbmax) THEN
               hup   = hnew
               hdown = hdown
             ELSEIF (newlist <= nbmin) THEN
               hup   = hup
               hdown = hnew
             ENDIF
!
! Kill the particle if the iteration fails
             IF (niter > itermax) THEN
               gxyzu(4,mlocal) = 123456789
               EXIT iterloop
             ENDIF
           ENDDO iterloop
!
           nb(1,mlocal) = newlist
           DO s=1,newlist
             nb(s+1,mlocal) = newnearl(s)
           ENDDO
           xyzhm(4,mlocal) = hnew 
         ELSE
           nb(1,mlocal) = nlist
           DO s=1,nlist
             nb(s+1,mlocal) = nearl(s)
           ENDDO
         ENDIF
      ENDDO partloop
!$OMP END DO
!$OMP END PARALLEL
!
!--Force synchronization
!
#ifdef MPI
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
!
!--Deallocate pointers
!
      NULLIFY(rx)
      NULLIFY(ry)
      NULLIFY(rz)
      NULLIFY(qxx)
      NULLIFY(qyy)
      NULLIFY(qzz)
      NULLIFY(qxy)
      NULLIFY(qyz)
      NULLIFY(qzx)
!
!--Transfer gxyzu across MPI processess
!
#ifdef MPI
      sentarray    = 0.0
      receivearray = 0.0
!$OMP PARALLEL DEFAULT(none) shared(sentarray,gxyzu,factor) &
!$OMP private(p,k)
!$OMP DO SCHEDULE(runtime)
      DO p = 1, factor
        DO k = 1, ndim+1
          sentarray((p-1)*(ndim+1) + k) = gxyzu(k,p)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
      CALL MPI_ALLGATHER(sentarray,(ndim+1)*factor,MPI_DOUBLE_PRECISION,&
                      receivearray,(ndim+1)*factor,MPI_DOUBLE_PRECISION,&
                      MPI_COMM_WORLD,ierr)
!
!$OMP PARALLEL DEFAULT(none) shared(receivearray,gxyzu,nbody,istep,istep0,step,nstep) &
!$OMP private(p,k)
!$OMP DO SCHEDULE(runtime)
      DO p = 1, nbody
        !IF ((nstep == 0) .OR. (istep0(p) + step(p) == istep)) THEN
          DO k = 1, ndim+1
            gxyzu(k,p) = receivearray((p-1)*(ndim+1) + k)
          ENDDO
        !ENDIF
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
#endif
!
!--Kill particles if necessary
!
      notdone = 0
!$OMP PARALLEL DEFAULT(none) shared(nbody,gxyzu,partype,istep0,istep,step,nstep) &
!$OMP private(p) reduction(+:notdone)
!$OMP DO SCHEDULE(runtime)
      DO p=1,nbody
        !IF ((nstep == 0) .OR. (istep0(p) + step(p) == istep)) THEN
          IF (gxyzu(4,p) == 123456789) THEN
            notdone    = notdone + 1
            partype(p) = 2
          ENDIF
        !ENDIF
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
!--Statistics
!
      nmin = 100000
      nmax = 0
      nvec = 0
      nsig = 0
!$OMP PARALLEL DEFAULT(none) shared(nb,rank,factor,nbody) private(m,p)  &
!$OMP reduction(+:nvec,nsig) reduction(MIN:nmin) reduction(MAX:nmax)
!$OMP DO SCHEDULE(runtime) 
      DO m = 1, factor
         p = m + rank*factor
         IF (p <= nbody) THEN
            nvec = nvec + nb(1,m)
            nsig = nsig + nb(1,m)*nb(1,m)
            nmin = MIN(nmin,nb(1,m))
            nmax = MAX(nmax,nb(1,m))
         ENDIF
      END DO
!$OMP END DO 
!$OMP END PARALLEL
!
#ifdef MPI
       CALL MPI_ALLREDUCE(nmin, globnmin, 1, MPI_INTEGER, MPI_MIN,      &
                          MPI_COMM_WORLD, ierr)
       CALL MPI_ALLREDUCE(nmax, globnmax, 1, MPI_INTEGER, MPI_MAX,      &
                          MPI_COMM_WORLD, ierr)
       CALL MPI_ALLREDUCE(nvec, globnvec, 1, MPI_INTEGER, MPI_SUM,      &
                          MPI_COMM_WORLD, ierr)
#else
       globnmin  = nmin
       globnmax  = nmax
       globnvec  = nvec
#endif
!
        IF (rank == MASTER) THEN
          globnvec = INT(FLOAT(globnvec)/FLOAT(nbody))
        ENDIF
!
      END SUBROUTINE treeconst
!
      SUBROUTINE grav_dist(difx,dify,difz,emn,eps,eps3,frx,fry,frz,pot)
!
!--Compute gravitational forces due to distant particles
!
      IMPLICIT NONE
!
!--I/O Variables
!      
      real, INTENT(IN)  :: difx, dify, difz, emn, eps, eps3 
      real, INTENT(OUT) :: frx, fry, frz, pot
!
!--Local variables
      real :: rr, rr05, fp, ga, u, uu2, u3 
!
      frx  = 0.0d0
      fry  = 0.0d0
      frz  = 0.0d0
      pot  = 0.0d0
      rr   = difx*difx+dify*dify+difz*difz
!
      IF (rr.EQ.0.0d0) GOTO 205
!
      rr05 = sqrt(rr)
      fp   = 1d0/rr05
      ga   = 1d0/rr05/rr
      u    = rr05/eps
      IF (u.GT.1d0.AND.u.LT.2d0) THEN
            uu2 = u*u
            u3  = uu2*u
            fp  = -1d0/(15d0*rr05)-((4d0/3d0)-u+.3d0*uu2-(.1d0/3d0)     &
                  *u3)*uu2/eps+1.6d0/eps
            ga  = (-1d0/15d0+u3*((8d0/3d0)-3d0*u+1.2d0*uu2-             &
                  (1d0/6d0)*u3))/rr05/rr
      ELSEIF (u.LT.1d0) THEN
            uu2 = u*u
            u3  = uu2*u
            fp  = (-((1d0/3d0)-.15d0*uu2+.05d0*u3)*2d0*uu2+1.4d0)/eps
            ga  = (4d0/3d0-1.2d0*uu2+.5d0*u3)/eps3
      ENDIF
!
      frx = emn*ga*difx
      fry = emn*ga*dify
      frz = emn*ga*difz
      pot = emn*fp
!
205   CONTINUE
      RETURN
!
      END SUBROUTINE grav_dist           
!
      SUBROUTINE grav_node(difx,dify,difz,emn,qxxn,qyyn,qzzn,qxyn,qzxn, &
                           qyzn,frx,fry,frz,pot)
!
!--Compute gravitational forces and potential due to single distant
!  node. Include quadrupole corrections.
!
      IMPLICIT NONE
!
!--I/O Variables
!
      real, INTENT(IN)  :: difx, dify, difz, emn, qxxn, qyyn, qzzn,  &
                              qxyn, qzxn, qyzn
      real, INTENT(OUT) :: frx, fry, frz, pot
!
!--Local variables
!
      real :: rr, rri, rri05, ff, fpr, fpprr, tx, ty, tz, rqr, fac,  &
                 phimono, phiquad
!
      frx  = 0.0d0
      fry  = 0.0d0
      frz  = 0.0d0
      pot  = 0.0d0
      rr   = difx*difx+dify*dify+difz*difz
!
      IF (rr.EQ.0.0d0) GOTO 215
!
      rri   = 1.0d0/rr
      rri05 = sqrt(rri)
      ff    = rri*rri05
      fpr   = (-3d0)*ff*rri
      fpprr = (-2.5d0)*fpr*rri
      tx    = qxxn*difx+qxyn*dify+qzxn*difz
      ty    = qxyn*difx+qyyn*dify+qyzn*difz
      tz    = qzxn*difx+qyzn*dify+qzzn*difz
      rqr   = difx*tx+dify*ty+difz*tz
      fac   = emn*ff+rqr*fpprr+0.5d0*(qxxn+qyyn+qzzn)*fpr
      frx   = fac*difx+fpr*tx
      fry   = fac*dify+fpr*ty
      frz   = fac*difz+fpr*tz
!
      phimono = emn*rri05
      phiquad = 0.5d0*(rqr*fpr-(qxxn+qyyn+qzzn)*ff)
      pot     = phimono + phiquad
!
215   CONTINUE
      RETURN
      END SUBROUTINE grav_node
