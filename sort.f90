subroutine sort
!=================================================================
!  This subroutine sorts particles by their individual time-steps
!  and to get rid of dead particles
!
!  Last revision: 9/March/2017
!=================================================================
!
!--Load modules
!
  use mod_parameters, only : ndim, nel
  use mod_commons
!
!--Force to declare EVERYTHING
!
  IMPLICIT NONE
!
!--Local variables
!
  real,    dimension(npart) :: temporal
  integer, dimension(npart) :: itemporal
  integer :: p, k
!
!--Sort particle list as a function of its integration time-step
!
  eps = 1.0e30
  do p = 1, npart
    eps = MIN(eps,xyzhm(4,p))
    ilist(p)  = p
    plist(p)  = p
  enddo
  eps  = 1.4d0*2.d0*eps
  eps3 = eps*eps*eps
!
  nCO   = 0
  nHe   = 0
  ndead = 0
  do p=1,npart
    if (partype(p) == 0) nCO   = nCO + 1       
    if (partype(p) == 1) nHe   = nHe + 1       
    if (partype(p) == 2) ndead = ndead + 1       
  end do
  nbody = npart - ndead
  call indexxi2(npart, ilist, partype, plist)
!
!--Readjust parallelisation
!
  factor = nbody/nprocs
  if (mod(nbody,factor) /= 0)  factor = factor + 1
!
!--Sort step
!
  do p=1,npart
    itemporal(p) = step(plist(p))
  enddo
  do p=1,npart
    step(p) = itemporal(p)
  enddo
!
!--Sort istep0
!
  do p=1,npart
    itemporal(p) = istep0(plist(p))
  enddo
  do p=1,npart
    istep0(p) = itemporal(p)
  enddo
!
!--Sort xyzhm
! 
  do k=1,ndim+2
    do p=1,npart
      temporal(p) = xyzhm(k,plist(p))
    enddo
    do p=1,npart
      xyzhm(k,p) = temporal(p)
    enddo
  enddo
!
  do k=1,ndim+2
    do p=1,npart
      temporal(p) = xyzhmp(k,plist(p))
    enddo
    do p=1,npart
      xyzhmp(k,p) = temporal(p)
    enddo
  enddo
!
!--Sort rho
!
  do p=1,npart
    temporal(p) = rho(plist(p))
  enddo
  do p=1,npart
    rho(p) = temporal(p)
  enddo
!
!--Sort vxyzut
!
  do k=1,ndim+2
    do p=1,npart
      temporal(p) = vxyzut(k,plist(p))
    enddo
    do p=1,npart
      vxyzut(k,p) = temporal(p)
    enddo
  enddo
!
  do k=1,ndim+2
    do p=1,npart
      temporal(p) = vxyzutp(k,plist(p))
    enddo
    do p=1,npart
      vxyzutp(k,p) = temporal(p)
    enddo
  enddo
!
!--Sort axyzut
!
  do k=1,ndim+2
    do p=1,npart
      temporal(p) = axyzut(k,plist(p))
    enddo
    do p=1,npart
      axyzut(k,p) = temporal(p)
    enddo
  enddo
!
  do k=1,ndim+2
    do p=1,npart
      temporal(p) = axyzutp(k,plist(p))
    enddo
    do p=1,npart
      axyzutp(k,p) = temporal(p)
    enddo
  enddo
!
!--Sort gxyzu
!
  do k=1,ndim+1
    do p=1,npart
      temporal(p) = gxyzu(k,plist(p))
    enddo
    do p=1,npart
      gxyzu(k,p) = temporal(p)
    enddo
  enddo
!
!--Sort div
!
  do p=1,npart
    temporal(p) = div(plist(p))
  enddo
  do p=1,npart
    div(p) = temporal(p)
  enddo
!
!--Sort divt
!
  do p=1,npart
    temporal(p) = divt(plist(p))
  enddo
  do p=1,npart
    divt(p) = temporal(p)
  enddo
!
!--Sort cur
!
  do p=1,npart
    temporal(p) = cur(plist(p))
  enddo
  do p=1,npart
    cur(p) = temporal(p)
  enddo
!
!--Sort css
!
  do p=1,npart
    temporal(p) = css(plist(p))
  enddo
  do p=1,npart
    css(p) = temporal(p)
  enddo
!
!--Sort press
!
  do p=1,npart
    temporal(p) = press(plist(p))
  enddo
  do p=1,npart
    press(p) = temporal(p)
  enddo
!
!--Sort cvs
!
  do p=1,npart
    temporal(p) = cvs(plist(p))
  enddo
  do p=1,npart
    cvs(p) = temporal(p)
  enddo
!
!--Sort dPdT
!
  do p=1,npart
    temporal(p) = dPdT(plist(p))
  enddo
  do p=1,npart
    dPdT(p) = temporal(p)
  enddo
!
!--Sort cps
!
  do p=1,npart
    temporal(p) = cps(plist(p))
  enddo
  do p=1,npart
    cps(p) = temporal(p)
  enddo
!
!--Sort ka1
!
  do p=1,npart
    temporal(p) = ka1(plist(p))
  enddo
  do p=1,npart
    ka1(p) = temporal(p)
  enddo
!
!--Sort fh
!
  do p=1,npart
    temporal(p) = fh(plist(p))
  enddo
  do p=1,npart
    fh(p) = temporal(p)
  enddo
!
!--Sort uintprev
!
  do p=1,npart
    temporal(p) = uintprev(plist(p))
  enddo
  do p=1,npart
    uintprev(p) = temporal(p)
  enddo
!
!--Sort tscdyn
!
  do p=1,npart
    temporal(p) = tscdyn(plist(p))
  enddo
  do p=1,npart
    tscdyn(p) = temporal(p)
  enddo
!
!--Sort tscnuc
!
  do p=1,npart
    temporal(p) = tscnuc(plist(p))
  enddo
  do p=1,npart
    tscnuc(p) = temporal(p)
  enddo
!
!--Sort enuc
!
  do p=1,npart
    temporal(p) = enuc(plist(p))
  enddo
  do p=1,npart
    enuc(p) = temporal(p)
  enddo
!
  do p=1,npart
    temporal(p) = enucp(plist(p))
  enddo
  do p=1,npart
    enucp(p) = temporal(p)
  enddo
!
!--Sort luminuc
!
  do p=1,npart
    temporal(p) = luminuc(plist(p))
  enddo
  do p=1,npart
    luminuc(p) = temporal(p)
  enddo
!
  do p=1,npart
    temporal(p) = luminucp(plist(p))
  enddo
  do p=1,npart
    luminucp(p) = temporal(p)
  enddo
!
!--Sort partype
!
  do p=1,npart
    itemporal(p) = partype(plist(p))
  enddo
  do p=1,npart
    partype(p) = itemporal(p)
  enddo
!
!--Sort star
!
  do p=1,npart
    itemporal(p) = star(plist(p))
  enddo
  do p=1,npart
    star(p) = itemporal(p)
  enddo
!
!--Sort xss
!
  do k=1,nel
    do p=1,npart
      temporal(p) = xss(k,plist(p))
    enddo
    do p=1,npart
      xss(k,p) = temporal(p)
    enddo
  enddo
!
!--Sort vsigmax
!
  do p=1,npart
    temporal(p) = vsigmax(plist(p))
  enddo
  do p=1,npart
    vsigmax(p) = temporal(p)
  enddo
!
end subroutine sort
