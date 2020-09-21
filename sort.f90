      SUBROUTINE sort
!============================================================
!  This subroutine sorts particles to get rid of dead particles
!
!  Last revision: 26/November/2015
!============================================================
!
!--Load modules
!
      USE mod_parameters, ONLY : ndim, nel
      USE mod_commons
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local variables
!
      REAL, DIMENSION(npart) :: temporal
      INTEGER, DIMENSION(npart) :: itemporal
      INTEGER :: p, k
!
!--Sort particle list to get rid of dead particles
!
      eps = 1.0e30
      DO p = 1, npart
         eps = MIN(eps,xyzhm(4,p))
         ilist(p)  = p
         plist(p)  = p
      ENDDO
      eps  = 1.4d0*2.d0*eps
      eps3 = eps*eps*eps
!
      nCO   = 0
      nHe   = 0
      ndead = 0
      DO p=1,npart
         IF (partype(p) == 0) nCO   = nCO + 1       
         IF (partype(p) == 1) nHe   = nHe + 1       
         IF (partype(p) == 2) ndead = ndead + 1       
      END DO
      nbody = npart - ndead
      CALL indexxi2(npart, ilist, partype, plist)
!
!--Readjust parallelisation
!
      factor = nbody/nprocs
      IF (MOD(nbody,factor) /= 0)  factor = factor + 1
!
!--Sort xyzhm
! 
      DO k=1,ndim+2
         DO p=1,npart
            temporal(p) = xyzhm(k,plist(p))
         ENDDO
         DO p=1,npart
            xyzhm(k,p) = temporal(p)
         ENDDO
      ENDDO
!
      DO k=1,ndim+2
         DO p=1,npart
            temporal(p) = xyzhmp(k,plist(p))
         ENDDO
         DO p=1,npart
            xyzhmp(k,p) = temporal(p)
         ENDDO
      ENDDO
!
!--Sort rho
!
      DO p=1,npart
         temporal(p) = rho(plist(p))
      ENDDO
      DO p=1,npart
         rho(p) = temporal(p)
      ENDDO
!
!--Sort vxyzut
!
      DO k=1,ndim+2
         DO p=1,npart
            temporal(p) = vxyzut(k,plist(p))
         ENDDO
         DO p=1,npart
            vxyzut(k,p) = temporal(p)
         ENDDO
      ENDDO
!
      DO k=1,ndim+2
         DO p=1,npart
            temporal(p) = vxyzutp(k,plist(p))
         ENDDO
         DO p=1,npart
            vxyzutp(k,p) = temporal(p)
         ENDDO
      ENDDO
!
!--Sort axyzut
!
      DO k=1,ndim+2
         DO p=1,npart
            temporal(p) = axyzut(k,plist(p))
         ENDDO
         DO p=1,npart
            axyzut(k,p) = temporal(p)
         ENDDO
      ENDDO
!
      DO k=1,ndim+2
         DO p=1,npart
            temporal(p) = axyzutp(k,plist(p))
         ENDDO
         DO p=1,npart
            axyzutp(k,p) = temporal(p)
         ENDDO
      ENDDO
!
!--Sort gxyzu
!
      DO k=1,ndim+1
         DO p=1,npart
            temporal(p) = gxyzu(k,plist(p))
         ENDDO
         DO p=1,npart
            gxyzu(k,p) = temporal(p)
         ENDDO
      ENDDO
!
!--Sort div
!
      DO p=1,npart
         temporal(p) = div(plist(p))
      ENDDO
      DO p=1,npart
         div(p) = temporal(p)
      ENDDO
!
!--Sort divt
!
      DO p=1,npart
         temporal(p) = divt(plist(p))
      ENDDO
      DO p=1,npart
         divt(p) = temporal(p)
      ENDDO
!
!--Sort cur
!
      DO p=1,npart
         temporal(p) = cur(plist(p))
      ENDDO
      DO p=1,npart
         cur(p) = temporal(p)
      ENDDO
!
!--Sort css
!
      DO p=1,npart
         temporal(p) = css(plist(p))
      ENDDO
      DO p=1,npart
         css(p) = temporal(p)
      ENDDO
!
!--Sort press
!
      DO p=1,npart
         temporal(p) = press(plist(p))
      ENDDO
      DO p=1,npart
         press(p) = temporal(p)
      ENDDO
!
!--Sort cvs
!
      DO p=1,npart
         temporal(p) = cvs(plist(p))
      ENDDO
      DO p=1,npart
         cvs(p) = temporal(p)
      ENDDO
!
!--Sort dPdT
!
      DO p=1,npart
         temporal(p) = dPdT(plist(p))
      ENDDO
      DO p=1,npart
         dPdT(p) = temporal(p)
      ENDDO
!
!--Sort cps
!
      DO p=1,npart
         temporal(p) = cps(plist(p))
      ENDDO
      DO p=1,npart
         cps(p) = temporal(p)
      ENDDO
!
!--Sort ka1
!
      DO p=1,npart
         temporal(p) = ka1(plist(p))
      ENDDO
      DO p=1,npart
         ka1(p) = temporal(p)
      ENDDO
!
!--Sort fh
!
      DO p=1,npart
         temporal(p) = fh(plist(p))
      ENDDO
      DO p=1,npart
         fh(p) = temporal(p)
      ENDDO
!
!--Sort uintprev
!
      DO p=1,npart
         temporal(p) = uintprev(plist(p))
      ENDDO
      DO p=1,npart
         uintprev(p) = temporal(p)
      ENDDO
!
!--Sort tscdyn
!
      DO p=1,npart
         temporal(p) = tscdyn(plist(p))
      ENDDO
      DO p=1,npart
         tscdyn(p) = temporal(p)
      ENDDO
!
!--Sort tscnuc
!
      DO p=1,npart
         temporal(p) = tscnuc(plist(p))
      ENDDO
      DO p=1,npart
         tscnuc(p) = temporal(p)
      ENDDO
!
!--Sort enuc
!
      DO p=1,npart
         temporal(p) = enuc(plist(p))
      ENDDO
      DO p=1,npart
         enuc(p) = temporal(p)
      ENDDO
!
      DO p=1,npart
         temporal(p) = enucp(plist(p))
      ENDDO
      DO p=1,npart
         enucp(p) = temporal(p)
      ENDDO
!
!--Sort luminuc
!
      DO p=1,npart
         temporal(p) = luminuc(plist(p))
      ENDDO
      DO p=1,npart
         luminuc(p) = temporal(p)
      ENDDO
!
      DO p=1,npart
         temporal(p) = luminucp(plist(p))
      ENDDO
      DO p=1,npart
         luminucp(p) = temporal(p)
      ENDDO
!
!--Sort dhdt
!
      DO p=1,npart
         temporal(p) = dhdt(plist(p))
      ENDDO
      DO p=1,npart
         dhdt(p) = temporal(p)
      ENDDO
!
      DO p=1,npart
         temporal(p) = dhdtp(plist(p))
      ENDDO
      DO p=1,npart
         dhdtp(p) = temporal(p)
      ENDDO
!
!--Sort partype
!
      DO p=1,npart
         itemporal(p) = partype(plist(p))
      ENDDO
      DO p=1,npart
         partype(p) = itemporal(p)
      ENDDO
!
!--Sort star
!
      DO p=1,npart
         itemporal(p) = star(plist(p))
      ENDDO
      DO p=1,npart
         star(p) = itemporal(p)
      ENDDO
!
!--Sort xss
!
      DO k=1,nel
         DO p=1,npart
            temporal(p) = xss(k,plist(p))
         ENDDO
         DO p=1,npart
            xss(k,p) = temporal(p)
         ENDDO
      ENDDO
!
!--Sort vsigmax
!
      DO p=1,npart
         temporal(p) = vsigmax(plist(p))
      ENDDO
      DO p=1,npart
         vsigmax(p) = temporal(p)
      ENDDO
!
      END SUBROUTINE sort
