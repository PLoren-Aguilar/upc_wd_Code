      SUBROUTINE corrector
!============================================================
!  This subroutine corrects the the integrated quantities 
!
!  Last revision: 15/March/2015
!============================================================
!
!--Load modules
!
      USE mod_parameters, ONLY : tmin, tmax, ndim
      USE mod_commons, ONLY : xyzhm, vxyzut, axyzut, xyzhmp, vxyzutp,   &
      axyzutp, dtmp_min, nbody, enuc, enucp, luminuc, luminucp, RELFLAG,&
      dhdt, dhdtp, npart
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local variables
!
      REAL, DIMENSION(npart) :: temporal
      REAL :: fa, dt1, dt2, dt3, denom, dtmpmin_old
      REAL, PARAMETER :: dtfact_b=0.3d0, dtfact_h=0.3d0, sigma=1.0d0,&
                            f1 = 1./256., f2 = 255./256.
      INTEGER, DIMENSION(npart) :: itemporal
      INTEGER :: p, k
!
!--Save dtmp_min for later use
!
      dtmpmin_old = dtmp_min  
!
!--Correct positions
!
      DO k=1,ndim
!$OMP PARALLEL DEFAULT(none) shared(nbody,xyzhm,xyzhmp,vxyzut,vxyzutp)  &
!$OMP shared(axyzut,axyzutp,dtmp_min,k) private(p)
!$OMP DO SCHEDULE(runtime)
         DO p=1,nbody
#ifdef Helium
            xyzhm(k,p) = xyzhmp(k,p) + (f1*vxyzutp(k,p) +               &
                         f2*vxyzut(k,p))*dtmp_min
#else
            xyzhm(k,p) = xyzhm(k,p) + 0.1667d0*(axyzut(k,p)-            &
                         axyzutp(k,p))*dtmp_min*dtmp_min
#endif
         ENDDO
!$OMP END DO
!$OMP END PARALLEL
      ENDDO
!
!--Correct velocities
!
      DO k=1,ndim
!$OMP PARALLEL DEFAULT(none) shared(vxyzut,vxyzutp,axyzut,axyzutp,dtmp_min,k)   &
!$OMP shared(nbody) private(p)
!$OMP DO SCHEDULE(runtime)
         DO p=1,nbody
#ifdef Helium
            vxyzut(k,p) = vxyzutp(k,p) + (f1*axyzutp(k,p) +             &
                         f2*axyzut(k,p))*dtmp_min
#else
            vxyzut(k,p) = vxyzut(k,p) + 0.5*(axyzut(k,p)-               &
                          axyzutp(k,p))*dtmp_min
#endif
         ENDDO
!$OMP END DO
!$OMP END PARALLEL
      ENDDO
!
!--Correct thermal energies and temperatures if necessary
!
!$OMP PARALLEL DEFAULT(none) shared(vxyzut,vxyzutp,axyzut,axyzutp,dtmp_min)     &
!$OMP shared(nbody,enuc,enucp,luminuc,luminucp,RELFLAG) private(p)
!$OMP DO SCHEDULE(runtime)
      DO p=1,nbody
#ifdef Helium
         vxyzut(4,p) = vxyzutp(4,p) + (f1*axyzutp(4,p) +                &
                       f2*axyzut(4,p))*dtmp_min + enuc(p) - enucp(p)
#else
         IF ((vxyzut(5,p) > tmin).AND.(vxyzut(5,p) < tmax)) THEN
            vxyzut(4,p) = vxyzut(4,p) + 0.5*(axyzut(4,p)-               &
                       axyzutp(4,p))*dtmp_min + enuc(p) - enucp(p)
         ENDIF
#endif

      ENDDO
!$OMP END DO
      IF (RELFLAG.EQV..false.) THEN
!$OMP DO SCHEDULE(runtime)
         DO p=1,nbody
#ifdef Helium
            vxyzut(5,p) = vxyzutp(5,p) + (f1*axyzutp(5,p) +             &
                       f2*axyzut(5,p))*dtmp_min + luminuc(p)-luminucp(p)
#else
            vxyzut(5,p) = vxyzut(5,p) + 0.5*(axyzut(5,p)-axyzutp(5,p))* &
                          dtmp_min + luminuc(p) - luminucp(p)
#endif

         ENDDO
!$OMP END DO
      ENDIF
!$OMP END PARALLEL
!
!--Correct He layer velocity
!
#ifdef Helium
         CALL norm_layer
         CALL layer
#endif
!
!--Save quantities for later use
!
      DO k=1,5 
!$OMP PARALLEL DEFAULT(none) shared(axyzut,axyzutp,vxyzut,vxyzutp,xyzhm,xyzhmp) &
!$OMP shared(k,nbody) private(p)
!$OMP DO SCHEDULE(runtime)
         DO p=1,nbody
            xyzhmp(k,p)  = xyzhm(k,p)
            vxyzutp(k,p) = vxyzut(k,p)
            axyzutp(k,p) = axyzut(k,p)
         ENDDO
!$OMP END DO
!$OMP END PARALLEL
      ENDDO
!
!$OMP PARALLEL DEFAULT(none) shared(enucp,enuc,luminucp,luminuc,dhdtp)  &
!$OMP shared(dhdt,nbody) private(p)
!$OMP DO SCHEDULE(runtime)
      DO p=1,nbody
         enucp(p)    = enuc(p)
         luminucp(p) = luminuc(p)
         dhdtp(p)    = dhdt(p)
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
      END SUBROUTINE corrector
