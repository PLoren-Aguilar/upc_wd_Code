      SUBROUTINE energy
!==================================================================
!     Subroutine to calculate system energy, momentums, etc 
!
!     Last revision: 15/March/2015
!==================================================================
!
!--Load modules
!
      USE mod_parameters, ONLY : ndim
      USE mod_commons, ONLY : xyzhm, vxyzut, enuc, gxyzu, nbody, nstep, &
      tnow
!
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local variables
!
      REAL, DIMENSION(3) :: cmpos, cmvel, amvec 
      REAL ::  mass, mtot, ektot, eptot
      REAL ::  eitot, etot, uitot, delta, enn, eneut
      REAL ::  x, y, z, vx, vy, vz 
      INTEGER ::  p, k
!
!--Initializations
!
      mtot  = 0.0
      eptot = 0.0
      uitot = 0.0
      ektot = 0.0
      enn   = 0.0
      eneut = 0.0
      amvec = 0.0
      cmpos = 0.0
      cmvel = 0.0
!
!--Energy calculation
!
!$OMP PARALLEL DEFAULT(none) shared(nbody,xyzhm,vxyzut,gxyzu,enuc) &
!$OMP private(p,mass) reduction(+:mtot,eptot,uitot,ektot,enn)
!$OMP DO SCHEDULE(runtime)
      DO 20 p=1,nbody
         mass  = xyzhm(5,p) 
         mtot  = mtot  + mass
         eptot = eptot + mass*gxyzu(4,p)
         uitot = uitot + mass*vxyzut(4,p)
         ektot = ektot + mass*(vxyzut(1,p)*vxyzut(1,p) +                &
                 vxyzut(2,p)*vxyzut(2,p) + vxyzut(3,p)*vxyzut(3,p))
         enn   = enn + mass*enuc(p)
         !eneut = eneut + xyzhm(5,p)*eneutr(p)
20    ENDDO
!$OMP END DO
!$OMP END PARALLEL
      ektot = 0.5*ektot
      eptot = 0.5*eptot
!
!--Angular momentum calculation
!
!$OMP PARALLEL DEFAULT(none) shared(nbody,xyzhm,vxyzut)  &
!$OMP private(mass,x,y,z,vx,vy,vz,p) reduction(+:amvec)
!$OMP DO SCHEDULE(runtime)
      DO 30 p=1,nbody
         mass = xyzhm(5,p)
         x  = xyzhm(1,p)
         y  = xyzhm(2,p)
         z  = xyzhm(3,p)
         vx = vxyzut(1,p)
         vy = vxyzut(2,p)
         vz = vxyzut(3,p)
!
         amvec(1) = amvec(1) + mass*(y*vz - z*vy)
         amvec(2) = amvec(2) - mass*(x*vz - z*vx)
         amvec(3) = amvec(3) + mass*(x*vy - y*vx)
30    ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
!--Center-of-mass calculations
!
      DO 40 k=1,ndim
!$OMP PARALLEL DEFAULT(none) shared(nbody,k,xyzhm,vxyzut)  &
!$OMP private(p,mass) reduction(+:cmpos,cmvel)
!$OMP DO SCHEDULE(runtime)
         DO 50 p=1,nbody
            mass = xyzhm(5,p)
            cmpos(k) = cmpos(k) + mass*xyzhm(k,p)
            cmvel(k) = cmvel(k) + mass*vxyzut(k,p)
50       ENDDO
!$OMP END DO
!$OMP END PARALLEL
40    ENDDO
      cmpos = cmpos/mtot
      cmvel = cmvel/mtot
!
!--Compute total system energy
!
      etot = ektot + eptot + uitot
!
!--Ouput diagnostics
!
      OPEN (1,FILE='energy.out',STATUS='unknown',access='append')
      OPEN (2,FILE='treelogs.out',STATUS='unknown',access='append')
!
      WRITE(1,'(i7,1x,7(1pe13.5))') nstep, tnow, etot, ektot,           &
      eptot, uitot, enn, eneut
      WRITE(2,'(10(1pe12.4))') tnow,(cmpos(k),k=1,ndim),                &
      SQRT(cmpos(1)**2+cmpos(2)**2+cmpos(3)**2),                       &
      (cmvel(k),k=1,ndim),SQRT(cmvel(1)**2+cmvel(2)**2+cmvel(3)**2),   &
      SQRT(amvec(1)**2+amvec(2)**2+amvec(3)**2)
!
      CLOSE(1)
      CLOSE(2)
!
      END SUBROUTINE energy
