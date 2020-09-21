subroutine energy
!==================================================================
!     Subroutine to calculate system energy, momentums, etc 
!
!     Last revision: 15/March/2015
!==================================================================
!
!--Load modules
!
  use mod_parameters, only : ndim
  use mod_commons,    only : xyzhmp, vxyzutp, enucp, gxyzu, nbody, &
                                 nstep, tnow
!
!--Force to declare EVERYTHING
!
  implicit none
!
!--Local variables
!
  real, dimension(3) :: cmpos, cmvel, amvec 
  real ::  mass, mtot, ektot, eptot
  real ::  eitot, etot, uitot, delta, enn, eneut
  real ::  x, y, z, vx, vy, vz 
  integer ::  p, k
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
!$omp parallel default(none) shared(nbody,xyzhmp,vxyzutp,gxyzu,enucp) &
!$omp private(p,mass) reduction(+:mtot,eptot,uitot,ektot,enn)
!$omp do schedule(runtime)
  do p=1,nbody
    mass  = xyzhmp(5,p) 
    mtot  = mtot  + mass
    eptot = eptot + mass*gxyzu(4,p)
    uitot = uitot + mass*vxyzutp(4,p)
    ektot = ektot + mass*(vxyzutp(1,p)*vxyzutp(1,p) +              &
            vxyzutp(2,p)*vxyzutp(2,p) + vxyzutp(3,p)*vxyzutp(3,p))
    enn   = enn + mass*enucp(p)
    !eneut = eneut + xyzhm(5,p)*eneutr(p)
  enddo
!$omp end do
!$omp end parallel
  ektot = 0.5*ektot
  eptot = 0.5*eptot
!
!--Angular momentum calculation
!
!$omp parallel default(none) shared(nbody,xyzhmp,vxyzutp)  &
!$omp private(mass,x,y,z,vx,vy,vz,p) reduction(+:amvec)
!$omp do schedule(runtime)
  do p=1,nbody
    mass = xyzhmp(5,p)
    x  = xyzhmp(1,p)
    y  = xyzhmp(2,p)
    z  = xyzhmp(3,p)
    vx = vxyzutp(1,p)
    vy = vxyzutp(2,p)
    vz = vxyzutp(3,p)
!
    amvec(1) = amvec(1) + mass*(y*vz - z*vy)
    amvec(2) = amvec(2) - mass*(x*vz - z*vx)
    amvec(3) = amvec(3) + mass*(x*vy - y*vx)
  enddo
!$omp end do
!$omp end parallel
!
!--Center-of-mass calculations
!
  do k=1,ndim
!$omp parallel default(none) shared(nbody,k,xyzhmp,vxyzutp)  &
!$omp private(p,mass) reduction(+:cmpos,cmvel)
!$omp do schedule(runtime)
    do p=1,nbody
      mass = xyzhmp(5,p)
      cmpos(k) = cmpos(k) + mass*xyzhmp(k,p)
      cmvel(k) = cmvel(k) + mass*vxyzutp(k,p)
    enddo
!$omp end do
!$omp end parallel
  enddo
  cmpos = cmpos/mtot
  cmvel = cmvel/mtot
!
!--Compute total system energy
!
  etot = ektot + eptot + uitot
!
!--Ouput diagnostics
!
  open (1,file='energy.out',status='unknown',access='append')
  open (2,file='treelogs.out',status='unknown',access='append')

  write(1,'(i7,1x,7(1pe13.5))') nstep, tnow, etot, ektot, eptot, uitot, enn, eneut
  write(2,'(10(1pe12.4))') tnow,(cmpos(k),k=1,ndim), sqrt(cmpos(1)**2+cmpos(2)**2+cmpos(3)**2), &
       (cmvel(k),k=1,ndim),sqrt(cmvel(1)**2+cmvel(2)**2+cmvel(3)**2), sqrt(amvec(1)**2+amvec(2)**2+amvec(3)**2)

  close(1)
  close(2)

end subroutine energy
