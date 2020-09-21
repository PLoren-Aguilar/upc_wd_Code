      SUBROUTINE indata
!===================================================================
!  This subroutine reads the initial data.
!
!  Last revision: 15/March/2015
!===================================================================
!
!--Load modules
!
      USE mod_parameters, ONLY : ndim, nel, MASTER
      USE mod_commons
!     
!--Force to declare EVERYTHING
!
      IMPLICIT NONE
!
!--Local variables
!
      INTEGER       :: p, k, m, nums
      INTEGER       :: AllocateStatus=0
      CHARACTER(30) :: filename1, filename2, numstep
!
!--Declare structures in order to read data
!
      TYPE :: outdata
        REAL(4), DIMENSION(ndim+2) :: ioxyzhm
        REAL(4), DIMENSION(ndim+2) :: iovxyzut
        REAL(4) :: iorho
        REAL(4) :: ioka1
        INTEGER(1) :: iopartype
        INTEGER(1) :: iostar
      END TYPE
      TYPE(outdata), ALLOCATABLE, DIMENSION(:) :: data_array
!
      TYPE :: coutdata
        REAL(4), DIMENSION(nel) :: comp
      END TYPE
      TYPE(coutdata), ALLOCATABLE, DIMENSION(:) :: cdata_array
!
!--Open data unit
!
      WRITE(numstep,'(i4.4)') nstep_infile
      filename1 = 'bodi'//trim(numstep)//'.out'
      filename2 = 'comp'//trim(numstep)//'.out'
      OPEN (1, FILE=filename1, FORM='unformatted')
      OPEN (2, FILE=filename2, FORM='unformatted')

      !WRITE(filename,'(a,i4.4,a)') 'bodi',nstep_ini,'.out'
      !OPEN (1, FILE=filename, FORM='unformatted')
      !WRITE(filename,'(a,i4.4,a)') 'comp',nstep_ini,'.out'
      !OPEN (2, FILE=filename, FORM='unformatted')
!
!--Read time and npart data
!
      READ(1) tnow, nstep, npart, nbody, ndead, nbody1, nbody2, Omega0
!
!--Allocate space for the bodi and comp data
!
      ALLOCATE (data_array(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for data_array!!!"
         STOP
      ENDIF
!
      ALLOCATE (cdata_array(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for cdata_array!!!"
         STOP
      ENDIF

      !READ(1) tnow, nums, data_array
!
!--Check if the number of particles is consistent
!
      !IF (nums /= npart) THEN
      !   IF (rank == MASTER) PRINT*, 'nparts NE nbody!!!'
      !   STOP
      !ENDIF
!
!--Read bodi and comp data
!
      READ(1) data_array
      READ(2) cdata_array
!
!--Close files
!
      CLOSE(UNIT=1)
      CLOSE(UNIT=2)
!
!--Allocate particle data arrays and vectors
!
      ALLOCATE (xyzhm(ndim+2,2*npart+2), xyzhmp(ndim+2,2*npart+2), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for xyzhm!!!"
         STOP
      ENDIF
!
      ALLOCATE (vxyzut(ndim+2,npart), vxyzutp(ndim+2,npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for vxyzut!!!"
         STOP
      ENDIF
!
      ALLOCATE (axyzut(ndim+2,npart), axyzutp(ndim+2,npart), STAT =  AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for axyzut!!!"
         STOP
      ENDIF
!
      ALLOCATE (gxyzu(ndim+1,npart), gxyzup(ndim+1,npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for axyzut!!!"
         STOP
      ENDIF
!
      ALLOCATE (rho(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for rho!!!"
         STOP
      ENDIF
!
      ALLOCATE (ka1(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ka1!!!"
         STOP
      ENDIF
!
      ALLOCATE (partype(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for partype!!!"
         STOP
      ENDIF
!
      ALLOCATE (eosflag(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for partype!!!"
         STOP
      ENDIF
!
      ALLOCATE (star(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for star!!!"
         STOP
      ENDIF
!
      ALLOCATE (rotforc(ndim,npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for rotforc!!!"
         STOP
      ENDIF
!
      ALLOCATE (xss(nel,npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for xss!!!"
         STOP
      ENDIF
!
!--Transfer data into the appropriate locations
!
      DO p = 1,npart
         xyzhm(1:ndim+2,p)  = data_array(p)%ioxyzhm(1:ndim+2)
         vxyzut(1:ndim+2,p) = data_array(p)%iovxyzut(1:ndim+2)
         rho(p)             = data_array(p)%iorho
         ka1(p)             = data_array(p)%ioka1
         partype(p)         = data_array(p)%iopartype
         star(p)            = data_array(p)%iostar
         xss(1:nel,p)       = cdata_array(p)%comp(1:nel)
      ENDDO
!
!--Start up the rest of variables
!
      ALLOCATE (vsigmax(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for vsigmax!!!"
         STOP 
      ENDIF
!     
      ALLOCATE (cur(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for cur!!!"
         STOP 
      ENDIF
!     
      ALLOCATE (div(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for div!!!"
         STOP
      ENDIF
!
      ALLOCATE (divt(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for divt!!!"
         STOP
      ENDIF
!
      ALLOCATE (luminuc(npart), luminucp(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for luminuc!!!"
         STOP
      ENDIF
!
      ALLOCATE (enuc(npart), enucp(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for enuc!!!"
         STOP
      ENDIF
!
      ALLOCATE (tscdyn(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for tscdyn!!!"
         STOP 
      ENDIF
!     
      ALLOCATE (tscnuc(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for tscnuc!!!"
         STOP 
      ENDIF
!
      ALLOCATE (ilist(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for ilist!!!"
         STOP
      ENDIF
!
      ALLOCATE (plist(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for plist!!!"
         STOP
      ENDIF
!
      ALLOCATE (css(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for css!!!"
         STOP
      ENDIF
!
      ALLOCATE (press(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for press!!!"
         STOP
      ENDIF
!
      ALLOCATE (cvs(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for cvs!!!"
         STOP
      ENDIF
!
      ALLOCATE (dPdT(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for dPdT!!!"
         STOP
      ENDIF
!
      ALLOCATE (cps(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for cps!!!"
         STOP
      ENDIF
!
      ALLOCATE (fh(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for fh!!!"
         STOP
      ENDIF
!
      ALLOCATE (uintprev(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for uintprev!!!"
         STOP
      ENDIF
!
      ALLOCATE (dhdt(npart), dhdtp(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for dhdt!!!"
         STOP
      ENDIF
!
      ALLOCATE (dtnuc(npart), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT*, "Not enough memory for dhdt!!!"
         STOP
      ENDIF
!
      nstep_ini = nstep
      vsigmax = 1.0
      cur     = 0.0
      div     = 0.0
      luminuc = 0.0
      enuc    = 0.0
      tscdyn  = 0.0
      tscnuc  = 0.0
      dtnuc   = 0.0
!
      eps = 1d30
      DO p = 1, npart
         eps = MIN(eps,xyzhm(4,p))
         ilist(p)  = p
         plist(p)  = p
      ENDDO
      eps  = 1.4d0*2.d0*eps
      eps3 = eps*eps*eps
!
!--Sort particle list to get rid of dead particles
!
      CALL sort
#ifdef debug
      IF (rank == MASTER) PRINT*, 'sort called'
#endif
!
!--Set the maximum number of neighbours and nodes in the tree and
!  allocate space for the neighbour list
!
      nbmax = 5000
      mmax  = 2*npart + 2
!
      !ALLOCATE (nb(nbmax,npart), STAT = AllocateStatus)
      !IF (AllocateStatus /= 0) THEN
      !   PRINT*, "Not enough memory for nb!!!"
      !   STOP
      !ENDIF
!
      END SUBROUTINE indata
