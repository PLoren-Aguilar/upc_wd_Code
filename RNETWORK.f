C  ------------------- NUCLEAR REACTION NETWORK ------------------------
C *********************************************************************         
      SUBROUTINE RNETWORK
C *********************************************************************         
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C---------------------------------------------------------------------
CC    READS THE ISOTOPES AND THE NETWORK
C---------------------------------------------------------------------
      parameter(NRE=29,NIS=14,NSP=NIS+1)
      PARAMETER( AA = 9.6368D17 )

      CHARACTER ONC(NSP)*5

      COMMON/CM/K1(NRE),K2(NRE),K3(NRE),K4(NRE),K5(NRE),K6(NRE),
     &K7(NRE),K8(NRE)
      COMMON/CNAME/ONC
      COMMON/CNETW/AN(NSP),ZN(NSP),Q(NRE),BE(NSP),BE0
      DIMENSION X(NSP)

C
C --  READING OF THE REACTING ISOTOPES
C

      OPEN(UNIT=4,FILE='ZHELI.10b',status='old')

      M1=0
      BE0=0.D0
      SUMX=0.D0
      AST=3.9534D+10
C
      do 10 i=1,NSP
      READ(4,8050) NNNZ,ONC(I),AN(I),ZN(I),X(I),BE(I)
      X(I)=X(I)*AN(I)/AST
 10   continue

      DO 15 I = 1 , NSP
   15    SUMX = SUMX + X(I)
      F = 1./SUMX
      DO 350 I=1,NIS
         X(I)= F * X(I)
  350 continue
      DO 20 I=1,NSP
      BE0 = BE0 + X(I)*BE(I)*AA/AN(I)
      M1=M1+1
      IF(M1.NE.2) GO TO 20
      K=I-1
      M1=0
      L=I
20    CONTINUE

C  ------------------- NUCLEAR REACTION NETWORK ------------------------
      M1=0
      DO 30 L=1,NRE
      M1=M1+1
      READ(4,8060) NNNR,NA,NB,ND,NC,Q(L),IA,IB,ID,IC
      K1(L)=IA
      K2(L)=NA
      K3(L)=IB
      K4(L)=NB
      K5(L)=IC
      K6(L)=NC
      K7(L)=ID
      K8(L)=ND
      IF((M1.NE.2).AND.(L.NE.NRE)) GO TO 30
      K=L-1
      M1=0
      I=L
30    CONTINUE

!         print*,'AN',AN
!         print*,'ZN',ZN
!         print*,'reacciones 1',K1
!         print*,'reacciones 2',K2
!         print*,'reacciones 3',K3
!         print*,'reacciones 4',K4
!         print*,'reacciones 5',K5
!         print*,'reacciones 6',K6
!         print*,'reacciones 7',K7
!         print*,'reacciones 8',K8

      CLOSE(4)
C
 8050 FORMAT(I4,1X,A5,1X,F4.0,1X,F4.0,1X,D11.4,F10.4,2X,D11.4,I5)
 8060 FORMAT(I5,2X,I2,9X,I2,8X,I2,9X,I2,10X,F7.3,9X,3(I4,1X),I4)
C
      RETURN
      END
