C ***************************************************************
C                         SUBROUTINES
C ***************************************************************
C======================================================================
      SUBROUTINE SNUC(T,DELTA,RHO2,TAU2,SUMY,DTOLD,NCAPA,iread)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C----------------------------------------------------------------------
CC    NUCLEOSYNTHESIS FOR THE TIME-STEP DELTA
C----------------------------------------------------------------------
      PARAMETER (LSHE=1)
      parameter(ngrid=60,NRE=29,NIS=14,NSP=NIS+1) 
      parameter( e = 9.6368D17, O = 1.D-100, XNAV = 6.023D23 )
C----------------------------------------------------------------------
      COMMON/CM/K1(NRE),K2(NRE),K3(NRE),K4(NRE),K5(NRE),K6(NRE),
     1K7(NRE),K8(NRE)
      COMMON/CNETW/AN(NSP),ZN(NSP),Q(NRE),BE(NSP),BE0
      COMMON/CPARM/PAR,PAMIN,DYMIN,DEMIN,YMIN,YTMIN,DEPS,ATEST,TES
      COMMON/CTEST/DEY2,KTEST2
      COMMON/CFACT/DF(0:6)
      COMMON/CSTEP/NSNUC,NREAL
      COMMON/CJORDI/IHENRI
      common/cvit/V2(NRE),tgrid(ngrid),aNegrid(nre,11)
      COMMON/CODER/AE2(NRE)
      COMMON/SNUC2/YY0(NSP),YY(NSP),YYT(NSP),DIV2(NSP)
      COMMON/CSNC2/YYB(NSP)
      COMMON/NUCLEOS/XY(NSP)
      COMMON/ABUND/XYTOT(NSP,LSHE)
      DIMENSION AMAT2(NSP,NSP)

      REAL ar(NIS),ar11,coefg,coeft
 
      LOGICAL C
C
      DELT=DELTA
      NSNUC=0
C.....COMPUTATION OF REACTION RATES AT TAU-RHO
      pme2=0.D0

      DO 2335 I=1,NSP
      XY(I)=XYTOT(I,NCAPA)
 2335 CONTINUE

      DO 5 I=1,NSP
      YY0(I)=XY(I)/AN(I)
      YYB(I)=0.D0
      YYT(I)=0.D0
      pme2=pme2+YY0(I)*ZN(I) 
   5  CONTINUE
      pme2=1.D0/pme2
c     WRITE(7,*) 'Cridem VIT', pme2
      CALL vit(TAU2,RHO2,pme2,iread)

c  incluyo screenings (no en las fotodesintegraciones! --> if)            
      CALL genpar(RHO2,TAU2/1.D9,AN,ZN,YY0,1.D0/pme2,
     &              ar,ar11,coefg,coeft)

      DO I=1,NRE
         IF(K3(I).NE.NSP) THEN
         V2(I)=V2(I)*scrng(ZN(K1(I)),ZN(K3(I)),AN(K1(I)),AN(K3(I)),
     &         ar(K1(I)),ar(K3(I)),coefg,coeft)
         END IF
      END DO
      V2(1)=V2(1)*scrng(ZN(1),ZN(1),AN(1),AN(1),ar(1),ar(1),
     &      coefg,coeft)
      V2(1)=V2(1)*scrng(ZN(1),4.0d0,AN(1),8.0d0,ar(1),ar11,
     &      coefg,coeft)

c
   10 C=.FALSE.
      NNNN=0
      NL=0
      PRO=TES
      PROHE=TES
 
110   DO 120 I=1,NSP
      YY(I)=YY0(I)
  120 CONTINUE
C.......................................................................
C     CONSTRUCTION AND SOLUTION OF THE LINEARIZED SYSTEM OF EQUATIONS
C     BY THE 2-STEP METHOD OF WAGONER(1969:AP.J.SUP.18,P.247)
C     CALCULATION OF THE MATRIX ELEMENTS (SUBR MATRIX)
C     SOLUTION OF THE LINEARIZED SYSTEM : A(I,J)*Y(J) = Y(I)(SUBR EIGEN
C.......................................................................
  130 DO 230 ISTEP=1,2
      IF(C) THEN
          DO 150 I=1,NIS
  150        YY(I)=YY0(I)
      ENDIF

      CALL MATRIX(DELT,V2,YY,AMAT2)
C
      DO 170 I=1,NIS
  170 YYT(I)=YY0(I)
      CALL EIGEN(YYT,YYB,AMAT2,DIV2)

      DO 210 I=1,NIS
  210 IF(YYB(I).LT.YMIN) YYB(I)=YMIN
c     CALL EQUIL(YYB,YMIN,V2)
      IF( ISTEP .LT. 2 ) THEN
          DO 2 I=1,NIS
             DYTT=(YYB(I)-YY0(I))/DELT
             IF(ABS(DYTT).LT.DYMIN) DYTT=DYMIN
             YY(I)=YY0(I)+DYTT*DELT
    2     CONTINUE
          C=.FALSE.
      ENDIF
  230 CONTINUE
C      do i=1,nis
C         print*, 'yy(i),yyb(i)',yy(i),yyb(i)
C      end do 
C         DETERMINATION OF THE NEW TIME-STEP, USING THE LARGEST
C                       ABUNDANCE VARIATION.
 250  NSNUC=NSNUC+1
      IF(NL.EQ.1) GO TO 280
      PRO1=PRO
 280  CALL STEST(YYB,YY0,PROHE,KTEST2,DEY2,TAU2)
      AKHE=PRO1/PROHE
      IF((AKHE.LE.ATEST)) GOTO 340
C       LARGEST ABUNDANCE VARIATION TOO LARGE : RESULTS CANCELLED
C       TIME-STEP REDUCED BY PAMIN,AND BACK TO 110 FOR NEW CALCULATION.
      DELT=DELT/PAMIN
      NL=1
      NNNN=NNNN+1
      IF(NNNN.GT.20) THEN
      WRITE(7,*) 'Problems with Largest Abundance Variation'
      WRITE(7,*) 'N. iter; AKHE; ATEST; DELT',NNNN,AKHE,ATEST,DELT
      WRITE(7,*) 'ICO, Shell no.:',NREAL,NCAPA
      GOTO 340
      ENDIF
      C=.TRUE.
      GO TO 110
C     340:       RESULTS ACCEPTABLE: PRECISION AND NEXT TIME-STEP
  340 NL=0
      DO 888 I=1,NRE
      AE2(I)=0.D0
  888 CONTINUE

      DO 350 M=1,NRE
      I1=K1(M)
      I2=K2(M)
      I3=K3(M)
      I4=K4(M)
      I21=I2-1
      I41=I4-1
      B1B=((YY(I1)+O)**I21)*((YY(I3)+O)**I4)*YYB(I1)*I2
      B2B=((YY(I3)+O)**I41)*((YY(I1)+O)**I2)*YYB(I3)*I4
      AE2(M)=e*Q(M)*V2(M)*(B1B+B2B)/(DF(I2)*DF(I4)*(I2+I4)) 
  350 CONTINUE

      SUMY=0.D0
      DO 360 I=1,NSP
      YY0(I)=YYB(I)
      XY(I)=AN(I)*YY0(I)
      SUMY=SUMY+XY(I)
  360 CONTINUE
c      CAMBIO!!! comento siguiente IF
c      IF((NREAL.EQ.1).OR.(IHENRI.EQ.1)) THEN
c       SUMY=1.D0
c      ENDIF
      SEHE=SUMY-1.D0
      IF((ABS(SEHE).GT.DEPS))  GO TO 430
C
      DTOLD=DELT
      PRO=PROHE
      DELTA=PAR*PRO*DELT

      RETURN


C FOR FIXED TIME-STEP DELTA, USE THE FOLLOWING :
C
  430 WRITE(7,*) 'ABS(SEHE).GT.DEPS',ABS(SEHE),DEPS
C      CALL EXIT                         
      RETURN
  555 WRITE(30,9998) NNNN
C      CALL EXIT                               
C
 9260 FORMAT(///,5X,'SUMX ERROR GREATER THAN :',1PD8.2,5X,'SUMX=',
     &0PF20.16,5X,'STEP NO :',I5,///,
     &5X,'CHECK :   1) BARYON NUMBER CONSERVATION IN THE NETWORK',/,
     &5X,'          2) NON-NEGATIVE REACTION RATES IN  VIT',/,
     &5X,'          3) MATRIX INVERSION ROUTINE (EIGEN OR OTHER)',///,
     &2X,6('*- PROGRAM ABORTED -*'))
 9995 FORMAT(2X,I2,2X,'SNUC, RATDX=',1PD10.3)
 6130 FORMAT(2X,I2,2X,1PD10.3,2X,I4)
 7000 FORMAT(2X,I3,6(2X,I2))
 7001 FORMAT(4(2X,1PD10.3))
 7002 FORMAT(3(2X,1PD10.3))
 7003 FORMAT(2X,I3,4(2X,1PD10.3))
 7004 FORMAT(2X,I2,2X,1PD10.3)
 7005 FORMAT(2(2X,1PD10.3))
 8003 FORMAT(2X,I3,2X,1PD12.5)
 2130 FORMAT(2X,'Reaction :',I3,2X,'V=',1PD12.5,2X,'AE=',1PD12.5)
 9996 FORMAT(///,132('='),///,2X,'NSNUC =',I5,' GREATER THAN NLAST =',
     &  I5)
 9998 FORMAT(///,132('='),///,2X,'NO CONVERGENCE',I5,' ITERATIONS')
 9911 FORMAT(2X,I3,2(2X,1PD12.5))
 700  RETURN
      END
C======================================================================
      SUBROUTINE MATRIX(DELT,V,Y,AT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C-----------------------------------------------------------------------
CC    COMPUTES MATRIX ELEMENTS FROM THE NETWORK
C-----------------------------------------------------------------------
      PARAMETER (NIS=14, NRE=29, NSP=NIS+1)
      PARAMETER(NDS=1,NDI=1,IJD=4,IDEL=IJD-1,JDEL=IJD-1)
      COMMON/CM/K1(NRE),K2(NRE),K3(NRE),K4(NRE),K5(NRE),K6(NRE),
     1K7(NRE),K8(NRE)
      COMMON/CFACT/DF(0:6)
      DIMENSION Y(NSP),V(NRE),AT(NSP,NSP)
 
      O=1.D-100
      AIN=0.D0
C.......................................................................
C     FOR EVERY ISOTOPE I,CALCULATE THE MATRIX LINE(A(I,J),J=I,NSP)
C     FOR EVERY REACTION M : N(II)*II+N(IJ)*IJ ---> N(IK)*IK+N(IL)*IL
C     THERE IS A CONTRIBUTION :
C      ---> TO AB(I, I),AB(I, J)  :  IF  I=II  OR  I=IJ
C      ---> TO AB(I,II),AB(I,IJ)  :  IF  I=IK  OR  I=IL
C      ---> = 0.                  :  IF  I . NE . II,IJ,IK,IL
C     ACCORDING TO THE LINEARISATION OF WAGONER (1969)
C.......................................................................
      DO 600 I=1,NSP
      DO 610 J=1,NSP
      AT(J,I)=AIN
 610  CONTINUE
 600  CONTINUE
      DO 50 M=1,NRE
      II   = K1(M)
      NII= K2(M)
      IJ   = K3(M)
      NIJ= K4(M)
      IK   = K5(M)
      NIK= K6(M)
      IL   = K7(M)
      NIL= K8(M)
      AMT=V(M) * DELT / (DF(NII)*DF(NIJ)*(NII+NIJ))
      NI1=NII-1
      NJ1=NIJ-1
      AM=AMT*NII
      AT(II,II)=AT(II,II)+AM*NII*((Y(II)+O)**NI1)*((Y(IJ)+O)**NIJ)
      AT(IJ,II)=AT(IJ,II)+AM*NIJ*((Y(IJ)+O)**NJ1)*((Y(II)+O)**NII)
      AM=AMT*NIJ
      AT(IJ,IJ)=AT(IJ,IJ)+AM*NIJ*((Y(IJ)+O)**NJ1)*((Y(II)+O)**NII)
      AT(II,IJ)=AT(II,IJ)+AM*NII*((Y(II)+O)**NI1)*((Y(IJ)+O)**NIJ)
      AM=AMT*NIK
      AT(II,IK)=AT(II,IK)-AM*NII*((Y(II)+O)**NI1)*((Y(IJ)+O)**NIJ)
      AT(IJ,IK)=AT(IJ,IK)-AM*NIJ*((Y(IJ)+O)**NJ1)*((Y(II)+O)**NII)
      AM=AMT*NIL
      AT(II,IL)=AT(II,IL)-AM*NII*((Y(II)+O)**NI1)*((Y(IJ)+O)**NIJ)
      AT(IJ,IL)=AT(IJ,IL)-AM*NIJ*((Y(IJ)+O)**NJ1)*((Y(II)+O)**NII)
   50 CONTINUE
C
      DO 150 I=1,NIS
      AT(I,I)=AT(I,I)+1.D0
 150  CONTINUE
C
      RETURN
      END
C======================================================================
      SUBROUTINE EIGEN(C,X,A,DIV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C-----------------------------------------------------------------------
CC    INVERSION OF THE MATRIX (SPECIAL SPARSE FORM)
C-----------------------------------------------------------------------
      PARAMETER (NIS=14, NSP=NIS+1, NRE=29)
      PARAMETER(NDS=1,NDI=1,IJD=4,IDEL=IJD-1,JDEL=IJD-1)
C  Warning ! : NDS & NDI are, respectively, the number of diagonals
C  containing non-zero elements, above and below the main diagonal.
C  IJD is the dimension of the upper left square matrix + 1.
      COMMON/CNETW/AN(NSP),ZN(NSP),Q(NRE),BE(NSP),BE0
      COMMON/AVHD/ATH(NSP,IDEL),ATV(JDEL,NSP),ATD(NDI+NDS+1,NSP)
      DIMENSION C(NSP),X(NSP),DIV(NSP),A(NSP,NSP),
     *CT(NSP),XSOM(NSP),C0(NSP),CTT(NSP)
      N=NSP-1
      SOMME=0.D0
      DO 6543 I=1,N
      XSOM(I)=C(I)
 6543 SOMME=SOMME+XSOM(I)*AN(I)
      DO 6544 I=1,N
 6544 C(I)=C(I)/SOMME
      DO 703 J=1,N
      CTT(J)=C(J)
703   C0(J)=0.D0
      DO 734 J=1,IJD-1
      DO 724 I=1,N
      ATH(I,J)=A(I,J)
724   C0(J)=C0(J)+A(I,J)*C(I)
734   C0(J)=C(J)-C0(J)
      DO 825 J=1,IJD-1
825   CONTINUE
      DO 735 J=IJD,N
      DO 725 I=1,IDEL
      ATV(I,J)=A(I,J)
725   C0(J)=C0(J)+A(I,J)*C(I)
      IMIN=J-NDI
      IMAX=J+NDS
      IF(IMIN.LT.IJD)IMIN=IJD
      IF(IMAX.GT.N)IMAX=N
      DO 726 I=IMIN,IMAX
      K=I-IMIN+1
      ATD(K,J)=A(I,J)
726   C0(J)=C0(J)+A(I,J)*C(I)
735   C0(J)=C(J)-C0(J)
      DO 866 J=IJD,N
866   CONTINUE
      DO 737 J=1,N
      C(J)=C0(J)
      CT(J)=C(J)
737   XSOM(J)=0.D0
C====> C A R E F U L L:  ITER=1,2
      DO 10000 ITER=1,2
C
C---BEGINNING OF THE GAUSSIAN ELIMINATION PROCEDURE
C-(1)-ELIMINATION OF THE LOWER DIAGONALS
C
      DO 1000 JBAL=IJD,N-1
      DIV(JBAL)=-1.D0/A(JBAL,JBAL)
      JMAX=JBAL+NDI
      IF(JMAX.GT.N)JMAX=N
      IMAX=JBAL+NDS
      IF(IMAX.GT.N)IMAX=N
      DO 1000 J=JBAL+1,JMAX
      IF(A(JBAL,J).EQ.0.D0)GOTO 1000
      DIVJ=DIV(JBAL)*A(JBAL,J)
      DO 10 I=1,IDEL
10    A(I,J)=DIVJ*A(I,JBAL)+A(I,J)
      DO 20 I=JBAL+1,IMAX
20    A(I,J)=DIVJ*A(I,JBAL)+A(I,J)
      C(J)=DIVJ*C(JBAL)+C(J)
1000  CONTINUE
      DIV(N)=-1.D0/A(N,N)
C
C-(2)-ELIMINATION OF THE UPPER DIAGONALS AND OF THE HORIZONTAL BAND
C
      DO 2000 JBAL=N,IJD+1,-1
      JMIN=JBAL-NDI
      IF(JMIN.LT.IJD)JMIN=IJD
      DO 200 J=JMIN,JBAL-1
      IF(A(JBAL,J).EQ.0.)GOTO 200
      DIVJ=DIV(JBAL)*A(JBAL,J)
      DO 30 I=1,IDEL
30    A(I,J)=DIVJ*A(I,JBAL)+A(I,J)
      C(J)=DIVJ*C(JBAL)+C(J)
200   CONTINUE
      DO 300 J=1,JDEL
      IF(A(JBAL,J).EQ.0.)GOTO 300
      DIVJ=DIV(JBAL)*A(JBAL,J)
      DO 40 I=1,IDEL
40    A(I,J)=DIVJ*A(I,JBAL)+A(I,J)
      C(J)=DIVJ*C(JBAL)+C(J)
300   CONTINUE
2000  CONTINUE
      DO 400 J=1,JDEL
      IF(A(IJD,J).EQ.0.D0)GOTO 400
      DIVJ=DIV(IJD)*A(IJD,J)
      DO 50 I=1,IDEL
50    A(I,J)=DIVJ*A(I,IJD)+A(I,J)
      C(J)=DIVJ*C(IJD)+C(J)
400   CONTINUE
C
C-(3)-GAUSSIAN ELIMINATION OF THE UPPER LEFT SQUARE MATRIX
C
      DO 3000 JBAL=1,IJD-2
      DIV(JBAL)=-1.D0/A(JBAL,JBAL)
      DO 3000 J=JBAL+1,IJD-1
      IF(A(JBAL,J).EQ.0.D0)GOTO 3000
      DIVJ=DIV(JBAL)*A(JBAL,J)
      DO 60 I=JBAL+1,IJD-1
60    A(I,J)=DIVJ*A(I,JBAL)+A(I,J)
      C(J)=DIVJ*C(JBAL)+C(J)
3000  CONTINUE
      DIV(IJD-1)=-1.D0/A(IJD-1,IJD-1)
      X(IJD-1)=-DIV(IJD-1)*C(IJD-1)
      DO 4000 JBAL=IJD-2,1,-1
      SOM=0.D0
      DO 70 I=IJD-1,JBAL+1,-1
  70  SOM=SOM+A(I,JBAL)*X(I)
4000  X(JBAL)=DIV(JBAL)*(SOM-C(JBAL))
      SOM=0.D0
      DO 80 I=1,IDEL
80    SOM=SOM+A(I,N)*X(I)
      X(N)=DIV(N)*(SOM-C(N))
      DO 5000 JBAL=N-1,IJD,-1
      SOM=0.D0
      DO 90 I=1,IDEL
90    SOM=SOM+A(I,JBAL)*X(I)
5000  X(JBAL)=DIV(JBAL)*(SOM-C(JBAL))
      DO 202 J=1,N
      XSOM(J)=XSOM(J)+X(J)
202   C(J)=0.D0
      IF(ITER.EQ.2) GO TO 10000
      DO 233 J=1,IJD-1
      DO 234 I=1,N
      A(I,J)=ATH(I,J)
234   C(J)=C(J)+A(I,J)*X(I)
233   C(J)=CT(J)-C(J)
      DO 235 J=IJD,N
      DO 236 I=1,IDEL
      A(I,J)=ATV(I,J)
236   C(J)=C(J)+A(I,J)*X(I)
      IMIN=J-NDI
      IMAX=J+NDS
      IF(IMIN.LT.IJD)IMIN=IJD
      IF(IMAX.GT.N)IMAX=N
      DO 237 I=IMIN,IMAX
      K=I-IMIN+1
      A(I,J)=ATD(K,J)
237   C(J)=C(J)+A(I,J)*X(I)
235   C(J)=CT(J)-C(J)
10000 CONTINUE
      DO 950 I=1,N
950   X(I)=XSOM(I)+CTT(I)
      RETURN
      END
C======================================================================
      SUBROUTINE STEST(YN,YO,PROV,KMM,DDDD,TAU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C-----------------------------------------------------------------------
C     CHOICE OF THE LARGEST ABUNDANCE VARIATION DDDD= DY/Y = (YN-Y0)/Y0
C     OR,EQUIVALENTLY,OF THE SMALLEST PRO=Y(I)/DELY(I).
C     IMPORTANT : ONLY SIGNIFICANT ISOTOPES (Y.GT.YTMIN) CONSIDERED.
C ----------------------------------------------------------------------
      PARAMETER (NIS=14,NSP=NIS+1)
      COMMON/CPARM/PAR,PAMIN,DYMIN,DEMIN,YMIN,YTMIN,DEPS,ATEST,TES
      COMMON/IHYDRO/INEUT,IPROT,IDEUT,ITRIT
      DIMENSION YN(NSP),YO(NSP),TEST(NSP)
C===> SPECIAL CARE FOR I=1,2,3 (N,P,A): FOLLOWED DOWN TO YLMIN
C     YLMIN=5.D-15
      DO 10 I=1,NIS
      if( (I.EQ.INEUT).OR.(I.EQ.IDEUT).OR.(I.EQ.ITRIT)) go to 8
      DELY=YN(I)-YO(I)
      IF(ABS(DELY).LT.DEMIN) DELY=DEMIN
   4  IF(ABS(YO(I))-YTMIN) 8,7,7
   7  TEST1= ABS(YN(I)/DELY)
      TEST2= ABS(YO(I)/DELY)
      TEST(I)=MIN(TEST1,TEST2)
      GO TO 10
   8  TEST(I)=TES
  10  CONTINUE
      PRO=TEST(1)
      KMM=1
      DO 20 I=1,NIS
      IF(PRO-TEST(I)) 20,20,15
  15  PRO=TEST(I)
      KMM=I
  20  CONTINUE
      PROV=PRO
      YYYY=YO(KMM)
      IF(YYYY.EQ.0.D0) GO TO 999
      DDDD=(YN(KMM)-YYYY)/YYYY
      RETURN
  999 DDDD=0.D0
C
      RETURN
      END

C======================================================================
      subroutine genpar(rho,t9,a,z,y,ye,ar,ar11,coefg,coeft)
c*************************************************************
c                                                            *
c  this subroutine computes the screening parameters for the *
c  Thielemann network.                                       *
c                                                            *
c*************************************************************
      implicit none

      integer i,NIS,NSP
      PARAMETER(NIS=14,NSP=NIS+1)
      real*8 y(NSP),a(NSP),z(NSP),ar(NIS),rho,t9,ye,ar11,coefg,coeft,
     $     coefa


      coefa=7.345889d-9*(rho*ye)**(-0.33333333333d0) 
      do 30 i=1,NIS
      ar(i)=coefa*z(i)**(0.3333333333d0)
   30 continue
      ar11=coefa*4.0d0**(0.3333333333d0)
      coefg=1.67100d-12/t9
      coeft=3.58039d4*coefg**(0.3333333333333d0)

      return
      end


      real*8 function scrng(z1,z2,a1,a2,ar1,ar2,coefg,coeft)

      implicit none
      real*8 z1,z2,a1,a2,ar1,ar2,coefg,coeft,gama,
     &      tauu,betta,fun

      gama=coefg*z1*z2*2.0d0/(ar1+ar2)
      tauu=coeft*(a1*a2/(a1+a2)*z1**2*z2**2)**(0.33333333333d0)

      betta=3.0d0*gama/tauu
      fun=0.0455d0*betta+0.348d0*betta**3+9.49d0*betta**6-
     &  0.123d0*betta**12+0.101d0*betta**13
      fun=fun/(1.0d0+100.0d0*betta**4+0.267d0*betta**12)
      scrng=dexp(1.25d0*gama-tauu*fun)

c      scrng=dexp(1.25d0*gamma-0.0975d0*tau*(3.0d0*gamma/tau)**2) 

      return
      end

C========================================================================
      subroutine vit(T,rho,pme,iread)

c     This routine provides the reaction rates v(i) (i=1,nreac) 
c     for temperature T and density RHO 
c     from a linear interpolation of the log of [the grid-point rates read
c     from the file 'vit.dat' containing the output from Netgen].
c
c     Note that the factorials accounting for identical particles have
c     already been included in v(i).
c     To obtain the evolution dYj/dt of species j, simply multiply v(i) by
c     the molar fractions of the reacting species j1 + j2: Yj1^n1 Yj2^n2,
c     where n1, n2 are the stoechiometric factors (stored in the vectors
c     with the same name).
c
c     The last lines of subroutine vit print all variables to standard output,
c     to check that it worked properly
c
c     The first call to vit should be done with 
c     iread = 0: to read the vit.dat file, and then store them in array vgrid
c
c     Subsequent calls may be done with
c     iread = 1: to compute rates at other temperatures, after vgrid
c                has been initialized by the first call with iread = 0     
 
c     The parameters 
c          NGRID = maximum number of temperature grid points 
c          NRE   = maximum number of reactions
c
c     may be lowered to better match your network.
c
c     The actual number of grid points and reaction rates are computed
c     automatically, and stored in variables IREAC and KGRID
c
c     PME = mean molecular weight per electron
c     T = temperature in K
c     RHO = density in g/cm3
c
c
c     The rates on grid points are stored in the array 
c         vgrid(ngrid,nre,11)
c     where the last index refers to the density grid for the weak rates
c     dependent upon electron number density:
c     vgrid(ngrid,nre,i) with i=1,11 
c             corresponds to the rate at different electron densities Ne 
c             as given by the variable aNegrid(nre,i)
c     For the other reactions (with flags < 0), 
c     only vgrid(ngrid,nre,1) is relevant 
c
          
c
c     Coding of reactions:
c
c       n1 to n4: stoechiometric factors
c       z1 to z4: chemical symbol
c       a1 to a4: atomic mass

c     The following table allows to convert proton number into
c     chemical symbol, if required:

c     CHARACTER*2 SYM(0:105)*2
c     DATA SYM/'n ',
c    & 'H ','HE','LI','BE','B ','C ','N ','O ','F ','NE',
c    1 'NA','MG','AL','SI','P ','S ','CL','AR','K ','CA',
c    2 'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',
c    3 'GA','GE','AS','SE','BR','KR','RB','SR','Y ','ZR', 
c    4 'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN',
c    5 'SB','TE','I ','XE','CS','BA','LA','CE','PR','ND',
c    6 'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB',
c    7 'LU','HF','TA','W ','RE','OS','IR','PT','AU','HG',
c    8 'TL','PB','BI','PO','AT','RN','FR','RA','AC','TH',
c    9 'PA','U ','NP','PU','AM','CM','BK','CF','ES','FM',
c    & 'MD','NO','LR','RF','HA'/
c
c     ====================================================================
c     Duplicate the following lines in the main program
c     
      implicit double precision (a-h,o-z)

      parameter(ngrid=60,NRE=29) 

      common/cvit/v(NRE),tgrid(ngrid),aNegrid(nre,11)
      common/cvgrid/vgrid(ngrid,nre,11),flag(nre),ndata(nre)
      common/network/ireac,kgrid,n1(nre),n2(nre),n3(nre),n4(nre),
     &               z1(nre),z2(nre),z3(nre),z4(nre),
     &               a1(nre),a2(nre),a3(nre),a4(nre),
     &               reaction(nre)
      common/Q/Qrad(nre),Qnu(nre) 

      character z1*6,z2*6,z3*6,z4*6
      character reaction*37

c     ======================================================================
      character longline(7)*132,dummy*1
      character aflag(10)*8
      dimension nn1(11),nn2(11),nn3(11),nn4(11)
      character zz1(11)*6,zz2(11)*6,zz3(11)*6,zz4(11)*6
      dimension QQrad(11),QQnu(11)
      dimension f(0:4)

c     factorials
C **
C **  Els factorials son tinguts en compte fora d'aquesta subrutina
C **
c      f(0) = 1.
c      f(1) = 1.
c      f(2) = 2.
c      f(3) = 6.
c      f(4) = 24.

 

      T8 = T * 1.E-08

c     ane is the electron density    [g / cm^3]
c     PME is the mean molecular weight per electron 
c         =  [ SUM (XZ/A) ]**(-1)    [Adimensional]
c     1amu is the atomic mass unit   [1 amu = 1.66053 E-24 g/particle]
c     1 / 1amu = 6.02E+23 particles/g
c     ane = rho / (PME * 1amu)       [particle / cm^3]

      ane   = rho * 6.02E+23 / PME 

c     if iread = 0, start by reading the data file, then compute the rates

      if(iread.eq.0) then

         OPEN(unit=10,file='VIT.HELI10b')
c        OPEN(unit=10,file='VIT.NOVAbr1')

c        read the data table 

         ireac = 0

         do ll=1,50000


c           initialize the string
            do i = 1,7
               do j=1,132
                  longline(i)(j:j)=' '
               enddo
            enddo

c           read the header
            read(10,1100,end=9999) (longline(j),j=1,7)
            ireac = ireac + 1

c           l = number of data records on line
            leng = index(longline(5),'         ') 
            if(leng.eq.0) leng = 132
            l = (leng - 11)/11
            ndata(ireac) = l


            read(longline(7),1200) (aflag(j),j=1,l)
            aNegrid(ireac,1) = 1.

c           l > 1: flag corresponds to electron densities
c                  read electron densities
            if(l.gt.1) then
                 read(longline(7),1210) (aNegrid(ireac,j),j=1,l)
                 flag(ireac) = aNegrid(ireac,1)
            endif

            read(longline(1),1201,end=1) (nn1(j),j=1,l)
 1          read(longline(1),1202,end=2) (zz1(j),j=1,l)
 2          read(longline(2),1201,end=3) (nn2(j),j=1,l)
 3          read(longline(2),1202,end=4) (zz2(j),j=1,l)
 4          read(longline(3),1201,end=5) (nn3(j),j=1,l)
 5          read(longline(3),1202,end=6) (zz3(j),j=1,l)
 6          read(longline(4),1201,end=7) (nn4(j),j=1,l)
 7          read(longline(4),1202,end=8) (zz4(j),j=1,l)
 8          read(longline(5),1205,end=11) (QQrad(j),j=1,l)
 11         read(longline(6),1205,err=20,end=12) (QQnu(j), j=1,l)
 12         continue

c           read the data

            read(10,*)
            read(10,*)

c           check the grid size
            if(ireac.eq.1) then 

               do i=1,ngrid
                  read(10,1206,err=21) 
     &                          dummy,tgrid(i),(vgrid(i,ireac,j),j=1,l)
                  if(dummy(1:1) .eq. '#') goto 13
                  kgrid = i
               enddo

            else 

               do i=1,kgrid
                  read(10,1206,end=13,err=21) dummy,
     &                                 tgrid(i),(vgrid(i,ireac,j),j=1,l)
               enddo
            endif

            read(10,*)
 13         continue


               n1(ireac) = nn1(1)
               n2(ireac) = nn2(1)
               n3(ireac) = nn3(1)
               n4(ireac) = nn4(1)
               z1(ireac) = zz1(1)
               z2(ireac) = zz2(1)
               z3(ireac) = zz3(1)
               z4(ireac) = zz4(1)
               Qrad(ireac)=QQrad(1)
               Qnu(ireac) =QQnu(1)

               if(aflag(1)(5:8).eq.'----') then
                   flag(ireac) = -200.
               else if(aflag(1)(6:8).eq.'---') then
                   flag(ireac) = -100.
               else if(aflag(1)(7:8).eq.'--') then
                   flag(ireac) = -10.
               else if(aflag(1)(7:8).eq.'+-') then
                   flag(ireac) = -27.
               else if(aflag(1)(6:8).eq.'+++') then
                  flag(ireac) = -14.
               else if(aflag(1)(7:8).eq.'++') then
                   flag(ireac) = -13.
               else if(aflag(1)(8:8).eq.'+') then
                   flag(ireac) = -11.
               endif        

                if(zz1(1)(1:4).eq.'NEUT'.or.zz1(1)(1:4).eq.'PROT') then 
                   a1(ireac) = 1.
                else if(zz1(1)(1:5).eq.'OOOOO') then
                   a1(ireac) = 0.
                else if(zz1(1)(1:4).eq.'DEUT') then
                   a1(ireac) = 2.
                else if(zz1(1)(1:4).eq.'TRIT') then
                   a1(ireac) = 3.
                else if(zz1(1)(1:5).eq.'AL 26') then
                   a1(ireac) = 26.
                else 
                   read(zz1(1)(3:5),'(i3)') i1
                   a1(ireac) = float(i1)
                endif
                
                if(zz2(1)(1:4).eq.'NEUT'.or.zz2(1)(1:4).eq.'PROT') then 
                   a2(ireac) = 1.
                else if(zz2(1)(1:5).eq.'OOOOO') then
                   a2(ireac) = 0.
                else if(zz2(1)(1:4).eq.'DEUT') then
                   a2(ireac) = 2.
                else if(zz2(1)(1:4).eq.'TRIT') then
                   a2(ireac) = 3.
                else if(zz2(1)(1:5).eq.'HE  4') then
                   a2(ireac) = 4.
                else if(zz2(1)(1:5).eq.'HE  3') then
                   a2(ireac) = 3.
                endif
                
                if(zz3(1)(1:4).eq.'NEUT'.or.zz3(1)(1:4).eq.'PROT') then 
                   a3(ireac) = 1.
                else if(zz3(1)(1:5).eq.'OOOOO') then
                   a3(ireac) = 0.
                else if(zz3(1)(1:4).eq.'DEUT') then
                   a3(ireac) = 2.
                else if(zz3(1)(1:4).eq.'TRIT') then
                   a3(ireac) = 3.
                else if(zz3(1)(1:5).eq.'HE  4') then
                   a3(ireac) = 4.
                else if(zz3(1)(1:5).eq.'HE  3') then
                   a3(ireac) = 3.
                endif
                
                if(zz4(1)(1:4).eq.'NEUT'.or.zz4(1)(1:4).eq.'PROT') then 
                   a4(ireac) = 1.
                else if(zz4(1)(1:5).eq.'OOOOO') then
                   a4(ireac) = 0.
                else if(zz4(1)(1:4).eq.'DEUT') then
                   a4(ireac) = 2.
                else if(zz4(1)(1:4).eq.'TRIT') then
                   a4(ireac) = 3.
                else if(zz4(1)(1:5).eq.'AL 26') then
                   a4(ireac) = 26.
                else 
                   read(zz4(1)(3:5),'(i3)') i4
                   a4(ireac) = float(i4)
                endif

                write(reaction(ireac),1204) 
     &                    nn1(1),zz1(1)(1:6),nn2(1),zz2(1)(1:5),
     &                    nn3(1),zz3(1)(1:5),nn4(1),zz4(1)(1:6)

           enddo


 1100    format(4(a132,//),3(a132,/))
 1101    format (i3,7(a132,/))
 1200    format(14x,11(a8,3x))
 1201    format(14x,11(i1,10x))
 1202    format(14x,11(2x,a6,3x))
 1203    format(14x,11(4x,i3,4x))
 1204    format(i1,1x,a6,'( ',i1,1x,a5,',',1x,i1,1x,a5,')',2x,i1,1x,a6)
 1205    format(14x,11(f8.3,3x))
 1206    format(a1,2x,f8.4,11(1x,e10.4))
 1210    format(14x,11(e8.1,3x))


         endif

 9999    continue
c        end of iread=0 loop
         close(10)


c            DO JJ=1,ngrid
c            WRITE(7,*) 'Escrivim algunes reaccions'
c            WRITE(7,*) 'Reaction 1'
c            WRITE(7,*) tgrid(JJ),vgrid(JJ,1,1)
c            WRITE(7,*) 'Reaction 10'
c            WRITE(7,*) tgrid(JJ),vgrid(JJ,10,1)
c            WRITE(7,*) 'Reaction 100'
c            WRITE(7,*) tgrid(JJ),vgrid(JJ,100,1)
c1714        CONTINUE

c        interpolate reaction rate

c     locate the position of current T within the grid
c     and compute step, a, b, c and d parameters for cubic spline interpolation

      klo = 1
      khi = kgrid
 100  if(khi-klo.gt.1) then
	k = (khi+klo)/2
	if(tgrid(k).gt.T8) then
	  khi=k
        else
	  klo=k
        endif
      goto 100
      endif
      h=tgrid(khi)-tgrid(klo)
      b=(T8 - tgrid(klo))/h    ! [b = (T - Tlow) / (Thi - Tlow)]


c     perform linear interpolation of log(v)


      do i=1,ireac

       if(vgrid(klo,i,1).GT.0.D0 .and. vgrid(khi,i,1).gt.0.D0) then 

	if(flag(i).le.0.) then 
c          reaction is not electron-density-dependent beta-decay

          v(i)= log10(vgrid(klo,i,1))  
     &        + b * (log10(vgrid(khi,i,1)) - log10(vgrid(klo,i,1))) 


c            test the reaction kind

             if (flag(i).eq.-14.) then
c               electron capture on Be7
	        v(i) = 10.**v(i) * RHO / PME
             else if (flag(i).eq.-13.) then
c               electron capture
	        v(i) = 10.**v(i) * RHO / PME
             else if (flag(i).eq.-27.) then
c               electron capture pep (special case)
	        v(i) = 10.**v(i) * RHO *RHO / PME
             else if (flag(i).eq.-11.) then
c               photodisintegration or beta-decay
                v(i) = 10.**v(i)
             else if (flag(i).eq.-10.) then
c               two-particle reaction
c               !if identical particles: factorials! 
c	        v(i) = 10.**v(i) * RHO /f(n1(i))/f(n2(i))
	        v(i) = 10.**v(i) * RHO  !Es divideix fora
	     else if (flag(i).eq.-100.) then
c               three-particle reactions
c               !if identical particles: factorials! 
c		v(i) = 10.**v(i) * RHO * RHO  /f(n1(i))/f(n2(i))
		v(i) = 10.**v(i) * RHO * RHO !Es divideix fora
	     else if (flag(i).eq.-200.) then
c               four-particle reactions
c               !if identical particles: factorials! 
c		v(i) = 10.**v(i) * RHO * RHO * RHO / f(n1(i))/f(n2(i))
		v(i) = 10.**v(i) * RHO * RHO * RHO !Es divideix fora
	     endif

	  else

c         beta-decay rate has to be interpolated in density

c         locate the position of current rho within the density grid aNegrid


            jlo = 1
            jhi = ndata(i)
 101        if(jhi-jlo.gt.1) then
	       j = (jhi+jlo)/2
	       if(aNegrid(i,j).gt.ane) then
	          jhi=j
               else
	          jlo=j
               endif
               goto 101
            endif
            deltagrid=aNegrid(i,jhi)-aNegrid(i,jlo)

c         linear interpolation in temperature for two density grid points
c                (note: log of beta decay rate is handled)    

            vbetalow    = log10(vgrid(klo,i,jlo))    
     &                  + b * 
     &                   (log10(vgrid(khi,i,jlo))  
     &                  - log10(vgrid(klo,i,jlo))) 
            vbetahigh   = log10(vgrid(klo,i,jhi))    
     &                  + b * 
     &                   (log10(vgrid(khi,i,jhi)) 
     &                  - log10(vgrid(klo,i,jhi))) 

            vbetalow    = 10.**vbetalow
	    vbetahigh   = 10.**vbetahigh

c            write(6,*) vbetalow,vbetahigh,jlo,jhi,ndata(i)

c           linear interpolation in density 

            v(i) = (vbetahigh-vbetalow)/deltagrid
     &              *(ane-aNegrid(i,jhi))
     &           +  vbetahigh

c           extrapolation often leads to negative rates 
c           check if so, and then put it to the closest density grid point
c                                     or to zero
c           (supplementary table values would be needed for those cases)

            if(v(i).lt.0.) then
               if(ane-aNegrid(i,jhi).lt.0.) then
c                 extrapolation at lower end of table
c                  v(i) = 0.
                 v(i) = vbetalow
               else
c                 extrapolation at upper end of table
c                  v(i) = 0.
                 v(i) = vbetahigh
               endif    
            endif

         endif

        else
           v(i) = 0.d0
        endif

       enddo

c        write reaction rates to standard output

c         do i = 1,ireac
c            write(22,6000) reaction(i),flag(i)
c            write(22,6002) n1(i),a1(i),z1(i),
c     &                 n2(i),a2(i),z2(i),
c     &                 n3(i),a3(i),z3(i),
c     &                 n4(i),a4(i),z4(i)
c            write(22,6001) Qrad(i),Qnu(i)
c            do l=1,ndata(i)
c               write(22,6003) (vgrid(k,i,l),k=1,kgrid)
c            enddo
c            write(22,6004) T8, rho,ane,v(i)
c         enddo


 6000    format(1x,a37,1x,'flag= ',f5.0)
 6001    format(1x,'Qrad (MeV) = ',f7.3,1x,'Qnu (MeV) = ',f7.3)
 6002    format(1x,'i1 = ',i1,1x,'a1 = ',f4.0,1x,'z1 = ',a6,/,
     &          1x,'i2 = ',i1,1x,'a2 = ',f4.0,1x,'z2 = ',a6,/,
     &          1x,'i3 = ',i1,1x,'a3 = ',f4.0,1x,'z3 = ',a6,/,
     &          1x,'i4 = ',i1,1x,'a4 = ',f4.0,1x,'z4 = ',a6)

 6004    format(' T8 = ',f5.2,1x,'rho =',e8.3,1x,
     &          'electron density =',e7.2,'cm-3',/,
     &          ' rate = ',e9.4,/)

 6003    format(11(1x,e9.4))


         return

c        branch here if NaN was detected in Qnu line
 20      continue

         write(6,7000) longline(1)(15:21),longline(2)(15:21),
     &                 longline(3)(15:21),longline(4)(15:21)

 7000    format('** WARNING: NaN have been detected among Qnu values',
     &          ' for reaction ',a7,' + ',a7,' = ',a7,' + ',a7,/,
     &          ' This may be because Qnu is density-dependent.',
     &          ' The Netgen log file provides a complete table',
     &          ' of Qnu(T,rho).',/,' Edit the vit.dat file to remove',
     &          ' NaNs and try again')  

         stop

c        branch here if NaN was detected among V
 21      continue

         write(6,7001) longline(1)(15:21),longline(2)(15:21),
     &                 longline(3)(15:21),longline(4)(15:21)

 7001    format('** WARNING: NaN has been detected in reaction rate',
     &          ' for reaction ',a7,' + ',a7,' = ',a7,' + ',a7,/,
     &          'This may be because your grid extends beyond the',
     &          ' valid data range for this reaction',
     &          ' (see Netgen log file)',/,
     &          ' Either change your temperature or density grid,',
     &          ' or remove NaN in vit.dat by extrapolating the',
     &          ' existing data (risky!)')  

         return
         end

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
C***********************************************************************
      SUBROUTINE RPARAM
C***********************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/CPARM/PAR,PAMIN,DYMIN,DEMIN,YMIN,YTMIN,DEPS,ATEST,TES
      COMMON/CFACT/DF(0:6)
      PAR=.15D0
      PAMIN=5.D0
      ATEST=5.D0
      TES=20.D0
      DYMIN=1.D-40
      DEMIN=1.D-40
      YMIN=1.D-35
      YTMIN=1.D-8
C--->  Test : lower YTMIN down to 10(-12)
      YEMINM=YTMIN
      DEPS=1.D-6
C--->  Test : lower DEPS down to 10(-10)
      DF(0)=1.D0
      DF(1)=1.D0
      DF(2)=2.D0
      DF(3)=6.D0
      DF(4)=24.D0
      DF(5)=120.D0
      DF(6)=720.D0
      RETURN
      END

