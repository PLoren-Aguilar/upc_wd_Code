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

