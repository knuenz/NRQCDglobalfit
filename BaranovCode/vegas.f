      SUBROUTINE VEGAS(FXN,AVGI,SD,CHI2A)
C====================================================================
C   SUBROUTINE PERFORMS N-DIMENSIONAL MONTE CARLO INTEG'N
C      - BY G.P. LEPAGE   SEPT 1976/(REV)APR 1978
C====================================================================
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      EXTERNAL FXN
      DIMENSION QRAN(10)
      COMMON/BVEG1/XL(10),XU(10),ACC
      COMMON/BVEGG/NDIM,NCALL,ITMX,NPRN
      COMMON/BVEG2/XI(50,10),SI,SI2,SWGT,SCHI
      COMMON/BVE22/NDO,IT
      SAVE
      DIMENSION D(50,10),DI(50,10),XIN(50),R(50),DX(10),DT(10),X(10)
     1   ,KG(10),IA(10)
      DATA NDMX/50/,ALPH/1.5D0/,QNE/1.D0/,MDS/1/
C
      NDO=1
      DO 1 J=1,NDIM
C      Print *, "vegas1"
 1    XI(1,J)=QNE
C
C      Print *, "vegas1a"
      ENTRY VEGAS1(FXN,AVGI,SD,CHI2A)
C         - INITIALIZES CUMMULATIVE VARIABLES, BUT NOT GRID
C      Print *, "vegas1b"
      IT=0
      SI=0.D0
      SI2=SI
      SWGT=SI
      SCHI=SI
C
C      Print *, "vegas1c"     
      ENTRY VEGAS2(FXN,AVGI,SD,CHI2A)
C         - NO INITIALIZATION
C     Print *, "vegas1d"     
      ND=NDMX
      NG=1
      IF(MDS.EQ.0) GO TO 2
      NG=INT((DFLOAT(NCALL)/2.D0)**(1.D0/DFLOAT(NDIM)))
      MDS=1
      IF((2*NG-NDMX).LT.0) GO TO 2
      MDS=-1
      NPG=NG/NDMX+1
      ND=NG/NPG
      NG=NPG*ND
C      Print *, "vegas2"
 2    K=NG**NDIM
C     Print *, "vegas2a"     
      NPG=NCALL/K
      IF(NPG.LT.2) NPG=2
      CALLS=NPG*K
      DXG=QNE/NG
      DV2G=(CALLS*DXG**NDIM)**2/NPG/NPG/(NPG-QNE)
      XND=ND
      NDM=ND-1
      DXG=DXG*XND
      XJAC=QNE/CALLS
      DO 3 J=1,NDIM
      DX(J)=XU(J)-XL(J)
 3    XJAC=XJAC*DX(J)
C     Print *, "vegas3"     
C
C   REBIN PRESERVING BIN DENSITY
C
      IF(ND.EQ.NDO) GO TO 8
      RC=NDO/XND
      DO 7 J=1,NDIM
      K=0
      XN=0.D0
      DR=XN
      I=K
 4    K=K+1
C     Print *, "vegas4"     
      DR=DR+QNE
      XO=XN
      XN=XI(K,J)
 5    IF(RC.GT.DR) GO TO 4
C     Print *, "vegas5"     
      I=I+1
      DR=DR-RC
      XIN(I)=XN-(XN-XO)*DR
      IF(I.LT.NDM) GO TO 5
      DO 6 I=1,NDM
 6    XI(I,J)=XIN(I)
 7    XI(ND,J)=QNE
      NDO=ND
C
 8    IF(NPRN.NE.0) WRITE(6,200) NDIM,CALLS,IT,ITMX,ACC,MDS,ND
     1                           ,(XL(J),XU(J),J=1,NDIM)
C
      ENTRY VEGAS3(FXN,AVGI,SD,CHI2A)
C         - MAIN INTEGRATION LOOP
 9    IT=IT+1
C     Print *, "vegas9"
      TI=0.D0
      TSI=TI
      DO 10 J=1,NDIM
      KG(J)=1
      DO 10 I=1,ND
      D(I,J)=TI
 10   DI(I,J)=TI
C
 11   FB=0.D0
      F2B=FB
      K=0
 12   K=K+1
      CALL ARAN9(QRAN,NDIM)
      WGT=XJAC
      DO 15 J=1,NDIM
      XN=(DFLOAT(KG(J))-QRAN(J))*DXG+QNE
      IA(J)=XN
      IF(IA(J).GT.1) GO TO 13
      XO=XI(IA(J),J)
      RC=(XN-IA(J))*XO
      GO TO 14
 13   XO=XI(IA(J),J)-XI(IA(J)-1,J)
      RC=XI(IA(J)-1,J)+(XN-IA(J))*XO
 14   X(J)=XL(J)+RC*DX(J)
 15   WGT=WGT*XO*XND
C     Print *, "vegas15"
C
      F=WGT
C     Print *, "vegas15a"
      F=F*FXN(X,WGT)
      F2=F*F
      FB=FB+F
C     Print *, "vegas15b"
      F2B=F2B+F2
      DO 16 J=1,NDIM
C     Print *, "vegas15c"
      DI(IA(J),J)=DI(IA(J),J)+F
C     Print *, "vegas15z"
 16   IF(MDS.GE.0) D(IA(J),J)=D(IA(J),J)+F2
      IF(K.LT.NPG) GO TO 12
C
      IF(FB.lt.1.d-20) FB = 0.D0
C      write(6,*) 'TI = ',ti,tsi,fb,fb2

      F2B=DSQRT(F2B*NPG)
      F2B=(F2B-FB)*(F2B+FB)
      TI=TI+FB
      TSI=TSI+F2B
      IF(MDS.GE.0) GO TO 18
      DO 17 J=1,NDIM
 17   D(IA(J),J)=D(IA(J),J)+F2B
C     Print *, "vegas17"
 18   K=NDIM
 19   KG(K)=MOD(KG(K),NG)+1
C     Print *, "vegas19"
      IF(KG(K).NE.1) GO TO 11
      K=K-1
      IF(K.GT.0) GO TO 19
C
C   FINAL RESULTS FOR THIS ITERATION
C
      TSI=TSI*DV2G
      TI2=TI*TI
      WGT=TI2/TSI
c      write(6,*) 'wgt ',wgt,ti2,tsi,ti
      SI=SI+TI*WGT
      SI2=SI2+TI2
      SWGT=SWGT+WGT
      SCHI=SCHI+TI2*WGT
      AVGI=SI/SWGT
      SD=SWGT*IT/SI2
      CHI2A=SD*(SCHI/SWGT-AVGI*AVGI)/(DFLOAT(IT)-.999D0)
      SD=DSQRT(QNE/SD)
C
      IF(NPRN.EQ.0) GO TO 21
      TSI=DSQRT(TSI)
      WRITE(6,201) IT,TI,TSI,AVGI,SD,CHI2A
c      WRITE(6,*) IT,TI,TSI,AVGI,SD,CHI2A
      IF(NPRN.GE.0) GO TO 21
      DO 20 J=1,NDIM
 20   WRITE(6,202) J,(XI(I,J),DI(I,J),D(I,J),I=1,ND)
C     Print *, "vegas20"
C
C   REFINE GRID
C
 21   DO 23 J=1,NDIM
      XO=D(1,J)
      XN=D(2,J)
      D(1,J)=(XO+XN)/2.D0
      DT(J)=D(1,J)
      DO 22 I=2,NDM
      D(I,J)=XO+XN
      XO=XN
      XN=D(I+1,J)
      D(I,J)=(D(I,J)+XN)/3.D0
 22   DT(J)=DT(J)+D(I,J)
      D(ND,J)=(XN+XO)/2.D0
 23   DT(J)=DT(J)+D(ND,J)
C
      DO 28 J=1,NDIM
      RC=0.D0
      DO 24 I=1,ND
      R(I)=0.
      IF(D(I,J).LE.0.) GO TO 24
      XO=DT(J)/D(I,J)
      R(I)=((XO-QNE)/XO/DLOG(XO))**ALPH
 24   RC=RC+R(I)
      RC=RC/XND
      K=0
      XN=0.D0
      DR=XN
      I=K
 25   K=K+1
      DR=DR+R(K)
      XO=XN
      XN=XI(K,J)
 26   IF(RC.GT.DR) GO TO 25
      I=I+1
      DR=DR-RC
      XIN(I)=XN-(XN-XO)*DR/R(K)
      IF(I.LT.NDM) GO TO 26
      DO 27 I=1,NDM
 27   XI(I,J)=XIN(I)
 28   XI(ND,J)=QNE
C
      IF(IT.LT.ITMX.AND.ACC*DABS(AVGI).LT.SD) GO TO 9
 200  FORMAT('0INPUT PARAMETERS FOR VEGAS:  NDIM=',I3,'  NCALL=',F8.0
     1    /28X,'  IT=',I5,'  ITMX=',I5/28X,'  ACC=',G9.3
     2    /28X,'  MDS=',I3,'   ND=',I4/28X,'  (XL,XU)=',
     3    (T40,'( ',G12.6,' , ',G12.6,' )'))
 201  FORMAT(///' INTEGRATION BY VEGAS' / '0ITERATION NO.',I3,
     1    ':   INTEGRAL =',G14.8/21X,'STD DEV  =',G10.4 /
     2    ' ACCUMULATED RESULTS:   INTEGRAL =',G14.8 /
     3    24X,'STD DEV  =',G10.4 / 24X,'CHI**2 PER IT''N =',G10.4)
 202  FORMAT('0DATA FOR AXIS',I2 / ' ',6X,'X',7X,'  DELT I  ',
     1    2X,' CONV''CE  ',11X,'X',7X,'  DELT I  ',2X,' CONV''CE  '
     2   ,11X,'X',7X,'  DELT I  ',2X,' CONV''CE  ' /
     2    (' ',3G12.4,5X,3G12.4,5X,3G12.4))
      RETURN
      END

      SUBROUTINE SAVE(NDIM)
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      COMMON/BVEG2/XI(50,10),SI,SI2,SWGT,SCHI
      COMMON/BVE22/NDO,IT
      SAVE
C
C   STORES VEGAS DATA (UNIT 7) FOR LATER RE-INITIALIZATION
C
      WRITE(7,200) NDO,IT,SI,SI2,SWGT,SCHI,
     1             ((XI(I,J),I=1,NDO),J=1,NDIM)
      RETURN
      ENTRY RESTR(NDIM)
C
C   ENTERS INITIALIZATION DATA FOR VEGAS
C
      READ(7,200) NDO,IT,SI,SI2,SWGT,SCHI,
     1            ((XI(I,J),I=1,NDO),J=1,NDIM)
 200  FORMAT(2I8,4Z16/(5Z16))
      RETURN
      END

      SUBROUTINE ARAN9(QRAN,NDIM)
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      COMMON/SEED/NUM
      SAVE
      DIMENSION QRAN(10)
      DO 1 I=1,NDIM
      QRAN(I)=RANDOM(NUM)
    1 CONTINUE
      RETURN
      END
      
      FUNCTION RANDOM(SEED)
* ----------------------------------------------------------------
* Ref.: K. Park and K.W. Miller, Comm. of the ACM 31 (1988) p.1192
* Use seed = 1 as first value.
* ----------------------------------------------------------------
      IMPLICIT INTEGER(A-Z)
      DOUBLE PRECISION MINV,RANDOM
      SAVE
      PARAMETER(M=2147483647,A=16807,Q=127773,R=2836)
      PARAMETER(MINV=0.46566128752458d-09)
      HI = SEED/Q
      LO = MOD(SEED,Q)
      SEED = A*LO - R*HI
      IF(SEED.LE.0) SEED = SEED + M
      RANDOM = SEED*MINV
      RETURN
      END
