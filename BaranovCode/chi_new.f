C###################################################################C
C                                                                   C
C     Routines to simulate the production of Quarkonia mesons       C
C                                                                   C
C    in the framework of Semihard Approach and/or Parton Model      C
C   via gluon-gluon, photon-gluon or photon-photon 2-->1 fusion     C
C                                                                   C
C###################################################################C
C
C     Complex arithmetics for chi_c helicities
C     
C     Clebsch-Gordan coefficients in the meson rest frame
C     all apply to amplitudes
C
C     Reduced version: no polarised beams (SYMSYM only)
C                      simple gluon/photon polarization: kt*kt
C
      DOUBLE PRECISION FUNCTION FXN(X,WGT)
C*******************************************************************C
C      The integrand expression for VEGAS  +  kinematics            C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      DIMENSION X(10), SYMSYM(-1:1,-1:1)
      COMMON/BVEGG/NDIM,NCALL,ITMX,NPRN
      COMMON/CONST/CME,PI,ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      COMMON/OUTPX/XChi0,XChi1(-1:1),XChi2(-2:2)
      COMMON/OUTPJ/SYMJ0,SYMJ1(-1:1),SYMJ2(-2:2)
      COMMON/TOTAL/CXTOT0,CXTOT1,CXTOT2
      COMMON/MOMEN/p1(4),p2(4),g1(4),g2(4),pj(4),a1(4),a2(4)
      COMMON/KINEM/x1,x2,g1T,g2T, pjT,yj,zj
      COMMON/XFACT/CXFACT,COLOR
      COMMON/TYPE/ITYPE,MODEL,IB1,IB2,NJ,NL
      COMMON/SEED/NUM
      SAVE
      JFXN=0
      FXN=0.D0
C... Random numbers within VEGAS
      g1T =dexp(X(1))!Transverse momentum of 1st parton
      g2T =dexp(X(2))!Transverse momentum of 2nd parton
       yj = X(3)     !Rapidity of the  final Quarkonium 
C... Random numbers other than VEGAS
      ph1 = RANDOM(NUM)*2.*PI !Azmt angle of 1st parton
      ph2 = RANDOM(NUM)*2.*PI !Azmt angle of 2nd parton
      cos1=dcos(ph1)
      sin1=dsin(ph1)
      cos2=dcos(ph2)
      sin2=dsin(ph2)
C
C... Particle and Parton momenta,  P(i) = (Beam, P_t, P_t, Energy)
      p1(1) = CME/2.
      p1(2) = 0.D0
      p1(3) = 0.D0
      p1(4) = CME/2.
      p2(1) =-CME/2.
      p2(2) = 0.D0
      p2(3) = 0.D0
      p2(4) = CME/2.
      g1(2) = g1T*cos1
      g1(3) = g1T*sin1
      g2(2) = g2T*cos2
      g2(3) = g2T*sin2
      pj(2) = g1(2) +g2(2) 
      pj(3) = g1(3) +g2(3)
      pjT= dsqrt(pj(2)**2+pj(3)**2)
      tmj= dsqrt(XJ2 +pjT*pjT)
      pj(1) = tmj*dsinh(yj)
      pj(4) = tmj*dcosh(yj)
      x1 = tmj*dexp(yj)/CME
      x2 = tmj*dexp(-yj)/CME
      IF(x1.GE..999D0 .OR. x2.GE..999D0) GOTO 99
      IF(x1.LE. 1.D-6 .OR. x2.LE. 1.D-6) GOTO 99
      y1 = 15.D0
      y2 =-15.D0
      DO 3 k=1,7  !Iterative solution
      IF(dabs(y1).GE.20.) GOTO 99
      IF(dabs(y2).GE.20.) GOTO 99
      arg1= CME*(1.-x1)-g2T*dexp(y2)
      arg2= CME*(1.-x2)-g1T*dexp(-y1)
      y1 = dlog(arg1/g1T)
      y2 =-dlog(arg2/g2T)
    3 CONTINUE
      IF(y1.LE. 0.D0) GOTO 99
      IF(y2.GE. 0.D0) GOTO 99
      g1(1) = CME/2. -g1T*dsinh(y1)
      g1(4) = CME/2. -g1T*dcosh(y1)
      g2(1) =-CME/2. -g2T*dsinh(y2)
      g2(4) = CME/2. -g2T*dcosh(y2)
C
C... Gluon/photon polarization vectors
      CALL GAUGEG(cos1,sin1,cos2,sin2)
C
C... Tests of momentum conservation
      ptest1=pj(1)-g1(1)-g2(1)
      ptest4=pj(4)-g1(4)-g2(4)
      IF(dabs(ptest1).GE. 1.D-5) GOTO 99
      IF(dabs(ptest4).GE. 1.D-5) GOTO 99
      zj = DOT(pj,p2)/DOT(g1,p2) -0.001
      IF(zj .LE. 0.D0) GOTO 99
      IF(zj .GT. 1.D0) GOTO 99
C
C... Phase Space boundary and integration Jacobian
      G1W = DOT(g1,g1)
      G2W = DOT(g2,g2)
      STOT= CME*CME
      SXX = XJ2
      IF(AF2(SXX,G1W,G2W).LE.0.D0) GOTO 99
      S12 = G1W +2.*DOT(g1,p2)
      S21 = G2W +2.*DOT(g2,p1)
      IF(GF(S12,G2W,XJ2,0.D0,G1W,0.D0).GE.0.D0) GOTO 99
      IF(GF(S21,G1W,XJ2,0.D0,G2W,0.D0).GE.0.D0) GOTO 99
C
      XJacob = PI*(4.*(g1T*g2T)**2)
C
C... Flux factor, Couplings, Parton distributions
      CALL FACTOR(CXFACT)
      CXFACT = CXFACT*XJacob
C
C... Partonic cross section
      CALL XSEC(SYMSYM)
C
      CXTOT = 0.D0
      DO 13 j=-NJ,NJ
      DO 13 l=-NL,NL
      CXTOT =CXTOT + CXFACT*SYMSYM(j,l)
   13 CONTINUE
C
      IF(NL.EQ.0 .AND. NJ.EQ.0) CXTOT0=CXTOT
      IF(NL.EQ.0 .AND. NJ.EQ.1) CXTOT1=CXTOT

      IF(NL.NE.1 .OR. NJ.NE.1) GOTO 19       
      XChi0     = CXFACT*SYMJ0
      DO 14 mj=-1,1
   14 XChi1(mj) = CXFACT*SYMJ1(mj)
      DO 15 mj=-2,2
   15 XChi2(mj) = CXFACT*SYMJ2(mj)
      CXTOT0= XChi0
      CXTOT1= Xchi1(-1) +Xchi1(0) +Xchi1(1)
      CXTOT2= Xchi2(-2) +Xchi2(-1)+Xchi2(0) +Xchi2(1) +Xchi2(2)
      CXTOT = CXTOT0 + CXTOT1 + CXTOT2
   19 CONTINUE
C
      IF(CXTOT.LT.0.D0) THEN
cw      write(6,*) 'Negative total cross section !!!!!!!!!'
        CXTOT=0.D0
      ENDIF
      IF(CXTOT.GE.1.D-20 .AND. CXTOT.LT.1.D20) JFXN=1
      IF(JFXN.EQ.0) RETURN
      FXN=CXTOT !*(1.+pjT**2) !to enlarge statistics at high pjT
C
C... Filling histograms
      IF(NPRN.NE.1) RETURN
      hwgt = sngl(WGT)/FLOAT(ITMX)
      CALL WRIOUT(hwgt)
   99 RETURN
      END

      SUBROUTINE XSEC(SYMSYM)
C*******************************************************************C
C      Partonic differential cross section                          C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SYMSYM(-1:1,-1:1)
      COMMON/TYPE/ITYPE,MODEL,IB1,IB2,NJ,NL
      COMMON/CONST/CME,PI,ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      COMMON/LOOPJ/AMPJ(2,-1:1,-1:1)
      COMMON/LOOPX/AMPX0(2),AMPX1(2,-1:1),AMPX2(2,-2:2)
      COMMON/OUTPJ/SYMJ0,SYMJ1(-1:1),SYMJ2(-2:2)
      COMMON/GMUNU/DF(4,4),DC(4),EP(4,4,4,4)
      COMMON/XFACT/CXFACT,COLOR
      SAVE
C
C...Preparing the full kinematics, namely:
C    Construct the reference ort system for chi_c decays,
C    Generate two-body radiative decay chi_c --> Jpsi + gamma,
C    Construct the reference ort system for J/psi decays,
C    Generate two-body leptonic decay  J/psi --> mu+ mu-
C
      IF((NL*NJ).EQ.1) CALL DECAYS
      IF((NL+NJ).EQ.1) CALL JDECAY
      IF((NL+NJ).EQ.0) CALL EDECAY
C
C...Chi-meson polarization vectors
      IF(NJ.EQ.1) CALL GAUGEJ
C
C...Evaluating the Matrix Elements
      CALL FEYNJ
      IF((NL*NJ).EQ.1) CALL CHIJ 
C
C...Color coefficients
      DO 1 j=-1,1
      DO 1 l=-1,1
      SYMSYM(j,l)=0.D0
   1  CONTINUE
      SYMJ0=0.D0
      DO 3 k=-1,1
   3  SYMJ1(k)=0.D0
      DO 5 k=-2,2
   5  SYMJ2(k)=0.D0
      COLOR=0.D0
      IF(ITYPE.EQ.10) COLOR=  8./3.      ! glu+glu->Eta(1)
      IF(ITYPE.EQ.12) COLOR=  8./3.      ! glu+glu->Chi(1)
      IF(ITYPE.EQ.30) COLOR= 12.*QC2*QC2 ! gam+gam->Eta(1)
      IF(ITYPE.EQ.32) COLOR= 12.*QC2*QC2 ! gam+gam->Chi(1)
      IF(ITYPE.EQ.40) COLOR=  5./6.      ! glu+glu->Eta(8)
      IF(ITYPE.EQ.41) COLOR=  3./2.      ! glu+glu->Psi(8)   
      IF(ITYPE.EQ.42) COLOR=  5./6.      ! glu+glu->Chi(8)
      IF(ITYPE.EQ.50) COLOR=  2.*QC2     ! gam+glu->Eta(8)
      IF(ITYPE.EQ.52) COLOR=  2.*QC2     ! gam+glu->Chi(8)
C
C...Chi_J production amplitudes squared
      DO 11 j=-NJ,NJ
      DO 11 l=-NL,NL
      SYMSYM(j,l)=SYMSYM(j,l) + AMPJ(1,j,l)*AMPJ(1,j,l)*COLOR
     .                        + AMPJ(2,j,l)*AMPJ(2,j,l)*COLOR
  11  CONTINUE
      IF(NL.NE.1 .OR. NJ.NE.1) GOTO 19 
      SYMJ0   = SYMJ0   +AMPX0(1)*AMPX0(1)     +AMPX0(2)*AMPX0(2)
      DO 13 k=-1,1
  13  SYMJ1(k)=SYMJ1(k) +AMPX1(1,k)*AMPX1(1,k) +AMPX1(2,k)*AMPX1(2,k)
      DO 15 k=-2,2
  15  SYMJ2(k)=SYMJ2(k) +AMPX2(1,k)*AMPX2(1,k) +AMPX2(2,k)*AMPX2(2,k)
  19  CONTINUE
C
      IF(NL.NE.1 .OR. NJ.NE.1) GOTO 29
      SYMJ0   = SYMJ0   *COLOR
      DO 23 k=-1,1
  23  SYMJ1(k)= SYMJ1(k)*COLOR
      DO 25 k=-2,2
  25  SYMJ2(k)= SYMJ2(k)*COLOR
  29  RETURN
      END

      SUBROUTINE GAUGEJ
C********************************************************************C
C    Polarization vectors of outgoing meson in a chosen ref.frame    C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/ORTS/OX1(4),OY1(4),OZ1(4),OX2(4),OY2(4),OZ2(4)
      COMMON/SPINJ/EJ(2,4,-1:1)
      SAVE
      sq2=sqrt(2.)
C... Longitudinal polarization,   zero helicity
      DO 1 j=1,3
      EJ(1,j,0)= OZ1(j)
   1  EJ(2,j,0)= 0.D0
      EJ(1,4,0)= 0.D0
      EJ(2,4,0)= 0.D0
C... Transverse polarization, negative helicity
      DO 2 j=1,3
      EJ(1,j,-1)=-OX1(j)/sq2
   2  EJ(2,j,-1)= OY1(j)/sq2
      EJ(1,4,-1)= 0.D0
      EJ(2,4,-1)= 0.D0
C... Transverse polarization, positive helicity
      DO 3 j=1,3
      EJ(1,j,1)= OX1(j)/sq2
   3  EJ(2,j,1)= OY1(j)/sq2
      EJ(1,4,1)= 0.D0
      EJ(2,4,1)= 0.D0
      RETURN
      END

      SUBROUTINE GAUGEG(cos1,sin1,cos2,sin2)
C*******************************************************************C
C      Photon and Gluon polarization vectors                        C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/MOMEN/p1(4),p2(4),g1(4),g2(4),pj(4),a1(4),a2(4)
      SAVE
      a1(1)=0.D0
      a1(2)=cos1
      a1(3)=sin1
      a1(4)=0.D0
      a2(1)=0.D0
      a2(2)=cos2
      a2(3)=sin2
      a2(4)=0.D0
C
      RETURN !Comment out this line for gauge invariance tests
      DO 9 i=1,4
      a1(i)=g1(i) 
      a2(i)=g2(i)
    9 Continue
      RETURN
      END

      SUBROUTINE FEYNJ
C*******************************************************************C
C      Matrix Elements for all Feynman diagrams                     C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION g3(4)
      COMMON/TYPE/ITYPE,MODEL,IB1,IB2,NJ,NL
      COMMON/CONST/CME,PI,ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      COMMON/FREST/p1(4),p2(4),g1(4),g2(4),pj(4),a1(4),a2(4)
      COMMON/LOOPJ/AMPJ(2,-1:1,-1:1)
      COMMON/LOOP1/AMP1(4)
      COMMON/GMUNU/DF(4,4),DC(4),EP(4,4,4,4)
      COMMON/SPINJ/EJ(2,4,-1:1)
      SAVE
      W1=DOT(g1,g1)
      W2=DOT(g2,g2)
      Z0=(W1+W2-XJ2)/2. ! = -DOT(g1,g2)
C...Clean main array
      DO 1 j=-1,1
      DO 1 l=-1,1
      AMPJ(1,j,l)= 0.D0
      AMPJ(2,j,l)= 0.D0
   1  CONTINUE
      DO 2 i=1,4
      AMP1(i)=0.D0
   2  g3(i) =0.5*(g1(i)-g2(i))
C
C... Eta = 1S0 state
      IF(NJ.NE.0 .OR. NL.NE.0) GOTO 110
      DO 105 m1=1,4
      DO 105 m2=1,4
      DO 105 k1=1,4
      DO 105 k2=1,4
      !Diagram << 1 >> = << 2 >>
      TERM0 = EP(k1,k2,m1,m2)*g1(k1)*g2(k2)*DC(k1)*DC(k2)/2./Z0
      AMPJ(1,0,0) = AMPJ(1,0,0) + TERM0*a1(m1)*a2(m2)*DC(m1)*DC(m2)
 105  CONTINUE
C
C... Psi = 3S1 state
 110  IF(NJ.NE.1 .OR. NL.NE.0) GOTO 120
      DO 115 m1=1,4
      DO 115 m2=1,4
      DO 115 mj=1,4      
      !Diagram << 1 >> = - << 2 >> = Diagram << 3 >>/2 
      TERM0 = DF(m1,m2)*g3(mj) - DF(m1,mj)*pj(m2) + DF(m2,mj)*pj(m1)
      TERM0 = TERM0*(1./Z0 +2./XJ2)*a1(m1)*a2(m2)*DC(m1)*DC(m2)
      AMP1(mj) = AMP1(mj) + TERM0*DC(mj)
      DO 115 j=-NJ,NJ
      AMPJ(1,j,0) = AMPJ(1,j,0) + TERM0*EJ(1,mj,j)*DC(mj)
      AMPJ(2,j,0) = AMPJ(2,j,0) + TERM0*EJ(2,mj,j)*DC(mj)
 115  CONTINUE
C
C... Chi = 3PJ state
 120  IF(NJ.NE.1 .OR. NL.NE.1) GOTO 199
      DO 125 m1=1,4
      DO 125 m2=1,4
      DO 125 mj=1,4
      DO 125 ml=1,4
      !Zeroth order in relative momentum, Diagram << 1 >> = - << 2 >>
      TERM0 = - DF(m1,m2)*g3(mj) + DF(m1,mj)*g1(m2) - DF(m2,mj)*g2(m1)
      !Linear terms in relative momentum, Diagram << 1 >> = << 2 >>
      TERM1 = - DF(m1,mj)*DF(m2,ml)*W1  + DF(m1,mj)*DF(m2,ml)*Z0
     .        + DF(m1,mj)*pj(m2)*g3(ml) - DF(m1,ml)*DF(m2,mj)*W2
     .        + DF(m1,ml)*DF(m2,mj)*Z0  - DF(m1,ml)*pj(m2)*g3(mj)
     .        - DF(m2,mj)*pj(m1)*g3(ml) + DF(m2,ml)*pj(m1)*g3(mj)
      Z1  = -2.*g3(ml) 
      SUM = TERM1/Z0 - 2*XC2*TERM0*Z1/(Z0)**2
      SUM = SUM*a1(m1)*a2(m2)*DC(m1)*DC(m2)
      DO 125 j=-NJ,NJ
      DO 125 l=-NL,NL
      AMPJ(1,j,l) = AMPJ(1,j,l)                           !real part
     .            + SUM*EJ(1,mj,j)*DC(mj)*EJ(1,ml,l)*DC(ml)
     .            - SUM*EJ(2,mj,j)*DC(mj)*EJ(2,ml,l)*DC(ml)
      AMPJ(2,j,l) = AMPJ(2,j,l)                           !imag part
     .            + SUM*EJ(1,mj,j)*DC(mj)*EJ(2,ml,l)*DC(ml)
     .            + SUM*EJ(2,mj,j)*DC(mj)*EJ(1,ml,l)*DC(ml)
 125  CONTINUE
 199  RETURN
      END

      SUBROUTINE BEAMS
      COMMON/TYPE/ITYPE,MODEL,IB1,IB2,NJ,NL
      SAVE
      IQ= MOD(ITYPE,10)
      IB= (ITYPE-IQ)/10
                               IB1=1
      IF(IB.EQ.1 .OR. IB.EQ.4) IB1=8
                               IB2=8
      IF(IB.EQ.3 .OR. IB.EQ.6) IB2=1
                  NJ=1
      IF(IQ.EQ.0) NJ=0
                  NL=0
      IF(IQ.EQ.2) NL=1
      RETURN
      END

      SUBROUTINE CHIJ
C*******************************************************************C
C     Converts (S,L) amplitudes into (J,mj) amplitudes              C
C       using the Clebsch-Gordan coefficients                       C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/LOOPJ/AMPJ(2,-1:1,-1:1)
      COMMON/LOOPX/AMPX0(2),AMPX1(2,-1:1),AMPX2(2,-2:2)
      SAVE
      sq2=sqrt(2.)
      sq3=sqrt(3.)
      sq6=sqrt(6.)
      DO 1  i=1,2 !real and imag parts
      AMPX2(i, 2)=AMPJ(i, 1, 1)
      AMPX2(i, 1)=AMPJ(i, 1, 0)/sq2 + AMPJ(i, 0, 1)/sq2
      AMPX2(i, 0)=AMPJ(i, 1,-1)/sq6 + AMPJ(i,-1, 1)/sq6
     .           +AMPJ(i, 0, 0)*sq2/sq3       !
      AMPX2(i,-1)=AMPJ(i,-1, 0)/sq2 + AMPJ(i, 0,-1)/sq2
      AMPX2(i,-2)=AMPJ(i,-1,-1)
      AMPX1(i, 1)=AMPJ(i, 1, 0)/sq2 - AMPJ(i, 0, 1)/sq2
      AMPX1(i, 0)=AMPJ(i, 1,-1)/sq2 - AMPJ(i,-1, 1)/sq2
      AMPX1(i,-1)=AMPJ(i,-1, 0)/sq2 - AMPJ(i, 0,-1)/sq2
      AMPX0(i)   =AMPJ(i, 1,-1)/sq3 + AMPJ(i,-1, 1)/sq3
     .           -AMPJ(i, 0, 0)/sq3           !
    1 CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION ALPHAS(Q2,XNF)
C*******************************************************************C
C     Strong coupling constant                                      C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      ALPHAS=0.D0
      ALAM2 =(0.2)**2
      Q2RUN =dabs(Q2)
      IF(Q2RUN.LE.ALAM2) RETURN
      ALPHA =12.*3.1415926/(33.-2.*XNF)/dlog(Q2RUN/ALAM2)
      IF(ALPHA.GE.1.D0) ALPHA=1.D0
    9 ALPHAS=ALPHA
      RETURN
      END

      DOUBLE PRECISION FUNCTION AF2(X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      AF2 = X*X + Y*Y + Z*Z - 2.*X*Y - 2.*X*Z - 2.*Y*Z
      RETURN
      END

      DOUBLE PRECISION FUNCTION AF(X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      AF = DSQRT(X*X + Y*Y + Z*Z - 2.*X*Y - 2.*X*Z - 2.*Y*Z)
      RETURN
      END

      DOUBLE PRECISION FUNCTION GF(X,Y,Z,U,V,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      GF =          X*Z*W + X*U*V + Y*Z*V + Y*U*W
     . -X*Y*(Z+U+V+W-X-Y) -Z*U*(X+Y+V+W-Z-U) -V*W*(X+Y+Z+U-V-W)
      RETURN
      END

      SUBROUTINE METRIC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/GMUNU/DF(4,4),DC(4),EP(4,4,4,4)
      DO 1 I=1,4
      DO 1 J=1,4
      DF(I,J)= 0.D0
      DO 1 M=1,4
      DO 1 N=1,4
      EP(I,J,M,N)= 0.D0
    1 CONTINUE !P(i) =(Beam, P_t, P_t, Energy)
      DO 2 I=1,3
      DF(I,I)=-1.D0
    2 DC(I)  =-1.D0
      DF(4,4)= 1.D0
      DC(4)  = 1.D0
      EP(1,2,3,4)= 1.D0
      EP(1,2,4,3)=-1.D0
      EP(1,3,2,4)=-1.D0
      EP(1,3,4,2)= 1.D0
      EP(1,4,2,3)= 1.D0
      EP(1,4,3,2)=-1.D0
      EP(2,1,3,4)=-1.D0
      EP(2,1,4,3)= 1.D0
      EP(2,3,1,4)= 1.D0
      EP(2,3,4,1)=-1.D0
      EP(2,4,1,3)=-1.D0
      EP(2,4,3,1)= 1.D0
      EP(3,1,2,4)= 1.D0
      EP(3,1,4,2)=-1.D0
      EP(3,2,1,4)=-1.D0
      EP(3,2,4,1)= 1.D0
      EP(3,4,1,2)= 1.D0
      EP(3,4,2,1)=-1.D0
      EP(4,1,2,3)=-1.D0
      EP(4,1,3,2)= 1.D0
      EP(4,2,1,3)= 1.D0
      EP(4,2,3,1)=-1.D0
      EP(4,3,1,2)=-1.D0
      EP(4,3,2,1)= 1.D0
      RETURN
      END

      DOUBLE PRECISION FUNCTION DOT(xx,yy)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION xx(4),yy(4)   !P(i) = (Beam, P_t, P_t, Energy)
      DOT= -xx(1)*yy(1) -xx(2)*yy(2) -xx(3)*yy(3) +xx(4)*yy(4)
      RETURN
      END

      DOUBLE PRECISION FUNCTION EPS(q1,q2,q3,q4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION q1(4),q2(4),q3(4),q4(4)
      G1 = q1(2)*q2(3)*q3(4) +q1(3)*q2(4)*q3(2) +q1(4)*q2(2)*q3(3)
     .    -q1(4)*q2(3)*q3(2) -q1(2)*q2(4)*q3(3) -q1(3)*q2(2)*q3(4)
      G2 = q1(1)*q2(3)*q3(4) +q1(3)*q2(4)*q3(1) +q1(4)*q2(1)*q3(3)
     .    -q1(4)*q2(3)*q3(1) -q1(1)*q2(4)*q3(3) -q1(3)*q2(1)*q3(4)
      G3 = q1(1)*q2(2)*q3(4) +q1(2)*q2(4)*q3(1) +q1(4)*q2(1)*q3(2)
     .    -q1(4)*q2(2)*q3(1) -q1(1)*q2(4)*q3(2) -q1(2)*q2(1)*q3(4)
      G4 = q1(1)*q2(2)*q3(3) +q1(2)*q2(3)*q3(1) +q1(3)*q2(1)*q3(2)
     .    -q1(3)*q2(2)*q3(1) -q1(1)*q2(3)*q3(2) -q1(2)*q2(1)*q3(3)
      EPS= -G1*q4(1) +G2*q4(2) -G3*q4(3) +G4*q4(4)
      RETURN
      END

      SUBROUTINE GAUGEX
C********************************************************************C
C       Tensor Polarization of outgoing meson, general case          C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION EX(2,4,4,-2:2)
      COMMON/LOOPX/AMPX0(2),AMPX1(2,-1:1),AMPX2(2,-2:2)
      COMMON/SPINX/EX1(2,4),EX2(2,4,4)
      COMMON/SPINJ/EJ(2,4,-1:1)
      COMMON/GMUNU/DF(4,4),DC(4),EP(4,4,4,4)
      SAVE
      sq2=sqrt(2.)
      sq3=sqrt(3.)
      sq6=sqrt(6.)
      DO 1 i=1,2
      DO 1 m=1,4
      EX1(i,m)=0.D0
      DO 1 n=1,4
      EX2(i,m,n)=0.D0
    1 CONTINUE
C    Polarization tensor
      DO 10 m=1,4
      DO 10 n=1,4
      EX(1,m,n, 2)= EJ(1,m, 1)*EJ(1,n, 1) -EJ(2,m, 1)*EJ(2,n, 1)
      EX(2,m,n, 2)= EJ(1,m, 1)*EJ(2,n, 1) +EJ(2,m, 1)*EJ(1,n, 1)
      EX(1,m,n, 1)=(EJ(1,m, 1)*EJ(1,n, 0) +EJ(1,m, 0)*EJ(1,n, 1))/sq2
      EX(2,m,n, 1)=(EJ(2,m, 1)*EJ(1,n, 0) +EJ(1,m, 0)*EJ(2,n, 1))/sq2
      EX(1,m,n, 0)=(EJ(1,m,-1)*EJ(1,n, 1) +EJ(1,m, 1)*EJ(1,n,-1))/sq6
     .            -(EJ(2,m,-1)*EJ(2,n, 1) +EJ(2,m, 1)*EJ(2,n,-1))/sq6
     .            + EJ(1,m, 0)*EJ(1,n, 0)*sq2/sq3  !
      EX(2,m,n, 0)=(EJ(2,m,-1)*EJ(1,n, 1) +EJ(2,m, 1)*EJ(1,n,-1))/sq6
     .            +(EJ(1,m,-1)*EJ(2,n, 1) +EJ(1,m, 1)*EJ(2,n,-1))/sq6
      EX(1,m,n,-1)=(EJ(1,m,-1)*EJ(1,n, 0) +EJ(1,m, 0)*EJ(1,n,-1))/sq2
      EX(2,m,n,-1)=(EJ(2,m,-1)*EJ(1,n, 0) +EJ(1,m, 0)*EJ(2,n,-1))/sq2
      EX(1,m,n,-2)= EJ(1,m,-1)*EJ(1,n,-1) -EJ(2,m,-1)*EJ(2,n,-1)
      EX(2,m,n,-2)= EJ(1,m,-1)*EJ(2,n,-1) +EJ(2,m,-1)*EJ(1,n,-1)
   10 CONTINUE
C    Vector meson production amplitude
      DO 11 j=-1,1
      DO 11 m=1,4
      EX1(1,m)= EX1(1,m) +EJ(1,m,j)*AMPX1(1,j) -EJ(2,m,j)*AMPX1(2,j)
      EX1(2,m)= EX1(2,m) +EJ(2,m,j)*AMPX1(1,j) +EJ(1,m,j)*AMPX1(2,j)
   11 CONTINUE
C    Tensor meson production amplitude
      DO 12 j=-2,2
      DO 12 m=1,4
      DO 12 n=1,4
      EX2(1,m,n)= EX2(1,m,n)
     .          + EX(1,m,n,j)*AMPX2(1,j) - EX(2,m,n,j)*AMPX2(2,j)
      EX2(2,m,n)= EX2(2,m,n)
     .          + EX(2,m,n,j)*AMPX2(1,j) + EX(1,m,n,j)*AMPX2(2,j)
   12 CONTINUE
      RETURN
      END

      SUBROUTINE VGAMMA(WGT1,WGT2)
C********************************************************************C
C    Radiative decay vertices Chi_J --> Jpsi + gamma                 C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION eg1(4),eg2(4),SUM5(4,4,4,4)
      COMMON/DECAY/xpsi,xmu2,vj(4),vg(4),v1(4),v2(4),csg,phg,csm,phm
      COMMON/FREST/p1(4),p2(4),g1(4),g2(4),pj(4),a1(4),a2(4)
      COMMON/CONST/CME,PI,ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      COMMON/GMUNU/DF(4,4),DC(4),EP(4,4,4,4)
      COMMON/TOTAL/CXTOT0,CXTOT1,CXTOT2
      COMMON/VERTG/SUM1(4,4,4,4),SUM2(4,4),SUMJ(4,4),SUME(4,4)
      COMMON/SPINX/EX1(2,4),EX2(2,4,4)
      COMMON/XFACT/CXFACT,COLOR
      SAVE
C
      xpsi2=xpsi*xpsi
      vg4=(XJ2-xpsi2)/XJ/2.
C
C... Photon polarization vectors
c     CALL GAUGAM(eg1,eg2)
C
C... J/psi and photon spin density matrices
      DO 1 m1=1,4
      DO 1 m2=1,4
      SUMJ(m1,m2)= v1(m1)*v2(m2) +v2(m1)*v1(m2) -xpsi2*DF(m1,m2)/2.
c     SUME(m1,m2)= eg1(m1)*eg1(m2) + eg2(m1)*eg2(m2)
      SUM2(m1,m2)=-DF(m1,m2) - vg(m1)*vg(m2)/vg4/vg4
    1 CONTINUE
C
C... Chi_1 and chi_2 decay spin density matrices
      DO 2 m1=1,4
      DO 2 m2=1,4
      DO 2 n1=1,4
      DO 2 n2=1,4
      SUM1(m1,m2,n1,n2) =
     . - DF(m1,m2)*vg(n1)*vg(n2) + DF(m1,n2)*vg(m2)*vg(n1)
     . - DF(n1,n2)*vg(m1)*vg(m2) + DF(m2,n1)*vg(m1)*vg(n2)
c     SUM5(m1,m2,n1,n2) =
c    .   vg(m1)*vg(m2)*(eg1(n1)*eg1(n2) +eg2(n1)*eg2(n2) +2*DF(n1,n2))
c    . - vg(m1)*vg(n2)*(eg1(m2)*eg1(n1) +eg2(m2)*eg2(n1) +2*DF(m2,n1))
c    . - vg(m2)*vg(n1)*(eg1(m1)*eg1(n2) +eg2(m1)*eg2(n2) +2*DF(m1,n2))
c    . + vg(n1)*vg(n2)*(eg1(m1)*eg1(m2) +eg2(m1)*eg2(m2) +2*DF(m1,m2))
    2 CONTINUE
C
C... Overall decay weights
      WGT1=0.D0
      WGT2=0.D0
      WGT4=0.D0
      WGT5=0.D0
      DO 3 i=1,2
      DO 3 m1=1,4
      DO 3 m2=1,4
      DO 3 n1=1,4
      DO 3 n2=1,4
c     WGT5= WGT5 -EX1(i,m1)*EX1(i,m2)*SUM5(m1,m2,n1,n2)*SUMJ(n1,n2)
      WGT1= WGT1 +EX1(i,m1)*EX1(i,m2)*SUM1(m1,m2,n1,n2)*SUMJ(n1,n2)
      WGT2= WGT2 +EX2(i,m1,n1)*EX2(i,m2,n2)*SUMJ(m1,m2)*SUM2(n1,n2)
c     WGT4= WGT4 +EX2(i,m1,n1)*EX2(i,m2,n2)*SUMJ(m1,m2)*SUME(n1,n2)
    3 CONTINUE
C
c     WGT5=WGT5*CXFACT*COLOR*(3./xpsi2)*(3./2.)/vg4/vg4*XJ2/(xpsi2+XJ2)
      WGT1=WGT1*CXFACT*COLOR*(3./xpsi2)*(3./2.)/vg4/vg4*XJ2/(xpsi2+XJ2)
      WGT2=WGT2*CXFACT*COLOR*(3./xpsi2)*(3./2.)
c     WGT4=WGT4*CXFACT*COLOR*(3./xpsi2)*(3./2.)
      RETURN
      END

      SUBROUTINE GAUGAM(eg1,eg2)
C********************************************************************C
C   Polarization vectors of the final-state photon, chi_J rest frame C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION eg1(4),eg2(4)
      COMMON/DECAY/xpsi,xmu2,vj(4),vg(4),v1(4),v2(4),csg,phg,csm,phm
      SAVE
      QJ=dsqrt(vg(3)**2 +vg(2)**2 +vg(1)**2)
      Qt=dsqrt(vg(3)**2 +vg(2)**2)
      CT=vg(1)/QJ
      ST=Qt/QJ
      CP=vg(2)/Qt
      SP=vg(3)/Qt
C... Transverse polarization 1
      eg1(1)=-ST
      eg1(2)= CT*CP
      eg1(3)= CT*SP
      eg1(4)= 0.D0
C... Transverse polarization 2
      eg2(1)= 0.D0
      eg2(2)=-SP
      eg2(3)= CP
      eg2(4)= 0.D0
  99  RETURN
      END

      SUBROUTINE DECAYS
C********************************************************************C
C     Construct the reference ort system for chi_c decays,           C
C     Generate two-body radiative decay chi_c --> Jpsi + gamma,      C
C     Store cos(theta) and phi v.r.t. orts in the chi_c rest frame,  C
C     Store J/psi and photon momenta in the Beam*Beam c.m.system     C
C      then                                                          C
C     Construct the reference ort system for J/psi decays,           C
C     Generate two-body leptonic decay  J/psi --> mu+ mu-,           C
C     Store cos(theta) and phi v.r.t. orts in the pj rest frame      C
C     Store lepton momenta v1 and v2 in the Beam*Beam c.m.system     C
C                                                                    C
C     /MOMEN/ all variables are in the C.M.S. frame                  C
C     /FREST/ all variables are in the Chi_J rest frame              C
C     /DECAY/ vj(4),vg(4),v1(4),v2(4),csg,phg  are in Chi_J restsyst C
C                                     csm,phm  are in J/psi restsyst C
C     /ORTS/  OX1(4),OY1(4),OZ1(4)             are in Chi_J restsyst C
C                OX2(4),OY2(4),OZ2(4)          are in J/psi restsyst C
C     /FJPSI/ f1(4),f2(4),vmu1(4),vmu2(4)      are in J/psi restsyst C
C     /FLABS/ qpsi(4),qgam(4),qmu1(4),qmu2(4)  are in the C.M.System C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/FJPSI/f1(4),f2(4),vmu1(4),vmu2(4)
      COMMON/FLABS/qpsi(4),qgam(4),qmu1(4),qmu2(4)
      COMMON/MOMEN/b1(4),b2(4),g1(4),g2(4),pj(4),a1(4),a2(4)
      COMMON/CONST/CME,PI,ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      COMMON/ORTS/OX1(4),OY1(4),OZ1(4),OX2(4),OY2(4),OZ2(4)
      COMMON/FREST/d1(4),d2(4),q1(4),q2(4),qj(4),e1(4),e2(4)
      COMMON/DECAY/xpsi,xmu2,vj(4),vg(4),v1(4),v2(4),csg,phg,csm,phm
      COMMON/FRAME/JF1,JF2
      COMMON/SEED/NUM
      SAVE
      xpsi2=xpsi**2
      xgam2=0.D0
C... Particle momenta in the chi_J rest frame
      CALL JREST(XJ,pj,b1,b2,d1,d2)
      CALL JREST(XJ,pj,g1,g2,q1,q2)
      CALL JREST(XJ,pj,a1,a2,e1,e2)
      DO 1 i=1,3
    1 qj(i)=0.D0
      qj(4)=XJ
C... Construct the reference ort system for chi_c decays
      CALL JFRAME(JF1,d1,d2,OX1,OY1,OZ1)
C
C... Generate two-body radiative decay chi_c --> Jpsi + gamma
      vj(4)= (XJ2+xpsi2-xgam2)/2./XJ
      vg(4)= (XJ2-xpsi2+xgam2)/2./XJ
      vtot = AF(XJ2,xpsi2,xgam2)/2./XJ
C    Two Euler's angles at random
      cosb=-1. + 2.*RANDOM(NUM)
      alph= 3.14159*RANDOM(NUM)*2.
      sinb= dabs(1.-cosb*cosb)
      sinb= dsqrt(sinb)
      cosa= dcos(alph)
      sina= dsin(alph)
      vj(1)= vtot*cosb
      vj(2)= vtot*sinb*cosa
      vj(3)= vtot*sinb*sina
      do 2 i=1,3
    2 vg(i)=-vj(i)
C     Angles with respect to the predefined orts
      csg = DOT3(vg,OZ1)/vtot
      cph = DOT3(vg,OX1)
      sph = DOT3(vg,OY1)
      phg = DATAN2(sph,cph)
C    Decay momenta in the beam-beam c.m.s.
      CALL CMSYS(XJ,pj,vj,vg,qpsi,qgam)
C
C Tests for the decay Chi_J --> J/psi+gamma ---------------
      GOTO 3
      vjtest=DOT(vj,vj)-xpsi*xpsi
      vgtest=DOT(vg,vg)
      test1=vj(1)+vg(1)-qj(1)
      test2=vj(2)+vg(2)-qj(2)
      test3=vj(3)+vg(3)-qj(3)
      test4=vj(4)+vg(4)-qj(4)
      IF(dabs(vjtest).GE.1.D-5) write(6,*) 'Bad J/psi mass'
      IF(dabs(vgtest).GE.1.D-5) write(6,*) 'Bad photo mass'
      IF(dabs(test1).GE.1.D-5) write(6,*) 'Bad component 1'
      IF(dabs(test2).GE.1.D-5) write(6,*) 'Bad component 2'
      IF(dabs(test3).GE.1.D-5) write(6,*) 'Bad component 3'
      IF(dabs(test4).GE.1.D-5) write(6,*) 'Bad component 4'
    3 CONTINUE !End of test section -----------------------
C
C... Construct the reference ort system for J/psi decays
      CALL JREST(xpsi,qpsi,b1,b2,f1,f2)
      CALL JFRAME(JF2,f1,f2,OX2,OY2,OZ2)
C
C... Generate two-body leptonic decay  J/psi --> mu+ mu-
      vmu1(4)=xpsi/2.
      vmu2(4)=xpsi/2.
      vtot=dsqrt(xpsi2/4.-xmu2)
c    Two Euler's angles at random
      cosb=-1. + 2.*RANDOM(NUM)
      alph= 3.14159*RANDOM(NUM)*2.
      sinb= dabs(1.-cosb*cosb)
      sinb= dsqrt(sinb)
      cosa= dcos(alph)
      sina= dsin(alph)
      vmu1(1)= vtot*cosb
      vmu1(2)= vtot*sinb*cosa
      vmu1(3)= vtot*sinb*sina
      do 5 i=1,3
    5 vmu2(i)=-vmu1(i)
C    Angles with respect to the predefined orts
      csm = DOT3(vmu1,OZ2)/vtot
      cph = DOT3(vmu1,OX2)
      sph = DOT3(vmu1,OY2)
      phm = DATAN2(sph,cph)
C... Decay momenta in the beam-beam c.m.s.
      CALL CMSYS(xpsi,qpsi, vmu1,vmu2,qmu1,qmu2)
C... Decay momenta in the chi_J rest frame
      CALL JREST(XJ,pj,qmu1,qmu2,v1,v2)
C
C Tests for the decay J/psi --> muon+muon -----------------
      GOTO 6
      v1test=DOT(v1,v1)-xmu2
      v2test=DOT(v2,v2)-xmu2
      test1=vj(1)-v1(1)-v2(1)
      test2=vj(2)-v1(2)-v2(2)
      test3=vj(3)-v1(3)-v2(3)
      test4=vj(4)-v1(4)-v2(4)
      IF(dabs(v1test).GE.1.D-5) write(6,*) 'Bad muon1 mass'
      IF(dabs(v2test).GE.1.D-5) write(6,*) 'Bad muon2 mass'
      IF(dabs(test1).GE.1.D-5) write(6,*) 'Bad component 1'
      IF(dabs(test2).GE.1.D-5) write(6,*) 'Bad component 2'
      IF(dabs(test3).GE.1.D-5) write(6,*) 'Bad component 3'
      IF(dabs(test4).GE.1.D-5) write(6,*) 'Bad component 4'
    6 CONTINUE !End of test section -----------------------
      RETURN
      END

      SUBROUTINE JMUON(WGTJ)
C********************************************************************C
C     Leptonic decay vertex J/psi --> mu+ + mu-                      C
C********************************************************************C
      DIMENSION SPJ(4,4)
      COMMON/CONST/CME,PI,ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      COMMON/FJPSI/f1(4),f2(4),vmu1(4),vmu2(4)
      COMMON/LOOP1/AMP1(4)
      COMMON/GMUNU/DF(4,4),DC(4),EP(4,4,4,4)
      COMMON/XFACT/CXFACT,COLOR
      SAVE
      WGTJ=0.D0
C... Jpsi decay spin density matrix
      DO 10 mj=1,4 
      DO 10 nj=1,4
      SPJ(mj,nj)= 3.*(vmu1(mj)*vmu2(nj) +vmu2(mj)*vmu1(nj) 
     .                                      -XJ2*DF(mj,nj)/2.)/XJ2
   10 CONTINUE
C... Convoluting with the production matrix elements
      DO 12 mj=1,4
      DO 12 nj=1,4
      WGTJ = WGTJ +AMP1(mj)*AMP1(nj)*SPJ(mj,nj)
   12 CONTINUE
      WGTJ = WGTJ*CXFACT*COLOR
      RETURN
      END

      SUBROUTINE JDECAY
C********************************************************************C
C     Generates two-body leptonic decay  J/psi --> mu+ mu-           C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/FJPSI/f1(4),f2(4),vmu1(4),vmu2(4)
      COMMON/FLABS/qpsi(4),qgam(4),qmu1(4),qmu2(4)
      COMMON/MOMEN/b1(4),b2(4),g1(4),g2(4),pj(4),a1(4),a2(4)
      COMMON/CONST/CME,PI,ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      COMMON/ORTS/OX1(4),OY1(4),OZ1(4),OX2(4),OY2(4),OZ2(4)
      COMMON/FREST/d1(4),d2(4),q1(4),q2(4),qj(4),e1(4),e2(4)
      COMMON/DECAY/xpsi,xmu2,vj(4),vg(4),v1(4),v2(4),csg,phg,csm,phm
      COMMON/FRAME/JF1,JF2
      COMMON/SEED/NUM
      SAVE
      xmu2=xmu*xmu
C... Particle momenta in the chi_J rest frame
      CALL JREST(XJ,pj,b1,b2,d1,d2)
      CALL JREST(XJ,pj,g1,g2,q1,q2)
      CALL JREST(XJ,pj,a1,a2,e1,e2)
      DO 1 i=1,3  
    1 qj(i)=0.D0
      qj(4)=XJ  
C... Construct the reference ort system for J/psi decays
      CALL JREST(XJ,pj,b1,b2,f1,f2)
      CALL JFRAME(JF2,f1,f2,OX2,OY2,OZ2)
C
C... Generate two-body leptonic decay  J/psi --> mu+ mu-
      vmu1(4)=XJ/2.
      vmu2(4)=XJ/2.
      vtot=dsqrt(XJ2/4.-xmu2)
c    Two Euler's angles at random
      cosb=-1. + 2.*RANDOM(NUM)
      alph=-PI + 2.*RANDOM(NUM)*PI
      sinb= dabs(1.-cosb*cosb)
      sinb= dsqrt(sinb)
      cosa= dcos(alph)
      sina= dsin(alph)
      DO 5 i=1,3
      vmu1(i)=(OX2(i)*sinb*cosa +OY2(i)*sinb*sina +OZ2(i)*cosb)*vtot
    5 vmu2(i)=-vmu1(i)
      csm = cosb
      phm = alph
C... Decay momenta in the beam-beam c.m.s.
      CALL CMSYS(XJ,pj, vmu1,vmu2,qmu1,qmu2)
C
C Tests for the decay J/psi --> muon+muon -----------------
c     GOTO 6
      v1test=DOT(qmu1,qmu1)-xmu2
      v2test=DOT(qmu2,qmu2)-xmu2
      test1=pj(1)-qmu1(1)-qmu2(1)
      test2=pj(2)-qmu1(2)-qmu2(2)
      test3=pj(3)-qmu1(3)-qmu2(3)
      test4=pj(4)-qmu1(4)-qmu2(4)
      IF(dabs(v1test).GE.1.D-5) write(6,*) 'Bad muon1 mass'
      IF(dabs(v2test).GE.1.D-5) write(6,*) 'Bad muon2 mass'
      IF(dabs(test1).GE.1.D-5) write(6,*) 'Bad component 1'
      IF(dabs(test2).GE.1.D-5) write(6,*) 'Bad component 2'
      IF(dabs(test3).GE.1.D-5) write(6,*) 'Bad component 3'
      IF(dabs(test4).GE.1.D-5) write(6,*) 'Bad component 4'
    6 CONTINUE !End of test section -----------------------
C
      RETURN
      END

      SUBROUTINE EDECAY
C********************************************************************C
C     Generates two-body leptonic decay  Eta --> mu+ mu-           C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/FJPSI/f1(4),f2(4),vmu1(4),vmu2(4)
      COMMON/FLABS/qpsi(4),qgam(4),qmu1(4),qmu2(4)
      COMMON/MOMEN/b1(4),b2(4),g1(4),g2(4),pj(4),a1(4),a2(4)
      COMMON/CONST/CME,PI,ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      COMMON/ORTS/OX1(4),OY1(4),OZ1(4),OX2(4),OY2(4),OZ2(4)
      COMMON/FREST/d1(4),d2(4),q1(4),q2(4),qj(4),e1(4),e2(4)
      COMMON/DECAY/xpsi,xmu2,vj(4),vg(4),v1(4),v2(4),csg,phg,csm,phm
      COMMON/FRAME/JF1,JF2
      COMMON/SEED/NUM 
      SAVE
      xmu2=xmu*xmu
C... Particle momenta in the chi_J rest frame
      CALL JREST(XJ,pj,b1,b2,d1,d2)
      CALL JREST(XJ,pj,g1,g2,q1,q2)
      CALL JREST(XJ,pj,a1,a2,e1,e2)
      DO 1 i=1,3
    1 qj(i)=0.D0
      qj(4)=XJ
      RETURN
      END

      SUBROUTINE JREST(XJ,pj,p1,p2,d1,d2)
C********************************************************************C
C      Recalculate momenta p1,p2 into pj rest frame making d1,d2     C
C ___________________________________________________________________C
C  The component notation: p(1, 2, 3, 4) = (Beam, p_T, p_T, Energy)  C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION pj(4),pjr(4),pjb(4),qj(4)
      DIMENSION p1(4),p1r(4),p1b(4),d1(4)
      DIMENSION p2(4),p2r(4),p2b(4),d2(4)
      DATA Tiny/1.D-5/
C
      pjtot = dsqrt(pj(2)**2 + pj(3)**2 + pj(1)**2)
      pjt   = dsqrt(pj(2)**2 + pj(3)**2)
C
C... Rotation around the beam axis (get pj(3) =0)
      cphi = pj(2)/pjt
      sphi = pj(3)/pjt
c     CALL ROTAT(cphi,sphi,2,3, pj,pjr)
      CALL ROTAT(cphi,sphi,2,3, p1,p1r)
      CALL ROTAT(cphi,sphi,2,3, p2,p2r)
c     IF(dabs(pjr(3)).GE.1.E-5) write(6,*) 'JREST: wrong pj(3)'
C
C... Rotation in the pj-production plane (get pj(2)= 0)
      cost = pj(1)/pjtot
      sint = pjt/pjtot
c     CALL ROTAT(cost,sint,1,2, pjr,pjb)
      CALL ROTAT(cost,sint,1,2, p1r,p1b)
      CALL ROTAT(cost,sint,1,2, p2r,p2b)
c     IF(dabs(pjb(2)).GE.1.E-5) write(6,*) 'JREST: wrong pj(2)'
C
C... Boost from the Lab system to pj rest (get pj(1)= 0)
      Gam   = pj(4)/XJ
      Gambet= pjtot/XJ
c     CALL BOOST(Gam,Gambet,1,4, pjb,qj)
      CALL BOOST(Gam,Gambet,1,4, p1b,d1)
      CALL BOOST(Gam,Gambet,1,4, p2b,d2)
c     IF(dabs(qj(1)).GE.1.E-5) write(6,*) 'JREST: wrong pj(1)'
C
      RETURN
      END

      SUBROUTINE JFRAME(JF,p1,p2,OX,OY,OZ)
C********************************************************************C
C                                                                    C
C   Define the reference coordinate system according to JF choice    C
C    p1  is assumed to be the Photon_1 momentum (in pJ rest frame)   C
C    p2  is assumed to be the Proton_2 momentum (in pJ rest frame)   C
C ___________________________________________________________________C
C  The component notation: p(1, 2, 3, 4) = (Beam, p_T, p_T, Energy)  C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION p1(4),p2(4), OX(4),OY(4),OZ(4)
C
      ptot1= DSQRT(p1(1)**2 + p1(2)**2 + p1(3)**2)
      ptot2= DSQRT(p2(1)**2 + p2(2)**2 + p2(3)**2)
C
      GOTO(10,20,30,40,50),JF !Choose the reference frame
C---  Recoil frame ----------
   10 DO 11 i=1,3
   11 OZ(i) =-p1(i)-p2(i)
      GOTO 90
C---  Gottfried-Jackson frame
   20 DO 21 i=1,3
   21 OZ(i) = p1(i)
      GOTO 90
C---  Target frame ----------
   30 DO 31 i=1,3
   31 OZ(i) =-p2(i)
      GOTO 90
C---  Collins-Soper frame ---
   40 DO 41 i=1,3
   41 OZ(i) = p1(i)/ptot1 -p2(i)/ptot2
      GOTO 90
C---  Perpendicular frame ---
   50 DO 51 i=1,3
   51 OZ(i) = p1(i)/ptot1 +p2(i)/ptot2
      GOTO 90
C... End of frame definitions
C    Construct the coordinate system orts
C
   90 OZ(4)= DSQRT(OZ(1)**2 + OZ(2)**2 + OZ(3)**2)
      CALL VEC3(p1,p2,OY)
      CALL VEC3(OY,OZ,OX)
      DO 91 i=1,3
      OZ(i) = OZ(i)/OZ(4)
      OY(i) = OY(i)/OY(4)
   91 OX(i) = OX(i)/OX(4)
      OZ(4) = 0.D0
      OY(4) = 0.D0
      OX(4) = 0.D0
C
C... Tests
c     GOTO 99
      testX  = DOT(oX,oX) +1.
      testY  = DOT(oY,oY) +1.
      testZ  = DOT(oZ,oZ) +1.
      testXY = DOT(oX,oY)
      testXZ = DOT(oX,oZ)
      testZY = DOT(oZ,oY)
      IF(dabs(oX(4)).GE.1.E-4) write(6,*) 'bad oX(4)', oX
      IF(dabs(oY(4)).GE.1.E-4) write(6,*) 'bad oY(4)', oY
      IF(dabs(oZ(4)).GE.1.E-4) write(6,*) 'bad oZ(4)', oZ
      IF(dabs(testX).GE.1.E-4) write(6,*) 'bad oX norm', testX
      IF(dabs(testY).GE.1.E-4) write(6,*) 'bad oY norm', testY
      IF(dabs(testZ).GE.1.E-4) write(6,*) 'bad oZ norm', testZ
      IF(dabs(testXY).GE.1.E-4) write(6,*) 'bad oXoY', testXY
      IF(dabs(testXZ).GE.1.E-4) write(6,*) 'bad oXoZ', testXZ
      IF(dabs(testZY).GE.1.E-4) write(6,*) 'bad oZoY', testZY
   99 RETURN
      END

      SUBROUTINE CMSYS(XJ,pj,d1,d2,p1,p2)
C********************************************************************C
C   Recalculate momenta d1,d2 from pj rest frame into c.m.s. p1,p2   C
C ___________________________________________________________________C
C  The component notation: p(1, 2, 3, 4) = (Beam, p_T, p_T, Energy)  C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION pj(4),pjr(4),pjb(4),qj(4)
      DIMENSION p1(4),p1r(4),p1b(4),d1(4)
      DIMENSION p2(4),p2r(4),p2b(4),d2(4)
      DATA Tiny/1.D-5/
C
      pjtot = dsqrt(pj(2)**2 + pj(3)**2 + pj(1)**2)
      pjt   = dsqrt(pj(2)**2 + pj(3)**2)
      qj(1) = 0.D0
      qj(2) = 0.D0
      qj(3) = 0.D0
      qj(4) = XJ
C
C... Boost from pj rest frame to c.m.s.  (restore nominal pj(4))
      Gam   = pj(4)/XJ
      Gambet=-pjtot/XJ
c     CALL BOOST(Gam,Gambet,1,4, qj,pjb)
      CALL BOOST(Gam,Gambet,1,4, d1,p1b)
      CALL BOOST(Gam,Gambet,1,4, d2,p2b)
c     IF(dabs(pjb(4)-pj(4)).GE.1.E-5) write(6,*) 'CMSYS: wrong pj(4)'
C
C... Rotation in the pj-production plane (restore nominal pj(1))
      cost = pj(1)/pjtot
      sint =-pjt/pjtot
c     CALL ROTAT(cost,sint,1,2, pjb,pjr)
      CALL ROTAT(cost,sint,1,2, p1b,p1r)
      CALL ROTAT(cost,sint,1,2, p2b,p2r)
c     IF(dabs(pjr(1)-pj(1)).GE.1.E-5) write(6,*) 'CMSYS: wrong pj(1)'
C
C... Rotation around the beam axis (restore nominal pj(2),pj(3))
      cphi = pj(2)/pjt
      sphi =-pj(3)/pjt
c     CALL ROTAT(cphi,sphi,2,3, pjr,qj)
      CALL ROTAT(cphi,sphi,2,3, p1r,p1)
      CALL ROTAT(cphi,sphi,2,3, p2r,p2)
c     IF(dabs(qj(2)-pj(2)).GE.1.E-5) write(6,*) 'CMSYS: wrong pj(2)'
c     IF(dabs(qj(3)-pj(3)).GE.1.E-5) write(6,*) 'CMSYS: wrong pj(3)'
C
      RETURN
      END

      SUBROUTINE ROTAT(Ct,St,nx,ny, Pin,Pnew)
C     Rotation in the "nx*ny" plane
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Pin(4),Pnew(4)
      Pnew(nx) = Pin(nx)*Ct + Pin(ny)*St
      Pnew(ny) = Pin(ny)*Ct - Pin(nx)*St
      DO 1 i=1,4
      IF(i.EQ.nx .OR. i.EQ.ny) GOTO 1
      Pnew(i) = Pin(i)
    1 CONTINUE
      RETURN
      END

      SUBROUTINE BOOST(Ga,GB,nz,nE, Pin,Pnew)
C     Lorentz boost along the "nz" axis
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Pin(4),Pnew(4)
      Pnew(nz) = Pin(nz)*Ga - Pin(nE)*GB
      Pnew(nE) = Pin(nE)*Ga - Pin(nz)*GB
      DO 1 i=1,4
      IF(i.EQ.nz .OR. i.EQ.nE) GOTO 1
      Pnew(i) = Pin(i)
    1 CONTINUE
      RETURN
      END
      
      SUBROUTINE VEC3(V1,V2,VR)
C     Vector product of the 3-vectors, VR = [V1 x V2]
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      DIMENSION V1(4),V2(4),VR(4)
      VR(1) = V1(2)*V2(3) - V1(3)*V2(2)
      VR(2) = V1(3)*V2(1) - V1(1)*V2(3)
      VR(3) = V1(1)*V2(2) - V1(2)*V2(1)
      VR(4) = DSQRT(VR(1)**2 +VR(2)**2 +VR(3)**2)
      RETURN
      END

      DOUBLE PRECISION FUNCTION DOT3(xx,yy)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION xx(4),yy(4)   !P(i) = (Beam, P_t, P_t, Energy)
      DOT3 = xx(1)*yy(1) + xx(2)*yy(2) + xx(3)*yy(3)
      RETURN
      END
