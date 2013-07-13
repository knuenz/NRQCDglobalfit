C###################################################################C
C                                                                   C
C       Routines to simulate the production of J/psi mesons         C
C                                                                   C
C    in the framework of Semihard Approach and/or Parton Model      C
C      via gluon-gluon, photon-gluon or photon-photon fusion        C
C                                                                   C
C###################################################################C
C
C     Real arithmetics for J/psi helicities
C
C     Options for helicity frames: Target, Recoil, Got-Jac, Col-Sop
C
C     Reduced version: no polarised beams (SYMSYM only)
C                      simple gluon/photon polarization: kt*kt
C
C     Take care of using COMMON/MOMEN/ or COMMON/FREST/  in FEYNJ
C
      DOUBLE PRECISION FUNCTION FXN(X,WGT)
C*******************************************************************C
C      The integrand expression for VEGAS  +  kinematics            C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      PARAMETER(PI=3.1415926D0)
      DIMENSION X(10), WH(4), XPQ1(-6:6),XPQ2(-6:6)
      COMMON/BVEG1/XL(10),XU(10),ACC
      COMMON/BVEGG/NDIM,NCALL,ITMX,NPRN
      COMMON/CONST/CME,ALP,ALS, XC,XC2,XJ,XJ2, QC,PSI
      COMMON/MOMEN/p1(4),p2(4),g1(4),g2(4),pj(4),g3(4)
      COMMON/KINEM/x1,x2,g1T,g2T, pjT,yj,zj
      COMMON/DELTA/dph,y3
      COMMON/OUTPX/CXTOT,CX000,FACTOR
      COMMON/TYPE/ITYPE,MODEL,IB1,IB2
      COMMON/SEED/NUM
      SAVE
      JFXN=0
      FXN=0.D0
C... Random numbers within VEGAS
       pjT=dexp(X(1))   !Transverse momentum of J/psi
       yj   =   X(2)    !Rapidity of J/psi
       y3   =   X(3)    !Rapidity of accompanying Gluon or Photon
      IF(NDIM.GE.4) g1T=dexp(X(4))  !pT of the 1st initial parton
      IF(NDIM.GE.5) g2T=dexp(X(5))  !pT of the 2nd initial parton
      IF(NDIM.LT.4) g1T = 0.001
      IF(NDIM.LT.5) g2T = 0.001
C
C... Random numbers other than VEGAS
      ph1 = 2.*PI*RANDOM(NUM)-PI !Azm.angle of the 1st parton
      ph2 = 2.*PI*RANDOM(NUM)-PI !Azm.angle of the 2nd parton
      phj = 2.*PI*RANDOM(NUM)-PI !Azimuthal angle of J/psi
      dph = ph2-ph1
      IF(dph.LT.0.D0) dph= dph + 2.*PI
      IF(dph.GE.PI)   dph= 2.*PI - dph
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
      g1(2) = g1T*dcos(ph1)
      g1(3) = g1T*dsin(ph1)
      g2(2) = g2T*dcos(ph2)
      g2(3) = g2T*dsin(ph2)
      pj(2) = pjT*dcos(phj)
      pj(3) = pjT*dsin(phj)
      tmj= dsqrt(XJ2 +pjT*pjT)
      g3(2) = g1(2) +g2(2) -pj(2)
      g3(3) = g1(3) +g2(3) -pj(3)
      tm3= dsqrt(g3(2)**2 +g3(3)**2)
      pj(1) = tmj*dsinh(yj)
      pj(4) = tmj*dcosh(yj)
      g3(1) = tm3*dsinh(y3)
      g3(4) = tm3*dcosh(y3)
      x1 = (tmj*dexp(yj) + tm3*dexp(y3))/CME
      x2 = (tmj*dexp(-yj)+ tm3*dexp(-y3))/CME
      IF(x1.GE..999D0 .OR. x2.GE..999D0) RETURN
      IF(x1.LE. 1.D-6 .OR. x2.LE. 1.D-6) RETURN
      y1 = 15.D0
      y2 =-15.D0
      DO 3 k=1,7  !Iterative solution
      IF(dabs(y1).GE.20.) RETURN
      IF(dabs(y2).GE.20.) RETURN
      arg1= CME*(1.-x1)-g2T*dexp(y2)
      arg2= CME*(1.-x2)-g1T*dexp(-y1)
      y1 = dlog(arg1/g1T)
      y2 =-dlog(arg2/g2T)
    3 CONTINUE
      IF(y1.LE. 0.D0) RETURN
      IF(y2.GE. 0.D0) RETURN
      g1(1) = CME/2. -g1T*dsinh(y1)
      g1(4) = CME/2. -g1T*dcosh(y1)
      g2(1) =-CME/2. -g2T*dsinh(y2)
      g2(4) = CME/2. -g2T*dcosh(y2)
       ptest1=pj(1)+g3(1)-g1(1)-g2(1)
       ptest4=pj(4)+g3(4)-g1(4)-g2(4)
      IF(dabs(ptest1).GE. 1.D-5) RETURN
      IF(dabs(ptest4).GE. 1.D-5) RETURN
C
C... Infrared cutoffs
c     RedCut= 0.25*XJ2         ![GeV2]
c     RedCut= 1.D0
c     RedCut= 0.D0             !Singlet only
c     CALL CUTOFF(RedCut, infra)
c     IF(infra .EQ. 0) RETURN
C... Cut on the inelasticity
c     zj = DOT(pj,p2)/DOT(g1,p2) -0.0001
c     IF(zj .LE. 0.30D0) RETURN
c     IF(zj .GT. 0.90D0) RETURN
C... Cut on the final state hadron mass
c     DO 4 i=1,4
c   4 WH(i) = g3(i)+p2(i)-g2(i)
c     WX=DOT(WH,WH)
c     IF(WX.LT.100.) RETURN
C
C... Phase Space boundary and flux factor
      G1W = DOT(g1,g1)
      G2W = DOT(g2,g2)
      STOT= CME*CME
       SXX = XJ2 +2.*DOT(pj,g3)
       SXg = G2W +2.*DOT(p1,g2)
      IF(AF2(SXX,G1W,G2W).LE.0.D0) RETURN
      T1  = G1W +XJ2 -2.*DOT(g1,pj)
      T2  = G2W      -2.*DOT(g2,g3)
      IF(GF(SXX,T1,0.D0,G1W,G2W,XJ2).GE.0.D0) RETURN
      IF(GF(SXX,T2,XJ2,G2W,G1W,0.D0).GE.0.D0) RETURN
C
       SXX = x1*x2*STOT
       SXg =    x2*STOT
      IF((IB1*IB2).EQ.64) DENOMR=    STOT*SXX/PI/PI/4.
      IF((IB1*IB2).EQ.8) DENOMR=     STOT*SXg*G1W**2 /PI
      IF((IB1*IB2).EQ.1) DENOMR= 4.*(STOT*G1W*G2W)**2
C
      IF(NDIM.EQ.3) FACTOR=        (2.*(pjT)**2)/DENOMR
      IF(NDIM.EQ.4) FACTOR=    (4.*(g1T*pjT)**2)/DENOMR
      IF(NDIM.EQ.5) FACTOR=(8.*(g1T*g2T*pjT)**2)/DENOMR
C
C... Couplings and other constants; averaging over spins and colors
c     qR2= SXX/4.  !Renormalization scale in the strong coupling
      qR2= tmj*tmj !
      XNF= 4.      ! Number of flavours
      ALS= ALPHAS(qR2,XNF)
c     AL1= ALPHAS(dabs(G1W),XNF)
c     AL2= ALPHAS(dabs(G2W),XNF)
      IF(ITYPE.EQ.0) FACTOR = FACTOR*PSI *ALP    *ALS**2 /4./64.
      IF(ITYPE.EQ.1) FACTOR = FACTOR*PSI         *ALS**3 /4./64.
c>    IF(ITYPE.EQ.1) FACTOR = FACTOR*PSI *AL1*AL2*ALS    /4./64.
      IF(ITYPE.EQ.2) FACTOR = FACTOR*PSI *ALP**2 *ALS**2 /4./8.
c>    IF(ITYPE.EQ.2) FACTOR = FACTOR*PSI *ALP**2 *ALS*AL2/4./8.
      IF(ITYPE.EQ.3) FACTOR = FACTOR*PSI *ALP**5         /4.
      IF(ITYPE.EQ.4) FACTOR = FACTOR*PSI *ALP    *ALS**2 /4./64.
      IF(ITYPE.EQ.5) FACTOR = FACTOR*PSI         *ALS**3 /4./64.
      IF(ITYPE.EQ.6) FACTOR = FACTOR*PSI *ALP**2 *ALS**2 /4./8.
      IF(ITYPE.EQ.7) FACTOR = FACTOR*PSI *ALP**2 *ALS**2 /4./8.
c     IF(ITYPE.EQ.7) FACTOR = FACTOR*PSI *ALP**2 *ALS*AL2/4./8.
      IF(ITYPE.EQ.8) FACTOR = FACTOR*PSI *ALP**4 *ALS    /4.
C
C     FACTOR=FACTOR*SXX/AF(SXX,G1W,G2W)
C
C... Gluon and Photon Distribution Functions
c     qF2= SXX/4.  !Factorization scale in parton distributions
      qF2= tmj*tmj !
      qF = tmj

c     qTpair2 = (g1(2)+g2(2))**2 + (g1(3)+g2(3))**2
c     qF2 = SXX + qTpair2
c     qF  = dsqrt(qF2)

      w1T = g1T*g1T
      w2T = g2T*g2T
C    
      IF(IB1.EQ.8) THEN
                   CALL cauniglu(2212,x1,w1T,qF,XPQ1)
                   PDF1 = XPQ1(0)/x1
                   ELSE
                   PDF1 = 1.D0
      ENDIF
      IF(IB2.EQ.8) THEN
                   CALL cauniglu(2212,x2,w2T,qF,XPQ2)
                   PDF2 = XPQ2(0)/x2
                   ELSE
                   PDF2 = 1.D0
      ENDIF
C
c     LFs= 0       !Factorization scheme: 0,1,2=LO,NLO(MS),NLO(DIS) 
c     IF(MODEL.EQ.0 .AND.IB1.EQ.8) PDF1 =GRVg(x1,qF2,LFs)
c     IF(MODEL.EQ.1 .AND.IB1.EQ.8) PDF1 =BLUEM(x1,w1T,qF2,LFs,ALS)
c     IF(MODEL.EQ.2 .AND.IB1.EQ.8) PDF1 =DGLAP(x1,w1T,qF2,LFs)
c                     IF(IB1.EQ.8) PDF1 =PDF1/x1
c                     IF(IB1.EQ.1) PDF1 =1.D0
c     IF(MODEL.EQ.0 .AND.IB2.EQ.8) PDF2 =GRVg(x2,qF2,LFs)
c     IF(MODEL.EQ.1 .AND.IB2.EQ.8) PDF2 =BLUEM(x2,w2T,qF2,LFs,ALS)
c     IF(MODEL.EQ.2 .AND.IB2.EQ.8) PDF2 =DGLAP(x2,w2T,qF2,LFs)
c                     IF(IB2.EQ.8) PDF2 =PDF2/x2
c                     IF(IB2.EQ.1) PDF2 =1.D0
C
      FACTOR =FACTOR*PDF1*PDF2*16.*XC2
      IF(PDF1.LT.0.D0) write(6,*) 'PDF1 is negative',PDF1
      IF(PDF2.LT.0.D0) write(6,*) 'PDF2 is negative',PDF2
      IF(PDF1.LT.0.D0) FACTOR=0.D0
      IF(PDF2.LT.0.D0) FACTOR=0.D0
C
C... Partonic cross section
      CALL XSEC(SYMSYM,SYM000)
      CXTOT = FACTOR*SYMSYM
      CX000 = FACTOR*SYM000
      IF(CXTOT.LT.0.D0) THEN
        write(6,*) 'Negative total cross section !!!!!!!!!'
        CXTOT=0.D0
      ENDIF
      IF(CXTOT.GE.1.D-20 .AND. CXTOT.LT.1.D20) JFXN=1
      IF(JFXN.EQ.0) RETURN
      FXN=CXTOT *(1.+pjT**2)**2 !to enlarge statistics at high pjT
C
C... Filling histograms
      IF(NPRN.NE.1) RETURN
c     hfxn = sngl(WGT*CXTOT)/FLOAT(ITMX)
      hwgt = sngl(WGT)/FLOAT(ITMX)
      CALL WRIOUT(hwgt)
      RETURN
      END

      SUBROUTINE XSEC(SYMSYM,SYM000)
C*******************************************************************C
C      Partonic differential cross section                          C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      Parameter(GeVmb=0.4D0, GeVmub=4.D2, GeVnb=4.D5, GeVpb=4.D8)
      COMMON/TYPE/ITYPE,MODEL,IB1,IB2
      COMMON/CONST/CME,ALP,ALS, XC,XC2,XJ,XJ2, QC,PSI
      COMMON/FREST/d1(4),d2(4),q1(4),q2(4),qj(4),q3(4),e1(4),e2(4)
      COMMON/MOMEN/p1(4),p2(4),g1(4),g2(4),pj(4),g3(4)
      COMMON/LOOP8/AMP8(4,4,9)
      COMMON/LOOPJ/AMPJ(4,4)
      COMMON/SPINJ/SPJ(4,4),SPS(4,4)
      COMMON/GMUNU/DF(4,4),DC(4)
      SAVE
C...Momentum dot-products
      CALL SCALAR
C...Gluon polarization vectors
      CALL GAUGEG
C...Go to the meson rest frame
c     CALL JREST         !
C...Define the reference axes
c     CALL JFRAME(d1,d2) !
C...Meson polarization vectors
c     CALL GAUGEF(1) !Explicit  = scalar + longitudinal
      CALL GAUGEF(2) !Invariant = g_mu_nu -(p_mu*p_nu)/XJ2
C...Evaluating the Matrix Elements
      CALL FEYNJ
      SYMSYM=0.D0
      SYM000=0.D0
      COLOR=0.D0    !Color and charge coefficients 
      IF(ITYPE.EQ.0) COLOR= 8./3. *QC    ! glu*glu -> Jpsi*gam
      IF(ITYPE.EQ.1) COLOR=10./9.        ! glu*glu -> Jpsi*glu
      IF(ITYPE.EQ.2) COLOR= 8./3. *QC    ! gam*glu -> Jpsi*glu
      IF(ITYPE.EQ.3) COLOR=   12. *QC**3 ! gam*gam -> Jpsi*gam
      IF(ITYPE.EQ.4) COLOR= 5./6. *QC    ! glu*glu -> Jpsi8*gam
      IF(ITYPE.EQ.5) COLOR= 1.           ! glu*glu -> Jpsi8*glu
      IF(ITYPE.EQ.6) COLOR= 2.    *QC**2 ! gam*glu -> Jpsi8*gam
      IF(ITYPE.EQ.7) COLOR= 5./6. *QC    ! gam*glu -> Jpsi8*glu
      IF(ITYPE.EQ.8) COLOR= 2.    *QC**2 ! gam*gam -> Jpsi8*glu
      DO 12  l=1,4
      DO 11 mj=1,4
      DO 11 nj=1,4
      IF(ITYPE.EQ.5) GOTO 5
      SYMSYM=SYMSYM-AMPJ(mj,l)*AMPJ(nj,l)*SPJ(mj,nj)*DC(l)
c     SYM000=SYM000-AMPJ(mj,l)*AMPJ(nj,l)*SPS(mj,nj)*DC(l)
      IF(ITYPE.NE.5) GOTO 10
   5  ANS0=
     .  AMP8(mj,l,1)*AMP8(nj,l,1)*72.
     . +AMP8(mj,l,2)*AMP8(nj,l,2)*200./9.
     . +AMP8(mj,l,3)*AMP8(nj,l,3)*72.
     . +AMP8(mj,l,4)*AMP8(nj,l,4)*200./9.
     . +AMP8(mj,l,5)*AMP8(nj,l,5)*72.
     . +AMP8(mj,l,6)*AMP8(nj,l,6)*200./9.
     . +AMP8(mj,l,7)*AMP8(nj,l,7)*64.
     . +AMP8(mj,l,8)*AMP8(nj,l,8)*64.
     . +AMP8(mj,l,9)*AMP8(nj,l,9)*64.
     . +AMP8(mj,l,7)*AMP8(nj,l,8)*16.
     . +AMP8(mj,l,7)*AMP8(nj,l,9)*16.
     . +AMP8(mj,l,8)*AMP8(nj,l,9)*16.
      ANS1=
     .  AMP8(mj,l,1)*AMP8(nj,l,3)*72.
     . -AMP8(mj,l,1)*AMP8(nj,l,4)*40.
     . +AMP8(mj,l,1)*AMP8(nj,l,5)*72.
     . +AMP8(mj,l,1)*AMP8(nj,l,6)*40.
     . +AMP8(mj,l,1)*AMP8(nj,l,7)*48.
     . -AMP8(mj,l,1)*AMP8(nj,l,9)*48.
     . -AMP8(mj,l,2)*AMP8(nj,l,3)*40.
     . -AMP8(mj,l,2)*AMP8(nj,l,4)*40./3.
     . +AMP8(mj,l,2)*AMP8(nj,l,5)*40.
     . -AMP8(mj,l,2)*AMP8(nj,l,6)*40./3.
     . +AMP8(mj,l,2)*AMP8(nj,l,7)*80./3.
     . +AMP8(mj,l,2)*AMP8(nj,l,9)*80./3.
      ANS2=
     . -AMP8(mj,l,3)*AMP8(nj,l,5)*72.
     . +AMP8(mj,l,3)*AMP8(nj,l,6)*40.
     . +AMP8(mj,l,3)*AMP8(nj,l,7)*48.
     . -AMP8(mj,l,3)*AMP8(nj,l,8)*48.
     . -AMP8(mj,l,4)*AMP8(nj,l,5)*40.
     . -AMP8(mj,l,4)*AMP8(nj,l,6)*40./3.
     . +AMP8(mj,l,4)*AMP8(nj,l,7)*80./3.
     . +AMP8(mj,l,4)*AMP8(nj,l,8)*80./3.
     . +AMP8(mj,l,5)*AMP8(nj,l,8)*48.
     . -AMP8(mj,l,5)*AMP8(nj,l,9)*48.
     . +AMP8(mj,l,6)*AMP8(nj,l,8)*80./3.
     . +AMP8(mj,l,6)*AMP8(nj,l,9)*80./3.
      SYMSYM=SYMSYM - (ANS0 +ANS1 +ANS2)*SPJ(mj,nj)*DC(l)
c     SYM000=SYM000 - (ANS0 +ANS1 +ANS2)*SPS(mj,nj)*DC(l)
  10  CONTINUE
  11  CONTINUE
  12  CONTINUE
      SYMSYM=SYMSYM*COLOR*GeVnb !Conversion to nanobarns
      SYM000=SYM000*COLOR*GeVnb
      RETURN
      END

      SUBROUTINE SCALAR
C*******************************************************************C
C      Dot-products of particle momenta                             C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      COMMON/MOMEN/p1(4),p2(4),g1(4),g2(4),pj(4),g3(4)
      COMMON/DOTPR/DG1P,DG2P,DG3P,DG12,DG13,DG23,G1W,G2W
      SAVE
      DG1P=DOT(g1,pj)
      DG2P=DOT(g2,pj)
      DG3P=DOT(g3,pj)
      DG12=DOT(g1,g2)
      DG13=DOT(g1,g3)
      DG23=DOT(g2,g3)
       G1W=DOT(g1,g1)
       G2W=DOT(g2,g2)
      RETURN
      END

      SUBROUTINE GAUGEF(Mode)
C********************************************************************C
C       Spin density of the outgoing J/psi                           C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      DIMENSION EJ(4,-1:1)
      COMMON/CONST/CME,ALP,ALS, XC,XC2,XJ,XJ2, QC,PSI
      COMMON/MOMEN/p1(4),p2(4),g1(4),g2(4),pj(4),g3(4)
c     COMMON/FREST/p1(4),p2(4),g1(4),g2(4),pj(4),g3(4),a1(4),a2(4)
      COMMON/ORTS/OX(4),OY(4),OZ(4)
      COMMON/SPINJ/SPJ(4,4),SPS(4,4)
      COMMON/GMUNU/DF(4,4),DC(4)
      SAVE
      GOTO(10,20),Mode
   10 sq2=sqrt(2.)
C... Longitudinal polarization, = zero helicity
      DO 11 j=1,3
   11 EJ(j, 0)= OZ(j)
      EJ(4, 0)= 0.D0
C... Transverse polarization, Y direction
      DO 12 j=1,3
   12 EJ(j,-1)= OY(j)
      EJ(4,-1)= 0.D0
C... Transverse polarization, X direction
      DO 13 j=1,3
   13 EJ(j, 1)= OX(j)
      EJ(4, 1)= 0.D0
C... Spin density matrix
      DO 14 mj=1,4
      DO 14 nj=1,4
      SPS(mj,nj)=EJ(mj, 0)*EJ(nj, 0)
      SPJ(mj,nj)=EJ(mj,-1)*EJ(nj,-1)
     .          +EJ(mj, 0)*EJ(nj, 0)
     .          +EJ(mj, 1)*EJ(nj, 1)
   14 CONTINUE
      GOTO 99
C... Invariant definition of Spin density matrix
   20 DO 24 mj=1,4
      DO 24 nj=1,4
      SPJ(mj,nj)= -(DF(mj,nj) - pj(mj)*pj(nj)/XJ2)
   24 CONTINUE
   99 RETURN
      END

      SUBROUTINE GAUGEJ
C********************************************************************C
C       Spin density of the outgoing J/psi via its Leptonic decay    C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      COMMON/CONST/CME,ALP,ALS, XC,XC2,XJ,XJ2, QC,PSI
      COMMON/MOMEN/p1(4),p2(4),g1(4),g2(4),pj(4),g3(4)
      COMMON/GMUNU/DF(4,4),DC(4)
      COMMON/MUONS/xmu,v1(4),v2(4),csm,phm
      COMMON/SPINJ/SPJ(4,4),SPS(4,4)
      SAVE
      CALL DECAY2(XJ,pj, xmu,v1, xmu,v2)
      DO 10 mj=1,4
      DO 10 nj=1,4
      SPJ(mj,nj)=3*(v1(mj)*v2(nj)+v2(mj)*v1(nj)-XJ2*DF(mj,nj)/2.)/XJ2
   10 CONTINUE
      RETURN
      END

      SUBROUTINE GAUGEG
C*******************************************************************C
C      Photon and Gluon polarization vectors                        C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      COMMON/MOMEN/p1(4),p2(4),g1(4),g2(4),pj(4),g3(4)
      COMMON/POLAR/            a1(4),a2(4)
      COMMON/KINEM/x1,x2,g1T,g2T, pjT,yj,zj
      COMMON/GMUNU/DF(4,4),DC(4)
      COMMON/TYPE/ITYPE,MODEL,IB1,IB2
      SAVE
      sq=sqrt(2.)
      a1(1)=0.D0
      a1(4)=0.D0
      a2(1)=0.D0
      a2(4)=0.D0
      IF(IB1.EQ.1) Fnorm1 = 2*sq/x1
      IF(IB1.EQ.8) Fnorm1 =   sq/g1T
      IF(IB2.EQ.1) Fnorm2 = 2*sq/x2
      IF(IB2.EQ.8) Fnorm2 =   sq/g2T
      DO 2 i=2,3
      a1(i)=g1(i)*Fnorm1
      a2(i)=g2(i)*Fnorm2
    2 Continue
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
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      COMMON/TYPE/ITYPE,MODEL,IB1,IB2
      COMMON/CONST/CME,ALP,ALS, XC,XC2,XJ,XJ2, QC,PSI
      COMMON/MOMEN/p1(4),p2(4),g1(4),g2(4),pj(4),g3(4)
      COMMON/POLAR/            a1(4),a2(4)
c     COMMON/FREST/p1(4),p2(4),g1(4),g2(4),pj(4),g3(4),a1(4),a2(4)
      COMMON/DOTPR/DG1P,DG2P,DG3P,DG12,DG13,DG23,W1,W2
      COMMON/LOOP8/AMP8(4,4,9)
      COMMON/LOOPJ/AMPJ(4,4)
      COMMON/GMUNU/DF(4,4),DC(4)
      SAVE
      W3=0.D0 !real final gluon
      SHAT=(W1+W2+2*DG12)
      THAT=(W2+W3-2*DG23)
      UHAT=(W1+W3-2*DG13)
C...Fill main array
      DO 1 mj=1,4
      DO 1 m3=1,4
      AMPJ(mj,m3)=0.D0
      DO 1 ic=1,9
      AMP8(mj,m3,ic)=0.D0
   1  CONTINUE
      DO 10 m1=1,4
      DO 10 m2=1,4
      DO 10 m3=1,4
      DO 10 mj=1,4
C...Diagram  =<< 1 >> = << 6 >>
      TERM01=  - DF(m1,m2)*DF(m3,mj)*DG13 - DF(m1,m2)*g1(m3)*g1(mj) +
     + DF(m1,m2)*g1(m3)*g3(mj) - DF(m1,m2)*g1(mj)*g2(m3) + DF(m1,m3)*
     + DF(m2,mj)*DG13 - DF(m1,m3)*g1(m2)*g3(mj) - DF(m1,m3)*g1(mj)*
     + g3(m2) - DF(m1,mj)*DF(m2,m3)*DG13 + DF(m1,mj)*g1(m2)*g1(m3) +
     + DF(m1,mj)*g1(m2)*g2(m3) + DF(m1,mj)*g1(m3)*g3(m2) + DF(m2,m3)*
     + g1(mj)*g3(m1) + DF(m2,m3)*g2(m1)*g3(mj) - DF(m2,m3)*g3(m1)*
     + g3(mj) - DF(m2,mj)*g1(m3)*g2(m1) - DF(m2,mj)*g2(m1)*g2(m3) +
     + DF(m2,mj)*g2(m3)*g3(m1) + DF(m3,mj)*g1(m2)*g3(m1) - DF(m3,mj)*
     + g2(m1)*g3(m2) + DF(m3,mj)*g3(m1)*g3(m2)
      TERM01=TERM01/((W1-DG1P)*(W3+DG3P))*a1(m1)*a2(m2)*DC(m1)*DC(m2)
C...Diagram  << 2 >> = << 5 >>
      TERM02= - DF(m1,m2)*DF(m3,mj)*DG12 + DF(m1,m2)*g1(m3)*g2(mj) +
     + DF(m1,m2)*g1(mj)*g2(m3) + DF(m1,m3)*DF(m2,mj)*DG12 - DF(m1,m3)*
     + g1(m2)*g1(mj) - DF(m1,m3)*g1(m2)*g2(mj) + DF(m1,m3)*g1(mj)*
     + g3(m2) + DF(m1,mj)*DF(m2,m3)*DG12 + DF(m1,mj)*g1(m2)*g1(m3) -
     + DF(m1,mj)*g1(m2)*g2(m3) - DF(m1,mj)*g1(m3)*g3(m2) - DF(m2,m3)*
     + g1(mj)*g2(m1) - DF(m2,m3)*g2(m1)*g2(mj) + DF(m2,m3)*g2(mj)*
     + g3(m1) - DF(m2,mj)*g1(m3)*g2(m1) + DF(m2,mj)*g2(m1)*g2(m3) -
     + DF(m2,mj)*g2(m3)*g3(m1) + DF(m3,mj)*g1(m2)*g3(m1) + DF(m3,mj)*
     + g2(m1)*g3(m2) - DF(m3,mj)*g3(m1)*g3(m2)
      TERM02=TERM02/((W1-DG1P)*(W2-DG2P))*a1(m1)*a2(m2)*DC(m1)*DC(m2)
C...Diagram  << 3 >> = << 4 >>
      TERM03= - DF(m1,m2)*DF(m3,mj)*DG23 - DF(m1,m2)*g1(m3)*g2(mj) -
     + DF(m1,m2)*g2(m3)*g2(mj) + DF(m1,m2)*g2(m3)*g3(mj) - DF(m1,m3)*
     + DF(m2,mj)*DG23 + DF(m1,m3)*g1(m2)*g3(mj) + DF(m1,m3)*g2(mj)*
     + g3(m2) - DF(m1,m3)*g3(m2)*g3(mj) + DF(m1,mj)*DF(m2,m3)*DG23 -
     + DF(m1,mj)*g1(m2)*g1(m3) - DF(m1,mj)*g1(m2)*g2(m3) + DF(m1,mj)*
     + g1(m3)*g3(m2) - DF(m2,m3)*g2(m1)*g3(mj) - DF(m2,m3)*g2(mj)*
     + g3(m1) + DF(m2,mj)*g1(m3)*g2(m1) + DF(m2,mj)*g2(m1)*g2(m3) +
     + DF(m2,mj)*g2(m3)*g3(m1) - DF(m3,mj)*g1(m2)*g3(m1) + DF(m3,mj)*
     + g2(m1)*g3(m2) + DF(m3,mj)*g3(m1)*g3(m2)
      TERM03=TERM03/((W2-DG2P)*(W3+DG3P))*a1(m1)*a2(m2)*DC(m1)*DC(m2)
C
      IF(ITYPE.NE.5) GOTO 5 !Color-octet section starts
C...Diagram << 21 >> = - << 22 >>
      TERM21= - DF(m1,m2)*DF(m3,mj)*DG13 + DF(m1,m2)*DF(m3,mj)*DG23 -
     + DF(m1,m2)*g1(m3)*g1(mj) + DF(m1,m2)*g1(m3)*g2(mj) + DF(m1,m2)*
     + g1(m3)*g3(mj) - DF(m1,m2)*g1(mj)*g2(m3) + DF(m1,m2)*g2(m3)*
     + g2(mj) - DF(m1,m2)*g2(m3)*g3(mj) - 2*DF(m1,m3)*g1(m2)*g3(mj) + 2
     + *DF(m1,mj)*g1(m2)*g1(m3) + 2*DF(m1,mj)*g1(m2)*g2(m3) + 2*
     + DF(m2,m3)*g2(m1)*g3(mj) - 2*DF(m2,mj)*g1(m3)*g2(m1) - 2*
     + DF(m2,mj)*g2(m1)*g2(m3) + 2*DF(m3,mj)*g1(m2)*g3(m1) - 2*
     + DF(m3,mj)*g2(m1)*g3(m2)
      TERM21=TERM21/((W3+DG3P)*SHAT)*a1(m1)*a2(m2)*DC(m1)*DC(m2)
C...Diagram << 23 >> = - << 24 >>
      TERM23= - 2*DF(m1,m2)*g1(mj)*g2(m3) - 2*DF(m1,m3)*g1(mj)*g3(m2) -
     + DF(m1,mj)*DF(m2,m3)*DG12 - DF(m1,mj)*DF(m2,m3)*DG13 + 2*
     + DF(m1,mj)*g1(m2)*g2(m3) + 2*DF(m1,mj)*g1(m3)*g3(m2) + DF(m2,m3)*
     + g1(mj)*g2(m1) + DF(m2,m3)*g1(mj)*g3(m1) + DF(m2,m3)*g2(m1)*
     + g2(mj) + DF(m2,m3)*g2(m1)*g3(mj) - DF(m2,m3)*g2(mj)*g3(m1) -
     + DF(m2,m3)*g3(m1)*g3(mj) - 2*DF(m2,mj)*g2(m1)*g2(m3) + 2*
     + DF(m2,mj)*g2(m3)*g3(m1) - 2*DF(m3,mj)*g2(m1)*g3(m2) + 2*
     + DF(m3,mj)*g3(m1)*g3(m2)
      TERM23=TERM23/((W1-DG1P)*THAT)*a1(m1)*a2(m2)*DC(m1)*DC(m2)
C...Diagram << 25 >> = - << 26 >>
      TERM25= - 2*DF(m1,m2)*g1(m3)*g2(mj) - DF(m1,m3)*DF(m2,mj)*DG12 -
     + DF(m1,m3)*DF(m2,mj)*DG23 + DF(m1,m3)*g1(m2)*g1(mj) + DF(m1,m3)*
     + g1(m2)*g2(mj) + DF(m1,m3)*g1(m2)*g3(mj) - DF(m1,m3)*g1(mj)*
     + g3(m2) + DF(m1,m3)*g2(mj)*g3(m2) - DF(m1,m3)*g3(m2)*g3(mj) - 2*
     + DF(m1,mj)*g1(m2)*g1(m3) + 2*DF(m1,mj)*g1(m3)*g3(m2) - 2*
     + DF(m2,m3)*g2(mj)*g3(m1) + 2*DF(m2,mj)*g1(m3)*g2(m1) + 2*
     + DF(m2,mj)*g2(m3)*g3(m1) - 2*DF(m3,mj)*g1(m2)*g3(m1) + 2*
     + DF(m3,mj)*g3(m1)*g3(m2)
      TERM25=TERM25/((W2-DG2P)*UHAT)*a1(m1)*a2(m2)*DC(m1)*DC(m2)
C...Diagram << 41 + 51 >>
      TERM41 =
     +  - DF(m1,m2)*DF(m3,mj)*W1 + DF(m1,m2)*DF(m3,mj)*W2 + 2*DF(m1,m2)
     + *DF(m3,mj)*DG13 - 2*DF(m1,m2)*DF(m3,mj)*DG23 + DF(m1,m2)*g1(m3)*
     + g1(mj) - 3*DF(m1,m2)*g1(m3)*g2(mj) - DF(m1,m2)*g1(m3)*g3(mj) + 3
     + *DF(m1,m2)*g1(mj)*g2(m3) - DF(m1,m2)*g2(m3)*g2(mj) + DF(m1,m2)*
     + g2(m3)*g3(mj) + 2*DF(m1,m3)*g1(m2)*g1(mj) + 2*DF(m1,m3)*g1(m2)*
     + g2(mj) + 2*DF(m1,m3)*g1(m2)*g3(mj) - 4*DF(m1,mj)*g1(m2)*g1(m3)
     +  - 4*DF(m1,mj)*g1(m2)*g2(m3) - 2*DF(m2,m3)*g1(mj)*g2(m1) - 2*
     + DF(m2,m3)*g2(m1)*g2(mj) - 2*DF(m2,m3)*g2(m1)*g3(mj) + 4*
     + DF(m2,mj)*g1(m3)*g2(m1) + 4*DF(m2,mj)*g2(m1)*g2(m3) - 4*
     + DF(m3,mj)*g1(m2)*g3(m1) + 4*DF(m3,mj)*g2(m1)*g3(m2)
      TERM41=(TERM41/SHAT +DF(m1,mj)*DF(m2,m3)-DF(m1,m3)*DF(m2,mj))/XJ2
     .      *a1(m1)*a2(m2)*DC(m1)*DC(m2)
C...Diagram << 42 + 52 >>
      TERM42 =
     + 2*DF(m1,m2)*g1(mj)*g2(m3) - 2*DF(m1,m2)*g2(m3)*g2(mj) + 2*
     + DF(m1,m2)*g2(m3)*g3(mj) + 2*DF(m1,m3)*g1(mj)*g3(m2) - 2*
     + DF(m1,m3)*g2(mj)*g3(m2) + 2*DF(m1,m3)*g3(m2)*g3(mj) + DF(m1,mj)*
     + DF(m2,m3)*W2 - DF(m1,mj)*DF(m2,m3)*W3 + 2*DF(m1,mj)*DF(m2,m3)*
     + DG12 + 2*DF(m1,mj)*DF(m2,m3)*DG13 - 4*DF(m1,mj)*g1(m2)*g2(m3) -
     + 4*DF(m1,mj)*g1(m3)*g3(m2) - DF(m2,m3)*g1(mj)*g2(m1) - DF(m2,m3)*
     + g1(mj)*g3(m1) - DF(m2,m3)*g2(m1)*g2(mj) - 3*DF(m2,m3)*g2(m1)*
     + g3(mj) + 3*DF(m2,m3)*g2(mj)*g3(m1) + DF(m2,m3)*g3(m1)*g3(mj) + 4
     + *DF(m2,mj)*g2(m1)*g2(m3) - 4*DF(m2,mj)*g2(m3)*g3(m1) + 4*
     + DF(m3,mj)*g2(m1)*g3(m2) - 4*DF(m3,mj)*g3(m1)*g3(m2)
      TERM42=(TERM42/THAT +DF(m1,m2)*DF(mj,m3)-DF(m1,m3)*DF(m2,mj))/XJ2
     .      *a1(m1)*a2(m2)*DC(m1)*DC(m2)
C...Diagram << 43 + 53 >>
      TERM43 =
     +  - 2*DF(m1,m2)*g1(m3)*g1(mj) + 2*DF(m1,m2)*g1(m3)*g2(mj) + 2*
     + DF(m1,m2)*g1(m3)*g3(mj) + DF(m1,m3)*DF(m2,mj)*W1 - DF(m1,m3)*
     + DF(m2,mj)*W3 + 2*DF(m1,m3)*DF(m2,mj)*DG12 + 2*DF(m1,m3)*
     + DF(m2,mj)*DG23 - DF(m1,m3)*g1(m2)*g1(mj) - DF(m1,m3)*g1(m2)*
     + g2(mj) - 3*DF(m1,m3)*g1(m2)*g3(mj) + 3*DF(m1,m3)*g1(mj)*g3(m2)
     +  - DF(m1,m3)*g2(mj)*g3(m2) + DF(m1,m3)*g3(m2)*g3(mj) + 4*
     + DF(m1,mj)*g1(m2)*g1(m3) - 4*DF(m1,mj)*g1(m3)*g3(m2) - 2*
     + DF(m2,m3)*g1(mj)*g3(m1) + 2*DF(m2,m3)*g2(mj)*g3(m1) + 2*
     + DF(m2,m3)*g3(m1)*g3(mj) - 4*DF(m2,mj)*g1(m3)*g2(m1) - 4*
     + DF(m2,mj)*g2(m3)*g3(m1) + 4*DF(m3,mj)*g1(m2)*g3(m1) - 4*
     + DF(m3,mj)*g3(m1)*g3(m2)
      TERM43=(TERM43/UHAT +DF(m1,m2)*DF(mj,m3)-DF(m1,mj)*DF(m2,m3))/XJ2
     .      *a1(m1)*a2(m2)*DC(m1)*DC(m2)
    5 CONTINUE !End of Color-octet section
C
      AMPJ(mj,m3)=AMPJ(mj,m3)+(TERM01+TERM02+TERM03)*DC(mj)
      IF(ITYPE.NE.5) GOTO 8
      AMP8(mj,m3,1)=AMP8(mj,m3,1)+
     + ( TERM01/8. +TERM23/4. -TERM42/4.)*DC(mj)
      AMP8(mj,m3,2)=AMP8(mj,m3,2) +TERM01*DC(mj)/8.
      AMP8(mj,m3,3)=AMP8(mj,m3,3) +
     + (-TERM02/8. +TERM25/4. -TERM43/4.)*DC(mj)
      AMP8(mj,m3,4)=AMP8(mj,m3,4) +TERM02*DC(mj)/8.
      AMP8(mj,m3,5)=AMP8(mj,m3,5) +
     + (-TERM03/8. +TERM21/4. -TERM41/4.)*DC(mj)
      AMP8(mj,m3,6)=AMP8(mj,m3,6) +TERM03*DC(mj)/8.
      AMP8(mj,m3,7)=AMP8(mj,m3,7) +TERM03*DC(mj)/12.
      AMP8(mj,m3,8)=AMP8(mj,m3,8) +TERM01*DC(mj)/12.
      AMP8(mj,m3,9)=AMP8(mj,m3,9) +TERM02*DC(mj)/12.
    8 CONTINUE 
   10 CONTINUE
      RETURN
      END

      SUBROUTINE CUTOFF(RedCut, infra)
C*******************************************************************C
C     Ifrared cutoffs for all denominators in all diagrams          C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      DIMENSION PROP(6)
      COMMON/CONST/CME,ALP,ALS, XC,XC2,XJ,XJ2, QC,PSI
      COMMON/DOTPR/DG1P,DG2P,DG3P,DG12,DG13,DG23,W1,W2
      SAVE
      infra=1
      IF(RedCut.LE.0.D0) RETURN
      CALL SCALAR
      W3=0.D0
      PROP(1)=(W1+W2+2*DG12) !S^hat, part of ZD21
      PROP(2)=(W2+W3-2*DG23) !T^hat, part of ZD23
      PROP(3)=(W1+W3-2*DG13) !U^hat, part of ZD25
      PROP(4)=(W1-DG1P)      !Part of ZD01,ZD02,ZD23
      PROP(5)=(W2-DG2P)      !Part of ZD02,ZD03,ZD25
      PROP(6)=(W3+DG3P)      !Part of ZD03,ZD01,ZD21
      DO 1 i=1,6
    1 IF(DABS(PROP(i)).LT.RedCut) infra=0 
      RETURN
      END

      SUBROUTINE BEAMS
      COMMON/TYPE/ITYPE,MODEL,IB1,IB2
      SAVE
                                                      IB1=8
      IF(ITYPE.EQ.2 .OR. ITYPE.EQ.3. .OR. ITYPE.GE.6) IB1=1
                                     IB2=8
      IF(ITYPE.EQ.3 .OR. ITYPE.EQ.8) IB2=1
      RETURN
      END

      DOUBLE PRECISION FUNCTION ALPHAS(Q2,XNF)
C*******************************************************************C
C     Strong coupling constant                                      C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      ALPHAS=0.D0
      ALAM2 =(0.232)**2
c     ALAM2 =(0.2)**2
      Q2RUN =dabs(Q2)
      IF(Q2RUN.LE.ALAM2) RETURN
      ALPHA =12.*3.1415926/(33.-2.*XNF)/dlog(Q2RUN/ALAM2)
      IF(ALPHA.GE.1.D0) ALPHA=1.D0
      ALPHAS=ALPHA
      RETURN
      END

      DOUBLE PRECISION FUNCTION AF2(X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      AF2 = X*X + Y*Y + Z*Z - 2.*X*Y - 2.*X*Z - 2.*Y*Z
      RETURN
      END

      DOUBLE PRECISION FUNCTION AF(X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      AF = DSQRT(X*X + Y*Y + Z*Z - 2.*X*Y - 2.*X*Z - 2.*Y*Z)
      RETURN
      END

      DOUBLE PRECISION FUNCTION GF(X,Y,Z,U,V,W)
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      GF =          X*Z*W + X*U*V + Y*Z*V + Y*U*W
     . -X*Y*(Z+U+V+W-X-Y) -Z*U*(X+Y+V+W-Z-U) -V*W*(X+Y+Z+U-V-W)
      RETURN
      END

      SUBROUTINE METRIC
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      COMMON/GMUNU/DF(4,4),DC(4)
      SAVE
      DO 1 I=1,4
      DO 1 J=1,4
      DF(I,J)= 0.D0
    1 CONTINUE !P(i) =(Beam, P_t, P_t, Energy)
      DO 2 I=1,3
      DF(I,I)=-1.D0
    2 DC(I)  =-1.D0
      DF(4,4)= 1.D0
      DC(4)  = 1.D0
      RETURN
      END

      DOUBLE PRECISION FUNCTION DOT(xx,yy)
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      DIMENSION xx(4),yy(4)   !P(i) = (Beam, P_t, P_t, Energy)
      DOT= -xx(1)*yy(1) -xx(2)*yy(2) -xx(3)*yy(3) +xx(4)*yy(4)
      RETURN
      END

      DOUBLE PRECISION FUNCTION EPS(q1,q2,q3,q4)
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
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

      SUBROUTINE JREST
C********************************************************************C
C      Recalculate momenta p1,p2 into pj rest frame making d1,d2     C
C                          g1,g2,g3                 making q1,q2,q3  C
C      gluon polarizations a1,a2                    making e1,e2     C
C ___________________________________________________________________C
C  The component notation: p(1, 2, 3, 4) = (Beam, p_T, p_T, Energy)  C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      DIMENSION p1r(4),p1b(4), p2r(4),p2b(4), pjr(4),pjb(4)
      DIMENSION g1r(4),g1b(4), g2r(4),g2b(4), g3r(4),g3b(4)
      DIMENSION a1r(4),a1b(4), a2r(4),a2b(4)
      COMMON/CONST/CME,ALP,ALS, XC,XC2,XJ,XJ2, QC,PSI
      COMMON/POLAR/            a1(4),a2(4)
      COMMON/MOMEN/p1(4),p2(4),g1(4),g2(4),pj(4),g3(4)
      COMMON/FREST/d1(4),d2(4),q1(4),q2(4),qj(4),q3(4),e1(4),e2(4)
      SAVE
C
      pjtot = dsqrt(pj(2)**2 + pj(3)**2 + pj(1)**2)
      pjt   = dsqrt(pj(2)**2 + pj(3)**2)
C
C... Rotation around the beam axis (get pj(3) =0)
      cphi = pj(2)/pjt
      sphi = pj(3)/pjt
      CALL ROTAT(cphi,sphi,2,3, pj,pjr)
      CALL ROTAT(cphi,sphi,2,3, p1,p1r)
      CALL ROTAT(cphi,sphi,2,3, p2,p2r)
      CALL ROTAT(cphi,sphi,2,3, g1,g1r)
      CALL ROTAT(cphi,sphi,2,3, g2,g2r)
      CALL ROTAT(cphi,sphi,2,3, g3,g3r)
      CALL ROTAT(cphi,sphi,2,3, a1,a1r)
      CALL ROTAT(cphi,sphi,2,3, a2,a2r)
c     IF(dabs(rj(3)).GE.1.E-5) write(6,*) 'JREST: wrong pj(3)'
C
C... Rotation in the pj-production plane (get pj(2)= 0)
      cost = pj(1)/pjtot
      sint = pjt/pjtot
      CALL ROTAT(cost,sint,1,2, pjr,pjb)
      CALL ROTAT(cost,sint,1,2, p1r,p1b)
      CALL ROTAT(cost,sint,1,2, p2r,p2b)
      CALL ROTAT(cost,sint,1,2, g1r,g1b)
      CALL ROTAT(cost,sint,1,2, g2r,g2b)
      CALL ROTAT(cost,sint,1,2, g3r,g3b)
      CALL ROTAT(cost,sint,1,2, a1r,a1b)
      CALL ROTAT(cost,sint,1,2, a2r,a2b)
c     IF(dabs(bj(2)).GE.1.E-5) write(6,*) 'JREST: wrong pj(2)'
C
C... Boost from the Lab system to pj rest (get pj(1)= 0)
      Gam   = pj(4)/XJ
      Gambet= pjtot/XJ
      CALL BOOST(Gam,Gambet,1,4, pjb,qj)
      CALL BOOST(Gam,Gambet,1,4, p1b,d1)
      CALL BOOST(Gam,Gambet,1,4, p2b,d2)
      CALL BOOST(Gam,Gambet,1,4, g1b,q1)
      CALL BOOST(Gam,Gambet,1,4, g2b,q2)
      CALL BOOST(Gam,Gambet,1,4, g3b,q3)
      CALL BOOST(Gam,Gambet,1,4, a1b,e1)
      CALL BOOST(Gam,Gambet,1,4, a2b,e2)
c     IF(dabs(qj(1)).GE.1.E-5) write(6,*) 'JREST: wrong pj(1)'
C
      RETURN
      END
