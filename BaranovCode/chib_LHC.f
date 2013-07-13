C
C     MAIN program; needs to link files:
C                      + chi_new.f  < The essential Model !!!!!!!!!!
C                      + vegas.f    < Integrating routine
C                      + pdf_grv.f  < Structure functions
C
C###################################################################C
C                                                                   C
C    A program to simulate the production of Quarkonia mesons       C
C                                                                   C
C    in the framework of Semihard Approach and/or Parton Model      C
C   via gluon-gluon, photon-gluon or photon-photon 2-->1 fusion     C
C                                                                   C
C###################################################################C
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      PARAMETER(hp=3.1415926)
      EXTERNAL FXN
      COMMON/BVEG1/XL(10),XU(10),ACC
      COMMON/BVEGG/NDIM,NCALL,ITMX,NPRN
      COMMON/CONST/CME,PI,ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      COMMON/DECAY/xpsi,xmu2,vj(4),vg(4),v1(4),v2(4),csg,phg,csm,phm
      COMMON/TYPE/ITYPE,MODEL,IB1,IB2,NJ,NL
      COMMON/FRAME/JF1,JF2
      COMMON/SEED/NUM
      COMMON/HISW/hw(9)
      COMMON/PAWC/HBOOK(1000000)
      Common/CAGLUON/Iglu
      COMMON/CASSHWR/IORDER,ITIMSHR,ICCFM
      SAVE
      CALL METRIC
      NUM=1
C
C... Theoretical model, interaction type, and cms total energy
      MODEL = 9    ! = Semihard Theory with CAGLUON settings
      Iccfm = 1
      Iglu  = 1010 !    see CAGLUON menu in  'cauniglu.F'
C
      ITYPE = 12        ! Eta = 1S0;  Psi = 3S1;  Chi = 3PJ
                        !  Colour-singlet final states
                        ! 10 = proton*proton [glu+glu->Eta(1)]
                        ! 12 = proton*proton [glu+glu->Chi(1)]
                        ! 30 = lepton*lepton [gam+gam->Eta(1)]
                        ! 32 = lepton*lepton [gam+gam->Chi(1)]
                        !  Colour-octet final states
                        ! 40 = proton*proton [glu+glu->Eta(8)]
                        ! 41 = proton*proton [glu+glu->Psi(8)]
                        ! 42 = proton*proton [glu+glu->Chi(8)] 
                        ! 50 = lepton*proton [gam+glu->Eta(8)]
                        ! 52 = lepton*proton [gam+glu->Chi(8)] 
C... Define reference frames:
      JF1 = 4           ! for decays chiJ --> Jpsi + gamma
      JF2 = 4           ! for decays Jpsi --> muon + muon
C                       ! 1 = Recoil=Helicity frame
C                       ! 2 = Gottfried-Jackson
C                       ! 3 = Target frame
C                       ! 4 = Collins-Soper
C                       ! 5 = Perpendicular
      CALL BEAMS
c     CME   =   7.      ! sqrt(s) [GeV], HERMES
c     CME   =  14.      ! sqrt(s) [GeV], COMPASS
c     CME   =  19.      ! sqrt(s) [GeV], SMC
c     CME   = 190.      ! sqrt(s) [GeV], HERA
c     CME   = 300.      ! sqrt(s) [GeV], HERA main ep
c     CME   = 1960.     ! sqrt(s) [GeV], TEVatron, pp
      CME   = 7000.     ! sqrt(s) [GeV], CERN LHC, pp
c     CME   =14000.     ! sqrt(s) [GeV], CERN LHC, pp
C... Parameters and constants
      pi  = 3.1415926   ! pi constant unchangeable
      ALP = 1./137.     ! electromagnetic constant
C...                                       Charmonium
c      XJ = 3.1         ! J/psi mass [GeV]
c      XJ = 3.7         ! psi'  mass [GeV]
c     PSI = 0.8/(4.*pi)/XJ ! 1S state [GeV^3 /XJ]
c     PSI = 0.5/(4.*pi)/XJ ! 2S state [GeV^3 /XJ]
c      XJ = 3.53        ! chi_c(1P) mass [GeV]
c      XJ = 3.7         ! chi_c(2P) mass [GeV]
c     PSI = .075/(4.*pi)/XJ ! 1P state [GeV^5 /XJ]
c     PSI = .102/(4.*pi)/XJ ! 2P state [GeV^5 /XJ]
C                                          Bottomonium
c      XJ = 9.46         ! 1S Upsilon mass [GeV]
c      XJ =10.02         ! 2S Upsilon mass [GeV]
c      XJ =10.36         ! 3S Upsilon mass [GeV]
c     PSI =6.48/(4.*pi)/XJ ! 1S state
c     PSI =3.23/(4.*pi)/XJ ! 2S state !? 2.96/...
c     PSI =2.47/(4.*pi)/XJ ! 3S state !? 2.14/...
        XJ = 9.90        ! chi_b(1P) mass [GeV]
c       XJ =10.26        ! chi_b(2P) mass [GeV]
c       XJ =10.52        ! chi_b(3P) mass [GeV]
       PSI=1.42/(4.*pi)/XJ  ! 1P state
c      PSI=1.65/(4.*pi)/XJ  ! 2P state
c      PSI=1.79/(4.*pi)/XJ  ! 3P state
C... Identities
      XC  = XJ/2.       ! Quark mass
      QC2 =(1./3.)**2   ! Quark charge, squared
      XJ2 = XJ*XJ
      XC2 = XC*XC
C
      xpsi= 9.46  ! Jpsi mass for decays chiJ --> Jpsi + gamma
      xmu = 0.105 ! muon mass for decays Jpsi --> muon + muon
      xmu2=xmu*xmu
C
C...Open output file.
      OPEN(UNIT=2,file='chib1p_cs.outpt')
C...Book histograms.
      CALL HLIMIT(1000000)
      CALL HROPEN(1,'BB','chib1p_cs.hbook','N',1024,ISTAT)
C
      hw(1)=50./50.
      CALL HBOOK1(900,'Pt(Ups from chi0)',50, 0.,50.,0.)
      CALL HBOOK1(901,'Pt(Ups from chi1)',50, 0.,50.,0.)
      CALL HBOOK1(902,'Pt(Ups from chi2)',50, 0.,50.,0.)
      hw(2)=45./4.5
      CALL HBOOK1(920,'Y(Ups from chi0)',45, 0.,4.5,0.)
      CALL HBOOK1(921,'Y(Ups from chi1)',45, 0.,4.5,0.)
      CALL HBOOK1(922,'Y(Ups from chi2)',45, 0.,4.5,0.)
C
C... Integration limits for VEGAS
      NDIM=3
      XL(1) =alog( 0.1)! min transverse momentum of 1st parton
      XU(1) =alog(60.0)! max...
      XL(2) =alog( 0.1)! min transverse momentum of 2nd parton
      XU(2) =alog(60.0)! max...
      XL(3) = 0.D0     ! min rapidity of Quarkonium
      XU(3) = 4.5      ! max...
C... Monte-Carlo integration by VEGAS
      NCALL = 2000 ! number of points per iteration
      ITMX  = 6      ! number of iterations
      NPRN  = 1      ! print/noprint statistics
      CALL VEGAS(FXN,AVGI,SD,CHI2A)
C
      NCALL = 200000 ! number of points per iteration
      ITMX  = 12     ! number of iterations
      NPRN  = 1      ! print statistics
c     CALL VEGAS1(FXN,AVGI,SD,CHI2A)
C
C... Write output: cross sections and histograms
      CALL HISTDO
      CALL HROUT(0,ICYCLE,' ')
      CALL HREND('BB')
      WRITE(6,*) 'Last random number=',NUM
      CLOSE(UNIT=2)
      STOP
      END

      SUBROUTINE WRIOUT(hwgt)
C*******************************************************************C
C      Filling output histograms                                    C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      COMMON/CONST/CME,PI,ALP,ALS, XC,XC2,XJ,XJ2, QC,PSI
      COMMON/FREST/d1(4),d2(4),q1(4),q2(4),qj(4),e1(4),e2(4)
      COMMON/MOMEN/p1(4),p2(4),g1(4),g2(4),pj(4),a1(4),a2(4)
      COMMON/KINEM/x1,x2,g1T,g2T, pjT,yj,Zj
      COMMON/OUTPX/XChi0,XChi1(-1:1),XChi2(-2:2)
      COMMON/FLABS/qpsi(4),qgam(4),qmu1(4),qmu2(4)
      COMMON/DECAY/xpsi,xmu2,vj(4),vg(4),v1(4),v2(4),csg,phg,csm,phm
      COMMON/TYPE/ITYPE,MODEL,IB1,IB2,NJ,NL
      COMMON/TOTAL/CXTOT0,CXTOT1,CXTOT2
      COMMON/HISW/hw(9)
      SAVE
c     DATA Br0,Br1,Br2 /0.013,  0.36,  0.20 / !Chi_c to Jpsi decays
c     DATA Br0,Br1,Br2 /0.06,   0.35,  0.22 / !Chi_b to Upsl decays
c     DATA Br0,Br1,Br2 /1.D0,   1.D0,  1.D0 / !Color octet model
c     DATA Br/0.0588/  !branching J/psi --> mu+mu-
c     DATA Br/0.0248/  !branching Upsilon --> mu+mu-
      xpsi2=xpsi**2
      tiny=1.D-5      
C
C... Kinematics of the event.....
C
       pT = DSQRT(qpsi(2)**2 + qpsi(3)**2)
       y  = 0.5*DLOG((qpsi(4)+qpsi(1))/(qpsi(4)-qpsi(1)))
C
C... Plotting the distributions
      IF(NJ.EQ.1 .AND. NL.EQ.1) THEN !3PJ decay option
      CALL GAUGEX
      CALL VGAMMA(WGT1,WGT2) !In unpolarized case use CXTOTi for WGTi
      CALL HF1(900,sngl(pT),sngl(CXTOT0)*hwgt*hw(1))
      CALL HF1(901,sngl(pT),sngl(WGT1)*hwgt*hw(1))
      CALL HF1(902,sngl(pT),sngl(WGT2)*hwgt*hw(1))
      CALL HF1(920,sngl(y),sngl(CXTOT0)*hwgt*hw(2))
      CALL HF1(921,sngl(y),sngl(WGT1)*hwgt*hw(2))
      CALL HF1(922,sngl(y),sngl(WGT2)*hwgt*hw(2))
      wdec0 = hwgt*CXTOT0
      wdec1 = hwgt*WGT1
      wdec2 = hwgt*WGT2
      write(2,200) pT,y,csg,phg,wdec0,wdec1,wdec2
  200 FORMAT(2E12.3,2F7.3,3E12.3)
      ENDIF
C
      IF(NJ.EQ.1 .AND. NL.EQ.0) THEN !3S1 decay option
      CALL JMUON(WGTJ)
      CALL HF1(901,sngl(pjT),sngl(WGTJ)*hwgt*hw(1))
      CALL HF1(902,sngl(pjT),sngl(CXTOT1)*hwgt*hw(1))
      CALL HF1(921,sngl(yj),sngl(WGTJ)*hwgt*hw(2))
      wdec1 = hwgt*WGTJ
      write(2,201) pjT,yj,csm,phm,wdec1
  201 FORMAT(2E12.3,2F7.3,E12.3)
      ENDIF
C
      IF(NJ.EQ.0 .AND. NL.EQ.0) THEN !1S0 decay option
      CALL HF1(900,sngl(pjT),sngl(CXTOT0)*hwgt*hw(1))
      CALL HF1(920,sngl(yj),sngl(CXTOT0)*hwgt*hw(2))
      wdec0 = hwgt*CXTOT0
      write(2,202) pjT,yj,wdec0
  202 FORMAT(2E12.3,2F7.3,E12.3)
      ENDIF
C
      RETURN
      END

      SUBROUTINE FACTOR(CXFACT)
C*******************************************************************C
C      Parton distrubutions, Couplings, Flux factor                 C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      Parameter(GeVmb=0.4D0, GeVmub=4.D2, GeVnb=4.D5, GeVpb=4.D8)
      DIMENSION XPQ1(-6:6),XPQ2(-6:6)
      COMMON/BVEG1/XL(10),XU(10),ACC
      COMMON/TYPE/ITYPE,MODEL,IB1,IB2,NJ,NL
      COMMON/CONST/CME,PI,ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      COMMON/MOMEN/Beam1(4),Beam2(4),g1(4),g2(4),pj(4),a1(4),a2(4)
      COMMON/KINEM/x1,x2,g1T,g2T, pjT,yj,zj
      SAVE
C
C... Phase Space factor
      STOT= CME*CME
      G1W = DOT(g1,g1)
      G2W = DOT(g2,g2)
      SXX  = XJ2
      tmj2 = XJ2+pjT*pjT
C
      FluxPS= 16.*PI*PI/SXX/SXX
C
c     IF((IB1*IB2).EQ.64) FluxPS=FluxPS*SXX/AF(SXX,G1W,G2W)
      IF((IB1*IB2).EQ.64) FluxPS=FluxPS*SXX/tmj2
C
C... Couplings and other constants; averaging over spins and colors
      qR2= tmj2    !Renormalization scale in the strong coupling
      XNF= 3.      ! Number of flavours
      ALS= ALPHAS(qR2,XNF)
c     ALS= 0.2
      IF((IB1*IB2).EQ.64) Couple = PSI      *ALS**2 /64.
      IF((IB1*IB2).EQ.8)  Couple = PSI *ALP *ALS    /8.
      IF((IB1*IB2).EQ.1)  Couple = PSI *ALP**2
C
C... Gluon and Photon Distribution Functions
      qF2= tmj2    !Factorization scale in parton distributions
      qF = DSQRT(qF2)
      LFs= 0       !Factorization scheme: 0,1,2=LO,NLO(MS),NLO(DIS)
      w1T = g1T*g1T
      w2T = g2T*g2T
      IF(IB1.EQ.8) THEN
                   CALL cauniglu(2212,x1,w1T,qF,XPQ1)
                   PDF1 = XPQ1(0)
                   ELSE
                   PDF1 = ALP/PI/w1T *(1.+(1.-x1)**2)/2.
      ENDIF
      IF(IB2.EQ.8) THEN
                   CALL cauniglu(2212,x2,w2T,qF,XPQ2)
                   PDF2 = XPQ2(0)
                   ELSE
                   PDF2 = ALP/PI/w2T *(1.+(1.-x2)**2)/2.
      ENDIF
C
      IF(PDF1.LT.0.D0) PDF1=0.D0
      IF(PDF2.LT.0.D0) PDF2=0.D0
C
C... Other corrections for the matrix elements and conversion to nb
      CXFACT = FluxPS*Couple *PDF1*PDF2 *GeVnb
      IF(NJ.EQ.0 .AND. NL.EQ.0) CXFACT = CXFACT*16.
      IF(NJ.EQ.1 .AND. NL.EQ.0) CXFACT = CXFACT*16.*XC2
      IF(NJ.EQ.1 .AND. NL.EQ.1) CXFACT = CXFACT*48./XJ2
C
      RETURN
      END
