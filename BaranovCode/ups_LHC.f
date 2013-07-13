C
C     MAIN program; needs to link files:
C                      + js_red.f   < The essential Model
C                      + js_aux.f   < Auxiliary  routines
C                      + vegas.f    < Integrating routine
C                      + pdf_grv.f  < Structure functions
C
C###################################################################C
C                                                                   C
C      A program to simulate the production of J/psi mesons         C
C                                                                   C
C    in the framework of Semihard Approach and/or Parton Model      C
C      via gluon-gluon, photon-gluon or photon-photon fusion        C
C                                                                   C
C             !!!!!!!! An updated version !!!!!!!!!                 C
C       Includes J/psi decays and full spin density matrix          C
C###################################################################C
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      PARAMETER(Ncos=15,Nphi=18)
      PARAMETER(pi=3.1415926D0)
      PARAMETER(hp=3.1415926)
      EXTERNAL FXN
      COMMON/BVEG1/XL(10),XU(10),ACC
      COMMON/BVEGG/NDIM,NCALL,ITMX,NPRN
      COMMON/CONST/CME,ALP,ALS, XC,XC2,XJ,XJ2, QC2,PSI
      COMMON/TYPE/ITYPE,MODEL,IB1,IB2
      COMMON/MUONS/xmu,v1(4),v2(4),csm,phm
      COMMON/FRAME/JF
      COMMON/SEED/NUM
      COMMON/HISW/hw(9)
      COMMON/PAWC/HBOOK(1000000)
      Common/CAGLUON/Iglu
      COMMON/CASSHWR/IORDER,ITIMSHR,ICCFM
      SAVE
      CALL METRIC
      NUM=1   ! Random number for the first event
C      Print *, "1"
C
C... Theoretical model and choice of Gluon distributions
      MODEL = 9    ! = Semihard Theory with CAGLUON settings
      Iccfm = 1
      Iglu  = 1010 !    see CAGLUON menu in  'cauniglu.F'
C
C... Interaction type 
      ITYPE = 1         !  Colour-singlet options, Jpsi=3S1(1)
                        ! 0 = proton*proton [glu+glu->Jpsi+gam]
                        ! 1 = proton*proton [glu+glu->Jpsi+glu]
                        ! 2 = lepton*proton [gam+glu->Jpsi+glu]
                        ! 3 = lepton*lepton [gam+gam->Jpsi+gam]
      CALL BEAMS
C... C.m.s. total energy
c     CME   =   7.      ! sqrt(s) [GeV], HERMES
c     CME   =  14.      ! sqrt(s) [GeV], COMPASS
c     CME   =  19.      ! sqrt(s) [GeV], SMC
c     CME   = 200.      ! sqrt(s) [GeV], RHIC pp
c     CME   = 190.      ! sqrt(s) [GeV], HERA
c     CME   = 300.      ! sqrt(s) [GeV], HERA main ep
c     CME   = 1960.     ! sqrt(s) [GeV], TEVatron, pp
      CME   = 7000.     ! sqrt(s) [GeV], CERN LHC, pp
c     CME   =14000.     ! sqrt(s) [GeV], CERN LHC, pp
C... Parameters and constants
      ALP = 1./137.     ! electromagnetic constant
c     ALS = 0.3         ! strong coupling constant
c     XJ  = 3.1         ! J/psi mass [GeV]
c     PSI = 0.038/XJ    ! J/psi wave function
c     PSI = 0.08/XJ
c     PSI = 0.8/(4.*pi)/XJ
       XJ = 9.46         ! 1S Upsilon mass [GeV]
c      XJ =10.02         ! 2S Upsilon mass [GeV]
c      XJ =10.36         ! 3S Upsilon mass [GeV]
c     PSI=0.4/XJ         !Upsilon wave function
       PSI=6.48/(4.*pi)/XJ ! 1S state
c      PSI=3.23/(4.*pi)/XJ ! 2S state
c      PSI=2.47/(4.*pi)/XJ ! 3S state
      XC  = XJ/2.       ! heavy quark mass
      QC2 =(1./3.)**2   ! heavy quark charge, squared
      XJ2 = XJ*XJ
      XC2 = XC*XC
      xmu = 0.105       ! muon mass used for decays
C
C... Define helicity frame
      JF = 1            ! 1 = Recoil frame  
                        ! 2 = Gottfried-Jackson
                        ! 3 = Target frame
                        ! 4 = Collins-Soper                        
                        ! 5 = Perpendicular
C
C     Print *, "2"
C...Open output file.
      OPEN(UNIT=2,file='ups1s_hx.outpt')
C     Print *, "3"
C...Book histograms.
      CALL HLIMIT(1000000)
C     Print *, "4"
      CALL HROPEN(1,'BB','ups1s_hx.hbook','N',1024,ISTAT)
C     Print *, "5"
C
      hw(1)=50./50.
      CALL HBOOK1(91,'Pt(Upsilon)',50, 0.,50.,0.)
      CALL HBOOK1(92,'Pt(gluJet)',50, 0.,50.,0.)
      CALL HBOOK1(93,'Y(Upsilon)',45, 0.,4.5,0.)
C     Print *, "6"
      hw(2)=45./9.
      CALL HBOOK1(94,'Y(gluJet)',45, -3.,6.,0.)
C
C... Integration limits for VEGAS
C     Print *, "7"
                                                   NDIM=5      
                       IF(IB1.EQ.1 .AND. IB2.EQ.1) NDIM=5
      IF(MODEL.EQ.0 .AND. IB1.EQ.1 .AND. IB2.EQ.8) NDIM=4
      IF(MODEL.EQ.0 .AND. IB1.EQ.8 .AND. IB2.EQ.8) NDIM=3
C
C     Print *, "8"
      XL(1) =alog( 0.1)! min transverse momentum of Jpsi
      XU(1) =alog(50.0)! max...
      XL(2) = 0.D0   ! min rapidity of Jpsi
      XU(2) = 4.5    ! max...
      XL(3) =-3.     ! min rapidity of accompanying quantum
      XU(3) = 6.     ! max...
      IF(NDIM.GE.4) THEN
      XL(4) =alog( 0.1)! min transverse momentum of 1st parton
      XU(4) =alog(60.0)! max...
      ENDIF
      IF(NDIM.GE.5) THEN
      XL(5) =alog( 0.1)! min transverse momentum of 2nd parton
      XU(5) =alog(60.0)! max...
      ENDIF
C
C     Print *, "9"
C... Monte-Carlo integration by VEGAS
      NCALL = 6000 ! number of points per iteration 600000 default
      ITMX  = 4     ! number of iterations
      NPRN  = 1     ! print/noprint results
C     Print *, "9a"
      CALL VEGAS(FXN,AVGI,SD,CHI2A)
C
C     Print *, "10"
c     NCALL = 60000 ! number of points per iteration
c     ITMX  = 4      ! number of iterations
c     NPRN  = 1      ! print statistics
c     CALL VEGAS1(FXN,AVGI,SD,CHI2A)
C
C     Print *, "11"
C... Write output: cross sections and histograms
      CALL HISTDO
      CALL HROUT(0,ICYCLE,' ')
      CALL HREND('BB')
      WRITE(6,*) 'Last random number=',NUM
      CLOSE(UNIT=2)
C
      STOP
      END

      SUBROUTINE WRIOUT(hwgt)
C*******************************************************************C
C      Filling output histograms                                    C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      Parameter(GeVmb=0.4D0, GeVmub=4.D2, GeVnb=4.D5, GeVpb=4.D8)
      PARAMETER(PI=3.1415926D0)
      DIMENSION hbin(6)
      DIMENSION qp1(4),qp2(4)
      COMMON/CONST/CME,ALP,ALS, XC,XC2,XJ,XJ2, QC,PSI
      COMMON/MOMEN/p1(4),p2(4),g1(4),g2(4),pj(4),g3(4)
      COMMON/KINEM/x1,x2,g1T,g2T, pjT,yj,zj
      COMMON/MUONS/xmu,v1(4),v2(4),csm,phm
      COMMON/OUTPX/CXTOT,CX000,FACTOR            !use with js_red.f 
      COMMON/HISW/hw(9)
      SAVE 
      DATA hbin/10., 12., 16., 20., 30., 50./
c     DATA Br/0.0588/  !branching J/psi --> mu+mu-
c     DATA Br/0.0248/  !branching Upsilon --> mu+mu-
      tiny=1.D-5
C
C... Kinematics of the event.....
C
C>     th1=datan2(pjT+tiny,pj(1))
C>     eta1=-dlog(tan(th1/2.))
       y3 = 0.5*DLOG((g3(4)+g3(1))/(g3(4)-g3(1)))
       g3T = DSQRT(g3(2)**2 + g3(3)**2)
       dyj = dabs(yj)
C
C... Plotting the distributions
      CALL HF1(91,sngl(pjT),sngl(CXTOT)*hwgt*hw(1))
      CALL HF1(92,sngl(g3T),sngl(CXTOT)*hwgt*hw(1))
      CALL HF1(93,sngl(yj),sngl(CXTOT)*hwgt*hw(1)*10.)
      CALL HF1(94,sngl(y3),sngl(CXTOT)*hwgt*hw(2))
C
C... Polarization analysis.......
      CALL JREST2(XJ,pj, p1,p2, qp1,qp2)
      CALL JFRAME(qp1,qp2)
      CALL DECAY(WGTDec)
C
      wdec = hwgt*WGTDec*FACTOR*GeVnb
C
      write(2,200) pjT,yj,csm,phm,wdec
  200 FORMAT(5E12.3)
    9 RETURN
      END

      SUBROUTINE XDEC(DGAM)
C********************************************************************C
C   Calculate the event weight DGAM                                  C
C   from Jpsi production matrix element and Jpsi decay spin density  C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      DIMENSION SPJ(4,4)
      COMMON/CONST/CME,ALP,ALS, XC,XC2,XJ,XJ2, QC,PSI
      COMMON/MUONS/xmu,v1(4),v2(4),csm,phm
      COMMON/GMUNU/DF(4,4),DC(4)
      COMMON/LOOPJ/AMPJ(4,4)
      SAVE
      DGAM=0.D0
C... Jpsi decay spin density matrix
      DO 10 mj=1,4
      DO 10 nj=1,4
      SPJ(mj,nj)= 
     .      3.*(v1(mj)*v2(nj) +v2(mj)*v1(nj) -XJ2*DF(mj,nj)/2.)/XJ2
   10 CONTINUE
C... Convoluting with the production matrix elements
      DO 12 mj=1,4
      DO 12 nj=1,4
      DO 12 l3=1,4
c     IF(ITYPE.EQ.5) GOTO 5!Color-octet is not supported
      DGAM = DGAM -AMPJ(mj,l3)*AMPJ(nj,l3)*SPJ(mj,nj)*DC(l3)
   12 CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION DOT3(xx,yy)
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      DIMENSION xx(4),yy(4)   !P(i) = (Beam, P_t, P_t, Energy)
      DOT3 = xx(1)*yy(1) + xx(2)*yy(2) + xx(3)*yy(3)
      RETURN
      END

      DOUBLE PRECISION FUNCTION DOT2(xx,yy)
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      DIMENSION xx(4),yy(4)   !P(i) = (Beam, P_t, P_t, Energy)
      DOT2 = xx(2)*yy(2) + xx(3)*yy(3)
      RETURN
      END

      SUBROUTINE DECAY(WGTDec)
C********************************************************************C
C     Generates two-body leptonic decay  J/psi --> mu+ mu-           C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(PI=3.1415926D0)
      DIMENSION vmu1(4),vmu2(4)
      COMMON/CONST/CME,ALP,ALS, XC,XC2,XJ,XJ2, QC,PSI
      COMMON/MOMEN/p1(4),p2(4),g1(4),g2(4),pj(4),g3(4)
      COMMON/MUONS/xmu,v1(4),v2(4),csm,phm
      COMMON/ORTS/OX(4),OY(4),OZ(4)
      COMMON/SEED/NUM
      SAVE
      xmu2=xmu*xmu
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
      vmu1(i)=(OX(i)*sinb*cosa + OY(i)*sinb*sina + OZ(i)*cosb)*vtot
    5 vmu2(i)=-vmu1(i)
      csm = cosb
      phm = alph
C... Decay momenta in the beam-beam c.m.s.
      CALL CMSYS(XJ,pj, vmu1,vmu2,v1,v2)
C
C Tests for the decay J/psi --> muon+muon -----------------
      GOTO 6
      v1test=DOT(v1,v1)-xmu2
      v2test=DOT(v2,v2)-xmu2
      test1=pj(1)-v1(1)-v2(1)   
      test2=pj(2)-v1(2)-v2(2)
      test3=pj(3)-v1(3)-v2(3)
      test4=pj(4)-v1(4)-v2(4)
      IF(dabs(v1test).GE.1.D-5) write(6,*) 'Bad muon1 mass' 
      IF(dabs(v2test).GE.1.D-5) write(6,*) 'Bad muon2 mass'
      IF(dabs(test1).GE.1.D-5) write(6,*) 'Bad component 1'
      IF(dabs(test2).GE.1.D-5) write(6,*) 'Bad component 2' 
      IF(dabs(test3).GE.1.D-5) write(6,*) 'Bad component 3'
      IF(dabs(test4).GE.1.D-5) write(6,*) 'Bad component 4'
    6 CONTINUE !End of test section -----------------------
C
      CALL XDEC(DGAM) 
      WGTDec=DGAM
    7 CONTINUE
      RETURN
      END
