      SUBROUTINE DECAYJ(xb,pb, xm,pm, xn,pn, alpha,costh)
************************************************************************
*     two body decay of a particle b --> m + n                         *
*     alpha is the asymmetry parameter dGamma = 1 + alpha*cos^2(theta) *
* ____________________________________________________________________ *
* Care! The component notations are different from those in PYHTIA     *
* Here:  p(1, 2, 3, 4) = (beam direction, p_T, p_T, Energy)            *
************************************************************************
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      dimension pb(4)  ! decaying particle momentum
      dimension pm(4),pm1(4),pm2(4),pm3(4) ! produced particle momentum
      dimension pn(4),pn1(4),pn2(4),pn3(4) ! ...  particle momentum
      dimension test(4)
      dimension pj(4),pj1(4),pj2(4),pj3(4),testj(4)
      COMMON/SEED/NUM
      SAVE
      xb2=xb*xb
      xm2=xm*xm
      xn2=xn*xn
      if(alpha.lt.0.D0) wmax=1.
      if(alpha.ge.0.D0) wmax=1.+alpha
      pj3(1)=0.
      pj3(2)=0.
      pj3(3)=0.
      pj3(4)=xb
C... Decay kinematics in the b-particle rest frame
      pbtot = dsqrt(pb(2)**2 + pb(3)**2 + pb(1)**2)
      pbt   = dsqrt(pb(2)**2 + pb(3)**2)
      pm3(4)= (xb2+xm2-xn2)/2./xb
      pn3(4)= (xb2+xn2-xm2)/2./xb
      pmn = AF(xb2,xm2,xn2)/2./xb
    1 cosb=-1. + 2.*RANDOM(NUM) !now take two Euler's angles at random
      alph= 3.14159*RANDOM(NUM)*2.
      weig= 1. + alpha*cosb*cosb
      wrnd= wmax*RANDOM(NUM)
      if(wrnd.ge.weig) goto 1
      costh=cosb
      sinb= 1.-cosb*cosb
      if(sinb.lt.0.D0) goto 1
      sinb= dsqrt(sinb)
      cosa= dcos(alph)
      sina= dsin(alph)
      pm3(1)= pmn*cosb
      pm3(2)= pmn*sinb*cosa
      pm3(3)= pmn*sinb*sina
      do 2 i=1,3
    2 pn3(i)=-pm3(i)
C
C... Boost from b-particle rest system to lab system
      Gam   = pb(4)/xb
      Gambet=-pbtot/xb
      CALL BOOST(Gam,Gambet,1,4, pm3,pm2)
      CALL BOOST(Gam,Gambet,1,4, pn3,pn2)
      CALL BOOST(Gam,Gambet,1,4, pj3,pj2)
C
C... Rotation to the b-production plane
      cost = pb(1)/pbtot
      sint =-pbt/pbtot
      CALL ROTAT(cost,sint,1,2, pm2,pm1)
      CALL ROTAT(cost,sint,1,2, pn2,pn1)
      CALL ROTAT(cost,sint,1,2, pj2,pj1)
C... Rotation around the beam axis
      cphi = pb(2)/pbt
      sphi =-pb(3)/pbt
      CALL ROTAT(cphi,sphi,2,3, pm1,pm)
      CALL ROTAT(cphi,sphi,2,3, pn1,pn)
      CALL ROTAT(cphi,sphi,2,3, pj1,pj)
C
C... Testing the mass and momentum conservation
      do 4 i=1,4
      testj(i)=pb(i)-pj(i)
    4 test(i)= pb(i)-pm(i)-pn(i)
      if(dabs(test(1)).ge.0.01) write(6,*) ' bad momentum 1',test(1)
      if(dabs(test(2)).ge.0.01) write(6,*) ' bad momentum 2',test(2)
      if(dabs(test(3)).ge.0.01) write(6,*) ' bad momentum 3',test(3)
      if(dabs(test(4)).ge.0.01) write(6,*) ' bad momentum 4',test(4)
      testm= pm(4)**2 -pm(3)**2 -pm(2)**2 -pm(1)**2 -xm2
      testn= pn(4)**2 -pn(3)**2 -pn(2)**2 -pn(1)**2 -xn2
      if(dabs(testm).ge.0.01) write(6,*) ' bad mass m',testm
      if(dabs(testn).ge.0.01) write(6,*) ' bad mass n',testn
C
      if(dabs(testj(1)).ge.0.001) write(6,*) ' bad j1',testj(1)
      if(dabs(testj(2)).ge.0.001) write(6,*) ' bad j2',testj(2)
      if(dabs(testj(3)).ge.0.001) write(6,*) ' bad j3',testj(3)
      if(dabs(testj(4)).ge.0.001) write(6,*) ' bad j4',testj(4)
C
   9  RETURN
      END

      SUBROUTINE DECAY2(xb,pb, xm,pm, xn,pn)
************************************************************************
*     two body decay of a particle b --> m + n   (fully at random)     *
* ____________________________________________________________________ *
* The component notation is: p(1, 2, 3, 4) = (Beam, p_T, p_T, Energy)  *
************************************************************************
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      dimension pb(4)                      ! decaying particle momentum
      dimension pm(4),pm1(4),pm2(4),pm3(4) ! produced particle momentum
      dimension pn(4),pn1(4),pn2(4),pn3(4) ! produced particle momentum
      dimension pj(4),pj1(4),pj2(4),pj3(4),test(4),testj(4)
      COMMON/ORTS/OX(4),OY(4),OZ(4)
      COMMON/ANGL/csm,phm
      COMMON/SEED/NUM
      SAVE
      xb2=xb*xb
      xm2=xm*xm
      xn2=xn*xn
      pj3(1)=0.
      pj3(2)=0.
      pj3(3)=0.
      pj3(4)=xb
C... Decay kinematics in the b-particle rest frame
      pbtot = dsqrt(pb(2)**2 + pb(3)**2 + pb(1)**2)
      pbt   = dsqrt(pb(2)**2 + pb(3)**2)
      pm3(4)= (xb2+xm2-xn2)/2./xb
      pn3(4)= (xb2+xn2-xm2)/2./xb
      pmn = AF(xb2,xm2,xn2)/2./xb
c    Two Euler's angles at random
      cosb=-1. + 2.*RANDOM(NUM) 
      alph= 3.14159*RANDOM(NUM)*2.
      sinb= dabs(1.-cosb*cosb)
      sinb= dsqrt(sinb)
      cosa= dcos(alph)
      sina= dsin(alph)
      pm3(1)= pmn*cosb
      pm3(2)= pmn*sinb*cosa
      pm3(3)= pmn*sinb*sina
      do 2 i=1,3
    2 pn3(i)=-pm3(i)
C    Angles with respect to the predefined orts
      csm = DOT3(pm3,OZ)/pmn
      cph = DOT3(pm3,OX)
      sph = DOT3(pm3,OY)
      phm = DATAN2(sph,cph) 
C
C... Boost from b-particle rest system to lab system
      Gam   = pb(4)/xb
      Gambet=-pbtot/xb
      CALL BOOST(Gam,Gambet,1,4, pm3,pm2)
      CALL BOOST(Gam,Gambet,1,4, pn3,pn2)
      CALL BOOST(Gam,Gambet,1,4, pj3,pj2)
C
C... Rotation to the b-production plane
      cost = pb(1)/pbtot
      sint =-pbt/pbtot
      CALL ROTAT(cost,sint,1,2, pm2,pm1)
      CALL ROTAT(cost,sint,1,2, pn2,pn1)
      CALL ROTAT(cost,sint,1,2, pj2,pj1)
C... Rotation around the beam axis
      cphi = pb(2)/pbt
      sphi =-pb(3)/pbt
      CALL ROTAT(cphi,sphi,2,3, pm1,pm)
      CALL ROTAT(cphi,sphi,2,3, pn1,pn)
      CALL ROTAT(cphi,sphi,2,3, pj1,pj)
C
C... Testing the mass and momentum conservation
      GOTO 9
      do 4 i=1,4
      testj(i)=pb(i)-pj(i)
    4 test(i)= pb(i)-pm(i)-pn(i)
      if(dabs(test(1)).ge.0.01) write(6,*) ' bad momentum 1',test(1)
      if(dabs(test(2)).ge.0.01) write(6,*) ' bad momentum 2',test(2)
      if(dabs(test(3)).ge.0.01) write(6,*) ' bad momentum 3',test(3)
      if(dabs(test(4)).ge.0.01) write(6,*) ' bad momentum 4',test(4)
      testm= pm(4)**2 -pm(3)**2 -pm(2)**2 -pm(1)**2 -xm2
      testn= pn(4)**2 -pn(3)**2 -pn(2)**2 -pn(1)**2 -xn2
      if(dabs(testm).ge.0.01) write(6,*) ' bad mass m',testm
      if(dabs(testn).ge.0.01) write(6,*) ' bad mass n',testn
C
      if(dabs(testj(1)).ge.0.001) write(6,*) ' bad j1',testj(1)
      if(dabs(testj(2)).ge.0.001) write(6,*) ' bad j2',testj(2)
      if(dabs(testj(3)).ge.0.001) write(6,*) ' bad j3',testj(3)
      if(dabs(testj(4)).ge.0.001) write(6,*) ' bad j4',testj(4)
    9 CONTINUE ! End of testing part
C
      RETURN
      END

      SUBROUTINE JREST2(XJ,pj, p1,p2, q1,q2)
C********************************************************************C
C    Recalculate the momenta p1,p2 into pj rest frame making q1,q2   C
C ___________________________________________________________________C
C  The component notation: p(1, 2, 3, 4) = (Beam, p_T, p_T, Energy)  C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      DIMENSION pj(4),rj(4),bj(4),qj(4)
      DIMENSION p1(4),r1(4),b1(4),q1(4)
      DIMENSION p2(4),r2(4),b2(4),q2(4)
C
      pjtot = dsqrt(pj(2)**2 + pj(3)**2 + pj(1)**2)
      pjt   = dsqrt(pj(2)**2 + pj(3)**2)
C
C... Rotation around the beam axis (get pj(3) =0)
      cphi = pj(2)/pjt
      sphi = pj(3)/pjt
      CALL ROTAT(cphi,sphi,2,3, pj,rj)
      CALL ROTAT(cphi,sphi,2,3, p1,r1)
      CALL ROTAT(cphi,sphi,2,3, p2,r2)

c     IF(dabs(rj(3)).GE.1.E-5) write(6,*) 'JREST2: wrong pj(3)'
C
C... Rotation in the pj-production plane (get pj(2)= 0)
      cost = pj(1)/pjtot
      sint = pjt/pjtot
      CALL ROTAT(cost,sint,1,2, rj,bj)
      CALL ROTAT(cost,sint,1,2, r1,b1)
      CALL ROTAT(cost,sint,1,2, r2,b2)
c     IF(dabs(bj(2)).GE.1.E-5) write(6,*) 'JREST2: wrong pj(2)'
C
C... Boost from the Lab system to pj rest (get pj(1)= 0)
      Gam   = pj(4)/XJ
      Gambet= pjtot/XJ
      CALL BOOST(Gam,Gambet,1,4, bj,qj)
      CALL BOOST(Gam,Gambet,1,4, b1,q1)
      CALL BOOST(Gam,Gambet,1,4, b2,q2)
c     IF(dabs(qj(1)).GE.1.E-5) write(6,*) 'JREST2: wrong pj(1)'
C
      RETURN
      END

      SUBROUTINE GPCMS(g1,p2, qPsi,qNew)
C********************************************************************C
C     Going to c.m.s. of the particles g1 and p2; qPsi goes to qNew  C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION (A-G,O-Z)
      DIMENSION rPsi(4),rGam(4),bGam(4)
      DIMENSION g1(4),p2(4),gp(4),qGam(4),qPro(4),qPsi(4),qNew(4)
      SAVE
C
C ... Going to the c.m.s.
      DO 1 i=1,4
    1 gp(i)=g1(i)+p2(i)
      Whadr=DOT(gp,gp)
      Whadr=dsqrt(Whadr)
      CALL JREST2(Whadr,gp, g1,p2, qGam,qPro)
      IF(dabs(qGam(1)+qPro(1)).GE.1.E-5) write(6,*) 'JPCMS: wrong cms'
      IF(dabs(qGam(2)+qPro(2)).GE.1.E-5) write(6,*) 'JPCMS: wrong cms'
      IF(dabs(qGam(3)+qPro(3)).GE.1.E-5) write(6,*) 'JPCMS: wrong cms'
C
      pjtot = dsqrt(qGam(2)**2 + qGam(3)**2 + qGam(1)**2)
      pjt   = dsqrt(qGam(2)**2 + qGam(3)**2)
C
C... Rotation to get qGam(3)=qPro(3)=0
      cphi = qGam(2)/pjt
      sphi = qGam(3)/pjt
      CALL ROTAT(cphi,sphi,2,3, qPsi,rPsi)
      CALL ROTAT(cphi,sphi,2,3, qGam,rGam)
      IF(dabs(rGam(3)).GE.1.E-5) write(6,*) 'GPCMS: wrong qGam(3)'
C
C... Rotation to get qGam(2)=qPro(2)=0
      cost = qGam(1)/pjtot
      sint = pjt/pjtot
      CALL ROTAT(cost,sint,1,2, rPsi,qNew)
      CALL ROTAT(cost,sint,1,2, rGam,bGam)
      IF(dabs(bGam(2)).GE.1.E-5) write(6,*) 'GPCMS: wrong qGam(2)'
C
      RETURN
      END

      SUBROUTINE JFRAME(p1,p2)
C********************************************************************C
C                                                                    C
C   Define the reference coordinate system according to JF choice    C
C      p1  is assumed to be the Photon momentum                      C
C      p2  is assumed to be the Proton momentum                      C
C ___________________________________________________________________C
C  The component notation: p(1, 2, 3, 4) = (Beam, p_T, p_T, Energy)  C
C********************************************************************C
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      DIMENSION p1(4),p2(4)
      COMMON/FRAME/JF
      COMMON/ORTS/OX(4),OY(4),OZ(4)
      SAVE
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
      GOTO 99
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
C
   99 RETURN
      END

      SUBROUTINE ROTAT(Ct,St,nx,ny, Pin,Pnew)
C     Rotation in the "nx*ny" plane
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
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
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
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
C     Vector product of the vectors, VR = [V1 x V2]
      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
      DIMENSION V1(4),V2(4),VR(4)
      VR(1) = V1(2)*V2(3) - V1(3)*V2(2)
      VR(2) = V1(3)*V2(1) - V1(1)*V2(3)
      VR(3) = V1(1)*V2(2) - V1(2)*V2(1)
      VR(4) = DSQRT(VR(1)**2 +VR(2)**2 +VR(3)**2)
      RETURN
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
      CALL BOOST(Gam,Gambet,1,4, qj,pjb)
      CALL BOOST(Gam,Gambet,1,4, d1,p1b)
      CALL BOOST(Gam,Gambet,1,4, d2,p2b)
      IF(dabs(pjb(4)-pj(4)).GE.1.E-5) write(6,*) 'CMSYS: wrong pj(4)'
C
C... Rotation in the pj-production plane (restore nominal pj(1))
      cost = pj(1)/pjtot
      sint =-pjt/pjtot
      CALL ROTAT(cost,sint,1,2, pjb,pjr)
      CALL ROTAT(cost,sint,1,2, p1b,p1r)
      CALL ROTAT(cost,sint,1,2, p2b,p2r)
      IF(dabs(pjr(1)-pj(1)).GE.1.E-5) write(6,*) 'CMSYS: wrong pj(1)'
C
C... Rotation around the beam axis (restore nominal pj(2),pj(3))
      cphi = pj(2)/pjt
      sphi =-pj(3)/pjt
      CALL ROTAT(cphi,sphi,2,3, pjr,qj)
      CALL ROTAT(cphi,sphi,2,3, p1r,p1)
      CALL ROTAT(cphi,sphi,2,3, p2r,p2)
      IF(dabs(qj(2)-pj(2)).GE.1.E-5) write(6,*) 'CMSYS: wrong pj(2)'
      IF(dabs(qj(3)-pj(3)).GE.1.E-5) write(6,*) 'CMSYS: wrong pj(3)'
C
      RETURN
      END
