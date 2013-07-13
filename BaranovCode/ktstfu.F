*CMZ :  1.00/04 22/09/2001  20.15.50  by  Hannes Jung
*CMZ :  1.00/02 12/07/2001  21.53.13  by  Hannes Jung
*CMZ :  1.00/00 12/06/2001  12.05.04  by  Hannes Jung
*-- Author : H. Jung
      SUBROUTINE KTSTFU(KF,X,Q2,XPQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL XPQ(-6:6),X,Q2
      DIMENSION XX(0:100),Q2X(0:100),XPD(0:100,0:100,-6:6)
      SAVE XX, Q2X, XPD   
      character*8 cna1,cna2,cna3
      Integer IGLU
      Common/cagluon/Iglu
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      DATA NPR/0/
      IF(FIRST) THEN
         i=0
         xold = 0.d0
         q2old = 0.d0
         nx = 0
         nq = 0
         If(Iglu.eq.4) then
            open(30,FILE='kms.dat', FORM='formatted',STATUS= 'OLD',
     +      IOSTAT=IRR,ERR=70 )
            write(6,*) ' Kwiecinski, Martin, Stasto set'
            write(6,*) ' read unintegrated gluon from file: kms.dat'
         elseif(Iglu.eq.6) then
            open(30,FILE='kmr.dat', FORM='formatted',STATUS= 'OLD',
     +      IOSTAT=IRR,ERR=70 )
            write(6,*) ' Kimber, Martin, Ryskin set'
            write(6,*) ' read unintegrated gluon from file: kmr.dat'
         else
            write(6,*) ' unintegrated gluon file not available'
         endif
         read(30,*) cna1,cna2,cna3
c         write(6,10001)
c10001    Format('    x    | kt**2  | xg(x,kt**2) ')
   10    continue
         i=i+1
         READ(30,*,END=20 ) RX,RQ2,RKMS
c         write(6,10002) RX,RQ2,RX*RKMS
c10002    Format(1X,e8.2,1X,e8.2,1X,e8.2 )
         if(q2old.ne.rq2) then
            nq = nq +1
            q2old = rq2
            nx = 1
         else
            if(xold.ne.rx) then
               nx = nx + 1
               xold = rx
            endif
         endif
c we have only gluons
         xx(nx) = rx
         q2x(nq) = rq2
         xpd(nx,nq,0) = max(0.d0,rkms/rq2)
         GOTO 10
   20    CONTINUE
         write(6,*) ' nr of points: nx = ',nx,' nq = ',nq
         DO 30 J=1,NQ
            DO 30 I=1,NX

   30    CONTINUE
         FIRST=.FALSE.
         write(6,*) '  parton densities read from file unit 30 '
      ENDIF


      XPRT = DBLE(X)
      Q2T= DBLE(Q2)
      if(dble(q2).lt.Q2x(1).and.npr.lt.5) then
         write(6,*) 'RES :  Q2 values outside grid '
         WRITE(6,*) ' Q2_min ',Q2X(1),' Q2_max ',Q2X(NQ), ' ' //
     +   'actual Q2 ',Q2T
         write(6,*) ' take parton densities at Q2_min = ',Q2X(1)
         npr = npr + 1
      endif
      IF(Q2T.LT.Q2X(1)) Q2T=Q2X(1)
      IF(XPRT.LT.XX(1).OR.XPRT.GT.XX(NX) .OR.Q2T.LT.Q2X(1) .OR.Q2T
     +.GT.Q2X(NQ)) THEN
         WRITE(6,*) 'RES : X or Q2 values outside grid '
         WRITE(6,*) ' X_min ',XX(1),' X_max ',XX(NX), ' actual '
     +   //'X ', XPRT
         WRITE(6,*) ' Q2_min ',Q2X(1),' Q2_max ',Q2X(NQ), ' ' //
     +   'actual Q2 ',Q2T
         IF(XPRT.LT.XX(1)) XPRT=XX(1)
         IF(XPRT.GT.XX(NX)) XPRT=XX(NX)
         IF(Q2T.LT.Q2X(1)) Q2T = Q2X(1)
         IF(Q2T.GT.Q2X(NQ)) Q2T = Q2X(NQ)
      ENDIF
      IX = -1
   40 IX = IX + 1
      IF(XPRT.GT.XX(IX+1)) GOTO 40
      IQ =  0
   50 IQ = IQ + 1
      IF(Q2T.GT.Q2X(IQ+1)) GOTO 50
      DO 60 IP=-6,6
         XPQ(IP)=0.0
         XD = (XPRT - XX(IX))/(XX(IX+1)-XX(IX))
         QD = (Q2T - Q2X(IQ))/(Q2X(IQ+1) - Q2X(IQ))
         X1P=(XPD(IX+1,IQ,IP)-XPD(IX,IQ,IP))*XD +XPD(IX,
     +   IQ,IP)
         X2P=(XPD(IX+1,IQ+1,IP)-XPD(IX,IQ+1,IP))*XD + XPD(
     +   IX,IQ+1,IP)
         XPQ(IP) = SNGL((X2P-X1P)*QD + X1P)
   60 CONTINUE
      return
   70 write(6,*) ' error in opening file '
      Print *, "kt"
      stop
      END
