*CMZ :  1.00/08 14/10/2001  16.41.42  by  Hannes Jung
*CMZ :  1.00/07 29/09/2001  11.25.59  by  Hannes Jung
*-- Author :
      SUBROUTINE RYSSHA(KF,X,Q2,XPQ)
      IMPLICIT None
      REAL XPQ(-6:6),X,Q2
      Integer KF
      DOUBLE PRECISION XXX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU,
     +VAL(20)
      LOGICAL PDFFIRST
      COMMON /W50516/PDFFIRST
      CHARACTER*20 PARM(20)
      Integer MST51,NPDF
      DATA NPDF/0/,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU/9*0.D0/
      real phi,Cq,k0x,q02
      real a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,alpha,beta,eps
      Data a1/163.2/,a2/-0.137/,a3/0.888/,a4/-0.423/
      Data a5/0.000187/,a6/1.0324/,a7/1.1038/,a8/0.000141/
      Data a9/-4.344/,a10/1.3846/,a11/-1.1098/
      Data alpha/11.94/,beta/1.789/,eps/0.00001/
      Real pi
      Data pi/3.1415927/
      Logical First
      Data First/.TRUE./

      PARM(1) = 'NPTYPE'
      VAL(1) = 1

      PARM(2) = 'NGROUP'
      VAL(2) = 5
      PARM(3) = 'NSET'
      VAL(3) = 3
      IF(FIRST) THEN
         CALL PDFSET(PARM,VAL)
         FIRST=.FALSE.
      ENDIF

*! for GRV 1990
      q02=4.
      k0x=2.+0.053*0.053*exp(3.56*sqrt(max(0.,log(1./(3.*x)))))
      if(q2.lt.k0x) then
         Cq=1.
      else
         Cq=k0x/q2
      endif

      phi=a1/(a2+a3+a4)*(a2+a3*q02/q2+a4*(q02/q2)**2 +
     + alpha*x+beta/(eps+log(1./x)))*Cq*sqrt(a5/(a5+x)) *
     + (1.-a6*x**a7 *log(q2/a8))*(1.+a11*x)*(1.-x)**(a9+a10*log(q2/a8))
      if(q2.gt.q02) then
         xpq(0) = phi/4./sqrt(2.)/pi**3
c     write(6,*) ' phi = ',phi
      else
         QQ=dble(sqrt(q02))
         XXX = dble(X)
         CALL STRUCTM(XXX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
c     write(6,*) ' glu = ',glu
         xpq(0) = sngl(glu)/q02
      endif
      return
      END
