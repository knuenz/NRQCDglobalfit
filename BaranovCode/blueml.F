*CMZ :  1.00/07 28/09/2001  18.19.09  by  Hannes Jung
*CMZ :  1.00/02 12/07/2001  21.53.12  by  Hannes Jung
*CMZ :  1.00/00 12/06/2001  12.05.04  by  Hannes Jung
*-- Author : S. Baranov
      DOUBLE PRECISION FUNCTION BLUEML(x,g2,Q2)
C*******************************************************************C
C     Gluon distribution by the method of J.Bluemlein               C
C      BLUEML = x*gluon(x,kt2,Q2)               DESY 94-121         C
C*******************************************************************C
      IMPLICIT DOUBLE PRECISION(A-G,P-Z)
      EXTERNAL BJ0,BI0
      COMMON/CONVOL/cx,cg2,cQ2
      COMMON/CONST/PI,ALP,ALS
C... Parameters and constants
*! pi constant
      PI  = 3.1415926d0
*! electromagnetic constant
      ALP = 1.d0/137.d0
*! strong coupling constant
      ALS = 0.2d0
      BLUEML=0.D0
      cx = x
      cg2= g2
      cQ2= dabs(Q2)
*!Precision of integration
      acc= 2.D-3
      IF(cg2.GE.cQ2) THEN
         BLUEML=DGAUSS(BI0, x, 1.D0, acc)
      ENDIF
      IF(cg2.LT.cQ2) THEN
         BLUEML=DGAUSS(BJ0, x, 1.D0, acc)
      ENDIF
      BLUEML=DMAX1(0.D0,BLUEML/cg2)
      RETURN
      END

      DOUBLE PRECISION FUNCTION BJ0(eta)
      IMPLICIT DOUBLE PRECISION(A-G,P-Z)
      COMMON/CONST/PI,ALP,ALS
      COMMON/CONVOL/cx,cg2,cQ2
	Double Precision XGLU
	External XGLU
      Alsb= ALS*3.d0/PI
      argj= 2.d0*dsqrt(-Alsb*dlog(eta)*dabs(dlog(cg2/cQ2)))
      argg= cx/eta
      BJ0 = Alsb/eta * DBESJ0(argj) * XGLU(argg,cQ2)
      RETURN
      END

      DOUBLE PRECISION FUNCTION BI0(eta)
      IMPLICIT DOUBLE PRECISION(A-G,P-Z)
      COMMON/CONST/PI,ALP,ALS
      COMMON/CONVOL/cx,cg2,cQ2
	Double Precision XGLU
	External XGLU
      Alsb= ALS*3.d0/PI
      argj= 2.d0*dsqrt(-Alsb*dlog(eta)*dabs(dlog(cg2/cQ2)))
      argg= cx/eta
      BI0 = Alsb/eta * DBESI0(argj) * XGLU(argg,cQ2)
      RETURN
      END
