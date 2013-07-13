*CMZ :  1.00/08 07/10/2001  13.08.51  by  Hannes Jung
*CMZ :  1.00/00 12/06/2001  12.05.04  by  Hannes Jung
*-- Author :
      FUNCTION XGX(QQ2)
	Implicit none
      DOUBLE PRECISION X,XGX,QQ2
      COMMON/gluon1/X
      DOUBLE PRECISION XXX,Q2,XGLU
	EXTERNAL XGLU
      Q2=QQ2
      XXX = X
      XGX = XGLU(XXX,Q2)
      RETURN
      END
