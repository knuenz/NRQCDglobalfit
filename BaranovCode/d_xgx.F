*CMZ :  1.00/08 07/10/2001  13.09.15  by  Hannes Jung
*CMZ :  1.00/00 12/06/2001  12.05.04  by  Hannes Jung
*-- Author :
      FUNCTION D_XGX(XX,QQ2)
      DOUBLE PRECISION XGX,X,XX,QQ2,DELTA,DXGX,RERR,D_XGX
      DOUBLE PRECISION DQ2
      COMMON/gluon1/X
      EXTERNAL XGX
      DELTA = 1.D-1*QQ2
      RERR = -9999.D0
      DXGX = -99999.D0
      DQ2 = QQ2
      X = XX
      CALL DFRIDR(XGX,DQ2,DELTA,DXGX,RERR)
      IF(RERR.EQ.-9999.D0) THEN
         write(6,*) ' D_XGX: error in calculating derivative of xg(x,'
     +   //'q2)'
         write(6,*) ' error = ',RERR
         DXGX = -99999.D0
      ENDIF
      IF(RERR/dxgx.GT.5.) THEN
         write(6,*) ' D_XGX: large error for derivative of xg(x,'
     +   //'q2)'
         write(6,*) ' error = ',RERR/dxgx,' larger than 5. '
         write(6,*) ' with RERR = ',RERR,' and dxgx = ',dxgx
	   write(6,*) ' results are not reliable'
	   write(6,*) ' check selected set of pdf"s'
	   write(6,*) ' are they valid for Q2 = ',DQ2,' and x = ',XX
	   write(6,*) ' stop the program; derivative set to 0.0000 '
         DXGX = 0.D0
      ENDIF
c      write(6,*) ' D_XGX: DFRIDR DXGX,DELTA,RERR ',DXGX,DELTA,RERR
      D_XGX = DXGX
      RETURN
      END
