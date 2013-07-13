 
C*********************************************************************
 
C...STRUCTM
C...Dummy routine, to be removed when PDFLIB is to be linked.
 
      SUBROUTINE STRUCTM(XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
 
C...Local variables
      DOUBLE PRECISION XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU
 
C...Stop program if this routine is ever called.
      WRITE(6,*) ' CASCADE: dummy routine structm called ' 
      WRITE(6,*) ' CASCADE: prgram stopped ' 
      Stop
 
      RETURN
      END
