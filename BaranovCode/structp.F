 
C*********************************************************************
 
C...STRUCTP
C...Dummy routine, to be removed when PDFLIB is to be linked.
 
      SUBROUTINE STRUCTP(XX,QQ2,P2,IP2,UPV,DNV,USEA,DSEA,STR,CHM,
     &BOT,TOP,GLU)
 
C...Local variables
      DOUBLE PRECISION XX,QQ2,P2,UPV,DNV,USEA,DSEA,STR,CHM,BOT,
     &TOP,GLU
 
C...Stop program if this routine is ever called.
      WRITE(6,*) ' CASCADE: dummy routine structp called ' 
      WRITE(6,*) ' CASCADE: prgram stopped ' 
      Stop
 
      RETURN
      END
