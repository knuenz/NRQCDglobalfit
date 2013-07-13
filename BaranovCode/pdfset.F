 
C*********************************************************************
 
C...PDFSET
C...Dummy routine, to be removed when PDFLIB is to be linked.
 
      SUBROUTINE PDFSET(PARM,VALUE)
 
C...Local arrays and character variables.
      CHARACTER*20 PARM(20)
      DOUBLE PRECISION VALUE(20)
 
C...Stop program if this routine is ever called.
      WRITE(6,*) ' CASCADE: dummy routine pdfset called ' 
      WRITE(6,*) ' CASCADE: prgram stopped ' 
      Stop
 
      RETURN
      END
