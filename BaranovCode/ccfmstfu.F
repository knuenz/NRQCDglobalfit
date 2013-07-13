      SUBROUTINE CCFMSTFU(KF,X,Q2,P,XPQ)
      Implicit None
      Integer n1,n2,n3
      REAL XPQ(-6:6),X,Q2,P
*! new
      Parameter (n1=51,n2=51,n3=51)
      Real Q2x,xx,px,xpd
      DIMENSION Q2X(0:n1),XX(0:n2),PX(0:n3),XPD(0:n2,0:n1,0:n3)
      Real amax
      Dimension amax(0:n2)
      Real kl,kt,xpqm
      REAL XA(3),f_grid(n1,n2,n3),A(N1+N2+N3)
      INTEGER NA(3)
      DATA NA/n1,n2,n3/
      REAL FINT
      INTEGER IORDER,ITIMSHR,ICCFM
      COMMON /CASSHWR/IORDER,ITIMSHR,ICCFM
	Integer  ikincut,Ipgg,ns_sel
	Real QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Integer IGLU
      Common/CAGLUON/Iglu
	Double Precision Qscal
	Integer Iqqbar
	Common/GLUDAT2/Qscal,Iqqbar
      Integer i,j,k,in,irr,igrid,ip,n2test
      Integer ipoint,kf,it
	Integer Iset,Iccfmset
      Real scal,rx,rq2,rp,rkms,glu,xpqmax,xm
      Real xl,pl
	Integer n1min,n1max,n2min,n2max,n3min,n3max,ncall
	common/caerrstf/n1min,n1max,n2min,n2max,n3min,n3max,ncall
	character *72 TXT
	character adum
      LOGICAL FIRST,afirst
      DATA FIRST/.TRUE./
      Data Afirst/.TRUE./
      SAVE   
      
	if(iabs(kf).eq.2212) then 
	else
         write(6,*) ' ccfmstfu: no valid kf = ',kf
         write(6,*) ' ccfmstfu: only applicable for proton' 	   
	   stop
	endif
      IF(FIRST) THEN
	   n1min = 0
	   n1max = 0
	   n2min = 0
	   n2max = 0
	   n3min = 0
	   n3max = 0
	   ncall = 0
         i=0
         scal = 1.0
         Ipgg = 0
	   ns_sel = -1
	   If(IGLU.GT.1000) then
	   ISET = Int(IGLU/1000)
	   ELSE
         ISET = IGLU
	   ENDIF
         if(iset.eq.1) then
            IF(ICCFM.EQ.0) THEN
C               Print *, "try to open dglap.dat"               open(30,FILE='dglap.dat', FORM='formatted',STATUS=
     +         'OLD', IOSTAT=IRR,ERR=80 )
               write(6,*) ' read unintegrated gluon from file: '
     +         //'dglap.dat'
            ELSEIF(ICCFM.EQ.1) THEN
		   if(iglu.lt.1000) then
               open(30,FILE='ccfm.dat', FORM='formatted',STATUS= 'OLD',
     +         IOSTAT=IRR,ERR=80 )
               write(6,*) ' read unintegrated gluon from file: '
     +         //'ccfm.dat'
		   else
	         Iccfmset =  MOD(IGLU,1000)
	         write(6,*) 'Iccfmset',Iccfmset	         write(6,*) 'iccfmset',iccfmset
	            if(iccfmset.eq.1) then
                    open(30,FILE='ccfm-ipgg0ns-1-0-all-lhc.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-ipgg0ns-1-0-all-lhc.dat'
	            elseif(iccfmset.eq.2) then
                    open(30,FILE='ccfm-ipgg1ns1-0-all-lhc.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-ipgg1ns1-0-all-lhc.dat'
	            elseif(iccfmset.eq.3) then
                    open(30,FILE='ccfm-ipgg2ns2-0-all-lhc.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-ipgg2ns2-0-all-lhc.dat'
		      elseif(iccfmset.eq.10) then
C                    Print *, "try to open ccfm-final-setA0.dat"                    open(30,FILE='ccfm-final-setA0.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-final-setA0.dat'
	            elseif(iccfmset.eq.11) then
                    open(30,FILE='ccfm-final-setA0+.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-final-setA0+.dat'
	            elseif(iccfmset.eq.12) then
                    open(30,FILE='ccfm-final-setA0-.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-final-setA0-.dat'
	            elseif(iccfmset.eq.13) then
                    open(30,FILE='ccfm-final-setA1.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-final-setA1.dat'
		      elseif(iccfmset.eq.20) then
                    open(30,FILE='ccfm-final-setB0.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-final-setB0.dat'
	            elseif(iccfmset.eq.21) then
                    open(30,FILE='ccfm-final-setB0+.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-final-setB0+.dat'
 	            elseif(iccfmset.eq.22) then
                    open(30,FILE='ccfm-final-setB0-.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-final-setB0-.dat'
	            elseif(iccfmset.eq.23) then
                    open(30,FILE='ccfm-final-setA1.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-final-setB1.dat'
	            elseif(iccfmset.eq.111) then
                    open(30,FILE='ccfm-test.dat', 
     +              FORM='formatted',STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
                    write(6,*) ' read: ccfm-test.dat'
                  endif		   
               endif		   
            ENDIF
         elseif(iset.eq.6) then
            open(30,FILE='kmr.dat', FORM='formatted',STATUS= 'OLD',
     +      IOSTAT=IRR,ERR=80 )
            write(6,*) ' read unintegrated gluon from file: '
     +      //'kmr.dat'
         else
            write(6,*) ' No valid IGLU selected for CCFMSTFU: IGLU = ',
     +      Iset
            STOP
         endif
c         read(30,10000) Qg0,ikincut
c10000    format(' Qg0 = ',f12.8,' ikincut= ',I6)
c         read(30,10100)
c10100    format('xg,  kt, p  xgx') 
200      Read(30,101) TXT
  101    Format(A72)
c         WRITE(6,101) '  line ',TXT
         If(TXT(1:4).EQ.'  Qg') then 
	      read(txt,1000) adum,Qg0,adum,ikincut
C1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
1000        format(A7,f12.8,A10,I6)
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Qg0') then 
	      read(txt,1000) adum,Qg0,adum,ikincut
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Ipg') then 
	      read(txt,1001) adum,Ipgg,adum,ns_sel
c1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
1001        format(A8,I4,A10,I4)
c         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:6).EQ.' Qscal') then 
	      read(txt,1002)  adum,Qscal, adum,Iqqbar
c1002        format(' Qscal = ',f7.3,' Iqqbar = ',I4)
1002        format(A9,f7.3,A10,I4)
C         WRITE(6,101) '2nd line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' ln(') then 
c         WRITE(6,101) '2 or 3rd line',TXT
         Endif
	   
         write(6,*) ' soft cut Q0 ',Qg0,' scale factor = ',scal
         write(6,*) ' kin cut ',ikincut,' Ipgg = ',Ipgg,
     &	   ' ns_sel = ',ns_sel
         write(6,*) ' type of splitting function and non-sudakov: '
	   If(Ipgg.eq.0) then
            write(6,*) ' Pgg = asb(k)/z*Delta_ns + asb(q)/(1-z) '
		elseif(Ipgg.eq.1) then
		write(6,*) 
     &	' Pgg = asb(k)*B*(1-z)/z*Delta_ns + asb(q)*B*z/(1-z)'
c		write(6,*) ' with BB = ',BB
            elseif(Ipgg.eq.2) then
	     write(6,*) ' Pgg = asb(q)[1/z*Delta_ns + 1/(1-z)]'
            elseif(Ipgg.eq.3) then
		write(6,*) 
     &	' Pgg = asb(q)*[B*(1-z)/z*Delta_ns + B*z/(1-z)]'
		else
c     Ipgg = 0 ! Pgg=(as(k)*1/z +as(q)*1/(1-z))   
c     Ipgg = 1 ! Pgg=(as(k)*B(1-z)/z +as(q)*B*z/(1-z)), Gavins full splitting
c     Ipgg = 2 ! Pgg=(as(q)*1/z +as(q)*1/(1-z)), as(q) eveywhere
c     Ipgg = 3 ! Pgg=(as(q)*B(1-z)/z +as(q)*B*z/(1-z)), as(q) eveywhere + full split
 		write(6,*) ' non - valid splitting fct selected ',Ipgg
		stop
	   Endif
         write(6,*) ' scal factor = ',Qscal,' fact. scale = ',Iqqbar
		
         do j=1,n1
            do i=1,n2
               do k=1,n3
                  READ(30,*,Err=90 ) RX,RQ2,RP,RKMS
                  xx(i) = rx
                  q2x(j) = rq2
                  px(k) = rp
                  xpd(i,j,k) = max(0.,rkms)*scal
               enddo
            enddo
         enddo
         write(6,*) ' end of file at ',j,i,k
c we have only gluons
         DO 10 J=1,N2
            DO 10 I=1,N1
               DO 10 k=1,N3

   10    CONTINUE
         IN=0
         DO I=1,n1
            IN=IN+1
            A(IN) = q2x(I)
         ENDDO
         DO I=1,N2
            IN=IN+1
            A(IN) = xx(I)
         ENDDO
         DO I=1,N3
            IN=IN+1
            A(IN) = px(I)
         ENDDO
         DO I=1,N2
            DO k=1,N3
               DO J=1,N1
                  f_grid(j,i,k) = xpd(i,j,k)
               ENDDO
            ENDDO
         ENDDO
         FIRST=.FALSE.
         write(6,*) '  parton densities read from file unit 30 '
         Close(30)
      ENDIF

      ncall = ncall +1 


c      XA(1) = Q2
c      XA(2) = X
c      XA(3) = P
      XA(1) = log(Q2)
      XA(2) = log(X)
      XA(3) = log(P)
      if(xa(3).lt.px(1)) then
	   n3min = n3min + 1
	  
	   if(n3min.lt.10) then 
c            write(6,*) ' p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' p out of range ',xa(3),px(1)
	   elseif(n3min.eq.10) then 
c            write(6,*) ' p out of range ',p,' min p ',exp(px(1))
c            write(6,*) ' last message printed: min p'
	   endif
         xa(3)=px(1)
      endif
      if(xa(3).gt.px(n3)) then
	   n3max = n3max + 1
	   if(n3max.lt.10) then 
c            write(6,*) ' p out of range ',p,' max p ',exp(px(n3))
	   elseif(n3max.eq.10) then 
c            write(6,*) ' p out of range ',p,' max p ',exp(px(n3))
c            write(6,*) ' last message printed: max p'
	   endif
         xa(3)=px(n3)
      endif
      if(xa(2).ge.xx(n2)) xa(2)=xx(n2)-0.0001
      if(xa(2).lt.xx(1)) then
	   n2min = n2min + 1
	   if(n2min.lt.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
	   elseif(n2min.eq.10) then 
c            write(6,*) '  x out of range ',x,' min ',exp(xx(1))
c            write(6,*) ' last message printed: min x'
	   endif
	   xa(2)=xx(1)
	   endif 
      if(xa(1).lt.q2x(1)) then
	   n1min = n1min + 1
	   if(n1min.lt.10) then 
c            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
	   elseif(n1min.eq.10) then 
c            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
c            write(6,*) ' last message printed: min k2'
	   endif
         xa(1)=q2x(1)
      endif
      if(xa(1).ge.q2x(n1)) then
	   n1max = n1max + 1
         if(n1max.lt.10) then 
c            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
         elseif(n1max.eq.10) then 
c            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
c            write(6,*) ' last message printed: max k2'
	   endif
         xa(1)=q2x(n1)-0.1
      endif
c check if interpolation or grid wanted
      igrid = 0
      if(igrid.eq.1) then
         if(xa(1).lt.q2x(1)) then
            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
            xa(1)=q2x(1)
         endif
         if(xa(1).ge.q2x(n1)) then
            write(6,*) '  k2 out of range: x = ', x,n1,q2x(n1)
            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
            xa(1)=q2x(n1)-0.1
         endif
         i=0
   20    i=i+1
         if(xa(1).gt.a(na(1))) then
c       write(6,*) ' q2  ',xa(1),a(na(1))
            write(6,*) ' q2 not found ',q2,a(na(1)),q2x(n1),xa(1)
            i=na(1)
         else
            if(xa(1).ge.A(i).and.xa(1).lt.a(i+1)) Then
            else
               if(i.le.na(1)) then
                  goto 20
               else
                  write(6,*) ' at 20: q2 not found ',i,q2
               endif
            endif
         endif
         j=0
   30    j=j+1
         if(xa(2).ge.A(na(1)+j).and.xa(2).lt.a(na(1)+j+1)) Then
         else
            if(j.le.na(2)) then
               goto 30
            else
               write(6,*) ' at 30: x not found ',x,xa(2),j
            endif
         endif
         k=0
   40    k=k+1
         if(xa(3).ge.a(na(1)+na(2)+na(3))) then
            k=na(3)
c       write(6,*) ' p  ',xa(3),a(na(1)+na(2)+na(3))
         else
            if(xa(3).ge.A(na(1)+na(2)+k).and. xa(3).lt.a(na(1)+na(2)+k+
     +      1)) Then
            else
               if(k.le.na(3)) then
                  goto 40
               else
                  write(6,*) ' at 40: p not found ',k,p
               endif
            endif
         endif
         glu = f_grid(i,j,k)/q2
      else
         glu = FINT(3,XA,NA,A,f_grid)/q2
      endif
      DO  IP=-6,6
         XPQ(IP)=0.0
      ENDDO
      if(glu.gt.0.) xpq(0) = glu
      if(glu.gt.1000000) then
         write(6,*) ' new ccfmstfu x,q2,p,xpq',x,q2,p,glu
         write(6,*) ' i,j,k ',i,j,k
      endif
      return
      ENTRY ccfmmax(x,xpqmax)

      xpqmax = 0.
c first calculate maxima and store it in vector amax
      if(Afirst) then
         write(6,*) ' calculate now maxima '
         DO I=1,N2
            n2test = I
            xpqmax = 0.
            do it=n2test,n2
               DO k=1,N3
                  DO J=1,N1
                     xm=f_grid(j,it,k)/exp(q2x(j))
c              write(6,*) f_grid(j,it,k),exp(q2x(j))
                     if(xm.gt.xpqmax) xpqmax=xm
                  Enddo
               Enddo
            Enddo
            amax(i)=xpqmax
         Enddo
         Afirst=.false.
         write(6,*) ' finished with maxima '
      endif
      xl = log(x)
      DO I=1,N2
         If(xx(I).gt.xl) then
            ipoint=i
            goto 50
         Endif
      Enddo
   50 continue
      xpqmax = amax(ipoint)

      return

      ENTRY ccfmmin(x,kt,p,xpqm)
      xpqm = 100000.
      xl=log(x)
      kl=log(kt)
      pl=log(p)

      if(kl.lt.q2x(1)) then
c         write(6,*) '  k2 out of range ',kt,' min ',exp(q2x(1))
         kl=q2x(1)
      endif
      if(kl.ge.q2x(n1)) then
c         write(6,*) ' ccfmmin: k2 out of range: p =',p,' x= ',x,n1
c         write(6,*) ' ccfmmin: k2 out of range ',kt,' max ',exp(q2x(n1))
         kl=q2x(n1)-0.1
      endif
      i=0
   60 i=i+1
      if(kl.gt.a(na(1))) then
c       write(6,*) ' q2  ',xa(1),a(na(1))
         write(6,*) ' q2 not found ',q2,a(na(1)),q2x(n1),xa(1)
         i=na(1)
      else
         if(kl.ge.A(i).and.kl.lt.a(i+1)) Then
         else
            if(i.le.na(1)) then
               goto 60
            else
               write(6,*) ' ccfmmin q2 not found ',i,kt
            endif
         endif
      endif
      j=0
   70 j=j+1
      if(xl.ge.A(na(1)+j).and.xl.lt.a(na(1)+j+1)) Then
      else
         if(j.le.na(2)) then
            goto 70
         else
            write(6,*) ' ccfmmin x not found: set j=na(2) ',x,xl,j
		j = na(2)
         endif
      endif
      DO K=1,N3
         if(pl.ge.px(k)) then
            xm=max(f_grid(i,j,k)/kt,0.005)
            if(xm.lt.xpqm) xpqm=xm
         endif
      Enddo
c        write(6,*) ' ccfm x,p,xpqm ',x,p,xpqm
      return


   80 write(6,*) ' error in opening file 'C      Print *, "ccu"
   90 write(6,*) ' end of file at ',i,j,k,RX,RQ2,RP,RKMS
      stop
      END
