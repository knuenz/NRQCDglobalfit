      SUBROUTINE CCFMSTFG(KF,X,Q2,P,XPQ)
      Implicit None
      INTEGER IORDER,ITIMSHR,ICCFM
      COMMON /CASSHWR/IORDER,ITIMSHR,ICCFM
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
	Integer  ikincut,Ipgg,ns_sel
	Real QG0
      COMMON /GLUDAT/QG0,ikincut,Ipgg,ns_sel
      Integer IGLU
      Common/CAGLUON/Iglu
      Integer i,j,k,in,irr,igrid,ip,n2test
      Integer ipoint,kf,it
	Integer Iset
      Real scal,rx,rq2,rp,rkms,glu,xpqmax,xm
      Real xl,pl
	character *72 TXT
      LOGICAL FIRST,afirst
      DATA FIRST/.TRUE./
      Data Afirst/.TRUE./
      SAVE   
      
	if(kf.eq.22.or.iabs(kf).eq.11) then
	else
         write(6,*) ' ccfmstfg: no valid kf = ',kf
         write(6,*) ' ccfmstfg: only applicable for gamma' 	   
	   stop
	endif
      IF(FIRST) THEN
         i=0
         scal = 1.0
c         scal = 0.318
	   Iset = Iglu
         if(iset.eq.1) then
            IF(ICCFM.EQ.0) THEN
               open(30,FILE='dglap-gam.dat', FORM='formatted',STATUS=
     +         'OLD', IOSTAT=IRR,ERR=80 )
               write(6,*) ' read unintegrated gluon from file: '
     +         //'dglap-gam.dat'
            ELSEIF(ICCFM.EQ.1) THEN
               open(30,FILE='ccfm-gam.dat', FORM='formatted',
     +         STATUS= 'OLD',IOSTAT=IRR,ERR=80 )
               write(6,*) ' read unintegrated gluon from file: '
     +         //'ccfm-gam.dat'
            ENDIF
         else
            write(6,*) ' No valid IGLU selected for ccfmstfg: IGLU = ',
     +      Iset
            STOP
         endif
c         read(30,10000) Qg0,ikincut
c10000    format(' Qg0 = ',f12.8,' ikincut= ',I6)
c         read(30,10100)
c10100    format('xg,  kt, p  xgx')
c         write(6,*) ' soft cut Q0 ',Qg0,' scale factor = ',scal
c         write(6,*) ' kin cut ',ikincut
200      Read(30,101) TXT
  101    Format(A72)
         If(TXT(1:4).EQ.'  Qg') then 
	      read(txt,1000) Qg0,ikincut
1000        format(' Qg0 = ',f12.8,' ikincut= ',I6)
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Qg0') then 
	      read(txt,1000) Qg0,ikincut
c         WRITE(6,101) ' 1st line ',TXT
            goto 200
         Endif
         If(TXT(1:4).EQ.' Ipg') then 
	      read(txt,1001) Ipgg,ns_sel
1001        format(' Ipgg = ',I4,' ns_sel = ',I4)
c         WRITE(6,101) '2nd line ',TXT
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



      XA(1) = Q2
      XA(2) = X
      XA(3) = P
      XA(1) = log(Q2)
      XA(2) = log(X)
      XA(3) = log(P)
      if(xa(3).lt.px(1)) then
c         write(6,*) ' p out of range ',p,' min p ',exp(px(1))
         xa(3)=px(1)
      endif
      if(xa(3).gt.px(n3)) then
         write(6,*) ' p out of range ',p,' max p ',exp(px(n3))
         xa(3)=px(n3)
      endif
      if(xa(2).ge.xx(n2)) xa(2)=xx(n2)-0.0001
      if(xa(1).lt.q2x(1)) then
c            write(6,*) '  k2 out of range ',q2,' min ',exp(q2x(1))
         xa(1)=q2x(1)
      endif
      if(xa(1).ge.q2x(n1)) then
c            write(6,*) '  k2 out of range ',q2,' max ',exp(q2x(n1))
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
                  write(6,*) ' q2 not found ',i,q2
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
               write(6,*) ' x not found ',x,xa(2),j
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
                  write(6,*) ' p not found ',k,p
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
         write(6,*) ' new ccfmstfg x,q2,p,xpq',x,q2,p,glu
         write(6,*) ' i,j,k ',i,j,k
      endif
      return
      ENTRY ccfmmaxg(x,xpqmax)

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

      ENTRY ccfmming(x,kt,p,xpqm)
      xpqm = 100000.
      xl=log(x)
      kl=log(kt)
      pl=log(p)

      if(kl.lt.q2x(1)) then
c         write(6,*) '  k2 out of range ',kt,' min ',exp(q2x(1))
         kl=q2x(1)
      endif
      if(kl.ge.q2x(n1)) then
         write(6,*) ' ccfmmin: k2 out of range: p =',p,' x= ',x,n1
         write(6,*) ' ccfmmin: k2 out of range ',kt,' max ',exp(q2x(n1))
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
            write(6,*) ' x not found ',x,xl,j
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


   80 write(6,*) ' error in opening file '
C      Print *, "ccg"
   90 write(6,*) ' end of file at ',i,j,k,RX,RQ2,RP,RKMS
      stop
      END
