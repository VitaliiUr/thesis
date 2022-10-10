      PROGRAM MULTIPOLES_AND_SIEGERT
c
c 08.02.2015: 
c AV18 deuteron wave function at 100 NDEUT points is used
c NPP is NOT equal NDEUT
c
      implicit none
C
      DOUBLE PRECISION FMM1
      PARAMETER (FMM1=197.327)
      integer JMAX,LCHMX2
      PARAMETER (JMAX=4,LCHMX2=2+4*JMAX) ! All, t=0 and t=1 channels
      INTEGER INDEX2,NINDEX2,IAL2,IAL22,L,IS,J,IT,IZ
      LOGICAL COUPL,NNNPD2
      COMMON /PWTB2/ INDEX2(LCHMX2),NINDEX2,COUPL(LCHMX2),NNNPD2(LCHMX2)
      double precision Edeut
c
      DOUBLE PRECISION Mp,Mn,Mmu
      PARAMETER (Mp=938.272046,Mn=939.565378,Mmu=105.658366) ! MeV
      DOUBLE PRECISION Md,Pi
      DOUBLE PRECISION QBIG
      INTEGER NP,NPP
      PARAMETER (NPP=46)
      double precision PP,WPP,PBAR ! pedy_p, gdzie PP(NPP)=P0
      COMMON/PUNKTYPP/ PP(NPP),WPP(NPP), PBAR  ! pedy_p, zerowy to p0
C
      CALL PRENFAC
C
      CALL ALFA2 ! we need only states with t=1 !
C
      CALL SET_DEUTERON(Edeut) ! wave function for the Bonn B potential
C
      CALL SET_ANGLES ! asimuthal and polar angles for integrals
C
      READ(15,*) QBIG,NP,PBAR,PP,WPP
      IF (NP.NE.NPP) STOP 'NPP !'
!     QBIG=1.0 ! 1/fm
C
      CALL DEUTERON_AND_HARMONICS(QBIG)
C
CC      CALL COULOMB_MULTIPOLES(QBIG)
C
      CALL MAGNETIC_MULTIPOLES(QBIG)
C
CC      CALL ELECTRIC_MULTIPOLES(QBIG)
C
      CALL ELECTRIC_MULTIPOLES_SIEGERT(QBIG)
C
      STOP
      END

      SUBROUTINE SET_DEUTERON(Edeut)
      implicit none
C     =============
      double precision Edeut
      INTEGER NDEUT,NDEUT1
      PARAMETER (NDEUT=100)
      INTEGER I,IP,ISP
      double precision S0,S2
      double precision PDEUT,WDEUT ! pedy_p i wagi w deuteronie
      double precision SWAV,DWAV ! fala S i fala D w deuteronie
      COMMON /DEUTERON/ PDEUT(NDEUT),WDEUT(NDEUT),
     &                  SWAV(NDEUT),DWAV(NDEUT) 
C
      double precision SPL(NDEUT)
      double precision PP1,PP2,PP3
      double precision PBAR
      integer NP1,NP2,NPP
!     PARAMETER (NP1=52,NP2=20,PP1=3.0,PP2=10.0,PP3=40.0)
      PARAMETER (NPP=46)
      double precision PP,WPP ! pedy_p, zerowy to p0
      COMMON/PUNKTYPP/ PP(NPP),WPP(NPP), PBAR  ! pedy_p, zerowy to p0
C
C**** Deuteron wave function is read in
      OPEN(10,FILE='deuwf',STATUS='OLD',FORM='UNFORMATTED')
      READ (10) PDEUT,Edeut
      READ (10) SWAV
      READ (10) DWAV
      CLOSE (10)
      WRITE(6,*) ' # EDEUT= ',EDEUT
      WRITE(6,*) ' # PDEUT(I),SWAV(I),DWAV(I) '
      do i=1,NDEUT
      WRITE(6,*) PDEUT(I),SWAV(I),DWAV(I)
      enddo
C
      RETURN
      END

      SUBROUTINE D01BCF(ITYPE,A,B,C,D,N,WEIGHT,ABSCIS,IFAIL)
      implicit none
      integer N
      integer ITYPE,IFAIL  ! now DUMMY
      double precision A,B,C,D ! now DUMMY
      double precision ABSCIS(N),WEIGHT(N)
      CALL GAULEG2(ABSCIS,WEIGHT,N)
      RETURN
      END
      SUBROUTINE GAULEG2(XX,WW,N)
      implicit none
      integer N
      double precision X(500),W(500)
      double precision XX(N),WW(N)
      double precision PI,EPS,Z,P1,P2,P3,PP,Z1
      integer M,I,J
C
      PI=DACOS(-1.D0)
      EPS=1.0D-9
      M=(N+1)/2
      DO 12 I=1,M
        Z=DCOS(PI*(I-.25D0)/(N+.5D0))
    1   CONTINUE
          P1=1.D0
          P2=0.D0
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
  11      CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
          Z=Z1-P1/PP
        IF (ABS(Z-Z1).GT.EPS) GOTO 1
        X(I)=-Z
        X(N+1-I)=Z
        W(I)=2.D0/((1.D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
   12 CONTINUE
C
        DO 99 I=1,N
        XX(I)=X(I)
   99   WW(I)=W(I)
C
      RETURN
      END

      SUBROUTINE TRANSFO(N,X1,GX1,A,B,X,GX)
      implicit none
      integer N
      double precision A,B
      double precision X1(N),GX1(N),X(N),GX(N)
      double precision C,D
      integer I
C
C*****MAKES TRANSFORMATION OF GAUSS INTEGRATION POINTS AND WEIGHTS
C     FOR INTERVAL (-1.,1.) (X1,GX1) TO THESE (X,GX) FOR (A,B)
C
      C=(B-A)/2.
      D=(B+A)/2.
      DO 1 I=1,N
      GX(I)=C*GX1(I)
    1 X(I)=C*X1(I)+D
      RETURN
      END

      SUBROUTINE TRNS(NP1,NP2,NP,P1,P2,P3,XP,AP)
C     =============
C
C     TRNS BELEGT DIE FELDER XP UND AP MIT TRANSFORMIERTEN
C     GAUSS-LEGENDRE-PUNKTEN UND GEWICHTEN
C
C     NP1 PUNKTE WERDEN UEBER DIE HYPERBOLISCHE TRANSFORMATION
C
C     X --> (1.+X) / (1./P1-(1./P1-2./P2)*X)
C
C     AUF DAS INTERVALL (0.;P2) ABGEBILDET, WOBEI
C     NP1/2 PUNKTE IN (0.;P1) UND
C     NP1/2 PUNKTE IN (P1;P2) LIEGEN
C
C     NP2 PUNKTE WERDEN UEBER DIE LINEARE TRANSFORMATION
C
C     X --> (P3+P2)/2. + (P3-P2)/2.*X
C
C     AUF DAS INTERVALL (P2;P3) ABGEBILDET
C
C     NP = NP1 + NP2
C
      DIMENSION XP1(96),AP1(96),XP2(96),AP2(96)
      DIMENSION XP(NP),AP(NP)
C
      CALL D01BCF(0.,-1.,1.,0.,0.,NP1,AP1,XP1,IF)
      DO 1 I=1,NP1
      X=XP1(I)
      A=AP1(I)
      XX=1./P1-(1./P1-2./P2)*X
      XP1(I)=(1.+X) / XX
    1 AP1(I)=(2./P1-2./P2)*A / XX**2
C
      IF(NP2 .NE. 0) THEN
      CALL D01BCF(0.,-1.,1.,0.,0.,NP2,AP2,XP2,IF)
      DO 2 I=1,NP2
      X=XP2(I)
      A=AP2(I)
      DELPH=(P3-P2)/2.
      XP2(I)=(P3+P2)/2. + DELPH*X
    2 AP2(I)=DELPH*A
      ENDIF
C
      DO 3 I=1,NP1
      XP(I)=XP1(I)
    3 AP(I)=AP1(I)
C
      IF(NP2 .NE. 0) THEN
      DO 4 I=1,NP2
      XP(I+NP1)=XP2(I)
    4 AP(I+NP1)=AP2(I)
      ENDIF
C
      RETURN
      END

      SUBROUTINE CARSPH(X,R,THTA,PHI)
      implicit double precision (a-h,o-z)
      DIMENSION X(3)
C
C*****TRANSITION FROM CARTESIAN TO SPHERICAL COORDINATES
C
      PI=DACOS(-1.0D0)
      R1=X(1)**2+X(2)**2
      R=SQRT(R1+X(3)**2)
      R1=SQRT(R1)
      IF (R.LT.1.E-15) THEN
      THTA=0.
      PHI=0.
      ELSE
      IF (R1.EQ.0.0) THEN
      PHI=0.0
      THTA=0.0
      IF (X(3).LT.0.0) THTA=PI
      ELSE
      THTA=DACOS(X(3)/R)
      PHI=DATAN2(X(2),X(1))
      ENDIF
      ENDIF
      RETURN
      END

      SUBROUTINE PRENFAC
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
C
      PARAMETER (NFAC=50)
      COMMON /FACUL  / F(0:NFAC),SF(0:NFAC)
      COMMON /HILF   / D(0:NFAC),SD(0:NFAC)
      COMMON /FKLT/ FK(0:NFAC),WF(0:NFAC)
      COMMON /FAKLOG/ FL(0:NFAC)
C
C**** n! and sqrt(n!) (/FACUL/); (2n+1) and sqrt(2n+1) (/HILF/)
      F(0)=1.
      DO 1 I=1,NFAC
    1 F(I)=F(I-1)*FLOAT(I)
      DO 2 I=0,NFAC
    2 SF(I)=SQRT(F(I))
      DO 3 I=0,NFAC
      D(I)=2.*FLOAT(I)+1.
    3 SD(I)=SQRT(D(I))
C**** n! and sqrt(n!) (/FKLT/)
      FK(0)=1.
      DO 4 I=1,NFAC
    4 FK(I)=I*FK(I-1)
      DO 5 I=0,NFAC
    5 WF(I)=SQRT(FK(I))
C**** ln(n!) (/FAKLOG/)
      FL(0)=0.0
      DO 6 I=1,NFAC
    6 FL(I)=ALOG(FLOAT(I))+FL(I-1)
C
      RETURN
      END

      FUNCTION plgndr(l,m,x)
      INTEGER l,m
      DOUBLE PRECISION plgndr,x
      INTEGER i,ll
      DOUBLE PRECISION fact,pll,pmm,pmmp1,somx2
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.) stop
     *'bad arguments in plgndr'
      pmm=1.
      if(m.gt.0) then
        somx2=sqrt((1.-x)*(1.+x))
        fact=1.
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software )&5*.

      FUNCTION YTAU(L,M,THP,PHP)
      IMPLICIT NONE
      DOUBLE PRECISION THP,PHP
      DOUBLE COMPLEX YTAU
      INTEGER L,M,NFAC
      DOUBLE PRECISION Pi,E
      DOUBLE PRECISION PLGNDR
      DOUBLE PRECISION F,SF
      PARAMETER (NFAC=50)
      COMMON /FACUL/ F(0:NFAC),SF(0:NFAC)
      DOUBLE COMPLEX I
      INTEGER MM
      Pi=DACOS(-1.0d0)
      E=EXP(1.0d0)
      I=(0.0d0,1.0d0)
c
      MM=ABS(M)
      IF (MM.GT.L) THEN
      YTAU=(0.0,0.0)
      ELSE
      YTAU=SQRT((2.0*L+1.0)/(4.0*Pi))
     &     *(SF(L-MM)/SF(L+MM))
c    &     *PLGNDR(L,MM,DCOS(THP))*E**(I*MM*PHP)
     &     *PLGNDR(L,MM,DCOS(THP))*CDEXP(I*MM*PHP)
      IF (M.lt.0) YTAU=(-1)**M*CONJG(YTAU)
      ENDIF
c
      RETURN
      END

      SUBROUTINE ALFA2
      implicit integer (i-l)
      implicit double precision(a-h,o-z)
      PARAMETER (JMAX=4,LCHMX2=2+4*JMAX) ! All, t=0 and t=1 channels
      LOGICAL COUPL,NNNPD2
      CHARACTER*4 TEXT(2)
      CHARACTER*4 TEXT2(2)
      COMMON /PWTB2/ INDEX2(LCHMX2),NINDEX2,COUPL(LCHMX2),NNNPD2(LCHMX2)
C=    CHARACTER*1 LC(0:10)
C=    DATA LC /'S','P','D','F','G','H','I','J','K','L','M'/
      DATA TEXT/'NODI','DIFF'/
      DATA TEXT2/'UNCO','COUP'/
      NINDEX2=0
      write (*,100)
  100 FORMAT(//,' Channels used for two-body interaction'//
     &5X,'L',4X,'S',4X,'J',4X,'T',4X,'N',4X,'INDEX(N)'/)
      DO 10 J=0,JMAX
       DO 10 IS=0,1
C=      ISS=2*IS+1
        IF (J. EQ. 0) THEN
         JE=IS
        ELSE
         JE=J
        ENDIF
c zmiana !!!        DO 10 L=ABS(J-IS),JE
        DO 10 L=ABS(J-IS),J+IS
c      IF(.NOT.(  (L.EQ.0.AND.IS.EQ.0.AND.J.EQ.0)
c    .  .OR.       (L.EQ.0.AND.IS.EQ.1.AND.J.EQ.1) ) ) GOTO 10
         IT=((-1)**(L+IS)+1)/2
         NINDEX2=NINDEX2+1
         NNNPD2(NINDEX2)=.FALSE.
         IOUT=1
C        IF(L.EQ.0.AND.IS.EQ.0.AND.J.EQ.0.AND.IT.EQ.1) THEN
C        IOUT=2
C        NNNPD2(NINDEX2)=.TRUE.
C        ENDIF
         CALL PACK1(4,L,IS,J,IT,INDEX2(NINDEX2),IZ,IZ,IZ,IZ)
         CALL UNPACK1(4,L1,IS1,J1,IT1,INDEX2(NINDEX2),IZ,IZ,IZ,IZ)
         IF ((L .EQ. J-1).or.((J.GT.0).and.(L .EQ. J+1))) THEN
          COUPL(NINDEX2)=.TRUE.
         write(*,101) L1,IS1,J1,IT1,NINDEX2,INDEX2(NINDEX2),TEXT2(2)
  101    FORMAT(1X,5I5,I11,1X,A4)
C         write (*,102) 'COUP',ISS,LC(L),J,'-',ISS,LC(L+2),J
C 102     FORMAT ('+',39X,A4,3X,I1,A1,I1,A1,I1,A1,I1)
         ELSE
          COUPL(NINDEX2)=.FALSE.
          write(*,101) L1,IS1,J1,IT1,NINDEX2,INDEX2(NINDEX2),TEXT2(1)
C         write (*,102) 'UNCO',ISS,LC(L),J
         END IF
   10    CONTINUE
      RETURN
      END

      SUBROUTINE PACK1(N,I1,I2,I3,I4,I5,I6,I7,I8,I9)
      implicit integer (i-l)
      implicit double precision(a-h,o-z)
C     ===============
C**** EACH OF THE MAXIMAL 8 PACKED (POSITIVE) QUANTUM NUMBERS
C**** HAS TO BE GREATER EQUAL 0 AND LESS THAN 2**MOD
      INTEGER PACKY
      DIMENSION IFELD(8)
      DATA MOD/6/
      IFELD(1)=I1
      IF(N.EQ.1) GOTO 1
      IFELD(2)=I2
      IF(N.EQ.2) GOTO 1
      IFELD(3)=I3
      IF(N.EQ.3) GOTO 1
      IFELD(4)=I4
      IF(N.EQ.4) GOTO 1
      IFELD(5)=I5
      IF(N.EQ.5) GOTO 1
      IFELD(6)=I6
      IF(N.EQ.6) GOTO 1
      IFELD(7)=I7
      IF(N.EQ.7) GOTO 1
      IFELD(8)=I8
    1 CONTINUE
      PACKY=0
      DO 2 ICOUNT=1,N
    2 PACKY=PACKY+IFELD(ICOUNT)*(2**(MOD*(ICOUNT-1)))
      GOTO (10,20,30,40,50,60,70,80),N
   10 I2=PACKY
      GOTO 100
   20 I3=PACKY
      GOTO 100
   30 I4=PACKY
      GOTO 100
   40 I5=PACKY
      GOTO 100
   50 I6=PACKY
      GOTO 100
   60 I7=PACKY
      GOTO 100
   70 I8=PACKY
      GOTO 100
   80 I9=PACKY
  100 RETURN
      END
      SUBROUTINE UNPACK1(N,I1,I2,I3,I4,I5,I6,I7,I8,I9)
      implicit integer (i-l)
      implicit double precision(a-h,o-z)
C     =================
C**** EACH OF THE MAXIMAL 8 PACKED (POSITIVE) QUANTUM NUMBERS
C**** HAS TO BE GREATER EQUAL 0 AND LESS THAN 2**MOD
      DIMENSION IFELD(8)
      DATA MOD/6/
      IFELD(1)=I1
      IREST=I2
      IF(N.EQ.1) GOTO 1
      IFELD(2)=I2
      IREST=I3
      IF(N.EQ.2) GOTO 1
      IFELD(3)=I3
      IREST=I4
      IF(N.EQ.3) GOTO 1
      IFELD(4)=I4
      IREST=I5
      IF(N.EQ.4) GOTO 1
      IFELD(5)=I5
      IREST=I6
      IF(N.EQ.5) GOTO 1
      IFELD(6)=I6
      IREST=I7
      IF(N.EQ.6) GOTO 1
      IFELD(7)=I7
      IREST=I8
      IF(N.EQ.7) GOTO 1
      IFELD(8)=I8
      IREST=I9
    1 CONTINUE
      DO 2 ICOUNT=1,N
      JCOUNT=N+1-ICOUNT
      IFELD(JCOUNT)=IREST/(2**(MOD*(JCOUNT-1)))
    2 IREST=IREST-IFELD(JCOUNT)*(2**(MOD*(JCOUNT-1)))
      I1=IFELD(1)
      IF(N.EQ.1) GOTO 3
      I2=IFELD(2)
      IF(N.EQ.2) GOTO 3
      I3=IFELD(3)
      IF(N.EQ.3) GOTO 3
      I4=IFELD(4)
      IF(N.EQ.4) GOTO 3
      I5=IFELD(5)
      IF(N.EQ.5) GOTO 3
      I6=IFELD(6)
      IF(N.EQ.6) GOTO 3
      I7=IFELD(7)
      IF(N.EQ.7) GOTO 3
      I8=IFELD(8)
    3 RETURN
      END
      SUBROUTINE SPLFAK(X,N)
      IMPLICIT NONE
      DOUBLE PRECISION X
      INTEGER N
      INTEGER NSPL
      PARAMETER (NSPL=100)
C     DIE 2.ABLEITUNGEN AN DEN RAENDERN MUESSEN = 0 SEIN !
C     MIT DEN N STUETZSTELLEN X WIRD DER SPLINE AUFGEBAUT
      DOUBLE PRECISION HI,U,FAK1,FAK2,FAK3,Q,C
      DOUBLE PRECISION AX,BX,CX,AL,AM,PI,H1,H2,H3
      INTEGER I,J,L,N1,N2,N3,K,J1
      DIMENSION X(N)
      DIMENSION HI(NSPL),U(NSPL)
      COMMON /FAKTOR/ FAK1(NSPL,NSPL),FAK2(NSPL,NSPL),FAK3(NSPL,NSPL),
     C                Q(NSPL,NSPL),C(NSPL,NSPL)
      U(1)=0.
      HI(2)=X(2)-X(1)
      N1=N-1
      DO 10 I=2,N1
      AX=X(I+1)-X(I)
      HI(I+1)=AX
      BX=X(I+1)-X(I-1)
      CX=X(I)-X(I-1)
      AL=AX/BX
      AM=1.-AL
      PI=1./(2.-AM*U(I-1))
      U(I)=AL*PI
      DO 20 J=1,N
      Q(1,J)=0.
      H1=0.
      H2=0.
      H3=0.
      IF (J.EQ.I-1) H1=1./(CX*BX)
      IF (J.EQ.I) H2=1./(CX*AX)
      IF (J.EQ.I+1) H3=1./(AX*BX)
      Q(I,J)=-PI*(AM*Q(I-1,J)-H1+H2-H3)
  20  CONTINUE
  10  CONTINUE
      N2=N+1
      N3=N+2
      DO 30 K=3,N2
      J1=N3-K
      H1=1./HI(J1+1)
      DO 40 L=1,N
      C(N,L)=0.
      C(J1,L)=Q(J1,L)-C(J1+1,L)*U(J1)
      FAK1(J1,L)=-HI(J1+1)*(2.*C(J1,L)+C(J1+1,L))
      IF (L.EQ.J1) FAK1(J1,L)=FAK1(J1,L)-H1
      IF (L.EQ.J1+1) FAK1(J1,L)=FAK1(J1,L)+H1
      FAK2(J1,L)=3.*C(J1,L)
      FAK3(J1,L)=(C(J1+1,L)-C(J1,L))*H1
      FAK1(N,L)=0.
      FAK2(N,L)=0.
      FAK3(N,L)=0.
  40  CONTINUE
  30  CONTINUE
      RETURN
      END
      SUBROUTINE SPLMOD(X,N,XA,SPL)
      IMPLICIT NONE
      DOUBLE PRECISION X,XA,SPL
      INTEGER N
      INTEGER NSPL
C     =================
      PARAMETER (NSPL=100)
C     MITTELS N STUETZSTELLEN X WERDEN DIE FUNKTIONSWERTE SPL
C     AN BELIEBIGEN STELLEN XA INTERPOLIERT
C     IM HP MUSS ZUR ENDGUELTIGEN RECHNUNG NOCH
C     SUMME UEBER J SPL(J)*Y(J) GEBILDET WERDEN
      DIMENSION X(N)
      DIMENSION SPL(N)
      DOUBLE PRECISION FAK1,FAK2,FAK3,QDUM,CDUM
      INTEGER I,J,I1
      DOUBLE PRECISION DX
      COMMON /FAKTOR/ FAK1(NSPL,NSPL),FAK2(NSPL,NSPL),FAK3(NSPL,NSPL),
     C                QDUM(NSPL,NSPL),CDUM(NSPL,NSPL)
   1  I=0
   2  I=I+1
      IF (I.GT.N) GO TO 3
      IF (XA.GE.X(I)) GO TO 2
   3  I1=I-1
      IF (I1.EQ.0) I1=1
      DX=XA-X(I1)
      DO 10 J=1,N
      SPL(J)=((FAK3(I1,J)*DX+FAK2(I1,J))*DX+FAK1(I1,J))*DX
      IF (J.EQ.I1) SPL(J)=SPL(J)+1.
  10  CONTINUE
      RETURN
      END
      FUNCTION CLEBSCH(IA2,IB2,IC2,ID2,IE2,IF2)
      implicit integer (i-l)
      implicit double precision(a-h,o-z)
      double precision CLEBSCH
      integer IA2,IB2,IC2,ID2,IE2,IF2
C     ================
C
      INTEGER NFAC
      PARAMETER (NFAC=50) 
      common /FACUL/ f(0:NFAC),sf(0:NFAC)
C
        IF (ID2+IE2.NE.IF2.OR.IABS(IA2-IB2).GT.IC2.OR.IC2.GT.IA2+IB2
     & .OR.IABS(ID2).GT.IA2.OR.IABS(IE2).GT.IB2.OR.IABS(IF2).GT.IC2)
     & THEN
       CLEBSCH=0.0
       RETURN
      ENDIF
C
      J12=IA2+ID2
      J22=IA2-ID2
      J32=IB2+IE2
      J42=IB2-IE2
      J52=IC2+IF2
      J62=IC2-IF2
      J72=IA2+IB2-IC2
      J82=IB2+IC2-IA2
      J92=IC2+IA2-IB2
      L11=(IC2-IB2+ID2)/2+1
      L21=(IC2-IA2-IE2)/2+1
      L31=(IA2+IB2-IC2)/2+1
      L41=(IA2-ID2)/2+1
      L51=(IB2+IE2)/2+1
      NY2=MIN(J12,J22,J32,J42,J52,J62,J72,J82,J92)
      NY1=NY2/2+1
      K11=(IA2+IB2-IC2)/2+1
      K21=(IB2+IC2-IA2)/2+1
      K31=(IC2+IA2-IB2)/2+1
      K41=(IA2+IB2+IC2)/2+1
      DEL=SF(K11-1)*SF(K21-1)*SF(K31-1)/SF(K41)
     & *(-1.)**((IA2-IB2+IF2)/2)
     & * SF(J12/2)*SF(J22/2)*SF(J32/2)*SF(J42/2)
     & * SF(J52/2)*SF(J62/2)
      S=0.0
      IZ=0
      DO 1 IT1=1,100
      IT=IT1-1
      A=(-1.)**IT
      IF(L11+IT.LT.1) GOTO 1
      IF(L21+IT.LT.1) GOTO 1
      IF(L31-IT.LT.1) GOTO 1
      IF(L41-IT.LT.1) GOTO 1
      IF(L51-IT.LT.1) GOTO 1
      B=F(IT1-1)*F(L11+IT-1)*F(L21+IT-1)*F(L31-IT-1)*F(L41-IT-1)
     &       *F(L51-IT-1)
      S=S+A/B
      IZ=IZ+1
      IF(IZ.EQ.NY1) GOTO 2
      IF(IT1.EQ.100) STOP' 100 ERREICHT'
    1 CONTINUE
    2 CONTINUE
      CLEBSCH=S*DEL*(-1.)**((IA2-IB2+IF2)/2)*DSQRT(DFLOAT(IC2+1))
      RETURN
      END

      FUNCTION CG (I,J,K,L,M,N)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
C**************************************************************************
C
C  CLEBSCH-GORDON   < I/2 J/2 L/2 M/2 / K/2 N/2 >
C
C**************************************************************************
      PARAMETER (NFAC=50)
      INTEGER T
      COMMON /FKLT/ F(0:NFAC),WF(0:NFAC)
      IF (L+M .NE. N .OR. ABS(I-J) .GT. K .OR. K .GT. I+J .OR.
     .        ABS(L) .GT. I .OR. ABS(M) .GT. J .OR. ABS(N) .GT. K) THEN
      CG=0.
      RETURN
      ENDIF
      I1=(I+L)/2
      I2=(I-L)/2
      I3=(J+M)/2
      I4=(J-M)/2
      I5=(K-N)/2
      I6=(K+N)/2
      I7=(I+J-K)/2
      I8=(J+K-I)/2
      I9=(K+I-J)/2
      I10=(I+J+K+2)/2
      XX=WF(I7)*WF(I8)*WF(I9)*WF(I1)*WF(I2)*WF(I3)*WF(I4)*WF(I5)*
     .       WF(I6)/WF(I10)
      J1=(K-J+L)/2
      J2=(K-I-M)/2
      NT=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9)+1
      IT=0
      SUM=0.
      T=-1
   10 T=T+1
      L1=J1+T
      L2=J2+T
      L3=I7-T
      L4=I2-T
      L5=I3-T
      IF (MIN(L1,L2,L3,L4,L5) .LT. 0) GOTO 10
       SUM=SUM+RM1H(T)/(F(T)*F(L1)*F(L2)*F(L3)*F(L4)*F(L5))
       IT=IT+1
       IF (IT .LT. NT) GOTO 10
      CG=XX*SUM*SQRT(K+1.)
      END
      FUNCTION C6J (I,J,K,L,M,N)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
C**************************************************************************
C
C  6-J-SYMBOL     ( I/2  J/2  K/2 )
C                 ( L/2  M/2  N/2 )
C
C**************************************************************************
      INTEGER T,G
      PARAMETER (NFAC=50)
      DIMENSION JF(2,3)
      COMMON /FAKLOG/ FL(0:NFAC)
C
      C6J=0.0
      CALL DELMIN(I,J,K,IJK,MIJK,DIJK)
      CALL DELMIN(I,M,N,IMN,MIMN,DIMN)
      CALL DELMIN(L,J,N,LJN,MLJN,DLJN)
      CALL DELMIN(L,M,K,LMK,MLMK,DLMK)
      NT=MIN0(MIJK,MIMN,MLJN,MLMK)
      IF (NT .LT. 0) RETURN
      IF (MIN0(I,J,K,L,M,N) .EQ. 0) THEN
       JF(1,1)=I
       JF(1,2)=J
       JF(1,3)=K
       JF(2,1)=L
       JF(2,2)=M
       JF(2,3)=N
       DO 10 JJ=1,3
      DO 10 II=1,2
       IF (JF(II,JJ) .EQ. 0) THEN
        I0=II
        J0=JJ
       ENDIF
   10  CONTINUE
       G=JF(3-I0,J0)
       IF (J0 .EQ. 3) THEN
      K0=1
       ELSE
      K0=3
       ENDIF
       JK=JF(1,K0)
       JG=JF(2,K0)
      C6J=(-1.)**((JK+JG+G)/2)/SQRT((FLOAT(JK)+1.0)*(FLOAT(JG)+1.0))
      ELSE
       NTMIN=MAX0(IJK,IMN,LJN,LMK)
       W=0.5*(DIJK+DIMN+DLJN+DLMK)
       IJLM=(I+J+L+M)/2
       JKMN=(J+K+M+N)/2
       KINL=(K+I+N+L)/2
       DO 30 T=NTMIN,NTMIN+NT
       SL=W+FL(T+1)-FL(T-IJK)-FL(T-IMN)-FL(T-LJN)-FL(T-LMK)-FL(IJLM-T)
     &         -FL(JKMN-T)-FL(KINL-T)
      S=(-1.)**T * EXP(SL)
   30 C6J=C6J+S
      ENDIF
      RETURN
      END
      FUNCTION CG0003J(IA,IB,IC)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
C**** Clebsch-Gordan coefficient  <a 0 b 0 | c 0 >
      PARAMETER (NFAC=50)
      COMMON /FACUL/ F(0:NFAC),SF(0:NFAC)
      COMMON/ HILF / D(0:NFAC),SD(0:NFAC)
c     CG0003J=0.
c     IF ( MOD(IA+IB+IC,2).EQ.0  .AND.
c    &     IC.GE.IABS(IA-IB) .AND. IC.LE.(IA+IB) .AND.
c    &     IA.GE.0 .AND. IB.GE.0 .AND. IC.GE.0 ) THEN
        IP=(IA+IB+IC)/2
        CG0003J = (-1.)**IP *
     &     SF(IA+IB-IC)*SF(IB+IC-IA)*SF(IC+IA-IB)/SF(IA+IB+IC+1) *
     &     F(IP) / (F(IP-IA)*F(IP-IB)*F(IP-IC)) *
     &     (-1.)**(-IA+IB) * SD(IC)
c     ENDIF
C
      RETURN
      END
      FUNCTION CG000 (I,J,K)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
C**************************************************************************
C
C  CLEBSCH-GORDON   < I J 0 0 / K 0 >
C
C  DREIECKSUNGLEICHUNG UND PHASE VON (I+J+K) WERDEN NICHT GETESTET.
C
C**************************************************************************
      PARAMETER (NFAC=50)
      COMMON /FKLT/ F(0:NFAC),WF(0:NFAC)
      IP=(I+J+K)/2
      CG000=(-1.)**(IP+I-J)*WF(I+J-K)*WF(J+K-I)*WF(K+I-J)/WF(I+J+K+1)
     &      *F(IP)/(F(IP-I)*F(IP-J)*F(IP-K))*SQRT(2.*K+1)
      RETURN
      END
      SUBROUTINE DELMIN (J1,J2,J3,MX,MNPER,DELTA)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      PARAMETER (NFAC=50)
      COMMON /FAKLOG/ FL(0:NFAC)
C
      MX=J1+J2+J3
      IF (MOD(MX,2) .NE. 0) THEN
       MNPER=-1
      ELSE
       MX=MX/2
       I=(J1+J2-J3)/2
       J=(J2+J3-J1)/2
       K=(J3+J1-J2)/2
       MNPER=MIN0(I,J,K)
       IF (MNPER .GE. 0) DELTA=FL(I)+FL(J)+FL(K)-FL(MX+1)
      ENDIF
      END
      FUNCTION RACAH (IA,IB,IE,ID,IC,IF)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
C     ==============
CCC  ALL  L*  MUST BE TWICE AS LARGE AS TRUE ARGUMENTS
      PARAMETER (NFAC=50)
      COMMON /FAKLOG/ FL(0:NFAC)
      DIMENSION LT(6),FACLOG(NFAC+1)
      EQUIVALENCE (FACLOG(1),FL(0))
C
      K1=IA+IB-IE
      K3=IC+ID-IE
      K5=IA+IC-IF
      K7=IB+ID-IF
      K2=IE-IABS(IA-IB)
      K4=IE-IABS(IC-ID)
      K6=IF-IABS(IA-IC)
      K8=IF-IABS(IB-ID)
      K9=MIN0(K1,K2,K3,K4,K5,K6,K7,K8)
      RACAH=0.
      IF(K9)4000,20,20
   20 K2=K1-2*(K1/2)
      K4=K3-2*(K3/2)
      K6=K5-2*(K5/2)
      K8=K7-2*(K7/2)
      IF(MAX0(K2,K4,K6,K8))4000,25,4000
   25 LTMIN=MIN0(IA,IB,IC,ID,IE,IF)
      IF(LTMIN)4000,30,150
   30 LT(1)=IA
      LT(2)=IB
      LT(3)=IC
      LT(4)=ID
      LT(5)=IE
      LT(6)=IF
      LTMIN=LT(1)
      KMIN=1
      DO 40 N=2,6
        IF(LT(N)-LTMIN)35,40,40
   35   LTMIN=LT(N)
        KMIN=N
   40 CONTINUE
      S1=1.000
      F1=IE
      F2=IF
      GOTO(55,55,55,55,45,50),KMIN
   45 F1=IA
      F2=IC
      S1=(-1.)**(K5/2)
      GOTO 55
   50 F1=IA
      F2=IB
      S1=(-1.)**(K1/2)
   55 RACAH=S1/SQRT((F1+1.0)*(F2+1.0))
      GOTO 4000
  150 IABEP=(IA+IB+IE)/2+1
      ICDEP=(IC+ID+IE)/2+1
      IACFP=(IA+IC+IF)/2+1
      IBDFP=(IB+ID+IF)/2+1
      IABE=IABEP-IE
      IEAB=IABEP-IB
      IBEA=IABEP-IA
      ICDE=ICDEP-IE
      IECD=ICDEP-ID
      IDEC=ICDEP-IC
      IACF=IACFP-IF
      IFAC=IACFP-IC
      ICFA=IACFP-IA
      IBDF=IBDFP-IF
      IFBD=IBDFP-ID
      IDFB=IBDFP-IB
      NZMAX=MIN0(IABE,ICDE,IACF,IBDF)
      IABCD1=(IA+IB+IC+ID+4)/2
      IEFMAD=(IE+IF-IA-ID)/2
      IEFMBC=(IE+IF-IB-IC)/2
      NZMI1=-IEFMAD
      NZMI2=-IEFMBC
      NZMIN=MAX0(0,NZMI1,NZMI2)+1
      SQLOG=0.500*(FACLOG(IABE)+FACLOG(IEAB)+FACLOG(IBEA)+FACLOG(ICDE)
     &              +FACLOG(IECD)+FACLOG(IACF)+FACLOG(IFAC)+FACLOG(IDEC)
     &              +FACLOG(ICFA)+FACLOG(IBDF)+FACLOG(IFBD)+FACLOG(IDFB)
     & -FACLOG(IABEP+1)-FACLOG(ICDEP+1)-FACLOG(IACFP+1)-FACLOG(IBDFP+1))
      DO 200 NZ=NZMIN,NZMAX
      NZM1=NZ-1
      K1=IABCD1-NZM1
      K2=IABE-NZM1
      K3=ICDE-NZM1
      K4=IACF-NZM1
      K5=IBDF-NZM1
      K6=NZ
      K7=IEFMAD+NZ
      K8=IEFMBC+NZ
      SSLOG=SQLOG+FACLOG(K1)-FACLOG(K2)-FACLOG(K3)-FACLOG(K4)
     &-FACLOG(K5)-FACLOG(K6)-FACLOG(K7)-FACLOG(K8)
      SSTERM=(-1.)**NZM1 * EXP(SSLOG)
      RACAH=RACAH+SSTERM
  200 CONTINUE
 4000 RACAH=RACAH*(-1.)**((IA+IB+IC+ID)/2)
      RETURN
      END
      FUNCTION COEF9J (J1,J2,J5,J3,J4,J6,J7,J8,J9)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
C     ===============
CCCCCC
CCCCCC          ALL L9 MUST BE TWICE AS LARGE AS TRUE ARGUMENTS
CCCCCC
CCCC  TAMURAS NUMBERING CONVENTION IS USED HERE.
CCCC                       1  2  5
CCCC                       3  4  6
CCCC                       7  8  9
      DIMENSION LT(9)
      U9=0.000
      LT(1)=J1
      LT(2)=J2
      LT(3)=J3
      LT(4)=J4
      LT(5)=J5
      LT(6)=J6
      LT(7)=J7
      LT(8)=J8
      LT(9)=J9
      LMIN=LT(1)
      IMIN=1
      DO20I=2,9
      IF(LT(I)-LMIN)15,20,20
   15 LMIN=LT(I)
      IMIN=I
   20 CONTINUE
      KEX=0
      GOTO(110,110,110,110,150,150,170,170,190),IMIN
  110 MM=(IMIN-1)/2+1
      M1=MM+MM-1
      M2=M1+1
      M3=MM+4
      L1=LT(7)
      LT(7)=LT(M1)
      LT(M1)=L1
      L1=LT(8)
      LT(8)=LT(M2)
      LT(M2)=L1
      L1=LT(9)
      LT(9)=LT(M3)
      LT(M3)=L1
      IMIN=IMIN+(7-M1)
      GOTO175
  150 KEX=1
      M1=7
      M2=8
      M3=IMIN+IMIN-9
      M4=M3+1
      GOTO 180
  170 KEX=1
  175 M1=5
      M2=6
      M3=IMIN-6
      M4=M3+2
  180 L1=LT(M1)
      L1=LT(M1)
      LT(M1)=LT(M3)
      LT(M3)=L1
      L1=LT(M2)
      LT(M2)=LT(M4)
      LT(M4)=L1
      L1=LT(9)
      LT(9)=LT(IMIN)
      LT(IMIN)=L1
  190 IF(LT(9))200,200,300
  200 IF(LT(5)-LT(6))1000,210,1000
  210 IF(LT(7)-LT(8))1000,220,1000
  220 RT=(LT(5)+1)*(LT(7)+1)
      K=(LT(5)+LT(7)+LT(2)+LT(3))/2
      RAC=RACAH (LT(1),LT(2),LT(5),LT(4),LT(3),LT(7))
      U9=(RAC/ SQRT(RT))*(-1.)**K
      GOTO370
  300 K1=IABS(LT(2)-LT(7))
      K2=IABS(LT(3)-LT(5))
      K3=IABS(LT(4)-LT(9))
      NMIN=MAX0(K1,K2,K3)
      K1=LT(2)+LT(7)
      K2=LT(3)+LT(5)
      K3=LT(4)+LT(9)
      NMAX=MIN0(K1,K2,K3)
      IF(NMIN-NMAX)320,320,1000
  320 DO350N=NMIN,NMAX,2
      W1=N+1
      RAC=RACAH (LT(2),LT(5),LT(1),LT(3),LT(7),N)
      IF(RAC)321,350,321
  321 W1=W1*RAC
      RAC=RACAH (LT(2),LT(4),LT(8),LT(9),LT(7),N)
      IF(RAC)322,350,322
  322 W1=W1*RAC
      RAC=RACAH (LT(3),LT(4),LT(6),LT(9),LT(5),N)
      IF(RAC)323,350,323
  323 U9=U9+W1*RAC
  350 CONTINUE
      U9=(-1.)**(LT(2)+LT(3)+LT(4)+LT(5)+LT(7)+LT(9)) *U9
  370 IF(KEX)400,1000,400
  400 KP=0
      DO410I=1,9
  410 KP=KP+LT(I)
      U9=U9*(-1.)**(KP/2)
 1000 COEF9J=U9
      RETURN
      END
      FUNCTION RM1H (I)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
C**************************************************************************
C
C  RM1H = (-1) ** I
C
C**************************************************************************
      RM1H=1-2*MOD(I,2)
      END
      FUNCTION PHASE(M1)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
C     ==============
      I1=2*(M1/2)
      PHASE=-1.
      IF(I1.EQ.M1)PHASE=1.
      RETURN
      END
      FUNCTION EIMIHO(K,FAC)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
C     ===============
      EIMIHO=FAC
      IF( (K/2) * 2 .NE. K) EIMIHO=-FAC
      RETURN
      END
      FUNCTION DAFA(N1,N2)
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
C     =============
      PARAMETER (NFAC=50)
      COMMON/FACUL/F(0:NFAC),SF(0:NFAC)
      N3=N1-N2
      A=F(2*N1+2)
      B=F(2*N2+1)
      C=F(2*N3+1)
      DAFA=SQRT(A/(B*C))
      RETURN
      END

      SUBROUTINE cmplxludcmp(a,n,NPP,indx,d)
      IMPLICIT NONE
      INTEGER n,NPP,indx(n), NMAX
      DOUBLE PRECISION d,TINY
      DOUBLE COMPLEX a(NPP,NPP)
      PARAMETER (NMAX=500,TINY=1.0e-15)
C
      INTEGER i,imax,j,k
      DOUBLE COMPLEX aamax,vv(NMAX),sum,dum
C
      d=1.
      DO i=1,n
      aamax=(0.,0.)
       DO j=1,n
       if (CDABS(a(i,j)).gt.CDABS(aamax)) aamax=CDABS(a(i,j))
       ENDDO 
      if(CDABS(aamax).eq.0.) stop ' singular matrix in ludcmp ! '
      vv(i)=1./aamax
      ENDDO
C
      DO j=1,n
      DO i=1,j-1
      sum=a(i,j)
      DO k=1,i-1
      sum=sum-a(i,k)*a(k,j)
      ENDDO
      a(i,j)=sum
      ENDDO
      aamax=(0.,0.)
      DO i=j,n
      sum=a(i,j)
      DO k=1,j-1
      sum=sum-a(i,k)*a(k,j)
      ENDDO
      a(i,j)=sum
      dum=vv(i)*CDABS(sum)
      if (CDABS(dum).gt.CDABS(aamax)) then
      imax=i
      aamax=dum
      ENDIF
      ENDDO
      if (j.ne.imax) then
      DO k=1,n
      dum=a(imax,k)
      a(imax,k)=a(j,k)
      a(j,k)=dum
      ENDDO
      d=-d
      vv(imax)=vv(j)
      ENDIF
      indx(j)=imax
      if(CDABS(a(j,j)).eq.0.) a(j,j)=TINY
      if(j.ne.n) then
      dum=1./a(j,j)
      DO i=j+1,n
      a(i,j)=a(i,j)*dum
      ENDDO
      ENDIF 
      ENDDO
C
      RETURN
      END

      SUBROUTINE cmplxlubksb(a,n,NPP,indx,b)
      IMPLICIT NONE
      INTEGER n,NPP,indx(n)
      DOUBLE COMPLEX a(NPP,NPP),b(n)
      INTEGER i,ii,j,ll
      DOUBLE COMPLEX sum
C
      ii=0
      DO i=1,n
      ll=indx(i)
      sum=b(ll)
      b(ll)=b(i)
      if (ii.ne.0) then
      DO j=ii,i-1
      sum=sum-a(i,j)*b(j)
      ENDDO
      else if (CDABS(sum).ne.0.) then
      ii=i
      ENDif
      b(i)=sum
      ENDDO
      DO i=n,1,-1
      sum=b(i)
      DO j=i+1,n
      sum=sum-a(i,j)*b(j)
      ENDDO
      b(i)=sum/a(i,i)
      ENDDO
C
      RETURN
      END

      SUBROUTINE odw13a(a,n,NPP)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
C     =========================
      INTEGER n,NPP,NMAX
      DOUBLE COMPLEX a(NPP,NPP)
      PARAMETER (NMAX=500)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      DOUBLE PRECISION big
      DOUBLE COMPLEX pivinv, dum
 
c
       if(n.gt.nmax.or.NPP.gt.nmax) then
         write(6,*) ' error in odw13a!'
         write(6,*) ' n,NPP,nmax=',n,NPP,nmax
         stop
       endif
c
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (CDABS(a(j,k)).ge.big)then
                  big=CDABS(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                pause 'singular matrix in gaussj'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (CDABS(a(icol,icol)).eq.0.) pause 'singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=(1.,0.)
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=(0.,0.)
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END

      SUBROUTINE SET_ANGLES
      implicit none

      integer NPHIP,NTHETAP
      parameter (NPHIP=15,NTHETAP=15)
      double precision POINTS1(NPHIP),WEIGHTS1(NPHIP)
      double precision POINTS2(NTHETAP),WEIGHTS2(NTHETAP)
      double precision PHIM,WPHIM
      double precision THETAM,WTHETAM
      COMMON /INTANGLES/ PHIM(NPHIP),WPHIM(NPHIP),
     &                   THETAM(NTHETAP),WTHETAM(NTHETAP)

      integer IFE,IPH,ITH
      double precision Pi
C
      Pi=dacos(-1.0d0)
c
      CALL D01BCF(0,-1.0d0,1.0d0,0.0d0,0.0d0,NPHIP,WEIGHTS1,POINTS1,IFE)
      CALL TRANSFO(NPHIP,POINTS1,WEIGHTS1,0.0d0,2*Pi,PHIM,WPHIM)
      CALL D01BCF(0,-1.0d0,1.0d0,0.0d0,0.0d0,NTHETAP,WEIGHTS2,
     &                                           POINTS2,IFE)
      CALL TRANSFO(NTHETAP,POINTS2,WEIGHTS2,0.0d0,Pi,THETAM,WTHETAM)
c
      write(6,*) ' PHIM(IPH), WPHIM(IPH) '
      do IPH=1,NPHIP
      write(6,*) PHIM(IPH), WPHIM(IPH)
      enddo
      write(6,*)
c
      write(6,*) ' THETAM(ITH), WTHETAM(ITH) '
      do ITH=1,NTHETAP
      write(6,*) THETAM(ITH), WTHETAM(ITH)
      enddo
      write(6,*)
c
      RETURN 
      END

      SUBROUTINE DEUTERON_AND_HARMONICS(QBIG)
      implicit none
      double precision QBIG

      integer JMAX,LCHMX2,LMAX
      PARAMETER (JMAX=4,LCHMX2=2+4*JMAX,LMAX=JMAX+2) ! All, t=0 and t=1 channels

      integer NDEUT
      PARAMETER (NDEUT=100)
      double precision PDEUT,WDEUT ! points and weights for the deuteron
      double precision SWAV,DWAV
      COMMON /DEUTERON/ PDEUT(NDEUT),WDEUT(NDEUT),
     &                  SWAV(NDEUT),DWAV(NDEUT) 

      integer NPHIP,NTHETAP,NPHIQ,NTHETAQ
      parameter (NPHIP=15,NTHETAP=15,NPHIQ=NPHIP,NTHETAQ=NTHETAP)
      double precision PHIM,WPHIM
      double precision THETAM,WTHETAM
      COMMON /INTANGLES/ PHIM(NPHIP),WPHIM(NPHIP),
     &                   THETAM(NTHETAP),WTHETAM(NTHETAP)

      integer NPP
      PARAMETER (NPP=46)
      double precision PP,WPP,PBAR ! pedy_p, zerowy to p0
      COMMON/PUNKTYPP/ PP(NPP),WPP(NPP), PBAR  ! pedy_p, zerowy to p0

      integer ithq,iphq
      integer ithp,iphp
      integer ipp,interp
      double precision thq,phq,p,php,thp,P1FM,SPL(NDEUT)
      double precision vecq(3),vec(3),theta,phi
      double precision wavinterp
      COMMON /DEUINTRP/ wavinterp(NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ,2)

      double complex YDEU
      COMMON /YDEUT/ YDEU(NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ,0:2,-2:2)

      double precision  VECP1 ! right 
      double precision VECP1P ! left
      COMMON /VEC1/ VECP1(3,NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ),
     &             VECP1P(3,NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ)

      double precision VECPBIG(3)

      double complex YLM
      COMMON /YPQ/ YLM(NTHETAP,NPHIP,0:LMAX,-LMAX:LMAX)

      LOGICAL COUPL,NNNPD2
      integer NINDEX2,INDEX2
      COMMON /PWTB2/ INDEX2(LCHMX2),NINDEX2,COUPL(LCHMX2),NNNPD2(LCHMX2)

      double complex YTAU
      external YTAU

      integer IAL2,LP,ISP,JP,ITP,IZ,MLP,LD,MLD

      VECPBIG(1)=0.0 ! LAB FRAME
      VECPBIG(2)=0.0 ! LAB FRAME
      VECPBIG(3)=0.0 ! LAB FRAME

      wavinterp=0.0d0
      YDEU=(0.0d0,0.0d0)
      YLM=(0.0d0,0.0d0)
 
      DO LP=0,LMAX
          do mlp=-LP,LP
            do ithp=1,NTHETAP
              thp=THETAM(ithp)
             do iphp=1,NPHIP
               php=PHIM(iphp)
        YLM(ithp,iphp,LP,MLP)=CONJG(YTAU(lp,mlp,thp,php))*
     &    Sin(thp)*WTHETAM(ithp)*WPHIM(iphp)
            enddo ! iphp
           enddo ! ithp
          enddo ! mlp
      ENDDO ! LP
C

      CALL SPLFAK(PDEUT,NDEUT)
      do ithq=1,NTHETAQ
       thq=THETAM(ithq)
        do iphq=1,NPHIQ
         phq=PHIM(iphq)
         vecq(1)= 0.5*QBIG*Sin(thq)*Cos(phq)
         vecq(2)= 0.5*QBIG*Sin(thq)*Sin(phq)
         vecq(3)= 0.5*QBIG*Cos(thq)
          do ipp=1,NPP
            p=PP(ipp)
            do ithp=1,NTHETAP
             thp=THETAM(ithp)
             do iphp=1,NPHIP
              php=PHIM(iphp)
              vec(1)=p*Sin(thp)*Cos(php) - vecq(1)
              vec(2)=p*Sin(thp)*Sin(php) - vecq(2)
              vec(3)=p*Cos(thp) - vecq(3)
               VECP1(1,ipp,ithp,iphp,ithq,iphq)=vec(1)
               VECP1(2,ipp,ithp,iphp,ithq,iphq)=vec(2)
               VECP1(3,ipp,ithp,iphp,ithq,iphq)=vec(3)
              VECP1P(1,ipp,ithp,iphp,ithq,iphq)=vec(1)+2*vecq(1)
              VECP1P(2,ipp,ithp,iphp,ithq,iphq)=vec(2)+2*vecq(2)
              VECP1P(3,ipp,ithp,iphp,ithq,iphq)=vec(3)+2*vecq(3)
              CALL CARSPH(vec,P1FM,theta,phi) ! angles in rad
!             P1FM= Sqrt( 4*p**2 + QBIG**2 - 4*p*QBIG*
!    &     (Cos(thp)*Cos(thq) + Cos(php - phq)*Sin(thp)*Sin(thq)) )/2.
              IF (P1FM.LT.PDEUT(NDEUT)) THEN
                CALL SPLMOD(PDEUT,NDEUT,P1FM,SPL)
              ELSE
                DO interp=1,NDEUT
                  SPL(interp)=0.0
                ENDDO
              ENDIF
             wavinterp(ipp,ithp,iphp,ithq,iphq,1)=0.0 
             wavinterp(ipp,ithp,iphp,ithq,iphq,2)=0.0 
               do interp=1,NDEUT
                wavinterp(ipp,ithp,iphp,ithq,iphq,1)=
     &         wavinterp(ipp,ithp,iphp,ithq,iphq,1) + 
     &              swav(interp)*SPL(interp)
                wavinterp(ipp,ithp,iphp,ithq,iphq,2)=
     &             wavinterp(ipp,ithp,iphp,ithq,iphq,2) + 
     &              dwav(interp)*SPL(interp)
               enddo ! interp

            do ld=0,2,2
            do mld=-ld,ld
             YDEU(ipp,ithp,iphp,ithq,iphq,ld,mld)=
     &          YTAU(ld,mld,theta,phi)
            enddo ! ld
            enddo ! mld

           enddo ! iphp
         enddo  ! ithp
       enddo ! ipp
       enddo ! iphq
       enddo ! ithq 
c
c
      RETURN
      END

      SUBROUTINE COULOMB_MULTIPOLES(QBIG)
      implicit none
      double precision QBIG

      integer JMAX,LCHMX2,LMAX
      PARAMETER (JMAX=4,LCHMX2=2+4*JMAX,LMAX=JMAX+2) ! All, t=0 and t=1 channels

      integer NDEUT
      PARAMETER (NDEUT=100)
      double precision PDEUT,WDEUT ! points and weights for the deuteron
      double precision SWAV,DWAV
      COMMON /DEUTERON/ PDEUT(NDEUT),WDEUT(NDEUT),
     &                  SWAV(NDEUT),DWAV(NDEUT) 

      integer NPHIP,NTHETAP,NPHIQ,NTHETAQ
      parameter (NPHIP=15,NTHETAP=15,NPHIQ=NPHIP,NTHETAQ=NTHETAP)
      double precision PHIM,WPHIM
      double precision THETAM,WTHETAM
      COMMON /INTANGLES/ PHIM(NPHIP),WPHIM(NPHIP),
     &                   THETAM(NTHETAP),WTHETAM(NTHETAP)

      integer NPP
      PARAMETER (NPP=46)
      double precision PP,WPP,PBAR ! pedy_p, zerowy to p0
      COMMON/PUNKTYPP/ PP(NPP),WPP(NPP), PBAR  ! pedy_p, zerowy to p0

      integer ithq,iphq
      integer ithp,iphp
      integer ipp,interp
      double precision thq,phq,p,php,thp,P1FM,SPL(NDEUT)
      double precision vecq(3),vec(3),theta,phi
      double precision wavinterp
      COMMON /DEUINTRP/ wavinterp(NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ,2)

      double complex YDEU
      COMMON /YDEUT/ YDEU(NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ,0:2,-2:2)

      double complex YLM
      COMMON /YPQ/ YLM(NTHETAP,NPHIP,0:LMAX,-LMAX:LMAX)

      double precision  VECP1 ! right 
      double precision VECP1P ! left
      COMMON /VEC1/ VECP1(3,NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ),
     &             VECP1P(3,NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ)

      LOGICAL COUPL,NNNPD2
      integer NINDEX2,INDEX2
      COMMON /PWTB2/ INDEX2(LCHMX2),NINDEX2,COUPL(LCHMX2),NNNPD2(LCHMX2)

      double complex TJM(NPP,LCHMX2,-JMAX:JMAX,-1:1,0:LMAX,-LMAX:LMAX)
      double complex DENSITY(NPP,LCHMX2,-JMAX:JMAX,-1:1)

      double precision CL1,CL2,CL3,CL4,PI

      double precision CLEBSCH
      external CLEBSCH

      double complex RHO
      external RHO

      integer IAL2,LP,ISP,JP,ITP,IZ
      integer MLD,MLP,IALD,IS,L,J,IT
      integer m1p,m1p2,m1d,m1d2,mjp,md
      integer JBIG,MBIG

      PI=DACOS(-1.0d0)

      DO IAL2=1,NINDEX2
       CALL UNPACK1(4,LP,ISP,JP,ITP,INDEX2(IAL2),IZ,IZ,IZ,IZ)
        do mjp=-JP,JP
          do mlp=-LP,LP
            IF (abs(mjp-mlp).le.isp) THEN
             CL1=CLEBSCH(2*lp,2*isp,2*jp,2*mlp,2*(mjp-mlp),2*mjp)
              DO IALD=1,2
               L=2*IALD-2
               IS=1
               J=1
               IT=0
                do md=-1,1
                 do mld=-L,L
                  IF (abs(md-mld).le.is) THEN
                   CL2=CL1*CLEBSCH(2*L,2,2,2*mld,2*(md-mld),2*md)
                    do m1p=1,2
                     m1p2=2*m1p-3
                      IF (ABS(2*(mjp-mlp)-m1p2).le.1) THEN
                        CL3=CL2*
     &           CLEBSCH(1,1,2*ISP,m1p2,2*(mjp-mlp)-m1p2,2*(mjp-mlp))
                       do m1d=1,2
                        m1d2=2*m1d-3
                           IF ((ABS(2*(md-mld)-m1d2).le.1).and.
     &            ((2*(mjp-mlp)-m1p2).eq.(2*(md-mld)-m1d2))) THEN
                         CL4=CL3*
     &           CLEBSCH(1,1,2,m1d2,2*(md-mld)-m1d2,2*(md-mld))
      DO JBIG=0,JMAX+1
       DO MBIG=-JBIG,JBIG
        do ithq=1,NTHETAQ
         do iphq=1,NPHIQ
          do ipp=1,NPP
           do ithp=1,NTHETAP
            do iphp=1,NPHIP
               TJM(ipp,IAL2,MJP,MD,JBIG,MBIG)=
     &             TJM(ipp,IAL2,MJP,MD,JBIG,MBIG) +
     &        CL4*YLM(ithp,iphp,LP,MLP)*
     &        YDEU(ipp,ithp,iphp,ithq,iphq,L,mld)*
     &        wavinterp(ipp,ithp,iphp,ithq,iphq,IALD)*
     &        DCONJG(YLM(ithq,iphq,JBIG,MBIG))*
     &        RHO(m1p2,m1d2,ITP,Qbig,
     &            VECP1P(1,ipp,ithp,iphp,ithq,iphq),
     &             VECP1(1,ipp,ithp,iphp,ithq,iphq))
            enddo ! iphp
           enddo ! ithp
          enddo ! ipp
         enddo ! iphq
        enddo ! ithq
       enddo ! MBIG
       enddo ! JBIG
c
c
                           ENDIF
                       enddo ! m1d
                      ENDIF 
                    enddo ! m1p
                   ENDIF 
              enddo ! mld
             enddo ! md
           ENDDO ! IALD
            ENDIF 
          enddo ! mlp
        enddo ! mjp
      ENDDO ! ial2
C
C
C
      WRITE(6,*) ' Coulomb multipoles: TJM '
      DO IAL2=1,NINDEX2
       CALL UNPACK1(4,LP,ISP,JP,ITP,INDEX2(IAL2),IZ,IZ,IZ,IZ)
        do mjp=-JP,JP
        do md=-1,1
         DO JBIG=0,JMAX+1
          DO MBIG=-JBIG,JBIG
            write(6,2000) IAL2,JP,MJP,MD,JBIG,MBIG
 2000  format(' IAL2,JP,MJP,MD,JBIG,MBIG= ',6I4)
             do ipp=1,NPP
               WRITE(6,2001) PP(ipp), TJM(ipp,IAL2,MJP,MD,JBIG,MBIG)
 2001  format(F14.8,2E16.8)
             enddo ! ipp
            write(6,*) 
          ENDDO ! MBIG
         ENDDO ! JBIG
        ENDDO ! md
       ENDDO ! mjp
      ENDDO ! IAL2
C
C
C
      DENSITY=(0.0d0,0.0d0)
      DO IAL2=1,NINDEX2
       CALL UNPACK1(4,LP,ISP,JP,ITP,INDEX2(IAL2),IZ,IZ,IZ,IZ)
        do mjp=-JP,JP
        do md=-1,1
         DO JBIG=abs(JP-1),JP+1
!         DO JBIG=0,JMAX+1
!          DO MBIG=-JBIG,JBIG
             do ipp=1,NPP
              DENSITY(ipp,IAL2,MJP,MD)=DENSITY(ipp,IAL2,MJP,MD) + 
     &         DSQRT((2*JBIG+1.0d0)/(4*Pi))*TJM(ipp,IAL2,MJP,MD,JBIG,0)
             enddo ! ipp
!          ENDDO ! MBIG
         ENDDO ! JBIG
        ENDDO ! md
       ENDDO ! mjp
      ENDDO ! IAL2
C
C
            write(6,*) 
            write(6,*) 
            write(6,*)
C
         write(6,*) ' DENSITY: general mjp ' 
      DO IAL2=1,NINDEX2
       CALL UNPACK1(4,LP,ISP,JP,ITP,INDEX2(IAL2),IZ,IZ,IZ,IZ)
        do mjp=-JP,JP
        do md=-1,1
          write(6,3000) IAL2,JP,MJP,MD
 3000  format(' IAL2,JP,MJP,MD= ',6I4)
             do ipp=1,NPP
              WRITE(6,3001) PP(ipp),DENSITY(ipp,IAL2,MJP,MD)
 3001  format(F14.8,2E16.8)
             enddo ! ipp
            write(6,*)
        ENDDO ! md
       ENDDO ! mjp
      ENDDO ! IAL2
C
C
            write(6,*)
            write(6,*)
            write(6,*)
C
         write(6,*) ' DENSITY: mjp=md ' 
      DO IAL2=1,NINDEX2
       CALL UNPACK1(4,LP,ISP,JP,ITP,INDEX2(IAL2),IZ,IZ,IZ,IZ)
!       do mjp=-JP,JP
        do md=-1,1
         mjp=md
          write(6,4000) IAL2,JP,MJP,MD
 4000  format(' IAL2,JP,MJP,MD= ',6I4)
             do ipp=1,NPP
              WRITE(6,4001) PP(ipp),DENSITY(ipp,IAL2,MJP,MD)
 4001  format(F14.8,2E16.8)
             enddo ! ipp
            write(6,*)
        ENDDO ! md
!      ENDDO ! mjp
      ENDDO ! IAL2
C
C
            write(6,*)
            write(6,*)
            write(6,*)
C
      RETURN
      END

      double complex function RHO(m1p2,m12,t,Qbig,p1p,p1)
      implicit none
      integer m1p2,m12,t
      double precision Qbig,p1p(3),p1(3)
      double precision GEP,GEN,GMP,GMN
      double precision PROTON(0:1),NEUTRON(0:1)
      double precision FMM1,MN,MP,MM,Mnuc
      PARAMETER (FMM1=197.327D0,
     &           MN=939.5653D0,MP=938.2720D0,MM=2.0D0*MP*MN/(MP+MN),
     &           Mnuc=MM/FMM1)
      integer cr1,cr2
      PARAMETER (cr1=0,cr2=0)
      PROTON(0)= 0.5
      PROTON(1)= 0.5
      NEUTRON(0)= 0.5
      NEUTRON(1)= -0.5
      GEP= 1.0
      GEN= 0.0
      GMP= 2.793
      GMN=-1.913
c
      if (m1p2.eq.-1 .and. m12.eq.-1) then
                RHO= GEN*NEUTRON(t) + GEP*PROTON(t) - 
     -   (cr1*Qbig**2*(GEN*NEUTRON(t) + GEP*PROTON(t)))/
     -    (8.*Mnuc**2) + (Complex(0,0.25)*cr2*
     -      (p1(2)*p1p(1) - p1(1)*p1p(2))*
     -      ((GEN - 2*GMN)*NEUTRON(t) + (GEP - 2*GMP)*PROTON(t))
     -      )/Mnuc**2
       else if (m1p2.eq.-1 .and. m12.eq.1) then
                RHO= -(cr2*
     -      (p1(3)*(p1p(1) + Complex(0,1)*p1p(2)) - 
     -        (p1(1) + Complex(0,1)*p1(2))*p1p(3))*
     -      ((GEN - 2*GMN)*NEUTRON(t) + (GEP - 2*GMP)*PROTON(t))
     -      )/(4.*Mnuc**2)
       else if (m1p2.eq.1 .and. m12.eq.-1) then
                RHO= (cr2*
     -     (p1(3)*(p1p(1) - Complex(0,1)*p1p(2)) - 
     -       (p1(1) - Complex(0,1)*p1(2))*p1p(3))*
     -     ((GEN - 2*GMN)*NEUTRON(t) + (GEP - 2*GMP)*PROTON(t)))
     -    /(4.*Mnuc**2)
       else if (m1p2.eq.1 .and. m12.eq.1) then
                RHO= GEN*NEUTRON(t) + GEP*PROTON(t) - 
     -   (cr1*Qbig**2*(GEN*NEUTRON(t) + GEP*PROTON(t)))/
     -    (8.*Mnuc**2) - (Complex(0,0.25)*cr2*
     -      (p1(2)*p1p(1) - p1(1)*p1p(2))*
     -      ((GEN - 2*GMN)*NEUTRON(t) + (GEP - 2*GMP)*PROTON(t))
     -      )/Mnuc**2
        else 
          stop 'in RHO'
        endif
      return
      end
 



      double complex function JEMPLUS(m1p2,m12,t,Qbig,p1p,p1)
      implicit none
      integer m1p2,m12,t
      double precision Qbig,p1p(3),p1(3)
      double precision GEP,GEN,GMP,GMN
      double precision PROTON(0:1),NEUTRON(0:1)
      double precision FMM1,MN,MP,MM,Mnuc
      PARAMETER (FMM1=197.327D0,
     &           MN=939.5653D0,MP=938.2720D0,MM=2.0D0*MP*MN/(MP+MN),
     &           Mnuc=MM/FMM1)
      PROTON(0)= 0.5
      PROTON(1)= 0.5
      NEUTRON(0)= 0.5
      NEUTRON(1)= -0.5
      GEP= 1.0
      GEN= 0.0
      GMP= 2.793
      GMN=-1.913
c
      if (m1p2.eq.-1 .and. m12.eq.-1) then
                JEMPLUS= 
     -  (NEUTRON(t)*(GMN*(p1(1) + Complex(0,1)*p1(2) - p1p(1) - 
     -           Complex(0,1)*p1p(2)) - 
     -        GEN*(p1(1) + Complex(0,1)*p1(2) + p1p(1) + 
     -           Complex(0,1)*p1p(2))) - 
     -     (GMP*(-p1(1) - Complex(0,1)*p1(2) + p1p(1) + 
     -           Complex(0,1)*p1p(2)) + 
     -        GEP*(p1(1) + Complex(0,1)*p1(2) + p1p(1) + 
     -           Complex(0,1)*p1p(2)))*PROTON(t))/
     -   (2.*DSQRT(2.0D0)*Mnuc)
       else if (m1p2.eq.-1 .and. m12.eq.1) then
        JEMPLUS= 0
       else if (m1p2.eq.1 .and. m12.eq.-1) then
                JEMPLUS= 
     -  ((p1(3) - p1p(3))*(GMN*NEUTRON(t) + GMP*PROTON(t)))/
     -   (DSQRT(2.0D0)*Mnuc)
       else if (m1p2.eq.1 .and. m12.eq.1) then
                JEMPLUS= 
     -  (-(NEUTRON(t)*(GMN*
     -           (p1(1) + Complex(0,1)*p1(2) - p1p(1) - 
     -             Complex(0,1)*p1p(2)) + 
     -          GEN*(p1(1) + Complex(0,1)*p1(2) + p1p(1) + 
     -             Complex(0,1)*p1p(2)))) - 
     -     (GMP*(p1(1) + Complex(0,1)*p1(2) - p1p(1) - 
     -           Complex(0,1)*p1p(2)) + 
     -        GEP*(p1(1) + Complex(0,1)*p1(2) + p1p(1) + 
     -           Complex(0,1)*p1p(2)))*PROTON(t))/
     -   (2.*DSQRT(2.0D0)*Mnuc)
        else 
          stop 'in JEMPLUS'
        endif
      return
      end
 



      double complex function JEMMINUS(m1p2,m12,t,Qbig,p1p,p1)
      implicit none
      integer m1p2,m12,t
      double precision Qbig,p1p(3),p1(3)
      double precision GEP,GEN,GMP,GMN
      double precision PROTON(0:1),NEUTRON(0:1)
      double precision FMM1,MN,MP,MM,Mnuc
      PARAMETER (FMM1=197.327D0,
     &           MN=939.5653D0,MP=938.2720D0,MM=2.0D0*MP*MN/(MP+MN),
     &           Mnuc=MM/FMM1)
      PROTON(0)= 0.5
      PROTON(1)= 0.5
      NEUTRON(0)= 0.5
      NEUTRON(1)= -0.5
      GEP= 1.0
      GEN= 0.0
      GMP= 2.793
      GMN=-1.913
c
      if (m1p2.eq.-1 .and. m12.eq.-1) then
                JEMMINUS= 
     -  (NEUTRON(t)*(GEN*(p1(1) - Complex(0,1)*p1(2) + p1p(1) - 
     -           Complex(0,1)*p1p(2)) + 
     -        GMN*(p1(1) - Complex(0,1)*p1(2) - p1p(1) + 
     -           Complex(0,1)*p1p(2))) + 
     -     (GEP*(p1(1) - Complex(0,1)*p1(2) + p1p(1) - 
     -           Complex(0,1)*p1p(2)) + 
     -        GMP*(p1(1) - Complex(0,1)*p1(2) - p1p(1) + 
     -           Complex(0,1)*p1p(2)))*PROTON(t))/
     -   (2.*DSQRT(2.0D0)*Mnuc)
       else if (m1p2.eq.-1 .and. m12.eq.1) then
                JEMMINUS= 
     -  ((p1(3) - p1p(3))*(GMN*NEUTRON(t) + GMP*PROTON(t)))/
     -   (DSQRT(2.0D0)*Mnuc)
       else if (m1p2.eq.1 .and. m12.eq.-1) then
        JEMMINUS= 0
       else if (m1p2.eq.1 .and. m12.eq.1) then
                JEMMINUS= 
     -  (NEUTRON(t)*(GEN*p1(1) - GMN*p1(1) - 
     -        Complex(0,1)*GEN*p1(2) + Complex(0,1)*GMN*p1(2) + 
     -        GEN*p1p(1) + GMN*p1p(1) - 
     -        Complex(0,1)*(GEN + GMN)*p1p(2)) + 
     -     (GEP*p1(1) - GMP*p1(1) - Complex(0,1)*GEP*p1(2) + 
     -        Complex(0,1)*GMP*p1(2) + GEP*p1p(1) + 
     -        GMP*p1p(1) - Complex(0,1)*(GEP + GMP)*p1p(2))*
     -      PROTON(t))/(2.*DSQRT(2.0D0)*Mnuc)
        else 
          stop 'in JEMMINUS'
        endif
      return
      end
 
      double complex function JEMZ(m1p2,m12,t,Qbig,p1p,p1)
      implicit none
      integer m1p2,m12,t
      double precision Qbig,p1p(3),p1(3)
      double precision GEP,GEN,GMP,GMN
      double precision PROTON(0:1),NEUTRON(0:1)
      double precision FMM1,MN,MP,MM,Mnuc
      PARAMETER (FMM1=197.327D0,
     &           MN=939.5653D0,MP=938.2720D0,MM=2.0D0*MP*MN/(MP+MN),
     &           Mnuc=MM/FMM1)
      PROTON(0)= 0.5
      PROTON(1)= 0.5
      NEUTRON(0)= 0.5
      NEUTRON(1)= -0.5
      GEP= 1.0
      GEN= 0.0
      GMP= 2.793
      GMN=-1.913
c
       if (m1p2.eq.-1 .and. m12.eq.-1) then
                JEMZ= ((p1(3) + p1p(3))*
     -     (GEN*NEUTRON(t) + GEP*PROTON(t)))/(2.*Mnuc)
       else if (m1p2.eq.-1 .and. m12.eq.1) then
                JEMZ= -((p1(1) + 
     -        Complex(0,1)*
     -         (p1(2) + Complex(0,1)*p1p(1) - p1p(2)))*
     -      (GMN*NEUTRON(t) + GMP*PROTON(t)))/(2.*Mnuc)
       else if (m1p2.eq.1 .and. m12.eq.-1) then
                JEMZ= ((p1(1) - Complex(0,1)*p1(2) - p1p(1) + 
     -       Complex(0,1)*p1p(2))*
     -     (GMN*NEUTRON(t) + GMP*PROTON(t)))/(2.*Mnuc)
       else if (m1p2.eq.1 .and. m12.eq.1) then
                JEMZ= ((p1(3) + p1p(3))*
     -     (GEN*NEUTRON(t) + GEP*PROTON(t)))/(2.*Mnuc)
        else 
          stop 'in JEMZ'
        endif
      return
      end

      SUBROUTINE MAGNETIC_MULTIPOLES(QBIG)
      implicit none
      double precision QBIG

      integer JMAX,LCHMX2,LMAX
      PARAMETER (JMAX=4,LCHMX2=2+4*JMAX,LMAX=JMAX+2) ! All, t=0 and t=1 channels

      integer NDEUT
      PARAMETER (NDEUT=100)
      double precision PDEUT,WDEUT ! points and weights for the deuteron
      double precision SWAV,DWAV
      COMMON /DEUTERON/ PDEUT(NDEUT),WDEUT(NDEUT),
     &                  SWAV(NDEUT),DWAV(NDEUT) 

      integer NPHIP,NTHETAP,NPHIQ,NTHETAQ
      parameter (NPHIP=15,NTHETAP=15,NPHIQ=NPHIP,NTHETAQ=NTHETAP)
      double precision PHIM,WPHIM
      double precision THETAM,WTHETAM
      COMMON /INTANGLES/ PHIM(NPHIP),WPHIM(NPHIP),
     &                   THETAM(NTHETAP),WTHETAM(NTHETAP)

      integer NPP
      PARAMETER (NPP=46)
      double precision PP,WPP,PBAR ! pedy_p, zerowy to p0
      COMMON/PUNKTYPP/ PP(NPP),WPP(NPP), PBAR  ! pedy_p, zerowy to p0

      integer ithq,iphq
      integer ithp,iphp
      integer ipp,interp
      double precision thq,phq,p,php,thp,P1FM,SPL(NDEUT)
      double precision vecq(3),vec(3),theta,phi
      double precision wavinterp
      COMMON /DEUINTRP/ wavinterp(NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ,2)

      double complex YDEU
      COMMON /YDEUT/ YDEU(NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ,0:2,-2:2)

      double complex YLM
      COMMON /YPQ/ YLM(NTHETAP,NPHIP,0:LMAX,-LMAX:LMAX)

      double precision  VECP1 ! right 
      double precision VECP1P ! left
      COMMON /VEC1/ VECP1(3,NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ),
     &             VECP1P(3,NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ)

      LOGICAL COUPL,NNNPD2
      integer NINDEX2,INDEX2
      COMMON /PWTB2/ INDEX2(LCHMX2),NINDEX2,COUPL(LCHMX2),NNNPD2(LCHMX2)

      double complex MAGJM
      COMMON /MAGNETIC/ MAGJM(NPP,LCHMX2,-JMAX:JMAX,-1:1,0:LMAX,1)

      double precision CL1,CL2,CL3,CL4,PI

      double precision CLEBSCH
      external CLEBSCH

      double complex JEMPLUS
      external JEMPLUS
      double complex JEMMINUS
      external JEMMINUS
      double complex JEMZ
      external JEMZ

      integer IAL2,LP,ISP,JP,ITP,IZ
      integer MLD,MLP,IALD,IS,L,J,IT
      integer m1p,m1p2,m1d,m1d2,mjp,md
      integer JBIG,MBIG

      PI=DACOS(-1.0d0)

      DO IAL2=1,NINDEX2
       CALL UNPACK1(4,LP,ISP,JP,ITP,INDEX2(IAL2),IZ,IZ,IZ,IZ)
        do mjp=-JP,JP
          do mlp=-LP,LP
            IF (abs(mjp-mlp).le.isp) THEN
             CL1=1.0d0/(4.0d0*Pi)*
     &         CLEBSCH(2*lp,2*isp,2*jp,2*mlp,2*(mjp-mlp),2*mjp)
              DO IALD=1,2
               L=2*IALD-2
               IS=1
               J=1
               IT=0
                do md=-1,1
                 do mld=-L,L
                  IF (abs(md-mld).le.is) THEN
                   CL2=CL1*CLEBSCH(2*L,2,2,2*mld,2*(md-mld),2*md)
                    do m1p=1,2
                     m1p2=2*m1p-3
                      IF (ABS(2*(mjp-mlp)-m1p2).le.1) THEN
                        CL3=CL2*
     &           CLEBSCH(1,1,2*ISP,m1p2,2*(mjp-mlp)-m1p2,2*(mjp-mlp))
                       do m1d=1,2
                        m1d2=2*m1d-3
                           IF ((ABS(2*(md-mld)-m1d2).le.1).and.
     &            ((2*(mjp-mlp)-m1p2).eq.(2*(md-mld)-m1d2))) THEN
                         CL4=CL3*
     &           CLEBSCH(1,1,2,m1d2,2*(md-mld)-m1d2,2*(md-mld))
      DO JBIG=1,JMAX+1
!      DO MBIG=-JBIG,JBIG
       DO MBIG=1,1
        do ithq=1,NTHETAQ
         do iphq=1,NPHIQ
          do ipp=1,NPP
           do ithp=1,NTHETAP
            do iphp=1,NPHIP
               MAGJM(ipp,IAL2,MJP,MD,JBIG,MBIG)=
     &             MAGJM(ipp,IAL2,MJP,MD,JBIG,MBIG) +
     &        CL4*YLM(ithp,iphp,LP,MLP)*
     &        YDEU(ipp,ithp,iphp,ithq,iphq,L,mld)*
     &        wavinterp(ipp,ithp,iphp,ithq,iphq,IALD)*
     &     (
     &        -1.0d0/DSQRT(2.0d0)*
     &        DCONJG(YLM(ithq,iphq,JBIG,0))*
     &        JEMPLUS(m1p2,m1d2,ITP,Qbig,
     &            VECP1P(1,ipp,ithp,iphp,ithq,iphq),
     &             VECP1(1,ipp,ithp,iphp,ithq,iphq))
     &       + 1.0d0/DSQRT(JBIG*(JBIG+1.0d0))*
     &        DCONJG(YLM(ithq,iphq,JBIG,1))*
     &        JEMZ(m1p2,m1d2,ITP,Qbig,
     &            VECP1P(1,ipp,ithp,iphp,ithq,iphq),
     &             VECP1(1,ipp,ithp,iphp,ithq,iphq))
     &       + DSQRT((JBIG**2+JBIG-2.0d0)/(2*JBIG**2+2*JBIG))*
     &        DCONJG(YLM(ithq,iphq,JBIG,2))*
     &        JEMMINUS(m1p2,m1d2,ITP,Qbig,
     &            VECP1P(1,ipp,ithp,iphp,ithq,iphq),
     &             VECP1(1,ipp,ithp,iphp,ithq,iphq))
     &     )
            enddo ! iphp
           enddo ! ithp
          enddo ! ipp
         enddo ! iphq
        enddo ! ithq
       enddo ! MBIG
       enddo ! JBIG
c
c
                           ENDIF
                       enddo ! m1d
                      ENDIF 
                    enddo ! m1p
                   ENDIF 
              enddo ! mld
             enddo ! md
           ENDDO ! IALD
            ENDIF 
          enddo ! mlp
        enddo ! mjp
      ENDDO ! ial2
C
C
      write(6,*) ' MAGJM for M=1 '
      DO IAL2=1,NINDEX2
       CALL UNPACK1(4,LP,ISP,JP,ITP,INDEX2(IAL2),IZ,IZ,IZ,IZ)
        do mjp=-JP,JP
        do md=-1,1
         DO JBIG=1,JMAX+1
!         DO MBIG=-JBIG,JBIG
          DO MBIG=1,1
            write(6,2000) IAL2,JP,MJP,MD,JBIG,MBIG
 2000  format(' IAL2,JP,MJP,MD,JBIG,MBIG= ',6I4)
             do ipp=1,NPP
               WRITE(6,2001) PP(ipp), MAGJM(ipp,IAL2,MJP,MD,JBIG,MBIG)
 2001  format(F14.8,2E16.8)
             enddo ! ipp
            write(6,*) 
          ENDDO ! MBIG
         ENDDO ! JBIG
        ENDDO ! md
       ENDDO ! mjp
      ENDDO ! IAL2
C
C
            write(6,*) 
            write(6,*) 
            write(6,*) 
C
      RETURN
      END

      SUBROUTINE ELECTRIC_MULTIPOLES(QBIG)
      implicit none
      double precision QBIG

      integer JMAX,LCHMX2,LMAX
      PARAMETER (JMAX=4,LCHMX2=2+4*JMAX,LMAX=JMAX+2) ! All, t=0 and t=1 channels

      integer NDEUT
      PARAMETER (NDEUT=100)
      double precision PDEUT,WDEUT ! points and weights for the deuteron
      double precision SWAV,DWAV
      COMMON /DEUTERON/ PDEUT(NDEUT),WDEUT(NDEUT),
     &                  SWAV(NDEUT),DWAV(NDEUT) 

      integer NPHIP,NTHETAP,NPHIQ,NTHETAQ
      parameter (NPHIP=15,NTHETAP=15,NPHIQ=NPHIP,NTHETAQ=NTHETAP)
      double precision PHIM,WPHIM
      double precision THETAM,WTHETAM
      COMMON /INTANGLES/ PHIM(NPHIP),WPHIM(NPHIP),
     &                   THETAM(NTHETAP),WTHETAM(NTHETAP)

      integer NPP
      PARAMETER (NPP=46)
      double precision PP,WPP,PBAR ! pedy_p, zerowy to p0
      COMMON/PUNKTYPP/ PP(NPP),WPP(NPP), PBAR  ! pedy_p, zerowy to p0

      integer ithq,iphq
      integer ithp,iphp
      integer ipp,interp
      double precision thq,phq,p,php,thp,P1FM,SPL(NDEUT)
      double precision vecq(3),vec(3),theta,phi
      double precision wavinterp
      COMMON /DEUINTRP/ wavinterp(NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ,2)

      double complex YDEU
      COMMON /YDEUT/ YDEU(NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ,0:2,-2:2)

      double complex YLM
      COMMON /YPQ/ YLM(NTHETAP,NPHIP,0:LMAX,-LMAX:LMAX)

      double precision  VECP1 ! right 
      double precision VECP1P ! left
      COMMON /VEC1/ VECP1(3,NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ),
     &             VECP1P(3,NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ)

      LOGICAL COUPL,NNNPD2
      integer NINDEX2,INDEX2
      COMMON /PWTB2/ INDEX2(LCHMX2),NINDEX2,COUPL(LCHMX2),NNNPD2(LCHMX2)

      double complex MAGJM
      COMMON /MAGNETIC/ MAGJM(NPP,LCHMX2,-JMAX:JMAX,-1:1,0:LMAX,1)

      double complex ELJM
      COMMON /ELECTRIC/ ELJM(NPP,LCHMX2,-JMAX:JMAX,-1:1,0:LMAX,1)
      
      double complex JPL(NPP,LCHMX2,-JMAX:JMAX,-1:1)

      double precision CL1,CL2,CL3,CL4,PI

      double precision CLEBSCH
      external CLEBSCH

      double complex JEMPLUS
      external JEMPLUS
      double complex JEMMINUS
      external JEMMINUS
      double complex JEMZ
      external JEMZ

      integer IAL2,LP,ISP,JP,ITP,IZ
      integer MLD,MLP,IALD,IS,L,J,IT
      integer m1p,m1p2,m1d,m1d2,mjp,md
      integer JBIG,MBIG

      PI=DACOS(-1.0d0)

      DO IAL2=1,NINDEX2
       CALL UNPACK1(4,LP,ISP,JP,ITP,INDEX2(IAL2),IZ,IZ,IZ,IZ)
        do mjp=-JP,JP
          do mlp=-LP,LP
            IF (abs(mjp-mlp).le.isp) THEN
             CL1=-1.0d0/(4.0d0*Pi)*
     &         CLEBSCH(2*lp,2*isp,2*jp,2*mlp,2*(mjp-mlp),2*mjp)
              DO IALD=1,2
               L=2*IALD-2
               IS=1
               J=1
               IT=0
                do md=-1,1
                 do mld=-L,L
                  IF (abs(md-mld).le.is) THEN
                   CL2=CL1*CLEBSCH(2*L,2,2,2*mld,2*(md-mld),2*md)
                    do m1p=1,2
                     m1p2=2*m1p-3
                      IF (ABS(2*(mjp-mlp)-m1p2).le.1) THEN
                        CL3=CL2*
     &           CLEBSCH(1,1,2*ISP,m1p2,2*(mjp-mlp)-m1p2,2*(mjp-mlp))
                       do m1d=1,2
                        m1d2=2*m1d-3
                           IF ((ABS(2*(md-mld)-m1d2).le.1).and.
     &            ((2*(mjp-mlp)-m1p2).eq.(2*(md-mld)-m1d2))) THEN
                         CL4=CL3*
     &           CLEBSCH(1,1,2,m1d2,2*(md-mld)-m1d2,2*(md-mld))
      DO JBIG=1,JMAX+1
!      DO MBIG=-JBIG,JBIG
       DO MBIG=1,1
        do ithq=1,NTHETAQ
         do iphq=1,NPHIQ
          do ipp=1,NPP
           do ithp=1,NTHETAP
            do iphp=1,NPHIP
               ELJM(ipp,IAL2,MJP,MD,JBIG,MBIG)=
     &             ELJM(ipp,IAL2,MJP,MD,JBIG,MBIG) +
     &        CL4*YLM(ithp,iphp,LP,MLP)*
     &        YDEU(ipp,ithp,iphp,ithq,iphq,L,mld)*
     &        wavinterp(ipp,ithp,iphp,ithq,iphq,IALD)*
     &   (
     &  (((1 + JBIG)*
     &     DCONJG(YLM(ithq,iphq,-1 + JBIG,
     &       0)))/Sqrt(-2.0d0 + 8*JBIG**2) + 
     &  (JBIG*DCONJG(YLM(ithq,iphq,1 + JBIG,
     &       0)))/
     &   Sqrt(6.0d0 + 16*JBIG + 8*JBIG**2))*
     &         JEMPLUS(m1p2,m1d2,ITP,Qbig,
     &            VECP1P(1,ipp,ithp,iphp,ithq,iphq),
     &             VECP1(1,ipp,ithp,iphp,ithq,iphq))
     &  +  ((1 + JBIG)*
     &   Sqrt((1.0d0 - JBIG)/
     &     (JBIG - 4*JBIG**3))*
     &   DCONJG(YLM(ithq,iphq,-1 + JBIG,1))
     &   - JBIG*Sqrt((2.0d0 + JBIG)/
     &     (3.0d0 + 11*JBIG + 12*JBIG**2 + 
     &       4*JBIG**3))*
     &   DCONJG(YLM(ithq,iphq,1 + JBIG,1)))*
     &        JEMZ(m1p2,m1d2,ITP,Qbig,
     &            VECP1P(1,ipp,ithp,iphp,ithq,iphq),
     &             VECP1(1,ipp,ithp,iphp,ithq,iphq))
     &  + ((Sqrt(((1.0d0 + JBIG)*
     &         (2.0d0 - 3*JBIG + JBIG**2))/
     &       (JBIG*(-1.0d0 + 4*JBIG**2)))*
     &     DCONJG(YLM(ithq,iphq,-1 + JBIG,
     &       2)))/Sqrt(2.0d0) + 
     &  Sqrt((JBIG*(6.0d0 + 5*JBIG + JBIG**2))/
     &     (6.0d0 + 22*JBIG + 24*JBIG**2 + 
     &       8*JBIG**3))*
     &   DCONJG(YLM(ithq,iphq,1 + JBIG,2)))*
     &        JEMMINUS(m1p2,m1d2,ITP,Qbig,
     &            VECP1P(1,ipp,ithp,iphp,ithq,iphq),
     &             VECP1(1,ipp,ithp,iphp,ithq,iphq))
     &  )
            enddo ! iphp
           enddo ! ithp
          enddo ! ipp
         enddo ! iphq
        enddo ! ithq
       enddo ! MBIG
       enddo ! JBIG
c
c
                           ENDIF
                       enddo ! m1d
                      ENDIF 
                    enddo ! m1p
                   ENDIF 
              enddo ! mld
             enddo ! md
           ENDDO ! IALD
            ENDIF 
          enddo ! mlp
        enddo ! mjp
      ENDDO ! ial2
C
C
      write(6,*) ' ELJM for M=1 '
      DO IAL2=1,NINDEX2
       CALL UNPACK1(4,LP,ISP,JP,ITP,INDEX2(IAL2),IZ,IZ,IZ,IZ)
        do mjp=-JP,JP
        do md=-1,1
         DO JBIG=1,JMAX+1
!         DO MBIG=-JBIG,JBIG
          DO MBIG=1,1
            write(6,2000) IAL2,JP,MJP,MD,JBIG,MBIG
 2000  format(' IAL2,JP,MJP,MD,JBIG,MBIG= ',6I4)
             do ipp=1,NPP
               WRITE(6,2001) PP(ipp), ELJM(ipp,IAL2,MJP,MD,JBIG,MBIG)
 2001  format(F14.8,2E16.8)
             enddo ! ipp
            write(6,*) 
          ENDDO ! MBIG
         ENDDO ! JBIG
        ENDDO ! md
       ENDDO ! mjp
      ENDDO ! IAL2
C
            write(6,*) 
            write(6,*) 
            write(6,*) 
C
C
      JPL=(0.0d0,0.0d0)
      DO IAL2=1,NINDEX2
       CALL UNPACK1(4,LP,ISP,JP,ITP,INDEX2(IAL2),IZ,IZ,IZ,IZ)
        do mjp=-JP,JP
        do md=-1,1
         DO JBIG=max(abs(JP-1),1),JP+1
!         DO JBIG=0,JMAX+1
!          DO MBIG=-JBIG,JBIG
             do ipp=1,NPP
              JPL(ipp,IAL2,MJP,MD)=JPL(ipp,IAL2,MJP,MD) - 
     &         DSQRT(2*Pi)*DSQRT(2*JBIG+1.0d0)*
     &   ( MAGJM(ipp,IAL2,MJP,MD,JBIG,1) +
     &     ELJM(ipp,IAL2,MJP,MD,JBIG,1) )
             enddo ! ipp
!          ENDDO ! MBIG
         ENDDO ! JBIG
        ENDDO ! md
       ENDDO ! mjp
      ENDDO ! IAL2
C
C
            write(6,*) 
            write(6,*) 
            write(6,*)
C
         write(6,*) ' JPL: general mjp ' 
      DO IAL2=1,NINDEX2
       CALL UNPACK1(4,LP,ISP,JP,ITP,INDEX2(IAL2),IZ,IZ,IZ,IZ)
        do mjp=-JP,JP
        do md=-1,1
          write(6,3000) IAL2,JP,MJP,MD
 3000  format(' IAL2,JP,MJP,MD= ',6I4)
             do ipp=1,NPP
              WRITE(6,3001) PP(ipp),JPL(ipp,IAL2,MJP,MD)
 3001  format(F14.8,2E16.8)
             enddo ! ipp
            write(6,*)
        ENDDO ! md
       ENDDO ! mjp
      ENDDO ! IAL2
C
C
            write(6,*)
            write(6,*)
            write(6,*)
C
         write(6,*) ' JPL: mjp=md+1 ' 
      DO IAL2=1,NINDEX2
       CALL UNPACK1(4,LP,ISP,JP,ITP,INDEX2(IAL2),IZ,IZ,IZ,IZ)
!       do mjp=-JP,JP
        do md=-1,1
         mjp=md+1
          write(6,4000) IAL2,JP,MJP,MD
 4000  format(' IAL2,JP,MJP,MD= ',6I4)
             do ipp=1,NPP
              WRITE(6,4001) PP(ipp),JPL(ipp,IAL2,MJP,MD)
 4001  format(F14.8,2E16.8)
             enddo ! ipp
            write(6,*)
        ENDDO ! md
!      ENDDO ! mjp
      ENDDO ! IAL2
C
C
            write(6,*)
            write(6,*)
            write(6,*)
C
      RETURN
      END

      SUBROUTINE ELECTRIC_MULTIPOLES_SIEGERT(QBIG)
      implicit none
      double precision QBIG

      integer JMAX,LCHMX2,LMAX
      PARAMETER (JMAX=4,LCHMX2=2+4*JMAX,LMAX=JMAX+2) ! All, t=0 and t=1 channels

      integer NDEUT
      PARAMETER (NDEUT=100)
      double precision PDEUT,WDEUT ! points and weights for the deuteron
      double precision SWAV,DWAV
      COMMON /DEUTERON/ PDEUT(NDEUT),WDEUT(NDEUT),
     &                  SWAV(NDEUT),DWAV(NDEUT) 

      integer NPHIP,NTHETAP,NPHIQ,NTHETAQ
      parameter (NPHIP=15,NTHETAP=15,NPHIQ=NPHIP,NTHETAQ=NTHETAP)
      double precision PHIM,WPHIM
      double precision THETAM,WTHETAM
      COMMON /INTANGLES/ PHIM(NPHIP),WPHIM(NPHIP),
     &                   THETAM(NTHETAP),WTHETAM(NTHETAP)

      integer NPP
      PARAMETER (NPP=46)
      double precision PP,WPP,PBAR ! pedy_p, zerowy to p0
      COMMON/PUNKTYPP/ PP(NPP),WPP(NPP), PBAR  ! pedy_p, zerowy to p0

      integer ithq,iphq
      integer ithp,iphp
      integer ipp,interp
      double precision thq,phq,p,php,thp,P1FM,SPL(NDEUT)
      double precision vecq(3),vec(3),theta,phi
      double precision wavinterp
      COMMON /DEUINTRP/ wavinterp(NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ,2)

      double complex YDEU
      COMMON /YDEUT/ YDEU(NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ,0:2,-2:2)

      double complex YLM
      COMMON /YPQ/ YLM(NTHETAP,NPHIP,0:LMAX,-LMAX:LMAX)

      double precision  VECP1 ! right 
      double precision VECP1P ! left
      COMMON /VEC1/ VECP1(3,NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ),
     &             VECP1P(3,NPP,NTHETAP,NPHIP,NTHETAQ,NPHIQ)

      LOGICAL COUPL,NNNPD2
      integer NINDEX2,INDEX2
      COMMON /PWTB2/ INDEX2(LCHMX2),NINDEX2,COUPL(LCHMX2),NNNPD2(LCHMX2)

      double complex MAGJM
      COMMON /MAGNETIC/ MAGJM(NPP,LCHMX2,-JMAX:JMAX,-1:1,0:LMAX,1)

      double complex ELJM
      COMMON /ELECTRIC2/ ELJM(NPP,LCHMX2,-JMAX:JMAX,-1:1,0:LMAX,1)
      
      double complex JPL(NPP,LCHMX2,-JMAX:JMAX,-1:1)
      double precision PLUS(NPP,LCHMX2,-1:1)

      double precision CL1,CL2,CL3,CL4,PI

      double precision CLEBSCH
      external CLEBSCH

      double complex RHO
      external RHO
      double complex JEMPLUS
      external JEMPLUS
      double complex JEMMINUS
      external JEMMINUS
      double complex JEMZ
      external JEMZ

      integer IAL2,LP,ISP,JP,ITP,IZ
      integer MLD,MLP,IALD,IS,L,J,IT
      integer m1p,m1p2,m1d,m1d2,mjp,md
      integer JBIG,MBIG

      PI=DACOS(-1.0d0)

      DO IAL2=1,NINDEX2
       CALL UNPACK1(4,LP,ISP,JP,ITP,INDEX2(IAL2),IZ,IZ,IZ,IZ)
        do mjp=-JP,JP
          do mlp=-LP,LP
            IF (abs(mjp-mlp).le.isp) THEN
             CL1=-1.0d0/(4.0d0*Pi)*
     &         CLEBSCH(2*lp,2*isp,2*jp,2*mlp,2*(mjp-mlp),2*mjp)
              DO IALD=1,2
               L=2*IALD-2
               IS=1
               J=1
               IT=0
                do md=-1,1
                 do mld=-L,L
                  IF (abs(md-mld).le.is) THEN
                   CL2=CL1*CLEBSCH(2*L,2,2,2*mld,2*(md-mld),2*md)
                    do m1p=1,2
                     m1p2=2*m1p-3
                      IF (ABS(2*(mjp-mlp)-m1p2).le.1) THEN
                        CL3=CL2*
     &           CLEBSCH(1,1,2*ISP,m1p2,2*(mjp-mlp)-m1p2,2*(mjp-mlp))
                       do m1d=1,2
                        m1d2=2*m1d-3
                           IF ((ABS(2*(md-mld)-m1d2).le.1).and.
     &            ((2*(mjp-mlp)-m1p2).eq.(2*(md-mld)-m1d2))) THEN
                         CL4=CL3*
     &           CLEBSCH(1,1,2,m1d2,2*(md-mld)-m1d2,2*(md-mld))
      DO JBIG=1,JMAX+1
!      DO MBIG=-JBIG,JBIG
       DO MBIG=1,1
        do ithq=1,NTHETAQ
         do iphq=1,NPHIQ
          do ipp=1,NPP
           do ithp=1,NTHETAP
            do iphp=1,NPHIP
               ELJM(ipp,IAL2,MJP,MD,JBIG,MBIG)=
     &             ELJM(ipp,IAL2,MJP,MD,JBIG,MBIG) +
     &        CL4*YLM(ithp,iphp,LP,MLP)*
     &        YDEU(ipp,ithp,iphp,ithq,iphq,L,mld)*
     &        wavinterp(ipp,ithp,iphp,ithq,iphq,IALD)*
     &   (
     &      Sqrt((1.0d0 + 2*JBIG)/(6.0d0 + 4*JBIG))*
     &      DCONJG(YLM(ithq,iphq,1 + JBIG,0))*
     &         JEMPLUS(m1p2,m1d2,ITP,Qbig,
     &            VECP1P(1,ipp,ithp,iphp,ithq,iphq),
     &             VECP1(1,ipp,ithp,iphp,ithq,iphq))
     &    -    Sqrt(((2.0d0 + JBIG)*(1.0d0 + 2*JBIG))/
     &    (3.0d0 + 5*JBIG + 2*JBIG**2))*
     &      DCONJG(YLM(ithq,iphq,1 + JBIG,1))*
     &         JEMZ(m1p2,m1d2,ITP,Qbig,
     &            VECP1P(1,ipp,ithp,iphp,ithq,iphq),
     &             VECP1(1,ipp,ithp,iphp,ithq,iphq))
     &  +  (Sqrt(((1.0d0 + 2*JBIG)*
     &      (6.0d0 + 5*JBIG + JBIG**2))/
     &    (JBIG*(3.0d0 + 5*JBIG + 2*JBIG**2)))/Sqrt(2.0d0))* 
     &      DCONJG(YLM(ithq,iphq,1 + JBIG,2))*
     &         JEMMINUS(m1p2,m1d2,ITP,Qbig,
     &            VECP1P(1,ipp,ithp,iphp,ithq,iphq),
     &             VECP1(1,ipp,ithp,iphp,ithq,iphq))
     &   + Sqrt((JBIG+1.0d0)/JBIG)*
     &      DCONJG(YLM(ithq,iphq,JBIG,1))*
     &         RHO(m1p2,m1d2,ITP,Qbig,
     &            VECP1P(1,ipp,ithp,iphp,ithq,iphq),
     &             VECP1(1,ipp,ithp,iphp,ithq,iphq))
     &  )
            enddo ! iphp
           enddo ! ithp
          enddo ! ipp
         enddo ! iphq
        enddo ! ithq
       enddo ! MBIG
       enddo ! JBIG
c
c
                           ENDIF
                       enddo ! m1d
                      ENDIF 
                    enddo ! m1p
                   ENDIF 
              enddo ! mld
             enddo ! md
           ENDDO ! IALD
            ENDIF 
          enddo ! mlp
        enddo ! mjp
      ENDDO ! ial2
C
C
      write(6,*) ' ELJM2 for M=1 '
      DO IAL2=1,NINDEX2
       CALL UNPACK1(4,LP,ISP,JP,ITP,INDEX2(IAL2),IZ,IZ,IZ,IZ)
        do mjp=-JP,JP
        do md=-1,1
!         DO JBIG=1,JMAX+1
!         DO MBIG=-JBIG,JBIG
          DO JBIG=max(abs(JP-1),1),JP+1
          DO MBIG=1,1
            write(6,2000) IAL2,JP,MJP,MD,JBIG,MBIG
 2000  format(' IAL2,JP,MJP,MD,JBIG,MBIG= ',6I4)
             do ipp=1,NPP
               WRITE(6,2001) PP(ipp), ELJM(ipp,IAL2,MJP,MD,JBIG,MBIG)
 2001  format(F14.8,2E16.8)
             enddo ! ipp
            write(6,*) 
          ENDDO ! MBIG
         ENDDO ! JBIG
        ENDDO ! md
       ENDDO ! mjp
      ENDDO ! IAL2
C
            write(6,*) 
            write(6,*) 
            write(6,*) 
C
C
      JPL=(0.0d0,0.0d0)
      DO IAL2=1,NINDEX2
       CALL UNPACK1(4,LP,ISP,JP,ITP,INDEX2(IAL2),IZ,IZ,IZ,IZ)
        do mjp=-JP,JP
        do md=-1,1
         DO JBIG=max(abs(JP-1),1),JP+1
!         DO JBIG=0,JMAX+1
!          DO MBIG=-JBIG,JBIG
             do ipp=1,NPP
              JPL(ipp,IAL2,MJP,MD)=JPL(ipp,IAL2,MJP,MD) - 
     &         DSQRT(2*Pi)*DSQRT(2*JBIG+1.0d0)*
     &   ( MAGJM(ipp,IAL2,MJP,MD,JBIG,1) +
     &     ELJM(ipp,IAL2,MJP,MD,JBIG,1) )
             enddo ! ipp
!          ENDDO ! MBIG
         ENDDO ! JBIG
        ENDDO ! md
       ENDDO ! mjp
      ENDDO ! IAL2
C
C
            write(6,*) 
            write(6,*) 
            write(6,*)
C
         write(6,*) ' JPL2: general mjp ' 
      DO IAL2=1,NINDEX2
       CALL UNPACK1(4,LP,ISP,JP,ITP,INDEX2(IAL2),IZ,IZ,IZ,IZ)
        do mjp=-JP,JP
        do md=-1,1
          write(6,3000) IAL2,JP,MJP,MD
 3000  format(' IAL2,JP,MJP,MD= ',6I4)
             do ipp=1,NPP
              WRITE(6,3001) PP(ipp),JPL(ipp,IAL2,MJP,MD)
 3001  format(F14.8,2E16.8)
             enddo ! ipp
            write(6,*)
        ENDDO ! md
       ENDDO ! mjp
      ENDDO ! IAL2
C
C
            write(6,*)
            write(6,*)
            write(6,*)
C
         write(6,*) ' JPL2: mjp=md+1 ' 
      DO IAL2=1,NINDEX2
       CALL UNPACK1(4,LP,ISP,JP,ITP,INDEX2(IAL2),IZ,IZ,IZ,IZ)
!       do mjp=-JP,JP
        do md=-1,1
         mjp=md+1
          write(6,4000) IAL2,JP,MJP,MD
 4000  format(' IAL2,JP,MJP,MD= ',6I4)
             do ipp=1,NPP
              WRITE(6,4001) PP(ipp),JPL(ipp,IAL2,MJP,MD)
CCCCCCCCC   FINAL RESULT !!!!   CCCCCCCCC
              PLUS(ipp,IAL2,MD)=2*JPL(ipp,IAL2,MJP,MD)
CCCCCCCCC   FINAL RESULT !!!!   CCCCCCCCC
 4001  format(F14.8,2E16.8)
             enddo ! ipp
            write(6,*)
        ENDDO ! md
!      ENDDO ! mjp
      ENDDO ! IAL2
C
C
            WRITE(150,*) PLUS
C
            write(6,*)
            write(6,*)
            write(6,*)
C
      RETURN
      END





