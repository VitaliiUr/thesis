      PROGRAM SCHROED_IMPROVED_CHIRAL
C FROM DR. WITALA
      PARAMETER (NP1=52,NP2=20,NP=NP1+NP2,NL=2,P1=3.,P2=10.,P3=50.,
     .           NDEUT=100,NPL=NP*NL)
      DIMENSION P(NP),WP(NP),SWAV(NDEUT),DWAV(NDEUT),QDEUT(NDEUT),
     .          SK(NP,NL,NP,NL),EWR(NPL),EWI(NPL),EVR(NP,NL),EVI(NP,NL),
     .          WS1(NPL),WS2(NPL),WS3(NPL),WS4(NPL+2,NPL),WS5(NPL),
     .          WS6(NPL),V(NL,NL),SN(NL)
      DIMENSION IS(NPL),IZ(NPL)
      REAL MN,MP,MM,M,FMM1,HQM
      PARAMETER (FMM1=197.327D0,
     &           MN=939.5653D0,MP=938.2720D0,MM=2.0D0*MP*MN/(MP+MN))
      character(len=2) FORCE
      integer OSTAT,CUTNUM
      common /pot/ OSTAT,CUTNUM,FORCE
C
      FORCE='np'
C     OSTAT=4 ! N4LO
c     CUTNUM=2 ! R=0.9 fm
      READ(1,*) OSTAT,CUTNUM
C
      HQM=FMM1**2/MM
      write(6,*) ' adjusted HQM= ',HQM
      ICTR=2 
      IF (ICTR .EQ. 1) THEN
        WRITE (6,*) 'MALFLIET-TJON POTENTIAL'
      ELSE
        WRITE (6,*) 'OTHER POTENTIAL'
      ENDIF
C
      CALL TRNS (NP1,NP2,NP,P1,P2,P3,P,WP)
      OPEN(12,FILE='potential.txt',STATUS='UNKNOWN',FORM='FORMATTED')
      write(12, *) "p, p', V11, V12, V21, V22"
      DO 10 I=1,NP
        DO 20 J=1,NP
          WPJ=WP(J)*P(J)**2
          IF (ICTR .EQ. 1) THEN
            CALL POTMAT1(P(I),P(J),0,1,1,V(1,1))
            V(1,2)=0.
            V(2,1)=0.
            V(2,2)=0.
          ELSE
           CALL POTLSJIEE3(P(I),P(J),0,1,1,V(1,1),V(1,2),V(2,1),V(2,2),
     .                   .TRUE.)
C
C  DIE FOLGENDEN ANWEISUNGEN NUR FUER PARIS- UND BONN-POTENTIAL
C
            PIJ=1./(P(I)*P(J))
            DO 30 L=1,NL
              DO 30 LS=1,NL
   30           V(L,LS)=V(L,LS)*PIJ
C
          ENDIF
          write(12, *) P(I),P(J),V(1,1),V(1,2),V(2,1),V(2,2)
          DO 20 L=1,NL
            DO 20 LS=1,NL
   20         SK(I,L,J,LS)=WPJ*V(L,LS)*HQM
        PM=P(I)**2*HQM
        DO 10 L=1,NL
 10     SK(I,L,I,L)=SK(I,L,I,L)+PM
C
*      Call F02BCF(SK,NPL,NPL,2.2,2.25,1,NEW,EWR,EWI,EVR,NPL,EVI,NPL,
*     .             WS1,WS2,WS3,WS4,NPL+2,WS5,WS6,IF)
C
      EPS1=1.E-15
      WRITE(6,*) ' EPS1= ',EPS1
      ITSMX=30
      CALL GNVEC(SK,NPL,-2.2,EVR,EVI,EPS1,IS,IZ,ITSMX,EWR(1))
      WRITE (6,98)
   98 FORMAT (1H1,///,10X,'PROGRAM SCHROED',///)
      WRITE (6,99) NP1,NP2,P1,P2,P3
   99 FORMAT (1X,'NP1= ',I3,3X,'NP2= ',I3,/,'P1= ',F5.2,3X,
     1           'P2= ',F5.2,3X,'P3= ',F5.2,//)
      WRITE (6,102) EWR(1)
  102 FORMAT (5X,'EIGENWERT= ',F20.13,2X,'MEV',//)
      SG=0.
      DO 40 L=1,NL
        SN(L)=0.
        DO 50 J=1,NP
   50     SN(L)=SN(L)+(EVR(J,L)*P(J))**2*WP(J)
        SN(L)=SN(L)
   40   SG=SG+SN(L)
      DO 45 L=1,NL
        FN=1./SG
        SN(L)=SN(L)*FN
        WRITE (6,104) 2*L-2,SN(L)*100
  104   FORMAT (/,5X,'L=',I1,' -WELLEN-BEIMISCHUNG',F16.8,'%',//)
        DO 45 J=1,NP
   45     EVR(J,L)=EVR(J,L)*SQRT(FN)
      DO 60 J=1,NP
c       WRITE (20)    P(J),EVR(J,1)*P(J),EVR(J,2)*P(J)
   60   WRITE (6,103) P(J),EVR(J,1),EVR(J,2)
  103   FORMAT (3(2X,E16.6))
      OPEN(11,FILE='deuwf.txt',STATUS='UNKNOWN',FORM='FORMATTED')
      WRITE(11, *) 'I,QD,S,D,D2'
      DO 210 I=1,NDEUT
        IF(I.EQ.1) QDEUT(I)=0.001
        IF(I.EQ.2) QDEUT(I)=0.01
        IF(I.GT.2.AND.I.LE.7) QDEUT(I)=0.19/5.*FLOAT(I-2)+0.01
        IF(I.GT.7.AND.I.LE.15) QDEUT(I)=0.8 /8.*FLOAT(I-7)+0.2
        IF(I.GT.15.AND.I.LE.60) QDEUT(I)=9./45. *FLOAT(I-15)+1.0
        IF(I.GT.60.AND.I.LE.80) QDEUT(I)=10./20.*FLOAT(I-60)+10.0
        IF(I.GT.80 ) QDEUT(I)=FLOAT(I-80)+20.0
        S0=0.
        S2=0.
        DO 220 J=1,NP
          IF (ICTR .EQ. 1) THEN
            CALL POTMAT1 (QDEUT(I),P(J),0,1,1,V(1,1))
            V(1,2)=0.
            V(2,1)=0.
            V(2,2)=0.
          ELSE
           CALL POTLSJIEE3(QDEUT(I),P(J),0,1,1,V(1,1),V(1,2),V(2,1),
     .                   V(2,2),.TRUE.)
C
C  DIE FOLGENDEN ANWEISUNGEN NUR FUER PARIS- UND BONN-POTENTIAL
C
            PIJ=1./(QDEUT(I)*P(J))
            DO 230 L=1,NL
              DO 230 LS=1,NL
  230           V(L,LS)=V(L,LS)*PIJ
C
          ENDIF
          S0=S0+P(J)**2*WP(J)*(V(1,1)*EVR(J,1)+V(1,2)*EVR(J,2))*HQM
  220     S2=S2+P(J)**2*WP(J)*(V(2,1)*EVR(J,1)+V(2,2)*EVR(J,2))*HQM
        HELP=1./(EWR(1)-QDEUT(I)**2*HQM)
        SWAV(I)=S0*HELP
        DWAV(I)=S2*HELP
        HELP=DWAV(I)/QDEUT(I)**2
        WRITE (6,676) I,QDEUT(I),SWAV(I),DWAV(I),HELP
  210   WRITE (11,676) I,QDEUT(I),SWAV(I),DWAV(I),HELP
  676   FORMAT(1X, I3,1X,4(E12.5,1X))
C
      OPEN(10,FILE='deuwf',STATUS='UNKNOWN',FORM='UNFORMATTED')
      WRITE(10) QDEUT,EWR(1)
      WRITE(10) SWAV
      WRITE(10) DWAV
      CLOSE (10)
C
      STOP
      END
      SUBROUTINE POTMAT(P,PS,L,IS,JJ,PNIJM)
C     FOR YUKAWA SEPARABLE
      REAL KAPPA
      DIMENSION KAPPA(2),BETA(2)
C     NP 1S0 INTERACTION
      DATA KAPPA/-0.03992,0.2317/, BETA/1.1771,1.4057/
C     PP 1S0 INTERACTION
C     DATA KAPPA/-0.05398,0.2317/, BETA/1.3224,1.4057/
      HELP=-2./ACOS(-1.)
      I=1
      IF( IS.EQ.1 ) I=2
      PNIJM=HELP*( BETA(I)+KAPPA(I) )**2*2.*BETA(I)
      PNIJM=PNIJM/( (P**2+BETA(I)**2)*(PS**2+BETA(I)**2) )
      RETURN
      END
      SUBROUTINE POTMAT1(P,PS,L,IS,JJ,PNIJM)
      DATA VA1,AMI1,VR1,RMI1/513.968,1.550,1438.720,3.110/
      DATA VA3,AMI3,VR3,RMI3/626.885,1.550,1438.720,3.110/
      IF(JJ) 1,1,2
    1 VA=VA1
      AMI=AMI1
      VR=VR1
      RMI=RMI1
      GO TO 3
    2 VA=VA3
      AMI=AMI3
      VR=VR3
      RMI=RMI3
    3 CONTINUE
      AMI2=AMI**2
      RMI2=RMI**2
      HELPA=2.0*VA/3.1415926
      HELPR=2.0*VR/3.1415926
      X1=(P+PS)**2
      X2=(P-PS)**2
      IF( P.EQ.0. .OR. PS.EQ.0. ) GO TO 20
      PNIJM= -0.25*HELPA*ALOG( (AMI2+X1)/(AMI2+X2) )
      PNIJM=PNIJM+0.25*HELPR*ALOG( (RMI2+X1)/(RMI2+X2) )
      PNIJM=PNIJM/(P*PS)
      GO TO 30
   20 PNIJM=-HELPA/(AMI2+X1)+HELPR/(RMI2+X1)
   30 PNIJM=PNIJM/41.467
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
CJG   CALL D01BCF(0.,-1.,1.,0.,0.,NP1,AP1,XP1,IF)
      CALL GAULEG(XP1,AP1,NP1) 
      DO 1 I=1,NP1
      X=XP1(I)
      A=AP1(I)
      XX=1./P1-(1./P1-2./P2)*X
      XP1(I)=(1.+X) / XX
    1 AP1(I)=(2./P1-2./P2)*A / XX**2
C
      IF(NP2 .NE. 0) THEN
CJG   CALL D01BCF(0.,-1.,1.,0.,0.,NP2,AP2,XP2,IF)
      CALL GAULEG(XP2,AP2,NP2) 
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
      SUBROUTINE GAULEG(XX,WW,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(100),W(100)
      REAL XX(N),WW(N)
C
      PI=DACOS(-1.D0)
C
C      Z=1.D-6
C   20 Z=Z*.1D0
C      IF (DABS(1.D0+Z).EQ.1.D0) GOTO 30
C      GOTO 20
C   30 EPS=Z
C      WRITE(*,*) 'EPS= ',EPS
C
      EPS=1.0D-14
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
C         WRITE(*,*) Z-Z1
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
      SUBROUTINE GNVEC(AK,N,EIG,VEC,VEC2,EPS1,IS,IZ,ITSMX,EIGENO)
C COMMON  -GNVEC   -S-00-000-H.KMD(920412)
      PARAMETER(EPS=1.E-10,MOD=1,IPIC=5)
      DIMENSION AK(N,N),VEC(N),VEC2(N)
      DIMENSION IS(N),IZ(N)
C
      IF(IPIC.GT.N) WRITE(6,*) 'EIGEN VALUE IS NOT GIVEN TRUE.'
C
        DO 2100 I=1,N
          AK(I,I)=AK(I,I)-EIG
 2100   CONTINUE
        CALL INVERT(AK,N,EPS,MOD,KRIT,IS,IZ,N)
        DO 2400 I=1,N
          VEC(I)=1.0
 2400   CONTINUE
C
        ITS=0
   10   CONTINUE
        ITS=ITS+1
        IF (ITS.GT.ITSMX) STOP 'ITS'
        DO 2500 J=1,N
          S=0.0
          DO 2410 I=1,N
           S=S + AK(J,I)*VEC(I)
 2410    CONTINUE
          VEC2(J)=S
 2500   CONTINUE
        EIGENO= EIG + VEC(IPIC)/VEC2(IPIC)
        S=0.0
        DO 2510 I=1,N
          S=S+VEC2(I)**2
 2510   CONTINUE
        S=1.0/SQRT(S)
        SS=0.0
        DO 2520 I=1,N
          VEC2(I)=VEC2(I)*S
          SS=SS+(ABS(VEC2(I))-ABS(VEC(I)) )**2
 2520   CONTINUE
        WRITE(6,*) (VEC2(I)/VEC(I),I=1,3) ,SS
        IF ( SS .GT.EPS1 ) THEN
          DO 2530 I=1,N
            VEC(I)=VEC2(I)
 2530     CONTINUE
          GOTO 10
        END IF
      RETURN
      END
      SUBROUTINE INVERT (A,N,EPS,MOD,KRIT,S,Z,M)
C     =================
      INTEGER S(N),Z(N)
      DIMENSION A(M,M)
C*****M=MATR.-DIM. IM CALLING PROGRAM, N=MATR.-DIM. IN INVERT*****
C*****PRESETTING OF BOOK-KEEPING ARRAYS*****
      DO 101 I=1,N
      Z(I)=-I
  101 S(I)=I
      MS=0
      MZ=0
      KRIT=0
      DO 999 L=1,N
C*****SEARCH FOR PIVOT*****
      IF (MOD) 102,105,102
  102 W=0.0
      DO 104 I=1,N
      IF ( Z(I) .GT. 0) GO TO 104
      DO 103 K=1,N
      IF ( S(K) .LT. 0) GO TO 103
      U=ABS(A(I,K))
      IF (W .GT. U) GO TO 103
      W=U
      MZ=I
      MS=K
  103 CONTINUE
  104 CONTINUE
      GO TO 106
  105 MS=MS+1
      MZ=MZ+1
C*****REPLACING*****
  106 H=A(MZ,MS)
      IF(ABS(H) .GE. EPS) GO TO 107
      H=SIGN(EPS,H)
      KRIT=1
  107 H=1.0/H
      A(MZ,MS)=H
      DO 109 MY=1,N
      IF (MY .EQ. MZ) GO TO 109
      G=-A(MY,MS)*H
      A(MY,MS)=G
      DO 108 NY=1,N
      IF (NY .EQ. MS) GO TO 108
      A(MY,NY)=G*A(MZ,NY)+A(MY,NY)
  108 CONTINUE
  109 CONTINUE
      DO 110 NY=1,N
      IF (NY .EQ. MS) GO TO 110
      A(MZ,NY)=A(MZ,NY)*H
  110 CONTINUE
C*****CHANGE ROWS AND COLOMNS*****
      IF (MOD) 111,999,111
  111 DO 112 I=1,N
      K=S(MS)
      U=A(MZ,I)
      A(MZ,I)=A(K,I)
  112 A(K,I)=U
      DO 113 I=1,N
      K=-Z(MZ)
      U=A(I,MS)
      A(I,MS)=A(I,K)
  113 A(I,K)=U
C*****CHANGE BOOK-KEEPING*****
      K=Z(MZ)
      Z(MZ)=S(MS)
      S(MS)=K
      K=Z(MZ)
      Z(MZ)=Z(K)
      Z(K)=K
      K=-S(MS)
      S(MS)=S(K)
      S(K)=-K
  999 CONTINUE
      RETURN
      END
      SUBROUTINE BALANC(A,N,NP)
      PARAMETER(RADIX=2.0,SQRDX=RADIX**2)
      DIMENSION A(NP,NP)
    1 CONTINUE
      LAST=1
      DO 14 I=1,N
        C=0.0
        R=0.0
        DO 11 J=1,N
          IF (J.NE.I) THEN
            C=C+ABS(A(J,I))
            R=R+ABS(A(I,J))
          END IF
   11   CONTINUE
        IF (C.NE.0.0 .AND. R.NE.0.0 ) THEN
          G=R/RADIX
          F=1.0
          S=C+R
    2     IF (C.LT.G) THEN
            F=F*RADIX
            C=C*SQRDX
            GO TO 2
          END IF
          G=R*RADIX
    3     IF (C.GT.G) THEN
            F=F/RADIX
            C=C/SQRDX
            GO TO 3
          END IF
          IF ( (C+R)/F.LT.0.95*S ) THEN
            LAST=0
            G=1.0/F
            DO 12 J=1,N
              A(I,J)=A(I,J)*G
   12       CONTINUE
            DO 13 J=1,N
              A(J,I)=A(J,I)*F
   13       CONTINUE
          END IF
        END IF
   14 CONTINUE
      IF (LAST.EQ.0) GO TO 1
      RETURN
      END
      SUBROUTINE ELMHES(A,N,NP)
      DIMENSION A(NP,NP)
      DO 17 M=2,N-1
        X=0.0
        I=M
        DO 11 J=M,N
          IF ( ABS( A(J,M-1) ).GT.ABS(X)) THEN
            X=A(J,M-1)
            I=J
          END IF
   11   CONTINUE
        IF(I.NE.M) THEN
          DO 12 J=M-1,N
            Y=A(I,J)
            A(I,J)=A(M,J)
            A(M,J)=Y
   12     CONTINUE
          DO 13 J=1,N
            Y=A(J,I)
            A(J,I)=A(J,M)
            A(J,M)=Y
   13     CONTINUE
        END IF
        IF (X.NE.0.0) THEN
          DO 16 I=M+1,N
            Y=A(I,M-1)
            IF (Y.NE.0.0) THEN
              Y=Y/X
              A(I,M-1)=Y
              DO 14 J=M,N
                A(I,J)=A(I,J)-Y*A(M,J)
   14         CONTINUE
              DO 15 J=1,N
                A(J,M)=A(J,M)+Y*A(J,I)
   15         CONTINUE
            END IF
   16     CONTINUE
        END IF
   17 CONTINUE
      RETURN
      END
      SUBROUTINE HQR(A,N,NP,WR,WI)
      DIMENSION A(NP,NP),WR(NP),WI(NP)
      ANORM=ABS(A(1,1))
      DO 12 I=2,N
        DO 11 J=I-1,N
          ANORM=ANORM + ABS(A(I,J))
   11   CONTINUE
   12 CONTINUE
      NN=N
      T=0.0
    1 IF ( NN.GE.1) THEN
        ITS=0
    2   DO 13 L=NN,2,-1
          S=ABS(A(L-1,L-1)) + ABS(A(L,L))
          IF (S.EQ.0.0) S=ANORM
          IF (ABS(A(L,L-1))+S.EQ.S) GO TO 3
   13   CONTINUE
        L=1
    3   X=A(NN,NN)
        IF (L.EQ.NN) THEN
          WR(NN)=X+T
          WI(NN)=0.0
          NN=NN-1
        ELSE
          Y=A(NN-1,NN-1)
          W=A(NN,NN-1)*A(NN-1,NN)
          IF (L.EQ.NN-1) THEN
            P=0.5*(Y-X)
            Q=P*P+W
            Z=SQRT(ABS(Q))
            X=X+T
            IF (Q.GE.0.0) THEN
              Z=P+SIGN(Z,P)
              WR(NN)=X+Z
              WR(NN-1)=WR(NN)
              IF (Z.NE.0.0) WR(NN)=X-W/Z
              WI(NN)=0.0
              WI(NN-1)=0.0
            ELSE
              WR(NN)=X+P
              WR(NN-1)=WR(NN)
              WI(NN)=Z
              WI(NN-1)=-Z
            END IF
            NN=NN-2
          ELSE
            IF (ITS.EQ.30) PAUSE 'too many iterations'
            IF (ITS.EQ.10 .OR. ITS.EQ.20 ) THEN
              T=T+X
              DO 14 I=1,NN
                A(I,I)=A(I,I)-X
   14         CONTINUE
              S=ABS( A(NN,NN-1) ) + ABS( A(NN-1,NN-2) )
              X=0.75*S
              Y=X
              W=-0.4375*S*S
            END IF
            ITS=ITS+1
            DO 15 M=NN-2,L,-1
              Z=A(M,M)
              R=X-Z
              S=Y-Z
              P=(R*S-W)/A(M+1,M)+A(M,M+1)
              Q=A(M+1,M+1)-(Z+R+S)
              R=A(M+2,M+1)
              S=ABS(P)+ABS(Q)+ABS(R)
              P=P/S
              Q=Q/S
              R=R/S
              IF (M.EQ.L) GO TO 4
              U=ABS(A(M,M-1))*(ABS(Q)+ABS(R))
              V=ABS(P)*(ABS(A(M-1,M-1))+ABS(Z)+ABS(A(M+1,M+1)))
              IF (U+V.EQ.V) GO TO 4
   15       CONTINUE
    4       DO 16 I=M+2,NN
              A(I,I-2)=0.0
              IF (I.NE.M+2) A(I,I-3)=0.0
   16       CONTINUE
            DO 19 K=M,NN-1
              IF (K.NE.M) THEN
                P=A(K,K-1)
                Q=A(K+1,K-1)
                R=0.0
                IF (K.NE.NN-1) R=A(K+2,K-1)
                X=ABS(P)+ABS(Q)+ABS(R)
                IF (X.NE.0.0) THEN
                  P=P/X
                  Q=Q/X
                  R=R/X
                END IF
              END IF
              S=SIGN(SQRT(P*P+Q*Q+R*R),P)
              IF(S.NE.0.0) THEN
                IF (K.EQ.M) THEN
                  IF(L.NE.M) A(K,K-1)=-A(K,K-1)
                ELSE
                  A(K,K-1)=-S*X
                END IF
                P=P+S
                X=P/S
                Y=Q/S
                Z=R/S
                Q=Q/P
                R=R/P
                DO 17 J=K,NN
                  P=A(K,J)+Q*A(K+1,J)
                  IF (K.NE.NN-1) THEN
                    P=P+R*A(K+2,J)
                    A(K+2,J)=A(K+2,J)-P*Z
                  END IF
                  A(K+1,J)=A(K+1,J)-P*Y
                  A(K,J)=A(K,J)-P*X
   17           CONTINUE
                DO 18 I=L,MIN(NN,K+3)
                  P=X*A(I,K)+Y*A(I,K+1)
                  IF (K.NE.NN-1) THEN
                    P=P+Z*A(I,K+2)
                    A(I,K+2)=A(I,K+2)-P*R
                  END IF
                  A(I,K+1)=A(I,K+1)-P*Q
                  A(I,K)=A(I,K)-P
   18           CONTINUE
              END IF
   19       CONTINUE
            GOTO 2
          END IF
        END IF
        GOTO 1
      END IF
      RETURN
      END
