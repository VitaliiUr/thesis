      PROGRAM DEUTERON_PHOTON_CM
C
C##   Deuteron photodisintegration with AV18 in the c.m. frame
C##   Made in April 2008 to compare with chiral calculations
C##   pi- and rho-like MEC consistent with the AV18 force can be included 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      REAL FMM1,MNEUT,MPROT,MM
      PARAMETER (FMM1=197.327D0,
     &           MNEUT=939.5653D0,MPROT=938.2720D0,
     &           MM=2.0D0*MPROT*MNEUT/(MPROT+MNEUT))
      PARAMETER (PBAR=30.0)
      PARAMETER (NP1=10,NP2=35,NP=NP1+NP2+1,NPS=NP1+NP2+1,NPI=NP1+NP2)
      COMMON /GAUSS/ AP(NP),P(NP)
      COMMON /GAUSS2/ GPNEW(NPI),PNEW(NPI)
c     COMMON /GAUSS/ APS,PS          in MEC part
c     COMMON /GAUSS2/ AP,P           in MEC part
c     COMMON /GAUSS2/ GPNEW,PNEW     in T-matrix part
      COMMON /KINEM/ OMEGA,QBIG,ELABNN,P0,THEQ
      DIMENSION AP1(NP1),XP1(NP1),AP2(NP1),XP2(NP1)
      DIMENSION AP3(NP2),XP3(NP2),AP4(NP2),XP4(NP2)
      PARAMETER (NDEUT=100)
      DIMENSION QDEUT(NDEUT),SFWAVE(NDEUT),DFWAVE(NDEUT)
      DIMENSION SPLD(NDEUT)
      COMMON /DEUTINT/ PHID(NP,2)
      double precision M,MD
!     COMMON /DETECT/ MTTWE1,MTTWE2
C     DATA P1/1.0/,P2/9.0/,P3/60.0/
C
      PI=DACOS(-1.0d0)
      AMNMEV=MM    ! MeV
      AMN=MM/FMM1  ! 1/fm
C
      WRITE(6,*)'INPUT  OMEGALAB '
      READ(5,*) OMEGALAB
!     WRITE(6,*)'CHOICE FIRST OUTGOING PARTICLE:1 -proton, -1 - neutron'
!     READ(5,*) MTTWE1 ! tu definiuje, ktora czatka ma wylatywac pierwsza 
!     IF (ABS(MTTWE1).NE.1) THEN
!       STOP
!       WRITE(6,*) ' MTTWE1= ',MTTWE1
!     ENDIF 
!     MTTWE2=-MTTWE1
C
      M=AMNMEV ! MeV
      absb2= 2.225 ! MeV
      MD=2*M-absb2 ! MeV
      OMEGACM= (-2*M*MD + Sqrt(2.)*Sqrt(-4*absb2*M**2*MD + 8*M**3*MD -
     &           2*M**2*MD**2 + 4*M**2*MD*OMEGALAB - M*MD*OMEGALAB**2))
     &            /(2.*M) ! non-relativistically
      QBIGCM=OMEGACM/FMM1 ! 1/fM
C     omegacm= (Sqrt(md)*omegalab)/Sqrt(md + 2*omegalab) ! relativistically
C
      WRITE(6,700) OMEGALAB,OMEGACM
  700 FORMAT(1X,'    DEUTERON PHOTO-DISINTEGRATION '
     2 /5X,'INCOMING LAB PHOTON ENERGY=',F10.4,' (MeV)'
     2 /5X,'CORRESPONDING CM PHOTON ENERGY=',F10.4,' (MeV)')
C
C**** Deuteron wave function is read in
      OPEN(10,FILE='deuwf',STATUS='OLD',FORM='UNFORMATTED')
      READ (10) QDEUT,EBDEUT
      READ (10) SFWAVE
      READ (10) DFWAVE
      CLOSE (10)
      WRITE(6,*) ' # EBDEUT= ',EBDEUT
      WRITE(6,*) ' # QDEUT(I),SFWAVE(I),DFWAVE(I) '
      do i=1,NDEUT
      WRITE(6,*) QDEUT(I),SFWAVE(I),DFWAVE(I)
      enddo
C
      ELABNN=OMEGALAB-ABS(EBDEUT)   ! MeV
      ECMNN=ELABNN-OMEGALAB**2/(4.0*AMNMEV)   ! MeV
      P0MEV=SQRT(ECMNN*AMNMEV) ! MeV
      P0=P0MEV/FMM1  ! 1/fm
      CORR=MD/(MD+OMEGACM) ! strumien 1/|v1-v2|
      CAPTURE=1.5*OMEGACM**2/P0MEV**2
C
C**** More kinematical quantities
      OMEGA=OMEGACM
      QBIG=QBIGCM
      WRITE(6,*) ' # OMEGALAB,QBIGLAB= ',OMEGALAB,OMEGALAB/FMM1
      WRITE(6,*) ' # OMEGACM,QBIGCM= ',OMEGACM,QBIGCM
      WRITE(6,*) ' # ELABNN,ECMNN,P0MEV,P0= ',ELABNN,ECMNN,P0MEV,P0
C
C
C**** Choice of P-points
C     CALL TRNS(NP1,NP2,NP1+NP2,P1,P2,P3,P,AP)
      CALL D01BCF(0,-1.0d0,1.0d0,0.0d0,0.0d0,NP1,AP1,XP1,IFE)
      CALL TRANSFO(NP1,XP1,AP1,0.0d0,P0,XP2,AP2)
      CALL D01BCF(0,-1.0d0,1.0d0,0.0d0,0.0d0,NP2,AP3,XP3,IFE)
      CALL TRANSFO(NP2,XP3,AP3,P0,PBAR,XP4,AP4)
      DO IP=1,NP1
      P(IP)=XP2(IP)
      AP(IP)=AP2(IP)
      ENDDO
      DO IP=1,NP2
      P(IP+NP1)=XP4(IP)
      AP(IP+NP1)=AP4(IP)
      ENDDO
      P(NP)=P0
      AP(NP)=0.0
C
C
C**** Deuteron wave function interpolation to p integral points
      CALL SPLFAK(QDEUT,NDEUT)
C
      DO 30 IP=1,NP
        IF( P(IP).LE.QDEUT(NDEUT) ) THEN
          CALL SPLMOD(QDEUT,NDEUT,P(IP),SPLD)
        ELSE
          DO 20 I=1,NDEUT
   20     SPLD(I)=0.
        ENDIF
        PHID(IP,1)=SKALPR(SPLD,SFWAVE,NDEUT)
        PHID(IP,2)=SKALPR(SPLD,DFWAVE,NDEUT)
   30 CONTINUE
C
C**** Deuteron wave function norm check
      S0=0.
      S2=0.
      DO 150 IP=1,NP-1
      S0= S0 + PHID(IP,1)**2*P(IP)**2*AP(IP)
  150 S2= S2 + PHID(IP,2)**2*P(IP)**2*AP(IP)
      WRITE(6,1000) S0*100.0,S2*100.0,(S0+S2)*100.
 1000 FORMAT(' == Norm check ==  S:',F12.5,'  D:',F12.5,'  T:',F12.5)
C
C
      WRITE(15,*) QBIGCM,NP,PBAR,P,AP
C
      STOP
      END
      SUBROUTINE GAULEG(X,W,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION X(N),W(N)
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
      EPS=1.0D-15
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
      RETURN
      END
      SUBROUTINE D01BCF(ITYPE,A,B,C,D,N,WEIGHT,ABSCIS,IFAIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      integer N
      integer ITYPE,IFAIL  ! now DUMMY
c     real A,B,C,D ! now DUMMY
      DIMENSION ABSCIS(N),WEIGHT(N)
      CALL GAULEG(ABSCIS,WEIGHT,N)
      RETURN
      END
      SUBROUTINE TRANSFO(N,X1,GX1,A,B,X,GX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION X1(N),GX1(N),X(N),GX(N)
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
      SUBROUTINE SPLFAK(X,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER ( NSPL=100 )
C     DIE 2.ABLEITUNGEN AN DEN RAENDERN MUESSEN = 0 SEIN !
C     MIT DEN N STUETZSTELLEN X WIRD DER SPLINE AUFGEBAUT
      DIMENSION X(N)
      DIMENSION HI(NSPL),U(NSPL)
      COMMON /FAKTOR/ FAK1(NSPL,NSPL),FAK2(NSPL,NSPL),FAK3(NSPL,NSPL),
     &                Q(NSPL,NSPL),C(NSPL,NSPL)
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
      IF(J.EQ.I-1) H1=1./(CX*BX)
      IF(J.EQ.I) H2=1./(CX*AX)
      IF(J.EQ.I+1) H3=1./(AX*BX)
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
      IF(L.EQ.J1) FAK1(J1,L)=FAK1(J1,L)-H1
      IF(L.EQ.J1+1) FAK1(J1,L)=FAK1(J1,L)+H1
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
C     =================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (NSPL=100)
C     MITTELS N STUETZSTELLEN X WERDEN DIE FUNKTIONSWERTE SPL
C     AN BELIEBIGEN STELLEN XA INTERPOLIERT
C     IM HP MUSS ZUR ENDGUELTIGEN RECHNUNG NOCH
C     SUMME UEBER J SPL(J)*Y(J) GEBILDET WERDEN
      DIMENSION X(N)
      DIMENSION SPL(N)
      COMMON /FAKTOR/ FAK1(NSPL,NSPL),FAK2(NSPL,NSPL),FAK3(NSPL,NSPL),
     &                QDUM(NSPL,NSPL),CDUM(NSPL,NSPL)
   1  I=0
   2  I=I+1
      IF(I.GT.N) GO TO 3
      IF(XA.GE.X(I)) GO TO 2
   3  I1=I-1
      IF(I1.EQ.0) I1=1
      DX=XA-X(I1)
      DO 10 J=1,N
      SPL(J)=((FAK3(I1,J)*DX+FAK2(I1,J))*DX+FAK1(I1,J))*DX
      IF(J.EQ.I1) SPL(J)=SPL(J)+1.
  10  CONTINUE
      RETURN
      END
      FUNCTION SKALPR(A1,A2,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION A1(N),A2(N)
      SUM=0.
      DO 1 I=1,N
    1 SUM=SUM+A1(I)*A2(I)
      SKALPR=SUM
      RETURN
      END
      DOUBLE COMPLEX FUNCTION ZSKALPR(R,Z,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DOUBLE PRECISION R(N)
      DOUBLE COMPLEX SUM,Z(N)
      SUM=(0.,0.)
      DO 1 I=1,N
    1 SUM=SUM+Z(I)*R(I)
      ZSKALPR=SUM
      RETURN
      END
