      PROGRAM DEUTERON_PHOTON_CM
C
C##   Deuteron photodisintegration with IEE in the c.m. frame
C##   Corrected on February 21, 2015 (kinematical factors)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      double precision FMM1,MN,MP,MM,Mnuc
      PARAMETER (FMM1=197.327D0,
     &           MN=939.5653D0,MP=938.2720D0,MM=2.0D0*MP*MN/(MP+MN))
      PARAMETER (NFAC=50,NDEUT=100,NX=40,NTHMAX=181)
      PARAMETER (PBAR=30.0)
      PARAMETER (NP1=10,NP2=35,NP=NP1+NP2+1,NPS=NP1+NP2+1,NPI=NP1+NP2)
      INTEGER RMAXT
      PARAMETER (JMAX=4,LCHMX2C=2+4*JMAX,RMAXT=4+JMAX)
      COMMON /GAUSS/ AP(NP),P(NP)
      COMMON /GAUSS2/ GPNEW(NPI),PNEW(NPI)
c     COMMON /GAUSS/ APS,PS          in MEC part
c     COMMON /GAUSS2/ AP,P           in MEC part
c     COMMON /GAUSS2/ GPNEW,PNEW     in T-matrix part
      COMMON /PPOW/ PPW(NP,0:5)
      COMMON /QBPOW/ QBPW(0:5)
      COMMON /KINEM/ OMEGA,QBIG,ELABNN,P0,THEQ
      COMMON /XSEC/ PHSFAC(NTHMAX),CORR,CAPTURE
      COMMON /NFFQB/ FP,FN,FPS,FNS
      DIMENSION QDEUT(NDEUT),SFWAVE(NDEUT),DFWAVE(NDEUT)
      DIMENSION SPLD(NDEUT)
      COMMON /DEUTINT/ PHID(NP-1,2)
      DIMENSION XP(NX), AXP(NX)
      DIMENSION PHIDTILD(NX,NP,2)
      DIMENSION PLX(NX,0:RMAXT)
      DIMENSION AP1(NP1),XP1(NP1),AP2(NP1),XP2(NP1)
      DIMENSION AP3(NP2),XP3(NP2),AP4(NP2),XP4(NP2)
      COMMON /GKALL/ GK(NP,0:RMAXT,2)
      COMMON /SGLE/ RHOM(NP,LCHMX2C,-1:1),
     &            CONVPL(NP,LCHMX2C,-1:1),
     &            SPINPL(NP,LCHMX2C,-1:1)
      COMMON /CURRENT/ CHARGE(NP,LCHMX2C,-1:1),
     &                   PLUS(NP,LCHMX2C,-1:1)
      double precision m,md
      COMMON /DETECT/ MTTWE1,MTTWE2
C     DATA P1/1.0/,P2/9.0/,P3/60.0/
C
      character(len=2) FORCE
      integer OSTAT,CUTNUM
      common /pot/ OSTAT,CUTNUM,FORCE
      FORCE='np'
!     OSTAT=4 ! N4LO
!     CUTNUM=2 ! R=0.9 fm
C
      READ(1,*) OSTAT,CUTNUM
C
      PI=DACOS(-1.0d0)
      AMNMEV=MM ! MeV
      AMN=AMNMEV/FMM1 ! 1/fm
      WRITE(6,*) ' (new) HQM= ',FMM1**2/MM
C
c      IRESC=1
      WRITE(6,*)'INPUT: OMEGALAB, IRESC, IMEC'
      READ(5,*) OMEGALAB, IRESC, IMEC
      WRITE(6,*)'CHOICE FIRST OUTGOING PARTICLE:1 -proton, -1 - neutron'
      READ(5,*) MTTWE1 ! tu definiuje, ktora czatka ma wylatywac pierwsza
      IF (ABS(MTTWE1).NE.1) THEN
        STOP
        WRITE(6,*) ' MTTWE1= ',MTTWE1
      ENDIF
      MTTWE2=-MTTWE1
C
      m=AMNMEV ! MeV
      absb2= 2.225 ! MeV
      md=2*m-absb2 ! MeV
      !OMEGA - momentum of photon
      !obtain the formula from kinematics
      OMEGACM= (-2*m*md + Sqrt(2.)*Sqrt(-4*absb2*m**2*md + 8*m**3*md -
     &           2*m**2*md**2 + 4*m**2*md*omegalab - m*md*omegalab**2))
     &            /(2.*m) ! non-relativistically
      QBIGCM=OMEGACM/FMM1 ! 1/fm
C     omegacm= (Sqrt(md)*omegalab)/Sqrt(md + 2*omegalab) ! relativistically
      THEQ= 0.0 ! deg
C
      WRITE(6,700) OMEGALAB,OMEGACM,IRESC,IMEC
  700 FORMAT(1X,'    DEUTERON PHOTO-DISINTEGRATION '
     2 /5X,'INCOMING LAB PHOTON ENERGY=',F10.4,' (MeV)'
     2 /5X,'CORRESPONDING CM PHOTON ENERGY=',F10.4,' (MeV)'
     5 /5X,'IRESC=',I4
     6 /5X,'IMEC=',I4)
C
C**** Deuteron wave function is read in
      OPEN(10,FILE='deuwf',STATUS='OLD',FORM='UNFORMATTED')
      READ(10) QDEUT,EBDEUT
      READ(10) SFWAVE
      READ(10) DFWAVE
      CLOSE (10)
      WRITE(6,*) ' # EBDEUT= ',EBDEUT
      WRITE(6,*) ' # QDEUT(I),SFWAVE(I),DFWAVE(I) '
      do i=1,NDEUT
      WRITE(6,*) QDEUT(I),SFWAVE(I),DFWAVE(I)
      enddo
C
      ! derive these
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
      WRITE(6,*) ' # CORR= ',CORR
      WRITE(6,*) ' # CAPTURE= ',CAPTURE
C
C**** Auxiliary quantities
      CALL AUXIL
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
C**** Filling COMMON/GAUSS2/
      DO IP=1,NP-1
      PNEW(IP)=P(IP)
      GPNEW(IP)=AP(IP)
      ENDDO
      write(6,*) ' p '
      write(6,*) p
C
C**** Powers of P-points
      DO 12 IP=1,NP
   12 PPW(IP,0)=1.
      DO 14 L=1,5
      DO 14 IP=1,NP
   14 PPW(IP,L)=PPW(IP,L-1)*P(IP)
C
C**** Powers of Q
      DO 41 L=0,5
      QBPW(L)=QBIG**L
  41  CONTINUE
C
C**** Deuteron wave function interpolation to p integral points
      !make with analytical function (e.g. polynomial)
      CALL SPLFAK(QDEUT,NDEUT)
C
      DO 30 IP=1,NP-1
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
      S0= S0 + PHID(IP,1)**2*PPW(IP,2)*AP(IP)
  150 S2= S2 + PHID(IP,2)**2*PPW(IP,2)*AP(IP)
      WRITE(6,1000) S0*100.0,S2*100.0,(S0+S2)*100.
 1000 FORMAT(' == Norm check ==  S:',F10.3,'  D:',F10.3,'  T:',F10.3)
C
C**** EM formfactors
c     INFSEL=4
c     CALL NUCLFFJG(QBIG,FP,FN,FPS,FNS,INFSEL)
c     QBIGGEV=QBIG*FMM1*0.001
c     CALL     hmd3(QBIGGEV**2,FP,FPS,FN,FNS) !!! GEp,GMp,GEn,GMn
      CALL     hmd3(0.0d0,FP,FPS,FN,FNS) !!! GEp,GMp,GEn,GMn
c     QQSQ=QBIGMEV**2 - OMEGA**2 ! MeV**2
c     T=0.000001*QQSQ ! GeV**2
c     CALL dipolfit(t,FP,FPS,FN,FNS)
C
C**** Deuteron wave function interpolation to \tilde{k}(p,Q,x) points
      CALL SPLFAK(QDEUT,NDEUT)
      IFEHL=0
      CALL D01BCF(0,-1.0d0,1.0d0,0.0d0,0.0d0,NX,AXP,XP,IFEHL)
      CALL LGDRE(XP,NX,PLX,RMAXT)
C
      QQ14=0.25*QBIG**2
      DO 290 IP=1,NP
      PP2=P(IP)**2
      PQ=P(IP)*QBIG
      DO 291 IX=1,NX
      PQX=PQ*XP(IX)
      PTILD=SQRT(QQ14+PP2+PQX)
        IF( PTILD.LE.QDEUT(NDEUT) ) THEN
          CALL SPLMOD(QDEUT,NDEUT,PTILD,SPLD)
        ELSE
          DO 21 I=1,NDEUT
   21     SPLD(I)=0.
        ENDIF
        PHIDTILD(IX,IP,1)=SKALPR(SPLD,SFWAVE,NDEUT)*AXP(IX)
        PHIDTILD(IX,IP,2)=SKALPR(SPLD,DFWAVE,NDEUT)/PTILD**2*AXP(IX)
  291 CONTINUE
  290 CONTINUE
C
        DO 481 K=0,RMAXT
          DO 481 IP=1,NP
              XINT1=0.
              XINT2=0.
C**** Integral over x - gk is calculated
              DO 471 IX=1,NX
              XINT1=XINT1 + PHIDTILD(IX,IP,1)*PLX(IX,K)
  471         XINT2=XINT2 + PHIDTILD(IX,IP,2)*PLX(IX,K)
          GK(IP,K,1) = XINT1
          GK(IP,K,2) = XINT2
  481 CONTINUE
C
C**** 2N channels are set
      CALL ALFA2CJG
C
C**** Isospin matrix elements are calculated
      CALL ISOSPIN
C
C**** Single-nucleon current matrix elements are calculated
      CALL RHO
      CALL CONV
      CALL SPIN
c     write(6,*) ' RHOM ', RHOM
c     write(6,*) ' CONVPL ', CONVPL
c     write(6,*) ' SPINPL ', SPINPL
C
!     IF (IMEC.EQ.1) THEN
C**** MEC contributions are added
!     CALL ADDMEC
!     ENDIF
      IF (IMEC.EQ.2) THEN
C**** Siegert theorem is used
      READ(150,*) PLUS
      ENDIF

C
C**** Kinematics for final proton-neutron state
      CALL KINEMATICS
C
C**** Plane wave amplitudes are calculated
      CALL PLANEWAVE
C
      IF (IRESC.EQ.1) THEN
C**** Rescattering amplitudes are calculated
      CALL RESCATTER
      ENDIF
C
C**** Cross section is calculated
      CALL CROSSSECTION
C
      STOP
      END
      SUBROUTINE ISOSPIN
C     =============
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DOUBLE PRECISION ISOCHG, ISOMAG
      COMMON /ISOPART/ ISOCHG(0:1), ISOMAG(0:1)
      COMMON /NFFQB/ FP,FN,FPS,FNS
      ISOCHG(0)= 0.5*(FP+FN)
      ISOCHG(1)=-0.5*(FP-FN)
      ISOMAG(0)= 0.5*(FPS+FNS)
      ISOMAG(1)=-0.5*(FPS-FNS)
C
      RETURN
      END
      SUBROUTINE RHO
C     =============
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (NFAC=50)
      DOUBLE PRECISION FMM1,MN,MP,MM
      PARAMETER (FMM1=197.327D0,
     &           MN=939.5653D0,MP=938.2720D0,MM=2.0D0*MP*MN/(MP+MN))
      PARAMETER (NP1=10,NP2=35,NP=NP1+NP2+1)
      INTEGER RMAXT
      PARAMETER (JMAX=4,LCHMX2C=2+4*JMAX,RMAXT=4+JMAX)
      COMMON /GAUSS/ AP(NP),P(NP)
      COMMON /PPOW/ PPW(NP,0:5)
      COMMON /QBPOW/ QBPW(0:5)
      COMMON /KINEM/ OMEGA,QBIG,ELABNN,P0,THEQ
      COMMON /NFFQB/ FP,FN,FPS,FNS
      DOUBLE PRECISION ISOCHG, ISOMAG
      COMMON /ISOPART/ ISOCHG(0:1), ISOMAG(0:1)
      COMMON /GKALL/ GK(NP,0:RMAXT,2)
      COMMON /FACUL/ F(0:NFAC),SF(0:NFAC)
      COMMON /HILF/ D(0:NFAC),SD(0:NFAC)
      COMMON /FKLT/ FK(0:NFAC),WF(0:NFAC)
      COMMON /FAKLOG/ FL(0:NFAC)
      LOGICAL COUPLC
      COMMON /PWTB2C/ INDX2C(LCHMX2C),NINDX2C,COUPLC(LCHMX2C)
      COMMON /SGLE/ RHOM(NP,LCHMX2C,-1:1),
     &            CONVPL(NP,LCHMX2C,-1:1),
     &            SPINPL(NP,LCHMX2C,-1:1)
      INTEGER FMIN,FMAX,FF,FF2
      INTEGER GMIN,GMAX,GG,GG2
      INTEGER RMIN,RMAX,RR,RR2
      INTEGER WMIN,WMAX,WW,WW2,WMIN1,WMAX1,WMIN2,WMAX2
C
      PI=DACOS(-1.0d0)
      AMN=MM/FMM1! 1/fm
      AMNINV=1/AMN ! fm
C
      DO MD=-1,1
      DO IALS=1,LCHMX2C
      DO IP=1,NP
      RHOM(IP,IALS,MD)=0.0
      ENDDO
      ENDDO
      ENDDO
C
      DO 1 IALS=1,NINDX2C
      CALL UNPACK1(4,LS,ISS,JS,ITS,INDX2C(IALS),IZ,IZ,IZ,IZ)
      IF (ISS.NE.1) GOTO 1
      LS2=2*LS
      ISS2=2*ISS
      JS2=2*JS
      WMIN1=IABS(JS-1)
      WMAX1=     JS+1
      S1=ISOCHG(ITS)*SD(JS)
C
      DO 2 L=0,2,2
      LTWO=2*L
      LI=L/2+1
      WMIN2=MAX(WMIN1, IABS(LS-L))
      WMAX2=MIN(WMAX1,      LS+L )
      S2=S1*SD(L)*SF(LTWO+1)
C
      DO 3 L1=0,L
      L2=L-L1
      L12=2*L1
      L22=2*L2
      RMIN=IABS(LS-L1)
      RMAX=     LS+L1
      S3=S2*0.5**L2*QBPW(L2)/(SF(L12)*SF(L22))
C
      DO 4 RR=RMIN,RMAX,2
      RR2=2*RR
      WMIN=MAX(WMIN2, IABS(L2-RR))
      WMAX=MIN(WMAX2,      L2+RR )
      IF (MOD(WMIN+L2+RR,2).NE.0) WMIN=WMIN+1
      S4=S3*D(RR)*(-1.)**RR*CG0003J(RR,L1,LS)
C
      DO 5 WW=WMIN,WMAX,2
      WW2=2*WW
      S5=S4*SD(WW)*CG0003J(RR,L2,WW)
     &      *C6J(WW2,LS2,LTWO,L12,L22,RR2)
     &      *C6J(WW2,LS2,LTWO,2,2,JS2)
C
      DO 6 MD=-1,1
      MD2=2*MD
      S6=S5*CLEBSCH(WW2,JS2,2,0,MD2,MD2)
C
      DO 7 IP=1,NP
      RHOM(IP,IALS,MD)=RHOM(IP,IALS,MD)
     & + S6*PPW(IP,L1)*GK(IP,RR,LI)
C
    7 CONTINUE
    6 CONTINUE
    5 CONTINUE
    4 CONTINUE
    3 CONTINUE
    2 CONTINUE
    1 CONTINUE
C
      RETURN
      END
      SUBROUTINE CONV
C     =============
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (NFAC=50)
      DOUBLE PRECISION FMM1,MN,MP,MM
      PARAMETER (FMM1=197.327D0,
     &           MN=939.5653D0,MP=938.2720D0,MM=2.0D0*MP*MN/(MP+MN))
      PARAMETER (NP1=10,NP2=35,NP=NP1+NP2+1)
      INTEGER RMAXT
      PARAMETER (JMAX=4,LCHMX2C=2+4*JMAX,RMAXT=4+JMAX)
      COMMON /GAUSS/ AP(NP),P(NP)
      COMMON /PPOW/ PPW(NP,0:5)
      COMMON /QBPOW/ QBPW(0:5)
      COMMON /KINEM/ OMEGA,QBIG,ELABNN,P0,THEQ
      COMMON /NFFQB/ FP,FN,FPS,FNS
      DOUBLE PRECISION ISOCHG, ISOMAG
      COMMON /ISOPART/ ISOCHG(0:1), ISOMAG(0:1)
      COMMON /GKALL/ GK(NP,0:RMAXT,2)
      COMMON /FACUL/ F(0:NFAC),SF(0:NFAC)
      COMMON /HILF/ D(0:NFAC),SD(0:NFAC)
      COMMON /FKLT/ FK(0:NFAC),WF(0:NFAC)
      COMMON /FAKLOG/ FL(0:NFAC)
      LOGICAL COUPLC
      COMMON /PWTB2C/ INDX2C(LCHMX2C),NINDX2C,COUPLC(LCHMX2C)
      COMMON /SGLE/ RHOM(NP,LCHMX2C,-1:1),
     &            CONVPL(NP,LCHMX2C,-1:1),
     &            SPINPL(NP,LCHMX2C,-1:1)
      INTEGER FMIN,FMAX,FF,FF2,FMIN1,FMAX1
      INTEGER GMIN,GMAX,GG,GG2
      INTEGER RMIN,RMAX,RR,RR2
      INTEGER WMIN,WMAX,WW,WW2,WMIN1,WMAX1,WMIN2,WMAX2
C
      PI=DACOS(-1.0d0)
      AMN=MM/FMM1 ! 1/fm
      AMNINV=1/AMN ! fm
C
      DO MD=-1,1
      DO IALS=1,LCHMX2C
      DO IP=1,NP
      CONVPL(IP,IALS,MD)=0.0
      ENDDO
      ENDDO
      ENDDO
C
      ITAU=1
      ITAU2=2
C
      DO 1 IALS=1,NINDX2C
      CALL UNPACK1(4,LS,ISS,JS,ITS,INDX2C(IALS),IZ,IZ,IZ,IZ)
      IF (ISS.NE.1) GOTO 1
      LS2=2*LS
      ISS2=2*ISS
      JS2=2*JS
      WMIN1=IABS(JS-1)
      WMAX1=     JS+1
      GMIN=IABS(LS-1)
      GMAX=     LS+1
      S1=ISOCHG(ITS)*AMNINV*SD(JS)*SD(LS)*(-1.)**(1+LS)
C
      DO 2 L=0,2,2
      LTWO=2*L
      LI=L/2+1
      WMIN2=MAX(WMIN1, IABS(LS-L))
      WMAX2=MIN(WMAX1,      LS+L )
      S2=S1*SD(L)*SF(LTWO+1)
C
      DO 3 L1=0,L
      L2=L-L1
      L12=2*L1
      L22=2*L2
      S3=S2*0.5**L2*QBPW(L2)/(SF(L12)*SF(L22))
C
      DO 35 GG=GMIN,GMAX,2
      GG2=2*GG
      RMIN=IABS(GG-L1)
      RMAX=     GG+L1
      FMIN1=IABS(GG-L)
      FMAX1=     GG+L
      S35=S3*CG0003J(LS,1,GG)
C
      DO 4 RR=RMIN,RMAX,2
      RR2=2*RR
      FMIN=MAX(FMIN1, IABS(L2-RR))
      FMAX=MIN(FMAX1,      L2+RR )
      IF (MOD(FMIN+L2+RR,2).NE.0) FMIN=FMIN+1
      S4=S35*D(RR)*(-1.)**RR*CG0003J(RR,L1,GG)
C
      DO 45 FF=FMIN,FMAX,2
      FF2=2*FF
      WMIN=MAX(WMIN2, IABS(FF-1))
      WMAX=MIN(WMAX2,      FF+1 )
      S45=S4*SD(FF)*(-1.)**FF*CG0003J(RR,L2,FF)
     &      *C6J(FF2,GG2,LTWO,L12,L22,RR2)
C
      DO 5 WW=WMIN,WMAX
      WW2=2*WW
      S5=S45*SD(WW)*C6J(LS2,2,GG2,FF2,LTWO,WW2)
     &      *C6J(WW2,LS2,LTWO,2,2,JS2)
     &      *CLEBSCH(2,FF2,WW2,-ITAU2,0,-ITAU2)
C
      DO 6 MD=-1,1
      MD2=2*MD
      S6=S5*CLEBSCH(WW2,JS2,2,-ITAU2,MD2+ITAU2,MD2)
C
      DO 7 IP=1,NP
      CONVPL(IP,IALS,MD)=CONVPL(IP,IALS,MD)
     & + S6*PPW(IP,L1+1)*GK(IP,RR,LI)
C
    7 CONTINUE
    6 CONTINUE
    5 CONTINUE
   45 CONTINUE
    4 CONTINUE
   35 CONTINUE
    3 CONTINUE
    2 CONTINUE
    1 CONTINUE
C
      RETURN
      END
      SUBROUTINE SPIN
C     =============
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (NFAC=50)
      DOUBLE PRECISION FMM1,MN,MP,MM
      PARAMETER (FMM1=197.327D0,
     &           MN=939.5653D0,MP=938.2720D0,MM=2.0D0*MP*MN/(MP+MN))
      PARAMETER (NP1=10,NP2=35,NP=NP1+NP2+1)
      INTEGER RMAXT
      PARAMETER (JMAX=4,LCHMX2C=2+4*JMAX,RMAXT=4+JMAX)
      COMMON /GAUSS/ AP(NP),P(NP)
      COMMON /PPOW/ PPW(NP,0:5)
      COMMON /QBPOW/ QBPW(0:5)
      COMMON /KINEM/ OMEGA,QBIG,ELABNN,P0,THEQ
      COMMON /NFFQB/ FP,FN,FPS,FNS
      DOUBLE PRECISION ISOCHG, ISOMAG
      COMMON /ISOPART/ ISOCHG(0:1), ISOMAG(0:1)
      COMMON /GKALL/ GK(NP,0:RMAXT,2)
      COMMON /FACUL/ F(0:NFAC),SF(0:NFAC)
      COMMON /HILF/ D(0:NFAC),SD(0:NFAC)
      COMMON /FKLT/ FK(0:NFAC),WF(0:NFAC)
      COMMON /FAKLOG/ FL(0:NFAC)
      LOGICAL COUPLC
      COMMON /PWTB2C/ INDX2C(LCHMX2C),NINDX2C,COUPLC(LCHMX2C)
      COMMON /SGLE/ RHOM(NP,LCHMX2C,-1:1),
     &            CONVPL(NP,LCHMX2C,-1:1),
     &            SPINPL(NP,LCHMX2C,-1:1)
      COMMON /CURRENT/ CHARGE(NP,LCHMX2C,-1:1),
     &                   PLUS(NP,LCHMX2C,-1:1)
      INTEGER FMIN,FMAX,FF,FF2,FMIN1,FMAX1
      INTEGER RMIN,RMAX,RR,RR2
      INTEGER WMIN,WMAX,WW,WW2,WMIN1,WMAX1
C
      PI=DACOS(-1.0d0)
      AMN=MM/FMM1 ! 1/fm
      AMNINV=1/AMN ! fm
C
      DO MD=-1,1
      DO IALS=1,LCHMX2C
      DO IP=1,NP
      SPINPL(IP,IALS,MD)=0.0
      ENDDO
      ENDDO
      ENDDO
C
      ITAU=1
      ITAU2=2
C
      DO 1 IALS=1,NINDX2C
      CALL UNPACK1(4,LS,ISS,JS,ITS,INDX2C(IALS),IZ,IZ,IZ,IZ)
      LS2=2*LS
      ISS2=2*ISS
      JS2=2*JS
      WMIN1=MAX(IABS(JS-1),IABS(LS-1))
      WMAX1=MIN(     JS+1,      LS+1 )
      S1=QBIG*AMNINV*ITAU*ISOMAG(ITS)*0.5*SQRT(18.)*SD(ISS)
     &  *(-1.)**(1+LS+ISS)*C6J(2,1,1,1,2,ISS2)
C
      DO 2 L=0,2,2
      LTWO=2*L
      LI=L/2+1
      FMIN1=IABS(LS-L)
      FMAX1=     LS+L
      S2=S1*SD(L)*SF(LTWO+1)
C
      DO 3 L1=0,L
      L2=L-L1
      L12=2*L1
      L22=2*L2
      RMIN=IABS(LS-L1)
      RMAX=     LS+L1
      S3=S2*0.5**L2*QBPW(L2)/(SF(L12)*SF(L22))
C
      DO 4 RR=RMIN,RMAX,2
      RR2=2*RR
      FMIN=MAX(FMIN1, IABS(L2-RR))
      FMAX=MIN(FMAX1,      L2+RR )
      IF (MOD(FMIN+L2+RR,2).NE.0) FMIN=FMIN+1
      S4=S3*D(RR)*(-1.)**RR*CG0003J(RR,L1,LS)
C
      DO 45 FF=FMIN,FMAX,2
      FF2=2*FF
      WMIN=MAX(WMIN1, IABS(FF-1))
      WMAX=MIN(WMAX1,      FF+1 )
      S45=S4*SD(FF)*CG0003J(RR,L2,FF)
     &      *C6J(FF2,LS2,LTWO,L12,L22,RR2)
C
      DO 5 WW=WMIN,WMAX
      WW2=2*WW
      S5=S45*D(WW)*(-1.)**WW*C6J(FF2,LS2,LTWO,2,2,WW2)
     &      *C6J(2,2,ISS2,LS2,JS2,WW2)
C
      DO 6 MD=-1,1
      MD2=2*MD
      S6=S5*CLEBSCH(FF2,WW2,2,0,MD2,MD2)
     &     *CLEBSCH(2,WW2,JS2,ITAU2,MD2,MD2+ITAU2)
C
      DO 7 IP=1,NP
      SPINPL(IP,IALS,MD)=SPINPL(IP,IALS,MD)
     & + S6*PPW(IP,L1)*GK(IP,RR,LI)
C
    7 CONTINUE
    6 CONTINUE
    5 CONTINUE
   45 CONTINUE
    4 CONTINUE
    3 CONTINUE
    2 CONTINUE
    1 CONTINUE
C
C
      DO MD=-1,1
      DO IALS=1,LCHMX2C
      DO IP=1,NP
      CHARGE(IP,IALS,MD)=RHOM(IP,IALS,MD)
      PLUS(IP,IALS,MD)=CONVPL(IP,IALS,MD)+SPINPL(IP,IALS,MD)
      ENDDO
      ENDDO
      ENDDO
C
C
      RETURN
      END
      SUBROUTINE ALFA2CJG
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (JMAX=4,LCHMX2C=2+4*JMAX)
      LOGICAL COUPL
      CHARACTER*4 TEXT
      COMMON /PWTB2C/ INDEX2(LCHMX2C),NINDEX2,COUPL(LCHMX2C)
C
      NINDEX2=0
      WRITE(6,100)
  100 FORMAT(//,'    Channels(C JG)  used for two-body interaction'/
     &5X,'L',4X,'S',4X,'J',4X,'T',4X,'N',4X,'INDEX(N)')
      DO 10 J=0,JMAX
       DO 11 IS=0,1
        DO 11 L=ABS(J-IS),J+IS
         NINDEX2=NINDEX2+1
         IT=((-1)**(L+IS)+1)/2
         CALL PACK1(4,L,IS,J,IT,INDEX2(NINDEX2),IZ,IZ,IZ,IZ)
         CALL UNPACK1(4,L1,IS1,J1,IT1,INDEX2(NINDEX2),IZ,IZ,IZ,IZ)
         IF (  L.NE.J .AND. IS.EQ.1 .AND. J.GT.0  ) THEN
       TEXT='COUP'
         COUPL(NINDEX2)=.TRUE.
         WRITE(6,101) L1,IS1,J1,IT1,NINDEX2,INDEX2(NINDEX2),TEXT
  101    FORMAT(1X,5I5,I11,1X,A4)
         ELSE
       TEXT='UNCO'
         COUPL(NINDEX2)=.FALSE.
         WRITE(6,101) L1,IS1,J1,IT1,NINDEX2,INDEX2(NINDEX2),TEXT
         END IF
   11    CONTINUE
CC         WRITE(6,*)
   10   CONTINUE
      RETURN
      END
      SUBROUTINE TRNS(NP1,NP2,NP,P1,P2,P3,XP,AP)
C     ===============
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION XP1(96),AP1(96),XP2(96),AP2(96)
      DIMENSION XP(NP),AP(NP)
      CALL D01BCF(0,-1.0d0,1.0d0,0.0d0,0.0d0,NP1,AP1,XP1,IF)
      DO 1 I=1,NP1
      X=XP1(I)
      A=AP1(I)
      XX=1./P1-(1./P1-2./P2)*X
      XP1(I)=(1.+X) / XX
    1 AP1(I)=(2./P1-2./P2)*A / XX**2
C
      IF(NP2 .NE. 0) THEN
      CALL D01BCF(0,-1.0d0,1.0d0,0.0d0,0.0d0,NP2,AP2,XP2,IF)
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
      SUBROUTINE TRNS1(NP1,NP2,NP,PMIN,P1,P2,P3,XP,AP)
C     =============
C
C     TRNS BELEGT DIE FELDER XP UND AP MIT TRANSFORMIERTEN
C     GAUSS-LEGENDRE-PUNKTEN UND GEWICHTEN
C
C     NP1 PUNKTE WERDEN UEBER DIE HYPERBOLISCHE TRANSFORMATION
C
C     X --> (A+B*X) / (1.+C*X)
C     AA=P1 ; BB=( P1*(PMIN+P2)-2.*PMIN*P2 )/( P2-PMIN ) ;
C     CC=( 2.*P1-PMIN-P2 )/( P2-PMIN )
C
C     AUF DAS INTERVALL (PMIN;P2) ABGEBILDET, WOBEI
C     NP1/2 PUNKTE IN (PMIN;P1) UND
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION XP1(96),AP1(96),XP2(96),AP2(96)
      DIMENSION XP(NP),AP(NP)
C
      AA=P1
      BB=( P1*(PMIN+P2)-2.*PMIN*P2 )/( P2-PMIN )
      CC=( 2.*P1-PMIN-P2 )/( P2-PMIN )
      IFE=0
      CALL D01BCF(0,-1.0d0,1.0d0,0.0d0,0.0d0,NP1,AP1,XP1,IFE)
CC    CALL GAULEG(XP1,AP1,NP1)
      DO 1 I=1,NP1
      X=XP1(I)
      A=AP1(I)
      XX1=AA+BB*X
      XX2=1.+CC*X
      XP1(I)=XX1/XX2
    1 AP1(I)=( BB*XX2-XX1*CC )*A/(XX2*XX2)
C
      IF(NP2 .NE. 0) THEN
      IFE=0
      CALL D01BCF(0,-1.0d0,1.0d0,0.0d0,0.0d0,NP2,AP2,XP2,IFE)
CC    CALL GAULEG(XP2,AP2,NP2)
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
      SUBROUTINE TRNSX(NX,P,XP,AP)
C     =============
C
C     TRNS BELEGT DIE FELDER XP UND AP MIT TRANSFORMIERTEN
C     GAUSS-LEGENDRE-PUNKTEN UND GEWICHTEN
C
C     NX PUNKTE WERDEN UEBER DIE HYPERBOLISCHE TRANSFORMATION
C
C     X --> (P+X) / (1.+P*X)
C
C     AUF DAS INTERVALL (-1.,1.) ABGEBILDET, WOBEI
C     NX/2 PUNKTE IN (-1.;P) UND
C     NX/2 PUNKTE IN (P;1.) LIEGEN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION XP1(96),AP1(96)
C     DIMENSION XP2(96),AP2(96)
      DIMENSION XP(NX),AP(NX)
C
      IFE=0
      CALL D01BCF(0,-1.0d0,1.0d0,0.0d0,0.0d0,NX,AP1,XP1,IFE)
CC    CALL GAULEG(XP1,AP1,NX)
      DO 1 I=1,NX
      X=XP1(I)
      A=AP1(I)
      XX=1.0+P*X
      XP1(I)=(P+X) / XX
    1 AP1(I)=(1.0-P**2 )*A / XX**2
C
      DO 3 I=1,NX
      XP(I)=XP1(I)
    3 AP(I)=AP1(I)
C
      RETURN
      END
      SUBROUTINE PACK1(N,I1,I2,I3,I4,I5,I6,I7,I8,I9)
C     ===============
C**** EACH OF THE MAXIMAL 8 PACKED (POSITIVE) QUANTUM NUMBERS
C**** HAS TO BE GREATER EQUAL 0 AND LESS THAN 2**MOD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
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
C     =================
C**** EACH OF THE MAXIMAL 8 PACKED (POSITIVE) QUANTUM NUMBERS
C**** HAS TO BE GREATER EQUAL 0 AND LESS THAN 2**MOD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
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
      SUBROUTINE LGDRE (X,NX,PL,KM)
C*******************************************************************************
C
C  LEGENDRE-POLYNOME P (K,X) FUER X(1)...X(NX) UND 0 <= K <= KM
C
C*******************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION X(NX),PL(NX,0:KM)
      DO 10 I=1,NX
       PL(I,0)=1.
   10  PL(I,1)=X(I)
      DO 20 K=1,KM-1
       DO 20 I=1,NX
   20 PL(I,K+1)=(2.*K+1)*X(I)/(K+1.)*PL(I,K)-K/(K+1.)*PL(I,K-1)
      RETURN
      END
      SUBROUTINE DELMIN (J1,J2,J3,MX,MNPER,DELTA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
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
      FUNCTION CG0003J(IA,IB,IC)
C**** Clebsch-Gordan coefficient  <a 0 b 0 | c 0 >
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (NFAC=50)
      COMMON /FACUL/ F(0:NFAC),SF(0:NFAC)
      COMMON/ HILF / D(0:NFAC),SD(0:NFAC)
      CG0003J=0.
      IF ( MOD(IA+IB+IC,2).EQ.0  .AND.
     &     IC.GE.IABS(IA-IB) .AND. IC.LE.(IA+IB) .AND.
     &     IA.GE.0 .AND. IB.GE.0 .AND. IC.GE.0 ) THEN
        IP=(IA+IB+IC)/2
        CG0003J = (-1.)**IP *
     &     SF(IA+IB-IC)*SF(IB+IC-IA)*SF(IC+IA-IB)/SF(IA+IB+IC+1) *
     &     F(IP) / (F(IP-IA)*F(IP-IB)*F(IP-IC)) *
     &     (-1.)**(-IA+IB) * SD(IC)
      ENDIF
C
      RETURN
      END
      FUNCTION CG000 (I,J,K)
C*******************************************************************************
C
C  CLEBSCH-GORDON   < I J 0 0 / K 0 >
C
C  DREIECKSUNGLEICHUNG UND PHASE VON (I+J+K) WERDEN NICHT GETESTET.
C
C*******************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (NFAC=50)
      COMMON /FKLT/ F(0:NFAC),WF(0:NFAC)
      IP=(I+J+K)/2
      CG000=(-1.)**(IP+I-J)*WF(I+J-K)*WF(J+K-I)*WF(K+I-J)/WF(I+J+K+1)
     &      *F(IP)/(F(IP-I)*F(IP-J)*F(IP-K))*SQRT(2.*K+1)
      RETURN
      END
      FUNCTION RACAH (IA,IB,IE,ID,IC,IF)
C     ==============
CCC  ALL  L*  MUST BE TWICE AS LARGE AS TRUE ARGUMENTS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
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
      FUNCTION CLEBSCH(IA2,IB2,IC2,ID2,IE2,IF2)
C     ================
C      THIS IS A CORRECTED ROUTINE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (NFAC=50)
      COMMON/FACUL/F(0:NFAC),SF(0:NFAC)
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
      CLEBSCH=S*DEL*(-1.)** ((IA2-IB2+IF2)/2)*SQRT(FLOAT(IC2+1))
      RETURN
      END
      FUNCTION COEF9J (J1,J2,J5,J3,J4,J6,J7,J8,J9)
C     ===============
CCCCCC
CCCCCC          ALL L9 MUST BE TWICE AS LARGE AS TRUE ARGUMENTS
CCCCCC
CCCC  TAMURAS NUMBERING CONVENTION IS USED HERE.
CCCC                       1  2  5
CCCC                       3  4  6
CCCC                       7  8  9
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
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
      FUNCTION C6J (I,J,K,L,M,N)
C***********************************************************************
C
C  6-J-SYMBOL     ( I/2  J/2  K/2 )
C                 ( L/2  M/2  N/2 )
C
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
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
      FUNCTION C6J2 (I,J,K,L,M)
C*******************************************************************************
C
C  6-J-SYMBOL   ( I/2  J/2  1/2 )
C               ( K/2  L/2  M/2 )
C
C*******************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER G
      JK=MIN(I,J)
      JG=MIN(K,L)
      IF ((J-I)*(L-K) .GE. 0) THEN
       G=M-1
       F=(1+(G+JK-JG)/2)*(1+(G-JK+JG)/2)
      ELSE
       G=M
       F=(1+(JK+JG-G)/2)*(2+(JK+JG+G)/2)
      ENDIF
      IF (ABS(JK-JG) .LE. G .AND. G .LE. JK+JG) THEN
       C6J2=(-1.)**(1+(JK+JG+G)/2)*SQRT(F/((JK+1)*(JK+2)*(JG+1)
     &        *(JG+2)))
      ELSE
       C6J2=0.
      ENDIF
      RETURN
      END
      FUNCTION C9J2 (I,J,K,L,M,N,O,P)
C*******************************************************************************
C
C               ( I/2  J/2  K/2 )
C  9-J-SYMBOL   ( L/2  1/2  M/2 )
C               ( N/2  O/2  P/2 )
C
C*******************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER O,P,G,GMIN,GMAX
      C9J2=0.
      GMIN=MAX(ABS(L-P),ABS(K-1),ABS(O-I))
      GMAX=MIN(L+P,K+1,O+I)
      DO 10 G=GMIN,GMAX,2
   10  C9J2=C9J2+(-1.)**G *(G+1)*C6J2(L,M,K,G,P)*C6J2(K,G,O,J,I)
     &      *C6J(N,O,P,G,L,I)
      C9J2=(-1.)**((I+J+K+L+1+M+N+O+P)/2)*C9J2
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
      SUBROUTINE AUXIL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (NFAC=50)
      COMMON /HILF/ D(0:NFAC),SD(0:NFAC)
      COMMON /FAKLOG/ FL(0:NFAC)
      COMMON /FACUL/ F(0:NFAC),SF(0:NFAC)
      COMMON /FKLT/ FK(0:NFAC),WF(0:NFAC)
C
C**** (2n+1) and sqrt(2n+1) (/HILF/)
      DO 1 I=0,NFAC
      D(I)=2.*FLOAT(I)+1.
    1 SD(I)=SQRT(D(I))
C
C**** ln(n!) (/FAKLOG/)
      FL(0)=0.0
      DO 2 I=1,NFAC
    2 FL(I)=ALOG(FLOAT(I))+FL(I-1)
C
C**** n! and sqrt(n!) (/FACUL/)
      F(0)=1.
      DO 3 I=1,NFAC
    3 F(I)=F(I-1)*FLOAT(I)
      DO 4 I=0,NFAC
    4 SF(I)=SQRT(F(I))
C
C**** n! and sqrt(n!) (/FKLT/)
      FK(0)=1.
      DO 5 I=1,NFAC
    5 FK(I)=I*FK(I-1)
      DO 6 I=0,NFAC
    6 WF(I)=SQRT(FK(I))
C
      RETURN
      END
      SUBROUTINE KINEMATICS
C     =============
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DOUBLE PRECISION FMM1,MN,MP,MM
      PARAMETER (FMM1=197.327D0,
     &           MN=939.5653D0,MP=938.2720D0,MM=2.0D0*MP*MN/(MP+MN))
      PARAMETER (NTHMAX=181)
      COMMON /KINEM/ OMEGA,QBIG,ELABNN,P0,THEQ
      COMMON /KINEM12/ PM1(NTHMAX),THM1(NTHMAX),PM2(NTHMAX),
     &                 THM2(NTHMAX),THPLM1(NTHMAX)
      COMMON /XSEC/ PHSFAC(NTHMAX),CORR,CAPTURE
      COMMON /WINKEL/ THP(NTHMAX),PHP(NTHMAX)  ! rad
      DIMENSION QQQ(3)
C
      PI=DACOS(-1.0d0)
      FF=PI/180.0
      AMNMEV=MM    ! MeV
      AMN=MM/FMM1  ! 1/fm
      QBIGMEV=QBIG*FMM1  ! MeV
      P0MEV=P0*FMM1
C
C**** Azimuthal CM angle of particle 1 in the system where Q || z
      PH1=0.0
      PH1R=PH1*FF
      CSPH1=COS(PH1R)
      SNPH1=SIN(PH1R)
C**** Loop over the proton CM angle with respect to the photon beam
      WRITE(6,*) ' TH1,P1,P2,TH2,THQQQ/FF,PH1R/FF '
      DO ITH=1,NTHMAX
      TH1=FLOAT(ITH-1) ! deg
      TH1R=TH1*FF ! rad
      THPLM1(ITH)=TH1R ! rad
      CSTH1=COS(TH1R)
      SNTH1=SIN(TH1R)

      P1= P0MEV   ! MeV
      PM1(ITH)=P1 ! MeV
      P2= P0MEV   ! MeV
      PM2(ITH)=P2 ! MeV
      TH2R= PI-TH1R ! rad
      TH2=TH2R/FF   ! deg
      THM1(ITH)=TH1R
      THM2(ITH)=TH2R
C**** p0=0.5*(p1-p2)=p1, bo jestesmy w CM
      THP(ITH)=TH1R
      PHP(ITH)=PH1R
      WRITE(6,100) TH1,P1,P2,TH2,THP(ITH)/FF,PH1
  100 FORMAT(7F12.3)
      P1FM=P1/FMM1
c     PHSFAC(ITH)=AMN*P1FM**2/ABS(2.*P1FM-QBIG*CSTH1)
      PHSFAC(ITH)=0.5*AMN*P1FM
      ENDDO
C
      RETURN
      END
      SUBROUTINE CARSPH(X,R,THTA,PHI)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION X(3)
C
C*****TRANSITION FROM CARTESIAN TO SPHERICAL COORDINATES
C
      PI=DACOS(-1.0d0)
      R1=X(1)**2+X(2)**2
      R=SQRT(R1+X(3)**2)
      R1=SQRT(R1)
      IF (R.LT.1.E-16) THEN
      THTA=0.
      PHI=0.
      ELSE
      IF (R1.EQ.0.0) THEN
      PHI=0.0
      THTA=0.0
      IF (X(3).LT.0.0) THTA=PI
      ELSE
      THTA=ACOS(X(3)/R)
      PHI=ATAN2(X(2),X(1))
      ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE YKQFN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (JMAX=4,LMAX=JMAX+1,NFAC=50,NTHMAX=181)
      DOUBLE COMPLEX YLMP
      COMMON /WINKEL/ TH(NTHMAX),PHI(NTHMAX)
      COMMON /YLM/ YLMP(NTHMAX,0:LMAX,-LMAX:LMAX)
      COMMON /PL/ X(NTHMAX),H(NTHMAX),P(NTHMAX,0:LMAX,0:LMAX)
      COMMON /FACUL/ F(0:NFAC),SF(0:NFAC)
C
      PI=DACOS(-1.0d0)
      DO 5 IW=1,NTHMAX
    5 X(IW)=COS(TH(IW))
      CALL LEGDRE
      DO 10 M=0,LMAX
        DO 10 L=M,LMAX
            DO 11 IW=1,NTHMAX
   11         YLMP(IW,L,M)=(-1)**M*SQRT((2*L+1)*F(L-M)/(4*PI*
     .        F(L+M)))*P(IW,L,M)*EXP(CMPLX(0.,M*PHI(IW)))
            IF (M.NE.0) THEN
              DO 12 IW=1,NTHMAX
   12           YLMP(IW,L,-M)=(-1)**M*CONJG(YLMP(IW,L,M))
            ENDIF
   10 CONTINUE
      END
      SUBROUTINE LEGDRE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (JMAX=4,LMAX=JMAX+1,NTHMAX=181)
      COMMON /PL/ X(NTHMAX),H(NTHMAX),P(NTHMAX,0:LMAX,0:LMAX)
      DO 10 IW=1,NTHMAX
        P(IW,0,0)=1.
   10   H(IW)=SQRT(1-X(IW)*X(IW))
      DO 20 M=1,LMAX
        DO 20 IW=1,NTHMAX
   20     P(IW,M,M)=(2*M-1.)*2.*M*0.5*H(IW)/M*P(IW,M-1,M-1)
      DO 30 M=0,LMAX-1
        DO 40 IW=1,NTHMAX
   40     P(IW,M+1,M)=(2*M+1)*X(IW)*P(IW,M,M)
        DO 30 L=M+1,LMAX-1
          DO 30 IW=1,NTHMAX
   30       P(IW,L+1,M)=((2*L+1)*X(IW)*P(IW,L,M)-(L+M)*P(IW,L-1,M))
     .                  /(L-M+1.)
      END
      SUBROUTINE PLANEWAVE
C     =============
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DOUBLE PRECISION FMM1,MNEUT,MPROT,MM
      PARAMETER (FMM1=197.327D0,
     &           MNEUT=939.5653D0,MPROT=938.2720D0,
     &           MM=2.0D0*MPROT*MNEUT/(MPROT+MNEUT))
      PARAMETER (NTHMAX=181)
      PARAMETER (NP1=10,NP2=35,NP=NP1+NP2+1)
      PARAMETER (JMAX=4,LCHMX2C=2+4*JMAX,LMAX=JMAX+1)
      COMMON /KINEM/ OMEGA,QBIG,ELABNN,P0,THEQ
      COMMON /CURRENT/ CHARGE(NP,LCHMX2C,-1:1),
     &                   PLUS(NP,LCHMX2C,-1:1)
      LOGICAL COUPLC
      COMMON /PWTB2C/ INDX2C(LCHMX2C),NINDX2C,COUPLC(LCHMX2C)
      DOUBLE COMPLEX YLMP
      COMMON /YLM/ YLMP(NTHMAX,0:LMAX,-LMAX:LMAX)
      DOUBLE COMPLEX TRAMPRHO, TRAMPPLUS, TRAMPMINUS
      DOUBLE COMPLEX TRAMPRHO2, TRAMPPLUS2, TRAMPMINUS2
c     COMMON /TRAMPPW/ TRAMPRHO(ith,im1,im2,md) etc.
      COMMON /TRAMPPW/ TRAMPRHO(NTHMAX,2,2,-1:1),
     &                TRAMPPLUS(NTHMAX,2,2,-1:1),
     &               TRAMPMINUS(NTHMAX,2,2,-1:1)
      COMMON /TRAMPFU/ TRAMPRHO2(NTHMAX,2,2,-1:1),
     &                TRAMPPLUS2(NTHMAX,2,2,-1:1),
     &               TRAMPMINUS2(NTHMAX,2,2,-1:1)
      COMMON /DETECT/ MTTWE1,MTTWE2
C
C**** Spherical harmonics for all angles p
      CALL YKQFN
c     write(6,*) ' YLMP ', YLMP
C
C**** TRAMPRHO
      do md=-1,1
      do im2=1,2
      do im1=1,2
      do ith=1,NTHMAX
      TRAMPRHO(ith,im1,im2,md)=(0.,0.)
      enddo
      enddo
      enddo
      enddo
      DO 1 IALS=1,NINDX2C
      CALL UNPACK1(4,LS,ISS,JS,ITS,INDX2C(IALS),IZ,IZ,IZ,IZ)
      LS2=2*LS
      ISS2=2*ISS
      JS2=2*JS
      ITS2=2*ITS
      DO 2 MD=-1,1
      MD2=2*MD
      MP=MD  ! m' for rho
      IF (ABS(MP).GT.JS) GOTO 2
      MP2=2*MP
      DO 3 im2=1,2
      M22=2*IM2-3
      DO 4 im1=1,2
      M12=2*IM1-3
      MUP2=MP2-M12-M22
      MUP=MUP2/2
      IF (ABS(MUP2).GT.LS2) GOTO 4
      IF (ABS(MP2-MUP2).GT.ISS2) GOTO 4
      CCC=CLEBSCH(1,1,ISS2,M12,M22,M12+M22)
     &   *CLEBSCH(1,1,ITS2,MTTWE1,MTTWE2,MTTWE1+MTTWE2)
     &   *CLEBSCH(LS2,ISS2,JS2,MUP2,MP2-MUP2,MP2)
     &   *CHARGE(NP,IALS,MD)
      DO 5 ITH=1,NTHMAX
      TRAMPRHO(ith,im1,im2,md)=TRAMPRHO(ith,im1,im2,md)
     &  + CCC*YLMP(ITH,LS,MUP)
    5 CONTINUE
    4 CONTINUE
    3 CONTINUE
    2 CONTINUE
    1 CONTINUE
C
C**** TRAMPPLUS
      do md=-1,1
      do im2=1,2
      do im1=1,2
      do ith=1,NTHMAX
      TRAMPPLUS(ith,im1,im2,md)=(0.,0.)
      enddo
      enddo
      enddo
      enddo
      DO 11 IALS=1,NINDX2C
      CALL UNPACK1(4,LS,ISS,JS,ITS,INDX2C(IALS),IZ,IZ,IZ,IZ)
      LS2=2*LS
      ISS2=2*ISS
      JS2=2*JS
      ITS2=2*ITS
      DO 13 MD=-1,1
      MD2=2*MD
      MP=MD+1  ! m' for j_+
      IF (ABS(MP).GT.JS) GOTO 13
      MP2=2*MP
      DO 14 im2=1,2
      M22=2*IM2-3
      DO 15 im1=1,2
      M12=2*IM1-3
      MUP2=MP2-M12-M22
      MUP=MUP2/2
      IF (ABS(MP-MUP).GT.ISS) GOTO 15
      CCC=CLEBSCH(1,1,ISS2,M12,M22,MP2-MUP2)
     &   *CLEBSCH(1,1,ITS2,MTTWE1,MTTWE2,0)
     &   *CLEBSCH(LS2,ISS2,JS2,MUP2,MP2-MUP2,MP2)
     &   *PLUS(NP,IALS,MD)
      DO 16 ITH=1,NTHMAX
      TRAMPPLUS(ith,im1,im2,md)=TRAMPPLUS(ith,im1,im2,md)
     &  + CCC*YLMP(ITH,LS,MUP)
   16 CONTINUE
   15 CONTINUE
   14 CONTINUE
   13 CONTINUE
   11 CONTINUE
C
C**** TRAMPMINUS
      do md=-1,1
      do im2=1,2
      do im1=1,2
      do ith=1,NTHMAX
      TRAMPMINUS(ith,im1,im2,md)=(0.,0.)
      enddo
      enddo
      enddo
      enddo
      DO 111 IALS=1,NINDX2C
      CALL UNPACK1(4,LS,ISS,JS,ITS,INDX2C(IALS),IZ,IZ,IZ,IZ)
      LS2=2*LS
      ISS2=2*ISS
      JS2=2*JS
      ITS2=2*ITS
      DO 113 MD=-1,1
      MD2=2*MD
      MP=MD-1  ! m'
      IF (ABS(MP).GT.JS) GOTO 113
      MP2=2*MP
      DO 114 im2=1,2
      M22=2*IM2-3
      DO 115 im1=1,2
      M12=2*IM1-3
      MUP2=MP2-M12-M22
      MUP=MUP2/2
      IF (ABS(MP-MUP).GT.ISS) GOTO 115
      CCC=CLEBSCH(1,1,ISS2,M12,M22,MP2-MUP2)
     &   *CLEBSCH(1,1,ITS2,MTTWE1,MTTWE2,0)
     &   *CLEBSCH(LS2,ISS2,JS2,MUP2,MP2-MUP2,MP2)
     &   *PLUS(NP,IALS,-MD)
     &   *(-1.)**(LS+JS-1)  !  symmetry properties are used !!!!
      DO 116 ITH=1,NTHMAX
      TRAMPMINUS(ith,im1,im2,md)=TRAMPMINUS(ith,im1,im2,md)
     &  + CCC*YLMP(ITH,LS,MUP)
  116 CONTINUE
  115 CONTINUE
  114 CONTINUE
  113 CONTINUE
  111 CONTINUE
C
c     write(6,*) ' TRAMPRHO ', TRAMPRHO
c     write(6,*) ' TRAMPPLUS ', TRAMPPLUS
c     write(6,*) ' TRAMPMINUS ', TRAMPMINUS
c     do ith=1,nthmax
c     write(6,333) ith,
c    &                real(TRAMPRHO(ith,1,1,-1)),
c    &                real(TRAMPRHO(ith,1,1, 0)),
c    &                real(TRAMPRHO(ith,1,1, 1)),
c    &                real(TRAMPRHO(ith,1,2,-1)),
c    &                real(TRAMPRHO(ith,1,2, 0)),
c    &                real(TRAMPRHO(ith,1,2, 1)),
c    &                real(TRAMPRHO(ith,2,1,-1)),
c    &                real(TRAMPRHO(ith,2,1, 0)),
c    &                real(TRAMPRHO(ith,2,1, 1)),
c    &                real(TRAMPRHO(ith,2,2,-1)),
c    &                real(TRAMPRHO(ith,2,2, 0)),
c    &                real(TRAMPRHO(ith,2,2, 1))
c     enddo
c 333 format(i5,12E10.2)
c     do ith=1,nthmax
c     rr1=0.0
c     rr2=0.0
c     do md=-1,1
c     do im2=1,2
c     do im1=1,2
c     rr1=rr1 + cabs(TRAMPRHO(ith,im1,im2,md))**2
c     rr2=rr2 + cabs(TRAMPPLUS(ith,im1,im2,md))**2
c    &        + cabs(TRAMPMINUS(ith,im1,im2,md))**2
c     enddo
c     enddo
c     enddo
c     write(6,*) ith, rr1,rr2
c     enddo
C
      do md=-1,1
      do im2=1,2
      do im1=1,2
      do ith=1,NTHMAX
      TRAMPRHO2(ith,im1,im2,md)=TRAMPRHO(ith,im1,im2,md)
      TRAMPPLUS2(ith,im1,im2,md)=TRAMPPLUS(ith,im1,im2,md)
      TRAMPMINUS2(ith,im1,im2,md)=TRAMPMINUS(ith,im1,im2,md)
      enddo
      enddo
      enddo
      enddo
C
C
      RETURN
      END

      SUBROUTINE CROSSSECTION
C     =============
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DOUBLE PRECISION FMM1,MNEUT,MPROT,MM
      PARAMETER (FMM1=197.327D0,
     &           MNEUT=939.5653D0,MPROT=938.2720D0,
     &           MM=2.0D0*MPROT*MNEUT/(MPROT+MNEUT))
      PARAMETER (NTHMAX=181)
      COMMON /KINEM/ OMEGA,QBIG,ELABNN,P0,THEQ
      COMMON /KINEM12/ PM1(NTHMAX),THM1(NTHMAX),PM2(NTHMAX),
     &                 THM2(NTHMAX),THPLM1(NTHMAX)
      COMMON /WINKEL/ THP(NTHMAX),PHP(NTHMAX)  ! rad
      COMMON /XSEC/ PHSFAC(NTHMAX),CORR,CAPTURE
      DOUBLE COMPLEX TRAMPRHO, TRAMPPLUS, TRAMPMINUS, TRAMPDIF,TRAMPSUM
c     COMMON /TRAMPPW/ TRAMPRHO(ith,im1,im2,md) etc.
      COMMON /TRAMPPW/ TRAMPRHO(NTHMAX,2,2,-1:1),
     &                TRAMPPLUS(NTHMAX,2,2,-1:1),
     &               TRAMPMINUS(NTHMAX,2,2,-1:1)
      DOUBLE COMPLEX TRAMPRHO2, TRAMPPLUS2, TRAMPMINUS2, TRAMPDIF2,
     &               TRAMPSUM2
      COMMON /TRAMPFU/ TRAMPRHO2(NTHMAX,2,2,-1:1),
     &                TRAMPPLUS2(NTHMAX,2,2,-1:1),
     &               TRAMPMINUS2(NTHMAX,2,2,-1:1)
      DIMENSION CROSS(NTHMAX), WL(NTHMAX), WT(NTHMAX), WTT(NTHMAX),
     &                        WTL(NTHMAX), WTP(NTHMAX), WTPLP(NTHMAX)
      DIMENSION PART1(NTHMAX),PART2(NTHMAX),PART3(NTHMAX),PART4(NTHMAX)
      DIMENSION CROSS2(NTHMAX), WL2(NTHMAX), WT2(NTHMAX), WTT2(NTHMAX),
     &                        WTL2(NTHMAX), WTP2(NTHMAX),
     &                        WTPLP2(NTHMAX)
      DIMENSION PART12(NTHMAX),PART22(NTHMAX),PART32(NTHMAX),
     &          PART42(NTHMAX)
C
C for deuteron analyzing powers
      DOUBLE COMPLEX TAUD(NTHMAX,-2:2,2), TAUD2(NTHMAX,-2:2,2) ! t(ith,q,k) <--> Tkq
      DIMENSION SIGMA(NTHMAX), SIGMA2(NTHMAX)
      INTEGER K,Q
C
C for photon analyzing powers
      DOUBLE PRECISION DENOM(NTHMAX),AX(NTHMAX),
     .                    AY(NTHMAX),AZ(NTHMAX)
      DOUBLE PRECISION DENOM2(NTHMAX),AX2(NTHMAX),
     .                    AY2(NTHMAX),AZ2(NTHMAX)
C
C for polarization of outgoing nucleon 1
      INTEGER QN,QN2,QI,QI2,QNS
      DOUBLE PRECISION DSMALL(NTHMAX,-1:1,-1:1,0:1)
      DOUBLE COMPLEX TKQ(NTHMAX,-1:1),TKQS(NTHMAX,-1:1,0:1)
      DOUBLE COMPLEX POLOUT(NTHMAX,-1:1,0:1)
      DOUBLE COMPLEX TKQ2(NTHMAX,-1:1),TKQS2(NTHMAX,-1:1,0:1)
      DOUBLE COMPLEX POLOUT2(NTHMAX,-1:1,0:1)
      DOUBLE PRECISION POLN(NTHMAX)
      DOUBLE PRECISION POLN2(NTHMAX)
C for np capture observables
      DOUBLE PRECISION D(NTHMAX,-1:1,-1:1)
      DOUBLE PRECISION D12(NTHMAX,2,2)
      DOUBLE COMPLEX TRAMPJP(NTHMAX,2,2,-1:1),
     &               TRAMPJM(NTHMAX,2,2,-1:1),
     &               TRAMPJP2(NTHMAX,2,2,-1:1),
     &               TRAMPJM2(NTHMAX,2,2,-1:1)
      DOUBLE COMPLEX T1,T2,T12,T22,T1P,T2P,T1P2,T2P2
      DOUBLE COMPLEX TAUN1CAP(NTHMAX,-1:1)
      DOUBLE COMPLEX TAUN1CAP2(NTHMAX,-1:1)
      DOUBLE PRECISION AYN, AYN2
C
      PI=DACOS(-1.0d0)
      FF=PI/180.0
      AMN=MM/FMM1 ! 1/fm
C
C**** For unpolarized photon beam and deuteron target
C**** Final loop - factors
C**** Factor  2. comes from the antisymmetrizing operator
C****                       1./SQRT(2.) * ( 1 - P_12 )
C**** Factor  1./3. comes from average over initial deuteron spin states
      FFF= 2.*PI**2/137.0359895/QBIG
      FFF= FFF*CORR
      WRITE(16,*)
     &   '# PHOTODISINTEGRATION CROSS SECTION in microbarns/sr '
      WRITE(16,*) ' # THCM,CROSS,CROSS2,P1,P2,PHIP '
      DO ITH=1,NTHMAX
      WT(ITH)=0.0
      CROSS(ITH)=0.0
      WT2(ITH)=0.0
      CROSS2(ITH)=0.0
      DO MD=-1,1
      DO IM2=1,2
      DO IM1=1,2
      WT(ITH)=WT(ITH) + CDABS(TRAMPPLUS(ith,im1,im2,md))**2
     &                + CDABS(TRAMPMINUS(ith,im1,im2,md))**2
      WT2(ITH)=WT2(ITH) + CDABS(TRAMPPLUS2(ith,im1,im2,md))**2
     &                + CDABS(TRAMPMINUS2(ith,im1,im2,md))**2
      ENDDO ! im1
      ENDDO ! im2
      ENDDO ! md
      CROSS(ITH)=  FFF*  WT(ITH)* 2. / 3.  * PHSFAC(ITH) *10000.0
      CROSS2(ITH)= FFF* WT2(ITH)* 2. / 3.  * PHSFAC(ITH) *10000.0
      WRITE(16,101) THPLM1(ITH),CROSS(ITH),CROSS2(ITH),
     &             PM1(ITH),PM2(ITH),PHP(ITH)/FF
  101 FORMAT(F10.4,2E12.4,3F8.2)
      ENDDO ! ith
C
      CROSS_INC= 0.0
      CROSS2_INC= 0.0
      DO ITH=1,NTHMAX
      CROSS_INC= CROSS_INC +
     &       0.5*CROSS(ITH+1)*SIN(THM1(ITH+1))*2*PI*PI/180.0 +
     &       0.5*CROSS(ITH)*SIN(THM1(ITH))*2*PI*PI/180.0
      CROSS2_INC= CROSS2_INC +
     &       0.5*CROSS2(ITH+1)*SIN(THM1(ITH+1))*2*PI*PI/180.0 +
     &       0.5*CROSS2(ITH)*SIN(THM1(ITH))*2*PI*PI/180.0
      ENDDO
C
      WRITE(16,*) ' # CM integrated CROSS SECTION in microbarns '
      WRITE(16,*) ' # CROSS=  ',CROSS_INC
      WRITE(16,*) ' # CROSS2= ',CROSS2_INC
C
      WRITE(16,*)
     & ' # Integrated capture CROSS SECTION in microbarns '
      WRITE(16,*) ' # capture CROSS=  ',CROSS_INC*CAPTURE
      WRITE(16,*) ' # capture CROSS2= ',CROSS2_INC*CAPTURE
C
C
      DO 2 I1=1,2
      DO 2 I2=-2,2
      DO 2 I3=1,NTHMAX
      TAUD(I3,I2,I1)=(0.,0.)
    2 TAUD2(I3,I2,I1)=(0.,0.)
      DO 3 I1=1,NTHMAX
      SIGMA(I1)=0.
    3 SIGMA2(I1)=0.
C
      DO 50 MD=-1,1
        MD2=2*MD
        DO 50 K=1,2
          DO 50 Q=-K,K
            MDP2=MD2-2*Q
            MDP=MDP2/2
            IF (MDP .LE. 1.AND.MDP.GE.-1 ) THEN
              XD=(-1)**(1-MD2/2)*CLEBSCH(2,2,2*K,MDP2,-MD2,-2*Q)
                DO 60 IM1=1,2
                  DO 60 IM2=1,2
         DO 60 I=1,NTHMAX
                    TAUD(I,Q,K)=TAUD(I,Q,K)+XD*
     . ( TRAMPPLUS(I,IM1,IM2,MD)*DCONJG(TRAMPPLUS(I,IM1,IM2,MDP))
     . + TRAMPMINUS(I,IM1,IM2,MD)*DCONJG(TRAMPMINUS(I,IM1,IM2,MDP)) )
                    TAUD2(I,Q,K)=TAUD2(I,Q,K)+XD*
     . ( TRAMPPLUS2(I,IM1,IM2,MD)*DCONJG(TRAMPPLUS2(I,IM1,IM2,MDP))
     . + TRAMPMINUS2(I,IM1,IM2,MD)*DCONJG(TRAMPMINUS2(I,IM1,IM2,MDP)) )
   60    CONTINUE
            END IF
   50       CONTINUE
C
        DO 100 MD=-1,1
          DO 100 IM1=1,2
            DO 100 IM2=1,2
        DO 100 I=1,NTHMAX
              SIGMA(I)=SIGMA(I)
     .                   + CDABS(  TRAMPPLUS(I,IM1,IM2,MD) )**2
     .                   + CDABS( TRAMPMINUS(I,IM1,IM2,MD) )**2
              SIGMA2(I)=SIGMA2(I)
     .                   + CDABS(  TRAMPPLUS2(I,IM1,IM2,MD) )**2
     .                   + CDABS( TRAMPMINUS2(I,IM1,IM2,MD) )**2
!             WRITE(6,*)'SIGMA2', I,SIGMA2(I)
  100       CONTINUE
      DO 300 K=1,2
        DO 300 Q=-K,K
        DO 300 I=1,NTHMAX
          TAUD(I,Q,K)=TAUD(I,Q,K)*SQRT(3.)*(-1)**Q/SIGMA(I)
          TAUD2(I,Q,K)=TAUD2(I,Q,K)*SQRT(3.)*(-1)**Q/SIGMA2(I)
  300       CONTINUE
C
      WRITE(26,*) ' # THCM,T11D,T20D,T21D,T22D,T11D2,T20D2,T21D2,T22D2'
      DO 1000 ITH=1,NTHMAX
       T11D=-DIMAG(TAUD(ITH,1,1))
       T20D=DREAL(TAUD(ITH,0,2))
       T21D=DREAL(TAUD(ITH,1,2))
       T22D=DREAL(TAUD(ITH,2,2))
       T11D2=-DIMAG(TAUD2(ITH,1,1))
       T20D2=DREAL(TAUD2(ITH,0,2))
       T21D2=DREAL(TAUD2(ITH,1,2))
       T22D2=DREAL(TAUD2(ITH,2,2))
       WRITE(26,600) THM1(ITH),T11D,T20D,T21D,T22D,
     &                T11D2,T20D2,T21D2,T22D2
  600  FORMAT(F10.4,9E12.4)
 1000  CONTINUE
C
C
C
        DO ITH=1,NTHMAX
        DENOM(ITH)=0.0
           AX(ITH)=0.0
           AY(ITH)=0.0
           AZ(ITH)=0.0
        DENOM2(ITH)=0.0
           AX2(ITH)=0.0
           AY2(ITH)=0.0
           AZ2(ITH)=0.0
        ENDDO
C
      DO MD=-1,1
      DO IM2=1,2
      DO IM1=1,2
      DO ITH=1,NTHMAX
        DENOM(ITH)=DENOM(ITH) +
     &  ( CDABS( TRAMPPLUS(ITH,IM1,IM2,MD) )**2
     &  + CDABS( TRAMPMINUS(ITH,IM1,IM2,MD) )**2 )
        AX(ITH)=AX(ITH) +
     &  2.0*DREAL( TRAMPPLUS(ITH,IM1,IM2,MD)*
     &      DCONJG( TRAMPMINUS(ITH,IM1,IM2,MD) ) )
        AY(ITH)=AY(ITH) +
     &  2.0*DIMAG( TRAMPPLUS(ITH,IM1,IM2,MD)*
     &      DCONJG( TRAMPMINUS(ITH,IM1,IM2,MD) ) )
        AZ(ITH)=AZ(ITH) +
     &  ( CDABS( TRAMPPLUS(ITH,IM1,IM2,MD) )**2
     &  - CDABS( TRAMPMINUS(ITH,IM1,IM2,MD) )**2 )
c
        DENOM2(ITH)=DENOM2(ITH) +
     &  ( CDABS( TRAMPPLUS2(ITH,IM1,IM2,MD) )**2
     &  + CDABS( TRAMPMINUS2(ITH,IM1,IM2,MD) )**2 )
        AX2(ITH)=AX2(ITH) +
     &  2.0*DREAL( TRAMPPLUS2(ITH,IM1,IM2,MD)*
     &      DCONJG( TRAMPMINUS2(ITH,IM1,IM2,MD) ) )
        AY2(ITH)=AY2(ITH) +
     &  2.0*DIMAG( TRAMPPLUS2(ITH,IM1,IM2,MD)*
     &      DCONJG( TRAMPMINUS2(ITH,IM1,IM2,MD) ) )
        AZ2(ITH)=AZ2(ITH) +
     &  ( CDABS( TRAMPPLUS2(ITH,IM1,IM2,MD) )**2
     &  - CDABS( TRAMPMINUS2(ITH,IM1,IM2,MD) )**2 )
      ENDDO ! ITH
      ENDDO ! IM1
      ENDDO ! IM2
      ENDDO ! MD
C
      WRITE(36,*) ' # THCM,AX,AY,AZ,AX2,AY2,AZ2 '
      DO ITH=1,NTHMAX
        AX(ITH)=AX(ITH)/DENOM(ITH)
        AY(ITH)=AY(ITH)/DENOM(ITH)
        AZ(ITH)=AZ(ITH)/DENOM(ITH)
        AX2(ITH)=AX2(ITH)/DENOM2(ITH)
        AY2(ITH)=AY2(ITH)/DENOM2(ITH)
        AZ2(ITH)=AZ2(ITH)/DENOM2(ITH)
        WRITE(36,600) THM1(ITH),AX(ITH),AY(ITH),AZ(ITH),
     &                         AX2(ITH),AY2(ITH),AZ2(ITH)
      ENDDO ! ITH
C
C
C
      DO ITH=1,NTHMAX
      THT=THM1(ITH) ! rad
      DSMALL(ITH, 0, 0,0)=1.
      DSMALL(ITH, 0, 0,1)=COS(THT)
      DSMALL(ITH, 0, 1,1)=SIN(THT)/SQRT(2.)
      DSMALL(ITH,-1, 0,1)= DSMALL(ITH, 0, 1,1)
      DSMALL(ITH, 0,-1,1)=-DSMALL(ITH, 0, 1,1)
      DSMALL(ITH, 1, 0,1)=-DSMALL(ITH, 0, 1,1)
      DSMALL(ITH, 1, 1,1)=COS(THT/2.)**2
      DSMALL(ITH,-1,-1,1)= DSMALL(ITH, 1, 1,1)
      DSMALL(ITH, 1,-1,1)=SIN(THT/2.)**2
      DSMALL(ITH,-1, 1,1)= DSMALL(ITH, 1,-1,1)
      ENDDO
C
      DO 10 KN=0,1
      KN2=2*KN
      DO 11 QN=-KN,KN
      QN2=2*QN
      DO 15 I=1,NTHMAX
      TKQ(I,QN)=(0.,0.)
   15 TKQ2(I,QN)=(0.,0.)
      DO 200 INO=1,2
      MNO2=2*(INO-1)-1
      MNO2S=MNO2+QN2
      INOS=(MNO2S+1)/2+1
      IF(INOS.LT.1.OR.INOS.GT.2) GOTO 200
      PHASE=(-1)**((1-MNO2)/2)*CLEBSCH(1,1,KN2,MNO2S,-MNO2,QN2)
      DO 25 MD=-1,1
      DO 25 IM2=1,2
      DO 25 I=1,NTHMAX
      TKQ(I,QN)=TKQ(I,QN)+PHASE*
     . ( TRAMPPLUS(I,INO,IM2,MD)*DCONJG(TRAMPPLUS(I,INOS,IM2,MD))
     . + TRAMPMINUS(I,INO,IM2,MD)*DCONJG(TRAMPMINUS(I,INOS,IM2,MD)) )
   25 TKQ2(I,QN)=TKQ2(I,QN)+PHASE*
     . ( TRAMPPLUS2(I,INO,IM2,MD)*DCONJG(TRAMPPLUS2(I,INOS,IM2,MD))
     . + TRAMPMINUS2(I,INO,IM2,MD)*DCONJG(TRAMPMINUS2(I,INOS,IM2,MD)) )
  200 CONTINUE
      PHASE=SQRT(2.)
      DO 30 I=1,NTHMAX
      TKQ(I,QN)=PHASE*TKQ(I,QN)/SIGMA(I)
   30 TKQ2(I,QN)=PHASE*TKQ2(I,QN)/SIGMA2(I)
   11 CONTINUE
      DO 12 QN=-KN,KN
      DO 13 I=1,NTHMAX
      TKQS(I,QN,KN)=(0.,0.)
   13 TKQS2(I,QN,KN)=(0.,0.)
      DO 12 QNS=-KN,KN
      DO 12 I=1,NTHMAX
      TKQS(I,QN,KN)=TKQS(I,QN,KN)+DSMALL(I,QNS,QN,KN)*TKQ(I,QNS)
   12 TKQS2(I,QN,KN)=TKQS2(I,QN,KN)+DSMALL(I,QNS,QN,KN)*TKQ2(I,QNS)
      DO 14 QN=-KN,KN
      DO 14 I=1,NTHMAX
      POLOUT(I,QN,KN)=TKQS(I,QN,KN)
   14 POLOUT2(I,QN,KN)=TKQS2(I,QN,KN)
   10 CONTINUE
C
      DO I=1,NTHMAX
        POLN(I)= -SQRT(2.)*DIMAG(POLOUT(I,1,1))
        POLN2(I)= -SQRT(2.)*DIMAG(POLOUT2(I,1,1))
      ENDDO
      WRITE(46,*) ' # THCM,POLNOUT(y), POLNOUT2(y) '
      DO ITH=1,NTHMAX
        WRITE(46,600) THM1(ITH),POLN(ITH),POLN2(ITH)
      ENDDO

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCC  np capture CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DO I1=-1,1
      DO ITH=1,NTHMAX
      TAUN1CAP(ITH,I1)=(0.,0.)   ! for nucleon 1
      TAUN1CAP2(ITH,I1)=(0.,0.)
      ENDDO
      ENDDO

      DO ITH=1,NTHMAX
      SIGMA(ITH)=0.0
      SIGMA2(ITH)=0.0
      ENDDO
C

      DO MD=-1,1
      DO IM2=1,2
      DO IM1=1,2
      DO ITH=1,NTHMAX
      TRAMPJP(ith,im1,im2,md)=(0.0,0.0)
      TRAMPJM(ith,im1,im2,md)=(0.0,0.0)
      TRAMPJP2(ith,im1,im2,md)=(0.0,0.0)
      TRAMPJM2(ith,im1,im2,md)=(0.0,0.0)
      ENDDO ! ith
      ENDDO ! im1
      ENDDO ! im2
      ENDDO ! md

C###  D-functions are real !
      DO 5 ITH=1,NTHMAX
      THG= PI-THM1(ITH) ! the angle between gamma and neutron
      D(ITH,  0, 0 )= COS(THG)
      D(ITH,  0, 1 )= SIN(THG)/SQRT(2.)
      D(ITH, -1, 0 )= D(ITH, 0, 1 )
      D(ITH,  0,-1 )=-D(ITH, 0, 1 )
      D(ITH,  1, 0 )=-D(ITH, 0, 1 )
      THG2=0.5*THG
      D(ITH,  1, 1 )= COS(THG2)**2
      D(ITH, -1,-1 )= D(ITH, 1, 1 )
      D(ITH,  1,-1 )= SIN(THG2)**2
      D(ITH, -1, 1 )= D(ITH, 1,-1 )
C
      D12(ITH, 1, 1 )= COS(0.5*THG)
      D12(ITH, 1, 2 )= SIN(0.5*THG)
      D12(ITH, 2, 1 )=-SIN(0.5*THG)
      D12(ITH, 2, 2 )= COS(0.5*THG)
    5 CONTINUE
C
C
C___Rotation of deuteron and nucleon magnetic quantum numbers
      DO 6 IM1=1,2
       DO 6 IM2=1,2
        DO 6 MD=-1,1
         PHASE1=(-1.)**((2*IM1-3 + 2*IM2-3 + 2 + 2*MD)/2)
         PHASE2=(-1.)**((2*IM1-3 + 2*IM2-3 - 2 + 2*MD)/2)
         DO 6 IM1S=1,2
          DO 6 IM2S=1,2
           DO 6 MDS=-1,1
            DO 6 ITH=1,NTHMAX
             TRAMPJP(ITH,IM1,IM2,MD)=
     &       TRAMPJP(ITH,IM1,IM2,MD)+
     &       TRAMPMINUS(ITH,IM1S,IM2S,MDS)*D(ITH,MDS,-MD)
     &                                     *D12(ITH,IM2S,3-IM2)
     &                                     *D12(ITH,IM1S,3-IM1)
     &                                     *PHASE1
             TRAMPJP2(ITH,IM1,IM2,MD)=
     &       TRAMPJP2(ITH,IM1,IM2,MD)+
     &       TRAMPMINUS2(ITH,IM1S,IM2S,MDS)*D(ITH,MDS,-MD)
     &                                     *D12(ITH,IM2S,3-IM2)
     &                                     *D12(ITH,IM1S,3-IM1)
     &                                     *PHASE1
             TRAMPJM(ITH,IM1,IM2,MD)=
     &       TRAMPJM(ITH,IM1,IM2,MD)+
     &       TRAMPPLUS(ITH,IM1S,IM2S,MDS)*D(ITH,MDS,-MD)
     &                                     *D12(ITH,IM2S,3-IM2)
     &                                     *D12(ITH,IM1S,3-IM1)
     &                                     *PHASE2
             TRAMPJM2(ITH,IM1,IM2,MD)=
     &       TRAMPJM2(ITH,IM1,IM2,MD)+
     &       TRAMPPLUS2(ITH,IM1S,IM2S,MDS)*D(ITH,MDS,-MD)
     &                                     *D12(ITH,IM2S,3-IM2)
     &                                     *D12(ITH,IM1S,3-IM1)
     &                                     *PHASE2
    6 CONTINUE
C
C___Polarization observables in the system,
C___where the quantization axis (z) is parallel to the initial neutron momentum
C
      DO 20 MN=1,2
        MN2=2*MN-3
        DO 20 Q=-1,1
          MNP2=MN2-2*Q
          MNP=(MNP2+3)/2
           IF (MNP .LE. 2.AND.MNP.GE.1  ) THEN
            XN=(-1)**((1-MN2)/2)*CLEBSCH(1,1,2,MNP2,-MN2,-2*Q)
            DO 26 MD=-1,1
              DO 26 IM2=1,2
                 DO 26 ITH=1,NTHMAX
                  T1=  TRAMPJP(ITH,MN,IM2,MD)
                  T1P= TRAMPJP(ITH,MNP,IM2,MD)
                  T2=  TRAMPJM(ITH,MN,IM2,MD)
                  T2P= TRAMPJM(ITH,MNP,IM2,MD)
          TAUN1CAP(ITH,Q)= TAUN1CAP(ITH,Q)
     &              + (T1*DCONJG(T1P)+T2*DCONJG(T2P))*XN
                  T12=  TRAMPJP2(ITH,MN,IM2,MD)
                  T1P2= TRAMPJP2(ITH,MNP,IM2,MD)
                  T22=  TRAMPJM2(ITH,MN,IM2,MD)
                  T2P2= TRAMPJM2(ITH,MNP,IM2,MD)
          TAUN1CAP2(ITH,Q)= TAUN1CAP2(ITH,Q)
     &              + (T12*DCONJG(T1P2)+T22*DCONJG(T2P2))*XN
   26     CONTINUE
          END IF
   20     CONTINUE


      DO 103 MN=1,2
        DO 103 MD=-1,1
          DO 103 IM2=1,2
           DO 103 ITH=1,NTHMAX
            T1= TRAMPJP(ITH,MN,IM2,MD)
            T2= TRAMPJM(ITH,MN,IM2,MD)
            SIGMA(ITH)=SIGMA(ITH) + CDABS(T1)**2 + CDABS(T2)**2
            T12= TRAMPJP2(ITH,MN,IM2,MD)
            T22= TRAMPJM2(ITH,MN,IM2,MD)
            SIGMA2(ITH)=SIGMA2(ITH) + CDABS(T12)**2 + CDABS(T22)**2
  103       CONTINUE


      DO 201 Q=-1,1
        DO 201 ITH=1,NTHMAX
        TAUN1CAP(ITH,Q)=TAUN1CAP(ITH,Q)*SQRT(2.)*(-1)**Q/SIGMA(ITH)
        TAUN1CAP2(ITH,Q)=TAUN1CAP2(ITH,Q)*SQRT(2.)*(-1)**Q/SIGMA2(ITH)
  201   CONTINUE

      WRITE(56,*) ' # THCM, CROSSCAP, CROSSCAP2, AYN1CAP, AYN1CAP2 '
       DO ITH=1,NTHMAX
       AYN=-SQRT(2.)*DIMAG(TAUN1CAP(ITH,1))
       AYN2=-SQRT(2.)*DIMAG(TAUN1CAP2(ITH,1))
        WRITE(56,600) THM1(ITH), CAPTURE*CROSS(ITH),
     &               CAPTURE*CROSS2(ITH), AYN, AYN2
      ENDDO

      RETURN
      END

      subroutine hmd3(t,gep,gmp,gen,gmn)
c
c   Date: Thu, 15 Jul 2004 17:29:29 -0700 (PDT)
c   From: Hans-Werner Hammer <hammer@pooh.phys.washington.edu>
c   To: Walter.Gloeckle@tp2.ruhr-uni-bochum.de
c   Subject: EM form factors
c
c   Lieber Herr Gloeckle,
c
c   beigefuegt sind die Daten-Files fuer unseren besten Fit (rote
c    durchgezogene Linie) in Fig. 1 von  hep-ph/0312081.
c    Die Files enthalten die folgenden Groessen:
c
c    genthf4.dat:  Q^2[GeV^2],  G_E^n
c    gepthnf4.dat:  Q^2[GeV^2],  G_E^p/G_dipole
c    gmnthnf4.dat:  Q^2[GeV^2],  G_M^n/G_dipole/mu_n
c    gmpthnf4.dat:  Q^2[GeV^2],  G_M^p/G_dipole/mu_p
c
c    where G_dipole = 1/(1+Q^2/(0.71 GeV^2))^2
c          mu_n = -1.913
c          mu_p = 2.793
c
c
c        Viele Gruesse,
c          Hans-Werner Hammer
c
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      dimension qq(42),gen1(42),gep1(42),gmn1(42),gmp1(42),spl(42)
c
      qq( 1)= 0.1000000050E-02
      qq( 2)= 0.9999999780E-02
      qq( 3)= 0.1122018420E-01
      qq( 4)= 0.1258925440E-01
      qq( 5)= 0.1412537510E-01
      qq( 6)= 0.1584893280E-01
      qq( 7)= 0.1778279430E-01
      qq( 8)= 0.1995262320E-01
      qq( 9)= 0.2238721210E-01
      qq(10)= 0.2511886510E-01
      qq(11)= 0.2818382900E-01
      qq(12)= 0.3162277860E-01
      qq(13)= 0.3548134120E-01
      qq(14)= 0.3981071710E-01
      qq(15)= 0.4466836150E-01
      qq(16)= 0.5011872570E-01
      qq(17)= 0.5623413250E-01
      qq(18)= 0.6309573350E-01
      qq(19)= 0.7079458240E-01
      qq(20)= 0.7943282280E-01
      qq(21)= 0.8912509680E-01
      qq(22)= 0.1000000010E+00
      qq(23)= 0.1122018470E+00
      qq(24)= 0.1258925500E+00
      qq(25)= 0.1412537540E+00
      qq(26)= 0.1584893320E+00
      qq(27)= 0.1778279540E+00
      qq(28)= 0.1995262350E+00
      qq(29)= 0.2238721250E+00
      qq(30)= 0.2511886660E+00
      qq(31)= 0.2818382980E+00
      qq(32)= 0.3162277940E+00
      qq(33)= 0.3548133970E+00
      qq(34)= 0.3981072010E+00
      qq(35)= 0.4466836150E+00
      qq(36)= 0.5011872650E+00
      qq(37)= 0.5623413320E+00
      qq(38)= 0.6309573650E+00
      qq(39)= 0.7079458240E+00
      qq(40)= 0.7943282720E+00
      qq(41)= 0.8912509680E+00
      qq(42)= 0.1000000120E+01
      gen1( 1)= 0.4887843970E-03
      gen1( 2)= 0.4651191179E-02
      gen1( 3)= 0.5184431560E-02
      gen1( 4)= 0.5774374586E-02
      gen1( 5)= 0.6425972562E-02
      gen1( 6)= 0.7144348696E-02
      gen1( 7)= 0.7934731431E-02
      gen1( 8)= 0.8802369237E-02
      gen1( 9)= 0.9752419777E-02
      gen1(10)= 0.1078981813E-01
      gen1(11)= 0.1191910356E-01
      gen1(12)= 0.1314422954E-01
      gen1(13)= 0.1446833462E-01
      gen1(14)= 0.1589348167E-01
      gen1(15)= 0.1742037758E-01
      gen1(16)= 0.1904806867E-01
      gen1(17)= 0.2077363804E-01
      gen1(18)= 0.2259188518E-01
      gen1(19)= 0.2449505590E-01
      gen1(20)= 0.2647260390E-01
      gen1(21)= 0.2851101197E-01
      gen1(22)= 0.3059372678E-01
      gen1(23)= 0.3270118311E-01
      gen1(24)= 0.3481098264E-01
      gen1(25)= 0.3689820692E-01
      gen1(26)= 0.3893588111E-01
      gen1(27)= 0.4089555889E-01
      gen1(28)= 0.4274802282E-01
      gen1(29)= 0.4446407408E-01
      gen1(30)= 0.4601530358E-01
      gen1(31)= 0.4737488925E-01
      gen1(32)= 0.4851827770E-01
      gen1(33)= 0.4942376539E-01
      gen1(34)= 0.5007296801E-01
      gen1(35)= 0.5045112222E-01
      gen1(36)= 0.5054729804E-01
      gen1(37)= 0.5035450682E-01
      gen1(38)= 0.4986979440E-01
      gen1(39)= 0.4909434170E-01
      gen1(40)= 0.4803357646E-01
      gen1(41)= 0.4669738188E-01
      gen1(42)= 0.4510026425E-01
      gep1( 1)= 0.9969341160E+00
      gep1( 2)= 0.9700371620E+00
      gep1( 3)= 0.9664841290E+00
      gep1( 4)= 0.9625232220E+00
      gep1( 5)= 0.9581110480E+00
      gep1( 6)= 0.9532004000E+00
      gep1( 7)= 0.9477399590E+00
      gep1( 8)= 0.9416745900E+00
      gep1( 9)= 0.9349449870E+00
      gep1(10)= 0.9274879690E+00
      gep1(11)= 0.9192365410E+00
      gep1(12)= 0.9101203680E+00
      gep1(13)= 0.9000660780E+00
      gep1(14)= 0.8889981510E+00
      gep1(15)= 0.8768398170E+00
      gep1(16)= 0.8635141250E+00
      gep1(17)= 0.8489455580E+00
      gep1(18)= 0.8330616950E+00
      gep1(19)= 0.8157955410E+00
      gep1(20)= 0.7970879080E+00
      gep1(21)= 0.7768900390E+00
      gep1(22)= 0.7551667690E+00
      gep1(23)= 0.7318998580E+00
      gep1(24)= 0.7070911530E+00
      gep1(25)= 0.6807661650E+00
      gep1(26)= 0.6529769300E+00
      gep1(27)= 0.6238050460E+00
      gep1(28)= 0.5933635230E+00
      gep1(29)= 0.5617982750E+00
      gep1(30)= 0.5292881130E+00
      gep1(31)= 0.4960436520E+00
      gep1(32)= 0.4623042640E+00
      gep1(33)= 0.4283341770E+00
      gep1(34)= 0.3944161240E+00
      gep1(35)= 0.3608441350E+00
      gep1(36)= 0.3279148340E+00
      gep1(37)= 0.2959181070E+00
      gep1(38)= 0.2651273610E+00
      gep1(39)= 0.2357904760E+00
      gep1(40)= 0.2081215380E+00
      gep1(41)= 0.1822939960E+00
      gep1(42)= 0.1584360750E+00
      gmn1( 1)=-0.1906855340E+01
      gmn1( 2)=-0.1852136850E+01
      gmn1( 3)=-0.1844980120E+01
      gmn1( 4)=-0.1837020400E+01
      gmn1( 5)=-0.1828176620E+01
      gmn1( 6)=-0.1818361160E+01
      gmn1( 7)=-0.1807480340E+01
      gmn1( 8)=-0.1795433880E+01
      gmn1( 9)=-0.1782116290E+01
      gmn1(10)=-0.1767416360E+01
      gmn1(11)=-0.1751218080E+01
      gmn1(12)=-0.1733401660E+01
      gmn1(13)=-0.1713844780E+01
      gmn1(14)=-0.1692423820E+01
      gmn1(15)=-0.1669015290E+01
      gmn1(16)=-0.1643499490E+01
      gmn1(17)=-0.1615760920E+01
      gmn1(18)=-0.1585692880E+01
      gmn1(19)=-0.1553199650E+01
      gmn1(20)=-0.1518200400E+01
      gmn1(21)=-0.1480632900E+01
      gmn1(22)=-0.1440457580E+01
      gmn1(23)=-0.1397662040E+01
      gmn1(24)=-0.1352264170E+01
      gmn1(25)=-0.1304317710E+01
      gmn1(26)=-0.1253915070E+01
      gmn1(27)=-0.1201191190E+01
      gmn1(28)=-0.1146325710E+01
      gmn1(29)=-0.1089545250E+01
      gmn1(30)=-0.1031123280E+01
      gmn1(31)=-0.9713792200E+00
      gmn1(32)=-0.9106752280E+00
      gmn1(33)=-0.8494116070E+00
      gmn1(34)=-0.7880197170E+00
      gmn1(35)=-0.7269531490E+00
      gmn1(36)=-0.6666770580E+00
      gmn1(37)=-0.6076566580E+00
      gmn1(38)=-0.5503436330E+00
      gmn1(39)=-0.4951642750E+00
      gmn1(40)=-0.4425060750E+00
      gmn1(41)=-0.3927070500E+00
      gmn1(42)=-0.3460458820E+00
      gmp1( 1)= 0.2784225700E+01
      gmp1( 2)= 0.2707937960E+01
      gmp1( 3)= 0.2697918420E+01
      gmp1( 4)= 0.2686764240E+01
      gmp1( 5)= 0.2674357180E+01
      gmp1( 6)= 0.2660570380E+01
      gmp1( 7)= 0.2645267010E+01
      gmp1( 8)= 0.2628300190E+01
      gmp1( 9)= 0.2609513040E+01
      gmp1(10)= 0.2588740830E+01
      gmp1(11)= 0.2565809010E+01
      gmp1(12)= 0.2540535690E+01
      gmp1(13)= 0.2512733940E+01
      gmp1(14)= 0.2482212310E+01
      gmp1(15)= 0.2448777910E+01
      gmp1(16)= 0.2412239070E+01
      gmp1(17)= 0.2372410060E+01
      gmp1(18)= 0.2329114680E+01
      gmp1(19)= 0.2282191040E+01
      gmp1(20)= 0.2231498000E+01
      gmp1(21)= 0.2176920890E+01
      gmp1(22)= 0.2118379120E+01
      gmp1(23)= 0.2055833100E+01
      gmp1(24)= 0.1989292030E+01
      gmp1(25)= 0.1918820740E+01
      gmp1(26)= 0.1844548340E+01
      gmp1(27)= 0.1766673210E+01
      gmp1(28)= 0.1685468790E+01
      gmp1(29)= 0.1601286530E+01
      gmp1(30)= 0.1514556880E+01
      gmp1(31)= 0.1425786970E+01
      gmp1(32)= 0.1335554960E+01
      gmp1(33)= 0.1244500880E+01
      gmp1(34)= 0.1153313400E+01
      gmp1(35)= 0.1062713740E+01
      gmp1(36)= 0.9734361170E+00
      gmp1(37)= 0.8862064480E+00
      gmp1(38)= 0.8017200830E+00
      gmp1(39)= 0.7206206320E+00
      gmp1(40)= 0.6434792280E+00
      gmp1(41)= 0.5707783700E+00
      gmp1(42)= 0.5028988120E+00
c
      CALL SPLFAK(qq,42)
      CALL SPLMOD(qq,42,t,spl)
      gep=0.0
      gmp=0.0
      gen=0.0
      gmn=0.0
      do i=1,42
      gen=gen+spl(i)*gen1(i)
      gep=gep+spl(i)*gep1(i)
      gmn=gmn+spl(i)*gmn1(i)
      gmp=gmp+spl(i)*gmp1(i)
      enddo
c
      return
      end
      SUBROUTINE RESCATTER
C     =============
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      CALL TMATRV_IEE_DEUTERON
      CALL INTEGRAL_P0
      CALL RESCAMP
C
      RETURN
      END
      SUBROUTINE TMATRV_IEE_DEUTERON
C     =============
C !!!! Be careful about LCHMX2=2+3*JMAX !!!!!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      LOGICAL SHELCON
      PARAMETER (NPS=72,NDEUT=100,JMAX=4,LCHMX2=2+3*JMAX,
     . LCHMX2C=2+4*JMAX,PBAR=30.0,
     . NPP1=52,NPP2=20,NPP=NPP1+NPP2,PP1=3.0,PP2=10.0,PP3=50.0,
     . NP1=10,NP2=35,NPI=NP1+NP2,
     . IPOT=2,SHELCON=.FALSE.,
     . COF1=1.0,COF2=1.0,COF3=1.0)
      INTEGER Z
      LOGICAL COUPL,COUPLC
      DOUBLE PRECISION FMM1,MN,MP,MM
      PARAMETER (FMM1=197.327D0,
     &           MN=939.5653D0,MP=938.2720D0,MM=2.0D0*MP*MN/(MP+MN))
      CHARACTER*4 TEXT(2)
      DIMENSION PD(NDEUT)
      COMMON/PWTB2/ INDEX2(LCHMX2),NINDEX2,COUPL(LCHMX2)
      COMMON /PWTB2C/ INDEX2C(LCHMX2C),NINDEX2C,COUPLC(LCHMX2C)
      COMMON /POINTS/ PPOINT(NPS),GPOINT(NPS)
      COMMON /GAUSS2/ GPNEW(NPI),PNEW(NPI)
      COMMON/COF/ COF(LCHMX2)
      COMMON /KINEM/ OMEGA,QBIG,ELABNN,P0,THEQ
      INTEGER LEFT, RIGHT
      COMMON /INDX/ LEFT(4,LCHMX2), RIGHT(4,LCHMX2)
      DOUBLE COMPLEX TM
      COMMON/ TMP0/ TM(LCHMX2C,LCHMX2C,NPI+1)
      DATA Z/0/
      DATA TEXT/'UNCO','COUP'/
C
C***** IPOT=1--M,TJON POTENTIAL; IPOT=2--OTHERS POTENTIALS
C***** NPS-NUMBER OF P-POINTS FOR INTEGRATION OF LIPMANN-SCHWINGER
C***** EQUATION AND DEUTERON PROBLEM
C***** NPI-NUMBER OF P-POINTS WITH WITCH MELTFI AND ITERAC WILL CALCULATE
C***** ITERMEL=1
C
      ICTRL=0
cjg   OPEN(80,FILE='tmnp',STATUS='UNKNOWN',FORM='UNFORMATTED')
      OPEN(14,FILE='deuwf',STATUS='OLD',FORM='UNFORMATTED')
      READ(14) PD,EBDEUT
      REWIND (14)
      IF( NPP.NE.NPS ) STOP ' NPP.NE.NPS'
C*****PPOINT SHOULD BE THE SAME P-POINTS WITH WITCH DEUTERON WAS CALCULATED
      CALL TRNS(NPP1,NPP2,NPP,PP1,PP2,PP3,PPOINT,GPOINT)
C     DO 20 I=1,NPS
C  20 IF( ABS(PD(I)-PPOINT(I)).GT.1.E-7)   STOP 'PD.NE.PPOINT'
C
      WRITE(6,700) OMEGA
  700 FORMAT(1X,'    T-MATRIX FOR DEUTERON PHOTO-DISINTEGRATION '
     2 /5X,'FOR INCOMING CM PHOTON ENERGY=',F10.4,' (MeV)')
C
        AMNMEV=MM ! MeV
        EBDEUT=ABS(EBDEUT)
        ECMNN=(P0*FMM1)**2/AMNMEV
        ETOT=ECMNN
        WRITE(6,*) ' ECMNN,P0= ',ECMNN,P0
C
C*****PNEW ARE P-POINTS WITH WITCH MELTFI AND ITERAC WILL BE CALCULATED
cjg   CALL PPTMT2(PNEW,GPNEW,QMAX,1)
cjg   see COMMON/GAUSS2/
C
      CALL ALFA2
C     READ(5,*) COF1,COF2,COF3
      DO 5555 J=1,NINDEX2
      COF(J)=1.0
      IF(J.EQ.2) COF(J)=COF1
      IF(J.EQ.5) COF(J)=COF2
      IF(J.EQ.7) COF(J)=COF3
 5555 WRITE(6,*) ' J=',J,' COF=',COF(J)
C
      CALL CALV(IPOT)
C     WRITE (80) ITERMEL,ELAB,IPOT,QMAX
C
      write(6,*) ' po CALV '
      E=ETOT
      DO 180 IALFA=1,NINDEX2
C****
      IF (.NOT.COUPL(IALFA)) THEN
        IAC=0
        DO I2=1,NINDEX2C
        IF (INDEX2(IALFA).EQ.INDEX2C(I2)) IAC=I2
        ENDDO
        IF (IAC.EQ.0) STOP 'IAC'
        LEFT(1,IALFA)=IAC
        RIGHT(1,IALFA)=IAC
        WRITE(6,*) 'IALFA,IAC=',IALFA,IAC
        WRITE(6,*) ' IALFA,LEFT(1,IALFA),RIGHT(1,IALFA)=',
     &               IALFA,LEFT(1,IALFA),RIGHT(1,IALFA)
      ELSE
        IAC=0
        DO I2=1,NINDEX2C
        IF (INDEX2(IALFA).EQ.INDEX2C(I2)) IAC=I2
        ENDDO
        IF (IAC.EQ.0) STOP 'IAC'
        LEFT(1,IALFA)=IAC
        RIGHT(1,IALFA)=IAC
        LEFT(2,IALFA)=IAC
        RIGHT(3,IALFA)=IAC
        WRITE(6,*) 'IALFA,IAC=',IALFA,IAC
        CALL UNPACK1(4,L,IS,J,IT,INDEX2(IALFA),IZ,IZ,IZ,IZ)
        L1=L+2
        CALL PACK1(4,L1,IS,J,IT,IDX,IZ,IZ,IZ,IZ)
        IAC=0
        DO I2=1,NINDEX2C
        IF (IDX.EQ.INDEX2C(I2)) IAC=I2
        ENDDO
        IF (IAC.EQ.0) STOP 'IAC'
        LEFT(3,IALFA)=IAC
        LEFT(4,IALFA)=IAC
        RIGHT(2,IALFA)=IAC
        RIGHT(4,IALFA)=IAC
        WRITE(6,*) 'IALFA,IAC=',IALFA,IAC
        WRITE(6,*) ' IALFA,LEFT(1,IALFA),RIGHT(1,IALFA)=',
     &               IALFA,LEFT(1,IALFA),RIGHT(1,IALFA)
        WRITE(6,*) ' IALFA,LEFT(2,IALFA),RIGHT(2,IALFA)=',
     &               IALFA,LEFT(2,IALFA),RIGHT(2,IALFA)
        WRITE(6,*) ' IALFA,LEFT(3,IALFA),RIGHT(3,IALFA)=',
     &               IALFA,LEFT(3,IALFA),RIGHT(3,IALFA)
        WRITE(6,*) ' IALFA,LEFT(4,IALFA),RIGHT(4,IALFA)=',
     &               IALFA,LEFT(4,IALFA),RIGHT(4,IALFA)
      ENDIF
C****
      IOUT=1
      IF(IALFA.GT.1) IOUT=3
      CALL TMTRX(E,INDEX2(IALFA),COUPL(IALFA),SHELCON,ICTRL,IOUT,EBDEUT,
     1           IPOT,IALFA)
  180 CONTINUE
C
cjg   CLOSE (80)
C
!      DO IA1=1,LCHMX2C
!      DO IA2=1,LCHMX2C
!      WRITE(6,*) ' IA1,IA2= ',IA1,IA2
!      DO IP=1,NPI+1
!      WRITE(6,*) TM(IA1,IA2,IP)
!      ENDDO
!      ENDDO
!      ENDDO
C
      RETURN
      END
      SUBROUTINE TMTRX(E,INDEX,KOPPEL,SHELCON,ICTRL,IOUT,EBDEUT,
     1                 IPOT,IALFA)
C
C*****IF ICTRL.NE.0----THEN ALSO SHOULD BE E=0.0 AND IN THIS CASE
C     QUANTIETY T(=K)LA,LAS(P,0.;0.)/(0.)**LAS IS CALCULATED AND WRITTEN
C     TO TAPE 80
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER(NPS=72,PP3=50.0,
     . NP1=10,NP2=35,NPI=NP1+NP2,NPSI=NPS+NPI,M5=2*NPSI,
     . M1=NPSI,M2=M1+1,M4=2*M2,M1I=NPI,M2I=M1I+1,M4I=2*M2I,
     . JMAX=4,LCHMX2=2+3*JMAX,LCHMX2C=2+4*JMAX)
      DOUBLE PRECISION FMM1,MN,MP,MM
      PARAMETER (FMM1=197.327D0,
     &           MN=939.5653D0,MP=938.2720D0,MM=2.0D0*MP*MN/(MP+MN))
      LOGICAL KOPPEL,ELOG,SHELCON,ALPHAD,LOG
      INTEGER BS1(M2),BS2(M4),BZ1(M2),BZ2(M4)
      DIMENSION T1(M2I,M2I),T2(M4I,M4I)
      DIMENSION TT1(M2I,M2I),RESI(3,NPI,NPI),PSS(NPSI),GPS(NPS),A(NPS)
      DOUBLE PRECISION KMAT0(4,NPI),V0(4,NPS)
      DOUBLE COMPLEX TT1,Z,Z1,Z2,Z3,Z4,TD1,TD2,TD3,TD4,TT2,     Z5
      COMMON/VKV/V1(M2,M2),V2(M4,M4)
      COMMON/TKT/T1I(M2,M2),T2I(M4,M4)
      COMMON/GAUSSK/ PMOM(M2),AS(M2)
      COMMON /POINTS/ PPOINT(NPS),GPOINT(NPS)
      COMMON /GAUSS2/ GPNEW(NPI),PNEW(NPI)
      COMMON/VCAL/ VV2(M5,M5,LCHMX2)
      COMMON/COF/ COF(LCHMX2)
      DIMENSION TT2(M4I,M4I)
      DOUBLE COMPLEX TM
      COMMON/ TMP0/ TM(LCHMX2C,LCHMX2C,NPI+1)
      INTEGER LEFT, RIGHT
      COMMON /INDX/ LEFT(4,LCHMX2), RIGHT(4,LCHMX2)
      EQUIVALENCE (TT1(1,1),T2I(1,1))
      EQUIVALENCE ( PSS(1),PMOM(1) ),( GPS(1),A(1) )
      CHARACTER *4 TEXT(2)
      DATA RM2/1./
      DATA TEXT/'COUP','UNCO'/
      HQM=FMM1**2/MM
      AMBDA=COF(IALFA)
      DO 1 I=1,NPS
      PSS(I)=PPOINT(I)
    1 GPS(I)=GPOINT(I)
      DO 400 I=NPS+1,NPSI
  400 PSS(I)=PNEW(I-NPS)
      ELOG=.FALSE.
      IF(E.LE.0.) ELOG=.TRUE.
      PI=3.141592654
      PIHALF=PI/2.
!     X2=41.467/RM2
      X2=HQM/RM2
      CALL UNPACK1(4,L,IS,JJ,IT,INDEX,I,I,I,I)
      IF(ABS(EIMIHO(L+IS+IT,1.0d0)+1.0d0).GE.1.D-9)
     &                                       STOP'TMTRX: PAULI RULE'
      IKAN2=2
      IF(KOPPEL) IKAN2=1
      ALPHAD=KOPPEL.AND.(L.EQ.0.AND.IS.EQ.1.AND.JJ.EQ.1.AND.IT.EQ.0)
      IF(IPOT.EQ.1) ALPHAD=L.EQ.0.AND.IS.EQ.1.AND.JJ.EQ.1.AND.IT.EQ.0
C     WRITE(6,3) L,IS,JJ,IT,TEXT(IKAN2)
cjg   WRITE (80) L,IS,JJ,IT,IKAN2
    3 FORMAT(1X,26HTMATRIX FOR L,S,J,T,COUPL=,4I5,4X,A4)
      NN1=NPSI+1
      IF(ELOG) NN1=NPSI
      NN12=2*NN1
      NN1I=NPI+1
      IF(ELOG) NN1I=NPI
      NN12I=2*NN1I
      IF(ABS(E+EBDEUT).LT.1.D-9.AND.ALPHAD) GO TO 5000
      IF(KOPPEL) GOTO 500
C*****SCHLEIFE UEBER ENERGIEN IM UNGEKOPPELTEN FALL*****
      DO 2 I=1,NPSI
      DO 2 J=1,NPSI
    2 V1(I,J)=VV2(I,J,IALFA)*AMBDA
!     P0=SQRT(RM2*ABS(E)/41.467)
      P0=SQRT(RM2*ABS(E)/HQM)
      P02=P0*P0
      PRVPL=0.5*P0*DLOG((PP3-P0)/(PP3+P0))
      PMOM(NPSI+1)=P0
      DO 10 I=1,NPSI
      LOG=I.LE.NPS
      P=PMOM(I)
      P2=P*P
      IF(ELOG) GO TO 30
      IF(LOG) THEN
      AS(I)=A(I)*P**2/(P**2-P0**2)
      ELSE
      AS(I)=0.0
      ENDIF
      GOTO(211,212),IPOT
  211 CALL POTMAT(P,P0,L,IS,JJ,PNIJM)
      GOTO 213
  212 IF( IT.EQ.1 ) THEN
      CALL POTLSJIEE3(P,P0,L,IS,JJ,PNIJM,V02,V20,V22,KOPPEL)
      ELSE
      CALL POTLSJIEE3(P,P0,L,IS,JJ,PNIJM,V02,V20,V22,KOPPEL)
      ENDIF
      PNIJM=PNIJM/(PMOM(I)*P0)
C*****FOR REID POT. CARD BEFORE SHOULD BE COMMENT CARD;
  213 V1(I,NN1)=PNIJM*AMBDA
      V1(NN1,I)=PNIJM*AMBDA
      GO TO 10
   30 IF(LOG) THEN
      IF(P02.LT.1.E-6) GO TO 31
      AS(I)=A(I)*P2/(P2+P02)
      GO TO 402
   31 AS(I)=A(I)
  402 CONTINUE
      ELSE
      AS(I)=0.0
      ENDIF
   10 CONTINUE
      IF(ELOG) GO TO 32
      GOTO(221,222),IPOT
  221 CALL POTMAT(P0,P0,L,IS,JJ,PNIJM)
      GOTO 223
  222 IF( IT.EQ.1 ) THEN
      CALL POTLSJIEE3(P0,P0,L,IS,JJ,PNIJM,V02,V20,V22,KOPPEL)
      ELSE
      CALL POTLSJIEE3(P0,P0,L,IS,JJ,PNIJM,V02,V20,V22,KOPPEL)
      ENDIF
      PNIJM=PNIJM/(P0*P0)
C*****FOR REID POT. CARD BEFORE SHOULD BE COMMENT CARD;
  223 V1(NN1,NN1)=PNIJM*AMBDA
      S=0.
      DO 11 K=1,NPS
      P=PMOM(K)
      S=S+P0**2*A(K)/(P**2-P0**2)
      IF(ABS(P0-P).LE.1.E-6) GOTO  9999
   11 CONTINUE
      AS(NN1)=-S+PRVPL
   32 CALL KMAT1(E,NN1,RM2,BS1,BZ1)
      DO 404 I=1,NPI
      PMOM(I)=PMOM(I+NPS)
      DO 404 J=1,NPI
  404 T1(I,J)=T1I(I+NPS,J+NPS)
      IF(.NOT.ELOG) THEN
      PMOM(NN1I)=PMOM(NN1)
      DO 406 I=1,NN1I
      T1(NN1I,I)=T1I(NN1,I+NPS)
  406 T1(I,NN1I)=T1I(I+NPS,NN1)
      ENDIF
      IF(ICTRL.NE.0) GO TO 1100
      IF(ELOG) GO TO 800
      Z1=CMPLX(0.,-PIHALF*PMOM(NN1I))/CMPLX(1.,PIHALF*PMOM(NN1I)
     1                                              *T1(NN1I,NN1I))
      DO 34 J=1,NN1I
      DO 34 I=1,NN1I
   34 TT1(I,J)=(T1(I,J)+Z1*T1(I,NN1I)*T1(NN1I,J))*X2
      GO TO 802
  800 DO 801 J=1,NN1I
      DO 801 I=1,NN1I
  801 TT1(I,J)=T1(I,J)*X2
  802 CONTINUE
      IF(.NOT.SHELCON) GO TO 38
C     Z1=1.-CMPLX(0.,PI*PMOM(NN1I))*TT1(NN1I,NN1I)/X2
C     EL=2.*E
C     DELTA=90.*ATAN2(AIMAG(Z1),REAL(Z1))/PI
C     WRITE(6,600) EL,TT1(NN1I,NN1I),DELTA
C 600 FORMAT(5X,24HON SHELL TMATRIX FOR E =,E12.5,4H MEV,
C    1       3H =(,E12.5,2X,E12.5,1H),5X,6HDELTA=,E12.5)
cjg   WRITE (80) (TT1(NN1I,J),J=1,NPI)
C     WRITE(6,5999)  (TT1(NN1I,J),J=1,NPI)
C5999 FORMAT(1X,5(E10.4,1X,E10.4,2X))
      GO TO 1000
   38 CONTINUE
cjg
cjg      IF(.NOT.ALPHAD) GO TO 55
cjg      DO 65 I=1,NPI
cjg      DO 65 J=1,NPI
cjg   65 TT1(I,J)=TT1(I,J)*(E+EBDEUT)
cjg   WRITE (80) ((TT1(I,J),I=1,NPI),J=1,NPI)
   55 CONTINUE
      DO IP=1,NPI+1
      TM(LEFT(1,IALFA),RIGHT(1,IALFA),IP)=TT1(NPI+1,IP)
      ENDDO
      GOTO 1000
  500 CONTINUE
C*****SCHLEIFE UEBER ENERGIEN IM GEKOPPELTEN FALL*****
      DO 22 I=1,NPSI
      DO 22 J=1,NPSI
      V2(I    ,J    )=VV2(I    ,J    ,IALFA)*AMBDA
      V2(I    ,J+NN1)=VV2(I    ,J+NPSI,IALFA)*AMBDA
      V2(I+NN1,J    )=VV2(I+NPSI,J    ,IALFA)*AMBDA
   22 V2(I+NN1,J+NN1)=VV2(I+NPSI,J+NPSI,IALFA)*AMBDA
!     P0=SQRT(RM2*ABS(E)/41.467)
      P0=SQRT(RM2*ABS(E)/HQM)
      P02=P0*P0
      PRVPL=0.5*P0*DLOG((PP3-P0)/(PP3+P0))
      PMOM(NPSI+1)=P0
      DO 19 I=1,NPSI
      LOG=I.LE.NPS
      P=PMOM(I)
      P2=P*P
      IF(ELOG) GO TO 300
      IF(LOG) THEN
      AS(I)=A(I)*P**2/(P**2-P0**2)
      ELSE
      AS(I)=0.0
      ENDIF
      IF( IT.EQ.1 ) THEN
      CALL POTLSJIEE3(P,P0,L,IS,JJ,P00,P02,P20,P22,KOPPEL)
      ELSE
      CALL POTLSJIEE3(P,P0,L,IS,JJ,P00,P02,P20,P22,KOPPEL)
      ENDIF
      P00=P00/(PMOM(I)*P0)
      P02=P02/(PMOM(I)*P0)
      P20=P20/(PMOM(I)*P0)
      P22=P22/(PMOM(I)*P0)
C*****FOR REID POT. CARDS BEFORE SHOULD BE COMMENT CARDS;
      V2(I,NN1)=P00*AMBDA
      V2(NN1,I)=P00*AMBDA
      V2(I,NN12)=P02*AMBDA
      V2(NN12,I)=P02*AMBDA
      V2(NN1,NN1+I)=P20*AMBDA
      V2(NN1+I,NN1)=P20*AMBDA
      V2(NN1+I,NN12)=P22*AMBDA
      V2(NN12,NN1+I)=P22*AMBDA
      GO TO 19
  300 IF(LOG) THEN
      IF(P02.LT.1.E-6) GO TO 310
      AS(I)=A(I)*P2/(P2+P02)
      GO TO 420
  310 AS(I)=A(I)
  420 CONTINUE
      ELSE
      AS(I)=0.0
      ENDIF
   19 CONTINUE
      IF(ELOG) GO TO 320
      IF( IT.EQ.1 ) THEN
      CALL POTLSJIEE3(P0,P0,L,IS,JJ,P00,P02,P20,P22,KOPPEL)
      ELSE
      CALL POTLSJIEE3(P0,P0,L,IS,JJ,P00,P02,P20,P22,KOPPEL)
      ENDIF
      P00=P00/(P0*P0)
      P02=P02/(P0*P0)
      P20=P20/(P0*P0)
      P22=P22/(P0*P0)
C*****FOR REID POT. CARDS BEFORE SHOULD BE COMMENT CARDS;
      V2(NN1,NN1)=P00*AMBDA
      V2(NN1,NN12)=P02*AMBDA
      V2(NN12,NN1)=P20*AMBDA
      V2(NN12,NN12)=P22*AMBDA
      S=0.
      DO 20 K=1,NPS
      P=PMOM(K)
      S=S+P0**2*A(K)/(P**2-P0**2)
      IF(ABS(P0-P).LE.1.E-6) GOTO  9999
   20 CONTINUE
      AS(NN1)=-S+PRVPL
  320 CALL KMAT2(E,NN1,NN12,RM2,BS2,BZ2)
      DO 424 I=1,NPI
      PMOM(I)=PMOM(I+NPS)
      DO 424 J=1,NPI
      T2(I      ,J      )=T2I(I+NPS     ,J+NPS     )
      T2(I      ,J+NN1I )=T2I(I+NPS     ,J+NPS+NN1 )
      T2(I+NN1I ,J      )=T2I(I+NPS+NN1 ,J+NPS     )
  424 T2(I+NN1I ,J+NN1I )=T2I(I+NPS+NN1 ,J+NPS+NN1 )
      IF(.NOT.ELOG) THEN
      PMOM(NN1I)=PMOM(NN1)
      DO 426 I=1,NN1I
      T2(NN1I   ,I      )=T2I(NN1       ,I+NPS     )
      T2(I      ,NN1I   )=T2I(I+NPS     ,NN1       )
      T2(NN12I  ,I      )=T2I(NN12      ,I+NPS     )
      T2(I+NN1I ,NN1I   )=T2I(I+NN1+NPS ,NN1       )
      T2(NN1I   ,I+NN1I )=T2I(NN1       ,I+NN1+NPS )
      T2(NN12I  ,I+NN1I )=T2I(NN12      ,I+NN1+NPS )
      T2(I      ,NN12I  )=T2I(I+NPS     ,NN12      )
  426 T2(I+NN1I ,NN12I  )=T2I(I+NN1+NPS ,NN12      )
      ENDIF
      IF(ICTRL.NE.0) GO TO 2100
      IF(ELOG) GO TO 8000
      Z=CMPLX(0.,PIHALF*PMOM(NN1I))
      Z1=Z*T2(NN1I,NN1I)
      Z4=Z*T2(NN12I,NN12I)
      Z3=Z*T2(NN12I,NN1I)
      Z2=Z*T2(NN1I,NN12I)
      Z5=1./((1.+Z1)*(1.+Z4)-Z2*Z3)
      DO 340 J=1,NN1I
      TD1=((1.+Z4)*T2(NN1I,J)-Z2*T2(NN12I,J))*Z5
      TD2=((1.+Z4)*T2(NN1I,NN1I+J)-Z2*T2(NN12I,NN1I+J))*Z5
      TD3=((1.+Z1)*T2(NN12I,J)-Z3*T2(NN1I,J))*Z5
      TD4=((1.+Z1)*T2(NN12I,NN1I+J)-Z3*T2(NN1I,NN1I+J))*Z5
      DO 340 I=1,NN1I
      TT2(I,J)=(T2(I,J)-Z*(TD1*T2(I,NN1I)+TD3*T2(I,NN12I)))*X2
      TT2(I,J+NN1I)=(T2(I,J+NN1I)-Z*(TD2*T2(I,NN1I)+TD4*T2(I,NN12I)))
     . *X2
      TT2(I+NN1I,J)=(T2(I+NN1I,J)-Z*(TD1*T2(I+NN1I,NN1I)+TD3
     1                                       *T2(I+NN1I,NN12I)))*X2
  340 TT2(I+NN1I,J+NN1I)=(T2(I+NN1I,J+NN1I)-Z*(TD2*T2(I+NN1I,NN1I)
     1                                   +TD4*T2(I+NN1I,NN12I)))*X2
      GO TO 8002
 8000 DO 8001 I=1,NN12I
      DO 8001 J=1,NN12I
 8001 TT2(I,J)=T2(I,J)*X2
 8002 CONTINUE
      IF(.NOT.SHELCON) GO TO 380
C     Z=CMPLX(0.,PMOM(NN1I)*PI)
C     Z1=1.-Z*TT2(NN1I,NN1I)/X2
C     Z2= Z*TT2(NN1I,NN12I)/X2
C     Z3=1.-Z*TT2(NN12I,NN12I)/X2
C     Z=2.*Z2/(Z1-Z3)
C     EPS=0.5*ATAN(REAL(Z))
C     Z=(Z1+Z3)*0.5
C     Z4=0.5*(Z1-Z3)/COS(2.*EPS)
C     Z1=Z+Z4
C     Z2=Z-Z4
C     DELTAM=90.*ATAN2(AIMAG(Z1),REAL(Z1))/PI
C     DELTAP=90.*ATAN2(AIMAG(Z2),REAL(Z2))/PI
C******STAPP PARAMETRIZATION:
C     DELTAM=DELTAM*PI/180.
C     DELTAP=DELTAP*PI/180.
C     EPSB=0.5*ASIN( SIN(2.*EPS)*SIN(DELTAM-DELTAP) )
C     BDEL=ASIN( TAN(2.*EPSB)/TAN(2.*EPS) )
C     DELTA2B=0.5*(DELTAM+DELTAP-BDEL)
C     DELTA1B=BDEL+DELTA2B
C     EPS=EPSB*180./PI
C     DELTAM=DELTA1B*180./PI
C     DELTAP=DELTA2B*180./PI
C     WRITE(6,601) E*2.
C 601 FORMAT(5X,24HON SHELL TMATRIX FOR EL=,E12.5,4H MEV)
C     WRITE(6,602) TT2(NN1I,NN1I)
C     WRITE(6,602) TT2(NN1I,NN12I)
C     WRITE(6,602) TT2(NN12I,NN1I)
C     WRITE(6,602) TT2(NN12I,NN12I)
C     WRITE(6,603) DELTAM,DELTAP,EPS
C 603 FORMAT(5X,18HDELTAM,DELTAP,EPS=,3(E12.5,2X))
C 602 FORMAT(5X,E12.5,2X,E12.5)
cjg   WRITE (80) (TT2(NN1I,J),J=1,NPI)
cjg   WRITE (80) (TT2(NN1I,J+NN1I),J=1,NPI)
cjg   WRITE (80) (TT2(NN12I,J   ),J=1,NPI)
cjg   WRITE (80) (TT2(NN12I,J+NN1I),J=1,NPI)
      GO TO 1000
  380 CONTINUE
cjg
cjg      IF(.NOT.ALPHAD) GO TO 550
cjg      XHELP=E+EBDEUT
cjg      DO 650 I=1,NPI
cjg      DO 650 J=1,NPI
cjg      TT2(I,J)=XHELP*TT2(I,J)
cjg      TT2(I,J+NN1I)=XHELP*TT2(I,J+NN1I)
cjg      TT2(I+NN1I,J)=XHELP*TT2(I+NN1I,J)
cjg  650 TT2(I+NN1I,J+NN1I)=XHELP*TT2(I+NN1I,J+NN1I)
cjg   WRITE (80) ((TT2(I,J),I=1,NPI),J=1,NPI)
  550 CONTINUE
      DO IP=1,NPI+1
      TM(LEFT(1,IALFA),RIGHT(1,IALFA),IP)=TT2(NPI+1,IP)
      ENDDO
      IF(IPOT.EQ.1) GOTO 18
cjg   WRITE (80) ((TT2(I,J+NN1I),I=1,NPI),J=1,NPI)
cjg   WRITE (80) ((TT2(I+NN1I,J),I=1,NPI),J=1,NPI)
cjg   WRITE (80) ((TT2(I+NN1I,J+NN1I),I=1,NPI),J=1,NPI)
      DO IP=1,NPI+1
      TM(LEFT(2,IALFA),RIGHT(2,IALFA),IP)=TT2(NPI+1,IP+NPI+1)
      ENDDO
      DO IP=1,NPI+1
      TM(LEFT(3,IALFA),RIGHT(3,IALFA),IP)=TT2(NPI+1+NPI+1,IP)
      ENDDO
      DO IP=1,NPI+1
      TM(LEFT(4,IALFA),RIGHT(4,IALFA),IP)=TT2(NPI+1+NPI+1,IP+NPI+1)
      ENDDO
   18 CONTINUE
 1111 CONTINUE
      GOTO 1000
 5000 CALL RESIDU(PNEW,RESI,EBDEUT1)
      DO 5005 I=1,NPI
      DO 5005 J=1,NPI
      TT2(I,J)=CMPLX(RESI(1,I,J),0.0)
      TT2(I,J+NN1I)=CMPLX(RESI(2,I,J),0.0)
      TT2(I+NN1I,J)=CMPLX(RESI(2,J,I),0.0)
 5005 TT2(I+NN1I,J+NN1I)=CMPLX(RESI(3,I,J),0.0)
      GO TO 550
 9999 CONTINUE
      CONTR=ABS(P-P0)
      WRITE(6,23) K,P,P0,CONTR
   23 FORMAT(4X,3HK= ,I2,2X,3HP= ,1PE12.6,2X,3HP0= ,1PE12.6,3X,
     &7HCONTR= ,1PE12.6)
      STOP 'P-P0'
 1100 DO 1115 I=1,NPS
      GOTO(231,232),IPOT
  231 CALL POTMAT(0.0d0,PPOINT(I),L,IS,JJ,P00)
      GOTO 1115
  232 IF( IT.EQ.1 ) THEN
      CALL POTLSJIEE3(0.0d0,PMOM(I),L,IS,JJ,P00,P02,P20,P22,KOPPEL)
      ELSE
      CALL POTLSJIEE3(0.0d0,PMOM(I),L,IS,JJ,P00,P02,P20,P22,KOPPEL)
      ENDIF
C*****FOR POTENTIALS OTHER THAN M.TJ. CHANGE POTLSJIEE3 FOR P OR PS =0.;
 1115 V0(1,I)=P00
      DO 1120 K=1,NPI
      GOTO(235,236),IPOT
  235 CALL POTMAT(0.0d0,PNEW(K),L,IS,JJ,P00)
      GOTO 238
  236 IF( IT.EQ.1 ) THEN
      CALL POTLSJIEE3(0.0d0,PNEW(K),L,IS,JJ,P00,P02,P20,P22,KOPPEL)
      ELSE
      CALL POTLSJIEE3(0.0d0,PNEW(K),L,IS,JJ,P00,P02,P20,P22,KOPPEL)
      ENDIF
C*****FOR POTENTIALS OTHER THAN M.TJ. CHANGE POTLSJIEE3 FOR P OR PS =0.;
  238 KMAT0(1,K)=P00
      DO 1120 I=1,NPS
 1120 KMAT0(1,K)=KMAT0(1,K)-      AS(I)*V0(1,I)*T1I(I,K+NPS)
cjg   WRITE (80) (KMAT0(1,K),K=1,NPI)
      GO TO 1000
 2100 DO 2115 I=1,NPS
      IF( IT.EQ.1 ) THEN
      CALL POTLSJIEE3(0.0d0,PPOINT(I),L,IS,JJ,P00,P02,P20,P22,KOPPEL)
      ELSE
      CALL POTLSJIEE3(0.0d0,PPOINT(I),L,IS,JJ,P00,P02,P20,P22,KOPPEL)
      ENDIF
C*****FOR POTENTIALS OTHER THAN M.TJ. CHANGE POTLSJIEE3 FOR P OR PS =0.;
      V0(1,I)=P00
      V0(2,I)=P02
      V0(3,I)=P20
 2115 V0(4,I)=P22
      DO 2120 K=1,NPI
      IF( IT.EQ.1 ) THEN
      CALL POTLSJIEE3(0.0d0,PNEW(K),L,IS,JJ,P00,P02,P20,P22,KOPPEL)
      ELSE
      CALL POTLSJIEE3(0.0d0,PNEW(K),L,IS,JJ,P00,P02,P20,P22,KOPPEL)
      ENDIF
C*****FOR POTENTIALS OTHER THAN M.TJ. CHANGE POTLSJIEE3 FOR P OR PS =0.;
      KMAT0(1,K)=P00
C
C*****NEXT CHANGING IS DUE TO THE FACT THAT K(L,LS)(P,0.)/0.**LS SHOULD
C     BE WRITTEN TO TAPE 80;
C
      KMAT0(3,K)=P02
      KMAT0(2,K)=P20
      KMAT0(4,K)=P22
      DO 2120 I=1,NPS
      KMAT0(1,K)=KMAT0(1,K)-     AS(I)*(V0(1,I)*T2I(I    ,K+NPS)+
     1                                  V0(2,I)*T2I(I+NN1,K+NPS))
      KMAT0(3,K)=KMAT0(3,K)-     AS(I)*(V0(2,I)*T2I(I+NN1,K+NN1+NPS)+
     1                                  V0(1,I)*T2I(I    ,K+NN1+NPS))
      KMAT0(2,K)=KMAT0(2,K)-     AS(I)*(V0(3,I)*T2I(I    ,K+NPS)+
     1                                  V0(4,I)*T2I(I+NN1,K+NPS))
 2120 KMAT0(4,K)=KMAT0(4,K)-     AS(I)*(V0(3,I)*T2I(I    ,K+NN1+NPS)+
     1                                  V0(4,I)*T2I(I+NN1,K+NN1+NPS))
cjg   WRITE (80) KMAT0
 1000 CONTINUE
      RETURN
      END
      SUBROUTINE CALV(IPOT)
C
C*****CALCULATE NEEDED V(P,PS)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (NPS=72,NP1=10,NP2=35,NPI=NP1+NP2,
     .           NPSI=NPS+NPI,M5=2*NPSI,
     .           JMAX=4,LCHMX2=2+3*JMAX)
      LOGICAL KOPPEL,COUPL
      DIMENSION PMOM(NPSI)
      COMMON /POINTS/ PPOINT(NPS),GPOINT(NPS)
      COMMON /GAUSS2/ GPNEW(NPI),PNEW(NPI)
      COMMON/VCAL/ VV2(M5,M5,LCHMX2)
      COMMON/PWTB2/ INDEX2(LCHMX2),NINDEX2,COUPL(LCHMX2)
      DO 1 I=1,NPS
    1 PMOM(I)=PPOINT(I)
      DO 2 I=NPS+1,NPSI
    2 PMOM(I)=PNEW(I-NPS)
      DO 1000 IALFA=1,NINDEX2
      KOPPEL=COUPL(IALFA)
      CALL UNPACK1(4,L,IS,JJ,IT,INDEX2(IALFA),I,I,I,I)
      IF(ABS(EIMIHO(L+IS+IT,1.0d0)+1.0d0).GE.1.D-9)
     &                               STOP'TMTRX: PAULI RULE'
      IF(KOPPEL) GOTO 500
      DO 5 I=1,NPSI
      P=PMOM(I)
      DO 4 J=1,I
      PS=PMOM(J)
      GOTO(201,202),IPOT
  201 CALL POTMAT(P,PS,L,IS,JJ,PNIJM)
      GOTO 4
  202 IF( IT.EQ.1 ) THEN
      CALL POTLSJIEE3(P,PS,L,IS,JJ,PNIJM,V02,V20,V22,KOPPEL)
      ELSE
      CALL POTLSJIEE3(P,PS,L,IS,JJ,PNIJM,V02,V20,V22,KOPPEL)
      ENDIF
      PNIJM=PNIJM/(PMOM(I)*PMOM(J))
C*****FOR REID POT. CARD BEFORE SHOULD BE COMMENT CARD;
    4 VV2(I,J,IALFA)=PNIJM
    5 CONTINUE
      DO 7 I=1,NPSI
      DO 6 J=1,I
    6 VV2(J,I,IALFA)=VV2(I,J,IALFA)
    7 CONTINUE
      GOTO 1000
  500 CONTINUE
      DO 14 I=1,NPSI
      P=PMOM(I)
      DO 13 J=1,I
      PS=PMOM(J)
      IF( IT.EQ.1 ) THEN
      CALL POTLSJIEE3(P,PS,L,IS,JJ,P00,P02,P20,P22,KOPPEL)
      ELSE
      CALL POTLSJIEE3(P,PS,L,IS,JJ,P00,P02,P20,P22,KOPPEL)
      ENDIF
      P00=P00/(PMOM(I)*PMOM(J))
      P02=P02/(PMOM(I)*PMOM(J))
      P20=P20/(PMOM(I)*PMOM(J))
      P22=P22/(PMOM(I)*PMOM(J))
C*****FOR REID POT. CARDS BEFORE SHOULD BE COMMENT CARDS;
      VV2(I    ,J    ,IALFA)=P00
      VV2(I    ,J+NPSI,IALFA)=P02
      VV2(I+NPSI,J    ,IALFA)=P20
   13 VV2(I+NPSI,J+NPSI,IALFA)=P22
   14 CONTINUE
      DO 16 I=1,NPSI
      DO 15 J=1,I
      VV2(J    ,I    ,IALFA)=VV2(I    ,J    ,IALFA)
      VV2(J+NPSI,I    ,IALFA)=VV2(I    ,J+NPSI,IALFA)
      VV2(J    ,I+NPSI,IALFA)=VV2(I+NPSI,J    ,IALFA)
   15 VV2(J+NPSI,I+NPSI,IALFA)=VV2(I+NPSI,J+NPSI,IALFA)
   16 CONTINUE
 1000 CONTINUE
      RETURN
      END
      SUBROUTINE POTMAT(P,PS,L,IS,JJ,PNIJM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
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
      PNIJM= -0.25*HELPA*DLOG( (AMI2+X1)/(AMI2+X2) )
      PNIJM=PNIJM+0.25*HELPR*DLOG( (RMI2+X1)/(RMI2+X2) )
      PNIJM=PNIJM/(P*PS)
      GO TO 30
   20 PNIJM=-HELPA/(AMI2+X1)+HELPR/(RMI2+X1)
   30 PNIJM=PNIJM/41.47
      RETURN
      END
      SUBROUTINE RESIDU(PMOM,RESI,E0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (NPS=72,NDEUT=100,NP1=10,NP2=35,NPI=NP1+NP2)
      DIMENSION PD(NDEUT),SD(NDEUT),DD(NDEUT),PMOM(NPI),SDEUT(NPI),
     1          DDEUT(NPI),SPL(NDEUT),RESI(3,NPI,NPI)
      DOUBLE PRECISION MN,MP,MM,M,FMM1,HQM
      PARAMETER (FMM1=197.327D0,
     &           MN=939.5653D0,MP=938.2720D0,MM=2.0D0*MP*MN/(MP+MN))
      HQM=FMM1**2/MM
      READ(14,*) PD,E0
      READ(14,*) SD
      READ(14,*) DD
      CLOSE (14,ERR=20)
      GO TO 30
   20 STOP 'RESIDU: ERROR IN CLOSING DEUWP1'
   30 CONTINUE
      DO 1234 IDEUT=1,NDEUT
 1234 DD(IDEUT)=DD(IDEUT)/PD(IDEUT)**2
      CALL SPLFAK(PD,NDEUT)
      DO 40 IP=1,NPI
      CALL SPLMOD(PD,NDEUT,PMOM(IP),SPL)
      SDEUT(IP)=0.0
      DDEUT(IP)=0.0
      IF(PMOM(IP).GT.PD(NDEUT)) GO TO 40
      DO 50 IPD=1,NDEUT
      SDEUT(IP)=SDEUT(IP)+SPL(IPD)*SD(IPD)
      DDEUT(IP)=DDEUT(IP)+SPL(IPD)*DD(IPD)
   50 CONTINUE
   40 CONTINUE
      DO 1235 IPI=1,NPI
 1235 DDEUT(IPI)=DDEUT(IPI)*PMOM(IPI)**2
      DO 6005 I=1,NPI
      DO 6005 J=1,NPI
!     X1=(-E0+41.467*PMOM(I)**2)*SDEUT(I)
!     X2=(-E0+41.467*PMOM(I)**2)*DDEUT(I)
!     X3=(-E0+41.467*PMOM(J)**2)*SDEUT(J)
!     X4=(-E0+41.467*PMOM(J)**2)*DDEUT(J)
      X1=(-E0 + HQM*PMOM(I)**2)*SDEUT(I)
      X2=(-E0 + HQM*PMOM(I)**2)*DDEUT(I)
      X3=(-E0 + HQM*PMOM(J)**2)*SDEUT(J)
      X4=(-E0 + HQM*PMOM(J)**2)*DDEUT(J)
      RESI(1,I,J)=X1*X3
      RESI(2,I,J)=X1*X4
      RESI(3,I,J)=X2*X4
 6005 CONTINUE
      RETURN
      END
      SUBROUTINE KMAT1(E,NN1,RM2,BS,BZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER(NPS=72,NP1=10,NP2=35,NPI=NP1+NP2,
     .          NPSI=NPS+NPI,M1=NPSI,M2=M1+1,M4=2*M2)
      INTEGER BS(M2),BZ(M2)
      COMMON/VKV/V1(M2,M2),V2(M4,M4)
      COMMON/TKT/T1(M2,M2),T2(M4,M4)
      DIMENSION H1(M2,M2),H2(M4,M4)
      COMMON/GAUSSK/PMOM(M2),AS(M2)
      DO 1 I=1,NN1
      DO 1 J=1,NN1
    1 H1(I,J)=RM2*V1(I,J)*AS(J)
      DO 2 I=1,NN1
    2 H1(I,I)=H1(I,I)+1.
      CALL INVERT(H1,NN1,1.D-6,1,KRIT,BS,BZ,M2)
      DO 3 I=1,NN1
      DO 4 J=1,NN1
      T1(I,J)=0.
      DO 5 K=1,NN1
    5 T1(I,J)=H1(I,K)*V1(K,J)+T1(I,J)
    4 CONTINUE
    3 CONTINUE
      RETURN
      END
      SUBROUTINE KMAT2(E,NN1,NN12,RM2,BS,BZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER(NPS=72,NP1=10,NP2=35,NPI=NP1+NP2,NPSI=NPS+NPI,
     .          M1=NPSI,M2=M1+1,M4=2*M2)
      INTEGER BS(M4),BZ(M4)
      COMMON/VKV/V1(M2,M2),V2(M4,M4)
      COMMON/TKT/T1(M2,M2),T2(M4,M4)
      DIMENSION H1(M2,M2),H2(M4,M4)
      COMMON/GAUSSK/ PMOM(M2),AS(M2)
      DO 1 I=1,NN1
      DO 2 J=1,NN1
      H2(I,J)=RM2*V2(I,J)*AS(J)
      H2(I+NN1,J+NN1)=RM2*V2(I+NN1,J+NN1)*AS(J)
      H2(I,J+NN1)=RM2*V2(I,J+NN1)*AS(J)
    2 H2(I+NN1,J)=RM2*V2(I+NN1,J)*AS(J)
    1 CONTINUE
      DO 3 I=1,NN12
    3 H2(I,I)=H2(I,I)+1.
      CALL INVERT(H2,NN12,1.D-6 ,1,KRIT,BS,BZ,M4)
      DO 4 I=1,NN12
      DO 5 J=1,NN12
      T2(I,J)=0.
      DO 6 K=1,NN12
    6 T2(I,J)=H2(I,K)*V2(K,J)+T2(I,J)
    5 CONTINUE
    4 CONTINUE
      RETURN
      END
      SUBROUTINE ALFA2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (JMAX=4,LCHMX2=2+3*JMAX)
      LOGICAL COUPL
      COMMON /PWTB2/ INDEX2(LCHMX2),NINDEX2,COUPL(LCHMX2)
      CHARACTER*1 LC(0:13)
      DATA LC /'S','P','D','F','G','H','I','J','K','L','M','N','O','P'/
      NINDEX2=0
      WRITE (*,100)
  100 FORMAT(//,' CHANNELS USED FOR TWO BODY INTERACTION'//
     .5X,'L',4X,'S',4X,'J',4X,'T',4X,'N',4X,'INDEX(N)'/)
      DO 10 J=0,JMAX
C     DO 10 J=0,1
       DO 10 IS=0,1
        ISS=2*IS+1
        IF (J. EQ. 0) THEN
         JE=IS
        ELSE
         JE=J
        ENDIF
        DO 10 L=ABS(J-IS),JE
C       IF( .NOT.((L.EQ.0.AND.IS.EQ.0.AND.J.EQ.0).OR.
C    .            (L.EQ.0.AND.IS.EQ.1.AND.J.EQ.1) ) ) GOTO 10
C        IF (.NOT. (L.EQ.1.AND.IS.EQ.1.AND.
C     .     (J.EQ.0.OR.J.EQ.1.OR.J.EQ.2)) )
C     .          GOTO 10
C        IF( L.EQ.1.AND.IS.EQ.0.AND.J.EQ.1) GOTO 10
C        IF( J.EQ.2.AND.(.NOT.(L.EQ.1.AND.IS.EQ.1.AND.J.EQ.2))) GOTO 10
         NINDEX2=NINDEX2+1
         IT=((-1)**(L+IS)+1)/2
         CALL PACK1(4,L,IS,J,IT,INDEX2(NINDEX2),IZ,IZ,IZ,IZ)
         CALL UNPACK1(4,L,IS,J,IT,INDEX2(NINDEX2),IZ,IZ,IZ,IZ)
         WRITE (*,101) L,IS,J,IT,NINDEX2,INDEX2(NINDEX2)
  101    FORMAT(1X,5I5,I11)
         IF (L .EQ. J-1) THEN
          COUPL(NINDEX2)=.TRUE.
          WRITE (*,102) 'COUP',ISS,LC(L),J,'-',ISS,LC(L+2),J
  102     FORMAT ('+',39X,A4,3X,I2,A1,I2,A1,I2,A1,I2)
         ELSE
          COUPL(NINDEX2)=.FALSE.
          WRITE (*,102) 'UNCO',ISS,LC(L),J
         END IF
   10    CONTINUE
      END

      FUNCTION EIMIHO(K,FAC)
C     ===============
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      EIMIHO=FAC
      IF( (K/2) * 2 .NE. K) EIMIHO=-FAC
      RETURN
      END
      SUBROUTINE INVERT (A,N,EPS,MOD,KRIT,S,Z,M)
C     =================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      INTEGER S(N),Z(N)
      DIMENSION A(M,M)
C*****M=MATR.-DIM. IM AUFRUFENDEN PROGRAMM,N=MATR.-DIM. IN INVERT*****
C*****VORBESETZEN BUCHHALTER*****
      DO 101 I=1,N
      Z(I)=-I
  101 S(I)=I
      MS=0
      MZ=0
      KRIT=0
      DO 999 L=1,N
C*****AUFSUCHEN PIVOT*****
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
C*****AUSTAUSCHSCHRITT*****
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
C*****ZEILEN UND SPALTENTAUSCH*****
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
C*****BUCHHALTERTAUSCH*******
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
      SUBROUTINE INTEGRAL_P0
C     =============
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (JMAX=4,LCHMX2=2+3*JMAX,LCHMX2C=2+4*JMAX,
     &           PBAR=30.0,NP1=10,NP2=35,NPI=NP1+NP2)
      DOUBLE PRECISION FMM1,MN,MP,MM
      PARAMETER (FMM1=197.327D0,
     &           MN=939.5653D0,MP=938.2720D0,MM=2.0D0*MP*MN/(MP+MN))
      DOUBLE COMPLEX TM
      COMMON/ TMP0/ TM(LCHMX2C,LCHMX2C,NPI+1)
      DOUBLE COMPLEX ZWWW
      COMMON /GAUSS/ AP(NPI+1),P(NPI+1)
      LOGICAL COUPLC
      COMMON /PWTB2C/ INDEX2C(LCHMX2C),NINDEX2C,COUPLC(LCHMX2C)
      COMMON /KINEM/ OMEGA,QBIG,ELABNN,P0,THEQ
      COMMON /CURRENT/ CHARGE(NPI+1,LCHMX2C,-1:1),
     &                   PLUS(NPI+1,LCHMX2C,-1:1)
      DOUBLE COMPLEX CHARGET, PLUST
      COMMON /RESCATT/ CHARGET(LCHMX2C,-1:1),
     &                   PLUST(LCHMX2C,-1:1)
      DIMENSION WW(NPI)
c
      AMNMEV=MM ! MeV
      HQM=FMM1**2/AMNMEV
      write(6,*) p0, p(npi+1)
      PI=DACOS(-1.0d0)
      DO IP=1,NPI
      WW(IP)=1./HQM*AP(IP)/(P(NPI+1)**2-P(IP)**2)
      ENDDO
      WWW1= 1./HQM*0.5*DLOG((PBAR+P0)/(PBAR-P0))/P0
      WWW2=-1./HQM*0.5*PI/P0
      ZWWW=CMPLX(WWW1,WWW2)
C
      DO IAL=1,NINDEX2C
      DO MD=-1,1
      CHARGET(IAL,MD)=(0.,0.)
      PLUST(IAL,MD)=(0.,0.)
C
      DO IALS=1,NINDEX2C
      DO IP=1,NPI
      CHARGET(IAL,MD)=CHARGET(IAL,MD) +
     & WW(IP)*( P(IP)**2*TM(IAL,IALS,IP)*CHARGE(IP,IALS,MD)
     &       - P0**2*TM(IAL,IALS,NPI+1)*CHARGE(NPI+1,IALS,MD) )
      PLUST(IAL,MD)=PLUST(IAL,MD) +
     & WW(IP)*( P(IP)**2*TM(IAL,IALS,IP)*PLUS(IP,IALS,MD)
     &       - P0**2*TM(IAL,IALS,NPI+1)*PLUS(NPI+1,IALS,MD) )
      ENDDO
      CHARGET(IAL,MD)=CHARGET(IAL,MD) +
     &       ZWWW*P0**2*TM(IAL,IALS,NPI+1)*CHARGE(NPI+1,IALS,MD)
      PLUST(IAL,MD)=PLUST(IAL,MD) +
     &       ZWWW*P0**2*TM(IAL,IALS,NPI+1)*PLUS(NPI+1,IALS,MD)
      ENDDO
C
      ENDDO
      ENDDO
C
      RETURN
      END
      SUBROUTINE RESCAMP
C     =============
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (NTHMAX=181,FMM1=197.327)
      PARAMETER (NP1=10,NP2=35,NP=NP1+NP2+1)
      PARAMETER (JMAX=4,LCHMX2C=2+4*JMAX,LMAX=JMAX+1)
      COMMON /KINEM/ OMEGA,QBIG,ELABNN,P0,THEQ
      DOUBLE COMPLEX CHARGET, PLUST
      COMMON /RESCATT/ CHARGET(LCHMX2C,-1:1),
     &                   PLUST(LCHMX2C,-1:1)
      LOGICAL COUPLC
      COMMON /PWTB2C/ INDX2C(LCHMX2C),NINDX2C,COUPLC(LCHMX2C)
      DOUBLE COMPLEX YLMP, CCC
      COMMON /YLM/ YLMP(NTHMAX,0:LMAX,-LMAX:LMAX)
      DOUBLE COMPLEX TRAMPRHO2, TRAMPPLUS2, TRAMPMINUS2
c     COMMON /TRAMPFU/ TRAMPRHO2(ith,im1,im2,md) etc.
      COMMON /TRAMPFU/ TRAMPRHO2(NTHMAX,2,2,-1:1),
     &                TRAMPPLUS2(NTHMAX,2,2,-1:1),
     &               TRAMPMINUS2(NTHMAX,2,2,-1:1)
      COMMON /DETECT/ MTTWE1,MTTWE2
C
C**** Spherical harmonics for all angles p
c     CALL YKQFN
c     write(6,*) ' YLMP ', YLMP
C
C**** TRAMPRHO2
c     do md=-1,1
c     do im2=1,2
c     do im1=1,2
c     do ith=1,NTHMAX
c     TRAMPRHO(ith,im1,im2,md)=(0.,0.)
c     enddo
c     enddo
c     enddo
c     enddo
      DO 1 IALS=1,NINDX2C
      CALL UNPACK1(4,LS,ISS,JS,ITS,INDX2C(IALS),IZ,IZ,IZ,IZ)
      LS2=2*LS
      ISS2=2*ISS
      JS2=2*JS
      ITS2=2*ITS
      DO 2 MD=-1,1
      MD2=2*MD
      MP=MD  ! m' for rho
      IF (ABS(MP).GT.JS) GOTO 2
      MP2=2*MP
      DO 3 im2=1,2
      M22=2*IM2-3
      DO 4 im1=1,2
      M12=2*IM1-3
      MUP2=MP2-M12-M22
      MUP=MUP2/2
      IF (ABS(MUP2).GT.LS2) GOTO 4
      IF (ABS(MP2-MUP2).GT.ISS2) GOTO 4
      CCC=CLEBSCH(1,1,ISS2,M12,M22,M12+M22)
     &   *CLEBSCH(1,1,ITS2,MTTWE1,MTTWE2,MTTWE1+MTTWE2)
     &   *CLEBSCH(LS2,ISS2,JS2,MUP2,MP2-MUP2,MP2)
     &   *CHARGET(IALS,MD)
      DO 5 ITH=1,NTHMAX
      TRAMPRHO2(ith,im1,im2,md)=TRAMPRHO2(ith,im1,im2,md)
     &  + CCC*YLMP(ITH,LS,MUP)
    5 CONTINUE
    4 CONTINUE
    3 CONTINUE
    2 CONTINUE
    1 CONTINUE
C
C**** TRAMPPLUS2
c     do md=-1,1
c     do im2=1,2
c     do im1=1,2
c     do ith=1,NTHMAX
c     TRAMPPLUS(ith,im1,im2,md)=(0.,0.)
c     enddo
c     enddo
c     enddo
c     enddo
      DO 11 IALS=1,NINDX2C
      CALL UNPACK1(4,LS,ISS,JS,ITS,INDX2C(IALS),IZ,IZ,IZ,IZ)
      LS2=2*LS
      ISS2=2*ISS
      JS2=2*JS
      ITS2=2*ITS
      DO 13 MD=-1,1
      MD2=2*MD
      MP=MD+1  ! m' for j_+
      IF (ABS(MP).GT.JS) GOTO 13
      MP2=2*MP
      DO 14 im2=1,2
      M22=2*IM2-3
      DO 15 im1=1,2
      M12=2*IM1-3
      MUP2=MP2-M12-M22
      MUP=MUP2/2
      IF (ABS(MP-MUP).GT.ISS) GOTO 15
      CCC=CLEBSCH(1,1,ISS2,M12,M22,MP2-MUP2)
     &   *CLEBSCH(1,1,ITS2,MTTWE1,MTTWE2,0)
     &   *CLEBSCH(LS2,ISS2,JS2,MUP2,MP2-MUP2,MP2)
     &   *PLUST(IALS,MD)
      DO 16 ITH=1,NTHMAX
      TRAMPPLUS2(ith,im1,im2,md)=TRAMPPLUS2(ith,im1,im2,md)
     &  + CCC*YLMP(ITH,LS,MUP)
   16 CONTINUE
   15 CONTINUE
   14 CONTINUE
   13 CONTINUE
   11 CONTINUE
C
C**** TRAMPMINUS2
c     do md=-1,1
c     do im2=1,2
c     do im1=1,2
c     do ith=1,NTHMAX
c     TRAMPMINUS(ith,im1,im2,md)=(0.,0.)
c     enddo
c     enddo
c     enddo
c     enddo
      DO 111 IALS=1,NINDX2C
      CALL UNPACK1(4,LS,ISS,JS,ITS,INDX2C(IALS),IZ,IZ,IZ,IZ)
      LS2=2*LS
      ISS2=2*ISS
      JS2=2*JS
      ITS2=2*ITS
      DO 113 MD=-1,1
      MD2=2*MD
      MP=MD-1  ! m'
      IF (ABS(MP).GT.JS) GOTO 113
      MP2=2*MP
      DO 114 im2=1,2
      M22=2*IM2-3
      DO 115 im1=1,2
      M12=2*IM1-3
      MUP2=MP2-M12-M22
      MUP=MUP2/2
      IF (ABS(MP-MUP).GT.ISS) GOTO 115
      CCC=CLEBSCH(1,1,ISS2,M12,M22,MP2-MUP2)
     &   *CLEBSCH(1,1,ITS2,MTTWE1,MTTWE2,0)
     &   *CLEBSCH(LS2,ISS2,JS2,MUP2,MP2-MUP2,MP2)
     &   *PLUST(IALS,-MD)
     &   *(-1.)**(LS+JS-1)  !  symmetry properties are used !!!!
      DO 116 ITH=1,NTHMAX
      TRAMPMINUS2(ith,im1,im2,md)=TRAMPMINUS2(ith,im1,im2,md)
     &  + CCC*YLMP(ITH,LS,MUP)
  116 CONTINUE
  115 CONTINUE
  114 CONTINUE
  113 CONTINUE
  111 CONTINUE
C
c     write(6,*) ' TRAMPRHO ', TRAMPRHO
c     write(6,*) ' TRAMPPLUS ', TRAMPPLUS
c     write(6,*) ' TRAMPMINUS ', TRAMPMINUS
c     do ith=1,nthmax
c     write(6,333) ith,
c    &                real(TRAMPRHO(ith,1,1,-1)),
c    &                real(TRAMPRHO(ith,1,1, 0)),
c    &                real(TRAMPRHO(ith,1,1, 1)),
c    &                real(TRAMPRHO(ith,1,2,-1)),
c    &                real(TRAMPRHO(ith,1,2, 0)),
c    &                real(TRAMPRHO(ith,1,2, 1)),
c    &                real(TRAMPRHO(ith,2,1,-1)),
c    &                real(TRAMPRHO(ith,2,1, 0)),
c    &                real(TRAMPRHO(ith,2,1, 1)),
c    &                real(TRAMPRHO(ith,2,2,-1)),
c    &                real(TRAMPRHO(ith,2,2, 0)),
c    &                real(TRAMPRHO(ith,2,2, 1))
c     enddo
c 333 format(i5,12E10.2)
c     do ith=1,nthmax
c     rr1=0.0
c     rr2=0.0
c     do md=-1,1
c     do im2=1,2
c     do im1=1,2
c     rr1=rr1 + cabs(TRAMPRHO(ith,im1,im2,md))**2
c     rr2=rr2 + cabs(TRAMPPLUS(ith,im1,im2,md))**2
c    &        + cabs(TRAMPMINUS(ith,im1,im2,md))**2
c     enddo
c     enddo
c     enddo
c     write(6,*) ith, rr1,rr2
c     enddo
C
C
C
      RETURN
      END

      SUBROUTINE ALFA2CJGBS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (JMAXBS=10,LCHMX2CBS=2)
      LOGICAL COUPL
      CHARACTER*4 TEXT
      COMMON /PWTB2CBS/ INDEX2(LCHMX2CBS),NINDEX2,COUPL(LCHMX2CBS)
C
      NINDEX2=0
      WRITE(6,100)
  100 FORMAT(//,
     &   '    Channels(C JG b.s.)  used for two-body interaction'/
     &5X,'L',4X,'S',4X,'J',4X,'T',4X,'N',4X,'INDEX(N)')
c     DO 10 J=0,JMAXBS
c      DO 11 IS=0,1
c       DO 11 L=ABS(J-IS),J+IS
      DO 10 J=1,1
       DO 11 IS=1,1
        DO 11 L=0,2,2
         NINDEX2=NINDEX2+1
         IT=((-1)**(L+IS)+1)/2
         CALL PACK1(4,L,IS,J,IT,INDEX2(NINDEX2),IZ,IZ,IZ,IZ)
         CALL UNPACK1(4,L1,IS1,J1,IT1,INDEX2(NINDEX2),IZ,IZ,IZ,IZ)
         IF (  L.NE.J .AND. IS.EQ.1 .AND. J.GT.0  ) THEN
       TEXT='COUP'
         COUPL(NINDEX2)=.TRUE.
         WRITE(6,101) L1,IS1,J1,IT1,NINDEX2,INDEX2(NINDEX2),TEXT
  101    FORMAT(1X,5I5,I11,1X,A4)
         ELSE
       TEXT='UNCO'
         COUPL(NINDEX2)=.FALSE.
         WRITE(6,101) L1,IS1,J1,IT1,NINDEX2,INDEX2(NINDEX2),TEXT
         END IF
   11    CONTINUE
CC         WRITE(6,*)
   10   CONTINUE
      RETURN
      END
