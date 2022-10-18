      subroutine potlsjiee3(p,ps,l,is,jmom,v00,v02,v20,v22,koup)
      use SMSCHIRAL
      implicit none
      real p,ps,v00,v02,v20,v22
      integer l,is,jmom
      logical koup
      real(kind(0.0D0)) potent,qpun1,qpun2,FACTOR,Pi
      dimension POTENT(6)
      REAL MN,MP,MM,M,FMM1
      PARAMETER (FMM1=197.327D0,
     &           MN=939.5653D0,MP=938.2720D0,MM=2.0D0*MP*MN/(MP+MN))
      character(len=2) FORCE
      integer OSTAT,CUTNUM
      common /pot/ OSTAT,CUTNUM,FORCE
      Pi=DACOS(-1.0D0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  M=938.926/FMM1  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      M=MM/FMM1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      QPUN1=p*FMM1/1000.0
      QPUN2=ps*FMM1/1000.0
c     FACTOR=p*ps*FMM1*938.919*2./Pi
c     FACTOR=FACTOR*Pi/2.0D0/1000000.0D0
      FACTOR=FMM1**2/1000000.0*p*ps*M
      call CHIRALMOMPWD(OSTAT,FORCE,QPUN1,QPUN2,JMOM,CUTNUM,POTENT)
      if(koup) goto 100
!
      if(l.eq.jmom .and. is.eq.0) v00=potent(1)*FACTOR
      if(l.eq.jmom .and. is.eq.1) v00=potent(2)*FACTOR

!     3P0

      if(l.eq.1.and.is.eq.1.and.jmom.eq.0) v00=potent(6)*FACTOR
      return
  100 continue
      v00=potent(3)*FACTOR
      v02=potent(4)*FACTOR
      v20=potent(5)*FACTOR
      v22=potent(6)*FACTOR
      return
      end
