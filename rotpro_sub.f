      subroutine rotpro(chi,phi,theta,rho,erot,esq,rhoprp,erotpr,erotsq,
     +                  jstop)
      implicit double precision(a-h,o-z)
      parameter(Pi=3.141592653589793238462643383279502884197d0)
      dimension rhoprp(0:23588100),erotpr(0:23588100),erotsq(0:23588100)

      ichi=int(chi)
      iphi=int(phi)
      itheta=int(theta)
      if(ichi.gt.360.or.ichi.lt.0) then
        write(6,*)'ichi out or range',ichi,chi
        ichi=0
        jstop=1
      endif
      if(iphi.gt.360.or.iphi.lt.0) then
        write(6,*)'iphi out or range',iphi,phi
        iphi=0
        jstop=1
      endif
      if(itheta.gt.180.or.itheta.lt.0) then
        write(6,*)'itheta out or range',itheta,theta
        itheta=0
        jstop=1
      endif

      indtpc=(itheta*361+iphi)*361+ichi
      rho0=rhoprp(indtpc)
      erot0=erotpr(indtpc)
      esq0=erotsq(indtpc)

c ... first to calculate propagator
      delchi=0.d0
      delphi=0.d0
      delthe=0.d0
      if(ichi.ne.360) then
c       delchi=rhoprp(ichi+1,iphi,itheta)-rhoprp(ichi,iphi,itheta)
        indtpc=(itheta*361+iphi)*361+ichi+1
        delchi=rhoprp(indtpc)-rho0
        delch2=erotpr(indtpc)-erot0
        delch3=erotsq(indtpc)-esq0
      endif
      if(iphi.ne.360) then
c       delphi=rhoprp(ichi,iphi+1,itheta)-rhoprp(ichi,iphi,itheta)
        indtpc=(itheta*361+iphi+1)*361+ichi
        delphi=rhoprp(indtpc)-rho0
        delph2=erotpr(indtpc)-erot0
        delph3=erotsq(indtpc)-esq0
      endif
      if(itheta.ne.180) then
c       delthe=rhoprp(ichi,iphi,itheta+1)-rhoprp(ichi,iphi,itheta)
        indtpc=((itheta+1)*361+iphi)*361+ichi
        delthe=rhoprp(indtpc)-rho0
        delth2=erotpr(indtpc)-erot0
        delth3=erotsq(indtpc)-esq0
      endif
      rho=rho0+delchi*(chi-dfloat(ichi))+
     +    delphi*(phi-dfloat(iphi))+delthe*(theta-dfloat(itheta))
      erot=erot0+delch2*(chi-dfloat(ichi))+
     +    delph2*(phi-dfloat(iphi))+delth2*(theta-dfloat(itheta))
      esq=esq0+delch3*(chi-dfloat(ichi))+
     +    delph3*(phi-dfloat(iphi))+delth3*(theta-dfloat(itheta))
      return

      end
