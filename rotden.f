      subroutine rotden(Eulan1,Eulan2,Eulrel,rho,erot,esq,rhoprp,erotpr,
     +                  erotsq,istop)
      implicit double precision(a-h,o-z)
      dimension Eulan1(3),Eulan2(3),Eulrel(3),rhoprp(0:23588100),
     +          erotpr(0:23588100),erotsq(0:23588100)
      parameter(pi=3.14159265358979323846d0,wno2k=0.6950356d0)

c     write(6,*)(Eulan1(i),Eulan2(i),i=1,3)
c     write(6,'(6f10.6)')(Eulan1(i),i=1,3)
c     write(6,'(6f10.6)')(Eulan2(i),i=1,3)
      call deleul(Eulan1,Eulan2,Eulrel,istop)

c ... Eulrel is in radian
c     write(6,*)(Eulrel(i)*180.d0/Pi,i=1,3)
      phi=Eulrel(1)*180.d0/Pi
      theta=Eulrel(2)*180.d0/Pi
      chi=Eulrel(3)*180.d0/Pi

      jstop=0
      call rotpro(chi,phi,theta,rho,erot,esq,rhoprp,erotpr,erotsq,jstop)
c     write(6,'(a,f10.5,a,f10.5)')'rho=',rho,' erot=',erot
c ... convert rotational energy from cm-1 to kelvin
      if(jstop.eq.1) then
        write(6,*)(Eulan1(i),i=1,3),(Eulan2(i),i=1,3),(Eulrel(i),i=1,3),
     +            phi,theta,chi
        istop=1
      endif
      erot=erot/wno2k
      esq=esq/(wno2k*wno2k)

      end
      subroutine deleul(Eulan1,Eulan2,Eulrel,istop)
c ... this program reads in two sets of Euler angles and return the set of
c ... Euler angles between the two WFF
c ... All Euler angles in Radia
      implicit double precision(a-h,o-z)

      dimension rotmat(3,3),
     +          Eulan1(3),Eulan2(3),
     +          rotinv(3,3),Eulrel(3),rotma1(3,3),
     +          rotmar(3,3),rotma2(3,3)
      parameter(pi=3.14159265358979323846d0,zero=0.d0,small=1.d-08,
     +          small1=1.d-02,small2=1.d-08)

      istop=0


c ... prepare rotational matricies
      call matpre(Eulan2,rotmat)
      call matpre(Eulan1,rotma1)

c ... directly extract relative euler angles from multiplied matrix
      do i=1,3
        do j=1,3
          rotma2(i,j)=0.d0
          do k=1,3
            rotma2(i,j)=rotma2(i,j)+rotma1(k,i)*rotmat(k,j)
          enddo
        enddo
      enddo

      cost=rotma2(3,3)
      call within(cost)
      theta2=acos(cost)

c     if(abs(theta-theta2).gt.small) then
c       write(6,*)'Warning for theta',theta,theta2
c     endif

      sint=sin(theta2)
      if(abs(1.d0-cost).lt.small) then
c ... theta=0
        phi2=0.d0
        cchi=rotma2(1,1)
        schi=rotma2(2,1)
        call within(cchi)
        call within(schi)
        if(schi.gt.zero) then
          chi2=acos(cchi)
        else
          chi2=2.0*Pi-acos(cchi)
        endif
      elseif(abs(1.d0+cost).lt.small) then
c ... theta=pi
        phi2=0.d0
        cchi=rotma2(2,2)
        schi=rotma2(1,2)
        call within(cchi)
        call within(schi)
        if(schi.gt.zero) then
          chi2=acos(cchi)
        else
          chi2=2.0*Pi-acos(cchi)
        endif
      else
c ... normal theta
        cphi=rotma2(1,3)/sint
        sphi=rotma2(2,3)/sint
        cchi=-rotma2(3,1)/sint
        schi=rotma2(3,2)/sint
        call within(cphi)
        call within(sphi)
        call within(cchi)
        call within(schi)
        if(sphi.gt.zero) then
          phi2=acos(cphi)
        else
          phi2=2.0*Pi-acos(cphi)
        endif
        if(schi.gt.zero) then
          chi2=acos(cchi)
        else
          chi2=2.0*Pi-acos(cchi)
        endif
      endif

      Eulrel(1)=phi2
      Eulrel(2)=theta2
      Eulrel(3)=chi2

c     if(abs(phi-phi2).gt.10.*small) then
c       write(6,*)'Warning for phi',phi,phi2
c     endif
c     if(abs(chi-chi2).gt.10.*small) then
c       write(6,*)'Warning for chi',chi,chi2
c     endif

c     itest=itest+1

c     if(itest.lt.ntest) goto 100

      return

      end
c-----------------------------------------------------------------------
      subroutine matpre(Eulang,rotmat)
      implicit double precision(a-h,o-z)

      dimension Eulang(3),rotmat(3,3)

      phi=Eulang(1)
      theta=Eulang(2)
      chi=Eulang(3)

      cp=cos(phi)
      sp=sin(phi)
      ct=cos(theta)
      st=sin(theta)
      ck=cos(chi)
      sk=sin(chi)

      rotmat(1,1)=cp*ct*ck-sp*sk
      rotmat(1,2)=-cp*ct*sk-sp*ck
      rotmat(1,3)=cp*st
      rotmat(2,1)=sp*ct*ck+cp*sk
      rotmat(2,2)=-sp*ct*sk+cp*ck
      rotmat(2,3)=sp*st
      rotmat(3,1)=-st*ck
      rotmat(3,2)=st*sk
      rotmat(3,3)=ct

      return
      end
c-----------------------------------------------------------------------
      subroutine rottrn(rotmat,rwf,rsf,rcom)
      implicit double precision(a-h,o-z)
      dimension rotmat(3,3),rwf(3),rsf(3),rcom(3)

      do i=1,3
        rsf(i)=rcom(i)
        do j=1,3
          rsf(i)=rsf(i)+rotmat(i,j)*rwf(j)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
c     subroutine rottrn(rotmat,rwf,rsf,rcom)
c     implicit double precision(a-h,o-z)
c     parameter(m=3,ip=3,n=1,alpha=1.d0,beta=1.d0)
c     dimension rotmat(3,3),rwf(3),rsf(3),rcom(3)

c ... use dgemm
c     do i=1,3
c       rsf(i)=rcom(i)
c     enddo
c     call dgemm('N','N',m,n,ip,alpha,rotmat,m,rwf,ip,beta,rsf,m)


c     return
c     end
c-----------------------------------------------------------------------
      subroutine invert(rotmat,rotinv)
c ... rotinv = rotmat transpose
      implicit double precision(a-h,o-z)
      dimension rotmat(3,3),rotinv(3,3)

      do i=1,3
        do j=1,3
          rotinv(i,j)=rotmat(j,i)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine within(value)
c ... this subroutine chops off the possible miscerllaneous value beyond
c ... 1 or -1
      implicit double precision(a-h,o-z)
      if(value.gt.1.d0)value=1.d0
      if(value.lt.-1.d0)value=-1.d0

      return
      end
c-----------------------------------------------------------------------
      subroutine rsrot(Eulan1,Eulan2,xrot,yrot,zrot,tauC,iodevn,eoff,
     +                  rho,erot)
      implicit double precision(a-h,o-z)
c ... calculate the rotational propagator and energy estimator using 
c ... RATTLE and SHAKE formula

      parameter(omass=15.994915d0,hmass=1.007826d0,zero=0.d0)
      parameter(wn2ha=219474.6313705d0,amu2me=1822.88839d0,
     +          a2bohr=1.8897162d0,wno2k=0.6950356d0)
      dimension Eulan1(3),Eulan2(3),rotma1(3,3),rotma2(3,3),ROwf(3),
     +          R1wf(3),R2wf(3),vec(3),rcom(3),rotm1t(3,3),uec(3),
     +          digrel(3),blist(3)
c ... WFF coordinates in Angs. They are derived from A0, B0, and C0
c ... and should not be misused as re.
c .. the following geometry is close to Valiron's geometry in his potential paper
      data ROwf/zero,zero,0.06562d0/,R1wf/0.7557d0,zero,-0.5223d0/,
     +     R2wf/-0.7557d0,zero,-0.5223d0/,rcom/0.d0,0.d0,0.d0/
c     data ROwf/zero,zero,0.065043d0/,
c    +     R1wf/0.7588977d0,zero,-0.516139d0/,
c    +     R2wf/-0.7588977d0,zero,-0.516139d0/,rcom/0.d0,0.d0,0.d0/
c ... the above H2O structure give is derived from A and B. It gives
c ... C=9.5483 cm-1, quite different from the C used in the Noya formula

c ... the intake tau is in the unit of K-1 and needs to be converted to cm
c ... tauC has the same address as MCRotTau in the C code and its value shouldn't
c ... be modified.

c     real*16 rho1

      tau=tauC/wno2k

c ... put the rotational constants in the list
      blist(1)=1.d0/xrot
      blist(2)=1.d0/yrot
      blist(3)=1.d0/zrot

c ... debug
c     write(6,'(7f10.5)')(Eulan1(i),i=1,3),(Eulan2(i),i=1,3),tau

      call matpre(Eulan1,rotma1)
      call matpre(Eulan2,rotma2)

c ... calculate the diagonal matrix elements for the rotational matrix of the relative euler angles
      do i=1,3
        digrel(i)=0.d0
        do j=1,3
          digrel(i)=digrel(i)+rotma1(j,i)*rotma2(j,i)
        enddo
      enddo

      sumaxs=(blist(1)-blist(2)-blist(3))*(1.d0-digrel(1))+
     +       (blist(2)-blist(3)-blist(1))*(1.d0-digrel(2))+
     +       (blist(3)-blist(1)-blist(2))*(1.d0-digrel(3))

      if(iodevn.eq.-1) then
c       rho=exp(sumaxs/(4.d0*tau))
c       rho=sumaxs/(4.d0*tau)
c       rhotemp=sumaxs/(4.d0*tau)
c ...   return rho without being divided by 4tau. the division will be done in the c-code
        rho=sumaxs
c       write(6,*)'in rsrot:',rhotemp,rho,4.d0*tau
c       rho=exp(rho)
        erot=sumaxs/(4.d0*tau*tau)
        erot=erot+1.5d0/tau+0.25d0*(xrot+yrot+zrot)
c ...   convert to K
        erot=erot/wno2k
      endif

      return

c ... loop over atoms to get the differential rotation term
c ... convert to atomic unit
      tauhar=tau*wn2ha
      sumatm=zero
      call rottrn(rotma2,R1wf,vec,rcom)
      call rottrn(rotma1,R1wf,uec,rcom)
      term11=hmass*dotprd(uec,vec)
      call rottrn(rotma1,R2wf,uec,rcom)
      term21=hmass*dotprd(uec,vec)
      term1=hmass*dotprd(R1wf,R1wf)

      call rottrn(rotma2,R2wf,vec,rcom)
      call rottrn(rotma1,R2wf,uec,rcom)
      term22=hmass*dotprd(uec,vec)
      call rottrn(rotma1,R1wf,uec,rcom)
      term12=hmass*dotprd(uec,vec)
      term2=hmass*dotprd(R2wf,R2wf)

      call rottrn(rotma2,ROwf,vec,rcom)
      call rottrn(rotma1,ROwf,uec,rcom)
      termoo=omass*dotprd(uec,vec)
      termo=omass*dotprd(ROwf,ROwf)

      if(iodevn.eq.-1) then
        sumatm=term11+term22+termoo-term1-term2-termo
        sumatm=sumatm*amu2me*a2bohr*a2bohr/tauhar
        rho=exp(sumatm)
        sumatm=sumatm*wn2ha/tauhar
        erot=sumatm
c       write(6,*)'erot=',erot
c ...   use 3/2(1/tau) as the energy reference
        erot=1.5d0/tau+erot
      else
        trmH12=term11+term22
        trmH21=term12+term21
c       trmH12=trmH12*amu2me*a2bohr*a2bohr/tauhar
c       trmH21=trmH21*amu2me*a2bohr*a2bohr/tauhar
c       trmxch=(exp(trmH12)+((-1)**iodevn)*exp(trmH21))*0.5d0
        sumatm=termoo+trmH12-term1-term2-termo
        sumat2=termoo+trmH21-term1-term2-termo
        sumatm=sumatm*amu2me*a2bohr*a2bohr/tauhar
        sumat2=sumat2*amu2me*a2bohr*a2bohr/tauhar
        rho=exp(sumatm)+((-1)**iodevn)*exp(sumat2)

        exp1=trmH12*amu2me*a2bohr*a2bohr/tauhar
        exp2=trmH21*amu2me*a2bohr*a2bohr/tauhar
        exp1=exp(exp1)
        exp2=exp(exp2)
        upper=exp1*trmH12+((-1)**iodevn)*exp2*trmH21
        denomi=exp1+((-1)**iodevn)*exp2
        xerot=upper/denomi
        erot=xerot+termoo-term1-term2-termo
        erot=erot*amu2me*a2bohr*a2bohr/tauhar
        erot=erot*wn2ha/tauhar
c ...   use 3/2(1/tau) as the energy reference
        erot=1.5d0/tau+erot
      endif

c ... offset erot by the read-in Noya energy, which has the unit of cm-1
      erot=eoff+erot
c ... convert erot from cm-1 to kelvin
      erot=erot/wno2k
c ... debug
c     write(6,'(2f10.5)')rho,erot
      return

      end
c-----------------------------------------------------------------------
      subroutine rsline(Brot,dprd,tauC,rho,erot)
      implicit double precision(a-h,o-z)
c ... rattle and shake propagator and energy estimator for linear rotor
      parameter(wn2ha=219474.6313705d0,amu2me=1822.88839d0,
     +          a2bohr=1.8897162d0,wno2k=0.6950356d0)
      parameter(pi=3.14159265358979323846d0)

      tau=tauC/wno2k

      rho=(1.d0-dprd)/(2.d0*brot*tau)
      erot=(1.d0-rho)/tau
c     rho=exp(-rho)/(4.d0*pi)
c ... do the exponent in the c++
      rho=-rho

c ... convert erot from cm-1 to kelvin
      erot=erot/wno2k

      return
      end
