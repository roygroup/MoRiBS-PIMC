      subroutine vcord(Eulang,RCOM,RpH2,vtable,nrgrd,nthgrd,nchgrd,
     +                 rvmax,rvmin,rvstep,vpot,
     +                 radret,theret,chiret,hatx,haty,hatz,ivcord)
c ... this program reads the coordinates of centre of mass of water and pH2,
c ... and three Euler angles of water, and returns the direction unit vectors
c ... for the three axes in water fixed frame and R, theta, chi for potential
c ... calculation.  ivcord determines whether the subroutine is run till the end.
c ... For only calculating hatx, haty, and hatz, the code does NOT need to go to the end.
c ... ivcord=0 goes to the end. ivcord=1 returns after hat(x,y,z) are calculated.
c ... NB: The angles are in RADIA, the coordinates are in ANGSTROM, and the potnetial are in KELVIN
      implicit double precision(a-h,o-z)

      dimension rotmat(3,3),
     +          hatz(3),haty(3),hatx(3),RCOM(3),RpH2(3),
     +          Eulang(3),RH2COM(3),unx(3),
     +          uny(3),unz(3),origin(3)
      dimension vtable(0:nrgrd*nthgrd*nchgrd-1)
      parameter (pi=3.14159265358979323846d0,zero=0.d0,small=1.d-08,
     +           bo2ang=0.529177249d0)
      data 
     +     unx/1.d0,zero,zero/,
     +     uny/zero,1.d0,zero/,unz/zero,zero,1.d0/,
     +     origin/zero,zero,zero/

c ... Eulang: three Euler angles, 1: phi, 2: theta: 3: chi
c ... rotmat: rotational matrix of the three Euler angles
c ... all lengths are in units of Angstrom, and all angles in radian

c     write(6,*)'in vcord',nrgrd,nthgrd,nchgrd,rvmax,rvmin,rvstep

      call matpre(Eulang,rotmat)

c ... obtain WFF orientation directly from rotation
      call rottrn(rotmat,unx,hatx,origin)
      call rottrn(rotmat,uny,haty,origin)
      call rottrn(rotmat,unz,hatz,origin)


      do i=1,3
        RH2COM(i)=RpH2(i)-RCOM(i)
      enddo


      if(ivcord.eq.1) return

c ... calculate thewff and chiwff
      thewff=dotang(hatz,RH2COM)
c     write(6,*)thewff*180d0/Pi
c     write(6,*)(hatx(i),i=1,3),(rh2com(i),i=1,3),dotprd(RH2COM,hatx)
      if(abs(dotprd(RH2COM,hatx)).lt.small) then
        chiwff=pi/2.d0
        goto 101
      endif
      tanchi=dotprd(RH2COM,haty)/dotprd(RH2COM,hatx)
c ... chiwff maps the chi into the first quadrant to have the C2v
c ... symmetrically unique chi
      chiwff=atan(abs(tanchi))
  101 continue
c     write(6,*)chiwff*180.d0/pi
c     write(6,*)chiwff*180d0/Pi
      radwff=dnorm(RH2COM)
c     write(6,*)radwff,radwff/bo2ang
c ... prepare r (Angs), theta (Radian), and chi (Radian) that are returned for density binning
      radret=radwff
      theret=thewff
      Rdotx=dotprd(RH2COM,hatx)
      Rdoty=dotprd(RH2COM,haty)
      if(Rdotx.ge.0.d0.and.Rdoty.ge.0.d0) then
        chiret=chiwff
      elseif(Rdotx.lt.0.d0.and.Rdoty.ge.0.d0) then
        chiret=Pi-chiwff
      elseif(Rdotx.lt.0.d0.and.Rdoty.lt.0.d0) then
        chiret=Pi+chiwff
      else
        chiret=2*Pi-chiwff
      endif
c ... get chiwff in the right range determined by nchgrd
      if(nchgrd.eq.91) then
c       do nothing because chiwff is in the first quadran already
      elseif(nchgrd.eq.181) then
        chiwff=chiret
        if(chiret.gt.Pi) chiwff=2*pi-chiret
      elseif(nchgrd.eq.361) then
        chiwff=chiret
      endif
c ... convert radwff to the unit of bohr
      radwff=radwff/bo2ang
c ... convert thewff and chiwff to degree
      thewff=thewff*180.d0/pi
      chiwff=chiwff*180.d0/pi
c     write(6,*)chiwff
c     write(6,*)radwff,thewff,chiwff,chiret*180.d0/pi

      call vcalc(radwff,thewff,chiwff,rvmin,rvmax,rvstep,nrgrd,nthgrd,
     +           nchgrd,vtable,vpot)
c     write(6,*)vpot

      end
c-----------------------------------------------------------------------
c     subroutine matpre(Eulang,rotmat)
c     implicit double precision(a-h,o-z)

c     dimension Eulang(3),rotmat(3,3)

c     phi=Eulang(1)
c     theta=Eulang(2)
c     chi=Eulang(3)

c     cp=cos(phi)
c     sp=sin(phi)
c     ct=cos(theta)
c     st=sin(theta)
c     ck=cos(chi)
c     sk=sin(chi)

c     rotmat(1,1)=cp*ct*ck-sp*sk
c     rotmat(1,2)=-cp*ct*sk-sp*ck
c     rotmat(1,3)=cp*st
c     rotmat(2,1)=sp*ct*ck+cp*sk
c     rotmat(2,2)=-sp*ct*sk+cp*ck
c     rotmat(2,3)=sp*st
c     rotmat(3,1)=-st*ck
c     rotmat(3,2)=st*sk
c     rotmat(3,3)=ct

c     return
c     end
c-----------------------------------------------------------------------
c     subroutine rottrn(rotmat,rwf,rsf,rcom)
c     implicit double precision(a-h,o-z)
c     dimension rotmat(3,3),rwf(3),rsf(3),rcom(3)

c     do i=1,3
c       rsf(i)=rcom(i)
c       do j=1,3
c         rsf(i)=rsf(i)+rotmat(i,j)*rwf(j)
c       enddo
c     enddo

c     return
c     end
c-----------------------------------------------------------------------
      subroutine crsprd(vec1,vec2,vec3)
      implicit double precision(a-h,o-z)
      dimension vec1(3),vec2(3),vec3(3)
c ... vec3 = vec1 cross vec2
      vec3(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
      vec3(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
      vec3(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)

      return
      end
c-----------------------------------------------------------------------
      double precision function dotprd(vec1,vec2)
      implicit double precision(a-h,o-z)
      dimension vec1(3),vec2(3)
c ... dot=vec1 dot vec2
      dotprd=0.d0
      do i=1,3
        dotprd=dotprd+vec1(i)*vec2(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      double precision function dnorm(vec)
      implicit double precision(a-h,o-z)
      dimension vec(3)

      dnorm=dotprd(vec,vec)
      dnorm=sqrt(dnorm)

      return
      end
c-----------------------------------------------------------------------
      double precision function dotang(vec1,vec2)
      implicit double precision(a-h,o-z)
      dimension vec1(3),vec2(3)
c ... returns the cross angle between two vectors in radian

      dotang=dotprd(vec1,vec2)/(dnorm(vec1)*dnorm(vec2))
c ... the following two lines eliminate the out of range case
      if(dotang.gt.1.d0) dotang=1.d0
      if(dotang.lt.-1.d0)dotang=-1.d0
      dotang=acos(dotang)

      return
      end
c-----------------------------------------------------------------------
      subroutine reflec(coord,rcom,hatx,haty,hatz)
      implicit double precision(a-h,o-z)

      dimension coord(3),rcom(3),hatx(3),haty(3),hatz(3),rmff(3),
     +          rpmff(3),rpsff(3),rrel(3),rsff(3),vec(3)

c     write(6,*)(coord(i),i=1,3),(rcom(i),i=1,3)
c     write(6,*)(hatx(i),i=1,3),(haty(i),i=1,3),(hatz(i),i=1,3)

c ... calculate the coordinate in MFF of HCOOCH3
      do i=1,3
        rrel(i)=coord(i)-rcom(i)
      enddo

      rmff(1)=dotprd(rrel,hatx)
      rmff(2)=dotprd(rrel,haty)
      rmff(3)=dotprd(rrel,hatz)

c ... reflect in the MFF of HCOOCH3

      dist=abs(2.0*rmff(2))
      rpmff(1)=rmff(1)
      rpmff(2)=-rmff(2)
      rpmff(3)=rmff(3)

c     dist1=0.0d0
c     do i=1,3
c       dist1=dist1+(rpmff(i)-rmff(i))**2.d0
c     enddo
c     dist1=sqrt(dist1)

c ... try converting rmff back to rsff and compare with coord
c     vec(1)=rmff(1)+dotprd(rcom,hatx)
c     vec(2)=rmff(2)+dotprd(rcom,haty)
c     vec(3)=rmff(3)+dotprd(rcom,hatz)

c     rsff(1)=vec(1)*hatx(1)+vec(2)*haty(1)+vec(3)*hatz(1)
c     rsff(2)=vec(1)*hatx(2)+vec(2)*haty(2)+vec(3)*hatz(2)
c     rsff(3)=vec(1)*hatx(3)+vec(2)*haty(3)+vec(3)*hatz(3)

c ... try converting rpmff back to rpsff
      vec(1)=rpmff(1)+dotprd(rcom,hatx)
      vec(2)=rpmff(2)+dotprd(rcom,haty)
      vec(3)=rpmff(3)+dotprd(rcom,hatz)

      rpsff(1)=vec(1)*hatx(1)+vec(2)*haty(1)+vec(3)*hatz(1)
      rpsff(2)=vec(1)*hatx(2)+vec(2)*haty(2)+vec(3)*hatz(2)
      rpsff(3)=vec(1)*hatx(3)+vec(2)*haty(3)+vec(3)*hatz(3)

      dist2=0.d0
      do i=1,3
        dist2=dist2+(rpsff(i)-coord(i))**2.d0
      enddo
      dist2=sqrt(dist2)

      if(abs(dist2-dist).gt.1.d-6) then
        write(6,*)(rpsff(i),i=1,3),(coord(i),i=1,3),dist,dist2
      endif

c .... replace the original coordinate
       do i=1,3
         coord(i)=rpsff(i)
       enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rflmfy(rcom,hatx,haty,hatz,eulang)
c ... reflect MF molecule wrt the SFF XZ plane and reflect wrt the MFF xz plane
c     implicit double precision(a-h,o-z)
      implicit none
      double precision zero,theta2,sint,sphi,small,schi,rotma2,pi,phi2,
     +                 hatzp,hatxp,cphi,cost,chi2,cchi,eulang,hatz,haty,
     +                 hatx,rcom,small1,small2

      integer i

      dimension rcom(3),hatx(3),haty(3),hatz(3),eulang(3),hatxp(3),
     +          hatzp(3),rotma2(3,3)
      parameter(pi=3.14159265358979323846d0,zero=0.d0,small=1.d-08,
     +          small1=1.d-02,small2=1.d-08)

c     write(6,*)'in F'
c     do i=1,3
c       write(6,*)rcom(i),hatx(i),haty(i),hatz(i),eulang(i)
c     enddo

c ... reflect MFF x and z
      hatx(2)=-hatx(2)
      hatz(2)=-hatz(2)

c ... calculate MFF y by cross product to maintain the right-handed coordinate
c     write(6,*)'before haty',(haty(i),i=1,3)
      call crsprd(hatz,hatx,haty)
c     write(6,*)'after haty',(haty(i),i=1,3)

c ... make the rotational matrix from the reflected x,y,z
      do i=1,3
        rotma2(i,1)=hatx(i)
        rotma2(i,2)=haty(i)
        rotma2(i,3)=hatz(i)
      enddo

c ... borrow code from deleul in rotden.f
      cost=rotma2(3,3)
      call within(cost)
      theta2=acos(cost)

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

c ... return angles in radia
      Eulang(1)=phi2
      Eulang(2)=theta2
      Eulang(3)=chi2
c     write(6,*)'in F:',Eulang(1),cos(Eulang(2)),Eulang(3)

      return
      end
c-----------------------------------------------------------------------
      subroutine rflmfx(rcom,hatx,haty,hatz,eulang)
c ... reflect MF molecule wrt the SFF XZ plane and reflect wrt the MFF yz plane
c     implicit double precision(a-h,o-z)
      implicit none
      double precision zero,theta2,sint,sphi,small,schi,rotma2,pi,phi2,
     +                 hatzp,hatxp,cphi,cost,chi2,cchi,eulang,hatz,haty,
     +                 hatx,rcom,small1,small2

      integer i

      dimension rcom(3),hatx(3),haty(3),hatz(3),eulang(3),hatxp(3),
     +          hatzp(3),rotma2(3,3)
      parameter(pi=3.14159265358979323846d0,zero=0.d0,small=1.d-08,
     +          small1=1.d-02,small2=1.d-08)

c     write(6,*)'in F'
c     do i=1,3
c       write(6,*)rcom(i),hatx(i),haty(i),hatz(i),eulang(i)
c     enddo

c ... reflect MFF y and z
      haty(2)=-haty(2)
      hatz(2)=-hatz(2)

c ... calculate MFF y by cross product to maintain the right-handed coordinate
c     write(6,*)'before haty',(haty(i),i=1,3)
      call crsprd(haty,hatz,hatx)
c     write(6,*)'after haty',(haty(i),i=1,3)

c ... make the rotational matrix from the reflected x,y,z
      do i=1,3
        rotma2(i,1)=hatx(i)
        rotma2(i,2)=haty(i)
        rotma2(i,3)=hatz(i)
      enddo

c ... borrow code from deleul in rotden.f
      cost=rotma2(3,3)
      call within(cost)
      theta2=acos(cost)

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

c ... return angles in radia
      Eulang(1)=phi2
      Eulang(2)=theta2
      Eulang(3)=chi2
c     write(6,*)'in F:',Eulang(1),cos(Eulang(2)),Eulang(3)

      return
      end
c-----------------------------------------------------------------------
      subroutine rflmfz(rcom,hatx,haty,hatz,eulang)
c ... reflect MF molecule wrt the SFF XZ plane and reflect wrt the MFF yz plane
c     implicit double precision(a-h,o-z)
      implicit none
      double precision zero,theta2,sint,sphi,small,schi,rotma2,pi,phi2,
     +                 hatzp,hatxp,cphi,cost,chi2,cchi,eulang,hatz,haty,
     +                 hatx,rcom,small1,small2

      integer i

      dimension rcom(3),hatx(3),haty(3),hatz(3),eulang(3),hatxp(3),
     +          hatzp(3),rotma2(3,3)
      parameter(pi=3.14159265358979323846d0,zero=0.d0,small=1.d-08,
     +          small1=1.d-02,small2=1.d-08)

c     write(6,*)'in F'
c     do i=1,3
c       write(6,*)rcom(i),hatx(i),haty(i),hatz(i),eulang(i)
c     enddo

c ... reflect MFF y and z
      hatx(2)=-hatx(2)
      haty(2)=-haty(2)

c ... calculate MFF y by cross product to maintain the right-handed coordinate
c     write(6,*)'before haty',(haty(i),i=1,3)
      call crsprd(hatx,haty,hatz)
c     write(6,*)'after haty',(haty(i),i=1,3)

c ... make the rotational matrix from the reflected x,y,z
      do i=1,3
        rotma2(i,1)=hatx(i)
        rotma2(i,2)=haty(i)
        rotma2(i,3)=hatz(i)
      enddo

c ... borrow code from deleul in rotden.f
      cost=rotma2(3,3)
      call within(cost)
      theta2=acos(cost)

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

c ... return angles in radia
      Eulang(1)=phi2
      Eulang(2)=theta2
      Eulang(3)=chi2
c     write(6,*)'in F:',Eulang(1),cos(Eulang(2)),Eulang(3)

      return
      end
