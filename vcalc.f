      subroutine vcalc(r,theta,chi,r0,rmax,rstep,nrgrd,nthgrd,nchgrd,
     +                 vtable,vpes)
      implicit double precision(a-h,o-z)

      dimension vtable(0:nrgrd*nthgrd*nchgrd-1)
      parameter(wno2k=0.6950356d0,wno2kj=1.196266d-2)
c     parameter(maxrpt=nrgrd-1,mxthpt=180,mxchpt=180)
c     character*30 argum

c ... the subroutine takes r / bohr, theta / deg, and chi / deg

      maxrpt=nrgrd-1
      mxthpt=nthgrd-1
      mxchpt=nchgrd-1

      if(r.lt.r0) r=r0
      if(r.gt.rmax)r=rmax

      ir=int((r-r0)/rstep)
      ith=int(theta)
      ich=int(chi)
      if(ich.gt.mxchpt)ich=mxchpt
      if(ich.lt.0)ich=0
      if(ith.gt.mxthpt)ith=mxthpt
      if(ith.lt.0)ith=0
c     write(6,*)ir,ith,ich

      indrtc=(ir*nthgrd+ith)*nchgrd+ich
      v0=vtable(indrtc)
c     write(6,*)v0,v0*wno2k

c ... linear interpolation
      if(ir.eq.maxrpt) then
        gradr=0.d0
        delr=0.d0
      else
        indrtc=((ir+1)*nthgrd+ith)*nchgrd+ich
        gradr=(vtable(indrtc)-v0)/rstep
        delr=r-(r0+ir*rstep)
      endif

      if(ith.eq.mxthpt) then
        gradth=0.d0
        delth=0.d0
      else
        indrtc=(ir*nthgrd+(ith+1))*nchgrd+ich
        gradth=(vtable(indrtc)-v0)
        delth=theta-dfloat(ith)
      endif

      if(ich.eq.mxchpt) then
        gradch=0.d0
        delch=0.d0
      else
        indrtc=(ir*nthgrd+ith)*nchgrd+ich+1
        gradch=(vtable(indrtc)-v0)
        delch=chi-dfloat(ich)
      endif

      vpes=v0+gradr*delr+gradth*delth+gradch*delch

c     write(6,*)vpes,vpes*wno2k,vpes*wno2k*wno2kj

      return
      end
