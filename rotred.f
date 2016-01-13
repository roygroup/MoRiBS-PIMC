      subroutine rotred(rhoprp,erotpr)
      implicit double precision(a-h,o-z)

      dimension rhoprp(0:23588100),erotpr(0:23588100)

      open(2,file='rho.den_rho')
      open(3,file='rho.den_eng')

      do il=0,23588100
        read(2,*)rhoprp(il)
        read(3,*)erotpr(il)
      enddo

      close(2,status='keep')
      close(3,status='keep')

      return
      end
