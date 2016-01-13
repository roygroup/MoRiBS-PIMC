      subroutine potred(vtable)
      implicit double precision(a-h,o-z)

      dimension vtable(0:501*181*181-1)

      open(2,file='pes.tab',status='old')

      do i=0,501*181*181-1
        read(2,*)vtable(i)
      enddo

      close(2,status='keep')

      return
      end
