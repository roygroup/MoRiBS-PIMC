      subroutine initconf(coords,angles,ntotal,nboson,indexp,indexr)
      implicit double precision(a-h,o-z)
      character argum*30

      dimension coords(0:3*ntotal-1),angles(0:3*ntotal-1),
     +          indexp(0:nboson-1),indexr(0:nboson-1)

c     write(6,*)'nboson=',nboson,(indexp(i),i=0,nboson-1)
c     write(6,*)'nboson=',nboson,(indexr(i),i=0,nboson-1)
      open(2,file='xyz.init',status='old')
      read(2,*)idump,(indexp(i),i=0,nboson-1)
      do i=0,nboson-1
        indexr(indexp(i))=i
      enddo
c     write(6,*)'idump=',idump,'ntotal=',ntotal
      read(2,*)argum
c     write(6,*)argum
      do i=0,ntotal-1
        read(2,*)argum,(coords(i*3+j),angles(i*3+j),j=0,2) 
c       write(6,'(A8,I7,6F12.6)')
c    +    argum,i,(coords(i*3+j),angles(i*3+j),j=0,2)
      enddo
      close(2,status='keep')


      return
      end
      subroutine prtper(indexp,nboson,nblock)
      implicit double precision(a-h,o-z)

      dimension indexp(0:nboson-1)
      character string*80
      integer*8 nblock

      open(2,file='permutation.tab',status='unknown')
   10 continue
      read(2,'(A80)',end=20)string
      goto 10
   20 continue
      write(2,'(A,I4,A,80I3)')
     +     'BLOCK',nblock,':',(indexp(i),i=0,nboson-1)

      close(2,status='keep')

      return
      end
