c     ==================================================
      subroutine caleng(com_1, com_2, E_2H2O,
                        Eulang_1, Eulang_2,)
c     ==================================================
c     Initially:
c     __________________________________________________
c     this subroutine calculates TIP4P potential between
c     two rigid waters given the coordinates
c     of their centres of mass and 
c     their respective Euler angles
c     e: com1, com2, rotmat1, rotmat2
c        ROwf, R1wf, R2wf, RMwf
c     s: E2H2O
c     __________________________________________________
c     GG (Dec. 16th 2011):
c     __________________________________________________
c     rotmat1, rotmat2 computed within the code with 
c     Eulang1, Eulang2
c     ROwf, R1wf, R2wf, RMwf put as data   
c     e: com_1, com_2, Eulang_1, Eulang_2
c     s: E_2H2O 
      implicit double precision(a-h,o-z)
      parameter(zero=0.d0)
      dimension ROwf(3), RH1wf(3), RH2wf(3), RMwf(3),
     +          com_1(3), com_2(3), Eulang_1(3),
     +          Eulang_2(3), RO_1_sf(3), RO_2_sf(3),
     +          RH1_1_sf(3), RH1_2_sf(3),
     +          RH2_1_sf(3), RH2_2_sf(3),
     +          RM_1_sf(3), RM_2_sf(3), vec(3),
     +          rotmat_2(3,3), rotmat_1(3,3), 
     +          crossA(3), crossB(3)
c     TIP4P parameters (L-J) & conversion factors: 
      parameter(epsoo=0.154875717017208413d0,sigoo=3.15365d0,
     +          qo=-1.040d0,qh=0.520d0,
     +          br2ang=0.52917721092d0,hr2kcl=627.509469d0)
      data ROwf/zero,zero,0.06562d0/,RH1wf/0.7557d0,zero,-0.5223d0/,
     +     RH2wf/-0.7557d0,zero,-0.5223d0/
     +     RMwf/0.d0,0.d0,ROwf(3)-.15d0/
c
c
      call matpre(Eulang_1, rotmat_1)
c      do i=1,3
c         RO1sf(i)=0.d0
c      enddo
c      call DGEMV ('N', 3, 3, 1.d0, rotmat_1, 3, ROwf, 1, 1.d0, RO_1_sf, 1 )
      call rottrn(rotmat_1, ROwf, RO_1_sf, com1)

c      do i=1,3
c         R11sf(i)=0.d0
c      enddo
c      call DGEMV ('N', 3, 3, 1.d0, rotmat_1, 3, R1wf, 1, 1.d0, R1_1_sf, 1 )
      call rottrn(rotmat_1, RH1wf, RH1_1_sf, com1)
c
c      do i=1,3
c         R21sf(i)=0.d0
c      enddo
c      call DGEMV ('N', 3, 3, 1.d0, rotmat1, 3, R2wf, 1, 1.d0, R21sf, 1 )
      call rottrn(rotmat_1, RH2wf, RH2_1_sf, com1)
c
      do i=1,3
         RM1sf(i)=0.d0
      enddo
      call DGEMV ('N', 3, 3, 1.d0, rotmat1, 3, RMwf, 1, 1.d0, RM1sf, 1 )
c
c ... prepare rotational matrix for water 2
c ... obtain the SFF coordinates for H, H, and O of water 2
c     call rottrn(rotmat2,ROwf,RO2sf,com2)
      do i=1,3
         RO2sf(i)=com2(i)
      enddo
      call DGEMV ('N', 3, 3, 1.d0, rotmat2, 3, ROwf, 1, 1.d0, RO2sf, 1 )
c
c     call rottrn(rotmat2,R1wf,R12sf,com2)
      do i=1,3
         R12sf(i)=com2(i)
      enddo
      call DGEMV ('N', 3, 3, 1.d0, rotmat2, 3, R1wf, 1, 1.d0, R12sf, 1 )
c
c     call rottrn(rotmat2,R2wf,R22sf,com2)
      do i=1,3
         R22sf(i)=com2(i)
      enddo
      call DGEMV ('N', 3, 3, 1.d0, rotmat2, 3, R2wf, 1, 1.d0, R22sf, 1 )
c
c     call rottrn(rotmat2,RMwf,RM2sf,com2)
      do i=1,3
         RM2sf(i)=com2(i)
      enddo
      call DGEMV ('N', 3, 3, 1.d0, rotmat2, 3, RMwf, 1, 1.d0, RM2sf, 1 )
c
c ... calculate water dimer energies through SPC/WF formula
      E2H2O=0.d0
c ... O-O interaction
      roo=0.0d0
      rMM=0.0d0
      do i=1,3
        roo=roo+(RO1sf(i)-RO2sf(i))*(RO1sf(i)-RO2sf(i))
        rMM=rMM+(RM1sf(i)-RM2sf(i))*(RM1sf(i)-RM2sf(i))
      enddo
      rMM=sqrt(rMM)
      roo=sqrt(roo)
      roo4=roo*roo
      roo6=roo4*roo
      roo12=roo6*roo6
c     AMBER values 
      A=5.99896595E+05
      B=6.09865468E+02
      o2lj=A/roo12-B/roo6
c ... H-O, H-H and O-O Columbic interaction
      rho1=0.0
      rho2=0.0
      rho3=0.0
      rho4=0.0
      rhh1=0.0
      rhh2=0.0
      rhh3=0.0
      rhh4=0.0
      do i=1,3
        rho1=rho1+(RM1sf(i)-R12sf(i))*(RM1sf(i)-R12sf(i))
        rho2=rho2+(RM1sf(i)-R22sf(i))*(RM1sf(i)-R22sf(i))
        rho3=rho3+(RM2sf(i)-R11sf(i))*(RM2sf(i)-R11sf(i))
        rho4=rho4+(RM2sf(i)-R21sf(i))*(RM2sf(i)-R21sf(i))
        rhh1=rhh1+(R11sf(i)-R12sf(i))*(R11sf(i)-R12sf(i))
        rhh2=rhh2+(R11sf(i)-R22sf(i))*(R11sf(i)-R22sf(i))
        rhh3=rhh3+(R21sf(i)-R12sf(i))*(R21sf(i)-R12sf(i))
        rhh4=rhh4+(R21sf(i)-R22sf(i))*(R21sf(i)-R22sf(i))
      enddo
      rho1=sqrt(rho1)
      rho2=sqrt(rho2)
      rho3=sqrt(rho3)
      rho4=sqrt(rho4)
      rhh1=sqrt(rhh1)
      rhh2=sqrt(rhh2)
      rhh3=sqrt(rhh3)
      rhh4=sqrt(rhh4)
c
c ... ohcolm is the coulumbic term between O and H from different H2O in the unit of Hartree
      ohcolm=qo*qh*(1./rho1+1./rho2+1./rho3+1./rho4)
c ... hhcolm is ... between H and H ...
      hhcolm=qh*qh*(1./rhh1+1./rhh2+1./rhh3+1./rhh4)
c ... oocolm is ... between O and O ...
      oocolm=qo*qo*(1./rMM)
c
      E2H2O=o2lj+(ohcolm+oocolm+hhcolm)*hr2kcl*br2ang
c
      return
      end
