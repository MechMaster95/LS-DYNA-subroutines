      subroutine umat43v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsv,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,nnpcrv,cma,qmat,elsizv,idelev,
     . reject)
c******************************************************************
c|  NMCAD LAB MAT_002 (UMAT)                                      |
c|  ------------------------------------------------------------  |
c|  Copyright @ 2017 NMCAD LAB IISc                               |
c|  All Rights Reserved                                           |
c******************************************************************
      include 'nlqparm'
      include 'iounits.inc'
      DIMENSION d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      DIMENSION sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      DIMENSION cm(*),epsps(*),hsv(nlq,*),dt1siz(*)
      DIMENSION temps(*),crv(101,2,*),cma(*),qmat(nlq,3,3),elsizv(*)
      integer nnpcrv(*),i,j,m,n,k,node(8)
	  REAL*8 :: t,R,sigglb,sigloc,siz
	  DIMENSION :: t(6,6),R(6,6),sigloc(6),sigglb(6),siz(3)
	  REAL*8 E1,E2,E3,v12,v21,v23,v32,v31,v13,G12,G23,G31
	  REAL*8 l1(nlq),l2(nlq),l3(nlq),m1(nlq),m2(nlq),m3(nlq)
	  REAL*8 n1(nlq),n2(nlq),n3(nlq),xyz(3,8)
	  REAL*8 ::x17(nlq),x28(nlq),x35(nlq),x46(nlq),y17(nlq),y28(nlq),
     1 y35(nlq),y46(nlq),z17(nlq),z28(nlq),z35(nlq),z46(nlq)
	  REAL*8 :: S,C11,C12,C13,C23,C22,C33,C44,C55,C66
	  REAL*8 :: G1,G2,G3,G4,G5,G6,S1,S2,S3,S4,S5,S6,detF
	  REAL*8 ::FS11,FS12,FS13,FS21,FS22,FS23,FS31,FS32,FS33
      INTEGER idelev(*),ie
      LOGICAL failels(*),reject,ierr,res1,res2,res3,res4,res5
      CHARACTER*5 etype
	  REAL*8 ::
     & x1(nlq),x2(nlq),x3(nlq),x4(nlq),
     & x5(nlq),x6(nlq),x7(nlq),x8(nlq),
     & y1(nlq),y2(nlq),y3(nlq),y4(nlq),
     & y5(nlq),y6(nlq),y7(nlq),y8(nlq),
     & z1(nlq),z2(nlq),z3(nlq),z4(nlq),
     & z5(nlq),z6(nlq),z7(nlq),z8(nlq)
	  common/nani/no_hsvs
c create same common block in urmathn subroutine
c
c     eps(1)=local x  strain increment
c     eps(2)=local y  strain increment
c     eps(3)=local z  strain increment
c     eps(4)=local xy strain increment
c     eps(5)=local yz strain increment
c     eps(6)=local zx strain increment
c
c     sig(1)=local x  stress
c     sig(2)=local y  stress
c     sig(3)=local z  stress
c     sig(4)=local xy stress
c     sig(5)=local yz stress
c     sig(6)=local zx stress
c
c     material properties for orthotropic elastic material 
c     If transverse elastic properties are not provided, the transverse
c     directions take on the same properties as the longitudinal direction,
c     i.e. behaves as an isotropic material
c
      E1 = cm(1)
      E2 = cm(2)
      E3 = cm(3)
      v21 = cm(4)
      v31 = cm(5)
      v32 = cm(6)
      G12 = abs(cm(7))
      G23 = abs(cm(8))
      G31 = abs(cm(9))
	  m = no_hsvs
   
      v12 = v21*E1/E2
      v23 = v32*E2/E3
      v13 = v31*E1/E3
c
c     stiffness matrix 
c
      S = (1.-v12*v21-v23*v32-v31*v13-2*v21*v32*v13)/(E1*E2*E3)
c
      C11=(1. - v23*v32)/(E2*E3*S)
      C12=(v21 + v31*v23)/(E2*E3*S)
      C13=(v31 + v21*v32)/(E2*E3*S)
      C22=(1. - v31*v13)/(E1*E3*S)
      C23=(v32 + v31*v12)/(E1*E3*S)
      C33=(1. - v12*v21)/(E1*E2*S)
      C44=G12
      C55=G23
      C66=G31
c	  call mrchek(E1,E2,E3,v12,v21,v23,v32,v31,v13,G12,G23,G31,
c     1  res1,res2,res3,res4,res5)
c	 
c      if((res1.eqv..true.).and.(res2.eqv..true.).and.(res3.eqv..true.) 
c     1    .and.(res4.eqv..true.).and.(res5.eqv..true.)) then 
c	    goto 10 
c		print*,res1,res2,res3,res4,res5
c      else 
c		write(6,*) 'Check the material constants'
c		print*,res1,res2,res3,res4,res5
c      goto 11
c      end if
c10    continue
	  do 20 i=lft,llt
	  call lqfnode(i,2,ie,node,xyz,ierr)
		 x1(i) = xyz(1,1)
		 x2(i) = xyz(1,2)
		 x3(i) = xyz(1,3)
		 x4(i) = xyz(1,4)
		 x5(i) = xyz(1,5)
		 x6(i) = xyz(1,6)
		 x7(i) = xyz(1,7)
		 x8(i) = xyz(1,8)
		 
		 y1(i) = xyz(2,1)
		 y2(i) = xyz(2,2)
		 y3(i) = xyz(2,3)
		 y4(i) = xyz(2,4)
		 y5(i) = xyz(2,5)
		 y6(i) = xyz(2,6)
		 y7(i) = xyz(2,7)
		 y8(i) = xyz(2,8)
		 
		 z1(i) = xyz(3,1)
		 z2(i) = xyz(3,2)
		 z3(i) = xyz(3,3)
		 z4(i) = xyz(3,4)
		 z5(i) = xyz(3,5)
		 z6(i) = xyz(3,6)
		 z7(i) = xyz(3,7)
		 z8(i) = xyz(3,8)
		 
		 siz(1) = dsqrt(((x2(i)-x1(i))**2) +
     1                 ((y2(i)-y1(i))**2) +
     1                 ((z2(i)-z1(i))**2))
     	 
	     siz(2) = dsqrt(((x4(i)-x1(i))**2) +
     1                 ((y4(i)-y1(i))**2) +
     1                 ((z4(i)-z1(i))**2))
	  
	     siz(3) = dsqrt(((x5(i)-x1(i))**2) +
     1                 ((y5(i)-y1(i))**2) +
     1                 ((z5(i)-z1(i))**2))
        l1(i) = (x2(i)-x1(i))/siz(1)
		m1(i) = (y2(i)-y1(i))/siz(1)
		n1(i) = (z2(i)-z1(i))/siz(1)
	
		l2(i) = (x4(i)-x1(i))/siz(2)
		m2(i) = (y4(i)-y1(i))/siz(2)
		n2(i) = (z4(i)-z1(i))/siz(2)
		
c		l3(i) = (x5(i)-x1(i))/siz(3)
c		m3(i) = (y5(i)-y1(i))/siz(3)
c		n3(i) = (z5(i)-z1(i))/siz(3)

	    l3(i)=m1(i)*n2(i)-m2(i)*n1(i)
        m3(i)=n1(i)*l2(i)-n2(i)*l1(i)
        n3(i)=l1(i)*m2(i)-l2(i)*m1(i)
		
		l2(i)=m3(i)*n1(i)-m1(i)*n3(i)
        m2(i)=n3(i)*l1(i)-n1(i)*l3(i)
        n2(i)=l3(i)*m1(i)-l1(i)*m3(i)
		
c	  print*,'ELEMENT no',i
c	  print*, 'Coord X',x1(i),x2(i),x3(i),x4(i),x5(i),x6(i),x7(i),x8(i)
c	  print*, 'Coord Y',y1(i),y2(i),y3(i),y4(i),y5(i),y6(i),y7(i),y8(i)
c	  print*, 'Coord Z',z1(i),z2(i),z3(i),z4(i),z5(i),z6(i),z7(i),z8(i)

c      print*,'Direction Cosines_a',l1(i),m1(i),n1(i)
c      print*,'Direction Cosines_b',l2(i),m2(i),n2(i)
c      print*,'Direction Cosines_c',l3(i),m3(i),n3(i)
 
c     calculate transformation matrix

      t(1,1)=l1(i)*l1(i)
      t(1,2)=m1(i)*m1(i)
      t(1,3)=n1(i)*n1(i)
      t(1,4)=2.0*l1(i)*m1(i)
      t(1,5)=2.0*m1(i)*n1(i)
      t(1,6)=2.0*n1(i)*l1(i)
      t(2,1)=l2(i)*l2(i)
      t(2,2)=m2(i)*m2(i)
      t(2,3)=n2(i)*n2(i)
      t(2,4)=2.0*l2(i)*m2(i)
      t(2,5)=2.0*m2(i)*n2(i)
      t(2,6)=2.0*n2(i)*l2(i)
      t(3,1)=l3(i)*l3(i)
      t(3,2)=m3(i)*m3(i)
      t(3,3)=n3(i)*n3(i)
      t(3,4)=2.0*l3(i)*m3(i)
      t(3,5)=2.0*m3(i)*n3(i)
      t(3,6)=2.0*l3(i)*n3(i)
      t(4,1)=l1(i)*l2(i)
      t(4,2)=m1(i)*m2(i)
      t(4,3)=n1(i)*n2(i)
      t(4,4)=l2(i)*m1(i)+l1(i)*m2(i)
      t(4,5)=m2(i)*n1(i)+m1(i)*n2(i)
      t(4,6)=l2(i)*n1(i)+l1(i)*n2(i)
      t(5,1)=l2(i)*l3(i)
      t(5,2)=m2(i)*m3(i)
      t(5,3)=n2(i)*n3(i)
      t(5,4)=l2(i)*m3(i)+l3(i)*m2(i)
      t(5,5)=m2(i)*n3(i)+m3(i)*n2(i)
      t(5,6)=n2(i)*l3(i)+n3(i)*l2(i)
      t(6,1)=l3(i)*l1(i)
      t(6,2)=m3(i)*m1(i)
      t(6,3)=n3(i)*n1(i)
      t(6,4)=l3(i)*m1(i)+l1(i)*m3(i)
      t(6,5)=m3(i)*n1(i)+m1(i)*n3(i)
      t(6,6)=n3(i)*l1(i)+n1(i)*l3(i)
c      print*,'T1 Matrix---', t(1,1),t(1,2),t(1,3),t(1,4),t(1,5),t(1,6)
c      print*,'-------------------------------------------'
c      print*,'T2 Matrix---', t(2,1),t(2,2),t(2,3),t(2,4),t(2,5),t(2,6)
c	   print*,'-------------------------------------------'
c      print*,'T3 Matrix---', t(3,1),t(3,2),t(3,3),t(3,4),t(3,5),t(3,6)
c      print*,'-------------------------------------------'
c      print*,'T4 Matrix---', t(4,1),t(4,2),t(4,3),t(4,4),t(4,5),t(4,6)
c      print*,'-------------------------------------------'
c      print*,'T5 Matrix---', t(5,1),t(5,2),t(5,3),t(5,4),t(5,5),t(5,6)
c      print*,'-------------------------------------------'
c      print*,'T6 Matrix---', t(6,1),t(6,2),t(6,3),t(6,4),t(6,5),t(6,6)
c      print*,'-------------------------------------------'
	  
      R = TRANSPOSE(t)
   
      if (etype.eq.'solid'.or.etype.eq.'shl_t'.or.
     1     etype.eq.'sld2d'.or.etype.eq.'tshel'.or.
     2     etype.eq.'sph  '.or.etype.eq.'sldax') then
c       Right Green-St. Venant Strain Tensor
c       G1-G6 will be the local strains in Ea,Eb,Ec directions
        G1 = 0.5*(hsv(i,m+1)*hsv(i,m+1)+hsv(i,m+2)*hsv(i,m+2)+
     .       hsv(i,m+3)*hsv(i,m+3)-1.)
        G2 = 0.5*(hsv(i,m+4)*hsv(i,m+4)+hsv(i,m+5)*hsv(i,m+5)+
     .       hsv(i,m+6)*hsv(i,m+6)-1.)
        G3 = 0.5*(hsv(i,m+7)*hsv(i,m+7)+hsv(i,m+8)*hsv(i,m+8)+
     .       hsv(i,m+9)*hsv(i,m+9)-1.)
        G4 = 0.5*(hsv(i,m+1)*hsv(i,m+4)+hsv(i,m+2)*hsv(i,m+5)+
     .       hsv(i,m+3)*hsv(i,m+6))
        G5 = 0.5*(hsv(i,m+4)*hsv(i,m+7)+hsv(i,m+5)*hsv(i,m+8)+
     .       hsv(i,m+6)*hsv(i,m+9))
        G6 = 0.5*(hsv(i,m+1)*hsv(i,m+7)+hsv(i,m+2)*hsv(i,m+8)+
     .       hsv(i,m+3)*hsv(i,m+9))
c	 
c  Second Piola-Kirchhoff Stress Tensor
c  S1-S6 will be local stresses in Ea,Eb,Ec orientations
        S1 = C11*G1 + C12*G2 + C13*G3
        S2 = C12*G1 + C22*G2 + C23*G3
        S3 = C13*G1 + C23*G2 + C33*G3
        S4 = 2.*C44*G4
        S5 = 2.*C55*G5
        S6 = 2.*C66*G6
c First m slots for history variables(hsv)are unused
c Jacobian of the deformation gradient matrix
        detF=hsv(i,m+1)*(hsv(i,m+5)*hsv(i,m+9)-hsv(i,m+6)*hsv(i,m+8))-
     .       hsv(i,m+2)*(hsv(i,m+4)*hsv(i,m+9)-hsv(i,m+6)*hsv(i,m+7))+
     .       hsv(i,m+3)*(hsv(i,m+4)*hsv(i,m+8)-hsv(i,m+5)*hsv(i,m+7))
c
c       Cauchy Stresses (F*S*F^T)
        FS11 = hsv(i,m+1)*S1 + hsv(i,m+4)*S4 + hsv(i,m+7)*S6
        FS12 = hsv(i,m+1)*S4 + hsv(i,m+4)*S2 + hsv(i,m+7)*S5
        FS13 = hsv(i,m+1)*S6 + hsv(i,m+4)*S5 + hsv(i,m+7)*S3
        FS21 = hsv(i,m+2)*S1 + hsv(i,m+5)*S4 + hsv(i,m+8)*S6
        FS22 = hsv(i,m+2)*S4 + hsv(i,m+5)*S2 + hsv(i,m+8)*S5
        FS23 = hsv(i,m+2)*S6 + hsv(i,m+5)*S5 + hsv(i,m+8)*S3
        FS31 = hsv(i,m+3)*S1 + hsv(i,m+6)*S4 + hsv(i,m+9)*S6
        FS32 = hsv(i,m+3)*S4 + hsv(i,m+6)*S2 + hsv(i,m+9)*S5
        FS33 = hsv(i,m+3)*S6 + hsv(i,m+6)*S5 + hsv(i,m+9)*S3
c
c    sig1(i) - sig6(i) will be global stresses 
      sig1(i)=1./detF*(FS11*hsv(i,m+1)+FS12*hsv(i,m+4)+FS13*hsv(i,m+7))
      sig2(i)=1./detF*(FS21*hsv(i,m+2)+FS22*hsv(i,m+5)+FS23*hsv(i,m+8))
      sig3(i)=1./detF*(FS31*hsv(i,m+3)+FS32*hsv(i,m+6)+FS33*hsv(i,m+9))
      sig4(i)=1./detF*(FS11*hsv(i,m+2)+FS12*hsv(i,m+5)+FS13*hsv(i,m+8))
      sig5(i)=1./detF*(FS21*hsv(i,m+3)+FS22*hsv(i,m+6)+FS23*hsv(i,m+9))
      sig6(i)=1./detF*(FS11*hsv(i,m+3)+FS12*hsv(i,m+6)+FS13*hsv(i,m+9))
		
		hsv(i,1) = G1
		hsv(i,2) = G2
		endif
c
c  Calculating the local stresses 
      sigglb(1)=sig1(i)
      sigglb(2)=sig2(i)
      sigglb(3)=sig3(i)
      sigglb(4)=sig4(i)
      sigglb(5)=sig5(i)
      sigglb(6)=sig6(i)
      do 7 n=1,6
      sigloc(n)=0.0
      do 8 k=1,6
      sigloc(n)=sigloc(n)+t(n,k)*sigglb(k)
   8  continue
   7  continue
c******************************************************************
      hsv(i,3) = sigloc(1)
	  hsv(i,4) = sigloc(2)
	  hsv(i,5) = sigloc(3)
	  hsv(i,6) = sigloc(4)
	  hsv(i,7) = sigloc(5)
	  hsv(i,8) = sigloc(6)
c  his(6) should be compared with his(12)
c  his(7) should be compared with his(13)
c  his(8) should be compared with his(14)
	  hsv(i,9)  = sigglb(1)
	  hsv(i,10) = sigglb(2)
	  hsv(i,11) = sigglb(3)
	  hsv(i,12) = sigglb(4)
	  hsv(i,13) = sigglb(5)
	  hsv(i,14) = sigglb(6)
c******************************************************************
20    continue
c11    continue
c       
      return
      end
c
      subroutine mrchek(E1,E2,E3,n12,n21,n23,n32,n31,n13,G12,G23,G31,
     1  res1,res2,res3,res4,res5)
c This subroutine checks the restrictions on the material constants
c  of the Orthotropic material entered for the analysis
      real*8  E1,E2,E3,n12,n23,n31,n21,n32,n13,G12,G23,G31
      logical  res1,res2,res3,res4,res5
      real*8 s11,s22,s33,s12,s23,s13
      real*8 r1,r2,r3,r4,r5,r6,r7,r8,r9,rcn,rpr,delt
      real*8 rm1,rm2,rm3,a1,a2,a3,p1,p2,p3,p4,p5,p6
c Calculating the compliance components
       s11 = 1.0/E1
       s22 = 1.0/E2
       s33 = 1.0/E3
       s12 = -n21/E2
       s23 = -n32/E3
       s13 = -n31/E3
c   Restrictions on the poissons ratio's
        r1 = 1.0-n23*n32
        r2 = 1.0-n13*n31
        r3 = 1.0-n12*n21
c
       delt = 1.0-n12*n21-n13*n31-n23*n32-n12*n23*n31-n21*n32*n13
c
       r4 = sqrt(E2/E1)
       r5 = sqrt(E3/E2)
       r6 = sqrt(E1/E3)
       r7 = sqrt(E1/E2)
       r8 = sqrt(E2/E3)
       r9 = sqrt(E3/E1)

c  Restrictions on the material constants and poissons ratio combined
       rm1 = sqrt(s22*s33)
       rm2 = sqrt(s11*s33)
       rm3 = sqrt(s11*s22)
c
       a1 = abs(s23)
       a2 = abs(s13)
       a3 = abs(s12)
c
       p1 = abs(n21)
       p2 = abs(n32)
       p3 = abs(n13)
       p4 = abs(n12)
       p5 = abs(n23)
       p6 = abs(n31)
c
       rcn = (1.0-n21**2*(E1/E2)-n32**2*(E2/E3)-n13**2*(E3/E1))/2.0
       rpr = n21*n32*n13
c
      if((r1.gt.0).and.(r2.gt.0).and.(r3.gt.0).and.(delt.gt.0)) then
       res1 = .true.
      else
       res1 = .false.
       end if
      if((p1.lt.r4).and.(p2.lt.r5).and.(p3.lt.r6).and.(p4.lt.r7).and.
     1  (p5.lt.r8).and.(p6.lt.r9)) then
       res2 = .true.
      else 
       res2 = .false.
      end if
      if((a1.lt.rm1).and.(a2.lt.rm2).and.(a3.lt.rm3)) then
       res3 = .true.
        else 
         res3 = .false.
        end if
       if((rpr.lt.rcn).and.(rcn.lt.0.5)) then 
        res4 = .true.
       else 
        res4 = .false.
       end if
       if((G12.gt.0).and.(G23.gt.0).and.(G31.gt.0)) then
       res5 = .true.
       else
       res5 = .false.
       end if
       print*,res1,res2,res3,res4,res5
       return
       end
c
	  subroutine lqfnode(iext,ityp,ie,node,xyz,ierr)	
c
c     give element external id and
c     return internal element id, node internal id's and 
c     coordinates of element connectivity
c
c Input:
c     iext  I        external element id
c     ityp  I        type
c                     2- solid
c                     3- beam
c                     4- shell
c                     5- thick shell
c Output:
c     ie    I        internal id of the element
c     node  I(*)     internal node id's 
c                solid-8 nodes beam-2 nodes shell-4 nodes tshell-8 nodes
c     xyz   r*8(3,*) coordinates of nodes
c     ierr  I        error flag  0-ok 1-error
c
      include 'bk13.inc'
      include 'memaia.inc'
c
c     dynamic memory allocation stuff.  This has to be a C style 
c     include because of the use of integer*8
c
      integer i_mem(0:1)
      real r_mem(0:1)
      real*8 r8_mem(0:1)
      real*8 real8_mem
      common/dynmem/real8_mem(0:1)
      equivalence (i_mem,r_mem,r8_mem,real8_mem)
      integer*8 mem_alloc, mems_alloc,mem_realloc, memh_ptr
      integer mem_length,memh_alloc,memh_realloc,memh_length
      integer memh_enum,memh_attach
      integer mem_newmark,mem_getmark,memh_getmark
      integer*8 memh_detach,memh_ralloc
      integer*8 m_to_m8
      external mem_alloc, mems_alloc, mem_realloc, mem_length
      external memh_alloc, memh_realloc, memh_length, memh_ptr
      external memh_enum,memh_attach,memh_detach,memh_ralloc
      external m_to_m8
      integer*8 dm_x, dm_v, dm_xms, dm_me1, dm_rots, dm_a
      integer*8 dm_x0,dm_v0,dm_xms0,dm_me10,         dm_a0
      integer*8 dm_disp,dm_xtpz
      common /dynmem1/ dm_x, dm_v, dm_xms, dm_me1, dm_rots, dm_a,
     &                 dm_x0,dm_v0,dm_xms0,dm_me10,         dm_a0,
     &                 dm_disp,dm_xtpz
c
      integer*8 mem_allocm,mems_allocm,mem_reallocm
      integer memh_allocm,memh_reallocm
      external mem_allocm,mems_allocm,mem_reallocm
      external memh_allocm,memh_reallocm
      integer*8 mem_allocmchk,mem_allocchk,mems_allocmchk,mems_allocchk
      integer*8 mem_reallocmchk,mem_reallocchk
      integer memh_allocmchk,memh_allocchk
      external mem_allocmchk,mem_allocchk,mems_allocmchk,mems_allocchk
      external mem_reallocmchk,mem_reallocchk
      external memh_allocmchk,memh_allocchk
c
      dimension node(*)
      real*8 xyz(3,*)
c
c     convert user ID to internal number
      ie=lqfint(iext,ityp,ierr)
      if (ierr.ne.0) go to 999
c
c     solid
      if (ityp.eq.2) then
        l=lc1h+(ie-1)*9
        do i=1,8
          node(i)=ia(l+i)
          call getdxyz(r_mem(dm_x),node(i),xyz(1,i))
        enddo
c
c     beam
      elseif (ityp.eq.3) then
        l=lc1b+(ie-1)*4
        do i=1,2
          node(i)=ia(l+i)
          call getdxyz(r_mem(dm_x),node(i),xyz(1,i))
        enddo
c
c     shell
      elseif (ityp.eq.4) then
        l=lc1s+(ie-1)*5
        do i=1,4
          node(i)=ia(l+i)
          call getdxyz(r_mem(dm_x),node(i),xyz(1,i))
        enddo
c
c     thick shell
      elseif (ityp.eq.5) then
        l=lc1t+(ie-1)*9
        do i=1,8
          node(i)=ia(l+i)
          call getdxyz(r_mem(dm_x),node(i),xyz(1,i))
        enddo
      endif
c
  999 continue
      return
      end
c
      subroutine getdxyz(x,n,xx)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     get nodal coordinate x
c
      real*8 x(3,*)
      real*8 xx(3)
	  real*8 vx(3)
	  real*8 ax(3)
c
      xx(1)=x(1,n)
      xx(2)=x(2,n)
      xx(3)=x(3,n) 
c  vx(1) = velocity in x
c  vx(2) = velocity in y
c  vx(3) = velocity in z
c  ax(1) = acceleration in x
c  ax(2) = acceleration in y
c  ax(3) = acceleration in z
      vx(1) = x(4,n)
      vx(2) = x(5,n)
      vx(3) = x(6,n)
      ax(1) = x(7,n)
      ax(2) = x(8,n)
      ax(3) = x(9,n)
c	write(6,*) 'vx=',vx(1)
c	write(6,*) 'vy=',vx(2)
c	write(6,*) 'vz=',vx(3)
c
      return
      end