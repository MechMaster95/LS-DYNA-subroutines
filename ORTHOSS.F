      subroutine umat43(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt)

c Orthotropic elastic material User subroutine
c Variables
c
c This subroutine is written by N.L.Vishnuvardhan Raju
c
c cm(1)=first material constant, here young's modulus in a direction
c cm(2)=second material constant, youngs modulus in b direction
c 
c 
c 
c cm(n)=nth material constant
c
c eps(1)=local x strain increment
c eps(2)=local y strain increment
c eps(3)=local z strain increment
c eps(4)=local xy strain increment
c eps(5)=local yz strain increment
c eps(6)=local zx strain increment
c
c sig(1)=local x stress
c sig(2)=local y stress
c sig(3)=local z stress
c sig(4)=local xy stress
c sig(5)=local yz stress
c sig(6)=local zx stress
c
c hsv(1)=1st history variable
c hsv(2)=2nd history variable
c .
c .
c .
c .
c hsv(n)=nth history variable
c
c dt1=current time step size
c capa=reduction factor for transverse shear
c etype:
c eq."solid" for solid elements
c eq."shell" for all other shell elements plus thick shell forms 1
c and 2
c
c tt=current problem time.
c
c cma=additional memory for material data defined by LMCA at
c 6th field of 2nd crad of *DATA_USER_DEFINED
c
c elsiz=characteristic element size
c
c idele=element id
c
c All transformations into the element local system are
c performed prior to entering this subroutine. Transformations
c back to the global system are performed after exiting this
c routine.
c
c All history variables are initialized to zero in the input
c phase. Initialization of history variables to nonzero values
c may be done during the first call to this subroutine for each
c element.
c
c
      DIMENSION :: cm(*),eps(*),sig(*),hsv(*)
      REAL      :: Ea,Eb,Ec,v12,v23,v31,v21,v13,v32,Gab,Gbc,Gca,tt
      REAL      :: prdr,prdr1,prdr2,prdr3,s1nu1,s1nu2,s1nu3,efs
      REAL      :: s2nu1,s2nu2,s2nu3,s3nu1,s3nu2,s3nu3,Gyz,Gzx
      REAL      :: Xc,Xt,Yt,Yc,Zc,Zt,Qtc,Rtc,Stc
      logical   :: rcn1,rcn2,rcn3,rcn4,rcn5
      CHARACTER*5  etype
c
c    Assigning the values
c

       Ea = cm(1)
       Eb = cm(2)
       Ec = cm(3)
       
       v12 = cm(4)
       v23 = cm(5)
       v31 = cm(6) 

       Gab = cm(7)
       Gbc = cm(8)
       Gca = cm(9)

c   Assigning the strength variables
	  
	  
      Xc  = cm(12)
	  Xt  = cm(13)
	  Yc  = cm(14)
	  Yt  = cm(15)
	  Zt  = cm(16)
	  Zc  = cm(17)
	  Qtc = cm(18)
	  Rtc = cm(19)
	  Stc = cm(20)
	  
	  

c     Computing the Poissons ratios and multiple's

       v21 =(v12*Eb)/Ea
       v32 =(v23*Ec)/Eb
       v13 =(v31*Ea)/Ec

       
        prdr  = 1.0-v12*v21-v13*v31-v23*v32-v12*v23*v31-v13*v21*v32
        prdr1 = Ea/prdr
        prdr2 = Eb/prdr
        prdr3 = Ec/prdr

        s1nu1 = 1.0-v23*v32
        s1nu2 = v21+v23*v32
        s1nu3 = v31+v21*v32

        s2nu1 = v12+v13*v32
        s2nu2 = 1.0-v13*v31
        s2nu3 = v32+v12*v31

        s3nu1 = v13+v12*v23
        s3nu2 = v23+v13*v21
        s3nu3 = 1.0-v12*v21
        
c     checking restrictions on the Material constants
      call mrchek(Ea,Eb,Ec,v12,v23,v31,v21,v32,v13,Gab,Gbc,
     1   Gca,rcn1,rcn2,rcn3,rcn4,rcn5)
     
        IF(rcn1.eqv..true..and.rcn2.eqv..true..and.rcn3.eqv..true..and.
     1   rcn4.eqv..true..and.rcn5.eqv..true.) THEN

        IF (etype.eq.'solid') THEN

        sig(1) = sig(1)+prdr1*(eps(1)*s1nu1+eps(2)*s1nu2+eps(3)*s1nu3)
        sig(2) = sig(2)+prdr2*(eps(1)*s2nu1+eps(2)*s2nu2+eps(3)*s2nu3)
        sig(3) = sig(3)+prdr3*(eps(1)*s3nu1+eps(2)*s3nu2+eps(3)*s3nu3)
        sig(4) = sig(4)+Gab*eps(4)
        sig(5) = sig(5)+Gbc*eps(5)
        sig(6) = sig(6)+Gca*eps(6)

        ELSE IF (etype.eq.'shell') THEN

        Gyz = capa*Gbc
        Gzx = capa*Gca

        eps(3) = (eps(2)*s3nu2+eps(1)*s3nu1)/(-s3nu3)        
        sig(1) = sig(1)+prdr1*(eps(1)*s1nu1+eps(2)*s1nu2+eps(3)*s1nu3)
        sig(2) = sig(2)+prdr2*(eps(1)*s2nu1+eps(2)*s2nu2+eps(3)*s2nu3)
        sig(3) = 0.0
        sig(4) = sig(4)+Gab*eps(4)
        sig(5) = sig(5)+Gyz*eps(5)
        sig(6) = sig(6)+Gzx*eps(6)

        END IF
        ELSE 
        WRITE(*,*) 'CHECK THE MATERIAL CONSTANTS'
        END IF
        efs = -(eps(1)+eps(2)+eps(3))/3.0
        hsv(1) = hsv(1)+eps
        hsv(2) = - (1.0/3.0) * (sig(1)+sig(2)+sig(3))
        RETURN
        END
        
        subroutine mrchek(E1,E2,E3,n12,n23,n31,n21,n32,n13,G12,G23,G31,
     1  res1,res2,res3,res4,res5)
	 
c      implicit none
c This subroutine checks the restrictions on the material constants
c  of the Orthotropic material entered for the analysis
      real  E1,E2,E3,n12,n23,n31,n21,n32,n13,G12,G23,G31
      logical  res1,res2,res3,res4,res5
      real s11,s22,s33,s12,s23,s31,s13
      real r1,r2,r3,r4,r5,r6,r7,r8,r9,rcn,rpr,delt
      real rm1,rm2,rm3,a1,a2,a3,p1,p2,p3,p4,p5,p6
c
c
c Calculating the compliance components
       s11 = 1.0/E1
       s22 = 1.0/E2
       s33 = 1.0/E3
       s12 = -n21/E2
       s23 = -n32/E3
       s31 = -n13/E1
c      s21 = -n12/E1
c      s32 = -n23/E3
       s13 = -n31/E3
c
c
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
       rcn = (1.0-n21**2*(E2/E1)-n32**2*(E2/E3)-n13**2*(E3/E1))/2.0
       rpr = n21*n32*n13
c
      if(r1.gt.0 .and.r2.gt.0 .and.r3.gt.0 .and. delt.gt.0) then
       res1 = .true.
      else
       res1 = .false.
       end if
      if(p1.lt.r4 .and. p2.lt.r5 .and. p3.lt.r6 .and. p4.lt.r7 .and.
     1  p5.lt.r8 .and. p6.lt.r9) then
	      res2 = .true.
      else 
       res2 = .false.
      end if
      if(a1.lt.rm1 .and. a2.lt.rm2 .and. a3.lt.rm3) then
       res3 = .true.
        else 
         res3 = .false.
        end if
       if(rpr.lt.rcn .and. rcn.lt. 0.5) then 
        res4 = .true.
       else 
        res4 = .false.
       end if
       if(G12.gt.0 .and. G23.gt.0 .and. G31.gt.0 ) then
       res5 = .true.
       else
       res5 = .false.
       end if
c
       end

