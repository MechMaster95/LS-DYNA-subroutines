       subroutine umat43 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt, 
     1 temper,failel,crv,nnpcrv,cma)
c
c******************************************************************
c|N.L.Vishnuvardhan Raju                                            |
c| ------------------------------------------------------------ |
c| Copyright @ NMCAD Lab 2017                                   |
c| All rights reserved                                          |
c******************************************************************
c
c Orthotropic Elasto-Plastic material Model using Non-Associative flow
c rule
c
c cm(1)= E11 Young's modulus in 1-direction
c cm(2)= E22 Young's modulus in 2-direction
c cm(3)= E33 Young's modulus in 3-direction
c cm(4)= Poisson's ratio v21
c cm(5)= Poisson's ratio v31
c cm(6)= Poisson's ratio v32
c cm(7)= Shear Modulus G12
c cm(8)= Shear Modulus G23
c cm(9)= Shear Modulus G31
c
c cm(10)= Bulk modulus
c cm(11)= Shear modulus
c
c cm(12)-cm(23) = 12 stress vs strain curves input id's.
c cm(24)-cm(32) = H11,H22,H33,H12,H23,H31,H44,H55,H66  (PLASTIC
c                                            POTENTIAL COEFFICIENTS)
c
c cm(33)-cm(40) = sig11t,sig22t,sig33t,sig11c,sig22c,sig33c,sig12,sig23
c cma(1)-cma(4) = sig31,sig4512t,sig4523t,sig4531t

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
c hsv(1)= effective plastic strain.
c hsv(2)= Total hydrostatic pressure.
c
c hsv(3)= stress state in 1-dir tensile load corresponding to epsp
c hsv(4)= stress state in 2-dir tensile load corresponding to epsp
c hsv(5)= stress state in 3-dir tensile load corresponding to epsp
c hsv(6)= stress state in 1-dir compressive load corresponding to epsp
c hsv(7)= stress state in 2-dir compressive load corresponding to epsp
c hsv(8)= stress state in 3-dir compressive load corresponding to epsp
c hsv(9)= stress state in 1-2 plane shear load corresponding to epsp
c hsv(10)= stress state in 2-3 plane shear load corresponding to epsp
c hsv(11)= stress state in 3-1 plane shear load corresponding to epsp
c hsv(12)= Sress state corresponding to 45deg specimen tensile load
c hsv(13)= Sress state corresponding to 45deg specimen tensile load
c hsv(14)= Sress state corresponding to 45deg specimen tensile load
c hsv(15)= Scalar plastic multiplier, (lda)
c 
c
c
c dt1  = current time step size
c capa = reduction factor for transverse shear
c etype .eq. "solid"  for solid elements
c 
c
c tt     = current problem time.
c
c temper = current temperature
c
c failel = flag for failure, set to .true. to fail an integration point,
c if .true. on input the corresponding integration point will fail.
c
c crv(lq1,2,*) = array representation of curves in keyword deck.
c
c nnpcrv(*)    = no. of discretization points per crv()
c
c cma  = additional memory for material data defined by LMCA in .k file
c
c elsiz=characteristic element size
c
c idele=element id
c
c 
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
c Energy calculations for the dyna3d energy balance are done
c outside this subroutine.
c
      DIMENSION :: cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      INTEGER   :: nnpcrv(*)
      REAL      :: capa,epsp,dt1,tt,temper
      
      REAL      :: Ea,Eb,Ec,v12,v23,v31,v21,v32,v13,Gab,Gbc,Gca
      REAL      :: prdr,s1nu1,s1nu2,s1nu3
      REAL      :: s2nu1,s2nu2,s2nu3,s3nu1,s3nu2,s3nu3
C      REAL      :: Xt,Xc,Yt,Yc,Zt,Zc,P12,Q23,R31,ST12,ST23,ST31
      REAL      :: cc,eid,deps,sigold,qst,lda,strss
      DIMENSION :: cc(6,6),eid(12),deps(6),sigold(6),qst(12),strss(6)
      REAL      :: H,QT,neps
      DIMENSION :: H(9),QT(12),neps(6)
      INTEGER   :: i,j
      logical   :: 
      logical   :: failel 
      CHARACTER*5  etype
      common /coeff/ cc,eid,QT,H,lq1

c
c    Assigning the values
c

       Ea = cm(1)
       Eb = cm(2)
       Ec = cm(3)
       
       v21 = cm(4)
       v31 = cm(5)
       v32 = cm(6) 

       Gab = cm(7)
       Gbc = cm(8)
       Gca = cm(9)
	  
c  Assigning the curve id's
       eid(1) = cm(12)
       eid(2) = cm(13)
       eid(3) = cm(14)
       eid(4) = cm(15)
       eid(5) = cm(16)
       eid(6) = cm(17)
       eid(7) = cm(18)
       eid(8) = cm(19)
       eid(9) = cm(20)
       eid(10) = cm(21)
       eid(11) = cm(22)
       eid(12) = cm(23)


c  Assigning the Hij coefficients
        H(1) = cm(24)  !H11
        H(2) = cm(25)  !H22
        H(3) = cm(26)  !H33
        H(4) = cm(27)  !H12
        H(5) = cm(28)  !H23
        H(6) = cm(29)  !H13
        H(7) = cm(30)  !H44
        H(8) = cm(31)  !H55
        H(9) = cm(32)  !H66
 
c  Assigning the strength variables
        QT(1)  = cm(33)  !Xt
        QT(2)  = cm(34)  !Yt
        QT(3)  = cm(35)  !Zt 
        QT(4)  = cm(36)  !Xc
        QT(5)  = cm(37)  !Yc
        QT(6)  = cm(38)  !Zc
        QT(7)  = cm(39)  !P
        QT(8)  = cm(40)  !Q
        QT(9)  = cm(41)  !R
        QT(10) = cm(42)  !S4512T
        QT(11) = cm(43)  !S4523T
        QT(12) = cm(44)  !S4531T

c     Computing the Poissons ratios and Formulating Cij Matrix.                              

       v12 =(v21*Ea)/Eb
       v13 =(v31*Ea)/Ec
       v23 =(v32*Eb)/Ec

       
        prdr  = 1.0-v12*v21-v13*v31-v12*v23*v31-v13*v21*v32-v23*v32

        s1nu1 = Ea*(1.0 - v23*v32)
        s1nu2 = Ea*(v21 + v23*v31)
        s1nu3 = Ea*(v31 + v21*v32)

        s2nu1 = Eb*(v12 + v13*v32)
        s2nu2 = Eb*(1.0 - v13*v31)
        s2nu3 = Eb*(v32 + v12*v31)

        s3nu1 = Ec*(v13 + v12*v23)
        s3nu2 = Ec*(v23 + v13*v21)
        s3nu3 = Ec*(1.0 - v12*v21)

        do 10, i=1,6
         do 20, j = 1,6
           cc(i,j) = 0.0
20      continue
10      continue

       cc(1,1) = s1nu1/prdr
       cc(1,2) = s1nu2/prdr
       cc(1,3) = s1nu3/prdr
       cc(2,1) = s2nu1/prdr
       cc(2,2) = s2nu2/prdr
       cc(2,3) = s2nu3/prdr
       cc(3,1) = s3nu1/prdr
       cc(3,2) = s3nu2/prdr
       cc(3,3) = s3nu3/prdr
       cc(4,4) = Gab
       cc(5,5) = Gbc
       cc(6,6) = Gca

c    Assign inputs for main code

      deps(1)   = eps(1)
      deps(2)   = eps(2)
      deps(3)   = eps(3)
      deps(4)   = eps(4)
      deps(5)   = eps(5)
      deps(6)   = eps(6)
c
      sigold(1) = sig(1)
      sigold(2) = sig(2)
      sigold(3) = sig(3)
      sigold(4) = sig(4)
      sigold(5) = sig(5)
      sigold(6) = sig(6)

       epn  = hsv(1)
       
       DO 22, i = 1,12
         qst(i) = hsv(i+2)
22     CONTINUE

       lda = hsv(15)
       
c *********************************************************************
c **     Compute Stress based on 3D Elasto plastic model            **
c *********************************************************************
c
      call elpl3d(sigold,deps,epn,crv,nnpcrv,strss,neps,lda,qst)
c
c *********************************************************************
c
c    Update Stresses back from the subroutine
c
      sig(1)   = strss(1)
      sig(2)   = strss(2)
      sig(3)   = strss(3)
      sig(4)   = strss(4)
      sig(5)   = strss(5)
      sig(6)   = strss(6)
c
c    Calculate Total Hydrostatic Pressure
c
      davg     = -(1.0/3.0) * (sig(1)+sig(2)+sig(3))
      epsp     = epn

c
c    Update History Variables
c
      hsv(1)  = epn
      hsv(2)  = davg
      do 21,i= 1,12
       hsv(i+2) = qst(i)
21    continue
c
      hsv(15) = lda

      end 

c******************************************************************

      subroutine elpl3d(sigold,deps,epn,crv,nnpcrv,strss,neps,lda,qst)

c*****************************************************************
c  Here we calculate the trail stress and the stress are returned to
c  the plastic state if crossed the yield value, using the radial
c  return technique for scaling and secant iteration technique for
c  calculating the plastic multiplier.

       REAL :: sigold, deps,strss,neps,qst
       dimension ::sigold(6),deps(6),strss(6),neps(6),qst(12)
       real :: dlambda,lda,crv,epsp,tw1,tw2
       dimension :: crv(lq1,2,*),epsp(6)
       integer :: nnpcrv(*),i,j,ii,jj
       real :: cc,eid,QT,H,dqst,hdr
       dimension :: cc(6,6),eid(12),QT(12),H(9),dqst(12),hdr(6)
       logical tfail,stfail

       common /coeff/ cc,eid,QT,H,lq1
c
       ep  = epn
       
       do i=1,6
          epsp(i)    = 0.0
          sigt(i)     = sigold(i)
       end do 

       do 23, i = 1,6
        do 24, j = 1,6
         sigt(i) = sigt(i) + cc(i,j)*deps(j)
24      continue
23      continue

       call  tsaiwu(sigt,QT,tfail,tw1)
       call  tsaiwu(sigt,qst,stfail,tw2)
       if (tw1 .gt. tw2) then
          if (tfail .eqv. .true.) then
             call hardfn(sigt,deps,crv,nnpcrv,QT,lda,hdr)
          end if
          do 50,i = 1,6
            do 51, j = 1,6
            sigt(i) = sigt(i)-lda*cc(i,j)*deps(j)
51          continue
50        continue

          call crvval(crv,nnpcrv,eid(1),lda,dqst(1),csp(1))
          call crvval(crv,nnpcrv,eid(2),lda,dqst(2),csp(2))
          call crvval(crv,nnpcrv,eid(3),lda,dqst(3),csp(3))
          call crvval(crv,nnpcrv,eid(4),lda,dqst(4),csp(4))
          call crvval(crv,nnpcrv,eid(5),lda,dqst(5),csp(5))
          call crvval(crv,nnpcrv,eid(6),lda,dqst(6),csp(6))
          call crvval(crv,nnpcrv,eid(7),lda,dqst(7),csp(7))
          call crvval(crv,nnpcrv,eid(8),lda,dqst(8),csp(8))
          call crvval(crv,nnpcrv,eid(9),lda,dqst(9),csp(9))
          call crvval(crv,nnpcrv,eid(10),lda,dqst(10),csp(10))
          call crvval(crv,nnpcrv,eid(11),lda,dqst(11),csp(11))
          call crvval(crv,nnpcrv,eid(12),lda,dqst(12),csp(12))

       else
c          if (tw2 .gt. tw1) then
            if (stfail .eqv. .true.) then
             call hardfn(sigt,deps,crv,nnpcrv,qst,lda,hdr)
            end if
             do 52,i = 1,6
              do 53, j = 1,6
                 sigt(i) = sigt(i)-lda*cc(i,j)*deps(j)
53            continue
52           continue

          call crvval(crv,nnpcrv,eid(1),lda,dqst(1),csp(1))
          call crvval(crv,nnpcrv,eid(2),lda,dqst(2),csp(2))
          call crvval(crv,nnpcrv,eid(3),lda,dqst(3),csp(3))
          call crvval(crv,nnpcrv,eid(4),lda,dqst(4),csp(4))
          call crvval(crv,nnpcrv,eid(5),lda,dqst(5),csp(5))
          call crvval(crv,nnpcrv,eid(6),lda,dqst(6),csp(6))
          call crvval(crv,nnpcrv,eid(7),lda,dqst(7),csp(7))
          call crvval(crv,nnpcrv,eid(8),lda,dqst(8),csp(8))
          call crvval(crv,nnpcrv,eid(9),lda,dqst(9),csp(9))
          call crvval(crv,nnpcrv,eid(10),lda,dqst(10),csp(10))
          call crvval(crv,nnpcrv,eid(11),lda,dqst(11),csp(11))
          call crvval(crv,nnpcrv,eid(12),lda,dqst(12),csp(12))
          
       end if

       do 60,ii = 1,6
         strss(ii) = sigt(ii)
60     continue
       do 61, jj =1,12
       qst(jj) = dqst(jj)
61     continue

       do 63,i = 1,6
         do 64,j = 1,6
          ep = lda*cc(i,j)*hdr(j)
64      continue
63      continue
        epn = epn + ep

        return
        end
c*********************************************************************
      subroutine tsaiwu(sig,qmat,tfail,tw)
      
c  subroutine to calculate initial yield using the strengths values
c    provided

      real :: sig,qmat
      dimension :: sig(6),qmat(12)
      real :: F1,F2,F3,F11,F22,F33,F12,F23,F13,F44,F55,F66
      REAL :: N1,N2,N3,N4,N5,N6,tw
c      real  ::tol1,tol2
      logical tfail

c     common /coeff/ cc,eid,qmat,H,lq1

c      tol1 =  2.0e-3
c      tol2 = -2.0e-3


      F1 = (1.0/qmat(1))-(1.0/qmat(4))
      F2 = (1.0/qmat(2))-(1.0/qmat(5))
      F3 = (1.0/qmat(3))-(1.0/qmat(6))
      F11 = 1.0/(qmat(1)*qmat(4))
      F22 = 1.0/(qmat(2)*qmat(5))
      F33 = 1.0/(qmat(3)*qmat(6))
      F44 = 1.0/qmat(7)**2
      F55 = 1.0/qmat(8)**2
      F66 = 1.0/qmat(9)**2
      F12 = (2.0/qmat(10)**2)-((F1+F2)/qmat(10))-0.5*(F11+F22+F44)
      F23 = (2.0/qmat(11)**2)-((F3+F2)/qmat(11))-0.5*(F33+F22+F55)
      F13 = (2.0/qmat(12)**2)-((F1+F3)/qmat(12))-0.5*(F11+F33+F66)

      N1 = F1*sig(1)+F2*sig(2)+F3*sig(3)-1.0
      N2 = F11*sig(1)**2+F22*sig(2)**2+F33*sig(3)**2
      N3 = 2.0*F12*sig(1)*sig(2)
      N4 = 2.0*F23*sig(2)*sig(3)
      N5 = 2.0*F13*sig(1)*sig(3)
      N6 = F44*sig(4)**2+F55*sig(5)**2+F66*sig(6)**2

      tw = N1+N2+N3+N4+N5+N6

      if(tw .gt. 0.0 ) then
        tfail = .true.
      else
         tfail = .false.
      end if

      end
c**********************************************************************
      subroutine hardfn(sig,eps,crv,nnpcrv,qmat,lda,dph)
c****************************************************************
c  This is the subroutine to calculate incremental plastic strain.

      real :: sig,crv,eid,H,qmat,lda,dlambda,dlambda1,df
      dimension sig(6),crv(lq1,2,*)eid(12),H(9),qmat(12),df(6)
      real :: eps,cc,eid,QT,csp1,csp2,dlda1,dlda2,dlda3
      dimension :: eps(6),cc(6,6),eid(12),QT(12),csp1(12),csp2(12)
      integer :: nnpcrv(*),i,j,ii,jj
      real :: F1,F2,F3,F11,F22,F33,F44,F55,F66,F12,F23,F13
      real :: tf1,tf2,tf3,tol1,tol2,dqmat,dlambda2,dlambda3,csp3
      dimension ::csp3(12),dqmat(12)
      real  :: dsts,nr1,dr1,dph
      dimension :: dsts(6),dph(6)
      real  ::m1,m2,m3,m4,m5,m6,m7,pph
      logical :: tfail

      common /coeff/ cc,eid,QT,H,lq1

      tol1 = 2.0E-03
      tol2 = -2.0E-03
      
c      do 25,i = 1,12
c       dqmat(i) = qmat(i)
c25    continue

c      do 35, i = 1,6
c      dsts(i) = sig(i)
c35    continue

      F1 = (1.0/qmat(1))-(1.0/qmat(4))
      F2 = (1.0/qmat(2))-(1.0/qmat(5))
      F3 = (1.0/qmat(3))-(1.0/qmat(6))
      F11 = 1.0/(qmat(1)*qmat(4))
      F22 = 1.0/(qmat(2)*qmat(5))
      F33 = 1.0/(qmat(3)*qmat(6))
      F44 = 1.0/qmat(7)**2
      F55 = 1.0/qmat(8)**2
      F66 = 1.0/qmat(9)**2
      F12 = (2.0/qmat(10)**2)-((F1+F2)/qmat(10))-0.5*(F11+F22+F44)
      F23 = (2.0/qmat(11)**2)-((F3+F2)/qmat(11))-0.5*(F33+F22+F55)
      F13 = (2.0/qmat(12)**2)-((F1+F3)/qmat(12))-0.5*(F11+F33+F66)

      df(1) = F1+2.0*F11*sig(1)+2.0*F12*sig(2)+2.0*F13*sig(3)
      df(2) = F2+2.0*F12*sig(1)+2.0*F22*sig(2)+2.0*F23*sig(3)
      df(3) = F3+2.0*F13*sig(1)+2.0*F23*sig(2)+2.0*F33*sig(3)
      df(4) = 2.0*F44*sig(4)
      df(5) = 2.0*F55*sig(5)
      df(6) = 2.0*F66*sig(6)

      m1 = H(1)*sig(1)**2+H(2)*sig(2)**2+H(3)*sig(3)**2
      m2 = 2.0*H(4)*sig(1)*sig(2)
      m3 = 2.0*H(5)*sig(2)*sig(3)
      m4 = 2.0*H(6)*sig(1)*sig(3)
      m5 = H(7)*sig(4)**2
      m6 = H(8)*sig(5)**2
      m7 = H(9)*sig(6)**2
      pph = sqrt(m1+m2+m3+m4+m5+m6+m7)
      
      dph(1)  = (H(1)*sig(1)+H(4)*sig(2)*H(6)*sig(3))/pph
      dph(2)  = (H(4)*sig(1)+H(2)*sig(2)*H(5)*sig(3))/pph
      dph(3)  = (H(6)*sig(1)+H(5)*sig(2)*H(3)*sig(3))/pph
      dph(4)  = (H(7)*sig(4))/pph
      dph(5)  = (H(8)*sig(5))/pph
      dph(6)  = (H(9)*sig(6))/pph

      nr1 = 0.0
      do 20, i = 1,6
       do 22, j = 1,6
       nr1 = nr1+df(i)*cc(i,j)*eps(j)
22     continue
20    continue

      dr1 = 0.0
      do 23,i=1,6
       do 24, j = 1,6
        dr1 = dr1+df(i)*cc(i,j)*dph(j)
24     continue
23     continue

      dlambda1 = 0.0
      dlambda2 = nr1/dr1

100   continue
      dlda1 = lda + dlambda1

      do 40,ii = 1,6
       do 41,jj = 1,6
      dsts(ii) = sig(ii)-dlda1*cc(ii,jj)*dph(jj)
41    continue
40    continue
      
      call crvval(crv,nnpcrv,eid(1),dlda1,dqmat(1),csp1(1))
      call crvval(crv,nnpcrv,eid(2),dlda1,dqmat(2),csp1(2))
      call crvval(crv,nnpcrv,eid(3),dlda1,dqmat(3),csp1(3))
      call crvval(crv,nnpcrv,eid(4),dlda1,dqmat(4),csp1(4))
      call crvval(crv,nnpcrv,eid(5),dlda1,dqmat(5),csp1(5))
      call crvval(crv,nnpcrv,eid(6),dlda1,dqmat(6),csp1(6))
      call crvval(crv,nnpcrv,eid(7),dlda1,dqmat(7),csp1(7))
      call crvval(crv,nnpcrv,eid(8),dlda1,dqmat(8),csp1(8))
      call crvval(crv,nnpcrv,eid(9),dlda1,dqmat(9),csp1(9))
      call crvval(crv,nnpcrv,eid(10),dlda1,dqmat(10),csp1(10))
      call crvval(crv,nnpcrv,eid(11),dlda1,dqmat(11),csp1(11))
      call crvval(crv,nnpcrv,eid(12),dlda1,dqmat(12),csp1(12))

      call tsaiwu(dsts,dqmat,tfail,tf1)

      dlda2 = lda + dlambda2

      do 42,ii = 1,6
       do 43,jj = 1,6
      dsts(ii) = sig(ii)-dlda2*cc(ii,jj)*dph(jj)
43    continue
42    continue
      call crvval(crv,nnpcrv,eid(1),dlda2,dqmat(1),csp2(1))
      call crvval(crv,nnpcrv,eid(2),dlda2,dqmat(2),csp2(2))
      call crvval(crv,nnpcrv,eid(3),dlda2,dqmat(3),csp2(3))
      call crvval(crv,nnpcrv,eid(4),dlda2,dqmat(4),csp2(4))
      call crvval(crv,nnpcrv,eid(5),dlda2,dqmat(5),csp2(5))
      call crvval(crv,nnpcrv,eid(6),dlda2,dqmat(6),csp2(6))
      call crvval(crv,nnpcrv,eid(7),dlda2,dqmat(7),csp2(7))
      call crvval(crv,nnpcrv,eid(8),dlda2,dqmat(8),csp2(8))
      call crvval(crv,nnpcrv,eid(9),dlda2,dqmat(9),csp2(9))
      call crvval(crv,nnpcrv,eid(10),dlda2,dqmat(10),csp2(10))
      call crvval(crv,nnpcrv,eid(11),dlda2,dqmat(11),csp2(11))
      call crvval(crv,nnpcrv,eid(12),dlda2,dqmat(12),csp2(12))

      call tsaiwu(dsts,dqmat,tfail,tf2)

      dlambda3 = dlambda2 - tf2*((dlambda2-dlambda1)/(tf2-tf1))
      dlda3 = lda + dlambda3
      do 44,ii = 1,6
       do 45,jj = 1,6
      dsts(ii) = sig(ii)-dlda3*cc(ii,jj)*dph(jj)
45    continue
44    continue
      call crvval(crv,nnpcrv,eid(1),dlda3,dqmat(1),csp3(1))
      call crvval(crv,nnpcrv,eid(2),dlda3,dqmat(2),csp3(2))
      call crvval(crv,nnpcrv,eid(3),dlda3,dqmat(3),csp3(3))
      call crvval(crv,nnpcrv,eid(4),dlda3,dqmat(4),csp3(4))
      call crvval(crv,nnpcrv,eid(5),dlda3,dqmat(5),csp3(5))
      call crvval(crv,nnpcrv,eid(6),dlda3,dqmat(6),csp3(6))
      call crvval(crv,nnpcrv,eid(7),dlda3,dqmat(7),csp3(7))
      call crvval(crv,nnpcrv,eid(8),dlda3,dqmat(8),csp3(8))
      call crvval(crv,nnpcrv,eid(9),dlda3,dqmat(9),csp3(9))
      call crvval(crv,nnpcrv,eid(10),dlda3,dqmat(10),csp3(10))
      call crvval(crv,nnpcrv,eid(11),dlda3,dqmat(11),csp3(11))
      call crvval(crv,nnpcrv,eid(12),dlda3,dqmat(12),csp3(12))

      call tsaiwu(dsts,dqmat,tfail,tf3)
      
      if(tf3 .gt. tol1) then
         dlambda1 = dlambda3
         dlambda2 = dlambda2
         go to 100
      else
         if(tf3 .lt. tol2) then
           dlambda1 = dlambda1
           dlambda2 = dlambda3
           go to 100
          
      else
         if(tf3 .lt. tol1 .and. tf3 .gt. tol2) then
           dlambda = dlambda3
           go to 110
c
c      end if
      
110    continue
      lda = lda + dlambda
      do 30, j = 1,12
       qmat(j) = dqmat(j)
30    continue
c
      return
      end 
