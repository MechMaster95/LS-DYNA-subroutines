      subroutine umat43 (cm,eps,sig,epsp2,hsv,dt1,capa,etype,tt, 
        1 temper,failel,crv)
c
c Written by K.L.BHARATH REDDY.
c Copywright  @  NMCAD Lab.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the redistributions of
c source code must retain the above copyright notice, this condition
c and the following disclaimer.
c
c E-Mail: bharathrlk@gmail.com
c Last Modified: 16/06/2017
c
c***********************************************************************
c*                                                                     *
c*****************          DISCLAIMER          ************************
c*                                                                     *
c***********************************************************************
c
c
c**********************************************************************
c von Mises - Plasticity UMAT for von Mises Yield Function for solid
c             elements with failure limit.
c             using the Cutting Plane algorithm for DYNA Explicit Code
c**********************************************************************
c
c    Written BY  : K.L.BHARATH REDDY
c    Date: 16/06/2017
c
c   Characteristics:
c           1. Implemented for solid elements considering volumetric
c               strain as effectrive strain for failure.
c           2. Formulation based on Incremental Deformation Theory
c           3. Stress integration based on Euler Backward Theory
c           4. Using Tangent Cutting Back Theory
c           5. Rate-Independent Formulation
c
c    Note: 1) hisv(1) & hisv(2) are used
c             Where:
c              hisv(1) = Equivalent (Effective) Plastic Strain, epn
c              hisv(2) = Total Hydrostatic Pressure, davg
c          2) For solid elements.
c
c**********************************************************************
c        Definitions:
c              epn    : Equivalent plastic strain
c              epsp   : Incremental Plastic Strain
c              sigold : old stress State
c              strss  : Updated Stress
c              toli   : Tolerance for Convergence Criteria
c              hflag  : Flag for type of hardening rule (1 or 2)
c**********************************************************************
c
c variables
c
c     cm(1)   = Young's Modulus, E
c     cm(2)   = Poisson's Ratio
c     cm(3)   = Flag for type of Hardening Law
c                 1  = Power Law Hardening
c     cm(4)   = 'K' Strength Coefficient in  Hollomon Power law
c     cm(5)   = 'n' Hardening Exponent in Hollomon
c     cm(6)   = Ultimate strain for failure.
c     cm(9)   = Bulk Modulus
c     cm(10)  = Shear Modulus
c     cm(11)  = Print Flag to print check message
c
c     eps(1)  = local x  strain (Incremental)
c     eps(2)  = local y  strain
c     eps(3)  = local z  strain
c     eps(4)  = local xy strain
c     eps(5)  = local yz strain
c     eps(6)  = local zx strain
c
c     sig(1)  = local x  stress (Accumulative)
c     sig(2)  = local y  stress
c     sig(3)  = local z  stress
c     sig(4)  = local xy stress
c     sig(5)  = local yz stress
c     sig(6)  = local xz stress
c
c     hisv(1) = Equivalent (Effective) Plastic Strain at t=n, epn
c     hisv(2) = Total Hydrostatic Pressure, davg
c
c     dt1  =current time step size
c     capa =reduction factor for transverse shear
c     time = current problem time.
c
c     all transformations into the element local system are
c     performed prior to entering this subroutine.  transformations
c     back to the global system are performed after exiting this
c     routine.
c
c     all history variables are initialized to zero in the input
c     phase.   initialization of history variables to nonzero values
c     may be done during the first call to this subroutine for each
c     element.
c
c     energy calculations for the dyna3d energy balance are done
c     outside this subroutine.
c
      IMPLICIT NONE
c
      include 'nlqparm'
      include 'iounits.inc'
      common/bk06/idmmy,iaddp,ifil,maxsiz,ncycle,time(2,30)
      REAL cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*), time
      INTEGER :: idmmy,iaddp,ifil,maxsiz,ncycle
      character*5 etype
      logical failel,eps
c
      INTEGER :: i, j, iter, hflag, PR
      REAL    :: epn, ams1, ams2, ams3, G, G1, G2
      REAL    :: capa, tt, temper, dt1, epsp2
      REAL    :: ym, nu, Kh, nh, davg,stn1
      REAL    :: toli cc, f1, f2, f3
      REAL, DIMENSION (6,6) :: cc
      REAL, DIMENSION (6)   :: Strss, sigold, deps, epsp
      logical failel,efs
c
      common/khard/ams1,ams2,ams3,ym,nu,hflag
c
cAssign Variables to Common Names and Recall Values of History Variables
c
      ym     = cm(1)
      nu     = cm(2)
      hflag  = NINT(cm(3))
      Kh     = cm(4)
      nh     = cm(5)
      PR     = NINT(cm(11))  c Message Print Flag
      stn1   = cm(6)
c
      toli   = 1.0E-6        c Iteration tolerance value
c
c    Initialize matrices and arrays
c
      do i=1,6
        do j=1,6
          cc(i,j) = 0.0
        enddo cj
      enddo ci
c
c    Compute shear modulus, G
c
      G     = ym / (2.0*(1.0 + nu))  c Modulus of Rigidity
      G1    = (ym*(1.0-nu)) / (1.0 - nu-2.0*nu*nu)   c E(1-v)/(1-v-2v²)
      G2    = (ym*nu) / (1.0 - nu-2.0*nu*nu)   c Ev/(1-v-2v²)
c
c    Assemble Elastic Moduli
c
      cc(1,1) = G1
      cc(1,2) = G2
      cc(1,3) = G2
      cc(2,1) = G2
      cc(2,2) = G1
      cc(2,3) = G2    
      cc(3,1) = G2
      cc(3,2) = G2
      cc(3,3) = G1
      cc(4,4) = G
      cc(5,5) = G
      cc(6,6) = G
c
c    Assign values for hardening coefficients
c    Depends on hflag-Value:
c         1 : Power Law Hardening
c
      if (hflag.eq.1) then
          ams1 = Kh
          ams2 = nh
          ams3 = 0.0
      else
       write(*,*)
       write(*,*)
       write(*,*) '***************************************************'
       write(*,*) '*                                                 *'
       write(*,*) '*                  E R R O R                      *'
       write(*,*) '*                                                 *'
       write(*,*) '***************************************************'
       write(*,*) '*                                                 *'
       write(*,*) '*                                                 *'
       write(*,*) '* FLAG ERROR: Check Flag Number For Hardening     *'
       write(*,*) '*             Model Type. In Input File, cm(3)    *'
       write(*,*) '*             Should be:                          *'
       write(*,*) '*                                                 *'
       write(*,*) '*             1: Power Law                        *'
       write(*,*) '*                                                 *'
       write(*,*) '***************************************************'
       write(*,*)
       write(*,*)
       STOP
      endif
c
c    Assign inputs for main code
c
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
      
c
      epn       = hsv(1)
c
c *********************************************************************
c **     Compute Stress based on Incremental Deformation Theory      **
c *********************************************************************
c
      call idefth(cc,deps,sigold,toli,epsp,epn,iter,Strss)
c
c *********************************************************************
c
c    Now to update the stress to the new updated values
c
c    Update Stresses for output
c
      sig(1)   = Strss(1)
      sig(2)   = Strss(2)
      sig(3)   = Strss(3)
      sig(4)   = Strss(4)
      sig(5)   = Strss(5)
      sig(6)   = Strss(6)
c
c    Calculate Total Hydrostatic Pressure
c
      davg     = - (1.0/3.0) * (sig(1)+sig(2)+sig(3))
c
c    Calculate Strain Increment in Normal Direction (z)
c
c      f1       = epsp(1) + epsp(2)
c      f2       = nu*(eps(1)+eps(2))
c      f3       = 2.0*nu*f1
c
c      eps(3)   = (f1+f2-f3)/(nu-1.0)
      epsp2    = epn                      c effective plastic strain
c                                         c Reported back to model
       
       call epsfail(epsp2,stn1,efs)
       if(efs.eqv..true.) then
       failel = 1
       sig = 0.0
       end if

c
c    Update History Variables
c
      hsv(1)  = epn
      hsv(2)  = davg
      hsv(3)  = iter
c
c*******************************************************************
c
c    Print a check statement to test code usage. IFF (PR=1)
c
      if (PR.ne.1) PR = 0
      if (PR.eq.1) then
       write(*,*)
       write(*,*)
       write(*,*)'***********************************'
       write(*,*)'*                                 *'
       write(*,*)'*      You are Using UMAT43       *'
       write(*,*)'*         Using von Mises         *'
       write(*,*)'*           LS-Dyna 971           *'
       write(*,*)'*                                 *'
       write(*,*)'***********************************'
       write(*,*)
       write(*,*) '  Toli           = ', toli
       write(*,*) '  Youngs Modulus = ', ym
       write(*,*) '  HFLAG          = ', hflag
       write(*,*) '  ams1           = ', ams1
       write(*,*) '  ams2           = ', ams2
       write(*,*) '  ams3           = ', ams3
       write(*,*)
       write(*,*)'*****************************'
       write(*,*)
       write(*,*)
       STOP
      endif
c
c*******************************************************************
c
      return
      end
c
c *********************************************************************
c *********************************************************************
c
      subroutine idefth(cc,deps,sigold,toli,epsp,epn,iter,Strss)
c
c *********************************************************************
c
c  BY  : K.L.Bharath Reddy
c  Date: 16/06/2017
c
c  ********************************************************************
c  This is the main program for correctly scaling back the stress to the
c  yield surface. Inputs for this program are the yield function and its
c  derivative and the hardening function. The Cutting Plane Algorithm is
c  used to scale back the stresses to the yield surface. It gives back
c  the value of the updated stress, the equivalent plastic strain, the
c  plastic strain value and the number of iterations to the main program.
c  ********************************************************************
c  The cutting plane algorithm is based on the following paper:
c  M. Ortiz and J.C. Simo, "An Analysis Of A New Class Of Integration
c  Algorithm For Elastoplastic Constitutive Relations", Inter. Journal
c  For Numerical Methods In Engineering, Vol.23, pp. 353-366, 1986.
c  ********************************************************************
c
      REAL, Dimension (6,6):: cc
      REAL, Dimension (6)  :: sigt, sigold, deps, epsp, epsp2
      REAL, Dimension (6)  :: Strss, dphi, epspdot
      REAL :: ep, epn, ams1, ams2, ams3
      REAL :: f1, f2, f3, f4, ff1, epdot
      REAL :: phi, hard2d, yldf, yld2, yldn, dlam, fac, kprim
      REAL :: toli
c
      INTEGER :: i, j, k, ii, jj, ik, hflag, iter
c
      common/khard/ams1,ams2,ams3,ym,nu,hflag
c
      ep    = epn
c
      do i=1,6
          epsp2(i)    = 0.0
          sigt(i)     = sigold(i)
      enddo ci
c
c    Compute updated total trial stress
c
      do i=1,6
          do j=1,6
              sigt(i) = sigt(i) + cc(i,j)*deps(j)
          enddo c j
      enddo c i
c
c    Calculate the trial yield function (at k)
c
      call vmises(sigt,phi,dphi)
c
c    Calculate Hardening part (center of yield surface)
c
      call hardfunc(ep,hard2d,kprim)
c
c    Now to check for yield condition : Elastic or Plastic ?
c
      yldf = phi - hard2d
c
      if (phi.gt.hard2d) then         c Plastic State check
          iter  = 0
          dlam  = 0.0
c
  100     fac   = kprim
          do ii=1,6
             do jj=1,6
                  fac = fac + dphi(ii)*cc(ii,jj)*dphi(jj)
             enddo cjj
          enddo cii
          yld2 = phi - hard2d
          dlam = yld2 / fac
          ep   = ep + dlam               c new equivalent plastic strain
          do ik = 1,6
               epspdot(ik) = dlam * dphi(ik)         ! Plastic Strain increment
               epsp2(ik)   = epsp2(ik) + epspdot(ik) ! Total Incremental plastic strain
          enddo cik
c
c     Evaluate new stress values (cutting back)
          do ii=1,6
             do jj=1,6
                  sigt(ii) = sigt(ii) - dlam*cc(ii,jj)*dphi(jj)
             enddo cjj
          enddo cii
c
c    Evaluate new yield function according to new updated stress
          call vmises(sigt,phi,dphi)
c    Calculate Hardening value from new equivalent plastic strain
          call hardfunc(ep,hard2d,kprim)
c
c           Check for Convergence
c
          yldn    = phi - hard2d
          iter    = iter+1
c
          if (iter.gt.50) go to 110
          if (abs(yldn).gt.toli) go to 100
          go to 120
c
  110     write(*,*)"*Warning* Plasticity Algorithm Failed To Converge"
          write(*,1000)'|yldn| =',abs(yldn),'toli= ',toli,'Iter# =',iter
c
  120     continue
c
      else                    c Elastic state
          iter    = 0
      endif                  c end for Plastic Check
c
c    Now to report the updated values for the stress, plastic strain,
c    equivalent plastic strain and iteration number
c
      epn   = ep
      do k = 1,6
           Strss(k) = sigt(k)
           epsp(k)  = epsp2(k)
      enddo ck
c
1000  Format (A8,f12.8,A10,f12.8,A12,I5)
      return
      end
c
c *********************************************************************
c ****     Sub-Function Calls for von-Mises and its derivatives    ****
c *********************************************************************
c
      subroutine  vmises(sig,phi,dphi)
c
      implicit none
c
      REAL, Dimension (6)  :: sig, dphi
      REAL :: t1, t2, t3,t4,t5,t6, phi, s12,s23,s31,s123
      INTEGER :: i
c
c This subfunction calculates the yield function
c and it derivative for the von Mises function
c
      do i=1,6
          dphi(i) = 0.0
      enddo
c
      t1      = sig(1)
      t2      = sig(2)
      t3      = sig(3)
      t4      = sig(4)
      t5      = sig(5)
      t6      = sig(6)
      
      s12     = (t1-t2)**2
      s23     = (t2-t3)**2
      s31     = (t3-t1)**2
      s123    = 6*(t4**2+t5**2+t6**2)
        
      phi     = sqrt((s12+s23+s31+s123)/2.0)
c
      dphi(1) = (2.0*sig(1)-(sig(2)+sig(3)))/(2.0*phi)
      dphi(2) = (2.0*sig(2)-(sig(1)+sig(3)))/(2.0*phi)
      dphi(3) = (2.0*sig(3)-(sig(1)+sig(2)))/(2.0*phi)
      dphi(4) = (3.0*sig(4)) / phi
      dphi(5) = (3.0*sig(5)) / phi
      dphi(6) = (3.0*sig(6)) / phi
c
      return
      end
c
c *********************************************************************
c ***************            Hardening Rules            ***************
c *********************************************************************
c
      subroutine  hardfunc(ep,hard2d,kprim)
c
      implicit none
c
      REAL :: ams1, ams2, ams3, fact, qh, ym, nu
      REAL :: ep, ee1, ee2, hard2d, kprim
c
      INTEGER :: hflag
c
      common/khard/ams1,ams2,ams3,ym,nu,hflag
c
c    This Calculates the hardening part; (center of yield surface),
c    using Hollomon Power Law
c
      if (hflag.eq.1) then
        fact     = 1.0/(ams2-1.0)
        qh       = (ym/ams1)**fact
        ee1      = (qh+ep)**ams2
        ee2      = (qh+ep)**(ams2-1.0)
        hard2d   = ams1 * ee1
        c Isotropic hardening part dK/dep
        kprim    = ams1 * ams2 * ee2
      endif
c
      return
      end
c
c *********************************************************************
c **********       End OF The Algorithm                  **********
c *********************************************************************
c
c
       subroutine epsfail(epn,stn1,efs)
c    For failing the material if strain is given
      implicit none
      real  epn,stn1,tole
      logical efs
c  Checking the tolerance limit on failure strain       
       tole = 1.0e-06
       IF((ABS(stn1)-ABS(epn)) .LT. TOL) THEN
       efs =.TRUE.
       ELSE
       efs = .FALSE.
       END IF
       END       
