      subroutine umat43 (cm,eps,sig,epsp2,hsv,dt1,capa,etype,tt, 
        1 temper,failel,crv)
!
! Copyright (C) 2005 NADER ABEDRABBO
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the redistributions of
! source code must retain the above copyright notice, this condition
! and the following disclaimer.
!
! E-Mail: aenader1@yahoo.com
! Last Modified: 09/24/2010
!
!***********************************************************************
!*                                                                     *
!*****************          DISCLAIMER          ************************
!*                                                                     *
!***********************************************************************
!*                                                                     *
!* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS *
!* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT   *
!* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS   *
!* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE      *
!* COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, *
!* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES            *
!* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR  *
!* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)  *
!* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, *
!* STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)       *
!* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED *
!* OF THE POSSIBILITY OF SUCH DAMAGE.                                  *
!*                                                                     *
!***********************************************************************
!
!**********************************************************************
! von Mises - Plasticity UMAT for von Mises Yield Function
!             using the Cutting Plane algorithm for DYNA Explicit Code
!**********************************************************************
!
!    Written BY  : Nader Abedrabbo, PhD
!    Date: September 22, 2005
!
!   Characteristics:
!           1. Plane Stress Assumption, i.e., sig13=sig23=sig33=0
!           2. Formulation based on Incremental Deformation Theory
!           3. Stress integration based on Euler Backward Theory
!           4. Using Tangent Cutting Back Theory
!           5. Rate-Independent Formulation
!
!    Note: 1) hisv(1) & hisv(2) are always fixed in all Ls-Dyna UMAT's
!             Where:
!              hisv(1) = Equivalent (Effective) Plastic Strain , epn
!              hisv(2) = Total Hydrostatic Pressure, davg
!          2) For shell elements in Dyna, sig13 and sig23 should be
!             supplied. They are calculated from Elasticity Theory.
!
!**********************************************************************
!        Definitions:
!              epn    : Equivalent plastic strain
!              epsp   : Incremental Plastic Strain
!              sigold : old stress State
!              strss  : Updated Stress
!              toli   : Tolerance for Convergence Criteria
!              hflag  : Flag for type of hardening rule (1 or 2)
!**********************************************************************
!
! variables
!
!     cm(1)   = Young's Modulus, E
!     cm(2)   = Poisson's Ratio
!     cm(3)   = Flag for type of Hardening Law
!                 1  = Power Law Hardening
!     cm(4)   = 'K' Strength Coefficient in  Hollomon Power law
!     cm(5)   = 'n' Hardening Exponent in Hollomon
!     cm(9)   = Bulk Modulus
!     cm(10)  = Shear Modulus
!     cm(11)  = Print Flag to print check message
!
!     eps(1)  = local x  strain (Incremental)
!     eps(2)  = local y  strain
!     eps(3)  = local z  strain
!     eps(4)  = local xy strain (Shear Strain eps(4) = 2*eps-xy)
!
!     sig(1)  = local x  stress (Accumulative)
!     sig(2)  = local y  stress
!     sig(3)  = local z  stress
!     sig(4)  = local xy stress
!     sig(5)  = local yz stress
!     sig(6)  = local xz stress
!
!     hisv(1) = Equivalent (Effective) Plastic Strain at t=n, epn
!     hisv(2) = Total Hydrostatic Pressure, davg
!
!     dt1  =current time step size
!     capa =reduction factor for transverse shear
!     time = current problem time.
!
!     all transformations into the element local system are
!     performed prior to entering this subroutine.  transformations
!     back to the global system are performed after exiting this
!     routine.
!
!     all history variables are initialized to zero in the input
!     phase.   initialization of history variables to nonzero values
!     may be done during the first call to this subroutine for each
!     element.
!
!     energy calculations for the dyna3d energy balance are done
!     outside this subroutine.
!
      IMPLICIT NONE
!
      include 'nlqparm'
      include 'iounits.inc'
      common/bk06/idmmy,iaddp,ifil,maxsiz,ncycle,time(2,30)
      REAL cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*), time
      INTEGER :: idmmy,iaddp,ifil,maxsiz,ncycle
      character*5 etype
      logical failel
!
      INTEGER :: i, j, iter, hflag, PR
      REAL    :: epn, ams1, ams2, ams3, G, G2
      REAL    :: capa, tt, temper, dt1, epsp2
      REAL    :: ym, nu, Kh, nh, davg
      REAL    :: toli, f1, f2, f3
      REAL, DIMENSION (3,3) :: cc
      REAL, DIMENSION (3)   :: Strss, sigold, deps, epsp
!
      common/khard/ams1,ams2,ams3,ym,nu,hflag
!
!    Assign Variables to Common Names and Recall Values of History Variables
!
      ym     = cm(1)
      nu     = cm(2)
      hflag  = NINT(cm(3))
      Kh     = cm(4)
      nh     = cm(5)
      PR     = NINT(cm(11))  ! Message Print Flag
!
      toli   = 1.0E-6        ! Iteration tolerance value
!
!    Initialize matrices and arrays
!
      do i=1,3
        do j=1,3
          cc(i,j) = 0.0
        enddo !j
      enddo !i
!
!    Compute shear modulus, G
!
      G     = ym / (2.0*(1.0 + nu))  ! Modulus of Rigidity
      G2    = ym / (1.0 - nu*nu)     ! E/(1-v²)
!
!    Assemble Elastic Moduli
!
      cc(1,1) = G2
      cc(1,2) = nu * G2
      cc(2,1) = nu * G2
      cc(2,2) = G2
      cc(3,3) = G       ! Because Shear strain is 2*eps_xy else cc(4,4)=2*G
!
!    Assign values for hardening coefficients
!    Depends on hflag-Value:
!         1 : Power Law Hardening
!
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
!
!    Assign inputs for main code
!
      deps(1)   = eps(1)
      deps(2)   = eps(2)
      deps(3)   = eps(4)
!
      sigold(1) = sig(1)
      sigold(2) = sig(2)
      sigold(3) = sig(4)
!
      epn       = hsv(1)
!
! *********************************************************************
! **     Compute Stress based on Incremental Deformation Theory      **
! *********************************************************************
!
      call jaman(cc,deps,sigold,toli,epsp,epn,iter,Strss)
!
! *********************************************************************
!
!    Now to update the stress to the new updated values
!
!    Update Stresses for output
!
      sig(1)   = Strss(1)
      sig(2)   = Strss(2)
      sig(3)   = 0.0
      sig(4)   = Strss(3)
      sig(5)   = sig(5) + capa*G*eps(5)
      sig(6)   = sig(6) + capa*G*eps(6)
!
!    Calculate Total Hydrostatic Pressure
!
      davg     = - (1.0/3.0) * (sig(1)+sig(2)+sig(3))
!
!    Calculate Strain Increment in Normal Direction (z)
!
      f1       = epsp(1) + epsp(2)
      f2       = nu*(eps(1)+eps(2))
      f3       = 2.0*nu*f1
!
      eps(3)   = (f1+f2-f3)/(nu-1.0)
      epsp2    = epn                      ! effective plastic strain
!                                         ! Reported back to model
!
!    Update History Variables
!
      hsv(1)  = epn
      hsv(2)  = davg
      hsv(3)  = iter
!
!*******************************************************************
!
!    Print a check statement to test code usage. IFF (PR=1)
!
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
!
!*******************************************************************
!
      return
      end
!
! *********************************************************************
! *********************************************************************
!
      subroutine jaman(cc,deps,sigold,toli,epsp,epn,iter,Strss)
!
! *********************************************************************
!
!  BY  : Nader Abedrabbo
!  Date: March 25, 2003
!
!  ********************************************************************
!  This is the main program for correctly scaling back the stress to the
!  yield surface. Inputs for this program are the yield function and its
!  derivative and the hardening function. The Cutting Plane Algorithm is
!  used to scale back the stresses to the yield surface. It gives back
!  the value of the updated stress, the equivalent plastic strain, the
!  plastic strain value and the number of iterations to the main program.
!  ********************************************************************
!  The cutting plane algorithm is based on the following paper:
!  M. Ortiz and J.C. Simo, "An Analysis Of A New Class Of Integration
!  Algorithm For Elastoplastic Constitutive Relations", Inter. Journal
!  For Numerical Methods In Engineering, Vol.23, pp. 353-366, 1986.
!  ********************************************************************
!
      REAL, Dimension (3,3):: cc
      REAL, Dimension (3)  :: sigt, sigold, deps, epsp, epsp2
      REAL, Dimension (3)  :: Strss, dphi, epspdot
      REAL :: ep, epn, ams1, ams2, ams3
      REAL :: f1, f2, f3, f4, ff1, epdot
      REAL :: phi, hard2d, yldf, yld2, yldn, dlam, fac, kprim
      REAL :: toli
!
      INTEGER :: i, j, k, ii, jj, ik, hflag, iter
!
      common/khard/ams1,ams2,ams3,ym,nu,hflag
!
      ep    = epn
!
      do i=1,3
          epsp2(i)    = 0.0
          sigt(i)     = sigold(i)
      enddo !i
!
!    Compute updated total trial stress
!
      do i=1,3
          do j=1,3
              sigt(i) = sigt(i) + cc(i,j)*deps(j)
          enddo ! j
      enddo ! i
!
!    Calculate the trial yield function (at k)
!
      call vmises(sigt,phi,dphi)
!
!    Calculate Hardening part (center of yield surface)
!
      call hardfunc(ep,hard2d,kprim)
!
!    Now to check for yield condition : Elastic or Plastic ?
!
      yldf = phi - hard2d
!
      if (phi.gt.hard2d) then         ! Plastic State check
          iter  = 0
          dlam  = 0.0
!
  100     fac   = kprim
          do ii=1,3
             do jj=1,3
                  fac = fac + dphi(ii)*cc(ii,jj)*dphi(jj)
             enddo !jj
          enddo !ii
          yld2 = phi - hard2d
          dlam = yld2 / fac
          ep   = ep + dlam                      ! new equivalent plastic strain
          do ik = 1,3
               epspdot(ik) = dlam * dphi(ik)         ! Plastic Strain increment
               epsp2(ik)   = epsp2(ik) + epspdot(ik) ! Total Incremental plastic strain
          enddo !ik
!
        ! Evaluate new stress values (cutting back)
          do ii=1,3
             do jj=1,3
                  sigt(ii) = sigt(ii) - dlam*cc(ii,jj)*dphi(jj)
             enddo !jj
          enddo !ii
!
        ! Evaluate new yield function according to new updated stress
          call vmises(sigt,phi,dphi)
        ! Calculate Hardening value from new equivalent plastic strain
          call hardfunc(ep,hard2d,kprim)
!
!           Check for Convergence
!
          yldn    = phi - hard2d
          iter    = iter+1
!
          if (iter.gt.50) go to 110
          if (abs(yldn).gt.toli) go to 100
          go to 120
!
  110     write(*,*)"*Warning* Plasticity Algorithm Failed To Converge"
          write(*,1000)'|yldn| =',abs(yldn),'toli= ',toli,'Iter# =',iter
!
  120     continue
!
      else                    ! Elastic state
          iter    = 0
      endif                  ! end for Plastic Check
!
!    Now to report the updated values for the stress, plastic strain,
!    equivalent plastic strain and iteration number
!
      epn   = ep
      do k = 1,3
           Strss(k) = sigt(k)
           epsp(k)  = epsp2(k)
      enddo !k
!
1000  Format (A8,f12.8,A10,f12.8,A12,I5)
      return
      end
!
! *********************************************************************
! ****     Sub-Function Calls for von-Mises and its derivatives    ****
! *********************************************************************
!
      subroutine  vmises(sig,phi,dphi)
!
      implicit none
!
      REAL, Dimension (3)  :: sig, dphi
      REAL :: t1, t2, t3, phi
      INTEGER :: i
!
! This subfunction calculates the yield function
! and it derivative for the von Mises function
!
      do i=1,3
          dphi(i) = 0.0
      enddo
!
      t1      = sig(1)
      t2      = sig(2)
      t3      = sig(3)
      phi     = sqrt(t1**2+t2**2-t1*t2+3.0*t3**2)
!
      dphi(1) = (2.0*sig(1) - sig(2)) / (2.0*phi)
      dphi(2) = (2.0*sig(2) - sig(1)) / (2.0*phi)
      dphi(3) = (3.0*sig(3)) / phi
!
      return
      end
!
! *********************************************************************
! ***************            Hardening Rules            ***************
! *********************************************************************
!
      subroutine  hardfunc(ep,hard2d,kprim)
!
      implicit none
!
      REAL :: ams1, ams2, ams3, fact, qh, ym, nu
      REAL :: ep, ee1, ee2, hard2d, kprim
!
      INTEGER :: hflag
!
      common/khard/ams1,ams2,ams3,ym,nu,hflag
!
!    This Calculates the hardening part; (center of yield surface),
!    using Hollomon Power Law
!
      if (hflag.eq.1) then
        fact     = 1.0/(ams2-1.0)
        qh       = (ym/ams1)**fact
        ee1      = (qh+ep)**ams2
        ee2      = (qh+ep)**(ams2-1.0)
        hard2d   = ams1 * ee1
        ! Isotropic hardening part dK/dep
        kprim    = ams1 * ams2 * ee2
      endif
!
      return
      end
!
! *********************************************************************
! **********       End OF NADER ABEDRABBO's PROGRAM          **********
! *********************************************************************
!
