      subroutine umat43(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt, 
     1 temper,failel,crv,nnpcrv)
c
c******************************************************************
c| N.L.Vishnuvardhan Raju                                            |
c| ------------------------------------------------------------ |
c| Copyright @ NMCAD Lab                                        |
c| All rights reserved                                          |
c******************************************************************
c
c
c Piecewise linear plasticity with load curve acess.
c Variables
c
c cm(1) = young's modulus.
c cm(2) = Poisson's ratio.
c cm(3) = initial yield limit.
c cm(4) = yield stress vs plastic strain curve id.
c 
c 
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
c
c etype:
c eq."solid" for solid elements
c eq."shell" for shell elements 
c 
c
c tt=current problem time.
c
c temper=current temperature
c
c failel=flag for failure, set to .true. to fail an integration point,
c if .true. on input the integration point has failed earlier
c
c crv=array representation of curves in keyword deck
c
c nnpcrv= no. of discretization points per crv()
c
c cma=additional memory for material data defined by LMCA at
c 6th field of 2nd crad of *DATA_USER_DEFINED
c
c elsiz=characteristic element size
c
c idele=element id
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
      dimension cm(*),eps(*),sig(*),hsv(*)
      dimension  crv(lq1,2,*)
      real :: cc,g,nu,ym,sigold,deps,strss,neps,eid
      dimension  cc(6,6),sigold(6),deps(6),neps(6),strss(6)
      real :: ylds,toli,epsp,etan
      integer  i,j
      integer  :: nnpcrv(*)
      logical failel
      character*5 etype
      common/coeff/ylds,ym,nu,g,etan,eid,lq1
c Assigning the variables   	  
       ym   = cm(1)
       nu   = cm(2)
       ylds = cm(3)
       etan = cm(4)
       eid  = cm(5)
       toli = 1.0E-06
       
c Populating the Constitutive matrix
        
          do 10, i=1,6
          do 20, j=1,6
            cc(i,j) = 0.0
20       continue		 
10       continue		 

       prdn  = 1.0-nu-2.0*nu*nu         
       prnr1 = ym*(1.0-nu)
       prnr2 = ym*nu
       g     = (0.5*ym)/(1.0+nu)
		  
       cc(1,1) = prnr1/prdn
       cc(2,2) = prnr1/prdn
       cc(3,3) = prnr1/prdn
       cc(1,2) = prnr2/prdn
       cc(1,3) = prnr2/prdn
       cc(2,1) = prnr2/prdn
       cc(2,3) = prnr2/prdn
       cc(3,1) = prnr2/prdn
       cc(3,2) = prnr2/prdn
       cc(4,4) = g
       cc(5,5) = g
       cc(6,6) = g
       if(etype .eq. 'solid') then
c
c   Assigining inputs for main code

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

      epn = hsv(1)
      
c *********************************************************************
c **     Compute Stress based on Piecewise Linear Plasticity         **
c *********************************************************************
c
      call plpuc(cc,deps,sigold,toli,crv,nnpcrv,epn,strss,neps)
c
c *********************************************************************
c
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
      davg     =  -(1.0/3.0)*(sig(1)+sig(2)+sig(3))
c

c   effective plastic strain Reported back to model 
       epsp    = epn
c    Update History Variables
c
      hsv(1)  = epn
      hsv(3) = hsv(3)+(1.0/3.0)*(eps(1)+eps(2)+eps(3))
      hsv(2)  = davg
c
      else
      write(*,*) 'element type mismatch'
      end if
       end

c *********************************************************************
c
       subroutine plpuc(cc,deps,sigold,toli,crv,nnpcrv,epn,strss,neps)
c
c *********************************************************************
       dimension crv(lq1,2,*)
       real ::   cc,deps,sigold,epn,strss,neps
       dimension cc(6,6),deps(6),sigold(6),neps(6),strss(6)
       real :: eid
       integer nnpcrv(*)
       real toli,ylds,ym,nu,phi, dphi, epsp2, epspdot,g,sigt
       dimension  dphi(6),epsp2(6),epspdot(6),sigt(6)
       integer i,j,ii,jj
       real :: cyst,plspe,etan
       common/coeff/ylds,ym,nu,g,etan,eid,lq1
       
       ep = epn

       do i=1,6
          epsp2(i)    = 0.0
          sigt(i)     = sigold(i)
      enddo  
c
c    Compute updated total trial stress
c
      do i=1,6
          do j=1,6
              sigt(i) = sigt(i) + cc(i,j)*deps(j)
          enddo            
      enddo                
c    Calculate the trial yield function (at k)
c
      call vmises(sigt,phi,dphi)
c
c    Calculate Hardening part

      call hardfunc(ep,crv,nnpcrv,cyst,plspe)
c
c    Now to check for yield condition : Elastic or Plastic ?
c
      yldf = ((1.0/3.0)*(phi)**2) - ((1.0/3.0)*cyst**2)
c
      if (yldf .gt. 0.0) then         ! Plastic State check
          dlam  = 0.0

100    continue

       dlam =  (phi - cyst)/(3*g+plspe)
       ep = ep + dlam
       do i = 1,6
               epspdot(i) = dlam * dphi(i)         
               epsp2(i)   = epsp2(i) + epspdot(i) 
          enddo  
c
c     Evaluate new stress values (cutting back)
          do ii=1,6
              sigt(ii) = (cyst/phi)*sigt(ii)
          enddo 
c
c    Evaluate new yield function according to new updated stress
          call vmises(sigt,phi,dphi)
c    Calculate Hardening value from new equivalent plastic strain
          call hardfunc(ep,crv,nnpcrv,cyst,plspe)
c
c           Check for Convergence
c
          yldn    = ((1.0/3.0)*(phi)**2) - ((1.0/3.0)*cyst**2)
c
          if (abs(yldn).gt.toli) go to 100
          else
          go to 120
c
120     continue
c
         endif                  ! end for Plastic Check
c
c    Now to report the updated values for the stress, plastic strain,
c    equivalent plastic strain and iteration number
c
      epn   = ep
      
      do i = 1,6
           strss(i) = sigt(i)
           neps(i)  = epsp2(i)
      enddo 
c
      return
      end
c
c *********************************************************************
c ****     Sub-Function Calls for von-Mises function               ****
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

c*********************************************************************
c     Calculating the hardening part.
c
c*********************************************************************
       subroutine   hardfunc(ep,crv,nnpcrv,cyst,plspe)
c
c
       dimension  crv(lq1,2,*)
       REAL :: g,etan,crv,eid 
       integer  nnpcrv(*)
       REAL :: ep, cyst,plspe,ylds,ym,nu
       common/coeff/ylds,ym,nu,g,etan,eid,lq1
c
c    This Calculates the hardening part;
c
       if(eid .eq. 0.0) then 
       plspe = (ym*etan)/(ym-etan)
       cyst = ylds + plspe*ep
       else 
       call crvval(crv,nnpcrv,eid,ep,cyst,plspe)
c       
       end if
      return
      end
