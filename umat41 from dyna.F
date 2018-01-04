      subroutine umat41 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
      1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
c
c
c isotropic elastic material (sample user subroutine)
c Variables
c
c cm(1)=first material constant, here young's modulus
c cm(2)=second material constant, here poisson's ratio
c .
c .
c .
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
c eq."sld2d" for shell forms 13, 14, and 15 (2D solids)
c eq."shl_t" for shell forms 25, 26, and 27 (shells with thickness
c stretch)
c eq."shell" for all other shell elements plus thick shell forms 1
c and 2
c eq."tshel" for thick shell forms 3 and 5
c eq."hbeam" for beam element forms 1 and 11
c eq."tbeam" for beam element form 3 (truss)
c eq."dbeam" for beam element form 6 (discrete)
c eq."beam " for all other beam elements
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
c nnpcrv=# of discretization points per crv()
c
c cma=additional memory for material data defined by LMCA at
c 6th field of 2nd crad of *DATA_USER_DEFINED
c
c elsiz=characteristic element size
c
c idele=element id
c
c reject (implicit only) = set to .true. if this implicit iterate is
c to be rejected for some reason
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
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      logical failel,reject
      character*5 etype
c
      if (ncycle.eq.1) then
        if (cm(16).ne.1234567) then
        call usermsg('mat41')
        endif
      endif
c
c compute shear modulus, g
c
      g2 =abs(cm(1))/(1.+cm(2))
      g =.5*g2
c
      if (etype.eq.'solid'.or.etype.eq.'shl_t'.or.
      1 etype.eq.'sld2d'.or.etype.eq.'tshel') then
      if (cm(16).eq.1234567) then
      call mitfail3d(cm,eps,sig,epsp,hsv,dt1,capa,failel,tt,crv)
      else
      if (.not.failel) then
      davg=(-eps(1)-eps(2)-eps(3))/3.
      p=-davg*abs(cm(1))/(1.-2.*cm(2))
      sig(1)=sig(1)+p+g2*(eps(1)+davg)
      sig(2)=sig(2)+p+g2*(eps(2)+davg)
      sig(3)=sig(3)+p+g2*(eps(3)+davg)
      sig(4)=sig(4)+g*eps(4)
      sig(5)=sig(5)+g*eps(5)
      sig(6)=sig(6)+g*eps(6)
      if (cm(1).lt.0.) then
        if (sig(1).gt.cm(5)) failel=.true.
        endif
      endif
      end if
c
      else if (etype.eq.'shell') then
      if (cm(16).eq.1234567) then
      call mitfailure(cm,eps,sig,epsp,hsv,dt1,capa,failel,tt,crv)
      else
      if (.not.failel) then
      gc =capa*g
      q1 =abs(cm(1))*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
      q3 =1./(q1+g2)
      eps(3)=-q1*(eps(1)+eps(2))*q3
      davg =(-eps(1)-eps(2)-eps(3))/3.
       p =-davg*abs(cm(1))/(1.-2.*cm(2))
      sig(1)=sig(1)+p+g2*(eps(1)+davg)
      sig(2)=sig(2)+p+g2*(eps(2)+davg)
      sig(3)=0.0
      sig(4)=sig(4)+g *eps(4)
      sig(5)=sig(5)+gc*eps(5)
      sig(6)=sig(6)+gc*eps(6)
      if (cm(1).lt.0.) then
        if (sig(1).gt.cm(5)) failel=.true.
        endif
      endif
      end if
      elseif (etype.eq.'beam ' ) then
      q1 =cm(1)*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
      q3 =q1+2.0*g
      gc =capa*g
      deti =1./(q3*q3-q1*q1)
      c22i = q3*deti
      c23i =-q1*deti
      fac =(c22i+c23i)*q1
       eps(2)=-eps(1)*fac-sig(2)*c22i-sig(3)*c23i
       eps(3)=-eps(1)*fac-sig(2)*c23i-sig(3)*c22i
       davg =(-eps(1)-eps(2)-eps(3))/3.
       p =-davg*cm(1)/(1.-2.*cm(2))
      sig(1)=sig(1)+p+g2*(eps(1)+davg)
      sig(2)=0.0
      sig(3)=0.0
      sig(4)=sig(4)+gc*eps(4)
      sig(5)=0.0
      sig(6)=sig(6)+gc*eps(6)
c
      elseif (etype.eq.'tbeam') then
      q1 =cm(1)*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
      q3 =q1+2.0*g
      deti =1./(q3*q3-q1*q1)
      c22i = q3*deti
      c23i =-q1*deti
      fac =(c22i+c23i)*q1
      eps(2)=-eps(1)*fac
      eps(3)=-eps(1)*fac
      davg =(-eps(1)-eps(2)-eps(3))/3.
      p =-davg*cm(1)/(1.-2.*cm(2))
      sig(1)=sig(1)+p+g2*(eps(1)+davg)
      sig(2)=0.0
      sig(3)=0.0
c
      else
c      write(iotty,10) etype
c      write(iohsp,10) etype
c      write(iomsg,10) etype
c      call adios(TC_ERROR)
      cerdat(1)=etype
      call lsmsg(3,MSG_SOL+1150,ioall,ierdat,rerdat,cerdat,0)
      endif
c
c10    format(/
c       1 ' *** Error element type ',a,' can not be',
c      2 ' run with the current material model.')
      return
      end
