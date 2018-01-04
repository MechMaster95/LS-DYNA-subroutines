      subroutine umat47(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt)
c isotropic elastic material (sample user subroutine)
c Variables
c WRITTEN BY Nadimpalli L Vishnuvardhan Raju
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
c hisv(1) = Equivalent (Effective) Plastic Strain at t=n, epn
c hisv(2) = Total Hydrostatic Pressure, davg
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
c
c eq."shell" for all other shell elements plus thick shell forms 1
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

      dimension cm(*),eps(*),sig(*),hsv(*)
      real tt,dt1,efs
      real capa,epsp
      character*5 etype
      
      
      g2 =abs(cm(1))/(1.+cm(2))
      g =.5*g2
      
      if (etype.eq.'solid') THEN
      
      davg=(-eps(1)-eps(2)-eps(3))/3.
      p=-davg*abs(cm(1))/(1.-2.*cm(2))
      sig(1)=sig(1)+p+g2*(eps(1)+davg)
      sig(2)=sig(2)+p+g2*(eps(2)+davg)
      sig(3)=sig(3)+p+g2*(eps(3)+davg)
      sig(4)=sig(4)+g*eps(4)
      sig(5)=sig(5)+g*eps(5)
      sig(6)=sig(6)+g*eps(6)
      write(*,*) 'current time =',tt
      write(*,*) 'increment time =',dt1
      efs = -(eps(1)+eps(2)+eps(3))
      hsv(2)=-(1.0/3.0)*(sig(1)+sig(2)+sig(3))
      hsv(1) = hsv(1)+efs
      hsv(3) = hsv(3)+1
      write(*,*) 'total calling',hsv(3)
      epsp = 0.0
      
      else if (etype.eq.'shell') then
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
      write(*,*) 'current time =',tt
      write(*,*) 'increment time =',dt1
      efs = -(eps(1)+eps(2)+eps(3))
      hsv(2)=-(1.0/3.0)*(sig(1)+sig(2)+sig(3))
      hsv(1) = hsv(1)+efs
      hsv(3) = hsv(3)+1
      write(*,*) 'total calling',hsv(3)
      epsp = 0.0
      else
      write(*,*) etype
      end if
      end
