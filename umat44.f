      subroutine umat44 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
c******************************************************************
c|  NMCAD LAB MAT_002 (UMAT)                                      |
c|  ------------------------------------------------------------  |
c|  Copyright @ 2017 NMCAD LAB IISc                               |
c|  All Rights Reserved                                           |
c******************************************************************
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      INTEGER nnpcrv(*)
	  REAL*8 :: capa,epsp,dt1,tt,temper,elsiz
      character*5 etype
      logical failel,reject
      INTEGER8 idele
	  INTEGER:: m
	  REAL*8 :: E1,E2,E3,v12,v21,v23,v32,v31,v13,G12,G23,G31
	  REAL*8 :: S,C11,C12,C13,C23,C22,C33,C44,C55,C66
	  REAL*8 :: G1,G2,G3,G4,G5,G6,S1,S2,S3,S4,S5,S6,detF
	  REAL*8 ::FS11,FS12,FS13,FS21,FS22,FS23,FS31,FS32,FS33
	  common/nani/no_hsvs
c create same common block in urmathn subroutine
c******************************************************************
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
c******************************************************************      
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
	  
c     stiffness matrix 
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
c******************************************************************	  
	  if (etype.eq.'solid'.or.etype.eq.'shl_t'.or.
     1     etype.eq.'sld2d'.or.etype.eq.'tshel'.or.
     2     etype.eq.'sph  '.or.etype.eq.'sldax') then
c  Right Green-St. Venant Strain Tensor
c  G1-G6 will be the local strains in Ea,Eb,Ec directions
        G1 = 0.5*(hsv(m+1)*hsv(m+1)+hsv(m+2)*hsv(m+2)+
     .       hsv(m+3)*hsv(m+3)-1.)
        G2 = 0.5*(hsv(m+4)*hsv(m+4)+hsv(m+5)*hsv(m+5)+
     .       hsv(m+6)*hsv(m+6)-1.)
        G3 = 0.5*(hsv(m+7)*hsv(m+7)+hsv(m+8)*hsv(m+8)+
     .       hsv(m+9)*hsv(m+9)-1.)
        G4 = 0.5*(hsv(m+1)*hsv(m+4)+hsv(m+2)*hsv(m+5)+
     .       hsv(m+3)*hsv(m+6))
        G5 = 0.5*(hsv(m+4)*hsv(m+7)+hsv(m+5)*hsv(m+8)+
     .       hsv(m+6)*hsv(m+9))
        G6 = 0.5*(hsv(m+1)*hsv(m+7)+hsv(m+2)*hsv(m+8)+
     .       hsv(m+3)*hsv(m+9))
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
        detF=hsv(m+1)*(hsv(m+5)*hsv(m+9)-hsv(m+6)*hsv(m+8))-
     .       hsv(m+2)*(hsv(m+4)*hsv(m+9)-hsv(m+6)*hsv(m+7))+
     .       hsv(m+3)*(hsv(m+4)*hsv(m+8)-hsv(m+5)*hsv(m+7))
c
c       Cauchy Stresses (F*S*F^T)
        FS11 = hsv(m+1)*S1 + hsv(m+4)*S4 + hsv(m+7)*S6
        FS12 = hsv(m+1)*S4 + hsv(m+4)*S2 + hsv(m+7)*S5
        FS13 = hsv(m+1)*S6 + hsv(m+4)*S5 + hsv(m+7)*S3
        FS21 = hsv(m+2)*S1 + hsv(m+5)*S4 + hsv(m+8)*S6
        FS22 = hsv(m+2)*S4 + hsv(m+5)*S2 + hsv(m+8)*S5
        FS23 = hsv(m+2)*S6 + hsv(m+5)*S5 + hsv(m+8)*S3
        FS31 = hsv(m+3)*S1 + hsv(m+6)*S4 + hsv(m+9)*S6
        FS32 = hsv(m+3)*S4 + hsv(m+6)*S2 + hsv(m+9)*S5
        FS33 = hsv(m+3)*S6 + hsv(m+6)*S5 + hsv(m+9)*S3

c sig(1) - sig(6) will be global stresses 
      sig(1)=1./detF*(FS11*hsv(m+1)+FS12*hsv(m+4)+FS13*hsv(m+7))
      sig(2)=1./detF*(FS21*hsv(m+2)+FS22*hsv(m+5)+FS23*hsv(m+8))
      sig(3)=1./detF*(FS31*hsv(m+3)+FS32*hsv(m+6)+FS33*hsv(m+9))
      sig(4)=1./detF*(FS11*hsv(m+2)+FS12*hsv(m+5)+FS13*hsv(m+8))
      sig(5)=1./detF*(FS21*hsv(m+3)+FS22*hsv(m+6)+FS23*hsv(m+9))
      sig(6)=1./detF*(FS11*hsv(m+3)+FS12*hsv(m+6)+FS13*hsv(m+9))
		
		hsv(1) = G1
		hsv(2) = G2
		endif
      return
      end