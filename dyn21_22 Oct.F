c
c $Id: dyn21.F 107077 2016-04-13 08:01:42Z tobias $
c
c
#include "define.inc"
#include "define2.inc"
#define MXSCND 5
      subroutine umat43v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsv,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,nnpcrv,cma,qmat,elsizv,idelev,
     . reject)
c******************************************************************
c|  Livermore Software Technology Corporation (LSTC)              |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2007 LSTC                                      |
c|  All Rights Reserved                                           |
c******************************************************************
      include 'nlqparm'
      include 'iounits.inc'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension cm(*),epsps(*),hsv(nlq,*),dt1siz(*)
      dimension temps(*),crv(101,2,*),cma(*),qmat(nlq,3,3),elsizv(*)
      integer nnpcrv(*)
      INTEGER8 idelev(*)
      logical failels(*),reject
      character*5 etype
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
      E2 = cm(12)
      E3 = cm(13)
      v21 = cm(14)
      v31=cm(15)
      v32=cm(16)
      G12 = cm(17)
      G23 = cm(18)
      G31 = cm(19)            
c      
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
c
c     write(59,*) 'Stiffness Matrix'
c     write(59,*) C11, '  ',C12,  '  ',  C13      
c     write(59,*) C22, '  ',C23,  '  ',  C33
c     write(59,*) C44, '  ',C55,  '  ',  C66
c     write(59,*) '      '
c
      do  10 i=lft,llt
      if (etype.eq.'solid'.or.etype.eq.'shl_t'.or.
     1     etype.eq.'sld2d'.or.etype.eq.'tshel'.or.
     2     etype.eq.'sph  '.or.etype.eq.'sldax') then
c
c       integrate 3-D strains (six components)
c       +++++++++++++++++++++++++++++++++++++++++++++++++
c       hsv(i,1)=hsv(i,1)+d1(i)     
c       hsv(i,6)=hsv(i,6)+d2(i)
c       hsv(i,7)=hsv(i,7)+d3(i)
c       hsv(i,8)=hsv(i,8)+d4(i)
c       hsv(i,9)=hsv(i,9)+d5(i)
c       hsv(i,10)=hsv(i,10)+d6(i)
c       +++++++++++++++++++++++++++++++++++++++++++++++++
c
c       Computes stresses 
c      
c       sig1(i) = C11*hsv(i,1)+C12*hsv(i,6)+C13*hsv(i,7)
c       sig2(i) = C12*hsv(i,1)+C22*hsv(i,6)+C23*hsv(i,7)
c       sig3(i) = C13*hsv(i,1)+C23*hsv(i,6)+C33*hsv(i,7)
c       sig4(i) = C44*hsv(i,8)
c       sig5(i) = C55*hsv(i,9)
c       sig6(i) = C66*hsv(i,10)
c
c       sig1(i) = sig1(i)+(C11*d1(i)+C12*d2(i)+C13*d3(i))
c       sig2(i) = sig2(i)+(C12*d1(i)+C22*d2(i)+C23*d3(i))
c       sig3(i) = sig3(i)+(C13*d1(i)+C23*d2(i)+C33*d3(i))
c       sig4(i) = sig4(i)+C44*d4(i)
c       sig5(i) = sig5(i)+C55*d5(i)
c       sig6(i) = sig6(i)+C66*d6(i)
c
c       S1 = C11*d1(i)+C12*d2(i)+C13*d3(i)
c       S2 = C12*d1(i)+C22*d2(i)+C23*d3(i)
c       S3 = C13*d1(i)+C23*d2(i)+C33*d3(i)
c       S4 = C44*d4(i)
c       S5 = C55*d5(i)
c       S6 = C66*d6(i)
c
c       hsv(i,1) = 2.0
c       hsv(i,2) = 0.9
c       hsv(i,3) = 1.2
c       hsv(i,4) = 0.5
c       hsv(i,5) = 3.0
c       hsv(i,6) = 5.0
c       hsv(i,7) = 0.7
c       hsv(i,8) = 3.4
c       hsv(i,9) = 8.0
c
c       Right Green-St. Venant Strain Tensor
        G1 = 0.5*(hsv(i,1)*hsv(i,1)+hsv(i,2)*hsv(i,2)+
     .       hsv(i,3)*hsv(i,3)-1.)
        G2 = 0.5*(hsv(i,4)*hsv(i,4)+hsv(i,5)*hsv(i,5)+
     .       hsv(i,6)*hsv(i,6)-1.)
        G3 = 0.5*(hsv(i,7)*hsv(i,7)+hsv(i,8)*hsv(i,8)+
     .       hsv(i,9)*hsv(i,9)-1.)
        G4 = 0.5*(hsv(i,1)*hsv(i,4)+hsv(i,2)*hsv(i,5)+
     .       hsv(i,3)*hsv(i,6))
        G5 = 0.5*(hsv(i,4)*hsv(i,7)+hsv(i,5)*hsv(i,8)+
     .       hsv(i,6)*hsv(i,9))
        G6 = 0.5*(hsv(i,1)*hsv(i,7)+hsv(i,2)*hsv(i,8)+
     .       hsv(i,3)*hsv(i,9))
c
c       Left Green-St. Venant Strain Tensor
c       G1 = 0.5*(hsv(i,1)*hsv(i,1)+hsv(i,4)*hsv(i,4)+
c     .      hsv(i,7)*hsv(i,7)-1.)
c       G2 = 0.5*(hsv(i,2)*hsv(i,2)+hsv(i,5)*hsv(i,5)+
c     .      hsv(i,8)*hsv(i,8)-1.)
c       G3 = 0.5*(hsv(i,3)*hsv(i,3)+hsv(i,6)*hsv(i,6)+
c     .      hsv(i,9)*hsv(i,9)-1.)
c       G4 = 0.5*(hsv(i,1)*hsv(i,2)+hsv(i,4)*hsv(i,5)+
c     .      hsv(i,7)*hsv(i,8))
c       G5 = 0.5*(hsv(i,2)*hsv(i,3)+hsv(i,5)*hsv(i,6)+
c     .      hsv(i,8)*hsv(i,9))
c       G6 = 0.5*(hsv(i,1)*hsv(i,3)+hsv(i,4)*hsv(i,6)+
c     .      hsv(i,7)*hsv(i,9))                     
c
c       write(59,*) 'Green-St. Venant Strain Tensor'
c       write(59,*) G1, '  ' ,G4,  '  ' ,G6      
c       write(59,*) G4, '  ' ,G2,  '  ' ,G5
c       write(59,*) G6, '  ' ,G5,  '  ' ,G3
c       write(59,*) '      '
c
c       Second Piola-Kirchhoff Stress Tensor
        S1 = C11*G1 + C12*G2 + C13*G3
        S2 = C12*G1 + C22*G2 + C23*G3
        S3 = C13*G1 + C23*G2 + C33*G3
        S4 = 2.*C44*G4
        S5 = 2.*C55*G5
        S6 = 2.*C66*G6
c
c       write(59,*) 'Second Piola-Kirchhoff Stress Tensor'
c       write(59,*) S1, '  ' ,S4,  '  ' ,S6      
c       write(59,*) S4, '  ' ,S2,  '  ' ,S5
c       write(59,*) S6, '  ' ,S5,  '  ' ,S3
c       write(59,*) '      '
c
c       Jacobian
        detF=hsv(i,1)*(hsv(i,5)*hsv(i,9)-hsv(i,6)*hsv(i,8))-
     .       hsv(i,2)*(hsv(i,4)*hsv(i,9)-hsv(i,6)*hsv(i,7))+
     .       hsv(i,3)*(hsv(i,4)*hsv(i,8)-hsv(i,5)*hsv(i,7))
c
c       write(59,*) 'Jacobian'
c       write(59,*) detF
c       write(59,*) '      '
c
c       Cauchy Stresses (F*S*F^T)
        FS11 = hsv(i,1)*S1 + hsv(i,4)*S4 + hsv(i,7)*S6
        FS12 = hsv(i,1)*S4 + hsv(i,4)*S2 + hsv(i,7)*S5
        FS13 = hsv(i,1)*S6 + hsv(i,4)*S5 + hsv(i,7)*S3
        FS21 = hsv(i,2)*S1 + hsv(i,5)*S4 + hsv(i,8)*S6
        FS22 = hsv(i,2)*S4 + hsv(i,5)*S2 + hsv(i,8)*S5
        FS23 = hsv(i,2)*S6 + hsv(i,5)*S5 + hsv(i,8)*S3
        FS31 = hsv(i,3)*S1 + hsv(i,6)*S4 + hsv(i,9)*S6
        FS32 = hsv(i,3)*S4 + hsv(i,6)*S2 + hsv(i,9)*S5
        FS33 = hsv(i,3)*S6 + hsv(i,6)*S5 + hsv(i,9)*S3
c
c       write(59,*) 'F*S'
c       write(59,*) FS11, '  ' ,FS12,  '  ' ,FS13      
c       write(59,*) FS21, '  ' ,FS22,  '  ' ,FS23
c       write(59,*) FS31, '  ' ,FS32,  '  ' ,FS33
c       write(59,*) '      '
c
        sig1(i) = 1./detF*(FS11*hsv(i,1)+FS12*hsv(i,4)+FS13*hsv(i,7))
        sig2(i) = 1./detF*(FS21*hsv(i,2)+FS22*hsv(i,5)+FS23*hsv(i,8))
        sig3(i) = 1./detF*(FS31*hsv(i,3)+FS32*hsv(i,6)+FS33*hsv(i,9))
        sig4(i) = 1./detF*(FS11*hsv(i,2)+FS12*hsv(i,5)+FS13*hsv(i,8))
        sig5(i) = 1./detF*(FS21*hsv(i,3)+FS22*hsv(i,6)+FS23*hsv(i,9))
        sig6(i) = 1./detF*(FS11*hsv(i,3)+FS12*hsv(i,6)+FS13*hsv(i,9))
c
c       write(59,*) 'Cauchy Stresses '
c       write(59,*) sig1(i), '  ' ,sig4(i),  '  ' ,sig6(i)    
c       write(59,*) sig4(i), '  ' ,sig2(i),  '  ' ,sig5(i)
c       write(59,*) sig6(i), '  ' ,sig5(i),  '  ' ,sig3(i)
c       write(59,*) '      '
c
c       write(59,*) '-----------------------------'
c
      endif
   10 continue
c       
      return
      end
      subroutine umat42v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,nnpcrv,cma,qmat,elsizv,idelev,
     . reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation (LSTC)              |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2007 LSTC                                      |
c|  All Rights Reserved                                           |
c******************************************************************
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension cm(*),epsps(*),hsvs(nlq,*),dt1siz(*)
      dimension temps(*),crv(lq1,2,*),cma(*),qmat(nlq,3,3),elsizv(*)
      integer nnpcrv(*)
      INTEGER8 idelev(*)
      logical failels(*),reject
      character*5 etype
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
      E2 = cm(12)
      E3 = cm(13)
      v21 = cm(14)
      v31=cm(15)
      v32=cm(16)
      G12 = cm(17)
      G23 = cm(18)
      G31 = cm(19)            
c      
      v12 = v21*E1/E2
      v23 = v32*E2/E3
      v13 = v31*E1/E3
c
c     stiffness matrix 
c
      S = 1.-v12*v21-v23*v32-v31*v13-2.*v21*v32*v13
c
      C11=(1. - v23*v32)*E1/S
      C12=(v21 + v31*v23)*E1/S
      C13=(v31 + v21*v32)*E1/S
c     C13=(v31 + v21*v23)/(E2*E3*S)
      C22=(1. - v31*v13)*E2/S
      C23=(v32 + v31*v12)*E2/S
      C33=(1. - v12*v21)*E3/S
      C44=G12
      C55=G23
      C66=G31
c
c     write(59,*) 'Stiffness Matrix'
c     write(59,*) C11, '  ',C12,  '  ',  C13      
c     write(59,*) C22, '  ',C23,  '  ',  C33
c     write(59,*) C44, '  ',C55,  '  ',  C66
c     write(59,*) '      '
c
      do  10 i=lft,llt
      if (etype.eq.'solid'.or.etype.eq.'shl_t'.or.
     1     etype.eq.'sld2d'.or.etype.eq.'tshel'.or.
     2     etype.eq.'sph  '.or.etype.eq.'sldax') then
c
c       Computes stresses 
c
        sig1(i) = sig1(i)+(C11*d1(i)+C12*d2(i)+C13*d3(i))
        sig2(i) = sig2(i)+(C12*d1(i)+C22*d2(i)+C23*d3(i))
        sig3(i) = sig3(i)+(C13*d1(i)+C23*d2(i)+C33*d3(i))
        sig4(i) = sig4(i)+C44*d4(i)
        sig5(i) = sig5(i)+C55*d5(i)
        sig6(i) = sig6(i)+C66*d6(i) 
c
c       sig1(i) = sig1(i)+(C11*d1(i)+C12*d2(i)+C13*d3(i))*dt1siz(i)
c       sig2(i) = sig2(i)+(C12*d1(i)+C22*d2(i)+C23*d3(i))*dt1siz(i)
c       sig3(i) = sig3(i)+(C13*d1(i)+C23*d2(i)+C33*d3(i))*dt1siz(i)
c       sig4(i) = sig4(i)+C44*d4(i)*dt1siz(i)
c       sig5(i) = sig5(i)+C55*d5(i)*dt1siz(i)
c       sig6(i) = sig6(i)+C66*d6(i)*dt1siz(i)
c
      else if (etype.eq.'shell') then
c       hsv(i,7)=-(hsv(i,4)*C13*hsv(i,1) + C23*hsv(i,6))/C33
        d3(i) = -(C13*d1(i) + C23*d2(i))/C33
        sig1(i) = sig1(i)+(C11*d1(i)+C12*d2(i)+C13*d3(i))
        sig2(i) = sig2(i)+(C12*d1(i)+C22*d2(i)+C23*d3(i))               
        sig3(i) = 0.0
        sig4(i) = sig4(i)+C44*d4(i)
        sig5(i) = sig5(i)+capa*C55*d5(i)
        sig6(i) = sig6(i)+capa*C66*d6(i)
         endif
   10 continue
c
      return
      end
      subroutine urmats(lft,llt,cm,capa,mt,crv,nnpcrv,
     . ipt,rcoor,scoor,tcoor,
     . nnm1,nconstp,nip,ipt_thk,npc,plc,eltype,nshbwp)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
#include "bk05.inc"
      include 'bk06.inc'
#include "bk07.inc"
#include "bk13.inc"
      include 'bk19.inc'
      include 'iounits.inc'
      include 'memaia.inc'
      include 'umatss.inc'
      include 'shlioc.inc'
#include "dynmem.inc"
#include "implicit1.inc"
C_TASKCOMMON (aux2loc)
C_TASKCOMMON (aux8cloc)
C_TASKCOMMON (aux14loc)
C_TASKCOMMON (aux19loc)
C_TASKCOMMON (aux33loc)
C_TASKCOMMON (bk36loc)
C_TASKCOMMON (bux13loc)
C_TASKCOMMON (failuloc)
C_TASKCOMMON (failcmloc)
C_TASKCOMMON (frozloc)
C_TASKCOMMON (hourgloc)
C_TASKCOMMON (shloptloc)
C_TASKCOMMON (sidesloc)
C_TASKCOMMON (soundloc)
C_TASKCOMMON (subtssloc)
      common/aux2loc/e1(nlq),e2(nlq),d3(nlq),e4(nlq),e5(nlq),e6(nlq),
     1 wzzdt(nlq),wyydt(nlq),wxxdt(nlq),einc(nlq)
      common /aux8cloc/
     & d1(nlq),d2(nlq),d4(nlq),d5(nlq),d6(nlq),
     & dfg11(nlq),dfg21(nlq),dfg31(nlq),dfg12(nlq),dfg22(nlq),
     & dfg32(nlq),dfg13(nlq),dfg23(nlq),dfg33(nlq),stg5(nlq),
     & stg6(nlq),a11(nlq),a12(nlq),a21(nlq),a22(nlq)
      common/aux14loc/
     & sig1(nlq),sig2(nlq),sig3(nlq),sig4(nlq),
     & sig5(nlq),sig6(nlq),epsps(nlq),hsvs(nlq,NHISVAR)
      common/aux19loc/
     1 sign0(nlq),sign1(nlq),sign2(nlq),sign3(nlq),sign4(nlq),
     2 sign5(nlq),sign6(nlq)
      common/aux33loc/
     1 ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq),ixs(nlq,4),mxt(nlq)
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     1 ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     2 grvity,idirgv,nodspc,nspcor,numelh10,numels8,numelsx,numnpin,
     3 numelhx,nmuxct,ndscmxh,ndscmxs,numnurb2d,nvcbox,lscmh,numnurb3d
      common/bk36loc/index
      common/bux13loc/eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      common/failcmloc/ifail(nlq)
      common/failuloc/sieu(nlq),fail(nlq),ifaili(nlq),nfipts(nlq)
      common/frozloc/yhat(nlq,3,4),qloc(nlq,3,3),qmat(nlq,3,3)
      common/hourgloc/ym(nlq),gm(nlq),ifsv(nlq)
      common/shloptloc/ibelyt
c
      common/bk01/itherm,itemp,ntmp0,ntmp1,itempan,itempdr,itmpe(44)
      common/bk26/begtim,nintcy
      common/bk28/summss,xke,xpe,tt,xte0,erodeke,erodeie,selie,selke,
     . erodehg,eintcnt,engjnt
      logical failur,failgl
      common/failcm/failur,failgl
      common/numcpu/ncpu,ncpua,ncpub,lenvec(8)
      integer istrn,istupd,ibelyts,miter,ipstpd,intsts,nodsts,intstn,
     1 nodstn,jstrn,nddet,nsout01,nsout02,nsout03,nsout04,
     2 nsout05,nsout06
      real wrpang
      common/shlopt/istrn,istupd,ibelytt,miter,wrpang,
     1 ipstpd,intsts,nodsts,intstn,nodstn,jstrn,nddet,nsout01,
     2 nsout02,nsout03,nsout04,nsout05,nsout06,mporlk,lporlk,
     3 istrn_t,istrn_p
      common/sidesloc/sidmn(nlq)
      common/soundloc/ss(nlq),sndsp(nlq),diagm(nlq),sarea(nlq),
     . dxl(nlq)
      common/subtssloc/dt1siz(nlq)
      common/thropt/itopaz(101)
c
      dimension cm(*),crv(lq1,2,*),nconstp(*),npc(*),plc(*)
      integer nnpcrv(*)
c     dimension eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      dimension qmats(3,3)
      character*(*)eltype
      dimension nshbwp(*)
      dimension thhsv(nlq,100),thhsvi(100)
      dimension elsizv(nlq)
      logical failel,failels(nlq),reject
      INTEGER8 lqfinv8
      INTEGER8 idelev(nlq),idele
      MINDEX ptr_usrcf
c
c     data temps/nlq*0.0/
c
c     processing elements
c
      mx=48*(mxt(lft)-1)
      gm(lft)=cm(mx+abs(ishrmp(mt)))
      bk=cm(mx+ibulkp(mt))
      pr=(3.*bk-2.*gm(lft))/(2.*(3.*bk+gm(lft)))
      ym(lft)=2.*(1.+pr)*gm(lft)
      ss(lft)=ym(lft)/(1.-pr*pr)
c
c     Determining number of requested history variables
c
      no_hsvs=nconstp(mxt(lft))
      if (ivumat(mt).ne.0) no_hsvs=no_hsvs-9-6
      if (lenvec(5) .eq.1) no_hsvs=no_hsvs-6
      if (iorien(mt).ne.0) no_hsvs=no_hsvs-2
      if (ihyper(mt).ne.0) no_hsvs=no_hsvs-9
      if (iabs(ihyper(mt)).eq.3) no_hsvs=no_hsvs-4
      if (iabs(ihyper(mt)).eq.4) no_hsvs=no_hsvs-18
      if (ifipts(mt).lt.0) no_hsvs=no_hsvs-1
      if (ioshl(67) .eq.1) no_hsvs=no_hsvs-1
c
      if (ioshl(881).ne.0) then
c contact forces are store under r_mem(3,numnp)
        ptr_usrcf = memh_ptr(ioshl(881))
      endif
c
c     initializing indeces
c
      if(ivumat(mt).ne.0) then

        if (iorien(mt).ne.0) then
          if (ihyper(mt).ne.0) then
            ind_defgrad=1+no_hsvs
            ind_ortho=ind_defgrad+9
          else
            ind_ortho=1+no_hsvs
          endif
          if (ifipts(mt).lt.0) then
            ind_fipts=ind_ortho
            ind_ortho=ind_ortho+1
          endif
          ind_q1=ind_ortho
          ind_q2=ind_ortho+1
          ind_last=ind_q2
        else
          if (ihyper(mt).ne.0) then
            ind_defgrad=1+no_hsvs
            ind_last=ind_defgrad+8
          else
            ind_last=no_hsvs
          endif
          if (ifipts(mt).lt.0) then
            ind_fipts=ind_last+1
            ind_last=ind_last+1
          endif             
        endif
        ind_vumat=ind_last+1
        ind_last=ind_vumat+9+6-1
c        ind_last=ind_last+1

      else

        if (iorien(mt).ne.0) then
          if (ihyper(mt).ne.0) then
            ind_defgrad=1+no_hsvs
            ind_ortho=ind_defgrad+9
            if (iabs(ihyper(mt)).eq.3) then
               ind_thick=ind_ortho
               ind_ortho=ind_ortho+4
            endif
            if (iabs(ihyper(mt)).eq.4) ind_ortho=ind_ortho+18
          else
            ind_ortho=1+no_hsvs
          endif
          if (ifipts(mt).lt.0) then
            ind_fipts=ind_ortho
            ind_ortho=ind_ortho+1
          endif
          ind_q1=ind_ortho
          ind_q2=ind_ortho+1
          ind_last=ind_q2
        else
          if (ihyper(mt).ne.0) then
            ind_defgrad=1+no_hsvs
            ind_last=ind_defgrad+8
            if (iabs(ihyper(mt)).eq.3) then
               ind_thick=ind_last+1
               ind_last=ind_last+4
            endif
            if (iabs(ihyper(mt)).eq.4) ind_last=ind_last+18
          else
            ind_last=no_hsvs
          endif
          if (ifipts(mt).lt.0) then
            ind_fipts=ind_last+1
            ind_last=ind_last+1
          endif        
        endif 

      endif
c
      if (ind_last.gt.NHISVAR) then
c       write(iotty,900) lqfmiv(mxt(lft))
c       write(iohsp,900) lqfmiv(mxt(lft))
c       write(iomsg,900) lqfmiv(mxt(lft))
c       call adios(TC_ERROR)
        ierdat(1)=lqfmiv(mxt(lft))
        call lsmsg(3,MSG_SOL+1148,ioall,ierdat,rerdat,cerdat,0)
      endif
c$$$     if (tt.eq.0.0.and.nintcy.eq.0) then
c$$$c
c$$$c      move transformation matrix in history variables array
c$$$c
c$$$       if (iorien(mt).ne.0) then
c$$$         do i=lft,llt
c$$$           hsvs(i,ind_q2)=hsvs(i,2)
c$$$           hsvs(i,ind_q1)=hsvs(i,1)
c$$$         enddo
c$$$         if (no_hsvs.gt.0) then
c$$$           do j=1,no_hsvs
c$$$             do i=lft,llt
c$$$               hsvs(i,j)=0.0
c$$$             enddo
c$$$           enddo
c$$$         endif
c$$$       endif
c$$$     endif
c
      if (ipt.eq.1.and.ifipts(mt).ne.0) then
        do i=lft,llt
          nfipts(i)=0
        enddo
      endif
c
      num_nods=4
      if (ibelyt.eq.0) then
        num_nods=8
        do 10 i=lft,llt
        sign0(i)=sign3(i)+ym(lft)*d3(i)
        sign3(i)=sign0(i)
 10     continue
      endif
c
c     transform stresses and strains to material system
c
      if (iorien(mt).ne.0.and.(ihyper(mt).ge.0.or.ihyper(mt).eq.-2))
     $     then
c     if (iorien(mt).ne.0.and.ihyper(mt).ge.0) then
        do 20 i=lft,llt
        stg5(i)=sig5(i)
        stg6(i)=sig6(i)
        if (ncycle.le.1) then
          rr = 1./sqrt(hsvs(i,ind_q1)**2+hsvs(i,ind_q2)**2)
          hsvs(i,ind_q1)=hsvs(i,ind_q1)*rr
          hsvs(i,ind_q2)=hsvs(i,ind_q2)*rr
        end if
        a11(i) =hsvs(i,ind_q1)*sig1(i)-hsvs(i,ind_q2)*sig4(i)
        a12(i) =hsvs(i,ind_q2)*sig1(i)+hsvs(i,ind_q1)*sig4(i)
        a21(i) =hsvs(i,ind_q1)*sig4(i)-hsvs(i,ind_q2)*sig2(i)
        a22(i) =hsvs(i,ind_q2)*sig4(i)+hsvs(i,ind_q1)*sig2(i)
        sig1(i)=hsvs(i,ind_q1)*a11(i)-hsvs(i,ind_q2)*a21(i)
        sig2(i)=hsvs(i,ind_q2)*a12(i)+hsvs(i,ind_q1)*a22(i)
        sig4(i)=hsvs(i,ind_q1)*a12(i)-hsvs(i,ind_q2)*a22(i)
        sig5(i)=hsvs(i,ind_q2)*stg6(i)+hsvs(i,ind_q1)*stg5(i)
        sig6(i)=hsvs(i,ind_q1)*stg6(i)-hsvs(i,ind_q2)*stg5(i)
 20     continue
c
        do 30 i=lft,llt
        d4(i)  =.5*e4(i)
        a11(i) =hsvs(i,ind_q1)*e1(i)-hsvs(i,ind_q2)*d4(i)
        a12(i) =hsvs(i,ind_q2)*e1(i)+hsvs(i,ind_q1)*d4(i)
        a21(i) =hsvs(i,ind_q1)*d4(i)-hsvs(i,ind_q2)*e2(i)
        a22(i) =hsvs(i,ind_q2)*d4(i)+hsvs(i,ind_q1)*e2(i)
        d1(i)  =hsvs(i,ind_q1)*a11(i)-hsvs(i,ind_q2)*a21(i)
        d2(i)  =hsvs(i,ind_q2)*a12(i)+hsvs(i,ind_q1)*a22(i)
        d4(i)  =2.*(hsvs(i,ind_q1)*a12(i)-hsvs(i,ind_q2)*a22(i))
        d5(i)  =hsvs(i,ind_q2)*e6(i)+hsvs(i,ind_q1)*e5(i)
        d6(i)  =hsvs(i,ind_q1)*e6(i)-hsvs(i,ind_q2)*e5(i)
 30     continue
c
c     element local strains
c
      else
        do i=lft,llt
          d1(i)  =e1(i)
          d2(i)  =e2(i)
          d4(i)  =e4(i)
          d5(i)  =e5(i)
          d6(i)  =e6(i)
        enddo
      endif
c
c     storage of deformation gradient
c
      if (ivumat(mt).ne.0) then
        ifo11=ind_vumat
        ifo21=ind_vumat+1
        ifo31=ind_vumat+2
        ifo12=ind_vumat+3
        ifo22=ind_vumat+4
        ifo32=ind_vumat+5
        ifo13=ind_vumat+6
        ifo23=ind_vumat+7
        ifo33=ind_vumat+8
      endif
      if (ihyper(mt).ne.0) then
        if11=ind_defgrad
        if21=ind_defgrad+1
        if31=ind_defgrad+2
        if12=ind_defgrad+3
        if22=ind_defgrad+4
        if32=ind_defgrad+5
        if13=ind_defgrad+6
        if23=ind_defgrad+7
        if33=ind_defgrad+8
c
c       deformation gradient in element system at time = 0
c
        if (tt.eq.begtim.and.nintcy.eq.0) then
          do i=lft,llt
            hsvs(i,if11)=1.0
            hsvs(i,if21)=0.0
            hsvs(i,if31)=0.0
            hsvs(i,if12)=0.0
            hsvs(i,if22)=1.0
            hsvs(i,if32)=0.0
            hsvs(i,if13)=0.0
            hsvs(i,if23)=0.0
            hsvs(i,if33)=1.0
          enddo
          if (ivumat(mt).ne.0) then
            do i=lft,llt
              hsvs(i,ifo11)=1.0
              hsvs(i,ifo21)=0.0
              hsvs(i,ifo31)=0.0
              hsvs(i,ifo12)=0.0
              hsvs(i,ifo22)=1.0
              hsvs(i,ifo32)=0.0
              hsvs(i,ifo13)=0.0
              hsvs(i,ifo23)=0.0
              hsvs(i,ifo33)=1.0
            enddo
          endif
        endif
c
c       store deformation gradient in element system
c
        if (iorien(mt).eq.0) then
          do i=lft,llt
            dfg11(i)=hsvs(i,if11)
            dfg21(i)=hsvs(i,if21)
            dfg31(i)=hsvs(i,if31)
            dfg12(i)=hsvs(i,if12)
            dfg22(i)=hsvs(i,if22)
            dfg32(i)=hsvs(i,if32)
            dfg13(i)=hsvs(i,if13)
            dfg23(i)=hsvs(i,if23)
            dfg33(i)=hsvs(i,if33)
           enddo
        endif
      endif
c
      if (iorien(mt).ne.0) then
c
c       store deformation gradient in material system
c       as F_bar=Q^T*F
c
      if (ihyper(mt).gt.0) then
c        if (ihyper(mt).eq.1) then
          do i=lft,llt
            dfg11(i)=hsvs(i,ind_q1)*hsvs(i,if11)-
     1               hsvs(i,ind_q2)*hsvs(i,if21)
            dfg21(i)=hsvs(i,ind_q2)*hsvs(i,if11)+
     1               hsvs(i,ind_q1)*hsvs(i,if21)
            dfg31(i)=hsvs(i,if31)
            dfg12(i)=hsvs(i,ind_q1)*hsvs(i,if12)-
     1               hsvs(i,ind_q2)*hsvs(i,if22)
            dfg22(i)=hsvs(i,ind_q2)*hsvs(i,if12)+
     1               hsvs(i,ind_q1)*hsvs(i,if22)
            dfg32(i)=hsvs(i,if32)
            dfg13(i)=hsvs(i,ind_q1)*hsvs(i,if13)-
     1               hsvs(i,ind_q2)*hsvs(i,if23)
            dfg23(i)=hsvs(i,ind_q2)*hsvs(i,if13)+
     1               hsvs(i,ind_q1)*hsvs(i,if23)
            dfg33(i)=hsvs(i,if33)
c
            hsvs(i,if11)=dfg11(i)
            hsvs(i,if21)=dfg21(i)
            hsvs(i,if12)=dfg12(i)
            hsvs(i,if22)=dfg22(i)
            hsvs(i,if13)=dfg13(i)
            hsvs(i,if23)=dfg23(i)
          enddo
c
c       store deformation gradient in material system
c       as F_bar=Q^T*F*Q
c
        else if (ihyper(mt).eq.-2) then
          do i=lft,llt
            dfg11(i)=hsvs(i,ind_q1)*hsvs(i,if11)-
     1               hsvs(i,ind_q2)*hsvs(i,if21)
            dfg21(i)=hsvs(i,ind_q2)*hsvs(i,if11)+
     1               hsvs(i,ind_q1)*hsvs(i,if21)
            dfg31(i)=hsvs(i,if31)
            dfg12(i)=hsvs(i,ind_q1)*hsvs(i,if12)-
     1               hsvs(i,ind_q2)*hsvs(i,if22)
            dfg22(i)=hsvs(i,ind_q2)*hsvs(i,if12)+
     1               hsvs(i,ind_q1)*hsvs(i,if22)
            dfg32(i)=hsvs(i,if32)
            dfg13(i)=hsvs(i,ind_q1)*hsvs(i,if13)-
     1               hsvs(i,ind_q2)*hsvs(i,if23)
            dfg23(i)=hsvs(i,ind_q2)*hsvs(i,if13)+
     1               hsvs(i,ind_q1)*hsvs(i,if23)
            dfg33(i)=hsvs(i,if33)
c
            hsvs(i,if11)=dfg11(i)*hsvs(i,ind_q1)-dfg12(i)*hsvs(i,ind_q2)
            hsvs(i,if12)=dfg11(i)*hsvs(i,ind_q2)+dfg12(i)*hsvs(i,ind_q1)
            hsvs(i,if13)=dfg13(i)
            hsvs(i,if21)=dfg21(i)*hsvs(i,ind_q1)-dfg22(i)*hsvs(i,ind_q2)
            hsvs(i,if22)=dfg21(i)*hsvs(i,ind_q2)+dfg22(i)*hsvs(i,ind_q1)
            hsvs(i,if23)=dfg23(i)
            hsvs(i,if31)=dfg31(i)*hsvs(i,ind_q1)-dfg32(i)*hsvs(i,ind_q2)
            hsvs(i,if32)=dfg31(i)*hsvs(i,ind_q2)+dfg32(i)*hsvs(i,ind_q1)
            hsvs(i,if33)=dfg33(i)
            
            dfg11(i)=hsvs(i,if11)
            dfg22(i)=hsvs(i,if22)
            dfg33(i)=hsvs(i,if33)
            dfg12(i)=hsvs(i,if12)
            dfg21(i)=hsvs(i,if21)
            dfg13(i)=hsvs(i,if13)
            dfg31(i)=hsvs(i,if31)
            dfg23(i)=hsvs(i,if23)
            dfg32(i)=hsvs(i,if32)

          enddo

c
c       store deformation gradient in material system
c       as F_bar=F*Q
c

c        else if (ihyper(mt).lt.0) then
        else if (ihyper(mt).eq.-1) then
          do i=lft,llt
            dfg11(i)=hsvs(i,ind_q1)*hsvs(i,if11)-
     1               hsvs(i,ind_q2)*hsvs(i,if12)
            dfg12(i)=hsvs(i,ind_q2)*hsvs(i,if11)+
     1               hsvs(i,ind_q1)*hsvs(i,if12)
            dfg13(i)=hsvs(i,if13)
            dfg21(i)=hsvs(i,ind_q1)*hsvs(i,if21)-
     1               hsvs(i,ind_q2)*hsvs(i,if22)
            dfg22(i)=hsvs(i,ind_q2)*hsvs(i,if21)+
     1               hsvs(i,ind_q1)*hsvs(i,if22)
            dfg23(i)=hsvs(i,if23)
            dfg31(i)=hsvs(i,ind_q1)*hsvs(i,if31)-
     1               hsvs(i,ind_q2)*hsvs(i,if32)
            dfg32(i)=hsvs(i,ind_q2)*hsvs(i,if31)+
     1               hsvs(i,ind_q1)*hsvs(i,if32)
            dfg33(i)=hsvs(i,if33)
c
            hsvs(i,if11)=dfg11(i)
            hsvs(i,if12)=dfg12(i)
            hsvs(i,if21)=dfg21(i)
            hsvs(i,if22)=dfg22(i)
            hsvs(i,if31)=dfg31(i)
            hsvs(i,if32)=dfg32(i)
          enddo
        endif
      endif
c
c     update first two rows of deformation gradient at time = n+1
c
      if (ihyper(mt).ne.0) then
        call integf_s(hsvs(1,if11),hsvs(1,if21),hsvs(1,if31),
     &   hsvs(1,if12),hsvs(1,if22),hsvs(1,if32),
     &   hsvs(1,if13),hsvs(1,if23),hsvs(1,if33),lft,llt)
      endif
c
      if (iabs(ihyper(mt)).eq.3.or.iabs(ihyper(mt)).eq.4) then
         idxshl=5
         if (numels8.gt.0) idxshl=9
         if (iabs(ihyper(mt)).eq.3) then
            do i=lft,llt
               d3(i)=0.
            enddo
            if (hsvs(lft,ind_thick).eq.0.) then               
               do i=lft,llt
                  hsvs(i,ind_thick  )=a(ns05+9*(nnm1+i-1)  )
                  hsvs(i,ind_thick+1)=a(ns05+9*(nnm1+i-1)+1)
                  hsvs(i,ind_thick+2)=a(ns05+9*(nnm1+i-1)+2)
                  hsvs(i,ind_thick+3)=a(ns05+9*(nnm1+i-1)+3)
               enddo
            endif
         else
            ind_thick=ind_last+1
            do i=lft,llt
               hsvs(i,ind_last+1)=a(ns05+9*(nnm1+i-1)  )
               hsvs(i,ind_last+2)=a(ns05+9*(nnm1+i-1)+1)
               hsvs(i,ind_last+3)=a(ns05+9*(nnm1+i-1)+2)
               hsvs(i,ind_last+4)=a(ns05+9*(nnm1+i-1)+3)
            enddo
         endif
         call usrshl_defgrad(hsvs(1,if11),ihyper(mt),
     1        rcoor,scoor,tcoor,
     2        r_mem(dm_x),r_mem(dm_rots),
     3        a(ns03+3*(idxshl-1)*numels+12*nnm1),
     4        a(ns03+3*(idxshl+3)*numels+ 4*nnm1),
     5        hsvs(1,ind_thick),
     6        lft,-llt)         
         do k=1,3
            do j=1,3               
               do i=lft,llt
                  qmat(i,j,k)=qloc(i,j,k)
               enddo
            enddo
         enddo
         if (iorien(mt).ne.0) then
            do j=1,3
               do i=lft,llt
                  q1=hsvs(i,ind_q1)*qmat(i,j,1)-
     1                 hsvs(i,ind_q2)*qmat(i,j,2)
                  q2=hsvs(i,ind_q2)*qmat(i,j,1)+
     1                 hsvs(i,ind_q1)*qmat(i,j,2)
                  qmat(i,j,1)=q1
                  qmat(i,j,2)=q2
               enddo
            enddo
         endif
      endif
c
c     if ishrmp(mt).lt.0 then the flag for temperature is active
c     temperatures are stored in the temps array.
c
      if (ishrmp(mt).lt.0) then
        call usr_temps (a(ntmp0+1),a(n19),itemp,temps,lft,llt,ix1,ix2,
     .   ix3,ix4,ixs,ixs(1,2),ixs(1,3),ixs(1,4),num_nods,itopaz,
     .   nnm1,rcoor,scoor,tcoor)
      else
        do i=lft,llt
          temps(i)=0.0
        enddo
      endif
c
c     Get thermal history variables from thermal
c     user material
c
      call get_thhsv(nthhsv,thhsv,itemp,lft,llt,
     1 num_nods,itopaz,nnm1,rcoor,scoor,tcoor)
c
c     energy calculations
c
      do 33 i=lft,llt
      einc(i)=d1(i)*sig1(i)+d2(i)*sig2(i)+d4(i)*sig4(i)+d5(i)*sig5(i)
     1       +d6(i)*sig6(i)
 33   continue
c
c     characteristic element size and element id
      do i=lft,llt
        elsizv(i)=sqrt(sidmn(i))
        ielm=nshbwp(nnm1+i)
        if (ielm.gt.0) idelev(i)=lqfinv8(ielm,4)
      enddo
c
c     memory pointer for extra memory data
      lcma=ncma+nmmat+ia(ncma+mxt(lft)-1)-1
c
c     vectorization flag
c
      ivect=ivectr(mt)
c
c     implicit rejection flag
c
      reject=.false.
c
      mte40=mt-40
c
c     scalar umat
c
      if (ivect.eq.0) then
        do 120 i=lft,llt
c
c       gathering of variables
c
        index=i
        dt1   =dt1siz(i)
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
        epsp=epsps(i)
        temper=temps(i)
        failel=.false.
        if (ifipts(mt).lt.0) failel=hsvs(i,ind_fipts).ne.0.
        if (no_hsvs.gt.0) then
          do j=1,no_hsvs
            hsv(j)=hsvs(i,j)
          enddo
        endif
        if (ihyper(mt).ne.0) then
          do j=1,9
            hsv(no_hsvs+j)=hsvs(i,ind_defgrad-1+j)
          enddo
          if (iabs(ihyper(mt)).eq.3) then
             do j=1,4
                hsv(no_hsvs+j+9)=hsvs(i,ind_thick-1+j)
             enddo
          endif
          if (iabs(ihyper(mt)).eq.4) then
          do j=1,18
            hsv(no_hsvs+j+9)=hsvs(i,ind_defgrad-1+j+9)
          enddo
          endif
          if (iabs(ihyper(mt)).eq.3.or.iabs(ihyper(mt)).eq.4) then
             do k=1,3
                do j=1,3
                   qmats(j,k)=qmat(i,j,k)
                enddo
             enddo
          endif
        endif
        if (nthhsv.gt.0) then
          do j=1,nthhsv
            thhsvi(j)=thhsv(i,j)
          enddo
        endif        
        elsiz=elsizv(i)
        idele=idelev(i)
c
c       call user developed subroutines here
c
        go to (41,42,43,44,45,46,47,48,49,50), mte40
   41   call umat41 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,nnpcrv,a(lcma),qmats,elsiz,idele,reject)
        go to 60
   42   call umat42 (cm(mx+1),eps,sig,hsv,dt1,capa,'shell',tt,crv,
     .   nnpcrv,
     .   a(lcma),temper,qmats,elsiz,idele,reject)
        go to 60
   43   call umat43 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,nnpcrv,a(lcma),qmats,elsiz,idele,reject)
        go to 60
   44   call umat44 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,nnpcrv,a(lcma),qmats,elsiz,idele,reject)
        go to 60
   45   call umat45 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,nnpcrv,a(lcma),qmats,elsiz,idele,ihyper(mt),
     .   reject)
        go to 60
   46   call umat46 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,nnpcrv,
     .   a(lcma),thhsvi,nthhsv,qmats,elsiz,idele,reject)
        go to 60
   47   call umat47 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,nnpcrv,a(lcma),qmats,elsiz,idele,reject)
        go to 60
   48   call umat48 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,nnpcrv,a(lcma),qmats,elsiz,idele,reject)
        go to 60
   49   call umat49 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,nnpcrv,a(lcma),qmats,elsiz,idele,reject)
        go to 60
   50   call umat50 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,failel,crv,nnpcrv,a(lcma),qmats,elsiz,idele,reject)
c
c       scattering of variables
c
 60     if (no_hsvs.gt.0) then
          do j=1,no_hsvs
            hsvs(i,j)=hsv(j)
          enddo
        endif
c
        if (ihyper(mt).ne.0) then
          do j=1,9
            hsvs(i,ind_defgrad-1+j)=hsv(no_hsvs+j)
          enddo
          if (iabs(ihyper(mt)).eq.3) then
          do j=1,4
            hsvs(i,ind_thick-1+j)=hsv(no_hsvs+j+9)
          enddo
          endif
          if (iabs(ihyper(mt)).eq.4) then
          do j=1,18
            hsvs(i,ind_defgrad-1+j+9)=hsv(no_hsvs+j+9)
          enddo
          endif
        endif
c
        d3(i)  =eps(3)
c
        if (ifipts(mt).eq.0) failel=.false.
c
        if (failel) then
          sig1(i)=0.
          sig2(i)=0.
          sig3(i)=0.
          sig4(i)=0.
          sig5(i)=0.
          sig6(i)=0.
          if (ifipts(mt).lt.0) hsvs(i,ind_fipts)=1.
          nfipts(i)=nfipts(i)+1
          if (ifipts(mt).ge.0) then
          failur=.true.
          ifail(i)=1
          fail(i)=0.
          elseif (nfipts(i).ge.nint(cm(mx-ifipts(mt)))) then
          failur=.true.
          ifail(i)=1
          fail(i)=0.            
          endif          
        else
          sig1(i)=sig(1)
          sig2(i)=sig(2)
          sig3(i)=sig(3)
          sig4(i)=sig(4)
          sig5(i)=sig(5)
          sig6(i)=sig(6)
        endif
        epsps(i)=epsp
 120    continue
c
c     vector umat
c
      else
c
c       assume no failed elements
c
        if (ifipts(mt).ge.0) then
        do i=lft,llt
          failels(i)=.false.
        enddo
        else
          do i=lft,llt
            failels(i)=hsvs(i,ind_fipts).ne.0.
          enddo
        endif
c
c       call user developed subroutines here
c       first history variable is the one requested
c
        if(ivumat(mt).eq.0) then
c
        if (mt.eq.41) then
          call umat41v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),qmat,
     .     elsizv,idelev,reject)
        elseif (mt.eq.42) then
          call umat42v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),qmat,
     .     elsizv,idelev,reject)
        elseif (mt.eq.43) then
          call umat43v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),qmat,
     .     elsizv,idelev,reject)
        elseif (mt.eq.44) then
          call umat44v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),qmat,
     .     elsizv,idelev,reject)
        elseif (mt.eq.45) then
          call umat45v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),qmat,
     .     elsizv,idelev,ihyper(mt),reject)
        elseif (mt.eq.46) then
          call umat46v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,nnpcrv,a(lcma),
     .     thhsv,nthhsv,qmat,elsizv,idelev,reject)
        elseif (mt.eq.47) then
          call umat47v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),qmat,elsizv,idelev,reject)
        elseif (mt.eq.48) then
          call umat48v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,ipt,a(n8),a(n9),a(lcma),qmat,elsizv,
     .     idelev,reject)
        elseif (mt.eq.49) then
          call umat49v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,ipt,temps,a(lcma),qmat,elsizv,idelev,
     .     reject)
        elseif (mt.eq.50) then
          call umat50v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),qmat,
     .     elsizv,idelev,reject)
c
        elseif (mt.eq.281) then
          call mat281s(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
        elseif (mt.eq.282) then
          call mat282s(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
c       elseif (mt.eq.283) then
c         call umat283s(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
c    .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
c    .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,
c    .     a(lcma))
        elseif (mt.eq.284) then
          call mat284s(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
        elseif (mt.eq.285) then
          call mat285s(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
        elseif (mt.eq.286) then
          call mat286s(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
        elseif (mt.eq.287) then
          call mat287s(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,failels,nlq,crv,
     .    ipt, npc, plc, nnm1,a(lcma))
        endif
c
        endif
c
        if (ivumat(mt).ne.0) then 
          ndi=3
c         Set nsh=3 to allow for updating shear components inside VUMAT
c         TODO: Should use a flag to choose between nsh=1 and nsh=3
c         nsh=1
          nsh=3
          nblock=llt-lft+1
          curveid=ivumat(mt)
          icid=int(curveid)
          if (icid.lt.0) then
            iicid=lcids(-icid)
            k2=npc(iicid)
            k3=npc(iicid+1)
            nprops=nint((k3-k2)/2.0)
          else
            nprops=40+icmadd(mt)
          endif

c          if (eltype.eq.'shell') if13=ind_q1
          call ext_vecumat(lft,llt,cm,a(lcma),bqs,mt,crv,npc,plc,nnm1,
     $         curveid,rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,
     $         nprops,ndi,nsh, if11,if22,if33,if12,if21,if23,if32,if13
     $         ,if31,no_hsvs, ifo11,ifo22,ifo33,ifo12,ifo21,ifo23,ifo32
     $         ,ifo13,ifo31,ind_q1)

c          call ext_vecumat(lft,llt,cm,bqs,mt,crv,npc,plc,nnm1,curveid,
c     $         rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,nprops,
c     $         ndi,nsh,if11,if22,if33,if12,if21,if23,if32,if13,if31)
        endif
c
c       check for failure
c
        if (ifipts(mt).ne.0) then
        do i=lft,llt
          if (failels(i)) then
            sig1(i)=0.
            sig2(i)=0.
            sig3(i)=0.
            sig4(i)=0.
            sig5(i)=0.
            sig6(i)=0.
            if (ifipts(mt).lt.0) hsvs(i,ind_fipts)=1.
            nfipts(i)=nfipts(i)+1
            if (ifipts(mt).ge.0) then
            failur=.true.
            ifail(i)=1
            fail(i)=0.
            elseif (nfipts(i).ge.nint(cm(mx-ifipts(mt)))) then
            failur=.true.
            ifail(i)=1
            fail(i)=0.
            endif
           endif
        enddo
      endif
      endif
c
      if (iabs(ihyper(mt)).eq.4) then
         do i=lft,llt
            d3(i)=0.
         enddo
      endif
c
c     energy calculations
c
      do 130 i=lft,llt
      einc(i)=d1(i)*sig1(i)+d2(i)*sig2(i)+d4(i)*sig4(i)+d5(i)*sig5(i)
     1       +d6(i)*sig6(i)+einc(i)
 130  continue
c
      if (iorien(mt).ne.0) then
c
c       transform stresses to element system
c
        if (ihyper(mt).ge.0) then
          do 140 i=lft,llt
          stg5(i)=sig5(i)
          stg6(i)=sig6(i)
          a11(i) = sig1(i)*hsvs(i,ind_q1)+sig4(i)*hsvs(i,ind_q2)
          a12(i) =-sig1(i)*hsvs(i,ind_q2)+sig4(i)*hsvs(i,ind_q1)
          a21(i) = sig4(i)*hsvs(i,ind_q1)+sig2(i)*hsvs(i,ind_q2)
          a22(i) =-sig4(i)*hsvs(i,ind_q2)+sig2(i)*hsvs(i,ind_q1)
          sig1(i)= hsvs(i,ind_q1)*a11(i)+hsvs(i,ind_q2)*a21(i)
          sig2(i)=-hsvs(i,ind_q2)*a12(i)+hsvs(i,ind_q1)*a22(i)
          sig4(i)= hsvs(i,ind_q1)*a12(i)+hsvs(i,ind_q2)*a22(i)
          sig5(i)=-hsvs(i,ind_q2)*stg6(i)+hsvs(i,ind_q1)*stg5(i)
          sig6(i)= hsvs(i,ind_q1)*stg6(i)+hsvs(i,ind_q2)*stg5(i)
 140      continue
        endif
c
c       transform deformation gradient back to element system
c
      if (ihyper(mt).gt.0) then
c        if (ihyper(mt).eq.1) then
          do i=lft,llt
            ftemp=hsvs(i,ind_q1)*hsvs(i,if11)+
     1            hsvs(i,ind_q2)*hsvs(i,if21)
            hsvs(i,if21)=-hsvs(i,ind_q2)*hsvs(i,if11)+
     1                    hsvs(i,ind_q1)*hsvs(i,if21)
            hsvs(i,if11)=ftemp
            ftemp=hsvs(i,ind_q1)*hsvs(i,if12)+
     1            hsvs(i,ind_q2)*hsvs(i,if22)
            hsvs(i,if22)=-hsvs(i,ind_q2)*hsvs(i,if12)+
     1                    hsvs(i,ind_q1)*hsvs(i,if22)
            hsvs(i,if12)=ftemp
            ftemp=hsvs(i,ind_q1)*hsvs(i,if13)+
     1            hsvs(i,ind_q2)*hsvs(i,if23)
            hsvs(i,if23)=-hsvs(i,ind_q2)*hsvs(i,if13)+
     1                    hsvs(i,ind_q1)*hsvs(i,if23)
            hsvs(i,if13)=ftemp
          enddo
c
c
c        transform deformation gradient back to element system
c
        else if (ihyper(mt).eq.-2) then
          do i=lft,llt
            f11=hsvs(i,ind_q1)*hsvs(i,if11)+
     1            hsvs(i,ind_q2)*hsvs(i,if21)
            f21=-hsvs(i,ind_q2)*hsvs(i,if11)+
     1                    hsvs(i,ind_q1)*hsvs(i,if21)
            f31=hsvs(i,if31)
c            hsvs(i,if11)=ftemp
            f12=hsvs(i,ind_q1)*hsvs(i,if12)+
     1            hsvs(i,ind_q2)*hsvs(i,if22)
            f22=-hsvs(i,ind_q2)*hsvs(i,if12)+
     1                    hsvs(i,ind_q1)*hsvs(i,if22)
            f32=hsvs(i,if32)
c            hsvs(i,if12)=ftemp
            f13=hsvs(i,ind_q1)*hsvs(i,if13)+
     1            hsvs(i,ind_q2)*hsvs(i,if23)
            f23=-hsvs(i,ind_q2)*hsvs(i,if13)+
     1                    hsvs(i,ind_q1)*hsvs(i,if23)
            f33=hsvs(i,if33)
c            hsvs(i,if13)=ftemp
            hsvs(i,if11)=f11*hsvs(i,ind_q1)+f12*hsvs(i,ind_q2)
            hsvs(i,if12)=-f11*hsvs(i,ind_q2)+f12*hsvs(i,ind_q1)
            hsvs(i,if13)=f13
            hsvs(i,if21)=f21*hsvs(i,ind_q1)+f22*hsvs(i,ind_q2)
            hsvs(i,if22)=-f21*hsvs(i,ind_q2)+f22*hsvs(i,ind_q1)
            hsvs(i,if23)=f23
            hsvs(i,if31)=f31*hsvs(i,ind_q1)+f32*hsvs(i,ind_q2)
            hsvs(i,if32)=-f31*hsvs(i,ind_q2)+f32*hsvs(i,ind_q1)
            hsvs(i,if33)=f33
          enddo
c
c
c        transform deformation gradient back to element system
c
c         elseif (ihyper(mt).lt.0) then
         elseif (ihyper(mt).eq.-1) then
           do i=lft,llt
            ftemp=hsvs(i,ind_q1)*hsvs(i,if11)+
     1            hsvs(i,ind_q2)*hsvs(i,if12)
            hsvs(i,if12)=-hsvs(i,ind_q2)*hsvs(i,if11)+
     1                    hsvs(i,ind_q1)*hsvs(i,if12)
            hsvs(i,if11)=ftemp
            ftemp=hsvs(i,ind_q1)*hsvs(i,if21)+
     1            hsvs(i,ind_q2)*hsvs(i,if22)
            hsvs(i,if22)=-hsvs(i,ind_q2)*hsvs(i,if21)+
     1                    hsvs(i,ind_q1)*hsvs(i,if22)
            hsvs(i,if21)=ftemp
            ftemp=hsvs(i,ind_q1)*hsvs(i,if31)+
     1            hsvs(i,ind_q2)*hsvs(i,if32)
            hsvs(i,if32)=-hsvs(i,ind_q2)*hsvs(i,if31)+
     1                    hsvs(i,ind_q1)*hsvs(i,if32)
            hsvs(i,if31)=ftemp
c
            enddo
         endif
      endif
c     
c     flag if material model rejected implicit iterate
      if (isolvr(1).ne.0.and.reject) isolvr(78)=5
c
      return
c900  format(/
c    1 ' *** Error User defined material in part',i12,
c    2 ' exceeds the current limit'/
c    3 10x,'of history variables for shells.')
      end
      subroutine usrshl_updatfs(fs,ts,ss,es)
      implicit none
      include 'nlqparm'
C_TASKCOMMON (bk36loc)
      integer index
      common/bk36loc/index
      real fs(3,3),ts(4),ss(4),es
      real f(nlq,3,3),t(nlq,4),s(nlq,4),e(nlq)
      integer i,j
      e(index)=es
      do i=1,4
         t(index,i)=ts(i)
         do j=1,nlq
            s(j,i)=1.
         enddo
      enddo
      call usrshl_updatfv(f,t,s,e,index,index)
      do i=1,4
         ss(i)=s(index,i)
      enddo
      do i=1,3
         do j=1,3
            fs(i,j)=f(index,i,j)
         enddo
      enddo
      return
      end
      subroutine usrshl_updatfv(f,t,s,e,lft,llt)
      implicit none
      include 'nlqparm'
      real f(nlq,3,3),t(nlq,4),s(nlq,4),e(nlq)
      integer lft,llt
      integer i,j
      real q
      do i=lft,llt
         e(i)=min(1.9999,e(i))
         q=(1.+.5*e(i))/(1.-.5*e(i))
         do j=1,4
            s(i,j)=max(t(i,j)*q,1.e-20)
         enddo
      enddo
      call usrshl_defgrad(f,3,
     1     0.,0.,0.,     
     2     0.d0,0.d0,
     3     0.,s,
     4     s,
     5     lft,llt)
      return
      end
      subroutine urmatb (lft,llt,cm,capa,mt,crv,nnpcrv,
     $     nnm1,nconstp,npc,plc,
     $     eltype,nbmbwp)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
#include "bk07.inc"
#include "bk13.inc"
      include 'bk19.inc'
      include 'memaia.inc'
      include 'umatss.inc'
#include "dynmem.inc"
#include "implicit1.inc"
C_TASKCOMMON (aux2loc)
C_TASKCOMMON (aux14loc)
C_TASKCOMMON (aux33loc)
C_TASKCOMMON (aux35loc)
C_TASKCOMMON (bux13loc)
C_TASKCOMMON (failcmloc)
C_TASKCOMMON (failuloc)
C_TASKCOMMON (soundloc)
C_TASKCOMMON (subtssloc)
      common/aux2loc/d1(nlq),d2(nlq),d3(nlq),d4(nlq),d5(nlq),d6(nlq),
     1 wzzdt(nlq),wyydt(nlq),wxxdt(nlq),einc(nlq)
      common/aux14loc/
     1 sig1(nlq),sig2(nlq),sig3(nlq),sig4(nlq),
     2 sig5(nlq),sig6(nlq),epsps(nlq),hsvs(nlq,NHISVAR)
      common/aux33loc/
     1 ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq),ix5(nlq),mxt(nlq)
      common/aux35loc/rhoa(nlq),cb(nlq),davg(nlq),p(nlq)
      common/bux13loc/eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      common/failcmloc/ifail(nlq)
      common/failuloc/sieu(nlq),fail(nlq)
      common/soundloc/ss(nlq),sndsp(nlq),diagm(nlq),sarea(nlq),
     . dxl(nlq)
      common/subtssloc/dt1siz(nlq)
c
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,nusa
      common/bk01/itherm,itemp,ntmp0,ntmp1,itempan,itempdr,itmpe(44)
      common/bk28/summss,xke,xpe,tt,xte0,erodeke,erodeie
      logical failur
      common/failcm/failur
c
      common/thropt/itopaz(101)
c
      dimension cm(*),crv(lq1,2,*),nconstp(*),npc(*),plc(*)
      integer nnpcrv(*)
      character*(*)eltype
      dimension nbmbwp(*)
c     dimension eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      dimension elsizv(nlq),idelev(nlq)
      logical failel,failels(nlq),reject
      INTEGER8 lqfinv8,idele
c
c     data temps/nlq*0.0/
      data failels/nlq*.false./
c
      failel=.false.
c
      reject=.false.
      ivect=ivectr(mt)
      mx=48*(mxt(lft)-1)
      gm=cm(mx+abs(ishrmp(mt)))
      bk=cm(mx+ibulkp(mt))
      pr=(3.*bk-2.*gm)/(2.*(3.*bk+gm))
      ym=2.*(1.+pr)*gm
      ss(lft)=ym

      do 10 i=lft,llt
      einc(i)=d1(i)*sig1(i)+d4(i)*sig4(i)+d6(i)*sig6(i)
      cb(i)=ss(lft)
   10 continue
c
c     Determining number of requested history variables
c
      no_hsvs=nconstp(mxt(lft))
c
c     if ishrmp(mt).lt.0 then the flag for temperature is active
c     temperatures are stored in the temps array.
c
      if (ishrmp(mt).lt.0) then
        num_nods=2
        call usr_temps (a(ntmp0+1),a(n19),itemp,temps,lft,llt,ix1,
     .   ix2,ix3,ix4,ix4,ix4,ix4,ix4,num_nods,itopaz,0,0.,0.,0.)
      else
        do i=lft,llt
          temps(i)=0.0
        enddo
      endif
c
c     characteristic element size and element id
      do i=lft,llt
        elsizv(i)=diagm(i)
        ielm=nbmbwp(nnm1+i)
        if (ielm.gt.0) idelev(i)=lqfinv8(ielm,3)
      enddo
c
c     memory pointer for extra memory data
      lcma=ncma+nmmat+ia(ncma+mxt(lft)-1)-1
c
      mte40=mt-40
c
      if (ivect.eq.0) then
        do 120 i=lft,llt
        dt1   =dt1siz(i)
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
        epsp=epsps(i)
        temper=temps(i)
        if (no_hsvs.gt.0) then
          do 20 j=1,no_hsvs
          hsv(j)=hsvs(i,j)
 20       continue
        endif
        elsiz=elsizv(i)
        idele=idelev(i)
c
c       call user developed subroutines here
c
        go to (41,42,43,44,45,46,47,48,49,50), mte40
   41   call umat41 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   42   call umat42 (cm(mx+1),eps,sig,hsv,dt1,capa,'beam ',tt,crv,
     .   nnpcrv,
     1   a(lcma),temper,0,elsiz,idele,reject)
        go to 60
   43   call umat43 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   44   call umat44 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   45   call umat45 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,0,reject)
        go to 60
   46   call umat46 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0.0,0,0,elsiz,idele,reject)
        go to 60
   47   call umat47 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   48   call umat48 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   49   call umat49 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   50   call umat50 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'beam ',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
 60     continue
c
        if (no_hsvs.gt.0) then
          do 100 j=1,no_hsvs
          hsvs(i,j)=hsv(j)
 100      continue
        endif
        d2(i)  =eps(2)
        d3(i)  =eps(3)
        epsps(i)=epsp
        if (failel) then
          sig1(i)=0.
          sig2(i)=0.
          sig3(i)=0.
          sig4(i)=0.
          sig5(i)=0.
          sig6(i)=0.
          failur=.true.
          ifail(i)=1
          fail(i)=0.
        else
          sig1(i)=sig(1)
          sig2(i)=sig(2)
          sig3(i)=sig(3)
          sig4(i)=sig(4)
          sig5(i)=sig(5)
          sig6(i)=sig(6)
        endif
 120    continue
c
c     vector umat
c
      else
c
        if (ivumat(mt).eq.0) then

        if (mt.eq.41) then
          call umat41v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.42) then
          call umat42v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.43) then
          call umat43v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.44) then
          call umat44v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .   sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .   'beam ',tt,temps,failels,nlq,crv,nnpcrv,
     .   a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.45) then
          call umat45v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,0,reject)
        elseif (mt.eq.46) then
          call umat46v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0.0,0,0,elsizv,
     .     idelev,reject)
        elseif (mt.eq.47) then
          call umat47v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.48) then
          call umat48v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,ipt,a(n8),a(n9),a(lcma),0,elsizv,idelev,
     .     reject)
        elseif (mt.eq.49) then
          call umat49v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,ipt,temps,a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.50) then
          call umat50v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'beam ',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,reject)
        endif
        endif
      
        if (ivumat(mt).ne.0) then
          ndi=1
          nsh=1
          nblock=llt-lft+1
          curveid=ivumat(mt)
          icid=int(curveid)
          if (icid.lt.0) then
            iicid=lcids(-icid)
            k2=npc(iicid)
            k3=npc(iicid+1)
            nprops=nint((k3-k2)/2.0)
          else
            nprops=40+icmadd(mt)
          endif
          if11=0.0
          if22=0.0
          if33=0.0
          if12=0.0
          if21=0.0
          if23=0.0
          if32=0.0
          if13=0.0
          if31=0.0
          ifo11=0.0
          ifo22=0.0
          ifo33=0.0
          ifo12=0.0
          ifo21=0.0
          ifo23=0.0
          ifo32=0.0
          ifo13=0.0
          ifo31=0.0

          call ext_vecumat(lft,llt,cm,a(lcma),bqs,mt,crv,npc,plc,nnm1,
     $         curveid,rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,
     $         nprops,ndi,nsh, if11,if22,if33,if12,if21,if23,if32,if13
     $         ,if31,no_hsvs, ifo11,ifo22,ifo33,ifo12,ifo21,ifo23,ifo32
     $         ,ifo13,ifo31,tmp)

c          call ext_vecumat(lft,llt,cm,bqs,mt,crv,npc,plc,nnm1,curveid,
c     $         rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,nprops,
c     $         ndi,nsh,if11,if22,if33,if12,if21,if23,if32,if13,if31)
        endif

      endif
c
      do 130 i=lft,llt
      einc(i)=d1(i)*sig1(i)+d4(i)*sig4(i)+d6(i)*sig6(i)+einc(i)
 130  continue
c
c     flag if material model rejected implicit iterate
      if (isolvr(1).ne.0.and.reject) isolvr(78)=5
c
      return
      end
      subroutine urmatd (lft,llt,cm,capa,mt,crv,nnpcrv,
     $     nnm1,nconstp,npc,plc,
     $     eltype,nbmbwp)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     For now, call only vectorized umat routines for discrte beam.
c     umat47v includes curves
c
      include 'nlqparm'
      include 'nhisparm.inc'
#include "bk07.inc"
#include "bk13.inc"
      include 'bk19.inc'
      include 'memaia.inc'
      include 'umatss.inc'
#include "dynmem.inc"
#include "implicit1.inc"
C_TASKCOMMON (aux2loc)
C_TASKCOMMON (aux14loc)
C_TASKCOMMON (aux33loc)
C_TASKCOMMON (aux35loc)
C_TASKCOMMON (bux13loc)
C_TASKCOMMON (soundloc)
C_TASKCOMMON (subtssloc)
      common/aux2loc/d1(nlq),d2(nlq),d3(nlq),d4(nlq),d5(nlq),d6(nlq),
     1 wzzdt(nlq),wyydt(nlq),wxxdt(nlq),einc(nlq)
      common/aux14loc/
     1 sig1(nlq),sig2(nlq),sig3(nlq),sig4(nlq),
     2 sig5(nlq),sig6(nlq),epsps(nlq),hsvs(nlq,NHISVAR)
      common/aux33loc/
     1 ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq),ix5(nlq),mxt(nlq)
      common/aux35loc/rhoa(nlq),cb(nlq),davg(nlq),p(nlq)
      common/bux13loc/ eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      common/soundloc/ss(nlq),sndsp(nlq),diagm(nlq),sarea(nlq),
     . dxl(nlq)
      common/subtssloc/dt1siz(nlq)
c
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,nusa
      common/bk01/itherm,itemp,ntmp0,ntmp1,itempan,itempdr,itmpe(44)
      common/bk28/summss,xke,xpe,tt,xte0,erodeke,erodeie
      include 'slcntr.inc'
      common/thropt/itopaz(101)
c
c     capa here is ies in conmd6 and bqs in usrmat (all are internal energy)
c     real capa
      dimension capa(*),cm(*),crv(lq1,2,*),nconstp(*),npc(*),plc(*)
      integer nnpcrv(*)
c     dimension eps(6),sig(6),hsv(6),temps(nlq)
      dimension nbmbwp(*)
      character*(*)eltype
      dimension elsizv(nlq)
      logical failel,failels(nlq),reject
      INTEGER8 lqfinv8
      INTEGER8 idelev(nlq),idele
c
c     data temps/nlq*0.0/
      data failels/nlq*.false./
c
      failel=.false.
      reject=.false.
      ivect=ivectr(mt)
      mx=48*(mxt(lft)-1)
c
c     Determining number of requested history variables
c
      no_hsvs=nconstp(mxt(lft))
c
c     if ishrmp(mt).lt.0 then the flag for temperature is active
c     temperatures are stored in the temps array.
c
      if (ishrmp(mt).lt.0) then
        num_nods=2
        call usr_temps (a(ntmp0+1),a(n19),itemp,temps,lft,llt,
     .   ix1,ix2,ix3,ix4,ix4,ix4,ix4,ix4,num_nods,itopaz,0,
     .   0.,0.,0.)
      else
        do i=lft,llt
          temps(i)=0.0
        enddo
      endif
c
c     characteristic element size and element id
      do i=lft,llt
        elsizv(i)=diagm(i)
        ielm=nbmbwp(nnm1+i)
        if (ielm.gt.0) idelev(i)=lqfinv8(ielm,3)
      enddo
c
c     memory pointer for extra memory data
      lcma=ncma+nmmat+ia(ncma+mxt(lft)-1)-1
c
      mte40=mt-40
c
      if (ivect.eq.0) then
        do 120 i=lft,llt
        dt1   =dt1siz(i)
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=wxxdt(i)
        eps(5)=wyydt(i)
        eps(6)=wzzdt(i)
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
        temper=temps(i)
        epsp=epsps(i)
        if (no_hsvs.gt.0) then
          do 20 j=1,no_hsvs
          hsv(j)=hsvs(i,j)
 20       continue
        endif
        elsiz=elsizv(i)
        idele=idelev(i)
c
c       call user developed subroutines here
c
        go to (41,42,43,44,45,46,47,48,49,50), mte40
   41   call umat41 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   42   call umat42 (cm(mx+1),eps,sig,hsv,dt1,capa,'dbeam',tt,crv,
     .   nnpcrv,
     1   a(lcma),temper,0,elsiz,idele,reject)
        go to 60
   43   call umat43 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   44   call umat44 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   45   call umat45 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,0,reject)
        go to 60
   46   call umat46 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,nnpcrv,0.0,0,a(lcma),0,elsiz,idele,reject)
        go to 60
   47   call umat47 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   48   call umat48 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   49   call umat49 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   50   call umat50 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'dbeam',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
 60     continue
c
        if (no_hsvs.gt.0) then
          do 100 j=1,no_hsvs
          hsvs(i,j)=hsv(j)
 100      continue
        endif
        sig1(i)=sig(1)
        sig2(i)=sig(2)
        sig3(i)=sig(3)
        sig4(i)=sig(4)
        sig5(i)=sig(5)
        sig6(i)=sig(6)
        epsps(i)=epsp
 120  continue
c
c     vector umat
c
      else

        if (ivumat(mt).eq.0) then

        if (mt.eq.41) then
          call umat41v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.42) then
          call umat42v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.43) then
          call umat43v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.44) then
          call umat44v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.45) then
          call umat45v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,0,reject)
        elseif (mt.eq.46) then
          call umat46v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0.0,0,0,elsizv,
     .     idelev,reject)
        elseif (mt.eq.47) then
          call umat47v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.48) then
          call umat48v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,ipt,a(n8),a(n9),a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.49) then
          call umat49v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,ipt,temps,a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.50) then
          call umat50v(cm(mx+1),d1,d2,d3,wxxdt,wyydt,wzzdt,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'dbeam',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,reject)
        endif
        endif

        if (ivumat(mt).ne.0) then
          ndi=1
          nsh=1
          nblock=llt-lft+1
          curveid=ivumat(mt)
          icid=int(curveid)
          if (icid.lt.0) then
            iicid=lcids(-icid)
            k2=npc(iicid)
            k3=npc(iicid+1)
            nprops=nint((k3-k2)/2.0)
          else
            nprops=40+icmadd(mt)
          endif
          call ext_vecumat(lft,llt,cm,a(lcma),bqs,mt,crv,npc,plc,nnm1,
     $         curveid,rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,
     $         nprops,ndi,nsh, if11,if22,if33,if12,if21,if23,if32,if13
     $         ,if31,no_hsvs, ifo11,ifo22,ifo33,ifo12,ifo21,ifo23,ifo32
     $         ,ifo13,ifo31,tmp)

c          call ext_vecumat(lft,llt,cm,bqs,mt,crv,npc,plc,nnm1,curveid,
c     $         rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,nprops,
c     $         ndi,nsh,if11,if22,if33,if12,if21,if23,if32,if13,if31)
        endif
      
      endif
c
      do 130 i=lft,llt
      einc(i)=d1(i)*sig1(i)+d4(i)*sig4(i)+d6(i)*sig6(i)+einc(i)
  130 continue
c
c     flag if material model rejected implicit iterate
      if (isolvr(1).ne.0.and.reject) isolvr(78)=5
c
      return
      end
      subroutine urmatt (lft,llt,cm,capa,mt,crv,nnpcrv,
     $     nnm1,nconstp,npc,plc,
     $     eltype,nbmbwp)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
#include "bk07.inc"
#include "bk13.inc"
      include 'bk19.inc'
      include 'memaia.inc'
      include 'umatss.inc'
#include "dynmem.inc"
#include "implicit1.inc"
C_TASKCOMMON (aux2loc)
C_TASKCOMMON (aux14loc)
C_TASKCOMMON (aux33loc)
C_TASKCOMMON (bux13loc)
C_TASKCOMMON (aux35loc)
C_TASKCOMMON (failcmloc)
C_TASKCOMMON (failuloc)
C_TASKCOMMON (soundloc)
C_TASKCOMMON (subtssloc)
      common/aux2loc/d1(nlq),d2(nlq),d3(nlq),d4(nlq),d5(nlq),d6(nlq),
     1 wzzdt(nlq),wyydt(nlq),wxxdt(nlq),einc(nlq)
      common/aux14loc/
     1 sig1(nlq),sig2(nlq),sig3(nlq),sig4(nlq),
     2 sig5(nlq),sig6(nlq),epsps(nlq),hsvs(nlq,NHISVAR)
      common/aux33loc/
     1 ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq),ix5(nlq),mxt(nlq)
      common/aux35loc/rhoa(nlq),cb(nlq),davg(nlq),p(nlq)
      common/bux13loc/ eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      common/failcmloc/ifail(nlq)
      common/failuloc/sieu(nlq),fail(nlq)
      common/soundloc/ss(nlq),sndsp(nlq),diagm(nlq),sarea(nlq),
     . dxl(nlq)
      common/subtssloc/dt1siz(nlq)
c
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,nusa
      common/bk01/itherm,itemp,ntmp0,ntmp1,itempan,itempdr,itmpe(44)
      common/bk28/summss,xke,xpe,tt,xte0,erodeke,erodeie
      logical failur
      common/failcm/failur
      common/thropt/itopaz(101)
      include 'slcntr.inc'
c
c     real capa
      dimension capa(*),cm(*),crv(lq1,2,*),nconstp(*),npc(*),plc(*)
      integer nnpcrv(*)
c     dimension eps(6),sig(6),hsv(6),temps(nlq)
      character*(*)eltype
      dimension nbmbwp(*)
      dimension elsizv(nlq)
      logical failel,failels(nlq),reject
      INTEGER8 lqfinv8
      INTEGER8 idelev(nlq),idele
c
c     data temps/nlq*0.0/
      data failels/nlq*.false./
c
      reject=.false.
      ivect=ivectr(mt)
      mx=48*(mxt(lft)-1)
      gm=cm(mx+abs(ishrmp(mt)))
      bk=cm(mx+ibulkp(mt))
      ym=9.*bk*gm/(3.*bk+gm)
      ss(lft)=ym

      do 10 i=lft,llt
      einc(i)=d1(i)*sig1(i)
      cb(i)=ss(lft)
      wxxdt(i)=ym
   10 continue
c
c     Determining number of requested history variables
c
      no_hsvs=nconstp(mxt(lft))
c
c     if ishrmp(mt).lt.0 then the flag for temperature is active
c     temperatures are stored in the temps array.
c
      if (ishrmp(mt).lt.0) then
        num_nods=2
        call usr_temps (a(ntmp0+1),a(n19),itemp,temps,lft,llt,
     .   ix1,ix2,ix3,ix4,ix4,ix4,ix4,ix4,num_nods,itopaz,
     .   0,0.,0.,0.)
      else
        do i=lft,llt
          temps(i)=0.0
        enddo
      endif
c
c     characteristic element size
      do i=lft,llt
        elsizv(i)=diagm(i)
        ielm=nbmbwp(nnm1+i)
        if (ielm.gt.0) idelev(i)=lqfinv8(ielm,3)
      enddo
c
c     memory pointer for extra memory data
      lcma=ncma+nmmat+ia(ncma+mxt(lft)-1)-1
c
      mte40=mt-40
c
      if (ivect.eq.0) then
        do 120 i=lft,llt
        dt1   =dt1siz(i)
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        epsp=epsps(i)
        temper=temps(i)
        failel=.false.
        if (no_hsvs.gt.0) then
          do 20 j=1,no_hsvs
          hsv(j)=hsvs(i,j)
 20       continue
        endif
        elsiz=elsizv(i)
        idele=idelev(i)
c
c       call user developed subroutines here
c
        go to (41,42,43,44,45,46,47,48,49,50), mte40
   41   call umat41 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   42   call umat42 (cm(mx+1),eps,sig,hsv,dt1,capa,'tbeam',tt,crv,
     .   nnpcrv,
     1   a(lcma),temper,0,elsiz,idele,reject)
        go to 60
   43   call umat43 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   44   call umat44 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   45   call umat45 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,0,reject)
        go to 60
   46   call umat46 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0.0,0,0,elsiz,idele,reject)
        go to 60
   47   call umat47 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   48   call umat48 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   49   call umat49 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
        go to 60
   50   call umat50 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'tbeam',tt,
     1   temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,reject)
 60   continue
c
      if (no_hsvs.gt.0) then
        do 100 j=1,no_hsvs
        hsvs(i,j)=hsv(j)
 100   continue
      endif
c
      if (failel) then
        einc(i)=0.
        sig1(i)=0.
        sig2(i)=0.
        sig3(i)=0.
        failur=.true.
        ifail(i)=1
        fail(i)=0.
      else
        sig1(i)=sig(1)
        sig2(i)=sig(2)
        sig3(i)=sig(3)
      endif
      d2(i)  =eps(2)
      d3(i)  =eps(3)
      epsps(i)=epsp
 120  continue
c
c     vector umat
c
      else

        if (ivumat(mt).eq.0) then

        if (mt.eq.41) then
          call umat41v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.42) then
          call umat42v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.43) then
          call umat43v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.44) then
          call umat44v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.45) then
          call umat45v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,0,reject)
        elseif (mt.eq.46) then
          call umat46v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0.0,0,0,elsizv,
     .     idelev,reject)
        elseif (mt.eq.47) then
          call umat47v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.48) then
          call umat48v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,ipt,a(n8),a(n9),a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.49) then
          call umat49v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,ipt,temps,a(lcma),0,elsizv,idelev,reject)
        elseif (mt.eq.50) then
          call umat50v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     .     'tbeam',tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),0,elsizv,idelev,reject)
        endif

        endif

        if (ivumat(mt).ne.0) then
          ndi=1
          nsh=1
          nblock=llt-lft+1
          curveid=ivumat(mt)
          icid=int(curveid)
          if (icid.lt.0) then
            iicid=lcids(-icid)
            k2=npc(iicid)
            k3=npc(iicid+1)
            nprops=nint((k3-k2)/2.0)
          else
            nprops=40+icmadd(mt)
          endif
          call ext_vecumat(lft,llt,cm,a(lcma),bqs,mt,crv,npc,plc,nnm1,
     $         curveid,rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,
     $         nprops,ndi,nsh, if11,if22,if33,if12,if21,if23,if32,if13
     $         ,if31,no_hsvs, ifo11,ifo22,ifo33,ifo12,ifo21,ifo23,ifo32
     $         ,ifo13,ifo31,tmp)

c          call ext_vecumat(lft,llt,cm,bqs,mt,crv,npc,plc,nnm1,curveid,
c     $         rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,nprops,
c     $         ndi,nsh,if11,if22,if33,if12,if21,if23,if32,if13,if31)
        endif
        
c
        do i=lft,llt
          if (failels(i)) then
            sig1(i)=0.
            sig2(i)=0.
            sig3(i)=0.
            fail(i)=0.
            ifail(i)=1
            failur=.true.
         endif
        enddo
      endif
c
      do 130 i=lft,llt
      einc(i)=d1(i)*sig1(i)+einc(i)
 130  continue
c
c     flag if material model rejected implicit iterate
      if (isolvr(1).ne.0.and.reject) isolvr(78)=5
c
      return
      end
      subroutine umat41 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     isotropic elastic material (sample user subroutine)
c
c     Variables
c
c     cm(1)=first material constant, here young's modulus
c     cm(2)=second material constant, here poisson's ratio
c        .
c        .
c        .
c     cm(n)=nth material constant
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
c     hsv(1)=1st history variable
c     hsv(2)=2nd history variable
c        .
c        .
c        .
c        .
c     hsv(n)=nth history variable
c
c     dt1=current time step size
c     capa=reduction factor for transverse shear
c     etype:
c        eq."solid" for solid elements
c        eq."sph  " for smoothed particle hydrodynamics
c        eq."sld2d" for shell forms 13 (2D solids - plane strain)
c        eq."sldax" for shell forms 14, and 15 (2D solids - axisymmetric)
c        eq."shl_t" for shell forms 25, 26, and 27 (shells with thickness stretch)
c        eq."shell" for all other shell elements plus thick shell forms 1 and 2
c        eq."tshel" for thick shell forms 3 and 5
c        eq."hbeam" for beam element forms 1 and 11
c        eq."tbeam" for beam element form 3 (truss)
c        eq."dbeam" for beam element form 6 (discrete)
c        eq."beam " for all other beam elements
c
c     tt=current problem time.
c
c     temper=current temperature
c
c     failel=flag for failure, set to .true. to fail an integration point,
c            if .true. on input the integration point has failed earlier
c
c     crv=array representation of curves in keyword deck
c
c     nnpcrv=# of discretization points per crv()
c
c     cma=additional memory for material data defined by LMCA at 
c       6th field of 2nd crad of *DATA_USER_DEFINED
c
c     elsiz=characteristic element size
c
c     idele=element id
c
c     reject (implicit only) = set to .true. if this implicit iterate is
c                              to be rejected for some reason
c
c     All transformations into the element local system are
c     performed prior to entering this subroutine.  Transformations
c     back to the global system are performed after exiting this
c     routine.
c
c     All history variables are initialized to zero in the input
c     phase. Initialization of history variables to nonzero values
c     may be done during the first call to this subroutine for each
c     element.
c
c     Energy calculations for the dyna3d energy balance are done
c     outside this subroutine.
c
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      integer nnpcrv(*)
      logical failel,reject
      character*5 etype
      INTEGER8 idele
c
      if (ncycle.eq.1) then
        if (cm(16).ne.1234567) then
          call usermsg('mat41')
        endif
      endif
c
c     compute shear modulus, g
c
      g2 =abs(cm(1))/(1.+cm(2))
      g  =.5*g2
c
      if (etype.eq.'solid'.or.etype.eq.'shl_t'.or.
     1     etype.eq.'sld2d'.or.etype.eq.'tshel'.or.
     2     etype.eq.'sph  '.or.etype.eq.'sldax') then
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
          gc    =capa*g
          q1    =abs(cm(1))*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
          q3    =1./(q1+g2)
          eps(3)=-q1*(eps(1)+eps(2))*q3
          davg  =(-eps(1)-eps(2)-eps(3))/3.
          p     =-davg*abs(cm(1))/(1.-2.*cm(2))
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
          q1    =cm(1)*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
          q3    =q1+2.0*g
          gc    =capa*g
          deti  =1./(q3*q3-q1*q1)
          c22i  = q3*deti
          c23i  =-q1*deti
          fac   =(c22i+c23i)*q1
          eps(2)=-eps(1)*fac-sig(2)*c22i-sig(3)*c23i
          eps(3)=-eps(1)*fac-sig(2)*c23i-sig(3)*c22i
          davg  =(-eps(1)-eps(2)-eps(3))/3.
          p     =-davg*cm(1)/(1.-2.*cm(2))
          sig(1)=sig(1)+p+g2*(eps(1)+davg)
          sig(2)=0.0
          sig(3)=0.0
          sig(4)=sig(4)+gc*eps(4)
          sig(5)=0.0
          sig(6)=sig(6)+gc*eps(6)

c
      elseif (etype.eq.'tbeam') then
        q1    =cm(1)*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
        q3    =q1+2.0*g
        deti  =1./(q3*q3-q1*q1)
        c22i  = q3*deti
        c23i  =-q1*deti
        fac   =(c22i+c23i)*q1
        eps(2)=-eps(1)*fac
        eps(3)=-eps(1)*fac
        davg  =(-eps(1)-eps(2)-eps(3))/3.
        p     =-davg*cm(1)/(1.-2.*cm(2))
        sig(1)=sig(1)+p+g2*(eps(1)+davg)
        sig(2)=0.0
        sig(3)=0.0
c
      else
c       write(iotty,10) etype
c       write(iohsp,10) etype
c       write(iomsg,10) etype
c       call adios(TC_ERROR)
        cerdat(1)=etype
        call lsmsg(3,MSG_SOL+1150,ioall,ierdat,rerdat,cerdat,0)
      endif
c
c10   format(/
c    1 ' *** Error element type ',a,' can not be',
c    2 '           run with the current material model.')
      return
      end
      subroutine umat45 (cm,eps,sig,epsp,hsv,dt1,capa,
     . etype,time,temp,failel,crv,nnpcrv,cma,qmat,elsiz,idele,ih,
     . reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     Neo-Hookean material (sample user subroutine)
c
c     Also serves as a test example for the *MODULE feature,
c     see comments at the beginning of the code body
c
c     Variables
c
c     cm(1)=first material constant, here young's modulus
c     cm(2)=second material constant, here poisson's ratio
c        .
c        .
c        .
c     cm(n)=nth material constant
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
c     hsv(1)=1st history variable
c     hsv(2)=2nd history variable
c        .
c        .
c        .
c        .
c     hsv(n)=nth history variable
c
c     dt1=current time step size
c     capa=reduction factor for transverse shear
c     etype:
c        eq."solid" for solid elements
c        eq."sph  " for smoothed particle hydrodynamics
c        eq."sld2d" for shell forms 13 (2D solids - plane strain)
c        eq."sldax" for shell forms 14, and 15 (2D solids)
c        eq."shl_t" for shell forms 25, 26, and 27 (shells with thickness stretch)
c        eq."shell" for all other shell elements plus thick shell forms 1 and 2
c        eq."tshel" for thick shell forms 3 and 5
c        eq."hbeam" for beam element forms 1 and 11
c        eq."tbeam" for beam element form 3 (truss)
c        eq."dbeam" for beam element form 6 (discrete)
c        eq."beam " for all other beam elements
c
c     time=current problem time.
c     temp=current temperature
c
c     cma=additional memory for material data defined by LMCA at 
c       6th field of 2nd crad of *DATA_USER_DEFINED
c
c     All transformations into the element local system are
c     performed prior to entering this subroutine.  Transformations
c     back to the global system are performed after exiting this
c     routine.
c
c     All history variables are initialized to zero in the input
c     phase.   Initialization of history variables to nonzero values
c     may be done during the first call to this subroutine for each
c     element.
c
c     Energy calculations for the dyna3d energy balance are done
c     outside this subroutine.
c
      include 'nlqparm'
      include 'iounits.inc'
      include 'bk06.inc'
      character*5 etype
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      integer nnpcrv(*)
      logical failel,reject
      INTEGER8 idele
      real newthick(4)
c
      logical usemodule,usemetal
c
c     As an example of the use of *MODULE, this routine can be executed
c     both as a hyperelastic Neo-Hookean material (this routine) and as
c     an elastic-plastic material (subroutines metal*). The purpose of
c     the *MODULE concept is to load several shared object files to a
c     single LS-DYNA executable and whence (in this particular case) 
c     allow multiple instances of subroutine umat45 be executed without
c     conflicts. This is important if independent user material distributors 
c     happen to use the same user material number for their respective
c     material implementations, and neither want to reveal the source code.
c
c     To use this example, perform the following steps in order
c
c     1) Set usemodule=.true. below
c     2) Compile the code to a shared object file and call this file 
c        hyperelastic.so
c     3) Set usemetal=.true. below
c     4) Recompile the code to a new shared object file and call this file
c        elasticplastic.so
c     5) In the keyword file, use *MODULE_LOAD for loading the two object
c        files and *MODULE_USE to associate hyperelastic materials to 
c        material 45 in hyperelastic.so and elastic-plastic materials to
c        material 45 in elasticplastic.so
c
c     Note that for implicit, utan45 needs to modified accordingly
c
c     For details, refer to Appendix A of the Keyword User's Manual
c     
      usemodule=.false.
      if (usemodule) then
        usemetal=.false.
        if (usemetal) then
          if (etype.eq.'solid') then
            call metalsld(cm,eps,sig,epsp,hsv,dt1,capa,etype,time,temp,
     .            failel,crv,nnpcrv)
          else
            call metalshl(cm,eps,sig,epsp,hsv,dt1,capa,etype,time,temp,
     .            failel,crv,nnpcrv)
          endif
          return
        endif
      endif
      if (.not.usemodule) then
      if (ncycle.eq.1) then
        call usermsg('mat45')
      endif
      endif
c
c     compute lame parameters
c
      xlambda=cm(1)*cm(2)/((1.+cm(2))*(1.-2.*cm(2)))
      xmu=.5*cm(1)/(1.+cm(2))
c
      if (etype.eq.'solid'.or.etype.eq.'shl_t'.or.
     1     etype.eq.'sld2d'.or.etype.eq.'tshel'.or.
     2     etype.eq.'sph  '.or.etype.eq.'sldax'.or.
     3     etype.eq.'shell'.and.iabs(ih).eq.3) then
c
        if (iabs(ih).eq.3) then           
c     sample routine for ihyper=3
c     initialize iteration counter and new thickness,
c     the four thickness variables are stored immediately
c     after the deformation gradient always
           iter=0
           newthick(1)=hsv(10)
           newthick(2)=hsv(11)
           newthick(3)=hsv(12)
           newthick(4)=hsv(13)
        endif
 1      if (iabs(ih).eq.3) then
c     sample routine for ihyper=3
c     transform deformation gradient to local system
           b1=qmat(1,1)*hsv(1)+qmat(2,1)*hsv(2)+qmat(3,1)*hsv(3)
           b2=qmat(1,2)*hsv(1)+qmat(2,2)*hsv(2)+qmat(3,2)*hsv(3)
           b3=qmat(1,3)*hsv(1)+qmat(2,3)*hsv(2)+qmat(3,3)*hsv(3)
           hsv(1)=b1
           hsv(2)=b2
           hsv(3)=b3
           b1=qmat(1,1)*hsv(4)+qmat(2,1)*hsv(5)+qmat(3,1)*hsv(6)
           b2=qmat(1,2)*hsv(4)+qmat(2,2)*hsv(5)+qmat(3,2)*hsv(6)
           b3=qmat(1,3)*hsv(4)+qmat(2,3)*hsv(5)+qmat(3,3)*hsv(6)
           hsv(4)=b1
           hsv(5)=b2
           hsv(6)=b3
           b1=qmat(1,1)*hsv(7)+qmat(2,1)*hsv(8)+qmat(3,1)*hsv(9)
           b2=qmat(1,2)*hsv(7)+qmat(2,2)*hsv(8)+qmat(3,2)*hsv(9)
           b3=qmat(1,3)*hsv(7)+qmat(2,3)*hsv(8)+qmat(3,3)*hsv(9)
           hsv(7)=b1
           hsv(8)=b2
           hsv(9)=b3
        endif
c
c       deformation gradient stored in hsv(1),...,hsv(9)
c
c       compute jacobian
c
        detf=hsv(1)*(hsv(5)*hsv(9)-hsv(6)*hsv(8))
     1      -hsv(2)*(hsv(4)*hsv(9)-hsv(6)*hsv(7))
     2      +hsv(3)*(hsv(4)*hsv(8)-hsv(5)*hsv(7))
c
c       compute left cauchy-green tensor
c
        b1=hsv(1)*hsv(1)+hsv(4)*hsv(4)+hsv(7)*hsv(7)
        b2=hsv(2)*hsv(2)+hsv(5)*hsv(5)+hsv(8)*hsv(8)
        b3=hsv(3)*hsv(3)+hsv(6)*hsv(6)+hsv(9)*hsv(9)
        b4=hsv(1)*hsv(2)+hsv(4)*hsv(5)+hsv(7)*hsv(8)
        b5=hsv(2)*hsv(3)+hsv(5)*hsv(6)+hsv(8)*hsv(9)
        b6=hsv(1)*hsv(3)+hsv(4)*hsv(6)+hsv(7)*hsv(9)
c
c       compute cauchy stress
c
        detf=max(detf,1.e-8)
        detfinv=1./detf
        dmu=xmu-xlambda*log(detf)
        sig(1)=detfinv*(xmu*b1-dmu)
        sig(2)=detfinv*(xmu*b2-dmu)
        sig(3)=detfinv*(xmu*b3-dmu)
        sig(4)=detfinv*xmu*b4
        sig(5)=detfinv*xmu*b5
        sig(6)=detfinv*xmu*b6
c
        if (iabs(ih).eq.3) then
c     sample routine for ihyper=3
           if (abs(sig(3)).le.1.e-4*max(1.e-4*cm(1),abs(sig(1)),
     1          abs(sig(2)),abs(sig(4)),abs(sig(5)),abs(sig(6)))) then
c     convergence in thickness iterations, update thickness and return
              hsv(10)=newthick(1)
              hsv(11)=newthick(2)
              hsv(12)=newthick(3)
              hsv(13)=newthick(4)
              return
           endif
           iter=iter+1
           if (iter.eq.100) then
c     divergence, update thickness and return
              hsv(10)=newthick(1)
              hsv(11)=newthick(2)
              hsv(12)=newthick(3)
              hsv(13)=newthick(4)
              return
           elseif (iter.eq.1) then
c     the first stress was evaluated for eps(3)=0., store this point
              sigold=sig(3)
              epsold=eps(3)
c     new thickness strain estimate
              eps(3)=-xlambda*(eps(1)+eps(2))/(xlambda+2.*xmu)
           elseif (abs(sig(3)-sigold).gt.0.0) then
c     secant increment of thickness strain
              deps=-(eps(3)-epsold)/(sig(3)-sigold)*sig(3)
c     store previous point
              sigold=sig(3)
              epsold=eps(3)
c     update strain
              eps(3)=eps(3)+deps
           endif
c     recalculate deformation gradient with new thickness strain increment
           call usrshl_updatfs(hsv,hsv(10),newthick,eps(3))
c     go back and reevaluate stress
           goto 1
        endif
c
      else if (etype.eq.'shell') then
c
c       deformation gradient stored in hsv(1),...,hsv(9)
c
c       compute part of left cauchy-green tensor
c       independent of thickness strain increment
c
        b1=hsv(1)*hsv(1)+hsv(4)*hsv(4)+hsv(7)*hsv(7)
        b2=hsv(2)*hsv(2)+hsv(5)*hsv(5)+hsv(8)*hsv(8)
        b4=hsv(1)*hsv(2)+hsv(4)*hsv(5)+hsv(7)*hsv(8)
c
c       secant iterations for zero normal stress
c
        do iter=1,5
c
c         first thickness strain increment initial guess
c         assuming Poisson's ratio different from zero
c
          if (iter.eq.1) then
            eps(3)=-xlambda*(eps(1)+eps(2))/(xlambda+2.*xmu)
c
c         second thickness strain increment initial guess
c
          else if (iter.eq.2) then
            sigold=sig(3)
            epsold=eps(3)
            eps(3)=0.
c
c         secant update of thickness strain increment
c
          else if (abs(sig(3)-sigold).gt.0.0) then
            deps=-(eps(3)-epsold)/(sig(3)-sigold)*sig(3)
            sigold=sig(3)
            epsold=eps(3)
            eps(3)=eps(3)+deps
          endif
c
c         compute last row of deformation gradient
c
          call compute_f3s(hsv(3),hsv(6),hsv(9),eps(3))
c
c         compute jacobian
c
          detf=hsv(1)*(hsv(5)*hsv(9)-hsv(6)*hsv(8))
     1        -hsv(2)*(hsv(4)*hsv(9)-hsv(6)*hsv(7))
     2        +hsv(3)*(hsv(4)*hsv(8)-hsv(5)*hsv(7))
c
c         compute normal component of left cauchy-green tensor
c
          b3=hsv(3)*hsv(3)+hsv(6)*hsv(6)+hsv(9)*hsv(9)
c
c         compute normal stress
c
          detf=max(detf,1.e-8)
          detfinv=1./detf
          dmu=xmu-xlambda*log(detf)
          sig(1)=detfinv*(xmu*b1-dmu)
          sig(2)=detfinv*(xmu*b2-dmu)
          sig(3)=detfinv*(xmu*b3-dmu)
          sig(4)=detfinv*xmu*b4
c
c         exit if normal stress is sufficiently small
c
          if (abs(sig(3)).le.1.e-5*
     1     (abs(sig(1))+abs(sig(2))+abs(sig(4)))) goto 10
         enddo
c
c       compute remaining components of left cauchy-green tensor
c
 10     b5=hsv(2)*hsv(3)+hsv(5)*hsv(6)+hsv(8)*hsv(9)
        b6=hsv(1)*hsv(3)+hsv(4)*hsv(6)+hsv(7)*hsv(9)
c
c       compute remaining stress components
c
        sig(5)=detfinv*xmu*b5
        sig(6)=detfinv*xmu*b6
c
c     material model only available for solids and shells
c
      else
c       write(iotty,20) etype
c       write(iohsp,20) etype
c       write(iomsg,20) etype
c       call adios(TC_ERROR)
        cerdat(1)=etype
        call lsmsg(3,MSG_SOL+1151,ioall,ierdat,rerdat,cerdat,0)
      endif
c
c20   format(/
c    1 ' *** Error element type ',a,' can not be',
c    2 '           run with the current material model.')
      return
      end
c	
      subroutine umat43(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt)
	  
c Orthotropic elastic material User subroutine
      DIMENSION :: cm(*),eps(*),sig(*),hsv(*)
      REAL      :: dt1,capa,tt,epsp
	  REAL      :: JAC,G1,G2,G3,G4,G5,G6
      REAL      :: E1,E2,E3,v12,v23,v31,v21,v13,v32,G12,G23,G31
      REAL*8    :: C11,C22,C33,C12,C13,C23,efs
	  REAL*8    :: C44,C55,C66,S,S1,S2,S3,S4,S5,S6
	  REAL*8   :: FS11,FS12,FS13,FS21,FS22,FS23,FS31,FS32,FS33
	  INTEGER   :: i,j
      CHARACTER*5  etype
	              
c    Assigning the values
      E1 = cm(1)  
      E2 = cm(2)  
      E3 = cm(3)  
      v21 = cm(4) 
      v31 = cm(5)
      v32 = cm(6)
      G12 = cm(7)
      G23 = cm(8)
      G31 = cm(9)            
c      
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
	  
	  if (etype .eq. 'solid') Then
	 
c       Right Green-St. Venant Strain Tensor
        G1 = 0.5*(hsv(3)*hsv(3)+hsv(4)*hsv(4)+
     .       hsv(5)*hsv(5)-1.)
        G2 = 0.5*(hsv(6)*hsv(6)+hsv(7)*hsv(7)+
     .       hsv(8)*hsv(8)-1.)
        G3 = 0.5*(hsv(9)*hsv(9)+hsv(10)*hsv(10)+
     .       hsv(11)*hsv(11)-1.)
        G4 = 0.5*(hsv(3)*hsv(6)+hsv(4)*hsv(7)+
     .       hsv(5)*hsv(8))
        G5 = 0.5*(hsv(6)*hsv(9)+hsv(7)*hsv(10)+
     .       hsv(8)*hsv(11))
        G6 = 0.5*(hsv(3)*hsv(9)+hsv(4)*hsv(10)+
     .       hsv(5)*hsv(11))

c       Second Piola-Kirchhoff Stress Tensor
        S1 = C11*eps(1) + C12*eps(2) + C13*eps(3)
        S2 = C12*eps(1) + C22*eps(2) + C23*eps(3)
        S3 = C13*eps(1) + C23*eps(2) + C33*eps(3)
        S4 = 2.*C44*eps(4)
        S5 = 2.*C55*eps(5)
        S6 = 2.*C66*eps(6)
		
c       Jacobian
        JAC=hsv(3)*(hsv(7)*hsv(11)-hsv(8)*hsv(10))-
     .       hsv(4)*(hsv(6)*hsv(11)-hsv(8)*hsv(9))+
     .       hsv(5)*(hsv(6)*hsv(10)-hsv(7)*hsv(9))
c       Cauchy Stresses (F*S*F^T)
        FS11 = hsv(3)*S1 + hsv(6)*S4 + hsv(9)*S6
        FS12 = hsv(3)*S4 + hsv(6)*S2 + hsv(9)*S5
        FS13 = hsv(3)*S6 + hsv(6)*S5 + hsv(9)*S3
        FS21 = hsv(4)*S1 + hsv(7)*S4 + hsv(10)*S6
        FS22 = hsv(4)*S4 + hsv(7)*S2 + hsv(10)*S5
        FS23 = hsv(4)*S6 + hsv(7)*S5 + hsv(10)*S3
        FS31 = hsv(5)*S1 + hsv(8)*S4 + hsv(11)*S6
        FS32 = hsv(5)*S4 + hsv(8)*S2 + hsv(11)*S5
        FS33 = hsv(5)*S6 + hsv(8)*S5 + hsv(11)*S3
c
		sig(1)=sig(1)+1./JAC*(FS11*hsv(3)+FS12*hsv(6)+FS13*hsv(9))
        sig(2)=sig(2)+1./JAC*(FS21*hsv(4)+FS22*hsv(7)+FS23*hsv(10))
        sig(3)=sig(3)+1./JAC*(FS31*hsv(5)+FS32*hsv(8)+FS33*hsv(11))
        sig(4)=sig(4)+1./JAC*(FS11*hsv(4)+FS12*hsv(7)+FS13*hsv(10))
        sig(5)=sig(5)+1./JAC*(FS21*hsv(5)+FS22*hsv(8)+FS23*hsv(11))
        sig(6)=sig(6)+1./JAC*(FS11*hsv(5)+FS12*hsv(8)+FS13*hsv(11))


      efs = (G1+G2+G3)
      hsv(1) = hsv(1)+efs
      hsv(2) = -(1.0/3.0) * (sig(1)+sig(2)+sig(3))
	  END IF
      RETURN
      END	  
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
      real*8 x(9,*)
      real*8 xx(3)
	  real*8 vx(3)
	  real*8 ax(3)
c
      xx(1)=x(1,n)
      xx(2)=x(2,n)
      xx(3)=x(3,n)
c
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
c	  write(6,*) 'vx=',vx(1)
c	  write(6,*) 'vy=',vx(2)
c	  write(6,*) 'vz=',vx(3)
      return
      end
c
      subroutine umat44 (cm,eps,sig,epsp2,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv)
!
!    Written BY  : AMITHKUMAR R ACHARI
!
!   Characteristics:
!           1. Fully Integrated Solid Element
!           2. Formulation based on Incremental Deformation Theory
!           3. Orthotropic Material Model Behaviour
!
!    Note: 1) hisv(1) & hisv(2) are always fixed in all Ls-Dyna UMAT's
!             Where:
!              hisv(1) = Equivalent (Effective) Plastic Strain , epn
!              hisv(2) = Total Hydrostatic Pressure, davg
!
!**********************************************************************
!        Definitions:
!              E    : Young's Modulus
!              G    : Shear Modulus
!              nu   : Poisson's Ratio
!              cc   : Stiffness Matix
!**********************************************************************
!
! variables
!
!     cm(1)   = Young's Modulus, E1
!     cm(2)   = Young's Modulus, E2
!     cm(3)   = Young's Modulus, E3
!     cm(4)   = Shear Modulus, G12
!     cm(5)   = Shear Modulus, G23
!     cm(6)   = Shear Modulus, G31
!     cm(7)   = Poisson's Ratio, nu12
!     cm(8)   = Poisson's Ratio, nu23
!     cm(9)   = Poisson's Ratio, nu31
!
!     eps(1)  = local x  strain (Incremental)
!     eps(2)  = local y  strain
!     eps(3)  = local z  strain
!     eps(4)  = local xy strain (Shear Strain eps(4) = 2*eps-xy)
!     eps(5)  = local yz strain
!     eps(6)  = local zx strain
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
      INTEGER :: i, j
      REAL    :: E1, E2, E3, G12, G23, G31, detF
	  REAL    :: G1, G2, G3, G4, G5, G6, S1, S2, S3, S4, S5, S6 
	  REAL    :: FS11, FS22, FS33, FS12, FS13, FS21, FS23, FS31, FS32
	  REAL    :: capa, tt, temper, dt1, epsp2, Crate1, FAIL_NUM
      REAL    :: v12,v23,v31,v21,v32,v13, delta
      REAL    :: C11,C22,C33,C12,C13,C31,C23,C32,C44,C55,C66,davg
      REAL, DIMENSION (6,6) :: cc
      REAL, DIMENSION (6)   :: PK2, ee, SRT, So, SR, SRR
	  REAL, DIMENSION (3,3) :: DG, DGT, PK2M, SIGMA
	  
!
!    Assign Variables to Common Names
!
      E1 = cm(1)
      E2 = cm(2)
      E3 = cm(3)
      v21 = cm(4)
      v31=cm(5)
      v32=cm(6)
      G12 = cm(7)
      G23 = cm(8)
      G31 = cm(9)            
!      
      v12 = v21*E1/E2
      v23 = v32*E2/E3
      v13 = v31*E1/E3
!
!     stiffness matrix 
!
      delta = (1.-v12*v21-v23*v32-v31*v13-2*v21*v32*v13)/(E1*E2*E3)
!
      C11=(1. - v23*v32)/(E2*E3*delta)
      C12=(v21 + v31*v23)/(E2*E3*delta)
      C13=(v31 + v21*v32)/(E2*E3*delta)
      C22=(1. - v31*v13)/(E1*E3*delta)
      C23=(v32 + v31*v12)/(E1*E3*delta)
      C33=(1. - v12*v21)/(E1*E2*delta)
      C44=G12
      C55=G23
      C66=G31
!
!
!       Right Green-St. Venant Strain Tensor
!
        G1 = 0.5*(hsv(1)**2+hsv(2)**2+hsv(3)**2-1)
        G2 = 0.5*(hsv(4)**2+hsv(5)**2+hsv(6)**2-1)
        G3 = 0.5*(hsv(7)**2+hsv(8)**2+hsv(9)**2-1.)
        G4 = 0.5*(hsv(1)*hsv(4)+hsv(2)*hsv(5)+hsv(3)*hsv(6))
        G5 = 0.5*(hsv(4)*hsv(7)+hsv(5)*hsv(8)+hsv(6)*hsv(9))
        G6 = 0.5*(hsv(1)*hsv(7)+hsv(2)*hsv(8)+hsv(3)*hsv(9))
!
!       Second Piola-Kirchhoff Stress Tensor
!
        S1 = C11*G1 + C12*G2 + C13*G3
        S2 = C12*G1 + C22*G2 + C23*G3
        S3 = C13*G1 + C23*G2 + C33*G3
        S4 = 2.*C44*G4
        S5 = 2.*C55*G5
        S6 = 2.*C66*G6
!
!       Jacobian
!
        detF=hsv(1)*(hsv(5)*hsv(9)-hsv(6)*hsv(8))-
     .       hsv(2)*(hsv(4)*hsv(9)-hsv(6)*hsv(7))+
     .       hsv(3)*(hsv(4)*hsv(8)-hsv(5)*hsv(7))
!
!       Cauchy Stresses (F*S*F^T)
!
        FS11 = hsv(1)*S1 + hsv(4)*S4 + hsv(7)*S6
        FS12 = hsv(1)*S4 + hsv(4)*S2 + hsv(7)*S5
        FS13 = hsv(1)*S6 + hsv(4)*S5 + hsv(7)*S3
        FS21 = hsv(2)*S1 + hsv(5)*S4 + hsv(8)*S6
        FS22 = hsv(2)*S4 + hsv(5)*S2 + hsv(8)*S5
        FS23 = hsv(2)*S6 + hsv(5)*S5 + hsv(8)*S3
        FS31 = hsv(3)*S1 + hsv(6)*S4 + hsv(9)*S6
        FS32 = hsv(3)*S4 + hsv(6)*S2 + hsv(9)*S5
        FS33 = hsv(3)*S6 + hsv(6)*S5 + hsv(9)*S3
!
        sig(1) = 1./detF*(FS11*hsv(1)+FS12*hsv(4)+FS13*hsv(7))
        sig(2) = 1./detF*(FS21*hsv(2)+FS22*hsv(5)+FS23*hsv(8))
        sig(3) = 1./detF*(FS31*hsv(3)+FS32*hsv(6)+FS33*hsv(9))
        sig(4) = 1./detF*(FS11*hsv(2)+FS12*hsv(5)+FS13*hsv(8))
        sig(5) = 1./detF*(FS21*hsv(3)+FS22*hsv(6)+FS23*hsv(9))
        sig(6) = 1./detF*(FS11*hsv(3)+FS12*hsv(6)+FS13*hsv(9))
!
!	Printing local stresses
!
		DO i=1,6
			print*, 'Sig(',i,')=',sig(i)
		ENDDO
!
      return
      end
!
!	 END OF MAT_002 UMAT
      subroutine metalsld (cm,eps,sig,epsp,hisv,dt,capa,elform,tt,temp,
     .     failel,crv,nnpcrv)
      implicit none
      include 'nlqparm'
      real cm(*),eps(*),sig(*),hisv(*)
      real crv(lq1,2,*)
      integer nnpcrv(*)
      character*5 elform
      logical failel
      real epsp,dt,capa,tt,temp
c     
c     Elastic plastic material routine for solids
c     
c     Material parameters:
c     
c     cm(1)=Young's modulus
c     cm(2)=Poisson's ratio
c     cm(3)=Bulk modulus
c     cm(4)=Shear modulus
c     cm(5)=Hardening rule
c           .eq.0.:Linear hardening 
c           .ne.0.:Hardening determined by load curve with ID cm(5)
c     cm(6)=Yield stress for linear hardening
c     cm(7)=Hardening modulus for linear hardening
c     cm(8)=Failure strain
c
c     History variables:
c
c     epsp   =Effective plastic strain
c     hisv(1)=Multiplier for plastic correction in implicit tangent
c
      real em,pr,bm,gm,hrd,epsf,g2,g3,sigyld,sighrd,dvol,davg,p,sigv,
     1     sig1d,sig2d,sig3d,sig456,sigeff,deps,coef
c     
c     Compute material parameters
c     
      em=cm(1)
      pr=cm(2)
      bm=cm(3)
      gm=cm(4)
      hrd=cm(5)
      epsf=cm(8)
c     
c     Some useful shear modulus parameters
c     
      g2=2.*gm
      g3=3.*gm
c     
c     Compute yield stress for current value of plastic strain
c     
      if (hrd.le.0.) then
         sigyld=cm(6)+cm(7)*epsp
         sighrd=cm(7)
      else
         call crvval(crv,nnpcrv,hrd,epsp,sigyld,sighrd)
      endif
c     
c     Compute volumetric strain increment
c     
      dvol=eps(1)+eps(2)+eps(3)
      davg=dvol/3.
c     
c     Compute increment in (negative) pressure
c     
      p=bm*dvol
c     
c     Compute elastic trial stress
c     
      sig(1)=sig(1)+p+g2*(eps(1)-davg)
      sig(2)=sig(2)+p+g2*(eps(2)-davg)
      sig(3)=sig(3)+p+g2*(eps(3)-davg)
      sig(4)=sig(4)+gm*eps(4)
      sig(5)=sig(5)+gm*eps(5)
      sig(6)=sig(6)+gm*eps(6)
c     
c     Compute effective stress in the sense of von Mises
c     
      sigv=(sig(1)+sig(2)+sig(3))/3.
      sig1d=sig(1)-sigv
      sig2d=sig(2)-sigv
      sig3d=sig(3)-sigv
      sig456=sig(4)**2+sig(5)**2+sig(6)**2
c     
      sigeff=sqrt(1.5*(sig1d**2+sig2d**2+sig3d**2+2.*sig456))
c
c     Indicate elastic update as default
c
      hisv(1)=0.
c     
c     Check if yielding
c     
      if (sigeff.gt.sigyld) then
c
c     Plastic update
c
         hisv(1)=.9
c     
c     Yield stress exceeded, perform elastic-plastic stress return
c     
c     Increment plastic strain
c     
         deps=(sigeff-sigyld)/(sighrd+g3)
         epsp=epsp+deps
c     
c     Return stress in the direction of deviatoric stress
c     
         coef=1.-g3*deps/sigeff
         sig(1)=sigv+sig1d*coef
         sig(2)=sigv+sig2d*coef
         sig(3)=sigv+sig3d*coef
         sig(4)=sig(4)*coef
         sig(5)=sig(5)*coef
         sig(6)=sig(6)*coef
c     
      endif
c     
c     Check if failure
c
      if (epsp.gt.epsf) failel=.true.
c     
      return
      end
      subroutine metalshl (cm,eps,sig,epsp,hisv,dt,capa,elform,tt,temp,
     .     failel,crv,nnpcrv)
      implicit none
      include 'nlqparm'
      real cm(*),eps(*),sig(*),hisv(*)
      real crv(lq1,2,*)
      integer nnpcrv(*)
      character*5 elform
      logical failel
      real epsp,dt,capa,tt,temp
c     
c     Elastic plastic material routine for shells
c     
c     Material parameters:
c     
c     cm(1)=Young's modulus
c     cm(2)=Poisson's ratio
c     cm(3)=Bulk modulus
c     cm(4)=Shear modulus
c     cm(5)=Hardening rule
c           .eq.0.:Linear hardening 
c           .ne.0.:Hardening determined by load curve with ID cm(5)
c     cm(6)=Yield stress for linear hardening
c     cm(7)=Hardening modulus for linear hardening
c     cm(8)=Failure strain
c     
c     History variables:
c
c     epsp   =Effective plastic strain
c     hisv(1)=Multiplier for plastic correction in implicit tangent
c
      real em,pr,bm,gm,hrd,epsf,g2,g3,gc,gt,sigyld,sighrd,eps12,
     1     eps3e,dvol,davg,p,sig1,sig2,sig4,sig5,sig6,sig56,sig456,
     2     sigv,sig1d,sig2d,sig3d,sigeff2,sigeff,quot,deps,coef,
     3     sig3e,eps3p,sig3p
c     
c     Compute material parameters
c     
      em=cm(1)
      pr=cm(2)
      bm=cm(3)
      gm=cm(4)
      hrd=cm(5)
      epsf=cm(8)
c     
c     Some useful shear modulus parameters
c     
      g2=2.*gm
      g3=3.*gm
      gc=capa*gm
      gt=4.5*(gc-gm)
c     
c     Compute yield stress for current value of plastic strain
c     
      if (hrd.le.0.) then
         sigyld=cm(6)+cm(7)*epsp
         sighrd=cm(7)
      else
         call crvval(crv,nnpcrv,hrd,epsp,sigyld,sighrd)
      endif
c     
c     Compute elastic through thickness strain increment
c     
      eps12=eps(1)+eps(2)
      eps3e=-pr*eps12/(1.-pr)
c     
c     Compute volumetric strain increment
c     
      dvol=eps12+eps3e
      davg=dvol/3.
c     
c     Compute increment in (negative) pressure
c     
      p=bm*dvol
c     
c     Compute elastic trial stress
c     
      sig1=sig(1)+p+g2*(eps(1)-davg)
      sig2=sig(2)+p+g2*(eps(2)-davg)
      sig4=sig(4)+gm*eps(4)
      sig5=sig(5)+gc*eps(5)
      sig6=sig(6)+gc*eps(6)
c
c     Compute some parameters of shear stress components
c     useful for further calculations
c     
      sig56=sig5**2+sig6**2
      sig456=sig4**2+sig56
c     
c     Compute effective stress in the sense of von Mises
c     
      sigv=(sig1+sig2)/3.
      sig1d=sig1-sigv
      sig2d=sig2-sigv
      sig3d=-sigv
c     
      sigeff2=1.5*(sig1d**2+sig2d**2+sig3d**2+2.*sig456)
      sigeff=sqrt(sigeff2)
c
c     Indicate elastic update as default
c
      hisv(1)=0.
c     
c     Check if yielding
c     
      if (sigeff.gt.sigyld) then
c     
c     Yield stress exceeded, perform elastic-plastic stress return
c     
c     Compute through thickness stress 
c     for elastic through thickness strain
c     
         quot=sig56/sigeff2
         deps=(sigeff-sigyld)/(sighrd+g3+gt*quot)
         coef=g3*deps/sigeff
         sig3e=-coef*sig3d
c     
c     First point in secant iteration algorithm: (eps3e,sig3e)
c     
c     Compute "plastic" through thickness strain increment
c     
         eps3p=-eps12
c     
c     Compute elastic trial stress
c     
         sig1=sig(1)+g2*eps(1)
         sig2=sig(2)+g2*eps(2)
         sig3p=g2*eps3p
c     
c     Compute effective stress in the sense of von Mises
c     
         sigv=(sig1+sig2+sig3p)/3.
         sig1d=sig1-sigv
         sig2d=sig2-sigv
         sig3d=sig3p-sigv
c     
         sigeff2=1.5*(sig1d**2+sig2d**2+sig3d**2+2.*sig456)
         sigeff=sqrt(sigeff2)
c     
c     Check if yielding
c     
         if (sigeff.gt.sigyld) then
c     
c     Yield stress exceeded, perform stress return
c     
            quot=sig56/sigeff2
            deps=(sigeff-sigyld)/(sighrd+g3+gt*quot)
            coef=g3*deps/sigeff
            sig3p=sig3p-coef*sig3d
c     
         endif
c     
c     Second point in secant iteration algorithm: (eps3p,sig3p)
c     
c     Determine through thickness strain based on 
c     interpolation from previous two stress updates
c     
         if (abs(sig3e-sig3p).gt.1.e-4*
     1        max(abs(sig3e),abs(sig3p))) then
            eps(3)=eps3e-sig3e*(eps3e-eps3p)/(sig3e-sig3p)
         else
            eps(3)=eps3p
         endif
c     
c     Compute volumetric strain increment
c     
         dvol=eps12+eps(3)
         davg=dvol/3.
c     
c     Compute increment in (negative) pressure
c     
         p=bm*dvol
c     
c     Compute elastic trial stress
c     
         sig(1)=sig(1)+p+g2*(eps(1)-davg)
         sig(2)=sig(2)+p+g2*(eps(2)-davg)
         sig(3)=p+g2*(eps(3)-davg)
c     
c     Compute effective stress in the sense of von Mises
c     
         sigv=(sig(1)+sig(2)+sig(3))/3.
         sig1d=sig(1)-sigv
         sig2d=sig(2)-sigv
         sig3d=sig(3)-sigv
c     
         sigeff2=1.5*(sig1d**2+sig2d**2+sig3d**2+2.*sig456)
         sigeff=sqrt(sigeff2)
c     
c     Check if yielding
c     
         if (sigeff.gt.sigyld) then
c     
c     Plastic update
c     
            hisv(1)=.9
c     
c     Yield stress exceeded, perform stress return
c     
c     Increment plastic strain
c     
            quot=sig56/sigeff2
            deps=(sigeff-sigyld)/(sighrd+g3+gt*quot)
            epsp=epsp+deps
c     
c     Return stress in the direction of deviatoric stress
c     and set through thickness stress to zero
c     
            coef=g3*deps/sigeff
            sig(1)=sigv+sig1d*(1.-coef)
            sig(2)=sigv+sig2d*(1.-coef)
            sig(3)=0.
            sig(4)=sig4*(1.-coef)
            sig(5)=sig5*(1.-capa*coef)
            sig(6)=sig6*(1.-capa*coef)
c     
         else
c
c     Yield stress not exceeded
c
            sig(3)=0.
            sig(4)=sig4
            sig(5)=sig5
            sig(6)=sig6
c
         endif
         
c     
c     Elastic-plastic stress update finished
c         
      else
c     
c     The update was elastic, just store the stress
c     and strain increment with no update of history
c     variables, i.e., plastic strain
c     
         sig(1)=sig1
         sig(2)=sig2
         sig(3)=0.
         sig(4)=sig4
         sig(5)=sig5
         sig(6)=sig6
c     
c     Store through thickness strain increment
c
         eps(3)=eps3e
c     
c     Elastic stress update finished
c
      endif
c
c     Check if failure
c     
      if (epsp.gt.epsf) failel=.true.
c     
      return
      end
      subroutine metaltan (cm,eps,sig,epsp,hisv,dt,capa,elform,tt,temp,
     1 failel,crv,nnpcrv,es)
      implicit none
      include 'nlqparm'
      real cm(*),eps(*),sig(*),hisv(*)
      real crv(lq1,2,*)
      integer nnpcrv(*)
      character*5 elform
      logical failel
      real epsp,dt,capa,tt,temp
      real es(6,6)
c     
c     Elastic plastic material routine associated with 
c     the metal* routines above, tangent is the same for 
c     solids and shells. For shells this is a code that
c     assumes capa=1.0, a minor detail.
c
      real em,pr,bm,gm,hrd,g3,bg23,bg43,sigv,sig1d,sig2d,sig3d,sig456,
     1     sigeff,sighrd,coef,cdfds(6)
      integer i,j
c     
c     Compute material parameters
c     
      em=cm(1)
      pr=cm(2)
      bm=cm(3)
      gm=cm(4)
      hrd=cm(5)
c     
c     Some useful modulus parameters
c     
      g3=3.*gm
      bg23=bm-2.*gm/3.
      bg43=bm+4.*gm/3.
c
c     Elastic tangent
c     
      es(1,1)=bg43
      es(2,2)=bg43
      es(3,3)=bg43
      es(2,1)=bg23
      es(3,1)=bg23
      es(3,2)=bg23
      es(1,2)=es(2,1)
      es(1,3)=es(3,1)
      es(2,3)=es(3,2)
      es(4,4)=gm
      es(5,5)=gm
      es(6,6)=gm
c     
      if (hisv(1).gt.0.) then
c     
c     Need plastic correction
c     
c     Compute effective stress in the sense of von Mises
c     
         sigv=(sig(1)+sig(2)+sig(3))/3.
         sig1d=sig(1)-sigv
         sig2d=sig(2)-sigv
         sig3d=sig(3)-sigv
         sig456=sig(4)**2+sig(5)**2+sig(6)**2
c     
         sigeff=sqrt(1.5*(sig1d**2+sig2d**2+sig3d**2+2.*sig456))
c
c     Compute hardening
c
         if (hrd.le.0.) then
            sighrd=cm(7)
         else
            call crvval(crv,nnpcrv,hrd,epsp,coef,sighrd)
         endif
c     
c     Compute C*df/ds
c     
         coef=g3/sigeff
         cdfds(1)=coef*sig1d
         cdfds(2)=coef*sig2d
         cdfds(3)=coef*sig3d
         cdfds(4)=coef*sig(4)
         cdfds(5)=coef*sig(5)
         cdfds(6)=coef*sig(6)
c
c     Update tangent
c     
         coef=hisv(1)/(g3+sighrd)
         do j=1,6
            do i=1,6
               es(i,j)=es(i,j)-coef*cdfds(i)*cdfds(j)
            enddo
         enddo
c     
      endif
c
      return
      end
c$$$      subroutine umat42 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
c$$$     1 temper,failel,crv,cma,temper)
c$$$      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
c$$$      character*5 etype
c$$$      logical failel
c$$$c
c$$$      return
c$$$      end
      subroutine umat46 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,thhsvi,nthhsv,qmat,elsiz,idele,
     2 reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),thhsvi(*),
     . cma(*),qmat(3,3)
      integer nnpcrv(*)
      character*5 etype
      logical failel,reject
      INTEGER8 idele
c
c     Isotropic elastic material with access to history variables
c     of a corresponding thermal user material. See also
c     subroutine umat41 for information about non-thermal
c     specific information. Access to thermal history variables
c     requieres IHVE.eq.1 on *MAT_THERMAL_USER_DEFINED and
c     is only possible for a selection of shell and solid elements.
c     More information including a list of supported element
c     types is given in the comments to subroutine thusrmat.
c
c     Read only variables
c
c     cm(1)=first material constant, here young's modulus
c     cm(2)=second material constant, here poisson's ratio
c     cm(3)=.eq.0 do nothing
c           .eq.1 print mechanical and thermal history variables.
c
c     nthhsv=Number of thermal history variables.
c     thhsvi=If IHVE.eq.1 on *MAT_THERMAL_USER_DEFINED then
c            thhsvi(i) i=1...nthhsv, contains history
c            variable data from thermal user material routine
c            interpolated to current integration point using
c            shape functions for the corresponding thermal element.
c
c     compute shear modulus, g
c
      if (ncycle.eq.1) then
        call usermsg('mat46')
      endif
c
      g2 =cm(1)/(1.+cm(2))
      g  =.5*g2
c
      if (etype.eq.'solid'.or.etype.eq.'shl_t'.or.
     1     etype.eq.'sld2d'.or.etype.eq.'tshel'.or.
     2     etype.eq.'sph  '.or.etype.eq.'sldax') then
        davg=(-eps(1)-eps(2)-eps(3))/3.
        p=-davg*cm(1)/(1.-2.*cm(2))
        sig(1)=sig(1)+p+g2*(eps(1)+davg)
        sig(2)=sig(2)+p+g2*(eps(2)+davg)
        sig(3)=sig(3)+p+g2*(eps(3)+davg)
        sig(4)=sig(4)+g*eps(4)
        sig(5)=sig(5)+g*eps(5)
        sig(6)=sig(6)+g*eps(6)
c
      else if (etype.eq.'shell') then
        gc    =capa*g
        q1    =cm(1)*cm(2)/((1.0+cm(2))*(1.0-2.0*cm(2)))
        q3    =1./(q1+g2)
        eps(3)=-q1*(eps(1)+eps(2))*q3
        davg  =(-eps(1)-eps(2)-eps(3))/3.
        p     =-davg*cm(1)/(1.-2.*cm(2))
        sig(1)=sig(1)+p+g2*(eps(1)+davg)
        sig(2)=sig(2)+p+g2*(eps(2)+davg)
        sig(3)=0.0
        sig(4)=sig(4)+g *eps(4)
        sig(5)=sig(5)+gc*eps(5)
        sig(6)=sig(6)+gc*eps(6)
c
      else
c       write(iotty,10) etype
c       write(iohsp,10) etype
c       write(iomsg,10) etype
c       call adios(TC_ERROR)
        cerdat(1)=etype
        call lsmsg(3,MSG_SOL+1152,ioall,ierdat,rerdat,cerdat,0)
      endif
c
c     print history variables?
c
      if ( nint(cm(3)).eq.1) then
        write(iotty,*) " epsp ",epsp," thhsvi(1) ",thhsvi(1),
     1   " thhsvi(5) ",thhsvi(5),nthhsv
      endif
c
c10   format(/
c    1 ' *** Error element type ',a,' can not be',
c    2 '           run with the current user material model #47.')
      return
      end
#ifndef AEROQA
      subroutine umat47 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      logical failel,reject
      character*5 etype
      INTEGER8 idele
c
      if (ncycle.eq.1) then
        call usermsg('mat47')
      endif
c
      return
      end
      subroutine umat48 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      logical failel,reject
      character*5 etype
      INTEGER8 idele
c
      if (ncycle.eq.1) then
        call usermsg('mat48')
      endif
c
      return
      end
#endif
      subroutine umat49 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      integer nnpcrv(*)
      logical failel,reject
      character*5 etype
      INTEGER8 idele
c
      if (ncycle.eq.1) then
        call usermsg('mat49')
      endif
c
      return
      end
      subroutine umat50 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      integer nnpcrv(*)
      logical failel,reject
      character*5 etype
      INTEGER8 idele
c
      if (ncycle.eq.1) then
        call usermsg('mat50')
      endif
c
      return
      end
      subroutine usermsg(msg)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     print warning message
c
      include 'iounits.inc'
      character*(*) msg
c
c     write(iotty,2100) msg
c     write(iohsp,2100) msg
c     write(iomsg,2100) msg
      cerdat(1)=msg
      call lsmsg(1,MSG_SOL+1147,ioall,ierdat,rerdat,cerdat,0)
c2100 format(/
c    1 ' *** Warning'
c    2 5x,'Default ',a,' is used',/,
c    3 5x,'Please check user material library')
c
      return
      end
      subroutine umat41v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,nnpcrv,cma,qmat,elsizv,idelev,
     . reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension cm(*),epsps(*),hsvs(nlq,*),dt1siz(*)
      dimension temps(*),crv(lq1,2,*),cma(*),qmat(nlq,3,3),elsizv(*)
      integer nnpcrv(*)
      INTEGER8 idelev(*)
      logical failels(*),reject
      character*5 etype
c
      dimension sig(6),eps(6),hsv(NHISVAR)
c
      do i=lft,llt
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
c
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
c
        call umat41(cm,eps,sig,epsps(i),hsv,dt1siz(i),capa,etype,tt,
     1   temps(i),failels(i),crv,nnpcrv,cma,0,elsizv(i),idelev(i),
     2   reject)
c
        sig1(i)=sig(1)
        sig2(i)=sig(2)
        sig3(i)=sig(3)
        sig4(i)=sig(4)
        sig5(i)=sig(5)
        sig6(i)=sig(6)
        d3(i)=eps(3)
      enddo
c
      return
      end
      subroutine umat45v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,nnpcrv,cma,qmat,elsizv,idelev,ih,
     . reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
C_TASKCOMMON (bk36loc)
      common/bk36loc/index
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension cm(*),epsps(*),hsvs(nlq,*),dt1siz(*)
      dimension temps(*),crv(lq1,2,*),cma(*),qmat(nlq,3,3),elsizv(*)
      integer nnpcrv(*)
      INTEGER8 idelev(*)
      logical failels(*),reject
      character*5 etype
c
      dimension sig(6),eps(6),hsv(NHISVAR),qmats(3,3)
c
      do i=lft,llt
        index=i
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
c
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
c
        hsv(1)=hsvs(i,1)
        hsv(2)=hsvs(i,2)
        hsv(3)=hsvs(i,3)
        hsv(4)=hsvs(i,4)
        hsv(5)=hsvs(i,5)
        hsv(6)=hsvs(i,6)
        hsv(7)=hsvs(i,7)
        hsv(8)=hsvs(i,8)
        hsv(9)=hsvs(i,9)
c
        if (iabs(ih).eq.3) then
           hsv(10)=hsvs(i,10)
           hsv(11)=hsvs(i,11)
           hsv(12)=hsvs(i,12)
           hsv(13)=hsvs(i,13)
           do k=1,3
              do j=1,3
                 qmats(j,k)=qmat(i,j,k)
              enddo
           enddo
        endif
c
        call umat45(cm,eps,sig,epsps(i),hsv,dt1siz(i),capa,etype,tt,
     1  temps(i),failels(i),crv,nnpcrv,cma,qmats,elsizv(i),idelev(i),ih,
     2  reject)
c
        sig1(i)=sig(1)
        sig2(i)=sig(2)
        sig3(i)=sig(3)
        sig4(i)=sig(4)
        sig5(i)=sig(5)
        sig6(i)=sig(6)
c
        hsvs(i,1)=hsv(1)
        hsvs(i,2)=hsv(2)
        hsvs(i,3)=hsv(3)
        hsvs(i,4)=hsv(4)
        hsvs(i,5)=hsv(5)
        hsvs(i,6)=hsv(6)
        hsvs(i,7)=hsv(7)
        hsvs(i,8)=hsv(8)
        hsvs(i,9)=hsv(9)
        d3(i)=eps(3)
c
        if (iabs(ih).eq.3) then
           hsvs(i,10)=hsv(10)
           hsvs(i,11)=hsv(11)
           hsvs(i,12)=hsv(12)
           hsvs(i,13)=hsv(13)
        endif
      enddo
c
      return
      end
      subroutine umat44v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,eps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,nnpcrv,cma,qmat,elsizv,idelev,
     . reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c***  isotropic plasticity with linear hardening
c
c***  updates only the deviatoric stress so that it can be used with
c     an equation of state
c
      parameter (third=1.0/3.0)
      include 'nlqparm'
c
C_TASKCOMMON (eosdloc)
      common/eosdloc/pc(nlq)
c
      dimension cm(*),d1(*),d2(*),d3(*),d4(*),d5(*),d6(*),
     & sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*),
     & eps(*),hsvs(nlqa,*),dt1siz(*),temps(*),crv(lq1,2,*),cma(*),
     & failels(*),qmat(nlq,3,3),elsizv(*)
      integer nnpcrv(*)
      INTEGER8 idelev(*)
c
      character*5 etype
      logical failels,reject
c
c***  shear modulus, initial yield stress, hardening, and pressure cut-off
      g   =cm(1)
      sy0 =cm(2)
      h   =cm(3)
      pcut=cm(4)
c
c***  plastic strain for failure
      epsfail=cm(5)
c
      ofac=1.0/(3.0*g+h)
      twog=2.0*g
c
      do i=lft,llt
c
c***    elastic deviatoric stress
        davg=third*(d1(i)+d2(i)+d3(i))
        savg=third*(sig1(i)+sig2(i)+sig3(i))
        sig1(i)=sig1(i)-savg+twog*(d1(i)-davg)
        sig2(i)=sig2(i)-savg+twog*(d2(i)-davg)
        sig3(i)=sig3(i)-savg+twog*(d3(i)-davg)
        sig4(i)=sig4(i)+g*d4(i)
        sig5(i)=sig5(i)+g*d5(i)
        sig6(i)=sig6(i)+g*d6(i)
c
c***    radial return
        aj2=sqrt(1.5*(sig1(i)**2+sig2(i)**2+sig3(i)**2)+
     &           3.0*(sig4(i)**2+sig5(i)**2+sig6(i)**2))
        sy=sy0+h*eps(i)
        eps(i)=eps(i)+ofac*max(0.0,aj2-sy)
        synew=sy0+h*eps(i)
        scale=synew/max(synew,aj2)
c
c***    scaling for radial return
        sig1(i)=scale*sig1(i)
        sig2(i)=scale*sig2(i)
        sig3(i)=scale*sig3(i)
        sig4(i)=scale*sig4(i)
        sig5(i)=scale*sig5(i)
        sig6(i)=scale*sig6(i)
c
c***    set pressure cut-off
        pc(i)=pcut
c
c***    failure due to plastic strain
        failels(i)=eps(i).gt.epsfail
c
      enddo
c
      return
      end
      subroutine umat46v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,nnpcrv,
     . cma,thhsv,nthhsv,qmat,elsizv,
     . idelev,reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension cm(*),epsps(*),hsvs(nlq,*),dt1siz(*),thhsv(nlq,*)
      dimension temps(*),crv(lq1,2,*),cma(*),qmat(nlq,3,3),elsizv(*)
      integer nnpcrv(*)
      INTEGER8 idelev(*)
      logical failels(*),reject
      character*5 etype
c
c     For more information see comments in subroutine umat46.
c
c     nthhsv=Number of thermal history variables.
c     thhsv=If IHVE.eq.1 on *MAT_THERMAL_USER_DEFINED then
c           thhsvi(i,j) j=1...nthhsv, i=1...nlq, contains history
c           variable data from thermal user material routine
c           interpolated to current integration point using
c           shape functions for the corresponding thermal element.
c
      dimension sig(6),eps(6),hsv(NHISVAR),thhsvi(100)
c
      do i=lft,llt
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
c
        do j=1,nthhsv
          thhsvi(j)=thhsv(i,j)
        enddo
c
        call umat46(cm,eps,sig,epsps(i),hsv,dt1siz(i),capa,etype,tt,
     1   temps(i),failels(i),crv,nnpcrv,cma,thhsvi,nthhsv,0,elsizv(i),
     2   idelev(i),reject)
c
        do j=1,nthhsv
          thhsv(i,j)=thhsvi(j)
        enddo
c
        sig1(i)=sig(1)
        sig2(i)=sig(2)
        sig3(i)=sig(3)
        sig4(i)=sig(4)
        sig5(i)=sig(5)
        sig6(i)=sig(6)
        d3(i)=eps(3)
      enddo
c
      return
      end
      subroutine umat47v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,nnpcrv,cma,qmat,elsizv,idelev,
     . reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension cm(*),epsps(*),hsvs(nlq,*),dt1siz(*)
      dimension temps(*),crv(lq1,2,*),cma(*),qmat(nlq,3,3),elsizv(*)
      integer nnpcrv(*)
      INTEGER8 idelev(*)
      logical failels(*),reject
      character*5 etype
c
      return
      end
c$$$      subroutine umat48v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
c$$$     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
c$$$     . etype,tt,temps,failels,nlqa,crv,cma)
c$$$      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
c$$$      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
c$$$      dimension cm(*),epsps(*),hsvs(nlq,*),dt1siz(*)
c$$$      dimension temps(*),crv(lq1,2,*),cma(*)
c$$$      logical failels(*)
c$$$      character*5 etype
c$$$c
c$$$      return
c$$$      end
c$$$      subroutine (cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
c$$$     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
c$$$     . etype,tt,temps,failels,nlqa,crv,cma)
c$$$      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
c$$$      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
c$$$      dimension cm(*),epsps(*),hsvs(nlq,*),dt1siz(*)
c$$$      dimension temps(*),crv(lq1,2,*),cma(*)
c$$$      logical failels(*)
c$$$      character*5 etype
c$$$c
c$$$      return
c$$$      end
      subroutine umat50v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,failels,nlqa,crv,nnpcrv,cma,qmat,elsizv,idelev,
     . reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension cm(*),epsps(*),hsvs(nlq,*),dt1siz(*)
      dimension temps(*),crv(lq1,2,*),cma(*),qmat(nlq,3,3),elsizv(*)
      integer nnpcrv(*)
      INTEGER8 idelev(*)
      logical failels(*),reject
      character*5 etype
c
      return
      end
      subroutine loadsetud(time,lft,llt,crv,iduls,parm,nod,nnm1)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     Input (not modifiable)
c       time : analysis time
c       x    : coordinate of node or element center
c       d    : displacement of node or lement center
c       v    : velocity of node or lement center
c       temp : temperature of node or element center
c       crv  : value of LCURV at time=time
c       iduls: id of user_loading_set
c       parm : parameters from user_loading
c       nod  : internal node numbers
c       nnm1 : offset for node block
c     Output (defined by user)
c       udl : user-defined load curve value
      include 'nlqparm'
C_TASKCOMMON (aux8loc)
      common/aux8loc/
     & x1(nlq),x2(nlq),x3(nlq),v1(nlq),
     & v2(nlq),v3(nlq),d1(nlq),d2(nlq),
     & d3(nlq),temp(nlq),udl(nlq),tmp(nlq,12)
c
      dimension parm(*),nod(*)
c
c     node numbers are given as
c     
c     for LTYPE.eq.FORCEN,BODYFN,TEMPTN
c     do i=lft,llt
c       n=i+nnm1
c       nodint=nod(n)
c       nodext=lqfinv(nodint,1)
c     enddo
c     
c     for LTYPE.eq.PRESSS
c     do i=lft,llt
c       n=i+nnm1-1
c       np=n*7
c       n1=nod(np+1)
c       n2=nod(np+2)
c       n3=nod(np+3)
c       n4=nod(np+4)
c       (only first two for 2D case)
c       (get external node numbers using lqfinv as above)
c     enddo
c
c     for LTYPE.eq.PRESSB
c     do i=lft,llt
c       n=i+nnm1
c       n1=nod(n)
c       n2=nod(nlq+n)
c       ...
c     enddo
c       
c     sample
c     if (iduls.eq.100) then
c       do i=lft,llt
cc        to define a pressure loading proportional to sin of the angle
cc        between unit vector v,=(vx,vy,vz), and element normal
cc        take 42 x 13 as normal
c         x42=x2(i)-x4(i)
c         y42=y2(i)-y4(i) 
c         z42=z2(i)-z4(i)
c         x13=x3(i)-x1(i)
c         y13=y3(i)-y1(i) 
c         z13=z3(i)-z1(i)
c         xn = y42*z13-y13*z42
c         yn = z42*x13-z13*x42
c         zn = x42*y13-x13*y42
cc        v x normal
c         vnx= vy*zn-yn*vz
c         vny= vz*xn-zn*vx
c         vnz= vx*yn-xn*vy
c         vn = sqrt(vnx*vnx+vny*vny+vnz*vnz)
cc        abs(sin)=abs(v x n )/|n|
c         rn = 1./sqrt(xn*xn+yn*yn+zn*zn)
c         sint = vn*rn
c         fact = sint
c         udl(i)=crv*fact
c       enddo
c     else if (iduls.eq.200) then
c       do i=lft,llt
c         udl(i)=....
c       enddo
c       ....
c     endif
      do i=lft,llt
       udl(i)= crv
      enddo
c
      return
      end
      subroutine loadud(fnod,dt1,time,ires,x,d,v,a,ixs,
     . numels,ixb,numelb,idrflg,tfail,isf,p,npc,fval,iob,iadd64,numelh,
     . ixh,nhex_del,nbeam_del,nshell_del,hexarray,hextim,bemarray,
     . bemtim,shlarray,shltim,parm,numnp,fnodr,dr,vr,ndof,xmst,xmsr)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     input arrays
c
c     fnod - global nodal forces
c     fnodr- global nodal moment
c     dt1 - current time step size
c     time - current problem time
c     ires - restart flag, ( 0=solution phase     )
c                          (-n=input n parameters )
c                          ( 2=restart            )
c                          ( 3=write data into dump file  )
c                          ( 4=read data from restart file)
c
c     When data is read, DUMMY arrays are passed in the call.
c     Data should be read into a local common block which is
c     written into the restart database.
c
c
c     d - displacements
c     v - velocities
c     a - acceleations
c     dr- rotational displacements
c     vr- rotational velocities
c     xmst(numnp) =reciprocal of nodal translational masses in solution phase
c     xmsr(numnp) =reciprocal of nodal rotational masses in solution phase.
c                  this array is defined if and only if ndof=6
c
c     ixs - shell element connectivities (ixs(1,*)=part ID)
c                                        (ixs(2,*)=node 1)
c                                        (ixs(3,*)=node 2)
c                                        (ixs(4,*)=node 3)
c                                        (ixs(5,*)=node 4)
c
c     ixb - beam  element connectivities (ixb(1,*)=part ID)
c                                        (ixb(2,*)=node 1)
c                                        (ixb(3,*)=node 2)
c                                        (ixb(4,*)=orientation node)
c
c     ixh - shell element connectivities (ixh(1,*)=part ID)
c                                         or
c                                        (ixh(1,*)=0 implies element deleted)
c                                        (ixh(2,*)=node 1)
c                                        (ixh(3,*)=node 2)
c                                        (ixh(4,*)=node 3)
c                                        (ixh(5,*)=node 4)
c                                        (ixh(6,*)=node 5)
c                                        (ixh(7,*)=node 6)
c                                        (ixh(8,*)=node 7)
c                                        (ixh(9,*)=node 8)
c     numnp  - number of nodal points
c     numels - number of shell elements
c     numelb - number of beam  elements
c     numelh - number of solid elements
c     isf    - shell element failure flag (1=on)
c     tfail  - shell element failure time (eq.0:okay)
c                                         (ne.0:failure time)
c
c     idrflg - nonzero if dynamic relaxation phase
c     p      - load curve data pairs (abcissa,ordinate)
c     npc    - pointer into p.  (p(npc(lc)) points to the beginning
c              of load curve ID lc.  npoints=npc(lc+1)-npc(lc)=
c              number of points in the load curve.
c     fval   - fval(lc) is the value of load curve lc at t=time
c     iob    - i/o buffer
c     ndof   - number of degrees of freedom per node in solution phase,
c            - 0 in initialization phase
c
c
c     ELEMENT DELETION FROM THE USER LOADING SUBROUTINE
c
c     to delete elements the original file needs to have a
c     *DEFINE_ELEMENT_DEATH definition for one or more elements
c     of the type to be deleted.  The deletion time can be set
c     to a value that greatly exceeds the termination time for the
c     run.  With this keyword the necessary arrays are created that
c     allows elements of the same type to deleted during the run.
c
c     nhex_del if >0 element deletion option is active for solids
c     nbem_del if >0 element deletion option is active for beams
c     nshl_del if >0 element deletion option is active for shells
c     hexarray defines time to delete solid element,
c              the value should be >time
c     hextim   we check for solid element deletion when the current
c              time is greater to or equal to hextim
c     bemarray defines time to delete beam element,
c              the value should be >time
c     bemtim   we check for beam element deletion when the current
c              time is greater to or equal to bemim
c     shlarray defines time to delete shell element,
c              the value should be >time
c     shltim   we check for shell element deletion when the current
c              time is greater to or equal to shltim
c     parm     user loading parameter if ires<0
c
      include 'iounits.inc'
      include 'bigprb.inc'
      include 'txtline.inc'
c
      parameter (NPARM=1000)
c     common/usrldv/parm(NPARM)
c
      integer*8 iadd64
      real*8 x
      real*8 d,dr
      dimension a(3,*),v(3,*),d(3,*),fnod(3,*),ixs(5,*),ixb(4,*),
     . x(3,*),tfail(*),p(*),npc(*),fval(*),iob(*),ixh(9,*),
     . hexarray(*),bemarray(*),shlarray(*),parm(*),fnodr(3,*),
     . vr(3,*),dr(3,*),xmst(*),xmsr(*)
c     character*80 txts,mssg
c
#ifdef MPP
c  MPP special variables
c
c  remove all "cmpp" for output to a single file in MPP for UMAT routine
c
cmpp  include 'mpif.h'
cmpp  integer mpp_iam
cmpp  common/mpp_usermat/mpp_iam
c
cmpp  integer ibuf,nitems
cmpp  parameter (ibuf=512)
cmpp  parameter (nitems=2)
cmpp  integer mpp_itmp(nitems,0:ibuf),mpp_icnt,MT_UMAT0,MT_UMAT1
cmpp  real    mpp_rtmp(nitems,0:ibuf)
cmpp  equivalence(mpp_itmp,mpp_rtmp)
cmpp  integer mpp_j,icnt,iready,status(MPP_STATUS_SIZE)
cmpp  integer  lqfinvf
cmpp  external lqfinvf
c
cmpp  MT_UMAT0=999100
cmpp  MT_UMAT1=999101
#endif
c
      if (ires.lt.0) then
        n=abs(ires)
        write(iohsp,1030)
        mssg='reading user loading subroutine'
        if (longs) then
          do 11 i=1,n,8
          call gttxsg (txts,lcount)
          read (txts,'(8e20.0)',err=400) (parm(j),j=i,min(i+3,n))
          write(iohsp,1040) (j,parm(j),j=i,min(i+3,n))
   11     continue
        else
          do 10 i=1,n,8
          call gttxsg (txts,lcount)
          read (txts,1020,err=400) (parm(j),j=i,min(i+7,n))
          write(iohsp,1040) (j,parm(j),j=i,min(i+7,n))
   10     continue
        endif
        write(iohsp,1050)
#ifdef MPP
cmpp    if (mpp_iam.eq.0) then
c         open output file for processor 0
cmpp      call lsopen(21,'d3umat','formatted','unknown',ierr,0,0)
cmpp    endif
#endif
        return
      endif
c
c     sample coding of gravitational force bet. two nodes
c
c     if (ndof.eq.0) return
c     grav=6.67e-11
c     do i=1,numnp
c      xmsi=1./xmst(i)
c      xi=x(1,i)
c      yi=x(2,i)
c      zi=x(3,i)
c      do j=i+1,numnp
c       xj=x(1,j)
c       yj=x(2,j)
c       zj=x(3,j)
c       xij=xj-xi
c       yij=yj-yi
c       zij=zj-zi
c       rij2i=1./max(xij*xij+yij*yij+zij*zij,1.e-16)
c       riji =sqrt(rij2i)
c       xij=xij*riji 
c       yij=yij*riji
c       zij=zij*riji
c       xmsj=1./xmst(j)
c       fij=grav*xmsi*xmsj*rij2i
c       fx=fij*xij
c       fy=fij*yij
c       fz=fij*zij 
c       fnod(1,i)=fnod(1,i)+fx
c       fnod(2,i)=fnod(2,i)+fy
c       fnod(3,i)=fnod(3,i)+fz
c       fnod(1,j)=fnod(1,i)-fx
c       fnod(2,j)=fnod(2,i)-fy
c       fnod(3,j)=fnod(3,i)-fz
c      enddo
c     enddo
c
c     if (ires.eq.3) then
c       nloddv=160
c       call wrabsf64(iob,parm,nloddv,iadd64)
c       iadd64=iadd64+nloddv
c       return
c     endif
c
c     if (ires.eq.4) then
c       nloddv=160
c       call rdabsf64(iob,parm,nloddv,iadd64,ioerr)
c       iadd64=iadd64+nloddv
#ifdef MPP
cmpp    if (mpp_iam.eq.0) then
c         open output file for processor 0
cmpp      call lsopen(21,'d3umat','formatted','unknown',ierr,0,0)
cmpp    endif
#endif
c       return
c     endif
c
c     for shells only
c
#ifdef MPP
cmpp  do 9020 mpp_j=0,numproc-1
c loop through all CPU and output necessary one-by-one
cmpp    if (mpp_iam.ne.mpp_j) goto 9015
cmpp    if (mpp_iam.ne.0) then
cmpp      call MPI_RECV(iready,1,MPP_INTEGER,0,MT_UMAT0,
cmpp .     MPP_COMM_WORLD,status,ierr)
cmpp      mpp_icnt=0
cmpp    endif
#endif
      do 20 i=1,numels
      ixs2i=ixs(2,i)
      ixs3i=ixs(3,i)
      ixs4i=ixs(4,i)
      ixs5i=ixs(5,i)
c
      xx11 =x(1,ixs2i)
      xx21 =x(2,ixs2i)
      xx31 =x(3,ixs2i)
c
      xx12 =x(1,ixs3i)
      xx22 =x(2,ixs3i)
      xx32 =x(3,ixs3i)
c
      xx13 =x(1,ixs4i)
      xx23 =x(2,ixs4i)
      xx33 =x(3,ixs4i)
c
      xx14 =x(1,ixs5i)
      xx24 =x(2,ixs5i)
      xx34 =x(3,ixs5i)
c
      fs1 =-xx11+xx12+xx13-xx14
      fs2 =-xx21+xx22+xx23-xx24
      fs3 =-xx31+xx32+xx33-xx34
c
      ft1 =-xx11-xx12+xx13+xx14
      ft2 =-xx21-xx22+xx23+xx24
      ft3 =-xx31-xx32+xx33+xx34
c
      e=fs1*fs1+fs2*fs2+fs3*fs3
      f=fs1*ft1+fs2*ft2+fs3*ft3
      g=ft1*ft1+ft2*ft2+ft3*ft3
      area=sqrt((e*g-f*f)/16.)
      tr1 =fs2*ft3-fs3*ft2
      tr2 =fs3*ft1-fs1*ft3
      tr3 =fs1*ft2-fs2*ft1
      xmg =.25/sqrt(tr1**2+tr2**2+tr3**2)
c
c     unit normal surface vector
c
c     pressure set by user here
c
      pressure=-5.e-04*fval(1)
c
      tr1 =xmg*tr1*area*pressure
      tr2 =xmg*tr2*area*pressure
      tr3 =xmg*tr3*area*pressure
c
      fnod(1,ixs2i)=fnod(1,ixs2i)-tr1
      fnod(2,ixs2i)=fnod(2,ixs2i)-tr2
      fnod(3,ixs2i)=fnod(3,ixs2i)-tr3
      fnod(1,ixs3i)=fnod(1,ixs3i)-tr1
      fnod(2,ixs3i)=fnod(2,ixs3i)-tr2
      fnod(3,ixs3i)=fnod(3,ixs3i)-tr3
      fnod(1,ixs4i)=fnod(1,ixs4i)-tr1
      fnod(2,ixs4i)=fnod(2,ixs4i)-tr2
      fnod(3,ixs4i)=fnod(3,ixs4i)-tr3
      fnod(1,ixs5i)=fnod(1,ixs5i)-tr1
      fnod(2,ixs5i)=fnod(2,ixs5i)-tr2
      fnod(3,ixs5i)=fnod(3,ixs5i)-tr3
#ifdef MPP
c output necessary data to d3umat file - for example element ID and pressure
cmpp  if (mpp_iam.eq.0) then
cmpp    write(21,'(i8,2x,e14.7)') nelmntid(i,2),pressure
cmpp  else
cmpp    mpp_icnt=mpp_icnt+1
cmpp    mpp_itmp(1,mpp_icnt)=nelmntid(i,2)
cmpp    mpp_rtmp(2,mpp_icnt)=pressure
cmpp    if (mpp_icnt.eq.ibuf) then
cmpp      mpp_itmp(1,0)=mpp_icnt
cmpp      call MPI_SEND(mpp_itmp,nitems*(ibuf+1),MPP_INTEGER,
cmpp &     0,MT_UMAT1,MPP_COMM_WORLD,status,ierr)
cmpp      mpp_icnt=0
cmpp    endif
cmpp  endif
#endif
   20 continue
c
#ifdef MPP
cmpp 9015   continue
cmpp    if (mpp_iam.eq.0) then
cmpp      if (mpp_j.ne.0) then
cmpp        iready=1
cmpp        call MPI_SEND(iready,1,MPP_INTEGER,
cmpp &       mpp_j,MT_UMAT0,MPP_COMM_WORLD,status,ierr)
c
cmpp 9011       call MPI_RECV(mpp_itmp,nitems*(ibuf+1),MPP_INTEGER,
cmpp &           mpp_j,MT_UMAT1,MPP_COMM_WORLD,status,ierr)
cmpp        icnt=mpp_itmp(1,0)
cmpp        if (icnt.eq.-1) goto 9020
cmpp        do j=1,icnt
cmpp          write(21,'(I8,2x,e14.7)') mpp_itmp(1,j),mpp_rtmp(2,j)
cmpp        enddo
cmpp        goto 9011
cmpp      endif
cmpp    elseif (mpp_iam.eq.mpp_j) then
cmpp      if (mpp_icnt.gt.0) then
cmpp        mpp_itmp(1,0)=mpp_icnt
cmpp        call MPI_SEND(mpp_itmp,nitems*(ibuf+1),MPP_INTEGER,
cmpp &       0,MT_UMAT1,MPP_COMM_WORLD,status,ierr)
cmpp      endif
c - end of message
cmpp      mpp_itmp(1,0)=-1
cmpp      call MPI_SEND(mpp_itmp,nitems*(ibuf+1),MPP_INTEGER,
cmpp &     0,MT_UMAT1,MPP_COMM_WORLD,status,ierr)
cmpp    endif
cmpp 9020 enddo
#endif
c
      return
  400 call termin (txts,mssg,lcount,1)
 1020 format(8e10.0)
 1030 format(//' u s e r   d e f i n e d   l o a d i n g',
     .           '   p a r a m e t e r s',/)
 1040 format(
     1 5x,'   parameter number  ',i4,'=', e15.8)
 1050 format(//)
      end
      subroutine ushout(iop,sig,nzgpt,rslt,ix,x,nmtcon,nnel,mx,mte)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     rslt(1:8) array of force/moment resultants to be written to
c               the LS-TAURUS database unless changed in this user
c               routine. Components 26-33 in LS-TAURUS.
c
c     rslt(9:9) current shell thickness
c
c     iop       shell formulation flag
c               1=Hughes-Liu
c               2=Belytschko-Tsay
c               3=DKT
c               4=C0
c               5=Membrane
c               6=Hughes-Liu with full integration
c               7=Corotational Hughes-Liu with full integration
c               8=Belytschko-Leviathan
c               9=Membrane with 2 x 2 integratio
c              10=Belytschko-Wong
c              11=Corotational Hughes-Liu
c
c     sig(1:6) Stress array xx, yy, zz, xy, yz, zx
c              if iop=1 or 6 the stresses are in the global
c                            system
c              else          the stresses reference the local element
c                            system
c
c     sig(7:-) Constitutive history variables
c
c     nzgpt    Number of thru thickness integration points
c
c     ix(1:4)  element connectivities
c
c     nmtcon   length of stress+history variables at an integration poitn
c
c     nnel     internal element number
c
c     mx       part ID
c
c     mte      material type
c
      real*8 x
      dimension sig(nmtcon,*),rslt(9),ix(*),x(3,*)
      return
      end
      subroutine airusr (rbu,rbv,rba,time,dt1,dt2,param,hist,itrnon,
     . rbug,rbvg,rbag,icnv)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     user subroutine to initiate the inflation of the airbag
c
c     variables
c
c     displacements are defined at time n+1 in local system
c     velocites are defined at time n+1/2 in local system
c     accelerations are defined at time n in local system
c
c         rbu(1-3) total displacements in the local xyz directions
c         rbu(3-6) total rotations about the local xyz axes
c         rbv(1-3) velocities in the local xyz directions
c         rbv(3-6) rotational velocities about the local xyz axes
c         rba(1-3) accelerations in the local xyz directions
c         rba(3-6) rotational accelerations about the local xyz axes
c         time is the current time
c         dt1 is time step size at n-1/2
c         dt2 is time step size at n+1/2
c         param is user defined input parameters
c         hist is user defined history variables
c         itrnon is a flag to turn on the airbag inflation
c         rbug,rbvg,rbag, are similar to rbu,rbv,rba but are defined
c                         globally.
c         icnv is the airbag ID
c
c     the user subroutine sets the variable itrnon to:
c
c           itrnon=0 bag is not inflated
c           itrnon=1 bag inflation begins and this subroutine is
c                    not called again
c
      include 'iounits.inc'
      dimension rbu(6),rbv(6),rba(6),param(25),hist(25),
     . rbug(6),rbvg(6),rbag(6)
c
      itrnon=0
      ra=sqrt(rba(1)**2+rba(2)**2+rba(3)**2)
      if (ra.gt.param(1)) then
        itrnon=1
        write(iotty,100) time
        write(iohsp,100) time
        write(iomsg,100) time
      endif
 100  format (' Airbag activated at time ',1pe10.3)
c
      return
      end
      subroutine usrdmp (ioflg,iobuf,jadd64,maxsiz64,lendhu)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     dump out or read files
c
c     ioflag -
c         =0 return the length of magusr to lendhu
c         =1 read user defined data
c         =2 dump user defined data
c         useful for plot files etc.
c
      include 'bk06.inc'
      parameter (magusr=1)
      common/usrcom/usr(magusr)
      integer   iobuf
      integer*8 jadd64,maxsiz64
c
      magusrr=magusr
      if (magusrr.eq.1) return
c     if (magusr.eq.1) return
c
c     return the lenght of magusr
c
      if (ioflg.eq.0) then
        lendhu=magusr
c
c     read user's data
c
      elseif (ioflg.eq.1) then
        call rdabsf64 (iobuf,usr,magusr,jadd64,ioerr)
        jadd64=jadd64+magusr
c
c     dump user's data
c
      elseif (ioflg.eq.2) then
        call wrabsf64 (iobuf,usr,magusr,jadd64)
        jadd64=jadd64+magusr
      endif
c
      return
      end
      subroutine uctrl1 (numnp,ndof,time,dt1,dt2,prtc,pltc,frci,prto,
     . plto,frco,vt,vr,at,ar,ut,ur,xmst,xmsr,irbody,rbdyn,usrhv,
     . messag,totalm,cycle,idrint,mtype,mxrb,nrba,rbcor,x,rbv,nrbn,
     . nrb,xrb,yrb,zrb,axrb,ayrb,azrb,dtx,nmmat,rba,fvalnew,fvalold,
     . fvalmid,fvalnxt)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     user subroutine for solution control and is called at the
c     beginning of time step n+1.   The time at n+
c
c
c     note:  ls-dyna3d uses an internal numbering system to
c            accomodate arbitrary node numbering.   to access
c            information for user node n, address array location m,
c            m=lqf8(n,1). to obtain user node number, n,
c            corresponding to array address m, set n=lqfinv8(m,1)
c
c     arguments:
c          numnp=number of nodal points
c          ndof =number of degrees of freedom per node
c          time =current solution time at n+1
c          dt1  =time step size between time n-1 and n
c          dt2  =time step size between time n and n+1
c          prtc =output interval for taurus time history data
c          pltc =output interval for taurus state data
c          frci =output interval for taurus interface force data
c          prto =output time for time history file
c          plto =output time for state data
c          frco =output time for force data
c          vt(3,numnp) =nodal translational velocity vector
c          vr(3,numnp) =nodal rotational velocity vector.  this
c                       array is defined if and only if ndof=6
c          at(3,numnp) =nodal translational acceleration vector
c          ar(3,numnp) =nodal rotational acceleration vector.  this
c                       array is defined if and only if ndof=6
c          ut(3,numnp) =nodal translational displacement vector
c          ur(3,numnp) =nodal rotational displacement vector.  this
c                       array is defined if and only if ndof=6
c          xmst(numnp) =reciprocal of nodal translational masses
c          xmsr(numnp) =reciprocal of nodal rotational masses.  this
c                       array is defined if and only if ndof=6
c          irbody      =0 if no rigid bodies
c          rbdyn(numnp)=flag for rigid body nodal points
c                       if deformable node then set to 1.0
c                       if rigid body node then set to 0.0
c                       defined if and only if rigid bodies are present
c                       i.e., irbody.ne.0 if no rigid bodies are
c                       present
c          usrhv(lenhv)=user defined history variables that are stored
c                       in the restart file.  lenhv=100+7*nummat where
c                       nummat is the # of materials in the problem.
c                       array usrhv is updated only in this subroutine.
c          messag      =flag for dyna3d which may be set to:
c                      'sw1.' ls-dyna3d terminates with restart file
c                      'sw3.' ls-dyna3d writes a restart file
c                     'sw4.' ls-dyna3d writes a plot state
c           totalm     =total mass in problem
c           cycle     =cycle number
c           idrint    =flag for dynamic relaxation phase
c                      .ne.0: dynamic relaxation in progress
c                      .eq.0: solution phase
c          mtype(*)   = material type for each part in the model
c          mxrb(*)    = master rigid body for each rigid body
c                       If part n is rigid and is a master or has not
c                       been merged then mxrb(n)==n.
c          nrba(*)    = starting index in nrb(*) of the nodes for each
c                       rigid body
c          rbcor(3,*) = rigid body cg coordinates
c          x(3,*)     = Node coordinate array
c          rbv(6,*)   = rigid body cg velocity
c          nrbn(*)    = # nodes in each rigid body
c          nrb(*)     = List of all rigid body nodes
c          xrb(*)     = RB scratch array as long as longest nrbn() value
c          yrb(*)     = RB scratch array as long as longest nrbn() value
c          zrb(*)     = RB scratch array as long as longest nrbn() value
c          axrb(*)    = RB scratch array as long as longest nrbn() value
c          ayrb(*)    = RB scratch array as long as longest nrbn() value
c          azrb(*)    = RB scratch array as long as longest nrbn() value
c          dtx        = (dt1+dt2)*0.5 except at time 0 when it is = dt2
c          nmmat      = number of parts in the model
c          rba(6,*)   = rigid body cg acceleration
c          fvalnew     =array for setting load curve values at time n+1
c                       only load curves with 0 input points may be user
c                       defined.  When the load curve is user set, the
c                       value at time n must be stored in array fvalold.
c          fvalold     =array for storing load curve values at time n
c          fvalmid     =array for predicting load curve values at time n+3/2
c          fvalnxt     =array for predicting load curve values at time n+2
c                       for some applications it is necessary to predict the
c                       load curve values at time n+2, a time that is not known,
c                       this for instance for boundary prescribed motion. In
c                       this case the load curve values at time n+3/2 need to
c                       be predicted in fvalmid, and fvalold should be set to
c                       fvalnew and fvalnew should be set to fvalnxt. See
c                       coding below.
c
c
      include 'ptimes.inc'
c
c     prtims(1-37)=output intervals for ascii files
c
c      ascii files:
c           ( 1)-cross section forces
c           ( 2)-rigid wall forces
c           ( 3)-nodal data
c           ( 4)-element data
c           ( 5)-global data
c           ( 6)-discrete elements
c           ( 7)-material energies
c           ( 8)-nodal interface forces
c           ( 9)-resultant interface forces
c           (10)-smug animator
c           (11)-spc reaction forces
c           (12)-nodal constraint resultant forces
c           (13)-airbag statistics
c           (14)-avs database
c           (15)-nodal force groups
c           (16)-output intervals for nodal boundary conditions
c           (17)-(32) unused at this time
c           (37)-auto tiebreak damage output
c
c     prtlst(32)=output times for ascii files above.  when solution time
c                exceeds the output time a print state is dumped.
c
      common/rbkeng/enrbdy,rbdyx,rbdyy,rbdyz
c
c      total rigid body energies and momentums:
c           enrbdy=rigid body kinetic energy
c           rbdyx =rigid body x-momentum
c           rbdyy =rigid body y-momentum
c           rbdyz =rigid body z-momentum
c
      common/swmke/swxmom,swymom,swzmom,swkeng
c
c      total stonewall energies and momentums:
c           swxmom=stonewall x-momentum
c           swymom=stonewall y-momentum
c           swzmom=stonewall z-momentum
c           swkeng=stonewall kinetic energy
c
      common/deengs/deeng
c
c      deeng=total discrete element energy
c
      common/bk03/ndthl,nsthl,nstsl,nstbl,nsttl,endtim
c
c      endtim = termination time
c
      common/bk28/summss,xke,xpe,tt,xte0,erodeke,erodeie,selie,selke,
     . erodehg
c
c      xpe  =total internal energy in the finite elements
c
      common/sprengs/spreng
c
c      spreng=total spr energy
c
      character*(*) messag
      integer cycle
      real*8 x
      dimension vt(3,*),vr(3,*),at(3,*),ar(3,*),
     . xmst(*),xmsr(*),rbdyn(*),usrhv(*),mtype(*),mxrb(*),nrba(*),
     . rbcor(3,1),x(*),rbv(6,*),nrbn(*),nrb(*),xrb(*),yrb(*),
     . zrb(*),axrb(*),ayrb(*),azrb(*),rba(6,*),fvalnew(*),
     . fvalold(*),fvalmid(*),fvalnxt(*)
      real*8 ut(3,*),ur(3,*)
c
c     MPP special variables
c
c     remove all "cmpp" for output to a single file in MPP for UMAT routine
c
#ifdef MPP
cMPP
cMPP  user lqfe to check node whether on the CPU (>0) or not (-1)
cMPP
cMPP  include 'mpif.h'
      integer mpp_iam
      common/mpp_usermat/mpp_iam
      character*80 fname
cMPP
cMPP  temp arrays for data communication
      dimension rjunk(20),ijunk(20)
      equivalence (rjunk,ijunk)
#endif
c
c     sample momentum and kinetic energy calculations
c
c     remove all comments in column 1 below to activate
      i=1
      if (i.eq.1) return
c     return
cc
cc
cc     initialize kinetic energy, xke, and x,y,z momentums.
cc
c      xke=2.*swkeng+2.*enrbdy
c      xm=swxmom+rbdyx
c      ym=swymom+rbdyy
c      zm=swzmom+rbdyz
cc
c      numnp2=numnp
c      if (ndof.eq.6) then
c        numnp2=numnp+numnp
c      endif
c      write(iotty,*)ndof
cc
cc
cc     no rigid bodies present
cc
c      if (irbody.eq.0) then
cc       note in blank comment vr follows vt.  this fact is used below.
c        do 10 n=1,numnp2
c        xmsn=1./xmst(n)
c        vn1=vt(1,n)
c        vn2=vt(2,n)
c        vn3=vt(3,n)
c        xm=xm+xmsn*vn1
c        ym=ym+xmsn*vn2
c        zm=zm+xmsn*vn3
c        xke=xke+xmsn*(vn1*vn1+vn2*vn2+vn3*vn3)
c   10   continue
cc
cc
cc     rigid bodies present
cc
c      else
cc       nodal accerations for rigid bodies
cc
c        do 12 n=1,nmmat
c        if (mtype(n).ne.20.or.mxrb(n).ne.n) go to 12
c        lrbn=nrba(n)
c        call stvlut(rbcor(1,n),x,vt,at,ar,vr,rbv(1,n),dt2,
c    .   nrbn(n),nrb(lrbn),xrb,yrb,zrb,axrb,ayrb,azrb,dtx)
c
c        rigid body nodal accelerations
c
c        if (ndof.eq.6) then
c          call rbnacc(nrbn(n),nrb(lrbn),rba(4,n),ar)
c        endif
c
c  12    continue
cc
c        do 20 n=1,numnp
c        xmsn=1./xmst(n)
c        vn1=rbdyn(n)*vt(1,n)
c        vn2=rbdyn(n)*vt(2,n)
c        vn3=rbdyn(n)*vt(3,n)
c        xm=xm+xmsn*vn1
c        ym=ym+xmsn*vn2
c        zm=zm+xmsn*vn3
c        xke=xke+xmsn*(vn1*vn1+vn2*vn2+vn3*vn3)
c   20   continue
c        if (ndof.eq.6) then
c          do 30 n=1,numnp
c          xmsn=1./xmsr(n)
c          vn1=rbdyn(n)*vr(1,n)
c          vn2=rbdyn(n)*vr(2,n)
c          vn3=rbdyn(n)*vr(3,n)
c          xm=xm+xmsn*vn1
c          ym=ym+xmsn*vn2
c          zm=zm+xmsn*vn3
c          xke=xke+xmsn*(vn1*vn1+vn2*vn2+vn3*vn3)
c   30     continue
c        endif
c      endif
cc
cc     total kinetic energy
c      xke=.5*xke
cc     total internal energy
c      xie=xpe+deeng+spreng
cc     total energy
c      xte=xke+xpe+deeng+spreng
cc     total x-rigid body velocity
c      xrbv=xm/totalm
cc     total y-rigid body velocity
c     yrbv=ym/totalm
cc     total z-rigid body velocity
c      zrbv=zm/totalm
#ifdef MPP
cMPP
cMPP  operation MPP_MAX, MPP_MIN, MPP_SUM
cMPP  data type MPP_REAL, MPP_INTEGER
c
c     ioper=MPP_SUM
c     idtyp=MPP_REAL
c
c All values in common blocks are only for each processor.  If you would
c like to monitor a GLOBAL value, you need to sum across processors
c
c     xke=.5*xke
c     xie=xpe+deeng+spreng
c     xte=xke+xpe+deeng+spreng
c     xrbv=xm/totalm
c     yrbv=ym/totalm
c     zrbv=zm/totalm
c
c     nitem=6
c     rjunk(1) = xke
c     rjunk(2) = xie
c     rjunk(3) = xte
c     rjunk(4) = xrbv
c     rjunk(5) = yrbv
c     rjunk(6) = zrbv
c
c     call MPI_ALLREDUCE(rjunk(1),rjunk(11),nitem,idtyp,ioper,
c    . MPP_COMM_WORLD,ierr)
cMPP
cMPP use rjunk(11-16) for sense switch control and DO NOT overwritten
cMPP common block data
cMPP
#endif
c
c     New source to control load curve 2 in model during solution
c
c     External DYNA model load curve id
      lidext=2
c     Obtain the internal load curve id
      lidint=lqfl(lidext)
c     Put old 'new' value into fvalold
      fvalold(lidint)=fvalnew(lidint)
      fvalnew(lidint)=fvalnxt(lidint)
c     Calculate value for fvalnxt
      fvalnxt(lidint)=100.0*sin(time*3.1415926)
      fvalmid(lidint)=.5*(fvalnxt(lidint)+fvalnew(lidint))
c
c     New source to check that we can correctly obtain the nodal
c     displacements and velocities for an arbitrary node.
c     Using node 3 in this case on part A of model.
c     Obtain the displacements and velocities at t=0.3s.
c
      if (abs(time-0.3).lt.1e-6) then
c        Time is 0.3s (time step is 2.5e-6s)
#ifdef MPP
cMPP
cmpp  filename needs to be unique for each processor to avoid conflict
c        write (fname,'("node3.out_",i4.4)') mpp_iam
         call appendid('node3.out_',fname,mpp_iam)
         call lsopen(85,fname,'formatted','new',ierr,0,0)
         nodext=3
         nodint=lqfe(nodext,1)
         if (nodint.gt.0) then
#else
         call lsopen(85,'node3.out','formatted','new',ierr,0,0)
c        External DYNA model node number
         nodext=3
c        Obtain the internal node id
         nodint=lqf8(nodext,1)
#endif
         write (85,'('' Time = '',E10.4E2,/)')time
c        Obtain and write out the nodal displacements
         write (85,'('' Node 3 x disp. = '',E10.4E2)')ut(1,nodint)
         write (85,'('' Node 3 y disp. = '',E10.4E2)')ut(2,nodint)
         write (85,'('' Node 3 z disp. = '',E10.4E2)')ut(3,nodint)
c        Obtain and write out the nodal velocities
         write (85,'('' Node 3 x velo. = '',E10.4E2)')vt(1,nodint)
         write (85,'('' Node 3 y velo. = '',E10.4E2)')vt(2,nodint)
         write (85,'('' Node 3 z velo. = '',E10.4E2)')vt(3,nodint)
         close (unit=85)
#ifdef MPP
cMPP
        endif
#endif
      endif
c
      return
      end
      subroutine rbnacc(nrbn,nrb,rba,ar)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      dimension nrb(*),rba(6),ar(3,*)
c
      do 20 i=1,nrbn
      ar(1,nrb(i))=rba(1)
      ar(2,nrb(i))=rba(2)
      ar(3,nrb(i))=rba(3)
   20 continue
c
      return
      end
      subroutine uctrl2(nsi,nty,time,cycle,msr,nmn,nsv,nsn,
     1 thmr,thsv,vt,xi,ut,iskip,idrint,numnp,dt2,ninput,ua,
     2 irectm,nrtm,irects,nrts)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     user subroutine for interface control
c
c     note:  ls-dyna3d uses an internal numbering system to
c            accomodate arbitrary node numbering.   to access
c            information for user node n, address array location m,
c            m=lqf8(n,1). to obtain user node number, n,
c            corresponding to array address m, set n=lqfinv8(m,1)
c
c     arguments:
c          nsi        =number of sliding interface
c          nty        =interface type.
c                      .eq.4:single surface
c                      .ne.4:surface to surface
c          time       =current solution time
c          cycle      =cycle number
c          msr(nmn)   =list of master nodes numbers in internal
c                      numbering scheme
c          nmn        =number of master nodes
c          nsv(nsn)   =list of slave nodes numbers in internal
c                      numbering scheme
c          nsn        =number of slave nodes
c          thmr(nmn)  =master node thickness
c          thsv(nsn)  =slave node thickness
c          vt(3,numnp)=nodal translational velocity vector
c          xi(3,numnp)=initial coordinates at time=0
c          ut(3,numnp)=nodal translational displacement vector
c          idrint     =flag for dynamic relaxation phase
c                      .ne.0: dynamic relaxation in progress
c                      .eq.0: solution phase
c          numnp      =number of nodal points
c          dt2        =time step size at n+1/2
c          ninput     =number of variables input into ua
c          ua(*)      =users' array, first ninput locations
c                      defined by user.  the length of this
c                      array is defined on control card 10.
c                      this array is unique to interface nsi.
c          irectm(4,*)=list of master segments in internal
c                      numbering scheme
c          nrtm       =number of master segments
c          irects(4,*)=list of slave segments in internal
c                      numbering scheme
c          nrts       =number of master segments
c
c     set flag for active contact
c       iskip=0 active
c       iskip=1 inactive
c
c*******************************************************************
c
      integer cycle
      real*8 ut
      real*8 xi
      dimension msr(*),nsv(*),thmr(*),thsv(*),vt(3,*),xi(3,*),
     .          ut(3,*),ua(*),irectm(4,*),irects(4,*)
c
c     the following sample of codeing is provided to illustrate how
c     this subroutine might be used.   here we check to see if the
c     surfaces in the surface to surface contact are separated.   if
c     so the iskip=1 and the contact treatment is skipped.
c
c     if (nty.eq.4) return
c     dt2hlf=dt2/2.
c     xmins= 1.e+16
c     xmaxs=-xmins
c     ymins= 1.e+16
c     ymaxs=-ymins
c     zmins= 1.e+16
c     zmaxs=-zmins
c     xminm= 1.e+16
c     xmaxm=-xminm
c     yminm= 1.e+16
c     ymaxm=-yminm
c     zminm= 1.e+16
c     zmaxm=-zminm
c     thks=0.0
c     thkm=0.0
c     do 10 i=1,nsn
c     dsp1=ut(1,nsv(i))+dt2hlf*vt(1,nsv(i))
c     dsp2=ut(2,nsv(i))+dt2hlf*vt(2,nsv(i))
c     dsp3=ut(3,nsv(i))+dt2hlf*vt(3,nsv(i))
c     x1=xi(1,nsv(i))+dsp1
c     x2=xi(2,nsv(i))+dsp2
c     x3=xi(3,nsv(i))+dsp3
c     thks = max(thsv(i),thks)
c     xmins=min(xmins,x1)
c     xmaxs=max(xmaxs,x1)
c     ymins=min(ymins,x2)
c     ymaxs=max(ymaxs,x2)
c     zmins=min(zmins,x3)
c     zmaxs=max(zmaxs,x3)
c  10 continue
c     do 20 i=1,nmn
c     dsp1=ut(1,msr(i))+dt2hlf*vt(1,msr(i))
c     dsp2=ut(2,msr(i))+dt2hlf*vt(2,msr(i))
c     dsp3=ut(3,msr(i))+dt2hlf*vt(3,msr(i))
c     x1=xi(1,msr(i))+dsp1
c     x2=xi(2,msr(i))+dsp2
c     x3=xi(3,msr(i))+dsp3
c     thkm = max(thmr(i),thks)
c     xminm=min(xminm,x1)
c     xmaxm=max(xmaxm,x1)
c     yminm=min(yminm,x2)
c     ymaxm=max(ymaxm,x2)
c     zminm=min(zminm,x3)
c     zmaxm=max(zmaxm,x3)
c  20 continue
c
c     if thks or thkm equal zero set them to some reasonable value
c
c     if (thks.eq.0.0) then
c       e1=(xi(1,irects(1,1))-xi(1,irects(3,1)))**2
c    .    +(xi(2,irects(1,1))-xi(2,irects(3,1)))**2
c    .    +(xi(3,irects(1,1))-xi(3,irects(3,1)))**2
c       e2=(xi(1,irects(2,1))-xi(1,irects(4,1)))**2
c    .    +(xi(2,irects(2,1))-xi(2,irects(4,1)))**2
c    .    +(xi(3,irects(2,1))-xi(3,irects(4,1)))**2
c       thks=.3*sqrt(max(e1,e2))
c     endif
c     if (thkm.eq.0.0) then
c       e1=(xi(1,irectm(1,1))-xi(1,irectm(3,1)))**2
c    .    +(xi(2,irectm(1,1))-xi(2,irectm(3,1)))**2
c    .    +(xi(3,irectm(1,1))-xi(3,irectm(3,1)))**2
c       e2=(xi(1,irectm(2,1))-xi(1,irectm(4,1)))**2
c    .    +(xi(2,irectm(2,1))-xi(2,irectm(4,1)))**2
c    .    +(xi(3,irectm(2,1))-xi(3,irectm(4,1)))**2
c       thkm=.3*sqrt(max(e1,e2))
c     endif
c
c     if (xmaxs+thks.lt.xminm-thkm) go to 40
c     if (ymaxs+thks.lt.yminm-thkm) go to 40
c     if (zmaxs+thks.lt.zminm-thkm) go to 40
c     if (xmaxm+thkm.lt.xmins-thks) go to 40
c     if (ymaxm+thkm.lt.ymins-thks) go to 40
c     if (zmaxm+thkm.lt.zmins-thks) go to 40
c     iskip=0
c
c     return
c  40 iskip=1
c
      return
      end
      subroutine matusr_24(lft,llt,ietype,cm,crv,nnpcrv,
     . ipt,ipt_thk,nip,
     . nelbwp,lfail)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     user-defined damage subroutine for materials:
c     24, 36, 114, 123, 124, 133, 155, 182, 225, 238, 243, 251, 255
c
c     Description of Input Arguments:
c
c     sig1 - stress in x  direction, local x for shells
c     sig2 - stress in y  direction, local y for shells
c     sig3 - stress in z  direction, local z for shells
c     sig4 - stress in xy direction, local for shells
c     sig5 - stress in yz direction, local for shells
c     sig6 - stress in zx direction, local for shells
c
c     dt1siz-element time step size
c
c     incremental strains are used for shell elements
c     strain rates are used for solid elements
c
c     d1 - strain rate/increment in x  direction, local x for shells
c     d2 - strain rate/increment in y  direction, local y for shells
c     d3 - strain rate/increment in z  direction, local z for shells
c     d4 - strain rate/increment in xy direction, local for shells
c     d5 - strain rate/increment in yz direction, local for shells
c     d6 - strain rate/increment in zx direction, local for shells
c
c     wxxdt - time-step*spin rate in x  direction, only for solds       
c     wyydt - time-step*spin rate in y  direction, only for solds       
c     wzzdt - time-step*spin rate in z  direction, only for solds       
c
c     The following strain components are only available when
c     INTOUT in *DATABASE_EXTENT_BINARY is set to "STRAIN" or "ALL"
c
c     eps1 - strain in x  direction, local x for shells
c     eps2 - strain in y  direction, local y for shells
c     eps3 - strain in z  direction, local z for shells
c     eps4 - strain in xy direction, local for shells
c     eps5 - strain in yz direction, local for shells
c     eps6 - strain in zx direction, local for shells
c
c     eint - strain energy density increment
c     ep - plastic strain
c     fail - failure flag SET TO ZERO IF ELEMENT FAILS
c     n - number of elements in vector block
c     ietype - 1 - beams, 2 - shells, 3 - solids.
c
c     lft: 1st element to be processed
c     llt: last element to be processed
c
c     cm: material properties
c       -1: young's modulus
c       -2: poissons's ratio
c
c     crv: curve array
c
c     nnpcrv=# of discretization points per crv()
c
c     ipt: integration point number
c     ipt_thk: through-thickness integration point number
c     nip: total number of integration points
c
c     nelbwp: element id array
c
c     - solids: idele=lqfinv8(nelbwp(i),2)
c     - beams : idele=lqfinv8(nelbwp(i),3)
c     - shells: idele=lqfinv8(nelbwp(i),4)
c
c     nintcy: current cycle
c     tt: current time
c
c     lfail: where failure flag is stored, defined by ls-dyna,
c           =0: failure is stored in fail0
c           >2: failure is stored in fail1
c
c     Be aware that in addition to failure criteria set forth in this routine,
c     failure also occurs if ep >= abs(FAIL) where FAIL is the 'failure flag'
c     provided in the mat24 input.
c
      include 'nlqparm'
C_TASKCOMMON (aux2loc)
C_TASKCOMMON (aux13loc)
C_TASKCOMMON (aux14loc)
C_TASKCOMMON (failuloc)
C_TASKCOMMON (strnintloc)
C_TASKCOMMON (subtssloc)
      common/aux2loc/d1(nlq),d2(nlq),d3(nlq),d4(nlq),d5(nlq),d6(nlq),
     1 wzzdt(nlq),wyydt(nlq),wxxdt(nlq),einc(nlq)
      common/aux13loc/zeta(nlq),thick(nlq)
      common/aux14loc/
     1 sig1(nlq),sig2(nlq),sig3(nlq),sig4(nlq),
     2 sig5(nlq),sig6(nlq),  ep(nlq),fail1(nlq,142)
      common/failuloc/tmp1(nlq),fail0(nlq)
      common/subtssloc/dt1siz(nlq)
      common/strnintloc/eps1(nlq),eps2(nlq),eps3(nlq),eps4(nlq),
     1 eps5(nlq),eps6(nlq)
      common/bk26/begtim,nintcy
      common/bk28/summss,xke,xpe,tt
      dimension cm(*),crv(lq1,2,*),nelbwp(*)
      integer nnpcrv(*)
c
c     loop over elements
c     if (lfail.eq.0) then
c       do i=lft,llt
c        if (failure condition is met) fail0(i)=0.0
c       enddo
c     else
c       do i=lft,llt
c        if (failure condition is met) fail1(i,lfail)=0.0
c       enddo
c     endif
      return
      end
      subroutine matusr_103 (sig1,sig2,sig3,sig4,sig5,sig6,dt1siz,
     . d1,d2,d3,d4,d5,d6,ep,fail,n,ietype,q1,q2,d1m,d2m,d3m,d4m,d5m,
     . d6m,failure_flag,s11,s12,s13,s21,s22,s23,s31,s32,s33,cm,crv,
     . nnpcrv,
     . parms,nparms,dmg,aux,naux,fracs,nfracs,clen,cang,dep)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     description of Input Arguments, global for solid element
c
c     sig1 - stress in x  direction, local x for shells
c     sig2 - stress in y  direction, local y for shells
c     sig3 - stress in z  direction, local z for shells
c     sig4 - stress in xy direction, local for shells
c     sig5 - stress in yz direction, local for shells
c     sig6 - stress in zx direction, local for shells
c
c     dt1siz-element time step size
c
c     incremental strains are used for shell elements
c     strain rates are used for solid elements
c
c     d1  - strain rate/increment in x  direction, local x for shells
c     d2  - strain rate/increment in y  direction, local y for shells
c     d3  - strain rate/increment in z  direction, local z for shells
c     d4  - strain rate/increment in xy direction, local for shells
c     d5  - strain rate/increment in yz direction, local for shells
c     d6  - strain rate/increment in zx direction, local for shells
c
c     d1-d6 are in the a,b,c material local system for solids
c
c     q1,q2  - direction cosines to obtain material direction (shells)
c
c     To obtain stresses in shell a,b coordinate system:
c
c           a11    =q1(i)*sig1(i)-q2(i)*sig4(i)
c           a12    =q2(i)*sig1(i)+q1(i)*sig4(i)
c           a21    =q1(i)*sig4(i)-q2(i)*sig2(i)
c           a22    =q2(i)*sig4(i)+q1(i)*sig2(i)
c           sig1(i)=q1(i)*a11-q2(i)*a21
c           sig2(i)=q2(i)*a12+q1(i)*a22
c           sig4(i)=q1(i)*a12-q2(i)*a22
c
c     To transform back to shell local
c
c           a11    = sig1(i)*q1(i)+sig4(i)*q2(i)
c           a12    =-sig1(i)*q2(i)+sig4(i)*q1(i)
c           a21    = sig4(i)*q1(i)+sig2(i)*q2(i)
c           a22    =-sig4(i)*q2(i)+sig2(i)*q1(i)
c           sig1(i)= q1(i)*a11+q2(i)*a21
c           sig2(i)=-q2(i)*a12+q1(i)*a22
c           sig4(i)= q1(i)*a12+q2(i)*a22
c
c
c     d1m - strain rate/increment in a  direction (shells)
c     d2m - strain rate/increment in b  direction (shells)
c     d3m - strain rate/increment in c  direction (shells)
c     d4m - strain rate/increment in ab direction (shells)
c     d5m - strain rate/increment in bc direction (shells)
c     d6m - strain rate/increment in ca direction (shells)
c
c     ep - plastic strain
c     fail - failure flag SET TO ZERO IF ELEMENT FAILS
c     n - number of elements in vector block
c     ietype  beams, 2 - shells, 3 - solids.
c
c     Be aware that in addition to failure criteria set forth in this routine,
c     failure also occurs if ep >= abs(FAIL) where FAIL is the 'failure flag'
c     provided in the mat24 input.
c
c     s11-s33 transformation matrix for SOLIDS only.  To transform the stresses
c     into the local system:
c
c     a11    =sig1(i)*s11(i)+sig4(i)*s12(i)+sig6(i)*s13(i)
c     a12    =sig1(i)*s21(i)+sig4(i)*s22(i)+sig6(i)*s23(i)
c     a13    =sig1(i)*s31(i)+sig4(i)*s32(i)+sig6(i)*s33(i)
c     a21    =sig4(i)*s11(i)+sig2(i)*s12(i)+sig5(i)*s13(i)
c     a22    =sig4(i)*s21(i)+sig2(i)*s22(i)+sig5(i)*s23(i)
c     a23    =sig4(i)*s31(i)+sig2(i)*s32(i)+sig5(i)*s33(i)
c     a31    =sig6(i)*s11(i)+sig5(i)*s12(i)+sig3(i)*s13(i)
c     a32    =sig6(i)*s21(i)+sig5(i)*s22(i)+sig3(i)*s23(i)
c     a33    =sig6(i)*s31(i)+sig5(i)*s32(i)+sig3(i)*s33(i)
c     sig1(i)=s11(i)*a11+s12(i)*a21+s13(i)*a31
c     sig2(i)=s21(i)*a12+s22(i)*a22+s23(i)*a32
c     sig3(i)=s31(i)*a13+s32(i)*a23+s33(i)*a33
c     sig4(i)=s11(i)*a12+s12(i)*a22+s13(i)*a32
c     sig5(i)=s21(i)*a13+s22(i)*a23+s23(i)*a33
c     sig6(i)=s11(i)*a13+s12(i)*a23+s13(i)*a33
c
c     warning: For shells the stresses must be in the shell local system
c              when returning from this subroutine.  For solids, the
c              stresses must be in the global system when returning.
c
c     Additional arguments for supporting user defined damage model
c
c     parms(nparms)      = damage parameters, input
c     dmg(i)             = damage, input/output
c     aux(i,naux)        = history variable, input/output
c     fracs(i,nfracs)    = fraction material of phases, input
c     clen(i)            = characteristic element length divided by shell
c                          thickness (for solids this is zero), input
c     cang(i)            = cosine of angle between loading and first material
c                          direction, input
c     dep(i)             = plastic strain rate, input
c
      include 'nlqparm'
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*),
     . dt1siz(*),d1(*),d2(*),d3(*),d4(*),d5(*),d6(*),ep(*),
     . fail(*),q1(*),q2(*),d1m(*),d2m(*),d3m(*),d4m(*),d5m(*),d6m(*),
     . s11(*),s12(*),s13(*),s21(*),s22(*),s23(*),s31(*),s32(*),s33(*),
     . cm(*),crv(lq1,2,*)
      integer nnpcrv(*)
      dimension parms(*),dmg(*),aux(nlq,*),fracs(nlq,*),
     .     clen(*),cang(*),dep(*)
c
c     loop over elements, Put new failure criteria here.
c
       do i=1,n
c      if (failure condition is met) fail(i)=0.0
c      example
       if (ep(i).gt.abs(failure_flag)) fail(i)=0.
c
       enddo
      return
      end
#ifdef MPP
      subroutine usrfrc(fstt,fdyn,uc,nc,prs,temp,v,uh,nh,crv,nnpcrv,
     . nosl,
     . ictype,side,time,ncycle,dt2,fric1,fric2,fric3,fric4,lsv,idele8,
     . sfac1,sfac2,insv,fni,areas,stfk,ix1,ix2,ix3,ix4,aream,
     . rn1,rn2,rn3,ue,ne)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     user subroutine for interface friction control (mpp)
c
c     Purpose: To define static and dynamic frictional
c              coefficients for a slave node
c
c     Variables:
c
c     fstt            = static frictional coefficient (output)
c     fdyn            = dynamic frictional coefficient (output)
c     uc(nc)          = user defined friction parameters (input)
c     nc              = number of user defined friction parameters (input)
c     prs             = interface pressure (input)
c     temp            = temperature (input, available for coupled analysis)
c                       ... not yet supported for SOFT=4
c     v               = magnitude of relative velocity (input)
c     uh(nh)          = user defined friction history variables (input/output)
c     nh              = number of user defined friction history variables (input)
c     crv             = curve array (input)
c     nnpcrv          = # of discretization points per crv()
c     nosl            = number of sliding interface (input)
c     ictype          = contact type (input)
c     side            = 'master' for 1st pass or 'slave' for 2nd pass (input)
c     time            = current solution time (input)
c     ncycle          = number of current cycle (input)
c     dt2             = time step size at n+1/2 (input)
c     fric1           = static friction coefficient FS from keyword (input)
c     fric2           = dynamic friction coefficient FD from keyword (input)
c     fric3           = decay constant DC from keyword (input)
c     fric4           = viscous friction coefficient VC from keyword (input)
c     lsv             = master segment number (input)
c     idele8          = external user element number of master segment (input)
c     sfac1           = Coulomb friction scale factor FSF from keyword (input)
c     sfac2           = Viscous friction scale factor VSF from keyword (input)
c     insv            = slave node user id (input) only for output purpose
c     fni             = normal force (input)
c     areas           = slave node area (input)
c     stfk            = penalty stiffness (input)
c     ix1,ix2,ix3,ix4 = master segment nodes - internal node numbers (input)
c                       (external user node numbers: ix1ext=lqfinv8(ix1,1) ...)
c     aream           = master segment area (input)
c     rn1,rn2,rn3     = x,y,z components of master segment normal (input)
c     ue(ne)          = user defined friction element data (input)
c                         1 - plastic strain (not yet supported for SOFT=4)
c     ne              = number of user defined friction element variables (input)
c
      implicit none
      include 'nlqparm'
      real fstt,fdyn
      integer lqfinv
      external lqfinv
      integer nc,nh,nosl,ictype,ncycle,lsv,insv,ix1,ix2,ix3,ix4
      integer ne
      real uc(nc),uh(nh),ue(ne)
      real prs,temp,v,time,dt2,fric1,fric2,fric3,fric4,sfac1,sfac2,
     .     fni,areas,stfk,aream,rn1,rn2,rn3
      real crv(lq1,2,*)
      integer nnpcrv(*)
      character*(*) side
      INTEGER8 idele8
c
c     set coefficients to keyword values
c
      fstt=fric1
      fdyn=fric2
c
c     define frictional coefficients as function of temperature via load curve
c
c     call crvval(crv,uc(1),temp,fstt,fdyn)
c     fdyn=fstt
c
      return
      end
#else
      subroutine usrfrc(nosl,time,ncycle,dt2,insv,areas,xs,ys,zs,
     . lsv,ix1,ix2,ix3,ix4,aream,xx1,xx2,xx3,stfn,stf,fni,
     . dx,dy,dz,fdt2,ninput,ua,side,iisv5,niisv5,rn1,rn2,rn3,fric1,
     . fric2,fric3,fric4,bignum,fdat,iseg,fxis,fyis,fzis,ss,tt,
     . ilbsv,stfk,frc,numnp,npc,pld,lcfst,lcfdt,temp,temp_bot,
     . temp_top,isurface,nhv,ufhist,neps,ufeps,sfacn3,thkn,tk)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     user subroutine for interface friction control (smp)
c
c     note:  ls-dyna3d uses an internal numbering system to
c            accomodate arbitrary node numbering.   to access
c            information for user node n, address array location m,
c            m=lqf8(n,1). to obtain user node number, n,
c            corresponding to array address m, set n=lqfinv8(m,1)
c
c     arguments:
c
c          nosl       =number of sliding interface
c          time       =current solution time
c          ncycle     =ncycle number
c          dt2        =time step size at n+1/2
c          insv       =slave node array where the nodes are stored
c                      in ls-dyna3d internal numbering.  User numbers
c                      are given by function: lqfinv8(insv(ii),1)
c                      for slave node ii.
c          areas(ii)  =slave node area (interface types 5&10 only) for
c                      slave node ii (forming contact, too)
c          xs(ii)     =x-coordinate slave node ii (projected)
c          ys(ii)     =y-coordinate slave node ii (projected)
c          zs(ii)     =z-coordinate slave node ii (projected)
c          lsv(ii)    =master segment number for slave node ii
c          ix1(ii), ix2(ii), ix3(ii), ix4(ii)
c                     =master segment nodes in ls-dyna3d internal
c                      numbering for slave node ii
c          aream(ii)  =master segment area for slave node ii.
c                      (not for forming contacts)
c          xx1(ii,4)  =x-coordinates master surface (projected) for
c                      slave node ii (not for forming contacts)
c          xx2(ii,4)  =y-coordinates master surface (projected) for
c                      slave node ii (not for forming contacts)
c          xx3(ii,4)  =z-coordinates master surface (projected) for
c                      slave node ii (not for forming contacts)
c          stfn       =slave node penalty stiffness
c          stf        =master segment penalty stiffness
c          fni        =normal force
c          dx,dy,dz   =relative x,y,z-displacement between slave node and
c                      master surface.  Multipling by fdt2 defines the
c                      relative velocity.
c          rn1,rn2,rn3=x,y, and z components of master segments normal
c                      vector
c          sfacn3     =scaling factor for implicit calculations
c          thkn(ii)   =slave node thickness
c          tk(ii)     =master segment thickness
c
c***********************************************************************
c       frictional coefficients defined for the contact interface
c
c          fric1      =static friction coefficient
c          fric2      =dynamic friction coefficient
c          fric3      =decay constant
c          fric4      =viscous friction coefficient (setting fric4=0
c                      turns this option off)
c
c***********************************************************************
c
c          bignum     =0.0 for one way surface to surface and
c                      for surface to surface, and 1.e+10 for nodes
c                      to surface contact
c                     =0.0 for forming contact
c          ninput     =number of variables input into ua
c          ua(*)      =users' array, first ninput locations
c                      defined by user.  the length of this
c                      array is defined on control card 10.
c                      this array is unique to interface nosl.
c
c          side       ='master' for first pass.  the master
c                       surface is the surface designated in the
c                       input
c                     ='slave' for second pass after slave and
c                      master surfaces have be switched for
c                      the type 3 symmetric interface treatment.
c
c          iisv5      =an array giving the pointers to the active nodes
c                      in the arrays.
c
c          niisv5     =number of active nodes
c
c          fdat       =contact history data array containing:
c                       1 - ss (isoparametric coord x)
c                       2 - tt (isoparametric coord y)
c                       3 - old friction x-force
c                       4 - old friction y-force
c                       5 - old friction z-force
c          iseg       =contact master segment from previous step.
c          fxis       =slave node force component in global x dir.
c                      to be updated to include friction
c          fyis       =slave node force component in global y dir.
c                      to be updated to include friction
c          fzis       =slave node force component in global z dir.
c                      to be updated to include friction
c          ss(ii)     =s contact point (-1 to 1) in parametric coordinates
c                      for slave node ii.
c          tt(ii)     =t contact point (-1 to 1) in parametric coordinates
c                      for slave node ii.
c          ilbsv(ii)  =pointer for node ii into global arrays.
c          stfk(ii)   =penalty stiffness for slave node ii which was used
c                      to compute normal interface force.
c          frc(1,lsv(ii))
c                    =Coulomb friction scale factor for segment lsv(ii)
c          frc(2,lsv(ii))
c                    =viscous friction scale factor for segment lsv(ii)
c
c***********************************************************************
c       parameters for a coupled thermal-mechanical contact
c
c          numnp      = number of nodal points in the model
c          npc        = load curve pointer
c          pld        = load curve (x,y) data
c          lcfst(nosl)= load curve number for static coefficient of
c                       friction versus temperture for contact
c                       surface nosl
c          lcfdt(nosl)= load curve number for dynamic coefficient of
c                       friction versus temperture for contact
c                       surface nosl
c          temp(j)    = temperature for node point j
c          temp_bot(j)= temparature for thick thermal shell bottom
c                       surface
c          temp_top(j)= temparature for thick thermal shell top
c                       surface
c          numsh12    = number of thick thermal shells
c          itopaz(1)  = 999 ==> thermal-mechanical analysis
c          isurface   = thick thermal shell surface pointer
c
c***********************************************************************
c
c          nhv = size of user friction history data array
c          ufhist(nhv,ilbsv(ii)) = user friction history data array
c          neps= size of element history data array
c          ufeps(neps,ilbsv(ii)) = element history data array
c                                  1 - plastic strain (history var #1)
c
c***********************************************************************
      include 'nlqparm'
#include "bk05.inc"
#include "bk07.inc"
#include "bk08.inc"
#include "bk13.inc"
      common/bmcntc/ibmcnt(101)
      common/thropt/itopaz(101),iprm_thrm(20),rprm_thrm(20)
      common/blk03/rdumy1(2),idumy6(25),numsh12,ndtot,nsl_th,lenhsv,
     . numel2,numel4,numel6,numel8,numel10,numel20
      dimension insv(*),areas(*),xs(*),ys(*),zs(*),lsv(*),ix1(*),
     . ix2(*),ix3(*),ix4(*),aream(*),xx1(nlq,*),xx2(nlq,*),xx3(nlq,*),
     . stfn(*),stf(*),fni(*),dx(*),dy(*),dz(*),ua(*),iisv5(*),
     . fdat(5,*),iseg(*),fxis(*),fyis(*),fzis(*),ss(*),tt(*),
     . ilbsv(*),stfk(*),frc(2,*),npc(*),pld(2,*),lcfst(*),lcfdt(*),
     . temp(*),temp_bot(*),temp_top(*),isurface(numnp,*),
     . ufhist(nhv,*),ufeps(neps,*),thkn(*),tk(*)
      real rn1(*),rn2(*),rn3(*)
      integer ncycle
      character*(*) side
      INTEGER8 lqfinv8
c
c     set coefficient of friction multipliers to 1.
      fstm=1.
      fdym=1.
c
c     loop on the number of active nodes niisv5 in contact
      do kk = 1, niisv5
c
c     get the location pointer ii for this active node
      ii = iisv5(kk)
c
c------------------------- begin thermal -------------------------------
c     recalculate coefficient of friction multipliers if
c     function of temperature
c
      if (itopaz(1) .eq. 999) then
c
c     get the node number jj
c     the user node number is obtained using node_user=lqfinv8(jj,1)
      jj = insv(ii)
c
c     get the temperature of the active node.
c     use the array isurface to distinquish between a bottom,
c     midplane, or top surface node for a thick thermal shell
c      if (ibmcnt(19).ne.0) then
      if (numsh12.ne.0) then
      if (isurface(jj,nosl) .lt. 0) then
        tnode=temp_bot(jj)
      elseif (isurface(jj,nosl) .gt. 0) then
        tnode=temp_top(jj)
      else
        tnode=temp(jj)
      endif
      else
        tnode=temp(jj)
      endif
c
c     get the shape functions for this active node projected
c     onto the surface it contacts
      h1=.25*(1.-ss(ii))*(1.-tt(ii))
      h2=.25*(1.+ss(ii))*(1.-tt(ii))
      h3=.25*(1.+ss(ii))*(1.+tt(ii))
      h4=.25*(1.-ss(ii))*(1.+tt(ii))
c     get the segment temperature at the interpolation point
c     use the array isurface to distinquish between a bottom,
c     midplane, or top surface node for a thick thermal shell
c      if (ibmcnt(19).ne.0) then
      if (numsh12.ne.0) then
      if (isurface(ix1(ii),nosl) .lt. 0) then
        tseg=h1*temp_bot(ix1(ii))+h2*temp_bot(ix2(ii))
     1   +h3*temp_bot(ix3(ii))+h4*temp_bot(ix4(ii))
      elseif (isurface(ix1(ii),nosl) .gt. 0) then
        tseg=h1*temp_top(ix1(ii))+h2*temp_top(ix2(ii))
     1   +h3*temp_top(ix3(ii))+h4*temp_top(ix4(ii))
      else
        tseg=h1*temp(ix1(ii))+h2*temp(ix2(ii))
     1   +h3*temp(ix3(ii))+h4*temp(ix4(ii))
      endif
      else
        tseg=h1*temp(ix1(ii))+h2*temp(ix2(ii))
     1   +h3*temp(ix3(ii))+h4*temp(ix4(ii))
      endif
c
c     get the average contact temperature
      tavg=.5*(tnode+tseg)
c
c     get static coefficient of friction temperature multiplier fstm
      call value (nosl,lcfst,npc,pld,tavg,tavg,fstm,valp,1)
c     get dynamic coefficient of friction temperature multiplier fdym
      call value (nosl,lcfdt,npc,pld,tavg,tavg,fdym,valp,1)
      endif
c------------------------- end thermal ---------------------------------
c---------------------- begin mechanical -------------------------------
c
c     current & accumulated sliding distance
      ddist=sqrt(dx(ii)**2+dy(ii)**2+dz(ii)**2)
      ufhist(1,ilbsv(ii))=ufhist(1,ilbsv(ii))+ddist
      dist=ufhist(1,ilbsv(ii))
c
c     current sliding velocity
      vel=fdt2*ddist
c
c     current contact pressure = normal force / slave node area
      press=-fni(ii)/areas(ii)
c
c     master and slave node area
c     (needed for frict. force limit, see parameter VC in the manual)
      aream(ii)=aream(ii)*fric4*frc(2,lsv(ii))+
     .          1.e+16*(.5-sign(.5,fric4-1.e-10))
      areas(ii)=areas(ii)*fric4*frc(2,lsv(ii))+
     .          1.e+16*(.5-sign(.5,fric4-1.e-10))+bignum
      areas(ii)=areas(ii)+1.e+16*(.5-sign(.5,areas(ii)-1.e-10))
c
c     new contact force = old + current penalty force
      fdat(3,ii)=fdat(3,ii)+stfk(ii)*dx(ii)*sfacn3
      fdat(4,ii)=fdat(4,ii)+stfk(ii)*dy(ii)*sfacn3
      fdat(5,ii)=fdat(5,ii)+stfk(ii)*dz(ii)*sfacn3
c
c     tangential part = total - normal part
      proj=fdat(3,ii)*rn1(ii)+fdat(4,ii)*rn2(ii)+fdat(5,ii)*rn3(ii)
      fdat(3,ii)=fdat(3,ii)-proj*rn1(ii)
      fdat(4,ii)=fdat(4,ii)-proj*rn2(ii)
      fdat(5,ii)=fdat(5,ii)-proj*rn3(ii)
c
c     scalar tangential force = func(fric.coeff.) * normal force
      fmax=-(fdym*fric2+(fstm*fric1-fdym*fric2)*exp(-fric3*vel))
     .     *fni(ii)*frc(1,lsv(ii))
c
c     limit the friction force
      fmax=min(fmax,aream(ii),areas(ii))
      fmag=sqrt(fdat(3,ii)**2+fdat(4,ii)**2+fdat(5,ii)**2)+1.e-06
c
c     tangential force due to friction
      if (fmag.gt.fmax .and. fmag.ne.0.0) then
        sclf=fmax/fmag
        fdat(3,ii)=sclf*fdat(3,ii)
        fdat(4,ii)=sclf*fdat(4,ii)
        fdat(5,ii)=sclf*fdat(5,ii)
      endif
c
c     slave node force components in global directions
      fxis(ii)=fxis(ii)+fdat(3,ii)
      fyis(ii)=fyis(ii)+fdat(4,ii)
      fzis(ii)=fzis(ii)+fdat(5,ii)
c
      iseg(ii)  =lsv(ii)
      fdat(1,ii)=ss(ii)
      fdat(2,ii)=tt(ii)
c----------------------- end mechanical --------------------------------
c
      enddo
c
      return
      end
#endif
      subroutine userwear(wdepth,wtype,cside,
     1     wparm,mwearu,
     2     whisv,lwearu,
     3     nx,ny,nz,
     4     dx,dy,dz,
     5     press,temp,stemp,time,dt,
     6     crv,nnpcrv,
     7     idnode,idcntct,isslave)
      implicit none
      include 'nlqparm'
      real wdepth
      integer wtype,cside
      integer mwearu,lwearu
      real wparm(*),whisv(*)
      real nx,ny,nz
      real dx,dy,dz
      real press,temp,stemp,time,dt
      real crv(lq1,2,*)
      integer nnpcrv(*)
      integer idnode,idcntct,isslave
c     
c     User defined wear routine, see *CONTACT_ADD_WEAR, for processing
c     wear on a contact node. This routine is called for wear interfaces
c     with WTYPE<0, this denotes the type. The purpose is to increment the
c     parameter wdepth, which is the wear depth in length units. The
c     following is a list of input parameters, after that is a
c     sample code for the simple Archard's law corresponding to
c     WTYPE=0 and positive as well as negative parameters.
c
c     outputs:
c     
c     wdepth        = wear depth to be incremented (length unit)
c     whisv(lwearu) = list of wear history variables, initialized to zero,
c                     to be updated in this routine
c     
c     inputs:
c
c     wtype         = wear type, this is the value of WTYPE from the input
c                     and can be used to distinguish wear laws, see sample routine below
c     cside         = 0 if contact node is on slave side of the contact, 1 otherwise
c     wparm(mwearu) = list of input wear parameters from keyword deck
c     mwearu        = number of input wear parameters requested, for this wear interface
c     lwearu        = number of input wear histories requested, for this contact node
c     nx,ny,nz      = components of normal in the contact, vector is normalized,
c                     normal is pointing in the direction defined by the orientation
c                     of underlying element connectivity
c     dx,dy,dz      = relative incremental displacement vector of the contact point,
c                     in both normal and tangential direction (length unit), vector
c                     is defined as the contact nodal incremental displacement minus
c                     the incremental displacement of the opposite contact point
c     press         = contact pressure in the contact node, positive (force per area unit)
c     temp          = temperature in the contact node, zero if no thermal calculation
c     stemp         = temperature on the contact segment on opposing side,
c                     zero if no thermal calculation
c     time          = current simulation time
c     dt            = current time step
c     crv(lq1,*,*)  = curve array, use subroutines crvval and tabval for evaluating data,
c                     see Appendix A
c     nnpcrv(*)     = data required for curve data
c     idnode        = global node ID (from keyword input) for the node
c     idcntct       = contact ID (from keyword input) for the underlying contact
c     isslave       = 1 if the contact node is on the slave side of the contact interface, 
c                     0 otherwise
c     
c     Sample routine follows, Archard's law in two forms, mwearu=3, lwearu=0, 
c
c     Wear parameters for wtype=-1 are
c
c     wparm(1)      = dimensionless coefficient
c     wparm(2)      = hardness of slave surface
c     wparm(3)      = hardness of master surface
c
c     Wear parameters for wtype=-2 are
c
c     wparm(1)      = dimensionless coefficient table id as function
c                     of pressure and sliding velocity
c     wparm(2)      = hardness of slave surface curve id as function 
c                     of temperature
c     wparm(3)      = hardness of master surface curve id as function
c                     of temperature
c
c     Declare variables
      real dn,tx,ty,tz,tt,c,h
      real vt,dcdvt,dcdpress,dhdtemp
c
c     First do stuff necessary for both wtype=-1 and wtype=-2
c
c     Calculate relative normal displacement
      dn=dx*nx+dy*ny+dz*nz
c     Calculate relative tangential displacement
      tx=dx-dn*nx
      ty=dy-dn*ny
      tz=dz-dn*nz
c     Calculate relative sliding distance
      tt=sqrt(tx**2+ty**2+tz**2)
c     For law with constant parameters
      if (wtype.eq.-1) then
c     Dimensionless coefficient
         c=wparm(1)
c     Make sure correct hardness is used
         if (cside.eq.0) then
            h=wparm(2)
         else
            h=wparm(3)
         endif
      endif
c     store opposite side temperature as first history variable
      if (lwearu.ge.1) whisv(1)=stemp
c     For law with dependent parameters
      if (wtype.eq.-2) then
c     Sliding velocity
         vt=tt/max(dt,1.e-20)
c     Evaluate coefficient from table
         call tabval(crv,nnpcrv,wparm(1),vt,c,dcdvt,press,dcdpress)
c     Make sure correct hardness is used
         if (cside.eq.0) then
            call crvval(crv,nnpcrv,wparm(2),temp,h,dhdtemp)
         else
            call crvval(crv,nnpcrv,wparm(3),temp,h,dhdtemp)
         endif
      endif
c     Increment wear depth by the Archard's law parameters
      wdepth=wdepth+c*tt*press/h
c     Done, wdepth can be post-processed in the intfor file, see DATABASE_EXTENT_INTFOR
      return
      end
c     
      subroutine useradap(eso,b,maxsiz64)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
#include "bk05.inc"
#include "bk07.inc"
#include "bk08.inc"
#include "bk13.inc"
#include "dynmem.inc"
#include "sorter.inc"
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,nusa
      integer*8 maxsiz64
c
c
      dimension eso(7,*),b(1)
c
      call userrefin (b(n1),b(n4f),r_mem(dm_x),i_mem(dm_smtnum),b(lc1s),
     +  b(ns06),b(n1+nmmat),b(ns13),b(ns14),eso,b(n4d+3*nmmat),maxsiz64)
c
      return
      end
      subroutine userrefin (mtype,csprop,x,mtnum,ixp,auxvec,
     + ishlfm,nshpnt,lochvs,eso,nconstp,maxsiz64)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'bk06.inc'
      include 'bk19.inc'
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,nusa
c     common/bk28/summss,xke,xpe,time
      common/daimlr/noddie,neldie,ncunk,npersv,ivtmem
      common/ellcm2/nelnt,nelet,nstart,nelbt
      common/xcntrl/kontrl(200)
c
      real*8 x
      dimension mtype(*),csprop(24,*),x(3,*),mtnum(*),ixp(5,*),
     +          auxvec(*),ishlfm(*),nshpnt(*),lochvs(*),eso(7,*),
     +          stress(7,20),nconstp(*)
      dimension output(7)
      integer   iob28
      integer*8 maxsiz64
      integer*8 izero64,jadd64
      data izero64/0/
      data init/0/
      data output/7*0.0/
c
c     Variables in the subroutine:
c
c        numnp=number of nodal points
c        numels-neldie=number of shell elements
c
c        x(1,i)=x-coordinate of nodal point i
c        x(2,i)=y-coordinate of nodal point i
c        x(3,i)=z-coordinate of nodal point i
c
c
c        ixp - shell element connectivities (ixp(1,*)=part ID)
c                                           (ixp(2,*)=node 1)
c                                           (ixp(3,*)=node 2)
c                                           (ixp(4,*)=node 3)
c                                           (ixp(5,*)=node 4)
c
c        ipt=integration point number
c        nip=number of integration points
c
c        stress(1,ipt)=local x  stress
c        stress(2,ipt)=local y  stress
c        stress(3,ipt)=local z  stress
c        stress(4,ipt)=local xy stress
c        stress(5,ipt)=local yz stress
c        stress(6,ipt)=local zx stress
c        stress(7,ipt)=effective plastic strain
c
c     Write file with data for useradaptivity:
c
c     structure of file:
c     line 1: Number of nodes, number of elements (shell)
c     part 2: Nodal points (numnp lines)
c     part 3: Material # + Element connectivities (numels lines)
c     part 4: Element materials and stresses (numels*nip lines)
c             Material #,
c             Integration point #
c             local sigma_xx
c             local sigma_yy
c             local sigma_zz
c             local sigma_xy
c             local sigma_yz
c             local sigma_zx
c             effective plastic strain
c
      call lsopen(21,'adapt.user','formatted','unknown',ierr,0,0)
      write(21,900) numnp,numels-neldie
      write(21,915) (x(1,n),x(2,n),x(3,n),n=1,numnp)
      write(21,905) (ixp(1,nshpnt(n)),ixp(2,nshpnt(n)),
     +               ixp(3,nshpnt(n)),ixp(4,nshpnt(n)),
     +               ixp(5,nshpnt(n)),n=1,numels-neldie)
c
      do 100 n=1,numels-neldie
      nel =nshpnt(n)
      lav =lochvs(nel)
      mx=ixp(1,nel)
      if (mx.eq.0) go to 100
      iop=ishlfm(mx)
      nip=csprop(2,mx)
      mte=mtype(mx)
      nmtcon=7+nconstp(mx)
c
      do 50 ipt=1,nip
      if (mte.eq.20) then
        stress(1,ipt)=0.
        stress(2,ipt)=0.
        stress(3,ipt)=0.
        stress(4,ipt)=0.
        stress(5,ipt)=0.
        stress(6,ipt)=0.
        stress(7,ipt)=0.
      else
        call blkcpy (auxvec(lav),stress(1,ipt),7)
      endif
c
      write(21,540) mx,ipt,(stress(j,ipt),j=1,7)
      lav=lav+nmtcon
   50 continue
  100 continue
      close(21)
c
c     External program for caclulating refinement indicators
c     can be called here
c
c     call my_system('External program')
c
c     Write adapt.err file for userspecified adaptivity
c
c     The value of enerden specify
c     which element which will be refined
c
c     if enerden < 1 no refinement of element
c     if enerden > 1 refinement of element
c
c
c
c     Or refinement indicators can be calculated
c     inside the next loop
c
c     call llsopen(28,'adapt.err','unformatted','unknown',ierr,0,0)
c     rewind 28
#ifdef MPP
      call rwabsf64(iob28,'adapt.err',0,0,izero64)
#else
      call rwabsf64(iob28,'adapt.err',28,2**9,maxsiz64)
#endif
      nn=numels-neldie-nelet
      jadd64=0
      call wrabsf64(iob28,nn,1,jadd64)
      jadd64=1
c
      enerden=0.
      do 200 n=1,numels-neldie-nelet
c
c     calculate refinement for each element based
c     on the data in the "adapt.user" file
c
c     calculate the value of enerden here
c     and the value of enerden will be written to
c     the file "adapt.err"
c
      output(1)=enerden
      if (init .eq. 0) then
c       write(28) enerden,0.0,0.0,enerden,0.0,0.0,0.0
        output(4)=enerden
      else
c       write(28) enerden,0.0,0.0,eso(1,n),0.0,0.0,0.0
        output(4)=eso(1,n)
      endif
      call wrabsf64(iob28,output,7,jadd64)
      jadd=jadd+7
  200 continue
      kontrl(42)=jadd
      init = init+1
c     close(28)
      jadd64=0
      call rwabsf64(iob28,'keep',0,0,0_8)
c
      return
  540 format(2i5,7(1pe15.7))
  900 format(10i8)
  905 format(5i10)
c 910 format(6(1pe12.5))
  915 format(3(1pe15.7))
      end
      subroutine usrmat (lft,llt,cm,bqs,capa,eltype,mt,ipt,
     . npc,plc,crv,nnpcrv,rcoor,scoor,tcoor,nnm1,nip,ipt_thk)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
#include "bk07.inc"
#include "bk08.inc"
      include 'memaia.inc'
      include 'umatss.inc'
      include 'sph.inc'
#include "sorter.inc"
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     1 ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     2 grvity,idirgv,nodspc,nspcor,numelh10,numels8,numelsx,numnpin,
     3 numelhx,nmuxct
      dimension cm(*),bqs(*),crv(lq1,2,*),npc(*),plc(*)
      integer nnpcrv(*)
      character*(*) eltype
c
      if (eltype.eq.'solid'.or.eltype.eq.'shl_t'.or.
     1     eltype.eq.'sld2d'.or.eltype.eq.'tshel'.or.
     2     eltype.eq.'sldax') then
        call urmathn (lft,llt,cm,bqs,mt,crv,nnpcrv,npc,plc,nnm1,
     1   rcoor,scoor,tcoor,ia(n4d+3*nmmat),nip,ipt,eltype,
     2   a(nh13+numelh),a(ns13+numels),a(nt13+numelt))
      elseif (eltype.eq.'sph  ') then
        call urmathn (lft,llt,cm,bqs,mt,crv,nnpcrv,npc,plc,nnm1,
     1   rcoor,scoor,tcoor,ia(n4d+3*nmmat),nip,ipt,eltype,
     2   a(ldsph00),a(ns13+numels),a(nt13+numelt))
      elseif (eltype.eq.'shell') then
        call urmats (lft,llt,cm,capa,mt,crv,nnpcrv,
     .   ipt,rcoor,scoor,tcoor,
     1   nnm1,ia(n4d+3*nmmat),nip,ipt_thk,npc,plc,eltype,
     2   a(ns13+numels))
      elseif (eltype.eq.'hbeam ') then
        call urmatb (lft,llt,cm,capa,mt,crv,nnpcrv,
     .   nnm1,ia(n4d+3*nmmat),
     1   npc,plc,eltype,a(nb13+numelb))
      elseif (eltype.eq.'dbeam') then
        call urmatd (lft,llt,cm,bqs,mt,crv,nnpcrv,nnm1,ia(n4d+3*nmmat),
     1   npc,plc,eltype,a(nb13+numelb))
      elseif (eltype.eq.'tbeam') then
        call urmatt (lft,llt,cm,bqs,mt,crv,nnpcrv,nnm1,ia(n4d+3*nmmat),
     1   npc,plc,eltype,a(nb13+numelb))
      endif
c
      return
      end
      subroutine urmathn(lft,llt,cm,bqs,mt,crv,nnpcrv,npc,plc,nnm1,
     1 rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nhxbwp,nshbwp,ntsbwp)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
#include "bk07.inc"
#include "bk13.inc"
      include 'bk19.inc'
      include 'iounits.inc'
      include 'memaia.inc'
      include 'umatss.inc'
#include "dynmem.inc"
#include "implicit1.inc"
C_TASKCOMMON (aux2loc)
C_TASKCOMMON (aux14loc)
C_TASKCOMMON (aux18loc)
C_TASKCOMMON (aux33loc)
C_TASKCOMMON (aux34loc)
C_TASKCOMMON (aux35loc)
C_TASKCOMMON (aux40loc)
C_TASKCOMMON (aux41loc)
C_TASKCOMMON (bk36loc)
C_TASKCOMMON (failcmloc)
C_TASKCOMMON (failuloc)
C_TASKCOMMON (soundloc)
C_TASKCOMMON (subtssloc)
C_TASKCOMMON (vect32loc)
C_TASKCOMMON (bux13loc)
      common/bux13loc/ eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      common/aux14loc/sig1(nlq),sig2(nlq),sig3(nlq),sig4(nlq),sig5(nlq),
     & sig6(nlq),epsps(nlq),hsvs(nlq,NHISVAR)
      common/aux18loc/dd(nlq),def(nlq),ddq(nlq)
      common/aux2loc/d1(nlq),d2(nlq),d3(nlq),d4(nlq),d5(nlq),d6(nlq),
     1 wzzdt(nlq),wyydt(nlq),wxxdt(nlq),einc(nlq)
      common/aux33loc/ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq),ix5(nlq),
     1 ix6(nlq),ix7(nlq),ix8(nlq),mxt(nlq)
      common/aux34loc/xioff(nlq),dx(nlq)
      common/aux35loc/rhoa(nlq),cb(nlq),davg(nlq),p(nlq)
      common/aux40loc/
     1 a11(nlq),a12(nlq),a13(nlq),a21(nlq),a22(nlq),a23(nlq),
     2 a31(nlq),a32(nlq),a33(nlq),z11(nlq),z12(nlq),z13(nlq),
     3 z21(nlq),z22(nlq),z23(nlq),z31(nlq),z32(nlq),z33(nlq)
      common/aux41loc/qq1(nlq),cbb(nlq),aj1(nlq)
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,nusa
      common/bk36loc/index
      common/soundloc/ss(nlq),sndsp(nlq),diagm(nlq),sarea(nlq),
     . dxl(nlq)
      common/subtssloc/dt1siz(nlq)
      common/vect32loc/g11(nlq),g21(nlq),g31(nlq),g12(nlq),g22(nlq),
     .     g32(nlq),g13(nlq),g23(nlq),g33(nlq),ifg(nlq)
c
      common/bk01/itherm,itemp,ntmp0,ntmp1,itempan,itempdr,itmpe(44)
      common/bk26/begtim,nintcy
      common/bk28/summss,xke,xpe,tt,xte0,erodeke,erodeie
      logical failur
      common/failcm/failur
      common/failcmloc/ifail(nlq)
      common/failuloc/sieu(nlq),fail(nlq),ifaili(nlq),nfipts(nlq)
      common/numcpu/ncpu,ncpua,ncpub,lenvec(8)
      common/shlopt/istrn,istupd,ibelyts,miter,wrpang,ipstpd,intsts,
     1 nodsts,intstn,nodstn,jstrn
      common/thropt/itopaz(101)
c
      dimension cm(*),bqs(*),crv(lq1,2,*),npc(*),plc(*),nconstp(*)
      integer nnpcrv(*)
      character*(*) eltype
      dimension nhxbwp(*),nshbwp(*),ntsbwp(*)
c     dimension eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      dimension q21(nlq),q22(nlq),q23(nlq)
      dimension thhsv(nlq,100),thhsvi(100)
      dimension elsizv(nlq)
      INTEGER8 idelev(nlq),idele
c
      real mlt1,mlt2
      logical failel,failels(nlq),reject
      INTEGER8 lqfinv8
c     data temps/nlq*0.0/
c
      capa=0.0
c
c     processing elements
c
      mx=48*(mxt(lft)-1)
      gm=cm(mx+abs(ishrmp(mt)))
      bk=cm(mx+ibulkp(mt))
      ss(lft)=bk+4.*gm/3.
c
c     energy calculations
c
      do 20 i=lft,llt
      cb(i)=ss(lft)
      einc(i)=(d1(i)*sig1(i)+d2(i)*sig2(i)+d3(i)*sig3(i)+
     .        d4(i)*sig4(i)+d5(i)*sig5(i)+
     .        d6(i)*sig6(i)+dd(i)*bqs(i))*dt1siz(i)
 20   continue
c
c     Determining number of requested history variables
c
      no_hsvs=nconstp(mxt(lft))
      
      if (ivumat(mt).ne.0) no_hsvs=no_hsvs-9-6
      if (lenvec(5) .eq.1) no_hsvs=no_hsvs-6
      if (iorien(mt).ne.0) no_hsvs=no_hsvs-6
      if (ihyper(mt).ne.0) no_hsvs=no_hsvs-9
      if (iorien(mt).ne.0.and.ihyper(mt).lt.0) no_hsvs=no_hsvs-6
      if (ifipts(mt).lt.0) no_hsvs=no_hsvs-1
      if (jstrn.ne.0) no_hsvs=no_hsvs-6*jstrn 
c
c     initialization of indeces
c
      if(ivumat(mt).ne.0) then

        if (iorien(mt).ne.0) then
          if (ihyper(mt).ne.0) then
            ind_defgrad=1+no_hsvs
            ind_ortho=ind_defgrad+9
          else
            ind_ortho=1+no_hsvs
          endif
          if (ifipts(mt).lt.0) then
            ind_fipts=ind_ortho
            ind_ortho=ind_ortho+1
          endif
          ind_q11=ind_ortho
          ind_q12=ind_ortho+1
          ind_q13=ind_ortho+2
          ind_q31=ind_ortho+3
          ind_q32=ind_ortho+4
          ind_q33=ind_ortho+5
          ind_last=ind_q33
          if (ihyper(mt).lt.0) then
            ind_qstore=ind_last+1
            ind_last=ind_last+6
          endif
        else
          if (ihyper(mt).ne.0) then
            ind_defgrad=1+no_hsvs
            ind_last=ind_defgrad+8
          else
            ind_last=no_hsvs
          endif
          if (ifipts(mt).lt.0) then
            ind_fipts=ind_last+1
            ind_last=ind_last+1
          endif
        endif
        ind_vumat=ind_last+1
        ind_last=ind_vumat+9+6-1       
c        ind_last=ind_last+1
c        write(*,'(A,i5,A,i5,A,i5,a,i5,a,i5,a,i5)')'NHSV:',no_hsvs,
c     $       'DF',ind_defgrad,'ORTH',ind_ortho,'VMAT',ind_vumat,
c     $       'LST',ind_last,'QSTO',ind_qstore

      else

        if (iorien(mt).ne.0) then
          if (ihyper(mt).ne.0) then
            ind_defgrad=1+no_hsvs
            ind_ortho=ind_defgrad+9
          else
            ind_ortho=1+no_hsvs
          endif
          if (ifipts(mt).lt.0) then
            ind_fipts=ind_ortho
            ind_ortho=ind_ortho+1
          endif
          ind_q11=ind_ortho
          ind_q12=ind_ortho+1
          ind_q13=ind_ortho+2
          ind_q31=ind_ortho+3
          ind_q32=ind_ortho+4
          ind_q33=ind_ortho+5
          ind_last=ind_q33
          if (ihyper(mt).lt.0) then
            ind_qstore=ind_last+1
            ind_last=ind_last+6
          endif
        else
          if (ihyper(mt).ne.0) then
            ind_defgrad=1+no_hsvs
            ind_last=ind_defgrad+8
          else
            ind_last=no_hsvs
          endif
          if (ifipts(mt).lt.0) then
            ind_fipts=ind_last+1
            ind_last=ind_last+1
          endif
        endif

      endif

c
      if (ind_last.gt.NHISVAR) then
c       write(iotty,900) lqfmiv(mxt(lft))
c       write(iohsp,900) lqfmiv(mxt(lft))
c       write(iomsg,900) lqfmiv(mxt(lft))
c       call adios(TC_ERROR)
        ierdat(1)=lqfmiv(mxt(lft))
        call lsmsg(3,MSG_SOL+1149,ioall,ierdat,rerdat,cerdat,0)
      endif
c
c$$$      if (tt.eq.begtim.and.nintcy.eq.0) then
c$$$c
c$$$c       move transformation matrix in history variables array
c$$$c
c$$$        if (iorien(mt).ne.0) then
c$$$          do i=lft,llt
c$$$            hsvs(i,ind_q33)=hsvs(i,6)
c$$$            hsvs(i,ind_q32)=hsvs(i,5)
c$$$            hsvs(i,ind_q31)=hsvs(i,4)
c$$$            hsvs(i,ind_q13)=hsvs(i,3)
c$$$            hsvs(i,ind_q12)=hsvs(i,2)
c$$$            hsvs(i,ind_q11)=hsvs(i,1)
c$$$          enddo
c$$$c
c$$$          if (no_hsvs.gt.0) then
c$$$            do j=1,no_hsvs
c$$$              do i=lft,llt
c$$$                hsvs(i,j)=0.0
c$$$              enddo
c$$$            enddo
c$$$          endif
c$$$        endif
c$$$      endif
c
      if (ipt.eq.1.and.ifipts(mt).ne.0) then
        do i=lft,llt
          nfipts(i)=0
        enddo
      endif
c
c     storage of deformation gradient
c
      if (ivumat(mt).ne.0) then
        ifo11=ind_vumat
        ifo21=ind_vumat+1
        ifo31=ind_vumat+2
        ifo12=ind_vumat+3
        ifo22=ind_vumat+4
        ifo32=ind_vumat+5
        ifo13=ind_vumat+6
        ifo23=ind_vumat+7
        ifo33=ind_vumat+8
      endif
      
      if (ihyper(mt).ne.0) then
        if11=ind_defgrad
        if21=ind_defgrad+1
        if31=ind_defgrad+2
        if12=ind_defgrad+3
        if22=ind_defgrad+4
        if32=ind_defgrad+5
        if13=ind_defgrad+6
        if23=ind_defgrad+7
        if33=ind_defgrad+8
c
c
c       deformation gradient in global system at time = 0
c
        if (tt.eq.begtim.and.nintcy.eq.0) then
          do i=lft,llt
            hsvs(i,if11)=1.0
            hsvs(i,if21)=0.0
            hsvs(i,if31)=0.0
            hsvs(i,if12)=0.0
            hsvs(i,if22)=1.0
            hsvs(i,if32)=0.0
            hsvs(i,if13)=0.0
            hsvs(i,if23)=0.0
            hsvs(i,if33)=1.0
          enddo
          if (ivumat(mt).ne.0) then
            do i=lft,llt
              hsvs(i,ifo11)=1.0
              hsvs(i,ifo21)=0.0
              hsvs(i,ifo31)=0.0
              hsvs(i,ifo12)=0.0
              hsvs(i,ifo22)=1.0
              hsvs(i,ifo32)=0.0
              hsvs(i,ifo13)=0.0
              hsvs(i,ifo23)=0.0
              hsvs(i,ifo33)=1.0
            enddo
          endif
        endif
c
c       deformation gradient calculated exactly
c
        if (ifg(lft).eq.1) then
          do i=lft,llt
            hsvs(i,if11)=g11(i)
            hsvs(i,if21)=g21(i)
            hsvs(i,if31)=g31(i)
            hsvs(i,if12)=g12(i)
            hsvs(i,if22)=g22(i)
            hsvs(i,if32)=g32(i)
            hsvs(i,if13)=g13(i)
            hsvs(i,if23)=g23(i)
            hsvs(i,if33)=g33(i)
          enddo
c
c
c       update deformation gradient incrementally in global system
c       at time = n+1
c
        else
          call integf(hsvs(1,if11),hsvs(1,if21),hsvs(1,if31),
     &     hsvs(1,if12),hsvs(1,if22),hsvs(1,if32),
     &     hsvs(1,if13),hsvs(1,if23),hsvs(1,if33),lft,llt)
        endif

      endif
c
      if (iorien(mt).ne.0) then
c
c       set up global to element transformation matrix a11,...,a33
c
        if (ihyper(mt).ge.0.and.eltype.ne.'tshel') then
c        if (ihyper(mt).eq.1 .or. ihyper(mt).eq.0) then
          if (eltype.ne.'sld2d'.and.eltype.ne.'sldax') then
            call lcsm22(lft,llt)
          else
            call cs2d_g2l(lft,llt)
          endif
        endif
c     
        do 10 i=lft,llt
          q21(i)=hsvs(i,ind_q32)*hsvs(i,ind_q13)-
     1         hsvs(i,ind_q12)*hsvs(i,ind_q33)
          q22(i)=hsvs(i,ind_q33)*hsvs(i,ind_q11)-
     1         hsvs(i,ind_q13)*hsvs(i,ind_q31)
          q23(i)=hsvs(i,ind_q31)*hsvs(i,ind_q12)-
     1         hsvs(i,ind_q11)*hsvs(i,ind_q32)
 10     continue

c        write(*,*)'Q21,22,23',q21(1),q22(1),q23(1)

c     
c       set up global to material transformation matrix z11,...,z33
c
        if (ihyper(mt).ge.0) then
c     if (ihyper(mt).eq.1 .or. ihyper(mt).eq.0) then
          if (eltype.ne.'tshel') then
            do i=lft,llt
              z11(i)=a11(i)*hsvs(i,ind_q11)+a21(i)*hsvs(i,ind_q12)+
     1             a31(i)*hsvs(i,ind_q13)
              z12(i)=a12(i)*hsvs(i,ind_q11)+a22(i)*hsvs(i,ind_q12)+
     1             a32(i)*hsvs(i,ind_q13)
              z13(i)=a13(i)*hsvs(i,ind_q11)+a23(i)*hsvs(i,ind_q12)+
     1             a33(i)*hsvs(i,ind_q13)
              z21(i)=a11(i)*q21(i)+a21(i)*q22(i)+a31(i)*q23(i)
              z22(i)=a12(i)*q21(i)+a22(i)*q22(i)+a32(i)*q23(i)
              z23(i)=a13(i)*q21(i)+a23(i)*q22(i)+a33(i)*q23(i)
              z31(i)=a11(i)*hsvs(i,ind_q31)+a21(i)*hsvs(i,ind_q32)+
     1             a31(i)*hsvs(i,ind_q33)
              z32(i)=a12(i)*hsvs(i,ind_q31)+a22(i)*hsvs(i,ind_q32)+
     1             a32(i)*hsvs(i,ind_q33)
              z33(i)=a13(i)*hsvs(i,ind_q31)+a23(i)*hsvs(i,ind_q32)+
     1             a33(i)*hsvs(i,ind_q33)
            enddo
          else
            do i=lft,llt
              z11(i)=hsvs(i,ind_q11)
              z12(i)=hsvs(i,ind_q12)
              z13(i)=hsvs(i,ind_q13)
              z21(i)=q21(i)
              z22(i)=q22(i)
              z23(i)=q23(i)
              z31(i)=hsvs(i,ind_q31)
              z32(i)=hsvs(i,ind_q32)
              z33(i)=hsvs(i,ind_q33)
            enddo
          endif
c     
c     transform stresses and strains to local system
c     
          call ldsm22(lft,llt)
          call ldem22(lft,llt)
c     
c     else
c        endif
      elseif (ihyper(mt).lt.0) then
c     if ( (ihyper(mt).lt.0) .or. (ihyper(mt).eq.2) ) then
c     
          if (tt.eq.begtim.and.nintcy.eq.0) then
c     
c           set up global to element transformation matrix a11,...,a33
c
            if (eltype.ne.'tshel') then
              if (eltype.ne.'sld2d'.and.eltype.ne.'sldax') then
                call lcsm22(lft,llt)
              else
                call cs2d_g2l(lft,llt)
              endif
c
c             set up global to material transformation matrix z11,...,z33
c
              do i=lft,llt
                z11(i)=a11(i)*hsvs(i,ind_q11)+a21(i)*hsvs(i,ind_q12)+
     1               a31(i)*hsvs(i,ind_q13)
                z12(i)=a12(i)*hsvs(i,ind_q11)+a22(i)*hsvs(i,ind_q12)+
     1               a32(i)*hsvs(i,ind_q13)
                z13(i)=a13(i)*hsvs(i,ind_q11)+a23(i)*hsvs(i,ind_q12)+
     1               a33(i)*hsvs(i,ind_q13)
                z21(i)=a11(i)*q21(i)+a21(i)*q22(i)+a31(i)*q23(i)
                z22(i)=a12(i)*q21(i)+a22(i)*q22(i)+a32(i)*q23(i)
                z23(i)=a13(i)*q21(i)+a23(i)*q22(i)+a33(i)*q23(i)
                z31(i)=a11(i)*hsvs(i,ind_q31)+a21(i)*hsvs(i,ind_q32)+
     1               a31(i)*hsvs(i,ind_q33)
                z32(i)=a12(i)*hsvs(i,ind_q31)+a22(i)*hsvs(i,ind_q32)+
     1               a32(i)*hsvs(i,ind_q33)
                z33(i)=a13(i)*hsvs(i,ind_q31)+a23(i)*hsvs(i,ind_q32)+
     1               a33(i)*hsvs(i,ind_q33)
              enddo
            else
              do i=lft,llt
                z11(i)=hsvs(i,ind_q11)
                z12(i)=hsvs(i,ind_q12)
                z13(i)=hsvs(i,ind_q13)
                z21(i)=q21(i)
                z22(i)=q22(i)
                z23(i)=q23(i)
                z31(i)=hsvs(i,ind_q31)
                z32(i)=hsvs(i,ind_q32)
                z33(i)=hsvs(i,ind_q33)
              enddo
            endif
c            write(*,*)'A1:',a11(1),a22(1),a33(1)
c            write(*,*)'Z1:',z11(1),z22(1),z33(1)
c     
c     Store correct transformation used for transforming
c     the deformation gradient
c     
            do i=lft,llt
              hsvs(i,ind_qstore  )=hsvs(i,ind_q11)
              hsvs(i,ind_qstore+1)=hsvs(i,ind_q12)
              hsvs(i,ind_qstore+2)=hsvs(i,ind_q13)
              hsvs(i,ind_qstore+3)=hsvs(i,ind_q31)
              hsvs(i,ind_qstore+4)=hsvs(i,ind_q32)
              hsvs(i,ind_qstore+5)=hsvs(i,ind_q33)
              hsvs(i,ind_q11)=z11(i)
              hsvs(i,ind_q12)=z12(i)
              hsvs(i,ind_q13)=z13(i)
              hsvs(i,ind_q31)=z31(i)
              hsvs(i,ind_q32)=z32(i)
              hsvs(i,ind_q33)=z33(i)
            enddo
c     
c     set up element to material transformation matrix z11,...,z33
c     
          else
            do i=lft,llt
              z11(i)=hsvs(i,ind_q11)
              z12(i)=hsvs(i,ind_q12)
              z13(i)=hsvs(i,ind_q13)
              z21(i)=q21(i)
              z22(i)=q22(i)
              z23(i)=q23(i)
              z31(i)=hsvs(i,ind_q31)
              z32(i)=hsvs(i,ind_q32)
              z33(i)=hsvs(i,ind_q33)
            enddo
          endif
        endif
c     
c        deformation gradient transformed to local system
c        as F_bar=Q^T*F
c
         if (ihyper(mt).gt.0) then
c        if (ihyper(mt).eq.1) then
          do i=lft,llt
            f1=z11(i)*hsvs(i,if11)+
     1           z12(i)*hsvs(i,if21)+z13(i)*hsvs(i,if31)
            f2=z21(i)*hsvs(i,if11)+
     1           z22(i)*hsvs(i,if21)+z23(i)*hsvs(i,if31)
            hsvs(i,if31)=z31(i)*hsvs(i,if11)+
     1           z32(i)*hsvs(i,if21)+z33(i)*hsvs(i,if31)
            hsvs(i,if11)=f1
            hsvs(i,if21)=f2
            f1=z11(i)*hsvs(i,if12)+
     1           z12(i)*hsvs(i,if22)+z13(i)*hsvs(i,if32)
            f2=z21(i)*hsvs(i,if12)+
     1           z22(i)*hsvs(i,if22)+z23(i)*hsvs(i,if32)
            hsvs(i,if32)=z31(i)*hsvs(i,if12)+
     1           z32(i)*hsvs(i,if22)+z33(i)*hsvs(i,if32)
            hsvs(i,if12)=f1
            hsvs(i,if22)=f2
            f1=z11(i)*hsvs(i,if13)+
     1           z12(i)*hsvs(i,if23)+z13(i)*hsvs(i,if33)
            f2=z21(i)*hsvs(i,if13)+
     1           z22(i)*hsvs(i,if23)+z23(i)*hsvs(i,if33)
            hsvs(i,if33)=z31(i)*hsvs(i,if13)+
     1           z32(i)*hsvs(i,if23)+z33(i)*hsvs(i,if33)
            hsvs(i,if13)=f1
            hsvs(i,if23)=f2
          enddo
c
c     deformation gradient transformed to co-rotated coordinate system
c     as F_bar=Q^T*F*Q (suggested by Richard Hamm, P&G 2011)
        else if (ihyper(mt).eq.-2) then
          do i=lft,llt
            f11=z11(i)*hsvs(i,if11)+
     1           z12(i)*hsvs(i,if21)+z13(i)*hsvs(i,if31)
            f21=z21(i)*hsvs(i,if11)+
     1           z22(i)*hsvs(i,if21)+z23(i)*hsvs(i,if31)
            f31=z31(i)*hsvs(i,if11)+
     1           z32(i)*hsvs(i,if21)+z33(i)*hsvs(i,if31)
            f12=z11(i)*hsvs(i,if12)+
     1           z12(i)*hsvs(i,if22)+z13(i)*hsvs(i,if32)
            f22=z21(i)*hsvs(i,if12)+
     1           z22(i)*hsvs(i,if22)+z23(i)*hsvs(i,if32)
            f32=z31(i)*hsvs(i,if12)+
     1           z32(i)*hsvs(i,if22)+z33(i)*hsvs(i,if32)
            f13=z11(i)*hsvs(i,if13)+
     1           z12(i)*hsvs(i,if23)+z13(i)*hsvs(i,if33)
            f23=z21(i)*hsvs(i,if13)+
     1           z22(i)*hsvs(i,if23)+z23(i)*hsvs(i,if33)
            f33=z31(i)*hsvs(i,if13)+
     1           z32(i)*hsvs(i,if23)+z33(i)*hsvs(i,if33)
            hsvs(i,if11)=f11*z11(i)+f12*z12(i)+f13*z13(i)
            hsvs(i,if12)=f11*z21(i)+f12*z22(i)+f13*z23(i)
            hsvs(i,if13)=f11*z31(i)+f12*z32(i)+f13*z33(i)
            hsvs(i,if21)=f21*z11(i)+f22*z12(i)+f23*z13(i)
            hsvs(i,if22)=f21*z21(i)+f22*z22(i)+f23*z23(i)
            hsvs(i,if23)=f21*z31(i)+f22*z32(i)+f23*z33(i)
            hsvs(i,if31)=f31*z11(i)+f32*z12(i)+f33*z13(i)
            hsvs(i,if32)=f31*z21(i)+f32*z22(i)+f33*z23(i)
            hsvs(i,if33)=f31*z31(i)+f32*z32(i)+f33*z33(i)
          enddo
c     
c
c        deformation gradient transformed to local system
c        as F_bar=F*Q
c
c        else if (ihyper(mt).lt.0) then
        else if (ihyper(mt).eq.-1) then
          do i=lft,llt
            f1=z11(i)*hsvs(i,if11)+
     1           z12(i)*hsvs(i,if12)+z13(i)*hsvs(i,if13)
            f2=z21(i)*hsvs(i,if11)+
     1           z22(i)*hsvs(i,if12)+z23(i)*hsvs(i,if13)
            hsvs(i,if13)=z31(i)*hsvs(i,if11)+
     1           z32(i)*hsvs(i,if12)+z33(i)*hsvs(i,if13)
            hsvs(i,if11)=f1
            hsvs(i,if12)=f2
            f1=z11(i)*hsvs(i,if21)+
     1           z12(i)*hsvs(i,if22)+z13(i)*hsvs(i,if23)
            f2=z21(i)*hsvs(i,if21)+
     1           z22(i)*hsvs(i,if22)+z23(i)*hsvs(i,if23)
            hsvs(i,if23)=z31(i)*hsvs(i,if21)+
     1           z32(i)*hsvs(i,if22)+z33(i)*hsvs(i,if23)
            hsvs(i,if21)=f1
            hsvs(i,if22)=f2
            f1=z11(i)*hsvs(i,if31)+
     1           z12(i)*hsvs(i,if32)+z13(i)*hsvs(i,if33)
            f2=z21(i)*hsvs(i,if31)+
     1           z22(i)*hsvs(i,if32)+z23(i)*hsvs(i,if33)
            hsvs(i,if33)=z31(i)*hsvs(i,if31)+
     1           z32(i)*hsvs(i,if32)+z33(i)*hsvs(i,if33)
            hsvs(i,if31)=f1
            hsvs(i,if32)=f2
          enddo
        endif
      endif
c     
c     if ishrmp(mt).lt.0 then the flag for temperature is active
c     temperatures are stored in the temps array.
c
c        write(*,*)'Z2:',z11(1),z22(1),z33(1)
c        write(*,*)'dyn21_2 F:',
c     $       hsvs(1,if11),hsvs(1,if22),hsvs(1,if33)
            
      num_nods=8
      if (eltype.eq.'shl_t') num_nods=4
      if (eltype.eq.'sld2d'.or.eltype.eq.'sldax') num_nods=4
      if (ishrmp(mt).lt.0) then
        call usr_temps (a(ntmp0+1),a(n19),itemp,temps,lft,llt,ix1,ix2,
     .   ix3,ix4,ix5,ix6,ix7,ix8,num_nods,itopaz,0,0.,0.,0.)
      else
        do i=lft,llt
          temps(i)=0.0
        enddo
      endif
c
c     convert rates to increments
c
      do i=lft,llt
        d1(i)=d1(i)*dt1siz(i)
        d2(i)=d2(i)*dt1siz(i)
        if (eltype.eq.'sld2d') then
          d3(i)=0.0
        else
          d3(i)=d3(i)*dt1siz(i)
        endif
        d4(i)=d4(i)*dt1siz(i)
        d5(i)=d5(i)*dt1siz(i)
        d6(i)=d6(i)*dt1siz(i)
      enddo
c
c     Get thermal history variables from thermal
c     user material
c
      call get_thhsv(nthhsv,thhsv,itemp,lft,llt,
     1 num_nods,itopaz,nnm1,rcoor,scoor,tcoor)
c
c     characteristic element size
      do i=lft,llt
        elsizv(i)=dx(i)
      enddo
c
c     element id
      if (eltype.eq.'shl_t'.or.eltype.eq.'sld2d'.or.
     &    eltype.eq.'sldax') then
        do i=lft,llt
          ielm=nshbwp(nnm1+i)
          if (ielm.gt.0) idelev(i)=lqfinv8(ielm,4)
        enddo
      elseif (eltype.eq.'tshel') then
        do i=lft,llt
          ielm=ntsbwp(nnm1+i)
          if (ielm.gt.0) idelev(i)=lqfinv8(ielm,5)
        enddo
      elseif (eltype.eq.'sph  ') then
        do i=lft,llt
          ielm=nhxbwp(nnm1+i)
          if (ielm.gt.0) idelev(i)=lqfinv8(ielm,1)
        enddo
      else
        do i=lft,llt
          ielm=nhxbwp(nnm1+i)
          if (ielm.gt.0) idelev(i)=lqfinv8(ielm,2)
        enddo
      endif
c
c     memory pointer for extra memory data
      lcma=ncma+nmmat+ia(ncma+mxt(lft)-1)-1
c
c     vector parameter
c
      ivect=ivectr(mt)
c
c     implicit rejection flag
c
      reject=.false.
c
      mte40=mt-40
c
c     scalar umat
c
      if (ivect.eq.0) then
        do 90 i=lft,llt
c
c       gathering of variables
c
        index =i
        dt1   =dt1siz(i)
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
        epsp=epsps(i)
        temper=temps(i)
        failel=.false.
        if (ifipts(mt).lt.0) failel=hsvs(i,ind_fipts).ne.0.
        if (no_hsvs.gt.0) then
          do 30 j=1,no_hsvs
          hsv(j)=hsvs(i,j)
 30       continue
        endif
        if (ihyper(mt).ne.0) then
          do 40 j=1,9
          hsv(no_hsvs+j)=hsvs(i,ind_defgrad-1+j)
 40       continue
        endif
        if (nthhsv.gt.0) then
          do j=1,nthhsv
            thhsvi(j)=thhsv(i,j)
          enddo
        endif
        elsiz=elsizv(i)
        idele=idelev(i)
c
c       call user developed subroutines here
c
        go to (41,42,43,44,45,46,47,48,49,50), mte40
 41     call umat41 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,
     2   reject)
        go to 60
 42     call umat42 (cm(mx+1),eps,sig,hsv,dt1,
     1   capa,eltype,tt,crv,nnpcrv,a(lcma),temper,0,elsiz,idele,reject)
        go to 60
 43     call umat43 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,
     2   reject)
        go to 60
 44     call umat44 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,
     2   reject)
        go to 60
 45     call umat45 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,
     2   0,reject)
        go to 60
 46     call umat46 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,failel,crv,nnpcrv,
     2   a(lcma),thhsvi,nthhsv,0,elsiz,idele,reject)
        go to 60
 47     call umat47 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,
     2   reject)
        go to 60
 48     call umat48 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,
     2   reject)
        go to 60
 49     call umat49 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,
     2   reject)
        go to 60
 50     call umat50 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,failel,crv,nnpcrv,a(lcma),0,elsiz,idele,
     2   reject)
c
c       scattering of variables
c
 60     if (no_hsvs.gt.0) then
          do j=1,no_hsvs
            hsvs(i,j)=hsv(j)
          enddo
        endif
        if (ihyper(mt).ne.0) then
          do j=1,9
            hsvs(i,ind_defgrad-1+j)=hsv(no_hsvs+j)
          enddo
        endif
c
        if (ifipts(mt).eq.0) failel=.false.
c
        if (failel) then
          sig1(i)=0.
          sig2(i)=0.
          sig3(i)=0.
          sig4(i)=0.
          sig5(i)=0.
          sig6(i)=0.
          if (ifipts(mt).lt.0) hsvs(i,ind_fipts)=1.
          nfipts(i)=nfipts(i)+1
          if (ifipts(mt).ge.0) then
          failur=.true.
          ifail(i)=1
          fail(i)=0.
          elseif (nfipts(i).ge.nint(cm(mx-ifipts(mt)))) then
          failur=.true.
          ifail(i)=1
          fail(i)=0.
          endif
        else
          sig1(i)=sig(1)
          sig2(i)=sig(2)
          sig3(i)=sig(3)
          sig4(i)=sig(4)
          sig5(i)=sig(5)
          sig6(i)=sig(6)
        endif
        epsps(i)=epsp
 90     continue
c
c     vector umat
c
      else
c
c       assume no failed elements
c
        if (ifipts(mt).ge.0) then
        do i=lft,llt
          failels(i)=.false.
        enddo
        else
          do i=lft,llt
            failels(i)=hsvs(i,ind_fipts).ne.0.
          enddo
        endif
c
c       call user defined subroutines here
c       first history variable is the one requested
c
        if(ivumat(mt).eq.0) then
c
        if (mt.eq.41) then
          call umat41v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),
     .     0,elsizv,idelev,reject)
        elseif (mt.eq.42) then
          call umat42v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),
     .     0,elsizv,idelev,reject)
        elseif (mt.eq.43) then
          call umat43v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),
     .     0,elsizv,idelev,reject)
        elseif (mt.eq.44) then
          call umat44v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),
     .     0,elsizv,idelev,reject)
        elseif (mt.eq.45) then
          call umat45v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),
     .     0,elsizv,idelev,0,reject)
        elseif (mt.eq.46) then
          call umat46v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),
     .     thhsv,nthhsv,0,elsizv,idelev,reject)
        elseif (mt.eq.47) then
          call umat47v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),
     .     0,elsizv,idelev,reject)
        elseif (mt.eq.48) then
c         ipt=1
          call umat48v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,ipt,a(n8),a(n9),a(lcma),0,
     .     elsizv,idelev,reject)
        elseif (mt.eq.49) then
c         ipt=1
          call umat49v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,ipt,temps,a(lcma),0,elsizv,
     .     idelev,reject)
        elseif (mt.eq.50) then
          call umat50v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,nnpcrv,
     .     a(lcma),
     .     0,elsizv,idelev,reject)
        elseif (mt.eq.281) then
c         ipt=1
          call mat281h(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
        elseif (mt.eq.282) then
c         ipt=1
          call mat282h(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
c       elseif (mt.eq.283) then
c         call umat283h(cm(mx+1),d1,d2,d3,d4,d5,d6,
c    .     sig1,sig2,sig3,sig4,
c    .     sig5,sig6,epsps,hsvs,
c    .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,
c    .     a(lcma))
        elseif (mt.eq.284) then
c         ipt=1
          call mat284h(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
c       elseif (mt.eq.285) then
c         call umat285h(cm(mx+1),d1,d2,d3,d4,d5,d6,
c    .     sig1,sig2,sig3,sig4,
c    .     sig5,sig6,epsps,hsvs,
c    .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv)
        elseif (mt.eq.286) then
c         ipt=1
          call mat286h(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
        elseif (mt.eq.287) then
c         ipt=1
          call mat287h(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,failels,nlq,crv,
     .     ipt, npc, plc, nnm1,a(lcma))
        endif
c
        endif

        if (ivumat(mt).ne.0) then 
          ndi=3
          nsh=3
          nblock=llt-lft+1
          curveid=ivumat(mt)
          icid=int(curveid)
          if (icid.lt.0) then
            iicid=lcids(-icid)
            k2=npc(iicid)
            k3=npc(iicid+1)
            nprops=nint((k3-k2)/2.0)
          else
            nprops=40+icmadd(mt)
          endif
c          write(*,*)'INDICES DF:',
c     $         if11,if22,if33,if12,if21,if23,if32,if13,if31
c          write(*,*)'INDICES VMAT:',
c     $         ifo11,ifo22,ifo33,ifo12,ifo21,ifo23,ifo32,ifo13,ifo31
          call ext_vecumat(lft,llt,cm,a(lcma),bqs,mt,crv,npc,plc,nnm1,
     $         curveid,rcoor,scoor,tcoor,nconstp,nip,ipt,eltype,nblock,
     $         nprops,ndi,nsh, if11,if22,if33,if12,if21,if23,if32,if13
     $         ,if31,no_hsvs, ifo11,ifo22,ifo33,ifo12,ifo21,ifo23,ifo32
     $         ,ifo13,ifo31,tmp)
        endif
c     
        if (ifipts(mt).ne.0) then
          do i=lft,llt
            if (failels(i)) then
              sig1(i)=0.
              sig2(i)=0.
              sig3(i)=0.
              sig4(i)=0.
              sig5(i)=0.
              sig6(i)=0.
              if (ifipts(mt).lt.0) hsvs(i,ind_fipts)=1.
              nfipts(i)=nfipts(i)+1
              if (ifipts(mt).ge.0) then
                failur=.true.
                ifail(i)=1
                fail(i)=0.
              elseif (nfipts(i).ge.nint(cm(mx-ifipts(mt)))) then
                failur=.true.
                ifail(i)=1
                fail(i)=0.
              endif
            endif
          enddo
        endif
      endif
c
c     energy calculations
c
      do 100 i=lft,llt
      einc(i)=(d1(i)*sig1(i)+d2(i)*sig2(i)+d3(i)*sig3(i)+
     . d4(i)*sig4(i)+d5(i)*sig5(i)+d6(i)*sig6(i))+einc(i)
 100  continue
c
      if (iorien(mt).ne.0) then
c
c       transform stresses to global system
c     
        if (ihyper(mt).ge.0) then
c     if ((ihyper(mt).eq.0) .or. (ihyper(mt).eq.1)) then
          call gblm22 (lft,llt)
        endif
c     
c     if ((ihyper(mt).eq.1) .or. (ihyper(mt).eq.0) ) then
        if (ihyper(mt).gt.0) then
c     
c     deformation gradient stored in global system
c
          if (ifg(lft).eq.1) then
            do i=lft,llt
              hsvs(i,if11)=g11(i)
              hsvs(i,if21)=g21(i)
              hsvs(i,if31)=g31(i)
              hsvs(i,if12)=g12(i)
              hsvs(i,if22)=g22(i)
              hsvs(i,if32)=g32(i)
              hsvs(i,if13)=g13(i)
              hsvs(i,if23)=g23(i)
              hsvs(i,if33)=g33(i)
            enddo
c     
c     deformation gradient transformed back to global system
c
          else
            do i=lft,llt
              f1=z11(i)*hsvs(i,if11)+
     1             z21(i)*hsvs(i,if21)+z31(i)*hsvs(i,if31)
              f2=z12(i)*hsvs(i,if11)+
     1             z22(i)*hsvs(i,if21)+z32(i)*hsvs(i,if31)
              hsvs(i,if31)=z13(i)*hsvs(i,if11)+
     1             z23(i)*hsvs(i,if21)+z33(i)*hsvs(i,if31)
              hsvs(i,if11)=f1
              hsvs(i,if21)=f2
              f1=z11(i)*hsvs(i,if12)+
     1             z21(i)*hsvs(i,if22)+z31(i)*hsvs(i,if32)
              f2=z12(i)*hsvs(i,if12)+
     1             z22(i)*hsvs(i,if22)+z32(i)*hsvs(i,if32)
              hsvs(i,if32)=z13(i)*hsvs(i,if12)+
     1             z23(i)*hsvs(i,if22)+z33(i)*hsvs(i,if32)
              hsvs(i,if12)=f1
              hsvs(i,if22)=f2
              f1=z11(i)*hsvs(i,if13)+
     1             z21(i)*hsvs(i,if23)+z31(i)*hsvs(i,if33)
              f2=z12(i)*hsvs(i,if13)+
     1             z22(i)*hsvs(i,if23)+z32(i)*hsvs(i,if33)
              hsvs(i,if33)=z13(i)*hsvs(i,if13)+
     1             z23(i)*hsvs(i,if23)+z33(i)*hsvs(i,if33)
              hsvs(i,if13)=f1
              hsvs(i,if23)=f2
            enddo
          endif
c     
        elseif (ihyper(mt).eq.-2) then
c     
c     deformation gradient stored in global system
c     
          if (ifg(lft).eq.1) then
            do i=lft,llt
              hsvs(i,if11)=g11(i)
              hsvs(i,if21)=g21(i)
              hsvs(i,if31)=g31(i)
              hsvs(i,if12)=g12(i)
              hsvs(i,if22)=g22(i)
              hsvs(i,if32)=g32(i)
              hsvs(i,if13)=g13(i)
              hsvs(i,if23)=g23(i)
              hsvs(i,if33)=g33(i)
            enddo
c     
c     deformation gradient transformed back to global system
c     
          else
            do i=lft,llt
              f11=z11(i)*hsvs(i,if11)+
     1             z21(i)*hsvs(i,if21)+z31(i)*hsvs(i,if31)
              f21=z12(i)*hsvs(i,if11)+
     1             z22(i)*hsvs(i,if21)+z32(i)*hsvs(i,if31)
              f31=z13(i)*hsvs(i,if11)+
     1             z23(i)*hsvs(i,if21)+z33(i)*hsvs(i,if31)
c     hsvs(i,if11)=f1
c     hsvs(i,if21)=f2
              f12=z11(i)*hsvs(i,if12)+
     1             z21(i)*hsvs(i,if22)+z31(i)*hsvs(i,if32)
              f22=z12(i)*hsvs(i,if12)+
     1             z22(i)*hsvs(i,if22)+z32(i)*hsvs(i,if32)
              f32=z13(i)*hsvs(i,if12)+
     1             z23(i)*hsvs(i,if22)+z33(i)*hsvs(i,if32)
c     hsvs(i,if12)=f1
c     hsvs(i,if22)=f2
              f13=z11(i)*hsvs(i,if13)+
     1             z21(i)*hsvs(i,if23)+z31(i)*hsvs(i,if33)
              f23=z12(i)*hsvs(i,if13)+
     1             z22(i)*hsvs(i,if23)+z32(i)*hsvs(i,if33)
              f33=z13(i)*hsvs(i,if13)+
     1             z23(i)*hsvs(i,if23)+z33(i)*hsvs(i,if33)
c     hsvs(i,if13)=f1
c     hsvs(i,if23)=f2
              hsvs(i,if11)=f11*z11(i)+f12*z21(i)+f13*z31(i)
              hsvs(i,if12)=f11*z12(i)+f12*z22(i)+f13*z32(i)
              hsvs(i,if13)=f11*z13(i)+f12*z23(i)+f13*z33(i)
              hsvs(i,if21)=f21*z11(i)+f22*z21(i)+f23*z31(i)
              hsvs(i,if22)=f21*z12(i)+f22*z22(i)+f23*z32(i)
              hsvs(i,if23)=f21*z13(i)+f22*z23(i)+f23*z33(i)
              hsvs(i,if31)=f31*z11(i)+f32*z21(i)+f33*z31(i)
              hsvs(i,if32)=f31*z12(i)+f32*z22(i)+f33*z32(i)
              hsvs(i,if33)=f31*z13(i)+f32*z23(i)+f33*z33(i)
            enddo
          endif
          
          
c     elseif (ihyper(mt).lt.0) then
        elseif (ihyper(mt).eq.-1) then
c     
c     deformation gradient stored in global system
c     
          if (ifg(lft).eq.1) then
            do i=lft,llt
              hsvs(i,if11)=g11(i)
              hsvs(i,if21)=g21(i)
              hsvs(i,if31)=g31(i)
              hsvs(i,if12)=g12(i)
              hsvs(i,if22)=g22(i)
              hsvs(i,if32)=g32(i)
              hsvs(i,if13)=g13(i)
              hsvs(i,if23)=g23(i)
              hsvs(i,if33)=g33(i)
            enddo
c     
c     deformation gradient transformed back to global system
c     
          else
            do i=lft,llt
              f1=z11(i)*hsvs(i,if11)+
     1             z21(i)*hsvs(i,if12)+z31(i)*hsvs(i,if13)
              f2=z12(i)*hsvs(i,if11)+
     1             z22(i)*hsvs(i,if12)+z32(i)*hsvs(i,if13)
              hsvs(i,if13)=z13(i)*hsvs(i,if11)+
     1             z23(i)*hsvs(i,if12)+z33(i)*hsvs(i,if13)
              hsvs(i,if11)=f1
              hsvs(i,if12)=f2
              f1=z11(i)*hsvs(i,if21)+
     1             z21(i)*hsvs(i,if22)+z31(i)*hsvs(i,if23)
              f2=z12(i)*hsvs(i,if21)+
     1             z22(i)*hsvs(i,if22)+z32(i)*hsvs(i,if23)
              hsvs(i,if23)=z13(i)*hsvs(i,if21)+
     1             z23(i)*hsvs(i,if22)+z33(i)*hsvs(i,if23)
              hsvs(i,if21)=f1
              hsvs(i,if22)=f2
              f1=z11(i)*hsvs(i,if31)+
     1             z21(i)*hsvs(i,if32)+z31(i)*hsvs(i,if33)
              f2=z12(i)*hsvs(i,if31)+
     1             z22(i)*hsvs(i,if32)+z32(i)*hsvs(i,if33)
              hsvs(i,if33)=z13(i)*hsvs(i,if31)+
     1             z23(i)*hsvs(i,if32)+z33(i)*hsvs(i,if33)
              hsvs(i,if31)=f1
              hsvs(i,if32)=f2
            enddo
          endif
        endif
      endif
c     
c     flag if material model rejected implicit iterate
      if (isolvr(1).ne.0.and.reject) isolvr(78)=5
c
      return
c     900  format(/
c    1 ' *** Error User defined material in part',i12,
c    2 ' exceeds the current limit'/
c    3 10x,'of history variables for solids.')
      end
      subroutine usr_temps (tnew,fval,itemp,ctmp,lft,llt,ix1,ix2,
     . ix3,ix4,ix5,ix6,ix7,ix8,num_nods,itopaz,nnm1,rcoor,scoor,tcoor)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      common/blk03/sigma,tolb,nummat,numnpt,numel,numsh4,iunits,iband,
     1 nbs,nsl,nslvt,nmsrt,numelst,nprof,irtyp,itmaxb,numeltt,iphase,
     2 igenl,igenm,igene,igend,isotr,iem,isoln,igeom,mcut,numsh12,
     3 ndtot,nsl_th,lenhsv,numel2,numel4,numel6,numel8,numel10,numel20
      dimension tnew(*),fval(*),ctmp(*),ix1(*),ix2(*),ix3(*),ix4(*),
     .          ix5(*),ix6(*),ix7(*),ix8(*),itopaz(*)
c
c    temperatures at element center
c
      if (itemp.gt.0) then
        tempr=fval(itemp)
        do i=lft,llt
          ctmp(i)=tempr
        enddo
c
      elseif (itemp.lt.0 .or. itopaz(1).eq.999) then
        if (num_nods.eq.8) then
          do i=lft,llt
            ctmp(i) = .125*(tnew(ix1(i))+tnew(ix2(i))+tnew(ix3(i))
     .                      +tnew(ix4(i))+tnew(ix5(i))+tnew(ix6(i))
     .                      +tnew(ix7(i))+tnew(ix8(i)))
          enddo
c
        elseif (num_nods.eq.4) then
          if (numsh12.ne.0) then
           call get_tmp_tpz (ctmp,lft,llt,nnm1,rcoor,scoor,tcoor)
          else
            do i=lft,llt
            ctmp(i)=.25*(tnew(ix1(i))+tnew(ix2(i))
     .                  +tnew(ix3(i))+tnew(ix4(i)))
            enddo
          endif
c
        else
          do i=lft,llt
            ctmp(i)=.50*(tnew(ix1(i))+tnew(ix2(i)))
          enddo
        endif
      endif
c
      return
      end
      subroutine urtanh(cm,mt,crv,nnpcrv,lft,llt,nconstp,nnm1,eltype)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     assume that F and all coordinate transformations have been
c     performed by this time, and everything is already loaded in
c     common blocks
c
      include 'nlqparm'
      include 'nhisparm.inc'
#include "bk07.inc"
#include "bk13.inc"
      include 'bk19.inc'
      include 'memaia.inc'
      include 'umatss.inc'
C_TASKCOMMON (aux2loc)
C_TASKCOMMON (aux14loc)
C_TASKCOMMON (aux18loc)
C_TASKCOMMON (aux33loc)
C_TASKCOMMON (aux35loc)
C_TASKCOMMON (aux40loc)
C_TASKCOMMON (aux41loc)
C_TASKCOMMON (bk36loc)
C_TASKCOMMON (bux13loc)
C_TASKCOMMON (soundloc)
C_TASKCOMMON (subtssloc)
C_TASKCOMMON (vect8loc)
      common/aux14loc/sig1(nlq),sig2(nlq),sig3(nlq),sig4(nlq),sig5(nlq),
     & sig6(nlq),epsps(nlq),hsvs(nlq,NHISVAR)
      common/aux18loc/dd(nlq),def(nlq),ddq(nlq)
      common/aux33loc/
     1 ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq),ixs(nlq,4),mxt(nlq)
      common/aux40loc/
     1 a11(nlq),a12(nlq),a13(nlq),a21(nlq),a22(nlq),a23(nlq),
     2 a31(nlq),a32(nlq),a33(nlq),z11(nlq),z12(nlq),z13(nlq),
     3 z21(nlq),z22(nlq),z23(nlq),z31(nlq),z32(nlq),z33(nlq)
      common/aux35loc/rhoa(nlq),cb(nlq),davg(nlq),p(nlq)
      common/aux41loc/qq1(nlq),cbb(nlq),aj1(nlq)
      common/bk36loc/index
      common/bux13loc/ eps(6),sig(6),hsv(NHISVAR),temps(nlq),es(6,6)
      common/soundloc/ss(nlq),sndsp(nlq),diagm(nlq),sarea(nlq),
     . dxl(nlq)
      common/subtssloc/dt1siz(nlq)
      common/vect8loc/dsave(nlq,6,6)
c
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,nusa
      common/bk01/itherm,itemp,ntmp0,ntmp1,itempan,itempdr,itmpe(44)
      common/bk28/summss,xke,xpe,tt,xte0,erodeke,erodeie
      common/aux2loc/d1(nlq),d2(nlq),d3(nlq),d4(nlq),d5(nlq),d6(nlq),
     1 wzzdt(nlq),wyydt(nlq),wxxdt(nlq),einc(nlq)
      common/bk26/begtim,nintcy
      common/numcpu/ncpu,ncpua,ncpub,lenvec(8)
      common/thropt/itopaz(101)
      common/tsbsis/h(8,5,5),pr(8,5,5),ps(8,5,5),pt(8,5,5),ipt,
     1 nip,wgts(10,10),szeta(10,10)
c
      dimension cm(*),crv(lq1,2,*),nconstp(*)
      integer nnpcrv(*)
      character*(*) eltype
c     dimension eps(6),sig(6),hsv(NHISVAR),temps(nlq),es(6,6)
      dimension sigmaold(nlq,6)
      dimension dfgold(nlq,9)
      logical failel,failels(nlq)
c     data temps/nlq*0.0/
c
      capa=0.0
      mx=48*(mxt(lft)-1)
c
c     Determining number of requested history variables
c
      no_hsvs=nconstp(mxt(lft))
      if (lenvec(5) .eq.1) no_hsvs=no_hsvs-6
      if (iorien(mt).ne.0) no_hsvs=no_hsvs-6
      if (ihyper(mt).ne.0) no_hsvs=no_hsvs-9
      if (iorien(mt).ne.0.and.ihyper(mt).lt.0) no_hsvs=no_hsvs-6
      if (ifipts(mt).lt.0) no_hsvs=no_hsvs-1
c
c     initialization of indeces
c
      if (iorien(mt).ne.0) then
        if (ihyper(mt).ne.0) then
          ind_defgrad=1+no_hsvs
          ind_ortho=ind_defgrad+9
        else
          ind_ortho=1+no_hsvs
        endif
        if (ifipts(mt).lt.0) then
          ind_fipts=ind_ortho
          ind_ortho=ind_ortho+1
        endif
        ind_q11=ind_ortho
        ind_q12=ind_ortho+1
        ind_q13=ind_ortho+2
        ind_q31=ind_ortho+3
        ind_q32=ind_ortho+4
        ind_q33=ind_ortho+5
        ind_last=ind_q33
        if (ihyper(mt).lt.0) then
          ind_qstore=ind_last+1
          ind_last=ind_last+6
        endif
      else
        if (ihyper(mt).ne.0) then
          ind_defgrad=1+no_hsvs
          ind_last=ind_defgrad+8
        else
          ind_last=no_hsvs
        endif
        if (ifipts(mt).lt.0) then
          ind_fipts=ind_last+1
          ind_last=ind_last+1
        endif
      endif
c
      if (iorien(mt).ne.0) then
c
c       store global stresses
c
        if (ihyper(mt).ge.0) then
          do i=lft,llt
            sigmaold(i,1)=sig1(i)
            sigmaold(i,2)=sig2(i)
            sigmaold(i,3)=sig3(i)
            sigmaold(i,4)=sig4(i)
            sigmaold(i,5)=sig5(i)
            sigmaold(i,6)=sig6(i)
          enddo
c
c         transform stresses to local system
c
          call ldsm22(lft,llt)
        endif
c
c       storage of deformation gradient
c
        if (ihyper(mt).gt.0) then
          if11=ind_defgrad
          if21=ind_defgrad+1
          if31=ind_defgrad+2
          if12=ind_defgrad+3
          if22=ind_defgrad+4
          if32=ind_defgrad+5
          if13=ind_defgrad+6
          if23=ind_defgrad+7
          if33=ind_defgrad+8
c
c         store deformation gradient in global system
c
          do i=lft,llt
            dfgold(i,1)=hsvs(i,if11)
            dfgold(i,2)=hsvs(i,if21)
            dfgold(i,3)=hsvs(i,if31)
            dfgold(i,4)=hsvs(i,if12)
            dfgold(i,5)=hsvs(i,if22)
            dfgold(i,6)=hsvs(i,if32)
            dfgold(i,7)=hsvs(i,if13)
            dfgold(i,8)=hsvs(i,if23)
            dfgold(i,9)=hsvs(i,if33)
          enddo
c
c         transform deformation gradient to local system
c         as F_bar=Q^T*F
c
          do i=lft,llt
            f1=z11(i)*hsvs(i,if11)+
     1         z12(i)*hsvs(i,if21)+z13(i)*hsvs(i,if31)
            f2=z21(i)*hsvs(i,if11)+
     1         z22(i)*hsvs(i,if21)+z23(i)*hsvs(i,if31)
            hsvs(i,if31)=z31(i)*hsvs(i,if11)+
     1         z32(i)*hsvs(i,if21)+z33(i)*hsvs(i,if31)
            hsvs(i,if11)=f1
            hsvs(i,if21)=f2
            f1=z11(i)*hsvs(i,if12)+
     1         z12(i)*hsvs(i,if22)+z13(i)*hsvs(i,if32)
            f2=z21(i)*hsvs(i,if12)+
     1         z22(i)*hsvs(i,if22)+z23(i)*hsvs(i,if32)
            hsvs(i,if32)=z31(i)*hsvs(i,if12)+
     1         z32(i)*hsvs(i,if22)+z33(i)*hsvs(i,if32)
            hsvs(i,if12)=f1
            hsvs(i,if22)=f2
            f1=z11(i)*hsvs(i,if13)+
     1         z12(i)*hsvs(i,if23)+z13(i)*hsvs(i,if33)
            f2=z21(i)*hsvs(i,if13)+
     1         z22(i)*hsvs(i,if23)+z23(i)*hsvs(i,if33)
            hsvs(i,if33)=z31(i)*hsvs(i,if13)+
     1         z32(i)*hsvs(i,if23)+z33(i)*hsvs(i,if33)
            hsvs(i,if13)=f1
            hsvs(i,if23)=f2
          enddo
c
c       storage of deformation gradient
c
        else if (ihyper(mt).lt.0) then
          if11=ind_defgrad
          if21=ind_defgrad+1
          if31=ind_defgrad+2
          if12=ind_defgrad+3
          if22=ind_defgrad+4
          if32=ind_defgrad+5
          if13=ind_defgrad+6
          if23=ind_defgrad+7
          if33=ind_defgrad+8
c
c         store deformation gradient in global system
c
          do i=lft,llt
            dfgold(i,1)=hsvs(i,if11)
            dfgold(i,2)=hsvs(i,if21)
            dfgold(i,3)=hsvs(i,if31)
            dfgold(i,4)=hsvs(i,if12)
            dfgold(i,5)=hsvs(i,if22)
            dfgold(i,6)=hsvs(i,if32)
            dfgold(i,7)=hsvs(i,if13)
            dfgold(i,8)=hsvs(i,if23)
            dfgold(i,9)=hsvs(i,if33)
          enddo
c
c         transform deformation gradient to local system
c         as F_bar=F*Q
c
          do i=lft,llt
            f1=z11(i)*hsvs(i,if11)+
     1         z12(i)*hsvs(i,if12)+z13(i)*hsvs(i,if13)
            f2=z21(i)*hsvs(i,if11)+
     1         z22(i)*hsvs(i,if12)+z23(i)*hsvs(i,if13)
            hsvs(i,if13)=z31(i)*hsvs(i,if11)+
     1         z32(i)*hsvs(i,if12)+z33(i)*hsvs(i,if13)
            hsvs(i,if11)=f1
            hsvs(i,if12)=f2
            f1=z11(i)*hsvs(i,if21)+
     1         z12(i)*hsvs(i,if22)+z13(i)*hsvs(i,if23)
            f2=z21(i)*hsvs(i,if21)+
     1         z22(i)*hsvs(i,if22)+z23(i)*hsvs(i,if23)
            hsvs(i,if23)=z31(i)*hsvs(i,if21)+
     1         z32(i)*hsvs(i,if22)+z33(i)*hsvs(i,if23)
            hsvs(i,if21)=f1
            hsvs(i,if22)=f2
            f1=z11(i)*hsvs(i,if31)+
     1         z12(i)*hsvs(i,if32)+z13(i)*hsvs(i,if33)
            f2=z21(i)*hsvs(i,if31)+
     1         z22(i)*hsvs(i,if32)+z23(i)*hsvs(i,if33)
            hsvs(i,if33)=z31(i)*hsvs(i,if31)+
     1         z32(i)*hsvs(i,if32)+z33(i)*hsvs(i,if33)
            hsvs(i,if31)=f1
            hsvs(i,if32)=f2
          enddo
        endif
      endif
c
c     if ishrmp(mt).lt.0 then the flag for temperature is active
c     temperatures are stored in the temps array.
c
      num_nods=8
      if (eltype.eq.'shl_t') num_nods=4
      if (eltype.eq.'sld2d'.or.eltype.eq.'sldax') num_nods=4
      if (ishrmp(mt).lt.0) then
        call usr_temps (a(ntmp0+1),a(n19),itemp,temps,lft,llt,ix1,ix2,
     .   ix3,ix4,ixs(1,1),ixs(1,2),ixs(1,3),ixs(1,4),num_nods,
     .   itopaz,0,0.,0.,0.)
      else
        do i=lft,llt
          temps(i)=0.0
        enddo
      endif
c
c     memory pointer for extra memory data
      lcma=ncma+nmmat+ia(ncma+mxt(lft)-1)-1
c
c     vector flag
c
      ivect=ivectr(mt)
c
      mte40=mt-40
c
c     scalar utan
c
      if (ivect.eq.0) then
        do 90 i=lft,llt
c
c       gathering of variables
c
        index =i
        dt1   =dt1siz(i)
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
        epsp=epsps(i)
        temper=temps(i)
        failel=.false.
        if (ifipts(mt).lt.0) failel=hsvs(i,ind_fipts).ne.0.
        if (no_hsvs.gt.0) then
          do 30 j=1,no_hsvs
          hsv(j)=hsvs(i,j)
 30       continue
        endif
        if (ihyper(mt).ne.0) then
          do 40 j=1,9
          hsv(no_hsvs+j)=hsvs(i,ind_defgrad-1+j)
 40       continue
        endif
c
c       setting material tangent modulus to = 0
c
        do k=1,6
          do j=1,6
            es(j,k)=0.
          enddo
        enddo
c
c       call user developed subroutines here
c
        go to (41,42,43,44,45,46,47,48,49,50), mte40
   41   call utan41 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,nnpcrv,failel,a(lcma),0)
        go to 60
   42   call utan42 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,nnpcrv,failel,a(lcma),0)
        go to 60
   43   call utan43 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,nnpcrv,failel,a(lcma),0)
        go to 60
   44   call utan44 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,nnpcrv,failel,a(lcma),0)
        go to 60
   45   call utan45 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,nnpcrv,failel,a(lcma),0)
        go to 60
   46   call utan46 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,nnpcrv,failel,a(lcma),0)
        go to 60
   47   call utan47 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,nnpcrv,failel,a(lcma),0)
        go to 60
   48   call utan48 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,nnpcrv,failel,a(lcma),0)
        go to 60
   49   call utan49 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,nnpcrv,failel,a(lcma),0)
        go to 60
 50     call utan50 (cm(mx+1),eps,sig,epsp,hsv,dt1,
     1   capa,eltype,tt,temper,es,crv,nnpcrv,failel,a(lcma),0)
c
c       scattering tangent modulus
c
 60     do k=1,6
          do j=1,6
           dsave(i,j,k)=es(j,k)
          enddo
        enddo
 90    continue
c
c
c     vector umat
c
      else
c
c       setting material tangent modulus to = 0
c
        do k=1,6
          do j=1,6
            do i=lft,llt
              dsave(i,j,k)=0.
            enddo
          enddo
        enddo
c
        if (ifipts(mt).ge.0) then
          do i=lft,llt
            failels(i)=.false.          
          enddo
        else
          do i=lft,llt
            failels(i)=hsvs(i,ind_fipts).ne.0.
          enddo
        endif
c
c       call user developed routines here
c       first history variable is the one requested
c
        if (mt.eq.41) then
          call utan41v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,nnpcrv,
     .     failels,
     .     a(lcma),0)
        elseif (mt.eq.42) then
          call utan42v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,nnpcrv,
     .     failels,
     .     a(lcma),0)
        elseif (mt.eq.43) then
          call utan43v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,nnpcrv,
     .     failels,
     .     a(lcma),0)
        elseif (mt.eq.44) then
          call utan44v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,nnpcrv,
     .     failels,
     .     a(lcma),0)
        elseif (mt.eq.45) then
          call utan45v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,nnpcrv,
     .     failels,
     .     a(lcma),0)
        elseif (mt.eq.46) then
          call utan46v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,nnpcrv,
     .     failels,
     .     a(lcma),0)
        elseif (mt.eq.47) then
          call utan47v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,nnpcrv,
     .     failels,
     .     a(lcma),0)
        elseif (mt.eq.48) then
          call utan48v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,nnpcrv,
     .     failels,
     .     a(lcma),0)
        elseif (mt.eq.49) then
          call utan49v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,nnpcrv,
     .     failels,
     .     a(lcma),0)
        elseif (mt.eq.50) then
          call utan50v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,nnpcrv,
     .     failels,
     .     a(lcma),0)
        elseif (mt.eq.286) then
          call utan286v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,failels,
     .     a(lcma))
        elseif (mt.eq.287) then
          call utan287v(cm(mx+1),d1,d2,d3,d4,d5,d6,
     .     sig1,sig2,sig3,sig4,
     .     sig5,sig6,epsps,hsvs,
     .     lft,llt,dt1siz,capa,eltype,tt,temps,dsave,nlq,crv,failels,
     .     a(lcma))
        endif
      endif
c
      if (iorien(mt).ne.0) then
c
c       transform tangent stiffness to global system
c
        if (ihyper(mt).ge.0) then
          call transform44(z11,z12,z13,z21,z22,z23,z31,z32,z33,
     1     lft,llt)
c
c         restore stresses in global system
c
          do i=lft,llt
            sig1(i)=sigmaold(i,1)
            sig2(i)=sigmaold(i,2)
            sig3(i)=sigmaold(i,3)
            sig4(i)=sigmaold(i,4)
            sig5(i)=sigmaold(i,5)
            sig6(i)=sigmaold(i,6)
          enddo
        endif
c
c       restore deformation gradient in global system
c
        if (ihyper(mt).ne.0) then
          do i=lft,llt
            hsvs(i,if11)=dfgold(i,1)
            hsvs(i,if21)=dfgold(i,2)
            hsvs(i,if31)=dfgold(i,3)
            hsvs(i,if12)=dfgold(i,4)
            hsvs(i,if22)=dfgold(i,5)
            hsvs(i,if32)=dfgold(i,6)
            hsvs(i,if13)=dfgold(i,7)
            hsvs(i,if23)=dfgold(i,8)
            hsvs(i,if33)=dfgold(i,9)
          enddo
        endif
      endif
c
      return
      end
      subroutine urtans(cm,mt,capa,lft,llt,nconstp,nnm1)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
#include "bk07.inc"
#include "bk13.inc"
      include 'bk19.inc'
      include 'memaia.inc'
      include 'umatss.inc'
#include "dynmem.inc"
#include "ez_mem.inc"
c
C_TASKCOMMON (aux2loc)
C_TASKCOMMON (aux8cloc)
C_TASKCOMMON (aux14loc)
C_TASKCOMMON (aux19loc)
C_TASKCOMMON (aux33loc)
C_TASKCOMMON (bk36loc)
C_TASKCOMMON (bux13loc)
C_TASKCOMMON (failuloc)
C_TASKCOMMON (hourgloc)
C_TASKCOMMON (shloptloc)
C_TASKCOMMON (soundloc)
C_TASKCOMMON (subtssloc)
C_TASKCOMMON (vect8loc)
C_TASKCOMMON (frozloc)
      common/aux2loc/e1(nlq),e2(nlq),d3(nlq),e4(nlq),e5(nlq),e6(nlq),
     1 wzzdt(nlq),wyydt(nlq),wxxdt(nlq),einc(nlq)
      common /aux8cloc/
     & d1(nlq),d2(nlq),d4(nlq),d5(nlq),d6(nlq),
     & dfg11(nlq),dfg21(nlq),dfg31(nlq),dfg12(nlq),dfg22(nlq),
     & dfg32(nlq),dfg13(nlq),dfg23(nlq),dfg33(nlq),stg5(nlq),
     & stg6(nlq),a11(nlq),a12(nlq),a21(nlq),a22(nlq)
      common/aux14loc/
     & sig1(nlq),sig2(nlq),sig3(nlq),sig4(nlq),
     & sig5(nlq),sig6(nlq),epsps(nlq),hsvs(nlq,NHISVAR)
      common/aux19loc/
     1 sign0(nlq),sign1(nlq),sign2(nlq),sign3(nlq),sign4(nlq),
     2 sign5(nlq),sign6(nlq)
      common/aux33loc/
     1 ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq),ixs(nlq,4),mxt(nlq)
      common/bk36loc/index
      common/bux13loc/ eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      common/failuloc/sieu(nlq),fail(nlq)
      common/frozloc/yhat(nlq,3,4),qloc(nlq,3,3),qmat(nlq,3,3)
      common/hourgloc/ym(nlq),gm(nlq),ifsv(nlq)
      common/shlopt/istrn,istupd,ibelyts,miter
      common/shloptloc/ibelyt
      common/soundloc/ss(nlq),sndsp(nlq),diagm(nlq),sarea(nlq),
     . dxl(nlq)
      common/vect8loc/dsave(nlq,6,6)
c
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,nusa
      common/bk01/itherm,itemp,ntmp0,ntmp1,itempan,itempdr,itmpe(44)
      common/bk26/begtim,nintcy
      common/bk28/summss,xke,xpe,tt,xte0,erodeke,erodeie
      common/numcpu/ncpu,ncpua,ncpub,lenvec(8)
      common/slcntr/islcnt(21)
      common/subtssloc/dt1siz(nlq)
      common/thropt/itopaz(101)
c
      dimension cm(*),nconstp(*)
c     dimension eps(6),sig(6),hsv(NHISVAR),temps(nlq)
      dimension es(6,6)
      dimension sigmaold(nlq,6)
      logical failel,failels(nlq)
      dimension qmats(3,3)
c     data temps/nlq*0.0/
c
c     Note that no temperatures are calculated here
c
      mx=48*(mxt(lft)-1)
c
c     Determining number of requested history variables
c
      no_hsvs=nconstp(mxt(lft))
      if (lenvec(5) .eq.1) no_hsvs=no_hsvs-6
      if (iorien(mt).ne.0) no_hsvs=no_hsvs-2
      if (ihyper(mt).ne.0) no_hsvs=no_hsvs-9
      if (ifipts(mt).lt.0) no_hsvs=no_hsvs-1
c
c     initializing indeces
c
      if (iorien(mt).ne.0) then
        if (ihyper(mt).ne.0) then
          ind_defgrad=1+no_hsvs
          ind_ortho=ind_defgrad+9
        else
          ind_ortho=1+no_hsvs
        endif
        if (ifipts(mt).lt.0) then
          ind_fipts=ind_ortho
          ind_ortho=ind_ortho+1
        endif
        ind_q1=ind_ortho
        ind_q2=ind_ortho+1
        ind_last=ind_q2
      else
        if (ihyper(mt).ne.0) then
          ind_defgrad=1+no_hsvs
          ind_last=ind_defgrad+8
        else
          ind_last=no_hsvs
        endif
        if (ifipts(mt).lt.0) then
          ind_fipts=ind_last+1
          ind_last=ind_last+1
        endif
      endif
c
      if (iorien(mt).ne.0) then
c
c       store element local stresses
c
        if (ihyper(mt).ge.0) then
          do i=lft,llt
            sigmaold(i,1)=sig1(i)
            sigmaold(i,2)=sig2(i)
            sigmaold(i,3)=sig3(i)
            sigmaold(i,4)=sig4(i)
            sigmaold(i,5)=sig5(i)
            sigmaold(i,6)=sig6(i)
          enddo
c
c         transform stresses to material system
c
          do 20 i=lft,llt
          stg5(i)=sig5(i)
          stg6(i)=sig6(i)
          a11(i) =hsvs(i,ind_q1)*sig1(i)-hsvs(i,ind_q2)*sig4(i)
          a12(i) =hsvs(i,ind_q2)*sig1(i)+hsvs(i,ind_q1)*sig4(i)
          a21(i) =hsvs(i,ind_q1)*sig4(i)-hsvs(i,ind_q2)*sig2(i)
          a22(i) =hsvs(i,ind_q2)*sig4(i)+hsvs(i,ind_q1)*sig2(i)
          sig1(i)=hsvs(i,ind_q1)*a11(i)-hsvs(i,ind_q2)*a21(i)
          sig2(i)=hsvs(i,ind_q2)*a12(i)+hsvs(i,ind_q1)*a22(i)
          sig4(i)=hsvs(i,ind_q1)*a12(i)-hsvs(i,ind_q2)*a22(i)
          sig5(i)=hsvs(i,ind_q2)*stg6(i)+hsvs(i,ind_q1)*stg5(i)
          sig6(i)=hsvs(i,ind_q1)*stg6(i)-hsvs(i,ind_q2)*stg5(i)
 20       continue
        endif
c
c       storage of deformation gradient
c
        if (ihyper(mt).gt.0) then
          if11=ind_defgrad
          if21=ind_defgrad+1
          if31=ind_defgrad+2
          if12=ind_defgrad+3
          if22=ind_defgrad+4
          if32=ind_defgrad+5
          if13=ind_defgrad+6
          if23=ind_defgrad+7
          if33=ind_defgrad+8
c
c         transform deformation gradient to material system
c         as F_bar=Q^T*F
c
          do i=lft,llt
c
c           store deformation gradient in element system
c
            dfg11(i)=hsvs(i,if11)
            dfg21(i)=hsvs(i,if21)
            dfg12(i)=hsvs(i,if12)
            dfg22(i)=hsvs(i,if22)
            dfg13(i)=hsvs(i,if13)
            dfg23(i)=hsvs(i,if23)
c
            ftemp=hsvs(i,ind_q1)*hsvs(i,if11)-
     1            hsvs(i,ind_q2)*hsvs(i,if21)
            hsvs(i,if21)=hsvs(i,ind_q2)*hsvs(i,if11)+
     1            hsvs(i,ind_q1)*hsvs(i,if21)
            hsvs(i,if11)=ftemp
            ftemp=hsvs(i,ind_q1)*hsvs(i,if12)-
     1            hsvs(i,ind_q2)*hsvs(i,if22)
            hsvs(i,if22)=hsvs(i,ind_q2)*hsvs(i,if12)+
     1            hsvs(i,ind_q1)*hsvs(i,if22)
            hsvs(i,if12)=ftemp
            ftemp=hsvs(i,ind_q1)*hsvs(i,if13)-
     1            hsvs(i,ind_q2)*hsvs(i,if23)
            hsvs(i,if23)=hsvs(i,ind_q2)*hsvs(i,if13)+
     1            hsvs(i,ind_q1)*hsvs(i,if23)
            hsvs(i,if13)=ftemp
          enddo
c
c       storage of deformation gradient
c
        else if (ihyper(mt).lt.0) then
          if11=ind_defgrad
          if21=ind_defgrad+1
          if31=ind_defgrad+2
          if12=ind_defgrad+3
          if22=ind_defgrad+4
          if32=ind_defgrad+5
          if13=ind_defgrad+6
          if23=ind_defgrad+7
          if33=ind_defgrad+8
c
c         transform deformation gradient to material system
c         as F_bar=F*Q
c
          do i=lft,llt
c
c           store deformation gradient in element system
c
            dfg11(i)=hsvs(i,if11)
            dfg21(i)=hsvs(i,if21)
            dfg31(i)=hsvs(i,if31)
            dfg12(i)=hsvs(i,if12)
            dfg22(i)=hsvs(i,if22)
            dfg32(i)=hsvs(i,if32)
c
            ftemp=hsvs(i,ind_q1)*hsvs(i,if11)-
     1            hsvs(i,ind_q2)*hsvs(i,if12)
            hsvs(i,if12)=hsvs(i,ind_q2)*hsvs(i,if11)+
     1            hsvs(i,ind_q1)*hsvs(i,if12)
            hsvs(i,if11)=ftemp
            ftemp=hsvs(i,ind_q1)*hsvs(i,if21)-
     1            hsvs(i,ind_q2)*hsvs(i,if22)
            hsvs(i,if22)=hsvs(i,ind_q2)*hsvs(i,if21)+
     1            hsvs(i,ind_q1)*hsvs(i,if22)
            hsvs(i,if21)=ftemp
            ftemp=hsvs(i,ind_q1)*hsvs(i,if31)-
     1            hsvs(i,ind_q2)*hsvs(i,if32)
            hsvs(i,if32)=hsvs(i,ind_q2)*hsvs(i,if31)+
     1            hsvs(i,ind_q1)*hsvs(i,if32)
            hsvs(i,if31)=ftemp
          enddo
        endif
      endif
c
c     memory pointer for extra memory data
      lcma=ncma+nmmat+ia(ncma+mxt(lft)-1)-1
c
c     vector parameter
c
      ivect=ivectr(mt)
c
      mte40=mt-40
c
      if (ivect.eq.0) then
c
c       scalar utan
c
        do 120 i=lft,llt
c
c       gathering of variables
c
        index=i
        dt1   =dt1siz(i)
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
        epsp=epsps(i)
        temper=temps(i)
        failel=.false.
        if (ifipts(mt).lt.0) failel=hsvs(i,ind_fipts).ne.0.
        if (no_hsvs.gt.0) then
          do j=1,no_hsvs
            hsv(j)=hsvs(i,j)
          enddo
        endif
        if (ihyper(mt).ne.0) then
          do j=1,9
            hsv(no_hsvs+j)=hsvs(i,ind_defgrad-1+j)
          enddo
          if (iabs(ihyper(mt)).eq.3) then
             do k=1,3
                do j=1,3
                   qmats(j,k)=qmat(i,j,k)
                enddo
             enddo
          endif
        endif
c
c       setting material tangent modulus to = 0
c
        do k=1,6
          do j=1,6
            es(j,k)=0.
          enddo
        enddo
c
c       call user developed subroutines here
c
        go to (41,42,43,44,45,46,47,48,49,50), mte40
 41     call utan41 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),i_mem(mezp(34)),failel,a(lcma),qmats)
        go to 60
 42     call utan42 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),i_mem(mezp(34)),failel,a(lcma),qmats)
        go to 60
 43     call utan43 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),i_mem(mezp(34)),failel,a(lcma),qmats)
        go to 60
 44     call utan44 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),i_mem(mezp(34)),failel,a(lcma),qmats)
        go to 60
 45     call utan45 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),i_mem(mezp(34)),failel,a(lcma),qmats)
        go to 60
 46     call utan46 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),i_mem(mezp(34)),failel,a(lcma),qmats)
        go to 60
 47     call utan47 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),i_mem(mezp(34)),failel,a(lcma),qmats)
        go to 60
 48     call utan48 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),i_mem(mezp(34)),failel,a(lcma),qmats)
        go to 60
 49     call utan49 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),i_mem(mezp(34)),failel,a(lcma),qmats)
        go to 60
 50     call utan50 (cm(mx+1),eps,sig,epsp,hsv,dt1,capa,'shell',tt,
     .   temper,es,a(islcnt(16)),i_mem(mezp(34)),failel,a(lcma),qmats)
c
c       scattering material tangent modulus
c
 60     do k=1,6
          do j=1,6
           dsave(i,j,k)=es(j,k)
         enddo
        enddo
 120    continue
c
c     vector utan
c
      else
c
c       setting material tangent modulus to = 0
c
        do k=1,6
          do j=1,6
            do i=lft,llt
              dsave(i,j,k)=0.
            enddo
          enddo
        enddo
c
        if (ifipts(mt).ge.0) then
          do i=lft,llt
            failels(i)=.false.          
          enddo
        else
          do i=lft,llt
            failels(i)=hsvs(i,ind_fipts).ne.0.
          enddo
        endif
c
c       call user developed routines here
c       first history variable is the one requested
c
        if (mt.eq.41) then
          call utan41v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),
     .     i_mem(mezp(34)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.42) then
          call utan42v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),
     .     i_mem(mezp(34)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.43) then
          call utan43v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),
     .     i_mem(mezp(34)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.44) then
          call utan44v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),
     .     i_mem(mezp(34)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.45) then
          call utan45v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),
     .     i_mem(mezp(34)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.46) then
          call utan46v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),
     .     i_mem(mezp(34)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.47) then
          call utan47v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),
     .     i_mem(mezp(34)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.48) then
          call utan48v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),
     .     i_mem(mezp(34)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.49) then
          call utan49v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),
     .     i_mem(mezp(34)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.50) then
          call utan50v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),
     .     i_mem(mezp(34)),failels,
     .     a(lcma),qmat)
        elseif (mt.eq.286) then
          call utan286v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),failels,
     .     a(lcma))
        elseif (mt.eq.287) then
          call utan287v(cm(mx+1),d1,d2,d3,d4,d5,d6,sig1,sig2,
     .     sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,
     .     dt1siz,capa,'shell',tt,temps,dsave,nlq,a(islcnt(16)),failels,
     .     a(lcma))
        endif
      endif
c
      if (iorien(mt).ne.0) then
c
c       transform tangent stiffness to element system
c
        if (ihyper(mt).ge.0) then
          call transform44s(hsvs(1,ind_q1),hsvs(1,ind_q2),lft,llt)
c
c         restore element local stresses
c
          do i=lft,llt
            sig1(i)=sigmaold(i,1)
            sig2(i)=sigmaold(i,2)
            sig3(i)=sigmaold(i,3)
            sig4(i)=sigmaold(i,4)
            sig5(i)=sigmaold(i,5)
            sig6(i)=sigmaold(i,6)
          enddo
        endif
c
c       restore deformation gradient in element system
c
        if (ihyper(mt).ne.0) then
          do i=lft,llt
            hsvs(i,if11)=dfg11(i)
            hsvs(i,if21)=dfg21(i)
            hsvs(i,if12)=dfg12(i)
            hsvs(i,if22)=dfg22(i)
            hsvs(i,if13)=dfg13(i)
            hsvs(i,if23)=dfg23(i)
          enddo
        endif
      endif
c
      return
      end
      subroutine utan50v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,nnpcrv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      logical failels(*)
      character*5 etype
c
      return
      end
      subroutine utan49v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,nnpcrv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      logical failels(*)
      character*5 etype
c
      return
      end
      subroutine utan48v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,nnpcrv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      logical failels(*)
      character*5 etype
c
      return
      end
      subroutine utan47v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,nnpcrv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      logical failels(*)
      character*5 etype
c
      return
      end
      subroutine utan46v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,nnpcrv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      logical failels(*)
      character*5 etype
c
      return
      end
      subroutine utan42v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,nnpcrv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      logical failels(*)
      character*5 etype
c
      return
      end
      subroutine utan44v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,nnpcrv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      logical failels(*)
      character*5 etype
c
      return
      end
      subroutine utan43v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,nnpcrv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      logical failels(*)
      character*5 etype
c
      return
      end
      subroutine utan45v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,nnpcrv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      logical failels(*)
      character*5 etype
c
      dimension sig(6),eps(6),hsv(NHISVAR),es(6,6)
      logical failel
c
      do i=lft,llt
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
c
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
c
        hsv(1)=hsvs(i,1)
        hsv(2)=hsvs(i,2)
        hsv(3)=hsvs(i,3)
        hsv(4)=hsvs(i,4)
        hsv(5)=hsvs(i,5)
        hsv(6)=hsvs(i,6)
        hsv(7)=hsvs(i,7)
        hsv(8)=hsvs(i,8)
        hsv(9)=hsvs(i,9)
c
        do k=1,6
          do j=1,6
            es(j,k)=0.
          enddo
        enddo
c
        failel=failels(i)
c
        call utan45(cm,eps,sig,epsps(i),hsv,dt1siz(i),capa,etype,tt,
     1   temps(i),es,crv,nnpcrv,failel,cma,0)
c
        do k=1,6
          do j=1,6
            dsave(i,j,k)=es(j,k)
          enddo
        enddo
      enddo
c
      return
      end
      subroutine utan41v(cm,d1,d2,d3,d4,d5,d6,sig1,sig2,
     . sig3,sig4,sig5,sig6,epsps,hsvs,lft,llt,dt1siz,capa,
     . etype,tt,temps,dsave,nlqa,crv,nnpcrv,failels,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      include 'nhisparm.inc'
      dimension d1(*),d2(*),d3(*),d4(*),d5(*),d6(*)
      dimension sig1(*),sig2(*),sig3(*),sig4(*),sig5(*),sig6(*)
      dimension epsps(*),hsvs(nlq,*),dt1siz(*),cm(*),qmat(nlq,3,3)
      dimension temps(*),dsave(nlq,6,*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      logical failels(*)
      character*5 etype
c
      dimension sig(6),eps(6),hsv(NHISVAR),es(6,6)
      logical failel
c
      do i=lft,llt
        sig(1)=sig1(i)
        sig(2)=sig2(i)
        sig(3)=sig3(i)
        sig(4)=sig4(i)
        sig(5)=sig5(i)
        sig(6)=sig6(i)
c
        eps(1)=d1(i)
        eps(2)=d2(i)
        eps(3)=d3(i)
        eps(4)=d4(i)
        eps(5)=d5(i)
        eps(6)=d6(i)
c
        do k=1,6
          do j=1,6
            es(j,k)=0.
          enddo
        enddo
c
        failel=failels(i)
c
        call utan41(cm,eps,sig,epsps(i),hsv,dt1siz(i),capa,etype,tt,
     1   temps(i),es,crv,nnpcrv,failel,cma,0)
c
        do k=1,6
          do j=1,6
            dsave(i,j,k)=es(j,k)
          enddo
        enddo
      enddo
c
      return
      end
      subroutine utan41(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,nnpcrv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
c
      factor=1.
      if (failel) factor=1.e-8
c
      g=factor*.5*abs(cm(1))/(1.+cm(2))
      b=factor*abs(cm(1))/3./(1.-2.*cm(2))
      bg23=b-2.*g/3.
      bg43=b+4.*g/3.
c
      es(1,1)=bg43
      es(2,2)=bg43
      es(3,3)=bg43
      es(2,1)=bg23
      es(3,1)=bg23
      es(3,2)=bg23
      es(1,2)=es(2,1)
      es(1,3)=es(3,1)
      es(2,3)=es(3,2)
      es(4,4)=g
      es(5,5)=g
      es(6,6)=g
c
      return
      end
      subroutine utan45(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,nnpcrv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*)
      dimension es(6,*),crv(lq1,2,*),cma(*),qmat(3,3)
      integer nnpcrv(*)
      logical failel
      character*5 etype
c
      logical usemodule,usemetal
c
c     See comments in subroutine umat45 for use of *MODULE
c
      usemodule=.false.
      if (usemodule) then
        usemetal=.false.
        if (usemetal) then
          call metaltan(cm,eps,sig,epsp,hsv,dt1,capa,etype,time,temp,
     .          failel,crv,nnpcrv,es)
          return
        endif
      endif
c
c     no history variables, NHV=0
c     deformation gradient stored in hsv(1),...,hsv(9)
c
c     compute jacobian
c
      detf=hsv(1)*(hsv(5)*hsv(9)-hsv(6)*hsv(8))
     1     -hsv(2)*(hsv(4)*hsv(9)-hsv(6)*hsv(7))
     2     +hsv(3)*(hsv(4)*hsv(8)-hsv(5)*hsv(7))
c
c     compute lame parameters
c
      xlambda=cm(1)*cm(2)/((1.+cm(2))*(1.-2.*cm(2)))
      xmu=.5*cm(1)/(1.+cm(2))
c
c     compute tangent stiffness
c     same for both shells and solids
c
      detf=max(detf,1.e-8)
      detfinv=1./detf
      dmu=xmu-xlambda*log(detf)
      es(1,1)=detfinv*(xlambda+2.*dmu)
      es(2,2)=detfinv*(xlambda+2.*dmu)
      es(3,3)=detfinv*(xlambda+2.*dmu)
      es(4,4)=detfinv*dmu
      es(5,5)=detfinv*dmu
      es(6,6)=detfinv*dmu
      es(2,1)=detfinv*xlambda
      es(3,2)=detfinv*xlambda
      es(3,1)=detfinv*xlambda
      es(1,2)=es(2,1)
      es(2,3)=es(3,2)
      es(1,3)=es(3,1)
c
      return
      end
      subroutine utan43(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,nnpcrv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
c
      return
      end
      subroutine utan44(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,nnpcrv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
c
      return
      end
      subroutine utan42(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,nnpcrv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
c
      return
      end
      subroutine utan46(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,nnpcrv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
c
      return
      end
      subroutine utan47(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,nnpcrv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
c
      return
      end
      subroutine utan48(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,nnpcrv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
c
      return
      end
      subroutine utan49(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,nnpcrv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
c
      return
      end
      subroutine utan50(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,nnpcrv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      integer nnpcrv(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
      return
      end
      subroutine uweldfail(idweld,strr,fail,fibl,cm,tt,lft,llt)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     spotweld failure routine
c
c***  local coordinates: x is tangent to beam, y & z are normal
c
c***  variables
c          idweld ---- weld id number
c          strr ------ stress resultants
c                      (1) axial (x direction) force
c                      (2) y shear force
c                      (3) z shear force
c                      (4) moment about z
c                      (5) moment about y
c                      (6) torsional resultant
c          fail ------ failure flag
c                      =0 not failed
c                      =1 fail element
c          fibl ------ fibl(1,i)=diameter of weld
c          cm -------- 6 constants supplied by user
c          tt -------- current simulation time
c          lft,llt --- do-loop range for strr
c
      INTEGER8 idweld
      dimension strr(8,*),fail(*),cm(*),fibl(5,*)
c
c***  demonstration failure model
c
c***  geometry
c
c***  failure constants
      if (cm(1).gt.0.) then
        oamax=1./cm(1)
      else
        oamax=0.
      endif
      if (cm(2).gt.0.) then
        osmax=1./cm(2)
      else
        osmax=0.
      endif
c
c***  evaluate failure criterion:
c       axial stress = axial force/area + moment/sectionmod
c       shear stress = shear force/area + torsion/(2.*sectionmod)
c       failure surface = (axial stress/cm(1))**2 +
c                         (shear stress/cm(2))**2
      do i=lft,llt
      Z       =3.14159265359*fibl(1,i)**3/32
      osctmod =1.0/Z
      hosctmod=.5*osctmod
      oarea   =4./(3.14159265359*fibl(1,i)**2)
        axial=max(0.,strr(1,i))*oarea+
     &          sqrt(strr(4,i)**2+strr(5,i)**2)*osctmod
        shear=sqrt(strr(2,i)**2+strr(3,i)**2)*oarea+
     &        abs(strr(6,i))*hosctmod
        if ((axial*oamax)**2+(shear*osmax)**2.gt.1.0) fail(i)=1.
      enddo
c
      return
      end
      subroutine uweldfail22(idweld,strr,rscale,fail,fibl,cm,tt,lft,llt)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     sample user defined spotweld failure routine with 46 parameters
c
c
c***  local coordinates: x is tangent to beam, y & z are normal
c
c***  variables
c          idweld ---- weld id number
c          strr stress resultants at beam ends 1 and 2
c                      (1) axial (x direction) force
c                      (2) y shear force
c                      (3) z shear force
c                      (4) moment about z
c                      (5) moment about y
c                      (6) torsional resultant
c          fail ------ failure flag
c                      =0 not failed
c                      =1 fail element
c          fibl ------ fibl(1,i)=diameter of weld
c                      fibl(2,i)=pointer to rscale data
c          rscale ---- resultant force scale factor for
c                      *CONTROL_SPOTWELD_BEAM effects
c                      (1) node 1 tensile factor
c                      (2) node 2 tensile factor
c                      (3) node 1 shear factor
c                      (4) node 2 shear factor
c          cm -------- up to 46 constants supplied by user
c          tt -------- current simulation time
c          lft,llt --- do-loop range for strr
c
c     *** NOTE *** for solid element welds, the rscale factors are not
c                  useful and are set to 1.0
c
      INTEGER8 idweld
      dimension strr(8,*),rscale(4,*),fail(*),cm(*),fibl(5,*)
c
c***  demonstration failure model
c
c***  failure constants
      if (cm(1).gt.0.) then
        oamax=1./cm(1)
      else
        oamax=0.
      endif
      if (cm(2).gt.0.) then
        osmax=1./cm(2)
      else
        osmax=0.
      endif
c
c***  evaluate failure criterion:
c       axial stress = axial force/area + moment/sectionmod
c       shear stress = shear force/area + torsion/(2.*sectionmod)
c       failure surface = (axial stress/cm(1))**2 +
c                         (shear stress/cm(2))**2
      do i=lft,llt
      Z       =3.14159265359*fibl(1,i)**3/32
      osctmod =1.0/Z
      hosctmod=.5*osctmod
      oarea   =4./(3.14159265359*fibl(1,i)**2)
      axial=max(0.,strr(1,i))*oarea+
     &        sqrt(strr(4,i)**2+strr(5,i)**2)*osctmod
      shear=sqrt(strr(2,i)**2+strr(3,i)**2)*oarea+
     &      abs(strr(6,i))*hosctmod
      j=nint(fibl(2,i))
      axial1=axial*rscale(1,j)
      axial2=axial*rscale(2,j)
      shear1=shear*rscale(3,j)
      shear2=shear*rscale(4,j)
      if ((axial1*oamax)**2+(shear1*osmax)**2.gt.1.0) fail(i)=1.
      if ((axial2*oamax)**2+(shear2*osmax)**2.gt.1.0) fail(i)=1.
      enddo
c
      return
      end
      subroutine uweldfail12(eps_to_fail,beta,dmg,uhisv,num_uhisv,
     . idweld,iwtype,strr,eps,esr,npt,diameter,cm,idmgopt,tt,ncycle,
     . lft,llt)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     Sample user defined spotweld damage and failure routine with 46
c     parameters; this routine can be used to damage or fail beam welds
c     and hex weld assemblies.
c
c     There are two ways to use it.  The first way is to do a damage
c     initiation calculation and set eps_to_fail and beta when damage
c     should initiate.  Based on the values set (see below) either
c     immediate failure, or linear or exponential damage and failure
c     based on effective plastic strain will occur.  The second way is
c     for the user to do their own dmaage calculation and update the
c     dmg variable thoughout the calculation.
c           
c    
c     subroutine output variables (may be updated by user)
c
c        eps_to_fail  Plastic strain from damage initiation until
c                     failure.  A zero value means that damage has not
c                     yet begun.  Any negative value will cause
c                     immediate failure.  In the cycle when eps_to_fail
c                     is set, the current effective plastic strain will
c                     be stored and used as the strain at which damage
c                     inititiates (eps_di).  Note that eps_to_fail and
c                     beta can only be set only one time for each
c                     integration point and subsequent changes will be
c                     ignored.
c
c        beta ....... Set to positive number for exponential
c                     damage.  Set to zero for linear damage.
c
c                     For beta=0, the decay function is
c
c                        dmg=(eps-eps_di)/eps_to_fail
c
c                      For beta>0, the decay function is
c
c                        dmg=1.-exp(-beta*(eps-eps_di)/eps_to_fail)
c
c                     For exponential damage, the damage will grow
c                     according to the exponential equation above until
c                     the plastic strain reaches eps_to_fail at which
c                     time the weld will be failed.
c
c        uhisv ...... Weld history values.  These values will be
c                     retained, but will not be used for the
c                     calculation.
c
c     subroutine input/output variable (may be either calculated by
c     ls-dyna or updated by the user)
c     
c
c        dmg ........ Current damage at the integration point.
c                     Initially, dmg=0.  When dmg=1.0, the integration
c                     point will fail.  When all integration points
c                     fail, the element or assembly will be deleted.
c                     If the user wants ls-dyna to calculate damage
c                     using eps_to_fail and beta, then damage should
c                     not be changed.  Optionally, the user may update
c                     the dmg variable throughout the calculation and
c                     the user's value will be applied.
c                    
c
c     subroutine input variables
c
c        num_uhisv .. number of weld history variables (will be
c                     equal to FVAL on card 3 of *MAT_SPOTWELD)
c        idweld ..... weld id number
c
c        iwtype ..... weld type
c                      1=beam element
c                      2=solid element assembly
c
c        strr ....... stress resultants
c                       (1) axial (x direction) force
c                       (2) y shear force
c                       (3) z shear force
c                       (4) moment about z
c                       (5) moment about y
c                       (6) torsional mement
c         
c        eps ........ effective plastic strain at the integration point
c                     (or hex weld of an assembly)
c
c        esr ........ effective strain rate at the integration point
c                     (or hex weld of an assembly)
c
c        npt ........ number of integration points in the weld; for hex
c                     weld assemblies, npt is the number of elements
c                     in the assembly since each element has 1 integration
c                     point 
c
c        diameter ... diameter of weld
c        cm ......... up to 46 constants supplied by user
c        idmgopt .... damage option parameter DMGOPT
c        tt ......... current simulation time
c        ncycle ..... current cycle number
c        lft,llt .... do-loop range
c
c
c     For hex assmeblies, some variables may be dumped to d3plot using
c     NEIPH on *DATABASE_EXTENT_BINARY.  The table below shows the
c     variable number for each.
c
c        variable name    extra history variable #
c
c             dmg                   8
c            eps_di                10
c          eps_to_fail             11
c             beta                 12
c             esr                  13
c
      INTEGER8 idweld
      dimension eps_to_fail(*),beta(*),dmg(npt,*),uhisv(num_uhisv,*),
     . strr(8,*),eps(npt,*),esr(npt,*),diameter(*),cm(*)
c
c***  demonstration failure model
c
c***  failure constants
      if (cm(1).gt.0.) then
        oamax=1./cm(1)
      else
        oamax=0.
      endif
      if (cm(2).gt.0.) then
        osmax=1./cm(2)
      else
        osmax=0.
      endif
c
c***  evaluate damage initiation criterion:
c       axial stress = axial force/area + moment/sectionmod
c       shear stress = shear force/area + torsion/(2.*sectionmod)
c       failure surface = (axial stress/cm(1))**2 +
c                         (shear stress/cm(2))**2
      do i=lft,llt
      if (eps_to_fail(i).eq.0.) then
        Z       =3.14159265359*diameter(i)**3/32
        osctmod =1.0/Z
        hosctmod=.5*osctmod
        oarea   =4./(3.14159265359*diameter(i)**2)
c
        axial=max(0.,strr(1,i))*oarea+
     &         sqrt(strr(4,i)**2+strr(5,i)**2)*osctmod
        shear=sqrt(strr(2,i)**2+strr(3,i)**2)*oarea+
     &         abs(strr(6,i))*hosctmod
c
        if ((axial*oamax)**2+(shear*osmax)**2.gt.1.0) then
          eps_to_fail(i)=cm(3)
          beta(i)=cm(4)
        endif
      endif
      enddo
      return
      end
      subroutine hdot(lft,llt,hmin,hmax,iactive,e,h,dt,const,s)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     user integration in time of the smoothing length
c
c     lft      : left indice of the loop
c     llt      : right indice of the loop
c     hmin     : Minimum smoothing length defined in the corresponding
c                *SECTION_SPH
c     hmax     : Maximum smoothing length defined in the corresponding
c                *SECTION_SPH
c     iactive  : variable that indicates if the particle is active or not
c                = 1, particle is active
c                = 0, particle is not active
c     e        : tensor of gradient of deformation
c     h        : smoothing length
c     dt       : current time step
c     const    : constant = dt/dimension of the problem
c     s        : stress tensor
c
      dimension h(1),iactive(1),e(9,1),ro(1),smthi(1),matsph(1)
      dimension s(6,1)
c
C_DOACROSS
C_SHARE1(lft,llt,iactive,e,h,const,dt1,hmin,hmax)
C_LOCAL(i,divv,hnew)
c
      do 200 i=lft,llt
      divv = e(1,i) + e(2,i) + e(3,i)
c
      hnew = (1.+const*divv)*h(i)
      h(i) = max(hmin,hnew)
      h(i) = min(h(i),hmax)
 200  continue
c
      return
      end
#ifdef MPP
      subroutine ulsdyanfunc()
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     loadable functions from dyna
c
      include 'xjobid.inc'
      INTEGER8 lqfinv8,lqfinvf8
      character*80 mess
      INTEGER8 msr8,lcid8,nn8
      real*8 fchngl
c
      ij=lqf(ix,ityp)
      ij=lqfb(ix,msr,n,ierr)
c     ij=lqfbel(ix,ityp,mess)
      ij=lqfcfd(ix,nodes,ncfdsg,ierr)
      ij=lqfc(i,narray)
      ij=lqfcinv(i)
      ij=lqfcsinv(i)
      ij=lqfdei(i,iopt)
      ij=lqfdem(i,iopt)
      ij=lqfm(ix)
      ij=lqfmiv(ix)
      ij=lqfsti(ix,ityp)
      ij=lqfl(ix)
c     ij=lqfibe(ix,ityp)
      ij=lqfint(iextn,ityp,ierr)
      ij=lqfinb(ix,msr,n)
      ij=lqfinv(ix,ityp)
      ij=lqfinvf(ix,ityp)
      ij=lqfinv4(ix,ityp)
      ij=lqfinvf4(ix,ityp)
      ij=lqfinv8(ix,ityp)
      ij=lqfinvf8(ix,ityp)
      ij=lqjnt(jntid)
      ij=lqfmrs(ix)
      ij=lqfsi(ix,msr,n)
      ij=lqftm(ix,msr,n)
      ij=lqf8(ix,ityp)
      ij=lqlcs(lcsid)
      ij=lqsec(idsec)
c
c     ij=lqfnam(i1,i2,cname,i3)
c     ij=lqfpid(ipextid,4,jerr)
c
      ij=funct(s,t)
      ij=s1dwrn(islno,ncard)
      ij=rlcid(flcid)
      ij=ilcid(lcid)
      ij=rcsid(flcid)
      ij=lcids(lcsid)
      ij=icsid(lcid)
      ij=ivsid(lcid,nn8)
      ij=arclqh(r,z,crv,n)
      ij=fchngl(i,j)
      ij=rlqfl(rix)
      ij=pchst(arg1,arg2)
c
c  dummy calls
      call jidname()
c     call crvval()
c     call crvval_v()
c     call tabval()
c     call tabval_v()
      call compute_f3s()
      call compute_f3()
      call push_forward_2s()
      call push_forward_2()
      call push_forward_4s()
      call push_forward_4() 
c
c MPI routines we might need to make available
c
#if defined AUTODOUBLE && !defined MPI_AUTODOUBLE
#define MPI_ALLGATHER           my_allgather
#define MPI_ALLGATHERV          my_allgatherv
#define MPI_ALLREDUCE           my_allreduce
#define MPI_ALLTOALL            my_alltoall
#define MPI_ALLTOALLV           my_alltoallv
#define MPI_BARRIER             my_barrier
#define MPI_BCAST               my_bcast
#define MPI_CANCEL              my_cancel
#define MPI_COMM_CREATE         my_comm_create
#define MPI_COMM_DUP            my_comm_dup
#define MPI_COMM_FREE           my_comm_free
#define MPI_COMM_GROUP          my_comm_group
#define MPI_COMM_RANK           my_comm_rank
#define MPI_COMM_SIZE           my_comm_size
#define MPI_FINALIZE            my_finalize
#define MPI_GATHER              my_gather
#define MPI_GATHERV             my_gatherv
#define MPI_GET_COUNT           my_get_count
#define MPI_GROUP_FREE          my_group_free
#define MPI_GROUP_INCL          my_group_incl
#define MPI_INIT                my_init
#define MPI_IRECV               my_irecv
#define MPI_ISEND               my_isend
#define MPI_RECV                my_recv
#define MPI_REDUCE              my_reduce
#define MPI_REDUCE_SCATTER      my_reduce_scatter
#define MPI_REQUEST_FREE        my_request_free
#define MPI_SEND                my_send
#define MPI_TEST                my_test
#define MPI_WAIT                my_wait
#define MPI_WAITALL             my_waitall
#define MPI_ABORT               my_abort
#define MPI_WTIME               my_wtime
#define MPI_PROBE               my_probe
#endif
      call MPI_ALLGATHER
      call MPI_ALLGATHERV
      call MPI_ALLREDUCE
      call MPI_ALLTOALL
      call MPI_ALLTOALLV
      call MPI_BARRIER
      call MPI_BCAST
      call MPI_CANCEL
      call MPI_COMM_CREATE
      call MPI_COMM_DUP
      call MPI_COMM_FREE
      call MPI_COMM_GROUP
      call MPI_COMM_RANK
      call MPI_COMM_SIZE
      call MPI_FINALIZE
      call MPI_GATHER
      call MPI_GATHERV
      call MPI_GET_COUNT
      call MPI_GROUP_FREE
      call MPI_GROUP_INCL
      call MPI_INIT
      call MPI_IRECV
      call MPI_ISEND
      call MPI_RECV
      call MPI_REDUCE
      call MPI_REDUCE_SCATTER
      call MPI_REQUEST_FREE
      call MPI_SEND
      call MPI_TEST
      call MPI_WAIT
      call MPI_WAITALL
      call MPI_ABORT
      call MPI_WTIME
      call MPI_PROBE
      return
      end
#endif
      subroutine usrhcon(h,uc,nc,d,p,xs,xm,vs,vm,ts,tm,
     . istria,imtria,crv,nnpcrv)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     User subroutine for interface conductivity control (mpp/smp)
c
c     To define thermal contact conductivity
c
c     Variables:
c
c     h      = thermal contact conductivity (output)
c     uc(nc) = array of user defined parameters (input)
c     nc     = number of user defined conductivity parameters (input)
c     d      = distance (gap) between slave and master segment (input)
c     p      = interface pressure (input)
c     xs     = slave segment coordinates (input)
c     xm     = master segment coordinates (input)
c     vs     = slave segment velocities (input)
c     vm     = master segment velocities (input)
c     ts     = slave segment temperatures (input)
c     tm     = master segment temperatures (input)
c     istria = 1 iff slave segment is a triangle (input)
c     imtria = 1 iff master segment is a triangle (input)
c     crv    = curve array (input)
c     nnpcrv = # of discretization points per crv()
c
c     Note that usrvrel can be called to compute relative sliding velocity
c
      include 'nlqparm'
      dimension uc(nc)
      dimension xs(3,4),xm(3,4)
      dimension vs(3,4),vm(3,4)
      dimension crv(lq1,2,*)
      integer nnpcrv(*)
c
c     Define conductivity as function of pressure via load curve
c
      call crvval(crv,nnpcrv,uc(1),p,h,dh)
c
      return
      end
      subroutine usrvrel(vrel,xs,vs,vm,istria,imtria)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     compute sliding velocity vrel and that can be called from usrhcon.
c
c     Variables:
c     vrel    = relative sliding velocity (output)
c     xs      = slave segment coordinates (input)
c     vs      = slave segment velocities (input)
c     vm      = master segment velocities (input)
c     istria  = 1 iff slave segment is a triangle (input)
c     imtria  = 1 iff master segment is a triangle (input)
c
      dimension xs(3,4),vs(3,4),vm(3,4)
      dimension xnrm(3),vr(3)
c
      jmax=4
      if (istria.eq.1) jmax=3
      do i=1,3
         vr(i)=0.
         do j=1,jmax
            vr(i)=vr(i)+vs(i,j)
         enddo
         vr(i)=vr(i)/jmax
      enddo
      jmax=4
      if (imtria.eq.1) jmax=3
      do i=1,3
         vrl=0.
         do j=1,jmax
            vrl=vrl+vm(i,j)
         enddo
         vr(i)=vr(i)-vrl/jmax
      enddo
      x31=xs(1,3)-xs(1,1)
      y31=xs(2,3)-xs(2,1)
      z31=xs(3,3)-xs(3,1)
      x42=xs(1,4)-xs(1,2)
      y42=xs(2,4)-xs(2,2)
      z42=xs(3,4)-xs(3,2)
      xnrm(1)=y31*z42-z31*y42
      xnrm(2)=z31*x42-x31*z42
      xnrm(3)=x31*y42-y31*x42
      xlen=sqrt(xnrm(1)**2+xnrm(2)**2+xnrm(3)**2)
      xnrm(1)=xnrm(1)/xlen
      xnrm(2)=xnrm(2)/xlen
      xnrm(3)=xnrm(3)/xlen
      vnrm=vr(1)*xnrm(1)+vr(2)*xnrm(2)+vr(3)*xnrm(3)
      vr(1)=vr(1)-vnrm*xnrm(1)
      vr(2)=vr(2)-vnrm*xnrm(2)
      vr(3)=vr(3)-vnrm*xnrm(3)
      vrel=sqrt(vr(1)**2+vr(2)**2+vr(3)**2)
      return
      end
      subroutine usrflux(fl,flp,x,tnpl,tnl,nodes,
     . alpha,atime,atemp,dt,time,fhsv,nfhsv,crv,nnpcrv,v)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     User subroutine for boundary thermal flux
c
c     To define thermal flux parameter (heat per surface area and time)
c
c     Variables:
c
c     fl              = flux intensity (output)
c     flp             = flux intensity derivative wrt atemp (output)
c     x(3,nodes)      = global segment coordinates (input)
c     v(3,nodes)      = global segment velocities (input)
c     tnpl(nodes)     = temperatures at time time (input)
c     tnl(nodes)      = temperatures at time time-dt (input)
c     nodes           = number of nodes in segment (3,4 or 6) (input)
c     alpha           = time integration parameter (input)
c     atime           = time+(alpha-1)*dt
c     atemp           = average segment temperature at time atime
c     dt              = time step size (input)
c     time            = time at which the new temperature is sought (input)
c     fhsv(nfhsv)     = flux history variables (input/output)
c     nfhsv           = number of flux history variables for this segment (input)
c     crv             = curve array (input)
c     nnpcrv          = # of discretization points per crv()
c
      include 'nlqparm'
      dimension x(3,*),tnpl(*),tnl(*),v(3,*)
      dimension fhsv(*),crv(lq1,2,*),func(9)
      integer nnpcrv(*)
c---------------------------------------------------------------------
c Example 1
c     Define flux by linear convection
c     that optionally decays (in an ad-hoc way) as power
c     dissipates from surface
c
c     fhsv(1) = convection coefficient
c     fhsv(2) = ambient temperature
c     fhsv(3) = total amount of energy per surface area available
c     fhsv(4) = dissipated energy per surface area at current
c
c---- uncomment following lines for example 1
c     hcon=fhsv(1)
c     tinf=fhsv(2)
c     flin=hcon*(tinf-atemp)
c     if (nfhsv.gt.2) then
c        q=(1.-fhsv(4)/fhsv(3))/
c    .        (1.+.5*dt*flin/fhsv(3))
c        flp=-q*hcon
c        if (q.gt.1.) then
c           q=1.
c           flp=-hcon
c        elseif (q.lt.0.) then
c           q=0.
c           flp=0.
c        endif
c        fl=q*flin
c        fhsv(4)=fhsv(4)+dt*.5*fl
c        fhsv(4)=min(fhsv(3),fhsv(4))
c     else
c        fl=flin
c        flp=-hcon
c     endif
c     fl=-fl
c     flp=-flp
c---------------------------------------------------------------------
c Example 2
c     Define flux using a function ID
c
c     fhsv(1) =  function ID
c
      fl=0.
      flp=0.
c-- get segment centroid coordinates
      func(1)=(x(1,1)+x(1,2)+x(1,3)+x(1,4))/4.
      func(2)=(x(2,1)+x(2,2)+x(2,3)+x(2,4))/4.
      func(3)=(x(3,1)+x(3,2)+x(3,3)+x(3,4))/4.
c-- get segment average velocity
      func(4)=(v(1,1)+v(1,2)+v(1,3)+v(1,4))/4.
      func(5)=(v(2,1)+v(2,2)+v(2,3)+v(2,4))/4.
      func(6)=(v(3,1)+v(3,2)+v(3,3)+v(3,4))/4.
c-- get segment average temperature
      func(7)=(tnl(1)+tnl(2)+tnl(3)+tnl(4))/4.
c-- T_infinity is not a typical input for flux, set it to zero
      func(8)=0.
c-- get time
      func(9)=atime
c-- the first history variable for this segment is a function id
c-- make it an integer
      idff=int(fhsv(1))
c-- convert user function ID to the internal numbering system
      internal_f=ind_gf(idff)
c-- check to see if internal_f is valid. If so, then evaluate function
      if (internal_f .gt. 0) then
        fl=eval_function(internal_f,func,ierr)
        flp=0
      endif
c
      return
      end
      subroutine ujntfrc(idfrc,idjnt,const,t,dt,statevar,disp,vel,frc)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c***  user subroutine for applying forces to joint degrees of freedom
c
c***  variables
c          idfrc ----- ID for the joint force
c          idjnt ----- ID for the joint
c          const ----- array of user supplied constants
c          t --------- current time
c          dt -------- time step size
c          statevar -- array of user integrated state variables
c          disp ------ array of relative joint displacements
c          vel ------- array of relative joint velocities
c          frc ------- array of forces and torques calculated by user
c
      dimension const(*),statevar(*),disp(*),vel(*),frc(*)
c
c***  example for a spring and damper for a revolute joint
c     torque =  c*angular_velocity + k*(angle-offset)
      frc(1)=const(1)*vel(1)+const(2)*(disp(1)-const(3))
c
c***  statevariable integration: statevar(1)=0.5*t**2
      statevar(1)=statevar(1)+dt*t
c
      return
      end
      subroutine ujntfdrv(njntf,idat,kdat,idatjf,kdatjf,const,
     &                    statevar,q,rbv,rbf,x)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'iounits.inc'
c
c***  driver for the user joint force subroutines
c
c***  variables
c          njntf ----- number of joint forces
c          idatjf ---- integer user joint force data
c                      1) ID
c                      2) joint ID
c                      3) number of state variables
c                      4) joint index
c                      5) rigid body 1
c                      6) rigid body 2
c                      7) pointer into state variables
c          idat ------ 1) joint type
c                      2:7) nodes
c          const ----- user constants
c          statevar -- user state variables
c          q --------- integrated displacements
c          rbv ------- rigid body velocities
c          rbf ------- rigid body forces
c
c***  overall structure
c
c       qdot = B*{rigid body velocities}
c       q    = q + dt*qdot
c       f    = user_function(q,qdot,t,dt)
c       f_rb = f_rb - B^T*f
c       ener = ener + dt*f*qdot
c
      common/bk02/iburn,dt1,dt2,isdo,tet10jtol
      common/bk28/summss,xke,xpe,tt,xte0
c
      real*8 x
      character*19 atype(14)
      dimension idat(kdat,*),idatjf(kdatjf,*),const(48,*),
     &          statevar(*),q(6,*),rbv(6,*),rbf(6,*),x(3,*)
c
      data atype/'spherical          ',
     2           'revolute           ',
     3           'cylindrical        ',
     4           'planar             ',
     5           'universal          ',
     6           'translational      ',
     7           'locking            ',
     8           'translational motor',
     9           'rotational motor   ',
     .           'gears              ',
     1           'rack and pinion    ',
     2           'constant velocity  ',
     3           'pulley             ',
     4           'screw              '/
c
c***  over joint forces
      do 1000 ijntf=1,njntf
c
c***    the basics
        idjntf=idatjf(1,ijntf)
        idjnt =idatjf(2,ijntf)
        indjnt=idatjf(4,ijntf)
        irb1  =idatjf(5,ijntf)
        irb2  =idatjf(6,ijntf)
        ipsv  =idatjf(7,ijntf)
        jnttyp=idat(indjnt,1)
c
c***    branch on joint type
        go to (  10,  20,  30,  40,  50,  60,  70,
     &           80,  90, 100, 110, 120, 130, 140), jnttyp
c
c***    spherical
   10   go to 911
c
c***    revolute
   20   call rjufrc(idjntf,idjnt,const(1,ijntf),tt,dt1,
     &              statevar(ipsv),q(1,ijntf),
     &              rbv(1,irb1),rbf(1,irb1),rbv(1,irb2),rbf(1,irb2),
     &              x,idat,kdat,indjnt)
        go to 1000
c
c***    cylindrical
   30   go to 911
c
c***    planar
   40   go to 911
c
c***    universal
   50   go to 911
c
c***    translational
   60   go to 911
c
c***    locking
   70   go to 911
c
c***    tranlational motor
   80   go to 911
c
c***    rotational motor
   90   go to 911
c
c***    gears
  100   go to 911
c
c***    rack and pinion
  110   go to 911
c
c***    constant velocity
  120   go to 911
c
c***    pulley
  130   go to 911
c
c***    screw
  140   go to 911
c
c***    error message -- should never get here
c 911   write(iohsp,1010) idjntf,atype(jnttyp),idjnt
c       write(iomsg,1010) idjntf,atype(jnttyp),idjnt
c       write(iotty,1010) idjntf,atype(jnttyp),idjnt
  911   continue
        ierdat(1)=idjntf
        cerdat(1)=atype(jnttyp)
        ierdat(2)=idjnt
        call lsmsg(3,MSG_SOL+1153,ioall,ierdat,rerdat,cerdat,0)
c
 1000 continue
c
      return
c1010 format(1x/
c    &  '*** Error user joint force ID',i12,' for ',a19,i12/
c    &  '          this type of joint is currently not supported',
c    &  '          for user joint forces.')
      end
      subroutine rjufrc(idfrc,idjnt,const,t,dt,statevar,q,
     &                  rbv1,rbf1,rbv2,rbf2,x,idat,kdat,indjnt)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c***  revolute joint user force driver
c
c***  variables
c          idfrc ----- id of user force
c          idjnt ----- id of joint
c          const ----- user constants
c          t --------- current time
c          dt -------- time step size
c          statevar -- user state variables
c          q --------- displacements
c          rbv1,2 ---- rigid body velocities
c          rbf1,2 ---- rigid body forces
c          x --------- coordinates
c          idat ------ list of joint nodes
c          indjnt ---- joint indext
c
      real*8 x
      dimension const(*),statevar(*),q(*),
     &          rbv1(*),rbf1(*),rbv2(*),rbf2(*),x(3,*),idat(kdat,*)
c
c***  nodes
      n1=idat(indjnt,2)
      n2=idat(indjnt,3)
      n3=idat(indjnt,4)
      n4=idat(indjnt,5)
c
c***  angular velocity
      dx11=x(1,n3)-x(1,n1)
      dx21=x(2,n3)-x(2,n1)
      dx31=x(3,n3)-x(3,n1)
      omag=1.0/sqrt(dx11**2+dx21**2+dx31**2)
      dx11=omag*dx11
      dx21=omag*dx21
      dx31=omag*dx31
c
      dx12=x(1,n4)-x(1,n2)
      dx22=x(2,n4)-x(2,n2)
      dx32=x(3,n4)-x(3,n2)
      omag=1.0/sqrt(dx12**2+dx22**2+dx32**2)
      dx12=omag*dx12
      dx22=omag*dx22
      dx32=omag*dx32
c
      angv=rbv2(4)*dx12+rbv2(5)*dx22+rbv2(6)*dx32
     &    -rbv1(4)*dx11-rbv1(5)*dx21-rbv1(6)*dx31
c
c***  angular displacement integrated in time
c     copy so user can't overwrite our history variable value
      q(1)=q(1)+dt*angv
      ang=q(1)
c
c***  call user routine
      call ujntfrc(idfrc,idjnt,const,t,dt,statevar,ang,angv,torq)
c
c***  add torq to rhs
      rbf2(4)=rbf2(4)-dx12*torq
      rbf2(5)=rbf2(5)-dx22*torq
      rbf2(6)=rbf2(6)-dx32*torq
      rbf1(4)=rbf1(4)+dx11*torq
      rbf1(5)=rbf1(5)+dx21*torq
      rbf1(6)=rbf1(6)+dx31*torq
c
      return
      end
      subroutine matfailusercontrol(user)
      dimension user(*)
c
c     user interface to for user failure subroutine
c     when iuserfail equals to 1, user has the option to access stress, 
c     and history information for each element
c     Then user can make any change with 'subroutine matuserfail'
c     
c
      iuserfail=0
      if(iuserfail.eq.1) then
      call matfailusercontrol1(user)
      end if
      return
      end
      subroutine matuserfail(stress,shisv,nip,nplane,nhis)
      dimension stres(8,4,9),shisv(90,4,9)
c
c     user can access stress and history information now 
c     so, user's failure model can be defined here.
c
c     no more than 9 integration point in the thickness
c
c     stress(1,*,*): zeta (local coordinate) in the thickness direction
c     stress(8,*,*): effective strain
c
c     shisv(n,*,*) the nth history variable in one integration
c
c     Stress and history can be obtained in the following format
c     do i = 1, nplane
c     do j = 1, nip
c       write(*,'(t1,8e10.3)') (stres(k,i,j),k=1,8)
c       do m = 1, nhis,8
c         n = min(nhis,m+7)
c         write(*,'(t1,8e10.3)') (shisv(k,i,j),k=m,n)
c       end do
c     end do
c     end do

      return
      end
      subroutine rwumat(iphase,txt,lmc,ipnt,lmca,ipnta,idfoff,
     . fctmas,fcttim,fctlen,fcttem,incout)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'memaia.inc'
      include 'iodeg.inc'
c
c     In this routine, material parameter cards are read from
c     *MAT_USER_DEFINED_MATERIAL_MODELS in phase 1. 
c     In phase 2, the parameters are written as real numbers
c     to the structured input (dyna.str).
c
c     Main variables:
c       iphase : phase flag
c       txt    : incoming text
c       lmc    : number of material constants
c       ipnt   : pointer to material constants
c       lmca   : number of additional material constants
c       ipnta  : pointer to additional material constants
c
c     Offset and transformation factors from *INCLUDE_TRANSFORM 
c     can be used in phase 1:
c       idfoff : func id offset number, i.e., load curves and tables
c       fctmas : mass transformation factor
c       fcttim : time transformation factor
c       fctlen : length transformation factor
c       fcttem : temperature transformation factor
c       incout : the transformed included file will be written to dyna.inc
c                =0, no output file
c                =1, dyna.inc includes the complete transformed input
c
      character*CLENG txt
      character*80 mssg
      character*10 out10(8)
      dimension prop(48)
c
      isexample=0
c
c---  PHASE 1: read keyword input
c
c     In phase 1, the read format can be defined arbitrarily.
c     Parameters can be stored as real a(...) or integer ia(...)
c     as shown in the example. 
c
      if (iphase.eq.1) then
        mssg='readin *MAT_USER_DEFINED_MATERIAL_MODELS'
c
c       default read material constants: LMC
        if (isexample.eq.0) then
          do i=1,lmc,8
            iend=min(i+7,lmc)
            call gttxln (txt,lcount,2,8)
            read (txt,'(8e10.0)',err=900) (a(ipnt+k),k=i,iend)
          enddo
c
c       example: read 5 parameters
c       p1  - real with no unit 
c       p2  - real with time unit
c       p3  - real with stress unit
c       ip4 - integer with no unit
c       ip5 - load curve id (stress vs. strain)
        elseif (isexample.eq.1) then
          call gttxln (txt,lcount,2,8)
          read (txt,'(3e10.0,2i10)',err=900) p1,p2,p3,ip4,ip5
          p2=p2*fcttim
          fctfor=fctmas*fctlen/(fcttim*fcttim)
          fctpre=fctfor/(fctlen*fctlen)
          p3=p3*fctpre
          if (ip5.ne.0) then
            ip5=ip5+idfoff
c           adtbmdu : Add TaBle or Curve to be MoDified User-defined
c           sfa = scale factor for abscissa value
c           sf0 = scale factor for ordinate value
            sfa=1.0
            sfo=fctpre
            call adtbmdu(ip5,sfa,sfo,ipnt,ipnta,lmca)
          endif
          a(ipnt+1)=p1
          a(ipnt+2)=p2
          a(ipnt+3)=p3
          ia(ipnt+4)=ip4
          ia(ipnt+5)=ip5
        endif
c
c       read extra material constants: LMCA
        if (lmca.gt.0) then
          do i=1,lmca,8
            iend=min(i+7,lmca)
            call gttxln (txt,lcount,2,8)
            read (txt,'(8e10.0)',err=900) (a(ipnta+k),k=i,iend)
          enddo
        endif
c
c---  PHASE 2: write structured input
c
c     In phase 2, parameters must be transformed to real numbers,
c     before writing them to the structured input. Therefore,
c     external load curve ids should be changed to internal ids
c     as shown in the example. 
c
      elseif (iphase.eq.2) then
c
c       default write 48 parameters (LMC max)
        if (isexample.eq.0) then
          do i=1,48,8
            ii=min(48-i+1,8)
            do j=1,ii
              call form10(out10(j),a(ipnt+i+j-1))
            enddo
            write(io,'(8a10)') (out10(j),j=1,ii)
            call wttxsg(1)
          enddo
c
c       example: write 5 parameters and 43 zeros (sum=48!)
c       p1  - real with no unit 
c       p2  - real with time unit
c       p3  - real with stress unit
c       ip4 - integer with no unit
c       ip5 - load curve id (stress vs. strain)
        elseif (isexample.eq.1) then
          do i=1,48
            prop(i)=a(ipnt+i)
          enddo
          prop(4)=FLOAT(ia(ipnt+4))
c         change load curve id to internal number with lcidi
          lc=ia(ipnt+5)
          prop(5)=FLOAT(lcidi(lc))
          do i=1,48,8
            ii=min(48-i+1,8)
            do j=1,ii
c             change to special output format with form10
              call form10(out10(j),prop(i+j-1))
            enddo
            write(io,'(8a10)') (out10(j),j=1,ii)
            call wttxsg(1)
          enddo
        endif
c
c       write extra material constants: LMCA
        if (lmca.gt.0) then
          do i=1,lmca,8
            ii=min(lmca-i+1,8)
            write(io,'(1p8e10.3)') (a(ipnta+i+j-1),j=1,ii)
            call wttxsg(1)
          enddo
        endif
c
      endif
c
      return
  900 call terminn(txt,mssg,lcount,1)
      end
      subroutine usermatfpert ( x, y, z, fail )
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     user specification of material failure
c
c     Inputs:
c     x, y, z: coordinates
c
c     Returned:
c     fail:  failure criterion (stress, yield)
c
      fail = x
c
      return
      end
      subroutine chknhue(nhsv)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'iounits.inc'
      include 'nhisparm.inc'
c
      if (nhsv.gt.NHISVUE) then
        write(iotty,100) NHISVUE
        write(iohsp,100) NHISVUE
        write(iomsg,100) NHISVUE
      endif
      nhsv=min(NHISVUE,nhsv)
  100 format(/' *** Warning number of user element history variables',
     .       /'             is limited by NHISVUE in nhisparm.inc:',i10)
c
      return
      end
      subroutine chkusercomm
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     check the length of common block size
c
      include 'nlqparm'
      include 'nhisparm.inc'
      include 'iounits.inc'
c
      parameter (NBEL8A=max(100,2*(NHISVUE+45)))
      parameter (NBEL7 =max(2601,
     . 1180+576*NXDOFUE+64*NXDOFUE*NXDOFUE+12*MXSCND))
      common/aux14loc/ax14zz1(nlq,7+NHISVAR)
      common/bel8loc/bel8lz1(nlq,2*(NHISVUE+45))
      common/bel8aloc/bel8lz2(nlq,NBEL8A)
      common/bel7loc/bel7lz(nlq,NBEL7)
      common/usercom/iusercomm(4)
C_TASKCOMMON (aux14loc)
C_TASKCOMMON (bel8loc)
C_TASKCOMMON (bel8aloc)
C_TASKCOMMON (bel7loc)
c
      l14=nlq*(7+NHISVAR)
      l8 =nlq*2*(NHISVUE+45)
      l8a=nlq*NBEL8A
      l7 =nlq*NBEL7
c
      if (l14.gt.iusercomm(1)) then
        ierdat(1)=l14
        ierdat(2)=iusercomm(1)
        cerdat(1)='/aux14loc/'
        call lsmsg(2,MSG_OTH+81,iotty,ierdat,rerdat,cerdat,0)
        call cstop('E r r o r   t e r m i n a t i o n')
      endif
      if (l8.gt.iusercomm(2)) then
        ierdat(1)=l8
        ierdat(2)=iusercomm(2)
        cerdat(1)='/bel8loc/'
        call lsmsg(2,MSG_OTH+81,iotty,ierdat,rerdat,cerdat,0)
        call cstop('E r r o r   t e r m i n a t i o n')
      endif
      if (l8a.gt.iusercomm(3)) then
        ierdat(1)=l8a
        ierdat(2)=iusercomm(3)
        cerdat(1)='/bel8aloc/'
        call lsmsg(2,MSG_OTH+81,iotty,ierdat,rerdat,cerdat,0)
        call cstop('E r r o r   t e r m i n a t i o n')
      endif
      if (l7.gt.iusercomm(4)) then
        ierdat(1)=l7
        ierdat(2)=iusercomm(4)
        cerdat(1)='/bel7loc/'
        call lsmsg(2,MSG_OTH+81,iotty,ierdat,rerdat,cerdat,0)
        call cstop('E r r o r   t e r m i n a t i o n')
      endif
c
      return
      end
c
c user defined debonding law in *CONSTRAINED_BEAM_IN_SOLID
c*CONSTRAINED_BEAM_IN_SOLID
c#   slave    master     sstyp     mstyp     ctype     empty     nquad      idir
c        2         1         1         1         0                   2         1
c#   start       end              debond
c    0.000     0.000                 999
c
      subroutine rebar_bondslip_get_nhisvar(nhisvar)
c
c here user needs to request number of history variables per coupling point
c note: hisvar arrays allocated as float. If need double precision when
c running single precision executable, need make nhisvar twice the original value
      nhisvar=12
      return
      end
c
      subroutine rebar_bondslip_get_force(slip,dl,force,hsv,
     .  userparm,lawid)
      real hsv
      dimension hsv(12),cm(8),userparm(*)
c
c in this subroutine user defines debonding properties and
c call his own debonding subroutine to get force
c     du=slip
c     elsiz=dl 
c     cm(1)=1.0 
c     cm(2)=6.0
c     cm(3)=2.4*(cm(2)/5.0)**0.75
c     cm(8)=0.
c
c     call Bond_Slip(du,cm,hsv,elsiz,T,F)      
c     force=F
      cm(1)=userparm(1)
      cm(2)=userparm(2)
      cm(3)=2.4*(cm(2)/5.0)**0.75
      cm(8)=0.
c
      pi=3.1415926
      d=0.175
      area=pi*0.25*d*d*dl
      pf=1.0
      if (slip.lt.0.25) then
        shear=slip*pf
      else
        shear=0.25*pf
      endif
      force=sign(1.0,slip)*shear*area
      return
      end
c
c *SET_POROUS_LAGRANGIAN: user defined Aij and Bij matrix
c
      subroutine lagpor_getab_userdef(ab,iel,iuserdef,
     .                                ncon,ix,
     .                                xloc,yloc,zloc,
     .                                auxvec,locaux,
     .                                ro,cm,
     .                                x,v,a,
     .                                xms,tt,ncycle,tmpnd)
      real*8 x
      dimension ab(18),ix(ncon,*),xloc(3),yloc(3),zloc(3),
     .          auxvec(*),locaux(*),ro(*),cm(*),
     .          x(3,*),v(3,*),a(3,*),xms(*),tmpnd(*)
c
c ab(j) for j=1,9  : Aij matrix
c                    Axx=ab(1), Axy=ab(2), ... , Azz=ab(9)
c
c ab(j) for j=10,18: Bij matrix
c                    Bxx=ab(10), Bxy=ab(11), ... , Bzz=ab(18)
c
c iel              : internal Lagrangian element ID (structure element) 
c                    get user ID with lqfinvf(iel,2) if ncon=9
c                                     lqfinvf(iel,3) if ncon=5
c                                     lqfinvf(iel,4) if ncon=4
c
c iuserdef         : 6th parameter on the 1st line of *SET_POROUS_LAGRANGIAN
c
c ncon             : part ID + number of nodes connected to the elements
c                    ncon=9 : solids
c                    ncon=5 : shells
c                    ncon=4 : beams
c
c ix(j,iel)        : ix(1,iel) : part ID for the element iel
c                    ix(2,iel) : 1st node of the element iel
c                    .
c                    .
c                    .
c                    ix(ncon,iel) : (ncon-1)th node of the element iel
c
c xloc(j)          : x-axis of the local coordinate system attached to iel 
c yloc(j)          : y-axis of the local coordinate system attached to iel 
c zloc(j)          : z-axis of the local coordinate system attached to iel 
c
c locaux(iel)           : pointer to the auxiliary variables in the vector auxvec
c auxvec(locaux(iel)+j) : auxiliary variables of the element iel          
c                         auxvec(locaux(iel)+1) : xx-stress
c                         auxvec(locaux(iel)+2) : yy-stress
c                         auxvec(locaux(iel)+3) : zz-stress
c                         auxvec(locaux(iel)+4) : xy-stress
c                         auxvec(locaux(iel)+5) : yz-stress
c                         auxvec(locaux(iel)+6) : xz-stress
c                         auxvec(locaux(iel)+7) : plastic strain
c The history variables depends on the material (see histlist.txt)
c
c ro(ix(1,iel))    : initial density
c
c cm               : list of keyword inputs from *NAT
c
c x(3,nd)          : position of the node nd (double precision)
c
c v(3,nd)          : velocity of the node nd
c
c a(3,nd)          : acceleration of the node nd
c
c xms(nd)          : inverse of the nodal mass
c
c tt               : time
c
c ncycle           : current cycle
c
c tmpnd(nd)        : nodal temperature
      do j=1,18
        ab(j)=0.0
      enddo
      return
      end
c
c *MAT_ALE_VISCOUS: user defined dynamic viscosity
c
      subroutine f3dm9ale_userdef1(xmu,iel,gp,
     .                             ncon,ix,
     .                             auxvec,locaux,
     .                             ro,cm,
     .                             x,v,a,
     .                             xms,tt,ncycle,tmpnd)
      real*8 x
      dimension ix(ncon,*),auxvec(*),locaux(*),ro(*),cm(*),
     .          x(3,*),v(3,*),a(3,*),xms(*),tmpnd(*)
c
c xmu              : dynamic viscosity
c
c iel              : internal Lagrangian element ID (structure element) 
c                    get user ID with lqfinvf(iel,2) if ncon=9
c                                     lqfinvf(iel,3) if ncon=5
c                                     lqfinvf(iel,4) if ncon=4
c
c gp               : material group rank (defined by *ALE_MULTI-MATERIAL_GROUP)
c
c ncon             : part ID + number of nodes connected to the elements
c                    ncon=9 : solids
c                    ncon=5 : shells
c                    ncon=4 : beams
c
c ix(j,iel)        : ix(1,iel) : part ID for the element iel
c                    ix(2,iel) : 1st node of the element iel
c                    .
c                    .
c                    .
c                    ix(ncon,iel) : (ncon-1)th node of the element iel
c
c locaux(iel)           : pointer to the auxiliary variables in the vector auxvec
c                         for the group gp
c auxvec(locaux(iel)+j) : auxiliary variables of the element iel          
c                         auxvec(locaux(iel)+1) : xx-stress for the group gp
c                         auxvec(locaux(iel)+2) : yy-stress for the group gp
c                         auxvec(locaux(iel)+3) : zz-stress for the group gp
c                         auxvec(locaux(iel)+4) : xy-stress for the group gp
c                         auxvec(locaux(iel)+5) : yz-stress for the group gp
c                         auxvec(locaux(iel)+6) : xz-stress for the group gp
c                         auxvec(locaux(iel)+7) : plastic strain for the group gp
c                         auxvec(locaux(iel)+8) : internal energy for the group gp
c                         auxvec(locaux(iel)+9) : bulk viscosity for the group gp
c                         auxvec(locaux(iel)+10): previous material volume for the group gp
c
c ro(ix(1,iel))    : initial density
c
c cm               : list of keyword inputs from *NAT_ALE_VISCOUS
c                    cm(1) : PC
c                    cm(2) : MULO
c                    cm(3) : -IUSERDEF
c                    cm(4) : RK
c                    cm(6) : RT
c                    cm(7) : RN
c
c x(3,nd)          : position of the node nd (double precision)
c
c v(3,nd)          : velocity of the node nd
c
c a(3,nd)          : acceleration of the node nd
c
c xms(nd)          : inverse of the nodal mass
c
c tt               : time
c
c ncycle           : current cycle
c
c tmpnd(nd)        : nodal temperature
      return
      end
      
