      subroutine umat47(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt)
c Isotropic bilinear material User Materia Subroutine
c Variables
c WRITTEN BY N.L.Vishnuvardhan Raju
c cm(1)=first material constant, here young's modulus
c cm(2)=second material constant, here poisson's ratio
c cm(3)=third material constant,yield stress.
C CM(4)= 4TH MATERIAL CONSTANT %REDUCTION IN YOUNGS MODULUS
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
      REAL :: CC,SIG,EPS,NCC
      DIMENSION CC(6,6),NCC(6,6)
      dimension cm(*),eps(*),sig(*),hsv(*)
      real tt,dt1,efs
      real capa,epsp
      REAL E,NU,G,PRD,NR1,NR2,D,NE,NG,YD,TOL
      LOGICAL YLD
      character*5 etype
       E = CM(1)
       NU = CM(2)
       YD = CM(3)
       D = CM(4)
       
        G =E/(2.0*(1.0+NU))
        PRD = 1-NU-2.0*NU**2
        NR1 = E*(1-NU)
        NR2 = E*NU

        TOL = 1.0E-06
        
        CC =  0.0
        
        CC(1,1) = NR1/PRD
        CC(2,2) = NR1/PRD
        CC(3,3) = NR1/PRD
        CC(1,2) = NR2/PRD
        CC(1,3) = NR2/PRD
        CC(2,1) = NR2/PRD
        CC(2,3) = NR2/PRD
        CC(3,1) = NR2/PRD
        CC(3,2) = NR2/PRD
        CC(4,4) = G
        CC(5,5) = G
        CC(6,6) = G
        
        CALL BLIN(E,NU,D,NCC,NE,NG)
        
      CALL YLDCC(SIG,YD,YLD,TOL)
       IF (YLD .EQV. .FALSE.) THEN
       DO 15,I=1,6
         DO 20,J=1,6
          SIG(I)=SIG(I) + CC(I,J)*EPS(J)
20      CONTINUE
15      CONTINUE

        ELSE
        DO 25,I=1,6
         DO 30,J=1,6
          SIG(I)=SIG(I) + CC(I,J)*EPS(J)
30      CONTINUE
25      CONTINUE
        EFS = -(EPS(1)+EPS(2)+EPS(3))
        HSV(1)= -(1.0/3.0)*(SIG(1)+SIG(2)+SIG(3))
        HSV(2)=  HSV(2)+EFS
        EPSP = HSV(2)
        END IF
        END
C
        SUBROUTINE BLIN(E,NU,D,NCC,NE,NG)
         IMPLICIT NONE
C       THIS BLIN IS USED AFTER THE YIELD STRESS IS CROSSED 
c       FOR CALCULATING THE UPDATED STIFFNESS MATRIX
         REAL :: NCC
         DIMENSION  NCC(6,6)
         REAL:: NE,NNU,NG,E,NU,D
         REAL:: NG2,NPRD,NNR1,NNR2
C         INTEGER :: I,J

        
        NE  = E*(1-D)
        NNU = NU !*(1+PNU)
        
        NG2 = NE/(1+NNU)
        NG  = NG2*0.5
        NCC = 0.0
        NPRD = 1-NNU-2.0*NNU**2
        NNR1 = NE*(1-NNU)
        NNR2 = NE*NNU
        
        NCC(1,1) = NNR1/NPRD
        NCC(2,2) = NNR1/NPRD
        NCC(3,3) = NNR1/NPRD
        NCC(1,2) = NNR2/NPRD
        NCC(1,3) = NNR2/NPRD
        NCC(2,1) = NNR2/NPRD
        NCC(2,3) = NNR2/NPRD
        NCC(3,1) = NNR2/NPRD
        NCC(3,2) = NNR2/NPRD
        NCC(4,4) = NG
        NCC(5,5) = NG
        NCC(6,6) = NG
        END

        SUBROUTINE YLDCC(SIG,YD,YLD,TOL)
       IMPLICIT NONE
       REAL SIG,YD,SIGV,TOL,S1,S2,S3,S4
       DIMENSION  SIG(6)
       LOGICAL YLD
       S1 = (SIG(1)-SIG(2))**2
       S2 = (SIG(2)-SIG(3))**2
       S3 = (SIG(3)-SIG(1))**2
       S4 = (SIG(4)**2+SIG(5)**2+SIG(6)**2)
       SIGV = SQRT(0.5*(S1+S2+S3+6.0*S4))
       
       IF((ABS(YD)-ABS(SIGV)) .LT. TOL) THEN
       YLD =.TRUE.
       ELSE
       YLD = .FALSE.
       END IF
       END 
