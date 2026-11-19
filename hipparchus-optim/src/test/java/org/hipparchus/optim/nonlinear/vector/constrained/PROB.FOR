C
C*********************************************************************
C     
C               TESTPROBLEMS FOR NONLINEAR OPTIMIZATION
C
C
C   Author:      K. Schittkowski, Department of Computer Science, 
C                University of Bayreuth, D-95440 Bayreuth, Germany
C
C
C   References:  K. Schittkowski (2010): An Updated Set of 306 Test 
C                Problems for Nonlinear Programming with Validated 
C                Optimal Solutions, Report, Department of Computer 
C                Science, University of Bayreuth, Germany
C
C                W. Hock, K. Schittkowski (1981): Test Examples for 
C                Nonlinear Programming Codes, Lecture Notes in Economics 
C                and Mathematical Systems, Vol. 187, Springer
C
C                K. Schittkowski (1987): More Test Examples for 
C                Nonlinear Programming, Lecture Notes in Economics and 
C                Mathematical Systems, Vol. 182, Springer
C
C
C   Last change: Oct 2011 - bounds on variables 
C      
C
C***********************************************************************
C
      SUBROUTINE TP1(MODE)      
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      GOTO (1,2,3,4,4),MODE     
1     N=2       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=0    
      X(1)=-2.D0
      X(2)=1.D0 
      LXL(1)=.FALSE.    
      LXL(2)=.TRUE.     
      LXU(1)=.FALSE.    
      LXU(2)=.FALSE.    
      XL(2)=-1.5D0      
      LEX=.TRUE.
      NEX=1     
      XEX(1)=1.D0       
      XEX(2)=1.D0       
      FEX=0.D0  
      RETURN    
2     FX=100.D0*(X(2)-X(1)**2)**2+(1.D0-X(1))**2
      RETURN    
3     GF(2)=200.D0*(X(2)-X(1)**2)       
      GF(1)=-2.D0*(X(1)*(GF(2)-1.D0)+1.D0)      
4     RETURN    
      END       
C
      SUBROUTINE TP2(MODE)      
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION     X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION W1
      GOTO (1,2,3,4,4),MODE     
1     N=2       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=0    
      X(1)=-2.D0
      X(2)=1.0D0 
      LXL(1)=.FALSE.    
      LXU(1)=.FALSE.    
      LXL(2)=.TRUE.     
      XL(2)=1.5D0       
      LXU(2)=.FALSE.    
      LEX=.TRUE.
      NEX=1     
      W1=DSQRT(598.D0/1200.D0)  
      XEX(1)=2.D0*W1*DCOS(DACOS(2.5D-3/W1**3)/3.D0)    
      XEX(2)=1.5D0      
      FEX=100.D0*(XEX(2)-XEX(1)**2)**2+(1.D0-XEX(1))**2 
      RETURN    
2     FX=100.D0*(X(2)-X(1)**2)**2+(1.D0-X(1))**2
      RETURN    
3     GF(2)=200.D0*(X(2)-X(1)**2)       
      GF(1)=-2.D0*(X(1)*(GF(2)-1.D0)+1.D0)      
4     RETURN    
      END       
C
      SUBROUTINE TP3(MODE)      
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      GOTO (1,2,3,4,4),MODE     
1     N=2       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=0    
      X(1)=10.D0
      X(2)=1.D0 
      LXL(1)=.FALSE.    
      LXL(2)=.TRUE.     
      LXU(1)=.FALSE.    
      LXU(2)=.FALSE.    
      XL(2)=0.D0
      LEX=.TRUE.
      XEX(1)=0.D0       
      XEX(2)=0.D0       
      FEX=0.D0  
      NEX=1     
      RETURN    
2     FX=X(2)+(X(2)-X(1))**2*1.D-5      
      RETURN    
3     GF(1)=-2.D0*(X(2)-X(1))*1.D-5     
      GF(2)=1.D0-GF(1)  
4     RETURN    
      END       
C
      SUBROUTINE TP4(MODE)      
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,4),MODE     
1     N=2       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=0    
      X(1)=1.125D0      
      X(2)=0.125D0      
      DO 6 I=1,2
      LXU(I)=.FALSE.    
6     LXL(I)=.TRUE.     
      XL(1)=1.D0
      XL(2)=0.D0
      LEX=.TRUE.
      NEX=1     
      XEX(1)=1.D0       
      XEX(2)=0.D0       
      FEX=8.D0/3.D0     
      GF(2)=1.D0
      RETURN    
2     FX=(X(1)+1.D0)**3/3.D0+X(2)       
      RETURN    
3     GF(1)=(X(1)+1.D0)**2      
4     RETURN    
      END       
C
      SUBROUTINE TP5(MODE)      
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION A,V1,V2      
      GOTO (1,2,3,4,4),MODE     
1     N=2       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=0    
      X(1)=0.D0 
      X(2)=0.D0 
      DO 6 I=1,2
      LXL(I)=.TRUE.     
6     LXU(I)=.TRUE.     
      XL(1)=-1.5D0      
      XL(2)=-3.D0       
      XU(1)=4.D0
      XU(2)=3.D0
      A=4.D0*DATAN(1.D0)
      LEX=.TRUE.
      NEX=1     
      XEX(1)=0.5D0-A/3.D0       
      XEX(2)=XEX(1)-1.D0
      FEX=-DSQRT(3.D0)/2.D0-A/3.D0      
      RETURN    
2     FX=DSIN(X(1)+X(2))+(X(1)-X(2))**2-1.5D0*X(1)+2.5D0*X(2)+1.D0      
      RETURN    
3     V1=DCOS(X(1)+X(2))
      V2=2.D0*(X(1)-X(2))       
      GF(1)=V1+V2-1.5D0 
      GF(2)=V1-V2+2.5D0 
4     RETURN    
      END       
C
      SUBROUTINE TP6(MODE)      
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=1    
      X(1)=-1.2D0       
      X(2)=1.0D0 
      DO 6 I=1,2
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      LEX=.TRUE.
      NEX=1     
      XEX(1)=1.D0       
      XEX(2)=1.D0       
      FEX=0.D0  
      GG(1,2)=10.D0     
      GF(2)=0.D0
      RETURN    
2     FX=(1.D0-X(1))**2 
      RETURN    
3     GF(1)=2.D0*X(1)-2.D0      
      RETURN    
4     IF (INDEX1(1)) G(1)=10.D0*(X(2)-X(1)**2)  
      RETURN    
5     IF (INDEX2(1)) GG(1,1)=-20.D0*X(1)
      RETURN    
      END       
C
      SUBROUTINE TP7(MODE)      
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=1    
      X(1)=2.D0 
      X(2)=2.D0 
      DO 6 I=1,2
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      LEX=.TRUE.
      FEX=-DSQRT(3.D0)  
      XEX(1)=0.D0       
      XEX(2)=-FEX       
      NEX=1     
      GF(2)=-1.D0       
      RETURN    
2     FX=DLOG(1.D0+X(1)**2)-X(2)
      RETURN    
3     GF(1)=2.D0*X(1)/(1.D0+X(1)**2)    
      RETURN    
4     IF (INDEX1(1)) G(1)=(1.D0+X(1)**2)**2+X(2)**2-4.D0
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=4.D0*X(1)*(1.D0+X(1)**2)  
      GG(1,2)=2.D0*X(2) 
7     RETURN    
      END       
C
      SUBROUTINE TP8(MODE)      
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION A,B
      GOTO (1,3,3,4,5),MODE     
1     N=2       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=2    
      X(1)=2.D0 
      X(2)=1.D0 
      A=DSQRT((25.D0+DSQRT(301.D0))/2.D0)       
      B=DSQRT((25.D0-DSQRT(301.D0))/2.D0)       
      DO 6 I=1,2
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      LEX=.TRUE.
      NEX=4     
      XEX(1)=A  
      XEX(2)=9.D0/A     
      XEX(5)=B  
      XEX(6)=9.D0/B     
      DO 30 I=3,7,4     
      DO 30 J=1,2       
30    XEX(I+J-1)=-XEX(I+J-3)    
      FEX=-1.D0 
      GF(1)=0.D0
      GF(2)=0.D0
      FX=-1.D0  
3     RETURN    
4     IF (INDEX1(1)) G(1)=X(1)**2+X(2)**2-25.D0 
      IF (INDEX1(2)) G(2)=X(1)*X(2)-9.D0
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=2.D0*X(1) 
      GG(1,2)=2.D0*X(2) 
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,1)=X(2)      
      GG(2,2)=X(1)      
8     RETURN    
      END       
C
      SUBROUTINE TP9(MODE)      
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION V,V1,V2,V3,V4
      V=4.D0*DATAN(1.D0)
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=0    
      NINL=0    
      NELI=1    
      NENL=0    
      X(1)=0.D0 
      X(2)=0.D0 
      DO 6 I=1,2
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      LEX=.TRUE.
      NEX=1    
      FEX=-0.5D0
      XEX(1)=-3.D0      
      XEX(2)=-4.D0      
      GG(1,1)=4.D0      
      GG(1,2)=-3.D0     
      RETURN    
2     FX=DSIN(V*X(1)/12.D0)*DCOS(V*X(2)/16.D0)  
      RETURN    
3     V3=V/12.D0
      V4=V/16.D0
      V1=V3*X(1)
      V2=V4*X(2)
      GF(1)=V3*DCOS(V1)*DCOS(V2)
      GF(2)=-V4*DSIN(V1)*DSIN(V2)       
      RETURN    
4     IF(INDEX1(1)) G(1)=4.D0*X(1)-3.D0*X(2)    
5     RETURN    
      END       
C
      SUBROUTINE TP10(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=0    
      NINL=1    
      NELI=0    
      NENL=0    
      X(1)=-10.D0       
      X(2)=10.D0
      DO 6 I=1,2
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      LEX=.TRUE.
      NEX=1     
      XEX(1)=0.D0       
      XEX(2)=1.D0       
      FEX=-1.D0 
      GF(1)=1.D0
      GF(2)=-1.D0       
      RETURN    
2     FX=X(1)-X(2)      
3     RETURN    
4     IF (INDEX1(1)) G(1)=-3.D0*X(1)**2+2.D0*X(1)*X(2)-X(2)**2+1.D0     
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=-6.D0*X(1)+2.D0*X(2)      
      GG(1,2)=2.D0*(X(1)-X(2))  
7     RETURN    
      END       
C
      SUBROUTINE TP11(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION AEX,AW,QAW     
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=0    
      NINL=1    
      NELI=0    
      NENL=0    
      DO 6 I=1,2
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      X(1)=4.9D0
      X(2)=0.1D0
      LEX=.TRUE.
      NEX=1     
      AEX=7.5D0*DSQRT(6.D0)     
      AW=(DSQRT(AEX**2+1.D0)+AEX)**(1.D0/3.D0)  
      QAW=AW**2 
      XEX(1)=(AW-1.D0/AW)/DSQRT(6.D0)   
      XEX(2)=(QAW-2.D0+1.D0/QAW)/6.D0   
      FEX=(XEX(1)-5.D0)**2+XEX(2)**2-25.D0      
      GG(1,2)=1.D0      
      RETURN    
2     FX=(X(1)-5.D0)**2+X(2)**2-25.D0   
      RETURN    
3     GF(1)=2.D0*(X(1)-5.D0)    
      GF(2)=2.D0*X(2)   
      RETURN    
4     IF (INDEX1(1)) G(1)=-X(1)**2+X(2) 
      RETURN    
5     IF (INDEX2(1)) GG(1,1)=-2.D0*X(1) 
      RETURN    
      END       
C
      SUBROUTINE TP12(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=0    
      NINL=1    
      NELI=0    
      NENL=0    
      X(1)=0.D0 
      X(2)=0.D0 
      DO 6 I=1,2
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      LEX=.TRUE.
      NEX=1     
      XEX(1)=2.D0       
      XEX(2)=3.D0       
      FEX=-30.D0
      RETURN    
2     FX=0.5D0*X(1)**2+X(2)**2-X(1)*X(2)-7.D0*X(1)-7.D0*X(2)    
      RETURN    
3     GF(1)=X(1)-X(2)-7.D0      
      GF(2)=2.D0*X(2)-X(1)-7.D0 
      RETURN    
4     IF (INDEX1(1)) G(1)=25.D0-4.D0*X(1)**2-X(2)**2    
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=-8.D0*X(1)
      GG(1,2)=-2.D0*X(2)
7     RETURN    
      END       
C
      SUBROUTINE TP13(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=0    
      NINL=1    
      NELI=0    
      NENL=0    
      X(1)=0.D0
      X(2)=0.D0
      DO 6 I=1,2
      LXU(I)=.FALSE.    
      LXL(I)=.TRUE.     
6     XL(I)=0.D0
      LEX=.TRUE.
      NEX=1     
      XEX(1)=1.D0       
      XEX(2)=0.D0       
      FEX=1.D0  
      GG(1,2)=-1.D0     
      RETURN    
2     FX=(X(1)-2.D0)**2+X(2)**2 
      RETURN    
3     GF(1)=2.D0*(X(1)-2.D0)    
      GF(2)=2.D0*X(2)   
      RETURN    
4     IF (INDEX1(1)) G(1)=(1.D0-X(1))**3-X(2)   
      RETURN    
5     IF (INDEX2(1)) GG(1,1)=-3.D0*(1.D0-X(1))**2       
      RETURN    
      END       
C
      SUBROUTINE TP14(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION W7   
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=0    
      NINL=1    
      NELI=1    
      NENL=0    
      X(1)=2.D0 
      X(2)=2.D0 
      DO 6 I=1,2
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      LEX=.TRUE.
      NEX=1     
      W7=DSQRT(7.D0)    
      XEX(1)=(W7-1.D0)*0.5D0    
      XEX(2)=(W7+1.D0)*0.25D0   
      FEX=9.D0-23.D0*W7/8.D0    
      GG(2,1)=1.D0      
      GG(2,2)=-2.D0     
      RETURN    
2     FX=(X(1)-2.D0)**2+(X(2)-1.D0)**2  
      RETURN    
3     GF(1)=2.D0*(X(1)-2.D0)    
      GF(2)=2.D0*(X(2)-1.D0)    
      RETURN    
4     IF (INDEX1(1)) G(1)=1.D0-(X(1)**2)*0.25D0-X(2)**2 
      IF (INDEX1(2)) G(2)=X(1)-2.D0*X(2)+1.D0   
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=-X(1)*0.5D0       
      GG(1,2)=-2.D0*X(2)
7     RETURN    
      END       
C
      SUBROUTINE TP15(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=0    
      NINL=2    
      NELI=0    
      NENL=0    
      X(1)=-2.D0
      X(2)=1.D0 
      LXL(1)=.FALSE.    
      LXL(2)=.FALSE.    
      LXU(1)=.TRUE.     
      LXU(2)=.FALSE.    
      XU(1)=0.5D0       
      LEX=.TRUE.
      XEX(1)=0.5D0      
      XEX(2)=2.0       
      FEX=3.065D0      
      NEX=1     
      GG(2,1)=1.D0      
      RETURN    
2     FX=(X(2)-X(1)**2)**2+0.01*(1.D0-X(1))**2
      RETURN    
3     GF(2)=2.0*(X(2)-X(1)**2)       
      GF(1)=-2.D-2*(X(1)*(GF(2)-1.D0)+1.D0)      
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)*X(2)-1.D0
      IF (INDEX1(2)) G(2)=X(2)**2+X(1)  
      RETURN    
5     IF(.NOT.INDEX2(1)) GOTO 7 
      GG(1,1)=X(2)      
      GG(1,2)=X(1)      
7     IF (INDEX2(2)) GG(2,2)=2.D0*X(2)  
      RETURN    
      END       
C
      SUBROUTINE TP16(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=0    
      NINL=2    
      NELI=0    
      NENL=0    
      X(1)=-2.D0
      X(2)=1.D0 
      LXL(1)=.TRUE.     
      LXL(2)=.FALSE.    
      LXU(1)=.TRUE.     
      LXU(2)=.TRUE.     
      XL(1)=-0.5D0      
      XU(1)=0.5D0       
      XU(2)=1.D0
      LEX=.TRUE.
      NEX=1     
      XEX(1)=0.5D0      
      XEX(2)=0.25D0     
      FEX=0.25D0
      GG(1,1)=1.D0      
      GG(2,2)=1.D0      
      RETURN    
2     FX=100.D0*(X(2)-X(1)**2)**2+(1.D0-X(1))**2
      RETURN    
3     GF(2)=200.D0*(X(2)-X(1)**2)       
      GF(1)=-2.D0*(X(1)*(GF(2)-1.D0)+1.D0)      
      RETURN    
4     IF (INDEX1(1)) G(1)=X(2)**2+X(1)  
      IF (INDEX1(2)) G(2)=X(1)**2+X(2)  
      RETURN    
5     IF (INDEX2(1)) GG(1,2)=2.D0*X(2)  
      IF (INDEX2(2)) GG(2,1)=2.D0*X(1)  
      RETURN    
      END       
C
      SUBROUTINE TP17(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=0    
      NINL=2    
      NELI=0    
      NENL=0    
      X(1)=-2.0
      X(2)=1.D0 
      LXL(1)=.TRUE.     
      LXL(2)=.FALSE.    
      LXU(1)=.TRUE.     
      LXU(2)=.TRUE.     
      XL(1)=-2.0      
      XU(1)=0.5D0       
      XU(2)=1.D0
      LEX=.TRUE.
      NEX=1     
      XEX(1)=0.D0       
      XEX(2)=0.D0       
      FEX=1.D0  
      GG(1,1)=-1.D0     
      GG(2,2)=-1.D0     
      RETURN    
2     FX=(100.D0*(X(2)-X(1)**2)**2+(1.D0-X(1))**2)
      RETURN    
3     GF(2)=200.D0*(X(2)-X(1)**2)      
      GF(1)=-2.D0*(X(1)*(GF(2)-1.D0)+1.D0)   
      RETURN    
4     IF (INDEX1(1)) G(1)=X(2)**2-X(1)  
      IF(INDEX1(2)) G(2)=X(1)**2-X(2)   
      RETURN    
5     IF (INDEX2(1)) GG(1,2)=2.D0*X(2)  
      IF (INDEX2(2)) GG(2,1)=2.D0*X(1)  
      RETURN    
      END       
C
      SUBROUTINE TP18(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=0    
      NINL=2    
      NELI=0    
      NENL=0    
      X(1)=2.D0 
      X(2)=2.D0 
      DO 6 I=1,2
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
6     XU(I)=50.D0       
      XL(1)=2.D0
      XL(2)=0.D0
      LEX=.TRUE.
      NEX=1     
      XEX(1)=DSQRT(250.D0)      
      XEX(2)=0.1D0*XEX(1)       
      FEX=5.D0  
      RETURN    
2     FX=0.01D0*X(1)**2+X(2)**2 
      RETURN    
3     GF(1)=0.02D0*X(1) 
      GF(2)=2.D0*X(2)   
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)*X(2)-25.D0       
      IF(INDEX1(2)) G(2)=X(1)**2+X(2)**2-25.D0  
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=X(2)      
      GG(1,2)=X(1)      
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,1)=2.D0*X(1) 
      GG(2,2)=2.D0*X(2) 
8     RETURN    
      END       
C
      SUBROUTINE TP19(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /L19/ AEX,SAEX 
      DOUBLEPRECISION AEX,SAEX
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=0    
      NINL=2    
      NELI=0    
      NENL=0    
      X(1)=20.1D0       
      X(2)=5.84D0       
      DO 6 I=1,2
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
6     XU(I)=100.D0      
      XL(1)=13.D0       
      XL(2)=0.D0
      LEX=.TRUE.
      NEX=1     
      SAEX=1.7280975D+1 
      AEX=DSQRT(SAEX)   
      XEX(1)=14.095D0   
      XEX(2)=5.D0-AEX   
      FEX=(4.095D0**3-(15.D0+AEX)**3)
      RETURN    
2     FX=((X(1)-10.D0)**3+(X(2)-20.D0)**3)
      RETURN    
3     GF(1)=(3.D0*(X(1)-10.D0)**2)
      GF(2)=(3.D0*(X(2)-20.D0)**2)
      RETURN    
4     IF (INDEX1(1)) G(1)=(X(1)-5.D0)**2+(X(2)-5.D0)**2-100.D0  
      IF (INDEX1(2)) G(2)=82.81D0-(X(1)-6.D0)**2-(X(2)-5.D0)**2 
      RETURN    
5     IF(.NOT.INDEX2(1)) GOTO 7 
      GG(1,1)=2.D0*(X(1)-5.D0)  
      GG(1,2)=2.D0*(X(2)-5.D0)  
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,1)=-2.D0*(X(1)-6.D0) 
      GG(2,2)=-2.D0*(X(2)-5.D0) 
8     RETURN    
      END      
C       
      SUBROUTINE TP20(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=0    
      NINL=3    
      NELI=0    
      NENL=0    
      X(1)=1.D-1
      X(2)=1.D0 
      LXL(1)=.TRUE.     
      LXL(2)=.FALSE.    
      LXU(1)=.TRUE.     
      LXU(2)=.FALSE.    
      XL(1)=-0.5D0      
      XU(1)=0.5D0       
      LEX=.TRUE.
      NEX=1     
      XEX(1)=0.5D0      
      XEX(2)=DSQRT(3.D0)*0.5D0  
      FEX=(81.5D0-25.D0*DSQRT(3.D0))
      GG(1,1)=1.D0      
      GG(2,2)=1.D0      
      RETURN    
2     FX=(100.D0*(X(2)-X(1)**2)**2+(1.D0-X(1))**2)
      RETURN    
3     GF(2)=200.D0*(X(2)-X(1)**2)
      GF(1)=-2.D0*(X(1)*(GF(2)-1.D0)+1.D0)
      RETURN    
4     IF (INDEX1(1)) G(1)=X(2)**2+X(1)  
      IF (INDEX1(2)) G(2)=X(1)**2+X(2)  
      IF  (INDEX1(3)) G(3)=X(1)**2+X(2)**2-1.D0 
      RETURN    
5     IF (INDEX2(1)) GG(1,2)=2.D0*X(2)  
      IF(INDEX2(2)) GG(2,1)=2.D0*X(1)   
      IF (.NOT.INDEX2(3)) GOTO 9
      GG(3,1)=2.D0*X(1) 
      GG(3,2)=2.D0*X(2) 
9     RETURN    
      END       
C
      SUBROUTINE TP21(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=1    
      NINL=0    
      NELI=0    
      NENL=0    
      X(1)=2.D0
      X(2)=-1.D0
      DO 6 I=1,2
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
6     XU(I)=50.D0       
      XL(1)=2.D0
      XL(2)=-50.D0      
      LEX=.TRUE.
      NEX=1     
      XEX(1)=2.D0       
      XEX(2)=0.D0       
      FEX=-99.96D0      
      GG(1,1)=10.D0     
      GG(1,2)=-1.D0     
      RETURN    
2     FX=(0.01*X(1)**2 + X(2)**2 - 100.D0)
      RETURN    
3     GF(1)=0.02*X(1)
      GF(2)=2.0*X(2)
      RETURN    
4     IF (INDEX1(1)) G(1)=10.D0*X(1)-X(2)-10.D0 
5     RETURN    
      END       
C
      SUBROUTINE TP22(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=1    
      NINL=1    
      NELI=0    
      NENL=0    
      X(1)=2.D0 
      X(2)=2.D0 
      DO 6 I=1,2
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      LEX=.TRUE.
      NEX=1     
      XEX(1)=1.D0       
      XEX(2)=1.D0       
      FEX=1.D0  
      GG(1,1)=-1.D0     
      GG(1,2)=-1.D0     
      GG(2,2)=1.D0      
      RETURN    
2     FX=(X(1)-2.D0)**2+(X(2)-1.D0)**2  
      RETURN    
3     GF(1)=2.D0*(X(1)-2.D0)    
      GF(2)=2.D0*(X(2)-1.D0)    
      RETURN    
4     IF (INDEX1(1)) G(1)=2.D0-X(1)-X(2)
      IF (INDEX1(2)) G(2)=X(2)-X(1)**2  
      RETURN    
5     IF (INDEX2(2)) GG(2,1)=-2.D0*X(1) 
      RETURN    
      END       
C
      SUBROUTINE TP23(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=1    
      NINL=4    
      NELI=0    
      NENL=0    
      X(1)=3.D0 
      X(2)=1.D0 
      DO 6 I=1,2
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      XL(I)=-50.D0      
6     XU(I)=50.D0       
      LEX=.TRUE.
      NEX=1     
      XEX(1)=1.D0       
      XEX(2)=1.D0       
      FEX=2.D0  
      GG(1,1)=1.D0      
      GG(1,2)=1.D0      
      GG(4,2)=-1.D0     
      GG(5,1)=-1.D0     
      RETURN    
2     FX=X(1)**2+X(2)**2
      RETURN    
3     GF(1)=2.D0*X(1)   
      GF(2)=2.D0*X(2)   
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)+X(2)-1.D0
      IF (INDEX1(2)) G(2)=X(1)**2+X(2)**2-1.D0  
      IF (INDEX1(3)) G(3)=9.D0*X(1)**2+X(2)**2-9.D0     
      IF (INDEX1(4)) G(4)=X(1)**2-X(2)  
      IF (INDEX1(5)) G(5)=X(2)**2-X(1)  
      RETURN    
5     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,1)=2.D0*X(1) 
      GG(2,2)=2.D0*X(2) 
8     IF (.NOT.INDEX2(3)) GOTO 9
      GG(3,1)=18.D0*X(1)
      GG(3,2)=2.D0*X(2) 
9     IF (INDEX2(4)) GG(4,1)=2.D0*X(1)  
      IF (INDEX2(5)) GG(5,2)=2.D0*X(2)  
      RETURN    
      END       
C
      SUBROUTINE TP24(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /L24/ A
      DOUBLEPRECISION A    
      INTEGER I
      A=DSQRT(3.0D0)     
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=3    
      NINL=0    
      NELI=0    
      NENL=0    
      X(1)=1.D0 
      X(2)=0.5D0
      DO 6 I=1,2
      LXL(I)=.TRUE.     
      LXU(I)=.FALSE.    
6     XL(I)=0.D0
      LEX=.TRUE.
      NEX=1     
      XEX(1)=3.D0       
      XEX(2)=A  
      FEX=-1.D0 
      GG(1,1)=1.D0/A    
      GG(1,2)=-1.D0     
      GG(2,1)=1.D0      
      GG(2,2)=A 
      GG(3,1)=-1.D0     
      GG(3,2)=-A
      RETURN    
2     FX=((X(1)-3.D0)**2-9.D0)*X(2)**3/(27.D0*A)
      RETURN    
3     GF(1)=2.D0*(X(1)-3.D0)*X(2)**3/(27.D0*A)  
      GF(2)=((X(1)-3.D0)**2-9.D0)*X(2)**2/(9.D0*A)      
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)/A-X(2)   
      IF (INDEX1(2)) G(2)=X(1)+X(2)*A   
      IF (INDEX1(3)) G(3)=6.D0-X(2)*A-X(1)      
5     RETURN    
      END       
C
      SUBROUTINE TP25(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D25/ A(99),B(99),U(99),DA(99,3),V1,T,S,V2,V11,V22,X13
      DOUBLEPRECISION A,B,U,DA,V1,T,S,V2,V11,V22,X13    
      INTEGER I,J
      GOTO (1,2,3,4,4),MODE     
1     N=3       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=0    
      X(1)=100.D0       
      X(2)=12.5D0       
      X(3)=3.D0 
      DO 6 I=1,3
      LXL(I)=.TRUE.     
6     LXU(I)=.TRUE.     
      XL(1)=0.1D0       
      XL(2)=1.D-5       
      XL(3)=1.D-5       
      XU(1)=100.D0      
      XU(2)=25.6D0      
      XU(3)=5.0D0
      LEX=.TRUE.
      NEX=1     
      XEX(1)=50.D0      
      XEX(2)=25.D0      
      XEX(3)=1.5D0      
      FEX=0.D0  
      RETURN    
  2   CONTINUE  
      DO 30 I=1,99      
      V1=2.D0/3.D0      
      U(I)=25.D0+(-50.D0*DLOG(0.01D0*DBLE(I)))**V1    
      V1=U(I)-X(2) 
      IF (V1.LT.0) GOTO 7       
      V11=-(V1**X(3))/X(1)
      B(I)=DEXP(V11)    
30    A(I)=B(I)-0.01D0*DBLE(I)
      T=0.D0    
      DO 31 I=1,99      
31    T=T+A(I)**2       
      FX=T
      RETURN    
7     S=0.D0    
      DO 8 I=1,3
8     S=S+(X(I)-5.D0)**2
      FX=S
      RETURN    
 3    CONTINUE  
      DO 36 I=1,99      
      V1=2.D0/3.D0      
      U(I)=25.D0-(50.D0*DLOG(0.01D0*DBLE(I)))**V1    
      V2=U(I)-X(2)      
      IF (V2.LE.0) GOTO 9       
      V22=-V2**X(3)/X(1)
      IF(V22.GT. -150.D0) GOTO42
      B(I)=0.D0 
      GOTO43    
   42 B(I)=DEXP(V22)    
   43 CONTINUE  
      A(I)=B(I)-0.01D0*DBLE(I)
      DA(I,1)=V2**X(3)/X(1)**2*B(I)     
      DA(I,2)=X(3)*V2**(X(3)-1.D0)/X(1)*B(I)    
36    DA(I,3)=-V2**X(3)/X(1)*DLOG(V2)*B(I)      
      DO 34 I=1,3       
      T=0.D0    
      DO 33 J=1,99      
33    T=T+2.D0*A(J)*DA(J,I)     
34    GF(I)=T   
      RETURN    
9     DO 10 I=1,3       
10    GF(I)=2.D0*(X(I)-5.D0)    
4     RETURN    
      END       
C
      SUBROUTINE TP26(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION A    
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=1    
      X(1)=-2.6D0       
      X(2)=2.D0 
      X(3)=2.D0 
      DO 6 I=1,3
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      LEX=.TRUE.
      NEX=1     
      XEX(1)=1.D0       
      XEX(2)=1.D0       
      XEX(3)=1.D0       
      A=DSQRT(139.D0/108.D0)    
      XEX(4)=(A-61.D0/54.D0)**(1.D0/3.D0)-(61.D0/54.D0+A)**(1.D0/3.D0)  
     /     -2.D0/3.D0   
      XEX(5)=XEX(4)     
      XEX(6)=XEX(4)     
      FEX=0.D0  
      RETURN    
2     FX=(X(1)-X(2))**2+(X(2)-X(3))**4  
      RETURN    
3     GF(1)=2.D0*(X(1)-X(2))    
      GF(3)=-4.D0*(X(2)-X(3))**3
      GF(2)=-GF(1)-GF(3)
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)*(1.D0+X(2)**2)+X(3)**4-3.D0      
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=1.D0+X(2)**2      
      GG(1,2)=2.D0*X(1)*X(2)    
      GG(1,3)=4.D0*X(3)**3      
7     RETURN    
      END       
C
      SUBROUTINE TP27(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=1    
      DO 6 I=1,3
      X(I)=2.D0 
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      LEX=.TRUE.
      NEX=1     
      XEX(1)=-1.D0      
      XEX(2)=1.D0       
      XEX(3)=0.D0       
      FEX=4.0
      GF(3)=0.D0
      GG(1,1)=1.D0      
      GG(1,2)=0.D0      
      RETURN    
2     FX=(X(1)-1.0D0)**2 + 100.0D0*(X(2)-X(1)**2)**2
      RETURN    
3     GF(1)=(X(1)-1.0D0)*2.0D0 - 400.0D0*(X(2)-X(1)**2)*X(1)
      GF(2)=2.0D2*(X(2)-X(1)**2) 
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)+X(3)**2+1.D0     
      RETURN    
5     IF (INDEX2(1)) GG(1,3)=2.D0*X(3)  
      RETURN    
      END       
C      
      SUBROUTINE TP28(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=0    
      NINL=0    
      NELI=1    
      NENL=0    
      X(1)=-4.D0
      X(2)=1.D0 
      X(3)=1.D0 
      DO 6 I=1,3
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      LEX=.TRUE.
      NEX=1     
      XEX(1)=0.5D0      
      XEX(2)=-0.5D0     
      XEX(3)=0.5D0      
      FEX=0.D0  
      GG(1,1)=1.D0      
      GG(1,2)=2.D0      
      GG(1,3)=3.D0      
      RETURN    
2     FX=(X(1)+X(2))**2+(X(2)+X(3))**2  
      RETURN    
3     GF(1)=2.D0*(X(1)+X(2))    
      GF(3)=2.D0*(X(2)+X(3))    
      GF(2)=GF(1)+GF(3) 
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)+2.D0*X(2)+3.D0*X(3)-1.D0 
5     RETURN    
      END       
C
      SUBROUTINE TP29(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=0    
      NINL=1    
      NELI=0    
      NENL=0    
      DO 6 I=1,3
      X(I)=1.D0 
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      LEX=.TRUE.
      NEX=4     
      XEX(1)=4.D0       
      XEX(2)=2.D0*DSQRT(2.D0)   
      XEX(3)=2.D0       
      XEX(4)=XEX(1)     
      XEX(5)=-XEX(2)    
      XEX(6)=-XEX(3)    
      XEX(7)=-XEX(1)    
      XEX(8)=XEX(2)     
      XEX(9)=-XEX(3)    
      XEX(10)=-XEX(1)   
      XEX(11)=-XEX(2)   
      XEX(12)=XEX(3)    
      FEX=-16.D0*DSQRT(2.D0)    
      RETURN    
2     FX=-X(1)*X(2)*X(3)
      RETURN    
3     GF(1)=-X(2)*X(3)  
      GF(2)=-X(1)*X(3)  
      GF(3)=-X(1)*X(2)  
      RETURN    
4     IF (INDEX1(1)) G(1)=48.D0-X(1)**2-2.D0*X(2)**2-4.D0*X(3)**2       
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=-2.D0*X(1)
      GG(1,2)=-4.D0*X(2)
      GG(1,3)=-8.D0*X(3)
7     RETURN    
      END       
C
      SUBROUTINE TP30(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=0    
      NINL=1    
      NELI=0    
      NENL=0    
      DO 6 I=1,3
      X(I)=1.D0 
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
6     XU(I)=10.D0       
      XL(1)=1.D0
      XL(2)=-10.D0      
      XL(3)=-10.D0      
      LEX=.TRUE.
      NEX=1     
      XEX(1)=1.D0       
      XEX(2)=0.D0       
      XEX(3)=0.D0       
      FEX=1.D0  
      GG(1,3)=0.D0      
      RETURN    
2     FX=X(1)**2+X(2)**2+X(3)**2
      RETURN    
3     GF(1)=2.D0*X(1)   
      GF(2)=2.D0*X(2)   
      GF(3)=2.D0*X(3)   
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)**2+X(2)**2-1.D0  
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=2.D0*X(1) 
      GG(1,2)=2.D0*X(2) 
7     RETURN    
      END       
C
      SUBROUTINE TP31(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=0    
      NINL=1    
      NELI=0    
      NENL=0    
      DO 6 I=1,3
      X(I)=1.D0 
      LXL(I)=.TRUE.     
6     LXU(I)=.TRUE.     
      XL(1)=-10.D0      
      XL(2)=1.D0
      XL(3)=-10.D0      
      XU(1)=10.D0       
      XU(2)=10.D0       
      XU(3)=1.D0
      LEX=.TRUE.
      NEX=1     
      XEX(1)=1.D0/DSQRT(3.D0)   
      XEX(2)=DSQRT(3.D0)
      XEX(3)=0.D0       
      FEX=6.D0  
      GG(1,3)=0.D0      
      RETURN    
2     FX=9.D0*X(1)**2+X(2)**2+9.D0*X(3)**2      
      RETURN    
3     GF(1)=18.D0*X(1)  
      GF(2)=2.D0*X(2)   
      GF(3)=18.D0*X(3)  
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)*X(2)-1.D0
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=X(2)      
      GG(1,2)=X(1)      
7     RETURN    
      END       
C
      SUBROUTINE TP32(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=0    
      NINL=1    
      NELI=1    
      NENL=0    
      X(1)=0.1D0
      X(2)=0.7D0
      X(3)=0.2D0
      DO 6 I=1,3
      LXL(I)=.TRUE.     
      LXU(I)=.FALSE.    
6     XL(I)=0.D0
      LEX=.TRUE.
      NEX=1     
      XEX(1)=0.D0       
      XEX(2)=0.D0       
      XEX(3)=1.D0       
      FEX=1.D0  
      GG(1,2)=6.D0      
      GG(1,3)=4.D0      
      GG(2,1)=-1.D0     
      GG(2,2)=-1.D0     
      GG(2,3)=-1.D0     
      RETURN    
2     FX=(X(1)+3.D0*X(2)+X(3))**2+4.D0*(X(1)-X(2))**2   
      RETURN    
3     GF(1)=10.D0*X(1)-2.D0*X(2)+2.D0*X(3)      
      GF(2)=-2.D0*X(1)+26.D0*X(2)+6.D0*X(3)     
      GF(3)=2.D0*(X(1)+3.D0*X(2)+X(3))  
      RETURN    
4     IF (INDEX1(1)) G(1)=-X(1)**3+6.D0*X(2)+4.D0*X(3)-3.D0     
      IF (INDEX1(2)) G(2)=1.D0-X(1)-X(2)-X(3)   
      RETURN    
5     IF (INDEX2(1)) GG(1,1)=-3.D0*X(1)**2      
      RETURN    
      END       
C
      SUBROUTINE TP33(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=0    
      NINL=2    
      NELI=0    
      NENL=0    
      X(1)=0.D0 
      X(2)=0.D0 
      X(3)=3.D0 
      DO 6 I=1,3
      LXL(I)=.TRUE.     
6     XL(I)=0.D0
      LXU(1)=.FALSE.    
      LXU(2)=.FALSE.    
      LXU(3)=.TRUE.     
      XU(3)=5.D0
      LEX=.TRUE.
      NEX=1     
      XEX(1)=0.D0       
      XEX(2)=DSQRT(2.D0)
      XEX(3)=DSQRT(2.D0)
      FEX=DSQRT(2.D0)-6.0      
      GF(2)=0.D0
      GF(3)=1.D0
      RETURN    
2     FX=(X(1)-1.D0)*(X(1)-2.D0)*(X(1)-3.D0)+X(3)       
      RETURN    
3     GF(1)=3.D0*X(1)**2-12.D0*X(1)+11.D0       
      RETURN    
4     IF (INDEX1(1)) G(1)=X(3)**2-X(1)**2-X(2)**2       
      IF (INDEX1(2)) G(2)=X(1)**2+X(2)**2+X(3)**2-4.D0  
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=-2.D0*X(1)
      GG(1,2)=-2.D0*X(2)
      GG(1,3)=2.D0*X(3) 
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,1)=2.D0*X(1) 
      GG(2,2)=2.D0*X(2) 
      GG(2,3)=2.D0*X(3) 
8     RETURN    
      END       
C
      SUBROUTINE TP34(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=0    
      NINL=2    
      NELI=0    
      NENL=0    
      X(1)=0.D0 
      X(2)=1.05D0       
      X(3)=2.9D0
      DO 6 I=1,3
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
6     XL(I)=0.D0
      XU(1)=100.D0      
      XU(2)=100.D0      
      XU(3)=10.D0       
      LEX=.TRUE.
      NEX=1     
      XEX(1)=DLOG(DLOG(10.D0))  
      XEX(2)=DLOG(10.D0)
      XEX(3)=10.D0      
      FEX=-XEX(1)       
      GF(1)=-1.D0       
      GF(2)=0.D0
      GF(3)=0.D0
      GG(1,2)=1.D0      
      GG(1,3)=0.D0      
      GG(2,1)=0.D0      
      GG(2,3)=1.D0      
      RETURN    
2     FX=-X(1)  
3     RETURN    
4     IF (INDEX1(1)) G(1)=X(2)-DEXP(X(1))       
      IF(INDEX1(2)) G(2)=X(3)-DEXP(X(2))
      RETURN    
5     IF (INDEX2(1)) GG(1,1)=-DEXP(X(1))
      IF (INDEX2(2)) GG(2,2)=-DEXP(X(2))
      RETURN    
      END       
C
      SUBROUTINE TP35(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=1    
      NINL=0    
      NELI=0    
      NENL=0    
      DO 6 I=1,3
      X(I)=0.5D0
      LXL(I)=.TRUE.     
      LXU(I)=.FALSE.    
6     XL(I)=0.D0
      LEX=.TRUE.
      NEX=1     
      XEX(1)=4.D0/3.D0  
      XEX(2)=7.D0/9.D0  
      XEX(3)=4.D0/9.D0  
      FEX=1.D0/9.D0     
      GG(1,1)=-1.D0     
      GG(1,2)=-1.D0     
      GG(1,3)=-2.D0     
      RETURN    
2     FX=9.D0-8.D0*X(1)-6.D0*X(2)-4.D0*X(3)+2.D0*X(1)**2+2.D0*X(2)**2   
     /     +X(3)**2+2.D0*X(1)*X(2)+2.D0*X(1)*X(3)    
      RETURN    
3     GF(1)=-8.D0+4.D0*X(1)+2.D0*X(2)+2.D0*X(3) 
      GF(2)=-6.D0+4.D0*X(2)+2.D0*X(1)   
      GF(3)=-4.D0+2.D0*X(3)+2.D0*X(1)   
      RETURN    
4     IF (INDEX1(1)) G(1)=-X(1)-X(2)-2.D0*X(3)+3.D0     
5     RETURN    
      END       
C
      SUBROUTINE TP36(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=1    
      NINL=0    
      NELI=0    
      NENL=0    
      DO 6 I=1,3
      X(I)=10.D0
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
6     XL(I)=0.D0
      XU(1)=20.D0       
      XU(2)=11.D0       
      XU(3)=42.D0       
      LEX=.TRUE.
      NEX=1     
      XEX(1)=20.D0      
      XEX(2)=11.D0      
      XEX(3)=15.D0      
      FEX=-3.3D+3       
      GG(1,1)=-1.D0     
      GG(1,2)=-2.D0     
      GG(1,3)=-2.D0     
      RETURN    
2     FX=-X(1)*X(2)*X(3)
      RETURN    
3     GF(1)=-X(2)*X(3)  
      GF(2)=-X(1)*X(3)  
      GF(3)=-X(1)*X(2)  
      RETURN    
4     IF (INDEX1(1)) G(1)=72.D0-X(1)-2.D0*X(2)-2.D0*X(3)
      RETURN    
5     RETURN    
      END       
C
      SUBROUTINE TP37(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=2    
      NINL=0    
      NELI=0    
      NENL=0    
      DO 6 I=1,3
      X(I)=10.D0
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      XU(I)=42.D0       
6     XL(I)=0.D0
      LEX=.TRUE.
      NEX=1     
      XEX(1)=24.D0      
      XEX(2)=12.D0      
      XEX(3)=12.D0      
      FEX=-3.456D+3     
      GG(1,1)=-1.D0     
      GG(1,2)=-2.D0     
      GG(1,3)=-2.D0     
      GG(2,1)=1.D0      
      GG(2,2)=2.D0      
      GG(2,3)=2.D0      
      RETURN    
2     FX=-X(1)*X(2)*X(3)
      RETURN    
3     GF(1)=-X(2)*X(3)  
      GF(2)=-X(1)*X(3)  
      GF(3)=-X(1)*X(2)  
      RETURN    
4     IF (INDEX1(1)) G(1)=72.D0-X(1)-2.D0*X(2)-2.D0*X(3)
      IF (INDEX1(2)) G(2)=X(1)+2.D0*X(2)+2.D0*X(3)      
5     RETURN    
      END       
C
      SUBROUTINE TP38(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,4),MODE     
1     N=4       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=0    
      X(1)=-3.D0
      X(2)=-1.D0
      X(3)=-3.D0
      X(4)=-1.D0
      DO 6 I=1,4
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      XL(I)=-10.D0      
6     XU(I)=10.D0       
      LEX=.TRUE.
      NEX=1     
      DO 30 I=1,4       
30    XEX(I)=1.D0       
      FEX=0.D0  
      RETURN    
2     FX=(100.D0*(X(2)-X(1)**2)**2+(1.D0-X(1))**2
     /   + 90.D0*(X(4)-X(3)**2)**2
     /   + (1.D0-X(3))**2+10.1D0*((X(2)-1.D0)**2+(X(4)-1.D0)**2)
     /   + 19.8D0*(X(2)-1.D0)*(X(4)-1.D0))
      RETURN    
3     GF(1)=(-400.D0*X(1)*(X(2)-X(1)**2)-2.D0*(1.D0-X(1)))
      GF(2)=(200.D0*(X(2)-X(1)**2)+20.2D0*(X(2)-1.D0)+19.8D0
     /                 *(X(4)-1.D0))
      GF(3)=(-360.D0*X(3)*(X(4)-X(3)**2)-2.D0*(1.D0-X(3)))
      GF(4)=(180.D0*(X(4)-X(3)**2)+20.2D0*(X(4)-1.D0)+19.8D0
     /                 *(X(2)-1.D0))
4     RETURN    
      END       
C
      SUBROUTINE TP39(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=4       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=2    
      DO 6 I=1,4
      X(I)=2.D0 
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      NEX=1     
      LEX=.TRUE.
      XEX(1)=1.D0       
      XEX(2)=1.D0       
      XEX(3)=0.D0       
      XEX(4)=0.D0       
      FEX=-1.D0 
      GF(1)=-1.D0       
      GF(2)=0.D0
      GF(3)=0.D0
      GF(4)=0.D0
      GG(1,2)=1.D0      
      GG(1,4)=0.D0      
      GG(2,2)=-1.D0     
      GG(2,3)=0.D0      
      RETURN    
2     FX=-X(1)  
3     RETURN    
4     IF (INDEX1(1)) G(1)=X(2)-X(1)**3-X(3)**2  
      IF (INDEX1(2)) G(2)=X(1)**2-X(2)-X(4)**2  
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=-3.D0*X(1)**2     
      GG(1,3)=-2.D0*X(3)
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,1)=2.D0*X(1) 
      GG(2,4)=-2.D0*X(4)
8     RETURN    
      END       
C
      SUBROUTINE TP40(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      GOTO (1,2,3,4,5),MODE     
1     N=4       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=3    
      DO 6 I=1,4
      X(I)=0.8D0
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      DO 15 I=1,3       
      DO 15 J=1,4       
15    GG(I,J)=0.D0      
      GG(2,3)=-1.D0     
      GG(3,2)=-1.D0     
      LEX=.TRUE.
      XEX(1)=2.D0**(-1.D0/3.D0) 
      XEX(2)=2.D0**(-0.5D0)     
      XEX(3)=2.D0**(-11.D0/12.D0)       
      XEX(4)=2.D0**(-0.25D0)    
      XEX(5)=XEX(1)     
      XEX(6)=XEX(2)     
      XEX(7)=-XEX(3)    
      XEX(8)=-XEX(4)    
      FEX=-0.25D0       
      NEX=2     
      RETURN    
2     FX=-X(1)*X(2)*X(3)*X(4)   
      RETURN    
3     GF(1)=-X(2)*X(3)*X(4)     
      GF(2)=-X(1)*X(3)*X(4)     
      GF(3)=-X(1)*X(2)*X(4)     
      GF(4)=-X(1)*X(2)*X(3)     
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)**3+X(2)**2-1.D0  
      IF (INDEX1(2)) G(2)=X(1)**2*X(4)-X(3)     
      IF (INDEX1(3)) G(3)=X(4)**2-X(2)  
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=3.D0*X(1)**2      
      GG(1,2)=2.D0*X(2) 
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,1)=2.D0*X(1)*X(4)    
      GG(2,4)=X(1)**2   
8     IF (.NOT.INDEX2(3)) GOTO 9
      GG(3,4)=2.D0*X(4) 
9     RETURN    
      END       
C
      SUBROUTINE TP41(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=4       
      NILI=0    
      NINL=0    
      NELI=1    
      NENL=0    
      DO 6 I=1,4
      X(I)=1.D0 
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
6     XL(I)=0.D0
      XU(1)=1.D0
      XU(2)=1.D0
      XU(3)=1.D0
      XU(4)=2.D0
      GF(4)=0.D0
      GG(1,1)=1.D0      
      GG(1,2)=2.D0      
      GG(1,3)=2.D0      
      GG(1,4)=-1.D0     
      LEX=.TRUE.
      NEX=1     
      XEX(1)=2.D0/3.D0  
      XEX(2)=1.D0/3.D0  
      XEX(3)=XEX(2)     
      XEX(4)=2.D0       
      FEX=52.D0/27.D0   
      RETURN    
2     FX=2.D0-X(1)*X(2)*X(3)    
      RETURN    
3     GF(1)=-X(2)*X(3)  
      GF(2)=-X(1)*X(3)  
      GF(3)=-X(1)*X(2)  
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)+2.D0*X(2)+2.D0*X(3)-X(4) 
5     RETURN    
      END       
C
      SUBROUTINE TP42(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=4       
      NILI=0    
      NINL=0    
      NELI=1    
      NENL=1    
      DO 6 I=1,4
      X(I)=1.D0 
      LXU(I)=.FALSE.    
6     LXL(I)=.FALSE.    
      GG(1,1)=1.D0      
      GG(1,2)=0.D0      
      GG(1,3)=0.D0      
      GG(1,4)=0.D0      
      GG(2,1)=0.D0      
      GG(2,2)=0.D0      
      LEX=.TRUE.
      NEX=1     
      XEX(1)=2.D0       
      XEX(2)=2.D0       
      XEX(3)=DSQRT(0.72D0)      
      XEX(4)=DSQRT(1.28D0)      
      FEX=28.D0-10.D0*DSQRT(2.D0)       
      RETURN    
2     FX=(X(1)-1.D0)**2+(X(2)-2.D0)**2+(X(3)-3.D0)**2+(X(4)-4.D0)**2    
      RETURN    
3     DO 21 I=1,4       
21    GF(I)=2.D0*(X(I)-DBLE(I))       
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)-2.D0     
      IF (INDEX1(2)) G(2)=X(3)**2+X(4)**2-2.D0  
      RETURN    
5     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,3)=2.D0*X(3) 
      GG(2,4)=2.D0*X(4) 
8     RETURN    
      END       
C
      SUBROUTINE TP43(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=4       
      NILI=0    
      NINL=3    
      NELI=0    
      NENL=0    
      DO 6 I=1,4
      X(I)=0.D0 
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      LEX=.TRUE.
      NEX=1     
      XEX(1)=0.D0       
      XEX(2)=1.D0       
      XEX(3)=2.D0       
      XEX(4)=-1.D0      
      FEX=-44.D0
      GG(3,4)=1.D0      
      RETURN    
2     FX=X(1)**2+X(2)**2+2.D0*X(3)**2+X(4)**2-5.D0*X(1)-5.D0*X(2)-21.D0 
     /               *X(3)+7.D0*X(4)   
      RETURN    
3     GF(1)=2.D0*X(1)-5.D0      
      GF(2)=2.D0*X(2)-5.D0      
      GF(3)=4.D0*X(3)-21.D0     
      GF(4)=2.D0*X(4)+7.D0      
      RETURN    
4     IF (INDEX1(1)) G(1)=-X(1)**2-X(2)**2-X(3)**2-X(4)**2-X(1)+X(2)
     /          -X(3)+X(4)+8.D0    
      IF (INDEX1(2)) G(2)=-X(1)**2-2.D0*X(2)**2-X(3)**2-2.D0*X(4)**2    
     /          +X(1)+X(4)+10.D0  
      IF (INDEX1(3)) G(3)=-2.D0*X(1)**2-X(2)**2-X(3)**2-2.D0*X(1)+X(2)  
     /          +X(4)+5.D0
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=-2.D0*X(1)-1.D0   
      GG(1,2)=-2.D0*X(2)+1.D0   
      GG(1,3)=-2.D0*X(3)-1.D0   
      GG(1,4)=-2.D0*X(4)+1.D0   
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,1)=-2.D0*X(1)+1.D0   
      GG(2,2)=-4.D0*X(2)
      GG(2,3)=-2.D0*X(3)
      GG(2,4)=-4.D0*X(4)+1.D0   
8     IF (.NOT.INDEX2(3)) GOTO 9
      GG(3,1)=-4.D0*X(1)-2.D0   
      GG(3,2)=-2.D0*X(2)+1.D0   
      GG(3,3)=-2.D0*X(3)
9     RETURN    
      END       
C
      SUBROUTINE TP44(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      GOTO (1,2,3,4,5),MODE     
1     N=4       
      NILI=6    
      NINL=0    
      NELI=0    
      NENL=0    
      DO 6 I=1,4
      X(I)=0.D0 
      XL(I)=0.D0
      LXL(I)=.TRUE.     
6     LXU(I)=.FALSE.    
      DO 15 I=1,6       
      DO 15 J=1,4       
15    GG(I,J)=0.D0      
      GG(1,1)=-1.D0     
      GG(1,2)=-2.D0     
      GG(2,1)=-4.D0     
      GG(2,2)=-1.D0     
      GG(3,1)=-3.D0     
      GG(3,2)=-4.D0     
      GG(4,3)=-2.D0     
      GG(4,4)=-1.D0     
      GG(5,3)=-1.D0     
      GG(5,4)=-2.D0     
      GG(6,3)=-1.D0     
      GG(6,4)=-1.D0     
      LEX=.TRUE.
      NEX=1     
      XEX(1)=0.D0       
      XEX(2)=3.D0       
      XEX(3)=0.D0       
      XEX(4)=4.D0       
      FEX=-15.D0
      RETURN    
2     FX=X(1)-X(2)-X(3)-X(1)*X(3)+X(1)*X(4)+X(2)*X(3)-X(2)*X(4) 
      RETURN    
3     GF(1)=1.D0-X(3)+X(4)      
      GF(2)=-1.D0+X(3)-X(4)     
      GF(3)=-1.D0-X(1)+X(2)     
      GF(4)=X(1)-X(2)   
      RETURN    
4     IF (INDEX1(1)) G(1)=8.D0-X(1)-2.D0*X(2)   
      IF (INDEX1(2)) G(2)=12.D0-4.D0*X(1)-X(2)  
      IF (INDEX1(3)) G(3)=12.D0-3.D0*X(1)-4.D0*X(2)     
      IF (INDEX1(4)) G(4)=8.D0-2.D0*X(3)-X(4)   
      IF (INDEX1(5)) G(5)=8.D0-X(3)-2.D0*X(4)   
      IF (INDEX1(6)) G(6)=5.D0-X(3)-X(4)
5     RETURN    
      END       
C
      SUBROUTINE TP45(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,4),MODE     
1     N=5       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=0    
      DO 6 I=1,5
      X(I)=2.D0 
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      XL(I)=0.D0
6     XU(I)=DBLE(I)   
      LEX=.TRUE.
      NEX=1     
      DO 30 I=1,5       
30    XEX(I)=DBLE(I)  
      FEX=1.D0  
      RETURN    
2     FX=2.D0-X(1)*X(2)*X(3)*X(4)*X(5)/120.D0   
      RETURN    
3     GF(1)=-X(2)*X(3)*X(4)*X(5)/120.D0 
      GF(2)=-X(1)*X(3)*X(4)*X(5)/120.D0 
      GF(3)=-X(1)*X(2)*X(4)*X(5)/120.D0 
      GF(4)=-X(1)*X(2)*X(3)*X(5)/120.D0 
      GF(5)=-X(1)*X(2)*X(3)*X(4)/120.D0 
4     RETURN    
      END       
C
      SUBROUTINE TP46(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=5       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=2    
      X(1)=0.5D0*DSQRT(2.D0)    
      X(2)=1.75D0       
      X(3)=0.5D0
      X(4)=2.D0 
      X(5)=2.D0 
      DO 6 I=1,5
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      GG(1,2)=0.D0      
      GG(1,3)=0.D0      
      GG(2,1)=0.D0      
      GG(2,2)=1.D0      
      GG(2,5)=0.D0      
      LEX=.TRUE.
      NEX=1     
      DO 30 I=1,5       
30    XEX(I)=1.D0       
      FEX=0.D0  
      RETURN    
2     FX=(X(1)-X(2))**2+(X(3)-1.D0)**2+(X(4)-1.D0)**4+(X(5)-1.D0)**6    
      RETURN    
3     GF(1)=2.D0*(X(1)-X(2))    
      GF(2)=-GF(1)      
      GF(3)=2.D0*(X(3)-1.D0)    
      GF(4)=4.D0*(X(4)-1.D0)**3 
      GF(5)=6.D0*(X(5)-1.D0)**5 
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)**2*X(4)+DSIN(X(4)-X(5))-1.D0     
      IF (INDEX1(2)) G(2)=X(2)+X(3)**4*X(4)**2-2.D0     
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=2.D0*X(1)*X(4)    
      GG(1,5)=-DCOS(X(4)-X(5))  
      GG(1,4)=X(1)**2-GG(1,5)   
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,3)=4.D0*X(3)**3*X(4)**2      
      GG(2,4)=2.D0*X(3)**4*X(4) 
8     RETURN    
      END       
C
      SUBROUTINE TP47(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX)
      DOUBLEPRECISION V1,V2,V3,V4  
      INTEGER I,J
      GOTO (1,2,3,4,5),MODE     
1     N=5       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=3    
      X(1)=2.D0 
      X(2)=DSQRT(2.D0)  
      X(3)=-1.D0
      X(4)=2.D0-X(2)    
      X(5)=0.5D0
      DO 6 I=1,5
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      DO 15 I=1,3       
      DO 15 J=1,5       
15    GG(I,J)=0.D0      
      GG(1,1)=1.D0      
      GG(2,2)=1.D0      
      GG(2,4)=1.D0      
      LEX=.TRUE.
      NEX=1     
      DO 30 I=1,5       
30    XEX(I)=1.D0       
      FEX=0.D0  
      RETURN    
2     FX=(X(1)-X(2))**2+(X(2)-X(3))**2+(X(3)-X(4))**4+(X(4)-X(5))**4    
      RETURN    
3     V1=2.D0*(X(1)-X(2))       
      V2=2.D0*(X(2)-X(3))       
      V3=4.D0*(X(3)-X(4))**3    
      V4=4.D0*(X(4)-X(5))**3    
      GF(1)=V1  
      GF(2)=-V1+V2      
      GF(3)=-V2+V3      
      GF(4)=-V3+V4      
      GF(5)=-V4 
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)+X(2)**2+X(3)**3-3.D0     
      IF (INDEX1(2)) G(2)=X(2)-X(3)**2+X(4)-1.D0
      IF (INDEX1(3)) G(3)=X(1)*X(5)-1.D0
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,2)=2.D0*X(2) 
      GG(1,3)=3.D0*X(3)**2      
7     IF (INDEX2(2)) GG(2,3)=-2.D0*X(3) 
8     IF (.NOT.INDEX2(3)) GOTO 9
      GG(3,1)=X(5)      
      GG(3,5)=X(1)      
9     RETURN    
      END       
C
      SUBROUTINE TP48(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=5       
      NILI=0    
      NINL=0    
      NELI=2    
      NENL=0    
      X(1)=3.D0 
      X(2)=5.D0 
      X(3)=-3.D0
      X(4)=2.D0 
      X(5)=-2.D0
      DO 6 I=1,5
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      GG(2,1)=0.D0      
      GG(2,2)=0.D0      
      DO 20 I=1,5       
20    GG(1,I)=1.D0      
      GG(2,3)=1.D0      
      GG(2,4)=-2.D0     
      GG(2,5)=-2.D0     
      LEX=.TRUE.
      NEX=1     
      DO 30 I=1,5       
30    XEX(I)=1.D0       
      FEX=0.D0  
      RETURN    
2     FX=(X(1)-1.D0)**2+(X(2)-X(3))**2+(X(4)-X(5))**2   
      RETURN    
3     GF(1)=2.D0*(X(1)-1.D0)    
      GF(2)=2.D0*(X(2)-X(3))    
      GF(3)=-GF(2)      
      GF(4)=2.D0*(X(4)-X(5))    
      GF(5)=-GF(4)      
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)+X(2)+X(3)+X(4)+X(5)-5.D0 
      IF (INDEX1(2)) G(2)=X(3)-2.D0*(X(4)+X(5))+3.D0    
5     RETURN    
      END       
C
      SUBROUTINE TP49(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=5       
      NILI=0    
      NINL=0    
      NELI=2    
      NENL=0    
      X(1)=10.D0
      X(2)=7.D0 
      X(3)=2.D0 
      X(4)=-3.D0
      X(5)=0.8D0
      DO 6 I=1,5
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      GG(1,5)=0.D0      
      GG(2,1)=0.D0      
      GG(2,2)=0.D0      
      GG(2,4)=0.D0      
      GG(1,1)=1.D0      
      GG(1,2)=1.D0      
      GG(1,3)=1.D0      
      GG(1,4)=4.D0      
      GG(2,3)=1.D0      
      GG(2,5)=5.D0      
      LEX=.TRUE.
      NEX=1     
      DO 30 I=1,5       
30    XEX(I)=1.D0       
      FEX=0.D0  
      RETURN    
2     FX=((X(1)-X(2))**2+(X(3)-1.D0)**2+(X(4)-1.D0)**4
     /                                  +(X(5)-1.D0)**6)
      RETURN    
3     GF(1)=2.D0*(X(1)-X(2))
      GF(2)=-GF(1)      
      GF(3)=2.D0*(X(3)-1.D0)
      GF(4)=4.D0*(X(4)-1.D0)**3
      GF(5)=6.D0*(X(5)-1.D0)**5
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)+X(2)+X(3)+4.D0*X(4)-7.D0 
      IF (INDEX1(2)) G(2)=X(3)+5.D0*X(5)-6.D0   
5     RETURN    
      END       
C
      SUBROUTINE TP50(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION  V1,V2,V3,V4
      GOTO (1,2,3,4,5),MODE     
1     N=5       
      NILI=0    
      NINL=0    
      NELI=3    
      NENL=0    
      X(1)=35.D0
      X(2)=-31.D0       
      X(3)=11.D0
      X(4)=5.D0 
      X(5)=-5.D0
      DO 6 I=1,5
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      DO 15 I=1,3       
      DO 15 J=1,5       
15    GG(I,J)=0.D0      
      GG(1,1)=1.D0      
      GG(1,2)=2.D0      
      GG(1,3)=3.D0      
      GG(2,2)=1.D0      
      GG(2,3)=2.D0      
      GG(2,4)=3.D0      
      GG(3,3)=1.D0      
      GG(3,4)=2.D0      
      GG(3,5)=3.D0      
      LEX=.TRUE.
      NEX=1     
      DO 30 I=1,5       
30    XEX(I)=1.D0       
      FEX=0.D0  
      RETURN    
2     FX=(X(1)-X(2))**2+(X(2)-X(3))**2+(X(3)-X(4))**4+(X(4)-X(5))**4    
      RETURN    
3     V1=2.D0*(X(1)-X(2))       
      V2=2.D0*(X(2)-X(3))       
      V3=4.D0*(X(3)-X(4))**3    
      V4=4.D0*(X(4)-X(5))**3    
      GF(1)=V1  
      GF(2)=-V1+V2      
      GF(3)=-V2+V3      
      GF(4)=-V3+V4      
      GF(5)=-V4 
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)+2.D0*X(2)+3.D0*X(3)-6.D0 
      IF (INDEX1(2)) G(2)=X(2)+2.D0*X(3)+3.D0*X(4)-6.D0 
      IF (INDEX1(3)) G(3)=X(3)+2.D0*X(4)+3.D0*X(5)-6.D0 
5     RETURN    
      END       
C
      SUBROUTINE TP51(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      GOTO (1,2,3,4,5),MODE     
1     N=5       
      NILI=0    
      NINL=0    
      NELI=3    
      NENL=0    
      X(1)=2.5D0
      X(2)=0.5D0
      X(3)=2.D0 
      X(4)=-1.D0
      X(5)=0.5D0
      DO 6 I=1,5
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      DO 15 I=1,3       
      DO 15 J=1,5       
15    GG(I,J)=0.D0      
      GG(1,1)=1.D0      
      GG(1,2)=3.D0      
      GG(2,3)=1.D0      
      GG(2,4)=1.D0      
      GG(2,5)=-2.D0     
      GG(3,2)=1.D0      
      GG(3,5)=-1.D0     
      LEX=.TRUE.
      NEX=1     
      DO 30 I=1,5       
30    XEX(I)=1.D0       
      FEX=0.D0  
      RETURN    
2     FX=(X(1)-X(2))**2+(X(2)+X(3)-2.D0)**2+(X(4)-1.D0)**2      
     1     +(X(5)-1.D0)**2      
      RETURN    
3     GF(1)=2.D0*(X(1)-X(2))    
      GF(3)=2.D0*(X(2)+X(3)-2.D0)       
      GF(2)=GF(3)-GF(1) 
      GF(4)=2.D0*(X(4)-1.D0)    
      GF(5)=2.D0*(X(5)-1.D0)    
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)+3.D0*X(2)-4.D0   
      IF (INDEX1(2)) G(2)=X(3)+X(4)-2.D0*X(5)   
      IF (INDEX1(3)) G(3)=X(2)-X(5)     
5     RETURN    
      END       
C
      SUBROUTINE TP52(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      GOTO (1,2,3,4,5),MODE     
1     N=5       
      NILI=0    
      NINL=0    
      NELI=3    
      NENL=0    
      DO 6 I=1,5
      X(I)=2.D0 
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      DO 15 I=1,3       
      DO 15 J=1,5       
15    GG(I,J)=0.D0      
      GG(1,1)=1.D0      
      GG(1,2)=3.D0      
      GG(2,3)=1.D0      
      GG(2,4)=1.D0      
      GG(2,5)=-2.D0     
      GG(3,2)=1.D0      
      GG(3,5)=-1.D0     
      LEX=.TRUE.
      NEX=1     
      XEX(1)=-33.D0/349.D0      
      XEX(2)=11.D0/349.D0       
      XEX(3)=180.D0/349.D0      
      XEX(4)=-158.D0/349.D0     
      XEX(5)=XEX(2)     
      FEX=1859.D0/349.D0
      RETURN    
2     FX=(4.D0*X(1)-X(2))**2+(X(2)+X(3)-2.D0)**2+(X(4)-1.D0)**2+(X(5)   
     /      -1.D0)**2    
      RETURN    
3     GF(1)=8.D0*(X(1)*4.D0-X(2))       
      GF(3)=2.D0*(X(2)+X(3)-2.D0)       
      GF(2)=-0.25D0*GF(1)+GF(3) 
      GF(4)=2.D0*(X(4)-1.D0)    
      GF(5)=2.D0*(X(5)-1.D0)    
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)+3.D0*X(2)
      IF (INDEX1(2)) G(2)=X(3)+X(4)-2.D0*X(5)   
      IF (INDEX1(3)) G(3)=X(2)-X(5)     
5     RETURN    
      END       
C
      SUBROUTINE TP53(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      GOTO (1,2,3,4,5),MODE     
1     N=5       
      NILI=0    
      NINL=0    
      NELI=3    
      NENL=0    
      DO 6 I=1,5
      X(I)=2.D0 
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      XL(I)=-10.D0      
6     XU(I)=10.D0       
      DO 15 I=1,3       
      DO 15 J=1,5       
15    GG(I,J)=0.D0      
      GG(1,1)=1.D0      
      GG(1,2)=3.D0      
      GG(2,3)=1.D0      
      GG(2,4)=1.D0      
      GG(2,5)=-2.D0     
      GG(3,2)=1.D0      
      GG(3,5)=-1.D0     
      LEX=.TRUE.
      NEX=1     
      XEX(1)=-33.D0/43.D0       
      XEX(2)=11.D0/43.D0
      XEX(3)=27.D0/43.D0
      XEX(4)=-5.D0/43.D0
      XEX(5)=11.D0/43.D0
      FEX=176.D0/43.D0  
      RETURN    
2     FX=(X(1)-X(2))**2+(X(2)+X(3)-2.D0)**2+(X(4)-1.D0)**2+(X(5)
     /     -1.D0)**2    
      RETURN    
3     GF(1)=2.D0*(X(1)-X(2))    
      GF(3)=2.D0*(X(2)+X(3)-2.D0)       
      GF(2)=GF(3)-GF(1) 
      GF(4)=2.D0*(X(4)-1.D0)    
      GF(5)=2.D0*(X(5)-1.D0)    
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)+3.D0*X(2)
      IF (INDEX1(2)) G(2)=X(3)+X(4)-2.D0*X(5)   
      IF (INDEX1(3)) G(3)=X(2)-X(5)     
5     RETURN    
      END       
C
      SUBROUTINE TP54(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      COMMON     /D54/DQ,V1,V2,V3,V4,V5,V6,V7,V8,V9,Q  
      DOUBLEPRECISION DQ,V1,V2,V3,V4,V5,V6,V7,V8,V9,Q       
      DIMENSION DQ(6)   
      GOTO (1,2,3,4,5),MODE     
1     N=6       
      NILI=0    
      NINL=0    
      NELI=1    
      NENL=0    
      X(1)=6.D+3
      X(2)=1.5D0
      X(3)=4.D+6
      X(4)=2.D0 
      X(5)=3.D-3
      X(6)=5.D+7
      DO 6 I=1,6
      LXL(I)=.TRUE.     
6     LXU(I)=.TRUE.     
      XL(1)=0.D0
      XL(2)=-10.D0      
      XL(3)=0.D0
      XL(4)=0.D0
      XL(5)=0.D0       
      XL(6)=0.D0
      XU(1)=2.D+4       
      XU(2)=10.D0       
      XU(3)=1.D+7       
      XU(4)=20.D0       
      XU(5)=1.D0
      XU(6)=2.D+8       
      GG(1,1)=1.D0    
      GG(1,2)=4.D+3
      DO 30 I=3,6       
30    GG(1,I)=0.0      
      LEX=.TRUE.
      NEX=1     
      XEX(1)=9.16D+4/7.D0       
      XEX(2)=79.D0/70.D0
      XEX(3)=2.D+6      
      XEX(4)=10.D0      
      XEX(5)=1.D-3      
      XEX(6)=1.D+8      
      FEX=-DEXP(-27.D0/280.D0)  
c      FEX=-0.90807476D0 
      RETURN    
2     V1=X(1)-1.D+4     
      V2=X(2)-1.D0      
      V3=X(3)-2.D+6     
      V4=X(4)-10.D0     
      V5=X(5)-1.D-3     
      V6=X(6)-1.D+8     
      V7=1.D0/0.96D0    
      V8=1.D0/4.9D+13   
      V9=1.D0/2.45D+13  
      Q = (1.5625D-8*V1**2 + 5.D-5*V1*V2 + V2**2)*V7 + V3**2*V8
     /   + 4.D-4*V4**2 + 4.D+2*V5**2 + 4.D-18*V6**2       
      FX = -DEXP(-0.5D0*Q)
      RETURN    
3     V1=X(1)-1.D+4     
      V2=X(2)-1.D0      
      V3=X(3)-2.D+6     
      V4=X(4)-10.D0     
      V5=X(5)-1.D-3     
      V6=X(6)-1.D+8     
      V7=1.D0/0.96D0    
      V8=1.D0/4.9D+13   
      V9=1.D0/2.45D+13  
      DQ(1)=(3.125D-8*V1+5.D-5*V2)*V7   
      DQ(2)=(5.D-5*V1+2.D0*V2)*V7       
      DQ(3)=V3*V9       
      DQ(4)=8.D-4*V4    
      DQ(5)=800.D0*V5   
      DQ(6)=8.D-18*V6   
      DO 31 I=1,6       
31    GF(I)=0.5D0*DEXP(-0.5D0*Q)*DQ(I)  
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)+4.D+3*X(2)-1.76D+4
5     RETURN    
      END       
C
      SUBROUTINE TP55(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION V1,X14    
      GOTO (1,2,3,4,5),MODE     
1     N=6       
      NILI=0    
      NINL=0    
      NELI=6    
      NENL=0    
      X(1)=1.D0 
      X(2)=2.D0 
      X(3)=0.D0 
      X(4)=0.D0 
      X(5)=0.D0 
      X(6)=2.D0 
      DO 6 I=1,6
      LXL(I)=.TRUE.     
      LXU(I)=.FALSE.    
6     XL(I)=0.D0
      LXU(1)=.TRUE.     
      LXU(4)=.TRUE.     
      XU(1)=1.D+0
      XU(4)=1.D+0
      DO 15 I=1,6       
      DO 15 J=1,6       
15    GG(I,J)=0.D0      
      GF(2)=2.D0
      GF(3)=0.D0
      GF(5)=4.D0
      GF(6)=0.D0
      GG(1,1)=1.D0      
      GG(1,2)=2.D0      
      GG(1,5)=5.D0      
      GG(2,1)=1.D0      
      GG(2,2)=1.D0      
      GG(2,3)=1.D0      
      GG(3,4)=1.D0      
      GG(3,5)=1.D0      
      GG(3,6)=1.D0      
      GG(4,1)=1.D0      
      GG(4,4)=1.D0      
      GG(5,2)=1.D0      
      GG(5,5)=1.D0      
      GG(6,3)=1.D0      
      GG(6,6)=1.D0      
      LEX=.TRUE.
      NEX=1     
      XEX(1)=0.D0       
      XEX(2)=4.D0/3.D0  
      XEX(3)=5.D0/3.D0  
      XEX(4)=1.D0       
      XEX(5)=2.D0/3.D0  
      XEX(6)=1.D0/3.D0  
      FEX=19.D0/3.D0    
      RETURN  
2     X14 = X(1)*X(4)
      IF (X14.GT.1.0D1) X14 = 1.0D1   
      FX=X(1)+2.D0*X(2)+4.D0*X(5)+DEXP(X14)       
      RETURN    
3     V1=DEXP(X(1)*X(4))
      GF(1)=1.D0+X(4)*V1
      GF(4)=X(1)*V1     
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)+2.D0*X(2)+5.D0*X(5)-6.D0 
      IF (INDEX1(2)) G(2)=X(1)+X(2)+X(3)-3.D0   
      IF (INDEX1(3)) G(3)=X(4)+X(5)+X(6)-2.D0   
      IF (INDEX1(4)) G(4)=X(1)+X(4)-1.D0
      IF (INDEX1(5)) G(5)=X(2)+X(5)-2.D0
      IF (INDEX1(6)) G(6)=X(3)+X(6)-2.D0
5     RETURN    
      END       
C
      SUBROUTINE TP56(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      GOTO (1,2,3,4,5),MODE     
1     N=7       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=4    
      X(1)=1.D0 
      X(2)=1.D0 
      X(3)=1.D0 
      DO 6 I=1,7
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      X(4)=DASIN(DSQRT(1.D0/4.2D0))    
      X(5)=X(4) 
      X(6)=X(4) 
      X(7)=DASIN(DSQRT(5.D0/7.2D0))    
      DO 30 I=4,7       
30    GF(I)=0.D0
      DO 15 I=1,4       
      DO 15 J=1,7       
15    GG(I,J)=0.D0      
      GG(1,1)=1.D0      
      GG(2,2)=1.D0      
      GG(3,3)=1.D0      
      GG(4,1)=1.D0      
      GG(4,2)=2.D0      
      GG(4,3)=2.D0      
      LEX=.TRUE.
      NEX=1    
      XEX(1)=2.4D0      
      XEX(2)=1.2D0      
      XEX(3)=1.2D0      
      XEX(4)=DASIN(DSQRT(4.D0/7.D0))   
      XEX(5)=DASIN(DSQRT(2.D0/7.D0))   
      XEX(6)=XEX(5)     
      XEX(7)=2.D0*DATAN(1.D0)   
      FEX=-3.456D0      
      RETURN    
2     FX=-X(1)*X(2)*X(3)
      RETURN    
3     GF(1)=-X(2)*X(3)  
      GF(2)=-X(1)*X(3)  
      GF(3)=-X(1)*X(2)  
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)-4.2D0*DSIN(X(4))**2      
      IF (INDEX1(2)) G(2)=X(2)-4.2D0*DSIN(X(5))**2      
      IF (INDEX1(3)) G(3)=X(3)-4.2D0*DSIN(X(6))**2      
      IF (INDEX1(4)) G(4)=X(1)+2.D0*X(2)+2.D0*X(3)-7.2D0*DSIN(X(7))**2  
      RETURN    
5     IF (INDEX2(1)) GG(1,4)=-8.4D0*DSIN(X(4))*DCOS(X(4))       
      IF (INDEX2(2)) GG(2,5)=-8.4D0*DSIN(X(5))*DCOS(X(5))       
      IF (INDEX2(3)) GG(3,6)=-8.4D0*DSIN(X(6))*DCOS(X(6))       
      IF (INDEX2(4)) GG(4,7)=-14.4D0*DSIN(X(7))*DCOS(X(7))      
      RETURN    
      END       
C
      SUBROUTINE TP57(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)       
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      COMMON     /D57/ S,A,B
      DIMENSION S(2),A(44),B(44) 
      DOUBLEPRECISION F,DF,S,A,B,T,V1       
      IF (MODE - 2) 1,18,18     
1     N=2       
      NILI=0    
      NINL=1    
      NELI=0    
      NENL=0    
      X(1)=0.42D0       
      X(2)=5.D0 
      DO 6 I=1,2
      LXL(I)=.TRUE.     
6     LXU(I)=.FALSE.    
      XL(1)=0.4D0       
      XL(2)=-4.D0       
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.419952674511D+00     
      XEX(2)=0.128484562930D+01     
      FEX=0.284596697213D-01
      RETURN    
   18 CONTINUE
      DO 20 I=1,2       
      A(I)=8.D0 
      A(16+I)=18.D0     
      A(30+I)=28.D0     
      A(35+I)=32.D0     
      A(38+I)=36.D0     
      A(40+I)=38.D0     
      B(I)=0.49D0       
      B(6+I)=0.46D0     
      B(11+I)=0.43D0    
      B(14+I)=0.43D0    
      B(18+I)=0.42D0    
      B(21+I)=0.41D0    
      B(25+I)=0.4D0     
      B(29+I)=0.41D0    
      B(36+I)=0.4D0     
      B(40+I)=0.4D0     
  20  B(42+I)=0.39D0    
      DO 21 I=1,3       
      A(10+I)=14.D0     
      A(13+I)=16.D0     
      A(18+I)=20.D0     
      A(21+I)=22.D0     
      A(24+I)=24.D0     
      A(27+I)=26.D0     
      A(32+I)=30.D0     
21    B(31+I)=0.4D0     
      DO  22 I=1,4      
      A(2+I)=10.D0      
22    A(6+I)=12.D0      
      A(38)=34.D0       
      A(43)=40.D0       
      A(44)=42.D0       
      B(3)=0.48D0       
      B(4)=0.47D0       
      B(5)=0.48D0       
      B(6)=0.47D0       
      B(9)=0.45D0       
      B(10)=0.43D0      
      B(11)=0.45D0      
      B(14)=0.44D0      
      B(17)=0.46D0      
      B(18)=0.45D0      
      B(21)=0.43D0      
      B(24)=0.4D0       
      B(25)=0.42D0      
      B(28)=0.41D0      
      B(29)=0.4D0       
      B(35)=0.38D0      
      B(36)=0.41D0      
      B(39)=0.41D0      
      B(40)=0.38D0      
      IF (MODE - 4) 17,4,5      
17    DO 30 I=1,44      
30    F(I)=B(I)-X(1)-(0.49D0-X(1))*DEXP(-X(2)*(A(I)-8.D0))      
      GOTO (1,2,3),MODE 
2     T=0.D0    
      DO 19 I=1,44      
19    T=T+F(I)**2       
      FX=T
      RETURN    
3     S(1)=0.D0 
      S(2)=0.D0 
      DO 31 I=1,44      
      V1=DEXP(-X(2)*(A(I)-8.D0))
      DF(I,1)=-1.D0+V1  
      DF(I,2)=(A(I)-8.D0)*(0.49D0-X(1))*V1      
      DO 32 J=1,2       
32    S(J)=S(J)+2.D0*F(I)*DF(I,J)       
31    CONTINUE  
      GF(1)=S(1)
      GF(2)=S(2)
      RETURN    
4     IF (INDEX1(1)) G(1)=-X(1)*X(2)+0.49D0*X(2)-0.09D0 
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=-X(2)     
      GG(1,2)=-X(1)+0.49D0      
7     RETURN    
      END       
C
      SUBROUTINE TP58(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=0    
      NINL=3    
      NELI=0    
      NENL=0    
      X(1)=-2.D0
      X(2)=1.D0 
      LXL(1)=.TRUE.     
      LXU(1)=.TRUE.     
      LXL(2)=.FALSE.    
      LXU(2)=.FALSE.    
C      XL(1)=-0.5D0      
      XL(1)=-2.0D0      
      XU(1)=0.5D0       
      GG(1,1)=-1.D0     
      GG(2,2)=-1.D0     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=-0.786150483331
      XEX(2)=0.618034533851
      FEX=3.19033354957
      RETURN    
2     FX=100.0D0*(X(2)-X(1)**2)**2+(1.0D0-X(1))**2
      RETURN    
3     GF(2)=200.0D0*(X(2)-X(1)**2)       
      GF(1)=-2.0D0*(X(1)*(GF(2)-1.0D0)+1.0D0)      
      RETURN    
4     IF (INDEX1(1)) G(1)=X(2)**2-X(1)  
      IF (INDEX1(2)) G(2)=X(1)**2-X(2)  
      IF (INDEX1(3)) G(3)=X(1)**2+X(2)**2-1.0D0  
      RETURN    
5     IF (INDEX2(1)) GG(1,2)=2.0D0*X(2)  
      IF (INDEX2(2)) GG(2,1)=2.0D0*X(1)  
      IF (.NOT.INDEX2(3)) GOTO 9
      GG(3,1)=2.0D0*X(1) 
      GG(3,2)=2.0D0*X(2) 
9     RETURN    
      END       
C
      SUBROUTINE TP59(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D59/ X11,X12,X13,X14,X21,X22,X23,X24,XX11,XX12,XX13,  
     /     XX21,XX31   
      DOUBLEPRECISION X11,X12,X13,X14,X21,X22,X23,X24,XX11,XX12,XX13, 
     /     XX21,XX31    
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=2       
      NILI=0    
      NINL=3    
      NELI=0    
      NENL=0    
      X(1)=90.D0
      X(2)=10.D0
      DO 6 I=1,2
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
6     XL(I)=0.D0
      XU(1)=75.D0       
      XU(2)=65.D0       
      GG(2,2)=1.D0      
      GG(3,1)=-5.D0     
      NEX=2
      LEX=.FALSE.       
      XEX(1)=13.5501042366D0  
      XEX(2)=51.6601812877D0  
      FEX=-7.80422632408D0    
c      XEX(3)=0.463995762710D+2
c      XEX(4)=0.522196899513D+2
C     FEX=-0.675456604292D+1  
      RETURN    
2     X11=X(1)  
      X12=X11*X11       
      X13=X12*X11       
      X14=X13*X11       
      X21=X(2)  
      X22=X21*X21       
      X23=X22*X21       
      X24=X23*X21       
      FX=-75.196D0+3.8112D0*X11-0.12694D0*X12+2.0567D-3*X13
     /    -1.0345D-5*X14      
     /    +6.8306D0*X21-3.0234D-2*X11*X21+1.28134D-3*X12*X21
     /    -3.5256D-5*X13*X21      
     /    +2.266D-7*X14*X21-0.25645D0*X22+3.4604D-3*X23
     /    -1.3514D-5*X24+28.106D0 
     /         /(X21+1.D0)+5.2375D-6*X12*X22
     /    +6.3D-8*X13*X22-7.D-10*X13*X23       
     /    -3.4054D-4*X11*X22+1.6638D-6*X11*X23
     /    +2.8673D0*DEXP(5.D-4*X11*X21) 
      RETURN    
3     X11=X(1)  
      X12=X11*X11       
      X13=X12*X11       
      X14=X13*X11       
      X21=X(2)  
      X22=X21*X21       
      X23=X22*X21       
      XX11=X11*X21      
      XX12=X11*X22      
      XX21=X12*X21      
      XX31=X13*X21      
      GF(1)=3.8112D0-0.25388D0*X11+6.1701D-3*X12-4.138D-5*X13-3.0234D-2 
     -*X21      
     -+2.56268D-3*XX11-1.05768D-4*XX21+9.064D-7*XX31+1.0475D-5*XX12     
     -+1.89D-7*X12*X22-2.1D-9*X12*X23-3.4054D-4*X22+1.6638D-6*X23       
     -+1.43365D-3*X21*DEXP(5.D-4*XX11)  
      GF(2)=6.8306D0-3.0234D-2*X11+1.28134D-3*X12-3.5256D-5*X13+2.266D-7
     -*X14-0.5129D0*X21+1.03812D-2*X22-5.4056D-5*X23-28.106D0/(X21      
     -+1.D0)**2 
     -+1.0475D-5*XX21+1.26D-7*XX31-2.1D-9*X13*X22-6.8108D-4*XX11
     -+4.9914D-6*XX12+1.43365D-3*X11*DEXP(5.D-4*XX11)   
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)*X(2)-700.D0      
      IF (INDEX1(2)) G(2)=X(2)-8.D-3*X(1)**2    
      IF (INDEX1(3)) G(3)=(X(2)-50.D0)**2-5.D0*(X(1)-55.D0)     
      RETURN    
5     IF(.NOT.INDEX2(1)) GOTO 7 
      GG(1,1)=X(2)      
      GG(1,2)=X(1)      
7     IF (INDEX2(2)) GG(2,1)=-1.6D-2*X(1)       
      IF (INDEX2(3)) GG(3,2)=2.D0*(X(2)-50.D0)  
      RETURN    
      END       
C
      SUBROUTINE TP60(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION V1 
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=1    
      DO 6 I=1,3
      X(I)=2.D0 
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      XL(I)=-10.D0      
6     XU(I)=10.D0       
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.110485902423D+01     
      XEX(2)=0.119667419413D+01     
      XEX(3)=0.153526225739D+01     
      FEX= 0.325682002513D-01 
      RETURN    
2     FX=(X(1)-1.D0)**2+(X(1)-X(2))**2+(X(2)-X(3))**4   
      RETURN    
3     V1=2.D0*(X(1)-X(2))       
      GF(1)=2.D0*(X(1)-1.D0)+V1 
      GF(3)=-4.D0*(X(2)-X(3))**3
      GF(2)=-GF(3)-V1   
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)*(1.D0+X(2)**2)+X(3)**4-4.D0-3.D0 
     1     *DSQRT(2.D0) 
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=1.D0+X(2)**2      
      GG(1,2)=2.D0*X(1)*X(2)    
      GG(1,3)=4.D0*X(3)**3      
7     RETURN    
      END       
C
      SUBROUTINE TP61(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=2    
      DO 6 I=1,3
      X(I)=0.0D0 
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.532677015744D+01     
      XEX(2)=-0.211899863998D+01     
      XEX(3)=0.321046423906D+01     
      FEX=-0.143646142201D+03 
      RETURN    
2     FX=4.D0*X(1)**2+2.D0*X(2)**2+2.D0*X(3)**2-33.D0*X(1)+16.D0*X(2)   
     1     -24.D0*X(3)  
      RETURN    
3     GF(1)=8.D0*X(1)-33.D0     
      GF(2)=4.D0*X(2)+16.D0     
      GF(3)=4.D0*X(3)-24.D0     
      RETURN    
4     IF (INDEX1(1)) G(1)=3.D0*X(1)-2.D0*X(2)**2-7.D0   
      IF (INDEX1(2)) G(2)=4.D0*X(1)-X(3)**2-11.D0       
      RETURN    
5     IF (INDEX2(1)) GG(1,2)=-4.D0*X(2) 
      IF (INDEX2(2)) GG(2,3)=-2.D0*X(3) 
      GG(1,1)=3.D0      
      GG(1,3)=0.D0      
      GG(2,1)=4.D0      
      GG(2,2)=0.D0      
      RETURN    
      END       
C
      SUBROUTINE TP62(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION B1,B2,B3,C1,C2,C3,V1,V2,V3,V4,V5,V6,V7,S,RB1,
     /     RB2,RB3,RC1,RC2,RC3  
      IF (MODE - 2) 1,17,17     
1     N=3       
      NILI=0    
      NINL=0    
      NELI=1    
      NENL=0    
      X(1)=0.7D0
      X(2)=0.2D0
      X(3)=0.1D0
      DO 6 I=1,3
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      XL(I)=0.D0
      XU(I)=1.D0
6     GG(1,I)=1.D0      
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.617813298210D+00     
      XEX(2)=0.328202155786D+00     
      XEX(3)=0.539845460119D-01     
      FEX=-0.262725144873D+05 
      RETURN    
17    IF (MODE - 4) 18,4,5      
18    B3=X(3)+0.03D0    
      C3=0.13D0*X(3)+0.03D0     
      B2=B3+X(2)
      C2=B3+0.07D0*X(2) 
      B1=B2+X(1)
      C1=B2+0.09D0*X(1) 
      GOTO (1,2,3),MODE 
2     V5=B1/C1  
      V6=B2/C2  
      V7=B3/C3  
      IF (V5.LE.0.D0.OR.V6.LE.0.D0.OR.V7.LE.0.D0) GOTO 7
      FX=-32.174D0*(255.D0*DLOG(V5)+280.D0*DLOG(V6)+290.D0*DLOG(V7))    
      RETURN    
7     S=0.D0    
      DO 8 I=1,3
8     S=S+(X(I)-5.D0)**2
      FX=S+1.D+3-2.67D+4
      RETURN    
3     RB1=1.D0/B1       
      RB2=1.D0/B2       
      RB3=1.D0/B3       
      RC1=1.D0/C1       
      RC2=1.D0/C2       
      RC3=1.D0/C3       
      V1=-32.174D0*255.D0       
      V2=-32.174D0*280.D0       
      V3=-32.174D0*290.D0       
      V4=V1*(RB1-RC1)   
      GF(1)=V1*(RB1-0.09D0*RC1) 
      GF(2)=V4+V2*(RB2-0.07D0*RC2)      
      GF(3)=V4+V2*(RB2-RC2)+V3*(RB3-0.13D0*RC3) 
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)+X(2)+X(3)-1.D0   
5     RETURN    
      END       
C
      SUBROUTINE TP63(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=0    
      NINL=0    
      NELI=1    
      NENL=1    
      DO 6 I=1,3
      X(I)=2.D0 
      LXL(I)=.TRUE.     
      LXU(I)=.FALSE.    
6     XL(I)=0.D0
      GG(1,1)=8.D0      
      GG(1,2)=14.D0     
      GG(1,3)=7.D0      
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.351211841492D+01     
      XEX(2)=0.216988174172D+00     
      XEX(3)=0.355217403459D+01     
      FEX=0.961715172127D+03 
      RETURN    
2     FX=1.D+3-X(1)**2-2.D0*X(2)**2-X(3)**2-X(1)*X(2)-X(1)*X(3) 
      RETURN    
3     GF(1)=-2.D0*X(1)-X(2)-X(3)
      GF(2)=-4.D0*X(2)-X(1)     
      GF(3)=-2.D0*X(3)-X(1)     
      RETURN    
4     IF (INDEX1(1)) G(1)=8.D0*X(1)+14.D0*X(2)+7.D0*X(3)-56.D0  
      IF (INDEX1(2)) G(2)=X(1)**2+X(2)**2+X(3)**2-25.D0 
      RETURN    
5     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,1)=2.D0*X(1) 
      GG(2,2)=2.D0*X(2) 
      GG(2,3)=2.D0*X(3) 
8     RETURN    
      END       
C
      SUBROUTINE TP64(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=0    
      NINL=1    
      NELI=0    
      NENL=0    
      DO 6 I=1,3
      X(I)=1.D0 
      LXL(I)=.TRUE.     
      LXU(I)=.FALSE.    
6     XL(I)=1.D-5       
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.108734717597D+03     
      XEX(2)=0.851261394257D+02     
      XEX(3)=0.204324707858D+03     
      FEX=0.629984242821D+04 
      RETURN    
2     FX=5.D0*X(1)+5.D+4/X(1)+20.D0*X(2)+7.2D+4/X(2)+10.D0*X(3) 
     1     +1.44D+5/X(3)
      RETURN    
3     GF(1)=5.D0-5.D+4/X(1)**2  
      GF(2)=20.D0-7.2D+4/X(2)**2
      GF(3)=10.D0-1.44D+5/X(3)**2       
      RETURN    
4     IF (INDEX1(1)) G(1)=1.D0-4.D0/X(1)-32.D0/X(2)-120.D0/X(3) 
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=4.D0/X(1)**2      
      GG(1,2)=32.D0/X(2)**2     
      GG(1,3)=120.D0/X(3)**2    
7     RETURN    
      END       
C
      SUBROUTINE TP65(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION V1,V2      
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=0    
      NINL=1    
      NELI=0    
      NENL=0    
      X(1)=-5.D0
      X(2)=5.D0 
      X(3)=0.D0 
      DO 6 I=1,3
      LXL(I)=.TRUE.     
6     LXU(I)=.TRUE.     
      XL(1)=-4.5D0      
      XL(2)=-4.5D0      
      XL(3)=-5.D0       
      XU(1)=4.5D0       
      XU(2)=4.5D0       
      XU(3)=5.D0
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.365046182158D+01     
      XEX(2)=0.365046168940D+01     
      XEX(3)=0.462041750754D+01     
      FEX=0.953528856757D+00 
      RETURN    
2     FX=(X(1)-X(2))**2+((X(1)+X(2)-10.D0)/3.D0)**2+(X(3)-5.D0)**2      
      RETURN    
3     V1=2.D0*(X(1)-X(2))       
      V2=2.D0*(X(1)+X(2)-10.D0)/9.D0    
      GF(1)=V1+V2       
      GF(2)=-V1+V2      
      GF(3)=2.D0*(X(3)-5.D0)    
      RETURN    
4     IF (INDEX1(1)) G(1)=48.D0-X(1)**2-X(2)**2-X(3)**2 
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=-2.D0*X(1)
      GG(1,2)=-2.D0*X(2)
      GG(1,3)=-2.D0*X(3)
7     RETURN    
      END       
C
      SUBROUTINE TP66(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=3       
      NILI=0    
      NINL=2    
      NELI=0    
      NENL=0    
      X(1)=0.D0 
      X(2)=1.05D0       
      X(3)=2.9D0
      DO 6 I=1,3
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      XL(I)=0.D0
6     XU(I)=100.D0      
      XU(3)=10.D0       
      GF(1)=-0.8D0      
      GF(2)=0.D0
      GF(3)=0.2D0       
      GG(1,2)=1.D0      
      GG(1,3)=0.D0      
      GG(2,1)=0.D0      
      GG(2,3)=1.D0      
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.184126487951D+00     
      XEX(2)=0.120216787321D+01     
      XEX(3)=0.332732232258D+01     
      FEX=0.518163274159D+00 
      RETURN    
2     FX=0.2D0*X(3)-0.8D0*X(1)  
3     RETURN    
4     IF (INDEX1(1)) G(1)=X(2)-DEXP(X(1))       
      IF (INDEX1(2)) G(2)=X(3)-DEXP(X(2))       
      RETURN    
5     IF (INDEX2(1)) GG(1,1)=-DEXP(X(1))
      IF (INDEX2(2)) GG(2,2)=-DEXP(X(2))
      RETURN    
      END       
C
      SUBROUTINE TP67(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D67/ A(3),Y(8),DY2C(3),DY4C(3),DY(8,3),
     /                 RX,V1,V2,Y2C,V3,Y4C  
      DOUBLEPRECISION A,Y,DY2C,DY4C,DY,RX,V1,V2,Y2C,V3,Y4C  
      INTEGER I,J,IREP
      IF (MODE - 2) 1,17,17     
1     N=3       
      NILI=0    
      NINL=14   
      NELI=0    
      NENL=0    
      X(1)=1.745D+3     
      X(2)=1.2D+4       
      X(3)=110.D0       
      DO 6 I=1,3
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
6     XL(I)=1.D-5       
      XU(1)=2.D+3       
      XU(2)=1.6D+4      
      XU(3)=120.D0      
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.172837128614D+04     
      XEX(2)=0.160000000000D+05     
      XEX(3)=0.981415140238D+02     
      FEX=-0.116203650728D+04 
      RETURN    
17    RX=1.D0/X(1)      
      Y(2)=1.6D0*X(1)   
      DY(2,1)=1.6D0     
      DY(2,2)=0.D0      
      DY(2,3)=0.D0  
      IREP=0    
100   Y(3)=1.22D0*Y(2)-X(1)     
      DY(3,1)=1.22D0*DY(2,1)-1.D0       
      DY(3,2)=1.22D0*DY(2,2)    
      DY(3,3)=1.22D0*DY(2,3)    
      Y(6)=(X(2)+Y(3))*RX       
      DY(6,1)=(X(1)*DY(3,1)-X(2)-Y(3))*RX**2    
      DY(6,2)=(1.D0+DY(3,2))*RX 
      DY(6,3)=DY(3,1)*RX
      V1=0.01D0*X(1)*(13.167D0-2.D0*0.6667D0*Y(6))      
      V2=(112.D0+(13.167D0-0.6667D0*Y(6))*Y(6))*0.01D0  
      Y2C=X(1)*V2       
      DY2C(1)=V2+V1*DY(6,1)     
      DY2C(2)=V1*DY(6,2)
      DY2C(3)=V1*DY(6,3)
      IF (DABS(Y2C-Y(2))-1.D-3) 102,102,101     
101   DO 103 I=1,3      
103   DY(2,I)=DY2C(I)   
      Y(2)=Y2C
      IREP=IREP+1
      IF (IREP.GT.100) GOTO 102  
      GOTO 100  
102   Y(4)=93.D0
      DO 104 I=1,3      
104   DY(4,I)=0.D0  
      IREP=0    
105   Y(5)=86.35D0+1.098D0*Y(6)-0.038D0*Y(6)**2+0.325D0*(Y(4)-89.D0)    
      Y(8)=-133.D0+3.D0*Y(5)    
      Y(7)=35.82D0-0.222D0*Y(8) 
      DO 106 I=1,3      
      DY(5,I)=1.098D0*DY(6,I)-0.076D0*Y(6)*DY(6,I)+0.325D0*DY(4,I)      
      DY(8,I)=3.D0*DY(5,I)      
106   DY(7,I)=-0.222D0*DY(8,I)  
      V3=1.D0/(Y(2)*Y(7)+1.D+3*X(3))    
      Y4C=9.8D+4*X(3)*V3
      DO 107 I=1,2      
107   DY4C(I)=-9.8D+4*X(3)*(Y(2)*DY(7,I)+Y(7)*DY(2,I))/(Y(2)*Y(7)+1.D+3*
     -X(3))**2  
      DY4C(3)=9.8D+4*(Y(2)*Y(7)-X(3)*(Y(2)*DY(7,3)+Y(7)*DY(2,3)))*V3**2 
      IF (DABS(Y4C-Y(4))-1.D-4) 109,109,108     
108   Y(4)=Y4C  
      DO 110 I=1,3      
110   DY(4,I)=DY4C(I)
      IREP=IREP+1
      IF (IREP.GT.100) GOTO 109   
      GOTO 105  
109   GOTO (1,2,3,4,5),MODE     
2     FX=-(0.063D0*Y(2)*Y(5)-5.04D0*X(1)-3.36D0*Y(3)-0.035D0*X(2)-10.D0 
     1     *X(3))       
      RETURN    
3     DO 120 I=1,3      
120   A(I)=-0.063D0*(DY(2,I)*Y(5)+DY(5,I)*Y(2))+3.36D0*DY(3,I)  
      GF(1)=A(1)+5.04D0 
      GF(2)=A(2)+0.035D0
      GF(3)=A(3)+10.D0  
      RETURN    
4     IF (INDEX1(1)) G(1)=Y(2)  
      IF (INDEX1(2)) G(2)=Y(3)  
      IF (INDEX1(3)) G(3)=Y(4)-85.D0    
      IF (INDEX1(4)) G(4)=Y(5)-90.D0    
      IF (INDEX1(5)) G(5)=Y(6)-3.D0     
      IF (INDEX1(6)) G(6)=Y(7)-0.01D0   
      IF (INDEX1(7)) G(7)=Y(8)-145.D0   
      IF (INDEX1(8)) G(8)=5.D+3-Y(2)    
      IF (INDEX1(9)) G(9)=2.D+3-Y(3)    
      IF (INDEX1(10)) G(10)=93.D0-Y(4)  
      IF (INDEX1(11)) G(11)=95.D0-Y(5)  
      IF (INDEX1(12)) G(12)=12.D0-Y(6)  
      IF (INDEX1(13)) G(13)=4.D0-Y(7)   
      IF (INDEX1(14)) G(14)=162.D0-Y(8) 
      RETURN    
5     DO 133 J=1,7      
      IF (.NOT.INDEX2(J)) GOTO 131      
      DO 130 I=1,3      
130   GG(J,I)=DY(J+1,I) 
131   IF (.NOT.INDEX2(J+7)) GOTO 133    
      DO 132 I=1,3      
132   GG(J+7,I)=-DY(J+1,I)      
133   CONTINUE  
      RETURN    
      END       
C
      SUBROUTINE TP68(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D68/ A(2),B(2),Z(2),D(2)     
      DOUBLEPRECISION A,B,Z,D,V1,V2,V3,V4,V5,H1,H2,H3,X1
      INTEGER I,KN1
      NEX=1
      LEX=.FALSE.
      KN1=1     
      XEX(1)=0.678587452312D-01     
      XEX(2)=0.364617174165D+01     
      XEX(3)=0.266175189694D-03     
      XEX(4)=0.894862212037D+00     
      FEX=-0.920425020704D+00 
C      GOTO 9    
C      ENTRY TP69(MODE)  
C      KN1=2     
C      XEX(1)=0.293714180830D-01     
C      XEX(2)=0.119025343488D+01     
C      XEX(3)=0.233946796758D+00     
C      XEX(4)=0.791667815694D+00     
C      FEX=-0.956712887064D+03 
C9     CONTINUE
      GOTO (1,2,3,4,5),MODE     
1     N=4       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=2    
      DO 6 I=1,2
      X(I)=1.D0 
      X(I+2)=1.D0       
      LXL(I)=.TRUE.     
      LXL(I+2)=.TRUE.   
      LXU(I)=.TRUE.     
      LXU(I+2)=.TRUE.   
      XL(I+2)=0.D0      
      XU(I)=100.D0      
6     XU(I+2)=2.D0      
      XL(1)=1.D-4 
      XU(1)=1.0D0
      XL(2)=0.D0
      XL(4)=1.0D-4
      A(1)=1.D-4
      A(2)=0.1D0
      B(1)=1.D0 
      B(2)=1.D+3
      D(1)=1.D0 
      D(2)=1.D0 
      Z(1)=24.D0
      Z(2)=4.D0 
      GG(1,1)=0.D0      
      GG(1,4)=0.D0      
      GG(1,3)=1.D0      
      GG(2,1)=0.D0      
      GG(2,3)=0.D0      
      GG(2,4)=1.D0      
      GF(2)=0.D0
      LEX=.FALSE.       
      NEX=0
      RETURN    
2     CONTINUE
      X1=DMIN1(DMAX1(1.0D-8,X(1)),1.0D1)
      V1=DEXP(X1)-1.D0
      FX=(A(KN1)*Z(KN1)-X(4)*(B(KN1)*V1-X(3))/(V1+X(4)))/X1   
      RETURN    
3     V1=DEXP(X(1))     
      V2=V1-1.D0
      V3=1.D0/(V2+X(4)) 
      V4=1.D0/X(1)      
      V5=(B(KN1)*V2-X(3))*V4    
      GF(1)=-((V1*(X(4)*B(KN1)+X(3))*V3-V5)*X(4)*V3+A(KN1)*Z(KN1)*V4)*V4
      GF(3)=X(4)*V4*V3  
      GF(4)=-V5*V2*V3**2
      RETURN          
4     IF (.NOT.INDEX1(1)) GOTO 30       
      CALL MDNORD(-X(2),H1)
      G(1)=X(3)-2.D0*H1 
  30  IF (.NOT.INDEX1(2)) RETURN
      CALL MDNORD(-X(2)+DSQRT(Z(KN1)),H1)
      CALL MDNORD(-X(2)-DSQRT(Z(KN1)),H2)
      G(2)=X(4)-H1-H2
      RETURN
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,2)=2.D0*DEXP(-0.5D0*X(2)**2)/DSQRT(8.D0*DATAN(1.D0)) 
7     IF (.NOT.INDEX2(2)) GOTO 8
      H1=-X(2)-D(KN1)*DSQRT(Z(KN1))     
      H2=-X(2)+D(KN1)*DSQRT(Z(KN1))     
      H3=1.D0/DSQRT(8.D0*DATAN(1.D0))   
      GG(2,2)=(DEXP(-0.5D0*H1**2)+DEXP(-0.5D0*H2**2))*H3
8     RETURN    
      END       
C
      SUBROUTINE TP69(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D68/ A(2),B(2),Z(2),D(2)     
      DOUBLEPRECISION A,B,Z,D,V1,V2,V3,V4,V5,H1,H2,H3
      INTEGER I,KN1
      NEX=1
      LEX=.FALSE.
      KN1=2     
      XEX(1)=0.293714180830D-01     
      XEX(2)=0.119025343488D+01     
      XEX(3)=0.233946796758D+00     
      XEX(4)=0.791667815694D+00     
      FEX=-0.956712887064D+03 
      GOTO (1,2,3,4,5),MODE     
1     N=4       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=2    
      DO 6 I=1,2
      X(I)=1.D0 
      X(I+2)=1.D0       
      LXL(I)=.TRUE.     
      LXL(I+2)=.TRUE.   
      LXU(I)=.TRUE.     
      LXU(I+2)=.TRUE.   
      XL(I+2)=0.D0      
      XU(I)=100.D0      
6     XU(I+2)=2.D0      
      XL(1)=1.D-4 
      XL(2)=0.D0
      A(1)=1.D-4
      A(2)=0.1D0
      B(1)=1.D0 
      B(2)=1.D+3
      D(1)=1.D0 
      D(2)=1.D0 
      Z(1)=24.D0
      Z(2)=4.D0 
      GG(1,1)=0.D0      
      GG(1,4)=0.D0      
      GG(1,3)=1.D0      
      GG(2,1)=0.D0      
      GG(2,3)=0.D0      
      GG(2,4)=1.D0      
      GF(2)=0.D0
      LEX=.FALSE.       
      NEX=0
      RETURN    
2     V1=DEXP(X(1))-1.D0
      FX=(A(KN1)*Z(KN1)-X(4)*(B(KN1)*V1-X(3))/(V1+X(4)))/X(1)   
      RETURN    
3     V1=DEXP(X(1))     
      V2=V1-1.D0
      V3=1.D0/(V2+X(4)) 
      V4=1.D0/X(1)      
      V5=(B(KN1)*V2-X(3))*V4    
      GF(1)=-((V1*(X(4)*B(KN1)+X(3))*V3-V5)*X(4)*V3+A(KN1)*Z(KN1)*V4)*V4
      GF(3)=X(4)*V4*V3  
      GF(4)=-V5*V2*V3**2
      RETURN          
4     IF (.NOT.INDEX1(1)) GOTO 30       
      CALL MDNORD(-X(2),H1)
      G(1)=X(3)-2.D0*H1 
  30  IF (.NOT.INDEX1(2)) RETURN
      CALL MDNORD(-X(2)+DSQRT(Z(KN1)),H1)
      CALL MDNORD(-X(2)-DSQRT(Z(KN1)),H2)
      G(2)=X(4)-H1-H2
      RETURN
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,2)=2.D0*DEXP(-0.5D0*X(2)**2)/DSQRT(8.D0*DATAN(1.D0)) 
7     IF (.NOT.INDEX2(2)) GOTO 8
      H1=-X(2)-D(KN1)*DSQRT(Z(KN1))     
      H2=-X(2)+D(KN1)*DSQRT(Z(KN1))     
      H3=1.D0/DSQRT(8.D0*DATAN(1.D0))   
      GG(2,2)=(DEXP(-0.5D0*H1**2)+DEXP(-0.5D0*H2**2))*H3
8     RETURN    
      END       
C
      SUBROUTINE TP70(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D70/ DF(19,4),V4(19),U1(19),U2(19),YC(19),C(19),
     /                 V3(19),V8(19),V9(19),F(19),YO(19)       
      DOUBLEPRECISION DF,V4,U1,U2,YC,C,V3,V8,V9,F,YO,B,H1,H2,H3,  
     /     H4,H5,H6,H7,H30,H40,V10,V11,V12,Z1,Z2,Z5,Z6,Z7,
     /     T,H8,H9,H10,H11,H12,H13,H14,H15,H16,H17,H18,H19,
     /     H20,H21,H22,H23,S,SUM   
      INTEGER I,J
      LOGICAL LOG
      LOG=.FALSE.       
      IF (MODE - 2) 1,17,17     
1     N=4       
      NILI=0    
      NINL=1    
      NELI=0    
      NENL=0    
      X(1)=2.D0 
      X(2)=4.D0 
      X(3)=0.04D0       
      X(4)=2.D0 
      DO 6 I=1,4
      LXU(I)=.TRUE.     
      LXL(I)=.TRUE.     
      XL(I)=1.D-5       
6     XU(I)=100.D0      
      XU(3)=1.D0
      GG(1,1)=0.D0      
      GG(1,2)=0.D0      
      XEX(1)=0.122769537557D+02     
      XEX(2)=0.463178796886D+01     
      XEX(3)=0.312862470961D+00     
      XEX(4)=0.202929031570D+01     
      FEX=0.749846356143D-02 
      LEX=.FALSE.       
      NEX=1
      RETURN    
   17 CONTINUE
      C(1)=0.1D0
      DO 20 I=2,19      
20    C(I)=DBLE(I-1)  
      YO(1)=1.89D-3     
      YO(2)=0.1038D0    
      YO(3)=0.268D0     
      YO(4)=0.506D0     
      YO(5)=0.577D0     
      YO(6)=0.604D0     
      YO(7)=0.725D0     
      YO(8)=0.898D0     
      YO(9)=0.947D0     
      YO(10)=0.845D0    
      YO(11)=0.702D0    
      YO(12)=0.528D0    
      YO(13)=0.385D0    
      YO(14)=0.257D0    
      YO(15)=0.159D0    
      YO(16)=0.0869D0   
      YO(17)=0.0453D0   
      YO(18)=0.01509D0  
      YO(19)=1.89D-3    
      IF (MODE - 4) 18,4,5      
18    B=X(3)+(1.D0-X(3))*X(4)   
      H1=X(1)-1.D0      
      H2=X(2)-1.D0      
      H3=1.D0/7.658D0   
      H5=B*H3   
      H4=H5/X(4)
      H6=12.D0*X(1)/(12.D0*X(1)+1.D0)   
      H7=12.D0*X(2)/(12.D0*X(2)+1.D0)   
      H30=0.D0  
      H40=0.D0  
      V10=X(2)/6.2832D0 
      V11=B/X(4)
      V12=X(1)/6.2832D0 
      IF(B.LT.0.D0.OR.V10.LT.0.D0.OR.V11.LT.0.D0.OR.V12.LT.0.D0)
     1     LOG=.TRUE.   
      IF(LOG .AND. MODE.EQ.2) GOTO 8    
      IF (LOG .AND. MODE.EQ.3) GOTO 9   
      Z1=X(3)*B**X(2)   
      Z2=V10**0.5D0     
      Z5=1.D0-X(3)      
      Z6=V11**X(1)      
      Z7=V12**0.5D0     
      DO 30 I=1,19      
      V3(I)=(C(I)*H3)**H2       
      V4(I)=DEXP(X(2)*(1.D0-C(I)*H5))   
      V8(I)=(C(I)*H3)**H1       
      V9(I)=DEXP(X(1)*(1.D0-C(I)*H4))   
      U1(I)=Z1*Z2*V3(I)*V4(I)*H7
      U2(I)=Z5*Z6*Z7*V8(I)*V9(I)*H6     
      YC(I)=U1(I)+U2(I) 
30    F(I)=YC(I)-YO(I)  
      GOTO (1,2,3),MODE 
2     T=0.D0    
      DO 31 I=1,19      
31    T=T+(F(I))**2     
      FX=T      
      RETURN    
3     H8=X(4)-1.D0      
      H9=X(3)-1.D0      
      H10=1.D0/B
      H11=1.D0/X(4)     
      H12=(0.5D0+1.D0/(12.D0*X(1)+1.D0))/X(1)+1.D0      
      H13=(0.5D0+1.D0/(12.D0*X(2)+1.D0))/X(2)+1.D0      
      H16=X(2)*H8       
      H17=X(1)*H8       
      H18=1.D0/X(3)-H16*H10     
      H19=1.D0/H9-H17*H10       
      H16=H16*H3
      H17=H17*H11*H3    
      H20=X(2)*H9       
      H21=X(1)*X(3)*H11**2      
      H22=H20*H3
      H23=H21*H3
      H20=H20*H10       
      H21=H21*H10*X(4)  
      DO 33 J=1,19      
      H14=H4*C(J)       
      H15=H5*C(J)       
      IF (H14.GT.0.D0) GOTO 10  
      GF(1)=2.D0*(X(1)-5.D0)    
      H30=1.D0  
      GOTO 11   
10    DF(J,1)=U2(J)*(H12-H14+DLOG(H14)) 
11    IF (H15.GT.0.D0) GOTO 12  
      GF(2)=2.D0*(X(2)-5.D0)    
      H40=1.D0  
      GOTO 13   
12    DF(J,2)=U1(J)*(H13-H15+DLOG(H15)) 
13    DF(J,3)=U1(J)*(H18+C(J)*H16)+U2(J)*(H19+C(J)*H17) 
33    DF(J,4)=U1(J)*(H22*C(J)-H20)+U2(J)*(H23*C(J)-H21) 
      IF (H30.EQ.1.D0.OR.H40.EQ.1.D0) GOTO 14   
      DO 37 I=1,4       
      S=0.D0    
      DO 36 J=1,19      
      S=S+2.D0*F(J)*DF(J,I)     
36    CONTINUE  
      GF(I)=S   
37    CONTINUE  
      RETURN    
14    DO 38 I=1,4       
      IF (I.EQ.1.AND.H30.EQ.1.D0) GOTO 38       
      IF (I.EQ.2.AND.H40.EQ..1) GOTO 38 
      S=0.D0    
      DO 39 J=1,19      
39    S=S+2.D0*F(J)*DF(J,I)     
      GF(I)=S   
38    CONTINUE  
      H30=0.D0  
      H40=0.D0  
      RETURN    
8     LOG=.FALSE.       
      SUM=0.D0  
      DO 40 I=1,4       
40    SUM=SUM+(X(I)-5.D0)**2    
      FX=SUM
      RETURN    
9     LOG=.FALSE.       
      DO 41 I=1,4       
41    GF(I)=2.D0*(X(I)-5.D0)
      RETURN    
4     IF (INDEX1(1)) G(1)=X(3)+(1.D0-X(3))*X(4) 
      RETURN    
5     IF(.NOT.INDEX2(1)) GOTO 7 
      GG(1,3)=1.D0-X(4) 
      GG(1,4)=-X(3)+1.D0
7     RETURN    
      END       
C
      SUBROUTINE TP71(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=4       
      NILI=0    
      NINL=1    
      NELI=0    
      NENL=1    
      DO 6 I=1,4
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      XL(I)=1.D0
6     XU(I)=5.D0
      X(1)=1.D0 
      X(2)=5.D0 
      X(3)=5.D0 
      X(4)=1.D0 
      XEX(1)=0.100000000000D+01     
      XEX(2)=0.474299937545D+01     
      XEX(3)=0.382115032617D+01     
      XEX(4)=0.137940824585D+01     
      FEX=0.170140172895D+02 
      LEX=.FALSE.       
      NEX=0
      RETURN    
2     FX=X(1)*X(4)*(X(1)+X(2)+X(3))+X(3)
      RETURN    
3     GF(1)=X(4)*(2.D0*X(1)+X(2)+X(3))  
      GF(2)=X(1)*X(4)   
      GF(3)=GF(2)+1.D0  
      GF(4)=X(1)*(X(1)+X(2)+X(3))       
      RETURN    
4     IF (INDEX1(1)) G(1)=(X(1)*X(2)*X(3)*X(4)-25.D0)/25.0D0
      IF (INDEX1(2)) G(2)=(X(1)**2+X(2)**2+X(3)**2+X(4)**2-40.D0)/40.0D0
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=X(2)*X(3)*X(4)/25.0D0    
      GG(1,2)=X(1)*X(3)*X(4)/25.0D0    
      GG(1,3)=X(1)*X(2)*X(4)/25.0D0    
      GG(1,4)=X(1)*X(2)*X(3)/25.0D0    
7     IF (.NOT.INDEX2(2)) GOTO 8
      DO 20 I=1,4       
20    GG(2,I)=2.D0*X(I)/40.0D0  
8     RETURN    
      END       
C
      SUBROUTINE TP72(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=4       
      NILI=0    
      NINL=2    
      NELI=0    
      NENL=0    
      DO 6 I=1,4
      X(I)=1.D0 
      XL(I)=1.D-3       
      XU(I)=1.D+5*(5.D0-DBLE(I))      
      LXL(I)=.TRUE.     
6     LXU(I)=.TRUE.     
      DO 30 I=1,4       
30    GF(I)=1.D0
      LEX=.FALSE.
      NEX=1       
      XEX(1)=0.193407050141D+03     
      XEX(2)=0.179547504555D+03     
      XEX(3)=0.185018587841D+03     
      XEX(4)=0.168706233485D+03     
C      FEX=0.727679376021D+03 
      FEX=0.72767936D+3
      RETURN    
2     FX=1.D0+X(1)+X(2)+X(3)+X(4)       
3     RETURN    
4     IF (INDEX1(1)) G(1)=-4.D0/X(1)-2.25D0/X(2)-1.D0/X(3)-0.25D0/X(4)  
     1     +0.0401D0    
      IF (INDEX1(2)) G(2)=-0.16D0/X(1)-0.36D0/X(2)-0.64D0*(1.D0/X(3)    
     1     +1.D0/X(4))+0.010085D0       
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=4.D0/X(1)**2      
      GG(1,2)=2.25D0/X(2)**2    
      GG(1,3)=1.D0/X(3)**2      
      GG(1,4)=0.25D0/X(4)**2    
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,1)=0.16D0/X(1)**2    
      GG(2,2)=0.36D0/X(2)**2    
      GG(2,3)=0.64D0/X(3)**2    
      GG(2,4)=0.64D0/X(4)**2    
8     RETURN    
      END       
C
      SUBROUTINE TP73(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DIMENSION A(4)    
      DOUBLEPRECISION A  
      GOTO (1,2,3,4,5),MODE     
1     N=4       
      NILI=1    
      NINL=1    
      NELI=1    
      NENL=0    
      DO 6 I=1,4
      X(I)=1.D0 
      XL(I)=0.D0
      LXL(I)=.TRUE.     
6     LXU(I)=.FALSE.    
      GF(1)=24.55D0     
      GF(2)=26.75D0     
      GF(3)=39.D0       
      GF(4)=40.5D0      
      GG(1,1)=2.3D0     
      GG(1,2)=5.6D0     
      GG(1,3)=11.1D0    
      GG(1,4)=1.3D0     
      DO 31 I=1,4       
31    GG(3,I)=1.D0      
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.635521568605D+00     
      XEX(2)=-0.117862273760D-11     
      XEX(3)=0.312701880754D+00     
      XEX(4)=0.517765506011D-01     
      FEX=0.298943781573D+02 
      RETURN    
2     FX=24.55D0*X(1)+26.75D0*X(2)+39.D0*X(3)+40.50D0*X(4)      
3     RETURN    
4     IF (INDEX1(1)) G(1)=2.3D0*X(1)+5.6D0*X(2)+11.1D0*X(3)+1.3D0*X(4)  
     1     -5.D0
      IF (INDEX1(2)) G(2)=12.0D0*X(1)+11.9D0*X(2)+41.8D0*X(3)+52.1D0    
     /   *X(4)-1.645D0*    
     /   (0.28D0*X(1)**2+0.19D0*X(2)**2+20.5D0*X(3)**2+0.62D0*X(4)**2) 
     /   **.5D0-21.D0      
      IF (INDEX1(3)) G(3)=X(1)+X(2)+X(3)+X(4)-1.D0      
      RETURN    
5     IF (.NOT.INDEX2(2)) GOTO 8
      DO 30 I=1,4       
30    A(I)=1.645D0*X(I)*(0.28D0*X(1)**2+0.19D0*X(2)**2+20.5D0*X(3)**2   
     /     +0.62D0*X(4)**2)**(-0.5D0)
      GG(2,1)=12.D0-0.28D0*A(1) 
      GG(2,2)=11.9D0-0.19D0*A(2)
      GG(2,3)=41.8D0-20.5D0*A(3)
      GG(2,4)=52.1D0-0.62D0*A(4)
8     RETURN    
      END       
C
      SUBROUTINE TP74(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /L74/ A(2)    
      DOUBLEPRECISION A,V1,V2  
      INTEGER I,KN1
      XEX(1)=0.679945319802D+03     
      XEX(2)=0.102606713256D+04     
      XEX(3)=0.118876364490D+00     
      XEX(4)=-0.396233553180D+00     
      FEX=0.512649810934D+04 
      KN1=1     
      GOTO 7    
      ENTRY TP75(MODE)  
      KN1=2     
      XEX(1)=0.776159220293D+03     
      XEX(2)=0.925194939196D+03     
      XEX(3)=0.511087936804D-01     
      XEX(4)=-0.428891137432D+00     
      FEX=0.517441288686D+04 
7     GOTO (1,2,3,4,5),MODE    
1     N=4       
      NILI=2    
      NINL=0    
      NELI=0    
      NENL=3    
      A(1)=0.55D0       
      A(2)=0.48D0       
      DO 6 I=1,4
      X(I)=0.D0 
      LXL(I)=.TRUE.     
6     LXU(I)=.TRUE.     
      XL(1)=0.D0
      XL(2)=0.D0
      XL(3)=-A(KN1)     
      XL(4)=-A(KN1)     
      XU(1)=1.2D+3      
      XU(2)=1.2D+3      
      XU(3)=A(KN1)      
      XU(4)=A(KN1)      
      GF(3)=0.D0
      GF(4)=0.D0
      GG(1,1)=0.D0      
      GG(1,2)=0.D0      
      GG(1,3)=-1.D0     
      GG(1,4)=1.D0      
      GG(2,1)=0.D0      
      GG(2,2)=0.D0      
      GG(2,3)=1.D0      
      GG(2,4)=-1.D0     
      GG(3,1)=-1.D0     
      GG(3,2)=0.D0      
      GG(4,1)=0.D0      
      GG(4,2)=-1.D0     
      GG(5,1)=0.D0      
      GG(5,2)=0.D0      
      LEX=.FALSE.
      NEX=0       
      RETURN    
2     FX=3.D0*X(1)+X(1)**3*1.D-6+2.D0*X(2)+2.D0*1.D-6/3.D0*X(2)**3      
      RETURN    
3     GF(1)=3.D0+3.D-6*X(1)**2  
      GF(2)=2.D0+2.D-6*X(2)**2  
      RETURN    
4     IF (INDEX1(1)) G(1)=X(4)-X(3)+A(KN1)      
      IF (INDEX1(2)) G(2)=X(3)-X(4)+A(KN1)      
      IF (INDEX1(3)) G(3)=1.D+3*(DSIN(-X(3)-0.25D0)+DSIN(-X(4)-0.25D0))+
     -894.8D0-X(1)      
      IF (INDEX1(4)) G(4)=1.D+3*(DSIN(X(3)-0.25D0)+DSIN(X(3)-X(4)       
     --0.25D0))+
     -894.8D0-X(2)      
      IF(INDEX1(5)) G(5)=1.D+3*(DSIN(X(4)-0.25D0)+DSIN(X(4)-X(3)
     --0.25D0))+
     -1.2948D+3 
      RETURN    
5     IF (.NOT.INDEX2(3)) GOTO 9
      GG(3,3)=-1.D+3*DCOS(-X(3)-0.25D0) 
      GG(3,4)=-1.D+3*DCOS(-X(4)-0.25D0) 
9     IF (.NOT.INDEX2(4)) GOTO 10       
      V1=DCOS(X(3)-X(4)-0.25D0) 
      GG(4,3)=1.D+3*(DCOS(X(3)-0.25D0)+V1)      
      GG(4,4)=-1.D+3*V1 
10    IF (.NOT.INDEX2(5)) GOTO 11       
      V2=DCOS(X(4)-X(3)-0.25D0) 
      GG(5,3)=-1.D+3*V2 
      GG(5,4)=1.D+3*(DCOS(X(4)-0.25D0)+V2)      
11    RETURN    
      END       
C
      SUBROUTINE TP76(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=4       
      NILI=3    
      NINL=0    
      NELI=0    
      NENL=0    
      DO 6 I=1,4
      LXL(I)=.TRUE.     
      LXU(I)=.FALSE.    
      X(I)=0.5D0
6     XL(I)=0.D0
      GG(1,1)=-1.D0     
      GG(1,2)=-2.D0     
      GG(1,3)=-1.D0     
      GG(1,4)=-1.D0     
      GG(2,1)=-3.D0     
      GG(2,2)=-1.D0     
      GG(2,3)=-2.D0     
      GG(2,4)=1.D0      
      GG(3,1)=0.D0      
      GG(3,2)=1.D0      
      GG(3,3)=4.D0      
      GG(3,4)=0.D0      
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.272727272717D+00     
      XEX(2)=0.209090909094D+01     
      XEX(3)=-0.263371889808D-10     
      XEX(4)=0.545454545496D+00     
      FEX=-0.468181818182D+01 
      RETURN    
2     FX=X(1)**2+X(3)**2+0.5D0*(X(2)**2+X(4)**2)-X(1)*X(3)+X(3)*X(4)-   
     -X(1)-3.D0*X(2)+X(3)-X(4)  
      RETURN    
3     GF(1)=2.D0*X(1)-X(3)-1.D0 
      GF(2)=X(2)-3.D0   
      GF(3)=2.D0*X(3)-X(1)+X(4)+1.D0    
      GF(4)=X(4)+X(3)-1.D0      
      RETURN    
4     IF (INDEX1(1)) G(1)=-X(1)-2.D0*X(2)-X(3)-X(4)+5.D0
      IF (INDEX1(2)) G(2)=-3.D0*X(1)-X(2)-2.D0*X(3)+X(4)+4.D0   
      IF (INDEX1(3)) G(3)=X(2)+4.D0*X(3)-1.5D0  
5     RETURN    
      END       
C
      SUBROUTINE TP77(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION V1 
      GOTO (1,2,3,4,5),MODE     
1     N=5       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=2    
      DO 6 I=1,5
      X(I)=2.D0 
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      GG(1,2)=0.D0      
      GG(1,3)=0.D0      
      GG(2,1)=0.D0      
      GG(2,2)=1.D0      
      GG(2,5)=0.D0      
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.116617219726D+01     
      XEX(2)=0.118211138813D+01     
      XEX(3)=0.138025704044D+01     
      XEX(4)=0.150603627961D+01     
      XEX(5)=0.610920257517D+00     
      FEX=0.241505128786D+00 
      RETURN    
2     FX=(X(1)-1.D0)**2+(X(1)-X(2))**2+(X(3)-1.D0)**2+(X(4)-1.D0)**4    
     -+(X(5)-1.D0)**6   
      RETURN    
3     GF(1)=2.D0*(2.D0*X(1)-X(2)-1.D0)  
      GF(2)=-2.D0*(X(1)-X(2))   
      GF(3)=2.D0*(X(3)-1.D0)    
      GF(4)=4.D0*(X(4)-1.D0)**3 
      GF(5)=6.D0*(X(5)-1.D0)**5 
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)**2*X(4)+DSIN(X(4)-X(5))-2.D0*DSQRT(2.D0) 
      IF (INDEX1(2)) G(2)=X(2)+X(3)**4*X(4)**2-8.D0-DSQRT(2.D0) 
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      V1=DCOS(X(4)-X(5))
      GG(1,1)=2.D0*X(1)*X(4)    
      GG(1,4)=X(1)**2+V1
      GG(1,5)=-V1       
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,3)=4.D0*X(3)**3*X(4)**2      
      GG(2,4)=2.D0*X(3)**4*X(4) 
8     RETURN    
      END       
C
      SUBROUTINE TP78(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=5       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=3    
      DO 6 I=1,5
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      X(1)=-2.D0
      X(2)=1.5D0
      X(3)=2.D0 
      X(4)=-1.D0
      X(5)=-1.D0
      GG(2,1)=0.D0      
      GG(3,3)=0.D0      
      GG(3,4)=0.D0      
      GG(3,5)=0.D0      
      LEX=.FALSE.       
      NEX=1
      XEX(1)=-0.171714234230D+01     
      XEX(2)=0.159570826805D+01     
      XEX(3)=0.182724803488D+01     
      XEX(4)=-0.763642946600D+00     
      XEX(5)=-0.763643482853D+00     
      FEX=-0.291970040911D+01 
      RETURN    
2     FX=X(1)*X(2)*X(3)*X(4)*X(5)       
      RETURN    
3     GF(1)=X(2)*X(3)*X(4)*X(5) 
      GF(2)=X(1)*X(3)*X(4)*X(5) 
      GF(3)=X(1)*X(2)*X(4)*X(5) 
      GF(4)=X(1)*X(2)*X(3)*X(5) 
      GF(5)=X(1)*X(2)*X(3)*X(4) 
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)**2+X(2)**2+X(3)**2+X(4)**2+X(5)**2-10.D0 
      IF (INDEX1(2)) G(2)=X(2)*X(3)-5.D0*X(4)*X(5)      
      IF (INDEX1(3)) G(3)=X(1)**3+X(2)**3+1.D0  
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      DO 30 I=1,5       
30    GG(1,I)=2.D0*X(I) 
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,2)=X(3)      
      GG(2,3)=X(2)      
      GG(2,4)=-5.D0*X(5)
      GG(2,5)=-5.D0*X(4)
8     IF (.NOT.INDEX2(3)) GOTO 9
      GG(3,1)=3.D0*X(1)**2      
      GG(3,2)=3.D0*X(2)**2      
9     RETURN    
      END       
C
      SUBROUTINE TP79(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION V1,V2,V3,V4
      GOTO (1,2,3,4,5),MODE     
1     N=5       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=3    
      DO 6 I=1,5
      X(I)=2.D0 
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      GG(1,1)=1.D0      
      GG(1,4)=0.D0      
      GG(1,5)=0.D0      
      GG(2,1)=0.D0      
      GG(2,2)=1.D0      
      GG(2,4)=1.D0      
      GG(2,5)=0.D0      
      GG(3,2)=0.D0      
      GG(3,3)=0.D0      
      GG(3,4)=0.D0      
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.119112745626D+01     
      XEX(2)=0.136260316492D+01     
      XEX(3)=0.147281793150D+01     
      XEX(4)=0.163501661894D+01     
      XEX(5)=0.167908143619D+01     
      FEX=0.787768208538D-01 
      RETURN    
2     FX=(X(1)-1.D0)**2+(X(1)-X(2))**2+(X(2)-X(3))**2+(X(3)-X(4))**4    
     1     +(X(4)-X(5))**4      
      RETURN    
3     V1=X(1)-X(2)      
      V2=X(2)-X(3)      
      V3=X(3)-X(4)      
      V4=X(4)-X(5)      
      GF(1)=2.D0*(X(1)-1.D0+V1) 
      GF(2)=2.D0*(V2-V1)
      GF(3)=-2.D0*V2+4.D0*V3**3 
      GF(4)=4.D0*(V4**3-V3**3)  
      GF(5)=-4.D0*V4**3 
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)+X(2)**2+X(3)**3-2.D0-3.D0*DSQRT(2.D0)    
      IF (INDEX1(2)) G(2)=X(2)-X(3)**2+X(4)+2.D0-2.D0*DSQRT(2.D0)       
      IF (INDEX1(3)) G(3)=X(1)*X(5)-2.D0
      RETURN    
5     IF(.NOT.INDEX2(1)) GOTO 7 
      GG(1,2)=2.D0*X(2) 
      GG(1,3)=3.D0*X(3)**2      
7     IF (INDEX2(2)) GG(2,3)=-2.D0*X(3) 
      IF (.NOT.INDEX2(3)) GOTO 9
      GG(3,1)=X(5)      
      GG(3,5)=X(1)      
9     RETURN    
      END       
C
      SUBROUTINE TP80(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION V1,V2,T
      GOTO (1,2,3,4,5),MODE     
1     N=5       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=3    
      DO 6 I=1,5
      LXL(I)=.TRUE.     
6     LXU(I)=.TRUE.     
      X(1)=-2.D0
      X(2)=2.D0 
      X(3)=2.D0 
      X(4)=-1.D0
      X(5)=-1.D0
      DO 20 I=1,2       
      XL(I)=-2.3D0      
20    XU(I)=2.3D0       
      DO 21 I=3,5       
      XL(I)=-3.2D0      
21    XU(I)=3.2D0       
      GG(2,1)=0.D0      
      GG(3,3)=0.D0      
      GG(3,4)=0.D0      
      GG(3,5)=0.D0      
      LEX=.FALSE.       
      NEX=1
      XEX(1)=-0.171714294417D+01     
      XEX(2)=0.159570896503D+01     
      XEX(3)=0.182724691654D+01     
      XEX(4)=-0.763641279311D+00     
      XEX(5)=-0.763645016315D+00     
      FEX=0.539498477624D-01 
      RETURN    
2     FX=DEXP(X(1)*X(2)*X(3)*X(4)*X(5)) 
      RETURN    
3     V1=X(4)*X(5)      
      V2=X(1)*X(2)      
      T=DEXP(V2*X(3)*V1)
      GF(1)=X(2)*X(3)*V1*T      
      GF(2)=X(1)*X(3)*V1*T      
      GF(3)=V1*V2*T     
      GF(4)=V2*X(3)*X(5)*T      
      GF(5)=V2*X(3)*X(4)*T      
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)**2+X(2)**2+X(3)**2+X(4)**2+X(5)**2-10.D0 
      IF (INDEX1(2)) G(2)=X(2)*X(3)-5.D0*X(4)*X(5)      
      IF (INDEX1(3)) G(3)=X(1)**3+X(2)**3+1.D0  
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      DO 30 I=1,5       
30    GG(1,I)=2.D0*X(I) 
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,2)=X(3)      
      GG(2,3)=X(2)      
      GG(2,4)=-5.D0*X(5)
      GG(2,5)=-5.D0*X(4)
8     IF (.NOT.INDEX2(3)) GOTO 9
      GG(3,1)=3.D0*X(1)**2      
      GG(3,2)=3.D0*X(2)**2      
9     RETURN    
      END       
C
      SUBROUTINE TP81(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION V1,V2,V3,T    
      GOTO (1,2,3,4,5),MODE     
1     N=5       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=3    
      DO 6 I=1,5
      LXL(I)=.TRUE.     
6     LXU(I)=.TRUE.     
      DO 10 I=1,2       
      XL(I)=-2.3D0      
10    XU(I)=2.3D0       
      DO 20 I=3,5       
      XL(I)=-3.2D0      
20    XU(I)=3.2D0       
      X(1)=-2.D0
      X(2)=2.D0 
      X(3)=2.D0 
      X(4)=-1.D0
      X(5)=-1.D0
      GG(2,1)=0.D0      
      GG(3,3)=0.D0      
      GG(3,4)=0.D0      
      GG(3,5)=0.D0      
      LEX=.FALSE.       
      NEX=1
      XEX(1)=-0.171714240091D+01     
      XEX(2)=0.159570833592D+01     
      XEX(3)=0.182724792592D+01     
      XEX(4)=-0.763647440817D+00     
      XEX(5)=-0.763638975604D+00     
      FEX=0.539498477749D-01 
      RETURN    
2     FX=DEXP(X(1)*X(2)*X(3)*X(4)*X(5))-0.5D0*(X(1)**3+X(2)**3+1.D0)**2 
      RETURN    
3     V1=X(1)**3+X(2)**3+1.D0   
      V2=X(1)*X(2)      
      V3=X(4)*X(5)      
      T=DEXP(V2*V3*X(3))
      GF(1)=X(2)*X(3)*V3*T-3.D0*X(1)**2*V1      
      GF(2)=X(1)*X(3)*V3*T-3.D0*X(2)**2*V1      
      GF(3)=V2*V3*T     
      GF(4)=V2*X(3)*X(5)*T      
      GF(5)=V2*X(3)*X(4)*T      
      RETURN    
4     IF (INDEX1(1)) G(1)=X(1)**2+X(2)**2+X(3)**2+X(4)**2+X(5)**2-10.D0 
      IF (INDEX1(2)) G(2)=X(2)*X(3)-5.D0*X(4)*X(5)      
      IF (INDEX1(3)) G(3)=X(1)**3+X(2)**3+1.D0  
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      DO 30 I=1,5       
30    GG(1,I)=2.D0*X(I) 
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,2)=X(3)      
      GG(2,3)=X(2)      
      GG(2,4)=-5.D0*X(5)
      GG(2,5)=-5.D0*X(4)
8     IF(.NOT.INDEX2(3)) GOTO 9 
      GG(3,1)=3.D0*X(1)**2      
      GG(3,2)=3.D0*X(2)**2      
9     RETURN    
      END       
C
      SUBROUTINE TP83(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      COMMON/DATA83/A,B,C,D,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,
     /              V1,V2,V3
      DOUBLEPRECISION A,B,C,D,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,
     /              V1,V2,V3    
      GOTO (1,2,3,4,5),MODE     
1     N=5       
      NILI=0    
      NINL=6    
      NELI=0    
      NENL=0    
      X(1)=78.0D0
      X(2)=33.0D0
      X(3)=27.0D0
      X(4)=27.0D0
      X(5)=27.0D0
      DO 6 I=1,5
      LXL(I)=.TRUE.     
6     LXU(I)=.TRUE.     
      XL(1)=78.0D0       
      XL(2)=33.0D0       
      XU(1)=102.0D0      
      XU(2)=45.D0       
      DO 31 I=3,5       
      XL(I)=27.0D0       
31    XU(I)=45.0D0       
      A=5.3578547D0     
      B=0.8356891D0     
      C=37.293239D0     
      D=4.0792141D+4    
      A1=85.334407D0    
      A2=5.6858D-3      
      A3=6.262D-4       
      A4=2.2053D-3      
      A5=80.51249D0     
      A6=7.1317D-3      
      A7=2.9955D-3      
      A8=2.1813D-3      
      A9=9.300961D0     
      A10=4.7026D-3     
      A11=1.2547D-3     
      A12=1.9085D-3     
      GF(2)=0.D0
      GF(4)=0.D0
      GG(2,4)=0.D0      
      GG(3,2)=0.D0      
      GG(5,4)=0.D0      
      GG(6,2)=0.D0      
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.780000000000D+02     
      XEX(2)=0.330000000000D+02     
      XEX(3)=0.299952560253D+02     
      XEX(4)=0.450000000000D+02     
      XEX(5)=0.367758129081D+02     
      FEX=-0.306655386717D+05
      RETURN    
2     FX=A*X(3)**2+B*X(1)*X(5)+C*X(1)-D 
      RETURN    
3     GF(1)=B*X(5)+C    
      GF(3)=2.D0*A*X(3) 
      GF(5)=B*X(1) 
      RETURN    
4     IF (.NOT.(INDEX1(1).OR.INDEX1(4))) GOTO 41
      V1=A1+A2*X(2)*X(5)+A3*X(1)*X(4)-A4*X(3)*X(5)      
      IF (INDEX1(1)) G(1)=V1    
      IF (INDEX1(4)) G(4)=92.D0-V1      
41    IF (.NOT.(INDEX1(2).OR.INDEX1(5))) GOTO 42
      V2=A5+A6*X(2)*X(5)+A7*X(1)*X(2)+A8*X(3)**2-90.D0  
      IF (INDEX1(2)) G(2)=V2    
      IF (INDEX1(5)) G(5)=20.D0-V2      
42    IF (.NOT.(INDEX1(3).OR.INDEX1(6))) GOTO 43
      V3=A9+A10*X(3)*X(5)+A11*X(1)*X(3)+A12*X(3)*X(4)-20.D0      
      IF (INDEX1(3)) G(3)=V3    
      IF (INDEX1(6)) G(6)=5.0D0-V3       
43    RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=A3*X(4)   
      GG(1,2)=A2*X(5)   
      GG(1,3)=-A4*X(5)  
      GG(1,4)=A3*X(1)   
      GG(1,5)=A2*X(2)-A4*X(3)   
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,1)=A7*X(2)   
      GG(2,2)=A6*X(5)+A7*X(1)   
      GG(2,3)=2.D0*A8*X(3)      
      GG(2,5)=A6*X(2)   
8     IF(.NOT.INDEX2(3)) GOTO 9 
      GG(3,1)=A11*X(3)  
      GG(3,3)=A10*X(5)+A11*X(1)+A12*X(4)
      GG(3,4)=A12*X(3)  
      GG(3,5)=A10*X(3)  
9     IF (.NOT.INDEX2(4)) GOTO 10       
      GG(4,1)=-A3*X(4)  
      GG(4,2)=-A2*X(5)  
      GG(4,3)=A4*X(5)   
      GG(4,4)=-A3*X(1)  
      GG(4,5)=-A2*X(2)+A4*X(3)  
10    IF (.NOT.INDEX2(5)) GOTO 11       
      GG(5,1)=-A7*X(2)  
      GG(5,2)=-A6*X(5)-A7*X(1)  
      GG(5,3)=-2.D0*A8*X(3)     
      GG(5,5)=-A6*X(2)  
11    IF (.NOT.INDEX2(6)) GOTO 12       
      GG(6,1)=-A11*X(3) 
      GG(6,3)=-A10*X(5)-A11*X(1)-A12*X(4)       
      GG(6,4)=-A12*X(3) 
      GG(6,5)=-A10*X(3) 
12    RETURN    
      END       
C
      SUBROUTINE TP84(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J,I1
      COMMON     /D84/ A(21),B(3)      
      DOUBLEPRECISION A,B,V1     
      A(1)=-2.4345D+4   
      A(2)=-8.720288849D+6      
      A(3)=1.505125253D+5       
      A(4)=-1.566950325D+2      
      A(5)=4.764703222D+5       
      A(6)=7.294828271D+5       
      A(7)=-1.45421402D+5       
      A(8)=2.9311506D+3 
      A(9)=-40.427932D0 
      A(10)=5.106192D+3 
      A(11)=1.571136D+4 
      A(12)=-1.550111084D+5     
      A(13)=4.36053352D+3       
      A(14)=12.9492344D0
      A(15)=1.0236884D+4
      A(16)=1.3176786D+4
      A(17)=-3.266695104D+5     
      A(18)=7.39068412D+3       
      A(19)=-27.8986976D0       
      A(20)=1.6643076D+4
      A(21)=3.0988146D+4
      B(1)=2.94D+5      
      B(2)=2.94D+5      
      B(3)=2.772D+5     
      GOTO (1,2,3,4,5),MODE     
1     N=5       
      NILI=0    
      NINL=6    
      NELI=0    
      NENL=0    
      X(1)=2.52D0       
      X(2)=2.D0 
      X(3)=37.5D0       
      X(4)=9.25D0       
      X(5)=6.8D0
      DO 6 I=1,5
      LXU(I)=.TRUE.     
6     LXL(I)=.TRUE.     
      XL(1)=0.D0
      XL(2)=1.2D0       
      XL(3)=20.D0       
      XL(4)=9.D0
      XL(5)=6.5D0       
      XU(1)=1.D+3       
      XU(2)=2.4D0       
      XU(3)=60.D0       
      XU(4)=9.3D0       
      XU(5)=7.D0
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.453743097466D+01     
      XEX(2)=0.240000000002D+01     
      XEX(3)=0.600000000000D+02     
      XEX(4)=0.929999999999D+01     
      XEX(5)=0.700000000000D+01     
      FEX=-0.528033513306D+07
      RETURN    
2     FX=-(A(1)+X(1)*(A(2)+A(3)*X(2)+A(4)*X(3)+A(5)*X(4)+A(6)*X(5)))  
      RETURN    
3     GF(1)=-(A(2)+A(3)*X(2)+A(4)*X(3)+A(5)*X(4)+A(6)*X(5))
      DO 30 I=2,5       
30    GF(I)=-A(1+I)*X(1)
      RETURN    
4     DO 80 I=1,3       
      IF (.NOT.(INDEX1(I).OR.INDEX1(I+3))) GOTO 80      
      I1=I*5    
      V1=X(1)*(A(I1+2)+A(I1+3)*X(2)+A(I1+4)*X(3)+A(I1+5)*X(4)    
     /              +A(I1+6)*X(5))    
      IF (INDEX1(I)) G(I)=V1
      IF (INDEX1(I+3)) G(I+3)=(B(I)-V1)
80    CONTINUE  
      RETURN    
5     DO 90 I=1,3       
      IF (.NOT.(INDEX2(I).OR.INDEX2(I+3))) GOTO 90      
      I1=5*I+1  
      IF (.NOT.INDEX2(I)) GOTO 95       
      GG(I,1)=(A(I1+1)+A(I1+2)*X(2)+A(I1+3)*X(3)+A(I1+4)*X(4)     
     /+A(I1+5)*X(5))
      DO 91 J=2,5       
91    GG(I,J)=A(I1+J)*X(1)
      IF (.NOT.INDEX2(I+3)) GOTO 90     
      DO 92 J=1,5       
92    GG(I+3,J)=-GG(I,J)
      GOTO 90   
95    GG(I+3,1)=-(A(I1+1)+A(I1+2)*X(2)+A(I1+3)*X(3)+A(I1+4)      
     /*X(4)+A(I1+5)*X(5))
      DO 96 J=2,5       
96    GG(I+3,J)=-A(I1+J)*X(1)
90    CONTINUE  
      RETURN    
      END       
C
      SUBROUTINE TP85(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      COMMON     /D85/C(17),Y(17),DC(17,5),DY(17,5),A(17),B(17)  
      DOUBLEPRECISION C,Y,DC,DY,A,B,V1,V2,V3,V4,V5,V6,V7
      IF (MODE - 2) 1,17,17     
1     N=5       
      NILI=3    
      NINL=35   
      NELI=0    
      NENL=0    
      X(1)=900.D0       
      X(2)=80.D0
      X(3)=115.D0       
      X(4)=267.D0       
      X(5)=27.D0
      DO 6 I=1,5
      LXL(I)=.TRUE.     
6     LXU(I)=.TRUE.     
      XL(1)=704.4148D0  
      XL(2)=68.6D0      
      XL(3)=0.D0
      XL(4)=193.D0      
      XL(5)=25.D0       
      XU(1)=906.3855D0  
      XU(2)=288.88D0    
      XU(3)=134.75D0    
      XU(4)=287.0966D0  
      XU(5)=84.1988D0   
      A(2)=17.505D0     
      A(3)=11.275D0     
      A(4)=214.228D0    
      A(5)=7.458D0      
      A(6)=0.961D0      
      A(7)=1.612D0      
      A(8)=0.146D0      
      A(9)=107.99D0     
      A(10)=922.693D0   
      A(11)=926.832D0   
      A(12)=18.766D0    
      A(13)=1.072163D+3 
      A(14)=8.961448D+3 
      A(15)=0.063D0     
      A(16)=7.108433D+4 
      A(17)=2.802713D+6 
      B(2)=1.0536667D+3 
      B(3)=35.03D0      
      B(4)=665.585D0    
      B(5)=584.463D0    
      B(6)=265.916D0    
      B(7)=7.046D0      
      B(8)=0.222D0      
      B(9)=273.366D0    
      B(10)=1.286105D+3 
      B(11)=1.444046D+3 
      B(12)=537.141D0   
      B(13)=3.247039D+3 
      B(14)=2.6844086D+4
      B(15)=0.386D0     
      B(16)=1.4D+5      
      B(17)=1.2146108D+7
      DO 61 I=1,5       
      GG(1,I)=0.D0      
      GG(2,I)=0.D0      
      GG(3,I)=0.D0      
      DC(1,I)=0.D0      
      DC(5,I)=0.D0      
      DC(10,I)=0.D0     
61    CONTINUE  
      GG(1,2)=1.5D0     
      GG(1,3)=-1.D0     
      GG(2,2)=1.D0      
      GG(2,3)=1.D0      
      GG(3,2)=-1.D0     
      GG(3,3)=-1.D0     
      DY(1,1)=0.D0      
      DY(1,2)=1.D0      
      DY(1,3)=1.D0      
      DY(1,4)=0.D0      
      DY(1,5)=0.D0      
      DC(1,4)=0.024D0   
      DC(5,2)=100.D0    
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.705180328772D+03     
      XEX(2)=0.686000529425D+02     
      XEX(3)=0.102900013236D+03     
      XEX(4)=0.282324998587D+03     
      XEX(5)=0.375850413432D+02     
      FEX=-0.19051553D+01
      RETURN    
17    Y(1)=X(2)+X(3)+41.6D0     
      C(1)=0.024D0*X(4)-4.62D0  
      Y(2)=12.5D0/C(1)+12.0D0   
      V3=Y(2)*X(1)      
      C(2)=(3.535D-4*X(1)+0.5311D0)*X(1)+0.08705D0*V3   
      C(3)=0.052D0*X(1)+78.D0+2.377D-3*V3       
      Y(3)=C(2)/C(3)    
      Y(4)=19.D0*Y(3)   
      V1=X(1)-Y(3)      
      C(4)=(0.1956D0*V1/X(2)+0.04782D0)*V1+0.6376D0*Y(4)+1.594D0*Y(3)  
      C(5)=100.D0*X(2)  
      C(6)=V1-Y(4)      
      C(7)=0.95D0-C(4)/C(5)     
      Y(5)=C(6)*C(7)    
      V2=Y(5)+Y(4)      
      Y(6)=V1-V2
      C(8)=0.995D0*V2   
      Y(7)=C(8)/Y(1)    
      Y(8)=C(8)/3.798D+3
      C(9)=Y(7)-0.0663D0*Y(7)/Y(8)-0.3153D0     
      Y(9)=96.82D0/C(9)+0.321D0*Y(1)    
      Y(10)=1.29D0*Y(5)+1.258D0*Y(4)+2.29D0*Y(3)+1.71D0*Y(6)    
      Y(11)=1.71D0*X(1)-0.452D0*Y(4)+0.58D0*Y(3)
      C(10)=12.3D0/752.3D0      
      C(11)=1.74125D0*V3
      C(12)=0.995D0*Y(10)+1.998D+3      
      Y(12)=C(10)*X(1)+C(11)/C(12)      
      Y(13)=C(12)-1.75D0*Y(2)   
      V4=Y(9)+X(5)      
      Y(14)=3.623D+3+64.4D0*X(2)+58.4D0*X(3)+1.46312D+5/V4      
      C(13)=0.995D0*Y(10)+60.8D0*X(2)+48.D0*X(4)-0.1121D0*Y(14) 
     /      -5.095D+3 
      Y(15)=Y(13)/C(13) 
      Y(16)=1.48D+5-3.31D+5*Y(15)+40.D0*Y(13)-61.D0*Y(15)*Y(13)
      C(14)=2.324D+3*Y(10)-2.874D+7*Y(2)
      Y(17)=1.413D+7-1.328D+3*Y(10)-531.D0*Y(11)+C(14)/C(12) 
      C(15)=Y(13)/Y(15)-Y(13)/0.52D0    
      C(16)=1.104D0-0.72D0*Y(15)
      C(17)=V4  
      IF (MODE.EQ.3.OR.MODE.EQ.5) GOTO 71       
      GOTO (1,2,3,4,5),MODE     
71    DO 30 I=1,5       
30    DY(2,I)=-12.5D0*DC(1,I)/C(1)**2   
      V5=Y(2)+X(1)*DY(2,1)      
      DC(2,1)=7.07D-4*X(1)+0.5311D0+0.08705D0*V5
      DC(3,1)=0.052D0+2.377D-3*V5       
      DO 32 I=2,5       
      V6=X(1)*DY(2,I)   
      DC(2,I)=0.08705D0*V6      
32    DC(3,I)=2.377D-3*V6       
      DO 33 I=1,5       
      DY(3,I)=(C(3)*DC(2,I)-C(2)*DC(3,I))/C(3)**2       
33    DY(4,I)=19.D0*DY(3,I)     
      DC(4,1)=(0.04782D0+0.3912D0*V1/X(2))*(1.D0-DY(3,1))+      
     -0.6376D0*DY(4,1)  
     -+1.594D0*DY(3,1)  
      DC(4,2)=-0.1956D0*V1*(V1+2.D0*X(2)*DY(3,2))/X(2)**2+1     
     -.54618D0  
     -*DY(3,2)+0.6376D0*DY(4,2) 
      DO 34 I=3,5       
34    DC(4,I)=(1.54618D0-0.3912D0*V1/X(2))*DY(3,I)+0.6376D0*D   
     -Y(4,I)    
      DC(6,1)=1.D0-DY(3,1)-DY(4,1)      
      DO 35 I=2,5       
35    DC(6,I)=-DY(3,I)-DY(4,I)  
      DO 36 I=1,5       
      DC(7,I)=-(C(5)*DC(4,I)-C(4)*DC(5,I))/C(5)**2      
36    DY(5,I)=C(6)*DC(7,I)+C(7)*DC(6,I) 
      DO 37 I=1,5       
37    DY(6,I)=-DY(5,I)-DY(4,I)-DY(3,I)  
      DY(6,1)=DY(6,1)+1.D0      
      DO 38 I=1,5       
      DC(8,I)=0.995D0*(DY(5,I)+DY(4,I)) 
      DY(7,I)=(Y(1)*DC(8,I)-C(8)*DY(1,I))/Y(1)**2       
      DY(8,I)=DC(8,I)/3.798D+3  
      DC(9,I)=DY(7,I)-0.0663D0*(Y(8)*DY(7,I)-Y(7)*DY(8,I)       
     -)/Y(8)**2 
      DY(9,I)=-96.82D0*DC(9,I)/C(9)**2+0.321D0*DY(1,I)  
38    DY(10,I)=1.29D0*DY(5,I)+1.258D0*DY(4,I)+2.29D0*DY(3,I)+   
     -1.71D0*DY(6,I)    
      DY(11,1)=1.71D0-0.452D0*DY(4,1)+0.58D0*DY(3,1)    
      DO 39 I=2,5       
39    DY(11,I)=-0.452D0*DY(4,I)+0.58D0*DY(3,I)  
      DC(11,1)=1.74125D0*(Y(2)+X(1)*DY(2,1))    
      DO 40 I=2,5       
40    DC(11,I)=1.74125D0*X(1)*DY(2,I)   
      DO 41 I=1,5       
41    DC(12,I)=0.995D0*DY(10,I) 
      DY(12,1)=C(10)+X(1)*DC(10,1)+(C(12)*DC(11,1)-C(11 
     -)*DC(12,1))/C(12) 
     -**2       
      DO 42 I=2,5       
42    DY(12,I)=(C(12)*DC(11,I)-C(11)*DC(12,I))/C(12)**2 
      DO 43 I=1,5       
43    DY(13,I)=DC(12,I)-1.75D0*DY(2,I)  
      V7=-1.46312D+5/V4**2      
      DY(14,1)=V7*DY(9,1)       
      DY(14,2)=64.4D0+V7*DY(9,2)
      DY(14,3)=58.4D0+V7*DY(9,3)
      DY(14,4)=V7*DY(9,4)       
      DY(14,5)=V7*(1.D0+DY(9,5))
      DO 44 I=1,5       
44    DC(13,I)=0.995D0*DY(10,I)-0.1121D0*DY(14,I)       
      DC(13,2)=DC(13,2)+60.8D0  
      DC(13,4)=DC(13,4)+48.D0   
      DO 45 I=1,5       
      DY(15,I)=(C(13)*DY(13,I)-Y(13)*DC(13,I))/C(13)**2 
      DY(16,I)=-3.31D+5*DY(15,I)+40.D0*DY(13,I)-61.D0*(Y(15)    
     -*DY(13,I)+Y(13)*  
     -DY(15,I)) 
      DC(14,I)=2.324D+3*DY(10,I)-2.874D+7*DY(2,I)       
      DY(17,I)=-1.328D+3*DY(10,I)-531.D0*DY(11,I)+(C(12)*D      
     -C(14,I)-C(14)     
     -*DC(12,I))/C(12)**2       
      DC(15,I)=(Y(15)*DY(13,I)-Y(13)*DY(15,I))/Y(15)**2 
     --DY(13,I)/0.52D0  
45    DC(16,I)=-0.72D0*DY(15,I) 
      DO 46 I=1,4       
46    DC(17,I)=DY(9,I)  
      DC(17,5)=DY(9,5)+1.D0     
      GOTO (1,2,3,4,5),MODE     
2     FX = -(5.843D-7*Y(17) - 1.17D-4*Y(14) - 0.1365D0 
     /      - 2.358D-5*Y(13) - 1.502D-6*Y(16) - 0.0321D0*Y(12)
     /      - 4.324D-3*Y(5) - 1.D-4*C(15)/C(16) - 37.48D0*Y(2)/C(12))
      RETURN    
3     DO 47 I=1,5       
      GF(I)=-5.843D-7*DY(17,I)+1.17D-4*DY(14,I)+2.358D-5*DY(13,I)
     /   + 1.502D-6*DY(16,I)+0.0321D0*DY(12,I)+4.324D-3*DY(5,I)
     /   + 1.D-4*(C(16)*DC(15,I)-C(15)*DC(16,I))/C(16)**2
     /   + 37.48D0*(C(12)*DY(2,I)-Y(2)*DC(12,I))/C(12)**2       
47    CONTINUE  
      RETURN    
4     IF (INDEX1(1)) G(1)=1.5D0*X(2)-X(3)       
      IF (INDEX1(2)) G(2)=Y(1)-213.1D0  
      IF (INDEX1(3)) G(3)=405.23D0-Y(1) 
      DO 50 I=1,16      
      IF (INDEX1(I+3)) G(I+3)=Y(I+1)-A(I+1)     
      IF (INDEX1(I+19)) G(I+19)=B(I+1)-Y(I+1)   
50    CONTINUE  
      IF (INDEX1(36)) G(36)=Y(4)-0.28D0*Y(5)/0.72D0     
      IF (INDEX1(37)) G(37)=21.0D0-3.496D+3*Y(2)/C(12)  
      IF (INDEX1(38)) G(38)=6.2212D+4/C(17)-110.6D0-Y(1)
      RETURN    
5     DO 54 I=1,16      
      IF (.NOT.INDEX2(I+3)) GOTO 52     
      DO 51 J=1,5       
51    GG(I+3,J)=DY(I+1,J)       
52    IF (.NOT.INDEX2(I+19)) GOTO 54    
      DO 53 J=1,5       
53    GG(I+19,J)=-DY(I+1,J)     
54    CONTINUE  
      IF (.NOT.INDEX2(36)) GOTO 56      
      DO 55 J=1,5       
55    GG(36,J)=DY(4,J)-0.28D0*DY(5,J)/0.72D0    
56    IF (.NOT.INDEX2(37)) GOTO 58      
      DO 57 J=1,5       
57    GG(37,J)=-3.496D+3*(C(12)*DY(2,J)-Y(2)*DC(12,J))/C(12)**2   
58    IF (.NOT.INDEX2(38)) GOTO 60      
      DO 59 J=1,5       
59    GG(38,J)=-6.2212D+4*DC(17,J)/C(17)**2-DY(1,J)     
60    RETURN    
      END       
C
      SUBROUTINE TP86(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J,I1
      COMMON     /D86/ E,D,B,C,A,T,T1 
      DIMENSION E(5),D(5),B(10),C(5,5),A(10,5)  
      DOUBLEPRECISION E,D,B,C,A,T,T1     
      GOTO (1,2,3,4,5),MODE     
1     N=5       
      NILI=10   
      NINL=0    
      NELI=0    
      NENL=0    
      DO 26 I=1,4       
26    X(I)=0.D0 
      X(5)=1.D0 
      DO 6 I=1,5
      LXL(I)=.TRUE.     
      LXU(I)=.FALSE.    
6     XL(I)=0.D0
      E(1)=-15.D0       
      E(2)=-27.D0       
      E(3)=-36.D0       
      E(4)=-18.D0       
      E(5)=-12.D0       
      C(1,1)=30.D0      
      C(1,2)=-20.D0     
      C(1,3)=-10.D0     
      C(1,4)=32.D0      
      C(1,5)=-10.D0     
      C(2,2)=39.D0      
      C(2,3)=-6.D0      
      C(2,4)=-31.D0     
      C(2,5)=32.D0      
      C(3,3)=10.D0      
      C(3,4)=-6.D0      
      C(3,5)=-10.D0     
      C(4,4)=39.D0      
      C(4,5)=-20.D0     
      C(5,5)=30.D0      
      DO 27 I=1,4       
      I1=I+1    
      DO 27 J=I1,5      
27    C(J,I)=C(I,J)     
      D(1)=4.D0 
      D(2)=8.D0 
      D(3)=10.D0
      D(4)=6.D0 
      D(5)=2.D0 
      A(1,1)=-16.D0     
      A(1,2)=2.D0       
      A(1,3)=0.D0       
      A(1,4)=1.D0       
      A(1,5)=0.D0       
      A(2,1)=0.D0       
      A(2,2)=-2.D0      
      A(2,3)=0.D0       
      A(2,4)=0.4D0      
      A(2,5)=2.D0       
      A(3,1)=-3.5D0     
      A(3,2)=0.D0       
      A(3,3)=2.D0       
      A(3,4)=0.D0       
      A(3,5)=0.D0       
      A(4,1)=0.D0       
      A(4,2)=-2.D0      
      A(4,3)=0.D0       
      A(4,4)=-4.D0      
      A(4,5)=-1.D0      
      A(5,1)=0.D0       
      A(5,2)=-9.D0      
      A(5,3)=-2.D0      
      A(5,4)=1.D0       
      A(5,5)=-2.8D0     
      A(6,1)=2.D0       
      A(6,2)=0.D0       
      A(6,3)=-4.D0      
      A(6,4)=0.D0       
      A(6,5)=0.D0       
      A(8,1)=-1.D0      
      A(8,2)=-2.D0      
      A(8,3)=-3.D0      
      A(8,4)=-2.D0      
      A(8,5)=-1.D0      
      DO 29 I=1,5       
      A(7,I)=-1.D0      
      A(9,I)=I  
29    A(10,I)=1.D0      
      B(1)=-40.D0       
      B(2)=-2.D0
      B(3)=-0.25D0      
      B(4)=-4.D0
      B(5)=-4.D0
      B(6)=-1.D0
      B(7)=-40.D0       
      B(8)=-60.D0       
      B(9)=5.D0 
      B(10)=1.D0
      DO 25 I=1,10      
      DO 25 J=1,5       
25    GG(I,J)=A(I,J)    
      LEX=.FALSE.
      NEX=1       
      XEX(1)=0.299999999948D+00     
      XEX(2)=0.333467606492D+00     
      XEX(3)=0.400000000107D+00     
      XEX(4)=0.428310104740D+00     
      XEX(5)=0.223964873676D+00     
      FEX=-0.323486789716D+02 
      RETURN    
2     T=0.D0    
      DO 20 I=1,5       
      T1=0.D0   
      DO 21 J=1,5       
21    T1=T1+C(J,I)*X(I)*X(J)    
20    T=T+E(I)*X(I)+D(I)*X(I)**3+T1     
      FX=T      
      RETURN    
3     DO 23 I=1,5       
      T=0.D0    
      DO 22 J=1,5       
22    T=T+(C(I,J)+C(J,I))*X(J)  
23    GF(I)=E(I)+T+3.D0*D(I)*X(I)**2    
      RETURN    
4     DO 24 I=1,10      
24    IF (INDEX1(I)) G(I)=A(I,1)*X(1)+A(I,2)*X(2)+A(I,3 
     -)*X(3)+A(I,4)*    
     -X(4)+A(I,5)*X(5)-B(I)     
5     RETURN    
      END       
C      
      SUBROUTINE TP87(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION A,B,C,D,E,F1,F2,V1,V2,V3,V4,V5   
      A=131.078D0       
      B=1.48477D0       
      C=0.90798D0       
      D=DCOS(1.47588D0) 
      E=DSIN(1.47588D0) 
      GOTO (1,2,3,4,5),MODE     
1     N=6       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=4    
      X(1)=390.D0       
      X(2)=1.D+3
      X(3)=419.5D0      
      X(4)=340.5D0      
      X(5)=198.175D0    
      X(6)=0.5D0
      DO 6 I=1,6
      LXL(I)=.TRUE.     
6     LXU(I)=.TRUE.     
      XL(1)=0.D0
      XL(2)=0.D0
      XL(3)=340.D0      
      XL(4)=340.D0      
      XL(5)=-1.D+3      
      XL(6)=0.D0
      XU(1)=400.D0      
      XU(2)=1.D+3       
      XU(3)=420.D0      
      XU(4)=420.D0      
      XU(5)=1.D+3       
      XU(6)=0.5236D0    
      DO 70 I=3,6       
70    GF(I)=0.D0
      GG(1,1)=-1.D0     
      GG(1,2)=0.D0      
      GG(1,5)=0.D0      
      GG(2,1)=0.D0      
      GG(2,2)=-1.D0     
      GG(2,5)=0.D0      
      GG(3,1)=0.D0      
      GG(3,2)=0.D0      
      GG(3,5)=-1.D0     
      GG(4,1)=0.D0      
      GG(4,2)=0.D0      
      GG(4,5)=0.D0      
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.107811937779D+03     
      XEX(2)=0.196318606955D+03     
      XEX(3)=0.373830728516D+03     
      XEX(4)=0.420000000000D+03     
      XEX(5)=0.213071293896D+02     
      XEX(6)=0.153291953422D+00     
      FEX=0.892759773493D+04 
      RETURN    
2     IF (X(1)-300.D0) 31,32,32 
31    F1=30.D0*X(1)     
      GOTO 33   
32    F1=31.D0*X(1)     
33    IF (X(2)-100.D0) 34,35,35 
34    F2=28.D0*X(2)     
      GOTO 46   
35    IF (X(2)-200.D0) 36,37,37 
36    F2=29.D0*X(2)     
      GOTO 46   
37    F2=30.D0*X(2)     
46    FX=F1+F2  
      RETURN    
3     IF (X(1)-300.D0) 38,39,39 
38    GF(1)=30.D0       
      GOTO 40   
39    GF(1)=31.D0       
40    IF (X(2)-100.D0) 41,42,42 
41    GF(2)=28.D0       
      GOTO 45   
42    IF(X(2)-200.D0) 43,44,44  
43    GF(2)=29.D0       
      GOTO 45   
44    GF(2)=30.D0       
45    RETURN    
4     IF (INDEX1(1)) G(1)=-X(1)+300.D0-X(3)*X(4)/A*
     /     DCOS(B-X(6))+C*X(3)**2/A*D      
      IF (INDEX1(2)) G(2)=-X(2)-X(3)*X(4)/A*DCOS(B+X(6))
     /     +C*X(4)**2/A*D    
      IF (INDEX1(3)) G(3)=-X(5)-X(3)*X(4)/A*DSIN(B+X(6))
     /     +C*X(4)**2/A*E    
      IF (INDEX1(4)) G(4)=200.D0-X(3)*X(4)/A*DSIN(B-X(6))+      
     /     C*X(3)**2/A*E     
      RETURN    
5     V1=1.D0/A 
      IF (.NOT.(INDEX2(1).OR.INDEX2(4))) GOTO 8 
      V2=B-X(6) 
      V3=DCOS(V2)*V1    
      V4=DSIN(V2)*V1    
      V5=2.D0*C*X(3)*V1 
      IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,3)=-X(4)*V3+V5*D     
      GG(1,4)=-X(3)*V3  
      GG(1,6)=-X(3)*X(4)*V4     
7     IF (.NOT.INDEX2(4)) GOTO 8
      GG(4,3)=-X(4)*V4+V5*E     
      GG(4,4)=-X(3)*V4  
      GG(4,6)=X(3)*X(4)*V3      
8     IF (.NOT.(INDEX2(2).OR.INDEX2(3))) GOTO 10
      V2=B+X(6) 
      V3=DCOS(V2)*V1    
      V4=DSIN(V2)*V1    
      V5=2.D0*C*X(4)*V1 
      IF (.NOT.INDEX2(2)) GOTO 9
      GG(2,3)=-X(4)*V3  
      GG(2,4)=-X(3)*V3+V5*D     
      GG(2,6)=X(3)*X(4)*V4      
9     IF (.NOT.INDEX2(3)) GOTO 10       
      GG(3,3)=-X(4)*V4  
      GG(3,4)=-X(3)*V4+V5*E     
      GG(3,6)=-X(3)*X(4)*V3     
10    RETURN    
      END       
C
      SUBROUTINE TP88(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION GLEICH,MUE,A,DCOSKO,RHO,DV,T,DZ,INTKO,PI,Z,  
     /           V1,V2,V3,U,EP1,A1,W    
      COMMON     /D88/ MUE(30),A(30),DCOSKO(30),RHO(30),DV(6),T(6),
     /                 DZ(6),INTKO,PI,Z,V1,V2,V3,U,EP1,A1,W     
      INTEGER I,J,KN1,N1
      KN1=1     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.107431872940D+01     
      XEX(2)=-0.456613707247D+00     
      FEX=0.136265680997D+01 
      GOTO 7    
      ENTRY TP89(MODE)  
      KN1=2     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.107431872754D+01     
      XEX(2)=-0.456613706239D+00     
      XEX(3)=0.300836097604D-10     
      FEX=0.136265680508D+01 
      GOTO 7    
      ENTRY TP90(MODE)  
      KN1=3     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.708479399007D+00     
      XEX(2)=0.237919269592D-04     
      XEX(3)=0.807599939006D+00     
      XEX(4)=-0.456613723294D+00     
      FEX=0.136265681317D+01 
      GOTO 7    
      ENTRY TP91(MODE)  
      KN1=4     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.701892928031D+00     
      XEX(2)=0.221084326516D-11     
      XEX(3)=0.813330836201D+00     
      XEX(4)=0.456613707134D+00     
      XEX(5)=0.899937588382D-11     
      FEX=0.136265680910D+01 
      GOTO 7    
      ENTRY TP92(MODE)  
      KN1=5     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.494144465323D+00     
      XEX(2)=-0.103530473697D-04     
      XEX(3)=0.614950839550D+00     
      XEX(4)=-0.242186612731D-05     
      XEX(5)=0.729258528936D+00     
      XEX(6)=-0.456613099133D+00     
      FEX=0.136265681213D+01 
7     IF (MODE - 2) 1,17,17     
    1 N=KN1+1   
      NILI=0    
      NINL=1    
      NELI=0    
      NENL=0    
      DO 11 I=1,3       
      X(2*I-1)=0.5D0    
11    X(2*I)=-0.5D0     
      DO 6 I=1,6
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      XL(I)=-10.D0      
6     XU(I)=10.D0
      XL(1)=0.1       
      PI=DATAN(1.D0)*4.D0       
      DO 10 I=1,30      
      Z=PI*DBLE(I-1)  
      MUE(I)=GLEICH(Z)  
      V1=DSIN(MUE(I))   
      V2=DCOS(MUE(I))   
      DCOSKO(I)=(V1/MUE(I)-V2)/MUE(I)**2
10    A(I)=2.D0*V1/(MUE(I)+V1*V2)       
      INTKO=2.D0/15.D0  
      RETURN    
17    IF(MODE - 4)  19,18,18    
18    N1=N-1    
      T(N)=X(N)**2      
      DO 8 I=1,N1       
8     T(N-I)=T(N-I+1)+X(N-I)**2 
      V1=0.D0   
      DO 13 J=1,30      
      V2=MUE(J) 
      V3=-V2**2 
      RHO(J)=DBLE((-1)**N) 
      DO 14 I=1,N1      
      EP1=0.D0  
      A1=V3*T(N+1-I)    
      IF(A1.GT.-100.D0) EP1=DEXP(A1)    
14    RHO(J)=RHO(J)+DBLE((-1)**(N-I))*2.D0*EP1     
      EP1=0.D0  
      A1=V3*T(1)
      IF(A1.GT.-100.D0) EP1=DEXP(A1)    
      RHO(J)=(RHO(J)+EP1)/V3    
13    V1=V1-V3*A(J)*RHO(J)*(V2*DSIN(V2)*RHO(J)-2.D0*DCOSKO(J))       
19    GOTO (1,2,3,4,5),MODE     
2     U=0.D0    
      DO 20 I=1,N       
20    U=U+X(I)**2       
      FX=U      
      RETURN    
3     DO 21 I=1,N       
21    GF(I)=2.D0*X(I)   
      RETURN    
4     G(1)=1.D-4-V1-INTKO       
      RETURN    
5     DO 22 I=1,N       
22    DV(I)=0.D0
      DO 25 J=1,30      
      W=MUE(J)  
      V1=W**2*A(J)*(W*DSIN(W)*RHO(J)-DCOSKO(J)) 
      EP1=0.D0  
      A1=-MUE(J)**2*T(1)
      IF(A1.GT.-100.D0) EP1=DEXP(A1)    
      DZ(1)=EP1 
      DV(1)=DV(1)+DZ(1)*V1      
      DO 23 I=2,N       
      EP1=0.D0  
      A1=-MUE(J)**2*T(I)
      IF(A1.GT.-100.D0) EP1=DEXP(A1)    
      DZ(I)=DZ(I-1)+DBLE((-1)**(I+1))*2.D0*EP1     
23    DV(I)=DV(I)+DZ(I)*V1      
25    CONTINUE  
      DO 24 I=1,N       
24    GG(1,I)=-4.D0*DV(I)*X(I)  
      RETURN    
      END       
C      
      SUBROUTINE TP93(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION V1,V2,V3,V4,V5,V6,V7,V8,V9 
      INTEGER I
      GOTO (1,2,3,4,5),MODE     
1     N=6       
      NILI=0    
      NINL=2    
      NELI=0    
      NENL=0    
      X(1)=5.54D0       
      X(2)=4.4D0
      X(3)=12.02D0      
      X(4)=11.82D0      
      X(5)=0.702D0      
      X(6)=0.852D0      
      DO 6 I=1,6
      LXL(I)=.TRUE.     
      LXU(I)=.FALSE.    
6     XL(I)=0.D0
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.533266639884D+01     
      XEX(2)=0.465674439073D+01     
      XEX(3)=0.104329901123D+02     
      XEX(4)=0.120823085893D+02     
      XEX(5)=0.752607369745D+00     
      XEX(6)=0.878650836850D+00     
      FEX=0.135075961229D+03 
      RETURN    
2     V1=X(1)+X(2)+X(3) 
      V2=X(1)+1.57D0*X(2)+X(4)  
      V3=X(1)*X(4)      
      V4=X(3)*X(2)      
      FX=0.0204D0*V3*V1+0.0187D0*V4*V2+0.0607D0*V3*V1*X(5)**2   
     -+0.0437D0*V4*V2*X(6)**2    
      RETURN    
3     V1=X(1)*X(4)      
      V2=X(2)*X(3)      
      V3=X(2)+X(3)      
      V4=X(1)+X(4)      
      V5=V3+X(1)
      V6=1.57D0*X(2)+V4 
      V7=X(4)*X(5)**2   
      V8=X(3)*X(6)**2   
      V9=0.0607D0*X(1)*V7       
      GF(1)=0.0408D0*V1+0.0204D0*X(4)*V3+0.0187D0*V2+2.D0*V9+   
     -0.0607D0*V7*V3    
     -+0.0437D0*X(2)*V8 
      GF(2)=0.0204D0*V1+0.058718D0*V2+0.0187D0*X(3)*V4+V9+      
     -0.137218D0*X(2)   
     -*V8+0.0437D0*V8*V4
      GF(3)=0.0204D0*V1+0.0187D0*X(2)*V6+V9+0.0437D0*X(2)       
     -*X(6)**2*V6       
      GF(4)=0.0204D0*X(1)*V5+0.0187D0*V2+0.0437D0*X(2)*V8+      
     -.0607D0*X(1)*X(5)**2      
     -*V5       
      GF(5)=0.1214D0*V1*X(5)*V5 
      GF(6)=0.0874D0*X(6)*V2*V6 
      RETURN    
4     IF (INDEX1(1)) G(1)=1.D-3*X(1)*X(2)*X(3)*X(4)*X(5)*X(6)-2.07D0 
      IF (INDEX1(2)) G(2)=1.D0-6.2D-4*X(1)*X(4)*X(5)**2*(X(1)+X(2)+X(3))
     /    -5.8D-4*X(2)*X(3)*X(6)**2*(X(1)+1.57D0*X(2)+X(4))  
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      V1=X(1)*X(2)*X(3)*1.D-3   
      V2=X(4)*X(5)*X(6)*1.D-3   
      GG(1,1)=X(2)*X(3)*V2      
      GG(1,2)=X(1)*X(3)*V2      
      GG(1,3)=X(1)*X(2)*V2      
      GG(1,4)=X(5)*X(6)*V1      
      GG(1,5)=X(4)*X(6)*V1      
      GG(1,6)=X(4)*X(5)*V1      
7     IF (.NOT.INDEX2(2)) GOTO 8
      V1=-X(5)**2*6.2D-4
      V2=-X(6)**2*5.8D-4
      V3=X(1)+X(2)+X(3) 
      V4=X(1)+1.57D0*X(2)+X(4)  
      V5=V1*V3  
      V6=V2*V4  
      V7=V1*X(1)*X(4)   
      V8=V2*X(2)*X(3)   
      GG(2,1)=V7+V5*X(4)+V8     
      GG(2,2)=V7+V6*X(3)+1.57D0*V8      
      GG(2,3)=V7+V6*X(2)
      GG(2,4)=V5*X(1)+V8
      GG(2,5)=-1.24D-3*X(1)*X(4)*X(5)*V3
      GG(2,6)=-1.16D-3*X(2)*X(3)*X(6)*V4
8     RETURN    
      END       
C
      SUBROUTINE TP95(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D95/ B(4)    
      DOUBLEPRECISION B  
      INTEGER I,KN1
      KN1=1     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=-0.476149332788D-11     
      XEX(2)=-0.355239427962D-10     
      XEX(3)=-0.702611041315D-10     
      XEX(4)=-0.171856469485D-10     
      XEX(5)=-0.748993551642D-10     
      XEX(6)=0.332330328254D-02     
      FEX=0.156195144282D-01 
      GOTO 11   
      ENTRY TP96(MODE)  
      KN1=2     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=-0.519722825686D-11     
      XEX(2)=-0.387748184662D-10     
      XEX(3)=-0.766908552858D-10     
      XEX(4)=-0.187583442974D-10     
      XEX(5)=-0.817535626869D-10     
      XEX(6)=0.332330328612D-02     
      FEX=0.156195134384D-01 
      GOTO 11   
      ENTRY TP97(MODE)  
      KN1=3     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.268564912352D+00     
      XEX(2)=0.0D0  
      XEX(3)=0.0D0  
      XEX(4)=0.0D0  
      XEX(5)=0.28D-01     
      XEX(6)=0.134D-01     
C      FEX=0.313580912311D+01 
      FEX=0.31358089D+1 
      GOTO 11   
      ENTRY TP98(MODE)  
      KN1=4     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.268564912323D+00     
      XEX(2)=0.0D0  
      XEX(3)=0.0D0  
      XEX(4)=0.0D0  
      XEX(5)=0.280000000000D-01     
      XEX(6)=0.134000000001D-01     
C      FEX=0.313580912299D+01 
      FEX=0.31358089D+1 
11    GOTO (1,2,3,4,5),MODE     
    1 N=6       
      NILI=0    
      NINL=4    
      NELI=0    
      NENL=0    
      DO 6 I=1,6
      X(I)=0.D0 
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
6     XL(I)=0.D0
      XU(1)=0.31D0      
      XU(2)=0.046D0     
      XU(3)=0.068D0     
      XU(4)=0.042D0     
      XU(5)=0.028D0     
      XU(6)=0.0134D0    
      GF(1)=4.3D0       
      GF(2)=31.8D0      
      GF(3)=63.3D0      
      GF(4)=15.8D0      
      GF(5)=68.5D0      
      GF(6)=4.7D0       
      GG(1,2)=38.2D0    
      GG(2,2)=36.8D0    
      GG(3,1)=0.D0      
      GG(3,2)=-273.D0   
      GG(3,3)=0.D0      
      GG(3,6)=0.D0      
      GG(4,2)=-311.D0   
      GG(4,3)=0.D0      
      GG(4,4)=587.D0    
      GG(4,5)=391.D0    
      IF (KN1-2) 20,20,21       
20    B(1)=4.97D0       
      B(2)=-1.88D0      
      GOTO 22   
21    B(1)=32.97D0      
      B(2)=25.12D0      
22    IF (KN1-3) 25,24,23       
23    B(3)=-124.08D0    
      B(4)=-173.02D0    
      GOTO 27   
24    B(3)=-29.08D0     
      B(4)=-78.02D0     
      GOTO 27   
25    IF (KN1-2) 24,26,26       
26    B(3)=-69.08D0     
      B(4)=-118.02D0    
27    CONTINUE  
      RETURN    
2     FX=4.3D0*X(1)+31.8D0*X(2)+63.3D0*X(3)+15.8D0*X(4)
     /   +68.5D0*X(5)+4.7D0*X(6)     
3     RETURN    
4     IF (INDEX1(1)) G(1)=17.1*X(1)+38.2*X(2)+204.2*X(3)
     /     +212.3*X(4)+623.4*X(5)+1.4955D+3*X(6)-169.0*X(1)*X(3)
     /     -3.58D+3*X(3)*X(5)-3.81D+3*X(4)*X(5)-1.85D+4*X(4)*X(6)
     /     -2.43D+4*X(5)*X(6)-B(1)    
      IF (INDEX1(2)) G(2)=17.9*X(1)+36.8*X(2)+113.9*X(3)
     /     +169.7*X(4)+337.8*X(5)+1.3852D+3*X(6)-139.0*X(1)*X(3)
     /     -2.45D+3*X(4)*X(5)-1.66D+4*X(4)*X(6)-1.72D+4*X(5)*X(6)-B(2) 
      IF (INDEX1(3)) G(3)=-273.0*X(2)-70.0*X(4)-819.0*X(5)
     /     +2.6D+4*X(4)*X(5)-B(3)      
      IF (INDEX1(4)) G(4)=159.9D0*X(1)-311.D0*X(2)+587.D0*X(4)
     /     +391.0*X(5)+2.198D+3*X(6)-1.4D+4*X(1)*X(6)-B(4)      
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=17.1-169.0*X(3)
      GG(1,3)=204.2-169.0*X(1)-3.58D+3*X(5)  
      GG(1,4)=212.3-3.81D+3*X(5)-1.85D+4*X(6) 
      GG(1,5)=623.4-3.58D+3*X(3)-3.81D+3*X(4)-2.43D+4*X(6)    
      GG(1,6)=1.4955D+3-1.85D+4*X(4)-2.43D+4*X(5)       
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,1)=17.90-139.0*X(3)
      GG(2,3)=113.90-139.0*X(1)       
      GG(2,4)=169.70-2.45D+3*X(5)-1.66D+4*X(6) 
      GG(2,5)=337.80-2.45D+3*X(4)-1.72D+4*X(6) 
      GG(2,6)=1.3852D+3-1.66D+4*X(4)-1.72D+4*X(5)       
8     IF(.NOT.INDEX2(3)) GOTO 9 
      GG(3,4)=-70.D0+2.6D+4*X(5)
      GG(3,5)=-819.D0+2.6D+4*X(4)       
9     IF (.NOT.INDEX2(4)) GOTO 10       
      GG(4,1)=159.9-1.4D+4*X(6)       
      GG(4,6)=2.198D+3-1.4D+4*X(1)      
10    RETURN    
      END       
C
      SUBROUTINE TP99(MODE)     
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D99/ R,S,P,Q,DP,DQ,DR,DS,A,T,V1,V2,V3,V4 
      DIMENSION P(8),Q(8),R(8),S(8),DP(8,7),DQ(8,7),DR(8,7),DS(8,7),
     /          A(8),T(8)
      DOUBLEPRECISION R,S,P,Q,DP,DQ,DR,DS,A,T,V1,V2,V3,V4      
      INTEGER I,J,I1
      IF (MODE - 2)1,18,17      
1     N=7       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=2    
      DO 6 I=1,7
      X(I)=0.5D0
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      XL(I)=0.D0
6     XU(I)=1.58D0      
      A(1)=0.D0 
      A(2)=50.D0
      A(3)=50.D0
      A(4)=75.D0
      A(5)=75.D0
      A(6)=75.D0
      A(7)=100.D0       
      A(8)=100.D0       
      T(1)=0.D0 
      T(2)=25.D0
      T(3)=50.D0
      T(4)=100.D0       
      T(5)=150.D0       
      T(6)=200.D0       
      T(7)=290.D0       
      T(8)=380.D0       
      P(1)=0.D0 
      Q(1)=0.D0 
      R(1)=0.D0 
      S(1)=0.D0 
      DO 31 J=1,7       
      DO 31 I=1,J       
      DP(I,J)=0.D0      
      DQ(I,J)=0.D0      
      DR(I,J)=0.D0      
31    DS(I,J)=0.D0      
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.542460319142D+00     
      XEX(2)=0.529015870015D+00     
      XEX(3)=0.508450583169D+00     
      XEX(4)=0.480269265187D+00     
      XEX(5)=0.451235157238D+00     
      XEX(6)=0.409187805755D+00     
      XEX(7)=0.352784693565D+00     
      FEX=-0.831079891516D+09
      RETURN    
17    IF (MODE - 4) 18,18,19    
18    DO 30 I=2,8       
      I1=I-1    
      V1=A(I)*DSIN(X(I1))-32.D0 
      V2=A(I)*DCOS(X(I1))       
      V3=T(I)-T(I1)     
      V4=0.5D0*V3**2    
      P(I)=V2*V4+V3*R(I1)+P(I1) 
      Q(I)=V1*V4+V3*S(I1)+Q(I1) 
      R(I)=V2*V3+R(I1)  
30    S(I)=V1*V3+S(I1)  
      IF (MODE - 3) 40,19,40    
19    DO 34 I=2,8       
      DO 34 J=1,7       
      IF (J-I+1) 33,32,34       
32    I1=I-1    
      V1=A(I)*DSIN(X(I1))       
      V2=A(I)*DCOS(X(I1))       
      V3=T(I)-T(I1)     
      V4=0.5D0*V3**2    
      DP(I,I1)=-V1*V4+V3*DR(I1,I1)+DP(I1,I1)    
      DQ(I,I1)=V2*V4+V3*DS(I1,I1)+DQ(I1,I1)     
      DR(I,I1)=-V1*V3+DR(I1,I1) 
      DS(I,I1)=V2*V3+DS(I1,I1)  
      GOTO 34   
33    I1=I-1    
      V1=T(I)-T(I1)     
      DP(I,J)=V1*DR(I1,J)+DP(I1,J)      
      DQ(I,J)=V1*DS(I1,J)+DQ(I1,J)      
      DR(I,J)=DR(I1,J)  
      DS(I,J)=DS(I1,J)  
34    CONTINUE  
40    GOTO (1,2,3,4,5),MODE     
2     FX=-R(8)**2
      RETURN    
3     DO 35 I=1,7       
35    GF(I)=-2.D0*R(8)*DR(8,I)
      RETURN    
4     IF (INDEX1(1)) G(1)=Q(8)-1.D+5    
      IF (INDEX1(2)) G(2)=S(8)-1.D+3    
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      DO 36 I=1,7       
36    GG(1,I)=DQ(8,I)   
7     IF (.NOT.INDEX2(2)) GOTO 8
      DO 37 I=1,7       
37    GG(2,I)=DS(8,I)   
8     RETURN    
      END       
C
      SUBROUTINE TP100(MODE)    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION V1,V2      
      GOTO (1,2,3,4,5),MODE     
1     N=7       
      NILI=0    
      NINL=4    
      NELI=0    
      NENL=0    
      DO 6 I=1,7
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      X(1)=1.D0 
      X(2)=2.D0 
      X(3)=0.D0 
      X(4)=4.D0 
      X(5)=0.D0 
      X(6)=1.D0 
      X(7)=1.D0 
      GG(1,3)=-1.D0     
      GG(1,5)=-5.D0     
      GG(1,6)=0.D0      
      GG(1,7)=0.D0      
      GG(2,1)=-7.D0     
      GG(2,2)=-3.D0     
      GG(2,4)=-1.D0     
      GG(2,5)=1.D0      
      GG(2,6)=0.D0      
      GG(2,7)=0.D0      
      GG(3,1)=-23.D0    
      GG(3,3)=0.D0      
      GG(3,4)=0.D0      
      GG(3,5)=0.D0      
      GG(3,7)=8.D0      
      GG(4,4)=0.D0      
      GG(4,5)=0.D0      
      GG(4,6)=-5.D0     
      GG(4,7)=11.D0     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.233049937431D+01     
      XEX(2)=0.195137237315D+01     
      XEX(3)=-0.477541392625D+00     
      XEX(4)=0.436572623462D+01     
      XEX(5)=-0.624486970475D+00     
      XEX(6)=0.103813101881D+01     
      XEX(7)=0.159422671137D+01     
      FEX=0.680630057275D+03
      RETURN    
2     FX=(X(1)-10.D0)**2+5.D0*(X(2)-12.D0)**2+X(3)**4
     /    +3.D0*(X(4)-11.D0)**2+10.D0*X(5)**6+7.D0*X(6)**2
     /    +X(7)**4-4.D0*X(6)*X(7)-10.D0*X(6)-8.D0*X(7)      
      RETURN    
3     GF(1)=2.D0*(X(1)-10.D0)   
      GF(2)=10.D0*(X(2)-12.D0)  
      GF(3)=4.D0*X(3)**3
      GF(4)=6.D0*(X(4)-11.D0)   
      GF(5)=60.D0*X(5)**5       
      GF(6)=14.D0*X(6)-4.D0*X(7)-10.D0  
      GF(7)=4.D0*X(7)**3-4.D0*X(6)-8.D0 
      RETURN    
4     V1=2.D0*X(1)**2   
      V2=X(2)**2
      IF (INDEX1(1)) G(1)=-V1-3.D0*V2**2-X(3)-4.D0*X(4)**2-     
     /                    5.D0*X(5)+127.D0  
      IF (INDEX1(2)) G(2)=-7.D0*X(1)-3.D0*X(2)-10.D0*X(3)**2-   
     /                    X(4)+X(5)+282.D0  
      IF (INDEX1(3)) G(3)=-23.D0*X(1)-V2-6.D0*X(6)**2
     /                    +8.D0*X(7)+196.D0 
      IF (INDEX1(4)) G(4)=-2.D0*V1-V2+3.D0*X(1)*X(2)
     /                     -2.D0*X(3)**2-5.D0*X(6)+11.D0*X(7)  
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=-4.D0*X(1)
      GG(1,2)=-12.D0*X(2)**3    
      GG(1,4)=-8.D0*X(4)
7     IF (INDEX2(2)) GG(2,3)=-20.D0*X(3)
      IF(.NOT.INDEX2(3)) GOTO 9 
      GG(3,2)=-2.D0*X(2)
      GG(3,6)=-12.D0*X(6)       
9     IF (.NOT.INDEX2(4)) GOTO 10       
      GG(4,1)=-8.D0*X(1)+3.D0*X(2)      
      GG(4,2)=-2.D0*X(2)+3.D0*X(1)      
      GG(4,3)=-4.D0*X(3)
10    RETURN    
      END       
C
      SUBROUTINE TP101(MODE)    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D101/ A(3),GV(7),FMIN(3)      
      DOUBLEPRECISION A,GV,FMIN,SUM,V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,
     /     V11,V12,V13,V14,V15,V16,V17,V18,V19,V20,V21,V22,V23,V24,V25, 
     /     V26,V27,V28,V29,V30,V31,V32,V33,V34,V35,V36,V37,V38, 
     /     V39,V40,V41,V42,V43,V44,V45
      INTEGER I,KN1,K,M
      KN1=1     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.285615855584D+01     
      XEX(2)=0.610823030755D+00     
      XEX(3)=0.215081256203D+01     
      XEX(4)=0.471287370945D+01     
      XEX(5)=0.999487540961D+00     
      XEX(6)=0.134750750498D+01     
      XEX(7)=0.316527664991D-01     
      FEX=0.180976476556D+04 
      GOTO 13   
      ENTRY TP102(MODE) 
      KN1=2     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.389625319099D+01     
      XEX(2)=0.809358760118D+00     
      XEX(3)=0.266438599373D+01     
      XEX(4)=0.430091287458D+01     
      XEX(5)=0.853554935267D+00     
      XEX(6)=0.109528744459D+01     
      XEX(7)=0.273104596581D-01     
      FEX=0.911880571336D+03 
      GOTO 13   
      ENTRY TP103(MODE) 
      KN1=3     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.439410451026D+01     
      XEX(2)=0.854468738817D+00     
      XEX(3)=0.284323031380D+01     
      XEX(4)=0.339997866779D+01     
      XEX(5)=0.722926133025D+00     
      XEX(6)=0.870406381840D+00     
      XEX(7)=0.246388263302D-01     
      FEX=0.543667958424D+03 
   13 CONTINUE 
      M=KN1      
      GOTO (1,2,3,4,5),MODE     
    1 N=7       
      A(1)=-0.25D0      
      A(2)=0.125D0      
      A(3)=0.5D0
      NILI=0    
      NINL=6    
      NELI=0    
      NENL=0    
      DO 6 I=1,7
      X(I)=6.D0 
      XL(I)=0.1D0       
      XU(I)=10.D0       
      LXL(I)=.TRUE.     
6     LXU(I)=.TRUE.     
      XL(7)=0.01D0      
      GG(1,5)=0.D0      
      GG(2,7)=0.D0      
      GG(3,4)=0.D0      
      GG(4,6)=0.D0      
      RETURN    
2     DO 200 K=1,7      
200   IF (X(K).LT..0D-8) GOTO 14 
      FX=10.D0*X(1)*X(4)**2*X(7)**A(M)/(X(2)*X(6)**3)
     /    +15.D0*X(3)*X(4)/(X(1)*X(2)**2*X(5)*X(7)**0.5D0)
     /    +20.D0*X(2)*X(6)/(X(1)**2*X(4)*X(5)**2)
     /    +25.D0*X(1)**2*X(2)**2*X(5)**0.5D0*X(7)/(X(3)*X(6)**2)    
      RETURN    
14    SUM=0.D0  
      DO 40 I=1,7       
40    SUM=SUM+(X(I)-5.D0)**2    
      FMIN(1)=1.8D+3    
      FMIN(2)=9.1D+2    
      FMIN(3)=5.4D+2    
      FX=SUM+1.D+3+FMIN(KN1)    
      RETURN    
3     DO 201 K=1,7      
201   IF (X(K).LT.1.0D-8) GOTO 15 
      V1=10.D0*X(4)**2  
      V2=X(7)**A(M)     
      V3=X(2)*X(6)**3   
      V4=15.D0*X(3)*X(4)
      V5=X(1)*X(2)**2*X(5)*X(7)**0.5D0  
      V6=20.D0*X(2)*X(6)
      V7=X(1)**2*X(4)*X(5)**2   
      V8=25.D0*X(1)*X(2)*X(5)**0.5D0*X(7)       
      V9=X(3)*X(6)**2   
      V10=12.5D0*X(1)**2*X(2)**2*X(7)   
      V11=X(5)**0.5D0   
      GF(1)=V1*V2/V3-V4/(X(1)*V5)-2.D0*V6/(X(1)*V7)+2.D0*X(2)*V8/V9  
      GF(2)=-V1*X(1)*V2/(X(2)*V3)-2.D0*V4/(X(2)*V5)+20.D0*X(6)/V7
     /      +2.D0*X(1)*V8/V9       
      GF(3)=15.D0*X(4)/V5-X(1)*X(2)*V8/(X(3)*V9)
      GF(4)=20.D0*X(1)*X(4)*V2/V3+15.D0*X(3)/V5-V6/(X(4)*V7) 
      GF(5)=-V4/(X(5)*V5)-2.D0*V6/(X(5)*V7)+V10/(V9*V11)
      GF(6)=-3.D0*V1*X(1)*V2/(X(6)*V3)+20.D0*X(2)/V7-4.D0*V10   
     /        *V11/(X(6)*V9)    
      GF(7)=A(M)*V1*X(1)*X(7)**(A(M)-1.D0)/V3-0.5D0*V4/(V5*     
     /        X(7))+V8*X(1)*X(2)/(X(7)*V9)       
      RETURN    
15    DO 50 I=1,7       
50    GF(I)=2.D0*(X(I)-5.D0)    
      RETURN    
4     DO K=1,7      
         IF (X(K).LT.1.0D-8) GOTO 16 
      ENDDO   
      IF (INDEX1(1)) G(1)=1.D0-0.5D0*DABS(X(1))**0.5D0*X(7)/
     /    (X(3)*X(6)**2)-0.7D0
     /     *X(1)**3*X(2)*X(6)*DABS(X(7))**0.5D0/X(3)**2
     /     -0.2D0*X(3)*DABS(X(6))**(2.D0/3.D0)*DABS(X(7))**0.25D0
     /     /(X(2)*DABS(X(4))**0.5D0) 
      IF (INDEX1(2)) G(2)=1.D0-1.3D0*X(2)*X(6)/(DABS(X(1))**    
     /     0.5D0*X(3)*X(5))  
     /     -0.8D0*X(3)*X(6)**2/(X(4)*X(5))-3.1D0*DABS(X(2))**0.5D0*  
     /      DABS(X(6))**(1.D0/3.D0)/(X(1)*X(4)**2*X(5))     
      IF (INDEX1(3)) G(3)=1.D0-2.D0*X(1)*X(5)*DABS(X(7))**(1.D0/3.D0)
     /   /(DABS(X(3))**1.5D0*X(6))-0.1D0*X(2)*X(5)/(DABS(X(3)*X(7))
     /   **0.5D0*X(6))-X(2)*DABS(X(3))**0.5D0*X(5)/X(1)
     /    -0.65D0*X(3)*X(5)*X(7)/(X(2)**2*X(6))       
      IF (INDEX1(4)) G(4)=1.D0-0.2D0*X(2)*DABS(X(5))**0.5D0
     /     *DABS(X(7))**(1.D0/3.D0)/(X(1)**2*X(4))
     /    -0.3D0*DABS(X(1))**0.5D0*X(2)**2*X(3)*DABS(X(4))
     /     **(1.D0/3.D0)*DABS(X(7))**0.25D0/DABS(X(5))**(2.D0/3.D0)
     /     -0.4D0*X(3)*X(5)*DABS(X(7))**0.75D0/(X(1)**3*X(2)**2)
     /     -0.5D0*X(4)*DABS(X(7))**0.5D0/X(3)**2
      IF (INDEX1(5)) G(5)=10.D0*X(1)*X(4)**2*DABS(X(7))**A(M)
     /    /(X(2)*X(6)**3)+15.D0*X(3)*X(4)/(X(1)*X(2)**2*X(5)
     /    *DABS(X(7))**0.5D0)+20.D0*X(2)*X(6)/
     /    (X(1)**2*X(4)*X(5)**2)+25.D0*X(1)**2*X(2)**2
     /    *DABS(X(5))**0.5D0*X(7)/(X(3)*X(6)**2)-100.D0       
      IF (INDEX1(6)) G(6)=-(10.D0*X(1)*X(4)**2*DABS(X(7))**A(M)
     /     /(X(2)*X(6)**3)+15.D0*X(3)*X(4)/(X(1)*X(2)**2
     /     *X(5)*DABS(X(7))**0.5D0)+20.D0*X(2)*X(6)     
     /     /(X(1)**2*X(4)*X(5)**2)+25.D0*X(1)**2*X(2)**2
     /     *DABS(X(5))**0.5D0*X(7)/(X(3)*X(6)**2))+3.D+3       
      RETURN    
16    DO I=1,6
         G(I)=0.0D0
      ENDDO
      RETURN         
5     IF(.NOT.INDEX2(1)) GOTO 7 
      V1=DABS(X(1))**0.5D0      
      V2=X(1)**3
      V4=X(3)**2
      V5=V4**2  
      V6=DABS(X(4))**0.5D0      
      V7=X(6)**2
      V8=DABS(X(6))**(2.D0/3.D0)
      V9=DABS(X(7))**0.5D0      
      V10=DABS(X(7))**0.25D0    
      GG(1,1)=-0.25D0*X(7)/(V1*X(3)*V7)-2.1D0*X(1)**2*X(2)*X(6)*V9/V4
      GG(1,2)=-0.7D0*V2*X(6)*V9/V4+0.2D0*X(3)*V8*V10/(X(2)**2*V6)    
      GG(1,3)=0.5D0*V1*X(7)/(V4*V7)+1.4D0*V2*X(2)*X(6)*V9/(X(3)*V4)
     /      -0.2D0*V8*V10/(X(2)*V6)      
      GG(1,4)=0.1D0*X(3)*V8*V10/(X(2)*X(4)*V6)  
      GG(1,6)=V1*X(7)/(X(3)*V7*X(6))-0.7D0*V2*X(2)*V9/V4-0.4D0
     /     /3.D0*X(3)*V10/(X(2)*V6*DABS(X(6))**(1.D0/3.D0))
      GG(1,7)=-0.5D0*V1/(X(3)*V7)-0.35D0*V2*X(2)*X(6)/(V4*V9)
     /       -0.05D0*X(3)*V8/(X(2)*V6*V9*V10)      
7     IF (.NOT.INDEX2(2)) GOTO 8
      V11=DABS(X(1))**0.5D0     
      V12=DABS(X(2))**0.5D0     
      V13=X(4)**2       
      V14=X(5)**2       
      V15=DABS(X(6))**(1.D0/3.D0)       
      V16=X(6)**2       
      GG(2,1)=0.65D0*X(2)*X(6)/(X(1)*V11*X(3)*X(5))+3.1D0*V12*V15/
     /     (X(1)**2*V13*X(5))       
      GG(2,2)=-1.3D0*X(6)/(V11*X(3)*X(5))-1.55D0*V15/(X(1)*     
     /     V12*V13*X(5))     
      GG(2,3)=1.3D0*X(2)*X(6)/(V11*X(3)**2*X(5))-0.8D0*V16/(X(4)*X(5))  
      GG(2,4)=0.8D0*X(3)*V16/(V13*X(5))+6.2D0*V12*V15/(X(1)     
     /     *V13*X(4)*X(5))   
      GG(2,5)=1.3D0*X(2)*X(6)/(V11*X(3)*V14)+0.8D0*X(3)*V16     
     /           /(X(4)*V14)+3.1D0*V12*V15/(X(1)*V13*V14)   
      GG(2,6)=-1.3D0*X(2)/(V11*X(3)*X(5))-1.6D0*X(3)*X(6)/(X(4)*X(5))-  
     /        3.1D0/3.D0*V12/(X(1)*V13*X(5)*V15**2)     
8     IF(.NOT.INDEX2(3)) GOTO 9 
      V17=X(2)**2       
      V18=DABS(X(3))**0.5D0     
      V19=V18*X(3)      
      V20=X(6)**2       
      V21=DABS(X(7))**(1.D0/3.D0)       
      V22=DABS(X(7))**0.5D0     
      GG(3,1)=-2.D0*X(5)*V21/(V19*X(6))+X(2)*V18*X(5)/X(1)**2      
      GG(3,2)=-V18*X(5)/X(1)+1.3D0*X(3)*X(5)*X(7)/(V17*X(2)*X(6))
     /       -0.1D0*X(5)/(V18*V22*X(6))       
      GG(3,3)=3.D0*X(1)*X(5)*V21/(X(3)*V19*X(6))+0.05D0*X(2)*X(5)
     /     /(V19*X(6)*V22)-0.5D0*X(2)*X(5)/(X(1)*V18)
     /     -0.65D0*X(5)*X(7)/(V17*X(6))    
      GG(3,5)=-2.D0*X(1)*V21/(V19*X(6))-0.1D0*X(2)/(V18*X(6)*V22)
     /     -X(2)*V18/X(1)-0.65D0*X(3)*X(7)/(V17*X(6))  
      GG(3,6)=2.D0*X(1)*X(5)*V21/(V19*V20)+0.1D0*X(2)*X(5)/     
     /      (V18*V20*V22)+0.65D0*X(3)*X(5)*X(7)/(V17*V20)   
      GG(3,7)=-2.D0/3.D0*X(1)*X(5)/(V19*X(6)*V21**2)+0.05D0*X(2)*X(5)
     /     /(V18*X(6)*V22*X(7))-0.65D0*X(3)*X(5)/(V17*X(6))   
9     IF (.NOT.INDEX2(4)) GOTO 10       
      V23=DABS(X(1))**0.5D0     
      V24=X(1)**2       
      V25=V24*X(1)      
      V26=X(2)**2       
      V27=X(3)**2       
      V28=DABS(X(4))**(1.D0/3.D0)       
      V29=DABS(X(5))**(2.D0/3.D0)       
      V30=DABS(X(5))**0.5D0     
      V31=DABS(X(7))**0.25D0    
      V32=V31**2
      V33=V31*V32       
      V34=DABS(X(7))**(1.D0/3.D0)       
      GG(4,1)=0.4D0*X(2)*V30*V34/(V25*X(4))-0.15D0*V26*X(3)*V28*V31/ 
     /     (V23*V29)+1.2D0*X(3)*X(5)*V33/(V24**2*V26)
      GG(4,2)=-0.2D0*V30*V34/(V24*X(4))-0.6D0*V23*X(2)*X(3)*V28*V31/V29
     /     +0.8D0*X(3)*X(5)*V33/(V25*V26*X(2))       
      GG(4,3)=-0.3D0*V23*V26*V28*V31/V29-0.4D0*X(5)*V33/(V25*V26)
     /      +X(4)*V32/(V27*X(3))     
      GG(4,4)=0.2D0*X(2)*V30*V34/(V24*X(4)**2)-0.1D0*V23*V26*X(3)*V31
     /      /(V28**2*V29)-0.5D0*V32/V27   
      GG(4,5)=-0.1D0*X(2)*V34/(V24*X(4)*V30)+0.2D0*V23*V26*     
     /       X(3)*V28*V31/(X(5)*V29)-0.4D0*X(3)*V33/(V25*V26)  
      GG(4,7)=-0.2D0/3.D0*X(2)*V30/(V24*X(4)*V34**2)-0.075D0*   
     /       V23*V26*X(3)*V28/(V29*V33)-0.3D0*X(3)*X(5)/(V25*V26*V31)
     /       -0.25D0*X(4)/(V27*V32) 
10    IF (.NOT.INDEX2(5).AND..NOT.INDEX2(6)) GOTO 12    
      V35=10.D0*X(4)**2 
      V36=DABS(X(7))**A(M)      
      V37=X(2)*X(6)**3  
      V38=15.D0*X(3)*X(4)       
      V39=X(1)*X(2)**2*X(5)*DABS(X(7))**0.5D0   
      V40=20.D0*X(2)*X(6)       
      V41=X(1)**2*X(4)*X(5)**2  
      V42=25.D0*X(1)*X(2)*DABS(X(5))**0.5D0*X(7)
      V43=X(3)*X(6)**2  
      V44=12.5D0*X(1)**2*X(2)**2*X(7)   
      V45=DABS(X(5))**0.5D0     
      GV(1)=V35*V36/V37-V38/(X(1)*V39)-2.D0*V40/(X(1)*V41)+2.D0*X(2) 
     /     *V42/V43  
      GV(2)=-V35*X(1)*V36/(X(2)*V37)-2.D0*V38/(X(2)*V39)+20.D0*X(6)     
     /    /V41+2.D0*X(1)*V42/V43    
      GV(3)=15.D0*X(4)/V39-X(1)*X(2)*V42/(X(3)*V43)     
      GV(4)=20.D0*X(1)*X(4)*V36/V37+15.D0*X(3)/V39-V40/(X(4)*V41)    
      GV(5)=-V38/(X(5)*V39)-2.D0*V40/(X(5)*V41)+V44/(V43*V45)      
      GV(6)=-3.D0*V35*X(1)*V36/(X(6)*V37)+20.D0*X(2)/V41-4.D0   
     /      *V44*V45/(X(6)*V43)     
      GV(7)=A(M)*V35*X(1)*DABS(X(7))**(A(M)-1.D0)/V37-0.5D0*    
     /         V38/(V39*X(7))+V42*X(1)*X(2)/(X(7)*V43) 
      IF (.NOT.INDEX2(5)) GOTO 11       
      DO 20 I=1,7       
20    GG(5,I)=GV(I)     
11    IF (.NOT.INDEX2(6)) GOTO 12       
      DO 30 I=1,7       
30    GG(6,I)=-GV(I)    
12    RETURN    
      END       
C
      SUBROUTINE TP104(MODE)    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D104/A  
      DOUBLEPRECISION A,BX,V1,V2,EPS  
      INTEGER I,J
      GOTO (1,2,3,4,5),MODE     
1     N=8       
      NILI=0    
      NINL=6    
      NELI=0    
      NENL=0    
      DO 33 I=1,8       
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      XL(I)=0.1D0       
      XU(I)=10.D0       
33    CONTINUE  
      A=0.0588D0
      X(1)=6.D0 
      X(2)=3.D0 
      X(3)=0.4D0
      X(4)=0.2D0
      X(5)=6.D0 
      X(6)=6.D0 
      X(7)=1.D0 
      X(8)=0.5D0
      GF(3)=0.D0
      GF(4)=0.D0
      GF(5)=0.D0
      GF(6)=0.D0
      DO 34 I=1,4       
      DO 34 J=1,8       
34    GG(I,J)=0.D0      
      GG(1,1)=-0.1D0    
      GG(2,1)=-0.1D0    
      GG(2,2)=-0.1D0    
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.646511402796D+01     
      XEX(2)=0.223270864907D+01     
      XEX(3)=0.667397491303D+00     
      XEX(4)=0.595756422907D+00     
      XEX(5)=0.593267567811D+01     
      XEX(6)=0.552723456506D+01     
      XEX(7)=0.101332200907D+01     
      XEX(8)=0.400668229166D+00     
      FEX=0.395116343955D+01 
      RETURN    
2     CONTINUE
      EPS = 1.0D-3
      X(1)=DMAX1(EPS,X(1))
      X(2)=DMAX1(EPS,X(2))
      X(7)=DMAX1(EPS,X(7))
      X(8)=DMAX1(EPS,X(8))
      FX = 0.4D0*(X(1)**0.67D0/X(7)**0.67D0 + X(2)**0.67D0/X(8)**0.67D0)
     /     + 10.D0 - X(1) - X(2) 
      RETURN    
3     CONTINUE
      GF(1)=0.268D0*X(1)**(-0.33D0)*X(7)**(-0.67D0)-1.D0
      GF(2)=0.268D0*X(2)**(-0.33D0)*X(8)**(-0.67D0)-1.D0
      GF(7)=-0.268D0*X(1)**0.67D0*X(7)**(-1.67D0)       
      GF(8)=-0.268D0*X(2)**0.67D0*X(8)**(-1.67D0)       
      RETURN    
4     BX = 0.4D0*(X(1)**0.67D0/X(7)**0.67D0 + X(2)**0.67D0/X(8)**0.67D0)
     /     + 10.D0 - X(1) - X(2)  
      IF (INDEX1(1)) G(1)=-A*X(5)*X(7)-0.1D0*X(1)+1.D0  
      IF (INDEX1(2)) G(2)=-A*X(6)*X(8)-0.1D0*X(1)-0.1D0*X(2)+1.D0    
      IF (INDEX1(3)) G(3)=(-4.D0*X(3)-2.D0*X(3)**(-0.71D0))/X(5)
     /                    -A*DABS(X(3))**(-1.3D0)*X(7)+1.D0    
      IF (INDEX1(4)) G(4)=(-4.D0*X(4)-2.D0*DABS(X(4))**(-0.71D0))/X(6)
     /                    -A*DABS(X(4))**(-1.3D0)*X(8)+1.D0    
      IF (INDEX1(5)) G(5)=BX-1.D0       
      IF (INDEX1(6)) G(6)=4.2D0-BX      
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,5)=-A*X(7)   
      GG(1,7)=-A*X(5)   
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,6)=-A*X(8)   
      GG(2,8)=-A*X(6)   
8     IF(.NOT.INDEX2(3)) GOTO 9 
      V1=X(5)**2
      GG(3,3)=(-4.D0+1.42D0*DABS(X(3))**(-1.71D0))/X(5)+1.3D0*
     /         A*DABS(X(3))**(-2.3D0)*X(7)      
      GG(3,5)=(4.D0*X(3)+2.D0*DABS(X(3))**(-0.71D0))/V1 
      GG(3,7)=-A*DABS(X(3))**(-1.3D0)   
9     IF (.NOT.INDEX2(4)) GOTO 10       
      V2=X(6)**2
      GG(4,4)=(-4.D0+1.42D0*DABS(X(4))**(-1.71D0))/X(6)+1.3D0*
     /        A*DABS(X(4))**(-2.3D0)*X(8)      
      GG(4,6)=(4.D0*X(4)+2.D0*DABS(X(4))**(-0.71D0))/V2 
      GG(4,8)=-A*DABS(X(4))**(-1.3D0)   
10    IF (.NOT.INDEX2(5).AND..NOT.INDEX2(6)) GOTO 12    
      GF(1)=0.268D0*DABS(X(1))**(-0.33D0)*DABS(X(7))**(-0.67D0)-1.D0  
      GF(2)=0.268D0*DABS(X(2))**(-0.33D0)*DABS(X(8))**(-0.67D0)-1.D0  
      GF(7)=-0.268D0*DABS(X(1))**0.67D0*DABS(X(7))**(-1.67D0)   
      GF(8)=-0.268D0*DABS(X(2))**0.67D0*DABS(X(8))**(-1.67D0)   
      IF (.NOT.INDEX2(5)) GOTO 11       
      DO 31 I=1,8       
31    GG(5,I)=GF(I)     
11    IF (.NOT.INDEX2(6)) GOTO 12       
      DO 32 I=1,8       
32    GG(6,I)=-GF(I)    
12    RETURN    
      END       
C      
      SUBROUTINE TP105(MODE)    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
     /           /D105/Y(235)
      DIMENSION A(235),B(235),C(235),DA(235,8),D 
     /          B(235,8),DC(235,8)
      DOUBLEPRECISION Y,A,B,C,DA,DB,DC,S,V,V0,V1,V2,V3,V4,V5,V6,
     /          V7,V8,V9,V10,V11,T1,SUM      
      INTEGER I,J
      IF (MODE-1) 1,1,20      
1     N=8       
      NILI=1    
      NINL=0    
      NELI=0    
      NENL=0    
      DO 59 I=1,8       
      LXL(I)=.TRUE.     
59    LXU(I)=.TRUE.     
      XL(1)=1.D-3       
      XL(2)=1.D-3       
      XL(3)=100.D0      
      XL(4)=130.D0      
      XL(5)=170.D0      
      XU(1)=0.499D0     
      XU(2)=0.499D0     
      XU(3)=180.D0      
      XU(4)=210.D0      
      XU(5)=240.D0      
      DO 62 I=6,8       
      XL(I)=5.D0
62    XU(I)=25.D0       
      X(1)=0.1D0
      X(2)=0.2D0
      X(3)=100.D0       
      X(4)=125.D0       
      X(5)=175.D0       
      X(6)=11.2D0       
      X(7)=13.2D0       
      X(8)=15.8D0       
      Y(1)=95.D0
      Y(2)=105.D0       
      DO 30 I=3,6       
30    Y(I)=110.D0       
      DO 31 I=7,10      
31    Y(I)=115.D0       
      DO 32 I=11,25     
32    Y(I)=120.D0       
      DO 33 I=26,40     
33    Y(I)=125.D0       
      DO 34 I=41,55     
34    Y(I)=130.D0       
      DO 35 I=56,68     
35    Y(I)=135.D0       
      DO 36 I=69,89     
36    Y(I)=140.D0       
      DO 37 I=90,101    
37    Y(I)=145.D0       
      DO 38 I=102,118   
38    Y(I)=150.D0       
      DO 39 I=119,122   
39    Y(I)=155.D0       
      DO 40 I=123,142   
40    Y(I)=160.D0       
      DO 41 I=143,150   
41    Y(I)=165.D0       
      DO 42 I=151,167   
42    Y(I)=170.D0       
      DO 43 I=168,175   
43    Y(I)=175.D0       
      DO 44 I=176,181   
44    Y(I)=180.D0       
      DO 45 I=182,187   
45    Y(I)=185.D0       
      DO 46 I=188,194   
46    Y(I)=190.D0       
      DO 47 I=195,198   
47    Y(I)=195.D0       
      DO 48 I=199,201   
48    Y(I)=200.D0       
      DO 49 I=202,204   
49    Y(I)=205.D0       
      DO 50 I=205,212   
50    Y(I)=210.D0       
      Y(213)=215.D0     
      DO 51 I=214,219   
51    Y(I)=220.D0       
      DO 52 I=220,224   
52    Y(I)=230.D0       
      Y(225)=235.D0     
      DO 53 I=226,232   
53    Y(I)=240.D0       
      Y(233)=245.D0     
      Y(234)=260.D0     
      Y(235)=260.D0     
      GG(1,1)=-1.D0     
      GG(1,2)=-1.D0     
      DO 58 I=3,8       
58    GG(1,I)=0.D0      
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.412892753597D+00     
      XEX(2)=0.403352658261D+00     
      XEX(3)=0.131261311486D+03     
      XEX(4)=0.164313514476D+03     
      XEX(5)=0.217422221771D+03     
      XEX(6)=0.122801780396D+02     
      XEX(7)=0.157716989473D+02     
      XEX(8)=0.207468249193D+02     
      FEX=0.113841623960D+04 
      RETURN    
20    IF (MODE - 4) 21,4,5      
21    S=0.D0    
      V=1.D0/DSQRT(8.D0*DATAN(1.D0))    
      V1=X(1)/X(6)      
      V2=X(2)/X(7)      
      V3=(1.D0-X(1)-X(2))/X(8)  
      V4=1.D0/(2.D0*X(6)**2)    
      V5=1.D0/(2.D0*X(7)**2)    
      V6=1.D0/(2.D0*X(8)**2)    
      DO 54 I=1,235    
      A(I)=V1*DEXP(DMAX1(-(Y(I)-X(3))**2*V4,-1.0D1))  
      B(I)=V2*DEXP(DMAX1(-(Y(I)-X(4))**2*V5,-1.0D1))  
      C(I)=V3*DEXP(DMAX1(-(Y(I)-X(5))**2*V6,-1.0D1))  
      V11=(A(I)+B(I)+C(I))*V    
      IF (V11.LE.1.0D-5) GOTO 70  
54    S=S+DLOG(V11)     
      IF (MODE.EQ.3) GOTO 3     
2     FX=-S     
      RETURN    
3     DO 60 I=1,235     
      DO 60 J=1,8       
      DA(I,J)=0.D0      
      DB(I,J)=0.D0      
60    DC(I,J)=0.D0      
      DO 55 I=1,235     
      V0=X(6)**2
      V2=Y(I)-X(3)      
      V1=DEXP(-V2**2/(2.D0*V0)) 
      DA(I,1)=V1/X(6)   
      DA(I,3)=X(1)*V2/X(6)**3*V1
      DA(I,6)=X(1)/V0*(V2**2/V0-1.D0)*V1
      V3=X(7)**2
      V4=Y(I)-X(4)      
      V5=DEXP(-V4**2/(2.D0*V3)) 
      DB(I,2)=V5/X(7)   
      DB(I,4)=X(2)*V4/X(7)**3*V5
      DB(I,7)=X(2)/V3*(V4**2/V3-1.D0)*V5
      V7=X(8)**2
      V9=Y(I)-X(5)      
      V8=DEXP(-V9**2/(2.D0*V7)) 
      V10=1.D0-X(1)-X(2)
      DC(I,1)=-V8/X(8)  
      DC(I,2)=DC(I,1)   
      DC(I,5)=V10*V9/X(8)**3*V8 
      DC(I,8)=V10/V7*(V9**2/V7-1.D0)*V8 
55    CONTINUE  
      DO 57 J=1,8       
      T1=0.D0   
      DO 56 I=1,235     
56    T1=T1+(DA(I,J)+DB(I,J)+DC(I,J))/(A(I)+B(I)+C(I))  
      GF(J)=-T1 
57    CONTINUE  
      RETURN    
70    DO 71 I=1,8       
      SUM=0.D0  
71    SUM=SUM+(X(I)-5.D0)**2    
      FX=SUM+2.09D+3    
      RETURN    
4     IF (INDEX1(1)) G(1)=1.D0-X(1)-X(2)
5     RETURN    
      END       
C
      SUBROUTINE TP106(MODE)    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      GOTO (1,2,3,4,5),MODE     
1     N=8       
      NILI=3    
      NINL=3    
      NELI=0    
      NENL=0    
      DO 23 I=1,3       
      X(I)=5.D+3
      XL(I)=1.D+3       
23    XU(I)=1.D+4       
      XL(1)=100.D0      
      DO 24 I=4,8       
      XL(I)=10.D0       
24    XU(I)=1.D+3       
      X(4)=200.D0       
      X(5)=350.D0       
      X(6)=150.D0       
      X(7)=225.D0       
      X(8)=425.D0       
      DO 6 I=1,8
      LXL(I)=.TRUE.     
6     LXU(I)=.TRUE.     
      DO 22 I=1,6       
      DO 22 J=1,8       
22    GG(I,J)=0.D0      
      GF(1)=1.D0
      GF(2)=1.D0
      GF(3)=1.D0
      DO 20 I=4,8       
20    GF(I)=0.D0
      GG(1,4)=-2.5D-3   
      GG(1,6)=-2.5D-3   
      GG(2,5)=-2.5D-3   
      GG(2,7)=-2.5D-3   
      GG(2,4)=2.5D-3    
      GG(3,5)=0.01D0    
      GG(3,8)=-0.01D0   
      GG(4,4)=-833.33252D0      
      GG(5,5)=-1.25D+3
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.57930668D+03     
      XEX(2)=0.13599707D+04     
      XEX(3)=0.51099707D+04      
      XEX(4)=0.18201770D+03     
      XEX(5)=0.29560117D+03     
      XEX(6)=0.21798230D+03      
      XEX(7)=0.28641653D+03      
      XEX(8)=0.39560117D+03     
C      FEX=0.704933092308D+04 
      FEX=0.70492480D+04
      RETURN    
2     FX=X(1)+X(2)+X(3) 
3     RETURN    
4     IF (INDEX1(1)) G(1)=-2.5D-3*(X(4)+X(6))+1.D0      
      IF (INDEX1(2)) G(2)=-2.5D-3*(X(5)+X(7)-X(4))+1.D0 
      IF (INDEX1(3)) G(3)=-0.01D0*(X(8)-X(5))+1.D0      
      IF (INDEX1(4)) G(4)=(-833.33252D0*X(4)-100.D0*X(1)+8.3333333D+4 
     -                         +X(1)*X(6))
      IF (INDEX1(5)) G(5)=(-1.25D+3*X(5)-X(2)*X(4)+1.25D+3       
     -                        *X(4)+X(2)*X(7))
      IF (INDEX1(6)) G(6)=(-1.25D+6-X(3)*X(5)+2.5D+3*X(5)+       
     -                         X(3)*X(8))
      RETURN    
5     IF (.NOT.INDEX2(4)) GOTO 10       
      GG(4,1)=(-100.D0+X(6))
      GG(4,6)=X(1)
10    IF (.NOT.INDEX2(5)) GOTO 11       
      GG(5,2)=(-X(4)+X(7))
      GG(5,4)=(-X(2)+1.25D+3)
      GG(5,7)=X(2)      
11    IF (.NOT.INDEX2(6)) GOTO 12      
      GG(6,3)=(-X(5)+X(8))
      GG(6,5)=(-X(3)+2.5D+3)
      GG(6,8)=X(3)
12    RETURN    
      END       
C
      SUBROUTINE TP107(MODE)    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D107/V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,
     /                 V15,A,B,C,D,Y1,Y2,Y3,Y4,Y5,Y6  
      DOUBLEPRECISION  V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,
     /                 V15,A,B,C,D,Y1,Y2,Y3,Y4,Y5,Y6  
      INTEGER I,J
      GOTO (1,2,3,4,4),MODE     
1     N=9       
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=6    
      X(1)=0.8D0
      X(2)=0.8D0
      X(3)=0.2D0
      X(4)=0.2D0
      X(5)=1.0454D0     
      X(6)=1.0454D0     
      X(7)=1.0454D0     
      X(8)=0.D0 
      X(9)=0.D0 
      V1=48.4D0/50.176D0
      C=V1*DSIN(0.25D0) 
      D=V1*DCOS(0.25D0) 
      DO 20 I=1,6       
      DO 20 J=1,9       
20    GG(I,J)=0.D0      
      DO 21 I=1,2       
      LXL(I)=.TRUE.     
      LXU(I)=.FALSE.    
21    XL(I)=0.D0
      DO 22 I=3,4       
      LXL(I)=.FALSE.    
22    LXU(I)=.FALSE.    
      DO 23 I=5,7       
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      XL(I)=0.90909D0   
23    XU(I)=1.0909D0    
      DO 24 I=8,9       
      LXL(I)=.FALSE.    
24    LXU(I)=.FALSE.    
      DO 25 I=3,9       
25    GF(I)=0.D0
      GG(1,1)=-1.D0     
      GG(2,2)=-1.D0     
      GG(4,3)=-1.D0     
      GG(5,4)=-1.D0     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.667009506909D+00     
      XEX(2)=0.102238816675D+01     
      XEX(3)=0.228287932605D+00     
      XEX(4)=0.184821729352D+00     
      XEX(5)=0.10909D+01     
      XEX(6)=0.10909D+01     
      XEX(7)=0.106903593236D+01     
      XEX(8)=0.106612642267D+00     
      XEX(9)=-0.338786658776D+00     
      FEX=0.505501180339D+04
      RETURN    
2     FX=3.D+3*X(1)+1.D+3*X(1)**3+2.D+3*X(2)+666.66666667D0*X(2)**3
      RETURN    
3     GF(1)=3.D+3+3.D+3*X(1)**2 
      GF(2)=2.D+3+2.000001D+3*X(2)**2   
      RETURN    
 4    Y1=DSIN(X(8))     
      Y2=DCOS(X(8))     
      Y3=DSIN(X(9))     
      Y4=DCOS(X(9))     
      Y5=DSIN(X(8)-X(9))
      Y6=DCOS(X(8)-X(9))
      IF (MODE.EQ.5) GOTO 5     
      IF (INDEX1(1)) G(1)=0.4D0-X(1)+2.D0*C*X(5)**2+X(5)*X(6)
     /    *(-D*Y1-C*Y2)+X(5)*X(7)*(-D*Y3-C*Y4)    
      IF (INDEX1(2)) G(2)=0.4D0-X(2)+2.D0*C*X(6)**2+X(5)*X(6)
     /    *(D*Y1-C*Y2)+X(6)*X(7)*(D*Y5-C*Y6)     
      IF (INDEX1(3)) G(3)=0.8D0+2.D0*C*X(7)**2+X(5)*X(7)
     /    *(D*Y3-C*Y4)+X(6)*X(7)*(-D*Y5-C*Y6)   
      IF (INDEX1(4)) G(4)=0.2D0-X(3)+2.D0*D*X(5)**2-X(5)
     /    *X(6)*(-C*Y1+D*Y2)-X(5)*X(7)*(-C*Y3+D*Y4)     
      IF (INDEX1(5)) G(5)=0.2D0-X(4)+2.D0*D*X(6)**2
     /    -X(5)*X(6)*(C*Y1+D*Y2)-X(6)*X(7)*(C*Y5+D*Y6)       
      IF (INDEX1(6)) G(6)=-0.337D0+2.D0*D*X(7)**2-X(5)*X(7)     
     /    *(C*Y3+D*Y4)-X(6)*X(7)*(-C*Y5+D*Y6)
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      V1=-D*Y1-C*Y2     
      V2=-D*Y3-C*Y4     
      GG(1,5)=4.D0*C*X(5)+X(6)*V1+X(7)*V2       
      GG(1,6)=X(5)*V1   
      GG(1,7)=X(5)*V2   
      GG(1,8)=X(5)*X(6)*(-D*Y2+C*Y1)    
      GG(1,9)=X(5)*X(7)*(-D*Y4+C*Y3)    
7     IF (.NOT.INDEX2(2)) GOTO 8
      V2=D*Y1-C*Y2      
      V3=D*Y6+C*Y5      
      V4=D*Y5-C*Y6      
      GG(2,5)=X(6)*V2   
      GG(2,6)=4.D0*C*X(6)+X(5)*V2+X(7)*V4       
      GG(2,7)=X(6)*V4   
      GG(2,8)=X(5)*X(6)*(D*Y2+C*Y1)+X(6)*X(7)*V3
      GG(2,9)=-X(6)*X(7)*V3     
8     IF(.NOT.INDEX2(3)) GOTO 9 
      V5=D*Y3-C*Y4      
      V6=-D*Y5-C*Y6     
      V7=-D*Y6+C*Y5     
      GG(3,5)=X(7)*V5   
      GG(3,6)=X(7)*V6   
      GG(3,7)=4.D0*C*X(7)+X(5)*V5+X(6)*V6       
      GG(3,8)=X(6)*X(7)*V7      
      GG(3,9)=X(5)*X(7)*(D*Y4+C*Y3)-X(6)*X(7)*V7
9     IF (.NOT.INDEX2(4)) GOTO 10       
      V8=-C*Y1+D*Y2     
      V9=-C*Y3+D*Y4     
      GG(4,5)=4.D0*D*X(5)-X(6)*V8-X(7)*V9       
      GG(4,6)=-X(5)*V8  
      GG(4,7)=-X(5)*V9  
      GG(4,8)=X(5)*X(6)*(C*Y2+D*Y1)     
      GG(4,9)=X(5)*X(7)*(C*Y4+D*Y3)     
10    IF (.NOT.INDEX2(5)) GOTO 11       
      V10=C*Y1+D*Y2     
      V11=C*Y5+D*Y6     
      V12=(C*Y6-D*Y5)*X(6)      
      GG(5,5)=-X(6)*V10 
      GG(5,6)=4.D0*D*X(6)-X(5)*V10-X(7)*V11     
      GG(5,7)=-X(6)*V11 
      GG(5,8)=-X(5)*X(6)*(C*Y2-D*Y1)-X(7)*V12   
      GG(5,9)=X(7)*V12  
11    IF (.NOT.INDEX2(6)) GOTO 12       
      V13=C*Y3+D*Y4     
      V14=-C*Y5+D*Y6    
      V15=(C*Y6+D*Y5)*X(6)*X(7) 
      GG(6,5)=-X(7)*V13 
      GG(6,6)=-X(7)*V14 
      GG(6,7)=4.D0*D*X(7)-X(5)*V13-X(6)*V14     
      GG(6,8)=V15       
      GG(6,9)=-X(5)*X(7)*(C*Y4-D*Y3)-V15
12    RETURN    
      END       
C
      SUBROUTINE TP108(MODE)    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      GOTO (1,2,3,4,5),MODE     
1     N=9       
      NILI=0    
      NINL=13   
      NELI=0    
      NENL=0    
      DO 21 J=1,9       
21    X(J)=1.D0 
      DO 22 I=1,9 
      XL(I)=0.0D0      
      LXL(I)=.TRUE.
C      LXL(I)=.FALSE.
      XU(I)=1.0D0    
22    LXU(I)=.TRUE.    
      LXL(9)=.TRUE.     
      XL(9)=0.D0
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.884129216724D+00     
      XEX(2)=0.467242472598D+00     
      XEX(3)=0.374207573677D-01     
      XEX(4)=0.999299598210D+00     
      XEX(5)=0.884129216724D+00     
      XEX(6)=0.467242472594D+00     
      XEX(7)=0.374207573618D-01     
      XEX(8)=0.999299598210D+00     
      XEX(9)=0.261984643608D-19     
      FEX=-0.866025403841D+00 
      RETURN    
2     FX=-0.5D0*(X(1)*X(4)-X(2)*X(3)+X(3)*X(9)-X(5)*X(9)+       
     -X(5)*X(8)-X(6)*   
     -X(7))     
      RETURN    
3     GF(1)=-0.5D0*X(4) 
      GF(2)=0.5D0*X(3)  
      GF(3)=0.5D0*(X(2)-X(9))   
      GF(4)=-0.5D0*X(1) 
      GF(5)=0.5D0*(X(9)-X(8))   
      GF(6)=0.5D0*X(7)  
      GF(7)=0.5D0*X(6)  
      GF(8)=-0.5D0*X(5) 
      GF(9)=-GF(8)-GF(2)
      RETURN    
4     IF (INDEX1(1)) G(1)=1.D0-X(3)**2-X(4)**2  
      IF (INDEX1(2)) G(2)=1.D0-X(9)**2  
      IF (INDEX1(3)) G(3)=1.D0-X(5)**2-X(6)**2  
      IF (INDEX1(4)) G(4)=1.D0-X(1)**2-(X(2)-X(9))**2   
      IF (INDEX1(5)) G(5)=1.D0-(X(1)-X(5))**2-(X(2)-X(6))**2       
      IF (INDEX1(6)) G(6)=1.D0-(X(1)-X(7))**2-(X(2)-X(8))**2       
      IF (INDEX1(7)) G(7)=1.D0-(X(3)-X(5))**2-(X(4)-X(6))**2       
      IF (INDEX1(8)) G(8)=1.D0-(X(3)-X(7))**2-(X(4)-X(8))**2       
      IF (INDEX1(9)) G(9)=1.D0-X(7)**2-(X(8)-X(9))**2   
      IF (INDEX1(10)) G(10)=X(1)*X(4)-X(2)*X(3) 
      IF (INDEX1(11))G(11)=X(3)*X(9)    
      IF (INDEX1(12)) G(12)=-X(5)*X(9)  
      IF (INDEX1(13)) G(13)=X(5)*X(8)-X(6)*X(7) 
      RETURN    
5     DO 23 I=1,13      
      DO 23 J=1,9       
23    GG(I,J)=0.D0      
      IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,3)=-2.D0*X(3)
      GG(1,4)=-2.D0*X(4)
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,9)=-2.D0*X(9)
8     IF(.NOT.INDEX2(3)) GOTO 9 
      GG(3,5)=-2.D0*X(5)
      GG(3,6)=-2.D0*X(6)
9     IF (.NOT.INDEX2(4)) GOTO 10       
      GG(4,1)=-2.D0*X(1)
      GG(4,2)=-2.D0*(X(2)-X(9)) 
      GG(4,9)=-GG(4,2)  
10    IF (.NOT.INDEX2(5)) GOTO 11       
      GG(5,1)=-2.D0*(X(1)-X(5)) 
      GG(5,2)=-2.D0*(X(2)-X(6)) 
      GG(5,5)=-GG(5,1)  
      GG(5,6)=-GG(5,2)  
11    IF (.NOT.INDEX2(6)) GOTO 12       
      GG(6,1)=-2.D0*(X(1)-X(7)) 
      GG(6,2)=-2.D0*(X(2)-X(8)) 
      GG(6,7)=-GG(6,1)  
      GG(6,8)=-GG(6,2)  
12    IF (.NOT.INDEX2(7)) GOTO 13       
      GG(7,3)=-2.D0*(X(3)-X(5)) 
      GG(7,4)=-2.D0*(X(4)-X(6)) 
      GG(7,5)=-GG(7,3)  
      GG(7,6)=-GG(7,4)  
13    IF (.NOT.INDEX2(8)) GOTO 14       
      GG(8,3)=-2.D0*(X(3)-X(7)) 
      GG(8,4)=-2.D0*(X(4)-X(8)) 
      GG(8,7)=-GG(8,3)  
      GG(8,8)=-GG(8,4)  
14    IF (.NOT.INDEX2(9)) GOTO 15       
      GG(9,7)=-2.D0*X(7)
      GG(9,8)=-2.D0*(X(8)-X(9)) 
      GG(9,9)=-GG(9,8)  
15    IF (.NOT.INDEX2(10)) GOTO 16      
      GG(10,1)=X(4)     
      GG(10,2)=-X(3)    
      GG(10,3)=-X(2)    
      GG(10,4)=X(1)     
16    IF (.NOT.INDEX2(11)) GOTO 17      
      GG(11,3)=X(9)     
      GG(11,9)=X(3)     
17    IF (.NOT.INDEX2(12)) GOTO 18      
      GG(12,5)=-X(9)    
      GG(12,9)=-X(5)    
18    IF (.NOT.INDEX2(13)) GOTO 19      
      GG(13,5)=X(8)     
      GG(13,6)=-X(7)    
      GG(13,7)=-X(6)    
      GG(13,8)=X(5)     
19    RETURN    
      END       
C
      SUBROUTINE TP109(MODE)    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D109/A,RA,B,C,HV1,V1,V2,V3,V4,V5,V6,V7,V8,V9,
     /                V10,V11,V12,V13,V14,V15,V16,V17,V18,V19,V20    
      DOUBLEPRECISION A,RA,B,C,HV1,V1,V2,V3,V4,V5,V6,V7,V8,
     /                V9,V10,V11,V12,V13,V14,V15,V16,V17,V18,V19,V20  
      INTEGER I,J
      GOTO (1,2,3,4,5),MODE     
1     N=9       
      NILI=2    
      NINL=2    
      NELI=0    
      NENL=6    
      A=50.176D0
      RA=1.D0/A 
      B=DSIN(0.25D0)    
      C=DCOS(0.25D0)    
      DO 20 I=1,2       
      XL(I)=0.D0
      LXL(I)=.TRUE.     
      LXU(I)=.FALSE.    
      XL(I+2)=-0.55D0   
      XU(2+I)=0.55D0    
      XL(7+I)=-400.D0   
      XU(7+I)=800.D0    
      XL(4+I)=196.D0    
      XU(4+I)=252.D0    
20    CONTINUE  
      XL(7)=196.D0      
      XU(7)=252.D0      
      DO 21 I=3,9       
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
21    CONTINUE  
      DO 22 J=1,9       
      X(J)=0.D0 
      DO 22 I=1,10      
      GG(I,J)=0.D0      
22    CONTINUE  
      DO 30 I=3,9       
30    GF(I)=0.D0
      X(5)=250.D0
      X(6)=250.D0
      X(7)=200.D0
      GG(1,3)=-1.D0     
      GG(1,4)=1.D0      
      GG(2,3)=1.D0      
      GG(2,4)=-1.D0     
      GG(5,1)=-1.D0     
      GG(6,2)=-1.D0     
      GG(8,8)=1.D0      
      GG(9,9)=1.D0      
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.674888100445D+03     
      XEX(2)=0.113417039470D+04     
      XEX(3)=0.133569060261D+00     
      XEX(4)=-0.371152592466D+00     
      XEX(5)=0.252D+03     
      XEX(6)=0.252D+03     
      XEX(7)=0.201464535316D+03     
      XEX(8)=0.426660777226D+03     
      XEX(9)=0.368494083867D+03     
      FEX=0.536206927538D+04
      RETURN    
2     FX=3.D0*X(1)+1.D-6*X(1)**3+0.522074D-6*X(2)**3+2.D0*X(2) 
      RETURN    
3     GF(1)=3.D-6*X(1)**2+3.D0  
      GF(2)=2.D0+1.566222D-6*X(2)**2   
      RETURN    
4     IF (INDEX1(1)) G(1)=X(4)-X(3)+0.55D0      
      IF (INDEX1(2)) G(2)=X(3)-X(4)+0.55D0      
      IF (INDEX1(3)) G(3)=(2.25D+6-X(1)**2-X(8)**2)
      IF (INDEX1(4)) G(4)=(2.25D+6-X(2)**2-X(9)**2)
      IF (INDEX1(5)) G(5)=(X(5)*X(6)*DSIN(-X(3)-0.25D0)
     /   +X(5)*X(7)*DSIN(-X(4)-0.25D0)+2.D0*X(5)**2*B)*RA+400.D0-X(1)  
      IF (INDEX1(6)) G(6)=(X(5)*X(6)*DSIN(X(3)-0.25D0)+X(6)*X(7)
     /    *DSIN(X(3)-X(4)-0.25D0)+2.D0*X(6)**2*B)*RA+400.D0-X(2)  
      IF (INDEX1(7)) G(7)=(X(5)*X(7)*DSIN(X(4)-0.25D0)+X(6)*X(7)
     /    *DSIN(X(4)-X(3)-0.25D0)+2.D0*X(7)**2*B)*RA+881.779D0     
      IF (INDEX1(8)) G(8)=X(8)+(X(5)*X(6)*DCOS(-X(3)-0.25D0)+X(5)*X(7)
     /    *DCOS(-X(4)-0.25D0)-2.D0*X(5)**2*C)*RA+0.7533D-3*X(5)**2
     /    -200.0D0  
      IF (INDEX1(9)) G(9)=X(9)+(X(5)*X(6)*DCOS(X(3)-0.25D0)
     /     +X(7)*X(6)*DCOS(X(3)-X(4)-0.25D0)-2.D0*X(6)**2*C)*RA
     /     +0.7533D-3*X(6)**2-200.0D0
      IF (INDEX1(10)) G(10)=(X(5)*X(7)*DCOS(X(4)-0.25D0)+X(6)*X(7)
     /     *DCOS(X(4)-X(3)-0.25D0)-2.D0*X(7)**2*C)*RA+0.7533D-3*X(7)**2
     /     -22.938D0
      RETURN    
5     CONTINUE  
8     IF (.NOT.INDEX2(3)) GOTO 9
      GG(3,1)=-2.D0*X(1)
      GG(3,8)=-2.D0*X(8)
9     IF (.NOT.INDEX2(4)) GOTO 10       
      GG(4,2)=-2.D0*X(2)
      GG(4,9)=-2.D0*X(9)
10    IF (.NOT.INDEX2(5)) GOTO 11       
      V1=DSIN(-X(3)-0.25D0)     
      V2=DSIN(-X(4)-0.25D0)     
      V3=X(5)*RA
      GG(5,3)=-X(6)*V3*DCOS(-X(3)-0.25D0)       
      GG(5,4)=-X(7)*V3*DCOS(-X(4)-0.25D0)       
      GG(5,5)=(X(6)*V1+X(7)*V2+4.D0*X(5)*B)*RA  
      GG(5,6)=V3*V1     
      GG(5,7)=V3*V2     
11    IF (.NOT.INDEX2(6)) GOTO 12       
      HV1=X(3)-X(4)-0.25D0      
      V3=DCOS(HV1)      
      V4=DSIN(X(3)-0.25D0)      
      V5=X(6)*RA
      V6=DSIN(HV1)      
      GG(6,3)=X(5)*V5*DCOS(X(3)-0.25D0)+X(7)*V5*V3      
      GG(6,4)=-X(7)*V5*V3       
      GG(6,5)=V5*V4     
      GG(6,6)=(X(5)*V4+X(7)*V6)*RA+4.D0*V5*B    
      GG(6,7)=V5*V6     
12    IF (.NOT.INDEX2(7)) GOTO 13       
      HV1=X(4)-X(3)-0.25D0      
      V7=X(7)*RA
      V8=DCOS(HV1)      
      V9=DSIN(X(4)-0.25D0)      
      V10=DSIN(HV1)     
      GG(7,3)=-X(6)*V7*V8       
      GG(7,4)=X(5)*V7*DCOS(X(4)-0.25D0)+X(6)*V7*V8      
      GG(7,5)=V7*V9     
      GG(7,6)=V7*V10    
      GG(7,7)=(X(5)*V9+X(6)*V10)*RA+4.D0*V7*B   
13    IF (.NOT.INDEX2(8)) GOTO 14       
      V11=X(5)*RA       
      V12=DCOS(-X(3)-0.25D0)*RA 
      V13=DCOS(-X(4)-0.25D0)*RA 
      GG(8,3)=X(6)*V11*DSIN(-X(3)-0.25D0)       
      GG(8,4)=X(7)*V11*DSIN(-X(4)-0.25D0)       
      GG(8,5)=X(6)*V12+X(7)*V13-4.D0*V11*C+1.5066D-3*X(5)       
      GG(8,6)=X(5)*V12  
      GG(8,7)=X(5)*V13  
14    IF (.NOT.INDEX2(9)) GOTO 15       
      HV1=X(3)-X(4)-0.25D0      
      V14=DSIN(HV1)*X(6)*RA     
      V15=DCOS(X(3)-0.25D0)*RA  
      V16=DCOS(HV1)*RA  
      GG(9,3)=-X(5)*X(6)*DSIN(X(3)-0.25D0)*RA-X(7)*V14  
      GG(9,4)=X(7)*V14  
      GG(9,5)=X(6)*V15  
      GG(9,6)=X(5)*V15+X(7)*V16-4.D0*X(6)*C*RA+1.5066D-3*X(6)      
      GG(9,7)=X(6)*V16  
15    IF (.NOT.INDEX2(10)) GOTO 16      
      HV1=X(4)-X(3)-0.25D0      
      V17=DSIN(HV1)*X(6)*RA     
      V18=DCOS(X(4)-0.25D0)*RA  
      V19=DCOS(HV1)*RA  
      V20=X(7)*RA       
      GG(10,3)=X(7)*V17 
      GG(10,4)=-X(5)*V20*DSIN(X(4)-0.25D0)-X(7)*V17     
      GG(10,5)=X(7)*V18 
      GG(10,6)=X(7)*V19 
      GG(10,7)=X(5)*V18+X(6)*V19-4.D0*V20*C+1.5066D-3*X(7) 
16    RETURN    
      END       
C
      SUBROUTINE TP110(MODE)    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION T,S,U,SUM  
      INTEGER I
      GOTO (1,2,2,4,4),MODE     
1     N=10      
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=0    
      DO 20 I=1,10      
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      XL(I)=2.001D0     
      XU(I)=9.999D0     
      X(I)=9.D0 
20    CONTINUE  
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.935025654733D+01     
      XEX(2)=0.935025654733D+01     
      XEX(3)=0.935025654733D+01     
      XEX(4)=0.935025654733D+01     
      XEX(5)=0.935025654733D+01     
      XEX(6)=0.935025654733D+01     
      XEX(7)=0.935025654733D+01     
      XEX(8)=0.935025654733D+01     
      XEX(9)=0.935025654733D+01     
      XEX(10)=0.935025654733D+01     
      FEX=-0.457784697153D+02 
      RETURN    
 2    T=1.D0    
      DO 30 I=1,10      
30    T=T*X(I)  
      S=(DABS(T))**0.2D0
      S=DSIGN(S,T)      
      IF (MODE.EQ.3) GOTO 3     
      U=0.D0    
      DO 31 I=1,10      
      IF ((X(I)-2.D0).LE.0.D0.OR.(10.D0-X(I)).LE.0.D0) GOTO 33  
31    U=U+(DLOG(X(I)-2.D0))**2+(DLOG(10.D0-X(I)))**2    
      FX=U-S    
      RETURN    
33    SUM=0.D0  
      DO 34 I=1,10      
34    SUM=SUM+(X(I)-5.D0)**2    
      FX=SUM+1.D+3-45.8D0       
      RETURN    
3     DO 32 I=1,10      
      IF ((X(I)-2.D0).LE.0.D0.OR.(10.D0-X(I)).LE.0.D0) GOTO 35  
      GF(I)=2.D0*(DLOG(X(I)-2.D0)/(X(I)-2.D0)-DLOG(10.D0-X(I))/ 
     /     (10.D0-X(I)))-S/X(I)*0.2D0   
      GOTO 32   
35    GF(I)=2.D0*(X(I)-5.D0)    
32    CONTINUE  
4     RETURN    
      END       
C
      SUBROUTINE TP111(MODE)    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D111/C(10)   
      DOUBLEPRECISION  C,T,S    
      INTEGER I,J
      GOTO (1,2,2,4,5),MODE     
1     N=10      
      NILI=0    
      NINL=0    
      NELI=0    
      NENL=3    
      DO 20 J=1,10      
      X(J)=-2.3D0       
      DO 20 I=1,3       
20    GG(I,J)=0.D0      
      C(1)=-6.089D0     
      C(2)=-17.164D0    
      C(3)=-34.054D0    
      C(4)=-5.914D0     
      C(5)=-24.721D0    
      C(6)=-14.986D0    
      C(7)=-24.1D0      
      C(8)=-10.708D0    
      C(9)=-26.662D0    
      C(10)=-22.179D0   
      DO 6 I=1,10       
      XL(I)=-100.D0     
      XU(I)=100.D0      
      LXL(I)=.TRUE.     
6     LXU(I)=.TRUE.     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=-0.320121253241D+01     
      XEX(2)=-0.191205959435D+01     
      XEX(3)=-0.244441308369D+00     
      XEX(4)=-0.653748856532D+01     
      XEX(5)=-0.723152425984D+00     
      XEX(6)=-0.726773826993D+01     
      XEX(7)=-0.359671064233D+01     
      XEX(8)=-0.401776873216D+01     
      XEX(9)=-0.328746169619D+01     
      XEX(10)=-0.233558183059D+01     
      FEX=-0.477610902637D+02 
      RETURN    
 2    T=0.D0    
      DO 30 I=1,10      
30    T=T+DEXP(X(I))    
      IF (MODE.EQ.3) GOTO 3     
      S=0.D0    
      DO 31 I=1,10      
31    S=S+DEXP(X(I))*(C(I)+X(I)-DLOG(T))
      FX=S      
      RETURN    
3     DO 33 I=1,10      
33    GF(I)=DEXP(X(I))*(C(I)+X(I)-DLOG(T))      
      RETURN    
4     IF (INDEX1(1)) G(1)=DEXP(X(1))+2.D0*DEXP(X(2))+2.D0*DEXP  
     -(X(3))+DEXP(X(6)  
     -)+DEXP(X(10))-2.D0
      IF (INDEX1(2)) G(2)=DEXP(X(4))+2.D0*DEXP(X(5))+DEXP(X(    
     -6))+DEXP(X(7))-1.D0       
      IF (INDEX1(3)) G(3)=DEXP(X(3))+DEXP(X(7))+DEXP(X(8))      
     -+2.D0*DEXP(X(9))+ 
     -DEXP(X(10))-1.D0  
      RETURN    
5     IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=DEXP(X(1))
      GG(1,2)=2.D0*DEXP(X(2))   
      GG(1,3)=2.D0*DEXP(X(3))   
      GG(1,6)=DEXP(X(6))
      GG(1,10)=DEXP(X(10))      
7     IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,4)=DEXP(X(4))
      GG(2,5)=2.D0*DEXP(X(5))   
      GG(2,6)=DEXP(X(6))
      GG(2,7)=DEXP(X(7))
8     IF (.NOT.INDEX2(3)) GOTO 9
      GG(3,3)=DEXP(X(3))
      GG(3,7)=DEXP(X(7))
      GG(3,8)=DEXP(X(8))
      GG(3,9)=2.D0*DEXP(X(9))   
      GG(3,10)=DEXP(X(10))      
9     RETURN    
      END       
C
      SUBROUTINE TP112(MODE)    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D112/C(10)   
      DOUBLEPRECISION C,T,DLOGT,S   
      INTEGER I,J
      GOTO (1,2,2,4,5),MODE     
1     N=10      
      NILI=0    
      NINL=0    
      NELI=3    
      NENL=0    
      DO 6 I=1,10       
      XL(I)=1.D-4       
      LXL(I)=.TRUE.     
6     LXU(I)=.FALSE.    
      C(1)=-6.089D0     
      C(2)=-17.164D0    
      C(3)=-34.054D0    
      C(4)=-5.914D0     
      C(5)=-24.721D0    
      C(6)=-14.986D0    
      C(7)=-24.1D0      
      C(8)=-10.708D0    
      C(9)=-26.662D0    
      C(10)=-22.179D0   
      DO 20 J=1,10      
      X(J)=0.1D0
      DO 20 I=1,3       
20    GG(I,J)=0.D0      
      GG(1,1)=1.D0      
      GG(1,2)=2.D0      
      GG(1,3)=2.D0      
      GG(1,6)=1.D0      
      GG(1,10)=1.D0     
      GG(2,4)=1.D0      
      GG(2,5)=2.D0      
      GG(2,6)=1.D0      
      GG(2,7)=1.D0      
      GG(3,3)=1.D0      
      GG(3,7)=1.D0      
      GG(3,8)=1.D0      
      GG(3,9)=2.D0      
      GG(3,10)=1.D0     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.177354776881D-01     
      XEX(2)=0.820018011109D-01     
      XEX(3)=0.882564558920D+00     
      XEX(4)=0.723325625629D-03     
      XEX(5)=0.490785079062D+00     
      XEX(6)=0.433546900325D-03     
      XEX(7)=0.172729773078D-01     
      XEX(8)=0.776563912291D-02     
      XEX(9)=0.198492864597D-01     
      XEX(10)=0.526982611793D-01     
      FEX=-0.47761086D+2  
      RETURN    
2     T=0.D0    
      DO 30 I=1,10      
30    T=T+X(I)  
      IF (MODE.EQ.3) GOTO 3     
      IF (T.LT.1.D-5) GOTO 34    
      DLOGT=DLOG(T)     
      S=0.D0    
      DO 31 I=1,10      
      IF (X(I).LT.0.D0) GOTO 34  
   31 S=S+X(I)*(C(I)+DLOG(X(I))-DLOGT)  
      FX=S
      RETURN    
   34 S=0.D0    
      DO 35 I=1,10      
   35 IF (X(I).LT.0.D0) S=S+(X(I)-5.D0)**2       
      FX=(S+1.D+3-47.8D0)
      RETURN    
    3 IF (T.LT.1.D-5) GOTO 36    
      DLOGT=DLOG(T)     
      DO 33 I=1,10      
      IF (X(I).LT.0.D0) GOTO 36  
   33 GF(I)=(C(I)+DLOG(X(I))-DLOGT)
      RETURN    
   36 DO 37 I=1,10      
      GF(I)=0.D0
   37 IF (X(I).LT.0.D0) GF(I)=2.D0*(X(I)-5.D0)
      RETURN    
    4 IF (INDEX1(1)) G(1)=X(1)+2.D0*X(2)+2.D0*X(3)+X(6)+X(10)-2.0   
      IF (INDEX1(2)) G(2)=X(4)+2.D0*X(5)+X(6)+X(7)-1.0 
      IF (INDEX1(3)) G(3)=X(3)+X(7)+X(8)+2.D0*X(9)+X(10)-1.0      
5     RETURN    
      END       
C
      SUBROUTINE TP113(MODE)    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      GOTO (1,2,3,4,5),MODE     
1     N=10      
      NILI=3    
      NINL=5    
      NELI=0    
      NENL=0    
      DO 20 I=1,8       
      DO 20 J=1,10      
20    GG(I,J)=0.D0      
      X(1)=2.D0 
      X(2)=3.D0 
      X(3)=5.D0 
      X(4)=5.D0 
      X(5)=1.D0 
      X(6)=2.D0 
      X(7)=7.D0 
      X(8)=3.D0 
      X(9)=6.D0 
      X(10)=10.D0       
      DO 6 I=1,10       
      LXL(I)=.FALSE.    
6     LXU(I)=.FALSE.    
      GG(1,1)=-4.D0     
      GG(1,2)=-5.D0     
      GG(1,7)=3.D0      
      GG(1,8)=-9.D0     
      GG(2,1)=-10.D0    
      GG(2,2)=8.D0      
      GG(2,7)=17.D0     
      GG(2,8)=-2.D0     
      GG(3,1)=8.D0      
      GG(3,2)=-2.D0     
      GG(3,9)=-5.D0     
      GG(3,10)=2.D0     
      GG(4,4)=7.D0      
      GG(5,2)=-8.D0     
      GG(5,4)=2.D0      
      GG(6,6)=1.D0      
      GG(7,5)=-14.D0    
      GG(7,6)=6.D0      
      GG(8,1)=3.D0      
      GG(8,2)=-6.D0     
      GG(8,10)=7.D0     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.217199637118D+01     
      XEX(2)=0.236368297378D+01     
      XEX(3)=0.877392573847D+01     
      XEX(4)=0.509598448797D+01     
      XEX(5)=0.990654764992D+00     
      XEX(6)=0.143057397893D+01     
      XEX(7)=0.132164420805D+01     
      XEX(8)=0.982872580801D+01     
      XEX(9)=0.828009167017D+01     
      XEX(10)=0.837592666387D+01     
      FEX=0.243062090641D+02 
      RETURN    
2     FX=X(1)**2+X(2)**2+X(1)*X(2)-14.D0*X(1)-16.D0*X(2)
     /     +(X(3)-10.D0)**2+4.D0*(X(4)-5.D0)**2+(X(5)-3.D0)**2
     /     +2.D0*(X(6)-1.D0)**2+5.D0*X(7)**2+7.D0*(X(8)-11.D0)**2
     /     +2.D0*(X(9)-10.D0)**2+(X(10)-7.D0)**2+45.D0      
      RETURN    
3     GF(1)=2.D0*X(1)+X(2)-14.D0
      GF(2)=2.D0*X(2)+X(1)-16.D0
      GF(3)=2.D0*(X(3)-10.D0)   
      GF(4)=8.D0*(X(4)-5.D0)    
      GF(5)=2.D0*(X(5)-3.D0)    
      GF(6)=4.D0*(X(6)-1.D0)    
      GF(7)=10.D0*X(7)  
      GF(8)=14.D0*(X(8)-11.D0)  
      GF(9)=4.D0*(X(9)-10.D0)   
      GF(10)=2.D0*(X(10)-7.D0)  
      RETURN    
4     IF (INDEX1(1)) G(1)=-4.D0*X(1)-5.D0*X(2)+3.D0*X(7)-9.D0*X(8)
     /     +105.D0
      IF (INDEX1(2)) G(2)=-10.D0*X(1)+8.D0*X(2)+17.D0*X(7)
     /    -2.D0*X(8)     
      IF (INDEX1(3)) G(3)=8.D0*X(1)-2.D0*X(2)-5.D0*X(9)
     /     +2.D0*X(10)+12.D0 
      IF (INDEX1(4)) G(4)=-3.D0*(X(1)-2.D0)**2-4.D0*(X(2)-3.D0)**2
     /     -2.D0*X(3)**2+7.D0*X(4)+120.D0      
      IF (INDEX1(5)) G(5)=-5.D0*X(1)**2-8.D0*X(2)-(X(3)-6.D0)**2
     /     +2.D0*X(4)+40.D0       
      IF (INDEX1(6)) G(6)=-0.5D0*(X(1)-8.D0)**2-2.D0*(X(2)
     /     -4.D0)**2-3.D0*X(5)**2+X(6)+30.D0 
      IF (INDEX1(7)) G(7)=-X(1)**2-2.D0*(X(2)-2.D0)**2
     /     +2.D0*X(1)*X(2)-14.D0*X(5)+6.D0*X(6)       
      IF (INDEX1(8)) G(8)=3.D0*X(1)-6.D0*X(2)-12.D0*(X(9)-8.D0)**2
     /     +7.D0*X(10)    
      RETURN    
5     IF (.NOT.INDEX2(4)) GOTO 10       
      GG(4,1)=-6.D0*(X(1)-2.D0) 
      GG(4,2)=-8.D0*(X(2)-3.D0) 
      GG(4,3)=-4.D0*X(3)
10    IF (.NOT.INDEX2(5)) GOTO 11       
      GG(5,1)=-10.D0*X(1)       
      GG(5,3)=-2.D0*(X(3)-6.D0) 
11    IF (.NOT.INDEX2(6)) GOTO 12       
      GG(6,1)=8.D0-X(1) 
      GG(6,2)=-4.D0*(X(2)-4.D0) 
      GG(6,5)=-6.D0*X(5)
12    IF (.NOT.INDEX2(7)) GOTO 13       
      GG(7,1)=-2.D0*X(1)+2.D0*X(2)      
      GG(7,2)=-4.D0*(X(2)-2.D0)+2.D0*X(1)       
13    IF (INDEX2(8)) GG(8,9)=-24.D0*(X(9)-8.D0) 
      RETURN    
      END       
C
      SUBROUTINE TP114(MODE)    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION V1,V2,V3   
      GOTO (1,2,3,4,5),MODE     
1     N=10      
      NILI=4    
      NINL=4    
      NELI=1    
      NENL=2    
      X(1)=1.745D+3     
      X(2)=1.2D+4       
      X(3)=110.D0       
      X(4)=3.048D+3     
      X(5)=1.974D+3     
      X(6)=89.2D0       
      X(7)=92.8D0       
      X(8)=8.D0 
      X(9)=3.6D0
      X(10)=145.D0      
      DO 25 I=1,5       
25    XL(I)=1.D-5       
      XL(6)=85.D0       
      XL(7)=90.D0       
      XL(8)=3.D0
      XL(9)=1.2D0       
      XL(10)=145.D0     
      XU(1)=2.D+3       
      XU(2)=1.6D+4      
      XU(3)=120.D0      
      XU(4)=5.D+3       
      XU(5)=2.D+3       
      XU(6)=93.D0       
      XU(7)=95.D0       
      XU(8)=12.D0       
      XU(9)=4.D0
      XU(10)=162.D0     
      DO 6 I=1,10       
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      DO 6 J=1,11       
6     GG(J,I)=0.D0      
      GF(1)=5.04D0*1.0D-4      
      GF(2)=0.035D0*1.0D-4     
      GF(3)=10.D0*1.0D-4    
      GF(4)=0.0D0   
      GF(5)=3.36D0*1.0D-4      
      GF(6)=0.D0
      GF(7)=0.0D0   
      GF(8)=0.D0
      GF(9)=0.D0
      GF(10)=0.D0       
      GG(1,9)=-0.9D0    
      GG(1,10)=-0.222D0 
      GG(2,7)=3.D0      
      GG(2,10)=-0.99D0  
      GG(3,9)=10.D0/9.D0
      GG(3,10)=0.222D0  
      GG(4,7)=-3.D0     
      GG(4,10)=100.D0/99.D0     
      GG(5,4)=-0.99D0   
      GG(6,6)=0.325D0   
      GG(6,7)=-0.99D0   
      GG(7,4)=100.D0/99.D0      
      GG(8,6)=-0.325D0  
      GG(8,7)=100.D0/99.D0      
      GG(9,1)=-1.D0     
      GG(9,4)=1.22D0    
      GG(9,5)=-1.D0     
      GG(10,6)=-1.D0    
      GG(11,8)=-1.D0    
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.169809564792D+04     
      XEX(2)=0.158187256985D+05     
      XEX(3)=0.541022827849D+02     
      XEX(4)=0.303122594099D+04     
      XEX(5)=0.200000000000D+04     
      XEX(6)=0.901153669236D+02     
      XEX(7)=0.950000000000D+02     
      XEX(8)=0.104933580864D+02     
      XEX(9)=0.156163636380D+01     
      XEX(10)=0.153535353535D+03     
      FEX=-0.176880696344D+04
      RETURN    
2     FX=5.04D0*X(1)+0.035D0*X(2)+10.D0*X(3)+3.36D0*X(5)
     /       -0.063D0*X(4)*X(7) 
3     GF(4)=-0.063D0*X(7)      
      GF(7)=-0.063D0*X(4)       
      RETURN    
4     IF (INDEX1(1)) G(1)=35.82D0-0.222D0*X(10)-0.9D0*X(9)      
      IF (INDEX1(2)) G(2)=-133.D0+3.D0*X(7)-0.99D0*X(10)
      IF (INDEX1(3)) G(3)=-35.82D0+0.222D0*X(10)+10.D0/9.D0*X(9) 
      IF (INDEX1(4)) G(4)=133.D0-3.D0*X(7)+X(10)/0.99D0 
      IF (INDEX1(5)) G(5)=1.12D0*X(1)+0.13167D0*X(1)*X(8)
     /       -6.67D-3*X(1)*X(8)**2-0.99D0*X(4)    
      IF (INDEX1(6)) G(6)=57.425D0+1.098D0*X(8)-0.038D0*X(8)**2
     /       +0.325D0*X(6)-0.99D0*X(7)      
      IF (INDEX1(7)) G(7)=-1.12D0*X(1)-0.13167D0*X(1)*X(8)+     
     /        6.67D-3*X(1)*X(8)**2+X(4)/0.99D0   
      IF (INDEX1(8)) G(8)=-57.425D0-1.098D0*X(8)+0.038D0*X(8)**2
     /        -0.325D0*X(6)+X(7)/0.99D0       
      IF (INDEX1(9)) G(9)=1.22D0*X(4)-X(1)-X(5) 
      IF (INDEX1(10)) G(10)=9.8D+4*X(3)/(X(4)*X(9)+1.D+3*X(3))-X(6)
      IF (INDEX1(11)) G(11)=(X(2)+X(5))/X(1)-X(8)       
      RETURN    
5     IF (.NOT.INDEX2(5)) GOTO 11       
      GG(5,1)=1.12D0+0.13167D0*X(8)-6.67D-3*X(8)**2     
      GG(5,8)=0.13167D0*X(1)-2.D0*6.67D-3*X(1)*X(8)     
11    IF (INDEX2(6)) GG(6,8)=1.098D0-0.076D0*X(8)       
      IF (.NOT.INDEX2(7)) GOTO 13       
      GG(7,8)=-0.13167D0*X(1)+2.D0*6.67D-3*X(1)*X(8)    
      GG(7,1)=-0.13167D0*X(8)+6.67D-3*X(8)**2-1.12D0    
13    IF (INDEX2(8)) GG(8,8)=-1.098D0+0.076D0*X(8)      
      IF (.NOT.INDEX2(10)) GOTO 16      
      V1=(X(4)*X(9)+1.D+3*X(3))**2      
      V2=9.8D+4*X(9)    
      V3=V2/V1  
      GG(10,3)=X(4)*V3  
      GG(10,4)=-X(3)*V3 
      GG(10,9)=-9.8D+4*X(3)*X(4)/V1     
16    IF (.NOT.INDEX2(11)) GOTO 17      
      GG(11,1)=-(X(2)+X(5))/X(1)**2     
      GG(11,2)=1.D0/X(1)
      GG(11,5)=GG(11,2) 
17    RETURN    
      END       
C
      SUBROUTINE TP116(MODE)    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      GOTO (1,2,3,4,5),MODE     
1     N=13      
      NILI=5    
      NINL=10   
      NELI=0    
      NENL=0    
      X(1)=0.5D0
      X(2)=0.8D0
      X(3)=0.9D0
      X(4)=0.1D0
      X(5)=0.14D0       
      X(6)=0.5D0
      X(7)=489.D0       
      X(8)=80.D0
      X(9)=650.D0       
      X(10)=450.D0      
      X(11)=150.D0      
      X(12)=150.D0      
      X(13)=150.D0      
      DO 6 I=1,10       
6     XL(I)=0.1D0       
      XL(4)=1.D-4       
      XL(9)=500.D0      
      XL(11)=1.D0       
      XL(12)=1.D-4      
      XL(13)=1.D-4      
      DO 7 I=1,3
      XU(I)=1.D0
      XU(I+6)=1.D+3     
7     XU(I+10)=150.D0   
      XU(4)=0.1D0       
      XU(5)=0.9D0       
      XU(6)=0.9D0       
      XU(10)=500.D0     
      DO 32 I=1,13      
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      DO 32 J=1,15      
32    GG(J,I)=0.D0      
      DO 30 I=1,10      
30    GF(I)=0.D0
      DO 31 I=11,13     
31    GF(I)=1.D0
      GG(1,2)=-1.D0     
      GG(1,3)=1.D0      
      GG(2,1)=-1.D0     
      GG(2,2)=1.D0      
      GG(3,7)=-2.D-3    
      GG(3,8)=2.D-3     
      GG(4,11)=1.D0     
      GG(4,12)=1.D0     
      GG(4,13)=1.D0     
      GG(5,11)=-1.D0    
      GG(5,12)=-1.D0    
      GG(5,13)=-1.D0    
      GG(6,13)=1.D0     
      GG(14,11)=1.D0    
      GG(15,12)=1.D0    
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.803770278595D+00     
      XEX(2)=0.899986033670D+00     
      XEX(3)=0.970972419495D+00     
      XEX(4)=0.999995162129D-01     
      XEX(5)=0.190815447786D+00     
      XEX(6)=0.460571745738D+00     
      XEX(7)=0.574080310673D+03     
      XEX(8)=0.740804261398D+02     
      XEX(9)=0.500016155317D+03     
      XEX(10)=0.999999999985D-01     
      XEX(11)=0.202341325935D+02     
      XEX(12)=0.773475459898D+02     
      XEX(13)=0.673039736648D-02     
      FEX=0.975884089805D+02 
      RETURN    
2     FX=X(11)+X(12)+X(13)      
3     RETURN    
4     IF (INDEX1(1)) G(1)=X(3)-X(2)     
      IF (INDEX1(2)) G(2)=X(2)-X(1)     
      IF (INDEX1(3)) G(3)=1.D0-2.D-3*(X(7)-X(8))
      IF (INDEX1(4)) G(4)=X(11)+X(12)+X(13)-50.0D0      
      IF (INDEX1(5)) G(5)=250.D0-X(11)-X(12)-X(13)      
      IF (INDEX1(6)) G(6)=X(13)-1.262626D0*X(10)+1.231059D0     
     /     *X(3)*X(10)       
      IF (INDEX1(7)) G(7)=X(5)-0.03475D0*X(2)-0.975D0*X(2)*     
     /     X(5)+9.75D-3*X(2)**2       
      IF (INDEX1(8)) G(8)=X(6)-0.03475D0*X(3)-0.975D0*X(3)
     /     *X(6)+9.75D-3*X(3)**2       
      IF (INDEX1(9)) G(9)=X(5)*X(7)-X(1)*X(8)-X(4)*X(7) 
     /     +X(4)*X(8)
      IF (INDEX1(10)) G(10)=-2.D-3*(X(2)*X(9)+X(5)*X(8) 
     /     -X(1)*X(8)-X(6)*X(9))-X(6)-X(5)+1.D0      
      IF (INDEX1(11)) G(11)=X(2)*X(9)-X(3)*X(10)-X(6)*X(9)
     /     -500.D0*(X(2)-X(6))+X(2)*X(10)    
      IF (INDEX1(12)) G(12)=X(2)-0.9D0-2.D-3*(X(2)*X(10)-       
     /      X(3)*X(10))       
      IF (INDEX1(13)) G(13)=X(4)-0.03475D0*X(1)-0.975D0*X(1)*X(4)
     /    +9.75D-3*X(1)**2   
      IF(INDEX1(14)) G(14)=X(11)-1.262626D0*X(8)+1.231059D0*X(1)*X(8)
      IF (INDEX1(15)) G(15)=X(12)-1.262626D0*X(9)+1.231059D0*X(2)*X(9)
      RETURN    
5     IF (.NOT.INDEX2(6)) GOTO 12       
      GG(6,3)=1.231059D0*X(10)  
      GG(6,10)=-1.262626D0+1.231059D0*X(3)      
12    IF (.NOT.INDEX2(7)) GOTO 13       
      GG(7,2)=-0.03475D0-0.975D0*X(5)+1.95D-2*X(2)      
      GG(7,5)=1.D0-0.975D0*X(2) 
13    IF (.NOT.INDEX2(8)) GOTO 14       
      GG(8,3)=-0.03475D0-0.975D0*X(6)+1.95D-2*X(3)      
      GG(8,6)=1.D0-0.975D0*X(3) 
14    IF (.NOT.INDEX2(9)) GOTO 15       
      GG(9,1)=-X(8)     
      GG(9,4)=-X(7)+X(8)
      GG(9,5)=X(7)      
      GG(9,7)=X(5)-X(4) 
      GG(9,8)=-X(1)+X(4)
15    IF (.NOT.INDEX2(10)) GOTO 16      
      GG(10,1)=2.D-3*X(8)       
      GG(10,2)=-2.D-3*X(9)      
      GG(10,5)=-2.D-3*X(8)-1.D0 
      GG(10,6)=-1.D0+2.D-3*X(9) 
      GG(10,8)=-2.D-3*(X(5)-X(1))       
      GG(10,9)=-2.D-3*(X(2)-X(6))       
16    IF (.NOT.INDEX2(11)) GOTO 17      
      GG(11,2)=X(9)-500.D0+X(10)
      GG(11,3)=-X(10)   
      GG(11,6)=-X(9)+500.D0     
      GG(11,9)=X(2)-X(6)
      GG(11,10)=-X(3)+X(2)      
17    IF (.NOT.INDEX2(12)) GOTO 18      
      GG(12,2)=1.D0-2.D-3*X(10) 
      GG(12,3)=2.D-3*X(10)      
      GG(12,10)=2.D-3*(X(3)-X(2))       
18    IF (.NOT.INDEX2(13)) GOTO 19      
      GG(13,1)=-0.03475D0-0.975D0*X(4)+1.95D-2*X(1)     
      GG(13,4)=1.D0-0.975D0*X(1)
19    IF (.NOT.INDEX2(14)) GOTO 20      
      GG(14,1)=1.231059D0*X(8)  
      GG(14,8)=1.231059D0*X(1)-1.262626D0       
20    IF (.NOT.INDEX2(15)) GOTO 21      
      GG(15,2)=1.231059D0*X(9)  
      GG(15,9)=-1.262626D0+1.231059D0*X(2)      
21    RETURN    
      END       
C
      SUBROUTINE TP117(MODE)    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION E,D,B,C,A,T1,T2,T3,T4,T5,T6
      COMMON     /D117/E(5),D(5),B(10),C(5,5),A(10,5),T4(5),T5(5),T6(5) 
      INTEGER I,J
      GOTO (1,2,3,2,5),MODE     
1     N=15      
      NILI=0    
      NINL=5    
      NELI=0    
      NENL=0    
      DO 20 I=1,15      
20    X(I)=1.D-3
      X(7)=60.D0
      DO 22 I=1,15      
      XL(I)=0.D0
      LXL(I)=.TRUE.     
      LXU(I)=.FALSE.    
      DO 22 J=1,5       
22    GG(J,I)=0.D0      
      E(1)=-15.D0       
      E(2)=-27.D0       
      E(3)=-36.D0       
      E(4)=-18.D0       
      E(5)=-12.D0       
      C(1,1)=30.D0      
      C(1,2)=-20.D0     
      C(1,3)=-10.D0     
      C(1,4)=32.D0      
      C(1,5)=-10.D0     
      C(2,2)=39.D0      
      C(2,3)=-6.D0      
      C(2,4)=-31.D0     
      C(2,5)=32.D0      
      C(3,3)=10.D0      
      C(3,4)=-6.D0      
      C(3,5)=-10.D0     
      C(4,4)=39.D0      
      C(4,5)=-20.D0     
      C(5,5)=30.D0      
      DO 70 I=1,5       
      DO 70 J=1,5       
70    C(J,I)=C(I,J)     
      D(1)=4.D0 
      D(2)=8.D0 
      D(3)=10.D0
      D(4)=6.D0 
      D(5)=2.D0 
      DO 72 I=1,6       
      DO 72 J=1,5       
72    A(I,J)=0.D0       
      A(1,1)=-16.D0     
      A(1,2)=2.D0       
      A(1,4)=1.D0       
      A(2,2)=-2.D0      
      A(2,4)=0.4D0      
      A(2,5)=2.D0       
      A(3,1)=-3.5D0     
      A(3,3)=2.D0       
      A(4,2)=-2.D0      
      A(4,4)=-4.D0      
      A(4,5)=-1.D0      
      A(5,2)=-9.D0      
      A(5,3)=-2.D0      
      A(5,4)=1.D0       
      A(5,5)=-2.8D0     
      A(6,1)=2.D0       
      A(6,3)=-4.D0      
      A(8,1)=-1.D0      
      A(8,2)=-2.D0      
      A(8,3)=-3.D0      
      A(8,4)=-2.D0      
      A(8,5)=-1.D0      
      DO 73 I=1,5       
      A(7,I)=-1.D0      
      A(9,I)=DBLE(I)  
73    A(10,I)=1.D0      
      B(1)=-40.D0       
      B(2)=-2.D0
      B(3)=-0.25D0      
      B(4)=-4.D0
      B(5)=-4.D0
      B(6)=-1.D0
      B(7)=-40.D0       
      B(8)=-60.D0       
      B(9)=5.D0 
      B(10)=1.D0
      DO 35 I=1,10      
35    GF(I)=-B(I)       
      DO 40 I=1,10      
      DO 40 J=1,5       
40    GG(J,I)=-A(I,J)   
      DO 41 I=1,5       
      DO 42 J=1,5       
      IF (I-J) 45,42,45 
45    GG(J,I+10)=2.D0*C(I,J)    
42    CONTINUE  
41    CONTINUE  
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.0D0  
      XEX(2)=0.0D0  
      XEX(3)=0.517413630680D+01     
      XEX(4)=0.0D0  
      XEX(5)=0.306109271525D+01     
      XEX(6)=0.118396760290D+02     
      XEX(7)=0.0D0  
      XEX(8)=0.0D0  
      XEX(9)=0.103907059194D+00     
      XEX(10)=0.0D0  
      XEX(11)=0.299992902601D+00     
      XEX(12)=0.333470928832D+00     
      XEX(13)=0.399990975915D+00     
      XEX(14)=0.428314541579D+00     
      XEX(15)=0.223960749729D+00     
      FEX=0.323486789791D+02 
      RETURN    
2     T1=0.D0   
      T2=0.D0   
      DO 30 J=1,5       
      T2=T2+D(J)*X(10+J)**3     
      DO 30 I=1,5       
      T1=T1+C(I,J)*X(10+I)*X(10+J)      
30    CONTINUE  
      T3=0.D0   
      DO 31 I=1,10      
31    T3=T3+B(I)*X(I)   
      DO 34 J=1,5       
      T4(J)=0.D0
      T5(J)=0.D0
      DO 32 I=1,5       
32    T4(J)=T4(J)+C(I,J)*X(10+I)
      DO 33 I=1,10      
33    T5(J)=T5(J)+A(I,J)*X(I)   
34    CONTINUE  
      IF (MODE.EQ.4) GOTO 4     
      FX=-(T3-T1-2.D0*T2)       
      RETURN    
3     DO 37 I=1,5       
      T6(I)=0.D0
      DO 37 J=1,5       
37    T6(I)=T6(I)+(C(I,J)+C(J,I))*X(10+J)       
      DO 36 I=1,5       
36    GF(10+I)=T6(I)+6.D0*D(I)*X(10+I)**2       
      RETURN    
4     DO 38 J=1,5       
      IF (INDEX1(J)) G(J)=2.D0*T4(J)+3.D0*D(J)*X(10+J)**2+E     
     -(J)-T5(J) 
38    CONTINUE  
      RETURN    
5     DO 39 J=1,5       
      IF (.NOT.INDEX2(J)) GOTO 39       
      GG(J,10+J)=2.D0*C(J,J)+6.D0*D(J)*X(10+J)  
39    CONTINUE  
      RETURN    
      END       
      SUBROUTINE TP118(MODE)    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J,K,M
      DOUBLEPRECISION T  
      GOTO (1,2,3,4,5),MODE     
1     N=15      
      NILI=29   
      NINL=0    
      NELI=0    
      NENL=0    
      DO 6 I=1,15       
6     X(I)=20.D0
      X(2)=55.D0
      X(3)=15.D0
      X(5)=60.D0
      X(8)=60.D0
      X(11)=60.D0       
      X(14)=60.D0       
      XL(1)=8.D0
      XL(2)=43.D0       
      XL(3)=3.D0
      XU(1)=21.D0       
      XU(2)=57.D0       
      XU(3)=16.D0       
      DO 22 I=1,4       
      XL(3*I+1)=0.D0    
      XL(3*I+2)=0.D0    
      XL(3*I+3)=0.D0    
      XU(3*I+1)=90.D0   
      XU(3*I+2)=120.D0  
      XU(3*I+3)=60.D0   
22    CONTINUE  
      DO 25 I=1,15      
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      DO 25 J=1,29      
25    GG(J,I)=0.D0      
      DO 20 K=1,4       
      DO 20 I=1,3       
      GG(K+4*I-4,3*K+I)=1.D0    
      GG(K+4*I-4,3*K+I-3)=-1.D0 
      GG(K+8+4*I,3*K+I)=-1.D0   
      GG(K+8+4*I,3*K+I-3)=1.D0  
20    CONTINUE  
      DO 21 K=1,5       
      DO 21 I=1,3       
21    GG(24+K,3*K-3+I)=1.D0     
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.800000000000D+01     
      XEX(2)=0.490000000000D+02     
      XEX(3)=0.300000000000D+01     
      XEX(4)=0.100000000000D+01     
      XEX(5)=0.560000000000D+02     
      XEX(6)=0.0D0  
      XEX(7)=0.999999999545D+00     
      XEX(8)=0.630000000000D+02     
      XEX(9)=0.600000000000D+01     
      XEX(10)=0.299999999965D+01     
      XEX(11)=0.700000000000D+02     
      XEX(12)=0.120000000000D+02     
      XEX(13)=0.499999999971D+01     
      XEX(14)=0.770000000000D+02     
      XEX(15)=0.180000000000D+02     
      FEX=0.664820449993D+03 
      RETURN    
2     T=0.D0    
      DO 30 M=1,5       
      I=M-1     
30    T=T+2.3D0*X(3*I+1)+1.D-4*X(3*I+1)**2+1.7D0*X(3*I+2)+1.D-4
     /    *X(3*I+2)**2+2.2D0*X(3*I+3)+1.5D-4*X(3*I+3)**2       
      FX=T      
      RETURN    
3     DO 31 I=1,5       
      GF(3*I-2)=2.3D0+2.D-4*X(3*I-2)    
      GF(3*I-1)=1.7D0+2.D-4*X(3*I-1)    
31    GF(3*I)=2.2D0+3.D-4*X(3*I)
      RETURN    
4     DO 32 I=1,4       
      IF (INDEX1(I)) G(I)=X(3*I+1)-X(3*I-2)+7.D0
      IF (INDEX1(I+4)) G(I+4)=X(3*I+2)-X(3*I-1)+7.D0    
      IF (INDEX1(I+8)) G(I+8)=X(3*I+3)-X(3*I)+7.D0      
      IF (INDEX1(I+12)) G(I+12)=X(3*I-2)-X(3*I+1)+6.D0  
      IF (INDEX1(I+16)) G(I+16)=X(3*I-1)-X(3*I+2)+7.D0  
32    IF (INDEX1(I+20)) G(I+20)=X(3*I)-X(3*I+3)+6.D0    
      IF (INDEX1(25)) G(25)=X(1)+X(2)+X(3)-60.D0
      IF (INDEX1(26)) G(26)=X(4)+X(5)+X(6)-50.D0
      IF (INDEX1(27)) G(27)=X(7)+X(8)+X(9)-70.D0
      IF (INDEX1(28)) G(28)=X(10)+X(11)+X(12)-85.D0     
      IF (INDEX1(29)) G(29)=X(13)+X(14)+X(15)-100.D0    
5     RETURN    
      END       
      SUBROUTINE TP119(MODE)    
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D119/A(16,16),B(8,16),C(8),S(16),T
      DOUBLEPRECISION  A,B,C,S,T  
      INTEGER I,J
      GOTO (1,2,3,4,5),MODE     
1     N=16      
      NILI=0    
      NINL=0    
      NELI=8    
      NENL=0    
      DO 80 I=1,16      
      DO 21 J=1,16      
21    A(I,J)=0.D0       
      A(I,I)=1.D0       
      DO 81 J=1,8       
81    B(J,I)=0.D0       
80    CONTINUE  
      DO 22 I=1,8       
22    B(I,8+I)=1.D0     
      A(1,4)=1.D0       
      A(1,7)=1.D0       
      A(1,8)=1.D0       
      A(1,16)=1.D0      
      A(2,3)=1.D0       
      A(2,7)=1.D0       
      A(2,10)=1.D0      
      A(3,7)=1.D0       
      A(3,9)=1.D0       
      A(3,10)=1.D0      
      A(3,14)=1.D0      
      A(4,7)=1.D0       
      A(4,11)=1.D0      
      A(4,15)=1.D0      
      A(5,6)=1.D0       
      A(5,10)=1.D0      
      A(5,12)=1.D0      
      A(5,16)=1.D0      
      A(6,8)=1.D0       
      A(6,15)=1.D0      
      A(7,11)=1.D0      
      A(7,13)=1.D0      
      A(8,10)=1.D0      
      A(8,15)=1.D0      
      A(9,12)=1.D0      
      A(9,16)=1.D0      
      A(10,14)=1.D0     
      A(11,13)=1.D0     
      A(12,14)=1.D0     
      A(13,14)=1.D0     
      B(1,1)=0.22D0     
      B(1,2)=0.2D0      
      B(1,3)=0.19D0     
      B(1,4)=0.25D0     
      B(1,5)=0.15D0     
      B(1,6)=0.11D0     
      B(1,7)=0.12D0     
      B(1,8)=0.13D0     
      B(2,1)=-1.46D0    
      B(2,3)=-1.3D0     
      B(2,4)=1.82D0     
      B(2,5)=-1.15D0    
      B(2,7)=0.8D0      
      B(3,1)=1.29D0     
      B(3,2)=-0.89D0    
      B(3,5)=-1.16D0    
      B(3,6)=-0.96D0    
      B(3,8)=-0.49D0    
      B(4,1)=-1.1D0     
      B(4,2)=-1.06D0    
      B(4,3)=0.95D0     
      B(4,4)=-0.54D0    
      B(4,6)=-1.78D0    
      B(4,7)=-0.41D0    
      B(5,4)=-1.43D0    
      B(5,5)=1.51D0     
      B(5,6)=0.59D0     
      B(5,7)=-0.33D0    
      B(5,8)=-0.43D0    
      B(6,2)=-1.72D0    
      B(6,3)=-0.33D0    
      B(6,5)=1.62D0     
      B(6,6)=1.24D0     
      B(6,7)=0.21D0     
      B(6,8)=-0.26D0    
      B(7,1)=1.12D0     
      B(7,4)=0.31D0     
      B(7,7)=1.12D0     
      B(7,9)=-0.36D0    
      B(8,2)=0.45D0     
      B(8,3)=0.26D0     
      B(8,4)=-1.1D0     
      B(8,5)=0.58D0     
      B(8,7)=-1.03D0    
      B(8,8)=0.1D0      
      C(1)=2.5D0
      C(2)=1.1D0
      C(3)=-3.1D0       
      C(4)=-3.5D0       
      C(5)=1.3D0
      C(6)=2.1D0
      C(7)=2.3D0
      C(8)=-1.5D0       
      DO 20 I=1,16      
      LXL(I)=.TRUE.     
      LXU(I)=.TRUE.     
      X(I)=10.D0
      XL(I)=0.D0
      XU(I)=5.D0
      DO 20 J=1,8       
20    GG(J,I)=B(J,I)    
      LEX=.FALSE.       
      NEX=1
      XEX(1)=0.398473514099D-01     
      XEX(2)=0.791983155694D+00     
      XEX(3)=0.202870330224D+00     
      XEX(4)=0.844357916347D+00     
      XEX(5)=0.126990645286D+01     
      XEX(6)=0.934738707827D+00     
      XEX(7)=0.168196196924D+01     
      XEX(8)=0.155300877490D+00     
      XEX(9)=0.156787033356D+01     
      XEX(10)=-0.359021173251D-11     
      XEX(11)=-0.612900888082D-11     
      XEX(12)=-0.886794857449D-12     
      XEX(13)=0.660204066000D+00     
      XEX(14)=-0.254340725727D-11     
      XEX(15)=0.674255926901D+00     
      XEX(16)=-0.110433723798D-10     
      FEX=0.244899697515D+03 
      RETURN    
2     T=0.D0    
      DO 30 I=1,16      
      DO 30 J=1,16      
30    T=T+A(I,J)*(X(I)**2+X(I)+1.D0)*(X(J)**2+X(J)+1.D0)
      FX=T      
      RETURN    
3     DO 31 I=1,16       
      S(I)=0.D0 
      DO 32 J=1,16      
32    S(I)=S(I)+(A(I,J)+A(J,I))*(X(J)**2+X(J)+1.D0)*(2.D0*X(I)+1.D0) 
31    GF(I)=S(I)
      RETURN    
4     DO 33 I=1,8       
      IF (.NOT.INDEX1(I)) GOTO 33       
      S(I)=0.D0 
      DO 34 J=1,16      
34    S(I)=S(I)+B(I,J)*X(J)     
      G(I)=S(I)-C(I)    
33    CONTINUE  
5     RETURN    
      END       
C
      SUBROUTINE TP201(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0 
      NENL=0
      X(1)=8.D+0
      X(2)=9.D+0
      DO  6 I=1,2 
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      XEX(1)=5.D+0
      XEX(2)=6.D+0
      RETURN 
    2 FX=4.D+0*(X(1)-5.D+0)**2 + (X(2)-6.D+0)**2
      RETURN
    3 GF(1)=8.D+0*(X(1)-5.D+0)
      GF(2)=2.D+0*(X(2)-6.D+0)
    4 RETURN
      END
C
      SUBROUTINE TP202(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF
      INTEGER I
      GOTO(1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0 
      NENL=0
      X(1)=15.D+0
      X(2)=-2.D+0
      DO 6 I=1,2 
      LXU(I)=.TRUE.
    6 LXL(I)=.TRUE.
      XU(1)=20.0
      XU(2)=5.0
      XL(1)=1.0
      XL(2)=-5.0
      LEX=.TRUE.
      NEX=2
      FEX=0.D+0
      XEX(1)=5.D+0
      XEX(2)=4.D+0
      XEX(3)=48.98425D+0
      XEX(4)=-0.89681D+0
      RETURN 
    2 F(1)=-13.D+0+X(1)-2.D+0*X(2)+5.D+0*X(2)**2-X(2)**3
      F(2)=-29.D+0+X(1)-14.D+0*X(2)+X(2)**2+X(2)**3
      FX=F(1)**2+F(2)**2
      RETURN
    3 F(1)=-13.D+0+X(1)-2.D+0*X(2)+5.D+0*X(2)**2-X(2)**3
      F(2)=-29.D+0+X(1)-14.D+0*X(2)+X(2)**2+X(2)**3     
      DF(1,1)=1.D+0
      DF(1,2)=-2.D+0+10.D+0*X(2)-3.D+0*X(2)**2
      DF(2,1)=1.D+0
      DF(2,2)=-14.D+0+2.D+0*X(2)+3.D+0*X(2)**2
      DO 7 I=1,2
    7 GF(I)=2.D+0*F(1)*DF(1,I)+2.D+0*F(2)*DF(2,I)
    4 RETURN
      END
C
      SUBROUTINE TP203(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)   
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)     
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION C(3),F,DF
      DATA C/1.5D+0,2.25D+0,2.625D+0/
      INTEGER I
      GOTO (1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0 
      NENL=0
      X(1)=2.D+0
      X(2)=0.2D+0
      DO 6 I=1,2 
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      XEX(1)=3.D+0
      XEX(2)=0.5D+0
      RETURN 
    2 FX=0.D+0
      DO 7 I=1,3
      F(I)=C(I)-X(1)*(1.D+0-X(2)**I) 
    7 FX=FX+F(I)**2
      RETURN
    3 DO 8 I=1,3
    8 F(I)=C(I)-X(1)*(1.D+0-X(2)**I) 
      DF(1,1)=-1.D+0+X(2)
      DF(1,2)=X(1)
      DF(2,1)=-1.D+0+X(2)**2
      DF(2,2)=2.D+0*X(1)*X(2)
      DF(3,1)=-1.D+0+X(2)**3
      DF(3,2)=3.D+0*X(1)*X(2)**2
      DO 9 I=1,2
    9 GF(I)=2.D+0*F(1)*DF(1,I)+2.D+0*F(2)*DF(2,I)+2.D+0*F(3)*DF(3,I)
    4 RETURN
      END
C
      SUBROUTINE TP204(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION PROD,A(3),D(3),H(3,2),F,DF
      DATA A/0.13294D+0,-0.244378D+0,0.325895D+0/
      DATA D/2.5074D+0,-1.36401D+0,1.02282D+0/
      DATA H/-0.564255D+0,-0.404979D+0,-0.0735084D+0,0.392417D+0,
     1       0.927589D+0,0.535493D+0/
C      DATA B/5.66598D+0,2.77141D+0,2.77141D+01,2.12413D+0/
      INTEGER I
      GOTO(1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0 
      NENL=0
      X(1)=0.1D+0
      X(2)=0.1D+0
      DO 6 I=1,2 
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.183601D+0
      XEX(1)=0.D+0
      XEX(2)=0.D+0
      RETURN
    2 DO 10 I=1,3
      PROD=H(I,1)*X(1)+H(I,2)*X(2)
   10 F(I)=A(I)+PROD+.5D+0*PROD**2*D(I)
      FX=F(1)**2+F(2)**2+F(3)**2
      RETURN
    3 DO 11 I=1,3
      PROD=H(I,1)*X(1)+H(I,2)*X(2)
   11 F(I)=A(I)+PROD+.5D+0*PROD**2*D(I)
      DO 7 I=1,3 
    7 DF(I,1)=H(I,1)+(H(I,1)*X(1)+H(I,2)*X(2))*H(I,1)*D(I)
      DO 8 I=1,3
    8 DF(I,2)=H(I,2)+(H(I,1)*X(1)+H(I,2)*X(2))*H(I,2)*D(I)
      DO 9 I=1,2
    9 GF(I)=2.D+0*(F(1)*DF(1,I)+F(2)*DF(2,I)+F(3)*DF(3,I))
    4 RETURN
      END
C
      SUBROUTINE TP205(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF
      INTEGER I
      GOTO(1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0 
      NENL=0
      X(1)=0.D+0
      X(2)=0.D+0
      DO 6 I=1,2 
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      XEX(1)=3.D+0
      XEX(2)=0.5D+0 
      RETURN 
    2 F(1)=1.5D+0-X(1)*(1.D+0-X(2))
      F(2)=2.25D+0-X(1)*(1.D+0-X(2)**2)
      F(3)=2.625D+0-X(1)*(1.D+0-X(2)**3) 
      FX=F(1)**2+F(2)**2+F(3)**2
      RETURN
    3 F(1)=1.5D+0-X(1)*(1.D+0-X(2))
      F(2)=2.25D+0-X(1)*(1.D+0-X(2)**2)
      F(3)=2.625D+0-X(1)*(1.D+0-X(2)**3)
      DF(1,1)=X(2)-1.D+0
      DF(1,2)=X(1)
      DF(2,1)=X(2)**2-1.D+0
      DF(2,2)=2.D+0*X(1)*X(2)
      DF(3,1)=X(2)**3-1.D+0
      DF(3,2)=3.D+0*X(1)*X(2)**2
      DO 7 I=1,2
    7 GF(I)=2.D+0*F(1)*DF(1,I)+2.D+0*F(2)*DF(2,I)+2.D+0*F(3)*DF(3,I)
    4 RETURN
      END
C
      SUBROUTINE TP206(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0 
      NENL=0
      X(1)=-1.2D+0
      X(2)=1.D+0
      DO 6 I=1,2 
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      XEX(1)=1.D+0
      XEX(2)=1.D+0
      RETURN 
    2 FX=(X(2)-X(1)**2)**2+100.D+0*(1.D+0-X(1))**2
      RETURN
    3 GF(1)=-2.D+0*(X(2)-X(1)**2)*2.D+0*X(1)-200.D+0*(1.D+0-X(1))
      GF(2)=2.D+0*(X(2)-X(1)**2)
    4 RETURN
      END
C
      SUBROUTINE TP207(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION C
      DATA C/1.D+0/
      GOTO(1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0 
      NENL=0
      X(1)=-1.2D+0
      X(2)=1.D+0
      DO 6 I=1,2 
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      XEX(1)=1.D+0
      XEX(2)=1.D+0
      RETURN 
    2 FX=C*(X(2)-X(1)**2)**2+(1.D+0-X(1))**2
      RETURN
    3 GF(1)=-4.D+0*C*(X(2)-X(1)**2)*X(1)-2.D+0*(1.D+0-X(1))
      GF(2)=2.D+0*C*(X(2)-X(1)**2)
    4 RETURN
      END
C
      SUBROUTINE TP208(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION C
      DATA C/100.D+0/
      INTEGER I
      GOTO(1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0 
      NENL=0
      X(1)=-1.2D+0
      X(2)=1.D+0
      DO 6 I=1,2 
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      XEX(1)=1.D+0
      XEX(2)=1.D+0
      RETURN 
    2 FX=C*(X(2)-X(1)**2)**2+(1.D+0-X(1))**2
      RETURN
    3 GF(1)=-4.D+0*C*(X(2)-X(1)**2)*X(1)-2.D+0*(1.D+0-X(1))
      GF(2)=2.D+0*C*(X(2)-X(1)**2)
    4 RETURN
      END
      SUBROUTINE TP209(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION C
      DATA C/10000.D+0/
      INTEGER I
      GOTO(1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0 
      NENL=0
      X(1)=-1.2D+0
      X(2)=1.D+0
      DO 6 I=1,2 
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      XEX(1)=1.D+0
      XEX(2)=1.D+0
      RETURN 
    2 FX=C*(X(2)-X(1)**2)**2+(1.D+0-X(1))**2
      RETURN
    3 GF(1)=-4.D+0*C*(X(2)-X(1)**2)*X(1)-2.D+0*(1.D+0-X(1))
      GF(2)=2.D+0*C*(X(2)-X(1)**2)
    4 RETURN
      END
      SUBROUTINE TP210(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION C
      DATA C/1000000.D+0/
      INTEGER I
      GOTO(1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0 
      NENL=0
      X(1)=-1.2D+0
      X(2)=1.D+0
      DO 6 I=1,2 
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      XEX(1)=1.D+0
      XEX(2)=1.D+0
      RETURN 
    2 FX=(C*(X(2)-X(1)**2)**2+(1.D+0-X(1))**2)/C
      RETURN
    3 GF(1)=(-4.D+0*C*(X(2)-X(1)**2)*X(1)-2.D+0*(1.D+0-X(1)))/C
      GF(2)=(2.D+0*C*(X(2)-X(1)**2))/C
    4 RETURN
      END
      SUBROUTINE TP211(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0 
      NENL=0
      X(1)=-1.2D+0
      X(2)=1.D+0
      DO 6 I=1,2 
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      XEX(1)=1.D+0
      XEX(2)=1.D+0
      RETURN 
    2 FX=100.D+0*(X(2)-X(1)**3)**2+(1.D+0-X(1))**2
      RETURN
    3 GF(1) =-200.D+0*(X(2)-X(1)**3)*3.D+0*X(1)**2-2.D+0*(1.D+0-X(1))
      GF(2)=200.D+0*(X(2)-X(1)**3)
    4 RETURN
      END
      SUBROUTINE TP212(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0 
      NENL=0
      X(1)=2.D+0
      X(2)=0.D+0
      DO 6 I=1,2 
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      XEX(1)=0.D+0
      XEX(2)=0.D+0
      RETURN 
    2 FX=(4.D+0*(X(1)+X(2)))**2+(4.D+0*(X(1)+X(2))+(X(1)-X(2))
     /  *((X(1)-2.D+0)**2+X(2)**2-1.D+0))**2
      RETURN
    3 GF(1)=32.D+0*(X(1)+X(2))+2.D+0*(4.D+0*(X(1)+X(2))+(X(1)-X(2))
     /    *((X(1)-2.D+0)**2+X(2)**2-1.D+0))*(4.D+0+((X(1)-2.D+0)**2
     /    +X(2)**2-1.D+0)+(X(1)-X(2))*2.D+0*(X(1)-2.D+0))      
      GF(2)=32.D+0*(X(1)+X(2))+2.D+0*(4.D+0*(X(1)+X(2))+(X(1)-X(2))
     /    *((X(1)-2.D+0)**2+X(2)**2-1.D+0))*(4.D+0-(X(1)-2.D+0)**2
     /    +X(2)**2-1.D+0+(X(1)-X(2))*2.D+0*X(2))      
    4 RETURN
      END
C
      SUBROUTINE TP213(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0 
      NENL=0
      X(1)=3.D+0
      X(2)=1.D+0
      DO 6 I=1,2 
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      XEX(1)=1.D+0
      XEX(2)=1.D+0
      RETURN 
    2 FX=(10.D+0*(X(1)-X(2))**2+(X(1)-1.D+0)**2)**4
      RETURN
    3 GF(1)=(4.D+0*(10.D+0*(X(1)-X(2))**2+(X(1)-1.D+0)**2)**3
     /      *(20.D+0*(X(1)-X(2))+2.D+0*(X(1)-1.D+0)))
      GF(2)=(4.D+0*(10.D+0*(X(1)-X(2))**2+(X(1)-1.D+0)**2)**3*20.D+0
     /      *(X(2)-X(1)))
    4 RETURN
      END
C
      SUBROUTINE TP214(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0 
      NENL=0
      X(1)=-1.2D+0
      X(2)=1.0D+0
      DO 6 I=1,2 
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      XEX(1)=1.D+0
      XEX(2)=1.D+0
      RETURN 
    2 FX=(10.0D0*(X(1)-X(2))**2 + (X(1)-1.0D0)**2)**0.25D0
      RETURN
    3 GF(1)=((0.25D+0/(10.0D0*(X(1)-X(2))**2+(X(1)-1.0D0)**2)**0.75D0)
     /           *(22.0D0*X(1)-20.D0*X(2)-2.0D0))
      GF(2)=((0.25D+0/(10.0D0*(X(1)-X(2))**2+(X(1)-1.0D0)**2)**0.75D0)
     /           *20.0D0*(X(2)-X(1)))
    4 RETURN
      END
C
      SUBROUTINE TP215(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=1
      NELI=0 
      NENL=0
      X(1)=1.D+0
      X(2)=1.D+0
      LXU(1)=.FALSE.
      LXU(2)=.FALSE.
      LXL(1)=.TRUE.
      LXL(2)=.FALSE.
      XL(1)=0.D+0
      GG(1,2)=1.D+0 
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      XEX(1)=0.D+0
      XEX(2)=0.D+0
      RETURN 
    2 FX=X(2)
      RETURN
    3 GF(1)=0.D+0
      GF(2)=1.D+0
      RETURN
    4 IF (INDEX1(1)) G(1)=X(2)-X(1)**2
      RETURN
    5 IF (INDEX2(1)) GG(1,1)=-2.D+0*X(1)
      RETURN
      END
C
      SUBROUTINE TP216(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0 
      NENL=1
      X(1)=-1.2D+0
      X(2)=1.D+0
      DO 6 I=1,2 
      XU(I)=10.0
      LXU(I)=.TRUE.
      XL(I)=-3.0
    6 LXL(I)=.TRUE.
      GG(1,2)=-2.D+0
      LEX=.TRUE.
      NEX=1
      FEX=1.0D+0
      XEX(1)=2.0D+01
      XEX(2)=4.0D+0
      RETURN 
    2 FX=100.D+0*(X(1)**2-X(2))**2+(X(1)-1.D+0)**2
      RETURN
    3 GF(1)=400.D+0*(X(1)**2-X(2))*X(1)+2.D+0*(X(1)-1.D+0)
      GF(2)=-200.D+0*(X(1)**2-X(2))
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)*(X(1)-4.D+0)-2.D+0*X(2)+12.D+0
      RETURN
    5 IF (INDEX2(1)) GG(1,1)=2.D+0*X(1)-4.D+0
      RETURN
      END
C
      SUBROUTINE TP217(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=1
      NINL=0
      NELI=0 
      NENL=1
      X(1)=10.D+0
      X(2)=10.D+0
      LXU(1)=.FALSE.
      LXU(2)=.FALSE.
      LXL(1)=.TRUE.
      LXL(2)=.FALSE.
      XL(1)=0.D+0
      GG(1,1)= 1.D+0
      GG(1,2)= -2.D+0 
      LEX=.TRUE.
      NEX=1
      FEX=-0.8D+0
      XEX(1)=0.6D+0
      XEX(2)=0.8D+0
      RETURN 
    2 FX=-X(2)
      RETURN
    3 GF(1)=0.D+0
      GF(2)=-1.D+0
      RETURN
    4 IF (INDEX1(1)) G(1)=1.D+0+X(1)-2.D+0*X(2)
      IF (INDEX1(2)) G(2)=X(1)**2+X(2)**2-1.D+0 
      RETURN
    5 IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,1)=2.D+0*X(1)
      GG(2,2)=2.D+0*X(2)   
    8 RETURN
      END
C
      SUBROUTINE TP218(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=1
      NELI=0 
      NENL=0
      X(1)=9.D+0
      X(2)=100.D+0
      LXU(1)=.FALSE.
      LXU(2)=.FALSE.
      LXL(1)=.FALSE.
      LXL(2)=.TRUE.
      XL(2)=0
      GG(1,2)=1.D+0 
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      XEX(1)=0.D+0
      XEX(2)=0.D+0
      RETURN 
    2 FX=X(2)
      RETURN
    3 GF(1)=0.D+0
      GF(2)=1.D+0
      RETURN
    4 IF (INDEX1(1)) G(1)=X(2)-X(1)**2
      RETURN
    5 IF (INDEX2(1)) GG(1,1)=-2.D+0*X(1)
      RETURN
      END
C
      SUBROUTINE TP219(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=4
      NILI=0
      NINL=0
      NELI=0 
      NENL=2
      DO 8 I=1,4
      X(I)=10.D+0
      LXU(I)=.FALSE.
    8 LXL(I)=.FALSE.
      GG(1,2)=1.D+0
      GG(1,4)=0.D+0
      GG(2,2)=-1.D+0
      GG(2,3)=0.D+0  
      LEX=.TRUE.
      NEX=1
      FEX=-1.D+0
      XEX(1)=1.D+0
      XEX(2)=1.D+0
      XEX(3)=0.D+0
      XEX(4)=0.D+0
      RETURN 
    2 FX=-X(1)
      RETURN
    3 GF(1)=-1.D+0
      GF(2)=0.D+0
      GF(3)=0.D+0
      GF(4)=0.D+0
      RETURN
    4 IF (INDEX1(1)) G(1)=X(2)-X(1)**3-X(3)**2
      IF (INDEX1(2)) G(2)=X(1)**2-X(2)-X(4)**2
      RETURN
    5 IF (.NOT.INDEX2(1)) GOTO 6
      GG(1,1)=-3.D+0*X(1)**2
      GG(1,3)=-2.D+0*X(3)
    6 IF (.NOT.INDEX2(2)) GOTO 7
      GG(2,1)=2.D+0*X(1)
      GG(2,4)=-2.D+0*X(4)  
    7 RETURN
      END
C
      SUBROUTINE TP220(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=1
      NELI=0
      NENL=0
      DO 6 I=1,2
      X(I)=.25D+5
      LXU(I)=.FALSE.
    6 LXL(I)=.TRUE.
      XL(1)=1.0D0
      XL(2)=0.0D0
      LEX=.TRUE.
      NEX=1
      FEX=1.D+0
      XEX(1)=1.0D0
      XEX(2)=0.0D0
      RETURN
    2 FX=X(1)
      RETURN
    3 GF(1)=1.0D0
      GF(2)=0.0D0
      RETURN
    4 IF (INDEX1(1)) G(1)=(X(1)-1.0D0)**3-X(2)
      RETURN
    5 IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=3.0D0*(X(1)-1.0D0)**2
      GG(1,2)=-1.0D0
    7 RETURN
      END
C
      SUBROUTINE TP221(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=1
      NELI=0
      NENL=0
      DO 6 I=1,2
      X(I)=0.25D+0
      LXU(I)=.TRUE.
      XU(I)=1.0
      LXL(I)=.TRUE.
    6 XL(I)=0.0
      GG(1,2)=-1.0
      LEX=.TRUE.
      NEX=1
      FEX=-1.D+0
      XEX(1)=1.D+0
      XEX(2)=0.D+0
      RETURN
    2 FX=-X(1)
      RETURN
    3 GF(1)=-1.D+0
      GF(2)=0.D+0
      RETURN
    4 IF (INDEX1(1)) G(1)=-(X(1)-1.D+0)**3-X(2)
      RETURN
    5 IF (INDEX2(1)) GG(1,1)=-3.D+0*(X(1)-1.D+0)**2
      RETURN
      END
C      
      SUBROUTINE TP222(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=1
      NELI=0
      NENL=0
      X(1)=1.3D+0
      X(2)=.2D+0
      DO 6 I=1,2
      LXU(I)=.FALSE.
      LXL(I)=.TRUE.
    6 XL(I)=0.D+0
      GG(1,2)=-1.D+0
      LEX=.TRUE.
      NEX=1
      FEX=-1.5D+0
      XEX(1)=1.5D+0
      XEX(2)=0.D+0
      RETURN
    2 FX=-X(1)
      RETURN
    3 GF(1)=-1.D+0
      GF(2)=0.D+0
      RETURN
    4 IF (INDEX1(1)) G(1)=.125D+0-(X(1)-1.D+0)**3-X(2)
      RETURN
    5 IF (INDEX2(1)) GG(1,1)=-3.D+0*(X(1)-1.D+0)**2
      RETURN
      END
C      
      SUBROUTINE TP223(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=2
      NELI=0
      NENL=0
      X(1)=.1D+0
      X(2)=3.3D+0
      do 6 I=1,2
      LXU(I)=.TRUE.
      LXL(I)=.TRUE.
      XU(I)=10.D+0
    6 XL(I)=0.D+0
      XU(1)=1.0
      GG(1,2)=0.D+0
      GG(2,2)=1.D+0
      LEX=.TRUE.
      NEX=1
      FEX=-DLOG(DLOG(10.D+0))
      XEX(1)=DLOG(DLOG(10.D+0))
      XEX(2)=10.D+0
      RETURN
    2 FX=-X(1)
      RETURN
    3 GF(1)=-1.D+0
      GF(2)=0.D+0
      RETURN
    4 IF (INDEX1(1)) G(1)=EXP(EXP(X(1)))
      IF (INDEX1(2)) G(2)=X(2)-EXP(EXP(X(1)))
      RETURN
    5 IF (INDEX2(1)) GG(1,1)=EXP(X(1))*EXP(EXP(X(1)))
      IF (INDEX2(2)) GG(2,1)=-EXP(X(1))*EXP(EXP(X(1)))
      RETURN
      END
C      
      SUBROUTINE TP224(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=4
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,2
      X(I)=.1D+0
      LXU(I)=.TRUE.
      LXL(I)=.TRUE.
      XU(I)=6.D+0
    6 XL(I)=0.D+0
      GG(1,1)=1.D+0
      GG(1,2)=3.D+0
      GG(2,1)=-1.D+0
      GG(2,2)=-3.D+0
      GG(3,1)=1.D+0
      GG(3,2)=1.D+0
      GG(4,1)=-1.D+0
      GG(4,2)=-1.D+0
      LEX=.TRUE.
      NEX=1
      FEX=-304.D+0
      XEX(1)=4.D+0
      XEX(2)=4.D+0
      RETURN
    2 FX=2.D+0*X(1)**2+X(2)**2-48.D+0*X(1)-40.D+0*X(2)
      RETURN
    3 GF(1)=4.D+0*X(1)-48.D+0
      GF(2)=2.D+0*X(2)-40.D+0
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)+3.D+0*X(2)
      IF (INDEX1(2)) G(2)=18.D+0-X(1)-3.D+0*X(2)
      IF (INDEX1(3)) G(3)=X(1)+X(2)
      IF (INDEX1(4)) G(4)=8.D+0-X(1)-X(2)
    5 RETURN
      END
C      
      SUBROUTINE TP225(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=5
      NELI=0
      NENL=0
      X(1)=3.D+0
      X(2)=1.D+0
      DO 6 I=1,2
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      GG(1,1)=1.D+0
      GG(1,2)=1.d+0
      GG(4,2)=-1.D+0
      GG(5,1)=-1.D+0
      LEX=.TRUE.
      NEX=1
      FEX=2.D+0
      XEX(1)=1.D+0
      XEX(2)=1.D+0
      RETURN
    2 FX=X(1)**2+X(2)**2
      RETURN
    3 GF(1)=2.D+0*X(1)
      GF(2)=2.D+0*X(2)
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)+ X(2)-1.D+0
      IF (INDEX1(2)) G(2)=X(1)**2+ X(2)**2-1.D+0
      IF (INDEX1(3)) G(3)=9.D+0*X(1)**2+ X(2)**2-9.D+0
      IF (INDEX1(4)) G(4)=X(1)**2.D+0- X(2)
      IF (INDEX1(5)) G(5)=X(2)**2.D+0-X(1)
      RETURN
    5 IF (.NOT.INDEX2(2)) GOTO 7
      GG(2,1)=2.D+0*X(1)
      GG(2,2)=2.D+0*X(2)
    7 IF (.NOT.INDEX2(3)) GOTO 8
      GG(3,1)=18.D+0*X(1)
      GG(3,2)=2.D+0*X(2)
    8 IF (INDEX2(4)) GG(4,1)=2.D+0*X(1)
      IF (INDEX2(5)) GG(5,2)=2.D+0*X(2)
      RETURN
      END
C      
      SUBROUTINE TP226(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=2
      NELI=0
      NENL=0
      X(1)=.8D+0
      X(2)=.05D+0
      DO 6 I=1,2
      LXU(I)=.FALSE.
      LXL(I)=.TRUE.
      XL(I)=0.D+0
    6 XEX(I)=1.D+0/DSQRT(2.D+0)
      NEX=1
      LEX=.TRUE.
      FEX=-.5D+0  
      RETURN
    2 FX=-X(1)*X(2)
      RETURN
    3 GF(1)=-X(2)
      GF(2)=-X(1)
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)**2+X(2)**2
      IF (INDEX1(2)) G(2)=1.D+0-X(1)**2-X(2)**2
      RETURN
    5 IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=2.D+0*X(1)
      GG(1,2)=2.D+0*X(2)
    7 IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,1)=-2.D+0*X(1)
      GG(2,2)=-2.D+0*X(2)
    8 RETURN
      END
C      
      SUBROUTINE TP227(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=2
      NELI=0
      NENL=0
      DO 6 I=1,2
      X(I)=.5D+0
      LXU(I)=.FALSE.
      LXL(I)=.FALSE.
    6 XEX(I)=1.D+0
      GG(1,2)=1.D+0
      GG(2,1)=1.d+0
      LEX=.TRUE.
      NEX=1
      FEX=1.D+0
      RETURN
    2 FX=(X(1)-2.D+0)**2+(X(2)-1.D+0)**2
      RETURN
    3 GF(1)=2.D+0*(X(1)-2.D+0)
      GF(2)=2.D+0*(X(2)-1.D+0)
      RETURN
    4 IF (INDEX1(1)) G(1)=-X(1)**2+X(2)
      IF (INDEX1(2)) G(2)=X(1)-X(2)**2
      RETURN
    5 IF (INDEX2(1)) GG(1,1)=-2.D+0*X(1)
      IF (INDEX2(2)) GG(2,2)=-2.D+0*X(2)
      RETURN
      END
C      
      SUBROUTINE TP228(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=1
      NINL=1
      NELI=0
      NENL=0
      DO 6 I=1,2
      X(I)=0.D+0
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      XEX(1)=0.D+0
      XEX(2)=-3.D+0
      GG(1,1)=-1.D+0
      GG(1,2)=-1.d+0
      LEX=.TRUE.
      NEX=1
      FEX=-3.D+0
      RETURN
    2 FX=X(1)**2+X(2)
      RETURN
    3 GF(1)=2.D+0*X(1)
      GF(2)=1.D+0
      RETURN
    4 IF (INDEX1(1)) G(1)=-X(1)-X(2)+1.D+0
      IF (INDEX1(2)) G(2)=-(X(1)**2+X(2)**2)+9.D+0
      RETURN
    5 IF (.NOT.INDEX2(2)) GOTO 7
      GG(2,1)=-2.D+0*X(1)
      GG(2,2)=-2.D+0*X(2)
    7 RETURN
      END
C      
      SUBROUTINE TP229(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      X(1)=-1.2D+0
      X(2)=1.D+0
      DO 6 I=1,2
      LXU(I)=.TRUE.
      LXL(I)=.TRUE.
      XU(I)=2.D+0
      XL(I)=-2.D+0
    6 XEX(I)=1.D+0
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      RETURN
    2 FX=100.D+0*(X(2)-X(1)**2)**2+(1.D+0-X(1))**2
      RETURN
    3 GF(1)=-400.D+0*X(1)*(X(2)-X(1)**2)-2.D+0*(1.D+0-X(1))
      GF(2)=200.D+0*(X(2)-X(1)**2)
    4 RETURN
      END
C      
      SUBROUTINE TP230(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=2
      NELI=0
      NENL=0
      DO 6 I=1,2
      X(I)=0.D+0
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      XEX(1)=.5D+0
      XEX(2)=.375D+0
      GG(1,2)=1.D+0
      GG(2,2)=1.d+0
      LEX=.TRUE.
      NEX=1
      FEX=.375D+0
      RETURN
    2 FX=X(2)
      RETURN
    3 GF(1)=0.D+0
      GF(2)=1.D+0
      RETURN
    4 IF (INDEX1(1)) G(1)=-2.D+0*X(1)**2+X(1)**3+X(2)
      IF (INDEX1(2)) G(2)=-2.D+0*(1.D+0-X(1))**2+(1.D+0-X(1))**3+X(2)
      RETURN
    5 IF (INDEX2(1)) GG(1,1)=-4.D+0*X(1)+3.D+0*X(1)**2
      IF (INDEX2(2)) GG(2,1)=4.D+0*(1.D+0-X(1))-3.D+0*(1.D+0-X(1))**2
      RETURN
      END
C      
      SUBROUTINE TP231(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=2
      NINL=0
      NELI=0
      NENL=0
      X(1)=-1.2D+0
      X(2)=1.D+0
      DO 6 I=1,2
      XEX(I)=1.D+0
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      GG(1,1)=1.D+0/3.D+0
      GG(1,2)=1.d+0
      GG(2,1)=-1.D+0/3.D+0
      GG(2,2)=1.D+0
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      RETURN
    2 FX=100.D+0*(X(2)-X(1)**2)**2+(1.D+0-X(1))**2
      RETURN
    3 GF(1)=-400.D+0*X(1)*(X(2)-x(1)**2)-2.D+0*(1.D+0-x(1))
      GF(2)=200.D+0*(x(2)-x(1)**2)
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)/3.D+0+X(2)+.1D+0
      IF (INDEX1(2)) G(2)=-X(1)/3.D+0+X(2)+.1D+0
    5 RETURN
      END
C      
      SUBROUTINE TP232(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION HV
      HV=DSQRT(3.D+0)
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=3
      NINL=0
      NELI=0
      NENL=0
      X(1)=2.D+0
      X(2)=.5D+0
      DO 6 I=1,2
      LXU(I)=.FALSE.
      LXL(I)=.TRUE.
    6 XL(I)=0.D+0
      GG(1,1)=1/HV
      GG(1,2)=-1.d+0
      GG(2,1)=1.D+0
      GG(2,2)=HV
      GG(3,1)=-1.D+0
      GG(3,2)=-HV
      LEX=.TRUE.
      NEX=1
      FEX=-1.
      XEX(1)=3.D+0
      XEX(2)=HV
    2 FX=-(9.D+0-(X(1)-3.D+0)**2)*X(2)**3/(27.D+0*HV)
      RETURN
    3 GF(1)=2.D+0*(X(1)-3.D+0)*X(2)**3.D+0/(27.D+0*HV)
      GF(2)=-(9.D+0-(X(1)-3.D+0)**2)*X(2)**2/(9.D+0*HV)
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)/HV-X(2)
      IF (INDEX1(2)) G(2)=X(1)+HV*X(2)
      IF (INDEX1(3)) G(3)=6.D+0-X(1)-HV*X(2)
    5 RETURN
      END
C      
      SUBROUTINE TP233(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=1
      NELI=0
      NENL=0
      X(1)=1.2D+0
      X(2)=1.D+0
      DO 6 I=1,2
      LXU(I)=.FALSE.
      LXL(I)=.FALSE.
    6 XEX(I)=1.D+0
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      RETURN
    2 FX=100.D+0*(X(2)-X(1)**2)**2+(1.D+0-X(1))**2
      RETURN
    3 GF(1)=-400.d+0*X(1)*(X(2)-X(1)**2)-2.D+0*(1.D+0-X(1))
      GF(2)=200.D+0*(X(2)-X(1)**2)
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)**2+X(2)**2-.25D+0
      RETURN
    5 IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=2.D+0*X(1)
      GG(1,2)=2.D+0*X(2)
    7 RETURN
      END
C      
      SUBROUTINE TP234(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=1
      NELI=0
      NENL=0
      DO 6 I=1,2
      X(I)=1.D+0
      LXU(I)=.TRUE.
      LXL(I)=.TRUE.
      XL(I)=0.2D+0
      XU(I)=2.D+0
    6 XEX(I)=.2D+0
      LEX=.TRUE.
      NEX=1
      FEX=-0.8D+0
      RETURN
    2 FX=(X(2)-X(1))**4-(1.D+0-X(1))
      RETURN
    3 GF(1)=-4.D+0*(X(2)-X(1))**3+1.D+0
      GF(2)=4.D+0*(X(2)-X(1))**3
      RETURN
    4 IF (INDEX1(1)) G(1)=-X(1)**2-X(2)**2+1.D+0
      RETURN
    5 IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=-2.D+0*X(1)
      GG(1,2)=-2.D+0*X(2)
    7 RETURN
      END
C      
      SUBROUTINE TP235(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=3
      NILI=0
      NINL=0
      NELI=0
      NENL=1
      X(1)=-2.D+0
      X(2)=3.D+0
      X(3)=1.D+0
      DO 6 I=1,2
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      GG(1,1)=1.D+0
      GG(1,2)=0.D+0
      XEX(1)=-1.D+0
      XEX(2)=1.D+0
      XEX(3)=0.D+0
      LEX=.TRUE.
      NEX=1
      FEX=0.04D0
      RETURN
    2 FX=(X(2)-X(1)**2)**2+0.01D0*(X(1)-1.0D0)**2
      RETURN
    3 GF(1)=-4.0D0*X(1)*(X(2)-X(1)**2)+0.02D0*(X(1)-1.0D0)
      GF(2)=2.0D0*(X(2)-X(1)**2)
      GF(3)=0.0D0
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)+X(3)**2+1.0D0
      RETURN
    5 IF (INDEX2(1)) GG(1,3)=2.0D0*X(3)
      RETURN
      END
C            
      SUBROUTINE TP236239(IMODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    IMODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION B(20)
      DATA B /75.1963666677D+0,-3.8112755343D+0,0.1269366345D+0,
     1        -2.0567665D-3,1.0345D-5,-6.8306567613D+0,3.02344793D-2,
     2        -1.2813448D-3,3.52559D-5,-2.266D-7,0.2564581253D+0,
     3        -3.460403D-3,1.35139D-5,-28.1064434908D+0,-5.2375D-6,
     4        -6.3D-9,7.D-10,3.405462D-4,-1.6638D-6,-2.8673112392D+0/
      GOTO (2,3),IMODE
    2 FX=B(1)+B(2)*X(1)+B(3)*X(1)**2+B(4)*X(1)**3 +
     1   B(5)*X(1)**4+B(6)*X(2)+B(7)*X(1)*X(2)+B(8)*X(1)**2*X(2) +
     2   B(9)*X(1)**3*X(2)+B(10)*X(1)**4*X(2)+B(11)*X(2)**2 +
     3   B(12)*X(2)**3+B(13)*X(2)**4+B(14)*(1.D+0/(X(2)+1.D+0)) +
     4   B(15)*X(1)**2*X(2)**2+B(16)*X(1)**3*X(2)**2 +
     5   B(17)*X(1)**3*X(2)**3+B(18)*X(1)*X(2)**2 +
     6   B(19)*X(1)*X(2)**3+B(20)*(DEXP(5.D-4*X(1)*X(2)))
      FX=-FX
      RETURN
    3 GF(1)=B(2)+B(3)*2.D+0*X(1)+B(4)*3.D+0*X(1)**2 +
     1      B(5)*4.D+0*X(1)**3+B(7)*X(2)+B(8)*2.D+0*X(1)*X(2) +
     2      B(9)*3.D+0*X(1)**2*X(2)+B(10)*4.D+0*X(1)**3*X(2) +
     3      B(15)*2.D+0*X(1)*X(2)**2+B(16)*3.D+0*X(1)**2*X(2)**2 +
     4      B(17)*3.D+0*X(1)**2*X(2)**3+B(18)*X(2)**2+B(19)*X(2)**3 +
     5      B(20)*(DEXP(5.D-4*X(1)*X(2)))*(5.D-4*X(2))
      GF(1)=-GF(1)
      GF(2)=B(6)+B(7)*X(1)+B(8)*X(1)**2+B(9)*X(1)**3 +
     1      B(10)*X(1)**4+B(11)*2.D+0*X(2)+B(12)*3.D+0*X(2)**2 +
     2      B(13)*4.D+0*X(2)**3+B(14)*(-1.D+0/(X(2)+1.D+0)**2) +
     3      B(15)*X(1)**2*2*X(2)+B(16)*X(1)**3*2*X(2) +
     4      B(17)*X(1)**3*3*X(2)**2+B(18)*X(1)*2.D+0*X(2)+
     5      B(19)*X(1)*3.D+0*X(2)**2 +
     6      B(20)*(DEXP(5.D-4*X(1)*X(2)))*(5.D-4*X(1))
      GF(2)=-GF(2)
      RETURN
      END
C
      SUBROUTINE TP236(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=2
      NELI=0
      NENL=0
      X(1)=1.D+1
      X(2)=1.D+1
      DO 6 I=1,2
      LXU(I)=.TRUE.
      LXL(I)=.TRUE.
    6 XL(I)=0.D+0
      XU(1)=75.D+0
      XU(2)=65.D+0
      GG(2,2)=1.D+0
      LEX=.TRUE.
      NEX=1
      FEX=-58.9034360D+0
      XEX(1)=75.D+0
      XEX(2)=65.D+0
      RETURN
    2 CALL TP236239(1)
      RETURN
    3 CALL TP236239(2)
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)*X(2)-7.D+2
      IF (INDEX1(2)) G(2)=X(2)-5.D+0*((X(1)/25.D+0)**2)
      RETURN
    5 IF (.NOT. INDEX2(1)) GOTO 7
      GG(1,1)=X(2)
      GG(1,2)=X(1)
    7 IF (INDEX2(2)) GG(2,1)=-10.D+0/625.D+0*X(1)
      RETURN
      END
C
      SUBROUTINE TP237(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      GOTO (1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=3
      NELI=0
      NENL=0
      X(1)=65.0D+0
      X(2)=10.0D+0
      LXU(1)=.TRUE.
      LXU(2)=.TRUE.
      LXL(1)=.TRUE.
      LXL(2)=.FALSE.
      XU(1)=75.D+0
      XU(2)=65.D+0
      XL(1)=54.D+0
      GG(2,2)=1.D+0
      GG(3,1)=-5.D+0
      LEX=.TRUE.
      NEX=1
      FEX=-58.9034360D+0
      XEX(1)=75.0D+0
      XEX(2)=65.D+0
      RETURN
    2 CALL TP236239(1)
      RETURN
    3 CALL TP236239(2)
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)*X(2)-7.D+2
      IF (INDEX1(2)) G(2)=X(2)-5.D+0*((X(1)/25.D+0)**2)
      IF (INDEX1(3)) G(3)=(X(2)-5.D+1)**2-5.D+0*(X(1)-55.D+0)
      RETURN
    5 IF (.NOT. INDEX2(1)) GOTO 6
      GG(1,1)=X(2)
      GG(1,2)=X(1)
    6 IF (INDEX2(2)) GG(2,1)=-10.D+0/625.D+0*X(1)
      IF (INDEX2(3)) GG(3,2)=2.D+0*X(2)-1.D+2
      RETURN
      END
C
      SUBROUTINE TP238(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=3
      NELI=0
      NENL=0
      X(1)=1.D+1
      X(2)=1.D+1
      DO 6 I=1,2
      LXU(I)=.TRUE.
    6 LXL(I)=.FALSE.
      XU(1)=75.D+0
      XU(2)=65.D+0
      GG(2,2)=1.D+0
      GG(3,1)=-5.D+0
      LEX=.TRUE.
      NEX=1
      FEX=-58.9034360D+0
      XEX(1)=75.D+0
      XEX(2)=65.D+0
      RETURN
    2 CALL TP236239(1)
      RETURN
    3 CALL TP236239(2)
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)*X(2)-7.D+2
      IF (INDEX1(2)) G(2)=X(2)-5.D+0*((X(1)/25.D+0)**2)
      IF (INDEX1(3)) G(3)=(X(2)-5.D+1)**2-5.D+0*(X(1)-55.D+0)
      RETURN
    5 IF (.NOT. INDEX2(1)) GOTO 7
      GG(1,1)=X(2)
      GG(1,2)=X(1)
    7 IF (INDEX2(2)) GG(2,1)=-10.D+0/625.D+0*X(1)
      IF (INDEX2(3)) GG(3,2)=2.D+0*X(2)-1.D+2
      RETURN
      END
C
      SUBROUTINE TP239(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=1
      NELI=0
      NENL=0
      X(1)=1.D+1
      X(2)=1.D+1
      DO 6 I=1,2
      LXU(I)=.TRUE.
      LXL(I)=.TRUE.
    6 XL(I)=0.D+0
      XU(1)=75.D+0
      XU(2)=65.D+0
      LEX=.TRUE.
      NEX=1
      FEX=-58.903436D0
      XEX(1)=75.0D0
      XEX(2)=65.0D0
      RETURN
    2 CALL TP236239(1)
      RETURN
    3 CALL TP236239(2)
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)*X(2)-7.D+2
      RETURN
    5 IF (.NOT. INDEX2(1)) GOTO 7
      GG(1,1)=X(2)
      GG(1,2)=X(1)
    7 RETURN
      END
C
      SUBROUTINE TP240(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,4),MODE
    1 N=3
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      X(1)=1.D+2
      X(2)=-1.D+0
      X(3)=2.5D+0
      DO 6 I=1,3
      LXL(I)=.FALSE.
      LXU(I)=.FALSE.
    6 XEX(I)=0.D+0
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      RETURN
    2 FX=(X(1)-X(2)+X(3))**2+(-X(1)+X(2)+X(3))**2+(X(1)+X(2)
     /       -X(3))**2
    3 GF(1)=2.D+0*(X(1)-X(2)+X(3))-2.D+0*(X(2)-X(1)+X(3))
     /      +2.D+0*(X(1)+X(2)-X(3))
      GF(2)=2.D+0*(X(2)-X(1)+X(3))-2.D+0*(X(1)-X(2)+X(3))
     /      +2.D+0*(X(1)+X(2)-X(3))
      GF(3)=2.D+0*(X(1)-X(2)+X(3))+2.D+0*(X(2)-X(1)+X(3))
     /      -2.D+0*(X(1)+X(2)-X(3))
    4 RETURN
      END
C      
      SUBROUTINE TP241(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF
      INTEGER I
      GOTO (1,2,2,4,4),MODE
    1 N=3
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      X(1)=1.D+0
      X(2)=2.D+0
      X(3)=0.D+0
      DO 6 I=1,3
      LXL(I)=.FALSE.
    6 LXU(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      XEX(1)=0.D+0
      XEX(2)=0.D+0
      XEX(3)=1.D+0
      RETURN
    2 F(1)=X(1)**2+X(2)**2+X(3)**2-1.D+0
      F(2)=X(1)**2+X(2)**2+(X(3)-2.D+0)**2-1.D+0
      F(3)=X(1)+X(2)+X(3)-1.D+0
      F(4)=X(1)+X(2)-X(3)+1.D+0
      F(5)=X(1)**3+3.D+0*X(2)**2+(5.D+0*X(3)-X(1)+1.D+0)**2-3.6D+1
      IF (MODE.EQ.3) GOTO 3
      FX=F(1)**2+F(2)**2+F(3)**2+F(4)**2+F(5)**2
      RETURN
    3 DO 7 I=1,3
      DF(1,I)=2.D+0*X(I)
    7 DF(3,I)=1.D+0
      DF(2,1)=DF(1,1)
      DF(2,2)=DF(1,2)
      DF(2,3)=2.D+0*(X(3)-2.D+0)
      DF(4,1)=1.D+0
      DF(4,2)=1.D+0
      DF(4,3)=-1.D+0
      DF(5,1)=3.D+0*X(1)**2-2.D+0*(5.D+0*X(3)-X(1)+1.D+0)
      DF(5,2)=6.D+0*X(2)
      DF(5,3)=1.D+1*(5.D+0*X(3)-X(1)+1.D+0)
      GF(1)=0.D+0
      GF(2)=0.D+0
      GF(3)=0.D+0
      DO 8 I=1,5
      GF(1)=GF(1)+F(I)*DF(I,1)*2.D+0
      GF(2)=GF(2)+F(I)*DF(I,2)*2.D+0
    8 GF(3)=GF(3)+F(I)*DF(I,3)*2.D+0
    4 RETURN
      END
C      
      SUBROUTINE TP242(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION TI,F,DF
      INTEGER I
      GOTO (1,2,2,4,4),MODE
    1 N=3
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      X(1)=2.5D+0
      X(2)=1.D+1
      X(3)=X(2)
      DO 6 I=1,3
      LXL(I)=.TRUE.
      LXU(I)=.TRUE.
      XL(I)=0.D+0
    6 XU(I)=1.D+1
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      XEX(1)=1.D+0
      XEX(2)=1.D+1
      XEX(3)=1.D+0
      RETURN
    2 FX=0.D+0
      DO 7 I=1,10
      TI=0.1D+0*DBLE(I)
    7 F(I)=DEXP(-X(1)*TI)-DEXP(-X(2)*TI)-X(3)*
     F      (DEXP(-TI)-DEXP(-1.D+1*TI))
      IF (MODE.EQ.3) GOTO 3
      DO 8 I=1,10
    8 FX=FX+F(I)**2
      RETURN
    3 GF(1)=0.D+0
      GF(2)=0.D+0
      GF(3)=0.D+0
      DO 9 I=1,10
      TI=0.1D+0*DBLE(I)
      F(I)=DEXP(-X(1)*TI)-DEXP(-X(2)*TI)-X(3)*
     F     (DEXP(-TI)-DEXP(-1.D+1*TI))
      DF(I,1)=-TI*DEXP(-X(1)*TI)
      DF(I,2)=TI*DEXP(-X(2)*TI)
      DF(I,3)=DEXP(-1.D+1*TI)-DEXP(-TI)
      GF(1)=GF(1)+F(I)*DF(I,1)*2.D+0
      GF(2)=GF(2)+F(I)*DF(I,2)*2.D+0
    9 GF(3)=GF(3)+F(I)*DF(I,3)*2.D+0
    4 RETURN
      END
C      
      SUBROUTINE TP243(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION XBX,DXBX(3),A(4),B(3,3),D(4),E(4,3),F,DF
      DATA A,D,B /0.14272D+0,-0.184981D+0,-0.521869D+0,-0.685306D+0,
     F             1.75168D+0,-1.35195D+0,-0.479048D+0,-0.3648D+0,
     F             2.95137D+0,4.87407D+0,-2.0506D+0,
     F             4.87407D+0,9.39321D+0,-3.93185D+0,
     F             -2.0506D+0,-3.93189D+0,2.64745D+0 /
      INTEGER I
      E(1,1)=-0.564255D+0
      E(1,2)=0.392417D+0
      E(1,3)=-0.404979D+0
      E(2,1)=0.927589D+0
      E(2,2)=-0.0735083D+0
      E(2,3)=0.535493D+0
      E(3,1)=0.658799D+0
      E(3,2)=-0.636666D+0
      E(3,3)=-0.681091D+0
      E(4,1)=-0.869487D+0
      E(4,2)=0.586387D+0
      E(4,3)=0.289826D+0
      GOTO (1,2,2,4,4),MODE
    1 N=3
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,3
      LXL(I)=.FALSE.
      LXU(I)=.FALSE.
      X(I)=0.1D+0
    6 XEX(I)=0.D+0
      LEX=.TRUE.
      NEX=1
      FEX=0.79657853D+00
      RETURN
    2 FX=0.D+0
      XBX=(X(1)*B(1,1)+X(2)*B(2,1)+X(3)*B(3,1))*X(1)+
     /    (X(1)*B(1,2)+X(2)*B(2,2)+X(3)*B(3,2))*X(2)+
     /    (X(1)*B(1,3)+X(2)*B(2,3)+X(3)*B(3,3))*X(3)
      DO 7 I=1,4
    7 F(I)=A(I)+E(I,1)*X(1)+E(I,2)*X(2)+E(I,3)*X(3)+
     /     0.5D+0*XBX*D(I)
      IF (MODE.EQ.3) GOTO 3
      DO 10 I=1,4
   10 FX=FX+F(I)**2
      RETURN
    3 DXBX(1)=(X(1)*B(1,1)+X(2)*B(2,1)+X(3)*B(3,1))*2.D+0
      DXBX(2)=(X(1)*B(1,2)+X(2)*B(2,2)+X(3)*B(3,2))*2.D+0
      DXBX(3)=(X(1)*B(1,3)+X(2)*B(2,3)+X(3)*B(3,3))*2.D+0
      DO 9 I=1,3
    9 GF(I)=0.D+0
      DO 8 I=1,4
      DF(I,1)=E(I,1)+DXBX(1)*D(I)*0.5D+0
      DF(I,2)=E(I,2)+DXBX(2)*D(I)*0.5D+0
      DF(I,3)=E(I,3)+DXBX(3)*D(I)*0.5D+0
      GF(1)=GF(1)+F(I)*DF(I,1)*2.D+0
      GF(2)=GF(2)+F(I)*DF(I,2)*2.D+0
    8 GF(3)=GF(3)+F(I)*DF(I,3)*2.D+0
    4 RETURN
      END
C
      SUBROUTINE TP244(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION YI,ZI,F,DF
      INTEGER I
      GOTO (1,2,2,4,4),MODE
    1 N=3
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      NEX=1
      X(1)=1.D+0
      X(2)=2.D+0
      X(3)=1.D+0
      DO 6 I=1,3
      XL(I)=0.0
      XU(I)=1.0D+10
      LXL(I)=.TRUE.
    6 LXU(I)=.TRUE.
      LEX=.TRUE.
      FEX=0.D+0
      XEX(1)=1.D+0
      XEX(2)=1.D+1
      XEX(3)=5.D+0
      RETURN
    2 FX=0.D+0
      DO 7 I=1,10
      ZI=0.1D+0*DBLE(I)
      YI=DEXP(-ZI)-5.D+0*DEXP(-1.D+1*ZI)
    7 F(I)=DEXP(-X(1)*ZI)-X(3)*DEXP(-X(2)*ZI)-YI
      IF (MODE.EQ.3) GOTO 3
      DO 10 I=1,8
   10 FX=FX+F(I)**2
      RETURN
    3 DO 8 I=1,3
    8 GF(I)=0.D+0
      DO 9 I=1,10
      ZI=0.1D+0*DBLE(I)
      YI=DEXP(-ZI)-5.D+0*DEXP(-1.D+1*ZI)
      DF(I,1)=-ZI*DEXP(-X(1)*ZI)
      DF(I,2)=ZI*X(3)*DEXP(-X(2)*ZI)
      DF(I,3)=-DEXP(-X(2)*ZI)
      GF(1)=GF(1)+F(I)*DF(I,1)*2.D+0
      GF(2)=GF(2)+F(I)*DF(I,2)*2.D+0
    9 GF(3)=GF(3)+F(I)*DF(I,3)*2.D+0
    4 RETURN
      END
C
      SUBROUTINE TP245(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF,DI
      INTEGER I
      GOTO (1,2,2,4,4),MODE
    1 N=3
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      NEX=1
      X(1)=0.D+0
      X(2)=1.D+1
      X(3)=2.D+1
      DO 6 I=1,3
      XL(I)=0.0
      LXL(I)=.TRUE.
      XU(I)=20.0
    6 LXU(I)=.TRUE.
      XU(1)=12.0D0
      XU(2)=12.0D0
      LEX=.TRUE.
      FEX=0.D+0
      XEX(1)=1.D+0
      XEX(2)=1.D+1
      XEX(3)=1.D+0
      RETURN
    2 FX=0.D+0
      DO 7 I=1,10
      DI=DBLE(I)
    7 F(I)=DEXP(-DI*X(1)/1.D+1)-DEXP(-DI*X(2)/1.D+1)-
     /     X(3)*(DEXP(-DI/1.D+1)-DEXP(-DI))
      IF (MODE.EQ.3) GOTO 3
      DO 10 I=1,10
   10 FX=FX+F(I)**2
      RETURN
    3 DO 8 I=1,3
    8 GF(I)=0.D+0
      DO 9 I=1,10
      DI=DBLE(I)
      DF(I,1)=-DI/1.D+1*DEXP(-DI*X(1)/1.D+1)
      DF(I,2)=DI/1.D+1*DEXP(-DI*X(2)/1.D+1)
      DF(I,3)=DEXP(-DI)-DEXP(-DI/1.D+1)
      GF(1)=GF(1)+F(I)*DF(I,1)*2.D+0
      GF(2)=GF(2)+F(I)*DF(I,2)*2.D+0
    9 GF(3)=GF(3)+F(I)*DF(I,3)*2.D+0
    4 RETURN
      END
C
      SUBROUTINE TP246(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF
      INTEGER I,J
      GOTO (1,2,2,4,4),MODE
    1 N=3
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      X(1)=-1.2D+0
      X(2)=2.D+0
      X(3)=0.D+0
      DO 6 I=1,3
      DO 7 J=1,3
    7 DF(I,J)=0.D+0
      LXL(I)=.FALSE.
      LXU(I)=.FALSE.
    6 XEX(I)=1.D+0
      DF(1,3)=10.D+0
      DF(2,1)=-1.D+0
      DF(3,2)=-1.D+0
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      RETURN
    2 F(1)=10.D+0*(X(3)-((X(1)+X(2))/2.D+0)**2)
      F(2)=1.D+0-X(1)
      F(3)=1.D+0-X(2)
      IF (MODE.EQ.3) GOTO 3
      FX=F(1)**2+F(2)**2+F(3)**2
      RETURN
    3 DF(1,1)=-10.D+0*(X(1)+X(2))
      DF(1,2)=-10.D+0*(X(1)+X(2))
      DO 8 I=1,3
      GF(I)=0.D+0
      DO 8 J=1,3
    8 GF(I)=GF(I)+2.D+0*F(J)*DF(J,I) 
    4 RETURN
      END
C
      SUBROUTINE TP247(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION THETA,DTHETA(3),XPI
      INTEGER I
      GOTO (1,2,3,4,4),MODE
    1 N=3
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      NEX=1
      DO 6 I=1,2
      LXL(I)=.FALSE.
    6 LXU(I)=.FALSE.
      LXL(3)=.TRUE.
      LXU(3)=.TRUE.
      LXL(1)=.TRUE.
      XL(1)=0.1D+0
      XL(3)=-2.5D+0
      XU(3)=7.5D+0
      LEX=.TRUE.
      X(1)=0.1
      X(2)=0.D+0
      X(3)=0.D+0
      FEX=0.D+0
      XEX(1)=1.D+0
      XEX(2)=0.D+0
      XEX(3)=0.D+0
      RETURN
    2 XPI=DASIN(1.D+0)*2.D+0
      THETA=1.D+0/(2.D+0*XPI)*DATAN(X(2)/X(1))
      IF (X(1).LT.0.D+0) THETA=THETA+0.5D+0
      FX=1.D+2*((X(3)-1.D+1*THETA)**2+
     /   (DSQRT(X(1)**2+X(2)**2)-1.D+0)**2)+X(3)**2
      RETURN
    3 XPI=DASIN(1.D+0)*2.D+0
      THETA=1.D+0/(2.D+0*XPI)*DATAN(X(2)/X(1))
      DTHETA(1)=-X(2)/((1.D+0+(X(2)/X(1))**2)*X(1)**2)
      DTHETA(2)=1.D+0/((1.D+0+(X(2)/X(1))**2)*X(1))
      DTHETA(3)=0.D+0
      IF (X(1).LT.0.D+0) THETA=THETA+0.5D+0
      GF(1)=1.D+2*(2.D+1*(X(3)-1.D+1*THETA)*DTHETA(1)+
     /      2.D+0*(DSQRT(X(1)**2+X(2)**2)-1.D+0)/(DSQRT(X(1)**2
     /      +X(2)**2))*X(1))
      GF(2)=1.D+2*(2.D+1*(X(3)-1.D+1*THETA)*DTHETA(2)+
     /      2.D+0*(DSQRT(X(1)**2+X(2)**2)-1.D+0)/(DSQRT(X(1)**2
     /      +X(2)**2))*X(2))
      GF(3)=1.D+2*(2.D+0*(X(3)-1.D+1*THETA))+2.D+0*X(3)
    4 RETURN
      END
C      
      SUBROUTINE TP248(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=3
      NILI=1
      NINL=0
      NELI=0
      NENL=1
      X(1)=-0.1D+0
      X(2)=-1.D+0
      X(3)=0.1D+0
      DO 6 I=1,3
      LXL(I)=.FALSE.
    6 LXU(I)=.FALSE.
      GG(1,1)=1.D+0
      GG(1,2)=-2.D+0
      GG(1,3)=0.D+0
      LEX=.TRUE.
      NEX=1
      FEX=-0.8D+0
      XEX(1)=0.6D+0
      XEX(2)=0.8D+0
      XEX(3)=0.D+0
      GF(1)=0.D+0
      GF(2)=-1.D+0
      GF(3)=0.D+0
      RETURN
    2 FX=-X(2)
    3 RETURN
    4 IF (INDEX1(1)) G(1)=1.D+0-2.D+0*X(2)+X(1)
      IF (INDEX1(2)) G(2)=X(1)**2+X(2)**2+X(3)**2-1.D+0
      RETURN
    5 IF (.NOT.INDEX2(2)) GOTO 8
      DO 9 I=1,3
    9 GG(2,I)=2.D+0*X(I)
    8 RETURN
      END
C      
      SUBROUTINE TP249(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=3
      NILI=0
      NINL=1
      NELI=0
      NENL=0
      DO 6 I=1,3
      X(I)=1.D+0
      LXL(I)=.FALSE.
    6 LXU(I)=.FALSE.
      LXL(1)=.TRUE.
      GG(1,3)=0.D+0
      XL(1)=1.D+0
      LEX=.TRUE.
      NEX=1
      FEX=1.D+0
      XEX(1)=1.D+0
      XEX(2)=0.D+0
      XEX(3)=0.D+0
      RETURN
    2 FX=X(1)**2+X(2)**2+X(3)**2
      RETURN
    3 DO 7 I=1,3
    7 GF(I)=2.D+0*X(I)
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)**2+X(2)**2-1.D+0
      RETURN
    5 IF (.NOT.INDEX2(1)) GOTO 8
      GG(1,1)=2.D+0*X(1)
      GG(1,2)=2.D+0*X(2)
    8 RETURN
      END
C      
      SUBROUTINE TP250(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=3
      NILI=2
      NINL=0
      NELI=0
      NENL=0
      NEX=1
      DO 6 I=1,3
      LXL(I)=.TRUE.
      LXU(I)=.TRUE.
      XL(I)=0.D+0
    6 X(I)=1.D+1
      XU(1)=2.D+1
      XU(2)=1.1D+1
      XU(3)=4.2D+1
      LEX=.TRUE.
      FEX=-3.3D+3
      XEX(1)=2.0D+1
      XEX(2)=1.1D+1
      XEX(3)=1.5D+1
      GG(1,1)=1.D+0
      GG(1,2)=2.D+0
      GG(1,3)=2.D+0
      GG(2,1)=-1.D+0
      GG(2,2)=-2.D+0
      GG(2,3)=-2.D+0
      RETURN
    2 FX=-X(1)*X(2)*X(3)
      RETURN
    3 GF(1)=-X(2)*X(3)
      GF(2)=-X(1)*X(3)
      GF(3)=-X(1)*X(2)
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)+2.D+0*X(2)+2.D+0*X(3)
      IF (INDEX1(2)) G(2)=7.2D+1-X(1)-2.D+0*X(2)-2.D+0*X(3)
    5 RETURN
      END
C      
      SUBROUTINE TP251(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=3
      NILI=1
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,3
      X(I)=1.D+1
      LXL(I)=.TRUE.
      LXU(I)=.TRUE.
      XL(I)=0.D+0
    6 XU(I)=4.2D+1
      GG(1,1)=-1.D+0
      GG(1,2)=-2.D+0
      GG(1,3)=-2.D+0
      LEX=.TRUE.
      NEX=1
      FEX=-3.456D+3
      XEX(1)=2.4D+1
      XEX(2)=1.2D+1
      XEX(3)=1.2D+1
      RETURN
    2 FX=-X(1)*X(2)*X(3)
      RETURN
    3 GF(1)=-X(2)*X(3)
      GF(2)=-X(1)*X(3)
      GF(3)=-X(1)*X(2)
      RETURN
    4 IF (INDEX1(1)) G(1)=7.2D+1-X(1)-2.D+0*X(2)-2.D+0*X(3)
    5 RETURN
      END
C      
      SUBROUTINE TP252(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX)
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=3
      NILI=0
      NINL=0
      NELI=0
      NENL=1
      X(1)=-1.D+0
      X(2)=2.D+0
      X(3)=2.D+0
      DO 6 I=1,3
      LXL(I)=.FALSE.
    6 LXU(I)=.FALSE.
      LXU(1)=.TRUE.
      XU(1)=-1.D+0
      GG(1,1)=1.D+0
      GG(1,2)=0.D+0
      LEX=.TRUE.
      NEX=1
      FEX=0.4D-1
      XEX(1)=-1.D+0
      XEX(2)=1.D+0
      XEX(3)=0.D+0
      GF(3)=0.D+0
      RETURN
    2 FX=0.1D-1*(X(1)-1.D+0)**2+(X(2)-X(1)**2)**2
      RETURN
    3 GF(1)=0.2D-1*(X(1)-1.D+0)-4.D+0*(X(2)-X(1)**2)*X(1)
      GF(2)=2.D+0*(X(2)-X(1)**2)
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)+X(3)**2+1.D+0
      RETURN
    5 IF (INDEX2(1)) GG(1,3)=2.D+0*X(3)
      RETURN
      END
C      
      SUBROUTINE TP253(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION A(3,8)
      DATA ((A(I,J),I=1,3),J=1,8)
     1     /3*0.D+0,10.D+0,2*0.D+0,2*10.D+0,2*0.D+0,10.D+0,3*0.D+0,
     2      2*10.D+0,0.D+0,4*10.D+0,0.D+0,2*10.D+0/
      INTEGER I,J
      GOTO (1,2,3,4,5),MODE
    1 N=3
      NILI=1
      NINL=0
      NELI=0
      NENL=0
      X(1)=0.D+0
      X(2)=2.D+0
      X(3)=0.D+0
      DO 6 I=1,3
      LXU(I)=.FALSE.
      LXL(I)=.TRUE.
    6 XL(I)=0.D+0
      GG(1,1)=-3.D+0
      GG(1,2)=0.D+0
      GG(1,3)=-3.D+0
      LEX=.TRUE.
      NEX=1
      FEX=0.69282032D+02
      DO 7 I=1,3
    7 XEX(I)=5.0D0
      RETURN
    2 FX=0.D+0
      DO 8 J=1,8
    8 FX=FX+DSQRT((A(1,J)-X(1))**2+(A(2,J)-X(2))**2 +
     /      (A(3,J)-X(3))**2)
      RETURN
    3 DO 9 I=1,3
      GF(I)=0.D+0
      DO 9 J=1,8
    9 GF(I)=GF(I)+(X(I)-A(I,J))/DSQRT((A(1,J)-X(1))**2 +
     /      (A(2,J)-X(2))**2+(A(3,J)-X(3))**2)
      RETURN
    4 IF (INDEX1(1)) G(1)=3.D+1-3.D+0*X(1)-3.D+0*X(3)
    5 RETURN
      END
C
      SUBROUTINE TP254(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I     
      GOTO (1,2,3,4,5),MODE
    1 N=3
      NILI=0
      NINL=0
      NELI=0
      NENL=2
      X(1)=1.D+0
      X(2)=1.D+0
      X(3)=1.D+0
      DO 6 I=1,2
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LXU(3)=.FALSE.
      LXL(3)=.TRUE.
      XL(3)=1.D+0
      GG(1,1)=0.D+0
      GG(2,2)=0.D+0
      GG(2,3)=1.D+0
      LEX=.TRUE.
      NEX=1
      FEX=-DSQRT(3.D+0)
      XEX(1)=0.D+0
      XEX(2)=DSQRT(3.D+0)
      XEX(3)=1.D+0
      GF(1)=0.D+0
      GF(2)=-1.D+0
      RETURN
    2 FX=DLOG10(X(3))-X(2)
      RETURN
    3 GF(3)=1.D+0/(X(3)*DLOG(1.D+1))
      RETURN
    4 IF (INDEX1(1)) G(1)=X(2)**2+X(3)**2-4.D+0
      IF (INDEX1(2)) G(2)=X(3)-1.D+0-X(1)**2
      RETURN
    5 IF (.NOT. INDEX2(1)) GOTO 7
      GG(1,2)=2.D+0*X(2)
      GG(1,3)=2.D+0*X(3)
    7 IF (INDEX2(2)) GG(2,1)=-2.D+0*X(1)
      RETURN
      END
C
      SUBROUTINE TP255(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,4),MODE
    1 N=4
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      X(1)=-3.D+0
      X(2)=1.D+0
      X(3)=-3.D+0
      X(4)=1.D+0
      DO 6 I=1,4
      XU(I)=10.0
      LXU(I)=.TRUE.
      XL(I)=-10.0
    6 LXL(I)=.TRUE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      DO 7 I=1,4
    7 XEX(I)=1.D+0
      RETURN
    2 FX=100.0*(X(2)-X(1)**2)+(1.0-X(1))**2 +
     1   90.0*(X(4)-X(3)**2)+(1.0-X(3))**2 +
     2   10.1*((X(2)-1.0)**2+(X(4)-1.0)**2) +
     3   19.8*(X(2)-1.0)*(X(4)-1.0)  
      FX=0.5*FX**2
      RETURN
    3 CONTINUE
      FX=100.0*(X(2)-X(1)**2)+(1.0-X(1))**2 +
     1   90.0*(X(4)-X(3)**2)+(1.0-X(3))**2 +
     2   10.1*((X(2)-1.0)**2+(X(4)-1.0)**2) +
     3   19.8*(X(2)-1.0)*(X(4)-1.0)  
      GF(1)=FX*(-198.D+0*X(1)-2.D+0)
      GF(2)=FX*(20.2D+0*X(2)+19.8D+0*X(4)+6.D+1)
      GF(3)=FX*(-178.D+0*X(3)-2.D+0)
      GF(4)=FX*(19.8D+0*X(2)+20.2D+0*X(4)+5.D+1)
    4 RETURN
      END
C
      SUBROUTINE TP256(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,4),MODE
    1 N=4
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      X(1)=3.D+0
      X(2)=-1.D+0
      X(3)=0.D+0
      X(4)=1.D+0
      DO 6 I=1,4
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      DO 7 I=1,4
    7 XEX(I)=0.D+0
      RETURN
    2 FX=((X(1)+1.D+1*X(2))**2+5.D+0*(X(3)-X(4))**2 +
     /   (X(2)-2.D+0*X(3))**4+1.D+1*(X(1)-X(4))**4)
      RETURN
    3 GF(1)=(2.D+0*(X(1)+1.D+1*X(2))+4.D+1*(X(1)-X(4))**3)
      GF(2)=(2.D+1*(X(1)+1.D+1*X(2))+4.D+0*(X(2)-2.D+0*X(3))**3)
      GF(3)=(1.D+1*(X(3)-X(4))-8.D+0*(X(2)-2.D+0*X(3))**3)
      GF(4)=(-1.D+1*(X(3)-X(4))-4.D+1*(X(1)-X(4))**3)
    4 RETURN
      END
C
      SUBROUTINE TP257(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,4),MODE
    1 N=4
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,4
      LXL(I)=.TRUE.
      XL(I)=0.0D0
      X(I)=0.0D0
    6 LXU(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      DO 7 I=1,4
    7 XEX(I)=1.D+0
      RETURN
    2 FX=1.D+2*(X(1)**2-X(2))**2+(X(1)-1.D+0)**2 +
     1   9.D+1*(X(3)**2-X(4))**2+(X(3)-1.D+0)**2 +
     2   10.1D+0*((X(2)-1.D+0)**2+(X(4)-1.D+0)**2) +
     3   19.8D+0*(X(1)-1.D+0)*(X(4)-1.D+0)  
      RETURN
    3 GF(1)=4.D+2*(X(1)**3-X(1)*X(2))+2.D+0*X(1) +
     /      19.8D+0*X(4)-21.8D+0
      GF(2)=-2.D+2*X(1)**2+220.2D+0*X(2)-20.2D+0
      GF(3)=36.D+1*(X(3)**3-X(3)*X(4))+2.D+0*X(3)-2.D+0
      GF(4)=-18.D+1*X(3)**2+200.2D+0*X(4)+19.8D+0*X(1)-4.D+1
    4 RETURN
      END
C
      SUBROUTINE TP258(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,4),MODE
    1 N=4
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      X(1)=-3.D+0
      X(2)=-1.D+0
      X(3)=-3.D+0
      X(4)=-1.D+0
      DO 6 I=1,4
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      DO 7 I=1,4
    7 XEX(I)=1.D+0
      RETURN
    2 FX=(1.D+2*(X(2)-X(1)**2)**2+(1.D+0-X(1))**2 +
     1   9.D+1*(X(4)-X(3)**2)**2+(1.D+0-X(3))**2 +
     2   10.1D+0*((X(2)-1.D+0)**2+(X(4)-1.D+0)**2) +
     3   19.8D+0*(X(2)-1.D+0)*(X(4)-1.D+0))
      RETURN
    3 GF(1)=(4.D+2*(X(1)**3-X(1)*X(2))+2.D+0*X(1)-2.D+0)
      GF(2)=(-2.D+2*X(1)**2+220.2D+0*X(2)+19.8D+0*X(4)-4.D+1)
      GF(3)=(36.D+1*(X(3)**3-X(3)*X(4))+2.D+0*X(3)-2.D+0)
      GF(4)=(-18.D+1*X(3)**2+200.2D+0*X(4)+19.8D+0*X(2)-4.D+1)
    4 RETURN
      END
C
      SUBROUTINE TP259(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,4),MODE
    1 N=4
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,4
      X(I)=0.D+0
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LXU(4)=.TRUE.
      XU(4)=1.D+0
      LEX=.FALSE.
      NEX=2
      FEX=-0.85446210D+01
      XEX(1)=0.14358451D+01
      XEX(2)=0.20631635D+01
      XEX(3)=0.69002268D-01 
      XEX(4)=-0.99963939D-01
c      DO 7 I=1,4
c    7 XEX(N+I)=1.D+0
      RETURN
    2 FX=1.D+2*(X(2)-X(1)**2)**2+(1.D+0-X(1))**2 +
     1   9.D+1*(X(4)-X(3)**2)**2+(1.D+0-X(3))**3 +
     2   10.1D+0*(X(2)-1.D+0)**2+(X(4)-1.D+0)**2 +
     3   19.8D+0*(X(2)-1.D+0)*(X(4)-1.D+0)  
      RETURN
    3 GF(1)=4.D+2*(X(1)**3-X(1)*X(2))+2.D+0*X(1)-2.D+0
      GF(2)=-2.D+2*X(1)**2+220.2D+0*X(2)+19.8D+0*X(4)-4.D+1
      GF(3)=36.D+1*(X(3)**3-X(3)*X(4))-3.D+0*(1.D+0-X(3))**2
      GF(4)=-18.D+1*X(3)**2+182.D+0*X(4)+19.8D+0*X(2)-21.8D+0
    4 RETURN
      END
C
      SUBROUTINE TP260(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF
      INTEGER I,J
      GOTO (1,2,3,4,4),MODE
    1 N=4
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      X(1)=-3.D+0
      X(2)=-1.D+0
      X(3)=-3.D+0
      X(4)=-1.D+0
      DO 6 I=1,4
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      DO 7 I=1,4
    7 XEX(I)=1.D+0
      DO 8 I=1,7
      DO 8 J=1,4
    8 DF(I,J)=0.D+0
      DF(1,2)=10.D+0
      DF(2,1)=-1.D+0
      DF(3,4)=DSQRT(9.D+1)
      DF(4,3)=-1.D+0
      DF(5,2)=DSQRT(9.9D+0)
      DF(5,4)=DSQRT(9.9D+0)
      DF(6,2)=DSQRT(.2D+0)
      DF(7,4)=DSQRT(.2D+0)
      RETURN
    2 F(1)=10.D+0*(X(2)-X(1)**2)
      F(2)=1.D+0-X(1)
      F(3)=DSQRT(9.D+1)*(X(4)-X(3)**2)
      F(4)=1.D+0-X(3)
      F(5)=DSQRT(9.9D+0)*((X(2)-1.D+0)+(X(4)-1.D+0))
      F(6)=DSQRT(.2D+0)*(X(2)-1.D+0)
      F(7)=DSQRT(.2D+0)*(X(4)-1.D+0)
      IF (MODE.EQ.3) GOTO 3
      FX=0.D+0
      DO 9 I=1,7
    9 FX=FX+F(I)**2
      RETURN
    3 DF(1,1)=-20.D+0*X(1)
      DF(3,3)=-DSQRT(9.D+1)*2.0*X(3)
      DO 10 I=1,4
      GF(I)=0.D+0
      DO 10 J=1,7
   10 GF(I)=GF(I)+2.D+0*F(J)*DF(J,I)      
    4 RETURN
      END
C
      SUBROUTINE TP261(MODE) 
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION  F,DF,A,B,C
      INTEGER I,J
      GOTO (1,2,2,4,4),MODE 
 1    N=4 
      NILI=0 
      NINL=0 
      NELI=0 
      NENL=0 
      DO 6 I=1,4 
      XL(I)=0.0
      LXL(I)=.TRUE.
      XU(I)=1.0D+1
      LXU(I)=.TRUE.
      X(I)=0.D+0 
      DO 6 J=1,5 
 6    DF(J,I)=0.D+0 
      DF(5,4)=0.1D+1 
      LEX=.TRUE.
      NEX=1 
      FEX=0.D+0 
      XEX(1)=0.D+0 
      DO 7 I=2,4 
 7    XEX(I)=0.1D+1 
      RETURN 
 2    CONTINUE
      F(1)=(DEXP(X(1))-X(2))**2 
      F(2)=0.1D+2*(X(2)-X(3))**3 
      F(3)=DTAN(X(3)-X(4))**2
      F(4)=X(1)**4
      F(5)=X(4)-0.1D+1
      IF (MODE.EQ.3) GOTO 3
      FX=0.D+0 
      DO 8 I=1,5 
 8    FX=FX+F(I)**2 
      RETURN 
 3    A=DEXP(X(1))-X(2) 
      B=DTAN(X(3)-X(4)) 
      C=B/(DCOS(X(3)-X(4)))**2
      DF(1,1)=0.2D+1*DEXP(X(1))*A 
      DF(1,2)=-0.2D+1*A 
      DF(2,2)=0.3D+2*(X(2)-X(3))**2
      DF(2,3)=-DF(2,2) 
      DF(3,3)=0.2D+1*C 
      DF(3,4)=-DF(3,3) 
      DF(4,1)=0.4D+1*X(1)**3
      GF(1)=0.4D+1*DEXP(X(1))*A**3+0.8D+1*X(1)**7 
      GF(2)=-0.4D+1*A**3+0.6D+3*(X(2)-X(3))**5
      GF(3)=0.4D+1*B*B*C-0.6D+3*(X(2)-X(3))**5 
      GF(4)=-0.4D+1*B*B*C+0.2D+1*(X(4)-0.1D+1) 
 4    RETURN 
      END 
C
      SUBROUTINE TP262(MODE) 
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION HGG(4,4) 
      DATA HGG/-1.D+0,-0.2D+0,-2.D+0,1.D+0,-1.D+0,-0.5D+0,-1.D+0,1.D+0,
     F     	2*-1.D+0,-0.5D+0,1.D+0,-1.D+0,-2.D+0,-0.2D+0,-2.D+0/
      GOTO (1,2,3,4,3),MODE 
 1    N=4 
      NILI=3 
      NINL=0 
      NELI=1 
      NENL=0 
      DO 6 I=1,4 
      X(I)=0.1D+1 
      LXL(I)=.TRUE.
      LXU(I)=.FALSE.
      XL(I) =0.D+0 
      DO 6 J=1,4
      GG(I,J)=HGG(I,J)
 6    CONTINUE
      DO 7 I=1,3,2
      GF(I)=-0.5D+0 
 7    GF(I+1)=-0.1D+1 
      LEX=.TRUE.
      NEX=1 
      FEX=-0.1D+2 
      XEX(1)=0.D+0 
      XEX(2)=0.26D+2/0.3D+1 
      XEX(3)=0.D+0 
      XEX(4)=0.4D+1/0.3D+1 
      RETURN 
 2    FX=-0.5D+0*X(1)-X(2)-0.5D+0*X(3)-X(4) 
 3    RETURN 
 4    IF (INDEX1(1)) G(1)=0.1D+2-X(1)-X(2)-X(3)-X(4) 
      IF (INDEX1(2)) G(2)=0.1D+2-0.2D+0*X(1)-0.5D+0*X(2)
     1                     -X(3)-0.2D+1*X(4) 
      IF (INDEX1(3)) G(3)=0.1D+2-0.2D+1*X(1)-X(2) 
     1                     -0.5D+0*X(3)-0.2D+0*X(4) 
      IF (INDEX1(4)) G(4)=X(1)+X(2)+X(3)-0.2D+1*X(4)-0.6D+1
      RETURN 
      END 
C
      SUBROUTINE TP263(MODE) 
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE 
 1    N=4 
      NILI=0 
      NINL=2 
      NELI=0 
      NENL=2 
      DO 6 I=1,4 
      X(I)=0.1D+2 
      LXL(I)=.FALSE.
 6    LXU(I)=.FALSE.
      GG(1,2)=0.1D+1 
      GG(1,3)=0.D+0 
      GG(2,2)=-0.1D+1 
      GG(2,3)=0.D+0 
      GG(3,2)=0.1D+1 
      GG(4,2)=-0.1D+1 
      GG(4,3)=0.D+0 
      GF(1)=-0.1D+1 
      DO 7 I=1,3  
      GG(I,4)=0.D+0 
 7    GF(I+1)=0.D+0 
      LEX=.TRUE. 
      NEX=1 
      FEX=-0.1D+1 
      DO 8 I=1,2 
      XEX(I)=0.1D+1 
 8    XEX(I+2)=0.D+0 
      RETURN 
 2    FX=-X(1) 
 3    RETURN 
 4    IF (INDEX1(1)) G(1)=X(2)-X(1)**3
      IF (INDEX1(2)) G(2)=X(1)**2-X(2) 
      IF (INDEX1(3)) G(3)=X(2)-X(1)**3-X(3)**2
      IF (INDEX1(4)) G(4)=X(1)**2-X(2)-X(4)**2
      RETURN 
 5    IF (INDEX2(1)) GG(1,1)=-0.3D+1*X(1)**2
      IF (INDEX2(2)) GG(2,1)=0.2D+1*X(1) 
      IF (.NOT.INDEX2(3)) GOTO 9 
      GG(3,1)=-0.3D+1*X(1)**2
      GG(3,3)=-0.2D+1*X(3) 
 9    IF (.NOT.INDEX2(4)) GOTO 10 
      GG(4,1)=0.2D+1*X(1) 
      GG(4,4)=-0.2D+1*X(4) 
10    RETURN 
      END 
C
      SUBROUTINE TP264(MODE) 
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE 
 1    N=4 
      NILI=0 
      NINL=3 
      NELI=0 
      NENL=0 
      DO 6 I=1,4 
      X(I)=0.D+0 
      LXL(I)=.FALSE.
 6    LXU(I)=.FALSE.
      GG(3,4)=0.1D+1 
      LEX=.TRUE.
      NEX=1 
      FEX=-0.44D+2 	 
      XEX(1)=0.D+0 
      XEX(2)=0.1D+1 
      XEX(3)=0.2D+1 
      XEX(4)=-0.1D+1 
      RETURN 
 2    FX=(X(1)**2+X(2)**2+0.2D+1*X(3)**2+X(4)**2-0.5D+1*X(1)
     /       -0.5D+1*X(2)-0.21D+2*X(3)+0.7D+1*X(4))
      RETURN 
 3    GF(1)=(0.2D+1*X(1)-0.5D+1)    	       
      GF(2)=(0.2D+1*X(2)-0.5D+1)      	       
      GF(3)=(0.4D+1*X(3)-0.21D+2)      	       
      GF(4)=(0.2D+1*X(4)+0.7D+1) 
      RETURN 
 4    IF (INDEX1(1)) G(1)=0.8D+1-X(1)**2-X(2)**2-X(3)**2 
     /    -X(4)**2-X(1)+X(2)-X(3)+X(4)  
      IF (INDEX1(2)) G(2)=0.9D+1-X(1)**2-0.2D+1*X(2)**2
     /    -X(3)**2-0.2D+1*X(4)**2+X(1)+X(4) 
      IF (INDEX1(3)) G(3)=0.5D+1-0.2D+1*X(1)**2-X(2)**2
     /    -X(3)**2-0.2D+1*X(1)+X(2)+X(4) 
      RETURN 
 5    IF (.NOT.INDEX2(1)) GOTO 9 
      GG(1,1)=-0.1D+1-0.2D+1*X(1) 
      GG(1,2)=0.1D+1-0.2D+1*X(2) 
      GG(1,3)=-0.1D+1-0.2D+1*X(3) 
      GG(1,4)=0.1D+1-0.2D+1*X(4) 
 9    IF (.NOT.INDEX2(2)) GOTO 10 
      GG(2,1)=0.1D+1-0.2D+1*X(1) 
      GG(2,2)=-0.4D+1*X(2) 
      GG(2,3)=-0.2D+1*X(3) 
      GG(2,4)=-0.1D+1-0.4D+1*X(4)   
10    IF (.NOT.INDEX2(3)) GOTO 11 
      GG(3,1)=-0.2D+1-0.4D+1*X(1) 
      GG(3,2)=0.1D+1-0.2D+1*X(2) 
      GG(3,3)=-0.2D+1*X(3) 
11    RETURN 
      END       
C
      SUBROUTINE TP265(MODE) 
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION HGG(2,4)
      DATA HGG /1.D+0,0.D+0,1.D+0,0.D+0,0.D+0,1.D+0,0.D+0,1.D+0/
      GOTO (1,2,3,4,5),MODE 
 1    N=4 
      NILI=0 
      NINL=0 
      NELI=2 
      NENL=0 
      DO 6 I=1,4 
      X(I)=0.D+0
      LXL(I)=.TRUE.
      XL(I)=0.D+0 
      LXU(I)=.FALSE.
      XU(I)=1.D+0
      DO 6 J=1,2
      GG(J,I)=HGG(J,I) 
 6    CONTINUE   
      LEX=.TRUE.
      NEX=2
      FEX=0.97474658D+0
      DO 7 I=1,3,2
      XEX(I)=1.D+0
 7    XEX(I+1)=0.D+0
      DO 8 I=5,8
 8    XEX(I)=XEX(9-I)
      RETURN 
 2    FX=2.D+0 
      DO 9 I=1,2 
 9    FX=FX-DEXP(-0.1D+2*X(I)*DEXP(-X(I+2))) 
      RETURN 
 3    DO 10 I=1,2 
      GF(I)=0.1D+2*DEXP(-0.1D+2*X(I)*DEXP(-X(I+2))-X(I+2)) 
10    GF(I+2)=-X(I)*GF(I) 
      RETURN 
 4    IF (INDEX1(1)) G(1)=X(1)+X(2)-0.1D+1 
      IF (INDEX1(2)) G(2)=X(3)+X(4)-0.1D+1 
 5    RETURN 
      END 
C
      SUBROUTINE TP266(MODE) 
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,K,L
      DOUBLEPRECISION F,DF,A(10),D(10),C(10,5),B(5,5),HF 
      DATA A /0.426149D-1,0.352053D-1,0.878058D-1,0.330812D-1,
     /        0.580924D-1,
     F        0.649704D0,0.344144D0,-0.627443D0,0.1828D-2,-0.224783D0/
      DATA D /0.234659D+1,0.284048D+1,0.113888D+1,0.302286D+1,
     F        0.172139D+1,0.153917D+0,0.290577D+0,-0.159378D+0,
     F        0.546910D+2,-0.444873D+0/
      DATA C /-0.564255D+0,0.535493D+0,0.586387D+0,0.608734D+0,
     /    0.774227D+0,
     F   -0.435033D+0,0.759468D+0,-0.152448D+0,-0.821772D+0,0.819831D+0,
     F   .0392417D+0,0.658799D+0,0.289826D+0,0.984915D+0,0.325421D+0,
     F   -0.688583D+0,-0.627795D+0,-0.546437D+0,-0.53412D0,-0.910632D0,
     F   -0.404979D0,-0.0636666D0,0.854402D0,0.375699D0,-0.151719D0,
     F   .0222278D+0,0.0403142D+0,0.484134D+0,-0.798498D+0,-0.480344D+0,
     F   0.927589D+0,-0.681091D+0,0.789312D+0,0.239547D+0,0.448051D+0,
     F   -0.524653D+0,0.724666D+0,0.353951D+0,-0.658572D+0,-0.871758D+0,
     F   -0.0735083D+0,-0.869487D+0,0.949721D+0,0.463136D+0,0.149926D+0,
     F   0.413248D+0,-0.0182537D+0,0.887866D+0,0.662362D+0,-0.978666D+0/
      DATA B /.354033D+0,-0.0230349D+0,-0.211938D+0,-0.0554288D+0,
     F   0.220429D+0,-0.0230349D+0,0.29135D+0,-0.00180333D0,-0.111141D0,
     F   0.0485461D+0,-0.211938D0,-0.00180333D0,-0.815808D0,-0.133538D0,
     F   -0.38067D+0,-0.0554288D+0,-0.111141D+0,-0.133538D+0,0.389198D0,
     F   -0.131586D+0,0.220429D+0,0.0485461D+0,-0.38067D+0,-0.131586D0, 
     F   0.534706D+0/
      GOTO (1,2,2,4,4),MODE 
 1    N=5 
      NILI=0 
      NINL=0 
      NELI=0 
      NENL=0 
      DO 6 I=1,5 
      X(I)=0.1D+0 
      LXU(I)=.FALSE.
      LXL(I)=.TRUE.
      XL(I)=0.0D0
 6    XEX(I)=0.D+0 
      XEX(3)=0.29297857D-1
      NEX=1 
      LEX=.TRUE.
      NEX=1
      FEX= 0.99597447D+0 
      RETURN 
 2    DO 20 I=1,10
      CALL TP266A(A(I),B,C,D(I),I,X,F(I))
 20   CONTINUE
      IF (MODE.EQ.3) GOTO 3 
      FX=0.D+0 
      DO 7 I=1,10 
 7    FX=FX+F(I)**2
      RETURN 
 3    DO 8 K=1,5 
      GF(K)=0.D+0 
      DO 8 I=1,10 
      HF=0.D+0 
      DO 10 L=1,5 
10    HF=HF+(B(K,L)+B(L,K))*X(L) 
      DF(I,K)=C(I,K)+.5D+0*D(I)*HF 
      GF(K)=GF(K)+0.2D+1*F(I)*DF(I,K) 
 8    CONTINUE
 4    RETURN 
      END 
C      
      SUBROUTINE TP266A(AI,B,C,DI,I,X,TP) 
C 
C  FUNKTION ZUR BERECHNUNG VON  
C      F(I,X)=(A+C*X+0.5*X'*B*X*D)(I) 
C  D.H. DAS I-TE ELEMENT DES DOUBLEPRECISION VEKTORS F(10) FUER I=1,..,10 
C 
      DOUBLE PRECISION AI,C(10,5),B(5,5),DI,X(5),HF,TP
      INTEGER I,K,L
      TP=AI 
      DO 10 K=1,5 
      HF=0.0D0 
      DO 20 L=1,5 
 20   HF=HF+B(K,L)*X(L) 
      TP=TP+X(K)*(C(I,K)+0.5D0*DI*HF) 
 10   CONTINUE 
      RETURN 
      END 
C
      SUBROUTINE TP267(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION H,F,DF
      GOTO (1,2,2,4,4),MODE
    1 N=5
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,5
      X(I)=2.D+0
      LXL(I)=.FALSE.
      XU(I)=15.0
    6 LXU(I)=.TRUE.
      LXL(1)=.TRUE.
      LXL(2)=.TRUE.
      LXL(5)=.TRUE.
      XL(1)=0.0
      XL(2)=0.0
      XL(5)=0.0
      LEX=.TRUE.
      NEX=2
      FEX=0.D+0
      XEX(1)=1.D+0
      XEX(2)=1.D+1
      XEX(3)=1.D+0
      XEX(4)=5.D+0
      XEX(5)=4.D+0
      XEX(6)=1.D+1
      XEX(7)=1.D+0
      XEX(8)=-5.D+0
      XEX(9)=-1.D+0
      XEX(10)=4.D+0
      RETURN
    2 DO 20 I=1,11
      H=.1D+0*DBLE(I)
   20 F(I)=X(3)*DEXP(-X(1)*H)-X(4)*DEXP(-X(2)*H)+3.D+0*DEXP(-X(5)*H)
     F      -(DEXP(-H)-5.D+0*DEXP(-1.D+1*H)+3.D+0*DEXP(-4.D+0*H))
      IF (MODE.EQ.3) GOTO 3 
      FX=0.D+0
      DO 7 I=1,11
    7 FX=FX+F(I)**2
      RETURN
    3 DO 8 J=1,11
      H=.1D+0*DBLE(J)
      DF(J,3)=DEXP(-X(1)*H)
      DF(J,4)=-DEXP(-X(2)*H)
      DF(J,1)=-H*X(3)*DF(J,3)
      DF(J,2)=-H*X(4)*DF(J,4)
    8 DF(J,5)=-H*3.D+0*DEXP(-X(5)*H)
      DO 13 I=1,5
      GF(I)=0.D+0
      DO 13 J=1,11
   13 GF(I)=GF(I)+2.D+0*F(J)*DF(J,I)
    4 RETURN
      END
C
      SUBROUTINE TP268(MODE) 
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION HF
      INTEGER HGG(5,5),DD(5,5),DDVEKT(5),DVDV 
      DATA HGG /-1,10,-8,8,-4,-1,10,1,-1,-2,-1,-3,-2,2,3,-1,5, 
     1          -5,5,-5,-1,4,3,-3,1/
C KONSTANTE DATEN DER ZIELFUNKTION:
C         DD=D'*D
C         DDVEKT=DVEKT'*D
C         DVDV=DVEKT'*DVEKT=14463
C MIT
C          -                      -
C          |  -74   80   18  -11  -4  |
C          |   14  -69   21   28   0  |
C     D= |   66  -72   -5    7   1  |
C          |  -12   66  -30  -23   3  |
C          |    3    8   -7   -4   1  |
C          |    4  -12    4    4   0  |
C          -                       -
C
C     DVEKT=( 51, -61, -56, 69, 10, -12 )'
      DATA DD /10197,-12454,-1013,1948,329,-12454,20909,-1733,-4914,
     1         -186,-1013,-1733,1755,1089,-174,1948,-4914,1089,1515,
     2         -22,329,-186,-174,-22,27/
      DATA DDVEKT,DVDV/-9170,17099,-2271,-4336,-43,14463/
      GOTO (1,2,3,4,5),MODE 
 1    N=5 
      NILI=5 
      NINL=0 
      NELI=0 
      NENL=0 
      DO 6 I=1,5 
      X(I)=0.1D+1 
      LXL(I)=.FALSE.
      LXU(I)=.FALSE.
      DO 6 J=1,5 
 6    GG(I,J)=DBLE(HGG(I,J)) 
      LEX=.TRUE.
      NEX=1 
      FEX=0.D+0 
      XEX(1)=0.1D+1 
      XEX(2)=0.2D+1 
      XEX(3)=-0.1D+1 
      XEX(4)=0.3D+1 
      XEX(5)=-0.4D+1 
      RETURN 
 2    FX=DVDV 
      DO 7 I=1,5 
      HF=0.D+0
      DO 8 J=1,5 
 8    HF=HF+DBLE(DD(I,J))*X(J) 
 7    FX=FX+X(I)*(HF-0.2D+1*DBLE(DDVEKT(I)))
      RETURN 
 3    DO 9 I=1,5 
      GF(I)=- 0.2D+1*DBLE(DDVEKT(I)) 
      DO 9 J=1,5 
 9    GF(I)=GF(I)+DBLE(DD(I,J)+DD(J,I))*X(J) 
      RETURN 
 4    IF (INDEX1(1)) G(1)=-X(1)-X(2)-X(3)-X(4)-X(5)+5.0D0 
      IF (INDEX1(2)) G(2)=0.1D+2*X(1)+0.1D+2*X(2) 
     1  -0.3D+1*X(3)+0.5D+1*X(4)+0.4D+1*X(5)-0.2D+2 
      IF (INDEX1(3)) G(3)=-0.8D+1*X(1)+X(2)-0.2D+1*X(3) 
     1  -0.5D+1*X(4)+0.3D+1*X(5)+0.4D+2 
      IF (INDEX1(4)) G(4)=0.8D+1*X(1)-X(2)+0.2D+1*X(3) 
     1  +0.5D+1*X(4)-0.3D+1*X(5)-0.11D+2 
      IF (INDEX1(5)) G(5)=-0.4D+1*X(1)-0.2D+1*X(2) 
     1  +0.3D+1*X(3)-0.5D+1*X(4)+X(5)+0.3D+2 
 5    RETURN 
      END 
C
      SUBROUTINE TP269(MODE) 
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION  F,DF,HGG(3,5),HDF(4,5) 
      DATA HGG/1.D+0,2*0.D+0,3.D+0,0.D+0,1.D+0,0.D+0,1.D+0,2*0.D+0,
     F         1.D+0,2*0.D+0,-2.D+0,-1.D+0/
      DATA HDF/1.D+0,3*0.D+0,-1.D+0,1.D+0,3*0.D+0,1.D+0,4*0.D+0,
     F         1.D+0,4*0.D+0,1.D+0/
      GOTO (1,2,2,4,5),MODE 
 1    N=5 
      NILI=0 
      NINL=0 
      NELI=3 
      NENL=0 
      DO 6 I=1,5 
      X(I)=2.D+0
      LXL(I)=.FALSE.
      LXU(I)=.FALSE.
      DO 7 J=1,3 
      GG(J,I)=HGG(J,I) 
 7    DF(J,I)=HDF(J,I) 
 6    DF(4,I)=HDF(4,I) 
      LEX=.TRUE. 
      NEX=1 
      FEX=0.176D+3/0.43D+2 
      XEX(1)=-0.33D+2/0.43D+2 
      XEX(2)=0.11D+2/0.43D+2 
      XEX(3)=0.27D+2/0.43D+2 
      XEX(4)=-0.5D+1/0.43D+2 
      XEX(5)=0.11D+2/0.43D+2 
      RETURN 
 2    F(1)=X(1)-X(2) 
      F(2)=X(2)+X(3)-0.2D+1 
      F(3)=X(4)-0.1D+1 
      F(4)=X(5)-0.1D+1 
      IF (MODE.EQ.3) GOTO 3
      FX=0.D+0 
      DO 8 I=1,4 
 8    FX=FX+F(I)**2
      RETURN 
 3    GF(1)=0.2D+1*(X(1)-X(2)) 
      GF(2)=0.2D+1*(0.2D+1*X(2)+X(3)-X(1)-0.2D+1) 
      GF(3)=0.2D+1*(X(2)+X(3)-0.2D+1) 
      GF(4)=0.2D+1*(X(4)-0.1D+1) 
      GF(5)=0.2D+1*(X(5)-0.1D+1) 
      RETURN 
 4    IF (INDEX1(1)) G(1)=X(1)+0.3D+1*X(2) 
      IF (INDEX1(2)) G(2)=X(3)+X(4)-0.2D+1*X(5) 
      IF (INDEX1(3)) G(3)=X(2)-X(5) 
 5    RETURN	      
      END 
      SUBROUTINE TP270(MODE) 
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE 
 1    N=5 
      NILI=0 
      NINL=1 
      NELI=0 
      NENL=0 
      DO 6 I=1,4 
      X(I)=DBLE(I)+0.1D+0 
      XEX(I)=DBLE(I) 
      LXL(I)=.TRUE. 
      LXU(I)=.FALSE. 
 6    XL(I)=DBLE(I) 
      LXL(5)=.FALSE. 
      LXU(5)=.FALSE. 
      X(5)=-0.1D+1 
      XEX(5)=0.2D+1 
      NEX=1 
      LEX=.TRUE. 
      FEX= -0.1D+1 
      RETURN 
 2    FX=X(1)*X(2)*X(3)*X(4)-3.D+0*X(1)*X(2)*X(4)
     F -4.D+0*X(1)*X(2)*X(3)+12.D+0*X(1)*X(2)-X(2)*X(3)*X(4)
     F +3.D+0*X(2)*X(4)+4.D+0*X(2)*X(3)-12.D+0*X(2)
     F -2.D+0*X(1)*X(3)*X(4) +6.D+0*X(1)*X(4)+8.D+0*X(1)*X(3)
     F -24.D+0*X(1)+2.D+0*X(3)*X(4)-6.D+0*X(4)-8.D+0*X(3)
     F +24.D+0+1.5D+0*X(5)**4-5.75D+0*X(5)**3+5.25D+0*X(5)**2
      RETURN 
 3    GF(1)=X(2)*X(3)*X(4)-3.D+0*X(2)*X(4)-4.D+0*X(2)*X(3)   
     F +12.D+0*X(2)-2.D+0*X(3)*X(4)+6.D+0*X(4)+8.D+0*X(3)-24.D+0  
      GF(2)=X(1)*X(3)*X(4)-3.D+0*X(1)*X(4)-4.D+0*X(1)*X(3) 
     F +12.D+0*X(1)-X(3)*X(4)+3.D+0*X(4)+4.D+0*X(3)-12.D+0
      GF(3)=X(1)*X(2)*X(4)-4.D+0*X(1)*X(2)-X(2)*X(4)+4.D+0*X(2)
     F -2.D+0*X(1)*X(4)+8.D+0*X(1)+2.D+0*X(4)-8.D+0  
      GF(4)=X(1)*X(2)*X(3)-3.D+0*X(1)*X(2)-X(2)*X(3)+3.D+0*X(2)
     F -2.D+0*X(1)*X(3)+6.D+0*X(1)+2.D+0*X(3)-6.D+0 
      GF(5)=10.5D+0*X(5)-17.25D+0*X(5)**2+6.D+0*X(5)**3 
      RETURN 
 4    IF (INDEX1(1)) G(1)=34.D+0-X(1)**2-X(2)**2-X(3)**2
     F             -X(4)**2-X(5)**2 
      RETURN 
 5    IF (.NOT.INDEX2(1)) GOTO 7 
      DO 8 I=1,5 
 8    GG(1,I)=-0.2D+1*X(I) 
 7    RETURN 
      END 
C      
      SUBROUTINE TP271(MODE) 
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION  F,DF
      INTEGER I,J,LSUM
      GOTO (1,2,2,4,4),MODE 
 1    N=6 
      LSUM=6 
      NILI=0 
      NINL=0 
      NELI=0 
      NENL=0 
      DO 6 I=1,6 
      LXL(I)=.FALSE. 
      LXU(I)=.FALSE. 
      X(I)=0.D+0 
 6    XEX(I)=0.1D+1 
      LEX=.TRUE. 
      NEX=1
      FEX=0.D+0 
      RETURN 
 2    DO 10 I=1,6
 10   F(I)=DSQRT(DBLE(10*(16-I)))*(X(I)-1.D+0)
      IF (MODE.EQ.3) GOTO 3
      FX=0.D+0 
      DO 7 I=1,6 
 7    FX=FX+F(I)**2
      RETURN 
 3    DO 8 I=1,6 
      GF(I)=DBLE(20*(16-I))*(X(I)-0.1D+1) 
      DO 9 J=1,6 
 9    DF(I,J)=0.D+0 
 8    DF(I,I)=DSQRT(DBLE(10*(16-I))) 
 4    RETURN 
      END
C      
      SUBROUTINE TP272(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION H,F,DF
      INTEGER I,J
      GOTO (1,2,2,4,4),MODE
    1 N=6
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,6
      X(I)=1.D+0
      XL(I)=0.0
      LXL(I)=.TRUE.
    6 LXU(I)=.FALSE.
      X(2)=2.D+0
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      XEX(1)=1.D+0
      XEX(2)=10.D+0
      XEX(3)=4.D+0
      XEX(4)=1.D+0
      XEX(5)=5.D+0
      XEX(6)=3.D+0
      RETURN
    2 DO 20 I=1,13
      H=.1D+0*DBLE(I)
   20 F(I)=X(4)*DEXP(-X(1)*H)-X(5)*DEXP(-X(2)*H)+X(6)*DEXP(-X(3)*H)
     F     -DEXP(-H)+5.D+0*DEXP(-1.D+1*H)-3.D+0*DEXP(-4.D+0*H)
      IF (MODE.EQ.3) GOTO 3
      FX=0.D+0
      DO 7 I=1,13
    7 FX=FX+F(I)**2
      RETURN
    3 DO 8 J=1,13
      H=.1D+0*DBLE(J)
      DF(J,4)=DEXP(-X(1)*H)
      DF(J,5)=-DEXP(-X(2)*H)
      DF(J,6)=DEXP(-X(3)*H)
      DO 8 I=1,3
    8 DF(J,I)=-H*X(I+3)*DF(J,I+3)
      DO 13 I=1,6
      GF(I)=0.D+0
      DO 13 J=1,13
   13 GF(I)=GF(I)+2.D+0*DF(J,I)*F(J)
    4 RETURN
      END
C      
      SUBROUTINE TP273(MODE) 
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION HX,TP273A
      INTEGER I
      GOTO (1,2,3,4,4)MODE 
 1    N=6 
      NILI=0 
      NINL=0 
      NELI=0 
      NENL=0 
      DO 6 I=1,6 
      X(I)=0.D+0 
      XEX(I)=0.1D+1 
      LXL(I)=.FALSE. 
6     LXU(I)=.FALSE. 
      LEX=.TRUE. 
      NEX=1 
      FEX=0.D+0 
      RETURN 
 2    HX=TP273A(X) 
      FX=0.1D+2*HX*(0.1D+1+HX) 
      RETURN 
 3    HX=TP273A(X) 
      DO 7 I=1,6 
 7    GF(I)=0.2D+2*(0.16D+2-DBLE(I))*(X(I)-0.1D+1) 
     1       *(0.1D+1+0.2D+1*HX) 
 4    RETURN 
      END       
C      
      DOUBLEPRECISION FUNCTION TP273A (X) 
      DOUBLEPRECISION X(6)
      INTEGER I
C 
C  BERECHNUNG VON H(X) IN TP273 
C  HX=SUM( (16-I)*(X(I)-1)**2,I=1,..,6)
C 
      TP273A=0 
      DO 10 I=1,6 
10    TP273A=TP273A+(0.16D+2-DBLE(I))*(X(I)-0.1D+1)**2 
      RETURN 
      END
C      
      SUBROUTINE TP274(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D274/A(6,6)
      DOUBLEPRECISION A
      INTEGER I,J
      N=2
      GOTO 10
      ENTRY TP275(MODE)
      N=4
      GOTO 10
      ENTRY TP276(MODE)
      N=6
   10 GOTO (1,2,3,4,4), MODE
    1 NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,N
      X(I)=-4.D+0/DBLE(I)
      DO 11 J=1,N
   11 A(I,J)=1.D+0/DBLE(I+J-1)
      XEX(I)=0.D+0
      LXL(I)=.FALSE.
    6 LXU(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      RETURN
    2 FX=0.D+0
      DO 7 I=1,N
      DO 7 J=1,N
    7 FX=FX+A(I,J)*X(I)*X(J)
      RETURN
    3 DO 8 I=1,N
      GF(I)=0.D+0
      DO 8 J=1,N
    8 GF(I)=GF(I)+X(J)*(A(I,J)+A(J,I))
    4 RETURN
      END
C      
      SUBROUTINE TP277(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION H
      N=4
      GOTO 10
      ENTRY TP278(MODE)
      N=6
      GOTO 10
      ENTRY TP279(MODE)
      N=8
      GOTO 10
      ENTRY TP280(MODE)
      N=10
   10 GOTO (1,2,3,4,5), MODE
    1 NILI=N
      NINL=0
      NELI=0
      NENL=0
      FEX=0.D+0
      DO 6 I=1,N
      X(I)=0.D+0
      DO 16 J=1,N
      FEX=FEX+1.D+0/DBLE(I+J-1)
   16 GG(J,I)=1.D+0/DBLE(I+J-1)
      LXL(I)=.TRUE.
      LXU(I)=.FALSE.
      XL(I)=0.D+0
    6 XEX(I)=1.D+0
      LEX=.TRUE.
      NEX=1
      RETURN
    2 FX=0.D+0
      DO 7 I=1,N
      H=0.D+0
      DO 17 J=1,N
   17 H=H+1.D+0/DBLE(I+J-1)
    7 FX=FX+H*X(I)
      RETURN
    3 DO 8 I=1,N
      H=0.D+0
      DO 18 J=1,N
   18 H=H+1.D+0/DBLE(I+J-1)
    8 GF(I)=H
      RETURN
    4 DO 9 I=1,N
      H=0.D+0
      DO 19 J=1,N
   19 H=H+(X(J)-1.D+0)/DBLE(I+J-1)
    9 IF (INDEX1(I)) G(I)=H
    5 RETURN
      END
C
      SUBROUTINE TP281(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLE PRECISION H
      GOTO (1,2,3,4,4),MODE
    1 N=10
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,10
      X(I)=0.D+0
      LXL(I)=.FALSE.
      LXU(I)=.FALSE.
    6 XEX(I)=1.D+0
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      RETURN
    2 FX=0.D+0
      DO 7 I=1,10
    7 FX=FX+DBLE(I**3)*(X(I)-1.0D0)**2
      FX=FX**(1.D+0/3.0D0)
      RETURN
    3 H=0.D+0
      DO 8 I=1,10
    8 H=H+DBLE(I**3)*(X(I)-1.D+0)**2
      DO 13 I=1,10
   13 GF(I)=(2.D+0/3.D+0)*DBLE(I**3)*(X(I)-1.D+0)
     1                                     *(H**(-2.D+0/3.D+0))
    4 RETURN
      END
C
      SUBROUTINE TP282(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION H,F,DF
      GOTO (1,2,2,4,4), MODE
    1 N=10
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,10
      X(I)=0.D+0
      LXL(I)=.FALSE.
      LXU(I)=.FALSE.
    6 XEX(I)=1.D+0
      X(1)=-1.2D+0
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      RETURN
    2 F(10)=X(1)-1.D+0
      F(11)=X(10)-1.D+0
      DO 20 I=1,9
   20 F(I)=DSQRT((1.D+1-DBLE(I))*1.D+1)*(X(I)**2-X(I+1))
      IF (MODE.EQ.3) GOTO 3
      FX=F(10)**2+F(11)**2
      DO 7 I=1,9
    7 FX=FX+F(I)**2
      RETURN
    3 DO 8 I=1,11
      DO 8 J=1,10
    8 DF(I,J)=0.D+0
      DO 13 I=1,9
      H=DSQRT((1.D+1-DBLE(I))*1.D+1)
      DF(I,I)=2.D+0*H*X(I)
   13 DF(I,I+1)=-H
      DF(10,1)=1.D+0
      DF(11,10)=1.D+0
      DO 18 I=1,10
      GF(I)=0.D+0
      DO 18 J=1,11
   18 GF(I)=GF(I)+2.D+0*DF(J,I)*F(J)
    4 RETURN
      END
C
      SUBROUTINE TP283(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION H
      GOTO (1,2,3,4,4),MODE
    1 N=10
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,10
      X(I)=0.D+0
      LXL(I)=.FALSE.
      LXU(I)=.FALSE.
    6 XEX(I)=1.D+0
      LEX=.TRUE.
      NEX=1
      FEX=0.D+0
      RETURN
    2 FX=0.D+0
      DO 7 I=1,10
    7 FX=FX+DBLE(I**3)*(X(I)-1.D+0)**2
      FX=FX**3
      RETURN
    3 H=0.D+0
      DO 8 I=1,10
    8 H=H+DBLE(I**3)*(X(I)-1.D+0)**2
      DO 13 I=1,10
   13 GF(I)=6.0*DBLE(I**3)*(X(I)-1.D+0)*H**2
    4 RETURN
      END
C
      SUBROUTINE TP284(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION SUM
      INTEGER A(10,15),B(10),C(15)
      DATA C/20,40,400,20,80,20,40,140,380,280,80,40,140,40,120/
      DATA B/385,470,560,565,645,430,485,455,390,460/
      DATA A/100,90,70,2*50,40,30,20,10,5,2*100,50,0,10,
     F 0,60,30,70,3*10,2*0,70,50,30,40,10,100,5,35,55,
     F 65,60,95,90,25,35,5,10,20,25,35,45,50,0,40,25,
     F 20,0,5,2*100,45,35,30,25,65,5,2*0,40,35,0,10,5,
     F 15,0,10,25,35,50,60,35,60,25,10,30,35,0,55,2*0,
     F 65,2*0,80,0,95,10,25,30,15,5,45,70,20,0,70,55,
     F 20,60,0,75,15,20,30,25,20,5,0,10,75,100,20,25,
     F 30,0,10,45,40,30,35,75,0,70,5,15,35,20,25,0,30,
     F 10,5,15,65,50,10,0,10,40,65,0,5,15,20,55,30/
      GOTO(1,2,3,4,5),MODE
    1 N=15
      NILI=0
      NINL=10
      NELI=0
      NENL=0
      FEX=.0D+0
      DO 6 I=1,15
      X(I)=.0D+0
      XEX(I)=.1D+1
      FEX=FEX-DBLE(C(I))
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      RETURN
    2 FX=.0D+0
      DO 7 I=1,15
    7 FX=FX-DBLE(C(I))*X(I)
      RETURN
    3 DO 8 I=1,15
    8 GF(I)=-DBLE(C(I))
      RETURN
    4 DO 10 I=1,10
      SUM=.0D+0
      DO 9 J=1,15
    9 SUM=SUM+(DBLE(A(I,J))*X(J)**2)
   10 IF (INDEX1(I)) G(I)=DBLE(B(I))-SUM
      RETURN
    5 DO 12 J=1,10
      DO 11 I=1,15
      IF (.NOT.INDEX2(J)) GOTO 12
   11 GG(J,I)=-.2D+1*DBLE(A(J,I))*X(I)
   12 CONTINUE
      RETURN
      END
C      
      SUBROUTINE TP285(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION SUM
      INTEGER A(10,15),B(10),C(15)
      DATA C/486,640,758,776,477,707,175,619,627,614,475,377,524,
     F       468,529/
      DATA B/385,470,560,565,645,430,485,455,390,460/
      DATA A/100,90,70,2*50,40,30,20,10,5,2*100,50,0,10,
     F 0,60,30,70,3*10,2*0,70,50,30,40,10,100,5,35,55,
     F 65,60,95,90,25,35,5,10,20,25,35,45,50,0,40,25,
     F 20,0,5,2*100,45,35,30,25,65,5,2*0,40,35,0,10,5,
     F 15,0,10,25,35,50,60,35,60,25,10,30,35,0,55,2*0,
     F 65,2*0,80,0,95,10,25,30,15,5,45,70,20,0,70,55,
     F 20,60,0,75,15,20,30,25,20,5,0,10,75,100,20,25,
     F 30,0,10,45,40,30,35,75,0,70,5,15,35,20,25,0,30,
     F 10,5,15,65,50,10,0,10,40,65,0,5,15,20,55,30/
      GOTO(1,2,3,4,5),MODE
    1 N=15
      NILI=0
      NINL=10
      NELI=0
      NENL=0
      FEX=.0D+0
      DO 6 I=1,15
      X(I)=.0D+0
      XEX(I)=.1D+1
      FEX=FEX-DBLE(C(I))
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      RETURN
    2 FX=.0D+0
      DO 7 I=1,15
    7 FX=FX-DBLE(C(I))*X(I)
      RETURN
    3 DO 8 I=1,15
    8 GF(I)=-DBLE(C(I))
      RETURN
    4 DO 10 I=1,10
      SUM=.0D+0
      DO 9 J=1,15
    9 SUM=SUM+(DBLE(A(I,J))*X(J)**2)
   10 IF (INDEX1(I)) G(I)=DBLE(B(I))-SUM
      RETURN
    5 DO 12 J=1,10
      DO 11 I=1,15
      IF (.NOT.INDEX2(J)) GOTO 12
   11 GG(J,I)=-.2D+1*DBLE(A(J,I))*X(I)
   12 CONTINUE
      RETURN
      END
C      
      SUBROUTINE TP286(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF
      INTEGER I,J
      GOTO(1,3,3,4,4),MODE
    1 N=20
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,20
      X(I)=-.12D+1
      XEX(I)=.1D+1	
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      DO 7 I=11,20
    7 X(I)=.1D+1
      LEX=.TRUE.
      NEX=1
      FEX=.0D+0
      RETURN
    2 FX=.0D+0
      DO 8 I=1,20
    8 FX=FX+F(I)**2
      RETURN
    3 DO 12 I=1,10
      F(I)=X(I)-.1D+1
   12 F(I+10)=.1D+2*(X(I)**2-X(I+10))
      IF (MODE.EQ.2) GOTO 2
      DO 9 I=1,10
      DO 9 J=1,20
      DF(I,J)=.0D+0
      IF (I.EQ.J) DF(I,J)=.1D+1
      DF(I+10,J)=.0D+0
      IF (I.EQ.J) DF(I+10,J)=.2D+2*X(I)
    9 IF (J.EQ.(I+10)) DF(I+10,J)=-.1D+2
      DO 11 J=1,20
      GF(J)=.0D+0
      DO 11 I=1,20
   11 GF(J)=GF(J)+.2D+1*F(I)*DF(I,J)
    4 RETURN
      END
C      
      SUBROUTINE TP287(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,4),MODE
    1 N=20
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,5
      X(I)=-.3D+1
      X(I+5)=-.1D+1
      X(I+10)=-.3D+1
   6  X(I+15)=-.1D+1
      DO 7 I=1,20
      LXU(I)=.FALSE.
      LXL(I)=.FALSE.
    7 XEX(I)=.1D+1
      LEX=.TRUE.
      NEX=1
      FEX=.0D+0
      RETURN
    2 FX=.0D+0
      DO 8 I=1,5
    8 FX=FX+.1D+3*(X(I)**2-X(I+5))**2+(X(I)-.1D+1)**2+
     F   .9D+2*(X(I+10)**2-X(I+15))**2+(X(I+10)-.1D+1)**2+
     F   .101D+2*((X(I+5)-.1D+1)**2+(X(I+15)-.1D+1)**2)+
     F   .198D+2*(X(I+5)-.1D+1)*(X(I+15)-.1D+1)
      FX=FX*1.0D-5
      RETURN
    3 DO 9 I=1,5
      GF(I)=.4D+3*(X(I)**2-X(I+5))*X(I)+.2D+1*(X(I)-.1D+1)
      GF(I)=GF(I)*1.0d-5
      GF(I+5)=-.2D+3*(X(I)**2-X(I+5))+.202D+2*(X(I+5)-.1D+1)+
     F        .198D+2*(X(I+15)-.1D+1)
      GF(I+5)=GF(I+5)*1.0d-5
      GF(I+10)=.36D+3*X(I+10)*(X(I+10)**2-X(I+15))+.2D+1*(X(I+10)-.1D+1)
      GF(I+10)=GF(I+10)
      GF(I+15)=-.18D+3*(X(I+10)**2-X(I+15))+.202D+2*(X(I+15)-.1D+1)+
     F         .198D+2*(X(I+5)-.1D+1)
    9 GF(I+15)=GF(I+15)*1.0d-5
    4 RETURN
      END
C      
      SUBROUTINE TP288(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF
      INTEGER I,J
      GOTO(1,3,3,4,4),MODE
    1 N=20
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,5
      X(I)=.3D+1
      X(I+5)=-.1D+1
      X(I+10)=.0D+0
    6 X(I+15)=.1D+1
      DO 7 I=1,20
      XEX(I)=.0D+0
      LXL(I)=.FALSE.
    7 LXU(I)=.FALSE.
      FEX=.0D+0
      LEX=.TRUE.
      NEX=1
      RETURN
    2 FX=.0D+0
      DO 9 I=1,20
    9 FX=FX+F(I)**2
      RETURN
    3 DO 8 I=1,5
      F(I)=X(I)+.1D+2*X(I+5)
      F(I+5)=DSQRT(.5D+1)*(X(I+10)-X(I+15))
      F(I+10)=(X(I+5)-.2D+1*X(I+10))**2
    8 F(I+15)=DSQRT(.1D+2)*(X(I)-X(I+15))**2
      IF (MODE.EQ.2) GOTO 2
      DO 11 I=1,5
      DO 11 J=1,20
      DF(I,J)=.0D+0
      IF (J.EQ.I) DF(I,J)=.1D+1
      IF (J.EQ.(I+5)) DF(I,J)=.1D+2
      DF(I+5,J)=.0D+0
      IF (J.EQ.(I+10)) DF(I+5,J)=DSQRT(.5D+1)
      IF (J.EQ.(I+15)) DF(I+5,J)=-DSQRT(.5D+1)
      DF(I+10,J)=.0D+0
      IF (J.EQ.(I+5)) DF(I+10,J)=.2D+1*(X(I+5)-.2D+1*X(I+10))
      IF (J.EQ.(I+10)) DF(I+10,J)=-.4D+1*(X(I+5)-.2D+1*X(I+10))
      DF(I+15,J)=.0D+0
      IF (J.EQ.I) DF(I+15,J)=DSQRT(.1D+2)*.2D+1*(X(I)-X(I+15))
   11 IF (J.EQ.(I+15)) DF(I+15,J)=-DSQRT(.1D+2)*.2D+1*(X(I)-X(I+15))
      DO 10 J=1,20
      GF(J)=.0D+0
      DO 10 I=1,20
   10 GF(J)=GF(J)+.2D+1*F(I)*DF(I,J)
    4 RETURN
      END
C      
      SUBROUTINE TP289(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,2,4,4),MODE
    1 N=30
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,30
      X(I)=(-.1D+1)**I*(.1D+1+DBLE(I)/.3D+2)
      XEX(I)=.0D+0
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=.0D+0
      RETURN
    2 FX=.0D+0
      DO 7 I=1,30
    7 FX=FX+X(I)**2
      FX=.1D+1-DEXP(-FX/.6D+2)
      IF (MODE.EQ.2) RETURN
      DO 8 I=1,30
    8 GF(I)=(FX-.1D+1)*(-.2D+1*X(I)/.6D+2)
    4 RETURN
      END  
C                    
      SUBROUTINE TP290(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION F,DF
      N=2
      GOTO 10
      ENTRY TP291(MODE)
      N=10
      GOTO 10
      ENTRY TP292(MODE)
      N=30
      GOTO 10
      ENTRY TP293(MODE)
      N=50
   10 GOTO(1,3,3,4,4),MODE
    1 NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,N
      X(I)=.1D+1
      XEX(I)=.0D+0
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      FEX=.0D+0
      LEX=.TRUE.
      NEX=1
      RETURN
    2 FX=F(1)**2
      RETURN
    3 F(1)=.0D+0
      DO 7 I=1,N
    7 F(1)=F(1)+DBLE(I)*X(I)**2
      IF (MODE.EQ.2) GOTO 2
      DO 8 I=1,N
    8 DF(1,I)=DBLE(I)*.2D+1*X(I)
      DO 9 J=1,N
    9 GF(J)=.2D+1*F(1)*DF(1,J)
    4 RETURN
      END
C      
      SUBROUTINE TP294(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF
      INTEGER I,J,K
      N=6
      GOTO 20
      ENTRY TP295(MODE)
      N=10
      GOTO 20
      ENTRY TP296(MODE)
      N=16
      GOTO 20
      ENTRY TP297(MODE)
      N=30
      GOTO 20
      ENTRY TP298(MODE)
      N=50
      GOTO 20
      ENTRY TP299(MODE)
      N=100
   20 K=N-1
      GOTO(1,3,3,4,4),MODE
    1 NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,N
      X(I)=-.12D+1
      XEX(I)=.1D+1
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      DO 7 I=1,(N/2)
    7 X(2*I)=.1D+1  
      FEX=.0D+0
      LEX=.TRUE.
      NEX=1
      RETURN
    2 FX=.0D+0
      DO 8 I=1,K
    8 FX=FX+F(I)**2+F(I+K)**2
      FX=FX*1.0D-4
      RETURN
    3 DO 12 I=1,K  
      F(I)=.1D+2*(X(I+1)-X(I)**2)
   12 F(I+K)=.1D+1-X(I)
      IF (MODE.EQ.2) GOTO 2
      DO 9 I=1,K
      DO 9 J=1,N
      DF(I,J)=.0D+0
      IF (J.EQ.I) DF(I,J)=-.2D+2*X(I)
      IF (J.EQ.(I+1)) DF(I,J)=.1D+2
      DF(I+K,J)=.0D+0
    9 IF (J.EQ.I) DF(I+K,J)=-.1D+1
      DO 11 J=1,N
      GF(J)=.0D+0
      DO 13 I=1,2*K
   13 GF(J)=GF(J)+.2D+1*F(I)*DF(I,J)
      GF(J)=GF(J)*1.0D-4
   11 CONTINUE
    4 RETURN
      END
C
      SUBROUTINE TP300(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF
      INTEGER I
      N=20
      FEX=-.2D+2
      GOTO 10
      ENTRY TP301(MODE)
      N=50
      FEX=-.5D+2
      GOTO 10
      ENTRY TP302(MODE)
      N=100
      FEX=-.1D+3
   10 GOTO(1,2,3,4,4),MODE
    1 NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,N
      X(I)=.0D+0
      XEX(I)=DBLE(N-I)+.1D+1
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      RETURN
    2 FX=X(1)**2-.2D+1*X(1)
      DO 7 I=2,N
    7 FX=FX+.2D+1*X(I)**2-.2D+1*X(I-1)*X(I)
      RETURN
    3 GF(1)=.2D+1*X(1)-.2D+1*X(2)-.2D+1
      DO 8 I=2,(N-1)
    8 GF(I)=.4D+1*X(I)-.2D+1*(X(I-1)+X(I+1))
      GF(N)=.4D+1*X(N)-.2D+1*X(N-1)
    4 RETURN                
      END
C
      SUBROUTINE TP303(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION POM
      N=20
      GOTO 10
      ENTRY TP304(MODE)
      N=50
      GOTO 10
      ENTRY TP305(MODE)
      N=100
   10 GOTO(1,2,3,4,4),MODE
    1 NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,N
      X(I)=.1D+0
      XEX(I)=.0D+0
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=.0D+0
      RETURN
    2 POM=.0D+0
      DO 7 I=1,N
    7 POM=POM+.5D+0*DBLE(I)*X(I)
      FX=POM**2+POM**4
      DO 8 I=1,N
    8 FX=FX+X(I)**2
      RETURN
    3 DO 9 I=1,N
    9 GF(I)=.2D+1*X(I)+POM*DBLE(I)+.2D+1*DBLE(I)*POM**3
    4 RETURN
      END
C
      SUBROUTINE TP306(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION A,B
      GOTO(1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,2
      X(I)=1.D+0
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=-.11036D+1
      XEX(1)=.0D+0
      XEX(2)=.1D+1
      RETURN
    2 FX=-DEXP(-X(1)-X(2))*(.2D+1*X(1)**2+.3D+1*X(2)**2)
      RETURN
    3 A=DEXP(-X(1)-X(2))
      B=.2D+1*X(1)**2+.3D+1*X(2)**2
      GF(1)=A*(B-.4D+1*X(1))
      GF(2)=A*(B-.6D+1*X(2))
    4 RETURN
      END
C
      SUBROUTINE TP307(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF,WI,XI1,XI2
      INTEGER I
      GOTO (1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      X(1)=0.3D0
      X(2)=0.4D0
      DO 6 I=1,2
      LXL(I)=.TRUE.
      XL(I)=0.0
      LXU(I)=.TRUE.
C      XU(I)=0.26
      XU(I)=1.0D+10
    6 XEX(I)=.25783D+0
      NEX=0
      LEX=.FALSE.
      FEX=0.12436D+3
      RETURN
    2 FX=0.0D0
      GOTO 9
    3 GF(1)=0.0D0
      GF(2)=0.0D0
    9 DO 7 I=1,10
      WI=DBLE(I)
      XI1=WI*X(1)
      IF (XI1.GT.20) XI1=0.0D0
      XI2=WI*X(2)
      IF (XI2.GT.20) XI2=0.0D0
      F(I)=2.0D0 + 2.0D0*WI - DEXP(XI1) - DEXP(XI2)
      IF (MODE.EQ.2) GOTO 8
      DF(I,1)=-WI*DEXP(XI1)
      DF(I,2)=-WI*DEXP(XI2)
      GF(1)=GF(1) + 2.0D0*F(I)*DF(I,1)
      GF(2)=GF(2) + 2.0D0*F(I)*DF(I,2)
      GOTO 7
    8 FX=FX+F(I)**2
    7 CONTINUE
    4 RETURN
      END
C
      SUBROUTINE TP308(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION F,DF
      GOTO (1,2,2,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      X(1)=.3D+1
      X(2)=.1D+0
      LXU(1)=.FALSE.
      LXU(2)=.FALSE.
      LXL(1)=.FALSE.
      LXL(2)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=.77319906D+0
      XEX(1)=-.15543724D+0
      XEX(2)=.69456378D+0
      RETURN
    2 F(1)=X(1)**2+X(2)**2+X(1)*X(2)
      F(2)=DSIN(X(1))
      F(3)=DCOS(X(2))
      IF (MODE.EQ.3) GOTO 3
      FX=0.D+0
      DO 5 I=1,3
    5 FX=FX+F(I)**2
      RETURN
    3 DF(1,1)=2.D+0*X(1)+X(2)
      DF(1,2)=2.D+0*X(2)+X(1)
      DF(2,1)=2.D+0*DSIN(X(1))*DCOS(X(1))
      DF(2,2)=0.D+0
      DF(3,1)=0.D+0
      DF(3,2)=-2.D+0*DCOS(X(2))*DSIN(X(2))
      GF(1)=0.D+0
      GF(2)=0.D+0
      DO 6 I=1,3
      GF(1)=GF(1)+2.D+0*F(I)*DF(I,1)
    6 GF(2)=GF(2)+2.D+0*F(I)*DF(I,2)
    4 RETURN
      END
C      
      SUBROUTINE TP309(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 5 I=1,2
      X(I)=.0D+0
      LXL(I)=.FALSE.
    5 LXU(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      XEX(1)=.34826826D+1
      XEX(2)=.39D+1
      FEX=-.39871708D+1
    4 RETURN
    2 FX=.141D+1*X(1)**4-.1276D+2*X(1)**3+.3991D+2*X(1)**2
     F   -.5193D+2*X(1)+.2437D+2+(X(2)-.39D+1)**2
      RETURN
    3 GF(1)=.564D+1*X(1)**3-.3828D+2*X(1)**2+.7982D+2*X(1)-.5193D+2
      GF(2)=.2D+1*X(2)-.78D+1
      RETURN
      END
C      
      SUBROUTINE TP310(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION A,B,C
      GOTO (1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      X(1)=-.12D+1
      X(2)=.1D+1
      DO 7 I=1,2
      LXU(I)=.FALSE.
    7 LXL(I)=.FALSE.
      LEX=.TRUE.
      NEX=1
      FEX=0.0D0
      XEX(1)=0.0D0
      XEX(2)=0.0D0
      RETURN
    2 FX=(X(1)*X(2))**2*(.1D+1-X(1))**2
     /    *(.1D+1-X(1)-X(2)*(.1D+1-X(1))**5)**2
      RETURN
    3 A=X(1)*X(2)
      B=.1D+1-X(1)
      C=B-X(2)*(B**5)
      GF(1)=2.D+0*A*B*C*(X(2)-1.D+0-5.D+0*X(2)*(B**4))
      GF(2)=2.D+0*A*B*C*(X(1)-(B**5))
    4 RETURN
      END
C      
      SUBROUTINE TP311(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX)
      INTEGER I
      GOTO (1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 5 I=1,2
      X(I)=.1D+1
      LXL(I)=.FALSE.
    5 LXU(I)=.FALSE.
      LEX=.TRUE.
      NEX=2
      XEX(1)=.3D+1
      XEX(2)=.2D+1
      XEX(3)=3.58443D+0
      XEX(4)=-1.84813D+0
      FEX=.0D+0
    4 RETURN
    2 FX=(X(1)**2+X(2)-.11D+2)**2+(X(1)+X(2)**2-.7D+1)**2
      RETURN
    3 GF(1)=.4D+1*X(1)*(X(1)**2+X(2)-.11D+2)
     F      +.2D+1*(X(1)+X(2)**2-.7D+1)
      GF(2)=.2D+1*(X(1)**2+X(2)-.11D+2)
     F      +.4D+1*X(2)*(X(1)+X(2)**2-.7D+1)
      RETURN
      END    
C      
      SUBROUTINE TP312(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION A,B
      INTEGER I
      GOTO (1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,2
      X(I)=.1D+1
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=0.D+0
      XEX(1)=-.21026652D+2
      XEX(2)=-.36760009D+2
      RETURN
    2 A=X(1)**2+.12D+2*X(2)-.1D+1
      B=.49D+2*(X(1)**2+X(2)**2)+.84D+2*X(1)+.2324D+4*X(2)-.681D+3
      FX=A**2+B**2
      RETURN 
    3 A=X(1)**2+.12D+2*X(2)-.1D+1
      B=.49D+2*(X(1)**2+X(2)**2)+.84D+2*X(1)+.2324D+4*X(2)-.681D+3
      GF(1)=.2D+1*(.2D+1*X(1)*A+B*(.98D+2*X(1)+.84D+2))
      GF(2)=.2D+1*(.12D+2*A+B*(.98D+2*X(2)+.2324D+4))
    4 RETURN
      END
C      
      SUBROUTINE TP313(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION XH
      INTEGER I
      GOTO (1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      X(1)=.0D+0
      X(2)=-.1D+1
      DO 5 I=1,2
      LXL(I)=.FALSE.
    5 LXU(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=.199786D+0
      XEX(1)=.3D+1
      XEX(2)=.2850214D+1
      RETURN
    2 FX=.1D-3*(X(1)-.3D+1)**2-(X(2)-X(1))+DEXP(.2D+2*(X(2)-X(1)))
      RETURN
    3 XH=.2D+2*DEXP(.2D+2*(X(2)-X(1)))
      GF(1)=.1D+1+.2D-3*(X(1)-.3D+1)-XH
      GF(2)=XH-.1D+1
    4 RETURN
      END  
C            
      SUBROUTINE TP314(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION G1,H1,A,B
      INTEGER I
      GOTO (1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,2
      X(I)=.2D+1
      LXU(I)=.FALSE.
      XL(I)=1.0D0      
    6 LXL(I)=.TRUE.
      LEX=.FALSE.
      NEX=2
C      FEX=-0.79712534D+00
C      XEX(1)=0.10326462D+01
C      XEX(2)=0.86695956D+00 
      FEX=.16904D+0
      XEX(1)=.1789039D+1
      XEX(2)=.13740024D+1
      RETURN
    2 A=X(1)-.2D+1
      B=X(2)-.1D+1
      G1=(X(1)**2)/(-.4D+1)-X(2)**2+.1D+1
      H1=X(1)-.2D+1*X(2)+.1D+1
      FX=A**2+B**2+.4D-1/G1+H1**2/.2D+0
      RETURN
    3 G1=(X(1)**2)/(-.4D+1)-X(2)**2+.1D+1 
      H1=X(1)-.2D+1*X(2)+.1D+1
      GF(1)=.2D+1*(X(1)-.2D+1+X(1)*.1D-1/G1**2+.5D+1*H1)
      GF(2)=.2D+1*(X(2)-.1D+1+X(2)*.4D-1/G1**2-.1D+2*H1)
    4 RETURN
      END
C      
      SUBROUTINE TP315(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=1
      NINL=2
      NELI=0
      NENL=0
      X(1)=-.1D+0
      X(2)=-.9D+0
      LXU(1)=.FALSE.
      LXU(2)=.FALSE.
      LXL(1)=.FALSE.
      LXL(2)=.FALSE.
      GG(1,1)=.1D+1
      GG(1,2)=-.2D+1
      GF(1)=.0D+0
      GF(2)=-.1D+1
      LEX=.TRUE.
      NEX=1
      FEX=-.8D+0
      XEX(1)=.6D+0
      XEX(2)=.8D+0
      RETURN
    2 FX=-X(2)
    3 RETURN
    4 IF (INDEX1(1)) G(1)=.1D+1-.2D+1*X(2)+X(1)
      IF (INDEX1(2)) G(2)=X(1)**2+X(2)**2
      IF (INDEX1(3)) G(3)=.1D+1-X(1)**2-X(2)**2
      RETURN
    5 IF (.NOT. INDEX2(2)) GOTO 6
      GG(2,1)=.2D+1*X(1)
      GG(2,2)=.2D+1*X(2)
    6 IF (.NOT. INDEX2(3)) GOTO 7
      GG(3,1)=-.2D+1*X(1)
      GG(3,2)=-.2D+1*X(2)
    7 RETURN
      END
C      
      SUBROUTINE TP316(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=1
      DO 6 I=1,2
      X(I)=.0D+0
      LXL(I)=.FALSE.
    6 LXU(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=0.33431458D+3
      XEX(1)=0.70710678D+1
      XEX(2)=-0.70710678D+1
      RETURN
    2 FX=(X(1)-.2D+2)**2+(X(2)+.2D+2)**2
      RETURN
    3 GF(1)=.2D+1*X(1)-.4D+2
      GF(2)=.2D+1*X(2)+.4D+2
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)**2*.1D-1+X(2)**2*.1D-1-.1D+1
      RETURN
    5 IF (.NOT.INDEX2(1)) RETURN
      GG(1,1)=.2D-1*X(1)
      GG(1,2)=.2D-1*X(2)
      RETURN
      END
C      
      SUBROUTINE TP317(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=1
      DO 6 I=1,2
      X(I)=.0D+0
      LXL(I)=.FALSE.
    6 LXU(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=.37246661D+3
      XEX(1)=.73519262D+1
      XEX(2)=-.5422866D+1
      RETURN
    2 FX=(X(1)-.2D+2)**2+(X(2)+.2D+2)**2
      RETURN
    3 GF(1)=.2D+1*X(1)-.4D+2
      GF(2)=.2D+1*X(2)+.4D+2
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)**2*.1D-1+X(2)**2/.64D+2-.1D+1
      RETURN
    5 IF (.NOT.INDEX2(1)) RETURN
      GG(1,1)=.2D-1*X(1)
      GG(1,2)=.2D+1*X(2)/.64D+2
      RETURN
      END
C      
      SUBROUTINE TP318(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=1
      DO 6 I=1,2
      X(I)=.0D+0
      LXL(I)=.FALSE.
    6 LXU(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=.41275005D+3
      XEX(1)=.78091266D+1
      XEX(2)=-.37478414D+1
      RETURN
    2 FX=(X(1)-.2D+2)**2+(X(2)+.2D+2)**2
      RETURN
    3 GF(1)=.2D+1*X(1)-.4D+2
      GF(2)=.2D+1*X(2)+.4D+2
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)**2*.1D-1+X(2)**2/.36D+2-.1D+1
      RETURN
    5 IF (.NOT.INDEX2(1)) RETURN
      GG(1,1)=.2D-1*X(1)
      GG(1,2)=.2D+1*X(2)/.36D+2
      RETURN
      END
C      
      SUBROUTINE TP319(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=1
      DO 6 I=1,2
      X(I)=.0D+0
      LXL(I)=.FALSE.
    6 LXU(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=.4524044D+3
      XEX(1)=.84922857D+1
      XEX(2)=-.21121017D+1
      RETURN
    2 FX=(X(1)-.2D+2)**2+(X(2)+.2D+2)**2
      RETURN
    3 GF(1)=.2D+1*X(1)-.4D+2
      GF(2)=.2D+1*X(2)+.4D+2
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)**2*.1D-1+X(2)**2/.16D+2-.1D+1
      RETURN
    5 IF (.NOT.INDEX2(1)) RETURN
      GG(1,1)=.2D-1*X(1)
      GG(1,2)=.2D+1*X(2)/.16D+2
      RETURN
      END
C      
      SUBROUTINE TP320(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=1
      DO 6 I=1,2
      X(I)=.0D+0
      LXL(I)=.FALSE.
    6 LXU(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=.48553146D+3
      XEX(1)=.939525D+1
      XEX(2)=-.68459019D+0
      RETURN
    2 FX=((X(1)-.2D+2)**2+(X(2)+.2D+2)**2)
      RETURN
    3 GF(1)=(0.2D+1*X(1)-0.4D+2)
      GF(2)=(0.2D+1*X(2)+0.4D+2)
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)**2*.1D-1+X(2)**2/.4D+1-.1D+1
      RETURN
    5 IF (.NOT.INDEX2(1)) RETURN
      GG(1,1)=.2D-1*X(1)
      GG(1,2)=.5D+0*X(2)
      RETURN
      END
C      
      SUBROUTINE TP321(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=1
      DO 6 I=1,2
      X(I)=.0D+0
      LXL(I)=.FALSE.
    6 LXU(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=.49611237D+3
      XEX(1)=.98160292D+1
      XEX(2)=-.19093377D+0
      RETURN
    2 FX=(X(1)-.2D+2)**2+(X(2)+.2D+2)**2
      RETURN
    3 GF(1)=.2D+1*X(1)-.4D+2
      GF(2)=.2D+1*X(2)+.4D+2
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)**2*.1D-1+X(2)**2-.1D+1
      RETURN
    5 IF (.NOT.INDEX2(1)) RETURN
      GG(1,1)=.2D-1*X(1)
      GG(1,2)=.2D+1*X(2)
      RETURN
      END
C
      SUBROUTINE TP322(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=1
      DO 6 I=1,2
      X(I)=0.0001
      LXL(I)=.FALSE.
    6 LXU(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=.49996001D+3
      XEX(1)=.99980018D+1
      XEX(2)=-.19990011D-2
      RETURN
    2 FX=(X(1)-20.0)**2+(X(2)+20.0)**2
      RETURN
    3 GF(1)=.2D+1*X(1)-.4D+2
      GF(2)=.2D+1*X(2)+.4D+2
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)**2*.1D-1+X(2)**2*.1D+3-.1D+1
      RETURN
    5 IF (.NOT.INDEX2(1)) RETURN
      GG(1,1)=.2D-1*X(1)
      GG(1,2)=.2D+3*X(2)
      RETURN
      END
C      
      SUBROUTINE TP323(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=1
      NINL=1
      NELI=0
      NENL=0
      X(1)=.0D+0
      X(2)=.1D+1
      LXL(1)=.TRUE.
      XL(1)=.0D+0
      LXL(2)=.TRUE.
      XL(2)=.0D+0
      LXU(1)=.FALSE.
      LXU(2)=.FALSE.
      GG(1,1)=.1D+1
      GG(1,2)=-.1D+1
      GG(2,2)=.1D+1
      LEX=.FALSE.
      NEX=1
      FEX=.37989446D+1
      XEX(1)=.55357378D+0
      XEX(2)=.13064439D+1
      RETURN
    2 FX=X(1)**2+X(2)**2-.4D+1*X(1)+.4D+1
      RETURN
    3 GF(1)=.2D+1*X(1)-.4D+1
      GF(2)=.2D+1*X(2)
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)-X(2)+.2D+1
      IF (INDEX1(2)) G(2)=-.1D+1*X(1)**2+X(2)-.1D+1
      RETURN
    5 IF (INDEX2(2)) GG(2,1)=-.2D+1*X(1)
      RETURN
      END
C      
      SUBROUTINE TP324(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=2
      NELI=0
      NENL=0
      DO 6 I=1,2
      X(I)=.2D+1
    6 LXU(I)=.FALSE.
      LXL(1)=.TRUE.
      LXL(2)=.FALSE.
      XL(1)=.2D+1
      LEX=.FALSE.
      NEX=1
      FEX=.5D+1
      XEX(1)=.15811389D+2
      XEX(2)=.15811387D+1
      RETURN
    2 FX=.1D-1*X(1)**2+X(2)**2
      RETURN
    3 GF(1)=.2D-1*X(1)
      GF(2)=.2D+1*X(2)
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)*X(2)-.25D+2
      IF (INDEX1(2)) G(2)=X(1)**2+X(2)**2-.25D+2
      RETURN
    5 IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=X(2)
      GG(1,2)=X(1)
    7 IF (.NOT.INDEX2(2)) RETURN
      GG(2,1)=.2D+1*X(1)
      GG(2,2)=.2D+1*X(2)
      RETURN
      END
C      
      SUBROUTINE TP325(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=2
      NILI=1
      NINL=1
      NELI=0
      NENL=1
      X(1)=-.3D+1
      X(2)=.0D+0
      DO 6 I=1,2
      GG(1,I)=-.1D+1
      LXL(I)=.FALSE.
    6 LXU(I)=.FALSE.
      GG(2,1)=-.1D+1
      GF(2)=.1D+1
      LEX=.FALSE.
      NEX=1
      FEX=0.37913414D+1 
      XEX(1)=-0.23722813D+1
      XEX(2)=-0.18363772D+1
      RETURN
    2 FX=X(1)**2+X(2)
      RETURN
    3 GF(1)=0.2D+1*X(1)
      RETURN
    4 IF (INDEX1(1)) G(1)=-(X(1)+X(2))+0.1D+1
      IF (INDEX1(2)) G(2)=-(X(1)+X(2)**2)+0.1D+1
      IF (INDEX1(3)) G(3)=X(1)**2+X(2)**2-0.9D+1
      RETURN
    5 IF (INDEX2(2)) GG(2,2)=-.2D+1*X(2)
      IF (.NOT.INDEX2(3)) RETURN
      GG(3,1)=.2D+1*X(1)
      GG(3,2)=.2D+1*X(2)
      RETURN
      END
C      
      SUBROUTINE TP326(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=2
      NELI=0
      NENL=0
      X(1)=.4D+1
      X(2)=.3D+1
      DO 6 I=1,2
      LXL(I)=.TRUE.
      LXU(I)=.TRUE.
      XU(I)=10.0
    6 XL(I)=.0D+0
      GG(1,2)=-.4D+1
      LEX=.FALSE.
      NEX=1
      FEX=-.79807821D+2
      XEX(1)=.52396091D+1
      XEX(2)=.37460378D+1
      RETURN
    2 FX=X(1)**2+X(2)**2-.16D+2*X(1)-.1D+2*X(2)
      RETURN
    3 GF(1)=.2D+1*X(1)-.16D+2
      GF(2)=.2D+1*X(2)-.1D+2
      RETURN
    4 IF (INDEX1(1)) G(1)=.11D+2-X(1)**2+.6D+1*X(1)-.4D+1*X(2)
      IF (INDEX1(2)) G(2)=X(1)*X(2)-.3D+1*X(2)-DEXP(X(1)-.3D+1)+.1D+1
      RETURN
    5 IF (INDEX2(1)) GG(1,1)=-.2D+1*X(1)+.6D+1
      IF (.NOT.INDEX2(2)) RETURN
      GG(2,1)=X(2)-DEXP(X(1)-.3D+1)
      GG(2,2)=X(1)-.3D+1
      RETURN
      END
C      
      SUBROUTINE TP327(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION Y(44),Z(44),F,DF
      INTEGER I
      DATA Y/2*.49D+0,.48D+0,.47D+0,.48D+0,.47D+0,2*.46D+0,.45D+0,
     F       .43D+0,.45D+0,2*.43D+0,.44D+0,2*.43D+0,.46D+0,.45D+0,
     F       2*.42D+0,.43D+0,2*.41D+0,.4D+0,.42D+0,2*.4D+0,.41D+0,
     F       .4D+0,2*.41D+0,3*.4D+0,.38D+0,.41D+0,2*.4D+0,.41D+0,
     F       .38D+0,2*.4D+0,2*.39D+0/
      DATA Z/2*.8D+1,4*.1D+2,4*.12D+2,3*.14D+2,3*.16D+2,2*.18D+2,
     F       3*.2D+2,3*.22D+2,3*.24D+2,3*.26D+2,2*.28D+2,3*.3D+2,
     F       2*.32D+2,.34D+2,2*.36D+2,2*.38D+2,.4D+2,.42D+2/
      GOTO(1,2,2,4,5),MODE
    1 N=2
      NILI=0
      NINL=1
      NELI=0
      NENL=0
      X(1)=.42D+0
      X(2)=.5D+1
      DO 12 I=1,2
      LXU(I)=.FALSE.
      LXL(I)=.TRUE.
   12 XL(I)=.4D+0
      LEX=.FALSE.
      NEX=2
      FEX=0.28459670D-01
      XEX(1)=0.41995259D+0
      XEX(2)=0.12848442D+01
C      FEX=.30646306D-1
c      XEX(3)=.42190424D+0
c      XEX(4)=.50000526D+1
      RETURN
    2 DO 6 I=1,44
    6 F(I)=Y(I)-X(1)-(.49D+0-X(1))*DEXP(-X(2)*(Z(I)-.8D+1))
      IF (MODE.EQ.3) GOTO 3
      FX=.0D+0
      DO 7 I=1,44
    7 FX=FX+F(I)**2
      RETURN
    3 GF(1)=.0D+0
      GF(2)=.0D+0
      DO 8 I=1,44
      DF(I,1)=-.1D+1+DEXP(-X(2)*(Z(I)-.8D+1))
      DF(I,2)=(.49D+0-X(1))*DEXP(-X(2)*(Z(I)-.8D+1))*(Z(I)-.8D+1)
      GF(1)=GF(1)+DF(I,1)*F(I)*.2D+1
    8 GF(2)=GF(2)+DF(I,2)*F(I)*.2D+1
      RETURN
    4 IF (INDEX1(1)) G(1)=-.09D+0-X(1)*X(2)+.49D+0*X(2)
      RETURN
    5 IF (.NOT.INDEX2(1)) RETURN
      GG(1,1)=-X(2)
      GG(1,2)=.49D+0-X(1)
      RETURN
      END
C      
      SUBROUTINE TP328(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION A,B
      GOTO(1,2,3,4,4),MODE
    1 N=2
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,2
      X(I)=.5D+0
      LXU(I)=.TRUE.
      XU(I)=.3D+1
      LXL(I)=.TRUE.
    6 XL(I)=.1D+0
      LEX=.FALSE.
      NEX=1
      FEX=.1744152D+1
      XEX(1)=.1743439D+1
      XEX(2)=.20297056D+1
      RETURN
    2 A=(.1D+1+X(2)**2)/X(1)**2
      B=((X(1)*X(2))**2+.1D+3)/(X(1)*X(2))**4
      FX=(.12D+2+X(1)**2+A+B)/.1D+2
      RETURN
    3 A=(.1D+1+X(2)**2)/X(1)**3
      B=.1D+1/(X(1)**3*X(2)**2)+.2D+3/(X(1)**5*X(2)**4)
      GF(1)=(X(1)-A-B)/.5D+1
      A=.1D+1/(X(1)**2*X(2)**3)+.2D+3/(X(1)**4*X(2)**5)
      GF(2)=(X(2)/X(1)**2-A)/.5D+1
    4 RETURN
      END
C      
      SUBROUTINE TP329(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      GOTO (1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=3
      NELI=0
      NENL=0
      LXL(1)=.TRUE.
      LXL(2)=.TRUE.
      LXU(1)=.TRUE.
      LXU(2)=.TRUE.
      X(1)=.1435D+2
      X(2)=.86D+1
      XL(1)=.13D+2
      XL(2)=.0D+0
      XU(1)=.16D+2
      XU(2)=.15D+2
      LEX=.FALSE.
      NEX=1
      FEX=-.69618139D+4
      XEX(1)=.14095D+2
      XEX(2)=.84296079D+0
      RETURN
    2 FX=(X(1)-.1D+2)**3+(X(2)-.2D+2)**3
      RETURN
    3 GF(1)=.3D+1*(X(1)-.1D+2)**2
      GF(2)=.3D+1*(X(2)-.2D+2)**2
      RETURN
    4 IF (INDEX1(1)) G(1)=(X(1)-.5D+1)**2+(X(2)-.5D+1)**2-.1D+3
      IF (INDEX1(2)) G(2)=(X(1)-.6D+1)**2+(X(2)-.5D+1)**2
      IF (INDEX1(3)) G(3)=.8281D+2-(X(1)-.6D+1)**2-(X(2)-.5D+1)**2
      RETURN
    5 IF (.NOT.INDEX2(1)) GOTO 6
      GG(1,1)=.2D+1*X(1)-.1D+2
      GG(1,2)=.2D+1*X(2)-.1D+2
    6 IF (.NOT.INDEX2(2)) GOTO 7
      GG(2,1)=.2D+1*X(1)-.12D+2
      GG(2,2)=.2D+1*X(2)-.1D+2
    7 IF (.NOT.INDEX2(3)) RETURN
      GG(3,1)=-.2D+1*X(1)+.12D+2
      GG(3,2)=-.2D+1*X(2)+.1D+2
      RETURN
      END
C      
      SUBROUTINE TP330(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=2
      NILI=0
      NINL=1
      NELI=0
      NENL=0
      DO 6 I=1,2
      X(I)=0.25D+1
      LXL(I)=.TRUE.
      LXU(I)=.TRUE.
      XL(I)=0.1D-4
    6 XU(I)=0.5D+1
      LEX=.FALSE.
      NEX=1
      FEX=.16205833D+1
      XEX(1)=0.12866773D+1
      XEX(2)=0.53046181D+0
      RETURN
    2 CONTINUE
      IF (X(1).LT.XL(1)) X(1)=XL(1)
      IF (X(2).LT.XL(2)) X(2)=XL(2)
      FX=0.44D-1*X(1)**3/X(2)**2+1.D+0/X(1)+0.592D-1*X(1)/X(2)**3
      RETURN
    3 GF(1)=.132D+0*X(1)**2/X(2)**2-X(1)**(-2)+.592D-1/X(2)**3
      GF(2)=-.88D-1*X(1)**3/X(2)**3-.1776D+0*X(1)/X(2)**4
      RETURN
    4 CONTINUE
      IF (X(1).LT.XL(1)) X(1)=XL(1)
      IF (INDEX1(1)) G(1)=.1D+1-.862D+1*X(2)**3/X(1)
      RETURN
    5 CONTINUE
      IF (X(1).LT.XL(1)) X(1)=XL(1)
      IF (.NOT.INDEX2(1)) RETURN
      GG(1,1)=.862D+1*X(2)**3/X(1)**2
      GG(1,2)=-.2586D+2*X(2)**2/X(1)
      RETURN
      END
C
      SUBROUTINE TP331(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF,A,B,C,X1
      GOTO(1,2,3,4,5),MODE
    1 N=2
      NILI=1
      NINL=0
      NELI=0
      NENL=0
      X(1)=0.5D0
      X(2)=0.1D0
      LXU(1)=.TRUE.
      LXU(2)=.TRUE.
      XU(1)=0.7D0
      XU(2)=0.2D0
      LXL(1)=.TRUE.
      XL(1)=0.3D0
      LXL(2)=.TRUE.
      XL(2)=0.1D0
      GG(1,1)=-0.1D+1
      GG(1,2)=-0.1D+1
      LEX=.FALSE.
      NEX=1
      FEX=0.4258D+1
      XEX(1)=0.6175D+0
      XEX(2)=0.1039D+0
      RETURN
    2 CONTINUE
      IF (X(1).LT.XL(1)) X(1)=XL(1)
      IF (X(2).LT.XL(2)) X(2)=XL(2)
      FX=(DLOG(DABS(2.0D0*DLOG(X(2))/DLOG(X(1)+X(2)))))/X(1)
      RETURN
    3 CONTINUE
      IF (X(1).LT.XL(1)) X(1)=XL(1)
      IF (X(2).LT.XL(2)) X(2)=XL(2)
      A=X(1)+X(2)
      B=DLOG(A)
      C=2.0D0*DLOG(X(2))
      GF(1)=(-0.1D+1/X(1))*(DLOG(DMAX1(C/B,1.0D-4))/X(1)+0.1D+1/(B*A))
      GF(2)=((0.2D+1*B)/X(2)-C/A)/(C*B*X(1))
      RETURN
    4 IF (INDEX1(1)) G(1)=1.0D0-X(1)-X(2)
    5 RETURN
      END
C
      SUBROUTINE TP332(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION A,B,C,PI,TR,XXX,YYY,PIM,PBIG,PANGLE
      INTEGER I      
      DATA PI/.31415926535D+1/
      GOTO(1,2,3,4,3),MODE
    1 N=2
      NILI=0
      NINL=2
      NELI=0
      NENL=0
      DO 6 I=1,2
      X(I)=.75D+0
      LXU(I)=.TRUE.
      XU(I)=.15D+1
      LXL(I)=.TRUE.
    6 XL(I)=.0D+0
      LEX=.FALSE.
      NEX=1
      FEX=.11495015D+3
      XEX(1)=.91139872D+0
      XEX(2)=.29280207D-1
    3 RETURN
    2 PIM=PI/.36D+1
      FX=.0D+0
      DO 7 I=1,100
      TR=PI*((.1D+1/.3D+1)+((DBLE(I)-.1D+1)/.18D+3))
      A=DLOG(TR)
      B=DSIN(TR)
      C=DCOS(TR)
      XXX=((A+X(2))*B+X(1)*C)
      YYY=((A+X(2))*C-X(1)*B)
    7 FX=FX+PIM*(XXX**2+YYY**2)
      RETURN
    4 PBIG=-.36D+3
      PIM=.18D+3/PI
      DO 8 I=1,100
      TR=PI*((.1D+1/.3D+1)+((DBLE(I)-.1D+1)/.18D+3))
      A=.1D+1/TR-X(1)
      B=DLOG(TR)+X(2)
      PANGLE=PIM*DATAN(DABS(A/B))
    8 IF (PANGLE.GT.PBIG) PBIG=PANGLE
      IF (INDEX1(1)) G(1)=.3D+2-PBIG
      IF (INDEX1(2)) G(2)=PBIG+.3D+2
      RETURN
      END
C
      SUBROUTINE TP333(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF,A(8),Y(8)
      DATA A/.4D+1,.575D+1,.75D+1,.24D+2,.32D+2,.48D+2,.72D+2,.96D+2/
      DATA Y/.721D+2,.656D+2,.559D+2,.171D+2,.98D+1,.45D+1,.13D+1,.6D+0/
      INTEGER I,LSUM
      GOTO(1,2,2,4,4),MODE
    1 N=3
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      X(1)=.3D+2
      X(2)=.04D+0
      X(3)=.3D+1
      DO 6 I=1,3
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LXL(2)=.TRUE.
      XL(2)=0.0
      LXU(2)=.TRUE.
      XU(2)=0.07D0
      LEX=.FALSE.
      NEX=1
      FEX=0.0432
      XEX(1)=.89902D+2
      XEX(2)=.06699D+0
      XEX(3)=.47809D+0
      LSUM=8
      RETURN
    2 DO 7 I=1,8
    7 F(I)=(Y(I)-X(1)*DEXP(-X(2)*A(I))-X(3))/Y(I)
      IF (MODE.EQ.3) GOTO 3
      FX=.0D+0
      DO 8 I=1,8
    8 FX=FX+F(I)**2
      RETURN
    3 GF(1)=.0D+0
      GF(2)=.0D+0
      GF(3)=.0D+0
      DO 9 I=1,8
      DF(I,1)=(-DEXP(-X(2)*A(I)))/Y(I)
      DF(I,2)=(X(1)*A(I)*DEXP(-X(2)*A(I)))/Y(I)
      DF(I,3)=-1.D+0/Y(I)
      GF(1)=GF(1)+DF(I,1)*F(I)*.2D+1
      GF(2)=GF(2)+DF(I,2)*F(I)*.2D+1
    9 GF(3)=GF(3)+DF(I,3)*F(I)*.2D+1
    4 RETURN
      END
C
      SUBROUTINE TP334(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,LSUM
      DOUBLEPRECISION F,DF,Y(15),UI,VI,WI
      DATA(Y(I),I=1,15)/0.14D+0,0.18D+0,0.22D+0,0.25D+0,0.29D+0,
     F          0.32D+0,0.35D+0,0.39D+0,0.37D+0,0.58D+0,0.73D+0,
     F          0.96D+0,0.134D+1,0.21D+1,0.439D+1/
      GOTO(1,2,2,4,4),MODE
 1    N=3
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      LSUM=15
      DO 6  I=1,3
      X(I)=0.1D+1
      LXL(I)=.FALSE.
 6    LXU(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      XEX(1)=0.82410544D-1
      XEX(2)=0.11330358D+1
      XEX(3)=0.23436955D+1
C      FEX=0.82149184D-2
      FEX=0.82148773D-2 
      RETURN
 2    DO 7  I=1,15
      UI=DBLE(I)
      VI=DBLE(16-I)
      WI=DMIN1(UI,VI)
 7    F(I)=Y(I)-(X(1)+I/(X(2)*VI+X(3)*WI))
      IF (MODE .EQ. 3) GOTO 3
      FX=0.D+0
      DO 10 I=1,15
 10   FX=FX+F(I)**2
      RETURN
 3    DO 8  I=1,3
 8    GF(I)=0.0D+0
      DO 9  I=1,15
      UI=DBLE(I)
      VI=DBLE(16-I)
      WI=DMIN1(UI,VI)
      DF(I,1)=-0.1D+1
      DF(I,2)=(UI*VI)/(X(2)*VI+X(3)*WI)**2
      DF(I,3)=(UI*WI)/(X(2)*VI+X(3)*WI)**2
      GF(1)=GF(1)+DF(I,1)*F(I)*.2D+1
      GF(2)=GF(2)+DF(I,2)*F(I)*.2D+1
 9    GF(3)=GF(3)+DF(I,3)*F(I)*.2D+1
 4    RETURN
      END             
C                         
      SUBROUTINE TP335(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
 1    N=3
      NILI=0
      NINL=0
      NELI=0
      NENL=2
      DO 6  I=1,3
      X(I)=0.1D+1
      LXL(I)=.FALSE.
 6    LXU(I)=.FALSE.
      GG(1,3)=-0.1D+1
      GG(2,3)= 0.1D+1
      LEX=.FALSE.
      NEX=1
      XEX(1)=.20309475D-05
      XEX(2)=.44721349D-02
      XEX(3)=.20000032D-02
      FEX=-.44721370D-02
      RETURN
 2    FX=-(0.1D-2*X(1)+X(2))
      RETURN
 3    GF(1)=-0.1D-2
      GF(2)=-0.1D+1
      GF(3)=0.0D+0
      RETURN
 4    IF (INDEX1(1))  G(1)=0.1D+4*X(1)**2+0.1D+3*X(2)**2-X(3)
      IF (INDEX1(2))  G(2)=0.1D+3*X(1)**2+0.4D+3*X(2)**2+X(3)-0.1D-1  
      RETURN
 5    IF (.NOT. INDEX2(1))  GOTO 7
      GG(1,1)=0.2D+4*X(1)
      GG(1,2)=0.2D+3*X(2)
      IF (.NOT. INDEX2(2))  GOTO 7
      GG(2,1)=0.2D+3*X(1)
      GG(2,2)=0.8D+3*X(2)
 7    RETURN
      END
C      
      SUBROUTINE TP336 (MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
 1    N=3
      NILI=0
      NINL=0
      NELI=1
      NENL=1
      DO 6 I=1,3
      X(I)=.0D+0
      LXU(I)=.FALSE.
 6    LXL(I)=.FALSE.
      GG(1,1)=.5D+1
      GG(1,2)=.5D+1
      GG(1,3)=-.3D+1
      LEX=.FALSE.
      NEX=1
      XEX(1)=.53459441D+00
      XEX(2)=.53397092D+00
      XEX(3) =-.21905778D+00
      FEX=-.33789573D+00
      RETURN
 2    FX=.7D+1*X(1)-.6D+1*X(2)+.4D+1*X(3)
      RETURN
 3    GF(1)=.7D+1
      GF(2)=-.6D+1
      GF(3)=.4D+1
      RETURN
 4    IF (INDEX1(1)) G(1)=.5D+1*X(1)+.5D+1*X(2)-.3D+1*X(3)-.6D+1
      IF (INDEX1(2)) G(2)=X(1)**2+.2D+1*X(2)**2+.3D+1*X(3)**2-.1D+1
      RETURN
 5    IF (.NOT.INDEX2(2)) GOTO 8
      GG(2,1)=.2D+1*X(1)
      GG(2,2)=.4D+1*X(2)
      GG(2,3)=.6D+1*X(3)
 8    RETURN
      END
C      
      SUBROUTINE TP337(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
 1    N=3
      NILI=0
      NINL=1
      NELI=0
      NENL=0
      DO 6 I=1,3
 6    X(I)=.1D+1
      LXU(1)=.FALSE.
      LXL(1)=.FALSE.
      LXU(2)=.FALSE.
      LXL(2)=.TRUE.
      LXU(3)=.TRUE.
      LXL(3)=.FALSE.
      XL(2)=.1D+1
      XU(3)=.1D+1
      GG(1,3)=.0D+0
      LEX=.FALSE.
      NEX=1
      XEX(1)=.57735194D+00
      XEX(2)=.17320458D+01
      XEX(3) =-.20256839D-05
      FEX=.6D+1
      RETURN 
 2    FX=.9D+1*X(1)**2+X(2)**2+.9D+1*X(3)**2
      RETURN
 3    GF(1)=.18D+2*X(1)
      GF(2)=.2D+1*X(2)
      GF(3)=.18D+2*X(3)
      RETURN
 4    IF (INDEX1(1)) G(1)=X(1)*X(2)-.1D+1
      RETURN
 5    IF (.NOT.INDEX2(1)) GOTO7
      GG(1,1)=X(2)
      GG(1,2)=X(1)
 7    RETURN
      END
C      
      SUBROUTINE TP338(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
 1    N=3
      NILI=0
      NINL=0
      NELI=1
      NENL=1
      DO 6  I=1,3
      X(I)=0.0D+0
      LXL(I)=.FALSE.
 6    LXU(I)=.FALSE.
      GG(1,1)=0.5D+0
      GG(1,2)=0.1D+1
      GG(1,3)=0.1D+1
      LEX=.FALSE.
      NEX=1
c      XEX(1)=.36689438D+00
c      XEX(2)=.22437202D+01
c      XEX(3)=-.14271674D+01
c      FEX=-.72056984D+01
      XEX(1)=-0.36653028D+00
      XEX(2)=-0.16620759D+01
      XEX(3)=0.28453410D+01 
      FEX=-0.10992806D+02
      RETURN
 2    FX=-(X(1)**2+X(2)**2+X(3)**2)
      RETURN
 3    DO 7  I=1,3
 7    GF(I)=-0.2D+1*X(I)
      RETURN
 4    IF (INDEX1(1)) G(1)=0.5D+0*X(1)+X(2)+X(3)-0.1D+1
      IF (INDEX1(2)) G(2)=X(1)**2+(0.2D+1/0.3D+1)*X(2)**2
     F                    +0.25D+0*X(3)**2-0.4D+1
      RETURN
 5    IF (.NOT.INDEX2(2))  GOTO 8
      GG(2,1)=0.2D+1*X(1)
      GG(2,2)=(0.4D+1/0.3D+1)*X(2)
      GG(2,3)=0.5D+0*X(3)
 8    RETURN
      END         
C                                          
      SUBROUTINE TP339 (MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
 1    N=3
      NILI=0
      NINL=1
      NELI=0
      NENL=0
      DO 6 I=1,3	
      X(I)=0.1D+1
      LXU(I)=.FALSE.
      LXL(I)=.TRUE.
 6    XL(I)=0.01D+0
      LEX=.FALSE.
      NEX=1
      XEX(1)=0.23797626D+01
      XEX(2)=0.31622787D+00
      XEX(3)=0.19429359D+01
      FEX=0.33616797D+01
      RETURN
 2    FX=0.2D+0/(X(1)*X(2)*X(3))+0.4D+1/X(1)+0.3D+1/X(3)
      RETURN
 3    GF(1)=-0.2D+0/(X(2)*X(3)*X(1)**2)-0.4D+1/X(1)**2
      GF(2)=-0.2D+0/(X(1)*X(3)*X(2)**2)
      GF(3)=-0.2D+0/(X(1)*X(2)*X(3)**2)-0.3D+1/X(3)**2
      RETURN
 4    IF (INDEX1(1)) G(1)=0.1D+2-0.2D+1*X(1)*X(3)-X(1)*X(2)
      RETURN
 5    IF (.NOT.INDEX2(1)) GOTO 8
      GG(1,1)=-0.2D+1*X(3)-X(2)
      GG(1,2)=-X(1)
      GG(1,3)=-0.2D+1*X(1)
 8    RETURN
      END
C      
      SUBROUTINE TP340(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
 1    N=3
      NILI=1
      NINL=0
      NELI=0
      NENL=0
      DO6 I=1,3
      X(I)=.1D+1
      LXU(I)=.TRUE.
C 6    LXL(I)=.FALSE.
      XL(I)=0.0D0
 6    LXL(I)=.TRUE.
      LXU(1)=.TRUE.
      XU(1)=.1D+1
      GG(1,1)=-.1D+1
      GG(1,2)=-.2D+1
      GG(1,3)=-.2D+1
      LEX=.TRUE.
      NEX=1
      XEX(1)=.6D+00
      XEX(2)=.3D+00
      XEX(3)=.3D+00
      FEX=-.54D-01
      RETURN
 2    FX=-X(1)*X(2)*X(3)
      RETURN
 3    GF(1)=-X(2)*X(3)
      GF(2)=-X(1)*X(3)
      GF(3)=-X(1)*X(2)
      RETURN
 4    IF (INDEX1(1)) G(1)=.18D+1-X(1)-.2D+1*X(2)-.2D+1*X(3)
 5    RETURN
      END
C      
      SUBROUTINE TP341(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
 1    N=3
      NILI=0
      NINL=1
      NELI=0
      NENL=0
      DO 6  I=1,3
      X(I)=0.1D+1
      LXL(I)=.TRUE.
      LXU(I)=.FALSE.
 6    XL(I)=0.0D+0
      LEX=.FALSE.
      NEX=1
      XEX(1)=.40000000D+01
      XEX(2)=.28284271D+01
      XEX(3)=.20000000D+01
      FEX=-.22627417D+02
      RETURN
 2    FX=-X(1)*X(2)*X(3)
      RETURN
 3    GF(1)=-X(2)*X(3)
      GF(2)=-X(1)*X(3)
      GF(3)=-X(1)*X(2)
      RETURN
 4    IF (INDEX1(1))  G(1)=-.1D+1*X(1)**2-0.2D+1*X(2)**2
     F                      -0.4D+1*X(3)**2+0.48D+2
      RETURN
 5    IF (.NOT.INDEX2(1))  GOTO 7
      GG(1,1)=-0.2D+1*X(1)
      GG(1,2)=-0.4D+1*X(2)
      GG(1,3)=-0.8D+1*X(3)
 7    RETURN
      END     
C          
      SUBROUTINE TP342 (MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
 1    N=3
      NILI=0
      NINL=1
      NELI=0
      NENL=0
      DO 6 I=1,3	
      X(I)=.1D+1
      LXU(I)=.FALSE.
      LXL(I)=.TRUE.
 6    XL(I)=.0D+0
      LEX=.FALSE.
      NEX=1
      XEX(1)=.40000000D+01
      XEX(2)=.28284271D+01
      XEX(3)=.20000000D+01
      FEX=-.22627417D+02
      RETURN
 2    FX=-X(1)*X(2)*X(3)
      RETURN
 3    GF(1)=-X(2)*X(3)
      GF(2)=-X(1)*X(3)
      GF(3)=-X(1)*X(2)
      RETURN
 4    IF (INDEX1(1)) G(1)=.48D+2-X(1)**2-.2D+1*X(2)**2-.4D+1*X(3)**2  
      RETURN
 5    IF (.NOT.INDEX2(1)) GOTO 8
      GG(1,1)=-.2D+1*X(1)
      GG(1,2)=-.4D+1*X(2)
      GG(1,3)=-.8D+1*X(3)
 8    RETURN
      END
C      
      SUBROUTINE TP343(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
 1    N=3
      NILI=0
      NINL=2
      NELI=0
      NENL=0
      DO 6 I=1,3
      LXU(I)=.TRUE.
      LXL(I)=.TRUE.
 6    XL(I)=.0D+0
      X(1)=.223D+2
      X(2)=.5D+0
      X(3)=.125D+3
      XU(1)=.36D+2
      XU(2)=.5D+1
      XU(3)=.125D+3
      GG(1,3)=.0D+0
      GG(2,2)=.0D+0
      LEX=.FALSE.
      NEX=1
      XEX(1)=.16508383D+02
      XEX(2)=.24768216D+01
      XEX(3)=.12399452D+03
      FEX=-.56847825D+01
      RETURN
 2    FX=-.201D-1*X(1)**4*X(2)*X(3)**2*.1D-6
      RETURN
 3    GF(1)=-.804D-1*X(1)**3*X(2)*X(3)**2*.1D-6
      GF(2)=-.201D-1*X(1)**4*X(3)**2*.1D-6
      GF(3)=-.402D-1*X(1)**4*X(2)*X(3)*.1D-6
      RETURN
 4    IF (INDEX1(1)) G(1)=.675D+3-X(1)**2*X(2)
      IF (INDEX1(2)) G(2)=.419D+0-X(1)**2*X(3)**2*.1D-6
      RETURN
 5    IF (.NOT.INDEX2(1)) GOTO7
      GG(1,1)=-.2D+1*X(1)*X(2)
      GG(1,2)=-X(1)**2
 7    IF (.NOT.INDEX2(2)) GOTO8
      GG(2,1)=-.2D-6*X(1)*X(3)**2
      GG(2,3)=-.2D-6*X(1)**2*X(3)
 8    RETURN
      END
C      
      SUBROUTINE TP344(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO(1,2,3,4,5),MODE
 1    N=3
      NILI=0
      NINL=0
      NELI=0
      NENL=1
      DO 6  I=1,3
      X(I)=0.2D+1
      LXL(I)=.FALSE.
 6    LXU(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      XEX(1)=.11048590D+01
      XEX(2)=.11966742D+01
      XEX(3)=.15352623D+01
      FEX=.32568200D-01
      RETURN
 2    FX=(X(1)-0.1D+1)**2+(X(1)-X(2))**2+(X(2)-X(3))**4
      RETURN
 3    GF(1)=0.2D+1*(X(1)-0.1D+1)+0.2D+1*(X(1)-X(2))
      GF(2)=-0.2D+1*(X(1)-X(2))+0.4D+1*(X(2)-X(3))**3
      GF(3)=-0.4D+1*(X(2)-X(3))**3
      RETURN
 4    IF (INDEX1(1))  G(1)=X(1)*(0.1D+1+X(2)**2)+X(3)**4
     F                     -0.4D+1-0.3D+1*DSQRT(0.2D+1)
      RETURN
 5    IF (.NOT. INDEX2(1))  GOTO 7
      GG(1,1)=0.1D+1+X(2)**2
      GG(1,2)=0.2D+1*X(1)*X(2)
      GG(1,3)=0.4D+1*X(3)**3
 7    RETURN
      END
C      
      SUBROUTINE TP345(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
 1    N=3
      NILI=0
      NINL=0
      NELI=0
      NENL=1
      DO 6 I=1,3	
      X(I)=.0D+0
      LXU(I)=.FALSE.
 6    LXL(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      XEX(1)=.11048584D+01
      XEX(2)=.11966752D+01
      XEX(3)=.15352622D+01
      FEX=.32568200D-01
      RETURN
 2    FX=(X(1)-.1D+1)**2+(X(1)-X(2))**2+(X(2)-X(3))**4
      RETURN
 3    GF(1)=.2D+1*(X(1)-1)+.2D+1*(X(1)-X(2))
      GF(2)=-.2D+1*(X(1)-X(2))+.4D+1*(X(2)-X(3))**3
      GF(3)=-.4D+1*(X(2)-X(3))**3
      RETURN
 4    IF (INDEX1(1)) G(1)=X(1)*(.1D+1+X(2)**2)+X(3)**4-.4D+1
     F                    -DSQRT(.18D+2)
      RETURN
 5    IF (.NOT.INDEX2(1)) GOTO 8
      GG(1,1)=.1D+1+X(2)**2
      GG(1,2)=.2D+1*X(1)*X(2)
      GG(1,3)=.4D+1*X(3)**3
 8    RETURN
      END
      
      SUBROUTINE TP346(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
 1    N=3
      NILI=0
      NINL=2
      NELI=0
      NENL=0
      DO 6 I=1,3
      LXU(I)=.TRUE.
      LXL(I)=.TRUE.
 6    XL(I)=.0D+0
      X(1)=.223D+2
      X(2)=.5D+0
      X(3)=.125D+3
      XU(1)=.36D+2
      XU(2)=.5D+1
      XU(3)=.125D+3
      GG(1,3)=.0D+0
      GG(2,2)=.0D+0
      LEX=.FALSE.
      NEX=1
      XEX(1)=.16508383D+02
      XEX(2)=.24768216D+01
      XEX(3)=.12399452D+03
      FEX=-.56847825D+01
      RETURN
 2    FX=-.201D-1*X(1)**4*X(2)*X(3)**2*.1D-6
      RETURN
 3    GF(1)=-.804D-1*X(1)**3*X(2)*X(3)**2*.1D-6
      GF(2)=-.201D-1*X(1)**4*X(3)**2*.1D-6
      GF(3)=-.402D-1*X(1)**4*X(2)*X(3)*.1D-6
      RETURN
 4    IF (INDEX1(1)) G(1)=.675D+3-X(1)**2*X(2)
      IF (INDEX1(2)) G(2)=.419D+0-X(1)**2*X(3)**2*.1D-6
      RETURN
 5    IF (.NOT.INDEX2(1)) GOTO7
      GG(1,1)=-.2D+1*X(1)*X(2)
      GG(1,2)=-X(1)**2
 7    IF (.NOT.INDEX2(2)) GOTO8
      GG(2,1)=-.2D-6*X(1)*X(3)**2
      GG(2,3)=-.2D-6*X(1)**2*X(3)
 8    RETURN
      END
C      
      SUBROUTINE TP347(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D347/H
      DOUBLEPRECISION H(8),A(3)
      DATA(A(I),I=1,3)/0.820437D+4,0.900872D+4,0.933046D+4/
      INTEGER I
      GOTO(1,2,2,4,5),MODE
 1    N=3
      NILI=0
      NINL=0
      NELI=1
      NENL=0
      X(1)=0.7D+0
      X(2)=0.2D+0
      X(3)=0.1D+0
      DO 7  I=1,3
      LXL(I)=.TRUE.
      LXU(I)=.TRUE.
      XL(I)=0.0D+0
      XU(I)=0.1D+1
 7    GG(1,I)=0.1D+1
      LEX=.TRUE.
      NEX=1
      XEX(1)=.00000000D+00
      XEX(2)=.00000000D+00
      XEX(3)=.10000000D+01
      FEX=.17374625D+05
      RETURN
 2    H(1)=X(1)+X(2)+X(3)+0.3D-1
      H(2)=0.9D-1*X(1)+X(2)+X(3)+0.3D-1
      H(3)=H(1)*H(2)
      H(4)=X(2)+X(3)+0.3D-1
      H(5)=0.7D-1*X(2)+X(3)+0.3D-1
      H(6)=H(4)*H(5)
      H(7)=X(3)+0.3D-1
      H(8)=0.13D+0*X(3)+0.3D-1
      IF (MODE .EQ. 3) GOTO 3
      FX = A(1)*DLOG(DMAX1(H(1)/H(2),1.0D-4))
     /   + A(2)*DLOG(DMAX1(H(4)/H(5),1.0D-4))
     /   + A(3)*DLOG(DMAX1(H(7)/H(8),1.0D-4))
      RETURN
 3    GF(1)=A(1)*(H(2)-0.9D-1*H(1))/H(3)
      GF(2)=A(1)*(H(2)-H(1))/H(3)+A(2)*(H(5)-0.7D-1*H(4))/H(6)
      GF(3)=A(1)*(H(2)-H(1))/H(3)+A(2)*(H(5)-H(4))/H(6)
     F        +A(3)*(H(8)-.13D+0*H(7))/(H(7)*H(8))
      RETURN
 4    IF (INDEX1(1)) G(1)=X(1)+X(2)+X(3)-0.1D+1
 5    RETURN
      END            
C     
      SUBROUTINE TP348 (MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D348/AF,AT,AC,GI,RE,XMDOT,DELP,H1,COSTM,COSTT,COSTF,
     /       HO,XVAL,ETAF,ETAS,HEF,Q
      DOUBLEPRECISION AF,AT,AC,GI,RE,XMDOT,DELP,H1,COSTM,COSTT,COSTF,
     /       HO,XVAL,ETAF,ETAS,HEF,Q,XX
      DOUBLEPRECISION RHO,XMU,CP,PR,PI,D,TIN,TSURF,H,W,RHOC,RHOA
      DATA RHO,XMU,CP,PR,PI,D,TIN,TSURF,H,W,RHOC,RHOA/.747D-1,
     / .443D-1,.240D+0,.709D+0,3.14159D+0,.525D+0,75.D+0,45.D+0,
     / 13.13D+0,3.166D+0,559.D+0,169.D+0/  
      INTEGER I
      GOTO (1,2,3,4,5),MODE
 1    N=3
      NILI=0
      NINL=0
      NELI=0
      NENL=1
      DO 6 I=1,3	   
      LXL(I)=.FALSE.
 6    LXU(I)=.TRUE.
      LXL(2)=.TRUE.
      XL(2)=.1313D+2
      XU(1)=.44D-1
      XU(2)=.24D+2
      XU(3)=.6D+3
      X(1)=.04D0
      X(2)=.18D+2
      X(3)=.144D+3
      LEX=.FALSE.
      NEX=1
      XEX(1)=0.44D-01
      XEX(2)=0.24D+02
      XEX(3)=0.85607576D+02
      FEX=0.36970840D+02
      RETURN
2     CONTINUE
      IF (X(1).LT.XL(1)) X(2)=XL(1)
      IF (X(2).LT.XL(2)) X(2)=XL(2)
      IF (X(3).LT.XL(3)) X(2)=XL(3)
      AF=X(2)/X(1)*0.2D+1*(W*H-0.3D+2*PI*D**2/.4D+1)/0.144D+3
      AT=0.3D+2*PI*D*X(2)/0.144D+3
      AC=(H*X(2)-0.1D+2*D*X(2)-X(2)/X(1)*0.6D-2*H)/0.144D+3
      IF (AC.EQ.0.0D0) AC=1.0D-20
      GI=(RHO*X(3)*(H*X(2))/(AC*0.144D+3))*0.6D+2
      RE=GI*.1083D+1/(.12D+2*XMU)
      IF (RE.LT..1D-9) RE=.1D-9
      HO=(.195D+0*GI*CP)/(PR**.67*RE**.35)
      XMDOT=RHO*X(3)*H*X(2)/0.144D+3*0.6D+2
      DELP=0.1833D-5/RHO*GI**2*0.3D+1*(AF/AC*RE**(-0.5)+0.1D+0*AT/AC)
      IF (HO.LT..1D-9) HO=.1D-9
      XVAL=.732D-1*DSQRT(HO)
      ETAF=DTANH(XVAL)/XVAL
      ETAS=.1D+1-AF/(AF+AT)*(0.1D+1-ETAF)
      XX=XMDOT*CP
      HEF=0.1D+1-DEXP(DMAX1(-ETAS*HO*(AF+AT)/XX,-1.0D2))
      Q=HEF*(TIN-TSURF)*XMDOT*CP
      IF (MODE .EQ. 4) GOTO 7
      H1=DELP/RHO*XMDOT/.198D+7
      IF (H1.LT.0.1D-9) H1=.1D-9
      COSTM=DSQRT(H1)/.718D-1+.4D+1
      COSTT=.101D+1*.3D+2*X(2)*PI/.4D+1*(D**2-(D-.36D-1)**2)
      COSTF=.47D+0*H*W*.6D-2*RHOA/.1728D+4*X(2)/X(1)
      COSTT=COSTT*RHOC/.1728D+4
      FX=COSTM+COSTT+COSTF
 3    RETURN
 4    IF (.NOT.INDEX1(1)) GOTO 5
      GOTO 2
 7    G(1)=6.D+3-Q
 5    RETURN
      END
C
      SUBROUTINE TP349(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D349/P1,XK1,XK2,V,C1,UT,ARGU,XLMTD,HEAT,AREA,ARE,HEA,
     /        DIA,PRESS,WATE,C2,C3,C4,C5,C6,C7,COST,VEST,C0,U,
     /        PHI,A11,A12,A13,A22,A23,A31,A32,A33,TEMP1,TEMP2,TEMP3
      DOUBLEPRECISION P1,XK1,XK2,V,C1,UT,ARGU,XLMTD,HEAT,AREA,ARE,HEA,
     /        DIA,PRESS,WATE,C2,C3,C4,C5,C6,C7,COST,VEST,C0,U,
     /        PHI(9),A11,A12,A13,A22,A23,A31,A32,A33,TEMP1,TEMP2,TEMP3
      DOUBLEPRECISION P2,C1F,C2F,H1,H2,E1,E2,CPP,P
      DATA    P2,C1F,C2F,H1,H2,E1,E2,CPP,P/0.1D+2,0.75D-1,0.25D-1,
     /        0.8D+4,0.8D+4,0.1D+4,0.1D+4,0.36938503D+1,0.2D+2/
      INTEGER I
      GOTO (1,2,3,4,3),MODE
    1 N=3
      NILI=0
      NINL=9
      NELI=0
      NENL=0
      X(1)=0.5D+4
      X(2)=0.2D+3
      X(3)=0.1D+3
      DO 6 I=1,2
      LXL(I)=.TRUE.
    6 LXU(I)=.TRUE.
      LXL(3)=.FALSE.
      LXU(3)=.FALSE.
      XL(1)=0.1D+4
      XL(2)=0.1D+3
      XU(1)=0.8D+4
      XU(2)=0.5D+3
      LEX=.FALSE.
      NEX=1
      FEX=-0.41489499D+1
      XEX(1)=0.78287954D+4
      XEX(2)=0.18881406D+3
      XEX(3)=0.11381406D+3
      RETURN
    2 P1=0.1D+3
      DO 35 I=1,3
   35 IF(X(I).LT.0.1D-5) X(I)=0.1D-5
      XK1=P1*DEXP(-E1/(0.46D+3+X(2)))
      XK2=P2*DEXP(-E2/(0.46D+3+X(2)))
      V=P*X(1)/(XK2*(X(1)*C2F-P))
      C1=(X(1)*C1F-P)/(X(1)+V*XK1)
      UT=0.43D+2+0.452D-1*X(2)
   39 ARGU=(X(2)-X(3)-0.75D+2)/(X(2)-0.1D+3)
      IF(ARGU.EQ.0.0D+0) GOTO 48
      XLMTD=(0.25D+2-X(3))/DLOG(DABS(ARGU))
      HEAT=X(1)*CPP*(0.1D+3-X(2))+XK1*(X(1)*C1F-P)*V*H1/(X(1)+V*XK1)
     /    +P*H2
      AREA=HEAT/(UT*XLMTD)
      ARE=DABS(AREA)
      HEA=DABS(HEAT)
      DIA=(V/0.1272D+2)**0.33333333
      IF(X(2).LT.0.2D+3) GOTO 40
      PRESS=0.236D+2+0.33D-5*(X(2)**3)
      GOTO 41
   48 X(2)=X(2)*0.10001D+1
      GOTO 39
   40 PRESS=0.5D+2
   41 WATE=(0.909D-1*(DIA**3)+0.482D+0*(DIA**2))*PRESS+0.366D+2*(DIA**2)
     /     +0.1605D+3*DIA
      C1 =0.48D+1*(WATE**0.782)
      IF(X(2).LT.0.2D+3) GOTO 42
      C2=(0.172D+2+0.133D-1*X(2))*DIA**2
      GOTO 43
   42 C2=0.0D+0
   43 IF(PRESS.LT.0.15D+3) GOTO 44
      C3=0.27D+3*(ARE**0.546)*(0.962D+0+0.168D-6*(X(2)**3))
      GOTO 45
   44 C3=0.27D+3*(ARE**0.546)
   45 C4=0.14D+4+0.14D+3*DIA
      C5=0.875D+3*((0.5D-1*V)**0.3)
      C6=0.812D+3*(((0.695D-3+0.459D-10*(X(2)**3))+X(1))**0.467)
      IF(X(2).LT.0.25D+3) GOTO 46
      C7=0.1291D+4*((0.298D+3*HEA/X(3))**0.467)
      GOTO 47
   46 C7=0.812D+3*((0.298D+3*HEA/X(3))**0.467)
   47 COST=C1+C2+C3+C4+C5+C6+C7
      VEST=0.5D+1*COST
      C0=0.22D+5+0.18D+0*VEST+0.31D+1*V+0.611D+2*((0.695D-3+0.459D-10
     F   *(X(2)**3))*X(1))+0.115D-2*HEAT+0.692D+1*HEAT+0.574D+3
     F   *X(1)*(C1F-C1)+0.1148D+6
      FX=(0.688D+6-C0)/(0.2D+1*VEST)*(-0.1D-2)
    3 RETURN
    4 P1=0.1D+3
      DO 37 I=1,3
   37 IF(X(I).LT.0.1D-5) X(I)=0.1D-5
      XK1=P1*DEXP(-E1/(0.46D+3+X(2)))
      XK2=P2*DEXP(-E2/(0.46D+3+X(2)))
      V=P*X(1)/(XK2*(X(1)*C2F-P))
      C1=(X(1)*C1F-P)/(X(1)+V*XK1)
      UT=0.43D+2+0.452D-1*X(2)
   36 ARGU=(X(2)-X(3)-0.75D+2)/(X(2)-0.1D+3)
      IF(ARGU.EQ.0.0D+0) GOTO 49
      XLMTD=(0.25D+2-X(3))/DLOG(DABS(ARGU))
      HEAT=X(1)*CPP*(0.1D+3-X(2))+XK1*(X(1)*C1F-P)*V*H1/(X(1)+V*XK1)
     /     +P*H2
      AREA=HEAT/(UT*XLMTD)
      DIA=(V/0.1272D+2)**0.33333333
      IF(X(2).LT.0.2D+3) GOTO 50
      PRESS=0.236D+2+0.33D-5*(X(2)**3)
      GOTO 51
   49 X(2)=X(2)*0.10001D+1
      GOTO 36
   50 PRESS=0.5D+2
   51 PHI(1)=DIA-0.125D+1
      PHI(2)=0.967D+1-DIA
      PHI(3)=AREA-0.5D+2
      PHI(4)=0.4D+4-AREA
      A11=XK1+X(1)/V
      A12=XK2
      A13=(X(1)*C1F-PRESS)*XK1*E1/((X(1)+V*XK1)*((X(2)+0.46D+3)**2))
     F    +PRESS*E2/(V*((X(2)+0.46D+3)**2))
      A22=XK2+X(1)/V
      A23=PRESS*E2/(V*((X(2)+0.46D+3)**2))
      A31=-H1*XK1/CPP
      A32=-H2*XK2/CPP
      A33=X(1)/V+UT*AREA/(V*CPP)-(X(1)*C1F-PRESS)*XK1*E1*H1/((X(1)+V*XK1
     F    )*CPP*((X(2)+0.46D+3)**2))-PRESS*E2*H2/(V*CPP*((X(2)
     F    +0.46D+3)**2))
      TEMP1=A11+A22+A33
      PHI(5)=TEMP1
      TEMP2=A11*A22+A22*A33+A33*A11-A13*A31-A23*A32
      PHI(6)=TEMP2
      TEMP3=A11*A22*A33+A12*A23*A31-A13*A31*A22-A23*A32*A11
      PHI(7)=TEMP3
      PHI(8)=TEMP1*TEMP2-TEMP3
      PHI(9)=HEAT
      DO 7 I=1,9
    7 IF (INDEX1(I)) G(I)=PHI(I)
      RETURN
      END
C      
      SUBROUTINE TP350(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J,LSUM
      DOUBLEPRECISION  F,DF,Y(11),U(11),H(11)
      DATA (Y(I),I=1,11)/0.1957D+0,0.1947D+0,0.1735D+0,0.1600D+0,
     F                   0.844D-1,0.627D-1,0.456D-1,0.342D-1,
     F                   0.323D-1,0.235D-1,0.246D-1/
     F     (U(I),I=1,11)/0.4D+1,0.2D+1,0.1D+1,0.5D+0,0.25D+0,0.167D+0,
     F                   0.125D+0,0.1D+0,0.833D-1,0.714D-1,0.625D-1/
      DO 12 I=1,11
   12 H(I)=U(I)**2+X(3)*U(I)+X(4)
      GOTO (1,2,2,4,4),MODE
    1 N=4
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      X(1)=0.25D+0
      X(2)=0.39D+0
      X(3)=0.415D+0
      X(4)=0.39D+0
      DO 6 I=1,4
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LSUM=11
      LEX=.FALSE.
      NEX=1
      FEX=0.30750560D-3
      XEX(1)=0.19280644D+0
      XEX(2)=0.19126279D+0
      XEX(3)=0.12305098D+0
      XEX(4)=0.13605235D+0
      RETURN
    2 DO 20 I=1,11
   20 F(I)=Y(I)-X(1)/H(I)*(U(I)**2+X(2)*U(I))
      IF (MODE.EQ.3) GOTO 3
      FX=0.0D+0
      DO 7 I=1,11
    7 FX=FX+F(I)**2
      RETURN
    3 DO 8 I=1,11
      DF(I,1)=(-U(I)**2-X(2)*U(I))/H(I)
      DF(I,2)=(-X(1)*U(I))/H(I)
      DF(I,3)=X(1)*U(I)*(U(I)**2+X(2)*U(I))/H(I)**2
    8 DF(I,4)=X(1)*(U(I)**2+X(2)*U(I))/H(I)**2
      DO 9 J=1,4 
      GF(J)=0.0D+0
      DO 10 I=1,11
   10 GF(J)=GF(J)+2.D+0*F(I)*DF(I,J)
    9 CONTINUE
    4 RETURN
      END
C
      SUBROUTINE TP351(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J,LSUM
      DOUBLEPRECISION XH1,XH2,XH3,XH4,XH5,F,DF,A(7),B(7)
      DATA(A(I),I=1,7)/0.0D+0,0.428D-3,0.1D-2,0.161D-2,0.209D-2,
     F                 0.348D-2,0.525D-2/
     F    (B(I),I=1,7)/0.7391D+1,0.1118D+2,0.1644D+2,0.162D+2,
     F                 0.222D+2,0.2402D+2,0.3132D+2/
      XH1=X(1)**2
      XH2=X(2)**2
      XH3=X(3)**2
      XH4=X(4)**2
      GOTO (1,2,2,4,4),MODE
    1 N=4
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      X(1)=0.27D+1
      X(2)=0.9D+2
      X(3)=0.15D+4
      X(4)=0.1D+2
      DO 6 I=1,4
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
C      FEX=0.31858082D+3
      FEX=0.31857175D+3  
      XEX(1)=0.27143661D+1
      XEX(2)=0.14043580D+3
      XEX(3)=0.17075155D+4
      XEX(4)=0.31512867D+2
      LSUM=7
      RETURN
    2 DO 20 I=1,7
   20 F(I)=(((XH1+A(I)*XH2+A(I)*A(I)*XH3)/(0.1D+1+A(I)*XH4))
     F      -B(I))/B(I)*1.D+2
      IF (MODE.EQ.3) GOTO 3
      FX=0.0D+0
      DO 7 I=1,7
    7 FX=FX+F(I)**2
      RETURN
    3 DO 8 J=1,4
    8 GF(J)=0.0D+0
      DO 10 I=1,7
      XH5=(0.1D+1+XH4*A(I))*B(I)
      DF(I,1)=0.2D+3*X(1)/XH5
      DF(I,2)=0.2D+3*X(2)*A(I)/XH5
      DF(I,3)=0.2D+3*X(3)*A(I)**2/XH5
      DF(I,4)=-0.2D+3*X(4)*A(I)*B(I)*(XH1+XH2*A(I)+XH3*A(I)**2)/XH5**2
      DO 10 J=1,4    
   10 GF(J)=GF(J)+2.D+0*F(I)*DF(I,J)
    4 RETURN
      END
C      
      SUBROUTINE TP352(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J,LSUM
      DOUBLEPRECISION F,DF,TI
      GOTO (1,2,2,4,4),MODE
    1 N=4
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      X(1)=0.25D+2
      X(2)=0.5D+1
      X(3)=-0.5D+1
      X(4)=-0.1D+1
      DO 6 I=1,4
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=0.90323433D+3
      XEX(1)=-0.10223574D+2
      XEX(2)=0.11908429D+2
      XEX(3)=-0.45804134D+0
      XEX(4)=0.58031996D+0
      LSUM=40
      RETURN
    2 DO 20 I=1,20
      TI=I*0.2D+0
      F(I)=X(1)+X(2)*TI-DEXP(TI)
   20 F(I+20)=X(3)+X(4)*DSIN(TI)-DCOS(TI)
      IF (MODE.EQ.3) GOTO 3
      FX=0.0D+0
      DO 7 I=1,20
    7 FX=FX+F(I)**2+F(I+20)**2
      RETURN
    3 DO 8 J=1,4
    8 GF(J)=0.0D+0
      DO 10 I=1,20
      TI=I*0.2D+0
      DF(I,1)=0.1D+1
      DF(I+20,1)=0.0D+0
      DF(I,2)=TI
      DF(I+20,2)=0.0D+0
      DF(I,3)=0.0D+0
      DF(I+20,3)=0.1D+1
      DF(I,4)=0.0D+0
      DF(I+20,4)=DSIN(TI)
      DO 10 J=1,4
  10  GF(J)=GF(J)+2.D+0*F(I)*DF(I,J)+2.D+0*F(I+20)*DF(I+20,J)
    4 RETURN
      END
C      
      SUBROUTINE TP353(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION Q
      GOTO (1,2,3,4,4),MODE
    1 N=4
      NILI=1
      NINL=1
      NELI=1
      NENL=0
      X(1)=0.0D+0
      X(2)=0.0D+0
      X(3)=0.4D+0
      X(4)=0.6D+0
      GG(1,1)=0.23D+1
      GG(1,2)=0.56D+1
      GG(1,3)=0.111D+2
      GG(1,4)=0.13D+1
      DO 6 I=1,4
      GG(3,I)=0.1D+1
      LXL(I)=.TRUE.
      LXU(I)=.FALSE.
    6 XL(I)=0.0D+0
      LEX=.FALSE.
      NEX=1
      FEX=-0.39933673D+2
      XEX(1)=0.D+0
      XEX(2)=0.D+0
      XEX(3)=0.37755102D+0
      XEX(4)=0.62244898D+0
      RETURN
    2 FX=0.2455D+2*X(1)+0.2675D+2*X(2)+0.39D+2*X(3)+0.405D+2*X(4)
      FX=-FX
      RETURN
    3 GF(1)=-0.2455D+2
      GF(2)=-0.2675D+2
      GF(3)=-0.39D+2
      GF(4)=-0.405D+2
      RETURN
    4 Q=(0.53D+0*X(1))**2+(0.44D+0*X(2))**2+(0.45D+1*X(3))**2
     F  +(0.79D+0*X(4))**2
      IF (MODE.EQ.5) GOTO 5
      IF (INDEX1(1)) G(1)=0.23D+1*X(1)+0.56D+1*X(2)+0.111D+2*X(3)
     F                    +0.13D+1*X(4)-0.5D+1
      IF (INDEX1(2)) G(2)=0.12D+2*X(1)+0.119D+2*X(2)+0.418D+2*X(3)
     F               +0.521D+2*X(4)-0.1645D+1*DSQRT(Q)-0.12D+2
      IF (INDEX1(3)) G(3)=X(1)+X(2)+X(3)+X(4)-0.1D+1
      RETURN
    5 IF (.NOT. INDEX2(2)) GOTO 7
      GG(2,1)=0.12D+2-0.1645D+1*0.53D+0**2*X(1)/DSQRT(Q)
      GG(2,2)=0.119D+2-0.1645D+1*0.44D+0**2*X(2)/DSQRT(Q)
      GG(2,3)=0.418D+2-0.1645D+1*0.45D+1**2*X(3)/DSQRT(Q)
      GG(2,4)=0.521D+2-0.1645D+1*0.79D+0**2*X(4)/DSQRT(Q)
    7 RETURN
      END
C      
      SUBROUTINE TP354(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=4  
      NILI=1
      NINL=0
      NELI=0
      NENL=0
      X(1)=0.3D+1
      X(2)=-0.1D+1
      X(3)=0.0D+0
      X(4)=0.1D+1
      DO 6 I=1,4
      GG(1,I)=0.1D+1
      LXU(I)=.TRUE.
      LXL(I)=.FALSE.
    6 XU(I)=0.2D+2
      LEX=.FALSE.
      NEX=1
      FEX=0.11378385D+0
      XEX(1)=0.50336521D+0
      XEX(2)=-0.45569070D-1
      XEX(3)=0.23580504D+0
      XEX(4)=0.30639882D+0
      RETURN
    2 FX=(X(1)+0.1D+2*X(2))**2+0.5D+1*(X(3)-X(4))**2+(X(2)-0.2D+1
     F   *X(3))**4+0.1D+2*(X(1)-X(4))**4
      RETURN
    3 GF(1)=0.2D+1*X(1)+0.2D+2*X(2)+0.4D+2*(X(1)-X(4))**3
      GF(2)=0.2D+2*X(1)+0.2D+3*X(2)+0.4D+1*(X(2)-0.2D+1*X(3))**3
      GF(3)=0.1D+2*X(3)-0.1D+2*X(4)-0.8D+1*(X(2)-0.2D+1*X(3))**3
      GF(4)=-0.1D+2*X(3)+0.1D+2*X(4)-0.4D+2*(X(1)-X(4))**3
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)+X(2)+X(3)+X(4)-0.1D+1
    5 RETURN
      END
C
      SUBROUTINE TP355(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION R(4),H0,H1,H2,H3,H4,H5,H6
      INTEGER I
      H0=X(2)*X(4)
C      R(1)=0.11D+2-X(1)*X(4)-H0+X(3)*X(4)
C      R(2)=X(1)+0.1D+2*X(2)-X(3)+X(4)+H0*(X(3)-X(1))
C      R(3)=0.11D+2-0.4D+1*X(1)*X(4)-0.4D+1*H0+X(3)*X(4)
C      R(4)=0.2D+1*X(1)+0.2D+2*X(2)-0.5D+0*X(3)+0.2D+1*X(4)+0.2D+1*H0*
C     F(X(3)-0.4D+1*X(1))
      R(1)=11.0D0-(X(1)+X(2)-X(3))*X(4)
      R(2)=X(1)+10.0D0*X(2)-X(3)+X(4)+H0*(X(3)-X(1))
      R(3)=11.0D0-(4.0D0*X(1)+4.0D0*X(2)- X(3))*X(4)
      R(4)=2.0D0*X(1)+20.0D0*X(2)-0.5D0*X(3)+2.0D0*X(4) 
     /     + 2.0D0*H0*(X(3)-4.0D0*X(1))
      H1=0.2D+1*R(1) 
      H2=0.2D+1*R(2)
      H3=0.2D+1*R(3)
      H4=0.2D+1*R(4)
      H5=H1*X(4)
      H6=H2*(0.1D+1-H0)
      GOTO(1,2,3,4,5),MODE
    1 N=4
      NILI=0
      NINL=0
      NELI=0
      NENL=1
      DO 6 I=1,4
      LXU(I)=.FALSE.
      LXL(I)=.TRUE.
    6 X(I)=0.0D+0
      X(1)=0.1D0
      X(2)=0.1D0
      XL(1)=0.1D+0
      XL(2)=0.1D+0
      XL(3)=0.0D+0
      XL(4)=0.0D+0
      LEX=.FALSE.
      NEX=1
      FEX=0.69675463D+2
      XEX(1)=0.19166330D+1
      XEX(2)=0.10000000D+0
      XEX(3)=0.00000000D+0
      XEX(4)=0.19718118D+1
      RETURN
    2 FX=R(1)*R(1)+R(2)*R(2)
      RETURN
    3 GF(1)=-H5+H6
      GF(2)=-H5+H2*(0.1D+2+X(4)*(X(3)-X(1)))
      GF(3)=H5-H6
      GF(4)=H1*(-X(1)-X(2)+X(3))+H2*(0.1D+1+X(2)*(X(3)-X(1)))
      RETURN
    4 IF (INDEX1(1)) G(1)=R(1)**2+R(2)**2-R(3)**2-R(4)**2
      RETURN   
    5 IF (.NOT.INDEX2(1)) GOTO 7
      GG(1,1)=-H5+H6-H3*(-0.4D+1)*X(4)-H4*(0.2D+1-0.8D+1*H0)
      GG(1,2)=-H5+H2*(0.1D+2+X(4)*(X(3)-X(1)))-H3*(-0.4D+1)*X(4)-
     F   H4*(0.2D+2+0.2D+1*X(4)*(X(3)-0.4D+1*X(1)))
      GG(1,3)=H5-H6-H3*X(4)-H4*(-0.5D+0+0.2D+1*H0)
      GG(1,4)=H1*(-X(1)-X(2)+X(3))+H2*(0.1D+1+X(2)*(X(3)-X(1)))-H3*
     F   (-0.4D+1*X(1)-0.4D+1*X(2)+X(3))-H4*(0.2D+1+0.2D+1*X(2)*(X(3)
     F   -0.4D+1*X(1)))
    7 RETURN
      END
C
      SUBROUTINE TP356(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION L,LOAD,TD,SIGD,FH,T1,M,R,J,T2,COSA,WP,
     /        T,SIG,E,EI,GH,GJ,EITC,EIDC,REITC,REIDC,PC,DEL,PHI(5)
      INTEGER I
      GOTO(1,2,3,4,3),MODE
    1 N=4
      NILI=1
      NINL=4
      NENL=0
      NENL=0
      X(1)=0.1D+1
      X(2)=0.7D+1
      X(3)=0.8D+1
      X(4)=0.1D+1
      DO 6 I=1,3
      LXU(I)=.FALSE.
    6 LXL(I)=.TRUE.
      LXL(4)=.FALSE.
      LXU(4)=.FALSE.
      XL(1)=0.125D+0
      XL(2)=0.0D+0
      XL(3)=0.0D+0
      LEX=.FALSE.
      NEX=1
      FEX=0.23811648D+1
      XEX(1)=0.24436898D+0
      XEX(2)=0.62187934D+1
      XEX(3)=0.82914714D+1
      XEX(4)=0.24436898D+0
      RETURN
    2 FX=0.110471D+1*X(1)*X(1)*X(2)+0.4811D-1*X(3)*X(4)*(0.14D+2
     /   +X(2))
    3 RETURN
    4 L=0.14D+2
      LOAD=0.6D+4
      TD=0.136D+5
      SIGD=0.3D+5
      FH=LOAD
      T1=FH/(0.1414D+1*X(1)*X(2))
      M=FH*(L+(X(2)/0.2D+1))
      R=DSQRT((X(2)*X(2)/0.4D+1)+((X(3)+X(1))/0.2D+1)**2)
      J=0.2D+1*(0.707D+0*X(1)*X(2)*((X(2)*X(2)/0.12D+2)+((X(3)
     /   +X(1))/0.2D+1)**2))
      T2=M*R/J
      COSA=X(2)/(0.2D+1*R)
      WP=DABS(T1*T1+0.2D+1*T1*T2*COSA+T2*T2)
      IF (WP.GT.0.0) THEN
         T=DSQRT(WP)
      ELSE
         T=0.0
      ENDIF      
      SIG=0.6D+1*FH*L/(X(4)*X(3)*X(3))
      PHI(1)=(TD-T)/0.1D+5
      PHI(2)=(SIGD-SIG)/0.1D+5
      PHI(3)=X(4)-X(1)
      E=0.3D+8
      EI=E*X(3)*X(4)*X(4)*X(4)/0.12D+2
      GH=0.12D+8
      GJ=GH*X(3)*X(4)*X(4)*X(4)/0.3D+1
      EITC=EI*GJ
      EIDC=EI/GJ
      IF (EITC.GT.0.0) THEN
         REITC=DSQRT(EITC)
      ELSE
         REITC=0.0
      ENDIF      
      IF (EIDC.GT.0.0) THEN
         REIDC=DSQRT(EIDC)
      ELSE
         REIDC=0.0
      ENDIF      
      PC=0.4013D+1*REITC*(0.1D+1-(X(3)/(0.2D+1*L))*REIDC)/(L*L)
      PHI(4)=(PC-FH)/0.1D+5
      DEL=0.4D+1*FH*L*L*L/(E*X(4)*X(3)*X(3)*X(3))
      PHI(5)=0.25D+0-DEL
      IF (INDEX1(1)) G(1)=PHI(3)
      IF (INDEX1(2)) G(2)=PHI(1)
      IF (INDEX1(3)) G(3)=PHI(2)
      IF (INDEX1(4)) G(4)=PHI(4)
      IF (INDEX1(5)) G(5)=PHI(5)
      RETURN
      END
C
      SUBROUTINE TP357(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      COMMON     /D357/DALPHA,SUM,P1,Q1,R1,S1,ALPHA,CA,SA,PI,QI,
     /           A,B,C,AABB,TEST,AJ,PH,SP,CP,RI,SI,TEST1,CALCX,CALCY,SQL
      DOUBLEPRECISION DALPHA,SUM,P1,Q1,R1,S1,ALPHA,CA,SA,PI,QI,
     /           A,B,C,AABB,TEST,AJ,PH,SP,CP,RI,SI,TEST1,CALCX,CALCY,SQL
      DOUBLEPRECISION XPT(36),YPT(36),P0,Q0,R0,S0
      DATA(XPT(KI),KI=1,36)/0.113D+3,0.1101D+3,0.1062D+3,0.1013D+3,
     /             0.954D+2,0.888D+2,0.816D+2,0.74D+2,0.661D+2,
     /             0.584D+2,0.51D+2,0.443D+2,0.387D+2,0.345D+2,
     /             0.324D+2,0.329D+2,0.364D+2,0.428D+2,0.509D+2,
     /             0.59D+2,0.658D+2,0.715D+2,0.765D+2,0.811D+2,
     /             0.856D+2,0.902D+2,0.946D+2,0.989D+2,0.103D+3,
     /             0.1067D+3,0.1099D+3,0.1125D+3,0.1144D+3,
     /             0.1155D+3,0.1157D+3,0.1149D+3/,
     /    (YPT(KJ),KJ=1,36)/0.402D+2,0.468D+2,0.533D+2,0.594D+2,
     /             0.65D+2,0.699D+2,0.739D+2,0.769D+2,0.789D+2,
     /             0.798D+2,0.797D+2,0.785D+2,0.765D+2,0.736D+2,
     /             0.702D+2,0.66D+2,0.609D+2,0.543D+2,0.458D+2,
     /             0.361D+2,0.265D+2,0.181D+2,0.114D+2,0.62D+1,
     /             0.26D+1,0.3D+0,-0.7D+0,-0.6D+0,0.7D+0,0.31D+1,
     /             0.64D+1,0.105D+2,0.155D+2,0.21D+2,0.271D+2,
     /             0.336D+2/
      DATA P0,Q0,R0,S0/0.9D+2,0.0D+0,0.0D+0,0.0D+0/
      INTEGER I,K
      GOTO(1,2,3,3,3),MODE
    1 N=4
      NILI=0
      NINL=0
      NENL=0
      NENL=0
      X(1)=0.136D+3
      X(2)=0.0D+0
      X(3)=0.748D+2
      X(4)=0.757D+2
      DO 6 I=1,4
      LXL(I)=.TRUE.
    6 LXU(I)=.TRUE.
      LEX=.FALSE.
      NEX=1
      FEX=0.35845660D+0
      XEX(1)=0.13600762D+3
      XEX(2)=0.31371415D-1
      XEX(3)=0.73594390D+2
      XEX(4)=0.72187426D+2
      DO 7 K=1,4
    7 XL(K)=0.0D+0      
      XU(1)=0.150D+3
      XU(2)=0.5D+2
      XU(3)=0.1D+3
      XU(4)=0.1D+3
      RETURN
    2 DALPHA=0.3141527D+1/0.18D+2
      SUM=0.0D+0
      P1=X(1)
      Q1=X(2)
      R1=X(3)
      S1=X(4)
      DO 54 I=2,36
      ALPHA=DALPHA*(I-1)
      CA=DCOS(ALPHA)
      SA=DSIN(ALPHA)
      PI=P1*CA-Q1*SA+P0*(0.1D+1-CA)+Q0*SA
      QI=P1*SA+Q1*CA+Q0*(0.1D+1-CA)-P0*SA
      A=R0*S1-S0*R1-Q1*R0+P1*S0+PI*Q1-P1*QI+QI*R1-PI*S1
      B=-R0*R1-S0*S1+P1*R0+Q1*S0-P1*PI-Q1*QI+PI*R1+QI*S1
      C=-R1*R0-S1*S0+PI*R0+QI*S0+P1*R1+Q1*S1-(P1*P1+Q1*Q1+PI*PI+QI
     F*QI)/0.2D+1
      AABB=A*A+B*B
      AJ=0.1D+1
      TEST=0.0D0
      IF (AABB.LT.0.1D-29) THEN 
         IF (DABS(A).LT.1.0D-29) A = 1.0D-29
         GOTO 50
      ENDIF   
      TEST=C/DSQRT(AABB)
      IF(DABS(TEST).GT.0.1D+1) GOTO 51
   52 PH=DASIN(TEST)-DATAN(B/A)
   55 SP=DSIN(PH)
      CP=DCOS(PH)
      RI=R1*CP-S1*SP+PI-P1*CP+Q1*SP
      SI=R1*SP+S1*CP+QI-P1*SP-Q1*CP
      TEST1=(R1-R0)**2+(S1-S0)**2
      IF(TEST1.LT.0.1D-9) TEST1=0.1D-9
      IF(DABS((TEST1-(RI-R0)**2-(SI-S0)**2)/TEST1).LT.0.1D-2)GOTO53
      IF(AJ.EQ.0.2D+1) GOTO 51
      TEST=-TEST
      AJ=0.2D+1
      GOTO 52
   50 PH=-DATAN(B/A)
      GOTO 55
   51 FX=0.1D+21
      RETURN
   53 CALCX=XPT(1)*CP-YPT(1)*SP+PI-P1*CP+Q1*SP
      CALCY=XPT(1)*SP+YPT(1)*CP+QI-P1*SP-Q1*CP
   54 SUM=SUM+(CALCX-XPT(I))**2+(CALCY-YPT(I))**2
      SQL=(R1-R0)**2+(S1-S0)**2+(R1-P1)**2+(S1-Q1)**2+(P1-P0)**2+
     F(Q1-Q0)**2
      FX=SUM/0.1D+3+SQL/0.625D+5
    3 RETURN
      END
C
      SUBROUTINE TP358(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION Y(33),TI,F,DF
      INTEGER I,J,LSUM
      DATA(Y(I),I=1,33)/0.844D+0,0.908D+0,0.932D+0,0.936D+0,
     F          0.925D+0,0.908D+0,0.881D+0,0.850D+0,0.818D+0,
     F          0.784D+0,0.751D+0,0.718D+0,0.685D+0,0.658D+0,
     F          0.628D+0,0.603D+0,0.580D+0,0.558D+0,0.538D+0,
     F          0.522D+0,0.506D+0,0.490D+0,0.478D+0,0.467D+0,
     F          0.457D+0,0.448D+0,0.438D+0,0.431D+0,0.424D+0,
     F          0.420D+0,0.414D+0,0.411D+0,0.406D+0/
      GOTO(1,2,2,4,4),MODE
    1 N=5
      NILI=0
      NINL=0
      NENL=0
      NENL=0
      X(1)=0.5D+0
      X(2)=0.15D+1
      X(3)=-0.1D+1
      X(4)=0.1D-1
      X(5)=0.2D-1
      DO 6 I=1,5
      LXU(I)=.TRUE.
    6 LXL(I)=.TRUE.
      XL(1)=-0.5D+0
      XU(1)=0.5D+0
      XL(2)=.15D+1
      XU(2)=.25D+1
      XL(3)=-0.2D+1
      XU(3)=-0.1D+1
      XL(4)=0.1D-2
      XU(4)=0.1D+0
      XL(5)=0.1D-2
      XU(5)=0.1D+0
      LEX=.FALSE.
      NEX=1
      FEX=0.546D-4
      XEX(1)=0.3754D+0
      XEX(2)=0.19358D+1
      XEX(3)=-0.14647D+1
      XEX(4)=0.1287D-1
      XEX(5)=0.2212D-1
      LSUM=33
      RETURN
    2 CONTINUE
      DO I=1,5
      IF (X(I).GT.XU(I)) X(I)=XU(I)
      IF (X(I).LT.XL(I)) X(I)=XL(I)
      ENDDO 
      DO 20 I=1,33
      TI=DBLE((I-1))*0.1D+2
   20 F(I)=Y(I)-(X(1)+X(2)*DEXP(-X(4)*TI)+X(3)*DEXP(-X(5)*TI))
      IF (MODE.EQ.3) GOTO 3
      FX=0.0D+0
      DO 7 I=1,33
    7 FX=FX+F(I)**2
      RETURN
    3 DO 8 J=1,5
    8 GF(J)=0.0D+0
      DO 10 I=1,33
      TI=DBLE((I-1))*0.1D+2
      DF(I,1)=-0.1D+1
      DF(I,2)=-DEXP(-X(4)*TI)
      DF(I,3)=-DEXP(-X(5)*TI)
      DF(I,4)=X(2)*DEXP(-X(4)*TI)*TI
      DF(I,5)=X(3)*DEXP(-X(5)*TI)*TI
      DO 10 J=1,5
   10 GF(J)=GF(J)+2.D+0*F(I)*DF(I,J)
    4 RETURN
      END
C
      SUBROUTINE TP359(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION A(5),B(5),C(5),D(5),H(6)            
      DATA (A(I),I=1,5)/-0.8720288849D+7,0.1505125253D+6,
     /                  -0.1566950325D+3,0.4764703222D+6,
     /                  0.7294828271D+6/
     /     (B(I),I=1,5)/-0.145421402D+6,0.29311506D+4,-0.40427932D+2,
     /                 0.5106192D+4,0.1571136D+5/
     /     (C(I),I=1,5)/-0.1550111084D+6,0.436053352D+4,0.129492344D+2
     /                  ,0.10236884D+5,0.13176786D+5/
     /     (D(I),I=1,5)/-0.3266695104D+6,0.739068412D+4,
     /                  -0.278986976D+2,0.16643076D+5,0.30988146D+5/
      INTEGER I,J
      GOTO (1,2,3,4,5),MODE
    1 N=5
      NILI=14
      NINL=0
      NELI=0
      NENL=0
      X(1)=0.252D+1
      X(2)=0.504D+1
      X(3)=0.945D+2
      X(4)=0.2331D+2
      X(5)=0.17136D+2
      DO 10 I=1,5
   10 GF(I)=-A(I)
      DO 6 J=1,8
      DO 6 I=1,5
    6 GG(J,I)=0.0D+0
      DO 7 I=1,4
      GG(2*I-1,I+1)=-0.1D+1
    7 GG(2*I,I+1)=0.1D+1
      GG(1,1)=0.24D+1
      GG(2,1)=-0.12D+1
      GG(3,1)=0.6D+2
      GG(4,1)=-0.2D+2
      GG(5,1)=0.93D+1
      GG(6,1)=-0.9D+1
      GG(7,1)=0.7D+1
      GG(8,1)=-0.65D+1
      DO 8 I=1,5
      GG(9,I)=B(I)
      GG(10,I)=C(I)
      GG(11,I)=D(I)
      GG(12,I)=-B(I)
      GG(13,I)=-C(I)
      GG(14,I)=-D(I)
      LXU(I)=.FALSE.
    8 LXL(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=-0.52804168D+7
      XEX(1)=0.45374D+1
      XEX(2)=0.10889D+2
      XEX(3)=0.27224D+3
      XEX(4)=0.42198D+2
      XEX(5)=0.31762D+2
      RETURN
    2 FX=-0.24345D+5
      DO 9 I=1,5
    9 FX=FX+A(I)*X(I)
      FX=-FX
    3 RETURN
    4 IF (INDEX1(1)) G(1)=0.24D+1*X(1)-X(2)
      IF (INDEX1(2)) G(2)=-0.12D+1*X(1)+X(2)
      IF (INDEX1(3)) G(3)=0.6D+2*X(1)-X(3)
      IF (INDEX1(4)) G(4)=-0.2D+2*X(1)+X(3)
      IF (INDEX1(5)) G(5)=0.93D+1*X(1)-X(4)
      IF (INDEX1(6)) G(6)=-0.9D+1*X(1)+X(4)
      IF (INDEX1(7)) G(7)=0.7D+1*X(1)-X(5)
      IF (INDEX1(8)) G(8)=-0.65D+1*X(1)+X(5)
      DO 11 I=1,3
   11 H(I)=0.0D+0 
      H(4)=0.294D+6 
      H(5)=0.294D+6
      H(6)=0.2772D+6
      DO 12 I=1,5
      H(1)=H(1)+B(I)*X(I)
      H(2)=H(2)+C(I)*X(I)
      H(3)=H(3)+D(I)*X(I)
      H(4)=H(4)-B(I)*X(I)
      H(5)=H(5)-C(I)*X(I)
   12 H(6)=H(6)-D(I)*X(I)
      DO 13 I=9,14
   13 IF (INDEX1(I)) G(I)=H(I-8)
    5 RETURN
      END
C
      SUBROUTINE TP360(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION  C(10),H,HH(5)
      DATA (C(I),I=1,10)/-0.8720288849D+7,0.1505125233D+6,
     F     -0.1566950325D+3,0.4764703222D+6,0.7294828271D+6,
     F     -0.3266695104D+6,0.739068412D+4,-0.278986976D+2,
     F     0.16643076D+5,0.30988146D+5/
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=5  
      NILI=0
      NINL=2
      NELI=0
      NENL=0
      X(1)=0.252D+1
      X(2)=0.2D+1
      X(3)=0.375D+2
      X(4)=0.925D+1
      X(5)=0.68D+1
      LXL(1)=.TRUE.
      LXU(1)=.FALSE.
      DO 6 I=2,5
      LXL(I)=.TRUE.
    6 LXU(I)=.TRUE.
      XL(1)=0.0D+0
      XL(2)=0.12D+1
      XL(3)=0.2D+2
      XL(4)=0.9D+1
      XL(5)=0.65D+1
      XU(2)=0.24D+1
      XU(3)=0.6D+2
      XU(4)=0.93D+1
      XU(5)=0.7D+1
      LEX=.FALSE.
      NEX=1
      FEX=-0.52803351D+7
c      FEX=-0.52803351D0
      XEX(1)=0.4537431D+1
      XEX(2)=0.24D+1
      XEX(3)=0.6D+2
      XEX(4)=0.93D+1
      XEX(5)=0.7D+1
      RETURN
    2 FX=(-C(1)-C(2)*X(2)-C(3)*X(3)-C(4)*X(4)-C(5)*X(5))*X(1)+
     F   0.24345D+5       
      RETURN
    3 GF(1)=-C(1)-C(2)*X(2)-C(3)*X(3)-C(4)*X(4)-C(5)*X(5)
      DO 7 I=2,5
      GF(I)=-C(I)*X(1)
    7 CONTINUE
      RETURN
    4 H=(C(6)+C(7)*X(2)+C(8)*X(3)+C(9)*X(4)+C(10)*X(5))*X(1)
      IF (INDEX1(1)) G(1)=H
      IF (INDEX1(2)) G(2)=0.2772D+6-H
      RETURN
    5 HH(1)=C(6)+C(7)*X(2)+C(8)*X(3)+C(9)*X(4)+C(10)*X(5)
      DO 8 I=2,5     
    8 HH(I)=C(I+5)*X(1)  
      IF (.NOT.INDEX2(1)) GOTO 11
      DO 9 I=1,5
    9 GG(1,I)=HH(I)
   11 IF (.NOT.INDEX2(2)) GOTO 12
      DO 10 I=1,5
   10 GG(2,I)=-HH(I)
   12 RETURN
      END
C      
      SUBROUTINE TP361(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION A(5),B(5),C(5),D(5),H(3),HH(3,5)    
      DATA (A(I),I=1,5)/-0.8720288849D+7,0.1505125253D+6,
     F                  -0.1566950325D+3,0.4764703222D+6,
     F                  0.7294828271D+6/
     F     (B(I),I=1,5)/-0.145421402D+6,0.29311506D+4,-0.40427932D+2,
     F                 0.5106192D+4,0.1571136D+5/
     F     (C(I),I=1,5)/-0.1550111084D+6,0.436053352D+4,0.129492344D+2
     F                  ,0.10236884D+5,0.13176786D+5/
     F     (D(I),I=1,5)/-0.3266695104D+6,0.739068412D+4,
     F                  -0.278986976D+2,0.16643076D+5,0.30988146D+5/
      GOTO (1,2,3,4,5),MODE
    1 N=5
      NILI=0
      NINL=6
      NELI=0
      NENL=0
      X(1)=0.252D+1
      X(2)=0.2D+1
      X(3)=0.375D+2
      X(4)=0.925D+1
      X(5)=0.68D+1
      LXL(5)=.FALSE.
      LXU(1)=.FALSE.
      DO 6 I=1,4
      LXL(I)=.TRUE.
    6 LXU(I+1)=.TRUE.
      XL(1)=0.0D+0
      XL(2)=0.12D+1
      XL(3)=0.2D+2
      XL(4)=0.9D+1
      XU(2)=0.24D+1
      XU(3)=0.6D+2
      XU(4)=0.93D+1
      XU(5)=0.7D+1
      LEX=.FALSE.
      NEX=1
      FEX=-0.77641212D+6
      XEX(1)=0.68128605D+0
      XEX(2)=0.24D+1
      XEX(3)=0.2D+2
      XEX(4)=0.93D+1
      XEX(5)=0.7D+1
      RETURN
    2 FX=A(1)
      DO 7 I=2,5
    7 FX=FX+A(I)*X(I)
      FX=X(1)*FX-0.24345D+5
      FX=-FX
      RETURN
    3 GF(1)=A(1)
      DO 8 I=2,5
      GF(1)=GF(1)+A(I)*X(I)
    8 GF(I)=A(I)*X(1)
      DO 20 I=1,5
   20 GF(I)=-GF(I)
      RETURN
    4 H(1)=B(1)
      H(2)=C(1)
      H(3)=D(1)
      DO 9 I=2,5
      H(1)=H(1)+B(I)*X(I)
      H(2)=H(2)+C(I)*X(I)
    9 H(3)=H(3)+D(I)*X(I)
      DO 10 I=1,3
   10 H(I)=X(1)*H(I)
      DO 11 I=1,3
   11 IF (INDEX1(I)) G(I)=H(I)
      IF (INDEX1(4)) G(4)=0.294D+5-H(1)
      IF (INDEX1(5)) G(5)=0.294D+5-H(2)
      IF (INDEX1(6)) G(6)=0.2772D+6-H(3)
      RETURN
    5 HH(1,1)=B(1)
      HH(2,1)=C(1)
      HH(3,1)=D(1)
      DO 12 I=2,5
      HH(1,1)=HH(1,1)+B(I)*X(I)
      HH(2,1)=HH(2,1)+C(I)*X(I)
      HH(3,1)=HH(3,1)+D(I)*X(I)
      HH(1,I)=B(I)*X(1)
      HH(2,I)=C(I)*X(1)
   12 HH(3,I)=D(I)*X(1)
      DO 13 J=1,3
      IF (.NOT.INDEX2(J)) GOTO 13
      DO 130 I=1,5
  130 GG(J,I)=HH(J,I)
   13 CONTINUE
      DO 14 J=4,6
      IF (.NOT.INDEX2(J)) GOTO 14
      DO 140 I=1,5
  140 GG(J,I)=-HH(J-3,I)
   14 CONTINUE
      RETURN
      END
C
      SUBROUTINE TP362(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION TP362A
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=5
      NILI=4
      NINL=0
      NELI=0
      NENL=0
      X(1)=15.0D0
      X(2)=9.05D0
      X(3)=6.14D0
      X(4)=4.55D0
      X(5)=3.61D0
      DO 14 I=1,5
      XL(I)=3.0D0
      XU(I)=20.0D0
      LXL(I)=.TRUE.
   14 LXU(I)=.TRUE.
      XL(1)=15.D+0
      XL(5)=2.0D0
      LEX=.FALSE.
      NEX=3
      FEX=0.418D-01
      XEX(1)=0.19960505D+02 
      XEX(2)=0.19960505D+02  
      XEX(3)=0.19960505D+02
      XEX(4)=0.19960505D+02
      XEX(5)=0.19960505D+02
C      FEX=0.26229998D0
c      XEX(6)=0.15050962D+2
c      XEX(7)=0.88751199D+1
c      XEX(8)=0.59088230D+1
c      XEX(9)=0.48604810D+1
c      XEX(10)=0.43992690D+1
C      FEX=0.92400000D-01
c      XEX(11)=0.18414858D+02
c      XEX(12)=0.15909146D+02 
c      XEX(13)=0.15909146D+02
c      XEX(14)=0.15909146D+02
c      XEX(15)=0.86690897D+01
      RETURN
    2 FX=TP362A(X)
    3 RETURN
    4 DO 41 I=1,4
   41 IF (INDEX1(I)) G(I)=X(I)-X(I+1)
    5 RETURN 
      END
C      
      DOUBLE PRECISION FUNCTION TP362A(X)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)                     
      DOUBLEPRECISION X(5),RPM,TORQUE,RAD,CON1,CON2,RPMIN,RPMAX,EI,VI,
     /                DT,VMAX,V0,TSHIFT,TMAX,ACC,FORCE,V,ACC0,T,TT
      DATA RAD,CON1,CON2,RPMIN,RPMAX,EI,VI,DT,VMAX,V0,TSHIFT,TMAX/
     F 1.085D+0,1.466667D+0,12.90842D+0,6.D+2,5.7D+3,.6D+0,98.D+0,
     F .1D-1,1.D+2,5.D+0,.25D+0,1.D+2/
      INTEGER I,IT
   13 IT=0
      ACC=0.D+0
      V=V0
      I=1
  302 FORCE=.0239D+0*V**2+31.2D+0
  301 RPM=V*CON2*X(I)
      IF (RPM.LT.RPMIN) GOTO 300
      IF (RPM.GE.RPMAX) GOTO 305
      IF (RPM.GE.6.D+2.AND.RPM.LT.1.9D+3)
     F TORQUE=.3846154D-7*RPM**3-.2108974359D-3*RPM**2
     F +.42455128205133D+0*RPM-1.8711538461540295D+2
      IF (RPM.GE.1.9D+3.AND.RPM.LT.3.D+3)
     F TORQUE=-.492424D-8*RPM**3+.1867424242D-4*RPM**2
     F +.1229545454547D-1*RPM+64.999999999986D+0
      IF (RPM.GE.3.D+3.AND.RPM.LT.4.5D+3) 
     F TORQUE=-.26667D-9*RPM**3+.3D-5*RPM**2
     F -.1263333333336D-1*RPM+1.5510000000002947D+2
      IF (RPM.GE.4.5D+3.AND.RPM.LT.5.6D+3)
     F TORQUE=-.664141D-8*RPM**3+.8337626263D-4*RPM**2
     F -.34351868688129D+0*RPM+5.973636363847145D+2
      IF (RPM.GE.5.6D+3.AND.RPM.LT.6.D+3)
     F TORQUE=-.2539683D-7*RPM**3+.38158730157D-3*RPM**2
     F -1.9223492062348D+0*RPM+3.38066666645715304D+3
      ACC0=ACC
      ACC=RAD*(X(I)*TORQUE-FORCE*RAD)/(EI*X(I)**2+VI)
      IT=IT+1
      T=DT*IT
      V=V+(ACC0+ACC)/2.D+0*DT/CON1
      IF (T.GT.TMAX) GOTO 311
      IF (V.GE.VMAX) GOTO 311
      GOTO 302
  300 TP362A=TMAX
      RETURN
  305 I=I+1
      IF (I.GT.5) GOTO 311
      IF (T.EQ.0.D+0) GOTO 301
      TT=T+TSHIFT
  306 ACC=-FORCE*RAD**2/VI
      IT=IT+1
      T=DT*IT
      V=V+ACC*DT/CON1
      IF (T.LT.TT) GOTO 307
      GOTO 302
  307 FORCE=.0293D+0*V**2+31.2D+0
      GOTO 306
  311 TP362A=T/100.D+0
      RETURN 
      END
C
      SUBROUTINE TP363(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION TP363A
      GOTO (1,2,3,4,5),MODE
    1 N=5
      NILI=0
      NINL=3
      NELI=0
      NENL=0
      X(1)=-.3359769D+0
      X(2)=-.1432398D+1
      X(3)=0.D+0
      X(4)=4.D+0
      X(5)=9.D+0
      DO 11 I=1,5
      LXL(I)=.TRUE.
      LXU(I)=.TRUE.
   11 XU(I)=1.D+1
      XL(1)=-1.D+1
      XL(2)=-1.D+1
      XL(3)=-1.D+1
      XL(4)=.1D+0
      XL(5)=1.D+0
      LEX=.FALSE.
      NEX=1
      XEX(1)=.19852438D+0
      XEX(2)=-3.01059794D+0
      XEX(3)=-0.0530266138D+0
      XEX(4)=2.83165094D+0
      XEX(5)=10.D+0
      FEX=-5.55840576
      RETURN
    2 FX=TP363A(X)
    3 RETURN
    4 CALL TP363B(INDEX1,X,G)
    5 RETURN
      END
C      
      DOUBLE PRECISION FUNCTION TP363A(X)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      COMMON/B363/XMU,RHO,THICK,W2,DTHICK,KKK
      COMMON/ST363/VALUE2       
      DOUBLEPRECISION X(5),XC(100),THICK(100),DTHICK(100)
      DOUBLEPRECISION XMU,RHO,X6,X11,W2,VALUE2,ALPHA,W,EPSI,YI,XI,
     /     EPS,DR
      INTEGER I,KKK,II
C      DATA XMU,ALPHA,W,EPSI,RHO,YI,XI,EPS,X6,X11
C     1     /0.3,1.000,628.0,0.0001,7.263D-4,0.0,1.0,0.1,1.0,1.0/
      XMU=0.3D0
      ALPHA=1000.0D0
      W=628.0D0
      EPSI=0.0001D0
      RHO=7.263D-4
      YI=0.0D0
      XI=1.0D0
      EPS=0.1D0
      X6=1.0D0
      X11=1.0D0
      II=99
      KKK=98
      W2=X6*1.D+2
      DR=X(5)-XI
      XC(1)=XI
      DO 60 I=1,KKK
60    XC(I+1)=XC(I)+DR/DBLE(KKK)
      CALL TP363C(XC,THICK,DTHICK,X,KKK,X11,XI)
      CALL TP363D(XC,TP363A)
      TP363A=-TP363A/1.D+6
      RETURN
      END
C      
      SUBROUTINE TP363B (INDEX1,C,G)
      DIMENSION Y1(100),Y(100),X(100),RST(100),TST(100)
      DIMENSION THICK(100),DTHICK(100),STOT(100),C(5),G(3)
      COMMON/B363/XMU,RHO,THICK,W2,DTHICK,KKK
      COMMON/TFN1363/TMAX
      COMMON/ST363/VALUE2
      LOGICAL INDEX1(3)
      INTEGER I,NN,KKK,II
      DOUBLEPRECISION XMU,ALPHA,W,EPSI,RHO,YI,XI,EPS,Y1,Y,X,RST,TST
      DOUBLEPRECISION THICK,DTHICK,STOT,C,G,DR,TMAX,VALUE2,FXL,XL,Y1I
      DOUBLEPRECISION XR,FXR,ROOT,SMAX,C6,W2,V
C      DATA XMU,ALPHA,W,EPSI,RHO,YI,XI,EPS,C6
C     1     /0.3,1000.0,628.0,0.0001,7.263D-4,0.0,1.0,0.1,1.0/
      XMU=0.3D0
      ALPHA=1000.0D0
      W=628.0D0
      EPSI=0.0001D0
      RHO=7.263D-4
      YI=0.0D0
      XI=1.0D0
      EPS=0.1D0
      C6=1.0D0
      II=99
      KKK=98
      W2=C6*100.0D0
      DR=C(5)-XI
      X(1)=XI
      DO 60 I=1,KKK
   60 X(I+1)=X(I)+DR/DBLE(KKK)
      CALL TP363C(X,THICK,DTHICK,C,KKK,C(4),XI)
    5 Y1I=1.D+5
    2 CALL TP363E(Y1I,YI,XI,C(5),II,Y,Y1,X)
      FXL=-Y(II)
      XL=Y1I
      YI=0.D+0
      Y1I=2503000.D+0
      CALL TP363E(Y1I,YI,XI,C(5),II,Y,Y1,X)
      FXR=-Y(II)
      XR=Y1I
      CALL TP363F(XL,XR,FXL,FXR,EPS,Y,Y1,X,XI,C(5),YI,ROOT,II)
      SMAX=0.D+0
      VALUE2=0.D+0
      DO 300 NN=1,II
      RST(NN)=Y(NN)/(THICK(NN)*X(NN))
      TST(NN)=(Y1(NN)+THICK(NN)*RHO*W2**2*X(NN)**2)/THICK(NN)
      STOT(NN)=DSQRT((RST(NN)-TST(NN))**2+RST(NN)**2+TST(NN)**2)
      IF (STOT(NN).GT.SMAX) SMAX=STOT(NN)
      VALUE2=VALUE2+(3.D+4-STOT(NN))**2
  300 CONTINUE
      VALUE2=DSQRT(VALUE2)
      IF(INDEX1(1)) G(1)=(3.D+4-SMAX)/1.D+3
      IF(INDEX1(2)) G(2)=5.D+0-TMAX
      CALL TP363G(V,THICK,KKK,X)
      IF(INDEX1(3)) G(3)=(625.D+0-V)/10.D+0
  900 RETURN
      END
      SUBROUTINE TP363C(X,THICK,DTHICK,C,KKK,C0,XI)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      DIMENSION X(100),THICK(100),DTHICK(100),C(5)
      COMMON/TFN1363/TMAX
      DOUBLEPRECISION TMAX,XL,X,THICK,DTHICK,C,C0,XI
      INTEGER I,NFST,LM,KKK,JKL
      THICK(1)=C0
      XL=X(KKK+1)-X(1)
      NFST=2
      TMAX=THICK(1)
      DO 10 I=1,KKK
      THICK(I+1)=C0+C(1)*(X(I+1)-XI)
      DO 9 LM=1,NFST
      JKL=LM+1
    9 THICK(I+1)=THICK(I+1)+C(JKL)*DSIN((2.0*DBLE(JKL)-3.0)
     1 *3.1415926535897932D+0*(X(I)-X(1))/XL)
      IF (THICK(I+1).GT.TMAX) TMAX=THICK(I+1)
   10 DTHICK(I)=(THICK(I+1)-THICK(I))/(X(I+1)-X(I))
      RETURN
      END
C      
      SUBROUTINE TP363D(X,XKE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      DIMENSION X(100),THICK(100),DTHICK(100)
      COMMON/B363/XMU,RHO,THICK,W,DTHICK,KKK
      DOUBLEPRECISION X,XKE,THICK,DTHICK,XMU,RHO,W,CONST
      INTEGER I,KKK
      CONST=3.1415926535897932D+0*RHO*W**2
      XKE=0.D+0
      DO 10 I=1,KKK
   10 XKE=XKE+X(I)**3*THICK(I)*(X(I+1)-X(I))
      XKE=XKE*CONST
      RETURN
      END
C      
      SUBROUTINE TP363E(Y1I,YI,XI,XF,II,Y,Y1,X)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      DIMENSION Y1(100),Y(100),X(100)
      DOUBLEPRECISION M0,M1,M2,M3,Y1I,YI,XI,XF,Y,Y1,X,H,XR,YR,Y1R,TP363H
      INTEGER J,II,KK,LL
      X(1)=XI  
      Y(1)=YI
      Y1(1)=Y1I
      H=(XF-XI)/DBLE(II-1)
      KK=II-1
      DO 10 J=1,KK
      LL=J
      XR=X(J)
      YR=Y(J)
      Y1R=Y1(J)
      M0=H*TP363H(XR,YR,Y1R,LL)
      XR=X(J)+H/2.D+0
      YR=Y(J)+H*Y1(J)/2.D+0
      Y1R=Y1(J)+M0/2.D+0
      M1=H*TP363H(XR,YR,Y1R,LL)
      YR=YR+H*M0/4.D+0
      Y1R=Y1(J)+M1/2.D+0
      M2=H*TP363H(XR,YR,Y1R,LL)
      XR=X(J)+H
      YR=Y(J)+H*Y1(J)+H*M1/2.D+0
      Y1R=Y1(J)+M2
      M3=H*TP363H(XR,YR,Y1R,LL)
      Y(J+1)=Y(J)+H*Y1(J)+H/6.D+0*(M0+M1+M2)
      Y1(J+1)=Y1(J)+(M0+2.D+0*M1+2.D+0*M2+M3)/6.D+0
   10 X(J+1)=X(J)+H
      RETURN
      END
C      
      SUBROUTINE TP363F(XL,XR,FXL,FXR,EPS,Y,Y1,X,XI,XF,YI,ROOT,II)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      DIMENSION Y1(100),Y(100),X(100)
      DOUBLEPRECISION XL,XR,FXL,FXR,EPS,Y,Y1,X,XI,XF,YI,ROOT,XAPP
      DOUBLEPRECISION XSAVE,FXAPP,VALUE
      INTEGER II
      XSAVE=XL
  105 XAPP=XL+(FXL*(XR-XL)/(FXL-FXR))
      CALL TP363E(XAPP,YI,XI,XF,II,Y,Y1,X)
      FXAPP=-Y(II)
      IF (DABS(XAPP-XSAVE)/XAPP.LE.EPS) GO TO 250
      VALUE=FXAPP*FXL
      IF (VALUE.LT.0) GO TO 110
      XL=XAPP
      XSAVE=XL
      FXL=FXAPP
      GO TO 105
  110 XR=XAPP
      XSAVE=XR
      FXR=FXAPP
      GO TO 105
  250 ROOT=XAPP
      RETURN
      END
C      
      SUBROUTINE TP363G(V,THICK,KKK,X)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      DIMENSION X(100),THICK(100)
      DOUBLEPRECISION V,THICK,X,PI,DELTX,R1,R2,R3
      INTEGER I,KKK,LMN
      V=0.D+0
      PI=3.141592654D+0
      DELTX=(X(KKK+1)-X(1))/KKK
      LMN=KKK-1
      DO 10 I=1,LMN,2
      R1=(X(I+1)+X(I))/2.D+0
      R2=(X(I+1)+X(I+2))/2.D+0
      R3=(R1+R2)/2.D+0
   10 V=V+2.D+0*PI*DELTX/3.D+0*(THICK(I)*R1+4.D+0*THICK
     1 (I+1)*R3+THICK(I+2)*R2)
      RETURN
      END
C      
      DOUBLE PRECISION FUNCTION TP363H(XR,YR,Y1R,I)
      DIMENSION THICK(100),DTHICK(100)
      INTEGER I
      COMMON/B363/XMU,RHO,THICK,W2,DTHICK,KKK
      INTEGER KKK
      DOUBLEPRECISION XR,YR,Y1R,THICK,DTHICK,XMU,RHO,W2
      TP363H=(1.D+0/THICK(I)*DTHICK(I)-1.D+0/XR)*Y1R+(1.D+0/
     1 (XR**2)-XMU/(XR*THICK(I))*DTHICK(I))*YR-
     2 (3.D+0+XMU)*RHO*W2**2*THICK(I)*XR
      RETURN
      END
C      
      SUBROUTINE TP364(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION XMU1,XMU2,TP364A
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=6
      NILI=2
      NINL=2
      NELI=0
      NENL=0      
      X(1)=1.D+0
      X(2)=4.5D+0
      X(3)=4.D+0
      X(4)=5.D+0
      X(5)=3.D+0
      X(6)=3.D+0
      DO 11 I=1,4
   11 LXL(I)=.TRUE.
      LXL(5)=.FALSE.
      LXL(6)=.FALSE.
      LXU(5)=.FALSE.
      LXU(6)=.FALSE.
      LXU(1)=.TRUE.
      LXU(4)=.TRUE.
      LXU(2)=.FALSE.
      LXU(3)=.FALSE.
      XL(1)=.5D+0
      XL(2)=0.01D+0
      XL(3)=0.D+0
      XL(4)=2.D+0
      XU(1)=3.D+0
      XU(4)=10.D+0
      LEX=.FALSE.
      NEX=1
      FEX=0.0606002
      XEX(1)=0.99616882
      XEX(2)=0.41960616D+01
      XEX(3)=0.29771652D+01
      XEX(4)=0.39631949D+01
      XEX(5)=0.16536702D+01
      XEX(6)=0.12543998D+01
      RETURN
    2 FX=TP364A(X)
    3 RETURN
    4 XMU1=.7853981633D+0
      XMU2=2.356194491D+0
      IF (INDEX1(1)) G(1)=-X(1)+X(2)+X(3)-X(4)
      IF (INDEX1(2)) G(2)=-X(1)-X(2)+X(3)+X(4)
      IF (INDEX1(3)) G(3)=-X(2)*X(2)-X(3)*X(3)+(X(4)-X(1))*
     F (X(4)-X(1))+2.D+0*X(2)*X(3)*DCOS(XMU1)
      IF (INDEX1(4)) G(4)=X(2)*X(2)+X(3)*X(3)-(X(4)+X(1))*
     F (X(4)+X(1))-2.D+0*X(2)*X(3)*DCOS(XMU2)
    5 RETURN
      END
C      
      DOUBLE PRECISION FUNCTION TP364A(X)
      DOUBLEPRECISION X(6),X1A(31),Y1A(31),PHI(31),X1(31),Y1(31),WP,
     F  COSS,COSY,SINS,SINY,XINC,PI
      INTEGER I
      PI=3.141592654D+0
      XINC=2.D+0*PI/30.D+0
      DO 1 I=1,31
    1 PHI(I)=XINC*DBLE(I-1)
      CALL TP364B(PHI,X1,Y1)
      TP364A=0.D+0
      DO 2 I=1,31
      CALL TP364C(X,PHI(I),COSS)
      WP=DABS(1.D+0-COSS*COSS)
      IF (WP.GT.0.0) THEN
         SINS=DSQRT(WP)
      ELSE
         SINS=0.0
      ENDIF      
      COSY=(X(4)+X(3)*COSS-X(1)*DCOS(PHI(I)))/X(2)
      SINY=(X(3)*SINS-X(1)*DSIN(PHI(I)))/X(2)
      X1A(I)=X(1)*DCOS(PHI(I))+X(5)*COSY-X(6)*SINY
      Y1A(I)=X(1)*DSIN(PHI(I))+X(5)*SINY+X(6)*COSY
    2 TP364A=TP364A+(X1A(I)-X1(I))**2+(Y1A(I)-Y1(I))**2
      WP=TP364A/31.D+0
      IF (WP.GT.0.0) THEN
         TP364A=DSQRT(WP)
      ELSE
         TP364A=0.0
      ENDIF      
      RETURN
      END
C
      SUBROUTINE TP364B(PHI,X1,Y1)
      DOUBLEPRECISION PHI(31),X1(31),Y1(31),PI
      INTEGER I
      PI=3.141592654D+0
      DO 1 I=1,31
      X1(I)=.4D+0+DSIN((2.D+0*PI)*((PI-PHI(I))/(2.D+0*PI)-.16D+0))
    1 Y1(I)=2.D+0+.9D+0*DSIN(PI-PHI(I))
      RETURN
      END
C      
      SUBROUTINE TP364C(X,PHI,W)
      DOUBLEPRECISION X(6),K,L,M,W,PHI,PI,A,B,C,TERM
      PI=3.141592654D+0
      M=2.D+0*X(1)*X(3)*DSIN(PHI)
      L=2.D+0*X(3)*X(4)-2.D+0*X(1)*X(3)*DCOS(PHI)
      K=X(1)*X(1)-X(2)*X(2)+X(3)*X(3)+X(4)*X(4)-
     F 2.D+0*X(4)*X(1)*DCOS(PHI)
      A=L*L+M*M
      B=2.D+0*K*L
      C=K*K-M*M
      TERM=DSQRT(DABS(B*B-4.D+0*A*C))
      IF ((PI-PHI).LT.0.D+0) TERM=-TERM
      W=(-B+TERM)/(2.D+0*A)
      RETURN
      END
C
      SUBROUTINE TP365(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION P,Q
      INTEGER I
      GOTO (1,2,3,4,5), MODE
    1 N=7
      NILI=0
      NINL=5 
      NELI=0
      NENL=0
      X(1)=3.D+0
      X(2)=0.D+0
      X(3)=2.D+0
      X(4)=-1.5D+0
      X(5)=1.5D+0
      X(6)=5.D+0
      X(7)=0.D+0
      DO 11 I=1,3
      LXL(I*2)=.FALSE.
   11 LXL(I*2-1)=.TRUE.
      LXL(7)=.TRUE.
      DO 12 I=1,7
   12 LXU(I)=.FALSE.
      XL(1)=0.D+0
      XL(3)=0.D+0
      XL(5)=1.D+0
      XL(7)=1.D+0
      LEX=.FALSE.
      NEX=1
      FEX=0.23313708D+2
      XEX(1)=0.48284266D+1
      XEX(2)=0.47529555D-5
      XEX(3)=0.48284276D+1
      XEX(4)=0.10000024D+1
      XEX(5)=0.24142144D+1
      XEX(6)=0.24142151D+1
      XEX(7)=0.10000000D+1
      RETURN
    2 FX=X(1)*X(3)
    3 RETURN
    4 P=DSQRT(X(2)**2+X(3)**2)
      Q=DSQRT(X(3)**2+(X(2)-X(1))**2)
      IF (INDEX1(1)) G(1)=(X(4)-X(6))**2+(X(5)-X(7))**2-4.D+0
      IF (INDEX1(2)) G(2)=(X(3)*X(4)-X(2)*X(5))/P-1.D+0
      IF (INDEX1(3)) G(3)=(X(3)*X(6)-X(2)*X(7))/P-1.D+0
      IF (INDEX1(4)) G(4)=(X(1)*X(3)+(X(2)-X(1))*X(5)-X(3)*X(4))
     F    /Q-1.D+0
      IF (INDEX1(5)) G(5)=(X(1)*X(3)+(X(2)-X(1))*X(7)-X(3)*X(6))
     F    /Q-1.D+0
    5 RETURN 
      END
C
      SUBROUTINE TP366(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION C(38)
      DATA C/.59553571D-3,.88392857D+0,-.11756250D+0,1.1088D+0,
     F .1303533D+0,-.0066033D+0,.66173269D-3,.17239878D-1,
     F -.56595559D-2,-.19120592D-1,.5685075D+2,1.08702D+0,
     F .32175D+0,-.03762D+0,.006198D+0,.24623121D+4,-.25125634D+2,
     F .16118996D+3,5.D+3,-.48951D+6,.44333333D+2,.33D+0,.022556D+0,
     F -.007595D+0,.00061D+0,-.5D-3,.819672D+0,.819672D+0,.245D+5,
     F -.25D+3,.10204082D-1,.12244898D-4,.625D-4,.625D-4,-.7625D-4,
     F 1.22D+0,1.D+0,-1.D+0/
      GOTO(1,2,3,4,5),MODE
    1 N=7
      NILI=0
      NINL=14
      NELI=0
      NENL=0
      X(1)=1745.D+0
      X(2)=110.D+0
      X(3)=3048.D+0
      X(4)=89.D+0
      X(5)=92.8D+0
      X(6)=8.D+0
      X(7)=145.D+0
      DO 11 I=1,7
      LXL(I)=.TRUE.
   11 LXU(I)=.TRUE.
      XL(1)=1.D+0
      XL(2)=1.D+0
      XL(3)=1.D+0
      XL(4)=85.D+0
      XL(5)=90.D+0
      XL(6)=3.D+0
      XL(7)=145.D+0
      XU(1)=2.D+3
      XU(2)=1.2D+2
      XU(3)=5.D+3
      XU(4)=93.D+0
      XU(5)=95.D+0
      XU(6)=12.D+0
      XU(7)=162.D+0
      LEX=.FALSE.
      NEX=1
      FEX=0.70430560D+03 
      XEX(1)=0.90540351D+3
      XEX(2)=0.36394998D+2
      XEX(3)=0.23814783D+4
      XEX(4)=0.88987691D+2
      XEX(5)=0.95D+2
      XEX(6)=0.12D+2
      XEX(7)=0.15353535D+3
      RETURN
    2 FX=(1.715D+0*X(1)+.035D+0*X(1)*X(6)+4.0565D+0*X(3)
     1 +10.D+0*X(2)+3000.D+0-.063D+0*X(3)*X(5))
    3 RETURN
    4 IF (INDEX1(1)) G(1)=1.D+0-C(1)*X(6)**2-C(2)*X(3)/X(1)-C(3)*X(6)
      IF (INDEX1(2)) G(2)=1.D+0-C(4)*X(1)/X(3)-C(5)*X(1)/X(3)*X(6)
     F -C(6)*X(1)/X(3)*X(6)**2
      IF (INDEX1(3)) G(3)=1.D+0-C(7)*X(6)**2-C(8)*X(5)-C(9)*X(4)
     F -C(10)*X(6)
      IF (INDEX1(4)) G(4)=1.D+0-C(11)/X(5)-C(12)/X(5)*X(6)-C(13)
     F *X(4)/X(5)-C(14)/X(5)*X(6)**2
      IF (INDEX1(5)) G(5)=1.D+0-C(15)*X(7)-C(16)*X(2)/X(3)/X(4)
     F -C(17)*X(2)/X(3)
      IF (INDEX1(6)) G(6)=1.D+0-C(18)/X(7)-C(19)*X(2)/X(3)/X(7)
     F -C(20)*X(2)/X(3)/X(4)/X(7)
      IF (INDEX1(7)) G(7)=1.D+0-C(21)/X(5)-C(22)*X(7)/X(5)
      IF (INDEX1(8)) G(8)=1.D+0-C(23)*X(5)-C(24)*X(7)
      IF (INDEX1(9)) G(9)=1.D+0-C(25)*X(3)-C(26)*X(1)
      IF (INDEX1(10)) G(10)=1.D+0-C(27)*X(1)/X(3)-C(28)/X(3)
      IF (INDEX1(11)) G(11)=1.D+0-C(29)*X(2)/X(3)/X(4)-C(30)*X(2)/X(3)
      IF (INDEX1(12)) G(12)=1.D+0-C(31)*X(4)-C(32)/X(2)*X(3)*X(4)
      IF (INDEX1(13)) G(13)=1.D+0-C(33)*X(1)*X(6)-C(34)*X(1)-C(35)*X(3)
      IF (INDEX1(14)) G(14)=1.D+0-C(36)/X(1)*X(3)-C(37)/X(1)-C(38)*X(6)
    5 RETURN 
      END
C
      SUBROUTINE TP367(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5), MODE
    1 N=7
      NILI=2
      NINL=1
      NELI=1
      NENL=1
      DO 11 I=1,7
      X(I)=.1D+0
      LXL(I)=.TRUE.
      LXU(I)=.FALSE.
   11 XL(I)=0.D+0
      LEX=.FALSE.
      NEX=1
      FEX=-0.37412960D+2
      XEX(1)=0.14688103D+1
      XEX(2)=0.19839711D+1
      XEX(3)=0.35187754D+0
      XEX(4)=0.11953411D+1
      XEX(5)=0.56940029D+0
      XEX(6)=0.78474478D+0
      XEX(7)=0.14121216D+1
      DO 12 I=1,7
   12 GG(1,I)=-1.D+0
      DO 13 I=1,4
   13 GG(2,I)=-1.D+0
      DO 14 I=5,7
   14 GG(2,I)=0.D+0
      GG(3,1)=-1.D+0
      GG(3,2)=0.D+0
      GG(3,3)=-1.D+0
      GG(3,4)=0.D+0
      GG(3,5)=-1.D+0
      GG(4,1)=0.D+0
      GG(4,2)=0.D+0
      GG(4,3)=0.D+0
      GG(4,4)=2.D+0
      GG(4,5)=1.D+0
      GG(4,6)=.8D+0
      GG(4,7)=1.D+0
      GG(5,1)=0.D+0
      GG(5,4)=0.D+0
      GG(5,7)=0.D+0
      GF(2)=-5.D+0
      GF(4)=-6.D+0
      RETURN
    2 FX=-5.D+0*X(1)-5.D+0*X(2)-4.D+0*X(3)-X(1)*X(3)-6.D+0*X(4)
     F -5.D+0*X(5)/(1.D+0+X(5))-8.D+0*X(6)/(1.D+0+X(6))
     F -10.D+0*(1.D+0-2.D+0*DEXP(-X(7))+DEXP(-2.D+0*X(7)))
      RETURN
    3 GF(1)=-5.D+0-X(3)
      GF(3)=-4.D+0-X(1)
      GF(5)=-5.D+0/(1.D+0+X(5))**2
      GF(6)=-8.D+0/(1.D+0+X(6))**2
      GF(7)=-20.D+0*(DEXP(-X(7))-DEXP(-2.D+0*X(7)))
      RETURN
    4 IF (INDEX1(1)) G(1)=10.D+0-X(1)-X(2)-X(3)-X(4)-X(5)-X(6)-X(7)
      IF (INDEX1(2)) G(2)=5.D+0-X(1)-X(2)-X(3)-X(4)
      IF (INDEX1(3)) G(3)=5.D+0-X(1)-X(3)-X(5)-X(6)**2-X(7)**2
      IF (INDEX1(4)) G(4)=2.D+0*X(4)+X(5)+.8D+0*X(6)+X(7)-5.D+0
      IF (INDEX1(5)) G(5)=X(2)**2+X(3)**2+X(5)**2+X(6)**2-5.D+0
      RETURN
    5 IF (.NOT.INDEX2(3)) GOTO 51
      GG(3,6)=-2.D+0*X(6)
      GG(3,7)=-2.D+0*X(7)
   51 IF (.NOT.INDEX2(5)) GOTO 52
      GG(5,2)=2.D+0*X(2)
      GG(5,3)=2.D+0*X(3)
      GG(5,5)=2.D+0*X(5)
      GG(5,6)=2.D+0*X(6)
   52 RETURN 
      END
C
      SUBROUTINE TP368(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      COMMON     /D368/S2,S3,S4
      DOUBLEPRECISION S2,S3,S4
      GOTO (1,2,2,4,4), MODE
    1 N=8
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 11 I=1,8
      X(I)=1.0-1.0/DBLE(I)
      LXL(I)=.TRUE.
      LXU(I)=.TRUE.
      XL(I)=0.D+0
   11 XU(I)=1.D+0
      X(7)=0.7D0
      X(8)=0.7D0
      NEX=3
      LEX=.TRUE.
      XEX(1)=1.0D+0
      XEX(2)=0.5D+0
      XEX(3)=0.5D+0
      XEX(4)=1.0D+0
      XEX(5)=1.0D+0
      XEX(6)=1.0D+0
      XEX(7)=0.5D+0
      XEX(8)=0.5D+0     
      FEX=-0.74997564D+0
c      XEX(9)=0.49834105D+0
c      XEX(10)=0.49977950D+0
c      XEX(11)=0.50201378D+0
c      XEX(12)=0.50378302D+0
c      XEX(13)=0.50263008D+0
c      XEX(14)=0.50232579D+0
c      XEX(15)=0.10000000D+1
c      XEX(16)=0.10000000D+1        
c      FEX=-0.75
c      XEX(17)=0.5
c      XEX(18)=0.5
c      XEX(19)=0.5
c      XEX(20)=0.5
c      XEX(21)=0.5
c      XEX(22)=0.5
c      XEX(23)=0.10000000D+1
c      XEX(24)=0.10000000D+1
c      FEX=-1.0D0
      RETURN
    2 S2=0.D+0    
      S3=0.D+0
      S4=0.D+0
      DO 10 I=1,8
      S2=S2+X(I)**2
      S3=S3+X(I)**3     
   10 S4=S4+X(I)**4
      IF (MODE .EQ.3) GOTO 3
      FX=-S2*S4+S3**2
      RETURN
    3 DO 31 I=1,8
   31 GF(I)=-2.D+0*X(I)*S4-4.D+0*X(I)**3*S2+6.D+0*X(I)**2*S3
    4 RETURN
      END
C      
      SUBROUTINE TP369(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION C(16)
      DATA C/833.33252D+0,100.D+0,-83333.333D+0,1250.D+0,1.D+0,-1250.D+0
     F ,1250000.D+0,1.D+0,-2500.D+0,2.5D-3,2.5D-3,2.5D-3,2.5D-3,-2.5D-3
     F ,1.D-2,-1.D-2/   
      GOTO (1,2,3,4,5),MODE
    1 N=8 
      NILI=3
      NINL=3 
      NELI=0 
      NENL=0
      X(1)=5.D+3
      X(2)=5.D+3
      X(3)=5.D+3
      X(4)=2.D+2
      X(5)=3.5D+2
      X(6)=1.5D+2
      X(7)=225.D+0
      X(8)=425.D+0        
      DO 11 I=1,8
      LXU(I)=.TRUE.
   11 LXL(I)=.TRUE.
      XL(1)=1.D+2
      XL(2)=1.D+3
      XL(3)=1.D+3
      DO 12 I=1,3
   12 XU(I)=1.D+4
      DO 13 I=4,8
      XL(I)=10.D+0
   13 XU(I)=1.D+3
      LEX=.FALSE.
      NEX=1
      FEX=0.70492480D+4
      XEX(1)=0.57930657D+3
      XEX(2)=0.13599705D+4
      XEX(3)=0.51099709D+4
      XEX(4)=0.18201769D+3
      XEX(5)=0.29560116D+3
      XEX(6)=0.21798231D+3
      XEX(7)=0.28641653D+3
      XEX(8)=0.39560116D+3
      RETURN
    2 FX=X(1)+X(2)+X(3)
    3 RETURN
    4 CONTINUE
      DO I=1,5
      IF (X(I).GT.XU(I)) X(I)=XU(I)
      IF (X(I).LT.XL(I)) X(I)=XL(I)
      ENDDO 
      IF (INDEX1(1)) G(1)=1.0-C(10)*X(4)-C(11)*X(6)
      IF (INDEX1(2)) G(2)=1.0-C(12)*X(5)-C(13)*X(7)-C(14)*X(4)
      IF (INDEX1(3)) G(3)=1.0-C(15)*X(8)-C(16)*X(5)
      IF (INDEX1(4)) G(4)=1.0-C(1)/X(1)*X(4)/X(6)-C(2)/X(6)
     F -C(3)/X(1)/X(6)
      IF (INDEX1(5)) G(5)=1.0-C(4)/X(2)*X(5)/X(7)-C(5)*X(4)/X(7)
     F -C(6)/X(2)*X(4)/X(7)
      IF (INDEX1(6)) G(6)=1.0-C(7)/X(3)/X(8)-C(8)*X(5)/X(8)            
     F -C(9)/X(3)*X(5)/X(8)
    5 RETURN
      END
C
      SUBROUTINE TP370(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J,LSUM
      DOUBLEPRECISION F,DF,SUM,SUM1,SUM2(31),BASIS
      N=6
      FEX=.228767005355D-2
      XEX(1)=-0.15663881D-01
      XEX(2)=0.10124222D+01
      XEX(3)=-0.23290211D+00
      XEX(4)=0.12604426D+01
      XEX(5)=-0.15138534D+01
      XEX(6)=0.99314584D+00
      GOTO 10
      ENTRY TP371(MODE)
      N=9
C      FEX=0.1399766D-5
      FEX=0.13997601D-5
      XEX(1)=-0.10630204D-03
      XEX(2)=0.99902770D+00
      XEX(3)=0.30998203D-01 
      XEX(4)=0.61696369D-01
      XEX(5)=0.11466839D+01
      XEX(6)=-0.25891812D+01
      XEX(7)=0.37452496D+01
      XEX(8)=-0.27632595D+01
      XEX(9)=0.92600597D+00
   10 GOTO (1,2,4,4,4),MODE
    1 NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 11 I=1,N
      X(I)=0.D+0
      LXL(I)=.FALSE.
   11 LXU(I)=.FALSE.
      LSUM=31
      LEX=.FALSE.
      NEX=1
      RETURN
    2 F(1)=X(1)
      F(2)=X(2)-X(1)**2-1.D+0
      DO 20 I=2,30
      BASIS=DBLE(I-1)/29.D+0
      SUM1=0.D+0
      DO 21 J=2,N
   21 SUM1=SUM1+X(J)*DBLE(J-1)*BASIS**(J-2)
      SUM=0.D+0
      DO 23 J=1,N
   23 SUM=SUM+X(J)*BASIS**(J-1)
      SUM2(I+1)=SUM
   20 F(I+1)=SUM1-SUM2(I+1)**2-1.D+0
      FX=0.D+0
      DO 22 I=1,31
   22 FX=FX+F(I)**2
    4 RETURN
      END
C
      SUBROUTINE TP372(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF
      INTEGER I,J,LSUM
      GOTO (1,2,3,4,5),MODE
    1 N=9 
      NILI=0
      NINL=12
      NELI=0
      NENL=0
      X(1)=3.D+2
      X(2)=-1.D+2
      X(3)=-.1997D+0
      X(4)=+127.D+0
      X(5)=+151.D+0
      X(6)=379.D+0
      X(7)=421.D+0
      X(8)=460.D+0
      X(9)=426.D+0
      DO 12 I=1,12
      DO 12 J=4,9
   12 GG(I,J)=0.D+0
      DO 18 I=1,6
      GG(I,1)=1.D+0
   18 GG(I,I+3)=1.D+0
      DO 13 I=7,12
      GG(I,1)=-1.D+0
   13 GG(I,I-3)=1.D+0
      GF(1)=0.D+0
      GF(2)=0.D+0
      GF(3)=0.D+0
      DO 14 I=4,9
      LXL(I)=.TRUE.
      LXU(I)=.FALSE.
   14 XL(I)=0.D+0
      DO 15 I=1,2
      LXL(I)=.FALSE.
   15 LXU(I)=.FALSE.
      XL(3)=-1.0
      XU(3)=0.0
      LXL(3)=.TRUE.
      LXU(3)=.TRUE.
      XEX(1)=0.52330555D+3
      XEX(2)=-0.15694787D+3
      XEX(3)=-0.19966457D+0
      XEX(4)=0.29608067D+2
      XEX(5)=0.86615521D+2
      XEX(6)=0.47326718D+2
      XEX(7)=0.26235604D+2
      XEX(8)=0.22915985D+2
      XEX(9)=0.39470742D+2
      LEX=.FALSE.
      NEX=1
      FEX=0.13390093D+5
      LSUM=6
      RETURN
    2 FX=0.D+0
      DO 21 I=4,9
      F(I-3)=X(I)
   21 FX=FX+X(I)**2 
      RETURN
    3 DO 35 I=4,9
   35 F(I-3)=X(I)
      DO 31 I=4,9
   31 GF(I)=2.D+0*X(I)
      DO 33 I=1,6
      DO 33 J=1,9
   33 DF(I,J)=0.D+0
      DO 34 I=1,6
   34 DF(I,I+3)=1.D+0
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)+X(2)*DEXP(-5.D+0*X(3))+X(4)-127.D+0
      IF (INDEX1(2)) G(2)=X(1)+X(2)*DEXP(-3.D+0*X(3))+X(5)-151.D+0
      IF (INDEX1(3)) G(3)=X(1)+X(2)*DEXP(-X(3))+X(6)-379.D+0
      IF (INDEX1(4)) G(4)=X(1)+X(2)*DEXP(X(3))+X(7)-421.D+0
      IF (INDEX1(5)) G(5)=X(1)+X(2)*DEXP(3.D+0*X(3))+X(8)-460.D+0
      IF (INDEX1(6)) G(6)=X(1)+X(2)*DEXP(5.D+0*X(3))+X(9)-426.D+0
      IF (INDEX1(7)) G(7)=-X(1)-X(2)*DEXP(-5.D+0*X(3))+X(4)+127.D+0
      IF (INDEX1(8)) G(8)=-X(1)-X(2)*DEXP(-3.D+0*X(3))+X(5)+151.D+0
      IF (INDEX1(9)) G(9)=-X(1)-X(2)*DEXP(-X(3))+X(6)+379.D+0
      IF (INDEX1(10)) G(10)=-X(1)-X(2)*DEXP(X(3))+X(7)+421.D+0
      IF (INDEX1(11)) G(11)=-X(1)-X(2)*DEXP(3.D+0*X(3))+X(8)+460.D+0
      IF (INDEX1(12)) G(12)=-X(1)-X(2)*DEXP(5.D+0*X(3))+X(9)+426.D+0
      RETURN
    5 DO 51 I=1,6
      IF (.NOT.INDEX2(I+6)) GOTO 51
      GG(I+6,2)=-DEXP(DBLE(I*2-7)*X(3))
      GG(I+6,3)=-X(2)*DBLE(I*2-7)*GG(I+6,2)
   51 CONTINUE
      DO 52 I=1,6
      IF (.NOT.INDEX2(I)) GOTO 52
      GG(I,2)=DEXP(DBLE(I*2-7)*X(3))
      GG(I,3)=X(2)*DBLE(I*2-7)*GG(I,2)
   52 CONTINUE
      RETURN 
      END
C
      SUBROUTINE TP373(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF
      INTEGER I,J,LSUM
      GOTO (1,2,3,4,5), MODE
    1 NILI=0
      NINL=0
      NELI=0
      NENL=6
      N=9 
      DO 11 I=1,9
      LXL(I)=.FALSE.
   11 LXU(I)=.FALSE.
      LXL(3)=.TRUE.
      LXU(3)=.TRUE.
      XL(3)=-1.0
      XU(3)=0.0
      X(1)=3.D+2
      X(2)=-1.D+2
      X(3)=-.1997D+0
      X(4)=-127.D+0
      X(5)=-151.D+0
      X(6)=379.D+0
      X(7)=421.D+0
      X(8)=460.D+0
      X(9)=426.D+0
      GF(1)=0.D+0
      GF(2)=0.D+0
      GF(3)=0.D+0
      DO 17 I=1,6
      DO 17 J=4,9
   17 GG(I,J)=0.D+0
      DO 18 I=1,6
      GG(I,1)=1.D+0
   18 GG(I,I+3)=1.D+0
      LEX=.TRUE.
      NEX=1
      FEX=0.13390093D+5
      XEX(1)=0.52330542D+3
      XEX(2)=-0.15694770D+3
      XEX(3)=-0.19966472D+0
      XEX(4)=0.29608061D+2 
      XEX(5)=-0.86615571D+2
      XEX(6)=0.47326669D+2
      XEX(7)=0.26235575D+2
      XEX(8)=0.22915982D+2
      XEX(9)=-0.39470718D+2
      LSUM=6
      RETURN
    2 FX=0.D+0
      DO 21 I=4,9
      F(I-3)=X(I)
   21 FX=FX+X(I)**2 
      RETURN
    3 DO 35 I=4,9
   35 F(I-3)=X(I)
      DO 31 I=4,9
   31 GF(I)=2.D+0*X(I)
      DO 33 I=1,6
      DO 33 J=1,9
   33 DF(I,J)=0.D+0
      DO 34 I=1,6
   34 DF(I,I+3)=1.D+0
      RETURN
    4 IF (INDEX1(1)) G(1)=X(1)+X(2)*DEXP(-5.D+0*X(3))+X(4)-127.D+0
      IF (INDEX1(2)) G(2)=X(1)+X(2)*DEXP(-3.D+0*X(3))+X(5)-151.D+0
      IF (INDEX1(3)) G(3)=X(1)+X(2)*DEXP(-X(3))+X(6)-379.D+0
      IF (INDEX1(4)) G(4)=X(1)+X(2)*DEXP(X(3))+X(7)-421.D+0
      IF (INDEX1(5)) G(5)=X(1)+X(2)*DEXP(3.D+0*X(3))+X(8)-460.D+0
      IF (INDEX1(6)) G(6)=X(1)+X(2)*DEXP(5.D+0*X(3))+X(9)-426.D+0
      RETURN
    5 DO 52 I=1,6
      IF (.NOT.INDEX2(I)) GOTO 52
      GG(I,2)=DEXP(DBLE(I*2-7)*X(3))
      GG(I,3)=X(2)*DBLE(I*2-7)*GG(I,2)
   52 CONTINUE
      RETURN 
      END
C
      SUBROUTINE TP374(MODE)	
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,K
      DOUBLEPRECISION Z,PI,TP374A,TP374B,TP374G
      GOTO(1,2,3,4,5)MODE
    1 N=10
      NILI=0
      NINL=35
      NELI=0
      NENL=0
      DO 6 I=1,10
      X(I)=0.1D+0
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.FALSE.
      NEX=2
      FEX=0.233264D+0
      XEX(1)=0.218212D+0
      XEX(2)=0.232640D+0
      XEX(3)=0.278457D+0
      XEX(4)=0.268125D+0
      XEX(5)=0.212010D+0
      XEX(6)=0.125918D+0
      XEX(7)=0.34102D-1
      XEX(8)=-0.26136D-1
      XEX(9)=-0.142233D+0
      XEX(10)=0.233264D+0
      XEX(11)=-0.142233D+0
      XEX(12)=-0.26136D-1
      XEX(13)=0.34102D-1
      XEX(14)=0.125918D+0
      XEX(15)=0.212010D+0
      XEX(16)=0.268125D+0
      XEX(17)=0.278457D+0
      XEX(18)=0.23264D+0
      XEX(19)=0.218212D+0
      XEX(20)=0.233264D+0
      DO 46 I=1,9
   46 GF(I)=0.0D+0
      GF(10)=0.1D+1
      RETURN
    2 FX=X(10)
    3 RETURN
    4 PI=0.4D+1*DATAN(0.1D+1)
      DO 8 I=1,10
      Z=PI/4.D+0*(DBLE(I-1)*0.1D+0)
    8 IF (INDEX1(I)) G(I)=TP374G(Z,X)-(1.D+0-X(10))**2
      DO 9 I=11,20
      Z=PI/4.D+0*(DBLE(I-11)*0.1D+0)
    9 IF(INDEX1(I)) G(I)=(1.D+0+X(10))**2-TP374G(Z,X)
      DO 10 I=21,35
      Z=PI/4.D+0*(1.2D+0+DBLE(I-21)*0.2D+0)
   10 IF(INDEX1(I)) G(I)=X(10)**2-TP374G(Z,X)
      RETURN
    5 PI=0.4D+1*DATAN(0.1D+1)
      DO 50 I=1,10
      IF (.NOT.INDEX2(I)) GOTO 50
      Z=PI/4.D+0*(DBLE(I-1)*0.1D+0)
      DO 51 K=1,9
      GG(I,K)=2.D+0*(TP374A(Z,X)*DCOS(K*Z)+TP374B(Z,X)*DSIN(K*Z))
   51 GG(I,10)=2.D+0*(1.D+0-X(10))
   50 CONTINUE
      DO 52 I=11,20
      IF (.NOT.INDEX2(I)) GOTO 52
      Z=PI/4.D+0*(DBLE(I-11)*0.1D+0)
      DO 53 K=1,9
      GG(I,K)=-2.D+0*(TP374A(Z,X)*DCOS(K*Z)+TP374B(Z,X)*DSIN(K*Z))
   53 GG(I,10)=2.D+0*(1.D+0+X(10))
   52 CONTINUE
      DO 54 I=21,35
      IF (.NOT.INDEX2(I)) GOTO 54
      Z=PI/4.D+0*(1.2D+0+DBLE(I-21)*0.2D+0)
      DO 55 K=1,9
      GG(I,K)=-2.D+0*(TP374A(Z,X)*DCOS(K*Z)+TP374B(Z,X)*DSIN(K*Z))
   55 GG(I,10)=2.D+0*X(10)
   54 CONTINUE
      RETURN
      END
C      
      DOUBLE PRECISION FUNCTION TP374G(A,X)
      DOUBLEPRECISION A,X(10),TP374A,TP374B
      TP374G=TP374A(A,X)**2+TP374B(A,X)**2
      RETURN
      END                       			
C      
      DOUBLE PRECISION FUNCTION TP374A(A,X)
      DOUBLEPRECISION X(10),A
      INTEGER K
      TP374A=0.D+0
      DO 10 K=1,9
   10 TP374A=TP374A+(X(K)*DCOS(K*A))
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION TP374B(A,X)
      DOUBLEPRECISION X(10),A
      INTEGER K
      TP374B=0.D+0
      DO 10 K=1,9
   10 TP374B=TP374B+(X(K)*DSIN(K*A))
      RETURN
      END
C
      SUBROUTINE TP375(MODE)	
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF,TP375A
      INTEGER I,J,LSUM
      GOTO(1,2,2,4,5),MODE
    1 N=10
      NILI=0
      NINL=0
      NELI=8
      NENL=1
      DO 6 I=1,10
      X(I)=0.1D+1
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=-0.15161D+2
      DO 60 I=1,8
   60 XEX(I)=0.1064812D+0
      XEX(9)=0.28438742D+1
      XEX(10)=-0.26424832D+1
      DO 7 J=1,8
      DO 7 I=1,10
    7 GG(J,I)=0.1D+1
      DO 8 I=1,8
    8 GG(I,I)=0.5D+0
      LSUM=10
      DO 17 I=1,10
      DO 17 J=1,10
      DF(I,J)=0.D+0
   17 DF(I,I)=-0.1D+1
      RETURN
    2 DO 16 I=1,10
   16 F(I)=-X(I)
      IF (MODE.EQ.3) GOTO 3
      FX=0.0D+0
      DO 9 I=1,10
    9 FX=FX-X(I)**2
      RETURN
    3 DO 10 I=1,10
   10 GF(I)=-0.2D+1*X(I)
    4 DO 11 J=1,8
      IF(.NOT.INDEX1(J)) GOTO 11
      G(J)=0.0D+0
      DO 12 I=1,10
   12 G(J)=G(J)+X(I)/TP375A(I,J)
      G(J)=G(J)-0.1D+1
   11 CONTINUE
      IF(.NOT.INDEX1(9)) GOTO 18
      G(9)=0.0D+0
      DO 13 I=1,10
   13 G(9)=G(9)+X(I)**2/(0.1D+1+DBLE(I-1)/0.3D+1)
      G(9)=G(9)-0.4D+1
   18 RETURN
    5 IF(.NOT.INDEX2(9)) GOTO 15
      DO 14 I=1,10
   14 GG(9,I)=0.2D+1*X(I)/(0.1D+1+DBLE(I-1)/0.3D+1)
   15 RETURN
      END
C      
      DOUBLE PRECISION FUNCTION TP375A(I,J)
      INTEGER I,J
      TP375A=0.1D+1
      IF(I.EQ.J) TP375A=0.2D+1
      RETURN
      END
C
      SUBROUTINE TP376(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      DOUBLEPRECISION F,DF
      INTEGER I
      GOTO(1,2,3,4,5),MODE
    1 N=10
      NILI=0
      NINL=14
      NELI=1
      NENL=0
      X(1)=1.0D-0
      X(2)=0.5D-2
      X(3)=0.81D-2
      X(4)=1.D+2
      X(5)=0.17D-2
      X(6)=0.13D-2
      X(7)=0.27D-2
      X(8)=0.2D-2
      X(9)=0.15
      X(10)=0.105
      DO 6 I=1,10
      LXU(I)=.TRUE.
    6 LXL(I)=.TRUE.
      XL(1)=0.D+0
      XL(2)=0.D+0	
      XL(3)=0.5D-4
      XL(4)=0.1D+2
      DO 7 I=5,10
    7 XL(I)=0.1D-2
      XU(1)=0.1D+2
      XU(2)=0.1D+0
      XU(3)=0.81D-2
      XU(4)=0.1D+4
      XU(5)=0.17D-2
      XU(6)=0.13D-2
      XU(7)=0.27D-2
      XU(8)=0.2D-2
      XU(9)=0.1D+1
      XU(10)=0.1D+1
      GG(1,1)=0.1D+1
      GG(1,2)=0.D+0
      DO 9 I=5,10
    9 GG(1,I)=0.D+0
      GG(2,1)=0.1D+1
      GG(2,2)=0.D+0
      GG(2,3)=0.D+0
      DO 10 I=6,8
   10 GG(2,I)=0.D+0
      GG(2,10)=0.D+0
      GG(3,1)=0.1D+1
      GG(3,2)=0.D+0
      GG(3,3)=0.D+0
      GG(3,5)=0.D+0
      DO 11 I=7,9
   11 GG(3,I)=0.D+0
      GG(4,1)=0.1D+1
      GG(4,2)=0.D+0
      GG(4,3)=0.D+0
      GG(4,5)=0.D+0
      GG(4,6)=0.D+0
      DO 12 I=8,10
   12 GG(4,I)=0.D+0
      GG(5,1)=0.1D+1
      GG(5,2)=0.D+0
      GG(5,3)=0.D+0
      DO 13 I=5,7
   13 GG(5,I)=0.D+0
      GG(5,9)=0.D+0
      GG(5,10)=0.D+0
      GG(6,1)=0.D+0	
      GG(6,2)=0.1D+5
      GG(6,3)=0.D+0
      DO 14 I=6,8
   14 GG(6,I)=0.D+0
      GG(6,10)=0.D+0
      GG(7,1)=0.D+0
      GG(7,2)=0.1D+5
      GG(7,3)=0.D+0
      GG(7,5)=0.D+0
      DO 15 I=7,9
   15 GG(7,I)=0.D+0
      GG(8,1)=0.D+0
      GG(8,2)=0.1D+5
      GG(8,3)=0.D+0
      GG(8,5)=0.D+0
      GG(8,6)=0.D+0
      DO 16 I=8,10
   16 GG(8,I)=0.D+0
      GG(9,1)=0.D+0
      GG(9,2)=0.1D+5
      GG(9,3)=0.D+0
      DO 17 I=5,7
   17 GG(9,I)=0.D+0
      GG(9,9)=0.D+0
      GG(9,10)=0.D+0
      GG(10,1)=0.D+0
      GG(10,2)=0.1D+5
      DO 18 I=5,10
   18 GG(10,I)=0.D+0
      GG(11,1)=0.D+0
      GG(11,2)=0.1D+5
      DO 19 I=5,10
   19 GG(11,I)=0.D+0
      GG(12,1)=0.D+0
      GG(12,2)=0.1D+5
      DO 20 I=5,10
   20 GG(12,I)=0.D+0
      GG(13,1)=0.D+0
      GG(13,2)=0.1D+5
      DO 21 I=5,10
   21 GG(13,I)=0.D+0
      GG(14,1)=0.D+0
      GG(14,2)=0.D+0
      GG(14,7)=0.D+0
      GG(14,9)=0.D+0
      GG(14,10)=0.0D+0
      DO 8 I=1,8
    8 GG(15,I)=0.D+0
      GG(15,9)=0.1D+1
      GG(15,10)=0.1D+1
      LEX=.FALSE.
      NEX=1
      FEX=-0.44300879D+04
      XEX(1)=.14727222
      XEX(2)=0.1
      XEX(3)=0.81D-2
      XEX(4)=0.62871731D+03
      XEX(5)=0.17D-2
      XEX(6)=0.11816143D-02
      XEX(7)=0.27D-2
      XEX(8)=0.135D-02
      XEX(9)=0.15740741
      XEX(10)=0.97592593D-01
      DO 45 I=3,10
   45 GF(I)=0.0D+0
      RETURN
    2 FX=0.2D+5*(0.15D+0*X(1)+0.14D+2*X(2)-0.6D-1)/(0.2D-2+X(1)+0.6D+2
     /    *X(2))
      FX=-FX
      RETURN
    3 CONTINUE
      GF(1)=((0.15D+0*(0.2D-2+X(1)+0.6D+2*X(2))-(0.15D+0*X(1)+0.14D+2
     /    *X(2)-0.6D-1))*0.2D+5)/(0.2D-2+X(1)+0.6D+2*X(2))**0.2D+1
      GF(2)=((0.14D+2*(0.2D-2+X(1)+0.6D+2*X(2))-((0.15D+0*X(1)+0.14D+2
     /    *X(2)-0.6D-1)*0.6D+2))*0.2D+5)/(0.2D-2+X(1)
     /      +0.6D+2*X(2))**0.2D+1
      GF(1)=-GF(1)
      GF(2)=-GF(2)
      RETURN
    4 CONTINUE
      DO I=1,10
        IF (X(I).LT.XL(I)) X(I)=XL(I)
      ENDDO
      IF (INDEX1(1)) G(1)=X(1)-0.75D+0/X(3)/X(4)
      IF (INDEX1(2)) G(2)=X(1)-X(9)/X(5)/X(4)
      IF (INDEX1(3)) G(3)=X(1)-X(10)/X(6)/X(4)-0.1D+2/X(4)
      IF (INDEX1(4)) G(4)=X(1)-0.19D+0/X(7)/X(4)-0.1D+2/X(4)
      IF (INDEX1(5)) G(5)=X(1)-0.125D+0/X(8)/X(4)
      IF (INDEX1(6)) G(6)=0.1D+5*X(2)-0.131D-2*X(9)*X(5)**0.666D+0
     1 *X(4)**1.5D+0
      IF (INDEX1(7)) G(7)=0.1D+5*X(2)-0.1038D-2*X(10)*X(6)**0.16D+1
     1 *X(4)**3
      IF (INDEX1(8)) G(8)=0.1D+5*X(2)-0.223D-3*X(7)**0.666D+0
     1 *X(4)**1.5D+0
      IF (INDEX1(9)) G(9)=0.1D+5*X(2)-0.76D-4*X(8)**3.55D+0
     1 *X(4)**5.66D+0
      IF (INDEX1(10)) G(10)=0.1D+5*X(2)-0.698D-3*X(3)**1.2D+0
     1 *X(4)**2
      IF (INDEX1(11)) G(11)=0.1D+5*X(2)-0.5D-4*X(3)**1.6D+0
     1 *X(4)**3.D+0
      IF (INDEX1(12)) G(12)=0.1D+5*X(2)-0.654D-5*X(3)**2.42D+0
     1 *X(4)**4.17D+0
      IF (INDEX1(13)) G(13)=0.1D+5*X(2)-0.257D-3*X(3)**0.666D+0
     1 *X(4)**1.5D+0
      IF (INDEX1(14)) G(14)=0.3D+2-0.2003D+1*X(5)*X(4)-0.1885D+1*X(6)
     1 *X(4)-0.184D+0*X(8)*X(4)-0.2D+1*X(3)**0.803D+0*X(4)
      IF (INDEX1(15)) G(15)=X(9)+X(10)-0.255D+0
      RETURN
    5 IF (.NOT.INDEX2(1)) GOTO 50
      GG(1,3)=0.75D+0/X(3)**0.2D+1/X(4)
      GG(1,4)=0.75D+0/X(3)/X(4)**0.2D+1
   50 IF (.NOT.INDEX2(2)) GOTO 51
      GG(2,4)=X(9)/X(5)/X(4)**0.2D+1
      GG(2,5)=X(9)/X(5)**0.2D+1/X(4)
      GG(2,9)=-0.1D+1/X(5)/X(4)
   51 IF (.NOT.INDEX2(3)) GOTO 52
      GG(3,4)=X(10)/X(6)/X(4)**0.2D+1+0.1D+2/X(4)**0.2D+1
      GG(3,6)=X(10)/X(6)**0.2D+1/X(4)
      GG(3,10)=-0.1D+1/X(6)/X(4)
   52 IF (.NOT.INDEX2(4)) GOTO 53
      GG(4,4)=0.19D+0/X(7)/X(4)**0.2D+1+0.1D+2/X(4)**0.2D+1
      GG(4,7)=0.19D+0/X(7)**0.2D+1/X(4)
   53 IF (.NOT.INDEX2(5)) GOTO 54
      GG(5,4)=0.125D+0/X(8)/X(4)**0.2D+1
      GG(5,8)=0.125D+0/X(8)**0.2D+1/X(4)
   54 IF (.NOT.INDEX2(6)) GOTO 55
      GG(6,4)=-0.15D+1*0.131D-2*X(9)*X(5)**0.666D+0*X(4)**0.5D+0
      GG(6,5)=-0.666D+0*0.131D-2*X(9)/X(5)**0.334D+0*X(4)**0.15D+1
      GG(6,9)=-0.131D-2*X(5)**0.666D+0*X(4)**0.15D+1
   55 IF (.NOT.INDEX2(7)) GOTO 56
      GG(7,4)=-0.3D+1*0.1038D-2*X(10)*X(6)**0.16D+1*X(4)**0.2D+1
      GG(7,6)=-0.16D+1*0.1038D-2*X(10)*X(6)**0.6D+0*X(4)**0.3D+1
      GG(7,10)=-0.1038D-2*X(6)**0.16D+1*X(4)**0.3D+1
   56 IF (.NOT.INDEX2(8)) GOTO 57
      GG(8,4)=-0.15D+1*0.223D-3*X(7)**0.666D+0*X(4)**0.5D+0
      GG(8,7)=-0.666D+0*0.223D-3/X(7)**0.334D+0*X(4)**0.15D+1
   57 IF (.NOT.INDEX2(9)) GOTO 58
      GG(9,4)=-0.566D+1*0.76D-4*X(8)**0.355D+1*X(4)**0.466D+1
      GG(9,8)=-0.355D+1*0.76D-4*X(8)**0.255D+1*X(4)**0.566D+1
   58 IF (.NOT.INDEX2(10)) GOTO 59
      GG(10,3)=-0.12D+1*0.698D-3*X(3)**0.2D+0*X(4)**0.2D+1
      GG(10,4)=-0.2D+1*0.698D-3*X(3)**0.12D+1*X(4)
   59 IF (.NOT.INDEX2(11)) GOTO 60
      GG(11,3)=-0.16D+1*0.5D-4*X(3)**0.6D+0*X(4)**0.3D+1
      GG(11,4)=-0.3D+1*0.5D-4*X(3)**0.16D+1*X(4)**0.2D+1
   60 IF (.NOT.INDEX2(12)) GOTO 61
      GG(12,3)=-0.242D+1*0.654D-5*X(3)**0.142D+1*X(4)**0.417D+1
      GG(12,4)=-0.417D+1*0.654D-5*X(3)**0.242D+1*X(4)**0.317D+1
   61 IF (.NOT.INDEX2(13)) GOTO 62
      GG(13,3)=-0.666D+0*0.257D-3/X(3)**0.334D+0*X(4)**0.15D+1
      GG(13,4)=-0.15D+1*0.257D-3*X(3)**0.666D+0*X(4)**0.5D+0
   62 IF (.NOT.INDEX2(14)) GOTO 63
      GG(14,3)=-0.803D+0*0.2D+1/X(3)**0.197D+0*X(4)
      GG(14,4)=-0.2003D+1*X(5)-0.1885D+1*X(6)-0.184D+0*X(8)
     1-0.2D+1*X(3)**0.803D+0
      GG(14,5)=-0.2003D+1*X(4)
      GG(14,6)=-0.1885D+1*X(4)
      GG(14,8)=-0.184D+0*X(4)
   63 RETURN
      END	
C
      SUBROUTINE TP377(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION A(10),SUM
      DATA A/-0.6089D+1,-0.17164D+2,-0.34054D+2,-0.5914D+1,
     1 -0.24721D+2,-0.14986D+2,-0.24100D+2,-0.10708D+2,
     1 -0.26662D+2,-0.22179D+2/	
      GOTO (1,2,3,4,5),MODE
    1 N=10
      NILI=0
      NINL=0
      NELI=3
      NENL=0
      DO 6 I=1,10
      X(I)=0.1D+0
      LXU(I)=.TRUE.
      LXL(I)=.TRUE.
      XL(I)=0.1D-3
    6 XU(I)=0.1D+2
      GG(1,1)=0.1D+1
      GG(1,2)=-0.2D+1
      GG(1,3)=0.2D+1
      GG(1,4)=0.D+0
      GG(1,5)=0.D+0
      GG(1,6)=0.1D+1
      DO 7 I=7,9
    7 GG(1,I)=0.D+0
      GG(1,10)=0.1D+1
      DO 8 I=1,3
    8 GG(2,I)=0.D+0
      GG(2,4)=0.1D+1
      GG(2,5)=-0.2D+1
      GG(2,6)=0.1D+1
      GG(2,7)=0.1D+1
      DO 9 I=8,10
    9 GG(2,I)=0.D+0
      GG(3,1)=0.D+0
      GG(3,2)=0.D+1
      GG(3,3)=0.1D+1
      DO 10 I=4,6
   10 GG(3,I)=0.D+0
      GG(3,7)=0.1D+1
      GG(3,8)=0.1D+1
      GG(3,9)=0.2D+1
      GG(3,10)=0.1D+1
      LEX=.FALSE.
      NEX=1
      FEX=-795.001
      XEX(1)=10.0
      XEX(2)=10.0
      XEX(3)=1.0
      XEX(4)=10.0      	        			  	
      XEX(5)=9.5
      XEX(6)=10.0
      XEX(7)=0.1D-3
      XEX(8)=0.1D-3
      XEX(9)=0.1D-3
      XEX(10)=0.1D-3
      RETURN
    2 FX=0.D+0
      SUM=0.D+0
      DO 11 I=1,10
   11 SUM=SUM+X(I)
      DO 12 I=1,10
   12 FX=FX+X(I)*(A(I)+DLOG(DMAX1(X(I)/SUM,1.0D-5)))
    3 RETURN
      SUM=0.D+0
      DO 45 I=1,10
   45 SUM=SUM+X(I)
      DO 46 I=1,10
   46 GF(I)=A(I)+DLOG(X(I)/SUM)
      RETURN
    4 IF(INDEX1(1)) G(1)=X(1)-2.D+0*X(2)+2.D+0*X(3)+X(6)+X(10)-2.D+0
      IF(INDEX1(2)) G(2)=X(4)-2.D+0*X(5)+X(6)+X(7)-1.D+0
      IF(INDEX1(3)) G(3)=X(3)+X(7)+X(8)+2.D+0*X(9)+X(10)-1.D+0
    5 RETURN
      END
C
      SUBROUTINE TP378(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION A(10),CON
      DATA A/-0.6089D+1,-0.17164D+2,-0.34054D+2,-0.5914D+1,
     1 -0.24721D+2,-0.14986D+2,-0.24100D+2,-0.10708D+2,
     1 -0.26662D+2,-0.22179D+2/	
      GOTO(1,2,3,4,5),MODE
    1 N=10
      NILI=0	
      NINL=0
      NELI=0
      NENL=3
      DO 6 I=1,10
      X(I)=-0.23D+1
      XU(I)=-0.1
      LXU(I)=.TRUE.
      XL(I)=-16.0
    6 LXL(I)=.TRUE.
      LEX=.TRUE.
      NEX=1
C      FEX=-0.47760D+2
      FEX=-0.47761091D+2 
      XEX(1)=-0.32024D+1
      XEX(2)=-0.19123D+1
      XEX(3)=-0.2444D+0
      XEX(4)=-0.15670D+2
      XEX(5)=-0.7217D+0
      XEX(6)=-0.72736D+1
      XEX(7)=-0.35965D+1
      XEX(8)=-0.40206D+1
      XEX(9)=-0.32885D+1
      XEX(10)=-0.23344D+1
      GG(1,4)=0.0D+0
      GG(1,5)=0.0D+0
      GG(1,7)=0.0D+0
      GG(1,8)=0.0D+0
      GG(1,9)=0.0D+0
      GG(2,1)=0.0D+0
      GG(2,2)=0.0D+0
      GG(2,3)=0.0D+0
      GG(2,8)=0.0D+0
      GG(2,9)=0.0D+0
      GG(2,10)=0.0D+0
      GG(3,1)=0.0D+0
      GG(3,2)=0.0D+0
      GG(3,4)=0.0D+0
      GG(3,5)=0.0D+0
      GG(3,6)=0.0D+0
      RETURN
    2 FX=0.0D+0
      CON=0.0D+0
      DO 7 J=1,10
    7 CON=CON+DEXP(X(J))
      CON=DLOG(CON)
      DO 8 I=1,10
    8 FX=FX+DEXP(X(I))*(A(I)+X(I)-CON)
    3 RETURN
      CON=0.0D+0
      DO 45 J=1,10
   45 CON=CON+DEXP(X(J))
      CON=DLOG(CON)
      DO 46 I=1,10
      GF(I)=DEXP(X(I))*(A(I)+X(I)-CON)
   46 CONTINUE 
      RETURN    
    4 IF(INDEX1(1)) G(1)=DEXP(X(1))+0.2D+1*DEXP(X(2))+0.2D+1*DEXP(X(3))+
     +DEXP(X(6))+DEXP(X(10))-0.2D+1
      IF(INDEX1(2)) G(2)=DEXP(X(4))+0.2D+1*DEXP(X(5))+DEXP(X(6))+
     +DEXP(X(7))-0.1D+1
      IF(INDEX1(3)) G(3)=DEXP(X(3))+DEXP(X(7))+DEXP(X(8))+0.2D+1*
     *DEXP(X(9))+DEXP(X(10))-0.1D+1
      RETURN
    5 IF(.NOT.INDEX2(1)) GOTO 47
      GG(1,1)=DEXP(X(1))
      GG(1,2)=DEXP(X(2))*0.2D+1
      GG(1,3)=DEXP(X(3))*0.2D+1
      GG(1,6)=DEXP(X(6))
      GG(1,10)=DEXP(X(10))
   47 IF(.NOT.INDEX2(2)) GOTO 48
      GG(2,4)=DEXP(X(4))
      GG(2,5)=DEXP(X(5))*0.2D+1
      GG(2,6)=DEXP(X(6))
      GG(2,7)=DEXP(X(7))
   48 IF(.NOT.INDEX2(3)) GOTO 49
      GG(3,3)=DEXP(X(3))
      GG(3,7)=DEXP(X(7))
      GG(3,8)=DEXP(X(8))
      GG(3,9)=DEXP(X(9))*0.2D+1
      GG(3,10)=DEXP(X(10))
   49 RETURN
      END
C      
      SUBROUTINE TP379(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX)
      INTEGER I,J,LSUM
      DOUBLEPRECISION F,DF,Y(65),T
      DATA Y/1.366D+0,1.191D+0,1.112D+0,1.013D+0,0.991D+0,
     1 0.885D+0,0.831D+0,0.847D+0,0.786D+0,0.725D+0,0.746D+0,
     2 0.679D+0,0.608D+0,0.655D+0,0.616D+0,0.606D+0,0.602D+0,
     3 0.626D+0,0.651D+0,0.724D+0,0.649D+0,0.649D+0,0.694D+0,
     4 0.644D+0,0.624D+0,0.661D+0,0.612D+0,0.558D+0,0.533D+0,
     5 0.495D+0,0.500D+0,0.423D+0,0.395D+0,0.375D+0,0.372D+0,
     6 0.391D+0,0.396D+0,0.405D+0,0.428D+0,0.429D+0,0.523D+0,
     7 0.562D+0,0.607D+0,0.653D+0,0.672D+0,0.708D+0,0.633D+0,
     8 0.668D+0,0.645D+0,0.632D+0,0.591D+0,0.559D+0,0.597D+0,
     9 0.625D+0,0.739D+0,0.710D+0,0.729D+0,0.720D+0,0.636D+0,
     A 0.581D+0,0.428D+0,0.292D+0,0.162D+0,0.098D+0,0.054D+0/
      GOTO(1,2,2,4,4),MODE
    1 N=11
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      LSUM=65
      X(1)=1.3D+0
      X(2)=0.65D+0
      X(3)=0.65D+0
      X(4)=0.7D+0
      X(5)=0.6D+0
      X(6)=3.0D+0
      X(7)=5.0D+0
      X(8)=7.0D+0
      X(9)=2.0D+0
      X(10)=4.5D+0
      X(11)=5.5D+0
      DO 6 I=1,11
      LXU(I)=.FALSE.
      XL(I)=0.0
    6 LXL(I)=.TRUE.
      LEX=.FALSE.
      NEX=1
      FEX=0.401377D-1
      XEX(1)=0.130997D+1
      XEX(2)=0.431554D+0
      XEX(3)=0.63366D+0
      XEX(4)=0.59943D+0
      XEX(5)=0.754183D+0
      XEX(6)=0.904286D+0
      XEX(7)=0.136581D+1
      XEX(8)=0.482369D+1
      XEX(9)=0.239868D+1
      XEX(10)=0.456887D+1
      XEX(11)=0.567534D+1
      RETURN
    2 CONTINUE
      DO I=1,11
      IF (X(I).GT.XU(I)) X(I)=XU(I)
      IF (X(I).LT.XL(I)) X(I)=XL(I)
      ENDDO 
      DO 7 I=1,65
      T=.1D+0*DBLE(I-1)
    7 F(I)=Y(I)-(X(1)*DEXP(-X(5)*T)+X(2)*DEXP(-X(6)*(T-X(9))**2)
     2 +X(3)*DEXP(-X(7)*(T-X(10))**2)+X(4)*DEXP(-X(8)*(T-X(11))**2))
      IF (MODE.EQ.3) GOTO 3
      FX=0.D+0
      DO 70 I=1,65
   70 FX=FX+F(I)**2
      RETURN
    3 DO 8 I=1,65
      T=.1D+0*DBLE(I-1)
      DF(I,1)=-DEXP(-X(5)*T)
      DF(I,2)=-DEXP(-X(6)*(T-X(9))**2)
      DF(I,3)=-DEXP(-X(7)*(T-X(10))**2)
      DF(I,4)=-DEXP(-X(8)*(T-X(11))**2)
      DF(I,5)=X(1)*T*DEXP(-X(5)*T)
      DF(I,6)=X(2)*(T-X(9))**2*DEXP(-X(6)*(T-X(9))**2)
      DF(I,7)=X(3)*(T-X(10))**2*DEXP(-X(7)*(T-X(10))**2)
      DF(I,8)=X(4)*(T-X(11))**2*DEXP(-X(8)*(T-X(11))**2)
      DF(I,9)=-X(2)*X(6)*2.D+0*(T-X(9))*DEXP(-X(6)*(T-X(9))**2)
      DF(I,10)=-X(3)*X(7)*2.D+0*(T-X(10))*DEXP(-X(7)*(T-X(10))**2)
    8 DF(I,11)=-X(4)*X(8)*2.D+0*(T-X(11))*DEXP(-X(8)*(T-X(11))**2)
      DO 19 I=1,11
   19 GF(I)=0.D+0
      DO 20 J=1,11
      DO 20 I=1,65
   20 GF(J)=GF(J)+2.D+0*F(I)*DF(I,J)
    4 RETURN
      END
C
      SUBROUTINE TP380(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION A(11),TEMP,C(30)
      DATA A/-0.00133172D+0,-0.002270927D+0,-0.00248546D+0,
     1 -0.467D+1,-0.4671973D+1,-0.00814D+0,-0.008092D+0,
     2 -0.005D+0,-0.000909D+0,-0.00088D+0,-0.00119D+0/
      DATA C/5.367373D-2,2.1863746D-2,9.7733533D-2,6.6940803D-3,
     1 1.0D-6,1.0D-5,1.0D-6,1.0D-10,1.0D-8,1.0D-2,1.0D-4,
     2 1.0898645D-1,1.6108052D-4,1.0D-23,1.9304541D-6,1.0D-3,
     3 1.0D-6,1.0D-5,1.0D-6,1.0D-9,1.0D-9,1.0D-3,1.0D-3,
     4 1.0898645D-1,1.6108052D-5,1.0D-23,1.9304541D-8,1.0D-5,
     5 1.1184059D-4,1.0D-4/
      GOTO (1,2,3,4,5),MODE
    1 N=12
      NILI=0
      NINL=3
      NELI=0
      NENL=0
      DO 6 I=1,12
      X(I)=0.4D+1
      LXU(I)=.TRUE.
      LXL(I)=.TRUE.
      XL(I)=.1D+0
    6 XU(I)=100.D+0
      GG(1,1)=-C(1)   
      GG(1,2)=-C(2)
      GG(1,3)=-C(3)
      DO 7 I=6,12
    7 GG(1,I)=0.D+0
      GG(2,1)=-C(5)
      GG(2,3)=-C(7)
      GG(2,8)=0.D+0
      GG(2,9)=0.D+0
      GG(2,11)=0.D+0
      GG(3,3)=-C(19)
      GG(3,6)=-C(22)
      GG(3,7)=0.D+0
      GG(3,8)=-C(23)
      GG(3,10)=0.D+0
      GG(3,11)=-C(30)
      GG(3,12)=0.D+0
      LEX=.FALSE.
      NEX=1
      FEX=0.31682215D+01
c      FEX=FEX*1.0D-6
      XEX(1)=0.26631947068D+1
      XEX(2)=0.4517277762D+1
      XEX(3)=0.7133802907D+1
      XEX(4)=0.2237268448D+1
      XEX(5)=0.407840382657D+1
      XEX(6)=0.131827569D+1
      XEX(7)=0.4125187034D+1
      XEX(8)=0.2856195978D+1
      XEX(9)=0.16765929748D+1
      XEX(10)=0.21789111052D+1
      XEX(11)=0.512343515D+1
      XEX(12)=0.6659338016D+1
      RETURN
c    2 FX=0.1D+6
    2 FX=0.1D+1
      DO 20 I=1,11
      TEMP=X(I)
      IF(X(I).LT.0.1D-14) TEMP=0.1D-14
   20 FX=FX*TEMP**A(I)
      FX=1.0D5*FX
      RETURN
    3 DO 10 I=1,11
      TEMP=X(I)
      IF(X(I).LT.0.1D-14) TEMP=0.1D-14
      GF(I)=1.0D5*FX*(A(I)/TEMP)
   10 CONTINUE   
      GF(12)=0.D+0
      RETURN
    4 IF(INDEX1(1)) G(1)=1.D+0-C(1)*X(1)-C(2)*X(2)-C(3)*X(3)-C(4)*X(4)
     1 *X(5)
      IF(INDEX1(2)) G(2)=1.D+0-C(5)*X(1)-C(6)*X(2)-C(7)*X(3)-C(8)
     1 *X(4)*X(12)-C(9)*X(5)/X(12)-C(10)*X(6)/X(12)-C(11)*X(7)
     2 *X(12)-C(12)*X(4)*X(5)-C(13)*X(2)*X(5)/X(12)-C(14)*X(2)
     3 *X(4)*X(5)-C(15)*X(2)/X(4)*X(5)/X(12)**2-C(16)*X(10)/X(12)
      IF(INDEX1(3)) G(3)=1.0D+0-C(17)*X(1)-C(18)*X(2)-C(19)*X(3)
     1 -C(20)*X(4)-C(21)*X(5)-C(22)*X(6)-C(23)*X(8)-C(24)*X(4)*X(5)
     2 -C(25)*X(2)*X(5)-C(26)*X(2)*X(4)*X(5)-C(27)*X(2)*X(5)/X(4)
     3 -C(28)*X(9)-C(29)*X(1)*X(9)-C(30)*X(11)
      RETURN
    5 IF(.NOT.INDEX2(1)) GOTO 50
      GG(1,4)=-C(4)*X(5)
      GG(1,5)=-C(4)*X(4)
   50 IF(.NOT.INDEX2(2)) GOTO 51
      GG(2,2)=-C(6)-C(13)*X(5)/X(12)-C(14)*X(4)*X(5)-C(15)/X(4)*X(5)/
     1 X(12)**2
      GG(2,4)=-C(8)*X(12)-C(12)*X(5)-C(14)*X(2)*X(5)+C(15)*X(2)/X(4)
     1 **2*X(5)/X(12)**2
      GG(2,5)=-C(9)/X(12)-C(12)*X(4)-C(13)*X(2)/X(12)-C(14)*X(2)*X(4)-
     1 C(15)*X(2)/X(4)/X(12)**2
      GG(2,6)=-C(10)/X(12)
      GG(2,7)=-C(11)*X(12)
      GG(2,10)=-C(16)/X(12)
      GG(2,12)=-C(8)*X(4)+C(9)*X(5)/X(12)**2+C(10)*X(6)/X(12)**2
     1 -C(11)*X(7)+C(13)*X(2)*X(5)/X(12)**2+2*C(15)*X(2)/X(4)
     2 *X(5)/X(12)**3.D+0+C(16)*X(10)/X(12)**2
   51 IF(.NOT.INDEX2(3)) GOTO 52
      GG(3,1)=-C(17)-C(29)*X(9)
      GG(3,2)=-C(18)-C(25)*X(5)-C(26)*X(4)*X(5)-C(27)*X(5)/X(4)
      GG(3,4)=-C(20)-C(24)*X(5)-C(26)*X(2)*X(5)+C(27)*X(2)*X(5)/
     1 X(4)**2
      GG(3,5)=-C(21)-C(24)*X(4)-C(25)*X(2)-C(26)*X(2)*X(4)-C(27)*X(2)
     1 /X(4)
      GG(3,9)=-C(28)-C(29)*X(1)
   52 RETURN
      END
C
      SUBROUTINE TP381(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION S(13),U(13),V(13),R(13)
      DATA R/0.8D+0,0.11D+1,0.85D+0,0.345D+1,0.2D+1,0.21D+1,
     1 0.3D+1,0.8D+0,0.45D+0,0.72D+0,0.18D+1,0.3D+1,0.6D+0/
      DATA S/0.116D+2,0.137D+2,0.95D+1,0.485D+2,0.319D+2,0.511D+2,
     1 0.655D+2,0.0D+0,0.D+0,0.0D+0,0.218D+2,0.469D+2,0.D+0/
      DATA U/0.5D-1,0.7D-1,0.D+0,0.33D+0,0.0D+0,0.127D+1,0.127D+1,
     1 0.2335D+2,0.3584D+2,0.81D+0,0.179D+1,0.734D+1,0.0D+0/
      DATA V/0.35D+0,0.37D+0,0.1D+0,0.62D+0,0.0D+0,0.103D+1,
     1 0.169D+1,0.1821D+2,0.1D-1,0.8D-1,0.31D+0,0.159D+1,0.2245D+2/	
      GOTO(1,2,3,4,5),MODE
    1 N=13
      NILI=3
      NINL=0
      NELI=1
      NENL=0
      DO 6 I=1,13
      X(I)=0.1D+0
      LXU(I)=.FALSE.
      LXL(I)=.TRUE.
    6 XL(I)=0.0D+0  	
      DO 7 I=1,13
      GG(1,I)=S(I)
      GG(2,I)=U(I)
    7 GG(4,I)=0.1D+1
      DO 8 I=1,13
    8 GG(3,I)=V(I)
      LEX=.FALSE.
      NEX=1
      FEX=1.0149D0
      XEX(1)=0.7864989D+0
      XEX(2)=0.0D+0
      XEX(3)=0.0D+0
      XEX(4)=0.0D+0
      XEX(5)=0.0D+0
      XEX(6)=0.17105694D+0
      XEX(7)=0.0D+0
      XEX(8)=0.0D+0
      XEX(9)=0.20676337D-1
      XEX(10)=0.0D+0
      XEX(11)=0.0D+0
      XEX(12)=0.0D+0
      XEX(13)=0.19883712D-1
      DO 11 I=1,13
   11 GF(I)=R(I)
      RETURN
    2 FX=0.0D+0
      DO 10 I=1,13
   10 FX=FX+R(I)*X(I)
    3 RETURN
    4 IF(.NOT.INDEX1(1)) GOTO 12
      G(1)=0.0D+0
      DO 13 I=1,13
   13 G(1)=G(1)+S(I)*X(I)
      G(1)=G(1)-0.18D+2
   12 IF(.NOT.INDEX1(2)) GOTO 14
      G(2)=0.0D+0
      DO 15 I=1,13
   15 G(2)=G(2)+U(I)*X(I)
      G(2)=G(2)-0.1D+1
   14 IF(.NOT.INDEX1(3)) GOTO 16
      G(3)=0.0D+0
      DO 17 I=1,13
   17 G(3)=G(3)+V(I)*X(I)
      G(3)=G(3)-0.9D+0
   16 IF(.NOT.INDEX1(4)) GOTO 5
      G(4)=0.0D+0
      DO 19 I=1,13
   19 G(4)=G(4)+X(I)
      G(4)=G(4)-0.1D+1
    5 RETURN
      END   
C           
      SUBROUTINE TP382(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION R(13),S(13),U(13),V(13),Z1(13),Z2(13),Z3(13),HELP
      DATA R/0.8D+0,1.1D+0,0.85D+0,3.45D+0,2.D+0,2.1D+0,3.D+0,0.8D+0,
     1 0.45D+0,0.72D+0,1.8D+0,3.D+0,0.6D+0/
      DATA S/11.6D+0,13.7D+0,9.5D+0,48.5D+0,31.9D+0,51.1D+0,65.5D+0,
     1 0.D+0,0.D+0,0.D+0,21.8D+0,46.9D+0,0.D+0/
      DATA Z1/0.4844D+0,0.3003D+0,0.1444D+0,0.0588D+0,4.9863D+0,
     1 0.0653D+0,21.0222D+0,0.D+0,0.D+0,0.D+0,0.2970D+0,9.2933D+0,
     2 0.D+0/
      DATA U/0.05D+0,0.07D+0,0.D+0,0.33D+0,0.D+0,1.27D+0,1.27D+0,
     1 23.35D+0,35.84D+0,0.81D+0,1.79D+0,7.34D+0,0.D+0/
      DATA Z2/0.0001D+0,0.D+0,0.D+0,0.D+0,0.D+0,0.0040D+0,0.1404D+0,
     1 1.3631D+0,0.5138D+0,0.0289D+0,0.0097D+0,0.3893D+0,0.D+0/
      DATA V/0.35D+0,0.37D+0,0.1D+0,0.62D+0,0.D+0,1.03D+0,1.69D+0,
     1 18.21D+0,0.01D+0,0.08D+0,0.31D+0,1.59D+0,22.45D+0/
      DATA Z3/0.001D+0,0.0009D+0,0.0001D+0,0.0005D+0,0.D+0,0.0021D+0,
     1 0.0825D+0,0.2073D+0,0.D+0,0.0004D+0,0.0005D+0,0.0107D+0,
     2 1.0206D+0/
      GOTO(1,2,3,4,5),MODE    
    1 N=13
      NILI=0
      NINL=3
      NELI=1
      NENL=0
      DO 6 I=1,13
      X(I)=0.1D+0
      LXU(I)=.FALSE.
      LXL(I)=.TRUE.
      XL(I)=0.0D+0
    6 GG(4,I)=1.0D+0
      LEX=.FALSE.
      NEX=1
      FEX=1.03831D+0
      XEX(1)=0.769578D+0
      XEX(2)=0.D+0
      XEX(3)=0.D+0
      XEX(4)=0.D+0
      XEX(5)=0.D+0
      XEX(6)=0.18880588D+0
      XEX(7)=0.D+0
      XEX(8)=0.D+0
      XEX(9)=0.20625416D+0
      XEX(10)=0.D+0
      XEX(11)=0.D+0
      XEX(12)=0.D+0
      XEX(13)=0.020990697D+0
      DO 8 I=1,13
    8 GF(I)=R(I)
      RETURN
    2 FX=0.D+0
      DO 7 I=1,13
    7 FX=FX+R(I)*X(I)
    3 RETURN
    4 IF(.NOT.INDEX1(1)) GOTO 9
      G(1)=0.D+0
      DO 10 I=1,13
   10 G(1)=G(1)+Z1(I)*X(I)**2
      G(1)=-DSQRT(G(1))*0.1645D+1-0.18D+2
      DO 11 I=1,13
   11 G(1)=G(1)+S(I)*X(I)
    9 IF(.NOT.INDEX1(2)) GOTO 12
      G(2)=0.D+0               	
      DO 13 I=1,13
   13 G(2)=G(2)+Z2(I)*X(I)**2
      G(2)=-DSQRT(G(2))*0.1645D+1-0.1D+1
      DO 14 I=1,13
   14 G(2)=G(2)+U(I)*X(I)
   12 IF(.NOT.INDEX1(3)) GOTO 15
      G(3)=0.0D+0
      DO 16 I=1,13
   16 G(3)=G(3)+Z3(I)*X(I)**2
      G(3)=-DSQRT(G(3))*1.645D+0-0.9D+0
      DO 17 I=1,13
   17 G(3)=G(3)+V(I)*X(I)
   15 IF(.NOT.INDEX1(4)) GOTO 18
      G(4)=-1.D+0
      DO 19 I=1,13
   19 G(4)=G(4)+X(I)
   18 RETURN
    5 IF(.NOT.INDEX2(1)) GOTO 27
      HELP=0.D+0	
      DO 21 I=1,13
   21 HELP=HELP+Z1(I)*X(I)**2
      HELP=-1.645D+0/2.D+0/DSQRT(HELP)
      DO 22 I=1,13
   22 GG(1,I)=S(I)+HELP*2.D+0*Z1(I)*X(I)
   27 IF(.NOT.INDEX2(2)) GOTO 28
      HELP=0.D+0
      DO 23 I=1,13
   23 HELP=HELP+Z2(I)*X(I)**2
      HELP=-1.645D+0/2.D+0/DSQRT(HELP)
      DO 24 I=1,13
   24 GG(2,I)=U(I)+HELP*2.D+0*Z2(I)*X(I)
   28 IF(.NOT.INDEX2(3)) GOTO 20
      HELP=0.D+0
      DO 25 I=1,13
   25 HELP=HELP+Z3(I)*X(I)**2
      HELP=-1.645D+0/2.D+0/DSQRT(HELP)
      DO 26 I=1,13
   26 GG(3,I)=V(I)+HELP*2.D+0*Z3(I)*X(I)
   20 RETURN
      END   
C      
      SUBROUTINE TP383(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DIMENSION A(14),B(14),C(14)
      DOUBLEPRECISION A,B,C
      DATA A/0.12842275D+5,0.63425D+3,0.63425D+3,0.634125D+3,
     1 0.1268D+4,0.633875D+3,0.63375D+3,0.1267D+4,0.76005D+3,
     2 0.63325D+3,0.126625D+4,0.632875D+3,0.39446D+3,0.940838D+3/
      DATA B/0.25D+2,0.26D+2,0.26D+2,0.27D+2,0.28D+2,0.29D+2,
     1 0.30D+2,0.32D+2,0.33D+2,0.34D+2,0.35D+2,0.37D+2,
     2 0.38D+2,0.36D+2/
      DATA C/0.547934D+1,0.83234D+0,0.94749D+0,0.111082D+1,
     1 0.264824D+1,0.155868D+1,0.173215D+1,0.390896D+1,
     2 0.274284D+1,0.260541D+1,0.596184D+1,0.329522D+1,
     3 0.183517D+1,0.281372D+1/
      GOTO (1,2,3,4,5),MODE
    1 N=14
      NILI=0
      NINL=0
      NELI=1
      NENL=0
      DO 6 I=1,14
      X(I)=0.1D-1
      LXU(I)=.TRUE.
      LXL(I)=.TRUE.
      XL(I)=0.1D-3
      XU(I)=0.1D+1/B(I)
    6 GG(1,I)=C(I)
      LEX=.FALSE.
      NEX=1
      FEX=0.728566D+6
      XEX(1)=0.4D-1
      XEX(2)=0.382D-1
      XEX(3)=0.358D-1
      XEX(4)=0.33D-1
      XEX(5)=0.303D-1
      XEX(6)=0.279D-1
      XEX(7)=0.265D-1
      XEX(8)=0.249D-1
      XEX(9)=0.23D-1
      XEX(10)=0.216D-1
      XEX(11)=0.202D-1
      XEX(12)=0.192D-1
      XEX(13)=0.203D-1
      XEX(14)=0.253D-1
      RETURN
    2 CONTINUE
      DO I=1,14
      IF (X(I).LT.XL(I)) X(I)=XL(I)
      ENDDO 
      FX=0.D+0
      DO 7 I=1,14
    7 FX=FX+A(I)/X(I)
      RETURN
    3 DO 8 I=1,14
    8 GF(I)=-A(I)/X(I)**2
      RETURN
    4 IF(.NOT.INDEX1(1)) GOTO 5
      G(1)=0.D+0
      DO 10 I=1,14
   10 G(1)=G(1)+(C(I)*X(I))
      G(1)=G(1)-0.1D+1
    5 RETURN
      END       	      	
C
      SUBROUTINE TP384(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION A(10,15),B(10),C,D(15)
      DATA B/3.85D+2,4.7D+2,5.6D+2,5.65D+2,6.45D+2,4.3D+2,4.85D+2,
     F 4.55D+2,3.9D+2,8.6D+2/
      DATA A/1.D+2,9.D+1,7.D+1,2*5.D+1,4.D+1,3.D+1,2.D+1,1.D+1,5.D+0,
     F 2*1.D+2,5.D+1,0.D+0,1.D+1,0.D+0,6.D+1,3.D+1,7.D+1,1.D+1,
     F 2*1.D+1,2*0.D+0,7.D+1,5.D+1,3.D+1,4.D+1,1.D+1,5.D+2,5.D+0,
     F 3.5D+1,5.5D+1,6.5D+1,6.D+1,9.5D+1,9.D+1,2.5D+1,3.5D+1,5.D+0,
     F 1.D+1,2.D+1,2.5D+1,3.5D+1,4.5D+1,5.D+1,0.D+0,4.D+1,2.5D+1,2.D+1,
     F 0.D+0,5.D+0,2*1.D+2,4.5D+1,3.5D+1,3.D+1,2.5D+1,6.5D+1,5.D+0,
     F 2*0.D+0,4.D+1,3.5D+1,0.D+0,1.D+1,5.D+0,1.5D+1,0.D+0,1.D+1,
     F 2.5D+1,3.5D+1,5.D+1,6.D+1,3.5D+1,6.D+1,2.5D+1,1.D+1,3.D+1,3.5D+1,
     F 0.D+0,5.5D+1,2*0.D+0,6.5D+1,2*0.D+0,8.D+1,0.D+0,9.5D+1,
     F 1.D+1,2.5D+1,3.D+1,1.5D+1,5.D+0,4.5D+1,7.D+1,2.D+1,0.D+0,7.D+1,
     F 5.5D+1,2.D+1,6.D+1,0.D+0,7.5D+1,1.5D+1,2.D+1,3.D+1,2.5D+1,2.D+1,
     F 5.D+0,0.D+0,1.D+1,7.5D+1,1.D+2,2.D+1,2.5D+1,3.D+1,0.D+0,1.D+1,
     F 4.5D+1,4.D+1,3.D+1,3.5D+1,7.5D+1,0.D+0,7.D+1,5.D+0,1.5D+1,3.5D+1,
     F 2.D+1,2.5D+1,0.D+0,3.D+1,1.D+1,5.D+0,1.5D+1,6.5D+1,5.D+1,1.D+1,
     F 0.D+0,1.D+1,4.D+1,6.5D+1,0.D+0,5.D+0,1.5D+1,2.D+1,5.5D+1,3.D+1/
      DATA D/4.86D+2,6.4D+2,7.58D+2,7.76D+2,4.77D+2,7.07D+2,1.75D+2,
     F 6.19D+2,6.27D+2,6.14D+2,4.75D+2,3.77D+2,5.24D+2,4.68D+2,5.29D+2/
      GOTO (1,2,3,4,5),MODE
    1 N=15
      NILI=0
      NINL=10
      NELI=0
      NENL=0
      DO 6 I=1,15
      X(I)=0.D+0
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=-0.83102590D+4
      XEX(1)=0.86095379D+0
      XEX(2)=0.91736139D+0
      XEX(3)=0.91973646D+0
      XEX(4)=0.89600562D+0
      XEX(5)=0.10372946D+1
      XEX(6)=0.97308908D+0
      XEX(7)=0.82243629D+0
      XEX(8)=0.11987219D+1
      XEX(9)=0.11563350D+1
      XEX(10)=0.11443868D+1
      XEX(11)=0.10305681D+1
      XEX(12)=0.90949479D+0
      XEX(13)=0.10820450D+1
      XEX(14)=0.84682383D+0
      XEX(15)=0.11723720D+1
      RETURN
    2 FX=0.D+0
      DO 15 I=1,15
   15 FX=FX-D(I)*X(I)
      RETURN
    3 DO 11 I=1,15
   11 GF(I)=-D(I) 
      RETURN
    4 DO 7 I=1,10
      IF (.NOT.INDEX1(I)) GOTO 7
      C=0.D+0
      DO 9 J=1,15
    9 C=C+A(I,J)*X(J)**2
      G(I)=B(I)-C
    7 CONTINUE
      RETURN
    5 DO 10 I=1,10
      IF (.NOT.INDEX2(I)) GOTO 10
      DO 13 J=1,15
   13 GG(I,J)=-2.D+0*A(I,J)*X(J)
   10 CONTINUE
      RETURN
      END
C      
      SUBROUTINE TP385(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION C,A(10,15),B(10),D(15)
      DATA B/3.85D+2,4.7D+2,5.6D+2,5.65D+2,6.45D+2,4.3D+2,4.85D+2,
     F 4.55D+2,8.9D+2,4.6D+2/
      DATA A/1.D+2,9.D+1,7.D+1,2*5.D+1,4.D+1,3.D+1,2.D+1,1.D+1,5.D+0,
     F 2*1.D+2,5.D+1,0.D+0,1.D+1,0.D+0,6.D+1,3.D+1,7.D+1,1.D+1,
     F 2*1.D+1,2*0.D+0,7.D+1,5.D+1,3.D+1,4.D+1,1.D+1,1.D+2,5.D+0,
     F 3.5D+1,5.5D+1,6.5D+1,6.D+1,9.5D+1,9.D+1,2.5D+1,3.5D+1,5.D+0,
     F 1.D+1,2.D+1,2.5D+1,3.5D+1,4.5D+1,5.D+1,0.D+0,4.D+1,2.5D+1,2.D+1,
     F 0.D+0,5.D+0,2*1.D+2,4.5D+1,3.5D+1,3.D+1,2.5D+1,6.5D+1,5.D+0,
     F 2*0.D+0,4.D+1,3.5D+1,0.D+0,1.D+1,5.D+0,1.5D+1,0.D+0,1.D+1,
     F 2.5D+1,3.5D+1,5.D+1,6.D+1,3.5D+1,6.D+1,2.5D+1,1.D+1,3.D+1,3.5D+1,
     F 0.D+0,5.5D+1,2*0.D+0,6.5D+1,2*0.D+0,8.D+1,5.D+2,9.5D+1,
     F 1.D+1,2.5D+1,3.D+1,1.5D+1,5.D+0,4.5D+1,7.D+1,2.D+1,0.D+0,7.D+1,
     F 5.5D+1,2.D+1,6.D+1,0.D+0,7.5D+1,1.5D+1,2.D+1,3.D+1,2.5D+1,2.D+1,
     F 5.D+0,0.D+0,1.D+1,7.5D+1,1.D+2,2.D+1,2.5D+1,3.D+1,0.D+0,1.D+1,
     F 4.5D+1,4.D+1,3.D+1,3.5D+1,7.5D+1,0.D+0,7.D+1,5.D+0,1.5D+1,3.5D+1,
     F 2.D+1,2.5D+1,0.D+0,3.D+1,1.D+1,5.D+0,1.5D+1,6.5D+1,5.D+1,1.D+1,
     F 0.D+0,1.D+1,4.D+1,6.5D+1,0.D+0,5.D+0,1.5D+1,2.D+1,5.5D+1,3.D+1/
      DATA D/4.86D+2,6.4D+2,7.58D+2,7.76D+2,4.77D+2,7.07D+2,1.75D+2,
     F 6.19D+2,6.27D+2,6.14D+2,4.75D+2,3.77D+2,5.24D+2,4.68D+2,5.29D+2/
      GOTO (1,2,3,4,5),MODE
    1 N=15
      NILI=0
      NINL=10
      NELI=0
      NENL=0
      DO 6 I=1,15
      X(I)=0.D+0
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=-0.83152859D+4
      XEX(1)=0.81347013D+00
      XEX(2)=0.11327964D+01
      XEX(3)=0.10861185D+01 
      XEX(4)=0.99832982D+00
      XEX(5)=0.10754861D+01
      XEX(6)=0.10688758D+01 
      XEX(7)=0.62781562D+00
      XEX(8)=0.10929981D+01
      XEX(9)=0.91363214D+00
      XEX(10)=0.86191234D+00
      XEX(11)=0.10047312D+01
      XEX(12)=0.87742923D+00
      XEX(13)=0.98671497D+00
      XEX(14)=0.10411268D+01
      XEX(15)=0.11860997D+01
      RETURN
    2 FX=0.D+0
      DO 15 I=1,15
   15 FX=FX-D(I)*X(I)
      RETURN
    3 DO 11 I=1,15
   11 GF(I)=-D(I)
      RETURN
    4 DO 7 I=1,10
      IF (.NOT.INDEX1(I)) GOTO 7
      C=0.D+0
      DO 9 J=1,15
    9 C=C+A(I,J)*X(J)**2
      G(I)=B(I)-C
    7 CONTINUE
      RETURN
    5 DO 10 I=1,10
      IF (.NOT.INDEX2(I)) GOTO 10
      DO 13 J=1,15
   13 GG(I,J)=-2.D+0*A(I,J)*X(J)
   10 CONTINUE
      RETURN
      END
C      
      SUBROUTINE TP386(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION C,A(10,15),B(10),D(15)
      DATA B/3.85D+2,4.7D+2,5.6D+2,5.65D+2,6.45D+2,4.3D+2,4.85D+2,
     F 4.55D+2,3.9D+2,4.6D+2/
      DATA A/1.D+2,9.D+1,7.D+1,2*5.D+1,4.D+1,3.D+1,2.D+1,1.D+1,5.D+0,
     F 2*1.D+2,5.D+1,0.D+0,1.D+1,0.D+0,6.D+1,3.D+1,7.D+1,1.D+1,
     F 2*1.D+1,2*0.D+0,7.D+1,5.D+1,3.D+1,4.D+1,1.D+1,1.D+2,5.D+0,
     F 3.5D+1,5.5D+1,6.5D+1,6.D+1,9.5D+1,9.D+1,2.5D+1,3.5D+1,5.D+0,
     F 1.D+1,2.D+1,2.5D+1,3.5D+1,4.5D+1,5.D+1,0.D+0,4.D+1,2.5D+1,2.D+1,
     F 0.D+0,5.D+0,2*1.D+2,4.5D+1,3.5D+1,3.D+1,2.5D+1,6.5D+1,5.D+0,
     F 2*0.D+0,4.D+1,3.5D+1,0.D+0,1.D+1,5.D+0,1.5D+1,0.D+0,1.D+1,
     F 2.5D+1,3.5D+1,5.D+1,6.D+1,3.5D+1,6.D+1,2.5D+1,1.D+1,3.D+1,3.5D+1,
     F 0.D+0,5.5D+1,2*0.D+0,6.5D+1,2*0.D+0,8.D+1,0.D+0,9.5D+1,
     F 1.D+1,2.5D+1,3.D+1,1.5D+1,5.D+0,4.5D+1,7.D+1,2.D+1,0.D+0,7.D+1,
     F 5.5D+1,2.D+1,6.D+1,0.D+0,7.5D+1,1.5D+1,2.D+1,3.D+1,2.5D+1,2.D+1,
     F 5.D+0,0.D+0,1.D+1,7.5D+1,1.D+2,2.D+1,2.5D+1,3.D+1,0.D+0,1.D+1,
     F 4.5D+1,4.D+1,3.D+1,3.5D+1,7.5D+1,0.D+0,7.D+1,5.D+0,1.5D+1,3.5D+1,
     F 2.D+1,2.5D+1,0.D+0,3.D+1,1.D+1,5.D+0,1.5D+1,6.5D+1,5.D+1,1.D+1,
     F 0.D+0,1.D+1,4.D+1,6.5D+1,0.D+0,5.D+0,1.5D+1,2.D+1,5.5D+1,3.D+1/
      DATA D/4.86D+2,6.4D+2,7.58D+2,7.76D+2,4.77D+2,7.07D+2,1.75D+2,
     F 6.19D+2,6.27D+2,6.14D+2,4.75D+2,3.77D+2,5.24D+2,4.68D+2,5.29D+2/
      GOTO (1,2,3,4,5),MODE 
    1 N=15
      NILI=0
      NINL=11
      NELI=0
      NENL=0
      DO 6 I=1,15
      X(I)=0.D+0
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=-0.81643688D+4
      XEX(1)=0.10042725D+1
      XEX(2)=0.10871174D+1
      XEX(3)=0.11033800D+1
      XEX(4)=0.10307192D+1
      XEX(5)=0.92857958D+0
      XEX(6)=0.12568055D+1
      XEX(7)=0.76058681D+0
      XEX(8)=0.85688931D+0
      XEX(9)=0.10897780D+1
      XEX(10)=0.98119425D+0
      XEX(11)=0.85106387D+0
      XEX(12)=0.96555941D+0
      XEX(13)=0.90644190D+0
      XEX(14)=0.83804049D+0
      XEX(15)=0.80932365D+0
      RETURN
    2 FX=0.D+0
      DO 20 I=1,15
   20 FX=FX-D(I)*X(I)
      RETURN
    3 DO 11 I=1,15
   11 GF(I)=-D(I)
      RETURN
    4 DO 7 I=1,10
      IF (.NOT.INDEX1(I)) GOTO 7
      C=0.D+0
      DO 9 J=1,15
    9 C=C+A(I,J)*X(J)**2
      G(I)=B(I)-C
    7 CONTINUE
      IF (.NOT.INDEX1(11)) GOTO 12
      C=0.D+0
      DO 14 J=1,15
   14 C=C+DBLE(J)*(X(J)-2.D+0)**2
      G(11)=C/2.D+0-7.D+1
   12 RETURN
    5 DO 10 I=1,10
      IF (.NOT.INDEX2(I)) GOTO 10
      DO 13 J=1,15
   13 GG(I,J)=-2.D+0*A(I,J)*X(J)
   10 CONTINUE
      IF (.NOT.INDEX2(11)) GOTO 15
      DO 16 J=1,15
   16 GG(11,J)=DBLE(J)*(X(J)-2.D+0)
   15 RETURN
      END
C      
      SUBROUTINE TP387(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION C,A(10,15),B(10),D(15)
      DATA B/3.85D+2,4.7D+2,5.6D+2,5.65D+2,6.45D+2,4.3D+2,4.85D+2,
     F 4.55D+2,3.9D+2,4.6D+2/
      DATA A/1.D+2,9.D+1,7.D+1,2*5.D+1,4.D+1,3.D+1,2.D+1,1.D+1,5.D+0,
     F 2*1.D+2,5.D+1,0.D+0,1.D+1,0.D+0,6.D+1,3.D+1,7.D+1,1.D+1,
     F 2*1.D+1,2*0.D+0,7.D+1,5.D+1,3.D+1,4.D+1,1.D+1,1.D+2,5.D+0,
     F 3.5D+1,5.5D+1,6.5D+1,6.D+1,9.5D+1,9.D+1,2.5D+1,3.5D+1,5.D+0,
     F 1.D+1,2.D+1,2.5D+1,3.5D+1,4.5D+1,5.D+1,0.D+0,4.D+1,2.5D+1,2.D+1,
     F 0.D+0,5.D+0,2*1.D+2,4.5D+1,3.5D+1,3.D+1,2.5D+1,6.5D+1,5.D+0,
     F 2*0.D+0,4.D+1,3.5D+1,0.D+0,1.D+1,5.D+0,1.5D+1,0.D+0,1.D+1,
     F 2.5D+1,3.5D+1,5.D+1,6.D+1,3.5D+1,6.D+1,2.5D+1,1.D+1,3.D+1,3.5D+1,
     F 0.D+0,5.5D+1,2*0.D+0,6.5D+1,2*0.D+0,8.D+1,0.D+0,9.5D+1,
     F 1.D+1,2.5D+1,3.D+1,1.5D+1,5.D+0,4.5D+1,7.D+1,2.D+1,0.D+0,7.D+1,
     F 5.5D+1,2.D+1,6.D+1,0.D+0,7.5D+1,1.5D+1,2.D+1,3.D+1,2.5D+1,2.D+1,
     F 5.D+0,0.D+0,1.D+1,7.5D+1,1.D+2,2.D+1,2.5D+1,3.D+1,0.D+0,1.D+1,
     F 4.5D+1,4.D+1,3.D+1,3.5D+1,7.5D+1,0.D+0,7.D+1,5.D+0,1.5D+1,3.5D+1,
     F 2.D+1,2.5D+1,0.D+0,3.D+1,1.D+1,5.D+0,1.5D+1,6.5D+1,5.D+1,1.D+1,
     F 0.D+0,1.D+1,4.D+1,6.5D+1,0.D+0,5.D+0,1.5D+1,2.D+1,5.5D+1,3.D+1/
      DATA D/4.86D+2,6.4D+2,7.58D+2,7.76D+2,4.77D+2,7.07D+2,1.75D+2,
     F 6.19D+2,6.27D+2,6.14D+2,4.75D+2,3.77D+2,5.24D+2,4.68D+2,5.29D+2/
      GOTO (1,2,3,4,5),MODE 
    1 N=15
      NILI=0
      NINL=11
      NELI=0
      NENL=0
      DO 6 I=1,15
      X(I)=0.D+0
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=-0.82501417D+4
      XEX(1)=0.10125415D+1
      XEX(2)=0.10158505D+1
      XEX(3)=0.10309039D+1
      XEX(4)=0.99697018D+0
      XEX(5)=0.98528372D+0
      XEX(6)=0.10368532D+1
      XEX(7)=0.99349349D+0
      XEX(8)=0.97201160D+0
      XEX(9)=0.99994095D+0
      XEX(10)=0.99547294D+0
      XEX(11)=0.96953850D+0
      XEX(12)=0.10080569D+1
      XEX(13)=0.98236999D+0
      XEX(14)=0.99057993D+0
      XEX(15)=0.97760168D+0
      RETURN
    2 FX=0.D+0
      DO 20 I=1,15
   20 FX=FX-D(I)*X(I)
      RETURN
    3 DO 11 I=1,15
   11 GF(I)=-D(I)
      RETURN
    4 DO 7 I=1,10
      IF (.NOT.INDEX1(I)) GOTO 7
      C=0.D+0
      DO 9 J=1,15
    9 C=C+A(I,J)*X(J)**2
      G(I)=B(I)-C
    7 CONTINUE
      IF (.NOT.INDEX1(11)) GOTO 12
      C=0.D+0
      DO 14 J=1,15
   14 C=C+DBLE(J)*(X(J)-2.D+0)**2
      G(11)=C/2.D+0-6.1D+1
   12 RETURN
    5 DO 10 I=1,10
      IF (.NOT.INDEX2(I)) GOTO 10
      DO 13 J=1,15
   13 GG(I,J)=-2.D+0*A(I,J)*X(J)
   10 CONTINUE
      IF (.NOT.INDEX2(11)) GOTO 15
      DO 16 J=1,15
   16 GG(11,J)=DBLE(J)*(X(J)-2.D+0)
   15 RETURN
      END
C      
      SUBROUTINE TP388(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J,L
      DOUBLEPRECISION C,A(10,15),A1(4,15),B(15),D(15)
      DATA B/3.85D+2,4.7D+2,5.6D+2,5.65D+2,6.45D+2,4.3D+2,4.85D+2,
     F 4.55D+2,3.9D+2,4.6D+2,0.D+0,7.D+1,3.61D+2,2.65D+2,3.95D+2/
      DATA A/1.D+2,9.D+1,7.D+1,2*5.D+1,4.D+1,3.D+1,2.D+1,1.D+1,5.D+0,
     F 2*1.D+2,5.D+1,0.D+0,1.D+1,0.D+0,6.D+1,3.D+1,7.D+1,1.D+1,
     F 2*1.D+1,2*0.D+0,7.D+1,5.D+1,3.D+1,4.D+1,1.D+1,1.D+2,5.D+0,
     F 3.5D+1,5.5D+1,6.5D+1,6.D+1,9.5D+1,9.D+1,2.5D+1,3.5D+1,5.D+0,
     F 1.D+1,2.D+1,2.5D+1,3.5D+1,4.5D+1,5.D+1,0.D+0,4.D+1,2.5D+1,2.D+1,
     F 0.D+0,5.D+0,2*1.D+2,4.5D+1,3.5D+1,3.D+1,2.5D+1,6.5D+1,5.D+0,
     F 2*0.D+0,4.D+1,3.5D+1,0.D+0,1.D+1,5.D+0,1.5D+1,0.D+0,1.D+1,
     F 2.5D+1,3.5D+1,5.D+1,6.D+1,3.5D+1,6.D+1,2.5D+1,1.D+1,3.D+1,3.5D+1,
     F 0.D+0,5.5D+1,2*0.D+0,6.5D+1,2*0.D+0,8.D+1,0.D+0,9.5D+1,
     F 1.D+1,2.5D+1,3.D+1,1.5D+1,5.D+0,4.5D+1,7.D+1,2.D+1,0.D+0,7.D+1,
     F 5.5D+1,2.D+1,6.D+1,0.D+0,7.5D+1,1.5D+1,2.D+1,3.D+1,2.5D+1,2.D+1,
     F 5.D+0,0.D+0,1.D+1,7.5D+1,1.D+2,2.D+1,2.5D+1,3.D+1,0.D+0,1.D+1,
     F 4.5D+1,4.D+1,3.D+1,3.5D+1,7.5D+1,0.D+0,7.D+1,5.D+0,1.5D+1,3.5D+1,
     F 2.D+1,2.5D+1,0.D+0,3.D+1,1.D+1,5.D+0,1.5D+1,6.5D+1,5.D+1,1.D+1,
     F 0.D+0,1.D+1,4.D+1,6.5D+1,0.D+0,5.D+0,1.5D+1,2.D+1,5.5D+1,3.D+1/
      DATA A1/1.D+0,4.5D+1,5.3D+1,1.2D+1,2.D+0,2.5D+1,7.4D+1,4.3D+1,
     F 3.D+0,3.5D+1,2.6D+1,5.1D+1,4.D+0,8.5D+1,1.7D+1,3.9D+1,5.D+0,
     F 4.D+1,2.5D+1,5.8D+1,6.D+0,7.3D+1,2.5D+1,4.2D+1,7.D+0,1.7D+1,
     F 2.6D+1,6.D+1,8.D+0,5.2D+1,2.4D+1,2.D+1,9.D+0,8.6D+1,8.5D+1,4.D+1,
     F 1.D+1,1.4D+1,3.5D+1,8.D+1,1.5D+1,3.D+1,1.4D+1,7.5D+1,1.6D+1,
     F 5.D+1,2.3D+1,8.5D+1,1.7D+1,4.D+1,3.7D+1,9.5D+1,1.8D+1,7.D+1,
     F 5.6D+1,2.3D+1,1.9D+1,6.D+1,1.D+1,6.7D+1/
      DATA D/4.86D+2,6.4D+2,7.58D+2,7.76D+2,4.77D+2,7.07D+2,1.75D+2,
     F 6.19D+2,6.27D+2,6.14D+2,4.75D+2,3.77D+2,5.24D+2,4.68D+2,5.29D+2/
      GOTO (1,2,3,4,5),MODE 
    1 N=15
      NILI=4
      NINL=11
      NELI=0
      NENL=0
      DO 6 I=1,15
      X(I)=0.D+0
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      DO 19 I=12,15
      DO 19 J=1,15
   19 GG(I-11,J)=-A1(I-11,J)
      LEX=.FALSE.
      NEX=1
      FEX=-0.58210842D+4
      XEX(1)=0.62683876D+0
      XEX(2)=0.14330999D+1
      XEX(3)=0.14625963D+1
      XEX(4)=0.73133338D+0
      XEX(5)=0.78614240D+0
      XEX(6)=0.12048598D+1
      XEX(7)=-0.11433978D+1
      XEX(8)=0.10611103D+1
      XEX(9)=-0.13389293D+0
      XEX(10)=0.11820107D+1
      XEX(11)=0.96917757D+0
      XEX(12)=-0.84501289D+0
      XEX(13)=0.48122454D+0
      XEX(14)=-0.33986164D+0
      XEX(15)=0.68589012D+0
      RETURN
    2 FX=0.D+0
      DO 20 I=1,15
   20 FX=FX-D(I)*X(I)
      RETURN
    3 DO 11 I=1,15
   11 GF(I)=-D(I)
      RETURN
    4 DO 7 I=12,15
      L=I-11
      IF (.NOT.INDEX1(L)) GOTO 7
      C=0.D+0
      DO 9 J=1,15
    9 C=C+A1(L,J)*X(J)
      G(L)=B(I)-C
    7 CONTINUE
      DO 14 I=1,10
      IF (.NOT.INDEX1(I+4)) GOTO 14
      C=0.D+0
      DO 16 J=1,15
   16 C=C+A(I,J)*X(J)**2
      G(I+4)=B(I)-C
   14 CONTINUE
      IF (.NOT.INDEX1(15)) GOTO 17
      C=0.D+0
      DO 18 J=1,15
   18 C=C+DBLE(J)*(X(J)-2.D+0)**2
      G(15)=C/2.D+0-193.121D+0
   17 RETURN
    5 DO 22 I=1,10
      IF (.NOT.INDEX2(I+4)) GOTO 22
      DO 24 J=1,15
   24 GG(I+4,J)=-2.D+0*A(I,J)*X(J)
   22 CONTINUE
      IF (.NOT.INDEX2(15)) GOTO 25
      DO 26 J=1,15
   26 GG(15,J)=DBLE(J)*(X(J)-2.D+0)
   25 RETURN
      END
C      
      SUBROUTINE TP389(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J,L
      DOUBLEPRECISION C,A(10,15),A1(4,15),B(15),D(15)
      DATA B/3.85D+2,4.7D+2,5.6D+2,5.65D+2,6.45D+2,4.3D+2,4.85D+2,
     F 4.55D+2,3.9D+2,4.6D+2,0.D+0,7.D+1,3.61D+2,2.65D+2,3.95D+2/
      DATA A/1.D+2,9.D+1,7.D+1,2*5.D+1,4.D+1,3.D+1,2.D+1,1.D+1,5.D+0,
     F 2*1.D+2,5.D+1,0.D+0,1.D+1,0.D+0,6.D+1,3.D+1,7.D+1,1.D+1,
     F 2*1.D+1,2*0.D+0,7.D+1,5.D+1,3.D+1,4.D+1,1.D+1,1.D+2,5.D+0,
     F 3.5D+1,5.5D+1,6.5D+1,6.D+1,9.5D+1,9.D+1,2.5D+1,3.5D+1,5.D+0,
     F 1.D+1,2.D+1,2.5D+1,3.5D+1,4.5D+1,5.D+1,0.D+0,4.D+1,2.5D+1,2.D+1,
     F 0.D+0,5.D+0,2*1.D+2,4.5D+1,3.5D+1,3.D+1,2.5D+1,6.5D+1,5.D+0,
     F 2*0.D+0,4.D+1,3.5D+1,0.D+0,1.D+1,5.D+0,1.5D+1,0.D+0,1.D+1,
     F 2.5D+1,3.5D+1,5.D+1,6.D+1,3.5D+1,6.D+1,2.5D+1,1.D+1,3.D+1,3.5D+1,
     F 0.D+0,5.5D+1,2*0.D+0,6.5D+1,2*0.D+0,8.D+1,0.D+0,9.5D+1,
     F 1.D+1,2.5D+1,3.D+1,1.5D+1,5.D+0,4.5D+1,7.D+1,2.D+1,0.D+0,7.D+1,
     F 5.5D+1,2.D+1,6.D+1,0.D+0,7.5D+1,1.5D+1,2.D+1,3.D+1,2.5D+1,2.D+1,
     F 5.D+0,0.D+0,1.D+1,7.5D+1,1.D+2,2.D+1,2.5D+1,3.D+1,0.D+0,1.D+1,
     F 4.5D+1,4.D+1,3.D+1,3.5D+1,7.5D+1,0.D+0,7.D+1,5.D+0,1.5D+1,3.5D+1,
     F 2.D+1,2.5D+1,0.D+0,3.D+1,1.D+1,5.D+0,1.5D+1,6.5D+1,5.D+1,1.D+1,
     F 0.D+0,1.D+1,4.D+1,6.5D+1,0.D+0,5.D+0,1.5D+1,2.D+1,5.5D+1,3.D+1/
      DATA A1/1.D+0,4.5D+1,5.3D+1,1.2D+1,2.D+0,2.5D+1,7.4D+1,4.3D+1,
     F 3.D+0,3.5D+1,2.6D+1,5.1D+1,4.D+0,8.5D+1,1.7D+1,3.9D+1,5.D+0,
     F 4.D+1,2.5D+1,5.8D+1,6.D+0,7.3D+1,2.5D+1,4.2D+1,7.D+0,1.7D+1,
     F 2.6D+1,6.D+1,8.D+0,5.2D+1,2.4D+1,2.D+1,9.D+0,8.6D+1,8.5D+1,4.D+1,
     F 1.D+1,1.4D+1,3.5D+1,8.D+1,1.5D+1,3.D+1,1.4D+1,7.5D+1,1.6D+1,
     F 5.D+1,2.3D+1,8.5D+1,1.7D+1,4.D+1,3.7D+1,9.5D+1,1.8D+1,7.D+1,
     F 5.6D+1,2.3D+1,1.9D+1,6.D+1,1.D+1,6.7D+1/
      DATA D/4.86D+2,6.4D+2,7.58D+2,7.76D+2,4.77D+2,7.07D+2,1.75D+2,
     F 6.19D+2,6.27D+2,6.14D+2,4.75D+2,3.77D+2,5.24D+2,4.68D+2,5.29D+2/
      GOTO (1,2,3,4,5),MODE 
    1 N=15
      NILI=4
      NINL=11
      NELI=0
      NENL=0
      DO 6 I=1,15
      X(I)=0.D+0
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      DO 19 I=12,15
      DO 19 J=1,15
   19 GG(I-11,J)=-A1(I-11,J)
      LEX=.FALSE.
      NEX=1
      FEX=-0.58097197D+4
      XEX(1)=0.67105172D+0
      XEX(2)=0.13885400D+1
      XEX(3)=0.14676761D+1
      XEX(4)=0.76023633D+0
      XEX(5)=0.82935674D+0
      XEX(6)=0.11638523D+1
      XEX(7)=-0.12578290D+1
      XEX(8)=0.98193399D+0
      XEX(9)=0.68416463D-1
      XEX(10)=0.11472773D+1
      XEX(11)=0.98662969D+0
      XEX(12)=-0.88834924D+0
      XEX(13)=0.56465631D+0
      XEX(14)=-0.58120082D+0
      XEX(15)=0.72096897D+0
      RETURN
    2 FX=0.D+0
      DO 20 I=1,15
   20 FX=FX-D(I)*X(I)
      RETURN
    3 DO 11 I=1,15
   11 GF(I)=-D(I)
      RETURN
    4 DO 7 I=12,15
      L=I-11
      IF (.NOT.INDEX1(L)) GOTO 7
      C=0.D+0
      DO 9 J=1,15
    9 C=C+A1(L,J)*X(J)
      G(L)=B(I)-C
    7 CONTINUE
      DO 14 I=1,10
      IF (.NOT.INDEX1(I+4)) GOTO 14
      C=0.D+0
      DO 16 J=1,15
   16 C=C+A(I,J)*X(J)**2
      G(I+4)=B(I)-C
   14 CONTINUE
      IF (.NOT.INDEX1(15)) GOTO 17
      C=0.D+0
      DO 18 J=1,15
   18 C=C+DBLE(J)*(X(J)-2.D+0)**2
      G(15)=C/2.D+0-2.D+2
   17 RETURN
    5 DO 22 I=1,10
      IF (.NOT.INDEX2(I+4)) GOTO 22
      DO 24 J=1,15
   24 GG(I+4,J)=-2.D+0*A(I,J)*X(J)
   22 CONTINUE
      IF (.NOT.INDEX2(15)) GOTO 25
      DO 26 J=1,15
   26 GG(15,J)=DBLE(J)*(X(J)-2.D+0)
   25 RETURN
      END
C      
      SUBROUTINE TP390(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION PSI(11),ZI1,ZI2,ZI3
      GOTO (1,2,3,4,3),MODE
    1 N=19
      NILI=1
      NINL=0
      NELI=0
      NENL=11
      X(1)=2.D-2
      X(2)=4.D+0
      X(3)=100.D+0
      X(4)=100.D+0
      X(5)=15.D+0
      X(6)=15.D+0
      X(7)=100.D+0
      X(8)=1000.D+0
      X(9)=1000.D+0
      X(10)=1000.D+0
      X(11)=9000.D+0
      X(12)=0.001D+0
      X(13)=0.001D+0
      X(14)=1.D+0
      X(15)=0.001D+0
      X(16)=0.001D+0
      X(17)=.1D+0
      X(18)=8000.D+0
      X(19)=0.001D+0
      DO 6 I=1,19
      LXU(I)=.TRUE.
      LXL(I)=.TRUE.
      XU(I)=1.D+5
    6 XL(I)=0.00001D+0
      DO 7 I=1,2
      XU(I)=5.D+1
    7 XU(I+15)=5.D+1
      DO 9 I=3,6
    9 XU(I)=1.D+2
      DO 10 I=12,15
   10 XU(I)=1.D+0
      LEX=.FALSE.
      NEX=1
      FEX=0.244724654D+2
      XEX(1)=0.004473667D+0
      XEX(2)=0.3441565D+1
      XEX(3)=0.9934824D+2
      XEX(4)=0.89130035D+2
      XEX(5)=0.15279316D+2
      XEX(6)=0.15279316D+2
      XEX(7)=0.94726127D+2
      XEX(8)=0.12304197D+5
      XEX(9)=0.12313263D+5
      XEX(10)=0.12313263D+5
      XEX(11)=0.95905631D+5
      XEX(12)=0.00001D+0
      XEX(13)=0.00001D+0
      XEX(14)=0.9999890D+0
      XEX(15)=0.00001D+0
      XEX(16)=0.00001D+0
      XEX(17)=0.1622235D+0
      XEX(18)=0.83051515D+4
      XEX(19)=0.0014797356D+0
      RETURN
    2 ZI1=25.D+0*(2268.D+0*X(16)*X(1))**0.827D0
      ZI2=1.75D+05*X(17)+3.65D+04*X(17)**.182D0
      ZI3=12.6D+0*X(18)+5.35D+0*10.D+0**3.378D0/X(18)**.126D0
      FX=1.4D+0*(ZI1+ZI2+ZI3+1.095D+04+1.15D+03*(X(1)*(X(13)-X(14))
     1 +X(2)*(1.D+0+X(12))-3.D+0*(1.D+0-X(19))))
    3 RETURN
    4 IF (INDEX1(1)) G(1)=1.D+0-X(13)-X(14)
      CALL TP390A(X,PSI)
      DO 8 I=2,12
      IF (.NOT.INDEX1(I)) GOTO 8
      G(I)=PSI(I-1)
    8 CONTINUE
      RETURN
      END
C      
      SUBROUTINE TP390A(X,PSI)
      DOUBLEPRECISION X(19),PSI(11),AK,XZ4,ZJ1,YZ4,ZJ2,ZJ3,ZJ4,ZJ5,
     /       QZ12,ZJ8,ZJ10,CK,TEST,ZK7
      AK=.0259D+0*25.D+0/20.D+0**.656D0
      XZ4=X(3)*DEXP(-AK*X(16))
      ZJ1=-(X(1)*X(13)*XZ4+300.D+0*X(19))
      PSI(1)=ZJ1+X(1)*X(3)-X(2)*X(5)*X(12)
      YZ4=X(7)+.5D+0*(X(3)-XZ4)
      ZJ2=-X(13)*X(1)*YZ4
      PSI(2)=ZJ2+X(1)*X(7)-X(2)*X(9)*X(12)
      ZJ3=-300.D+0*(1.D+0-X(19))+3.D+0*X(6)*(1.D+0-X(19))
     F    -X(1)*X(14)*XZ4
      PSI(3)=ZJ3+X(2)*(X(4)-X(6))+X(1)*X(6)*X(14)
      ZJ4=3.D+0*X(11)*(1.D+0-X(19))+X(1)*X(14)*(X(11)-YZ4)
      PSI(4)=ZJ4+X(2)*(X(8)-X(11))
      ZJ5=X(17)*(.48D+0*X(5)*X(9)/(100.D+0+X(5)))
      PSI(5)=-2.D+0*ZJ5+X(2)*(X(4)-X(5))
      PSI(6)=ZJ5+X(2)*(X(8)-X(9))-.048D+0*X(9)*X(17)
      ZK7=X(1)*(1.D+0-X(13)-X(14))
      QZ12=X(1)*(1.D+0-X(13)-X(14))+X(2)*(1.D+0-X(12))
      PSI(7)=-ZK7*XZ4+X(6)*QZ12-X(2)*X(5)*(1.D+0-X(12))
      ZJ8=X(10)*QZ12-ZK7*YZ4
      PSI(8)=ZJ8-X(2)*X(9)*(1.D+0-X(12))
      PSI(9)=6.D+0*(1.D+0-X(15))*(20.D+0-X(6))+X(11)*(X(2)-3.D+0
     F  *(1.D+0-X(15))-X(1)*X(14))+3.D+0*X(19)*X(11)-X(10)*QZ12
      CK=7.4D+0*2.D+0*1.2D+0**4/2.31D+04
      TEST=-CK*X(18)/QZ12
      IF (TEST.GT.99) ZJ10=-2.1D0*DSQRT(DABS(X(10)))*DEXP(99.D+0)
      IF (TEST.LT.99) ZJ10=-2.1D0*DSQRT(DABS(X(10)))
     /                              *DEXP(-CK*X(18)/QZ12)
      PSI(10)=ZJ10+2.D+0*(20.D+0-X(6))
      PSI(11)=(1.D+0-X(13))*X(1)-X(12)*X(2)-3.D+0*X(19)
      RETURN
      END
C      
      SUBROUTINE TP391(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J
      DOUBLEPRECISION SUM,WURZ
      GOTO (1,2,3,3,3),MODE
    1 N=30
      NILI=0
      NINL=0
      NELI=0
      NENL=0
      DO 6 I=1,30
      SUM=0.D+0
      DO 60 J=1,30
      IF (J .EQ. I) GOTO 60
      WURZ=DSQRT(DBLE(I)/DBLE(J))
      SUM=SUM+WURZ*((DSIN(DLOG(WURZ)))**5+(DCOS(DLOG(WURZ)))**5)
   60 CONTINUE
      X(I)=-2.8742711D+0*(DBLE((I-15)**3)+SUM)
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=0.D+0
      XEX(1)=0.64449021D+01 
      XEX(2)=0.51379199D+01
      XEX(3)=0.40135613D+01
      XEX(4)=0.30596619D+01
      XEX(5)=0.22646038D+01
      XEX(6)=0.16151497D+01
      XEX(7)=0.10956969D+01
      XEX(8)=0.69059536D+00
      XEX(9)=0.38552462D+00
      XEX(10)=0.16651219D+00
      XEX(11)=0.19318427D-01
      XEX(12)=-0.70414236D-01
      XEX(13)=-0.11703457D+00
      XEX(14)=-0.13486644D+00
      XEX(15)=-0.13822151D+00
      XEX(16)=-0.14140568D+00
      XEX(17)=-0.15872230D+00 
      XEX(18)=-0.20446919D+00
      XEX(19)=-0.29293764D+00
      XEX(20)=-0.43841769D+00
      XEX(21)=-0.65524697D+00
      XEX(22)=-0.95789935D+00
      XEX(23)=-0.13608416D+01
      XEX(24)=-0.18778557D+01
      XEX(25)=-0.25219169D+01
      XEX(26)=-0.33067694D+01
      XEX(27)=-0.42480867D+01
      XEX(28)=-0.53623746D+01
      XEX(29)=-0.66652372D+01
      XEX(30)=-0.81710133D+01
      RETURN
    2 FX=0.D+0
      DO 7 I=1,30
      SUM=0.D+0
      DO 70 J=1,30
      IF (J .EQ. I) GOTO 70
      WURZ=DSQRT(X(J)**2+DBLE(I)/DBLE(J))
      SUM=SUM+WURZ*((DSIN(DLOG(WURZ)))**5+(DCOS(DLOG(WURZ)))**5)
   70 CONTINUE
    7 FX=FX+(4.2D+2*X(I)+DBLE((I-15)**3)+SUM)**2
    3 RETURN
      END      
C
      SUBROUTINE TP392 (MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L15/ F(200)    
     /           /L16/ DF(200,200)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I,J,K,L
      DOUBLE PRECISION F,DF,SUM1,SUM,R1(3,5),R2(3,5),KA(3,5),K1(3,5),
     F       KP(3,5),K3(3,5),KL1(3,5),KL2(3,5),H(3,5),B(3,5),T(3,3)
      DATA R1/1.D+3,5.2D+2,9.1D+2,1.D+3,5.2D+2,9.1D+2,1.D+3,5.2D+2,
     F   1.D+3,1.1D+3,6.D+2,1.D+3,1.1D+3,6.D+2,1.D+3/
      DATA R2/0.3D+0,0.1D+0,0.2D+0,0.3D+0,0.1D+0,0.2D+0,0.3D+0,0.1D+0,
     F   0.2D+0,0.3D+0,0.1D+0,0.2D+0,0.3D+0,0.1D+0,0.2D+0/
      DATA KA/1.2D+2,6.5D+1,1.05D+2,1.5D+2,6.5D+1,1.05D+2,1.5D+2,8.D+1,
     F   1.2D+2,1.7D+2,8.D+1,1.2D+2,1.7D+2,8.D+1,1.2D+2/
      DATA K1/1.5D+2,7.5D+1,1.4D+2,1.5D+2,7.5D+1,1.4D+2,1.5D+2,7.5D+1,
     F   1.4D+2,1.7D+2,9.D+1,1.5D+2,1.7D+2,9.D+1,1.5D+2/
      DATA KP/1.6D+2,7.5D+1,1.4D+2,1.6D+2,7.5D+1,1.4D+2,1.6D+2,7.5D+1,
     F   1.4D+2,1.8D+2,9.D+1,1.5D+2,1.8D+2,9.D+1,1.5D+2/
      DATA K3/.2D-1,.1D-1,.15D-1,.2D+0,.1D+0,.15D+0,.25D+0,.1D+0,.15D+0,
     F   .25D+0,2*.15D+0,.25D+0,2*.15D+0/
      DATA KL1/3*.5D-2,3*.5D-1,9*.6D-1/
      DATA KL2/8.D+1,4.5D+1,7.5D+1,8.D+1,4.5D+1,7.5D+1,1.D+2,4.5D+1,
     F   9.D+1,1.D+2,5.D+1,9.D+1,1.D+2,5.D+1,9.D+1/
      DATA H/1.D+2,2.8D+2,5.2D+2,1.8D+2,2*4.D+2,2.2D+2,4.5D+2,5.D+2,
     F   1.5D+2,4.5D+2,6.3D+2,1.D+2,4.D+2,6.D+2/
      DATA T/.6D+0,.3D+0,.36D+0,.4D+0,.1D+0,.8D-1,.1D+0,.12D+0,.6D-1/
      DATA B/2*1.7D+2,1.8D+2,2*1.7D+2,1.8D+2,2*1.7D+2,1.8D+2,2*1.7D+2,
     F   1.8D+2,2*1.7D+2,1.8D+2/
      GOTO (1,2,3,4,5),MODE
    1 N=30
      NILI=45
      NINL=0
      NELI=0
      NENL=0
      X(1)=80.D+0
      X(2)=100.D+0
      X(3)=400.D+0
      X(4)=100.D+0
      X(5)=200.D+0
      X(6)=200.D+0
      X(7)=100.D+0
      X(8)=250.D+0
      X(9)=400.D+0
      X(10)=50.D+0
      X(11)=200.D+0
      X(12)=500.D+0
      X(13)=50.D+0
      X(14)=200.D+0
      X(15)=500.D+0
      X(16)=100.D+0
      X(17)=120.D+0
      X(18)=410.D+0
      X(19)=120.D+0
      X(20)=250.D+0
      X(21)=250.D+0
      X(22)=150.D+0
      X(23)=300.D+0
      X(24)=410.D+0
      X(25)=600.D+0
      X(26)=250.D+0
      X(27)=510.D+0
      X(28)=100.D+0
      X(29)=250.D+0
      X(30)=510.D+0
      DO 6 I=1,30
      LXU(I)=.FALSE.
      LXL(I)=.TRUE.
    6 XL(I)=0.D+0
      LEX=.FALSE.
      NEX=1
      FEX=-1.6960671D+6
      XEX(1)=100.00D+0
      XEX(2)=145.00D+0
      XEX(3)=520.00D+0
      XEX(4)=171.43D+0
      XEX(5)=67.86D+0
      XEX(6)=400.00D+0
      XEX(7)=177.37+0
      XEX(8)=4.74D+0
      XEX(9)=500.D+0
      XEX(10)=150.00D+0
      XEX(11)=71.70D+0
      XEX(12)=630.00D+0
      XEX(13)=100.00D+0
      XEX(14)=125.00D+0
      XEX(15)=600.00D+0
      XEX(16)=100.00D+0
      XEX(17)=145.00D+0
      XEX(18)=520.00D+0
      XEX(19)=171.43D+0
      XEX(20)=67.86D+0
      XEX(21)=400.00D+0
      XEX(22)=177.37D+0
      XEX(23)=33.94D+0
      XEX(24)=500.00D+0
      XEX(25)=150.00D+0
      XEX(26)=42.50D+0
      XEX(27)=630.00D+0
      XEX(28)=100.00D+0
      XEX(29)=125.00D+0
      XEX(30)=600.00D+0
      RETURN
    2 FX=0.D+0
      DO 70 I=1,5
      SUM=0.D+0
      DO 71 J=1,3
      SUM1=0.D+0
      DO 72 K=1,I
   72 SUM1=SUM1+X(12+J+3*K)-X(J-3+3*K)
   71 SUM=SUM+X(3*(I-1)+J)*(R1(J,I)-KA(J,I))-X(3*(I-1)+J)**2*R2(J,I)
     F    -X(12+3*I+J)*(K1(J,I)+KP(J,I))-(X(12+3*I+J)-X(J+3*I-3))**2
     F    *(K3(J,I)+KL1(J,I))-KL2(J,I)*SUM1
   70 FX=FX-SUM
    3 RETURN
    4 DO 8 I=1,5
      DO 8 J=1,3
      L=3*(I-1)+J
    8 IF(INDEX1(L)) G(L)=H(J,I)-X(L)
      DO 9 I=1,5
      DO 9 J=1,3
      L=3*(I-1)+J+15
      IF (.NOT.INDEX1(L)) GOTO 9
      G(L)=B(J,I)
      DO 10 K=1,3
   10 G(L)=G(L)-T(J,K)*X(12+3*I+K)
    9 CONTINUE
      DO 11 I=1,5
      DO 11 J=1,3
      L=3*(I-1)+J+30
      IF (.NOT. INDEX1(L)) GOTO 11
      G(L)=0.D+0
      DO 12 K=1,I
   12 G(L)=G(L)+X(12+3*K+J)-X(J-3+3*K)
   11 CONTINUE
    5 RETURN
      END           
C
      SUBROUTINE TP393(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      DOUBLEPRECISION PHI,C,E
      GOTO (1,2,3,4,5),MODE
    1 N=48
      NILI=0
      NINL=1
      NELI=2
      NENL=0
      DO 6 I=1,24
    6 X(I)=1.D+0
      DO 7 I=25,30
    7 X(I)=1.3D+0
      DO 8 I=31,48
    8 X(I)=1.D+0
      DO 9 I=1,48
      LXL(I)=.TRUE.
    9 XL(I)=.002D+0
      DO 10 I=1,24
      LXU(I)=.TRUE.
   10 XU(I)=2.D+0
      DO 11 I=25,48
   11 LXU(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=0.86337998D+0
      XEX(1)=2.D+0
      XEX(2)=.002D+0
      XEX(3)=2.D+0
      XEX(4)=.0339797D+0
      XEX(5)=.01657455D+0
      XEX(6)=2.D+0
      XEX(7)=1.8945347D+0
      XEX(8)=.002D+0
      XEX(9)=2.D+0
      XEX(10)=.03424074D+0
      XEX(11)=.016670308D+0
      XEX(12)=2.D+0
      XEX(13)=2.D+0
      XEX(14)=.002D+0
      XEX(15)=2.D+0
      XEX(16)=.002D+0
      XEX(17)=.002D+0
      XEX(18)=1.988000D+0
      XEX(19)=2.D+0
      XEX(20)=.002D+0
      XEX(21)=2.D+0
      XEX(22)=.002D+0
      XEX(23)=.002D+0
      XEX(24)=2.D+0
      XEX(25)=1.0159886D+0
      XEX(26)=.002D+0
      XEX(27)=1.003163D+0
      XEX(28)=.002D+0
      XEX(29)=.002D+0
      XEX(30)=.999691944D+0
      XEX(31)=1.11272844D+0
      XEX(32)=.002D+0
      XEX(33)=1.1024463D+0
      XEX(34)=.002D+0
      XEX(35)=.002D+0
      XEX(36)=1.1030764D+0
      XEX(37)=.92326572D+0
      XEX(38)=.9343325D+0
      XEX(39)=.92947437D+0
      XEX(40)=.91383802D+0
      XEX(41)=.90517162D+0
      XEX(42)=.89452569D+0
      XEX(43)=1.174573D+0
      XEX(44)=.002D+0
      XEX(45)=1.12080408D+0
      XEX(46)=.002D+0
      XEX(47)=.002D+0
      XEX(48)=1.1163321536D+0
      RETURN
    2 E=0.D+0
      DO 100 I=1,12
      C=1.D+0-X(I)
  100 E=E+10.D+0*C*C
      DO 120 I=25,36
      C=X(I)-1.D+0
  120 E=E+1000.D+0*(.1D+0+2.D+0*C*(C+DSQRT(.1D+0+C*C)))/4.D+0
      DO 140 I=37,42
      C=X(I)-1.D+0
  140 E=E+2000.D+0*(.1D+0+2.D+0*C*(C+DSQRT(.1D+0+C*C)))/4.D+0
      DO 160 I=43,48
  160 E=E+100.D+0*X(I)
      FX=E/1000.D+0
    3 RETURN
    4 IF (.NOT.INDEX1(1)) GOTO 12
      CALL TP393B(X,PHI)
      G(1)=PHI
   12 IF (.NOT.INDEX1(2)) GOTO 14
      G(2)=12.D+0
      DO 13 I=1,12
   13 G(2)=G(2)-X(I)
   14 IF (.NOT.INDEX1(3)) GOTO 5
      G(3)=12.D+0
      DO 15 I=1,12
   15 G(3)=G(3)-X(I+12)
    5 RETURN
      END
C
      SUBROUTINE TP393B(X,PHI)
      DOUBLEPRECISION X,PHI,A,ALP,U,SUM,R
      DIMENSION A(18),U(18),X(48)
      INTEGER I,K1,K2,K3
      DATA (A(I),I=1,18)/.9D+0,.8D+0,1.1D+0,1.D+0,.7D+0,1.1D+0,
     F     1.D+0,1.D+0,1.1D+0,.9D+0,.8D+0,1.2D+0,.9D+0,1.2D+0,
     F     1.2D+0,1.D+0,1.D+0,.9D+0/
C     1ST TIER OF GASFIERS
      DO 20 I=1,6
      K1=I+24
      K2=I+42
      K3=I+12
      ALP=X(K1)*X(K1)*A(I)*2.D+0*X(K2)/(1.D+0+X(K2))*X(K3)
   20 U(I)=X(I)*X(I)/(X(I)+ALP)
C     2ND TIER OF GASFIERS
      DO 40 I=7,12
      K1=I+24
      K2=I+36
      K3=I+12
      ALP=X(K1)*X(K1)*A(I)*2.D+0*X(K2)/(1.D+0+X(K2))*X(K3)
      SUM=X(I)+U(I-6)
   40 U(I)=SUM*SUM/(SUM+ALP)
C     1ST TIER OF METHANATORS
      DO 60 I=13,15
      K1=2*(I-10)+1
      K2=I+24
      ALP=X(K2)*X(K2)*A(I)
      SUM=U(K1)+U(K1+1)
   60 U(I)=SUM*SUM/(SUM+ALP)
C     2ND TIER OF METHANATORS
      DO 80 I=16,18
      K1=I+24
      ALP=X(K1)*X(K1)*A(I)
      SUM=U(I-3)
   80 U(I)=SUM*SUM/(SUM+ALP)
      R=U(16)+U(17)+U(18)
      PHI=1.5D+0-R
      RETURN
      END
C      
      SUBROUTINE TP394(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=20
      NILI=0
      NINL=0
      NELI=0
      NENL=1
      DO 6 I=1,20
      X(I)=2.D+0
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=1.9166667D+0
      XEX(1)=0.91287160D+0
      XEX(2)=0.40824680D+0
      XEX(3)=-0.16746493D-4
      XEX(4)=-0.54074613D-5
      XEX(5)=0.19606096D-5
      XEX(6)=-0.88626385D-5
      XEX(7)=0.81697576D-5
      XEX(8)=-0.14386551D-4
      XEX(9)=0.21831200D-4
      XEX(10)=-0.13873341D-4
      XEX(11)=0.13498048D-4
      XEX(12)=-0.39814429D-5
      XEX(13)=-0.11023953D-4
      XEX(14)=-0.12809830D-4
      XEX(15)=0.79408513D-5
      XEX(16)=0.20458900D-4
      XEX(17)=0.45644559D-5
      XEX(18)=-0.94429887D-5
      XEX(19)=-0.10142804D-4
      XEX(20)=-0.13788343D-5
      RETURN
    2 FX=0.D+0
      DO 8 I=1,20
    8 FX=FX+DBLE(I)*(X(I)**2+X(I)**4)
      RETURN
    3 DO 9 I=1,20
    9 GF(I)=DBLE(I)*(2.D+0*X(I)+4.D+0*X(I)**3)
      RETURN
    4 IF(.NOT.INDEX1(1))GOTO 10
      G(1)=0.D+0
      DO 11 I=1,20
   11 G(1)=G(1)+X(I)**2
      G(1)=G(1)-1.D+0
   10 RETURN
    5 IF(.NOT.INDEX2(1))GOTO 12
      DO 13 I=1,20
   13 GG(1,I)=2.D+0*X(I)
   12 RETURN
      END
C      
      SUBROUTINE TP395(MODE)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)  
      INTEGER    MODE,NMAX,MMAX,N,NILI,NINL,NELI,NENL,NEX   
      DOUBLEPRECISION X,G,GF,GG,FX,XL,XU,FEX,XEX 
      LOGICAL    INDEX1,INDEX2,LXL,LXU,LEX    
      PARAMETER (NMAX=101, MMAX=50)
      COMMON     /L1/  N,NILI,NINL,NELI,NENL 
     /           /L2/  X(NMAX)                              
     /           /L3/  G(MMAX)           
     /           /L4/  GF(NMAX)          
     /           /L5/  GG(MMAX,NMAX)        
     /           /L6/  FX
     /           /L9/  INDEX1(MMAX)     
     /           /L10/ INDEX2(MMAX)
     /           /L11/ LXL(NMAX)       
     /           /L12/ LXU(NMAX)    
     /           /L13/ XL(NMAX)     
     /           /L14/ XU(NMAX)    
     /           /L20/ LEX,NEX,FEX,XEX(NMAX) 
      INTEGER I
      GOTO (1,2,3,4,5),MODE
    1 N=50
      NILI=0
      NINL=0
      NELI=0
      NENL=1
      DO 6 I=1,50
      X(I)=2.D+0
      LXU(I)=.FALSE.
    6 LXL(I)=.FALSE.
      LEX=.FALSE.
      NEX=1
      FEX=0.19166668D+1
      XEX(1)=0.91285206D+0
      XEX(2)=0.40829045D+0
      XEX(3)=-0.64969989D-5
      XEX(4)=-0.99096716D-4
      XEX(5)=0.11891290D-3
      XEX(6)=-0.46486687D-4
      XEX(7)=0.57605078D-4
      XEX(8)=-0.48016383D-4
      XEX(9)=0.25691371D-4
      XEX(10)=0.11670144D-4
      XEX(11)=-0.30881321D-4
      XEX(12)=0.87202482D-5
      XEX(13)=0.19980370D-4
      XEX(14)=-0.12338706D-4
      XEX(15)=-0.16390153D-4
      XEX(16)=0.73383634D-5
      XEX(17)=0.16862980D-4
      XEX(18)=0.43922807D-5
      XEX(19)=-0.58623189D-5
      XEX(20)=-0.25188987D-5
      XEX(21)=0.45980202D-5
      XEX(22)=0.32507205D-5
      XEX(23)=-0.66596023D-5
      XEX(24)=-0.14419491D-4
      XEX(25)=-0.12164937D-4
      XEX(26)=-0.39129061D-5
      XEX(27)=0.98985037D-6
      XEX(28)=0.14776535D-6
      XEX(29)=-0.68312704D-6
      XEX(30)=0.24242977D-5
      XEX(31)=0.53892372D-5
      XEX(32)=0.26662956D-5
      XEX(33)=-0.29282090D-5
      XEX(34)=-0.38338271D-5
      XEX(35)=0.61198364D-6
      XEX(36)=0.43671860D-5
      XEX(37)=0.41104627D-5
      XEX(38)=0.14549012D-5
      XEX(39)=-0.12562117D-5
      XEX(40)=-0.30092086D-5
      XEX(41)=-0.38620459D-5
      XEX(42)=-0.42627256D-5
      XEX(43)=-0.45080325D-5
      XEX(44)=-0.44852099D-5
      XEX(45)=-0.37953194D-5
      XEX(46)=-0.23440318D-5
      XEX(47)=-0.74816106D-6
      XEX(48)=-0.54626804D-7
      XEX(49)=-0.10972677D-5
      XEX(50)=-0.21312770D-5
      RETURN
    2 FX=0.D+0
      DO 8 I=1,50
    8 FX=FX+DBLE(I)*(X(I)**2+X(I)**4)
      RETURN
    3 DO 9 I=1,50
      GF(I)=DBLE(I)*(2.D+0*X(I)+4.D+0*X(I)**3)
    9 GF(I)=GF(I)
      RETURN
    4 IF(.NOT.INDEX1(1))GOTO 10
      G(1)=0.D+0
      DO 11 I=1,50
   11 G(1)=G(1)+X(I)**2
      G(1)=G(1)-1.D+0
   10 RETURN
    5 IF(.NOT.INDEX2(1))GOTO 12
      DO 13 I=1,50
   13 GG(1,I)=2.D+0*X(I)
   12 RETURN
      END
C            
      DOUBLE PRECISION FUNCTION GLEICH(P)   
      DOUBLEPRECISION F,P,EPS,EPS2,Y,A         
      EPS=1.D-5
      Y=P+1.D+0                             
    2 F=Y-P-DATAN(1.D+0/Y)                  
      IF (DABS(F).LE.EPS) GOTO 1            
      A=Y*Y+1.D+0                           
      A=(A+1.D+0)/A                         
      Y=Y-F/A                               
      GOTO 2                                
    1 CONTINUE
      EPS2=EPS**2
      IF (Y.GT.EPS2) THEN                              
         GLEICH=Y                              
      ELSE
         GLEICH=EPS2
      ENDIF      
      RETURN                                
      END 
C      
      SUBROUTINE MDNORD(A,B)
      DOUBLEPRECISION A,B,NORINT
      B=NORINT(A)
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION NORINT(X)
C
C   COMPUTES THE GAUSSIAN NORMAL DISTRIBUTION INTEGRAL
C   PRECISION ABOUT 16 DIGITS
C
      IMPLICIT NONE
      DOUBLEPRECISION X
      DOUBLEPRECISION P1(0:8),Q1(0:9),P2(0:5),Q2(0:6)
      DOUBLEPRECISION SQRT2,RSQRTPI,ARG,ARG2,XABS,ERF,ERFC
      DATA P1/.37235079815548067D4, .71136632469540499D4,
     F        .67582169641104859D4, .40322670108300497D4,
     F        .16317602687537147D4, .45626145870609263D3,
     F        .86082762211948595D2, .10064858974909542D2,
     F        .56418958676181361D0/
      DATA Q1/.37235079815548065D4, .11315192081854405D5,
     F        .15802535999402043D5, .13349346561284457D5,
     F        .75424795101934758D4, .29680049014823087D4,
     F        .81762238630454408D3, .15307771075036222D3,
     F        .17839498439139557D2, .1D1/
      DATA P2/.29788656263939929D1, .74097406059647418D1,
     F        .61602098531096305D1, .50190497267842675D1,
     F        .12753666447299660D1, .56418958354775507D0/
      DATA Q2/.33690752069827528D1, .96089653271927879D1,
     F        .17081440747466004D2, .12048951927855129D2,
     F        .93960340162350542D1, .22605285207673270D1, .1D1/
      DATA SQRT2/1.41421356237390505D0/
      DATA RSQRTPI/.56418958354775629D0/
      XABS=DABS(X)
      IF (XABS.GT.0.5D0) THEN
      IF (XABS.GT.8.D0) THEN
      IF (XABS.GT.100.D0) THEN
         ERFC=0.D0
      ELSE
         ARG=XABS/SQRT2
         ERFC=(((((P2(5)*ARG+P2(4))*ARG+P2(3))*ARG+P2(2))*ARG+P2(1)
     F             )*ARG+P2(0))/
     F           ((((((ARG+Q2(5))*ARG+Q2(4))*ARG+Q2(3))*ARG+Q2(2))*ARG
     F                 +Q2(1))*ARG+Q2(0))*DEXP(-ARG**2)
      ENDIF
      ELSE
        ARG=XABS/SQRT2
        ERFC=((((((((P1(8)*ARG+P1(7))*ARG+P1(6))*ARG+P1(5))*ARG
     F                 +P1(4))*ARG+P1(3))*ARG+P1(2))*ARG+P1(1))*ARG
     F                 +P1(0))/
     F         (((((((((ARG+Q1(8))*ARG+Q1(7))*ARG+Q1(6))*ARG+Q1(5))
     F                 *ARG+Q1(4))*ARG+Q1(3))*ARG+Q1(2))*ARG+Q1(1))
     F                 *ARG+Q1(0))*DEXP(-ARG**2)
      ENDIF
      IF (X.LT.0.D0) THEN
         NORINT=ERFC*.5D0
         RETURN
      ELSE
         NORINT=(2.D0-ERFC)*.5D0
         RETURN
      ENDIF
      ELSE
         ARG=XABS/SQRT2
         ARG2=ARG**2
         ERF=ARG*2.D0*RSQRTPI*((((((((((ARG2/210.D0-1.D0/19.D0)*ARG2
     F     /9.D0+1.D0/17.D0)*ARG2/8.D0-1.D0/15.D0)*ARG2/7.D0+1.D0/13.D0)
     F     *ARG2/6.D0-1.D0/11.D0)*ARG2/5.D0+1.D0/9.D0)*ARG2/4.D0
     F     -1.D0/7.D0)*ARG2/3.D0+1.D0/5.D0)*ARG2/2.D0-1.D0/3.D0)*ARG2
     F     +1.D0)
         IF ( X .GE. 0.D0 ) THEN
            NORINT=(1.D0+ERF)*.5D0
         RETURN
      ELSE
         NORINT=(1.D0-ERF)*.5D0
         RETURN
      ENDIF
      ENDIF
      END
C     
