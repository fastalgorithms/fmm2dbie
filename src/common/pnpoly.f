cccc
c     Original Version: Randolph Franklin
c     Modifications (minor): Travis Askham
cccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Copyright notice from original version
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Copyright (c) 1970-2003, Wm. Randolph Franklin
      
c     Permission is hereby granted, free of charge, to any person 
c     obtaining a copy of this software and associated documentation 
c     files (the "Software"), to deal in the Software without restriction, 
c     including without limitation the rights to use, copy, modify, merge, 
c     publish, distribute, sublicense, and/or sell copies of the Software, 
c     and to permit persons to whom the Software is furnished to do so, 
c     subject to the following conditions:

c     Redistributions of source code must retain the above copyright 
c     notice, this list of conditions and the following disclaimers.

c     Redistributions in binary form must reproduce the above copyright 
c     notice in the documentation and/or other materials provided with 
c     the distribution.

c     The name of W. Randolph Franklin may not be used to endorse or 
c     promote products derived from this Software without specific prior
c     written permission.

c     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
c     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
c     OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
c     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
c     BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
c     ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
c     CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
c     SOFTWARE.


C>>>PNP1                                                                
C                                                                       
C     ..................................................................
C                                                                       
C        SUBROUTINE PNPOLY                                              
C                                                                       
C        PURPOSE                                                        
C           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON            
C                                                                       
C        USAGE                                                          
C           CALL PNPOLY1 (PX, PY, XX, YY, N, INOUT, IER )                     
C                                                                       
C        DESCRIPTION OF THE PARAMETERS                                  
C           PX      - X-COORDINATE OF POINT IN QUESTION.                
C           PY      - Y-COORDINATE OF POINT IN QUESTION.                
C           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         
C                     VERTICES OF POLYGON.                              
C           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           
C     VERTICES OF POLYGON.                              
C     N       - NUMBER OF VERTICES IN THE POLYGON.                
C     INOUT   - THE SIGNAL RETURNED:                              
C     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        
C     0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,     
C     1 IF THE POINT IS INSIDE OF THE POLYGON.
c     IER     - error flag: 0 means regular operation, 1
c               means the routine wasn't able to allocate
c               two real arrays of size N
C     
C     REMARKS                       
C     THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.      
C     THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY           
C     OPTIONALLY BE INCREASED BY 1.                               
C     THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING      
C     OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX    
C     OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING   
C     N, THESE FIRST VERTICES MUST BE COUNTED TWICE.              
C     INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.         
C     THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM      
C     WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.   
C     
C     SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED  
C     NONE                                                        
C     
C     METHOD             
C     A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT  
C     CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE        
C     POINT IS INSIDE OF THE POLYGON.                             
C     
C     ......................................................
C     
      SUBROUTINE PNPOLY(PX,PY,XX,YY,N,INOUT,IER) 
      IMPLICIT NONE
      REAL*8 XX(*),YY(*), PX, PY
      REAL*8, ALLOCATABLE :: X(:), Y(:)
      LOGICAL MX,MY,NX,NY            
      INTEGER :: N, INOUT, IER, NLOC, IER1
      PARAMETER (NLOC=200)
      REAL *8 :: TX(NLOC), TY(NLOC)

      IER=0
      if (N .gt. NLOC) then
         allocate(X(N),Y(N),stat=IER1)
         if (IER1 .ne. 0) then
            IER=1
            return
         endif
         call PNPOLY0(PX,PY,XX,YY,X,Y,N,INOUT)
      else
         call PNPOLY0(PX,PY,XX,YY,TX,TY,N,INOUT) 
      endif
      
      RETURN
      END
      
      SUBROUTINE PNPOLY0(PX,PY,XX,YY,X,Y,N,INOUT)      
      IMPLICIT NONE
      REAL *8 :: XX(*),YY(*), PX, PY, X(*), Y(*)
      LOGICAL :: MX,MY,NX,NY 
      INTEGER :: I, J, N, INOUT
      
      DO I=1,N 
         X(I)=XX(I)-PX
         Y(I)=YY(I)-PY         
      ENDDO

      INOUT=1 

      DO 2 I=1,N   
         J=1+MOD(I,N)      
         MX=X(I).GE.0.0
         NX=X(J).GE.0.0 
         MY=Y(I).GE.0.0 
         NY=Y(J).GE.0.0 
         IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX))
     1        GO TO 2       
         IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX)))
     1        GO TO 3  
         INOUT=-INOUT 
         GO TO 2  
 3       IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5 
 4       INOUT=0 
 5       INOUT=-INOUT  
 2    CONTINUE
      
      INOUT= -INOUT
      RETURN
      END   

      SUBROUTINE ZPNPOLY(PT,ZZ,N,INOUT,IER) 
      IMPLICIT NONE
      REAL*8 ZZ(2,*), PT(2)
      REAL*8, ALLOCATABLE :: Z(:,:)
      INTEGER N, INOUT, IER, NLOC, IER1
      PARAMETER (NLOC=200)
      REAL *8 :: TZ(2,NLOC)

      IER=0
      if (N .gt. NLOC) then
         allocate(Z(2,N),stat=IER1)
         if (IER1 .ne. 0) then
            IER=1
            return
         endif
         call ZPNPOLY0(PT,ZZ,Z,N,INOUT)
      else
         call ZPNPOLY0(PT,ZZ,TZ,N,INOUT)         
      endif
      
      RETURN
      END
      
      SUBROUTINE ZPNPOLY0(PT,ZZ,Z,N,INOUT)      
      IMPLICIT NONE
      REAL *8 :: ZZ(2,*), PT(2), Z(2,*)
      LOGICAL :: MX,MY,NX,NY 
      INTEGER :: I, J, N, INOUT
      
      DO I=1,N 
         Z(1,I)=ZZ(1,I)-PT(1)
         Z(2,I)=ZZ(2,I)-PT(2)      
      ENDDO

      INOUT=1 

      DO 2 I=1,N   
         J=1+MOD(I,N)      
         MX=Z(1,I).GE.0.0
         NX=Z(1,J).GE.0.0 
         MY=Z(2,I).GE.0.0 
         NY=Z(2,J).GE.0.0 
         IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX))
     1        GO TO 2       
         IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX)))
     1        GO TO 3  
         INOUT=-INOUT 
         GO TO 2  
 3       IF((Z(2,I)*Z(1,J)-Z(1,I)*Z(2,J))/(Z(1,J)-Z(1,I))) 2,4,5 
 4       INOUT=0 
 5       INOUT=-INOUT  
 2    CONTINUE
      
      INOUT= -INOUT
      RETURN
      END   
      
