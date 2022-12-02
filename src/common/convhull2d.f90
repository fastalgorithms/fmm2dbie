
! CONVEX HULL FINDING ROUTINE IN TWO DIMENSIONS
!
! This routine was accessed on 8/11/2022 from:
! https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
!
! It is an implementation of Andrew's monotone chain algorithm
!
! User-callable routines:
!
! -convhull2d(n,pts,nh,hull)
!   Order(n) routine for finding the convex hull (pts array must be
!   sorted on the x coordinate on input)


FUNCTION Conv2dCross(v1,v2,v3)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------
! INPUT VARIABLES
REAL *8,INTENT(IN) :: v1(2)    !< input vector 1
REAL *8,INTENT(IN) :: v2(2)    !< input vector 2
REAL *8,INTENT(IN) :: v3(2)    !< input vector 3
!-----------------------------------------------
! OUTPUT VARIABLES
REAL *8            :: Conv2dCross    !< cross product
!-----------------------------------------------
! LOCAL VARIABLES
!===============================================
Conv2dCross=(v2(1)-v1(1))*(v3(2)-v1(2))-(v2(2)-v1(2))*(v3(1)-v1(1))
END FUNCTION Conv2dCross

SUBROUTINE ConvHull2d(nPoints,Points,nHull,Hull)
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE 
!------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: nPoints
REAL *8,INTENT(IN)     :: Points(2,0:nPoints-1)
!------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT) :: nHull
! NOTE: allocate Hull always one point greater than Points, because we save the first value twice
REAL *8,INTENT(OUT)    :: Hull(2,0:nPoints)
!------------------------------------------------
! LOCAL VARIABLES
REAL *8                :: Lower(2,0:nPoints-1)
REAL *8                :: Upper(2,0:nPoints-1)
REAL *8                :: CONV2DCROSS
INTEGER             :: i,iLower,iUpper
!================================================
IF(nPoints.LE.1)THEN
  Hull  = Points
  nHull = nPoints
ELSE
  iLower = 0
  Lower  = -HUGE(1.)
  DO i=0,nPoints-1
    DO WHILE(iLower.GE.2.AND.Conv2dCross(Lower(:,iLower-2),Lower(:,iLower-1),Points(:,i)).LE.0.)
      Lower(:,iLower) = -HUGE(1.)
      iLower          = iLower - 1
    END DO
    Lower(:,iLower) = Points(:,i)
    iLower = iLower + 1
  END DO

  iUpper = 0
  Upper  = HUGE(1.)
  DO i=nPoints-1,0,-1
    DO WHILE(iUpper.GE.2.AND.Conv2dCross(Upper(:,iUpper-2),Upper(:,iUpper-1),Points(:,i)).LE.0.)
      Upper(:,iUpper) = HUGE(1.)
      iUpper          = iUpper - 1
    END DO
    Upper(:,iUpper) = Points(:,i)
    iUpper = iUpper + 1
  END DO

  iLower = iLower-1
  iUpper = iUpper-1
  nHull  = iLower+iUpper+1
  
  ! NOTE: Initialize Hull with zeros
  Hull   = 0.

  ! NOTE: save values in Hull
  Hull(:,0     :iLower       -1) = Lower(:,0:iLower-1)
  Hull(:,iLower:iLower+iUpper-1) = Upper(:,0:iUpper-1)

  ! NOTE: save first value twice
  Hull(:,       iLower+iUpper  ) = Hull (:,0         )
END IF

END SUBROUTINE ConvHull2d
