! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE date_conversions_mod
!-----------------------------------------------------------------------
!
! Purpose and description:
!
! These conversion algorithms are based on the algorithm published
! in a letter to the editor of Communications of the ACM (CACM, volume 1
! number 10, October 1968, p.657) by Henry F. Fliegel and
! Thomas Van Flandern
!
! This algorithm is valid only for dates from
! 1/3/-4900 G onward when converting from a Julian day number to a date,
! or from 1/3/-4800 when converting from a date to a Julian day number.
! It should be noted that these algorithms are valid only in the
! Gregorian Calendar and the Proleptic Gregorian Calendar (after the
! dates given above). They do not handle dates in the Julian Calendar.
!
!-----------------------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
!
IMPLICIT NONE

CONTAINS

LOGICAL FUNCTION isaleap(iy)
!
! Returns .TRUE. if IY is a Leap year
! Returns .FALSE. if IY is not a Leap year
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: iy

IF (iy/4*4  /=  iy) THEN    ! Divide by 4
  isaleap=.FALSE.
ELSE
  IF (iy/400*400  ==  iy) THEN  ! Century check
    isaleap=.TRUE.
  ELSE
    IF (iy/100*100  ==  iy) THEN   ! Century qualifier
      isaleap=.FALSE.
    ELSE
      isaleap=.TRUE.
    END IF
  END IF
END IF
END FUNCTION isaleap

! -------------------------------------------------

SUBROUTINE date13 (icd, id, im, iny)
!
! DAY, MONTH, YEAR FROM DAYS SINCE 1.1.1900
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: icd
INTEGER, INTENT(OUT) :: id, im, iny

! LOCAL VARIABLES
INTEGER :: idy, iy
INTEGER :: k,kd,ke,ky,i,k1x
INTEGER :: months(12)
INTEGER :: days_in_feb

k=icd
ke = 0
IF (k  >=  366) THEN  ! these allow for the non-leap years 1900
  k = k + 1
  IF (k  >=  73416) THEN         !2100, ...
    k = k + 1
    IF (k  >=  109941) THEN     !2200,
      k = k + 1
      IF (k  >=  146466) THEN  !2300 ...
        k = k + 1
      END IF
    END IF
  END IF
END IF
IF (k  <=  -36159) THEN   ! and 1800 respectively
  k = k - 1
END IF

ky = k/1461*4
kd = k - k/1461*1461
IF (kd  <   0) THEN
  kd = kd + 1461
  ky = ky - 4
END IF
ky = ky + 1900
IF (kd  >   366) THEN
  kd = kd - 1
  ke = kd/365
  kd = kd - kd/365*365
END IF
IF (kd  ==  0) THEN
  ke = ke - 1
  kd = 365
END IF
iny = ky + ke
idy = kd
iy  = iny

IF (isaleap(iy)) THEN
  days_in_feb = 29
ELSE
  days_in_feb = 28
END IF

months = (/31,days_in_feb,31,30,31,30,31,31,30,31,30,31/)

k1x = idy

DO i=1,12
  k1x = k1x - months(i)
  IF (k1x  >   0) THEN
    CYCLE
  ELSE
    id = k1x + months(i)
    im = i
    iny = iy
  END IF
  EXIT
END DO

END SUBROUTINE date13

! -------------------------------------------------

SUBROUTINE date31 (id, im, iy, icd)
!
! DAYS SINCE 1.1.1900 FROM DAY, MONTH, YEAR
!
IMPLICIT NONE

INTEGER, INTENT(IN) :: id, im, iy
INTEGER, INTENT(OUT) :: icd

! LOCAL VARIABLES
INTEGER :: idy, iny
INTEGER :: k,iyn, days_in_feb
INTEGER :: months(12)

IF (isaleap(iy)) THEN
  days_in_feb = 29
ELSE
  days_in_feb = 28
END IF

months = (/31,days_in_feb,31,30,31,30,31,31,30,31,30,31/)

k = SUM(months(1:(im-1)))    ! use array sections and intrinsics

idy = k + id
iny = iy
iyn = iny - 1900
IF (iyn  >   0) THEN
  icd = idy + iyn*365 + (iyn-1)/4 - (iyn-1)/100 + (iyn+299)/400
ELSE
  icd = idy + iyn*365 + iyn/4 - iyn/100
END IF

END SUBROUTINE date31

END MODULE date_conversions_mod
