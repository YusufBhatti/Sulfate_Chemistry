! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE setperlen_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SETPERLEN_MOD'

CONTAINS
!
!    SUBROUTINE SETPERLEN---------------------------------------------
!
!    Purpose:
!             Return length of current meaning period using mean-level
!             (0, 1, 2 or 3) & current date. Days_in_month is declared
!             and provided by module cdaydata_mod.
!
!    Method: where this routine is called at the end of a month, the
!            month will already have been incremented, which is why
!            i_month-1 is used in setting the period length. Every
!            100 years is not leap unless year is divisible by 400.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc
!
!  Code description:
!   FORTRAN 77 + common extensions also in Fortran 90
!
!
!   Programming standard :  UMDP 3 Version 7  (rev. 6/10/94)
!
!   Logical components covered :
!
!   Project task :
!
!   External documentation:
!
!  -----------------------------------------------------------------
!   Arguments:------------------------------------------------------
SUBROUTINE setperlen (meanlev,i_month,i_year,periodlen)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE cdaydata_mod, ONLY: days_per_4c, days_per_c, days_per_4y,     &
                        days_per_y, days_in_month, days_to_month
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
IMPLICIT NONE


INTEGER ::                                                        &
  meanlev,                                                        &
                   ! IN - Mean level indicator, e.g. 0, 1, 2 or 3
  i_month,                                                        &
                   ! IN - model time (month)
  i_year,                                                         &
                   ! IN - model time (year)
  periodlen        ! OUT - length of current meaning period (days)

!-----------------------------------------------------------------------
! Workspace usage:------------------------------------------------------
! NONE
!-----------------------------------------------------------------------
! External subroutines called:------------------------------------------
! NONE
! ----------------------------------------------------------------------
! Define local variables:-----------------------------------------------

LOGICAL :: l_leap     ! Leap year indicator

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SETPERLEN'

!-----------------------------------------------------------------------
! End of standard header info

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (MOD(i_year,4)  ==  0 .AND.                                    &
                                        ! is this a leap year?
    (MOD(i_year,400)  ==  0 .OR. MOD(i_year,100)  /=  0)) THEN
  l_leap = .TRUE.
ELSE
  l_leap = .FALSE.
END IF

IF (meanlev  ==  0) THEN       ! instantaneous data (e.g. part-
  periodlen = 1                ! way through a monthly mean)
ELSE IF (meanlev  ==  1) THEN   ! end of monthly mean
  IF (l_leap .AND. (i_month  ==  3)) THEN  ! Is it leap year Feb?
    periodlen = days_in_month(i_month-1) + 1
  ELSE IF (i_month  ==  1) THEN  ! end of Dec, so can't use
    periodlen = 31               ! days_in_month(i_month-1)
  ELSE
    periodlen = days_in_month(i_month-1)
  END IF
ELSE IF (meanlev  ==  2) THEN  ! seasonal mean

  !      find season length using current month as pointer
  IF (l_leap) THEN  ! do leap year seasons
    IF (i_month  ==  5) THEN  ! season=FebMarApr
      periodlen = 90
    ELSE IF ((i_month  ==  3) .OR. (i_month  ==  4) .OR.          &
            (i_month  ==  7) .OR. (i_month  ==  12)) THEN
      periodlen = 91          ! for DJF, JFM, AMJ or SON
    ELSE
      periodlen = 92
    END IF
  ELSE              ! do non-leap year seasons
    IF (i_month  ==  5) THEN  ! season=FebMarApr
      periodlen = 89
    ELSE IF ((i_month  ==  3) .OR. (i_month  ==  4)) THEN
      periodlen = 90  ! for DJF and JFM
    ELSE IF ((i_month  ==  7) .OR. (i_month  ==  12)) THEN
      periodlen = 91  ! for AMJ and SON
    ELSE
      periodlen = 92  ! for all other seasons
    END IF
  END IF  ! end of IF test of L_LEAP, and end of seasons.
ELSE IF (meanlev  ==  3) THEN  ! annual mean

  ! Bear in mind period 3 may be 366 days if _previous_ year was leap, and
  ! may not always be 366 days even if current year is a leap year, since
  ! annual means are often not for calendar years
  IF (l_leap .AND. (i_month  >=  3)) THEN
    periodlen = 366
  ELSE IF (MOD(i_year-1,4)  ==  0 .AND.                           &
                                             ! was last year leap?
    (MOD(i_year-1,400)  ==  0 .OR. MOD(i_year-1,100)  /=  0)) THEN
    IF (i_month  <=  2) THEN   !is it no later than end of Jan?
      periodlen = 366          !Then this year had an extra day!
    ELSE
      periodlen = 365
    END IF                 !end of Jan
  ELSE
    periodlen = 365       ! Normal year
  END IF                   !was last year a leap?
ELSE                      ! meanlev has unexpected value
  periodlen = 1           ! so set weighting factor=1
  WRITE(umMessage,*)'SETPERLEN: MEANLEV not in allowed range of 0 to 3'
  CALL umPrint(umMessage,src='setperlen')
END IF  ! end of IF tests on meanlev

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE setperlen
END MODULE setperlen_mod
