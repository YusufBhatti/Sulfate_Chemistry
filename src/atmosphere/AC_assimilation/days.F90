! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!    Purpose : Read from ACOBS Files,reformat and place OBS header
!              details in COMOBS. The bulk of the required OBS data
!              is put into dynamic work array OBS for transmission via
!              argument list to GETOBS. OBS is written out to a cache
!              file for subsequent reading at later timesteps.
!              Thus reread of ACOBS files only required intermittently
!              (The routine DAYS does a dd/mm/yy to dayno)
!
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Logical components covered:
!
!    Project Task : P3
!
!    External documentation:
!
!
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: AC Assimilation
MODULE days_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DAYS_MOD'

CONTAINS

SUBROUTINE days (kday,kmonth,kyear,kdays)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER :: kday
INTEGER :: kmonth
INTEGER :: kyear
INTEGER :: kdays
INTEGER :: ileap

!            FIND ELAPSED DAYS SINCE 1 JAN 1980.
INTEGER :: idays(12,2)


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DAYS'

DATA idays/0,31,59,90,120,151,181,212,243,273,304,334,            &
           0,31,60,91,121,152,182,213,244,274,305,335/

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ileap=1
IF (MOD(kyear,4) == 0)ileap=2
kdays=(kyear-1980)*365+idays(kmonth,ileap)+kday+(kyear-1977)/4

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE days
END MODULE days_mod
