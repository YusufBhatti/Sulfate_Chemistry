! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc
MODULE cdaydata_mod

IMPLICIT NONE
! CDATDATA start
!
! Constants needed by routines to calculate day/month/year from
! incremental day number since calendar zero point, and vice-versa,
! when using Gregorian calendar

INTEGER,PARAMETER:: days_per_4c = 146097
INTEGER,PARAMETER:: days_per_c  = 36524
INTEGER,PARAMETER:: days_per_4y = 1461
INTEGER,PARAMETER:: days_per_y  = 365

INTEGER,PARAMETER:: days_in_month(12) =                           &
  (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

INTEGER,PARAMETER:: days_to_month(12) =                           &
  (/0, 31, 59, 90,120,151,181,212,243,273,304,334/)

! CDAYDATA end

END MODULE cdaydata_mod
