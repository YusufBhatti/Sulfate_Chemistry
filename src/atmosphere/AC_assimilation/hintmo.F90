! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE HINTMO -------------------------------------------------
!
!    Purpose :
!
!       Performs bilinear horizontal interpolation of a field
!       on the model grid to observation locations.
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!
!
!
!    Arguments:---------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE hintmo_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HINTMO_MOD'

CONTAINS

SUBROUTINE hintmo(field,cf1pt,cf2pt,cf3pt,cf4pt,                  &
                  np1pt,np2pt,np3pt,np4pt,                        &
                  lenfld,nlev,lenobt,back,                        &
                  icode,cmessage)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Atmos_Max_Sizes
USE UM_ParParams
USE comobs_mod, ONLY: nobtypmx
USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER :: lenfld,nlev,lenobt
REAL :: field(lenfld*nlev),cf1pt(lenobt),cf2pt(lenobt),           &
        cf3pt(lenobt),cf4pt(lenobt),back(lenobt)
INTEGER :: np1pt(lenobt),np2pt(lenobt),np3pt(lenobt),np4pt(lenobt)
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
!
!-INTENT=IN--------------------------------------------------------
!
!     FIELD        - model field to be interpolated
!     LENFLD       - length of one level of field
!     NLEV         - no of levels in field
!     LENOBT       - no of observations
!     CF1,2,3,4PT  - interpolation coefficients for 4 field points
!     NP1,2,3,4PT  - pointers to 4 field points for interpolation
!
!-INTENT=OUT-----------------------------------------------------
!
!     BACK            -  values of field interpolated to obs locations
!     ICODE,CMESSAGE  -  error code and message
!
!----------------------------------------------------------------------
!     Workspace usage
!-----------------------------------------------------------------------
!     NONE
!
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!     Define local variables
!----------------------------------------------------------------------
INTEGER :: job

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='HINTMO'

!
!     JOB   -  Loop counter in loop over obs
!-----------------------------------------------------------------------
!
! --- 1.     HORIZONTAL INTERPOLATION OF MODEL
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (ltimer_ac) CALL timer('HINTMO  ',3)
!
!     LOOP OVER OBSERVATIONS
!     ORDER OF POINTS IS NEAREST POINT (I,J), (I+/-1,J)
!                                      (I,J+/-1), (I+/-1,J+/-1)
DO job=1,lenobt
  !
  back(job) = cf1pt(job) *  field ( np1pt(job) ) +          &
      cf2pt(job) *  field ( np2pt(job) ) +                  &
      cf3pt(job) *  field ( np3pt(job) ) +                  &
      cf4pt(job) *  field ( np4pt(job) )

END DO
!
IF (ltimer_ac) CALL timer('HINTMO  ',4)
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE hintmo
END MODULE hintmo_mod
