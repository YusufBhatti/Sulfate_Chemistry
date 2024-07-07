! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE HMRTORH ------------------------------------------------
!
!    Purpose : Convert Humidity Mixing Ratio to Relative Humidity
!              and vice versa.  Option to re-initialise cloud water
!              from rh and temp.
!
!              KRMODE = 1 - For HMR to RH
!              KRMODE = 2 - For RH  to HMR
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!
!
!     Arguments
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE hmrtorh_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HMRTORH_MOD'

CONTAINS

SUBROUTINE hmrtorh (krmode, exner, pressure, theta, rh,           &
                    p_field, icode, cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Atmos_Max_Sizes
USE UM_ParParams
USE comobs_mod, ONLY: nobtypmx
USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength
USE nlsizes_namelist_mod, ONLY: model_levels

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_new         => qsat,                           &
                    l_new_qsat_acassim !Currently defaults to FALSE

IMPLICIT NONE
!-----------------------------------------------------------------------
INTEGER :: krmode
INTEGER :: p_field
REAL :: exner (p_field,model_levels)
REAL :: pressure (p_field,model_levels)
REAL :: theta (p_field,model_levels)
REAL :: rh    (p_field,model_levels)

INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
!=======================================================================
!     Dynamic allocation
REAL :: temp     (p_field)
REAL :: smr      (p_field)
!=======================================================================
!     Local variables
INTEGER :: jlev,j

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='HMRTORH'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO jlev=1,model_levels

  !       Convert Theta to Temperature
  !       ----------------------------
  DO j=1,p_field
    temp(j) = theta(j,jlev) * exner(j,jlev)
  END DO


  !       Obtain Saturated Mixing Ratio
  !       -----------------------------
  IF ( l_new_qsat_acassim ) THEN
    CALL qsat_new(smr,temp,pressure(:,jlev),p_field)
  ELSE
    ! DEPENDS ON: qsat
    CALL qsat (smr,temp,pressure(1,jlev),p_field)
  END IF

  IF (krmode == 1) THEN

    !         Convert Mixing Ratio to Relative Humidity
    !         -----------------------------------------
    DO j=1,p_field
      rh(j,jlev) = rh(j,jlev)/smr(j)
      rh(j,jlev) = rh(j,jlev)*100.0
    END DO

  ELSE IF (krmode == 2) THEN

    !         Convert Relative Humidity to Mixing Ratio
    !         -----------------------------------------
    DO j=1,p_field
      rh(j,jlev) = rh(j,jlev)*smr(j)
      rh(j,jlev) = rh(j,jlev)*0.01
    END DO


  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE hmrtorh
END MODULE hmrtorh_mod
