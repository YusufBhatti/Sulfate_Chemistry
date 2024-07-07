! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for SCM namelist : INGEOFOR

MODULE s_ingeofor

USE scm_cntl_mod, ONLY: scm_nml

USE s_maxdim, ONLY:                                   &
  mx_rw_lng, mx_rw, mx_nobs, mx_mod_lv

USE scm_utils, ONLY:                                  &
  rmdi, rw_lng, rw, nobs

USE s_main_force, ONLY:                               &
  ug_opt, vg_opt                                      &
, scm_ug => ug                                        &
, scm_vg => vg

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE nlsizes_namelist_mod, ONLY: model_levels

USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE


!=============================================================================
!
! Description:
!   Allows SCM to read in INGEOFOR namelist from forcing file, scm_nml.
!   INGEOFOR contains information related to geostrophic wind forcing.
!
! Method:
!   Namelist INGEOFOR is defined in this module and read in by contained
!   subroutine read_ingeofor.  Scalar variables are read directly to
!   s_main_force.  Array variables are declared as private single rank arrays
!   of fixed length (arrays cannot be varying length in namelists) and
!   reshaped before being transferred to arrays of the correct size/shape in
!   s_main_force.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!
!=============================================================================
!
! Declarations:

  !--------------------------------------------------
  ! Fixed arrays to read in namelist before reshaping
  !--------------------------------------------------
REAL, PRIVATE ::                                &
  ug (mx_rw_lng*mx_rw*mx_nobs*mx_mod_lv) = rmdi &! Geostrophic u-wind (m/s)
, vg (mx_rw_lng*mx_rw*mx_nobs*mx_mod_lv) = rmdi  ! Geostrophic v-wind (m/s)


!---------------------------------------------------------------------------
! Define namelist
!---------------------------------------------------------------------------
NAMELIST/ingeofor/         &
  ug, vg, ug_opt, vg_opt

PRIVATE :: ingeofor

!=============================================================================
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='S_INGEOFOR'

CONTAINS

SUBROUTINE read_ingeofor

IMPLICIT NONE


! Dr Hook
!==============================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

INTEGER :: istatus
INTEGER :: icode

INTEGER :: i, j, j2, k

INTEGER, PARAMETER ::       &
  constant_value       = 0  &
, constant_profile     = 1  &
, time_varying_profile = 2  &
, time_varying_value   = 3
CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_INGEOFOR'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), ACTION='READ', IOSTAT=istatus,         &
     IOMSG=iomessage)

IF (istatus /= 0) THEN
  icode = 500
  WRITE(umMessage,*) ' Error opening file on unit 10 from '//routinename
  CALL umPrint(umMessage,src='s_ingeofor')
  WRITE(umMessage,*) ' Filename = '//TRIM(ADJUSTL(scm_nml))
  CALL umPrint(umMessage,src='s_ingeofor')
  WRITE(umMessage,*) ' IOstat =', istatus
  CALL umPrint(umMessage,src='s_ingeofor')
  WRITE(umMessage,*) ' Error: ', TRIM(iomessage)
  CALL umPrint(umMessage,src='s_ingeofor')
END IF

READ  (10, ingeofor)
CLOSE (10)


! Geostrophic - Zonal wind
!---------------------------------------
IF (ug(1) /= rmdi) THEN
  SELECT CASE (ug_opt)

  CASE (constant_value)
    ! Single value ug applied everywhere
    DO k=1, model_levels
      DO j2=1, nobs
        DO j=1, rw
          DO i=1, rw_lng
            scm_ug(i,j,j2,k) = ug(1)
          END DO
        END DO
      END DO
    END DO

  CASE (constant_profile)
    ! Single profile, non-time varying
    DO k=1, model_levels
      DO j2=1, nobs
        DO j=1, rw
          DO i=1, rw_lng
            scm_ug(i,j,j2,k) = ug(k)
          END DO
        END DO
      END DO
    END DO

  CASE (time_varying_profile)
    ! Time varying profile
    scm_ug = RESHAPE(ug ,(/rw_lng, rw, nobs, model_levels/))

  CASE (time_varying_value)
    ! Time varying value ug applied everywhere
    DO k=1, model_levels
      DO j2=1, nobs
        DO j=1, rw
          DO i=1, rw_lng
            scm_ug(i,j,j2,k) = ug(j2)
          END DO
        END DO
      END DO
    END DO

  END SELECT
END IF


! Geostrophic - Meridional wind
!---------------------------------------
IF (vg(1) /= rmdi) THEN
  SELECT CASE (vg_opt)

  CASE (constant_value)
    ! Single value vg applied everywhere
    DO k=1, model_levels
      DO j2=1, nobs
        DO j=1, rw
          DO i=1, rw_lng
            scm_vg(i,j,j2,k) = vg(1)
          END DO
        END DO
      END DO
    END DO

  CASE (constant_profile)
    ! Single profile, non-time varying
    DO k=1, model_levels
      DO j2=1, nobs
        DO j=1, rw
          DO i=1, rw_lng
            scm_vg(i,j,j2,k) = vg(k)
          END DO
        END DO
      END DO
    END DO

  CASE (time_varying_profile)
    ! Time varying profile
    scm_vg = RESHAPE(vg ,(/rw_lng, rw, nobs, model_levels/))

  CASE (time_varying_value)
    ! Time varying value vg applied everywhere
    DO k=1, model_levels
      DO j2=1, nobs
        DO j=1, rw
          DO i=1, rw_lng
            scm_vg(i,j,j2,k) = vg(j2)
          END DO
        END DO
      END DO
    END DO

  END SELECT
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE read_ingeofor

!=============================================================================
END MODULE s_ingeofor
