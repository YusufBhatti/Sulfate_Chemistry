! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE post_interp_transform_mod

! Description:
!   Perform checks and transformations on the interpolated
!   LBC fields
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.
USE field_mod, ONLY: field_type, ASSIGNMENT(=)

USE umPrintMgr, ONLY:  umPrint, umMessage, PrintStatus, PrStatus_Oper
! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'POST_INTERP_TRANSFORM_MOD'

CONTAINS

SUBROUTINE post_interp_transform(field, q_min)

USE um_stashcode_mod, ONLY:                                       &
     stashcode_lbc_w,                stashcode_lbc_w_adv,         &
     stashcode_lbc_qcf,              stashcode_lbc_qcl,           &
     stashcode_lbc_qcf2,             stashcode_lbc_qrain,         &
     stashcode_lbc_qgraup,           stashcode_lbc_cf_bulk,       &
     stashcode_lbc_cf_liquid,        stashcode_lbc_cf_frozen,     &
     stashcode_lbc_dust1_mmr,        stashcode_lbc_dust2_mmr,     &
     stashcode_lbc_dust3_mmr,        stashcode_lbc_dust4_mmr,     &
     stashcode_lbc_dust5_mmr,        stashcode_lbc_dust6_mmr,     &
     stashcode_lbc_so2,              stashcode_lbc_dms,           &
     stashcode_lbc_so4_aitken,       stashcode_lbc_so4_accu,      &
     stashcode_lbc_nh3,              stashcode_lbc_soot_new,      &
     stashcode_lbc_soot_agd,         stashcode_lbc_soot_cld,      &
     stashcode_lbc_bmass_new,        stashcode_lbc_bmass_agd,     &
     stashcode_lbc_bmass_cld,        stashcode_lbc_ocff_new,      &
     stashcode_lbc_ocff_agd,         stashcode_lbc_ocff_cld,      &
     stashcode_lbc_nitr_acc,         stashcode_lbc_nitr_diss,     &
     stashcode_lbc_so4_diss,         stashcode_lbc_ukca_150,      &
     stashcode_lbc_murk,             stashcode_lbc_free_tracer_1, &
     stashcode_lbc_free_tracer_150,  stashcode_lbc_ukca_1,        &
     stashcode_lbc_q
     

IMPLICIT NONE

TYPE(field_type), INTENT(INOUT) :: field
REAL, INTENT(IN) :: q_min
CHARACTER(LEN=*), PARAMETER :: routinename = 'POST_INTERP_TRANSFORM'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE(field%quantity_ident)

! Set bottom and top level vertical wind to zero
CASE (stashcode_lbc_w, stashcode_lbc_w_adv)
  field%lbc_rdata(:,1) = 0.0
  field%lbc_rdata(:,field%get_num_levels()) = 0.0

! Ensure no negative values
CASE (stashcode_lbc_qcf,              stashcode_lbc_qcl,         &
      stashcode_lbc_qcf2,             stashcode_lbc_qrain,       &
      stashcode_lbc_qgraup,           stashcode_lbc_so4_diss,    &
      stashcode_lbc_dust1_mmr,        stashcode_lbc_dust2_mmr,   &
      stashcode_lbc_dust3_mmr,        stashcode_lbc_dust4_mmr,   &
      stashcode_lbc_dust5_mmr,        stashcode_lbc_dust6_mmr,   &
      stashcode_lbc_so2,              stashcode_lbc_dms,         &
      stashcode_lbc_so4_aitken,       stashcode_lbc_so4_accu,    &
      stashcode_lbc_nh3,              stashcode_lbc_soot_new,    &
      stashcode_lbc_soot_agd,         stashcode_lbc_soot_cld,    &
      stashcode_lbc_bmass_new,        stashcode_lbc_bmass_agd,   &
      stashcode_lbc_bmass_cld,        stashcode_lbc_ocff_new,    &
      stashcode_lbc_ocff_agd,         stashcode_lbc_ocff_cld,    &
      stashcode_lbc_nitr_acc,         stashcode_lbc_nitr_diss,   &
      stashcode_lbc_free_tracer_1:stashcode_lbc_free_tracer_150, &
      stashcode_lbc_ukca_1:stashcode_lbc_ukca_150 )
  CALL field_set_to_min(field, 0.0)

! Ratios should be between 0 and 1
CASE (stashcode_lbc_cf_bulk,          stashcode_lbc_cf_liquid, &
      stashcode_lbc_cf_frozen)
  CALL field_set_to_max(field, 1.0)
  CALL field_set_to_min(field, 0.0)

! Reset q to minimum value specifed from namelist
CASE (stashcode_lbc_q)
  CALL field_set_to_min(field, q_min)
    
CASE (stashcode_lbc_murk)
  CALL field_set_to_min(field, 0.1)
       
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE post_interp_transform

SUBROUTINE field_set_to_max(field, max_value)
IMPLICIT NONE
TYPE(field_type), INTENT(INOUT) :: field
REAL, INTENT(IN) :: max_value
INTEGER :: i, k, counter, lbc_level_size
CHARACTER(LEN=*), PARAMETER :: routinename = 'FIELD_SET_TO_MAX'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

lbc_level_size = SIZE(field%lbc_rdata, 1)
DO k = 1, field%grid%vert_grid%num_levels
  counter = 0
  DO i = 1, lbc_level_size
    IF (field%lbc_rdata(i, k) > max_value) THEN
      field%lbc_rdata(i, k) = max_value
      counter = counter + 1
    END IF
  END DO
  IF (counter > 0 .AND. PrintStatus >= PrStatus_Oper) THEN
    WRITE(umMessage,'(A,I0,A,F5.2,A,I0,A,I0)') "Field STASHcode: ",      &
         field%quantity_ident, " reset to maximum value of ", max_value, &
         " for ", counter, " grid points on level ", k
    CALL umPrint(umMessage,src='rcf_post_interp_transform_mod')
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE field_set_to_max

SUBROUTINE field_set_to_min(field, min_value)
IMPLICIT NONE
TYPE(field_type), INTENT(INOUT) :: field
REAL, INTENT(IN) :: min_value
INTEGER :: i, k, counter, lbc_level_size
CHARACTER(LEN=*), PARAMETER :: routinename = 'FIELD_SET_TO_MIN'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

lbc_level_size = SIZE(field%lbc_rdata, 1)
DO k = 1, field%grid%vert_grid%num_levels
  counter = 0
  DO i = 1, lbc_level_size
    IF (field%lbc_rdata(i, k) < min_value) THEN
      field%lbc_rdata(i, k) = min_value
      counter = counter + 1
    END IF
  END DO
  IF (counter > 0 .AND. PrintStatus >= PrStatus_Oper) THEN
    WRITE(umMessage,'(A,I0,A,F5.2,A,I0,A,I0)') "Field STASHcode: ",       &
         field%quantity_ident,  " reset to minimum value of ", min_value, &
         " for ", counter, " grid points on level ", k
    CALL umPrint(umMessage,src='rcf_post_interp_transform_mod')
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE field_set_to_min

END MODULE post_interp_transform_mod
