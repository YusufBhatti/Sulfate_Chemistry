! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   magic numbers and code for setting up field interpolation flags

MODULE Rcf_Set_Interp_Flags_Mod

IMPLICIT NONE

!  Subroutine Rcf_Set_Interp_Flags - set field interpolation flags
!
! Description:
! This Module contains magic numbers and code for setting up
! field interpolation flags
!
! Method:
!  Based on the v_int_active and h_int_active logical switches.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


!------------------------------------------------------------------
! Magic Numbers
!------------------------------------------------------------------
INTEGER, PARAMETER     :: interp_all    = 1    ! h and v to be done
INTEGER, PARAMETER     :: interp_h_only = 2    ! Do h - copy for v
INTEGER, PARAMETER     :: interp_v_only = 3    ! Do v - copy for h
INTEGER, PARAMETER     :: interp_copy   = 4    ! copy for h and v
INTEGER, PARAMETER     :: interp_no_op  = 5    ! No interp. operations
INTEGER, PARAMETER     :: interp_done   = 6    ! interp. completed
INTEGER, PARAMETER     :: interp_copied = 7    ! copy completed

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SET_INTERP_FLAGS_MOD'

CONTAINS

!--------------------------------------------------------------------
! Subroutine to set interp flags for intput fields based on
! h_int_active, v_int_active and internal rules
!--------------------------------------------------------------------
SUBROUTINE Rcf_Set_Interp_Flags( fields_in, fields_out, field_count_in,&
                                 field_count_out, data_source )

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_Data_Source_Mod, ONLY: &
    Data_Source_Type,           &
    Other_Field

USE items_nml_mod, ONLY: &
    Input_Dump

USE Rcf_V_Int_Ctl_Mod, ONLY: &
    v_int_active,            &
    v_int_active_soil

USE Rcf_Interp_Weights_Mod, ONLY: &
    h_int_active,                  &
    h_int_active_u,                &
    h_int_active_v,                &
    h_int_method,                  &
    nearest_neighbour

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE um_stashcode_mod, ONLY: &
    stashcode_ccb,             &
    stashcode_cct,             &
    stashcode_rho,             &
    stashcode_exner,           &
    stashcode_land_frac,       &
    stashcode_p,               &
    stashcode_prog_sec

USE Rcf_locate_alt_field_mod, ONLY: &
    rcf_locate_alt_field

USE stparam_mod, ONLY: &
    st_levels_single,  &
    st_levels_deep_soil

USE cppxref_mod, ONLY: &
    ppx_atm_ozone,     &
    ppx_atm_cuall,     &
    ppx_atm_cvall

USE rcf_nlist_recon_technical_mod, ONLY: &
    l_basic_interp

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER       :: fields_in(:)
TYPE( field_type ), POINTER       :: fields_out(:)
TYPE( data_source_type ), POINTER :: data_source(:)
INTEGER, INTENT(IN)               :: field_count_in
INTEGER, INTENT(IN)               :: field_count_out

! Local variables
INTEGER                           :: pos
INTEGER                           :: i
LOGICAL                           :: vertical
LOGICAL                           :: horizontal

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_SET_INTERP_FLAGS'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------
! Find the output fields that should be sourced from the input dump
!------------------------------------------------------------------
DO i = 1, field_count_out
  IF ( data_source(i) % Source == Input_Dump ) THEN
    CALL Rcf_Locate( fields_out(i) % stashmaster % section,          &
                     fields_out(i) % stashmaster % item,             &
                     fields_in, field_count_in, pos )
  ELSE IF ( data_source(i) % Source == Other_Field ) THEN
    CALL Rcf_Locate_Alt_Field( fields_out(i),                        &
                               fields_in, field_count_in, pos)
  ELSE
    CYCLE
  END IF



  ! Interpolation decisions based on standard interpolation flags
  vertical   = .FALSE.
  horizontal = .FALSE.
  IF (v_int_active) THEN
    vertical = .TRUE.
  END IF

  IF (h_int_active) THEN
    horizontal = .TRUE.
  END IF

  !------------------------------------------------------------------
  ! Override standard decisions
  !------------------------------------------------------------------
    ! Only do h interpolation if single level (not ccb or cct)
  IF ( fields_in(pos) % stashmaster % lv_code ==              &
                                      st_levels_single ) THEN
    IF ( .NOT.                                                        &
        ( fields_in(pos) % stashmaster % section ==                   &
          stashcode_prog_sec .AND.                                    &
          ( fields_in(pos) % stashmaster % item == stashcode_ccb .OR. &
            fields_in(pos) % stashmaster % item == stashcode_cct ) ) ) THEN

      vertical = .FALSE.

    END IF
  END IF

  ! Soil levels need special treatment
  IF ( fields_in(pos) % stashmaster % lv_code == st_levels_deep_soil) THEN
    IF ( v_int_active_soil ) THEN

      vertical = .TRUE.

    ELSE

      vertical = .FALSE.

    END IF
  END IF

  ! We could have P grid the same but the U or V grid are different.  So lets
  ! allow this to occur to minimise interpolation (e.g. ND LAM -> EG LAM)
  IF ( fields_in(pos) % stashmaster % grid_type == ppx_atm_cuall) THEN
    horizontal = h_int_active_u
  END IF
  IF ( fields_in(pos) % stashmaster % grid_type == ppx_atm_cvall) THEN
    horizontal = h_int_active_v
  END IF

  ! Treat ozone due to it using its own grid.
  IF ( fields_in(pos) % stashmaster % grid_type == ppx_atm_ozone) THEN

    IF ( fields_in(pos) % levels /= fields_out(i) % levels ) THEN     

      ! If Ozone grid and levels differ then turn on vertical
      ! interpolation whatever the main interpolation flag
      vertical = .TRUE.
    END IF

    IF ( fields_in(pos) % glob_row_len /= fields_out(i) % glob_row_len ) THEN
      ! We might be moving from full field to zonal (or vice versa).
      horizontal = .TRUE.
    END IF
  END IF


  !------------------------------------------------------------------
  ! Set the field interp flag
  !------------------------------------------------------------------
  IF (vertical .AND. horizontal) THEN
    fields_in(pos) % interp = interp_all

  ELSE IF (vertical) THEN
    fields_in(pos) % interp = interp_v_only

  ELSE IF (horizontal) THEN
    fields_in(pos) % interp = interp_h_only

  ELSE
    fields_in(pos) % interp = interp_copy

  END IF

  ! For nearest neighbour we dont care about rebalancing atmosphere.
  ! Also do not rebalance when performing basic interpolation
  IF (h_int_method /= nearest_neighbour .AND. &
      .NOT. l_basic_interp ) THEN

    ! Special treatment for rho, exner and p as need calculation rather
    ! than interpolation if would have been interpolated
    IF ( fields_in(pos) % stashmaster % section == stashcode_prog_sec .AND. &
         (fields_in(pos) % stashmaster % item == stashcode_rho        .OR.  &
         fields_in(pos) % stashmaster % item == stashcode_exner      .OR.   &
         fields_in(pos) % stashmaster % item == stashcode_p    )            &
         .AND. ( horizontal .OR. vertical )             ) THEN

      fields_in(pos) % interp = interp_no_op

    END IF
  END IF

  IF ( fields_in(pos) % stashmaster % section == stashcode_prog_sec .AND. &
       fields_in(pos) % stashmaster % item == stashcode_land_frac   .AND. &
       ( horizontal .OR. vertical )             ) THEN
    IF (vertical .AND. horizontal) THEN
      fields_in(pos) % interp = interp_v_only
    ELSE IF (vertical) THEN
      fields_in(pos) % interp = interp_v_only
    ELSE
      fields_in(pos) % interp = interp_no_op
    END IF
  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Set_Interp_Flags
END MODULE Rcf_Set_Interp_Flags_Mod
