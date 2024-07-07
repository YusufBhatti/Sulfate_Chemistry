! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Perform transformations etc before a field is interpolated

MODULE Rcf_Pre_Interp_Transform_Mod

!  Subroutine Rcf_Pre_Interp_Transform
!
! Description:
!   Wrapper to perform pre-interpolation transforms to those fields
!   that require it.
!
! Method:
!   Choice of transform is based on stashcode.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_PRE_INTERP_TRANSFORM_MOD'

CONTAINS

SUBROUTINE Rcf_Pre_Interp_Transform( input_field, fields_in,     &
                                     field_count_in, hdr_in,     &
                                     fields_out, field_count_out,&
                                     hdr_out, data_source,       &
                                     orog_in )

USE UM_ParCore, ONLY: &
    mype

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_V_Int_Ctl_Mod, ONLY: &
    v_int_active_soil

USE Rcf_Init_Tile_T_Mod, ONLY: &
    Rcf_Init_Tile_T

USE Rcf_Exner_P_Convs_Mod, ONLY: &
    Rcf_Conv_Exner_P

USE umPrintMgr, ONLY:       &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE decomp_params, ONLY: &
    decomp_rcf_input

USE Rcf_Grid_Type_Mod, ONLY: &
    Input_Grid

USE Rcf_theta_t_convs_mod, ONLY: &
    Rcf_conv_theta_t

USE Rcf_Generate_Heights_Mod, ONLY: &
    Rcf_Generate_Heights

USE Rcf_Lsm_Mod, ONLY: &
    local_lsm_in

USE Rcf_Set_Interp_Flags_Mod, ONLY: &
    interp_h_only,                   &
    interp_v_only,                   &
    interp_all

USE um_stashcode_mod, ONLY: &
    stashcode_theta,           &
    stashcode_exner,           &
    stashcode_prog_sec,        &
    stashcode_soil_moist,      &
    stashcode_tstar_tile

USE rcf_nlist_recon_science_mod, ONLY: &
    l_use_zero_frac_tile_temp,         &
    use_smc_stress

USE Rcf_Data_Source_Mod, ONLY: &
    data_source_type

USE Rcf_smc_stress_Mod,ONLY: &
    Rcf_smc_stress

USE rcf_smc_conc_mod,ONLY: &
    rcf_smc_conc

USE Rcf_Exppx_Mod, ONLY: &
    Rcf_Exppx

USE Submodel_Mod, ONLY:  &
    atmos_im

USE cppxref_mod, ONLY: &
    ppx_atm_tall

USE stparam_mod, ONLY:     &
    st_levels_model_theta, &
    st_levels_model_rho

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), INTENT(INOUT)     :: input_field
TYPE( field_type ), POINTER           :: fields_in(:)
TYPE( field_type ), POINTER           :: fields_out(:)
TYPE( field_type ), INTENT(IN)        :: orog_in
TYPE( um_header_type ), INTENT(IN)    :: hdr_in
TYPE( um_header_type ), INTENT(IN)    :: hdr_out
TYPE (data_source_type), POINTER      :: data_source(:)
INTEGER,            INTENT(IN)        :: field_count_in
INTEGER,            INTENT(IN)        :: field_count_out

! Local variables
LOGICAL                       :: interpolate
LOGICAL                       :: l_init_tile_t_zerofrac
REAL                          :: theta_heights(                     &
                                    input_grid % loc_p_field,       &
                                    0 : input_grid % model_levels + 1)
                                    ! Height of theta levels
REAL                          :: rho_heights(                       &
                                    input_grid % loc_p_field,       &
                                    0 : input_grid % model_levels + 1)
                                    ! Height of rho levels
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_PRE_INTERP_TRANSFORM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

interpolate = ( input_field % interp == interp_h_only .OR. &
                input_field % interp == interp_v_only .OR. &
                input_field % interp == interp_all )
!---------------------------------------------------------------
! Only do transforms if interpolation is switched on
!---------------------------------------------------------------
IF ( interpolate ) THEN

  SELECT CASE( input_field % stashmaster % section )

  CASE (stashcode_prog_sec)
    ! Which fields do we wish to apply transforms to?
    SELECT CASE( input_field % stashmaster % item )

    CASE ( stashcode_exner )
      ! convert exner to P
      IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
        WRITE(umMessage,'(a25)') 'Converting exner to P'
        CALL umPrint(umMessage,src='rcf_pre_interp_transform_mod')
      END IF
      CALL Rcf_Conv_Exner_P( input_field )
      ! Reset stashmaster back to original.
      input_field % stashmaster => Rcf_Exppx( atmos_im, stashcode_prog_sec, &
                                              stashcode_exner )

    CASE ( stashcode_theta )
      ! convert theta to T
      IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
        WRITE(umMessage,'(a25)') 'Converting Theta to T'
        CALL umPrint(umMessage,src='rcf_pre_interp_transform_mod')
      END IF

      CALL rcf_generate_heights( input_grid, orog_in,                 &
                                 ppx_atm_tall, st_levels_model_theta, &
                                 theta_heights,                       &
                                 input_field % level_size)

      CALL rcf_generate_heights( input_grid, orog_in,                 &
                                 ppx_atm_tall, st_levels_model_rho,   &
                                 rho_heights,                         &
                                 input_field % level_size)

      CALL Rcf_Conv_Theta_T( input_field, fields_in,field_count_in, &
                             hdr_in, decomp_rcf_input, rho_heights, &
                             theta_heights )
      ! Reset STASHmaster back to correct one (only transforming data).
      input_field % stashmaster => Rcf_Exppx( atmos_im, stashcode_prog_sec, &
                                              stashcode_theta )

    END SELECT
  END SELECT

END IF

! Perform soil moisture processing either if:
! (i) there is a change in the horizontal grid: h_int_active
! (ii) the horizontal grid is unchanged but the soil properties have changed or
! (iii) there is a change in the vertical soil levels: v_int_active_soil
!
! (ii) or (iii) may be true even when interpolate=.false. so this is
! outside the if test on that logical
SELECT CASE( input_field % stashmaster % section )
CASE (stashcode_prog_sec)
  ! Which fields do we wish to apply transforms to?
  SELECT CASE( input_field % stashmaster % item )
  CASE ( stashcode_soil_moist )
    IF ( use_smc_stress ) THEN
      ! If user requests soil stress calculations, then apply 
      ! tests (i-iii) above in smc_stress routine
      CALL Rcf_smc_stress( input_field, fields_in,               &
                             field_count_in, hdr_in,             &
                             fields_out, field_count_out,        &
                             hdr_out, data_source)
    ELSE
      ! Test (iii) above
      IF (v_int_active_soil) THEN
        CALL rcf_smc_conc(   input_field )
      END IF
    END IF
  END SELECT
END SELECT


! For tstar_tile, overwrite zero fractions with sensible values before
! interpolation
IF ( .NOT. l_use_zero_frac_tile_temp .AND.                              &
     input_field % stashmaster % section == stashcode_prog_sec .AND.    &
     input_field % stashmaster % item == stashcode_tstar_tile ) THEN
  l_init_tile_t_zerofrac = .TRUE.
  CALL Rcf_Init_Tile_T(fields_in, field_count_in, hdr_in,               &
                       decomp_rcf_input, local_lsm_in,                  &
                       input_field, l_init_tile_t_zerofrac)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Pre_Interp_Transform
END MODULE Rcf_Pre_Interp_Transform_Mod
