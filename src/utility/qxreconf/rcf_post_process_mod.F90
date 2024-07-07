! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Post `main loop' processing of data

MODULE Rcf_Post_Process_Mod

!  Subroutine Rcf_Post_Process - process fields after `main loop'
!
! Description:
! This subroutine performs post-processing of data using transforms
! etc that require a number of output-grid fields. This involves
! some transforms to original data types. Performed after the main
! field data creation loop.
!
! Method:
!  The relevant processing is just run through - this may need to be
!  be made more modular in future.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_POST_PROCESS_MOD'

CONTAINS

SUBROUTINE Rcf_Post_Process_Atmos( fields_in, field_count_in, orog_in,&
                                   fields_out, field_count_out, orog, &
                                   hdr_in, hdr_out, data_source)

USE Rcf_Adjust_Tsoil_Mod, ONLY: &
     Rcf_Adjust_Tsoil

USE Rcf_Field_Equals_Mod, ONLY: &
    Rcf_Field_Equals

USE Rcf_Exppx_Mod, ONLY: &
    Rcf_Exppx

USE Submodel_Mod, ONLY:  &
    atmos_im

USE Rcf_Calc_Rho_mod, ONLY: &
    Rcf_Calc_Rho

USE Rcf_Set_Interp_Flags_Mod, ONLY: &
    interp_done,                      &
    interp_no_op

USE Rcf_Freeze_Soil_Mod, ONLY: &
    Rcf_Freeze_Soil

USE Rcf_Data_Source_Mod, ONLY: &
    Data_Source_Type

USE items_nml_mod, ONLY:        &
    Ancillary_File,             &
    Input_Dump

USE Rcf_Conv_Cld_Chk_Mod, ONLY: &
    Rcf_Conv_Cld_Chk

USE Rcf_Cloud_Frac_Chk_Mod, ONLY: &
    Rcf_Cloud_Frac_Chk

USE Rcf_Snow_Amount_Chk_Mod, ONLY: &
    Rcf_Snow_Amount_Chk

USE Rcf_Sea_Ice_Frac_Chk_Mod, ONLY: &
    Rcf_Sea_Ice_Frac_Chk

USE Rcf_Soil_Moist_Chk_Mod, ONLY: &
    Rcf_Soil_Moist_Chk

USE Rcf_WriteUMhdr_Mod, ONLY: &
    Rcf_WriteUMhdr

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_UMhead_Mod, ONLY:&
    um_header_type

USE Rcf_Grid_Type_Mod, ONLY: &
    input_grid,              &
    output_grid

USE decomp_params, ONLY: &
    decomp_rcf_output

USE Rcf_HeadAddress_Mod, ONLY: &
    RC_PressureTop,         &
    RC_AtmMoist,            &
    RC_AtmMass,             &
    RC_AtmEnergy,           &
    RC_EnergyCorr,          &
    LDC_RHCrit,             &
    FH_GridStagger,         &
    FH_GridStagger_Endgame

USE Rcf_V_Int_Ctl_Mod, ONLY: &
    v_int_active,             &
    v_int_active_soil

USE Rcf_Interp_Weights_Mod, ONLY: &
    h_int_active

USE Rcf_Locate_mod, ONLY: &
    Rcf_Locate

USE Rcf_Alloc_Field_mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

USE Rcf_generate_heights_mod, ONLY: &
    Rcf_generate_heights

USE rcf_theta_t_convs_mod, ONLY: &
    Rcf_conv_theta_t,             &
    Rcf_conv_t_theta

USE Rcf_read_field_mod, ONLY: &
    Rcf_read_field

USE Rcf_Calc_Output_Exner_Mod, ONLY: &
    Rcf_Calc_Output_Exner

USE Rcf_Write_field_mod, ONLY: &
    Rcf_write_field

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_Exner_P_Convs_Mod, ONLY: &
    Rcf_Conv_Exner_P,             &
    Rcf_Conv_P_Exner

USE nlstcall_mod, ONLY:   &
    LTimer

USE um_stashcode_mod, ONLY:&
    stashcode_theta,          stashcode_q,       &
    stashcode_rho,            stashcode_exner,        &
    stashcode_p,              stashcode_soil_moist,   &
    stashcode_soil_temp,      stashcode_orog,         &
    stashcode_prog_sec

USE Rcf_Soilstress_To_Soilmoist_Mod,ONLY:&
    Rcf_Soilstress_to_soilmoist

USE rcf_soilconc_to_soilmoist_mod,ONLY:&
    rcf_soilconc_to_soilmoist

USE rcf_nlist_recon_science_mod, ONLY: &
    l_adj_t_soil,                      &
    polar_check,                       &
    use_smc_stress

USE Rcf_Polar_Rows_Chk_Mod, ONLY:  &
    Rcf_Polar_Rows_Chk

USE Rcf_Lsm_Mod, ONLY: &
    local_lsm_out

USE cppxref_mod, ONLY: &
    ppx_atm_tall

USE stparam_mod, ONLY:     &
    st_levels_model_theta, &
    st_levels_model_rho

USE rcf_lsh_land_ice_chk_mod, ONLY: &
    rcf_lsh_land_ice_chk

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER           :: fields_in(:)
TYPE( field_type ), POINTER           :: fields_out(:)
TYPE( field_type ), INTENT(INOUT)     :: orog_in
TYPE( field_type ), INTENT(IN)        :: orog
TYPE( data_source_type ), POINTER     :: data_source(:)
TYPE( um_header_type ), INTENT(IN)    :: hdr_in
TYPE( um_header_type ), INTENT(INOUT) :: hdr_out
INTEGER, INTENT(IN)                   :: field_count_in
INTEGER, INTENT(IN)                   :: field_count_out

! Local variables
INTEGER                       :: pos
INTEGER                       :: pos_st   !} positions in array of
INTEGER                       :: pos_sm   !} soil temp and moisture
INTEGER                       :: pos_orog !} and of orography
INTEGER                       :: ErrorStatus
CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'RCF_POST_PROCESS_ATMOS'
CHARACTER (LEN=errormessagelength)            :: Cmessage

REAL                          :: theta_heights(                     &
                                    output_grid % loc_p_field,      &
                                    0 : output_grid % model_levels + 1)
                                    ! Height of theta levels
REAL                          :: rho_heights(                       &
                                    output_grid % loc_p_field,      &
                                    0 : output_grid % model_levels + 1)
                                    ! Height of rho levels
LOGICAL                       :: exner_out ! is exner in the output?
LOGICAL                       :: l_soil_change   ! is soil changed

TYPE( field_type ), POINTER   :: exner
TYPE( field_type ), POINTER   :: theta
TYPE( field_type ), POINTER   :: q
TYPE( field_type ), POINTER   :: rho
TYPE( field_type ), POINTER   :: p
TYPE( field_type ), TARGET    :: pressure  ! exner OR p depending on
                                           ! circumstance
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (LTimer) CALL Timer( RoutineName, 3)

l_soil_change = .FALSE.     ! initialise soil moisture correction flag

CALL Rcf_Locate( stashcode_prog_sec, stashcode_soil_moist,           &
                 fields_out, field_count_out, pos_sm, zero_ok_arg = .TRUE. )

IF (pos_sm /= 0) THEN
  IF (data_source( pos_sm ) % source == input_dump) THEN
    IF (use_smc_stress) THEN
      CALL rcf_soilstress_to_soilmoist( fields_out, field_count_out,           &
                                        output_grid, decomp_rcf_output, hdr_out)
    ELSE
      IF (v_int_active_soil) THEN
        CALL rcf_soilconc_to_soilmoist( fields_out, field_count_out,           &
                                        output_grid, decomp_rcf_output, hdr_out)
      END IF
    END IF
  END IF
END IF

!-------------------------------------------------------------------
! Find and setup Theta (will hold T if interpolated)
!-------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_theta,                &
                 fields_out, field_count_out, pos)
theta => fields_out( pos )
CALL Rcf_Alloc_Field( theta )
CALL Rcf_Read_Field( theta, hdr_out, decomp_rcf_output )

IF ( h_int_active .OR. v_int_active ) THEN
  !-------------------------------------------------------------------
  ! Calculate required heights
  !-------------------------------------------------------------------
  CALL rcf_generate_heights( output_grid, orog,                   &
                             ppx_atm_tall, st_levels_model_theta, &
                             theta_heights, theta % level_size )

  CALL rcf_generate_heights( output_grid, orog,                   &
                             ppx_atm_tall, st_levels_model_rho,   &
                             rho_heights,  theta % level_size )
END IF

!--------------------------------------------------------------------
! Find and read Q
!--------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_q,                    &
                 fields_out, field_count_out, pos )
q => fields_out( pos )
CALL Rcf_Alloc_Field( q )
CALL Rcf_Read_Field( q, hdr_out, decomp_rcf_output )

!-------------------------------------------------------------------
! Find exner and P (as required)
! Note assumption that will have exner *OR* P available
!-------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_exner,                &
                 fields_out, field_count_out, pos, zero_ok_arg = .TRUE. )

IF (pos /= 0) THEN
  exner => fields_out( pos )
  CALL Rcf_Alloc_Field( exner )
  exner_out = .TRUE.

  ! Also need space for P
  p => pressure
  CALL Rcf_Field_Equals( p, exner )
  CALL Rcf_Alloc_Field( p )
  p % stashmaster => Rcf_Exppx( atmos_im, 0, stashcode_p )

ELSE    ! exner not in output dump

  CALL Rcf_Locate( stashcode_prog_sec, stashcode_p,                  &
                   fields_out, field_count_out, pos )
  p => fields_out( pos )
  CALL Rcf_Alloc_Field( p )
  exner_out = .FALSE.

  ! Also need space for exner
  exner => pressure
  CALL Rcf_Field_Equals( exner, p )
  CALL Rcf_Alloc_Field( exner )
  exner % stashmaster => Rcf_Exppx( atmos_im, 0, stashcode_exner )

END IF

!--------------------------------------------------------------------
! Either read or calculate exner (as appropriate)
!--------------------------------------------------------------------
IF ( exner % interp /= interp_no_op ) THEN
  IF ( exner_out ) THEN
    CALL Rcf_Read_Field( exner, hdr_out, decomp_rcf_output )
  ELSE
    CALL Rcf_Read_Field( p, hdr_out, decomp_rcf_output )
    exner % DATA(:,:) = p % DATA(:,:)
    exner % stashmaster => Rcf_Exppx( atmos_im, 0, stashcode_p )
    CALL Rcf_Conv_P_Exner( exner )
  END IF

ELSE
  IF ( theta % interp /= interp_done ) THEN
    ErrorStatus = 10
    Cmessage = 'Only have Theta available - need T to calculate exner'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  CALL Rcf_Locate( stashcode_prog_sec, stashcode_orog,               &
                   fields_out, field_count_out, pos_orog )
  CALL Rcf_Calc_Output_Exner( fields_in, field_count_in, orog_in,     &
                              hdr_in, orog, theta, q, exner,          &
                              Data_Source( pos_orog ) % source,       &
                              rho_heights, theta_heights )

  IF (exner_out) THEN
    CALL Rcf_Write_Field( exner, hdr_out, decomp_rcf_output )
  END IF

  ! Also calculate P - needed internally. Note need to temporarily
  ! fool with stashcode
  p % DATA(:,:) = exner % DATA(:,:)
  p % stashmaster => Rcf_Exppx( atmos_im, 0, stashcode_exner )
  CALL Rcf_Conv_Exner_P( p )

  ! Write this newly calcuated field to dump
  IF (.NOT. exner_out) THEN
    CALL Rcf_Write_Field( p, hdr_out, decomp_rcf_output )
  END IF

END IF

! *******************************************************************
! Convert T (in dump) back to theta
! *******************************************************************
IF ( theta % interp == interp_done ) THEN
  CALL Rcf_Conv_T_Theta( theta, fields_out, field_count_out, hdr_out, &
                         decomp_rcf_output, rho_heights, theta_heights)

  CALL Rcf_Write_Field( theta, hdr_out, decomp_rcf_output )
END IF


! ******************************************************************
! Calculate rho
! ******************************************************************
CALL Rcf_Locate( stashcode_prog_sec, stashcode_rho,                  &
                 fields_out, field_count_out, pos)
rho => fields_out(pos)
IF ( rho % interp == interp_no_op ) THEN

  ! No interpolation for rho need to calculate
  CALL Rcf_Alloc_Field( rho )

  CALL Rcf_Calc_Rho( theta, q, exner, p, theta_heights, rho_heights, &
                     rho)

  CALL Rcf_Write_Field( rho, hdr_out, decomp_rcf_output )
  CALL Rcf_DeAlloc_Field( rho )

END IF

!--------------------------------------------------------------------
! Tidy up
!--------------------------------------------------------------------
CALL Rcf_DeAlloc_Field( p )
CALL Rcf_DeAlloc_Field( q )
CALL Rcf_DeAlloc_Field( exner )
CALL Rcf_DeAlloc_Field( theta )

! *******************************************************************
! Check consistency between soil moisture and land use/soil type
! *******************************************************************
IF (h_int_active .OR. use_smc_stress .OR. v_int_active_soil) THEN
  CALL Rcf_Soil_Moist_Chk( fields_out, field_count_out, output_grid, &
                           decomp_rcf_output, hdr_out, l_soil_change)
END IF
! *******************************************************************
! Check some LSH fields for consistency with land-ice.
! *******************************************************************
! We should always call this to make sure any pre-existing dumps are consistent
! even when not interpolating.
CALL rcf_lsh_land_ice_chk( fields_out, field_count_out,              &
                           decomp_rcf_output, hdr_out)

! *******************************************************************
! Adjust soil temperatures to orography lapse rate if interpolated
! *******************************************************************
IF ( l_adj_t_soil .AND. h_int_active) THEN
  CALL Rcf_Adjust_Tsoil( fields_out, field_count_out,                   &
                         decomp_rcf_output, hdr_out,                    &
                         orog_in, orog,                                 &
                         input_grid, output_grid )
END IF

! *******************************************************************
! Adjust frozen/unfrozen soil if updated from ancillary
! *******************************************************************
CALL Rcf_Locate( stashcode_prog_sec, stashcode_soil_temp,            &
                 fields_out, field_count_out, pos_st, zero_ok_arg = .TRUE. )


! Only need to correct if soil moisture/temperature are both in dump
IF (pos_st /= 0 .AND. pos_sm /= 0) THEN
  IF ( Data_Source( pos_st ) % Source == Ancillary_File  .OR. &
       Data_Source( pos_sm ) % Source == Ancillary_File  .OR. &
       l_soil_change ) THEN

    CALL Rcf_Freeze_Soil( fields_out, field_count_out, hdr_out, &
                          output_grid )
  END IF
END IF

! *******************************************************************
! Check convective cloud base and top are sensible
! *******************************************************************
CALL Rcf_Conv_Cld_Chk( fields_out, field_count_out, output_grid,  &
                       decomp_rcf_output, hdr_out )

! *******************************************************************
! Check that cloud fraction fields are consistent with q fields.
! *******************************************************************
CALL Rcf_Cloud_Frac_Chk( fields_out, field_count_out, output_grid,  &
                         decomp_rcf_output, hdr_out )

! *******************************************************************
! Check for negative snow amounts
! *******************************************************************
IF (h_int_active) THEN
  CALL Rcf_Snow_Amount_Chk( fields_in, field_count_in, hdr_in,       &
                            fields_out, field_count_out, hdr_out )
END IF

! *******************************************************************
! Check for small sea-ice fractions if interpolated
! *******************************************************************
IF (h_int_active) THEN
  CALL Rcf_Sea_Ice_Frac_Chk( fields_out, field_count_out,            &
                             decomp_rcf_output, hdr_out, data_source)
END IF

! *******************************************************************
! If appropriate, ensure scalar fields in the polar rows are uniform
! *******************************************************************
IF ( output_grid % global .AND. polar_check ) THEN
  IF ( Hdr_Out % Fixhd(FH_GridStagger) == FH_GridStagger_Endgame ) THEN
    ErrorStatus = -20
    WRITE(cmessage,*) 'Polar row checking activated for an invalid grid ',  &
                      '- disabling test'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  ELSE
    CALL Rcf_Polar_Rows_Chk( fields_out, field_count_out, output_grid,      &
                             decomp_rcf_output, hdr_out, local_lsm_out,     &
                             data_source )
  END IF
END IF

!----------------------------------------
! Perform some modifications to headers.
!----------------------------------------
! This would be better earlier in the code when the headers are setup
! but there is logic after the headers are setup to turn on interpolation
! due to change in orography ancillary.  This is the next best place to
! perform this task.
IF ( .NOT. v_int_active ) THEN
  IF ( .NOT. h_int_active ) THEN
    ! Want to keep energy correction information if no interpolation
    hdr_out % realc( rc_atmmoist   ) = hdr_in % realc( rc_atmmoist   )
    hdr_out % realc( rc_atmmass    ) = hdr_in % realc( rc_atmmass    )
    hdr_out % realc( rc_atmenergy  ) = hdr_in % realc( rc_atmenergy  )
    hdr_out % realc( rc_energycorr ) = hdr_in % realc( rc_energycorr )
  END IF
  ! We want to copy the rhcrit values across if not vertically interpolating.
  hdr_out % levdepc(:,ldc_rhcrit) = hdr_in % levdepc(:,ldc_rhcrit)
  !---------------------------------------------------------------
  ! Write out the header here again due to changes above.
  !---------------------------------------------------------------
  CALL rcf_writeumhdr( hdr_out )
END IF

IF (LTimer) CALL Timer( RoutineName, 4)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Post_Process_Atmos
END MODULE Rcf_Post_Process_Mod
