! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Performs Source=9 field initialisation calculations for some
!   EG prognostics which require special handling

MODULE Rcf_Field_Dependent_Calcs_Mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE
!  Subroutine Rcf_Field_Dependent_Calcs - Field initialisation calculations
!  requiring other fields to have been evaluated first.
!
! Description:
!   Some fields may have hard-coded Source=9 initialisation. These fields
!   are initialised by the calculations in this routine.
!
! Method:
!   For these prognostics - the order the calculation is done is important,
!   which is why they can't be done in the same way as all other prognostics.
!   Choice of method applied determined by stashcode.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_FIELD_DEPENDENT_CALCS_MOD'

CONTAINS

SUBROUTINE Rcf_Field_Dependent_Calcs( fields_in, fields_out, field_count_in, &
                                      field_count_out, data_source, hdr_in,  &
                                      hdr_out )

USE Rcf_Field_Type_Mod, ONLY:                             &
    field_type

USE Rcf_Data_Source_Mod, ONLY:                            &
    data_source_type

USE items_nml_mod, ONLY: &
    Field_Dependent_Calcs

USE Rcf_Locate_mod, ONLY:                                 &
    Rcf_Locate

USE Rcf_UMhead_Mod, ONLY:                                 &
    um_header_type

USE um_stashcode_mod, ONLY:                               &
    stashcode_prog_sec,                                   &
    stashcode_thetavd,        stashcode_dry_rho,          &
    stashcode_exner_surf,     stashcode_etadot,           &
    stashcode_snowdep_grd_tile,                           &
    stashcode_snowpack_bk_dens,                           &
    stashcode_snow_laythk_tiles

USE Rcf_Grid_Type_Mod, ONLY:                              &
    Input_Grid, Output_Grid

USE Rcf_HeadAddress_Mod, ONLY:                            &
    Fh_GridStagger_Endgame

USE rcf_filter_exner_mod, ONLY:                           &
    rcf_filter_exner

USE Rcf_Recompute_Wet_Rho_Mod, ONLY:                      &
    Rcf_Recompute_Wet_Rho

USE Rcf_Derv_thetavd_mod, ONLY:                           &
    Rcf_Derv_thetavd

USE Rcf_Derv_Dry_rho_mod, ONLY:                           &
    Rcf_Derv_Dry_rho

USE Rcf_Derv_Etadot_mod, ONLY:                            &
    Rcf_Derv_Etadot

USE Rcf_Derv_Exner_Surf_mod, ONLY:                        &
    Rcf_Derv_Exner_Surf

USE Rcf_Alloc_Field_mod, ONLY:                            &
    Rcf_Alloc_Field,                                      &
    Rcf_DeAlloc_Field

USE Rcf_Write_Field_Mod, ONLY:                            &
    Rcf_Write_Field

USE Rcf_Eg_Poles_mod, ONLY:                               &
    Rcf_Eg_Poles

USE rcf_nlist_recon_technical_mod, ONLY:                  &
    tstmsk_to_decide,                                     &
    select_output_fields

USE decomp_params, ONLY:                                  &
    decomp_rcf_output

USE rcf_calc_tiles_mod, ONLY:                             &
    rcf_calc_tiles

USE rcf_snowstores_mod, ONLY:                             &
    rcf_snowstores

USE rcf_ml_snowpack_mod, ONLY:                            &
    rcf_ml_snowpack

USE nlstcall_mod, ONLY:                                   &
    LTimer

USE um_parvars, ONLY:                                     &
    bound

USE um_parparams, ONLY:                                   &
    bc_cyclic

USE errormessagelength_mod, ONLY:                         &
    errormessagelength

USE ereport_mod, ONLY :                                   &
    ereport

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER,            INTENT(IN)       :: field_count_in
INTEGER,            INTENT(IN)       :: field_count_out
TYPE( field_type ), POINTER          :: fields_in( : )
TYPE( field_type ), POINTER          :: fields_out( : )
TYPE( data_source_type ), POINTER    :: data_source( : )
TYPE( um_header_type ), INTENT(IN)   :: hdr_in
TYPE( um_header_type ), INTENT(IN)   :: hdr_out

! Local variables
INTEGER                              :: pos
INTEGER                              :: errorstatus
CHARACTER (LEN=errormessagelength)   :: Cmessage
CHARACTER (LEN=*), PARAMETER         :: RoutineName='RCF_FIELD_DEPENDENT_CALCS'

LOGICAL       :: recompute_wet_rho   = .FALSE.
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (LTimer) CALL Timer( 'Field_Depend_Calcs', 3)

! process any fields with dependencies

CALL Rcf_Locate (stashcode_prog_sec, stashcode_thetavd,       &
                 fields_out, field_count_out, pos ,zero_ok_arg = .TRUE.)

IF ( pos /= 0 ) THEN ! Thetavd is required in the output dump
  IF ( data_source( pos ) % Source == Field_Dependent_Calcs ) THEN
  
    CALL Rcf_Alloc_Field( fields_out( pos ) )
    ! Rcf_Derv_Thetavd depends on Theta and mv already being in the output
    ! dump and the Rcf_Locate calls within will abort if they're not found.
    CALL Rcf_Derv_Thetavd( fields_in, field_count_in,           &
                           fields_out, field_count_out,         &
                           hdr_in, hdr_out, fields_out( pos ) )
    CALL Rcf_Write_Field( fields_out( pos ), Hdr_out, decomp_rcf_output )
    CALL Rcf_DeAlloc_Field( fields_out( pos ) )
  
  END IF
END IF

!
! compute exner surf
!

CALL Rcf_Locate (stashcode_prog_sec, stashcode_exner_surf,    &
                 fields_out, field_count_out, pos ,zero_ok_arg = .TRUE.)

IF ( pos /= 0 ) THEN ! Exner_Surf is required in the output dump
  IF ( data_source( pos ) % Source == Field_Dependent_Calcs ) THEN
  
    CALL Rcf_Alloc_Field( fields_out( pos ) )
    ! Rcf_Derv_Exner_Surf depends on Orography, Exner and Rho already being in
    ! the output dump and the Rcf_Locate calls within will abort if they're not
    ! found.
    CALL Rcf_Derv_Exner_Surf( fields_in, field_count_in,        &
                           fields_out, field_count_out,         &
                           hdr_in, hdr_out, fields_out( pos ) )
    CALL Rcf_Write_Field( fields_out( pos ), Hdr_out, decomp_rcf_output )
    CALL Rcf_DeAlloc_Field( fields_out( pos ) )
  
  END IF
END IF


! when simply interpolating input data one cannot guarantee input contents
! thus ND2EG assumptions must be protected. Only perform if producing a UM
! dump and fields will be guarenteed by tstmsk routine.
IF ( select_output_fields == tstmsk_to_decide) THEN
  ! filtering of exner & exner surf and v_at_poles
  IF (input_grid % grid_stagger  /=  FH_Gridstagger_Endgame) THEN

    IF (output_grid  % grid_stagger  ==  FH_Gridstagger_Endgame) THEN

      IF (output_grid % global) THEN
        ! rcf_filter_exner depends on Exner_Surf and Exner already being in
        ! the output dump and the Rcf_Locate calls within will abort if they're
        ! not found. This call must follow the call to Rcf_Derv_Exner_Surf
        CALL rcf_filter_exner(fields_out, field_count_out, hdr_out)
        CALL rcf_EG_Poles    (fields_out, field_count_out, hdr_out)
      END IF

    END IF

  END IF

END IF

!
! compute dry Rho
!

CALL Rcf_Locate (stashcode_prog_sec, stashcode_dry_rho,       &
                 fields_out, field_count_out, pos ,zero_ok_arg = .TRUE.)

IF ( pos /= 0 ) THEN ! Dry_Rho is required in the output dump
  IF ( data_source( pos ) % Source == Field_Dependent_Calcs ) THEN
  
    CALL Rcf_Alloc_Field( fields_out( pos ) )
    ! Rcf_Derv_Dry_Rho depends on Thetavd and Exner already being in
    ! the output dump and the Rcf_Locate calls within will abort if they're not
    ! found. This call must follow the calls to rcf_filter_exner and
    ! Rcf_Derv_Thetavd
    CALL Rcf_Derv_Dry_Rho( fields_in, field_count_in,           &
                           fields_out, field_count_out,         &
                           hdr_in, hdr_out, fields_out( pos ) )
    CALL Rcf_Write_Field( fields_out( pos ), Hdr_out, decomp_rcf_output )
    CALL Rcf_DeAlloc_Field( fields_out( pos ) )
  
    recompute_wet_rho = .TRUE.
  
  END IF
END IF

!
! compute Etadot
!

CALL Rcf_Locate (stashcode_prog_sec, stashcode_etadot,        &
                 fields_out, field_count_out, pos ,zero_ok_arg = .TRUE.)

IF ( pos /= 0 ) THEN ! Eta_Dot is required in the output dump
  IF ( data_source( pos ) % Source == Field_Dependent_Calcs ) THEN
  
    CALL Rcf_Alloc_Field( fields_out( pos ) )
    ! Rcf_Derv_Etadot depends on u,v,w and Dry_Rho already being in
    ! the output dump and the Rcf_Locate calls within will abort if they're not
    ! found. This call must follow the call to Rcf_Derv_Dry_Rho
    CALL Rcf_Derv_Etadot( fields_in, field_count_in,            &
                           fields_out, field_count_out,         &
                           hdr_in, hdr_out, fields_out( pos ) )
    CALL Rcf_Write_Field( fields_out( pos ), Hdr_out, decomp_rcf_output )
    CALL Rcf_DeAlloc_Field( fields_out( pos ) )
  
  END IF
END IF

! now that we have a new dry rho, compute new wet rho that is
! consistent with this!

IF (recompute_wet_rho) THEN
  CALL Rcf_Recompute_Wet_Rho(fields_out, field_count_out, hdr_out)
END IF

!
! Snow Fields
!
! If the snow depth on tiles is required in the output dump, it may
! be appropriate to adjust the canopy and ground snow stores first. 
! The bulk density can be set only thereafter.
CALL Rcf_Locate (stashcode_prog_sec, stashcode_snowdep_grd_tile,        &
                 fields_out, field_count_out, pos ,zero_ok_arg = .TRUE.)

IF ( pos /= 0 ) THEN ! Snow depth is required in the output dump.
  IF ( data_source( pos ) % Source == Field_Dependent_Calcs ) THEN
    CALL rcf_snowstores( fields_in, field_count_in, hdr_in,    &
                         fields_out, field_count_out, hdr_out, &
                         data_source )
    CALL Rcf_Alloc_Field( fields_out( pos ) )
    CALL rcf_calc_tiles( fields_in, field_count_in, hdr_in,    &
                         fields_out, field_count_out, hdr_out, &
                         fields_out(pos),                      &
                         fields_out(pos) % stashmaster % item, &
                         data_source )
    CALL Rcf_Write_Field( fields_out( pos ), Hdr_out, decomp_rcf_output )
    CALL Rcf_DeAlloc_Field( fields_out( pos ) )
  END IF
END IF

! Repeat for the bulk density.
CALL Rcf_Locate (stashcode_prog_sec, stashcode_snowpack_bk_dens,        &
                 fields_out, field_count_out, pos ,zero_ok_arg = .TRUE.)

IF ( pos /= 0 ) THEN ! Bulk density is required in the output dump.
  IF ( data_source( pos ) % Source == Field_Dependent_Calcs ) THEN
    CALL Rcf_Alloc_Field( fields_out( pos ) )
    CALL rcf_calc_tiles( fields_in, field_count_in, hdr_in,    &
                         fields_out, field_count_out, hdr_out, &
                         fields_out(pos),                      &
                         fields_out(pos) % stashmaster % item, &
                         data_source )
    CALL Rcf_Write_Field( fields_out( pos ), Hdr_out, decomp_rcf_output )
    CALL Rcf_DeAlloc_Field( fields_out( pos ) )
  END IF
END IF

!
! Multilayer Snow Fields
!

! Use snow layer thickness to trigger setting of snow fields
CALL Rcf_Locate (stashcode_prog_sec, stashcode_snow_laythk_tiles,        &
                 fields_out, field_count_out, pos ,zero_ok_arg = .TRUE.)

IF ( pos /= 0 ) THEN ! Multilayer snow fields are required in the output dump
  IF ( data_source( pos ) % Source == Field_Dependent_Calcs ) THEN
    CALL rcf_ml_snowpack(  fields_in, field_count_in, hdr_in,   &
                           fields_out, field_count_out, hdr_out,&
                           data_source )
  END IF
END IF

IF (LTimer) CALL Timer( 'Field_Depend_Calcs', 4)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Field_Dependent_Calcs
END MODULE Rcf_Field_Dependent_Calcs_Mod
