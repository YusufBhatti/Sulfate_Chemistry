! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Performs Source=8 field initialisation calculations

MODULE Rcf_Field_Calcs_Mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE
!  Subroutine Rcf_Field_Calcs - field initialisation calculations.
!
! Description:
!   Some fields may have hard-coded Source=8 initialisation. These fields
!   are initialised by the calculations in this routine.
!
! Method:
!   Choice of method applied determined by stashcode.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_FIELD_CALCS_MOD'

CONTAINS

SUBROUTINE Rcf_Field_Calcs( fields_in, fields_out, field_count_in, &
                            field_count_out, data_source, hdr_in,  &
                            hdr_out )

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE rcf_data_source_Mod, ONLY:  &
    data_source_type

USE items_nml_mod, ONLY: &
    Field_Calcs

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Submodel_Mod, ONLY: atmos_im

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE decomp_params, ONLY: &
    decomp_rcf_output

USE Rcf_Alloc_Field_mod, ONLY:  &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

USE Rcf_Write_Field_Mod, ONLY: &
    Rcf_Write_Field

USE Rcf_Derv_Ice_Temp_Mod, ONLY: &
    Rcf_Derv_Ice_Temp

USE Rcf_Derv_2D_CCA_Mod, ONLY: &
    Rcf_Derv_2D_CCA

USE Rcf_Derv_Ice_Thick_Mod, ONLY: &
    Rcf_Derv_Ice_Thick

USE Rcf_Derv_Sea_Ice_Temp_Mod, ONLY:&
    Rcf_Derv_Sea_Ice_Temp

USE nlstcall_mod, ONLY:   &
    LTimer

USE rcf_init_soil_temp_mod, ONLY: rcf_init_soil_temp

USE rcf_init_tile_frac_mod, ONLY: rcf_init_tile_frac

USE Rcf_Init_Flake_Mod, ONLY: &
    Rcf_Init_Flake

USE Rcf_Calc_Gamtot_Mod, ONLY: &
    Rcf_Calc_Gamtot

USE Rcf_Est_Zw_Mod, ONLY:      &
    Rcf_Est_Zw

USE Rcf_Est_Sthzw_Mod, ONLY:      &
    Rcf_Est_Sthzw

USE Rcf_Fit_Fsat_Mod, ONLY: &
    Rcf_Fit_Fsat

USE Rcf_Calc_Fsat_Mod, ONLY:    &
    Rcf_Calc_Fsat

USE Rcf_Derv_Ice_Cat_Thick_Mod, ONLY:    &
    Rcf_Derv_Ice_Cat_Thick

USE Rcf_Derv_Adv_Winds_Mod, ONLY: &
    Rcf_Derv_Adv_Winds

USE Rcf_Change_Dust_Bins_Mod, ONLY: &
    Rcf_Change_Dust_Bins

USE Rcf_Mixing_Ratios_Mod, ONLY: &
    Rcf_Mixing_Ratios, moisture_field_sum


USE rcf_nlist_recon_science_mod, ONLY: &
    q_min

USE Rcf_Post_Interp_Transform_Mod, ONLY: &
    field_set_to_min

USE rcf_calc_tiles_mod, ONLY: &
    rcf_calc_tiles

USE rcf_set_dominant_tile_mod, ONLY:       &
    i_dominant

USE rcf_init_urban_macdonald_mod, ONLY:                   &
    rcf_init_urban_macdonald

USE jules_surface_mod, ONLY:                              &
    l_urban2t

USE switches_urban, ONLY :                                &
    l_moruses_macdonald

USE um_stashcode_mod, ONLY: &
    stashcode_prog_sec,                                    &
    stashcode_ice_conc_cat,    stashcode_ice_thick_cat,    &
    stashcode_ice_temp_cat,    stashcode_ice_snow_depth,   &
    stashcode_icethick,                                    &
    stashcode_cca,             stashcode_can_water_tile,   &
    stashcode_sea_ice_temp,                                &
    stashcode_ice_surf_temp_cat,                           &
    stashcode_tstar_ice_cat_cpl,                           &
    stashcode_catch_snow,      stashcode_snow_grnd,        &
    stashcode_rgrain,                                      &
    stashcode_snow_tile,       stashcode_tstar_tile,       &
    stashcode_tsurf_elev_surft,                            &
    stashcode_gamtot,          stashcode_zw,               &
    stashcode_sthzw,           stashcode_fsat,             &
    stashcode_fwetl,                                       &
    stashcode_a_fsat,          stashcode_c_fsat,           &
    stashcode_a_fwet,          stashcode_c_fwet,           &
    stashcode_u_adv,           stashcode_v_adv,            &
    stashcode_dust1_mmr,       stashcode_dust2_mmr,        &
    stashcode_dust3_mmr,       stashcode_dust4_mmr,        &
    stashcode_dust5_mmr,       stashcode_dust6_mmr,        &
    stashcode_flake_t_mean,    stashcode_flake_t_mxl,      &
    stashcode_flake_t_ice,     stashcode_flake_h_mxl,      &
    stashcode_flake_h_ice,     stashcode_mv,               &
    stashcode_mcl,             stashcode_mcf,              &
    stashcode_mr,              stashcode_mgr,              &
    stashcode_mcf2,                                        &
    stashcode_urbztm,          stashcode_urbdisp,          &
    stashcode_soil_temp,       stashcode_frac_surf_type

USE Rcf_Grid_Type_Mod, ONLY:  Input_Grid, Output_Grid

USE errormessagelength_mod, ONLY: errormessagelength

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
INTEGER                            :: i
INTEGER                            :: ErrorStatus
INTEGER                            :: pos
CHARACTER (LEN=*), PARAMETER       :: RoutineName='RCF_FIELD_CALCS'
CHARACTER (LEN=errormessagelength) :: Cmessage
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (LTimer) CALL Timer( RoutineName, 3)

!-----------------------------------------------------------------
! Loop around output fields/sources until a source=field_calcs
! is found
!-----------------------------------------------------------------
DO i = 1, field_count_out
  IF ( data_source( i ) % Source == Field_Calcs ) THEN

    CALL Rcf_Alloc_Field( fields_out( i ) )

    SELECT CASE( fields_out( i ) % stashmaster % model )

    CASE ( atmos_im )

      ! -----------
      ! Atmos items
      ! -----------

      SELECT CASE( fields_out( i ) % stashmaster % section )

        ! For the moment only section zero fields can be used here.
      CASE ( stashcode_prog_sec )

        SELECT CASE( fields_out(i) % stashmaster % item )

        CASE ( stashcode_cca )
          CALL Rcf_Derv_2D_CCA( fields_in, fields_out, field_count_in, &
                                field_count_out, hdr_in, hdr_out ,     &
                                fields_out( i ) )

        CASE ( stashcode_icethick )
          CALL Rcf_Derv_Ice_Thick( fields_out, field_count_out,        &
                                   hdr_out, fields_out(i) )

        CASE ( stashcode_sea_ice_temp )
          CALL Rcf_Derv_Sea_Ice_Temp( fields_out, field_count_out,  &
                                      fields_out( i ), hdr_out )

        CASE ( stashcode_soil_temp )
           CALL rcf_init_soil_temp( fields_out, field_count_out,        &
                                    decomp_rcf_output, hdr_out,         &
                                    fields_out(i) )

        CASE ( stashcode_frac_surf_type )
           CALL rcf_init_tile_frac( fields_out(i) )

        CASE ( stashcode_catch_snow         &
              ,stashcode_rgrain             &
              ,stashcode_can_water_tile     &
              ,stashcode_tstar_tile         &
              ,stashcode_tsurf_elev_surft   &
              ,stashcode_snow_tile          &
              ,stashcode_snow_grnd )
          CALL rcf_calc_tiles( fields_in, field_count_in, hdr_in,    &
                               fields_out, field_count_out, hdr_out, &
                               fields_out(i),                        &
                               fields_out(i) % stashmaster % item,   &
                               data_source )

        CASE ( stashcode_flake_t_mean     &
              ,stashcode_flake_t_mxl      &
              ,stashcode_flake_t_ice      &
              ,stashcode_flake_h_mxl      &
              ,stashcode_flake_h_ice        )
          CALL Rcf_Init_Flake( fields_out, field_count_out,       &
                               hdr_out, fields_out( i ),          &
                               fields_out(i) % stashmaster % item )

        CASE ( stashcode_ice_temp_cat, stashcode_ice_surf_temp_cat, &
                                 stashcode_tstar_ice_cat_cpl )
          CALL Rcf_Derv_Ice_Temp( fields_out,field_count_out, hdr_out, &
                                 fields_out( i) )


          ! NB The fields for LSH have quite a few inter-dependencies.  It seems
          ! that due to the dump being in STASHcode order the dependencies should be okay
          ! but please check code that if another field require calculations that it isnt
          ! used to calculate another field unless you are really sure.
        CASE ( stashcode_Gamtot )
          CALL Rcf_Calc_Gamtot( fields_out, field_count_out,  &
                                    i, hdr_out )
          ! If calculation required for the LSH fitting variables,
          ! must first allocate space for them.

        CASE ( stashcode_a_fsat, stashcode_c_fsat, &
               stashcode_a_fwet, stashcode_c_fwet)

          CALL Rcf_Fit_Fsat ( fields_out, field_count_out,      &
                                   i, hdr_out )

        CASE ( stashcode_Sthzw )
          CALL Rcf_Est_Sthzw( fields_out, field_count_out,  &
                                   fields_out( i ), hdr_out )

        CASE ( stashcode_Zw )
          CALL Rcf_Est_Zw( fields_out, field_count_out,  &
                                   fields_out( i ), hdr_out )

        CASE ( stashcode_Fsat, stashcode_fwetl )

          CALL Rcf_Calc_Fsat( fields_out, field_count_out,  &
                              i, hdr_out )

        CASE ( stashcode_ice_conc_cat, stashcode_ice_thick_cat )
          CALL Rcf_Derv_Ice_Cat_Thick (                               &
                                fields_out( i ) % stashmaster % item, &
                                fields_out, field_count_out,          &
                                hdr_out, fields_out( i ) )

        CASE ( stashcode_u_adv, stashcode_v_adv )
          CALL Rcf_Derv_Adv_Winds(                                   &
                               fields_out( i ) % stashmaster % item, &
                               fields_out, field_count_out,          &
                               hdr_out, fields_out( i ) )

        CASE ( stashcode_dust1_mmr, stashcode_dust2_mmr,                       &
               stashcode_dust3_mmr, stashcode_dust4_mmr,                       &
               stashcode_dust5_mmr, stashcode_dust6_mmr )
          CALL Rcf_Change_Dust_Bins( fields_in, field_count_in,                &
                                     fields_out, field_count_out,              &
                                     fields_out( i ) % stashmaster % item,     &
                                     hdr_in, hdr_out, fields_out( i ) )



        CASE ( stashcode_mv,   stashcode_mcl,                          &
               stashcode_mcf,  stashcode_mr,                           &
               stashcode_mgr,  stashcode_mcf2 )

          CALL Rcf_Mixing_Ratios( fields_out, field_count_out,         &
                                  fields_out( i ) % stashmaster % item,&
                                  hdr_out, fields_out( i ), data_source)
          ! apply qc checks consistent with those in applied to fields
          ! after interpolation  in rcf_post_interp_transform
          SELECT CASE( fields_out( i ) % stashmaster % item )
          CASE ( stashcode_mv )
            CALL field_set_to_min( "Vapour Mixing Ratio (mv)",         &
                                  fields_out( i ), field_min=q_min )

          CASE ( stashcode_mcl )
            CALL field_set_to_min( "Cld Liq Mixing Ratio (mcl)",       &
                                  fields_out( i ) )

          CASE ( stashcode_mcf )
            CALL field_set_to_min( "Cld Ice Mixing Ratio (mcf)",       &
                                  fields_out( i ) )

          CASE ( stashcode_mr )
            CALL field_set_to_min( "Rain Mixing Ratio (mr)",           &
                                  fields_out( i ) )

          CASE ( stashcode_mgr )
            CALL field_set_to_min( "Graupel Mixing Ratio (mg)",        &
                                  fields_out( i ) )

          CASE ( stashcode_mcf2 )
            CALL field_set_to_min( "Ice Cry Mixing Rat. (mcf2)",       &
                                  fields_out( i ) )

          END SELECT

        CASE ( stashcode_urbztm, stashcode_urbdisp )
          ! MORUSES MacDonald roughness and displacement height. Recalculate
          ! to ensure conformity with either the new ancillary or interpolated
          ! morphology.
          IF ( l_moruses_macdonald ) THEN
            IF ( l_urban2t ) THEN
              CALL rcf_init_urban_macdonald( fields_out, field_count_out, &
                                             hdr_out, data_source )
            ELSE
              errorstatus = 20
              WRITE(cmessage, '(A)') &
                 'l_moruses_macdonald = T, but l_urban2t = F. ' //            &
                 'Please set l_urban2t = T as the required stash will '//     &
                 'not be present.'
              CALL ereport( routinename, errorstatus, cmessage )
            END IF
          END IF

        CASE DEFAULT
          ErrorStatus = 30
          WRITE (Cmessage, '(A, I3, A, I5)')                         &
                'No Field Calculations specified for section ',      &
                fields_out( i ) % stashmaster % section, ' item ',   &
                fields_out( i ) % stashmaster % item
          CALL Ereport( RoutineName, ErrorStatus, Cmessage )

        END SELECT

      CASE DEFAULT
        ErrorStatus = 30
        WRITE (Cmessage, '(A, I3)')                                  &
              'No Field Calculations specified for section ',        &
              fields_out( i ) % stashmaster % section
        CALL Ereport( RoutineName, ErrorStatus, Cmessage )

      END SELECT

    CASE DEFAULT                          ! Selection on Internal Model
      ErrorStatus = 40
      WRITE (Cmessage, '(A, I2)')                                    &
            'No Field Calculations specified for Submodel ',         &
            fields_out( i ) % stashmaster % model
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )

    END SELECT                            ! Selection on Internal Model
    !-----------------------------------------------------------------
    ! Clean up and write out
    !-----------------------------------------------------------------
    CALL Rcf_Write_Field( fields_out( i ), Hdr_out, decomp_rcf_output )
    CALL Rcf_DeAlloc_Field( fields_out( i ) )
  END IF

END DO

IF ( ALLOCATED( i_dominant ) ) DEALLOCATE( i_dominant )
IF (LTimer) CALL Timer( RoutineName, 4)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Field_Calcs
END MODULE Rcf_Field_Calcs_Mod
