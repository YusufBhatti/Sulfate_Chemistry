! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Sets the data_source array corresponding to the fields array.

MODULE Rcf_Set_Data_Source_Mod

USE ereport_mod,            ONLY: ereport
USE rcf_field_type_mod,     ONLY: field_type
USE rcf_umhead_mod,         ONLY: um_header_type
USE rcf_locate_mod,         ONLY: rcf_locate
USE conversions_mod,        ONLY: zerodegc   

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

PRIVATE
PUBLIC :: rcf_initialise_data_source, rcf_reset_data_source

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SET_DATA_SOURCE_MOD'

CONTAINS

SUBROUTINE rcf_initialise_data_source(data_source, fields_in, fields_out,   &
                                      field_count_in, field_count_out,      &
                                      hdr_in, hdr_out )
! Description:
!   Allocates the data_source array and sets the elements according
!   to how the data in the the associated fields array should be
!   initialised based on choices made from ITEMS namelist.  Any items
!   not set in the ITEMS namelist are set to source=1 (input_dump). After
!   the ancil processing has been performed a second routine 
!   rcf_reset_data_source will be called to allow the reconfiguration 
!   logic to determine the data_source for certain fields.
!
!   If using the option to interpolate ALL fields which are found in 
!   the input dump (interp_all_fields) or using nearest neighbour interpolation
!   then  ALL data source values are overwritten by contents list of input_dump
!
! Method:
!   Uses data from the ITEMS namelist.

USE rcf_items_mod, ONLY:                &
     source_array,          upas_array, &
     sctn_array,            upaa_array, &
     item_array,            uprc_array, &
     area_array,            num_items,  &
     upaf_array                                      

USE rcf_data_source_mod,    ONLY: &
    data_source_type

USE items_nml_mod, ONLY: &
    input_dump,          &
    whole_grid 

USE rcf_interp_weights_mod, ONLY: h_int_method, nearest_neighbour   
USE umprintmgr,             ONLY: umPrint, umMessage
USE rcf_nlist_recon_technical_mod, ONLY: &
    interp_all_fields,                   &
    select_output_fields

IMPLICIT NONE

! Arguments
TYPE( data_source_type ), POINTER  :: data_source(:)
TYPE( field_type),   POINTER       :: fields_in(:)
TYPE( field_type ) , POINTER       :: fields_out(:)
TYPE( um_header_type ), INTENT(IN) :: hdr_in
TYPE( um_header_type ), INTENT(IN) :: hdr_out
INTEGER                            :: field_count_in
INTEGER                            :: field_count_out
! Local variables
INTEGER :: item_index
INTEGER :: i,j,pos
INTEGER :: ErrorStatus
CHARACTER (LEN=*), PARAMETER      :: RoutineName='RCF_INITIALISE_DATA_SOURCE'
CHARACTER (LEN=160)                :: Cmessage
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Check some of the arrays needed have been set up
IF ( ASSOCIATED( data_source ) ) THEN
  ErrorStatus = -10
  Cmessage = 'data_source is already set up!'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF
IF (.NOT. ASSOCIATED( fields_out )  ) THEN
  ErrorStatus = 20
  Cmessage = 'Fields have not yet been set up'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! Allocate the space
ALLOCATE( data_source( field_count_out ) )

! Loop over all fields which are due to be written to output dump
! and find any matching stash codes from the ITEMS namelist and 
! set the value of data_source
DO i = 1, field_count_out
  item_index = 0
  DO j = 1, num_items
    IF ( fields_out( i ) % stashmaster % item == item_array( j ) .AND.         &
         fields_out( i ) % stashmaster % section == sctn_array( j ) ) THEN
      item_index = j
      EXIT
    END IF
  END DO
  
  IF ( item_index > 0 ) THEN
    ! Have found the item - fill in gaps
    data_source( i ) % source      = source_array( item_index )
    data_source( i ) % Domain      = area_array  ( item_index )
    data_source( i ) % Ancil_SctnC = upas_array  ( item_index )
    data_source( i ) % Ancil_ItemC = upaa_array  ( item_index )
    data_source( i ) % RConst      = uprc_array  ( item_index )
    data_source( i ) % Ancil_File  = upaf_array  ( item_index )
  ELSE   ! No item found in namelist - use defaults
    data_source( i ) % source      = Input_Dump
    data_source( i ) % Domain      = Whole_Grid
    data_source( i ) % Ancil_SctnC = 0
    data_source( i ) % Ancil_ItemC = 0
    data_source( i ) % RConst      = 0.0
    data_source( i ) % Ancil_File  = ' '
  END IF

  ! Determine location of field in input dump (pos)
  CALL Rcf_Locate( fields_out( i ) % stashmaster % section,                    &
       fields_out( i ) % stashmaster % item,                                   &
       fields_in, field_count_in, pos, zero_ok_arg = .TRUE.)

  IF (select_output_fields == interp_all_fields) THEN
    IF (pos == 0) THEN
      WRITE(umMessage,'(A,I4,I3)') "Error processing STASHcode: ",             &
           fields_out( i ) % stashmaster % section,                            &
           fields_out( i ) % stashmaster % item
      CALL umPrint(umMessage,src='rcf_set_data_source_mod')
      ErrorStatus = 20
      Cmessage = 'Should be interpolating ALL input dump fields but input' // &
           ' dump is missing the field above.'
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    ELSE
      IF (data_source ( i ) % source /= Input_Dump ) THEN
        data_source( i ) % source = Input_Dump
      END IF
    END IF
  END IF
END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rcf_initialise_data_source

SUBROUTINE rcf_reset_data_source( data_source, fields_in, fields_out,        &
    field_count_in, field_count_out,                                         &
    hdr_in, hdr_out )

! Description:
!   Resets the data_source for certain fields from "input_dump" to one of 
!   several options (field_calcs, set_to_mdi, set_to_zero, set_to_const).
!   The data_source of the fields is reset based on one of three scenarios:
!   i)   Action required if field cannot be found in input dump
!   ii)  Action required if field IS found in input dump
!   iii) Action required regardless of whether field is input dump or not

USE items_nml_mod, ONLY:                                                   &
    Input_Dump,                                                            &
    Ancillary_File,                                                        &
    Set_To_Zero,                                                           &
    Set_To_MDI,                                                            &
    Set_To_Const,                                                          &
    Field_Calcs,                                                           &
    Field_Dependent_Calcs

USE Rcf_Data_Source_Mod, ONLY:                                             &
    Already_Processed,                                                     &
    Other_Field,                                                           &
    data_source_type

USE Rcf_Interp_Weights_Mod, ONLY:                                          &
    h_int_method,                                                          &
    nearest_neighbour,                                                     &
    h_int_active

USE Rcf_V_Int_Ctl_Mod, ONLY:                                               &
    v_int_active

USE UM_ParCore, ONLY:                                                      &
    mype

USE umPrintMgr, ONLY:                                                      &
    umPrint,                                                               &
    umMessage,                                                             &
    PrintStatus,                                                           &
    PrStatus_Normal,                                                       &
    PrStatus_Oper

USE Submodel_Mod, ONLY:                                                    &
    atmos_im

USE rcf_nlist_recon_science_mod, ONLY: &
    l_force_relayer,                   &
    l_snow_tile_gbm_ancil_fix

USE rcf_nlist_recon_technical_mod, ONLY: &
    grib2ff_input_dump,                  &
    grib_input_dump,                     &
    input_dump_type,                     &
    l_rcf_init_flexi

USE jules_sea_seaice_mod, ONLY:                                            &
    nice

USE um_stashcode_mod   ! So many of these are used the whole lot
                       ! have been included. Anything starting
                       ! with stashcode_ came from here.

USE jules_snow_mod, ONLY: nsmax, r0

USE jules_surface_mod, ONLY: l_aggregate

USE lake_mod, ONLY:                                                        &
    lake_depth_0,                                                          &
    lake_fetch_0,                                                          &
    lake_shape_0,                                                          &
    g_dt_0

USE rcf_lsh_field_checks_mod, ONLY:                                        &
    rcf_lsh_field_checks

USE rcf_compare_tiles_mod, ONLY:                                           &
    rcf_compare_tiles

USE Lookup_addresses

USE cppxref_mod, ONLY:                                                     &
    ppx_atm_lbc_theta,                                                     &
    ppx_atm_lbc_u,                                                         &
    ppx_atm_lbc_v

USE Rcf_HeadAddress_Mod, ONLY:                                             &
    FH_GridStagger, FH_GridStagger_Endgame, FH_GridStagger_C

USE missing_data_mod,     ONLY: rmdi
USE decomp_params,        ONLY: decomp_rcf_input

USE rcf_alloc_field_mod, ONLY:                                             &
    rcf_alloc_field,                                                       &
    rcf_dealloc_field
USE rcf_read_field_mod, ONLY:                                              &
    rcf_read_field

USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_lam

USE switches_urban, ONLY : l_moruses_macdonald

USE rcf_nlist_recon_idealised_mod, ONLY:                                   &
    initial_profile,                                                       &
    l_init_idealised,                                                      &
    l_reset_mixing,                                                        &
    qprofile_number,                                                       &
    surface_type,                                                          &
    tprofile_number

USE idealise_run_mod, ONLY:                                                &
    l_spec_z0,                                                             &
    roughlen_z0m

USE rcf_ideal_surface_mod, ONLY:                                           &
    surface_dump,                                                          &
    surface_ancil,                                                         &
    surface_zero

USE rcf_ideal_qprofile_mod, ONLY: &
    qp_dump

USE rcf_ideal_tprofile_mod, ONLY: &
    tp_dump

USE rcf_ideal_initial_profiles_mod, ONLY: &
    baro_inst,                            &
    deep_baro_inst,                       &
    no_preset,                            &
    rot_solid_body

USE dynamics_testing_mod, ONLY: &
    l_dry,                      &
    problem_number

USE problem_mod, ONLY: &
    dynamical_core,    &
    idealised_problem, &
    idealised_planet

USE missing_data_mod, ONLY: &
    imdi

IMPLICIT NONE

! Arguments
TYPE( data_source_type ), POINTER  :: data_source(:)
TYPE( field_type), POINTER         :: fields_in(:)
TYPE( field_type ) , POINTER       :: fields_out(:)
TYPE( um_header_type ), INTENT(IN) :: hdr_in
TYPE( um_header_type ), INTENT(IN) :: hdr_out
INTEGER                            :: field_count_in
INTEGER                            :: field_count_out

! Local variables
LOGICAL                            :: l_match=.TRUE.
LOGICAL                            :: l_rcf_compare_tiles
INTEGER                            :: i
INTEGER                            :: j
INTEGER                            :: pos
INTEGER                            :: pos_in
INTEGER                            :: pos_out
INTEGER                            :: ErrorStatus
CHARACTER (LEN=*), PARAMETER       :: RoutineName='RCF_RESET_DATA_SOURCE'
CHARACTER (LEN=errormessagelength)        :: Cmessage

CHARACTER (LEN=*), PARAMETER    :: FORM                                    &
    ="(2(i5,' '),4(i4,' '),e11.4,' ',50a)"
CHARACTER (LEN=*), PARAMETER    :: form_c=                                 &
    "(2(a5,' '),4(a4,' '),a10,  ' ',50a)"

LOGICAL                            :: eg2nd
TYPE( field_type ), POINTER        :: dryrho_in
REAL :: rmdi_tol
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! --------------------------------------------------------------
! ENDGame to New Dynamics RCF requires special treatment
! as adv wind fields are not used in EG while required for ND
! set up logical to identify this case.
! --------------------------------------------------------------
eg2nd = .FALSE.

IF ( (Hdr_In % FixHd(FH_GridStagger) == FH_GridStagger_Endgame)            &
.AND. ( Hdr_Out % FixHd(FH_GridStagger) == FH_GridStagger_C )) THEN
  eg2nd = .TRUE.
END IF

! Loop over all fields due to be written to the output dump
DO i = 1, field_count_out

  CALL Rcf_Locate( fields_out( i ) % stashmaster % section,                &
      fields_out( i ) % stashmaster % item,                                &
      fields_in, field_count_in, pos, zero_ok_arg = .TRUE.)

  ! --------------------------------------------------------------------
  ! Some fields may need special treatment whether they are in the
  ! input dump or not - these are checked first before their data
  ! sources are potentially altered later in the routine.
  ! --------------------------------------------------------------------

  IF ( data_source( i ) % source == Input_Dump ) THEN

    SELECT CASE( fields_out( i ) % stashmaster % model )

    CASE ( atmos_im )                 ! Atmosphere Sub-Model

      ! ----------------
      ! Atmosphere items
      ! ----------------

      SELECT CASE( fields_out( i ) % stashmaster % section )

      CASE ( stashcode_prog_sec )

        SELECT CASE( fields_out( i ) % stashmaster % item )

        CASE (stashcode_orog_x_grad,                                       &
              stashcode_orog_y_grad,                                       &
              stashcode_unfilt_orog,                                       &
              stashcode_sil_orog_rough,                                    &
              stashcode_hlf_pk_to_trf_ht,                                  &
              stashcode_orog_var,                                          &
              stashcode_orog_gdxx,                                         &
              stashcode_orog_gdxy,                                         &
              stashcode_orog_gdyy)       ! Orography-derived fields

          ! Check if idealised initialisation
          IF (l_init_idealised) THEN
            ! Select on idealised surface type
            SELECT CASE (surface_type)

            ! Orography from dump or ancil
            CASE (surface_dump,                                            &
                  surface_ancil)

              ! Do nothing

            ! Flat orography
            CASE (surface_zero)

              data_source( i ) % source = set_to_zero

            ! Any other prescribed ideal orography
            CASE DEFAULT

              data_source( i ) % source = set_to_mdi

            END SELECT                   ! based on surface type
          END IF                         ! if idealised


        CASE (stashcode_u)
          ! Check if idealised initialisation
          IF (l_init_idealised) THEN
            ! Check type of idealised problem
            SELECT CASE(problem_number)
            CASE (dynamical_core)
              IF (initial_profile == no_preset) THEN
                IF (tprofile_number /= imdi) THEN
                  data_source( i ) % source = set_to_zero
                END IF
              END IF
            CASE (idealised_planet)
              data_source( i ) % source = set_to_zero
            END SELECT
          END IF

        CASE (stashcode_v)
          ! Check if idealised initialisation
          IF (l_init_idealised) THEN
            ! Check type of idealised problem
            SELECT CASE(problem_number)
            CASE (dynamical_core)
              SELECT CASE (initial_profile)
              CASE (baro_inst)
                data_source( i ) % source = set_to_zero
              CASE (no_preset)
                IF (tprofile_number /= imdi) THEN
                  data_source( i ) % source = set_to_zero
                END IF
              END SELECT
            CASE (idealised_planet)
              data_source( i ) % source = set_to_zero
            END SELECT
          END IF

        CASE (stashcode_w)
          ! Check if idealised initialisation
          IF (l_init_idealised) THEN
            ! Check type of idealised problem
            SELECT CASE(problem_number)
            CASE (dynamical_core)
              SELECT CASE (initial_profile)
              CASE (baro_inst, rot_solid_body, deep_baro_inst)
                data_source( i ) % source = set_to_zero
              CASE (no_preset)
                IF (tprofile_number /= imdi) THEN
                  data_source( i ) % source = set_to_zero
                END IF
              END SELECT
            CASE (idealised_problem)
              data_source( i ) % source = set_to_zero
            CASE (idealised_planet)
              data_source( i ) % source = set_to_zero
            END SELECT
          END IF

        CASE (stashcode_u_adv)
          ! Check if idealised initialisation
          IF (l_init_idealised) THEN
            ! Check type of idealised problem
            SELECT CASE(problem_number)
            CASE (dynamical_core)
              IF (initial_profile == no_preset) THEN
                IF (tprofile_number /= imdi) THEN
                  data_source( i ) % source = set_to_zero
                END IF
              END IF
            CASE (idealised_planet)
              data_source( i ) % source = set_to_zero
            END SELECT
          END IF

        CASE (stashcode_v_adv)
          ! Check if idealised initialisation
          IF (l_init_idealised) THEN
            ! Check type of idealised problem
            SELECT CASE(problem_number)
            CASE (dynamical_core)
              SELECT CASE (initial_profile)
              CASE (baro_inst)
                data_source( i ) % source = set_to_zero
              CASE (no_preset)
                IF (tprofile_number /= imdi) THEN
                  data_source( i ) % source = set_to_zero
                END IF
              END SELECT
            CASE (idealised_planet)
              data_source( i ) % source = set_to_zero
            END SELECT
          END IF
        CASE (stashcode_w_adv)
          ! Check if idealised initialisation
          IF (l_init_idealised) THEN
            ! Check type of idealised problem
            SELECT CASE(problem_number)
            CASE (dynamical_core)
              SELECT CASE (initial_profile)
              CASE (baro_inst, rot_solid_body, deep_baro_inst)
                 data_source( i ) % source = set_to_zero
              CASE (no_preset)
                IF (tprofile_number /= imdi) THEN
                   data_source( i ) % source = set_to_zero
                END IF
              END SELECT
            CASE (idealised_planet)
              data_source( i ) % source = set_to_zero
            END SELECT
          END IF

        CASE (stashcode_exner_surf)
          ! Check if idealised initialisation
          IF (l_init_idealised) THEN
            IF (tprofile_number /= tp_dump) THEN
              ! Idealised initialisation begins with calc. of exner_surf
              data_source( i ) % source = Field_Dependent_Calcs
            END IF
          END IF

        CASE (stashcode_mv,                                                &
              stashcode_mcl,                                               &
              stashcode_mcf,                                               &
              stashcode_mr,                                                &
              stashcode_mgr,                                               &
              stashcode_mcf2)
          ! Check if idealised initialisation
          IF (l_init_idealised) THEN
            ! Check type of idealised problem
            SELECT CASE(problem_number)
            CASE (idealised_problem)
              IF (l_reset_mixing .OR. l_dry) THEN
                data_source( i ) % source = set_to_zero
              END IF
            CASE (idealised_planet)
              IF (qprofile_number /= qp_dump) THEN
                data_source( i ) % source = set_to_zero
              END IF
            END SELECT
          END IF

        CASE (stashcode_q,   &
              stashcode_qcl, &
              stashcode_qcf, &
              stashcode_qrain,&
              stashcode_qgraup,&
              stashcode_qcf2)
          ! Check if idealised initialisation
          IF (l_init_idealised) THEN
            ! Check type of idealised problem
            SELECT CASE(problem_number)
            CASE (idealised_problem)
              IF (l_reset_mixing .OR. l_dry) THEN
                data_source( i ) % source = set_to_mdi
              END IF
            CASE (idealised_planet)
              IF (qprofile_number /= qp_dump) THEN
                data_source( i ) % source = set_to_mdi
              END IF
            END SELECT
          END IF

        CASE (stashcode_cca,                                               &
              stashcode_ccb,                                               &
              stashcode_cct,                                               &
              stashcode_cc_lwp,                                            &
              stashcode_surf_z_curr,                                       &
              stashcode_surf_m_curr,                                       &
              stashcode_n_turb_mixlvs,                                     &
              stashcode_lvl_bse_dp_sc,                                     &
              stashcode_lvl_top_dp_sc,                                     &
              stashcode_bl_conv_flag,                                      &
              stashcode_turb_temp,                                         &
              stashcode_turb_humid,                                        &
              stashcode_area_cf,                                           &
              stashcode_bulk_cf,                                           &
              stashcode_liquid_cf,                                         &
              stashcode_frozen_cf,                                         &
              stashcode_sfc_zonal_cur,                                     &
              stashcode_sfc_merid_cur )
          IF (l_init_idealised) THEN
            data_source( i ) % source = set_to_zero
          END IF

        CASE (stashcode_bl_depth)
          IF (l_init_idealised) THEN
            data_source( i ) % source = set_to_const
            data_source( i ) % rconst = 100.0
          END IF

        CASE (stashcode_rough_length)
          IF (l_init_idealised) THEN
            IF (l_spec_z0) THEN
              data_source( i ) % source = set_to_const
              data_source( i ) % rconst = roughlen_z0m
            END IF
          END IF

        END SELECT                       ! based on STASH item code

      END SELECT                         ! based on STASH section code

    CASE DEFAULT                 ! Couldn't find a proper Internal Model

      ! The Model type didn't match one of the three defined types.
      WRITE (Cmessage, '(2A, I2, A, I3, A, I5)')                           &
          " Couldn't find a Sub-Model ID type for : ",                     &
          " Model ", fields_out( i ) % stashmaster % model,                &
          " Section ", fields_out( i ) % stashmaster % section,            &
          " Item ",  fields_out( i ) % stashmaster % item
      ErrorStatus = 80
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )

    END SELECT                   ! select by Internal model

  END IF

  ! Setting Source for items required from input dump but not
  ! available (and known about with relevant code written)
  IF ( pos == 0 .AND. data_source( i ) % source == Input_Dump ) THEN

    ! ----------------------------------
    ! This item is not in the input dump
    ! ----------------------------------

    SELECT CASE( fields_out( i ) % stashmaster % model )

    CASE ( atmos_im )                 ! Atmosphere Sub-Model

      ! ----------------
      ! Atmosphere items
      ! ----------------

      SELECT CASE( fields_out( i ) % stashmaster % section )

      CASE ( stashcode_prog_sec )

        SELECT CASE( fields_out( i ) % stashmaster % item )

        CASE ( stashcode_3d_cca )           ! 3D Convective cloud

          data_source( i ) % source = set_to_zero

        CASE ( stashcode_cca )              ! 2D Convective cloud

          data_source( i ) % source = Field_Calcs

        CASE ( stashcode_npp_pft_acc,     & ! Carbon accumulation fields
            stashcode_g_lf_pft_acc,                                        &
            stashcode_g_ph_lf_pft_acc,                                     &
            stashcode_rsp_w_pft_acc,                                       &
            stashcode_rsp_s_acc,                                           &
            stashcode_catch_snow,      & ! NLT canopy snow capacity
            stashcode_catch_tile,      & ! Tiled canopy capacity
            stashcode_z0_tile,         & ! Tiled roughness length
            stashcode_z0h_tile         & ! Tiled thermal roughness length
            )
          data_source( i ) % source = set_to_const
          data_source( i ) % rconst = -1.0

        CASE ( stashcode_can_water_tile,  & ! Vegetation fields
            stashcode_tstar_tile,                                          &
            stashcode_tsurf_elev_surft,                                    &
            stashcode_ice_temp_cat,                                        &
            stashcode_ice_conc_cat,      & ! ice conc (fraction) cats
            stashcode_ice_thick_cat,     & ! ice thickness categories
            stashcode_Gamtot,            & !\ large-scale hydrology
            stashcode_Zw,                & !/ fields
            stashcode_Sthzw,             & !/
            stashcode_Fsat,              & !\.
            stashcode_Fwetl,             & !/.
            stashcode_a_fsat,            & !\.
            stashcode_c_fsat,            & !/.
            stashcode_a_fwet,            & !\.
            stashcode_c_fwet,            & !/.
            stashcode_snow_tile,         &
            stashcode_flake_t_mean,      & ! FLake prognostics
            stashcode_flake_t_mxl,       & ! "
            stashcode_flake_t_ice,       & ! "
            stashcode_flake_h_mxl,       & ! "
            stashcode_flake_h_ice        & ! "
            )

          data_source( i ) % source = Field_Calcs

        CASE ( stashcode_snowdep_grd_tile,  & ! JULES snowdepth on tiles
            stashcode_snowpack_bk_dens,  & ! JULES snowpack bulk density
            stashcode_nsnow_layrs_tiles, & ! JULES number of snow layers
            stashcode_snow_laythk_tiles, & ! JULES snow layer thicknesses
            stashcode_snow_ice_tile,     & ! JULES snow layer solid mass
            stashcode_snow_liq_tile,     & ! JULES snow layer liquid mass
            stashcode_snow_T_tile,       & ! JULES snow layer temps.
            stashcode_snow_laydns_tiles  & ! JULES snow layer density
            )

          data_source( i ) % source = Field_Dependent_Calcs
!
        CASE ( stashcode_flake_depth )  ! FLake lake depth

          data_source( i ) % source = set_to_const
          data_source( i ) % rconst = lake_depth_0

        CASE ( stashcode_flake_fetch )  ! FLake lake fetch

          data_source( i ) % source = set_to_const
          data_source( i ) % rconst = lake_fetch_0

        CASE ( stashcode_flake_shape )  ! FLake thermocline shape factor

          data_source( i ) % source = set_to_const
          data_source( i ) % rconst = lake_shape_0

        CASE ( stashcode_flake_g_over_dt ) ! (ht flux/delta Tsurf from FLake)

          data_source( i ) % source = set_to_const
          data_source( i ) % rconst = g_dt_0

        CASE ( stashcode_ice_surf_temp_cat,   &
               stashcode_tstar_ice_cat_cpl )  ! Sea ice surface temperature 
          !                                               !on categories 
          data_source( i ) % source = Field_calcs 

        CASE ( stashcode_fexp          & ! LSH conductivity decay scale
            )
          data_source( i ) % source = set_to_const
          data_source( i ) % rconst = 1.0

        CASE ( stashcode_infil_max_tile,                                   &
            stashcode_snow_grnd,       & ! Snow beneath NLT canopy
            stashcode_snow_on_ice,                                         &
            stashcode_ice_snow_depth,                                      &
            stashcode_sw_down_tile,                                        &
            stashcode_sw_down,                                             &
            stashcode_lw_up_diff,                                          &
            stashcode_mmr_smoke_fr,                                        &
            stashcode_mmr_smoke_ag,                                        &
            stashcode_mmr_smoke_cl,                                        &
            stashcode_biom_surf_em,                                        &
            stashcode_biom_elev_em,                                        &
            stashcode_biom_elev_em_h1, & ! min inject height for biom_elev_em
            stashcode_biom_elev_em_h2, & ! max inject height for biom_elev_em
            stashcode_dust1_mmr,                                           &
            stashcode_dust2_mmr,                                           &
            stashcode_dust3_mmr,                                           &
            stashcode_dust4_mmr,                                           &
            stashcode_dust5_mmr,                                           &
            stashcode_dust6_mmr,                                           &
            stashcode_ddmfx,                                               &
            stashcode_albedo_sice,                                         &
            stashcode_albedo_land,                                         &
            stashcode_soilcarb_dpm,    & ! RothC soil C prognostic
            stashcode_soilcarb_rpm,    & ! RothC soil C prognostic
            stashcode_soilcarb_bio,    & ! RothC soil C prognostic
            stashcode_soilcarb_hum,    & ! RothC soil C prognostic
            stashcode_soilnitro_dpm,   & ! RothC soil N prognostic
            stashcode_soilnitro_rpm,   & ! RothC soil N prognostic
            stashcode_soilnitro_bio,   & ! RothC soil N prognostic
            stashcode_soilnitro_hum,   & ! RothC soil N prognostic
            stashcode_soil_inorgnit,   & ! RothC soil N prognostic
            stashcode_qcf2,            & ! second cloud ice prognostic
            stashcode_qrain,           & ! prognostic rain
            stashcode_qgraup,          & ! prognostic graupel
            stashcode_lcbase,          & ! lowest convective cloud base
            stashcode_3d_ccw,          & ! CCW profile sent to radiation
            stashcode_deep_conv_flag,  & ! deep convection flag
            stashcode_past_conv_precip,& ! past convective precipitation
            stashcode_past_conv_depth, & ! past convective depth
            stashcode_cca_dp,          & ! CCA from deep convection
            stashcode_cca_md,          & ! CCA from mid-level convection
            stashcode_cca_sh,          & ! CCA from shallow convection
            stashcode_total_precip,    & ! Total precipitation rate
            stashcode_flash_pot,       & ! Lightning flash potential
            stashcode_blwvariance,     & ! w-variance from bl scheme
            stashcode_cloud_number,    & ! CASIM Cloud Number
            stashcode_rain_number,     & ! CASIM Rain Number
            stashcode_rain_3mom,       & ! CASIM Rain 3rd moment
            stashcode_ice_number,      & ! CASIM Ice Number
            stashcode_snow_number,     & ! CASIM Snow Number
            stashcode_snow_3mom,       & ! CASIM Snow 3rd moment
            stashcode_graup_number,    & ! CASIM Graupel Number
            stashcode_graup_3mom,      & ! CASIM Graupel 3rd Moment
            stashcode_activesol_liquid,& ! CASIM Activated soluable 
                                         ! aerosol in liquid
            stashcode_activesol_rain,  & ! CASIM Activated soluable 
                                         ! aerosol in rain
            stashcode_active_insol_ice,& ! CASIM Activated 
                                         ! insoluable aerosol in ice
            stashcode_active_sol_ice,  & ! CASIM Activated soluable 
                                         ! aerosol in ice
            stashcode_active_insol_liq,& ! CASIM Activated insoluable 
                                         ! aerosol in liquid
            stashcode_active_sol_num,  & ! CASIM Activated soluable 
                                         ! aerosol number
            stashcode_active_insol_num,& ! CASIM Activated insoluable 
                                         ! aerosol number
            stashcode_ice_subl_cat,    & ! Sea ice sublimation on categories 
            stashcode_ice_surf_cond_cat, & ! Sea ice surface layer 
                                           ! conductivity on categories 
            stashcode_pond_frac_cat, &  ! Sea ice meltpond fraction on 
                                        ! categories
            stashcode_pond_depth_cat, &  ! Sea ice meltpond depth
            stashcode_dPV_rad,         & ! PV-tracer scheme: dPV for rad
            stashcode_dPV_sw,          & ! dPV for SW rad
            stashcode_dPV_lw,          & ! dPV for LW rad
            stashcode_dPV_mic,         & ! dPV for microphysics
            stashcode_dPV_gwd,         & ! dPV for GW drag
            stashcode_dPV_ph1,         & ! dPV for slow physics 
            stashcode_dPV_conv,        & ! dPV for convection
            stashcode_dPV_bl,          & ! dPV for Boundary Layer
            stashcode_dPV_stph,        & ! dPV for stochastic physics
            stashcode_dPV_cld,         & ! dPV for cloud re-balance
            stashcode_dPV_iau,         & ! dPV for IAU
            stashcode_dPV_nud,         & ! dPV for nudging
            stashcode_dPV_tot,         & ! Total dPV per timestep
            stashcode_dPV_adv,         & ! dPV advection
            stashcode_dPV_sol,         & ! dPV solver
            stashcode_dPV_mass,        & ! Mass update
            stashcode_adv_only_PV,     & ! Advection only dPV
                                       ! End of PV-tracer progs
            stashcode_dtheta_0,        & ! dtheta_0
            stashcode_dtheta_bl,       & ! dtheta_bl
            stashcode_dtheta_bl_mix,   & ! dtheta_bl_mix
            stashcode_dtheta_bl_LH,    & ! dtheta_bl_LH
            stashcode_dtheta_conv,     & ! dtheta_conv
            stashcode_dtheta_mic,      & ! dtheta_mic
            stashcode_dtheta_rad,      & ! dtheta_rad
            stashcode_dtheta_SW,       & ! dtheta_SW
            stashcode_dtheta_LW,       & ! dtheta_LW
            stashcode_dtheta_slow,     & ! dtheta_slow
            stashcode_dtheta_cld,      & ! dtheta_cld            
                                         ! End of diabatic-tracer progs
            stashcode_bl_pert_rand_fld, & ! Stoch phys random fld for BL pert
            stashcode_bl_pert_flag,     & ! Stoch phys flag for BL pert
            stashcode_ux_ccp,          & ! ccp x cmpt frnt speed vector sum
            stashcode_uy_ccp,          & ! ccp y cmpt frnt speed vector sum
            stashcode_um_ccp,          & ! ccp front speed scalar sum
            stashcode_g_ccp,           & ! ccp gridbox reduced gravity
            stashcode_h_ccp,           & ! ccp gridbox depth
            stashcode_riso_ccp,        & ! ccp remain counter (isotropic)
            stashcode_rdir_ccp         & ! ccp remain counter (directed)
            )

          data_source( i ) % source = set_to_zero

        CASE ( stashcode_rgrain,                                           &
            stashcode_snow_grnsiz_tiles & ! JULES snow layer grain size
            )

          data_source( i ) % source = set_to_const
          data_source( i ) % rconst = r0

        CASE ( stashcode_tstar_land,                                       &
               stashcode_tstar_sea,                                        &
               stashcode_tstar_anom,                                       &
               stashcode_tstar_sice )

          data_source( i ) % source = Already_Processed

        CASE ( stashcode_dctemp_tile,                                      &
               stashcode_dctemp_ssi,                                       &
               stashcode_tm_trans,                                         &
               stashcode_e_trb,                                            &
               stashcode_tsq_trb,                                          &
               stashcode_qsq_trb,                                          &
               stashcode_cov_trb,                                          &
               stashcode_zhpar_shcu )

          data_source( i ) % source = set_to_mdi

        CASE ( stashcode_sstfrz )
          data_source( i ) % source = set_to_const
          data_source( i ) % rconst = zerodegc

        CASE ( stashcode_MV,   stashcode_MCL,                              &
               stashcode_MCF,  stashcode_MR,                               &
               stashcode_MGR,  stashcode_MCF2 )

          data_source( i ) % source = Field_Calcs

        CASE ( stashcode_THETAVD, stashcode_DRY_RHO,                       &
               stashcode_exner_surf )

          data_source( i ) % source = Field_Dependent_Calcs

        CASE ( stashcode_etadot )     
          IF (model_type == mt_lam) THEN
            data_source( i ) % source = Set_To_Zero
          ELSE
            data_source( i ) % source = Field_Dependent_Calcs
          END IF

        CASE ( stashcode_psiw_lid,                                         &
               stashcode_psiw_surf )

          data_source( i ) % source = Set_To_Mdi

        CASE( stashcode_urbdisp,        & ! URBAN displacement height
              stashcode_urbztm )          ! URBAN effective roughness length

          ! l_urban2t should be true if these stash are required
          IF ( l_moruses_macdonald ) THEN
            data_source( i ) % source = Field_Calcs
          ELSE
            errorstatus = 10
            WRITE (Cmessage, '(A, 2I5)')                                   &
               'l_moruses_macdonald = F: Data source not set for',         &
               fields_out( i ) % stashmaster % item
            CALL Ereport( RoutineName, ErrorStatus, Cmessage )
          END IF

        CASE( stashcode_urbalbwl )       ! URBAN wall albedo

          data_source( i ) % source = set_to_const
          data_source( i ) % rconst = 0.375

        CASE( stashcode_urbalbrd )       ! URBAN road albedo

          data_source( i ) % source = set_to_const
          data_source( i ) % rconst = 0.08

        CASE( stashcode_urbemisw )       ! URBAN wall emissivity

          data_source( i ) % source = set_to_const
          data_source( i ) % rconst = 0.875

        CASE( stashcode_urbemisr )       ! URBAN road emissivity

          data_source( i ) % source = set_to_const
          data_source( i ) % rconst = 0.95

        END SELECT                       ! based on STASH item code

     ! GLOMAP_CLIM related items. 
     CASE ( stashcode_glomap_clim_sec ) 
        ! Select on item number
        SELECT CASE( (fields_out( i ) % stashmaster % section)*1000 +  &
                      fields_out( i ) % stashmaster % item )
           
        CASE (stashcode_gc_nd_nuc_sol,                                 &
              stashcode_gc_nuc_sol_su,                                 &
              stashcode_gc_nuc_sol_oc,                                 &
              stashcode_gc_nd_ait_sol,                                 &
              stashcode_gc_ait_sol_su,                                 &
              stashcode_gc_ait_sol_bc,                                 &
              stashcode_gc_ait_sol_oc,                                 &
              stashcode_gc_nd_acc_sol,                                 &
              stashcode_gc_acc_sol_su,                                 &
              stashcode_gc_acc_sol_bc,                                 &
              stashcode_gc_acc_sol_oc,                                 &
              stashcode_gc_acc_sol_ss,                                 &
              stashcode_gc_nd_cor_sol,                                 &
              stashcode_gc_cor_sol_su,                                 &
              stashcode_gc_cor_sol_bc,                                 &
              stashcode_gc_cor_sol_oc,                                 &
              stashcode_gc_cor_sol_ss,                                 &
              stashcode_gc_nd_ait_ins,                                 &
              stashcode_gc_ait_ins_bc,                                 &
              stashcode_gc_ait_ins_oc,                                 &
              stashcode_gc_dryd_ait_sol,                               &
              stashcode_gc_dryd_acc_sol,                               &
              stashcode_gc_dryd_cor_sol,                               &
              stashcode_gc_dryd_ait_ins,                               &
              stashcode_gc_dryd_acc_ins,                               &
              stashcode_gc_dryd_cor_ins,                               &
              stashcode_gc_wetd_ait_sol,                               &
              stashcode_gc_wetd_acc_sol,                               &
              stashcode_gc_wetd_cor_sol,                               &
              stashcode_gc_rho_ait_sol,                                &
              stashcode_gc_rho_acc_sol,                                &
              stashcode_gc_rho_cor_sol,                                &
              stashcode_gc_rho_ait_ins,                                &
              stashcode_gc_rho_acc_ins,                                &
              stashcode_gc_rho_cor_ins,                                &
              stashcode_gc_pvol_ait_su_sol,                            &
              stashcode_gc_pvol_ait_bc_sol,                            &
              stashcode_gc_pvol_ait_oc_sol,                            &
              stashcode_gc_pvol_ait_so_sol,                            &
              stashcode_gc_pvol_ait_h2o_sol,                           &
              stashcode_gc_pvol_acc_su_sol,                            &
              stashcode_gc_pvol_acc_bc_sol,                            &
              stashcode_gc_pvol_acc_oc_sol,                            &
              stashcode_gc_pvol_acc_ss_sol,                            &
              stashcode_gc_pvol_acc_no3_sol,                           &
              stashcode_gc_pvol_acc_du_sol,                            &
              stashcode_gc_pvol_acc_so_sol,                            &
              stashcode_gc_pvol_acc_h2o_sol,                           &
              stashcode_gc_pvol_cor_su_sol,                            &
              stashcode_gc_pvol_cor_bc_sol,                            &
              stashcode_gc_pvol_cor_oc_sol,                            &
              stashcode_gc_pvol_cor_ss_sol,                            &
              stashcode_gc_pvol_cor_no3_sol,                           &
              stashcode_gc_pvol_cor_du_sol,                            &
              stashcode_gc_pvol_cor_so_sol,                            &
              stashcode_gc_pvol_cor_h2o_sol,                           &
              stashcode_gc_pvol_ait_bc_ins,                            &
              stashcode_gc_pvol_ait_oc_ins,                            &
              stashcode_gc_pvol_acc_du_ins,                            &
              stashcode_gc_pvol_cor_du_ins,                            &
              stashcode_gc_cdnc3,                                      &
              stashcode_gc_cdnc )
           
          
           ! These are all set to zero if there is no other data available
           data_source( i ) % source = set_to_zero
           
        END SELECT                       ! based on STASH item code
         
      ! UKCA related items. 
      CASE ( stashcode_ukca_sec )

        ! Select on item number
        SELECT CASE( (fields_out( i ) % stashmaster % section)*1000 +  &
                      fields_out( i ) % stashmaster % item )

        CASE (stashcode_surfarea,                                      &
                stashcode_cdnc3,                                       &
                stashcode_cdnc,                                        &
                stashcode_ho2s,                                        &
                stashcode_ohs,                                         &
                stashcode_o1ds,                                        &
                stashcode_o3ps,                                        &
                stashcode_het_ho2,                                     &
                stashcode_het_n2o5,                                    &
                stashcode_tolp1,                                       &
                stashcode_hoipo2,                                      &
                stashcode_homvko2,                                     &
                stashcode_memald1,                                     &
                stashcode_oxyl1,                                       &
                stashcode_hoc3h6o2,                                    &
                stashcode_hoc2h4o2,                                    &
                stashcode_meko2,                                       &
                stashcode_mecoch2oo,                                   &
                stashcode_mecoc2oo,                                    &
                stashcode_etco3,                                       &
                stashcode_iproo,                                       &
                stashcode_sbuoo,                                       &
                stashcode_nproo,                                       &
                stashcode_meco3,                                       &
                stashcode_etoo,                                        &
                stashcode_meoo,                                        &
                stashcode_hcl_unlmp,                                   &
                stashcode_ho2_ntp,                                     &
                stashcode_bro_unlmp,                                   &
                stashcode_oh_ntp,                                      &
                stashcode_no2_unlmp,                                   &
                stashcode_o1d_ntp,                                     &
                stashcode_o3p_ntp,                                     &
        ! RADAER
                stashcode_dryd_ait_sol,                                &
                stashcode_dryd_acc_sol,                                &
                stashcode_dryd_cor_sol,                                &
                stashcode_dryd_ait_insol,                              &
                stashcode_dryd_acc_insol,                              &
                stashcode_dryd_cor_insol,                              &
                stashcode_wetd_ait_sol,                                &
                stashcode_wetd_acc_sol,                                &
                stashcode_wetd_cor_sol,                                &
                stashcode_rho_ait_sol,                                 &
                stashcode_rho_acc_sol,                                 &
                stashcode_rho_cor_sol,                                 &
                stashcode_rho_ait_insol,                               &
                stashcode_rho_acc_insol,                               &
                stashcode_rho_cor_insol,                               &
                stashcode_pvol_ait_su_sol,                             &
                stashcode_pvol_ait_bc_sol,                             &
                stashcode_pvol_ait_oc_sol,                             &
                stashcode_pvol_ait_so_sol,                             &
                stashcode_pvol_ait_h2o_sol,                            &
                stashcode_pvol_acc_su_sol,                             &
                stashcode_pvol_acc_bc_sol,                             &
                stashcode_pvol_acc_oc_sol,                             &
                stashcode_pvol_acc_ss_sol,                             &
                stashcode_pvol_acc_so_sol,                             &
                stashcode_pvol_acc_du_sol,                             &
                stashcode_pvol_acc_no3_sol,                            &
                stashcode_pvol_acc_h2o_sol,                            &
                stashcode_pvol_cor_su_sol,                             &
                stashcode_pvol_cor_bc_sol,                             &
                stashcode_pvol_cor_oc_sol,                             &
                stashcode_pvol_cor_ss_sol,                             &
                stashcode_pvol_cor_so_sol,                             &
                stashcode_pvol_cor_du_sol,                             &
                stashcode_pvol_cor_no3_sol,                            &
                stashcode_pvol_cor_h2o_sol,                            &
                stashcode_pvol_ait_bc_insol,                           &
                stashcode_pvol_ait_oc_insol,                           &
                stashcode_pvol_acc_du_insol,                           &
                stashcode_pvol_cor_du_insol )

          ! These are all set to zero if there is no other data available
          data_source( i ) % source = set_to_zero

        END SELECT                       ! based on STASH item code

      END SELECT                         ! based on STASH section code

      ! adv winds may not be included in EG dumps as they are not used.
      ! So when going EG --> ND we need to initialise them

      IF (eg2nd) THEN

        SELECT CASE( fields_out( i ) % stashmaster % section )

        CASE ( stashcode_prog_sec )
          SELECT CASE( fields_out( i ) % stashmaster % item )

          CASE ( stashcode_u_adv,                                          &
                 stashcode_v_adv )

            data_source( i ) % source = Field_Calcs

          CASE ( stashcode_w_adv )

            data_source( i ) % source = Set_To_Zero

          END SELECT
        END SELECT
      END IF


      !For files not found when reading from GRIB data
      IF ( input_dump_type == grib_input_dump .OR.                         &
           input_dump_type == grib2ff_input_dump ) THEN

        SELECT CASE( fields_out( i ) % stashmaster % section )

        CASE ( stashcode_prog_sec )

          SELECT CASE( fields_out( i ) % stashmaster % item )

          CASE ( stashcode_sea_ice_temp,                                   &
                 stashcode_icethick,                                       &
                 stashcode_u_adv,                                          &
                 stashcode_v_adv,                                          &
                 stashcode_bulk_cf,                                        &
                 stashcode_liquid_cf,                                      &
                 stashcode_frozen_cf )

            data_source( i ) % source = Field_Calcs

          CASE ( stashcode_w_adv,                                          &
                 stashcode_qcf,                                            &
                 stashcode_cc_lwp,                                         &
                 stashcode_unfrozen_soil,                                  &
                 stashcode_frozen_soil,                                    &
                 stashcode_qcl,                                            &
                 stashcode_n_turb_mixlvs,                                  &
                 stashcode_lvl_bse_dp_sc,                                  &
                 stashcode_lvl_top_dp_sc,                                  &
                 stashcode_bl_conv_flag,                                   &
                 stashcode_turb_temp,                                      &
                 stashcode_turb_humid,                                     &
                 stashcode_area_cf,                                        &
                 stashcode_sfc_zonal_cur,                                  &
                 stashcode_sfc_merid_cur,                                  &
                 stashcode_3d_cca,           & ! has to overide source=8
                 stashcode_can_water_tile,                                 &
                 stashcode_rho,              & ! rho calc'd after interp
                 stashcode_exner,                                          &
                 stashcode_can_conduct,                                    &
                 stashcode_ccb,                                            &
                 stashcode_cct,                                            &
                 stashcode_mean_canopyw,                                   &
                 stashcode_surf_z_curr,                                    &
                 stashcode_surf_m_curr,                                    &
                 stashcode_w )

            data_source( i ) % source = Set_To_Zero

          CASE ( stashcode_bl_depth )

            data_source( i ) % source = set_to_const
            data_source( i ) % rconst = 500.000

          CASE ( stashcode_z0 )

            data_source( i ) % source = set_to_const
            data_source( i ) % rconst = 0.500
          CASE ( stashcode_theta )
            data_source( i ) % source = other_field

          END SELECT                 ! select by item code

        END SELECT                   ! select by section code

      END IF                       ! If GRIB section

    CASE DEFAULT                 ! Couldn't find a proper Internal Model

      ! The Model type didn't match one of the three defined types.
      WRITE (Cmessage, '(2A, I2, A, I3, A, I5)')                           &
          " Couldn't find a Sub-Model ID type for : ",                     &
          " Model ", fields_out( i ) % stashmaster % model,                &
          " Section ", fields_out( i ) % stashmaster % section,            &
          " Item ",  fields_out( i ) % stashmaster % item
      ErrorStatus = 25
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )

    END SELECT                   ! select by Internal model

    ! Check that Source is now set correctly otherwise, fail
    IF ( data_source( i ) % source == Input_Dump ) THEN
      WRITE ( Cmessage, '(A, I3, A, I5, A)')                               &
          'Section ',                                                      &
          fields_out( i ) % stashmaster % section,                         &
          ' Item ',                                                        &
          fields_out( i ) % stashmaster % item ,                           &
          ' : Required field is not in input dump!'
      ErrorStatus = 30
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF

  END IF  ! If (pos == 0 ..)

  ! ---------------------------------------------------------------
  ! Some fields, even when found in the input dump, still need
  ! extra work. Mostly for fields that may have differing numbers
  ! of pseudo levels. But also fields needing to be recalculated to
  ! ensure consistency with others.
  ! ---------------------------------------------------------------

  IF ( pos /= 0 .AND. data_source( i ) % source == Input_Dump .AND.        &
      h_int_method /= nearest_neighbour) THEN

    ! ----------------------------------
    ! This item is in the input dump
    ! ----------------------------------

    ! Check to make sure that land and vegetation surface types have been
    ! checked for a change in configuration; will happily run thinking all is
    ! hunky-dory if check is not performed. Catch all for new stashcodes being
    ! added.
    l_rcf_compare_tiles = .FALSE.
    IF ( fields_out( i ) % stashmaster % pt_code == 9 .AND. &
       .NOT. l_aggregate ) l_rcf_compare_tiles = .TRUE.

    SELECT CASE( fields_out( i ) % stashmaster % model )

      ! -----------
      ! Atmos items
      ! -----------

    CASE ( atmos_im )

      SELECT CASE( fields_out( i ) % stashmaster % section )

      CASE ( stashcode_prog_sec )

        SELECT CASE( fields_out( i ) % stashmaster % item )

          ! adv winds exist in some EG dumps but are incorrect
          ! as they are not used by EG
          ! So when going EG --> ND we need to ensure we initialise them

        CASE ( stashcode_u_adv,                                            &
               stashcode_v_adv )

          IF (eg2nd) data_source( i ) % source = Field_Calcs

        CASE ( stashcode_w_adv )

          IF (eg2nd) data_source( i ) % source = Set_To_Zero

          ! some EG prognostics should always be calculated
          ! rather than interpolated. These are only in EG dumps.
        CASE ( stashcode_mv,                                               &
               stashcode_mcl,                                              &
               stashcode_mcf,                                              &
               stashcode_mr,                                               &
               stashcode_mgr,                                              &
               stashcode_mcf2 )

          IF (h_int_active .OR. v_int_active) THEN
            data_source( i ) % source = Field_Calcs
          END IF

        CASE (stashcode_thetavd,                                           &
              stashcode_dry_rho,                                           &
              stashcode_exner_surf )

          IF (h_int_active .OR. v_int_active) THEN
            data_source( i ) % source = Field_Dependent_Calcs
          END IF

        CASE ( stashcode_etadot )

          IF (h_int_active .OR. v_int_active) THEN
            IF (model_type == mt_lam) THEN
              data_source( i ) % source = Set_To_Zero
            ELSE
              data_source( i ) % source = Field_Dependent_Calcs
            END IF
          END IF

        CASE ( stashcode_psiw_lid,                                         &
               stashcode_psiw_surf )

          IF (h_int_active .OR. v_int_active) THEN
            data_source( i ) % source = Set_To_Mdi
          END IF

        CASE ( stashcode_icethick )

          IF (h_int_active) THEN
            data_source( i ) % source = Field_Calcs
          END IF

        CASE ( stashcode_ice_temp_cat,                                     &
               stashcode_ice_conc_cat,                                     &
               stashcode_ice_thick_cat )

          !-------------------------------------------------------------
          ! Check number of pseudo levels in input dump is same as
          ! number of sea ice categories (nice) in output dump.
          ! If not -recalculate
          !-------------------------------------------------------------

          IF (fields_in (pos) % levels /= nice) THEN
            data_source( i ) % source = Field_Calcs
            IF ( mype == 0 .AND.                                           &
                Printstatus >= PrStatus_Normal ) THEN
              WRITE(umMessage,'(A,I7,A)') &
                  "Setting source of stashcode item ",    &
                  fields_out( i ) % stashmaster % item,                    &
                  " to calculations due to difference in NICE"
              CALL umPrint(umMessage,src='rcf_set_data_source_mod')
            END IF
          END IF

        CASE ( stashcode_ice_snow_depth )

          !-------------------------------------------------------------
          ! Check number of pseudo levels in input dump is same as
          ! number of sea ice categories (nice) in output dump.
          ! If not -set to zero
          !-------------------------------------------------------------

          IF (fields_in (pos) % levels /= nice) THEN
            data_source( i ) % source = set_to_zero
            IF ( mype == 0 .AND.                                           &
                Printstatus >= PrStatus_Normal ) THEN
              WRITE(umMessage,'(A,I7,A)')                                  &
                  "Setting source of stashcode item ",                     &
                  fields_out( i ) % stashmaster % item,                    &
                  " to -set-to-zero- due to difference in NICE"
              CALL umPrint(umMessage,src='rcf_set_data_source_mod')
            END IF
          END IF

        CASE ( stashcode_infil_max_tile,                                   &
            stashcode_sw_down_tile )

          !-------------------------------------------------------------
          ! Check that tiles in input dump are the same as
          ! tiles in output dump.
          ! If not -set to zero
          !-------------------------------------------------------------

          CALL rcf_compare_tiles( fields_in, field_count_in, hdr_in,       &
                                  fields_out( i ), hdr_out,                &
                                  fields_out(i) % stashmaster % item,      &
                                  l_match )
          l_rcf_compare_tiles = .FALSE.
          IF ( .NOT. l_match ) THEN
            data_source( i ) % source = set_to_zero
            IF ( mype == 0 .AND.                                           &
                 Printstatus >= PrStatus_Normal ) THEN
              WRITE(umMessage,'(A,I7,A)')                                  &
                  "Setting source of stashcode item ",                     &
                  fields_out( i ) % stashmaster % item,                    &
                  " to -set-to-zero- due to difference in tiles"
              CALL umPrint(umMessage,src='rcf_set_data_source_mod')
            END IF
          END IF

        CASE ( stashcode_catch_tile,                                       &
               stashcode_z0_tile,                                          &
               stashcode_z0h_tile,                                         &
               stashcode_npp_pft_acc,                                      &
               stashcode_g_lf_pft_acc,                                     &
               stashcode_g_ph_lf_pft_acc,                                  &
               stashcode_rsp_w_pft_acc,                                    &
               stashcode_rsp_s_acc )

          !-------------------------------------------------------------
          ! Check that tiles in input dump are the same as
          ! tiles in output dump.
          ! If not -set to user constant -1
          !-------------------------------------------------------------

          CALL rcf_compare_tiles( fields_in, field_count_in, hdr_in,       &
                                  fields_out( i ), hdr_out,                &
                                  fields_out(i) % stashmaster % item,      &
                                  l_match )
          l_rcf_compare_tiles = .FALSE.
          IF ( .NOT. l_match ) THEN
            data_source( i ) % source = set_to_const
            data_source( i ) % rconst = -1.000
            IF ( mype == 0 .AND.                                           &
                 Printstatus >= PrStatus_Normal ) THEN
              WRITE(umMessage,'(A,I7,A)')                                  &
                  "Setting source of stashcode item ",                     &
                  fields_out( i ) % stashmaster % item,                    &
                  " to -1 due to difference in tiles"
              CALL umPrint(umMessage,src='rcf_set_data_source_mod')
            END IF
          END IF

        CASE ( stashcode_can_water_tile,                                   &
               stashcode_tstar_tile,                                       &
               stashcode_tsurf_elev_surft,                                 &
               stashcode_snow_tile,                                        &
               stashcode_catch_snow,                                       &
               stashcode_snow_grnd,                                        &
               stashcode_rgrain,                                           &
               stashcode_frac_surf_type,                                   &
               stashcode_lai,                                              &
               stashcode_canopy_height )

          !-------------------------------------------------------------
          ! Check that tiles in input dump are the same as
          ! tiles in output dump.
          ! If not -set to Field Calcs
          !-------------------------------------------------------------

          CALL rcf_compare_tiles( fields_in, field_count_in, hdr_in,       &
                                  fields_out( i ), hdr_out,                &
                                  fields_out(i) % stashmaster % item,      &
                                  l_match )
          l_rcf_compare_tiles = .FALSE.
          IF ( .NOT. l_match ) THEN
            data_source( i ) % source = Field_Calcs
            IF ( mype == 0 .AND.                                           &
                 Printstatus >= PrStatus_Normal ) THEN
              WRITE(umMessage,'(A,I7,A)')                                  &
                  "Setting source of stashcode item ",                     &
                  fields_out( i ) % stashmaster % item,                    &
                  " to -FieldCalcs- due to difference in tiles"
              CALL umPrint(umMessage,src='rcf_set_data_source_mod')
            END IF
          ELSE IF ( l_snow_tile_gbm_ancil_fix .AND. l_aggregate .AND. &
             fields_out( i ) % stashmaster % item == &
             stashcode_snow_tile ) THEN
            ! Need to keep snow_tile consistent with mean_snow, which result
            ! in inconsistencies when initialising ML snow scheme
            CALL rcf_locate( stashcode_prog_sec, stashcode_mean_snow,          &
               fields_out, field_count_out, pos_out )
            IF ( data_source( pos_out ) % source == Ancillary_File ) THEN
              data_source( i ) % source = Field_Calcs
              IF ( mype == 0 .AND. Printstatus >= PrStatus_Normal ) THEN
                WRITE(umMessage,'(A,I7,A,I7,A)')                            &
                   "Setting source of stashcode item ",                     &
                   fields_out( i ) % stashmaster % item,                    &
                   " to -FieldCalcs- due to source of ",                    &
                   stashcode_mean_snow, " being -AncillaryFile-"
                CALL umPrint(umMessage,src='rcf_set_data_source_mod')
              END IF
            END IF
          END IF

        CASE ( stashcode_snowdep_grd_tile,                                 &
               stashcode_snowpack_bk_dens )

          !-------------------------------------------------------------
          ! Check that tiles in input dump are the same as
          ! tiles in output dump.
          ! This will be set to field_dependent_calcs anyway.
          !-------------------------------------------------------------

          CALL rcf_compare_tiles( fields_in, field_count_in, hdr_in,       &
                                  fields_out( i ), hdr_out,                &
                                  fields_out(i) % stashmaster % item,      &
                                  l_match )
          l_rcf_compare_tiles = .FALSE.
          IF ( .NOT. l_match ) THEN
            data_source( i ) % source = Field_Dependent_Calcs
            IF ( mype == 0 .AND.                                           &
                Printstatus >= PrStatus_Normal ) THEN
              WRITE(umMessage,'(A,I7,A)') &
                  "Setting source of stashcode item ",    &
                  fields_out( i ) % stashmaster % item,                    &
                  " to Field_Dependent_Calcs due to difference in tiles"
              CALL umPrint(umMessage,src='rcf_set_data_source_mod')
            END IF
          END IF
!         If starting a run with multilayer snow from a dump with
!         the zero-layer scheme, we need to tidy fields and remove
!         negative snow or snow depths. The input dump may contain
!         the snow depth, so we want to reset this if initializing,
!         but process it in the standard way otherwise. Use the
!         thickness of snow layers to check whether we have a
!         dump genuinely containing multilayer snow.
          CALL Rcf_Locate( stashcode_prog_sec,                             &
            stashcode_snow_laythk_tiles ,                                  &
            fields_in, field_count_in, pos_in, zero_ok_arg = .TRUE.)
          IF ( (pos_in == 0) .AND. (nsmax > 0) ) THEN
            data_source( i ) % source = Field_Dependent_Calcs
          END IF

        CASE ( stashcode_nsnow_layrs_tiles , &
               stashcode_snow_laythk_tiles , &
               stashcode_snow_ice_tile     , &
               stashcode_snow_liq_tile     , &
               stashcode_snow_t_tile       , &
               stashcode_snow_laydns_tiles , &
               stashcode_snow_grnsiz_tiles )

          !-------------------------------------------------------------
          ! Check that tiles in input dump are the same as
          ! tiles in output dump.
          ! This will be set to field_dependent_calcs anyway.
          !-------------------------------------------------------------

          CALL rcf_compare_tiles( fields_in, field_count_in, hdr_in,       &
                                  fields_out( i ), hdr_out,                &
                                  fields_out(i) % stashmaster % item,      &
                                  l_match )
          l_rcf_compare_tiles = .FALSE.
          IF ( .NOT. l_match ) THEN
            data_source( i ) % source = Field_Dependent_Calcs
            IF ( mype == 0 .AND.                                           &
                Printstatus >= PrStatus_Normal ) THEN
              WRITE(umMessage,'(A,I7,A)') &
                  "Setting source of stashcode item ",    &
                  fields_out( i ) % stashmaster % item,                    &
                  " to Field_Dependent_Calcs due to difference in tiles"
              CALL umPrint(umMessage,src='rcf_set_data_source_mod')
            END IF
          ELSE
            IF (l_force_relayer) THEN
              data_source( i ) % source = Field_Dependent_Calcs
              IF ( mype == 0 .AND.                                         &
                  Printstatus >= PrStatus_Normal ) THEN
                WRITE(umMessage,'(A,I7,A)') &
                    "Setting source of stashcode item ",                   &
                    fields_out( i ) % stashmaster % item,                  &
                    " to Field_Dependent_Calcs because relayering of" //   &
                    " snowpack explicitly requested."
                CALL umPrint(umMessage,src='rcf_set_data_source_mod')
              END IF
            END IF
          END IF

          IF ( data_source( i ) % source == Field_Dependent_Calcs  .AND.   &
             l_rcf_init_flexi ) THEN
            errorstatus = 60
            WRITE(cmessage, '(A)') &
               'Multilayer snow recon (rcf_ml_snowpack) is not ' //        &
               'compatible with l_rcf_init_flexi. Please use '   //        &
               'tile_map_ids instead.'
            CALL ereport( routinename, errorstatus, cmessage )
          END IF

        CASE ( stashcode_dctemp_tile )

          !-------------------------------------------------------------
          ! Check that tiles in input dump are the same as
          ! tiles in output dump.
          ! If not -set to Missing Data Indicator
          !-------------------------------------------------------------

          CALL rcf_compare_tiles( fields_in, field_count_in, hdr_in,       &
                                  fields_out( i ), hdr_out,                &
                                  fields_out(i) % stashmaster % item,      &
                                  l_match )
          l_rcf_compare_tiles = .FALSE.
          IF ( h_int_active .OR. .NOT. l_match ) THEN
            data_source( i ) % source = set_to_mdi
            IF ( mype == 0 .AND.                                           &
                 Printstatus >= PrStatus_Normal ) THEN
              WRITE(umMessage,'(A,I7,A)')                                  &
                 "Setting source of stashcode item ",                      &
                 fields_out( i ) % stashmaster % item,                     &
                 " to -MDI- due to difference in tiles or resolution"
              CALL umPrint(umMessage,src='rcf_set_data_source_mod')
            END IF
          END IF

        CASE (stashcode_dctemp_ssi,                                        &
            stashcode_tm_trans)
          IF (h_int_active .OR.                                            &
              fields_in (pos) % levels /=                                  &
              fields_out( i ) % levels ) THEN
            data_source( i ) % source = set_to_mdi
            IF ( mype == 0 .AND.                                           &
                 Printstatus >= PrStatus_Normal ) THEN
              WRITE(umMessage,'(A,I7,A)')                                  &
                  "Setting source of stashcode item ",                     &
                  fields_out( i ) % stashmaster % item,                    &
                  " to -MDI- due to difference in no. of "//               &
                  " tiles or resolution"
              CALL umPrint(umMessage,src='rcf_set_data_source_mod')
            END IF
          END IF

        CASE( stashcode_urbdisp,        & ! URBAN displacement height
              stashcode_urbztm )          ! URBAN effective roughness length

          ! l_urban2t should be true if these stash are required
          IF ( l_moruses_macdonald ) THEN
            data_source( i ) % source = Field_Calcs
          ELSE
            errorstatus = 60
            WRITE (Cmessage, '(A, 2I5)')                                   &
               'l_moruses_macdonald = F: Data source not set for',         &
               stashcode_urbdisp, stashcode_urbztm
            CALL Ereport( RoutineName, ErrorStatus, Cmessage )
          END IF

        END SELECT                      ! Select on item number
      END SELECT                          ! Select on section number

      !For fields found when reading from GRIB data
      IF ( input_dump_type == grib_input_dump .OR.                         &
           input_dump_type == grib2ff_input_dump ) THEN

        SELECT CASE( fields_out( i ) % stashmaster % section )

        CASE ( stashcode_prog_sec )

          SELECT CASE( fields_out( i ) % stashmaster % item )

          CASE ( stashcode_qrain,     &
                 stashcode_qgraup )

            ! These stash codes have potentially been hijacked to get the
            ! qrain and qsnow moist species from ECMWF GRIB in for the
            ! hybrid height calculation. However there is no known
            ! correspondance between these fields and the UM fields.
            ! Therefore setting both fields to zero explicitly
            ! If anything changes in this regard then obviously this
            ! setting needs to be removed!
            data_source( i ) % source = Set_To_Zero

          END SELECT                      ! Select on item number
        END SELECT                          ! Select on section number
      END IF

    CASE DEFAULT                 ! Couldn't find a proper Internal Model

      ! The Model type didn't match one of the three defined types.
      WRITE (Cmessage, '(2A, I2, A, I3, A, I5)')                           &
          " Couldn't find a Sub-Model ID type for : ",                     &
          " Model ", fields_out( i ) % stashmaster % model,                &
          " Section ", fields_out( i ) % stashmaster % section,            &
          " Item ",  fields_out( i ) % stashmaster % item
      ErrorStatus = 70
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )

    END SELECT                            ! Select on Internal model

    IF ( l_rcf_compare_tiles ) THEN
      WRITE(umMessage,'(A, I5)')                                          &
         "Please make provision for the unreferenced tiled prognostic: ", &
         fields_out( i ) % stashmaster % item
      CALL umPrint(umMessage,src='rcf_set_data_source_mod')
      WRITE (Cmessage,'(2A, I2, A, I3, A, I5)')                           &
         " Tiled prognostic was not checked by rcf_compare_tiles:",       &
         " Model ", fields_out( i ) % stashmaster % model,                &
         " Section ", fields_out( i ) % stashmaster % section,            &
         " Item ",  fields_out( i ) % stashmaster % item
      ErrorStatus = 71
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF

  END IF                                  ! Item was in Dump
END DO ! Loop over output fields

! ---------------------------------------------------------------
! Some fields still need extra work regardless of whether they're
! in the input dump or not.  Need to restart loop since we need to make sure
! all fields are set with some default values above.
! ---------------------------------------------------------------
DO i = 1, field_count_out
  SELECT CASE( fields_out( i ) % stashmaster % model )

    ! -----------
    ! Atmos items
    ! -----------

  CASE ( atmos_im )

    SELECT CASE( fields_out( i ) % stashmaster % section )

    CASE ( stashcode_prog_sec )

      SELECT CASE( fields_out( i ) % stashmaster % item )

        ! VN8.5 restart dumps contain a dry_rho set to RMDI. Here we pick it up
        ! and instruct the reconfiguration to populate it

      CASE ( stashcode_dry_rho )

        CALL Rcf_Locate( fields_out( i ) % stashmaster % section,          &
                         fields_out( i ) % stashmaster % item,             &
                         fields_in, field_count_in, pos, zero_ok_arg = .TRUE.)
        IF (pos/=0) THEN
          dryrho_in => fields_in(pos)
          CALL rcf_alloc_field( dryrho_in )
          CALL rcf_read_field(  dryrho_in, hdr_in, decomp_rcf_input )

          rmdi_tol  = ABS (rmdi) * 0.0001

          IF ( dryrho_in % DATA (1,1) > rmdi-rmdi_tol .AND.                &
               dryrho_in % DATA (1,1) < rmdi+rmdi_tol ) THEN  
            data_source( i ) % source = Field_Dependent_Calcs
            WRITE(umMessage,'(A,I7,A,A)')                                  &
                    "Setting source of stashcode item ",                   &
                    fields_out( i ) % stashmaster % item,                  &
                    " to -Field_Dependent_Calcs- as RMDI detected.",       &
                    " Should only occur  with VN8.5 input dumps"
            CALL umPrint(umMessage,src='rcf_set_data_source_mod')
          END IF

          CALL rcf_dealloc_field( dryrho_in )
        END IF

        ! Field in output dump is one of the dust bins.
      CASE ( stashcode_dust1_mmr,                                          &
             stashcode_dust2_mmr,                                          &
             stashcode_dust3_mmr,                                          &
             stashcode_dust4_mmr,                                          &
             stashcode_dust5_mmr,                                          &
             stashcode_dust6_mmr)

        !-------------------------------------------------------------
        ! Check for presence of a dust bin in both input dump
        ! and output dump. pos and pos_out will both be non-zero.
        ! This means both dumps are using a dust scheme.
        !-------------------------------------------------------------
        CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust1_mmr,          &
                         fields_in, field_count_in, pos,zero_ok_arg =  .TRUE.)
        CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust1_mmr,          &
                         fields_out, field_count_out, pos_out,              &
                          zero_ok_arg = .TRUE.)
        IF ( pos /= 0 .AND. pos_out /= 0 ) THEN
          !-------------------------------------------------------------
          ! Check for presence of 3rd dust bin in both input and output
          ! dumps.
          ! If dust3 is in input dump (pos != 0) but not in output (pos_out=0)
          ! OR dust3 not in input dump (pos=0) but is in output (pos_out != 0)
          ! Then changing scheme means dust need fieldcalcs
          !-------------------------------------------------------------

          CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust3_mmr,        &
                           fields_in, field_count_in, pos,                  &
                           zero_ok_arg = .TRUE.)
          CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust3_mmr,        &
                           fields_out, field_count_out, pos_out,           &
                           zero_ok_arg = .TRUE.)

          IF ( ( pos /= 0 .AND. pos_out == 0 ) .OR.                        &
              ( pos == 0 .AND. pos_out /= 0 ) ) THEN
            data_source( i ) % source = Field_Calcs
            IF ( mype == 0 .AND.                                           &
                 Printstatus >= PrStatus_Normal ) THEN
              WRITE(umMessage,'(A,I7,A)') &
                  "Setting source of stashcode item ",    &
                  fields_out( i ) % stashmaster % item,                    &
                  " to -FieldCalcs- due to difference in no. of dust bins"
              CALL umPrint(umMessage,src='rcf_set_data_source_mod')
            END IF
          END IF ! If dust bin 3 only in one of input or output dumps.

        END IF ! If dust bin 1 present in both input and output dumps.
      CASE ( stashcode_zw,                                                 &
             stashcode_fexp,                                               &
             stashcode_sthzw,                                              &
             stashcode_fsat,                                               &
             stashcode_fwetl,                                              &
             stashcode_gamtot,                                             &
             stashcode_a_fsat,                                             &
             stashcode_c_fsat,                                             &
             stashcode_a_fwet,                                             &
             stashcode_c_fwet )
        ! NB. Some of these may have had data_source set above.
        ! For LSH we can use subroutine to check consistency of data sources.
        CALL rcf_lsh_field_checks ( data_source, fields_out,               &
                                    field_count_out, i )

        ! LSM is set independently of everything else since we
        ! need LSM for land-packed fields.  This was done before this routine
        ! was called.
      CASE (stashcode_lsm)
        ! LSM are performed already.
        IF (data_source( i ) % source == input_dump) THEN
          data_source( i ) % source = already_processed
        END IF

      END SELECT                      ! Select on item number

    END SELECT                          ! Select on section number

  CASE DEFAULT                 ! Couldn't find a proper Internal Model

    ! The Model type didn't match one of the three defined types.
    WRITE (Cmessage, '(2A, I2, A, I3, A, I5)')                             &
        " Couldn't find a Sub-Model ID type for : ",                       &
        " Model ", fields_out( i ) % stashmaster % model,                  &
        " Section ", fields_out( i ) % stashmaster % section,              &
        " Item ",  fields_out( i ) % stashmaster % item
    ErrorStatus = 80
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )

  END SELECT                            ! Select on Internal model

  ! Check on Rimwidtha if an atmos LBC file for copying
  IF ( data_source( i ) % source == Input_Dump .AND.                       &
      (fields_out( i ) % stashmaster % grid_type ==ppx_atm_lbc_theta .OR.  &
       fields_out( i ) % stashmaster % grid_type == ppx_atm_lbc_u   .OR.   &
       fields_out( i ) % stashmaster % grid_type == ppx_atm_lbc_v) ) THEN

    IF ( hdr_in  % Lookup( lbrow, fields_out( i ) % dump_pos ) /=          &
         hdr_out % Lookup( lbrow, fields_out( i ) % dump_pos) ) THEN
      Cmessage = 'Rimwidth needs to be equal for input ' //                &
          'and output grids if lbcs are copied'
      ErrorStatus = 90
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF
  END IF

END DO ! Loop over output fields

!------------------------------------------------------------------
! Print out some diagnostics if required....
!------------------------------------------------------------------
IF (mype == 0 .AND. PrintStatus >= PrStatus_Oper ) THEN
  ! form   ="(2(i5,' '),4(i4,' '),e10.4,' ',50a)"
  ! form_c ="(2(a5,' '),4(a4,' '),a10,  ' ',50a)"
  CALL umPrint( '',src='rcf_set_data_source_mod')
  CALL umPrint( 'Data source (after ancil processing):',                   &
                 src='rcf_set_data_source_mod')
  WRITE(umMessage,form_c) 'Sect', 'Item', 'src', 'Dom', 'AncS', 'AncI',    &
                          'Real Const', 'Ancil or NetCDF File'
  CALL umPrint(umMessage,src='rcf_set_data_source_mod')

  DO i = 1, field_count_out
    WRITE(umMessage,FORM) fields_out( i ) % stashmaster % section,         &
          fields_out(  i ) % stashmaster % item,                           &
          data_source( i ) % source,                                       &
          data_source( i ) % Domain,                                       &
          data_source( i ) % Ancil_SctnC,                                  &
          data_source( i ) % Ancil_ItemC,                                  &
          data_source( i ) % RConst,                                       &
     TRIM(data_source( i ) % Ancil_File)
    CALL umPrint(umMessage,src='rcf_set_data_source_mod')
  END DO
  CALL umPrint( '',src='rcf_set_data_source_mod')

END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Reset_Data_Source
END MODULE Rcf_Set_Data_Source_Mod
