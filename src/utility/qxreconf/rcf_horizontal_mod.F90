! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Performs horizontal interpolation and related tasks

MODULE Rcf_horizontal_Mod
 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

!  Subroutine Rcf_horizontal - horizontal interpolation
!
! Description:
!   This module contains a wrapper subroutine for horizontal
!   interpolation. The interpolation is done level-by-level
!   but on a vertically decomposed grid. Thus we need to gather
!   data by levels onto the compute PEs and then rescatter it onto
!   the output grids. This won't be the fastest routine in the world,
!   but should be reliable and be easily verifiable for
!   bit-reproducibility etc.
!
! Method:
!   Note that land compressed fields are uncompressed and recompressed,
!   copying is done if that is all that is required, coastal
!   adjustment is performed and polar row averaging is done.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_HORIZONTAL_MOD'

CONTAINS

SUBROUTINE Rcf_horizontal( field_in, field_out, grid_in, grid_out )

USE Ereport_mod, ONLY: &
    Ereport

USE mask_compression, ONLY: compress_to_mask, expand_from_mask

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE UM_ParVars, ONLY:  &
    gc_all_proc_group,  &
    change_decomposition

USE UM_ParCore, ONLY: &
    nproc,            &
    mype

USE UM_ParParams, ONLY: &
    halo_type_single

USE decomp_params, ONLY: &
    decomp_rcf_input,     &
    decomp_rcf_output

USE Rcf_select_weights_mod, ONLY: &
    Rcf_select_weights

USE Rcf_average_polar_mod, ONLY: &
    Rcf_average_polar

USE Rcf_H_Int_Ctl_Mod, ONLY: &
    Rcf_H_Int_Ctl

USE umPrintMgr, ONLY:       &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Diag,          &
    newline

USE nlstcall_mod, ONLY:   &
    LTimer

USE Rcf_Lsm_Mod, ONLY: &
    n_coastal_points,                n_land_points_unres, &
    n_sea_points_unres,              coast_index_in,      &
    coast_index_out,                 index_targ_land,     &
    index_targ_sea,                  land_unres_index,    &
    sea_unres_index,                 local_lsm_in,        &
    local_lsm_out,                   glob_lsm_out,        &
    cyclic,                          l_lsm_out_present,   &
    land_unres_constrain_index,                           &
    land_unres_notconstrain_index,                        &
    sea_unres_notconstrain_index

USE rcf_global_to_local_mod, ONLY: rcf_get_fld_type

USE rcf_nlist_recon_science_mod, ONLY: &
    coast_adj_circle_method,           &
    coast_adj_method,                  &
    coast_adj_spiral_method,           &
    coast_adj_standard_method

USE jules_snow_mod, ONLY: nsmax

USE Rcf_Set_Interp_Flags_Mod, ONLY: &
    interp_h_only,                   interp_v_only,       &
    interp_all,                      interp_copy,         &
    interp_no_op

USE um_stashcode_mod, ONLY:      &
    stashcode_prog_sec,          &
    stashcode_tstar,             &
    stashcode_z0,                &
    stashcode_vol_smc_wilt,      &
    stashcode_vol_smc_cri,       &
    stashcode_vol_smc_sat,       &
    stashcode_mean_snow,         &
    stashcode_snow_tile,         &
    stashcode_snow_grnd,         &
    stashcode_snowdep_grd_tile,  &
    stashcode_snowpack_bk_dens,  &
    stashcode_nsnow_layrs_tiles, &
    stashcode_snow_laythk_tiles, &
    stashcode_snow_ice_tile,     &
    stashcode_snow_liq_tile,     &
    stashcode_snow_T_tile,       &
    stashcode_snow_laydns_tiles, &
    stashcode_snow_grnsiz_tiles

USE Rcf_Scatter_Zonal_Field_Mod, ONLY: &
    Rcf_Scatter_Zonal_Field

USE Rcf_Gather_Zonal_Field_Mod, ONLY: &
    Rcf_Gather_Zonal_Field

USE Rcf_Interp_Weights_Mod  ! All of it

USE Rcf_HeadAddress_Mod, ONLY: &
    FH_GridStagger

USE cppxref_mod, ONLY:                     &
    ppx_type_real,      ppx_type_int,       &
    ppx_type_log,       ppx_atm_compressed, &
    ppx_atm_ozone,      ppx_atm_tall,       &
    ppx_atm_tzonal,     ppx_atm_tsea

USE science_fixes_mod, ONLY: l_roughnesslength_fix

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Derived Type Arguments
TYPE (field_type), TARGET, INTENT(INOUT) :: field_in  ! Input data field
TYPE (field_type), TARGET, INTENT(INOUT) :: field_out !Output data field
TYPE (grid_type), INTENT(IN)     :: grid_in   ! Input grid sizes
TYPE (grid_type), INTENT(IN)     :: grid_out  ! Output grid sizes

! Local data
INTEGER             :: div         ! for calculation of partition
INTEGER             :: rem         ! for calculation of partition
INTEGER             :: pe          ! PE from which to gather/scatter
INTEGER             :: i,j,k       ! Looping
INTEGER             :: fld_type    ! P, U or V field?
INTEGER             :: mask_points ! check value for lsm
INTEGER             :: stat        ! GCOM status
INTEGER             :: field_averaged ! 1/0 status of field
INTEGER             :: orig_h_int_method ! Stores original interp method
                                         ! since soil can use different one.
INTEGER             :: sea_points_unres_tmp  ! tmp stuff for coast aj
INTEGER             :: land_points_unres_tmp ! tmp stuff for coast aj
INTEGER             :: lsm_tmp( grid_out % glob_p_rows * &
                                grid_out % glob_p_row_length )
INTEGER             :: Index_targ_land_tmp( grid_out % glob_p_rows * &
                                          grid_out % glob_p_row_length )
INTEGER             :: Index_targ_sea_tmp( grid_out % glob_p_rows * &
                                          grid_out % glob_p_row_length )
INTEGER             :: maxdim      ! a size parameter for coast aj
INTEGER             :: ErrorStatus
LOGICAL             :: averaged    ! return from polar_average
LOGICAL             :: l_print_coastal_warning ! Flag to print warning if unable
                                               ! to perform coastal adjustment
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_HORIZONTAL'
CHARACTER (LEN=errormessagelength)  :: Cmessage

                                       ! Single level after scatter
REAL, ALLOCATABLE   :: level_field_in( : )

                                        ! Single level after interp.
REAL, ALLOCATABLE   :: level_field_out( : )

TYPE (field_type), TARGET     :: field_in_tmp
TYPE (field_type), TARGET     :: field_out_tmp
TYPE (field_type), POINTER    :: ptr_field_in
TYPE (field_type), POINTER    :: ptr_field_out

! Pointers for choice of data etc
INTEGER, POINTER              :: ptr_bl_index_b_l(:)
INTEGER, POINTER              :: ptr_bl_index_b_r(:)
INTEGER, POINTER              :: ptr_bl_index_t_l(:)
INTEGER, POINTER              :: ptr_bl_index_t_r(:)
INTEGER, POINTER              :: ptr_aw_index_targ_lhs(:)
INTEGER, POINTER              :: ptr_aw_index_targ_top(:)
REAL, POINTER                 :: ptr_aw_colat_t(:)
REAL, POINTER                 :: ptr_aw_long_l(:)
REAL, POINTER                 :: ptr_weight_b_l(:)
REAL, POINTER                 :: ptr_weight_b_r(:)
REAL, POINTER                 :: ptr_weight_t_l(:)
REAL, POINTER                 :: ptr_weight_t_r(:)
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


averaged = .FALSE.

IF (LTimer) CALL Timer( RoutineName, 3)

! If input and output datatypes are different, we have a problem.
IF (  field_in % stashmaster % data_type /=                     &
     field_out % stashmaster % data_type ) THEN
  Cmessage = 'Input and Output field datatypes differ!'
  ErrorStatus = 10
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! Make sure that field_averaged and coastal warning flag are initialised
field_averaged = 0
l_print_coastal_warning = .FALSE.

!-------------------------------------------------------------------
! Is interpolation activated? If not, copy data across and that's
! all we will do.
!-------------------------------------------------------------------
SELECT CASE( field_in % interp )
CASE ( interp_copy, interp_v_only )     ! copy data

  ! Sizes should be the same, but will check...
  IF ( field_in % level_size /= field_out % level_size ) THEN
    WRITE(umMessage,*) "Aborting due to mismatch in local datasizes "
    CALL umPrint(umMessage,src='rcf_horizontal_mod')
    WRITE(umMessage,*) "Input dump data size =",field_in % level_size
    CALL umPrint(umMessage,src='rcf_horizontal_mod')
    WRITE(umMessage,*) "Output dump data size =",field_out % level_size
    CALL umPrint(umMessage,src='rcf_horizontal_mod')
    Cmessage = 'No interpolation required, but input and output '//&
        'data fields are different sizes'
    ErrorStatus = 20
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  SELECT CASE( field_in % stashmaster % data_type )
  CASE ( ppx_type_real )
    field_out % DATA(:,:) = field_in % DATA(:,:)

  CASE ( ppx_type_int )
    IF ( ALLOCATED( field_in %  data ) ) THEN
      field_out % data(:,:) = field_in % data(:,:)
    ELSE
      field_out % Data_Int(:,:) = field_in % Data_Int(:,:)
    END IF

  CASE ( ppx_type_log )
    IF ( ALLOCATED( field_in % data ) ) THEN
      field_out % data(:,:) = field_in % data(:,:)
    ELSE
      field_out % Data_Log(:,:) = field_in % Data_Log(:,:)
    END IF

  CASE DEFAULT
    Cmessage = 'Unsupported Data-Type'
    ErrorStatus = -30
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )

  END SELECT

CASE ( interp_all, interp_h_only )

  !---------------------------------------------------------------
  ! If data is compressed onto land points, need to expand it
  ! To do this neatly requires some fiddling with pointers.
  ! ptr_field_in and ptr_field_out are the ones to work with. Set
  ! up here.
  !--------------------------------------------------------------
  IF ( field_in % stashmaster % grid_type == ppx_atm_compressed) THEN
    ALLOCATE( field_in_tmp % DATA( grid_in % loc_p_rows *           &
                             grid_in % loc_p_row_length ,           &
                             field_in % levels ) )

    ALLOCATE( field_out_tmp % DATA( grid_out % loc_p_rows *         &
                              grid_out % loc_p_row_length,          &
                              field_out % levels ) )

    field_in_tmp % levels       = field_in % levels
    field_in_tmp % rows         = grid_in % loc_p_rows
    field_in_tmp % row_len      = grid_in % loc_p_row_length
    field_in_tmp % level_size   = field_in_tmp % rows *              &
                                  field_in_tmp % row_len
    field_in_tmp % glob_rows    = grid_in % glob_p_rows
    field_in_tmp % glob_row_len = grid_in % glob_p_row_length
    field_in_tmp % glob_level_size = field_in_tmp % glob_rows *      &
                                     field_in_tmp % glob_row_len
    field_in_tmp % stashmaster => field_in % stashmaster

    field_out_tmp % levels       = field_out % levels
    field_out_tmp % rows         = grid_out % loc_p_rows
    field_out_tmp % row_len      = grid_out % loc_p_row_length
    field_out_tmp % level_size   = field_out_tmp % rows *              &
                                   field_out_tmp % row_len
    field_out_tmp % glob_rows    = grid_out % glob_p_rows
    field_out_tmp % glob_row_len = grid_out % glob_p_row_length
    field_out_tmp % glob_level_size = field_out_tmp % glob_rows *      &
                                      field_out_tmp % glob_row_len
    field_out_tmp % stashmaster => field_out % stashmaster

    ! Expand level by level
    DO i = 1, field_in % levels
      CALL expand_from_mask( field_in_tmp % DATA(1:,i),                   &
                             field_in % DATA(1:,i),                    &
                             local_lsm_in,  field_in_tmp % level_size,&
                             field_in % level_size )
    END DO

    ptr_field_in  => field_in_tmp
    ptr_field_out => field_out_tmp
  ELSE IF ( field_in % stashmaster % grid_type == ppx_atm_ozone .AND. &
           ( field_in % glob_row_len  /= 1                      .OR. &
             field_out % glob_row_len /= 1) ) THEN
    ! We still want to make zonal to zonal work.  This is to fix anything to
    ! do with full field ozone.
    ALLOCATE( field_in_tmp % DATA( grid_in % loc_p_rows *           &
                             grid_in % loc_p_row_length ,           &
                             field_in % levels ) )

    ALLOCATE( field_out_tmp % DATA( grid_out % loc_p_rows *         &
                              grid_out % loc_p_row_length,          &
                              field_out % levels ) )

    field_in_tmp % levels       = field_in % levels
    field_in_tmp % rows         = grid_in % loc_p_rows
    field_in_tmp % row_len      = grid_in % loc_p_row_length
    field_in_tmp % level_size   = field_in_tmp % rows *              &
                                  field_in_tmp % row_len
    field_in_tmp % glob_rows    = grid_in % glob_p_rows
    field_in_tmp % glob_row_len = grid_in % glob_p_row_length
    field_in_tmp % glob_level_size = field_in_tmp % glob_rows *      &
                                     field_in_tmp % glob_row_len
    ! Create new stashmaster entry to change from ozone to full field.
    ALLOCATE(field_in_tmp % stashmaster)
    field_in_tmp % stashmaster   = field_in % stashmaster
    field_in_tmp % stashmaster % grid_type = ppx_atm_tall

    field_out_tmp % levels       = field_out % levels
    field_out_tmp % rows         = grid_out % loc_p_rows
    field_out_tmp % row_len      = grid_out % loc_p_row_length
    field_out_tmp % level_size   = field_out_tmp % rows *              &
                                   field_out_tmp % row_len
    field_out_tmp % glob_rows    = grid_out % glob_p_rows
    field_out_tmp % glob_row_len = grid_out % glob_p_row_length
    field_out_tmp % glob_level_size = field_out_tmp % glob_rows *      &
                                      field_out_tmp % glob_row_len
    ! Create new stashmaster entry to change from ozone to full field.
    ALLOCATE(field_out_tmp % stashmaster)
    field_out_tmp % stashmaster   = field_out % stashmaster
    field_out_tmp % stashmaster % grid_type = ppx_atm_tall

    ! Expand level by level
    IF (field_in % glob_row_len == 1) THEN
      ! Zonal ozone
      DO k = 1, field_in_tmp % levels
        DO j = 1, field_in_tmp % rows
          DO i = 1, field_in_tmp % row_len
            field_in_tmp % DATA((j-1)*field_in_tmp%row_len+i,k) =          &
                                                        field_in % DATA(j,k)
          END DO
        END DO
      END DO
    ELSE
      ! Full field ozone
      field_in_tmp % DATA(:,:) = field_in % DATA(:,:)
    END IF

    ptr_field_in  => field_in_tmp
    ptr_field_out => field_out_tmp

  ELSE
    ptr_field_in  => field_in
    ptr_field_out => field_out
  END IF

  !------------------------------------------------------------------
  ! Can now allocate space for local levels for interpolation
  !------------------------------------------------------------------
  ALLOCATE( level_field_in( ptr_field_in % glob_level_size ) )
  ALLOCATE( level_field_out( ptr_field_out % glob_level_size ) )

  !-------------------------------------------------------------------
  ! Get the field type for comms
  !-------------------------------------------------------------------
  fld_type = rcf_get_fld_type(field_out % stashmaster % grid_type)
  !-------------------------------------------------------------------
  ! Need to work out which weights we wish to use. This depends on
  ! the type of field that is being interpolated.
  !-------------------------------------------------------------------

  ! For various land surface fields, nearest-neighbour interpolation
  ! is required, despite what may be selected more generally.
  ! Enforce this condition, but allow the original option to be
  ! restored.

  orig_h_int_method = h_int_method

! Soil Moisture
  IF ( ptr_field_in % stashmaster % section  == stashcode_prog_sec      .AND.  &
       ( ptr_field_in % stashmaster % item   == stashcode_vol_smc_wilt  .OR.   &
         ptr_field_in % stashmaster % item   == stashcode_vol_smc_cri   .OR.   &
         ptr_field_in % stashmaster % item   == stashcode_vol_smc_sat ) .AND.  &
       smcp_int_nearest_neighbour ) THEN
    h_int_method = nearest_neighbour
  END IF
!
! Multilayer Snow: Multilayer fields cannot be interpolated continuously 
! and must be done using nearest neighbour methods. Snow amounts and depths 
! could in theory be done continously, but it is much better to interpolate 
! consistently.
  IF ( ptr_field_in % stashmaster % section ==  &
             stashcode_prog_sec          .AND.  &
       ( ptr_field_in % stashmaster % item  ==  &
             stashcode_nsnow_layrs_tiles .OR.   &
         ptr_field_in % stashmaster % item  ==  &
             stashcode_snow_laythk_tiles .OR.   &
         ptr_field_in % stashmaster % item  ==  &
             stashcode_snow_ice_tile     .OR.   &
         ptr_field_in % stashmaster % item  ==  &
             stashcode_snow_liq_tile     .OR.   &
         ptr_field_in % stashmaster % item  ==  &
             stashcode_snow_T_tile       .OR.   &
         ptr_field_in % stashmaster % item  ==  &
             stashcode_snow_laydns_tiles .OR.   &
         ptr_field_in % stashmaster % item  ==  &
             stashcode_snow_grnsiz_tiles )      &
       ) THEN
    h_int_method = nearest_neighbour
  END IF
  IF (nsmax > 0) THEN
    IF ( ptr_field_in % stashmaster % section == &
               stashcode_prog_sec         .AND.  &
         ( ptr_field_in % stashmaster % item  == &
               stashcode_mean_snow        .OR.   &
           ptr_field_in % stashmaster % item  == &
               stashcode_snow_tile        .OR.   &
           ptr_field_in % stashmaster % item  == &
               stashcode_snow_grnd        .OR.   &
           ptr_field_in % stashmaster % item  == &
               stashcode_snowdep_grd_tile .OR.   &
           ptr_field_in % stashmaster % item  == &
               stashcode_snowpack_bk_dens )      &
         ) THEN
      h_int_method = nearest_neighbour
    END IF
  END IF

  CALL Rcf_select_weights( ptr_bl_index_b_l, ptr_bl_index_b_r, &
                           ptr_bl_index_t_l, ptr_bl_index_t_r, &
                           ptr_weight_b_l, ptr_weight_b_r,     &
                           ptr_weight_t_l, ptr_weight_t_r,     &
                           ptr_aw_index_targ_lhs,              &
                           ptr_aw_index_targ_top,              &
                           ptr_aw_colat_t, ptr_aw_long_l,      &
                           ptr_field_in % stashmaster % grid_type, &
                           ptr_field_in % stashmaster % section,   &
                           ptr_field_in % stashmaster % item )

  !--------------------------------------------------------------------
  ! We need to gather <nproc> levels onto pes.
  !--------------------------------------------------------------------

  div = ptr_field_in % levels / nproc
  rem = MOD( ptr_field_in % levels, nproc )
  pe = 0


  DO i = 1, div

    CALL Change_Decomposition( decomp_rcf_input )

    DO j = ((i-1) * nproc) + 1, i * nproc
      ! Will gather level j on pe

      IF (ptr_field_in % glob_row_len == 1) THEN ! Zonal Data
        CALL Rcf_Gather_Zonal_Field( ptr_field_in % DATA(:,j),         &
                                     level_field_in,                   &
                                     ptr_field_in % level_size,        &
                                     ptr_field_in % glob_level_size, 1,&
                                     ppx_atm_tzonal, pe )

      ELSE

!DEPENDS ON: gather_field
        CALL Gather_Field( ptr_field_in % DATA(:,j),              &
                           level_field_in,                        &
                           ptr_field_in % row_len,                &
                           ptr_field_in % rows,                   &
                           ptr_field_in % glob_row_len,           &
                           ptr_field_in % glob_rows,              &
                           fld_type,                              &
                           halo_type_single,                      &
                           pe,                                    &
                           gc_all_proc_group )
      END IF

      pe = pe + 1
      IF (pe == nproc) pe = 0
    END DO

    !-------------------------------------------------------------------
    ! All PEs are currently full of level data - so do the interpolation
    !-------------------------------------------------------------------
    CALL Rcf_H_Int_Ctl( ptr_field_out % glob_level_size,                 &
                        ptr_field_in % glob_row_len,                     &
                        ptr_field_out % glob_row_len,                    &
                        ptr_field_in % glob_rows,                        &
                        ptr_field_out % glob_rows,                       &
                        grid_out % global, ptr_aw_index_targ_lhs,        &
                        ptr_aw_index_targ_top, ptr_bl_index_b_l,         &
                        ptr_bl_index_b_r, ptr_bl_index_t_l,              &
                        ptr_bl_index_t_r, ptr_aw_colat_t, ptr_aw_long_l, &
                        level_field_in, ptr_weight_t_r, ptr_weight_b_r,  &
                        ptr_weight_t_l, ptr_weight_b_l, level_field_out )

    !-------------------------------------------------------------------
    ! Coastal Adjustment for land-only, sea-only and T*  fields
    !-------------------------------------------------------------------
    IF (LTimer) CALL Timer( 'Coastal_Adjustment', 103)

    IF (ptr_field_out % stashmaster % grid_type == ppx_atm_compressed .OR.   &
         ptr_field_out % stashmaster % grid_type == ppx_atm_tsea .OR.        &
         ( ptr_field_out % stashmaster % section == stashcode_prog_sec .AND. &
         (ptr_field_out % stashmaster % item    == stashcode_tstar .OR.      &
         (ptr_field_out % stashmaster % item    == stashcode_z0 .AND.        &
          l_roughnesslength_fix )))) THEN
      ! Only able to perform coastal adjustment if land sea mask is present 
      ! in output dump
      IF ( l_lsm_out_present ) THEN 
        DO j = 1, n_coastal_points
          level_field_out( coast_index_out( j ) ) = &
                           level_field_in( coast_index_in( j ) )
        END DO

        IF (coast_adj_method == coast_adj_spiral_method) THEN 
                       ! Spiral adjustment
          IF (PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
            WRITE(umMessage,'(A)') "Spiral method is used"
            CALL umPrint(umMessage,src='replanca-rcf_replanca')
          END IF

          maxdim = MIN(grid_out % glob_p_rows, grid_out % glob_p_row_length)
          DO j = 1, grid_out % glob_p_rows * grid_out % glob_p_row_length
            index_targ_sea_tmp( j )  = index_targ_sea( j )
            index_targ_land_tmp( j ) = index_targ_land( j )
            IF ( glob_lsm_out( j ) ) THEN
              lsm_tmp( j ) = 1
            ELSE
              lsm_tmp( j ) = 0
            END IF
          END DO

          sea_points_unres_tmp  = n_sea_points_unres
          land_points_unres_tmp = n_land_points_unres

          IF (ptr_field_out % stashmaster % grid_type /=                   &
                                            ppx_atm_compressed) THEN
            ! Only do the sea coastal points if the field isn't land only.
            ! DEPENDS ON: intf_coast_aj
            CALL Intf_Coast_AJ( lsm_tmp, index_targ_sea_tmp,               &
                 sea_points_unres_tmp, grid_out % glob_p_rows,             &
                 grid_out % glob_p_row_length, level_field_out, 0, cyclic, &
                 maxdim )
          END IF

          ! DEPENDS ON: intf_coast_aj
          CALL Intf_Coast_AJ( lsm_tmp, index_targ_land_tmp,              &
               land_points_unres_tmp, grid_out % glob_p_rows,            &
               grid_out % glob_p_row_length, level_field_out, 1, cyclic, &
               maxdim )

        ELSE IF (coast_adj_method == coast_adj_standard_method) THEN
                       ! Non-spiral adjustment
          IF (PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
            WRITE(umMessage,'(A)') "Standard method is used"
            CALL umPrint(umMessage,src='replanca-rcf_replanca')
          END IF

          DO j = 1, n_land_points_unres
            level_field_out( index_targ_land( j ) ) = &
                             level_field_out( land_unres_index( j ) )
          END DO

          DO j = 1, n_sea_points_unres
            level_field_out( index_targ_sea( j ) ) = &
                             level_field_out( sea_unres_index( j ) )
          END DO

        ELSE IF (coast_adj_method == coast_adj_circle_method) THEN
                       ! Spiral circle method
          IF (PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
            WRITE(umMessage,'(A)') "Spiral circle method is used"
            CALL umPrint(umMessage,src='replanca-rcf_replanca')
          END IF

          ! SPIRAL CIRCLE METHOD

          IF (ptr_field_out % stashmaster % grid_type ==                   &
                                            ppx_atm_compressed) THEN
            ! If land only do not use the 200km constraint
            DO j = 1, n_land_points_unres
              level_field_out( index_targ_land( j ) ) = &
                      level_field_out( land_unres_notconstrain_index( j ) )
            END DO
          ELSE
            ! Not land only
            DO j = 1, n_land_points_unres
              level_field_out( index_targ_land( j ) ) = &
                      level_field_out( land_unres_constrain_index( j ) )
            END DO

            ! Do sea points
            DO j = 1, n_sea_points_unres
              level_field_out( index_targ_sea( j ) ) = &
                    level_field_out( sea_unres_notconstrain_index( j ) )
            END DO
          END IF
        ELSE ! Not a valid coastal_adjustment method
          Cmessage = 'Not a valid coastal adjustment method'
          ErrorStatus = 160
          CALL Ereport( RoutineName, ErrorStatus, Cmessage )
        END IF ! Spiral
      ELSE
        ! Unable to perform coastal adjustment due to lack of land sea mask.
        ! Print warning later in routine, only want to print once per field.
        l_print_coastal_warning = .TRUE.
      END IF
    END IF

    IF (LTimer) CALL Timer( 'Coastal_Adjustment', 104)

    !-----------------------------------------------------------------
    ! Average the polar rows
    !-----------------------------------------------------------------
    IF ( ptr_field_out % stashmaster % grid_type <= 3 .OR.               &
         ptr_field_out % stashmaster % grid_type == 21 ) THEN
      IF ( ptr_field_out % glob_rows > grid_out % glob_v_rows ) THEN
        CALL Rcf_Average_polar(level_field_out, ptr_field_out % glob_rows, &
                             ptr_field_out % glob_row_len,                 &
                             grid_out % global, averaged )
      END IF

      IF ( averaged .AND. PrintStatus >= PrStatus_Normal) THEN
        field_averaged = 1
      END IF
    END IF
    IF ( field_out % stashmaster % grid_type == ppx_atm_ozone ) THEN
      IF (field_out % glob_row_len == 1 .AND. field_in % glob_row_len /= 1) THEN
        ! Need to use field_out since we have change stashmaster and row_len to
        ! full field.
        DO j = 1, field_out % rows
          ! We can reuse polar rows and average over 1 row
          CALL Rcf_Average_polar(                                        &
                      level_field_out((j-1)*ptr_field_out % glob_row_len:&
                                      j*ptr_field_out % glob_row_len),   &
                                  1,                                     &
                                  ptr_field_out % glob_row_len,          &
                                  grid_out % global, averaged )

        END DO
      END IF

    END IF

    !-------------------------------------------------------------------
    ! And re-scatter the data back to original PEs
    !------------------------------------------------------------------
    CALL Change_Decomposition( decomp_rcf_output )
    DO j = ((i-1) * nproc) + 1, i * nproc

      IF (ptr_field_out % glob_row_len == 1) THEN ! Zonal Data
        CALL Rcf_Scatter_Zonal_Field( ptr_field_out % DATA(:,j),        &
                                  level_field_out,                      &
                                  ptr_field_out % level_size,           &
                                  ptr_field_out % glob_level_size, 1,   &
                                  ppx_atm_tzonal, pe )
      ELSE
!DEPENDS ON: scatter_field
        CALL Scatter_Field( ptr_field_out % DATA(:,j),              &
                            level_field_out,                        &
                            ptr_field_out % row_len,                &
                            ptr_field_out % rows,                   &
                            ptr_field_out % glob_row_len,           &
                            ptr_field_out % glob_rows,              &
                            fld_type,                               &
                            halo_type_single,                       &
                            pe,                                     & 
                            gc_all_proc_group )
      END IF
      pe = pe + 1
      IF (pe == nproc) pe = 0
    END DO
  END DO

  !-------------------------------------------------------------------
  ! There are rem levels left to process. Will do these now.
  !-------------------------------------------------------------------
  CALL Change_Decomposition( decomp_rcf_input )
  pe = 0
  DO i = 1, rem
    j = nproc * div + i

    IF (ptr_field_in % glob_row_len == 1) THEN ! Zonal Data
      CALL Rcf_Gather_Zonal_Field( ptr_field_in % DATA(:,j),            &
                                   level_field_in,                      &
                                   ptr_field_in % level_size,           &
                                   ptr_field_in % glob_level_size, 1,   &
                                   ppx_atm_tzonal, pe )
    ELSE
!DEPENDS ON: gather_field
      CALL Gather_Field( ptr_field_in % DATA(:,j),                  &
                         level_field_in,                            &
                         ptr_field_in % row_len,                    &
                         ptr_field_in % rows,                       &
                         ptr_field_in % glob_row_len,               &
                         ptr_field_in % glob_rows,                  &
                         fld_type,                                  &
                         halo_type_single,                          &
                         pe,                                        &
                         gc_all_proc_group )
    END IF

    pe = pe + 1
  END DO

  IF (mype < pe) THEN
    CALL Rcf_H_Int_Ctl( ptr_field_out % glob_level_size,                 &
                        ptr_field_in % glob_row_len,                     &
                        ptr_field_out % glob_row_len,                    &
                        ptr_field_in % glob_rows,                        &
                        ptr_field_out % glob_rows,                       &
                        grid_out % global, ptr_aw_index_targ_lhs,        &
                        ptr_aw_index_targ_top, ptr_bl_index_b_l,         &
                        ptr_bl_index_b_r, ptr_bl_index_t_l,              &
                        ptr_bl_index_t_r, ptr_aw_colat_t, ptr_aw_long_l, &
                        level_field_in, ptr_weight_t_r, ptr_weight_b_r,  &
                        ptr_weight_t_l, ptr_weight_b_l, level_field_out )

    !-------------------------------------------------------------------
    ! Coastal Adjustment for land-only, sea-only and T*  fields
    !-------------------------------------------------------------------
    IF (LTimer) CALL Timer( 'Coastal_Adjustment', 103)
    IF (ptr_field_out % stashmaster % grid_type == ppx_atm_compressed .OR.   &
         ptr_field_out % stashmaster % grid_type == ppx_atm_tsea .OR.        &
         ( ptr_field_out % stashmaster % section == stashcode_prog_sec .AND. &
         (ptr_field_out % stashmaster % item    == stashcode_tstar .OR.     &
         (ptr_field_out % stashmaster % item    == stashcode_z0 .AND.       &
          l_roughnesslength_fix )))) THEN
      ! Only able to perform coastal adjustment if land sea mask is present 
      ! in output dump
      IF ( l_lsm_out_present ) THEN 
        DO j = 1, n_coastal_points
          level_field_out( coast_index_out( j ) ) = &
                           level_field_in( coast_index_in( j ) )
        END DO

        IF (coast_adj_method == coast_adj_spiral_method) THEN
                       ! Spiral adjustment
          IF (PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
            WRITE(umMessage,'(A)') "Spiral method is used"
            CALL umPrint(umMessage,src='replanca-rcf_replanca')
          END IF

          maxdim = MIN( grid_out % glob_p_rows,grid_out % glob_p_row_length)
          DO j = 1, grid_out % glob_p_rows * grid_out % glob_p_row_length
            index_targ_sea_tmp( j )  = index_targ_sea( j )
            index_targ_land_tmp( j ) = index_targ_land( j )
            IF ( glob_lsm_out( j ) ) THEN
              lsm_tmp( j ) = 1
            ELSE
              lsm_tmp( j ) = 0
            END IF
          END DO

          sea_points_unres_tmp  = n_sea_points_unres
          land_points_unres_tmp = n_land_points_unres

          IF (ptr_field_out % stashmaster % grid_type /=                   &
                                            ppx_atm_compressed) THEN
            ! Only do the sea coastal points if the field isn't land only.
            ! DEPENDS ON: intf_coast_aj
            CALL Intf_Coast_AJ( lsm_tmp, index_targ_sea_tmp,               &
                 sea_points_unres_tmp, grid_out % glob_p_rows,             &
                 grid_out % glob_p_row_length, level_field_out, 0, cyclic, &
                 maxdim )
          END IF

          ! DEPENDS ON: intf_coast_aj
          CALL Intf_Coast_AJ( lsm_tmp, index_targ_land_tmp,              &
               land_points_unres_tmp, grid_out % glob_p_rows,            &
               grid_out % glob_p_row_length, level_field_out, 1, cyclic, &
               maxdim )

        ELSE IF (coast_adj_method == coast_adj_standard_method) THEN
                       ! Non-spiral adjustment
          IF (PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
            WRITE(umMessage,'(A)') "Standard method is used"
            CALL umPrint(umMessage,src='replanca-rcf_replanca')
          END IF

          DO j = 1, n_land_points_unres
            level_field_out( index_targ_land( j ) ) = &
                             level_field_out( land_unres_index( j ) )
          END DO

          DO j = 1, n_sea_points_unres
            level_field_out( index_targ_sea( j ) ) = &
                             level_field_out( sea_unres_index( j ) )
          END DO

        ELSE IF (coast_adj_method == coast_adj_circle_method) THEN
                       ! Spiral circle method
          IF (PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
            WRITE(umMessage,'(A)') "Spiral circle method is used"
            CALL umPrint(umMessage,src='replanca-rcf_replanca')
          END IF

          ! SPIRAL CIRCLE METHOD

          IF (ptr_field_out % stashmaster % grid_type ==                   &
                                            ppx_atm_compressed) THEN
            ! If land only do not use the 200km constraint
            DO j = 1, n_land_points_unres
              level_field_out( index_targ_land( j ) ) = &
                      level_field_out( land_unres_notconstrain_index( j ) )
            END DO
          ELSE
            ! Not land only
            DO j = 1, n_land_points_unres
              level_field_out( index_targ_land( j ) ) = &
                      level_field_out( land_unres_constrain_index( j ) )
            END DO

            ! Do sea points
            DO j = 1, n_sea_points_unres
              level_field_out( index_targ_sea( j ) ) = &
                    level_field_out( sea_unres_notconstrain_index( j ) )
            END DO
          END IF
        ELSE ! Not a valid coastal_adjustment method
          Cmessage = 'Not a valid coastal adjustment method'
          ErrorStatus = 160
          CALL Ereport( RoutineName, ErrorStatus, Cmessage )
        END IF ! Spiral
      ELSE
        ! Unable to perform coastal adjustment due to lack of land sea mask.
        ! Print warning later in routine, only want to print once per field.
        l_print_coastal_warning = .TRUE.
      END IF
    END IF

    IF (LTimer) CALL Timer( 'Coastal_Adjustment', 104)

    !-----------------------------------------------------------------
    ! Average the polar rows
    !-----------------------------------------------------------------
    IF ( ptr_field_out % stashmaster % grid_type <= 3 .OR.               &
         ptr_field_out % stashmaster % grid_type == 21 ) THEN
      ! Might need to do something here for ENDGAME grid... theta points not at poles
      ! anymore.
      IF ( ptr_field_out % glob_rows > grid_out % glob_v_rows ) THEN
        CALL Rcf_Average_polar( level_field_out, ptr_field_out % glob_rows, &
                                ptr_field_out % glob_row_len,               &
                                grid_out % global, averaged )
      END IF

      IF ( averaged .AND. PrintStatus >= PrStatus_Normal) THEN
        field_averaged = 1
      END IF
    END IF
    IF ( field_out % stashmaster % grid_type == ppx_atm_ozone ) THEN
      IF (field_out % glob_row_len == 1 .AND. field_in % glob_row_len /= 1) THEN
        ! Need to use field_out since we have change stashmaster and row_len to
        ! full field.
        DO j = 1, field_out % rows
          ! We can reuse polar rows and average over 1 row
          CALL Rcf_Average_polar(                                        &
                      level_field_out((j-1)*ptr_field_out % glob_row_len:&
                                      j*ptr_field_out % glob_row_len),   &
                                  1,                                     &
                                  ptr_field_out % glob_row_len,          &
                                  grid_out % global, averaged )

        END DO
      END IF
    END IF

  END IF

  !-------------------------------------------------------------------
  ! And re-scatter data
  !-------------------------------------------------------------------

  CALL Change_Decomposition( decomp_rcf_output )
  pe = 0
  DO i = 1, rem
    j = nproc * div + i
    IF (ptr_field_out % glob_row_len == 1) THEN ! Zonal Data
      CALL Rcf_Scatter_Zonal_Field( ptr_field_out % DATA(:,j),        &
                                level_field_out,                      &
                                ptr_field_out % level_size,           &
                                ptr_field_out % glob_level_size, 1,   &
                                ppx_atm_tzonal , pe )
    ELSE
! DEPENDS ON: scatter_field
      CALL Scatter_Field( ptr_field_out % DATA(:,j),              &
                          level_field_out,                        &
                          ptr_field_out % row_len,                &
                          ptr_field_out % rows,                   &
                          ptr_field_out % glob_row_len,           &
                          ptr_field_out % glob_rows,              &
                          fld_type,                               &
                          halo_type_single,                       &
                          pe,                                     &
                          gc_all_proc_group )
    END IF
    pe = pe + 1
  END DO

  !----------------------------------------------------------------
  ! Print out a message if unable to perform coastal adjustment
  !----------------------------------------------------------------

  IF ( l_print_coastal_warning .AND. mype == 0 .AND.                           &
       PrintStatus >= PrStatus_Normal ) THEN
    WRITE(cmessage,'(A,I2,A,I3,A)')                                   newline//&
      'Unable to perform coastal adjustment for STASHcode section: ',          &
      ptr_field_out % stashmaster % section, ', item: ',                       &
      ptr_field_out % stashmaster % item,                             newline//&
      'as land sea mask is not present in output dump'
    errorstatus = -40
    CALL ereport(routinename, errorstatus, cmessage)
  END IF

  !----------------------------------------------------------------
  ! Print out a message if required for polar averaged fields
  !----------------------------------------------------------------

  CALL gc_imax( 1, nproc, stat, field_averaged )

  IF ( field_averaged == 1 .AND. mype == 0 .AND.                       &
       PrintStatus >= PrStatus_Normal ) THEN
    WRITE(umMessage,*) 'Interpolation has required averaging of polar rows ', &
                'for section: ', ptr_field_out % stashmaster % section,&
                ', item: ', ptr_field_out % stashmaster % item
    CALL umPrint(umMessage,src='rcf_horizontal_mod')
  END IF

  !----------------------------------------------------------------
  ! Deallocate levels
  !----------------------------------------------------------------
  DEALLOCATE( level_field_in )
  DEALLOCATE( level_field_out )

  !-----------------------------------------------------------------
  ! Reverse process above, take fields back to land points
  !-----------------------------------------------------------------

  IF (field_out % stashmaster % grid_type == ppx_atm_compressed) THEN
    DO i = 1, field_out % levels
      CALL compress_to_mask( field_out_tmp % DATA(1:,i),                 &
                           field_out % DATA(1:,i),                     &
                           local_lsm_out, field_out_tmp % level_size, &
                           mask_points )

      IF (mask_points /= field_out % level_size ) THEN
        Cmessage = 'Recompression onto land points - sizes mismatch'
        ErrorStatus = 60
        CALL Ereport( RoutineName, ErrorStatus, Cmessage )
      END IF

    END DO

    ! Release temp. memory
    DEALLOCATE( field_out_tmp % DATA )
    DEALLOCATE( field_in_tmp % DATA )
  ELSE IF ( field_in % stashmaster % grid_type == ppx_atm_ozone .AND. &
           (field_in % glob_row_len  /= 1                        .OR. &
            field_out % glob_row_len /= 1 )) THEN
    IF (field_out % glob_row_len == 1) THEN
      ! Zonal ozone
      DO k = 1, field_out % levels
        DO j = 1, field_out % rows
          ! This has been averaged earlier.
          field_out % DATA(j,k) =                                     &
                field_out_tmp % DATA(1+(j-1)*field_out_tmp%row_len,k)
        END DO
      END DO
    ELSE
      ! Full field
      field_out % DATA(:,:) = field_out_tmp % DATA(:,:)
    END IF
    ! Release temporary stashmasters
    DEALLOCATE( field_out_tmp % stashmaster)
    DEALLOCATE( field_in_tmp % stashmaster)
    ! Release temp. memory
    DEALLOCATE( field_out_tmp % DATA )
    DEALLOCATE( field_in_tmp % DATA )
  END IF

  NULLIFY( ptr_field_out )
  NULLIFY( ptr_field_in  )

  ! Earlier check on land fields might have changed this.  Set it back for
  ! future fields.
  h_int_method = orig_h_int_method

CASE ( interp_no_op)
  ! do nothing

END SELECT

IF (LTimer) CALL Timer( RoutineName, 4)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_horizontal

END MODULE Rcf_horizontal_Mod
