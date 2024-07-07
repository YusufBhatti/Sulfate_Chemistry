! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Checks that snow amount fields contains no negative values.

MODULE Rcf_Snow_Amount_Chk_Mod

!  Subroutine Rcf_Snow_Amount_Chk
!
! Description:
!   Ensures that snow amount field contains no negative values.
!
! Method:
!   Reset any negative values to zero.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SNOW_AMOUNT_CHK_MOD'

CONTAINS

SUBROUTINE Rcf_Snow_Amount_Chk( fields_in, field_count_in, hdr_in, &
                                fields_out, field_count_out, hdr_out )

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Write_Field_Mod, ONLY: &
    Rcf_Write_Field

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

USE decomp_params, ONLY: &
    decomp_rcf_input, decomp_rcf_output

USE um_stashcode_mod, ONLY: &
    stashcode_mean_snow,       &
    stashcode_frac_surf_type,  &
    stashcode_prog_sec

USE UM_ParCore, ONLY: &
    nproc,                  &
    mype

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,                &
    PrStatus_Normal

USE Rcf_LSM_Mod, ONLY: &
    local_lsm_out, local_land_out

USE mask_compression, ONLY: compress_to_mask, expand_from_mask

USE rcf_nlist_recon_science_mod, ONLY: &
    l_regularize_landice,              &
    l_regularize_landice_all_soil_pts, &
    l_regularize_landice_all_lice_pts, &
    snow_landice_min, snow_icefree_max

USE Rcf_Get_Regridded_Tile_Fractions_Mod, ONLY: &
    Rcf_Get_Regridded_Tile_Fractions

USE jules_snow_mod, ONLY: nsmax

USE jules_surface_types_mod, ONLY: ice

USE Ereport_Mod, ONLY: &
    Ereport

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER        :: fields_in(:), fields_out(:)
TYPE( um_header_type ), INTENT(IN) :: hdr_in, hdr_out
INTEGER, INTENT(IN)                :: field_count_in, field_count_out

! Local variables
INTEGER              :: pos_snowa
INTEGER              :: snowa_changed

INTEGER              :: i
INTEGER              :: istat, errorstatus
INTEGER              :: npts, soil_pts_out_snow_excess,            &
                        lice_pts_out_snow_deficit,                 &
                        soil_pts_out_snow_excess_remain,           &
                        lice_pts_out_snow_deficit_remain
INTEGER              :: mask_points            ! For mask compression
INTEGER, ALLOCATABLE :: l_lice_point_igo(:), & ! Integer flag for points on the
                                               ! output grid with a non-zero
                                               ! land-ice fraction obtained
                                               ! from direct interpolation of
                                               ! the input land-ice fraction
                        l_lice_point_out(:), & ! Integer flag for points on the
                                               ! the output grid that actually
                                               ! have a non-zero land-ice
                                               ! fraction
                        l_lice_point_diff(:)   ! Land-ice points new/gone,
                                               ! obtained from differencing the
                                               ! above

REAL, ALLOCATABLE    :: snow_amount_landpts(:,:) ! On land points only

TYPE( field_type ), POINTER  :: snow_amount
TYPE( field_type ), POINTER  :: tile_fractions_out
TYPE( field_type )           :: input_tilefrac_out_grid
!                                              ! Tile fractions of the input
!                                              ! dump interpolated to the
!                                              ! output grid

CHARACTER (LEN=*), PARAMETER       :: RoutineName = 'RCF_SNOW_AMOUNT_CHK'
CHARACTER (LEN=errormessagelength) :: Cmessage
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-------------------------------------------------------------------
! Loacte snow amount in output dump
!-------------------------------------------------------------------

CALL Rcf_Locate ( stashcode_prog_sec, stashcode_mean_snow, &
                  fields_out, field_count_out, pos_snowa, zero_ok_arg = .TRUE. )

!--------------------------------------------------------------------
! Reset snow amount to zero where appropriate
!--------------------------------------------------------------------

snowa_changed = 0

IF (pos_snowa /= 0 ) THEN

  snow_amount => fields_out( pos_snowa )
  CALL Rcf_Alloc_Field( snow_amount )
  CALL Rcf_Read_Field( snow_amount, hdr_out, decomp_rcf_output )

  npts = 0
  DO i = 1, snow_amount % level_size

    ! If snow amount is negative, reset it to 0
    IF ( snow_amount % DATA (i,1) <  0.0 ) THEN

      snow_amount % DATA (i,1) = 0.0
      snowa_changed = 1
      npts = npts + 1

    END IF

    ! If snow amount is positive over sea, reset it to 0
    IF ( snow_amount % DATA (i,1) > 0.0 .AND.  &
                    .NOT. local_lsm_out(i) ) THEN

      snow_amount % DATA (i,1) = 0.0
      snowa_changed = 1
      npts = npts + 1

    END IF

  END DO

  CALL gc_isum (1, nproc, istat, npts)
  IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
    WRITE(umMessage,'(A, I6)') &
        'Snow Amount : No of inappropriate values reset to zero ', npts
    CALL umPrint(umMessage,src='rcf_snow_amount_chk_mod')
  END IF

  ! Adjust the snow amount when changing between ice and ice-free points.
  !
  ! This check is designed only for simple checking with a single ice tile
  ! (as required by the use of the variable ice) and not with the elevated
  ! land ice tiles that are used in some versions of the Earth System Model,
  ! where the multilayer snow scheme should be used anyway. Checks for the
  ! multilayer scheme are carried out in rcf_ml_snowpack.
  !
  ! There seems to be a number of ways to check for presence of land-ice,
  ! but this is consistent with rcf_lsh_land_ice_chk.
  ! Restricting snow on points that have changed from ice to ice-free seems
  ! to miss a large number of points that have excess snow, and likewise for
  ! land-ice points, therefore logicals have been added so that the limits can
  ! be applied to all soil points or all land-ice points instead.

  IF ( l_regularize_landice .AND. nsmax == 0 .AND. ice > 0 ) THEN

    ! Interpolate fractions of surface types on input grid to output grid in
    ! order to identify land-ice points that have changed. On exit,
    ! input_tilefrac_out_grid contains the tile fractions from the input dump
    ! interpolated to the output grid. The tile fractions on the output
    ! grid (tile_fractions_out) are retrieved by the routine.
    CALL Rcf_Get_Regridded_Tile_Fractions(                                &
                                    fields_in, field_count_in, hdr_in,    &
                                    fields_out, field_count_out, hdr_out, &
                                    input_tilefrac_out_grid,              &
                                    tile_fractions_out )

    ! To perform the check we populate two arrays for points on
    ! the output grid: l_lice_point_igo is set to 1 if direct horizontal
    ! interpolation of the tile fractions from the input grid would yield
    ! a non-zero land-ice fraction at this point; l_lice_point_out is set
    ! to 1 if there really is a non-zero land-ice fraction at the point
    ! in the output dump. Points at which these arrays differ indicate
    ! whether land-ice has appeared in a previously ice-free grid-box or
    ! vice versa. Given the large differences in snow loading between 
    ! ice and and ice-free areas, this allows us to determine where the
    ! directly reconfigured snow fields are reasonable and where 
    ! they should be regularized to match the output surface type.

    ! Set up arrays and initialise if needed
    ALLOCATE( l_lice_point_igo  ( local_land_out ) )
    ALLOCATE( l_lice_point_out ( local_land_out ) )
    ALLOCATE( l_lice_point_diff( local_land_out ) )

    l_lice_point_igo(:)  = 0
    l_lice_point_out(:) = 0

    ! Locate land-ice points on the output grid according to 
    ! the directly interpolated tile fractions
    WHERE ( input_tilefrac_out_grid % DATA(:,ice) > 0.0 )
      l_lice_point_igo(:) = 1
    END WHERE
    ! Locate actual land-ice points on the output grid
    WHERE ( tile_fractions_out % DATA(:,ice) > 0.0 )
      l_lice_point_out(:) = 1
    END WHERE
    ! New land ice points have been recreated where l_lice_point_diff
    ! is +1, while old land-ice points have been removed where it is -1.
    ! Otherwise, the value is 0 and no further action is required.
    l_lice_point_diff(:)   = l_lice_point_out(:) - l_lice_point_igo(:)

    ! Need to apply lsm to snow_amount
    ! NB. Size of snow_amount_landpts is still local_lsm_out not local_land_out
    ALLOCATE( snow_amount_landpts( snow_amount % level_size, 1) )
    CALL compress_to_mask( snow_amount % DATA, snow_amount_landpts,       &
                           local_lsm_out, snow_amount % level_size,       &
                           mask_points )

    ! Initialize counters for the numbers of points where changes have
    ! been made.
    soil_pts_out_snow_excess         = 0
    lice_pts_out_snow_deficit        = 0
    soil_pts_out_snow_excess_remain  = 0
    lice_pts_out_snow_deficit_remain = 0
    DO i = 1, local_land_out
      ! Consider first points where there is no land-ice in the output
      ! dump (a soil point) where the snow loading appears excessive.
      ! Regularize the amount if either direct
      ! interpolation of the input grid indicates that there is land-ice,
      ! when direct interpolation of the snow amount is likely to be
      ! excessive, or when the option to regularize all soil points is
      ! selected.
      IF ( l_lice_point_out(i) == 0 .AND. &
         snow_amount_landpts(i,1) > snow_icefree_max ) THEN
        IF ( l_lice_point_diff(i) < 0 .OR. &
           l_regularize_landice_all_soil_pts ) THEN
          snow_amount_landpts(i,1) = snow_icefree_max
          soil_pts_out_snow_excess = soil_pts_out_snow_excess + 1
        ELSE
          soil_pts_out_snow_excess_remain = &
             soil_pts_out_snow_excess_remain + 1
        END IF
      END IF
      ! Next determine points with land-ice in the output dump where the
      ! snow loading appears insufficient. Regularize the amount if either
      ! direct interpolation of the input grid indicates that there isn't
      ! land-ice, when direct interpolation of the snow amount is
      ! likely to yield an inadequate loading, or when the option to do
      ! this for all land-ice points is selected.
      IF ( l_lice_point_out(i) == 1 .AND. &
         snow_amount_landpts(i,1) < snow_landice_min ) THEN
        IF ( l_lice_point_diff(i) > 0 .OR. &
           l_regularize_landice_all_lice_pts ) THEN
          snow_amount_landpts(i,1) = snow_landice_min
          lice_pts_out_snow_deficit = lice_pts_out_snow_deficit + 1
        ELSE
          lice_pts_out_snow_deficit_remain = &
             lice_pts_out_snow_deficit_remain + 1
        END IF
      END IF
    END DO
    CALL gc_isum (1, nproc, istat, soil_pts_out_snow_excess)
    CALL gc_isum (1, nproc, istat, lice_pts_out_snow_deficit)
    CALL gc_isum (1, nproc, istat, soil_pts_out_snow_excess_remain)
    CALL gc_isum (1, nproc, istat, lice_pts_out_snow_deficit_remain)

    IF ( soil_pts_out_snow_excess > 0 ) THEN
      WRITE(umMessage,'(A, I6)') &
         'Global number of soil points with snow > snow_icefree_max reset: ', &
         soil_pts_out_snow_excess
      CALL umPrint(umMessage,src='rcf_snow_amount_chk_mod')
      snowa_changed = 1
    END IF

    IF ( soil_pts_out_snow_excess_remain > 0 ) THEN
      ErrorStatus = -10
      WRITE(cMessage,'(A, I6)') &
         'Global: Excess snow still lingering on soil points ', &
         soil_pts_out_snow_excess_remain
      CALL Ereport ( RoutineName, ErrorStatus, Cmessage )
    END IF

    IF ( lice_pts_out_snow_deficit > 0 ) THEN
      WRITE(umMessage,'(A, I6)') &
         'Global number of land-ice points with ' // &
         'snow < snow_landice_min reset: ', lice_pts_out_snow_deficit
      CALL umPrint(umMessage,src='rcf_snow_amount_chk_mod')
      snowa_changed = 1
    END IF

    IF ( lice_pts_out_snow_deficit_remain > 0 ) THEN
      ErrorStatus = -20
      WRITE(cMessage,'(A, I6)') &
         'Global: Snow deficit remaining on land-ice points ', &
         lice_pts_out_snow_deficit_remain
      CALL Ereport ( RoutineName, ErrorStatus, Cmessage )
    END IF

    CALL Expand_From_Mask( snow_amount % DATA, snow_amount_landpts,       &
                           local_lsm_out, snow_amount % level_size,       &
                           mask_points )

    ! Expand from mask sets sea points to missing data
    WHERE ( .NOT. local_lsm_out(:) ) snow_amount % DATA (:,1) = 0.0

    DEALLOCATE( l_lice_point_igo )
    DEALLOCATE( l_lice_point_out )
    DEALLOCATE( l_lice_point_diff )
    DEALLOCATE( snow_amount_landpts )
    CALL Rcf_Dealloc_Field( tile_fractions_out )
    CALL Rcf_Dealloc_Field( input_tilefrac_out_grid )

  END IF
END IF

!---------------------------------------------------------------------
! Synchronise `changed' flag
!---------------------------------------------------------------------
CALL GC_Imax( 1, nproc, istat, snowa_changed )

!---------------------------------------------------------------------
! Write out changed field
!---------------------------------------------------------------------

IF (snowa_changed == 1) THEN
  CALL Rcf_Write_Field( snow_amount, hdr_out, decomp_rcf_output )
END IF

IF (pos_snowa /= 0) THEN
  CALL Rcf_Dealloc_Field( snow_amount )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Snow_Amount_Chk
END MODULE Rcf_Snow_Amount_Chk_Mod
