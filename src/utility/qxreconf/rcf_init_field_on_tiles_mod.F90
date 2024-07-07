! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Initialisation of a field on tiles from the same field, not tiled.

MODULE Rcf_Init_Field_On_Tiles_Mod

! Subroutine Rcf_Init_Field_On_Tiles
!
! Description:
!   Initialises a field on land points AND tiles from the same field in the
!   input dump, which can either be a field just on land points (1-to-n) OR
!   a tiled field with a different tile configuration (m-to-n). It also
!   initialises an aggregate tile from a field on multiple tiles (n-to-1).
!
! Method:
!   Uses the tile mapping array to copy the input field into the output field.
!   For 1-to-n case each element of the mapping array is 1. In the case of the
!   aggregate tile the dominant tile is used.
!   (This routine was modelled on Rcf_Init_Canopy_Water.)
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Reconfiguration

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_INIT_FIELD_ON_TILES_MOD'

CONTAINS

SUBROUTINE Rcf_Init_Field_on_tiles( fields_in,  field_count_in,  hdr_in, &
                                  ntiles_in, fld_tiles, field_stashcode )

USE Rcf_Field_Type_Mod, ONLY: &
    Field_Type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE um_stashcode_mod, ONLY: &
    stashcode_prog_sec, stashcode_tstar_tile

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal

USE UM_ParCore, ONLY: &
    mype

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE rcf_nlist_recon_science_mod, ONLY: &
    l_use_zero_frac_tile_temp

USE Rcf_Interpolate_Mod, ONLY: &
    Rcf_Interpolate

USE Rcf_Init_Tile_T_Mod, ONLY: &
    Rcf_Init_Tile_T

USE Rcf_Field_Equals_Mod, ONLY: &
    Rcf_Field_Equals

USE nlsizes_namelist_mod, ONLY: &
    ntiles

USE decomp_params, ONLY: &
    decomp_rcf_input

USE Rcf_Grid_Type_Mod, ONLY: &
    Input_Grid,               &
    Output_Grid

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Interp_Weights_Mod, ONLY: &
    h_int_active

USE Rcf_Set_Interp_Flags_Mod, ONLY: &
    interp_h_only,                   &
    interp_copy

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Lsm_Mod, ONLY: &
    local_lsm_in

USE Rcf_Set_Dominant_Tile_Mod, ONLY:       &
   i_dominant

USE um_stashcode_mod, ONLY:       &
   stashcode_snow_laythk_tiles,   &
   stashcode_snow_ice_tile,       &
   stashcode_snow_liq_tile,       &
   stashcode_snow_t_tile,         &
   stashcode_snow_laydns_tiles,   &
   stashcode_snow_grnsiz_tiles

USE jules_snow_mod, ONLY : nsmax

USE jules_surface_types_mod, ONLY: tile_map_pslevs

USE ereport_mod, ONLY:                     &
    ereport

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER          :: fields_in(:)
TYPE( field_type ), INTENT(INOUT)  :: fld_tiles
TYPE( um_header_type ), INTENT(IN) :: hdr_in

INTEGER, INTENT(IN)                :: field_count_in
INTEGER, INTENT(IN)                :: field_stashcode
INTEGER, INTENT(IN)                :: ntiles_in ! Number of tiles in input dump

! Local variables
TYPE (field_type), POINTER         :: fld_tiles_in
TYPE (field_type)                  :: fld_tiles_out
TYPE (field_type)                  :: dummy  ! pretend orography
INTEGER                            :: tile_map_pslevs_tmp(ntiles)
INTEGER                            :: pos    ! field position
INTEGER                            :: i,j    ! looper

LOGICAL :: l_init_tile_t_zerofrac

INTEGER                            :: errorstatus
CHARACTER (LEN=*), PARAMETER       :: RoutineName='RCF_INIT_FIELD_ON_TILES'
CHARACTER (LEN=errormessagelength) :: cmessage

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  IF ( MAXVAL(tile_map_pslevs) < 0 ) THEN
    IF ( ALLOCATED(i_dominant) ) THEN
      WRITE(umMessage,'(A,A,I6)') 'Setting single-layer field on tiles from ', &
         'dominant tile of input multi-layer field: stashcode = ',&
         field_stashcode
    ELSE
      WRITE(umMessage,'(A,A,I6)') 'Initialising field on tiles from ', &
         'input single-layer field: stashcode = ',field_stashcode
    END IF
  ELSE
    WRITE(umMessage,'(A,A,I6)') 'Mapping field on tiles from ',   &
       'input multi-layer field: stashcode = ',field_stashcode
  END IF
  CALL umPrint(umMessage,src='rcf_init_field_on_tiles_mod')
END IF

!----------------------------------------------------------------------
! Find single-layer field in input fields and read it in
!----------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, field_stashcode,         &
                 fields_in, field_count_in, pos)
fld_tiles_in => fields_in(pos)
CALL Rcf_Alloc_Field( fld_tiles_in )
CALL Rcf_Read_Field( fld_tiles_in, hdr_in, decomp_rcf_input )

!----------------------------------------------------------------------
! Initialise fld_tiles_out - and allocate space
!----------------------------------------------------------------------
CALL Rcf_Field_Equals( fld_tiles_out, fld_tiles_in )
fld_tiles_out % rows            = fld_tiles % rows
fld_tiles_out % row_len         = fld_tiles % row_len
fld_tiles_out % level_size      = fld_tiles % level_size
fld_tiles_out % glob_rows       = fld_tiles % glob_rows
fld_tiles_out % glob_row_len    = fld_tiles % glob_row_len
fld_tiles_out % glob_level_size = fld_tiles % glob_level_size

IF (h_int_active) THEN
  fld_tiles_in % interp = interp_h_only
ELSE
  fld_tiles_in % interp = interp_copy
END IF

CALL Rcf_Alloc_Field( fld_tiles_out )

!----------------------------------------------------------------------
! Interpolate input field to dummy single level output field
! and then copy result to all ntiles pseudolevels.
!----------------------------------------------------------------------

! For tstar_tile, overwrite zero fractions with sensible values before
! interpolation
IF ( .NOT. l_use_zero_frac_tile_temp .AND.                              &
     fld_tiles_in % stashmaster % section == stashcode_prog_sec .AND.   &
     fld_tiles_in % stashmaster % item == stashcode_tstar_tile ) THEN

  l_init_tile_t_zerofrac = .TRUE.

  CALL Rcf_Init_Tile_T(fields_in, field_count_in, hdr_in,               &
                       decomp_rcf_input, local_lsm_in,                  &
                       fld_tiles_in, l_init_tile_t_zerofrac)

END IF

CALL Rcf_Interpolate( fld_tiles_in, fld_tiles_out, input_grid, &
                      output_grid, dummy, dummy )

! Compatability with m-to-n reconfiguration:
! tile_map_pslevs not set therefore called from either 1-to-n or n-to-1
! rcf_calc_tiles. Set tile_map_pslevs to map output tiles to single tile in
! input. ( Not needed if called from n-to-1, but does no harm )
tile_map_pslevs_tmp(:)=tile_map_pslevs(1:ntiles)
IF ( MAXVAL(tile_map_pslevs) < 0 ) tile_map_pslevs_tmp(:) = 1

IF ( ntiles > 1 ) THEN
  DO i = 1, ntiles
    fld_tiles % DATA(:,i) = fld_tiles_out % DATA(:,tile_map_pslevs_tmp(i))
    SELECT CASE ( field_stashcode )
    CASE ( stashcode_snow_laythk_tiles,   &
           stashcode_snow_ice_tile,       &
           stashcode_snow_liq_tile,       &
           stashcode_snow_t_tile,         &
           stashcode_snow_laydns_tiles,   &
           stashcode_snow_grnsiz_tiles )
      IF ( nsmax > 1 ) THEN
        DO j = 2, nsmax
          fld_tiles % DATA(:,i + (j-1) * ntiles) = &
             fld_tiles_out % DATA(:,tile_map_pslevs_tmp(i) + (j-1) * ntiles_in)
        END DO
      END IF
    END SELECT
  END DO
ELSE IF ( ALLOCATED(i_dominant) ) THEN
  DO j = 1, fld_tiles % levels ! 1 or nsmax as ntiles = 1
    DO i = 1, fld_tiles % level_size ! Grid
      fld_tiles % DATA(i,j) =  &
         fld_tiles_out % DATA(i,i_dominant(i) + (j-1) * ntiles_in)
    END DO
  END DO
ELSE
  errorstatus = 6173
  WRITE (cmessage, '(A, I5)')                                  &
     'l_aggregate although dominant tile array is not allocated for: ', &
     field_stashcode
  CALL ereport( routinename, errorstatus, cmessage )
END IF

!----------------------------------------------------------------------
! Tidy up
!----------------------------------------------------------------------
CALL Rcf_Dealloc_Field( fld_tiles_out )
CALL Rcf_Dealloc_Field( fld_tiles_in )
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Rcf_Init_Field_On_Tiles

END MODULE Rcf_Init_Field_On_Tiles_Mod
