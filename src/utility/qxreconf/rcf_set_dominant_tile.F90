! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Location of dominant tile type for changing to 1 tile.

MODULE Rcf_Set_Dominant_Tile_Mod

! Subroutine Rcf_Set_Dominant_Tile
!
! Description:
!   Returns an array of the dominant tile in each grid box.
!
! Method:
!   The dominant tile at points on the output grid is selected, provided
!   that the same tile has a non-zero fraction on the input grid. This is 
!   to avoid taking garbage from the input dump for tiles that have a
!   zero fraction. (It would perhaps be better to leave the dominant tile
!   from the output dump, but instead when copying the data should find
!   the neareast neighbour where the tile type exists.)
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Reconfiguration

IMPLICIT NONE

INTEGER, ALLOCATABLE, SAVE :: i_dominant(:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SET_DOMINANT_TILE_MOD'

CONTAINS

SUBROUTINE Rcf_Set_Dominant_Tile( fields_in, field_count_in, hdr_in, &
                                  fields_out, field_count_out, hdr_out )

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE decomp_params, ONLY: &
    decomp_rcf_input, decomp_rcf_output

USE Rcf_Field_Type_Mod, ONLY: &
    Field_Type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE um_stashcode_mod, ONLY:   &
    stashcode_frac_surf_type, &
    stashcode_prog_sec

USE umPrintMgr, ONLY:       &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal

USE UM_ParCore, ONLY: &
    mype

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Alloc_Field_Mod, ONLY:  &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Grid_Type_Mod, ONLY:  &
    Input_Grid,               &
    Output_Grid

USE Rcf_Interp_Weights_Mod, ONLY: &
    h_int_active

USE Rcf_Set_Interp_Flags_Mod, ONLY:  &
    interp_h_only,                   &
    interp_copy

USE Rcf_Interpolate_Mod, ONLY:  &
    Rcf_Interpolate

USE Rcf_Field_Equals_Mod, ONLY: &
    Rcf_Field_Equals

USE Rcf_Get_Regridded_Tile_Fractions_Mod, ONLY: &
    Rcf_Get_Regridded_Tile_Fractions

USE missing_data_mod, ONLY: imdi

USE ereport_mod, ONLY: ereport

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
INTEGER, INTENT(IN)                :: field_count_in, field_count_out
TYPE( um_header_type ), INTENT(IN) :: hdr_in, hdr_out

! Local variables
TYPE( field_type ), POINTER        :: tile_fractions_out
TYPE( field_type )                 :: input_tilefrac_out_grid
!                                            ! Tile fractions of the input
!                                            ! dump interpolated to the
!                                            ! output grid
TYPE (field_type)                  :: dummy  ! pretend orography
INTEGER                            :: pos    ! field position
INTEGER                            :: orig_h_int_method
!                                            ! Original method of interpolation
INTEGER                            :: i      ! looper
INTEGER                            :: n      ! Number of dominant tile fails
REAL, ALLOCATABLE                  :: frac_dominant(:)

INTEGER                            :: errorstatus
CHARACTER (LEN=*), PARAMETER       :: RoutineName='RCF_SET_DOMINANT_TILE'
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
  WRITE(umMessage,'(A)') 'Locating dominant tile for setting to 1 tile.'
  CALL umPrint(umMessage,src='rcf_set_dominant_tile_mod')
END IF

! Set up the tile fractions for the output (tile_fractions_out) and
! also set up a field of tile fractions from the input dump interpolated
! to the output grid (input_tilefrac_out_grid); this field retains the tile
! types from the input dump.
CALL Rcf_Get_Regridded_Tile_Fractions( fields_in, field_count_in, hdr_in,    &
                                fields_out, field_count_out, hdr_out, &
                                input_tilefrac_out_grid, tile_fractions_out )

!----------------------------------------------------------------------
! Initialize the arrays
!----------------------------------------------------------------------

ALLOCATE(frac_dominant(SIZE(tile_fractions_out % DATA(:, 1))))

i_dominant(:)    = imdi
frac_dominant(:) = 0.0
DO i = 1, tile_fractions_out % levels
  WHERE ( tile_fractions_out % DATA(:, i) > frac_dominant .AND. &
     input_tilefrac_out_grid % DATA(:, i) > 0.0 )
    i_dominant(:)    = i
    frac_dominant(:) = tile_fractions_out % DATA(:, i)
  END WHERE
END DO

n = COUNT ( i_dominant(:) == imdi )
IF ( n > 0 ) THEN
  WHERE ( i_dominant(:) == imdi )
    i_dominant(:) = MAXLOC( input_tilefrac_out_grid % DATA(:,:), DIM=2 )
  END WHERE
  errorstatus = -6175
  WRITE (cmessage, '(A, I8)')                                                 &
     'Cases of dominant tiles set from input rather than output dump = ', n
  CALL ereport( routinename, errorstatus, cmessage )
END IF

!----------------------------------------------------------------------
! Tidy Up
!----------------------------------------------------------------------

DEALLOCATE(frac_dominant)
CALL Rcf_Dealloc_Field( input_tilefrac_out_grid )
CALL Rcf_Dealloc_Field( tile_fractions_out )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE Rcf_Set_Dominant_Tile

END MODULE Rcf_Set_Dominant_Tile_Mod
