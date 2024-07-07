! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Interpolate fractions of surface types from input to output grid

MODULE Rcf_Get_Regridded_Tile_Fractions_Mod

! Subroutine Rcf_Get_Regridded_Tile_Fractions
!
! Description:
!   This routine retrieves the field of tile fractions in the output
!   dump (tile_fractions_out) and interpolates the tile fractions
!   from the input to the output grid (input_tilefrac_out_grid).
!
! Method:
!   Self explanatory
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Reconfiguration

IMPLICIT NONE

CONTAINS

SUBROUTINE Rcf_Get_Regridded_Tile_Fractions(                                &
                                      fields_in, field_count_in, hdr_in,    &
                                      fields_out, field_count_out, hdr_out, &
                                      input_tilefrac_out_grid,              &
                                      tile_fractions_out )
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
    PrintStatus,            & ! Not used
    PrStatus_Normal           ! Not used

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field, Rcf_Dealloc_Field

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

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER        :: fields_in(:), fields_out(:)
INTEGER, INTENT(IN)                :: field_count_in, field_count_out
TYPE( um_header_type ), INTENT(IN) :: hdr_in, hdr_out

TYPE( field_type ), POINTER        :: tile_fractions_out
!                                            ! Tile fractions of the output
!                                            ! dump: the routine merely
!                                            ! populates this field
TYPE( field_type )                 :: input_tilefrac_out_grid
!                                            ! Tile fractions of the 
!                                            ! input dump interpolated
!                                            ! to the output grid

! Local variables
TYPE (field_type)                  :: dummy  ! pretend orography
INTEGER                            :: pos    ! field position
INTEGER                            :: orig_h_int_method
!                                            ! Original method of interpolation

CHARACTER (LEN=*), PARAMETER       :: routinename =                    &
                                      'Rcf_Get_Regridded_Tile_Fractions'
CHARACTER (LEN=errormessagelength) :: cmessage

!----------------------------------------------------------------------
! Find surface type fractions in output & input fields and read them in
!----------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_frac_surf_type,         &
                 fields_out, field_count_out, pos)
tile_fractions_out => fields_out(pos)
CALL Rcf_Alloc_Field( tile_fractions_out )
CALL Rcf_Read_Field( tile_fractions_out, hdr_out, decomp_rcf_output )

CALL Rcf_Locate( stashcode_prog_sec, stashcode_frac_surf_type,         &
                 fields_in, field_count_in, pos)
CALL Rcf_Alloc_Field( fields_in(pos) )
CALL Rcf_Read_Field( fields_in(pos), hdr_in, decomp_rcf_input )

!---------------------------------------------------------------------
! Initialise input_tilefrac_out_grid - and allocate space
!----------------------------------------------------------------------
CALL Rcf_Field_Equals( input_tilefrac_out_grid, fields_in(pos) )
input_tilefrac_out_grid % rows         = tile_fractions_out % rows
input_tilefrac_out_grid % row_len      = tile_fractions_out % row_len
input_tilefrac_out_grid % level_size   = tile_fractions_out % level_size
input_tilefrac_out_grid % glob_rows    = tile_fractions_out % glob_rows
input_tilefrac_out_grid % glob_row_len = tile_fractions_out % glob_row_len
input_tilefrac_out_grid % glob_level_size =                                   &
                                          tile_fractions_out % glob_level_size

! Store the original method of interpolation to restore later, then
! set the working method.
orig_h_int_method = fields_in(pos) % interp
IF (h_int_active) THEN
  fields_in(pos) % interp = interp_h_only
ELSE
  fields_in(pos) % interp = interp_copy
END IF

CALL Rcf_Alloc_Field( input_tilefrac_out_grid )

!----------------------------------------------------------------------
! Interpolate input_tilefrac_out_grid field
!----------------------------------------------------------------------
CALL Rcf_Interpolate( fields_in(pos), input_tilefrac_out_grid, input_grid, &
                      output_grid, dummy, dummy )

!----------------------------------------------------------------------
! Tidy Up
!----------------------------------------------------------------------

! Restore the original method of interpolation.
fields_in(pos) % interp = orig_h_int_method

END SUBROUTINE Rcf_Get_Regridded_Tile_Fractions

END MODULE Rcf_Get_Regridded_Tile_Fractions_Mod

