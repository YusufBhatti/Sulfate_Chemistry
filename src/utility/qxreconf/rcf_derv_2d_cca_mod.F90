! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Deriving 2d convective cloud amounts from 2d.

MODULE Rcf_Derv_2D_CCA_Mod

!  Subroutine Rcf_Derv_2D_CCA
!
! Description:
!    Top level routine obtaining variables (and some conversions)
!    for the Rcf_Calc_2D_CCA routine.
!    To avoid a 3D interpolation, calculation is done with
!    input grid variables and the result is interpolated.
!    Some correction of ccb and cct is then required.
!
! Method:
!
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Reconfiguration
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_DERV_2D_CCA_MOD'

CONTAINS

SUBROUTINE Rcf_Derv_2D_CCA( fields_in, fields_out, field_count_in, &
                            field_count_out, hdr_in, hdr_out, cca_2d)

USE Rcf_Field_Type_Mod, ONLY: &
    Field_Type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Grid_Type_Mod, ONLY: &
    Input_Grid,           &
    Output_Grid

USE Rcf_Alloc_Field_mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE decomp_params, ONLY: &
    decomp_rcf_output,   &
    decomp_rcf_input

USE Rcf_Set_Interp_Flags_Mod, ONLY: &
    interp_h_only,                   &
    interp_copy

USE Rcf_Interp_Weights_Mod, ONLY: &
    h_int_active

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Write_Field_Mod, ONLY: &
    Rcf_Write_Field

USE Rcf_Field_Equals_Mod, ONLY: &
    Rcf_Field_Equals

USE Rcf_Interpolate_Mod, ONLY: &
    Rcf_Interpolate

USE Rcf_Calc_2d_CCA_Mod, ONLY: &
    Rcf_Calc_2d_CCA

USE UM_ParCore, ONLY: &
    mype

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal

USE um_stashcode_mod, ONLY: &
    stashcode_3d_cca,          &
    stashcode_ccb,             &
    stashcode_cct,             &
    stashcode_prog_sec

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER        :: fields_in(:)
TYPE( field_type ), POINTER        :: fields_out(:)
TYPE( field_type ), INTENT(INOUT)  :: cca_2d
TYPE( um_header_type ), INTENT(IN) :: hdr_in
TYPE( um_header_type ), INTENT(IN) :: hdr_out
INTEGER, INTENT(IN)                :: field_count_in
INTEGER, INTENT(IN)                :: field_count_out

! Local variables
INTEGER                       :: i           ! looper
INTEGER                       :: pos         ! field position

TYPE( field_type ), POINTER   :: cca_3d
TYPE( field_type ), POINTER   :: ccb
TYPE( field_type ), POINTER   :: cct
TYPE( field_type ), POINTER   :: ccb_in
TYPE( field_type ), POINTER   :: cct_in
TYPE( field_type )            :: dummy
TYPE( field_type )            :: cca_2d_in

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_DERV_2D_CCA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!--------------------------------------------------------------
! Write out out action if appropriate
!--------------------------------------------------------------
IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  WRITE(umMessage,*) 'Initialising 2D CCA'
  CALL umPrint(umMessage,src='rcf_derv_2d_cca_mod')
END IF

!--------------------------------------------------------------
! Find and setup 3D CCA from input
!--------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_3d_cca,               &
                 fields_in, field_count_in, pos )
cca_3d => fields_in(pos)
CALL Rcf_Alloc_Field( cca_3d )
CALL Rcf_Read_Field( cca_3d, hdr_in, decomp_rcf_input )

!---------------------------------------------------------------
! Set up the cca_2d_in field
!---------------------------------------------------------------
CALL rcf_field_equals( cca_2d_in, cca_2d )
cca_2d_in % rows            = cca_3d % rows
cca_2d_in % row_len         = cca_3d % row_len
cca_2d_in % level_size      = cca_3d % level_size
cca_2d_in % glob_rows       = cca_3d % glob_rows
cca_2d_in % glob_row_len    = cca_3d % glob_row_len
cca_2d_in % glob_level_size = cca_3d % glob_level_size

CALL Rcf_Alloc_Field( cca_2d_in )
!---------------------------------------------------------------
! Read conv cloud base and top as these are required
!---------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_ccb,                  &
                 fields_in, field_count_in, pos )
ccb_in => fields_in(pos)
CALL Rcf_Alloc_Field( ccb_in )
CALL Rcf_Read_Field( ccb_in, hdr_in, decomp_rcf_input )

CALL Rcf_Locate( stashcode_prog_sec, stashcode_cct,                  &
                 fields_in, field_count_in, pos )
cct_in => fields_in(pos)
CALL Rcf_Alloc_Field( cct_in )
CALL Rcf_Read_Field( cct_in, hdr_in, decomp_rcf_input )

!-----------------------------------------------------------------
! Have all fields we require - call the calculation routine
!-----------------------------------------------------------------
CALL Rcf_Calc_2D_CCA( cca_3d, ccb_in, cct_in, cca_2d_in )

!--------------------------------------------------------------
! Get rid of the fields we no longer require
!--------------------------------------------------------------
CALL Rcf_Dealloc_Field( cca_3d )
CALL Rcf_Dealloc_Field( ccb_in )
CALL Rcf_Dealloc_Field( cct_in )


!--------------------------------------------------------------
! We now need to interpolate cca_2d_in to the output grid
!--------------------------------------------------------------
IF (h_int_active) THEN
  cca_2d_in % interp = interp_h_only
ELSE
  cca_2d_in % interp = interp_copy
END IF

CALL Rcf_Interpolate( cca_2d_in, cca_2d, input_grid, output_grid, &
                      dummy, dummy)

!----------------------------------------------------------------
! Interpolation can mean that ccb, cct and cca_2d are not
! consistent. A quick `clean up' should prevent later crashes.
!----------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_ccb,                  &
                 fields_out, field_count_out, pos )
ccb => fields_out(pos)
CALL Rcf_Alloc_Field( ccb )
CALL Rcf_Read_Field( ccb, hdr_out, decomp_rcf_output )

CALL Rcf_Locate( stashcode_prog_sec, stashcode_cct,                  &
                 fields_out, field_count_out, pos )
cct => fields_out(pos)
CALL Rcf_Alloc_Field( cct )
CALL Rcf_Read_Field( cct, hdr_out, decomp_rcf_output )

DO i = 1, cca_2d % level_size
  ! Make sure cca is +ive
  IF (cca_2d % DATA(i,1) < 0.0 ) THEN
    cca_2d % DATA(i,1) = 0.0
  END IF

  ! Make sure cca only exists where ccb and cct are both non-zero
  IF ( ccb % Data_Int(i,1) == 0 .AND.     &
       cct % Data_Int(i,1) == 0 .AND.     &
       cca_2d % DATA(i,1) > 0.0 ) THEN
    cca_2d % DATA(i,1) = 0.0
  END IF
END DO

!-------------------------------------------------------------
! Write out the required fields and tidy up memory
!-------------------------------------------------------------
CALL Rcf_Write_Field( ccb, hdr_out, decomp_rcf_output )
CALL Rcf_Write_Field( cct, hdr_out, decomp_rcf_output )

CALL Rcf_Dealloc_Field( cca_2d_in )
CALL Rcf_Dealloc_Field( ccb )
CALL Rcf_Dealloc_Field( cct )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Derv_2D_CCA
END MODULE Rcf_Derv_2D_CCA_Mod
