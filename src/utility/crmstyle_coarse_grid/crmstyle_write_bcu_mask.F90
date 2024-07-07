! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  write out coarse grid means and SD to ff

MODULE crmstyle_write_bcu_mask_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_WRITE_BCU_MASK_MOD'

CONTAINS

SUBROUTINE crmstyle_write_bcu_mask(icall_type, maxfldsout, example_hdr,       &
                                 ref_modlev_field,ref_single_field,           &
                                 bcu_mask_hdr,  error_code  )

USE crmstyle_sample_arrays_mod, ONLY:                                       &
  bcu_mask_w 

USE hires_data_mod, ONLY:                                                   &
  iclass_col

USE crmstyle_grid_info_mod, ONLY:                                          &
  local_row_len,local_rows,mask_full_x, mask_full_y

USE crmstyle_cntl_mod, ONLY:                                                &
  num_x,num_y, mlevs, l_class_col, l_bcu_mask

USE word_sizes_mod, ONLY: iwp,wp

USE missing_data_mod, ONLY: rmdi, imdi

USE io, ONLY: file_close

USE io_constants, ONLY: ioNoDelete

USE filenamelength_mod, ONLY: filenamelength

USE file_manager, ONLY: release_file_unit

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type,          &
  UM_Header_type,         &
  lenfixhd

USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal,            &
  EndofFile

USE gather_cp_to_pp_struct_mod, ONLY: gather_cp_to_pp_struct
USE ereport_mod, ONLY: ereport, ereport_finalise

USE umPrintMgr                              ! for writing output
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!  Writeout all the fields averaged for the given partition
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! ------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::              &
  icall_type                        & ! call type
 ,MaxFldsout                          ! max output fields

TYPE(UM_Header_type), INTENT(IN) :: &
  example_hdr                         ! Header to copy

TYPE(PP_Field_type), INTENT(INOUT)  :: &
  ref_modlev_field(mlevs)             ! pp fields on model levels
TYPE(PP_Field_type), INTENT(INOUT)  :: &
  ref_single_field(1)                 ! pp fields single


TYPE(UM_Header_type), INTENT(INOUT)     :: &
  bcu_mask_hdr                                    ! Header to for file

INTEGER, INTENT(OUT) ::    &
  error_code                 ! return code
!-------------------------------------------------------------------------
! Local variables
INTEGER ::               &
  i,j                      ! loop counters

INTEGER ::               &
  ErrorStatus             ! Error code from operations on file

REAL(wp), ALLOCATABLE :: &
  class_col(:,:)           ! real version of integer array


CHARACTER(LEN=*), PARAMETER :: RoutineName = "CRMSTYLE_WRITE_BCU_MASK"

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------


!=============================================================================
error_code = 0

ErrorStatus = 0

IF (icall_type == 1) THEN
  ! Opens and sets up header
  ! Setups new header and allocates lookup table
  ! DEPENDS ON: new_umhdr
  CALL new_umhdr( example_hdr, MaxFldsOut, bcu_mask_hdr, ErrorStatus )

  WRITE(umMessage,'(A,3I6)') ' setup lookup header ',bcu_mask_hdr % NumFlds,  &
                     bcu_mask_hdr % Len2Lookup, ErrorStatus
  CALL umPrint(umMessage,src=RoutineName)
  IF ( ErrorStatus /= StatusOK ) THEN
    ErrorStatus = StatusFatal
    CALL EReport( "crmstyle_write_bcu_mask", ErrorStatus,    &
               "Error Setting up Output File header." )
  END IF

END IF


! Copy field to structure and set stashcode and processing code before
! writing to output fields file.

WRITE(umMessage,'(A,I6,2A)') ' Unit ',bcu_mask_hdr %UnitNum, ' File ',    &
                        TRIM(bcu_mask_hdr % FileName)
CALL umPrint(umMessage,src=RoutineName)

IF (l_bcu_mask) THEN
  CALL gather_cp_to_pp_struct(150, 0, local_row_len,local_rows,              &
                            mask_full_x,mask_full_y,mlevs,maxfldsout,        &
                       bcu_mask_w,ref_modlev_field, bcu_mask_hdr, ErrorStatus )

END IF

! column classification need to first convert to a real array

IF (l_class_col) THEN
  ALLOCATE(class_col(local_row_len,local_rows))
  DO j=1,local_rows
    DO i=1,local_row_len
      class_col(i,j) = REAL(iclass_col(i,j))
    END DO
  END DO
  ! Using stascode  3 476 from combined BL type to indicate column class 
  CALL gather_cp_to_pp_struct(3476, 0, local_row_len,local_rows,             &
                            mask_full_x,mask_full_y, 1,maxfldsout,           &
                       class_col,ref_single_field, bcu_mask_hdr, ErrorStatus )

  DEALLOCATE(class_col)
END IF

! Close output file

IF (icall_type == 3) THEN
  ! Write out updated lookup table
  ! DEPENDS ON: writelookup
  CALL WriteLookup( bcu_mask_hdr, ErrorStatus )
  WRITE(umMessage,'(I7,A)') bcu_mask_hdr % NumFlds, " fields written."
  CALL umPrint(umMessage,src=RoutineName)

  IF ( ASSOCIATED( bcu_mask_hdr % Lookup ) ) THEN
    CALL File_Close ( bcu_mask_hdr % UnitNum,               &
                      bcu_mask_hdr % FileName,              &
                      filenamelength,                       &
                      delete=ioNoDelete, error=ErrorStatus)

    CALL release_file_unit(bcu_mask_hdr % UnitNum, handler="portio")

  END IF
END IF
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_write_bcu_mask

END MODULE crmstyle_write_bcu_mask_mod
