! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  write out plume info to ff

MODULE crmstyle_write_plume_ff_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_WRITE_PLUME_FF_MOD'

CONTAINS

SUBROUTINE crmstyle_write_plume_ff(icall_type,maxfldsout, example_hdr,       &
                                   ref_modlev_field,                         &
                                   n_plume, plume_size, plume_diam_pdf,      &
                                   crm_hdr, error_code )

USE crmstyle_grid_info_mod, ONLY:                                       &
  local_new_x,local_new_y

USE crmstyle_cntl_mod, ONLY:                                              &
  num_x,num_y, mlevs

USE missing_data_mod, ONLY: rmdi, imdi
USE word_sizes_mod, ONLY: iwp,wp

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
!  Writeout all the fields averaged for the gievn partition
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.5.

! ------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::    &
  icall_type              & ! call type
 ,MaxFldsout                ! max output fields


TYPE(UM_Header_type), INTENT(IN)     :: &
  example_hdr                             ! Header to copy
TYPE(PP_Field_type), INTENT(INOUT) ::   &
  ref_modlev_field(mlevs)                 ! pp fields on model levels


REAL(wp), INTENT(IN) ::                         &
  n_plume(local_new_x,local_new_y,mlevs)        & ! number of plumes
 ,plume_size(local_new_x,local_new_y,mlevs)     & ! mean size of a plume
 ,plume_diam_pdf(local_new_x,local_new_y,mlevs)   ! pdf of plume sizes
                                                  ! NOT USED at present

TYPE(UM_Header_type), INTENT(INOUT)     :: &
  crm_hdr                                   ! Header for output file

INTEGER, INTENT(OUT) ::    &
  error_code                 ! return code
  ! error_code is intent out, but is never checked in parent routine.
  ! Should it be deleted ?

!-------------------------------------------------------------------------
! Local variables
INTEGER ::               &
  i, j, k, ij,ii,jj                  ! loop counters

INTEGER ::               &
  ErrorStatus             ! Error code from operations on file


CHARACTER(LEN=*), PARAMETER :: RoutineName = "CRMSTYLE_WRITE_PLUME_FF"

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
! Write already open and setup
ErrorStatus = 0

error_code = 0 ! setting error_code prevents compiler warnings, but it isn't
               ! used. Should it be deleted ?

! Copy field to structure and set stashcode and processing code before
! writing to output fields file.

WRITE(umMessage,'(A,I4,2A)') ' Unit ',crm_hdr %UnitNum,                    &
                             ' File ',TRIM(crm_hdr % FileName)
CALL umPrint(umMessage,src=RoutineName)

! Note using non-standard stashcodes

CALL gather_cp_to_pp_struct(5003, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            n_plume,ref_modlev_field, crm_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(5004, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            plume_size,ref_modlev_field, crm_hdr, ErrorStatus )

! How do I writeout pdf's for each level - want an extra dimension
! Or do I just do this for selected levels ?


!------------------------------------------------------------------------------
! Close output file if completed all processing
!------------------------------------------------------------------------------

IF (icall_type == 3) THEN
  ! Write out updated lookup table
  ! DEPENDS ON: writelookup
  CALL WriteLookup( crm_hdr, ErrorStatus )
  WRITE(umMessage,'(I7,A)') crm_hdr % NumFlds, " fields written."
  CALL umPrint(umMessage,src=RoutineName)

  IF ( ASSOCIATED( crm_hdr % Lookup ) ) THEN
    CALL File_Close ( crm_hdr % UnitNum,               &
                      crm_hdr % FileName,              &
                      filenamelength,                  &
                      delete=ioNoDelete, error=ErrorStatus )

    CALL release_file_unit(crm_hdr % UnitNum, handler="portio")
  END IF

END IF
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_write_plume_ff

END MODULE crmstyle_write_plume_ff_mod

