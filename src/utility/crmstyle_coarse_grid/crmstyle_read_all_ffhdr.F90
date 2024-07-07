! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  open and read in all fieldsfile headers

MODULE crmstyle_read_all_ffhdr_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_READ_ALL_FFHDR_MOD'

CONTAINS

SUBROUTINE crmstyle_read_all_ffhdr(num_ff,ff_hdr)


USE crmstyle_pp_data_mod, ONLY:                                             &
  iyear,imon,iday,ihour,imin,isec,isyear,ismon,isday,ishour,ismin,issec,    &
  bdx,bdy

USE missing_data_mod, ONLY: rmdi

USE word_sizes_mod, ONLY: iwp,wp    ! Allows use of 4 byte words
USE crmstyle_filenames_mod, ONLY:                                           &
  pp_file, ff_env_name

USE file_manager, ONLY: assign_file_unit

USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal,            &
  EndofFile

USE ereport_mod, ONLY: ereport, ereport_finalise

USE IO_Mod, ONLY: UM_Header_type,PP_Header_type

USE umPrintMgr                              ! for writing output

USE get_env_var_mod, ONLY: get_env_var

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   Open and read in all fieldsfile headers.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.5

! ------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------

INTEGER,INTENT(IN) ::  &
  num_ff                 ! number of fieldsfiles to open and read

TYPE(UM_Header_type), INTENT(INOUT) :: ff_hdr(num_ff)  ! UM Headers: fieldsfile

!-------------------------------------------------------------------------
! Local variables
INTEGER ::               &
  i                        ! loop counters

INTEGER ::               &
  unit_in                & ! unit number
 ,errorstatus              ! return code

LOGICAL :: l_append = .FALSE.

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'CRMSTYLE_READ_ALL_FFHDR'

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
! Do I need to open all the fieldsfiles at the start or have at least 2 files
! open at any time?
! Units numbers from 30 onwards

DO i = 1,num_ff

  CALL get_env_var(ff_env_name(i), ff_hdr(i) % FileName)

  CALL assign_file_unit(ff_hdr(i) % FileName, ff_hdr(i) % UnitNum, &
                        handler="portio")

  WRITE(umMessage,'(A,I3,A)') ' reading file header on unit ',               &
                 ff_hdr(i) % UnitNum , TRIM(ff_hdr(i) % FileName)
  CALL umPrint(umMessage,src=RoutineName)
  ! initialise to ok
  ErrorStatus = StatusOK
  ! DEPENDS ON: read_umhdr
  CALL read_umhdr( ff_hdr(i), l_append, ErrorStatus )

  IF (ErrorStatus /= StatusOK) THEN
    CALL ereport(RoutineName, errorstatus, &
        "Problem reading header from file "//TRIM(ff_hdr(i) % FileName))
  END IF

END DO

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_read_all_ffhdr

END MODULE crmstyle_read_all_ffhdr_mod
