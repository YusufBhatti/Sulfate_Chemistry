! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  read in orog orog
MODULE crmstyle_read_orog_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_READ_OROG_MOD'

CONTAINS

SUBROUTINE crmstyle_read_orog(in_cols, in_rows, mype, all_proc_group, &
                              local_row_len,local_rows)

USE hires_data_mod , ONLY:                                              &
    orog
USE crmwork_arrays_mod , ONLY:                                          &
    orog_full


USE io

USE file_manager, ONLY: assign_file_unit

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type,          &
  UM_Header_type

USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal,            &
  EndofFile

USE UM_ParParams, ONLY: halo_type_no_halo,halo_type_single, fld_type_p

USE ereport_mod, ONLY: ereport, ereport_finalise

USE get_env_var_mod, ONLY: get_env_var

! Subroutine
USE get_anc_flds_mod, ONLY: get_anc_flds

USE umPrintMgr                              ! for writing output

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   Read Land sea orog ancillary file
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.6

! ------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::   &
  in_cols                & ! Number of columns
 ,in_rows                & ! Number of rows
 ,mype                   & ! Processor number
 ,all_proc_group         & ! Group id for all processors
 ,local_row_len          & ! local row length
 ,local_rows               ! local rows

!-------------------------------------------------------------------------
! Local variables
INTEGER ::               &
  i, j, k                  ! loop counters

INTEGER ::               &
  ErrorStatus            & ! Error code from operations on file
 ,numflds                & ! number of fields
 ,MaxFlds

INTEGER ::       &
  Store(1)       &  !
 ,STCode(1)         ! stashcode to search for

REAL ::    &
  bzy, bzx, bdy, bdx

REAL :: &
  orog_temp(local_row_len,local_rows)  ! orog for this PE as 64 bit field
REAL :: &
  orog_full_temp(in_cols,in_rows)  ! orog for this PE as 64 bit field

LOGICAL :: PPHdrMod = .FALSE.
LOGICAL :: l_append = .FALSE.

TYPE(PP_Field_type)  :: orogfield(1)      ! orography
TYPE(UM_Header_type) :: orog_hdr       ! UM Headers:  orogaphy file

CHARACTER(LEN=*), PARAMETER :: RoutineName = "CRMSTYLE_READ_OROG"

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------

! read header record
! Not doing any special checking assuming the file being pointed to is correct

! Locate and read in land sea orog field
CALL get_env_var("OROG_FILE", orog_hdr % FileName)

CALL assign_file_unit(orog_hdr % FileName, orog_hdr % UnitNum, handler="portio")

! DEPENDS ON: read_umhdr
CALL read_umhdr( orog_hdr, l_append, ErrorStatus )

IF ( ErrorStatus == StatusOK ) THEN

  ! read in orog

  maxflds = 1
  NumFlds       =       1
  Store     (1) =       1
  STCode    (1) =      33 ! code for orography

  CALL set_unit_bcast_flag(orog_hdr % UnitNum)
  IF (mype == 0) THEN
    CALL get_anc_flds( NumFlds, MaxFlds, STCode,Store, PPHdrMod, orog_hdr,   &
                      bzy,bzx,bdy,bdx,                                       &
                     orogfield, ErrorStatus )
    IF ( ErrorStatus /= StatusOK ) THEN
      ErrorStatus = StatusFatal
      CALL EReport( RoutineName, ErrorStatus,              &
                 "Error Reading orography field from ancillary" )
    END IF

    ! copy field into where I want to hold it
    DO j= 1,in_rows
      DO i= 1,in_cols
        orog_full(i,j) = orogField(1) % RData(i,j)        ! 32 bit copy
        orog_full_temp(i,j) = orogField(1) % RData(i,j)   ! 64 bit for scatter
      END DO
    END DO
    ! Release space
    DEALLOCATE( orogfield(1) % RData )
    NULLIFY( orogField(1) % RData )

  END IF   ! mype == 0
  CALL clear_unit_bcast_flag(orog_hdr % UnitNum)

  ! Full orography on PE 0 - Want to scatter back
  ! Also keep for use when calculating weights for use when reading in other
  ! fields on PE 0.

  ! DEPENDS ON: scatter_field
  CALL scatter_field(orog_temp,orog_full_temp,                    &
                  local_row_len,local_rows,                       &
                  in_cols,in_rows,                                &
                  fld_type_p,halo_type_no_halo,                   &
                  0,all_proc_group)

  ! 64 bit to 32 bit

  DO j= 1,local_rows
    DO i= 1, local_row_len
      orog(i,j) =  orog_temp(i,j)
    END DO
  END DO

ELSE ! problems reading file

  ErrorStatus = StatusFatal
  CALL EReport( RoutineName, ErrorStatus,              &    ! error - no input
                 "Error Reading orography file Header" )
END IF

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_read_orog

END MODULE crmstyle_read_orog_mod
