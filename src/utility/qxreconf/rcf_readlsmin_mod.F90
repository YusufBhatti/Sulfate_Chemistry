! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Read in the input grid land-sea masks

MODULE Rcf_ReadLSMIn_Mod
IMPLICIT NONE

!  Subroutine Rcf_ReadLSMin - read land-sea masks
!
! Description:
! This module is designed to read-in and setup an input land-sea
! mask for multiple PEs. IE, both local and global LSMs are read
! and appropriate sizes calculated. It is assumed that space has
! already been allocated for them, but this is checked.
!
! Method:
!  Local LSM is read in on each PE. This is gathered on PE 0 and
!  broadcast so that all pes have a copy of the global LSM.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_READLSMIN_MOD'

CONTAINS

SUBROUTINE Rcf_ReadLSMIn(Dump_hdr, fields, field_count )

USE Rcf_Lsm_Mod, ONLY: &
    glob_lsm_in,    &
    local_lsm_in,   &
    local_land_in,  &
    glob_land_in

USE Rcf_UMhead_Mod, ONLY: &
    UM_header_type

USE Ereport_Mod, ONLY: &
    Ereport

USE UM_ParVars, ONLY: &
    lasize,             &
    glsizep,            &
    gc_all_proc_group,  &
    current_decomp_type,&
    change_decomposition

USE UM_ParCore, ONLY: &
    nproc

USE decomp_params, ONLY: &
    Decomp_RCF_input

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE umPrintMgr, ONLY:      &
    umPrint,                &
    PrStatus_Normal

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Gather_Field_Mod, ONLY: &
    Rcf_Gather_Field

USE um_stashcode_mod, ONLY: &
    stashcode_lsm,             &
    stashcode_prog_sec

USE Field_types

USE UM_ParParams

USE cppxref_mod, ONLY: &
    ppx_atm_compressed

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)               :: field_count
TYPE (field_type), INTENT(INOUT)  :: fields( field_count )
TYPE (UM_header_type), INTENT(IN) :: Dump_hdr

! Local Data
INTEGER                           :: i
INTEGER                           :: msg
INTEGER                           :: Pos
INTEGER                           :: orig_decomp
INTEGER                           :: ErrorStatus
CHARACTER (LEN=errormessagelength)                :: Cmessage
CHARACTER (LEN=*), PARAMETER      :: RoutineName = 'RCF_READLSMIN'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------
! Set decomposition to be input - communications later require it
!------------------------------------------------------------------
orig_decomp = current_decomp_type
IF ( orig_decomp /= decomp_rcf_input ) THEN
  CALL Change_Decomposition( decomp_rcf_input )
END IF

!-------------------------------------------------------------------
! Check that relevant memory is available
!-------------------------------------------------------------------
IF ( .NOT. ALLOCATED( Local_LSM_In) ) THEN
  Cmessage = 'Local LSM space not allocated!'
  ErrorStatus = 20
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

IF (.NOT. ALLOCATED( Glob_LSM_In ) ) THEN
  Cmessage = 'global LSM space not allocated!'
  ErrorStatus = 30
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

!-------------------------------------------------------------------
! Find Land-Sea Mask amongst input fields
! Allow return code of pos = 0
!------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_lsm,                  &
                 fields, field_count, pos, zero_ok_arg = .TRUE. )

IF (pos == 0) THEN
  ! Default - set land sea mask to .FALSE. everywhere
  Local_LSM_In(:) = .FALSE.
  Glob_LSM_In(:)  = .FALSE.

ELSE
  !------------------------------------------------------------------
  ! Read Local LSM on all PEs
  !------------------------------------------------------------------
    ! allocate space in data type and copy data to Local_Lsm_In
  CALL Rcf_Alloc_Field( fields(pos) )
  CALL Rcf_Read_Field( fields(pos), Dump_hdr, decomp_rcf_input )
  Local_Lsm_In(:) = fields(pos) % Data_Log(:,1)
  CALL Rcf_Dealloc_Field( fields(pos) )

  !--------------------------------------------------------------------
  ! And need global LSM
  !-------------------------------------------------------------------
  CALL Rcf_Gather_Field( local_lsm_in,          glob_lsm_in,        &
                         fields(pos) % row_len, fields(pos) % rows, &
                         fields(pos) % glob_row_len,                &
                         fields(pos) % glob_rows, 0,                &
                         gc_all_proc_group )

  ! Note that this communication should ideally use GC_BBcast and
  ! use Kind to determine number of bytes, but we'll not bother.
  ! After all Gather_Field (above) and the I/O assume that a
  ! Logical is the same size as a Real as is an Integer. We will do
  ! likewise.

  msg = 802
  CALL GC_IBcast( msg, fields(pos) % glob_level_size, 0, nproc,  &
                  ErrorStatus, glob_lsm_in )

  IF ( ErrorStatus /= 0 ) THEN
    Cmessage = 'Problem broadcasting global land-sea mask from PE0'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

END IF
!-------------------------------------------------------------------
! Count the local number of land points for addressing
!-------------------------------------------------------------------

local_land_in = 0

DO i = 1, lasize(1,fld_type_p,halo_type_single) * &
    lasize(2,fld_type_p,halo_type_single)
  IF ( Local_Lsm_In(i) ) THEN
    local_land_in = local_land_in + 1
  END IF
END DO

glob_land_in = 0
DO i = 1, glsizep(1) * glsizep(2)
  IF ( Glob_Lsm_In(i) ) THEN
    glob_land_in = glob_land_in + 1
  END IF
END DO

!-----------------------------------------------------------------
! Adjust the land-only fields sizes to match the above sizes
!-----------------------------------------------------------------
DO i = 1, field_count
  IF ( fields( i ) % stashmaster % grid_type == ppx_atm_compressed) THEN
    IF (pos == 0) THEN       ! have land only points but no lsm to match
      ErrorStatus = 40
      Cmessage = 'Land-Only fields exist in input dump but *NO* '//&
        'Land-Sea Mask'
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF

    fields( i ) % level_size = local_land_in
  END IF
END DO

!-----------------------------------------------------------------
! Reset to original decomposition
!-----------------------------------------------------------------
IF ( orig_decomp /= current_decomp_type ) THEN
  CALL Change_Decomposition( orig_decomp )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_ReadLSMIn

END MODULE Rcf_ReadLSMIn_Mod
