! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Wrapper for WRITFLDS

MODULE Rcf_Write_Field_Mod
IMPLICIT NONE

!  Subroutine Rcf_Write_Field
!
! Description:
!    Wrapper for writflds - utilising fields data types
!
! Method:
!    Deals seperately with real, integer and logical data.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_WRITE_FIELD_MOD'

CONTAINS

SUBROUTINE Rcf_Write_Field( field, hdr, decomp, free )


USE Rcf_Alloc_Field_mod, ONLY: &
    Rcf_Dealloc_Field

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Ereport_mod, ONLY: &
    Ereport

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE UM_ParVars, ONLY:   &
    current_decomp_type, &
    change_decomposition

USE cppxref_mod, ONLY: &
    ppx_type_real,      &
    ppx_type_int,       &
    ppx_type_log

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), INTENT(INOUT)    :: field
TYPE( um_header_type), INTENT(IN)    :: hdr
INTEGER, INTENT(IN)                  :: decomp
LOGICAL, OPTIONAL, INTENT(IN)        :: free

! Local vars
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_WRITE_FIELD'
CHARACTER (LEN=errormessagelength)  :: Cmessage
INTEGER                      :: ErrorStatus
INTEGER                      :: i,k
LOGICAL                      :: free_data
INTEGER                      :: orig_decomp ! decomposition when
                                            ! routine is entered
LOGICAL                      :: multi_pe_in = .TRUE.
                                            ! flag to confirm 'gather' req'd
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------
! Change decomposition if required not same as current
!------------------------------------------------------------------
orig_decomp = current_decomp_type
IF (orig_decomp /= decomp) THEN
  CALL Change_Decomposition( decomp )
END IF

!------------------------------------------------------------------
! First set the free_data logical appropriately
!------------------------------------------------------------------
IF ( PRESENT( free ) ) THEN
  free_data = free
ELSE
  free_data = .FALSE.
END IF

!------------------------------------------------------------------
! Now write out the data - need 3 cases real, logical, integer
!------------------------------------------------------------------
SELECT CASE( field % stashmaster % data_type )
CASE (ppx_type_real)
  ! DEPENDS ON: rcf_writflds
  CALL Rcf_WritFlds( hdr % UnitNum, field % levels, field % dump_pos,&
                 hdr % Lookup, hdr % Len1Lookup, field % DATA,       &
                 field % level_size, hdr % FixHd,                    &
                 ErrorStatus, Cmessage , multi_pe_in)

CASE (ppx_type_int)
  ! DEPENDS ON: rcf_writflds
  CALL Rcf_WritFlds( hdr % UnitNum, field % levels, field % dump_pos,&
                 hdr % Lookup, hdr % Len1Lookup, field % Data_Int,   &
                 field % level_size, hdr % FixHd,                    &
                 ErrorStatus, Cmessage , multi_pe_in)

CASE (ppx_type_log)
  ! Convert logical to integer rather than logical to rcf_writflds and 
  ! rely on non-standard conversion of logicals via argument list
  ALLOCATE( field % DATA_int(SIZE(field%data_log,1), SIZE(field%data_log,2)))
  DO k = 1, SIZE(field % DATA_int,2)
    DO i = 1, SIZE(field % DATA_int,1)
      IF (field % data_log(i,k)) THEN
        field % DATA_int (i,k) = 1
      ELSE
        field % DATA_int (i,k) = 0
      END IF
    END DO
  END DO
  ! DEPENDS ON: rcf_writflds
  CALL Rcf_WritFlds( hdr % UnitNum, field % levels, field % dump_pos,&
                 hdr % Lookup, hdr % Len1Lookup, field % DATA_int,   &
                 field % level_size, hdr % FixHd,                    &
                 ErrorStatus, Cmessage , multi_pe_in)
  DEALLOCATE(field % DATA_int)
CASE DEFAULT
  ErrorStatus = -10
  Cmessage = 'Unable to write out field as datatype cannot be '//&
      'determined'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )

END SELECT

!-------------------------------------------------------------------
! Check the returned error conditions
!-------------------------------------------------------------------
IF ( ErrorStatus /= 0 ) THEN
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

!-------------------------------------------------------------------
! Free the allocated data-space if required to
!-------------------------------------------------------------------
IF (free_data) THEN
  CALL Rcf_Dealloc_Field( field )
END IF

!--------------------------------------------------------------------
! Change the decomposition back if required
!---------------------------------------------------------------------
IF (current_decomp_type /= orig_decomp) THEN
  CALL Change_Decomposition( orig_decomp )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_Write_Field
END MODULE Rcf_Write_Field_Mod
