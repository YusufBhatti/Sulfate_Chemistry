! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Module to allocate and deallocate space for field data

MODULE Rcf_Alloc_Field_Mod
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

PRIVATE :: lhook, dr_hook, jpim, jprb
!  Subroutine Rcf_Alloc_Field       - allocates space
!  Subroutine Rcf_Dealloc_Field     - frees up space
!
! Description:
!  Allocates and deallocates space for field data
!
! Method:
!  Space is allocated to data pointers based on the data type from
!  stashmaster
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_ALLOC_FIELD_MOD'

CONTAINS

SUBROUTINE Rcf_Alloc_Field( field )

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Ereport_mod, ONLY: &
    Ereport

USE cppxref_mod, ONLY: &
    ppx_type_real,     &
    ppx_type_int,      &
    ppx_type_log

USE missing_data_mod, ONLY: &
    imdi,                   &
    rmdi  

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Arguments
TYPE (field_type), INTENT(INOUT)   :: field

! Local data
INTEGER                      :: ErrorStatus
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_ALLOC_FIELD'
CHARACTER (LEN=errormessagelength)           :: Cmessage

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!---------------------------------------------------------------
! Space allocated is local - however must be of the right `type'
!---------------------------------------------------------------

SELECT CASE (field % stashmaster % data_type)
CASE (ppx_type_real)
  IF ( .NOT. ALLOCATED( field % data) ) THEN  
    ALLOCATE( field % data( field % level_size, field % levels ) )
    field % data = rmdi
  END IF

CASE (ppx_type_int)
  IF ( .NOT. ALLOCATED( field % data_int) ) THEN
    ALLOCATE( field % data_int( field % level_size, field % levels ) )
    field % data_int = imdi
  END IF

CASE (ppx_type_log)
  IF ( .NOT. ALLOCATED( field % data_log) ) THEN
    ALLOCATE( field % data_log( field % level_size, field % levels ) )
    field % data_log = .FALSE.
  END IF

CASE DEFAULT
  Cmessage = 'Field data space cannot be allocated as datatype is '//&
      'not recognised'
  ErrorStatus = -10
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_Alloc_Field

!--------------------------------------------------------------------
! Dealloc_field frees up *any* allocated data
!--------------------------------------------------------------------

SUBROUTINE Rcf_Dealloc_Field( field )

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE ereport_mod, ONLY: ereport
IMPLICIT NONE

TYPE (field_type), INTENT(INOUT)    :: field
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_DEALLOC_FIELD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF ( ALLOCATED( field % data ) ) THEN
  DEALLOCATE( field % data )
END IF

IF ( ALLOCATED( field % data_int ) ) THEN
  DEALLOCATE( field % data_int )
END IF

IF ( ALLOCATED( field % data_log ) ) THEN
  DEALLOCATE( field % data_log )
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE Rcf_Dealloc_Field
END MODULE Rcf_Alloc_Field_mod
