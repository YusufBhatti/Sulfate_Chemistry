! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Extract a STASHmaster record

MODULE Rcf_Exppx_Mod

!  Subroutine Rcf_Exppx   - extract a STASHmaster record.
!
! Description:
!   Given a model, section and item this routine returns a pointer
!   to the specified STASHmaster record
!
! Method:
!   The relevant row in the record array (STM_record) is found
!   by the ppxptr index.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_EXPPX_MOD'

CONTAINS

FUNCTION Rcf_Exppx( Im_ident, section, item, NoFindArg )

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_Ppx_Info_mod, ONLY: &
    STM_record_type,     &
    STM_record,          &
    ppxptr

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Function arguments:
INTEGER, INTENT(IN) :: Im_ident    ! Internal model identifier
INTEGER, INTENT(IN) :: section     ! STASH section no.
INTEGER, INTENT(IN) :: item        ! STASH item no.
LOGICAL, OPTIONAL, INTENT(IN) :: NoFindArg ! Warning with no entry found

! Function value (out)
TYPE (STM_record_type), POINTER :: Rcf_Exppx


! Local variables
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_EXPPX'
CHARACTER (LEN=errormessagelength)     :: Cmessage
INTEGER                      :: row         ! Row no. in PPXI array
INTEGER                      :: ErrorStatus !+ve = fatal error
LOGICAL                      :: NoFind  !local value of NoFindArg
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( PRESENT (NoFindArg)) THEN
  NoFind = NoFindArg
ELSE
  NoFind = .FALSE.
END IF

NULLIFY( Rcf_Exppx )

ErrorStatus = 0
IF ( Im_ident <= 0 .OR. section < 0 .OR. item <= 0) THEN
  IF (Im_ident <= 0) THEN
    WRITE(umMessage,*) 'Exppxi: INVALID Im_ident'
    CALL umPrint(umMessage,src='rcf_exppx_mod')
  END IF
  IF (section  <  0) THEN
    WRITE(umMessage,*) 'Exppxi: INVALID SECTION NO.'
    CALL umPrint(umMessage,src='rcf_exppx_mod')
  END IF
  IF (item     <=  0) THEN
    WRITE(umMessage,*) 'Exppxi: INVALID ITEM NO.'
    CALL umPrint(umMessage,src='rcf_exppx_mod')
  END IF
  WRITE(umMessage,*) &
      'Im_ident ',Im_ident,' section ',section,' item ',item
  CALL umPrint(umMessage,src='rcf_exppx_mod')
  ErrorStatus=1
  Cmessage='ERROR Exppx: INVALID STASH RECORD ID'

  CALL Ereport( RoutineName, ErrorStatus, Cmessage )

ELSE

  ! Obtain row no. in PPXI array
  row = ppxptr(Im_ident,section,item)

  ! Obtain required data value
  IF ( row > 0 ) THEN
    Rcf_Exppx => STM_record( row )
  ELSE

    WRITE (Cmessage, '(A, I5, A, I3, A, I2, A)')                       &
          'Cant find required STASH item ', item,                      &
          ' section ', section,' model ',  Im_ident, ' in STASHmaster'
    IF ( NoFind ) THEN
      ErrorStatus = -2
    ELSE
      ErrorStatus = 2
    END IF

    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END FUNCTION Rcf_Exppx

END MODULE Rcf_Exppx_Mod
