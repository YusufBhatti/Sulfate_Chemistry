! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

MODULE Rcf_Grib_Dest_List_Mod
IMPLICIT NONE

! Description: Routine to DeAllocate all entries in a linked list
!              and nullify the pointers which locate the 'ends'.
!
! Method:
!         For lists with more than 1 member-
!           Step through the list deallocating the previous entry
!           Deallocate the end
!           Nullify the pointers to the end
!         For lists with only 1 member-
!           Deallocate the member.
!           Nullify the pointers.
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GRIB_DEST_LIST_MOD'

CONTAINS

SUBROUTINE Rcf_Grib_Dest_List(List_Header)

USE Rcf_GRIB_Block_Params_Mod, ONLY:                                      &
    List_Marker,                                                           &
    Grib_Record

USE EReport_Mod, ONLY:                                                    &
    EReport

USE UM_ParCore, ONLY:                                                     &
    mype

USE umPrintMgr, ONLY:                                                     &
    umPrint,                                                               &
    umMessage,                                                             &
    PrintStatus,                                                           &
    PrStatus_Diag                   ! =4 Extra Diagnostic output

USE lookup_addresses

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Subroutine arguments

!< Array  arguments with intent(InOut):>
TYPE (List_Marker)               :: List_Header

! Local constants
CHARACTER (LEN=*), PARAMETER     :: RoutineName='RCF_GRIB_DEST_LIST'

! Local variables

TYPE (Grib_Record), POINTER      :: Current

CHARACTER (LEN=errormessagelength) :: Cmessage    ! used for EReport
INTEGER                          :: ErrorStatus   ! used for EReport

INTEGER                          :: cnter
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!=======================================================================
!  Routine Code Start :
!=======================================================================

! Check list has _some_ members.
IF (ASSOCIATED (List_Header % Begin )) THEN

  cnter = 0

  ! Lists with more than 1 member
  IF ( List_Header % LstCount  > 1 ) THEN
    Current => List_Header % Begin

    DO WHILE (ASSOCIATED (Current % Next))
      ! Deallocate array attached to pointer
      IF (ASSOCIATED(Current % VertCoords)) THEN
        DEALLOCATE(Current % VertCoords)
      END IF
      Current => Current % Next
      DEALLOCATE (Current % Prev)
      cnter = cnter + 1
    END DO
  END IF

  ! do for all lists
  NULLIFY (Current)
  DEALLOCATE (List_Header % END)

  IF ( cnter < List_Header % LstCount -1 ) THEN
    WRITE(Cmessage,'(A)') 'Destroyed less list entries than I '//   &
      'thought existed'
    ErrorStatus = 10
    CALL EReport( RoutineName, ErrorStatus, Cmessage )
  ELSE
    IF ( PrintStatus >= PrStatus_Diag .AND. mype == 0) THEN
      WRITE(umMessage,*) "Destroyed a list containing ",cnter +1, " members."
      CALL umPrint(umMessage,src='rcf_grib_dest_list_mod')
    END IF
  END IF

END IF ! list had members

NULLIFY (List_Header % Begin)
NULLIFY (List_Header % END)
List_Header % LstCount = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_Grib_Dest_List
END MODULE Rcf_Grib_Dest_List_Mod
