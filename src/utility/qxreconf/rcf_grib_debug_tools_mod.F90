! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Diagnostic tools used with GRIB data

MODULE Rcf_Grib_Debug_Tools_Mod
!
! Description:
!   Small selection of routines allowing me to print information from
!   the linked lists set up to handle incoming GRIB data
!
! *** NOTE *** - Routines containrd herin have unprotected write
!                statements, therefore calls to this routine should be
!                protected by normal procedures for write statements
!
!-----------------------------------------------------------------------
!  Variables to do specifically with the GRIB record
!-----------------------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

USE Rcf_GRIB_Block_Params_Mod

USE Rcf_GRIB_FldSort_Mod, ONLY:  &
  Rcf_GRIB_FldSort

USE um_stashcode_mod

USE Rcf_GRIB_Lookups_Mod               ! Contains cross ref table for
                                       ! Stash => ECMWF parameter Id's
                                       ! and appropriate params
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

PRIVATE :: lhook, dr_hook, jpim, jprb

TYPE (Grib_Record),POINTER       :: Current       !\ Pointer to current
                                                  !/ grib record

! Function definitions
!=======================================================================
! *** Actual Subroutines Involved                                  ***
!=======================================================================

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GRIB_DEBUG_TOOLS_MOD'

CONTAINS
!=======================================================================
! **********************************************************************
!  Print some basic information on the entries in                      *
!  each list (except the 'unknown parameters' one                      *
! **********************************************************************
!=======================================================================

SUBROUTINE Grib_Debug_Print_Basics(Lists)

USE umPrintMgr, ONLY: &
    umPrint,           &
    umMessage

IMPLICIT NONE

!An array of pointer pairs to form the head and tail of the field lists
TYPE (List_Marker)               :: Lists(0:grib_max_fields)

INTEGER                          :: i
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GRIB_DEBUG_PRINT_BASICS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!=======================================================================
! Loop through lists displaying basic info              ! debug
!=======================================================================

DO i = 1, grib_max_fields                    ! Loop through all lists

  IF (.NOT. ASSOCIATED(Lists(i)%Begin)) THEN    ! check list has a start

    WRITE(umMessage,*) " List is Unassociated and contains no entrys"
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')

  ELSE                                        ! List _has_ a start entry

    WRITE(umMessage,'(3A,I0,A)') &
               "List ", TRIM(Lists(i) % Begin % Desc) , " contains ",  &
                Lists( i) % LstCount, " members"
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')

    Current => Lists(i) % Begin
    WRITE(umMessage,'(A15,1x,":")') "Parameter ID is"
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
    DO WHILE (ASSOCIATED(Current))
      WRITE(umMessage,'(2x,I0)') Current % Block_1 ( p_Param_ID )
      CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
      Current => Current % Next
    END DO   ! while associated current
    WRITE(umMessage,*)
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')

    Current => Lists(i) % Begin
    WRITE(umMessage,'(A15,1x,":")') "Type of Levl is"
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
    DO WHILE (ASSOCIATED(Current))
      WRITE(umMessage,'(2x,I0)') Current % Block_1 ( p_Lvl_Type )
      CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
      Current => Current % Next
    END DO   ! while associated current
    WRITE(umMessage,*)
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')

    Current => Lists(i) % Begin
    WRITE(umMessage,'(A15,1x,":")') "1st Lvl Desc is"
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
    DO WHILE (ASSOCIATED(Current))
      WRITE(umMessage,'(1x,I0)') Current % Block_1(p_Lvl_Desc_1)
      CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
      Current => Current % Next
    END DO   ! while associated current
    WRITE(umMessage,*)
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')

    Current => Lists(i) % Begin
    WRITE(umMessage,'(A15,1x,":")') "2nd Lvl Desc is"
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
    DO WHILE (ASSOCIATED(Current))
      WRITE(umMessage,'(1x,I0)') Current % Block_1(p_Lvl_Desc_2)
      CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
      Current => Current % Next
    END DO   ! while associated current
    WRITE(umMessage,*)
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')

    Current => Lists(i) % Begin
    WRITE(umMessage,'(A15,1x,":")') "Start point is "
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
    DO WHILE (ASSOCIATED(Current))
      WRITE(umMessage,'(1x,I0)') Current % Start_pos
      CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
      Current => Current % Next
    END DO   ! while associated current
    WRITE(umMessage,*)
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')

    Current => Lists(i) % Begin
    WRITE(umMessage,'(A15,1x,":")') "Record Length  "
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
    DO WHILE (ASSOCIATED(Current))
      WRITE(umMessage,'(1x,I0)') Current % Block_0(p_Mes_Len)
      CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
      Current => Current % Next
    END DO   ! while associated current
    WRITE(umMessage,*)
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')

    WRITE(umMessage,*)
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')

  END IF ! associated list begin
END DO ! loop over all lists

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Grib_Debug_Print_Basics

!=======================================================================
! **********************************************************************
! Print a count of entries for each list                ! debug
! **********************************************************************
!=======================================================================

SUBROUTINE Grib_Debug_ListCounts(Lists)

USE umPrintMgr, ONLY: &
    umPrint,           &
    umMessage
IMPLICIT NONE

!An array of pointer pairs to form the head and tail of the field lists
TYPE (List_Marker)               :: Lists(0:grib_max_fields)

INTEGER                          :: i
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GRIB_DEBUG_LISTCOUNTS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO i = 1, grib_max_fields
  IF (ASSOCIATED( Lists( i) % Begin )) THEN
    WRITE(umMessage,'(A20,A,I3,A)')                                      &
              ADJUSTR(Lists( i) % Begin % Desc) , " contains ",   &
              Lists( i) % LstCount, " members"
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
  ELSE
    WRITE(umMessage,'(A20,A,I3,A)')                                      &
              "** No Description **" , " contains ",              &
              Lists( i) % LstCount, " members"
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
  END IF
END DO

WRITE(umMessage,'(A,I3,A)') "     grib_Misc_field contains ",            &
            Lists(grib_misc_field)% LstCount,  " members"
CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Grib_Debug_ListCounts

!=======================================================================
! **********************************************************************
!  Print contents of the blocks for a Record passsed in      ! debug   *
! **********************************************************************
!=======================================================================

SUBROUTINE Grib_Debug_Print_Blocks(Record)

USE umPrintMgr, ONLY: &
    umPrint,           &
    umMessage
IMPLICIT NONE

TYPE (Grib_Record), POINTER      :: Record

INTEGER                          :: i
CHARACTER(LEN=50)                :: cFormat1,cFormat2
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GRIB_DEBUG_PRINT_BLOCKS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

cFormat1 = "(A,I2,A,I10)"
cFormat2 = "(A)"

! Loop through the blocks of the entry passed in.
IF (ASSOCIATED(Record)) THEN
  ! Block_0
  WRITE(umMessage,cFormat2) "** Block_0 **"
  CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
  DO i = 1,Len_Block_0
    WRITE(umMessage,cFormat1) "Octet ",i," reads ", Record % Block_0(i)
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
  END DO

  ! Block_1
  WRITE(umMessage,cFormat2) "** Block_1 **  -- first 50 elements only"
  CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
  DO i = 1,50    ! The rest is _normally_ blank
    WRITE(umMessage,cFormat1) "Octet ",i," reads ", Record % Block_1(i)
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
  END DO

  ! Block_2
  WRITE(umMessage,cFormat2) "** Block_2 **"
  CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
  DO i = 1,Len_Block_2
    WRITE(umMessage,cFormat1) "Octet ",i," reads ", Record % Block_2(i)
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
  END DO

  ! Block_3
  WRITE(umMessage,cFormat2) "** Block_3 **"
  CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
  DO i = 1,Len_Block_3
    WRITE(umMessage,cFormat1) "Octet ",i," reads ", Record % Block_3(i)
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
  END DO

  ! Block_4
  WRITE(umMessage,cFormat2) "** Block_4 **"
  CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
  DO i = 1,Len_Block_4
    WRITE(umMessage,cFormat1) "Octet ",i," reads ", Record % Block_4(i)
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
  END DO

  ! Block_R
  WRITE(umMessage,cFormat2) "** Block_R **"
  CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
  DO i = 1,Len_Block_R
    WRITE(umMessage,*) "Octet ",i," reads ", Record % Block_R(i)
    CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
  END DO

ELSE  ! record was null on calling

  WRITE(umMessage,cFormat2) "rcf_Grib_Debug_Print_Blocks was passed an " // &
              "unassociated record"
  CALL umPrint(umMessage,src='rcf_grib_debug_tools_mod')
END IF     ! Associated (Record)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Grib_Debug_Print_Blocks

END MODULE Rcf_Grib_Debug_Tools_Mod
