! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Assign 'GRIB record' to the correct linked list

MODULE Rcf_Grib_Assign_Mod

! SUBROUTINE Rcf_Grib_Assign  - File 'current' data into relevant list
!
! Description: This routine files the information in derived type
!              'Current' into one of the predefined 'lists' used for
!              each parameter. (as specified in rcf_grib_lookups.F90)
!
! Method:  Using the center ID no. (Block_1, octet 1) compare the
!          parameter ID of the data within the lookup table.
!          Use the lookup table to find the stash code with which to
!          file the data in the correct 'list'
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GRIB_ASSIGN_MOD'

CONTAINS
SUBROUTINE Rcf_Grib_Assign(Current, Lists)

USE Rcf_GRIB_Lookups_Mod               ! Contains cross ref table for
                                       ! Stash => ECMWF parameter Id's
                                       ! and appropriate params

USE Rcf_GRIB_Block_Params_Mod          ! provides type def for Current

USE EReport_Mod, ONLY:     &
    EReport

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Min,           &       ! =1 Minimum output
    PrStatus_Normal,        &       ! =2 Short informative output
    PrStatus_Oper,          &       ! =3 Full informative output
    PrStatus_Diag                   ! =4 Extra Diagnostic output

USE UM_ParCore, ONLY: &
    mype

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
TYPE (Grib_Record),POINTER       :: Current       !\ Pointer to current
                                                  !/ grib record
TYPE (List_Marker),INTENT(INOUT) :: Lists(0:grib_max_fields)

! Local variables

CHARACTER (LEN=*), PARAMETER     :: RoutineName='RCF_GRIB_ASSIGN'
CHARACTER (LEN=errormessagelength)   :: Cmessage  ! used for EReport
INTEGER                          :: ErrorStatus   ! used for EReport
INTEGER                          :: i             ! Loop counter
INTEGER                          :: LookupCol     ! Column in 'Table A'
                                                  ! of rcf_grib_lookups
INTEGER                          :: ListNo        ! List to store data

LOGICAL                          :: L_Error       ! Error Flag

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!=======================================================================
!  Main Routine
!=======================================================================

!=======================================================================
!  Test record in order to assign to correct list
!=======================================================================

! First - Check which originating center data came from
!         And set the lookup column accordingingly
SELECT CASE (Current % Block_1(p_Orig_cntr))

CASE (GrbOrigECMWF)   ! Grib came from ECMWF
  LookupCol = p_ECMWF_IDCol

CASE (GrbOrigUKMO)   ! Grib came from UKMO
  LookupCol = p_UKMO_IDCol

CASE DEFAULT          ! Unidentified center code
  WRITE (Cmessage,'(A,I3,1X,A)')                                    &
                           'Originating Center code number :',      &
                           Current % Block_1(p_Orig_cntr)           &
                           ,'not recognised.'
  ErrorStatus = 10
  CALL EReport( RoutineName, ErrorStatus, Cmessage )

END SELECT

! Second - Use the lookup to get the Stash Code and List number
ListNo = 0                     ! The default list for 'misc' fields
Current % StashCode = -1
Current % Desc      = "Unknown Parameter   "

IF ( PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
  WRITE(umMessage,'(A,I3)') "Looking for Parameter :",                       &
               Current % Block_1(p_Param_ID)
  CALL umPrint(umMessage,src='rcf_grib_assign_mod')
END IF

L_Error = .TRUE.
FILE: DO i = 1, p_Max_Rows
  IF (Current % Block_1(p_Param_ID) ==                                &
              Param_ID_CrossRef( i) % CrossRefIDs(LookupCol) .AND.    &
      Current % Block_0(p_Tbl_Vers_No) ==                             &
              Param_ID_CrossRef( i) % Table_No(LookupCol) ) THEN

    Current % StashCode =                                             &
                  Param_ID_CrossRef(i) % CrossRefIDs(p_STASH_IDCol)
    ListNo = Param_ID_CrossRef( i) % List_No
    Current % Desc = Param_ID_CrossRef(i) % DescText
    L_Error = .FALSE.

    IF ( PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
      WRITE(umMessage,'(2A)') "Found Parameter :",                           &
                                     Param_ID_CrossRef(i) % DescText
      CALL umPrint(umMessage,src='rcf_grib_assign_mod')
    END IF
    EXIT FILE                       ! found what I want - exit loop

  END IF
END DO FILE

! double check that new record was found in lookup table.
IF ( L_Error ) THEN
  IF ( PrintStatus >= PrStatus_Normal .AND. mype == 0 ) THEN
    WRITE (Cmessage,'(A,I3,A)') 'Parameter ',                         &
                         Current % Block_1(p_Param_ID),               &
                         ' in GRIB file not found in lookup table'
    ErrorStatus = -10    ! Warning only - unknown data filed in misc
    CALL EReport( RoutineName, ErrorStatus, Cmessage )
  END IF
END IF

!=======================================================================
!  Now put field in the selected list
!=======================================================================

Current % Prev   => Lists(ListNo) % END
                               ! Point Prev pointer at end of
                               ! current list.(Null if first entry)
IF (ASSOCIATED(Lists(ListNo) % END)) THEN
                               ! If current end of list is a
                               ! valid record (Not first entry)
  Lists(ListNo) % END % Next  => Current
                               ! Point 'next' for previous entry
                               ! at current entry

ELSE                           ! Else : must be 1st entry
  Lists(ListNo) % Begin  => Current
                               ! Point begining of List at Current
END IF

Lists(ListNo) % END      => Current
                               ! Point End of List at (now complete)
                               ! Current Entry
NULLIFY(Current % Next)        ! Ensure 'Next' is not associated

Lists(ListNo) % LstCount = Lists(ListNo) % LstCount + 1
                               ! Add one to count of list size

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_Grib_Assign
END MODULE Rcf_Grib_Assign_Mod
