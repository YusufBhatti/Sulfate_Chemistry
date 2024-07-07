! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE CHK_LOOK---------------------------------------
!
!
!    Purpose: Cross checks pointers in PP LOOKUP records with
!             model parameters
!
!    Programming standard: Unified Model Documentation Paper No 3
!
!    Documentation: Unified Model Documentation Paper No F3
!
!    Code Owner: Please refer to the UM file CodeOwners.txt
!    This file belongs in section: Dump I/O

SUBROUTINE chk_look(fixhd,lookup,len1_lookup,                     &
                                                  ! Intent (In)
                    len_data,                                     &
                                                  !
                    icode,cmessage)               ! Intent (Out)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE lookup_addresses

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE errormessagelength_mod, ONLY: errormessagelength
USE pr_look_mod, ONLY: pr_look
USE application_description, ONLY: exe_type, exe_convpp

IMPLICIT NONE


INTEGER ::                                                        &
 len_data                                                         &
                 !IN Length of model data
,len1_lookup                                                      &
                 !IN First dimension for LOOKUP
,lookup(:,:)                                            &
                 !IN Integer equivalence of PP LOOKUP
,fixhd(*)                                                         &
                 !IN Fixed length header
,icode          !OUT Return code; successful=0
                !                 error > 0

CHARACTER(LEN=errormessagelength) ::                              &
 cmessage       !OUT Error message if ICODE > 0

! -------------------------------------------------------------
! Workspace usage:---------------------------------------------
! None
!--------------------------------------------------------------
! Local variables:---------------------------------------------
INTEGER ::                                                        &
 k                                                                &
         ! Loop count
,len_d                                                            &
         ! Cumulative length of data
,n1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHK_LOOK'
!--------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
icode=0
cmessage=' '

! CHK_LOOK falls over with Boundary Datasets if pre-3.4
IF (fixhd(5) == 5 .AND. fixhd(12) <  304) THEN
  WRITE(umMessage,*) ' CHK_LOOK skipped for Boundary Dataset (Pre 3.4)'
  CALL umPrint(umMessage,src='chk_look')
ELSE

  len_d=0

  DO k=1,fixhd(152)

    IF ((exe_type == exe_convpp) .AND. (lookup(1,k) == -99)) EXIT

    ! Check that data_type is valid no: 1 to 3 or -1 to -3
    IF ((lookup(data_type,k) >= 1 .AND. lookup(data_type,k) <= 3)   &
   .OR. (lookup(data_type,k) <= -1 .AND. lookup(data_type,k) >= -3))&
    THEN
      len_d=len_d+lookup(lblrec,k)
    ELSE
      WRITE(umMessage,'(A,i9,A,i4,A)')' *ERROR* Wrong value of',            &
      lookup(data_type,k),' IN LOOKUP(DATA_TYPE',k,')'
      CALL umPrint(umMessage,src='chk_look')

      CALL pr_look(lookup,k)

      icode=1
      cmessage='CHK_LOOK: Consistency check'
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF

  END DO

END IF  !  Check on pre-3.4 Boundary Datasets
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE chk_look

