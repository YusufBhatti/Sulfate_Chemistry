! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Interface and arguments: ------------------------------------------
!  Routine: LOGICAL_TO_REAL ------------------------------------------
!
!  Purpose: To convert logical data within FIELD to real data.
!  the data in FIELD.
!
!  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
!  Logical components covered: ...
!
!  Project task: ...
!
!  External documentation:
!
!  -------------------------------------------------------------------
!  Interface and arguments: ------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
SUBROUTINE LOGICAL_TO_REAL_ffread1a(IDIM,logical_field,field,     &
                           nvals,ilabel,icode,cmessage)

USE errormessagelength_mod, ONLY: errormessagelength
USE lookup_addresses
IMPLICIT NONE
INTEGER ::                                                        &
     IDIM                                                         &
                          !IN  The full unpacked size of a field
    ,ilabel(45)                                                   &
                          !OUT holds integer part of LOOKUP
    ,icode                !OUT Non zero for any error
REAL ::                                                           &
     field(IDIM)          !OUT On Input contains Real data.
LOGICAL ::                                                        &
     logical_field(IDIM)  !INOUT On Input contains logical data.
!                               ! contains the un-packed data.
CHARACTER(LEN=errormessagelength) :: cmessage 
                                !OUT Will contain any error messages.
!
!     LOCAL  VARIABLES
INTEGER ::                                                        &
     i                                                            &
                            ! Loop counter
    ,nvals                  ! IN no of values in an input field
!
!
DO  i=1,nvals
  IF (logical_field(i)) THEN
    field(i)=1.0
  ELSE
    field(i)=0.0
  END IF
END DO
ilabel(data_type)=1     ! The data type must now be real
icode=0
RETURN
END SUBROUTINE LOGICAL_TO_REAL_ffread1a
