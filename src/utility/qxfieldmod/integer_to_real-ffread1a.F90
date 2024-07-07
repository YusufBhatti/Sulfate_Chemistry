! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine: INTEGER_TO_REAL ------------------------------------------
!
! Purpose: To convert logical data within FIELD to real data.
! the data in FIELD.
!
! Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
! Logical components covered: ...
!
! Project task: ...
!
! External documentation:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
SUBROUTINE INTEGER_TO_REAL_ffread1a(IDIM,integer_field,field,     &
                           nvals,ilabel,icode,cmessage)
USE errormessagelength_mod, ONLY: errormessagelength
USE lookup_addresses
IMPLICIT NONE
INTEGER ::                                                        &
     IDIM                                                         &
                          !IN  The full unpacked size of a field
    ,ilabel(45)                                                   &
                          !OUT holds integer part of LOOKUP
    ,icode                                                        &
                          !OUT Non zero for any error
    ,integer_field(IDIM)                                          &
                          !IN  On input contains integer data.
    ,nvals                !IN no of values in an input field

REAL ::                                                           &
     field(IDIM)          !OUT On Input contains Real data.
!                               ! contains the un-packed data.
CHARACTER(LEN=errormessagelength) :: cmessage    
                                !OUT Will contain any error mesages.
!
!     LOCAL  VARIABLES
INTEGER ::                                                        &
     i                    ! Loop counter
!
!
DO  i=1,nvals
  field(i)=integer_field(i)
END DO
ilabel(data_type)=1     ! The data type must now be real
icode=0
RETURN
END SUBROUTINE INTEGER_TO_REAL_ffread1a
