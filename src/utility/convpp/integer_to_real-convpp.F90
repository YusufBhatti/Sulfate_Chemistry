! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Routine: INTEGER_TO_REAL
!
!    Purpose: To convert logical data within FIELD to real data.
!    the data in FIELD.
!
!    Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
!    Logical components covered:
!
!    Project task:
!
!    External documentation:
!
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Small execs
SUBROUTINE INTEGER_TO_REAL_convpp(npts,integer_field,field,     &
                           nvals,ilabel,icode,cmessage)
USE lookup_addresses
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!     arguments
CHARACTER ::                                                      &
     cmessage*(errormessagelength)         !OUT error mesages.
INTEGER ::                                                        &
     npts                                                         &
                          !IN full unpacked size of a field
    ,nvals                                                        &
                          !IN no of values in an input field
    ,integer_field(npts)                                          &
                          !IN contains integer data.
    ,ilabel(44)                                                   &
                          !OUT integer part of LOOKUP
    ,icode                !OUT error code
REAL ::                                                           &
     field(npts)          !OUT contains Real data.
!     Local variables
INTEGER ::                                                        &
     i                    ! loop counter
!
!

DO  i=1,nvals
  field(i)=integer_field(i)
END DO
ilabel(data_type)=1       ! The data type must now be real
icode=0
RETURN
END SUBROUTINE INTEGER_TO_REAL_convpp
