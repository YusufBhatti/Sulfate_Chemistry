! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Subroutine: LOGICAL_TO_REAL ------------------------------------------
!
!    Purpose: To convert logical data within FIELD to real data.
!    the data in FIELD.
!
!    Tested under compiler:   cft77
!    Tested under OS version: UNICOS 5.1
!
!    Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
!    Logical components covered:
!
!    Project task:
!
!    External documentation:
!
!    -------------------------------------------------------------------
!
!    Code Owner: Please refer to the UM file CodeOwners.txt
!    This file belongs in section: Small execs
SUBROUTINE LOGICAL_TO_REAL_convpp(npts,logical_field,field,     &
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
    ,ilabel(44)                                                   &
                          !OUT integer part of LOOKUP
    ,icode                !OUT error code
REAL ::                                                           &
     field(npts)          !OUT contains Real data.
LOGICAL ::                                                        &
     logical_field(npts)  !IN contains logical data.
!                               ! contains the un-packed data.
!     Local variables
INTEGER ::                                                        &
     i                    ! loop counter
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
END SUBROUTINE LOGICAL_TO_REAL_convpp
