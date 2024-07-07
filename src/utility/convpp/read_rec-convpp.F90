! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Subroutine: READ_REC
!
!    Purpose: To read a data record from a  pp file
!
!    Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
!    Logical components covered: ...
!
!    Project task: ...
!
!    External documentation:
!
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Small execs
SUBROUTINE read_rec_convpp(field,num_cray_words,iwa,ppunit,       &
                    icode,cmessage)
USE io
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE
!     arguments
CHARACTER(LEN=errormessagelength) :: cmessage      !OUT error message
INTEGER ::                                                        &
     num_cray_words                                               &
                            !IN  No of CRAY words holding the data
    ,ppunit                                                       &
                            !IN  unit no of the PP FILE
    ,iwa                                                          &
                            !IN  WORD address of field to be read
    ,icode                  !OUT error code
REAL ::                                                           &
     field(num_cray_words)  !OUT array holding data
!     arguments for called routines
INTEGER ::                                                        &
     len_io                 ! length of data read by BUFFIN
REAL ::                                                           &
     a_io                   ! return code from BUFFIN
!    LOCAL VARIABLES
INTEGER ::                                                        &
     i                                                            &
                            ! local counter
    ,j                                                            &
                            ! local counter
    ,ix                     ! used in the UNIT command

CALL setpos(ppunit,iwa,icode)
CALL buffin(ppunit,field,num_cray_words,len_io,a_io)

RETURN
END SUBROUTINE read_rec_convpp

