! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Small execs
SUBROUTINE READ_REC_ffread1a(field,num_cray_words,iwa,ppunit,     &
                    icode,cmessage)
USE io
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE
CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER ::                                                        &
     icode                                                        &
                            !OUT return code
    ,num_cray_words                                               &
                            !IN  No of CRAY words holding the data
    ,ppunit                                                       &
                            !IN  FT no of the PP FILE
    ,iwa                    !IN  WORD address of field to be read
REAL ::                                                           &
     field(num_cray_words)  !OUT array holding data
!    LOCAL VARIABLES
INTEGER ::                                                        &
     i                                                            &
                            ! local counter
    ,j                                                            &
                            ! local counter
    ,ix                                                           &
                            ! used in the UNIT command
    ,len_io                 ! used for call to BUFFIN
REAL ::                                                           &
     a_io                   ! used for call to BUFFIN
CALL setpos(ppunit,iwa,icode) ! C coded routine
CALL buffin(ppunit,field,num_cray_words,len_io,a_io)

IF (icode  /=  0) CALL ereport("READ_REC", icode, cmessage)
RETURN
END SUBROUTINE READ_REC_ffread1a
