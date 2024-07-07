! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Subroutine ABORT_IO-----------------------------------------------
!
!    Purpose:  Prints out message and stops execution of program.
!              Called if ICODE  /=  0
!
!    Code Owner: Please refer to the UM file CodeOwners.txt
!    This file belongs in section: Dump I/O

SUBROUTINE abort_io(string,cmessage,icode,nft)


USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE UM_ParCore, ONLY: mype, nproc
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER ::                                                        &
 icode                                                            &
         !IN Code returned by UM routines
,nft     !IN Unit no being processed

CHARACTER(LEN=*) ::                                               &
 string           
        !IN Subroutine name and position

CHARACTER(LEN=errormessagelength) ::                              &
 cmessage!IN Message returned by UM routines

INTEGER, PARAMETER :: stdout = 0  ! Standard output unit number.
INTEGER, PARAMETER :: stderr = 6  ! Standard error unit number.
!----------------------------------------------------------------------

WRITE(stdout,*) 'Processor ',mype,' calling ABORT'
WRITE(stderr,*) 'Processor ',mype,' calling ABORT'

WRITE(stdout,'('' Error detected in subroutine '',A)') string
WRITE(stderr,'('' Error detected in subroutine '',A)') string

IF (nft /= 0) THEN
  WRITE(stdout,'('' while doing i/o on unit'',I3)') nft
  WRITE(stderr,'('' while doing i/o on unit'',I3)') nft
END IF
WRITE(stdout,'(A)') TRIM(cmessage)
WRITE(stderr,'(A)') TRIM(cmessage)

WRITE(stdout,'('' ICODE='',I6)') icode

CALL gc_abort(mype,nproc,cmessage)

END SUBROUTINE abort_io

