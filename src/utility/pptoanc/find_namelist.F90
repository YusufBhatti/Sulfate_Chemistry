! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

FUNCTION find_namelist(iunit, namelist_name)

USE umPrintMgr, ONLY:      &
          umPrint,           &
          umMessage

IMPLICIT NONE
!
! Description:
!  This routine searches the input stream given by 'iunit'
!  to find the NAMELIST given by 'namelist_name'.  The
!  input file is then correctly positioned to let F90
!  library routines read the namelist.
!
!  The namelist is assumed to be contained within the
!  first 24 characters of the record
!
!  Return values:
!
!  -1  error - could not find the namelist - file is now
!      at end-of-file.
!   0  namelist ready to be processed.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
!
! Declarations:
! 1.0 Subroutine arguments
!   1.1 Scalar arguments with intent(in):
INTEGER :: iunit

!   1.2 Scalar arguments with intent(out):
CHARACTER(LEN=*) :: namelist_name

INTEGER :: find_namelist

! 2.0 Local scalars:
INTEGER :: i, j

CHARACTER(LEN=24) :: chvar

!- End of header

1000  CONTINUE
READ(iunit, '(a)', END=9000, err=9000) chvar

!  Check for leading '&' for namelist
i=INDEX(chvar, '&')
!  Found '&' - check for the name we want
IF (i /= 0) THEN
  j=INDEX(chvar, namelist_name)
  !  Not the name we want - print skipped message
  IF (j == 0) THEN
    IF (INDEX('endEndeNdenDENdEnDeNDEND',                                     &
        chvar(i+1:i+3)) == 0) THEN
      WRITE(umMessage,*)'- Skipped record named: ',                          &
          chvar(i+1:),' On Unit:',iunit,                                     &
          ' - f90 version'
      CALL umPrint(umMessage,src='find_namelist')
    END IF
    GO TO 1000
  END IF
  GO TO 1100
END IF
GO TO 1000

!  Found the namelist we want - backspace to position correctly
1100  CONTINUE
BACKSPACE iunit
find_namelist=0
RETURN

!  Cannot find the namelist we want
9000  CONTINUE
find_namelist=1
RETURN
END FUNCTION find_namelist
