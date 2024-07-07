! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module dealing with some common/simple manipulations of the pp header

! Description:
!   pp files can have different pp headers to the corresponding lookup
!   elements in a fields file.  The module provides routines for converting
!   the field file lookup table elements into pp headers.

!   See the documentation for the individual subroutines for more details

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs

! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.

MODULE pp_header_manips

USE lookup_addresses

USE umPrintMgr, ONLY:      &
    umPrint,               &
    umMessage,             &
    newline
IMPLICIT NONE


INTEGER, PARAMETER :: max_len_ilabel=45

! Module Common variables
LOGICAL  :: lzero   = .FALSE.
PRIVATE  :: lzero
CONTAINS

LOGICAL FUNCTION is_gregorian(ilabel)
IMPLICIT NONE


INTEGER, INTENT(IN) ::  ilabel(max_len_ilabel)
is_gregorian = (MOD(ilabel(lbtim), 10) == 1)
RETURN
END FUNCTION is_gregorian

LOGICAL FUNCTION valid_start_year(ilabel)
IMPLICIT NONE


INTEGER, INTENT(IN) ::  ilabel(max_len_ilabel)
valid_start_year = (ilabel(lbyr) >  0)
RETURN
END FUNCTION valid_start_year

SUBROUTINE set_zero(zero)
IMPLICIT NONE


LOGICAL, INTENT(IN)    :: zero
lzero = zero

END SUBROUTINE set_zero

SUBROUTINE header_manip(ilabel)  ! change names

IMPLICIT NONE


! Description:
!   coordinates the calling of any pp-header changes that
!   need to be made
INTEGER :: ilabel(max_len_ilabel)

IF (lzero) THEN
  CALL zero_da_fields(ilabel)
END IF

END SUBROUTINE header_manip


SUBROUTINE zero_da_fields(ilabel)

IMPLICIT NONE

! Description:
! Some fields in the pp headers are only applicable to direct access data sets.
! This routine sets those fields to zero.

INTEGER, INTENT(INOUT)  :: ilabel(max_len_ilabel)

ilabel(lbnrec) = 0
ilabel(lbegin) = 0

END SUBROUTINE zero_da_fields

END MODULE pp_header_manips
