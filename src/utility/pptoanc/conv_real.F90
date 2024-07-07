! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
! Subroutine interface:
!
! Subroutine interface:
!
! Subroutine interface:
!
! Subroutine interface:
!
! Subroutine interface:
SUBROUTINE conv_real(rlookup,lookup_all,len2_lookup)

IMPLICIT NONE
!
! Description:
!             Convert's the real part of the lookup header (rlookup)
!             into integer's so that it can be represented as one
!             array(lookup_all)
!
! Method:
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
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!

! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER :: len2_lookup                !IN no. of fields

!   Array  arguments with intent(in):
REAL :: rlookup(19,len2_lookup)       !IN real part of lookup table

! Array  arguments with intent(out):
REAL :: lookup_all(64,len2_lookup)    !OUT whole lookup table

! Local scalar
INTEGER :: i                          ! loop counter

!- End of header

DO i = 1,len2_lookup
  lookup_all(46:64,i) = rlookup(1:19,i)
END DO

RETURN
END SUBROUTINE conv_real
!
! Subroutine interface:
! Purpose: Works out the lookup tables for the dump/ancillary    *
!           file header from the pp fields                       *
!
! Subroutine interface:

!  Skip namelists in f90 compiled UM code removing need for
!  assign -f 77 g:sf  in script
!
! Subroutine Interface:
