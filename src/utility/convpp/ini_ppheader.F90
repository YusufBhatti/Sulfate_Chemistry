! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE ini_ppheader_mod

IMPLICIT NONE

LOGICAL :: ibm

CONTAINS

! Deals with namelist input to set up initialisation of pp header manipluations
SUBROUTINE ini_ppheader(ichan)

! Description:
!   Defines and reas a namelist containing switches that determine the pp
!   header manipluations to be performed.

! Method:
!   Standard namelist defines and reads.

! See also:
!   pp_header_manips

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Small execs

! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.

USE pp_header_manips, ONLY: set_zero

IMPLICIT NONE

INTEGER, INTENT(IN)   :: ichan ! opened fortran unit to read namelist from

LOGICAL :: zero_diraccess

NAMELIST /ppheader/ ibm, zero_diraccess

ibm=.FALSE.
zero_diraccess=.FALSE.
READ(ichan, NML=ppheader)

CALL set_zero(zero_diraccess)

END SUBROUTINE ini_ppheader

END MODULE ini_ppheader_mod
