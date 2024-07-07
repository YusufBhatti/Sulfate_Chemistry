! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! 
MODULE cancila_mod

! Description:
!   Specifies control parameters and arrays for ancillary updating
!
! Code Description:
!  Language: FORTRAN 95.
!  This code is written to UMDP3 programming standards.
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Ancillaries
IMPLICIT NONE

INTEGER, PARAMETER :: nancil_fields = 201 ! maximum number of ancillary fields

! Position of ancillary field in lookup tables.
INTEGER :: nlookup     (nancil_fields) 

! Interval between PP Headers referring to the same 
!  ancillary fields at different times
INTEGER :: lookup_step (nancil_fields) 

! Number of levels of data in each ancillary field
INTEGER :: levels      (nancil_fields) 

! Address of ancillary field in main data block
INTEGER :: d1_anciladd (nancil_fields)    

! Updating intervals
INTEGER :: steps(nancil_fields)

! Updating switches in replanca
LOGICAL :: update(nancil_fields)

END MODULE cancila_mod
