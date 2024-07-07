! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  32 bit pp field

MODULE pp_field32_mod

USE pp_int_head32_mod
USE pp_real_head32_mod
USE word_sizes_mod, ONLY: FourByteReal



IMPLICIT NONE

! Description:
! PP header for 32 bit fields.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utilty - crmstyle_coarse_grid
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

INTEGER, PARAMETER :: max_data_size = 3000*3000

TYPE pp_field32

  TYPE(pp_int_header32)   :: ipphead

  TYPE(pp_real_header32)  :: rpphead

  REAL(KIND=FourByteReal) :: DATA(max_data_size)        !

END TYPE pp_field32

END MODULE pp_field32_mod
