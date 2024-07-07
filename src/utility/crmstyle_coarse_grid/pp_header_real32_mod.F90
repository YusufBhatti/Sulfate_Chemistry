! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  real part of pp header

MODULE pp_real_head32_mod

! Description:
! Real PP header for 32 bit fields.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utilty - crmstyle_coarse_grid
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

USE word_sizes_mod, ONLY: FourByteReal

IMPLICIT NONE

TYPE pp_real_header32

  ! real part of header  46-64

  REAL(KIND=FourByteReal)  :: rsvd(4)  ! Unused
  REAL(KIND=FourByteReal)  :: datum    ! constant value subtracted from field
  REAL(KIND=FourByteReal)  :: acc      ! packing accuracy
  REAL(KIND=FourByteReal)  :: lev      ! level
  REAL(KIND=FourByteReal)  :: rlev     ! reference level
  REAL(KIND=FourByteReal)  :: hlev     ! hybrid levels a(k)
  REAL(KIND=FourByteReal)  :: hrlev    ! hybrid reference level
  REAL(KIND=FourByteReal)  :: plat     ! REAL*4 latitude of pseudo N pole
  REAL(KIND=FourByteReal)  :: plon     ! REAL*4 longitude of pseudo N pole
  REAL(KIND=FourByteReal)  :: gor      ! Not used lat/long grids
  REAL(KIND=FourByteReal)  :: zy       ! Zeroth latitude
  REAL(KIND=FourByteReal)  :: dy       ! latitude interval
  REAL(KIND=FourByteReal)  :: zx       ! Zeroth longitude
  REAL(KIND=FourByteReal)  :: dx       ! longitude interval
  REAL(KIND=FourByteReal)  :: mdi      ! missing data value
  REAL(KIND=FourByteReal)  :: mks      ! scaling factor


END TYPE pp_real_header32

END MODULE pp_real_head32_mod

