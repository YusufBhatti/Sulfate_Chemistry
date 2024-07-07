! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!-------------------------------------------------------------------------
! integer part of pp header

MODULE pp_int_head32_mod

! Description:
! Integer pp header for 32 bit fields.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utilty - crmstyle_coarse_grid
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

USE word_sizes_mod, ONLY: FourByteInt

IMPLICIT NONE

TYPE pp_int_header32

  ! integer part of pp header

    ! Validity time or start of averaging period
  INTEGER(KIND=FourByteInt)  :: yr     ! year
  INTEGER(KIND=FourByteInt)  :: mon    ! month
  INTEGER(KIND=FourByteInt)  :: dat    ! day of month
  INTEGER(KIND=FourByteInt)  :: hr     ! hour
  INTEGER(KIND=FourByteInt)  :: minute ! minute
  INTEGER(KIND=FourByteInt)  :: day    ! day number
  ! Data time or end averaging period
  INTEGER(KIND=FourByteInt)  :: yrd    ! year
  INTEGER(KIND=FourByteInt)  :: mond   ! month
  INTEGER(KIND=FourByteInt)  :: datd   ! day of month
  INTEGER(KIND=FourByteInt)  :: hrd    ! hour
  INTEGER(KIND=FourByteInt)  :: mind   ! minute
  INTEGER(KIND=FourByteInt)  :: dayd   ! day number

  INTEGER(KIND=FourByteInt)  :: tim    ! Time indicator
  INTEGER(KIND=FourByteInt)  :: ft     ! Forecast period
  INTEGER(KIND=FourByteInt)  :: lrec   ! length of record in words including
                                       ! any extra data
  INTEGER(KIND=FourByteInt)  :: code   ! grid code
  INTEGER(KIND=FourByteInt)  :: hem    ! hemisphere indicator
  INTEGER(KIND=FourByteInt)  :: row    ! number of rows
  INTEGER(KIND=FourByteInt)  :: npt    ! number of points per row (columns)
  INTEGER(KIND=FourByteInt)  :: ext    ! length of extra data
  INTEGER(KIND=FourByteInt)  :: ipack  ! packing method indicator
  INTEGER(KIND=FourByteInt)  :: rel    ! Header release number (=2)
  INTEGER(KIND=FourByteInt)  :: fc     ! pp field code
  INTEGER(KIND=FourByteInt)  :: cfc    ! 2nd pp field code
  INTEGER(KIND=FourByteInt)  :: proc   ! processing code
  INTEGER(KIND=FourByteInt)  :: vc     ! vertical coordinate type
  INTEGER(KIND=FourByteInt)  :: rvc    ! vertical coord type for reference level
  INTEGER(KIND=FourByteInt)  :: expn   ! Experiment number
  INTEGER(KIND=FourByteInt)  :: begin  ! Direct access field only start address
  INTEGER(KIND=FourByteInt)  :: nrec   ! Direct access field only no. of records
  INTEGER(KIND=FourByteInt)  :: Met08_proj   ! projection number
  INTEGER(KIND=FourByteInt)  :: met08_typ    !  Met O 8 field code
  INTEGER(KIND=FourByteInt)  :: met08_lev    !  Met O 8 level code
  INTEGER(KIND=FourByteInt)  :: unset(4)    ! Future use
  INTEGER(KIND=FourByteInt)  :: srce    ! If 1111 then user items set
  INTEGER(KIND=FourByteInt)  :: usera(3) !

  INTEGER(KIND=FourByteInt)  :: stashcode ! UM stashcode
  INTEGER(KIND=FourByteInt)  :: userb(3) !

END TYPE pp_int_header32

END MODULE pp_int_head32_mod
!
