! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Data module defineing data held by trans namelist

MODULE Rcf_trans_Mod

! Description:
!   Defines arrays to hold data from the trans namelist
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


IMPLICIT NONE

INTEGER, ALLOCATABLE      :: ItemC_array(:)   ! item code array
INTEGER, ALLOCATABLE      :: SctnC_array(:)   ! section code array
INTEGER, ALLOCATABLE      :: Lev1_array(:)    ! bottom level
INTEGER, ALLOCATABLE      :: Lev2_array(:)    ! top level
INTEGER, ALLOCATABLE      :: Col1_array(:)    ! 1st column
INTEGER, ALLOCATABLE      :: Col2_array(:)    ! 2nd column
INTEGER, ALLOCATABLE      :: Row1_array(:)    ! 1st row
INTEGER, ALLOCATABLE      :: Row2_array(:)    ! 2nd row

INTEGER                   :: Num_Trans     ! number of transplant items

END MODULE Rcf_trans_Mod
