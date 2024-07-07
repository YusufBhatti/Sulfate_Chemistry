! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Purpose random number generation routine
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model
!
!    Programming standards :
!      Fortran 95 , UMDP3
!----------------------------------------------------------------
MODULE random_num_var

IMPLICIT NONE

! generic variables used by the random generator routines.
! constitutes the 'memory' of the timeseries 

INTEGER, PARAMETER :: tabs = 32
INTEGER :: drive
INTEGER, SAVE :: iv(tabs), iy 
INTEGER :: jr                  ! index

END MODULE random_num_var
