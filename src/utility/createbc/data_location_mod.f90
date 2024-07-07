! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE data_location_mod

IMPLICIT NONE

! Description:
!   A class to store the locations of data on disk. 
!   This information is related to the fieldsfile class, and is used to 
!   contain information on where fields are stored on disk. Each field object
!   should have a corresponding data_location object when they used within the
!   fieldsfile class.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.


TYPE, PUBLIC :: data_location_type

  INTEGER, ALLOCATABLE :: lookup_position(:)
  INTEGER, ALLOCATABLE :: data_start(:)
  INTEGER, ALLOCATABLE :: data_length(:)
  INTEGER, ALLOCATABLE :: disk_length(:)

END TYPE data_location_type

END MODULE data_location_mod
