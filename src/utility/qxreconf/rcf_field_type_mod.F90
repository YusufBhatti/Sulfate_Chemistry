! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Defines the field_type data-type

MODULE Rcf_Field_Type_Mod

! Description:
!   Data module defining the field_type data-type. This should contain
!   all information commonly required of a field such as sizes,
!   stashmaster, dump position, etc
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE Rcf_Ppx_Info_Mod, ONLY: &
    STM_record_type

IMPLICIT NONE

TYPE field_type

  ! data field - two dimensions x*y and z
  ! For some situations (interpolation of integer field for e.g.),
  ! more than one of these fields may be used at a given time.
  REAL, ALLOCATABLE :: data( :, : )
  INTEGER, ALLOCATABLE :: data_int( :, : )
  LOGICAL, ALLOCATABLE :: data_log( :, : )

  ! Sizes
  INTEGER          :: levels
  ! Range according to STASHmaster (do we need to distinguish pseudo-levels?)
  INTEGER          :: bottom_level
  INTEGER          :: top_level

  ! Local sizes are default
  INTEGER          :: rows
  INTEGER          :: row_len
  INTEGER          :: level_size           ! Size of single level

  ! Global sizes are also required
  INTEGER          :: glob_rows
  INTEGER          :: glob_row_len
  INTEGER          :: glob_level_size

  ! Information about the sort of field it is
  INTEGER          :: dump_pos             ! start position in dump
  INTEGER          :: interp               ! Should it be/has it been
                                           ! interpolated/copied
  ! Pointer to relevant stashmaster field
  TYPE (STM_record_type), POINTER :: stashmaster => NULL() ! stashmaster

END TYPE field_type

END MODULE Rcf_Field_Type_Mod

