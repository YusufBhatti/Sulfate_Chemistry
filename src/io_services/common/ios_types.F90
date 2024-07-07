! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: C96

! Derived datatypes used by IOS

MODULE IOS_types
USE UM_types
USE IOS_Constants
IMPLICIT NONE
!-------------------------------------------------------
! IO Server metadata type definition
!-------------------------------------------------------

TYPE IOS_metadata_type
  SEQUENCE
  ! Sequence and order of definition are important in comms.
  ! Action
  INTEGER                       :: ACTION
  ! Common to all
  INTEGER                       :: UNIT
  ! Open/Close
  INTEGER                       :: name_length
  INTEGER                       :: subtype
  INTEGER                       :: delete
  ! Read/Write
  INTEGER                       :: data_size
  ! Setpos
  INTEGER                       :: address
  INTEGER                       :: client
  INTEGER                       :: handle
  INTEGER                       :: Originating_Slot

  ! Multi-purpose String
  CHARACTER (LEN=IOS_string_max):: string
END TYPE IOS_metadata_type

! Elements of the request member of IOS_async_object
INTEGER, PARAMETER                  :: cl_as_requests     = 3
INTEGER, PARAMETER                  :: cl_as_request_md   = 1
INTEGER, PARAMETER                  :: cl_as_request_str  = 2
INTEGER, PARAMETER                  :: cl_as_request_payl = 3

TYPE IOS_async_object
  SEQUENCE
  INTEGER                           :: seq_id
  INTEGER                           :: request(cl_as_requests)
  INTEGER                           :: state
  INTEGER                           :: pe
  TYPE(IOS_metadata_type)           :: md
  INTEGER(KIND=IOS_tukind), POINTER :: payload(:)
END TYPE IOS_async_object

TYPE IOS_status
  SEQUENCE
  INTEGER                       :: UNIT
  INTEGER                       :: Owner   ! The rank that opened the file
  INTEGER                       :: POSITION! The byte position in the file
  INTEGER                       :: extent  ! The byte length of the file
  LOGICAL                       :: isOpen
  LOGICAL                       :: isReadOnly
  CHARACTER (LEN=IOS_string_max):: filename
END TYPE IOS_status
INTEGER, PARAMETER              :: IOS_state_size=6

END MODULE IOS_types
