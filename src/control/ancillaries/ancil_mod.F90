! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE Ancil_Mod

!
! Description:
! Define structures to hold Ancillary related information

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Ancillaries
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
USE filenamelength_mod, ONLY:          &
    filenamelength

USE missing_data_mod, ONLY:  &
    imdi

IMPLICIT NONE

INTEGER , SAVE :: AncF_UnitNo ! Unit No for ancillary files, used in recon only

! Allocatable work arrays for ancillary processing
INTEGER, ALLOCATABLE :: nlookup(:)
INTEGER, ALLOCATABLE :: lookup_step(:)
INTEGER, ALLOCATABLE :: levels(:)
INTEGER, ALLOCATABLE :: ancil_add(:)

! -----------------------------------------------------------------------
INTEGER  :: num_ancil_requests  = 0 ! total number of ancillary
                                    ! requests can be larger than 
                                    ! num_items read in from RECONA

INTEGER  :: num_ancil_files = 0     ! number of ancillary files to open

INTEGER, PARAMETER  :: max_items = 1000     ! maximum number of ITEM requests
INTEGER, PARAMETER  :: stash_num_max=25     ! maximum requests to a single
                                            ! ancillary file

! Define a type to hold information about each ancil file
TYPE ancil_file_type
  INTEGER :: stashcodes(stash_num_max) ! stashcode
  INTEGER :: num_stash = 0             ! number of stashcodes due to be 
                                       ! requested from file
  INTEGER :: unit_num = imdi
  ! The default set here cannot be same as the default items namelist filename 
  ! value otherwise it will cause problems when searching for matches
  CHARACTER (LEN=filenamelength) :: filename = 'no filename set'
END TYPE ancil_file_type

! Define structure to hold ancillary request information
TYPE ancil_req_type
  INTEGER :: stashcode ! STASHcode
  INTEGER :: section   ! STASH section number
  INTEGER :: item      ! STASH item number
  INTEGER :: period    ! Updating period (UM only)
  INTEGER :: interval  ! Updating interval (UM only)
  CHARACTER (LEN=36) :: stash_name
  ! Store the position in the ancil file array of the ancil file associated
  ! with this ancil request
  INTEGER :: ancil_file_num 
END TYPE  ancil_req_type

! the list of fields that need to be updated
TYPE (ancil_req_type), SAVE   ::   ancil_requests(max_items)
! the list of ancillary files to be read in
TYPE (ancil_file_type), SAVE  ::   ancil_files(max_items)

CONTAINS

INTEGER FUNCTION find_ancil_req_by_stash(stashcode)
IMPLICIT NONE
INTEGER, INTENT(IN) :: stashcode
INTEGER :: i
DO i = 1, num_ancil_requests
  IF (ancil_requests(i)%stashcode == stashcode) THEN
    find_ancil_req_by_stash = i
    RETURN
  END IF
END DO
find_ancil_req_by_stash = imdi
RETURN 
END FUNCTION find_ancil_req_by_stash

INTEGER FUNCTION find_ancil_file_by_stash(stashcode)
IMPLICIT NONE
INTEGER, INTENT(IN) :: stashcode
INTEGER :: i
DO i = 1, num_ancil_requests
  IF (ancil_requests(i)%stashcode == stashcode) THEN
    find_ancil_file_by_stash = ancil_requests(i) % ancil_file_num
    RETURN
  END IF
END DO
find_ancil_file_by_stash = imdi
RETURN 
END FUNCTION find_ancil_file_by_stash

SUBROUTINE add_ancil_request(stashcode,      section,       &
                             item,           period,        &
                             interval,       anc_filename  )
! Add ancil request to both the ancil_requests and ancil_files arrays 
! contained in the ancil_mod module.
IMPLICIT NONE
INTEGER, INTENT(IN) :: stashcode
INTEGER, INTENT(IN) :: section
INTEGER, INTENT(IN) :: item
INTEGER, INTENT(IN) :: period
INTEGER, INTENT(IN) :: interval
CHARACTER(LEN=filenamelength), INTENT(IN) :: anc_filename
INTEGER             :: k 

! Add to ancil_requests type
num_ancil_requests = num_ancil_requests + 1
ancil_requests(num_ancil_requests) % stashcode = stashcode 
ancil_requests(num_ancil_requests) % section   = section 
ancil_requests(num_ancil_requests) % item      = item
ancil_requests(num_ancil_requests) % period    = period 
ancil_requests(num_ancil_requests) % interval  = interval 
! Add file to ancil_files type
IF ( .NOT. ANY( ancil_files%filename == anc_filename )) THEN
  ! If the filename has not been encountered before then add to 
  ! list of ancil files
  num_ancil_files = num_ancil_files + 1
  ancil_files(num_ancil_files) % filename = TRIM(anc_filename)
  ! Add stashcode to the list of stash to be read from the ancil file
  ancil_files(num_ancil_files) % num_stash = &
       ancil_files (num_ancil_files) % num_stash + 1
  ancil_files(num_ancil_files) % stashcodes(ancil_files(num_ancil_files) &
       % num_stash) = stashcode 
  ! Point the ancil_request to the file that contains the stashcode
  ancil_requests(num_ancil_requests) % ancil_file_num = num_ancil_files
ELSE ! Filename has been encountered before
  DO k = 1, num_ancil_files
    IF (ancil_files(k)%filename == TRIM(anc_filename) ) THEN
      ! Locate file that matches and add stash to list
      ancil_files(k) % num_stash = ancil_files(k) % num_stash + 1
      ancil_files(k) % stashcodes(ancil_files(k)%num_stash) = & 
           stashcode
      ! Point the ancil_request to the file that contains the stashcode
      ancil_requests(num_ancil_requests) % ancil_file_num = k
      EXIT
    END IF
  END DO ! End loop over files
END IF
END SUBROUTINE add_ancil_request
END MODULE Ancil_Mod

