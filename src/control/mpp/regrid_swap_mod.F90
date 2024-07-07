! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Module: REGRID_SWAP
!
!  Purpose: Provides grid field data swapping routines to facilite
!  parallel regridding for both area average and map max method

!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: MPP

MODULE regrid_swap

USE regrid_utils, ONLY: atmos_grid, riv_grid,                         &
     global_to_local_gridpt, error_check_mpl

IMPLICIT NONE

INTEGER, PARAMETER :: swap_atmos_grid = 1 , swap_riv_grid = 2

CONTAINS

  !! sends and receives the required grid data to/from each
  !! process. Concerns passed to this routine must contain send
  !! receive grid points required as well as proc numbers
  !! concerned
  !!
  !! For both sendConcerns and recvConcerns arrays all space required
  !! by array element members
SUBROUTINE swap_regrid_data(send_concern, recv_concern, field,        &
     src_grid, row_length, rows, error)

USE regrid_types
USE mpl
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

TYPE(concern), INTENT(IN)    :: recv_concern(:)
TYPE(concern), INTENT(INOUT) :: send_concern(:)
! send concern contains data that will be sent to other processes
! to regrid their subdomain field, and recv_concern contain data
! to be set from data from other processes needed by this process
! to regrid its subdomain src field

INTEGER, INTENT(IN) :: row_length, rows
! row_length, rows of src grid on this process

REAL, INTENT(IN) :: field(row_length, rows)
! this process local field subdomain

INTEGER, INTENT(IN) ::  src_grid
! grid regridding from

INTEGER, INTENT(OUT) :: error
! set to -1 if an error occurs, otherwise 0

! local variables

INTEGER :: i, j, send_field_size, x, y, my_comm

INTEGER :: s_concern_size, r_concern_size, request_size,               &
     send_data_space, recv_data_space
! extent of send and recev concerns passed

INTEGER, ALLOCATABLE :: requests(:), statuses(:,:), ierror(:)
INTEGER :: error1, error2
CHARACTER(LEN=errormessagelength) :: cmessage

error = 0
s_concern_size = SIZE(send_concern,1)
r_concern_size = SIZE(recv_concern,1)

CALL gc_get_communicator(my_comm, error)

! prepare data to send

!! populate send concerns with data required by other processes
DO i=1, s_concern_size
  send_field_size = send_concern(i)%SIZE
  DO j=1, send_field_size

    x = send_concern(i)%x(j)
    y = send_concern(i)%y(j)

    CALL global_to_local_gridpt(x, y, src_grid)

    send_concern(i)%field(j) = field(x,y)
  END DO
END DO

! swap data
request_size = s_concern_size + r_concern_size

!! allocate requests for sends and recvs
ALLOCATE(requests(request_size))
ALLOCATE(ierror(request_size))

requests = 0
ierror = 0
error1 = 0
error2 = 0

!! send field point required by other procs
DO i=1, s_concern_size

  CALL mpl_isend(send_concern(i)%field, send_concern(i)%SIZE,   &
       mpl_real, send_concern(i)%proc_num, 0,  my_comm,         &
       requests(i), ierror(i))

END DO

!! recv field points you need
DO i=1, r_concern_size

  CALL mpl_irecv(recv_concern(i)%field, recv_concern(i)%SIZE,   &
       mpl_real, recv_concern(i)%proc_num, 0, my_comm,          &
       requests(i+s_concern_size), ierror(i+s_concern_size))

END DO

ALLOCATE(statuses(mpl_status_size, request_size))

statuses = 0

! check mpl errors of recv and send
CALL error_check_mpl(ierror, request_size, error1)

CALL mpl_waitall(request_size, requests, statuses, i)

! check mpl errors of waitall

! request_size = 0 for 1x1 decomposition
IF (request_size /= 0) THEN
  ierror(1) = i
  CALL error_check_mpl(ierror, 1, error2)
END IF

IF (error1 < 0 .OR. error2 < 0) THEN
  error = -1
  cmessage = "REGRID_SWAP_DATA: MPL error"
END IF

! release mem resources
DEALLOCATE(requests)
DEALLOCATE(statuses)
DEALLOCATE(ierror)

RETURN
END SUBROUTINE swap_regrid_data

!! This retrieves and sends target data to and from other
!! other processes needed for backwards regridding
!!
SUBROUTINE swap_regrid_data_max(src_field, l_src_row_length,          &
     l_src_rows, src_grid, send_concern_max, recv_concern_max, error, &
     cmessage)

USE regrid_types
USE mpl
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

TYPE(concern_max), INTENT(IN)    :: recv_concern_max(:)
TYPE(concern_max), INTENT(INOUT) :: send_concern_max(:)
! send contains information other processes need to regrid their
! subdomain and recv concern max will be set with field data
! needed by this process to regrid its subdomain
INTEGER, INTENT(IN) :: l_src_row_length, l_src_rows, src_grid
! src field dimension extents
REAL, INTENT(IN) :: src_field(l_src_row_length, l_src_rows)
! src field to be regridded to target grid
INTEGER, INTENT(OUT) :: error
! if error occurs this is set to -1, otherwise 0
CHARACTER(LEN=errormessagelength), INTENT(OUT) :: cmessage
! is set to value indicating cause of error if error = -1

! the purpose of this routine is to recv grid points this
! processor needs to send

! local variables

INTEGER :: i, j, error1, error2, error_i(1), my_comm
INTEGER :: s_concern_size, r_concern_size, request_size
INTEGER :: send_field_size, x, y
! extent of send and recev concerns passed

INTEGER, ALLOCATABLE :: requests(:), statuses(:,:), ierror(:)

s_concern_size = SIZE(send_concern_max,1)
r_concern_size = SIZE(recv_concern_max,1)
request_size = s_concern_size + r_concern_size

! allocate requests for sends and recvs
ALLOCATE(requests(request_size))
ALLOCATE(ierror(request_size))

requests = 0
ierror = 0
error = 0
error1 = 0
error2 = 0

CALL gc_get_communicator(my_comm, error)

!! populate send concerns with data required by other processes
DO i=1, s_concern_size
  send_field_size = send_concern_max(i)%SIZE
  DO j=1, send_field_size

    x = send_concern_max(i)%x(j)
    y = send_concern_max(i)%y(j)

    CALL global_to_local_gridpt(x, y, src_grid)

    send_concern_max(i)%field(j) = src_field(x,y)
  END DO
END DO


! send lambda/phi point components telling appropriate
! proc what grid points you need
DO i=1, s_concern_size

  CALL mpl_isend(send_concern_max(i)%field,                        &
        send_concern_max(i)%SIZE, mpl_real,                        &
        send_concern_max(i)%proc_num, 0, my_comm,                  &
        requests(i), ierror(i))

END DO

! recv lambda/phi point components which tells you
! what grid points other procs need
DO i=1, r_concern_size

  CALL mpl_irecv(recv_concern_max(i)%field,                        &
       recv_concern_max(i)%SIZE, mpl_real,                         &
       recv_concern_max(i)%proc_num, 0, my_comm,                   &
       requests(i+s_concern_size), ierror(i+s_concern_size))

END DO

! check for mpl errors
CALL error_check_mpl(ierror, request_size,  &
     error1)

ALLOCATE(statuses(mpl_status_size, request_size))
statuses = 0

CALL mpl_waitall(request_size, requests, statuses, i)

! check for mpl errors
! request_size = 0 for 1x1 decomposition
IF (request_size /= 0) THEN
  error_i(1) = i
  CALL error_check_mpl(error_i, 1, error2)
END IF

IF (error1 < 0 .OR. error2 < 0) THEN
  error = -1
  cmessage = "REGRID_SWAP_DATA_MAX: MPL error"
END IF

! release request mem resource
DEALLOCATE(requests)
DEALLOCATE(statuses)
DEALLOCATE(ierror)

RETURN
END SUBROUTINE swap_regrid_data_max

END MODULE regrid_swap
