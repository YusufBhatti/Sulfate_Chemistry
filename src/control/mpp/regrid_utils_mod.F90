! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Module: RIV_ROUTE_UTILS --------------------------------------------
!
!  Purpose: Provides convenience routines and functions used for
!  parallel regridding of fields

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

MODULE regrid_utils

USE UM_ParVars, ONLY: g_datastart, g_datastartr, g_blsize
USE field_types, ONLY: fld_type_p, fld_type_r
USE ereport_mod, ONLY: ereport

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
IMPLICIT NONE


INTEGER, PARAMETER :: atmos_grid = 1, riv_grid = 2

CONTAINS

  ! returns what processor a grid point lie in for
  ! a given pixel and grid type, returns -1 if point
  ! could not be found
INTEGER FUNCTION get_proc_for_gridpt(x, y, grid)

USE UM_ParVars, ONLY: nproc_x, nproc_y, g_pe_index_EW,            &
     g_pe_index_NS
IMPLICIT NONE

INTEGER, INTENT(IN) :: grid
! grid type (e.g. atmso, river ...)
INTEGER, INTENT(IN) :: x, y
! x is lambda position and y phi position

! locals
INTEGER :: proc_num
INTEGER :: slam, sphi
INTEGER :: i, j
LOGICAL :: found


! search down lambda
found = .FALSE.

IF (grid == atmos_grid) THEN

  i = g_pe_index_EW(x)
  j = g_pe_index_NS(y)
  proc_num = j*nproc_x + i

ELSE IF (grid == riv_grid) THEN

  slam = 0
  sphi = 0
  i = 1
  j = 1

  DO WHILE((.NOT. found) .AND. (i <= nproc_x))

    slam = slam + g_blsize(1, fld_type_r, i)

    IF (x <= slam) THEN
      found = .TRUE.
    END IF

    IF (.NOT. found) i = i + 1
  END DO

  ! the search down phi
  found = .FALSE.

  DO WHILE((.NOT. found) .AND. (j <= nproc_y))

    sphi = sphi + g_blsize(2, fld_type_r, j)

    IF (y <= sphi) THEN
      found = .TRUE.
    END IF

    IF (.NOT. found) j = j + 1
  END DO

  proc_num = (j-1)*nproc_x + i

  ! should always find proc either yourself or another proc
  ! is condition below is true then something is wrong
  IF ((i>nproc_x) .OR. (j>nproc_y)) THEN
    proc_num = -1
  END IF

ELSE
  proc_num  = -1
END IF

get_proc_for_gridpt = proc_num
END FUNCTION get_proc_for_gridpt

! returns the global coordinates as process local subdomain
! coordinates (local index values)
SUBROUTINE global_to_local_gridpt(x, y, grid)

USE UM_ParCore, ONLY: mype
IMPLICIT NONE

INTEGER, INTENT(INOUT) :: x, y
! the global lambda and phi coordinate to convert
! to local coordinates

INTEGER, INTENT(IN) :: grid
! set to choose atmosphere or river grid
INTEGER :: err_code     ! error code to pass to ereport.

IF (grid == atmos_grid) THEN

  x = x - g_datastart(1, mype) + 1
  y = y - g_datastart(2, mype) + 1

ELSE IF (grid == riv_grid) THEN
  x = x - g_datastartr(1, mype) + 1
  y = y - g_datastartr(2, mype) + 1

ELSE

  err_code = -1
  CALL ereport("GLOBAL_TO_LOCAL_GRIDPT", err_code,                &
       "incorrect grid number")

END IF

END SUBROUTINE global_to_local_gridpt

SUBROUTINE local_to_global_gridpt(x, y, grid)

USE UM_ParCore, ONLY: mype
IMPLICIT NONE

INTEGER, INTENT(INOUT) :: x, y
! the local lambda and phi coordinate to convert
! to global coordinates
INTEGER, INTENT(IN) :: grid
! set to choose atmosphere or river grid
INTEGER :: err_code     ! error code to pass to ereport.

IF (grid == atmos_grid) THEN
  x = x + g_datastart(1, mype) - 1
  y = y + g_datastart(2, mype) - 1

ELSE IF (grid == riv_grid) THEN
  x = x + g_datastartr(1, mype) - 1
  y = y + g_datastartr(2, mype) - 1
ELSE
  err_code = -2
  CALL ereport("GLOBAL_TO_LOCAL_GRIDPT", err_code,                &
       "incorrect grid number")
END IF

END SUBROUTINE local_to_global_gridpt

! determine if global grid point is within subdomain of calling
! process
LOGICAL FUNCTION gridpt_outside_proc_domain(x, y, grid)

USE UM_ParCore, ONLY: mype
IMPLICIT NONE

INTEGER, INTENT(INOUT) :: x, y
! the global lambda and phi grid point
INTEGER, INTENT(IN) :: grid
! set to choose type of grid
INTEGER :: err_code     ! error code to pass to ereport.

LOGICAL :: ok

ok = .FALSE.

IF (grid == atmos_grid) THEN

  IF ((x < g_datastart(1,mype)) .OR. (x > (g_datastart(1, mype) +    &
       g_blsize(1, fld_type_p, mype)-1))) ok = .TRUE.

  IF ((y < g_datastart(2,mype)) .OR. (y > (g_datastart(2, mype) +    &
       g_blsize(2, fld_type_p, mype)-1))) ok = .TRUE.

ELSE
  err_code = -3
  CALL ereport("GLOBAL_TO_LOCAL_GRIDPT", err_code,                   &
  "GRID TYPE NOT SUPPORTED")
END IF

gridpt_outside_proc_domain = ok

END FUNCTION gridpt_outside_proc_domain

! Convenience routine which scans mpl
! ierror array to determine if an error
! has occured. If an element of ierror is
! is not equal to zero, error is set to -1 otherwise to 0
SUBROUTINE error_check_mpl(ierror, lenl, error)

IMPLICIT NONE
INTEGER, INTENT(IN) :: lenl, ierror(lenl)
! array of mpl error codes

INTEGER, INTENT(OUT) :: error
! result of error check, -1=present of error
! 0 equals no error

! local variable

INTEGER :: i

error = 0

DO i=1, lenl
  IF (ierror(i) /= 0) THEN
    error = 1
    RETURN
  END IF
END DO

RETURN
END SUBROUTINE error_check_mpl


!! searches recvConcern for a particular global grid
!! point. If found the value for the point is
!! retrieved and stored in val, and param found is set
!! to true on return otherwise false
SUBROUTINE find_value(x, y, recv_size, recv_concern_v, val, found)

USE regrid_types
IMPLICIT NONE

INTEGER, INTENT(IN) :: x, y
! lambda and phi position

LOGICAL, INTENT(OUT) :: found
! set to true if value found, false otherwise

REAL, INTENT(OUT) :: val
! the value found

INTEGER, INTENT(IN) :: recv_size
! size of recv_concern array

TYPE (concern), INTENT(IN) :: recv_concern_v(recv_size)
! concern that is searched for grid point x,y

! local
INTEGER :: i, j, k, LEN

found = .FALSE.
val = 0.0

! search all processors and for
! grid pt val

j = 1
DO WHILE((j <= recv_size) .AND. (.NOT. found))
  LEN = recv_concern_v(j)%SIZE

  DO i=1, LEN
    IF ((recv_concern_v(j)%y(i)==y) .AND. (recv_concern_v(j)%x(i)    &
         == x)) THEN
      found = .TRUE.
      val = recv_concern_v(j)%field(i)
    END IF
  END DO
  j = j + 1
END DO

END SUBROUTINE find_value

! provided to retrieve the field value at (xpos, ypos) in concern_max
! passed returns -1 if not found
REAL FUNCTION get_val_from_concern(recv_concern_max, recv_size,    &
     xpos, ypos, proc, error)

USE regrid_types
IMPLICIT NONE

INTEGER, INTENT(IN) :: recv_size
TYPE(concern_max), INTENT(IN) :: recv_concern_max(recv_size)
!
INTEGER, INTENT(IN) :: xpos, ypos, proc
INTEGER, INTENT(OUT) :: error

! local variable
INTEGER :: i, j, k, npts

error = 0

! search through all processors
DO i=1, recv_size

  IF (recv_concern_max(i)%proc_num == proc) THEN
    npts = recv_concern_max(i)%SIZE
    DO j=1, npts
      IF ((recv_concern_max(i)%x(j) == xpos) .AND. (               &
             recv_concern_max(i)%y(j) &
              == ypos) .AND. recv_concern_max(i)%contribute(j)) THEN
        get_val_from_concern = recv_concern_max(i)%field(j)
        RETURN
      END IF
    END DO
  END IF

END DO
error = -1
get_val_from_concern = 0.0
RETURN
END FUNCTION get_val_from_concern


!
! sorts each element of contribution info array
! that is weights are now in ascending order of their
! global 1-D index
!
SUBROUTINE sort_contributors(contribution, lenl, global_row_length,   &
     error)

USE regrid_types
IMPLICIT NONE


INTEGER, INTENT(IN) :: global_row_length, lenl
TYPE(contribution_info), INTENT(INOUT) :: contribution(lenl)
! the target grid local rows
INTEGER, INTENT(OUT) :: error
! set to 0 if not error occurs, to -1 if an error occurs

! local variables
INTEGER :: i, j, size_c, largest_size, rank, cg_index, g_index, k
INTEGER, ALLOCATABLE :: xtemp(:), ytemp(:), redirect(:), procTemp(:)
REAL, ALLOCATABLE :: weight_temp(:)

! set error flag
error = 0

size_c = 0
largest_size = 0

! get largest contribution info and allocate with
! that so we don't have to worry about src points with larger
! number of contributions

DO i=1, lenl

  size_c = contribution(i)%SIZE
  IF (size_c > largest_size) largest_size = size_c

END DO

! contributions cannot be less than zero
IF (size_c < 0) THEN
  WRITE(umMessage,*) "Error in SORT_CONTRIBUTION, SIZE < 0",             &
       ", Exiting!"
  CALL umPrint(umMessage,src='regrid_utils_mod')
  WRITE(umMessage,*) "SIZE: ", size_c
  CALL umPrint(umMessage,src='regrid_utils_mod')
  error = -1
  RETURN
END IF


ALLOCATE(xtemp(largest_size))
ALLOCATE(ytemp(largest_size))
ALLOCATE(weight_temp(largest_size))
ALLOCATE(procTemp(largest_size))
ALLOCATE(redirect(largest_size))

! iterate over each contribution info and arrange in
! ascending order (globalindex). Its really comparison sort
DO i=1, lenl

  size_c = contribution(i)%SIZE
  DO j=1, size_c
    g_index = (contribution(i)%y(j)-1)*global_row_length +          &
         contribution(i)%x(j)

    rank = 1
    ! determine target's rank
    DO k=1, size_c

      ! don't compare to yourself
      IF (k /= j) THEN
        cg_index = (contribution(i)%y(k)-1)*global_row_length +     &
             contribution(i)%x(k)

        IF (g_index > cg_index) rank = rank + 1

        ! should not be possible to have target appear twice
        IF (g_index == cg_index) THEN
          error = -1
          RETURN
        END IF
      END IF

    END DO

    redirect(j) = rank
  END DO

  ! now sort by applying redirection
  DO j=1, size_c
    xtemp(redirect(j)) = contribution(i)%x(j)
    ytemp(redirect(j)) = contribution(i)%y(j)
    weight_temp(redirect(j)) = contribution(i)%weight(j)
    procTemp(redirect(j)) = contribution(i)%contrib_proc(j)
  END DO

  DO j=1, size_c
    contribution(i)%x(j) = xtemp(j)
    contribution(i)%y(j) = ytemp(j)
    contribution(i)%weight(j) = weight_temp(j)
    contribution(i)%contrib_proc(j) = procTemp(j)
  END DO

END DO

! release resources used for sort
DEALLOCATE(xtemp)
DEALLOCATE(ytemp)
DEALLOCATE(weight_temp)
DEALLOCATE(redirect)
DEALLOCATE(procTemp)

RETURN
END SUBROUTINE sort_contributors

END MODULE regrid_utils


