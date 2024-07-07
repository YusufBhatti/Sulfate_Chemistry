! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module GLOBAL_2D_SUMS_MOD
MODULE global_2d_sums_mod


! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! Purpose:
!    This module contains a number of methods for performing global
!    2D (and by limitation of dimension 1D) sums across a group of
!    processors.
!
! Code Description:
!   Language: FORTRAN 90

USE UM_ParVars
USE missing_data_mod, ONLY: imdi
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Parameters for the different possible methods,
INTEGER, PARAMETER :: global_sum_reprod_orig = 1  ! Original method
INTEGER, PARAMETER :: global_sum_reprod_dd   = 2  ! "double double" method
                                                  ! in GCOM
INTEGER, PARAMETER :: global_sum_fast        = 3  ! Fast non-reprod method

INTEGER :: global_sum_method = imdi

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GLOBAL_2D_SUMS_MOD'

CONTAINS
SUBROUTINE global_2d_sums_sp(field, row_length, rows, off_x, off_y, &
     levels, global_sum, gid, global_sum_method_override)

  ! Purpose: this routine provides the global sums according to some
  !          the options chosen. Group ID is an optional argument but
  !          is ignored for the original method.


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE um_types, ONLY: real32, real64
IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)  :: row_length
INTEGER, INTENT(IN)  :: rows
INTEGER, INTENT(IN)  :: off_x
INTEGER, INTENT(IN)  :: off_y
INTEGER, INTENT(IN)  :: levels

REAL (KIND=real32),    INTENT(IN)  :: field(1-off_x : row_length+off_x,  &
                            1-off_y : rows+off_y,        &
                            levels)

REAL(KIND=real32),    INTENT(OUT) :: global_sum(levels)
INTEGER, INTENT(IN), OPTIONAL  :: gid   ! Group ID
INTEGER, INTENT(IN), OPTIONAL :: global_sum_method_override


! Local variables
INTEGER                       :: i,j,k
INTEGER                       :: istat
INTEGER                       :: my_gid
INTEGER                       :: global_sum_method_use
REAL (KIND=real32)            :: local_sum(MAX(levels,rows))
REAL (KIND=real64), ALLOCATABLE  :: global_sum64(:)   ! resultant double
                                                      ! precision sum
REAL (KIND=real64), ALLOCATABLE  :: field64(:, :, :)  ! dp copy of field

CHARACTER (LEN=errormessagelength)           :: cmessage
CHARACTER (LEN=*), PARAMETER  :: RoutineName='GLOBAL_2D_SUMS_SP'

INTEGER                       :: err_code

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Deal with optional group
IF (PRESENT(gid)) THEN
  my_gid = gid
ELSE
  ! Default position is to use the current global communicator
  CALL gc_get_communicator(my_gid, istat)
  IF (istat /= 0) THEN
    WRITE(cmessage,*) 'Problem in gc_get_communicator'
    err_code = 50
    CALL ereport(routinename, err_code, cmessage  )
  END IF

END IF

! If overriding namelist setting for global_sum_method then do it here
IF (PRESENT(global_sum_method_override)) THEN
  global_sum_method_use = global_sum_method_override
ELSE
  global_sum_method_use = global_sum_method
END IF

SELECT CASE (global_sum_method_use)
CASE (global_sum_reprod_orig)
  WRITE(cmessage,'(A)')       &
  'Original Reproducible global sum NOT implemented in 32bit, compile ' // &
  'with C_DP_HLM CPP or choose double-double precision '                // &
  'or fast, non-reproducible method.'
  err_code = 100
  CALL ereport(routinename, err_code, cmessage  )

CASE (global_sum_reprod_dd)
  ! As implementing the global 2d sums using single precision is difficult
  ! due to need to re-work GCOM, we promote to 64-bit here to do the sum
  ! noting that the majority of benefit will be found from the improved 
  ! cache use in the solver.
  ALLOCATE(field64(1-off_x : row_length+off_x, 1-off_y : rows+off_y, levels))
  ALLOCATE(global_sum64(levels))

  field64(:,:,:) = field(:,:,:)
  CALL gcg_r2darrsum(field64, row_length, rows, off_x, off_y, levels, &
                     my_gid, global_sum64, istat)
  global_sum(:) = global_sum64(:)

  DEALLOCATE(global_sum64)
  DEALLOCATE(field64)

CASE (global_sum_fast)
  ! Create local partial sums
  DO k = 1, levels
    local_sum(k) = 0.0

    DO j = 1, rows
      DO i = 1, row_length
        local_sum(k) = local_sum(k) + field(i,j,k)
      END DO
    END DO
  END DO

  ! Sum up partial sums
  CALL gcg_rsum_kind32(levels, my_gid, istat, local_sum)
  IF (istat /= 0) THEN
    WRITE(cmessage,*) 'Problem in gcg_rsum'
    err_code = 400
    CALL ereport(routinename, err_code, cmessage  )
  END IF
  global_sum(1:levels) = local_sum(1:levels)

CASE DEFAULT
  WRITE(cmessage,'(A,I6)') 'Unrecognised global 2D sum method ', &
                     global_sum_method_use
  err_code = 500
  CALL ereport(routinename, err_code, cmessage  )

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE global_2d_sums_sp

SUBROUTINE global_2d_sums(field, row_length, rows, off_x, off_y, &
                          levels, global_sum, gid)

! Purpose: this routine provides the global sums according to some
!          the options chosen. Group ID is an optional argument but
!          is ignored for the original method.


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)  :: row_length
INTEGER, INTENT(IN)  :: rows
INTEGER, INTENT(IN)  :: off_x
INTEGER, INTENT(IN)  :: off_y
INTEGER, INTENT(IN)  :: levels
INTEGER, INTENT(IN), OPTIONAL  :: gid   ! Group ID

REAL,    INTENT(IN)  :: field(1-off_x : row_length+off_x,  &
                              1-off_y : rows+off_y,        &
                              levels)

REAL,    INTENT(OUT) :: global_sum(levels)

! Local variables
INTEGER                       :: i,j,k
INTEGER                       :: istat
INTEGER                       :: my_gid
REAL                          :: local_sum(MAX(levels,rows))
CHARACTER (LEN=errormessagelength)           :: cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName='GLOBAL_2D_SUMS'

INTEGER :: err_code           ! error code used in call to ereport.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Deal with optional group
IF (PRESENT(gid)) THEN
  my_gid = gid
ELSE
  ! Default position is to use the current global communicator
  CALL gc_get_communicator(my_gid, istat)
  IF (istat /= 0) THEN
    WRITE(cmessage,*) 'Problem in gc_get_communicator'
    err_code = 50
    CALL ereport(routinename, err_code, cmessage  )
  END IF

END IF


SELECT CASE (global_sum_method)
CASE (global_sum_reprod_orig)

  ! Loop over levels since we could have an offset for rows so summation is
  ! not in memory order which GCOM expects since the offset is between
  ! levels.
  DO k = 1, levels

    ! Sum first along rows
    CALL gcg_rvecsumr(row_length+2*off_x, row_length, 1+off_x,       &
                      rows, field(:,1:,k), gc_proc_row_group, &
                      istat, local_sum)
    IF (istat /= 0) THEN
      WRITE(cmessage,*) 'Problem in gcg_rvecsumr'
      err_code = 100
      CALL ereport(routinename, err_code, cmessage  )
    END IF

    ! If we have more than one row, sum up the rows
    IF (rows > 1) THEN
      CALL gcg_rvecsumr(rows, rows, 1, 1, local_sum,            &
                        gc_proc_col_group, istat, global_sum(k))
      IF (istat /= 0) THEN
        WRITE(cmessage,*) 'Problem in 2nd gcg_rvecsumr'
        err_code = 200
        CALL ereport(routinename, err_code, cmessage  )
      END IF
    ELSE
      global_sum(k) = local_sum(1)
    END IF
  END DO

CASE (global_sum_reprod_dd)
  ! Use the GCOM routine to do this
  CALL gcg_r2darrsum(field, row_length, rows, off_x, off_y, levels, &
                     my_gid, global_sum, istat)
  IF (istat /= 0) THEN
    WRITE(cmessage,*) 'Problem in gcg_r2darrsum'
    err_code = 300
    CALL ereport(routinename, err_code, cmessage  )
  END IF

CASE (global_sum_fast)
  ! Create local partial sums
  DO k = 1, levels
    local_sum(k) = 0.0

    DO j = 1, rows
      DO i = 1, row_length
        local_sum(k) = local_sum(k) + field(i,j,k)
      END DO
    END DO
  END DO

  ! Sum up partial sums
  CALL gcg_rsum(levels, my_gid, istat, local_sum)
  IF (istat /= 0) THEN
    WRITE(cmessage,*) 'Problem in gcg_rsum'
    err_code = 400
    CALL ereport(routinename, err_code, cmessage  )
  END IF

  global_sum(1:levels) = local_sum(1:levels)

CASE DEFAULT
  WRITE(cmessage,'(A,I6)') 'Unrecognised global 2D sum method ', &
                           global_sum_method
  err_code = 500
  CALL ereport(routinename, err_code, cmessage  )

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE global_2d_sums

END MODULE global_2d_sums_mod
