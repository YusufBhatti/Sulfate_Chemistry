! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE generate_frame_mod
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

PRIVATE
PUBLIC :: generate_frame

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'GENERATE_FRAME_MOD'

CONTAINS

SUBROUTINE generate_frame(interp_weights, output_field)
!
! Description: Performs a cutout of the input data using the frame mask
!  contained in the interp_weights_type.
!
USE field_mod,              ONLY: field_type, ASSIGNMENT(=)
USE interp_weights_mod,     ONLY: interp_weights_type
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE missing_data_mod,       ONLY: rmdi
IMPLICIT NONE
! Arguments
TYPE(interp_weights_type), INTENT(IN) :: interp_weights
TYPE(field_type), INTENT(INOUT) :: output_field
! Local variables
REAL, ALLOCATABLE :: cutout_data(:,:,:)
INTEGER :: row,level
INTEGER :: num_levels, num_rows, num_cols, start_row, start_col, input_row_index, input_col_index
CHARACTER(LEN=*), PARAMETER :: routinename = 'GENERATE_FRAME'
CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER :: icode
INTEGER :: col_max

! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check that grid and data array sizes match
IF (output_field%get_num_cols() /= SIZE(output_field%rdata,1)) THEN
  icode = 10
  WRITE(cmessage, '(A,I0,A,I0)') "Number of columns of output field grid: ",    &
       output_field%get_num_cols(),                                             &
       " does not match the number of column in the output field data array: ", &
       SIZE(output_field%rdata,1)
  CALL ereport(routinename, icode, cmessage)
END IF
IF (output_field%get_num_rows() /= SIZE(output_field%rdata,2)) THEN
  icode = 20
  WRITE(cmessage, '(A,I0,A,I0)') "Number of rows of output field grid: ",       &
       output_field%get_num_rows(),                                             &
       " does not match the number of rows in the output field data array: ",   &
       SIZE(output_field%rdata,2)
  CALL ereport(routinename, icode, cmessage)
END IF
IF (.NOT. ALLOCATED(interp_weights%frame_mask)) THEN
  icode = 30
  WRITE(cmessage, '(A)') "Frames logical mask array is not allocted"
  CALL ereport(routinename, icode, cmessage)
END IF

! Perform cutout
CALL interp_weights%get_frames_domain(start_row, start_col, num_rows, num_cols)
ALLOCATE(cutout_data(num_cols, num_rows, output_field%get_num_levels()))

! The field will be copied from column "start_col" across, until it reaches
! either:
!   * the last column of the field, in which case the copy must cross the
!     meridian back to column 1 to continue
!   * or the required column to ensure "num_cols" are copied.

col_max = MIN(output_field%get_num_cols(),start_col+num_cols-1)

! Avoid using the function procedure inside the PARALLEL region to avoid a
! gfortran compiler bug (tested up to gfortran v6.3.0), potential race
! conditions with ifort and other possible compiler issues
num_levels = output_field%get_num_levels()

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& PRIVATE(level, row, input_row_index, input_col_index)                   &
!$OMP& SHARED(output_field, interp_weights, start_row, start_col, cutout_data, &
!$OMP&        num_rows, num_cols, col_max, num_levels)                         &
!$OMP& COLLAPSE(2)
DO level = 1, num_levels

  DO row = 1, num_rows

    input_row_index = start_row+row-1

    DO input_col_index = start_col, col_max

      ! Apply mask
      IF (interp_weights%frame_mask(input_col_index,input_row_index)) THEN
        cutout_data(input_col_index-start_col+1,row,level) =                   &
                       output_field%rdata(input_col_index,input_row_index,level)
      ELSE
        cutout_data(input_col_index-start_col+1,row,level) = rmdi
      END IF

    END DO

    ! Crossed meridian

    DO input_col_index = 1, num_cols - (col_max - start_col + 1)

      ! Apply mask
      IF (interp_weights%frame_mask(input_col_index,input_row_index)) THEN
        cutout_data(col_max+input_col_index-start_col+1,row,level) =           &
                       output_field%rdata(input_col_index,input_row_index,level)
      ELSE
        cutout_data(col_max+input_col_index-start_col+1,row,level) = rmdi
      END IF

    END DO

  END DO

END DO
!$OMP END PARALLEL DO

CALL MOVE_ALLOC(cutout_data, output_field%rdata)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE generate_frame

END MODULE generate_frame_mod
