! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE rad3d_inp_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'RAD3D_INP_MOD'
CONTAINS

SUBROUTINE rad3d_inp(                                             &
  L_complete_North, L_complete_South, L_complete_deg,             &
  row_length, rows, off_x, off_y, first_row, last_row,            &
  first_data_interp, ES_space_interp,                             &
  levels,                                                         &
  RadField                                                        &
  )

!  Subroutine to complete 3D radiation fields.

USE mpp_conf_mod,          ONLY: swap_field_is_scalar
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Field_Types
USE interp_field_mod, ONLY: interp_field
IMPLICIT NONE
!
! Description: At points where full radiation calculations were not
!   done, radiation fields are completed by copying or interpolation.
!   This routine operates on 3D horizontal fields.
!
! Method: Radiation fields are calculated only at the first point on
!   the northern or southern polar rows in a full global configuration:
!   they must be copied to the other points along the row. When spatial
!   degradation is in operation, fields are calculated only at
!   alternate points and must be completed by interpolation.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: radiation control
!
! Declarations:
!
! Global variables:
!
! Arguments of the subroutine
!
LOGICAL, INTENT(IN) :: L_complete_North
!                              Completion required at North pole
LOGICAL, INTENT(IN) :: L_complete_South
!                              Completion required at South pole
LOGICAL, INTENT(IN) :: L_complete_deg
!                              Completion by interpolation is
!                              required at points where radiation
!                              was not calculated because of
!                              spatial degradation
!
INTEGER, INTENT(IN) :: row_length
!                              Number of points along a row
!                              in the domain
INTEGER, INTENT(IN) :: rows
!                              Number of rows in the domain
INTEGER, INTENT(IN) :: levels
!                              Number of vertical levels in the domain
INTEGER, INTENT(IN) :: off_x
!                              Size of small halo in i-direction
INTEGER, INTENT(IN) :: off_y
!                              Size of small halo in j-direction
!
INTEGER, INTENT(IN) :: first_row
!                              First row requiring interpolation
!                              under spatial degradation
INTEGER, INTENT(IN) :: last_row
!                              Last row requiring interpolation
!                              under spatial degradation
INTEGER, INTENT(IN) :: first_data_interp
!                              The first data point, in data
!                              coordinates, which needs needs to be
!                              interpolated in the local domain
!
REAL, INTENT(IN)    :: ES_space_interp(4, row_length, rows)
!                              Coefficients for spatial interpolation
!                              of radiation quantities
REAL, INTENT(INOUT) :: RadField(row_length, rows, levels)
!                              Radiation field to be completed
!
!
!     Local Variables:
!
!
INTEGER :: i
!                              Loop variable
INTEGER :: j
!                              Loop variable
INTEGER :: k
!                              Loop variable
REAL    :: RadField_inp(1-off_x:row_length+off_x,                 &
                     1-off_y:rows+off_y, levels)
!                              Field for interpolation including halos

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RAD3D_INP'

!
!- End of header
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!
!     Fill in at the north and south poles as required.
IF ( L_complete_North ) THEN
  DO k = 1, levels
    DO i = 2, row_length
      RadField(i,rows,k) = RadField(1,rows,k)
    END DO
  END DO
END IF
IF ( L_complete_South ) THEN
  DO k = 1, levels
    DO i = 2, row_length
      RadField(i,1,k) = RadField(1,1,k)
    END DO
  END DO
END IF
!
!
!     Fill in fields when spatial degradation is in operation.
IF ( L_complete_deg  ) THEN
  !
  DO k = 1, levels
    DO j = 1, rows
      DO i = 1, row_length
        RadField_inp(i,j,k) = RadField(i,j,k)
      END DO
    END DO
  END DO
  ! DEPENDS ON: swap_bounds
  CALL Swap_Bounds(RadField_inp, row_length, rows,                &
    levels, off_x, off_y, fld_type_p, swap_field_is_scalar)

  ! DEPENDS ON: set_external_halos
  CALL set_external_halos(RadField_inp,row_length, rows,          &
                          levels,off_x,off_y,0.0)

  CALL Interp_Field(ES_space_interp, RadField_inp,                &
    levels, first_data_interp, row_length, rows,                  &
    first_row, last_row, off_x, off_y)
  !
  DO k = 1, levels
    DO j = 1, rows
      DO i = 1, row_length
        RadField(i,j,k) = RadField_inp(i,j,k)
      END DO
    END DO
  END DO
  !
END IF
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rad3d_inp
END MODULE rad3d_inp_mod
