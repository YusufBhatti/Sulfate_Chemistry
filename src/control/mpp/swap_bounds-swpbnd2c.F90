#if defined(C96_1C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose:
!   This subroutine takes care of all boundary swapping and
!   extending of arrays at the global boundaries. It is a wrapper,
!   used when 'USE halo_exchange' is absent from the calling routine.
!   Core functionality now in halo_exchange.F90
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: MPP
!
!   Attention, this routine is not thread-safe 
!   when l_swap_bounds_ddt is .TRUE.

SUBROUTINE swap_bounds(                &
    field, row_length, rows, levels,   & ! field
    halo_x, halo_y,                    & ! halos
    field_type, l_vector               & ! supporting information
    )

USE mpl, ONLY:                      &
    mpl_real,                        &
    mpl_status_size
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Halo_Exchange, ONLY: &
    Swap_bounds_hub,       &
    he_initialised
USE ereport_mod,   ONLY: ereport
USE um_parvars,    ONLY: sb_model_type
USE model_domain_mod, ONLY: mt_lam

IMPLICIT NONE

INTEGER, INTENT(IN) ::  row_length ! number of points on a row
                                   ! (not including halos)
INTEGER, INTENT(IN) ::  rows       ! number of rows in a theta field
                                   ! (not including halos)
INTEGER, INTENT(IN) ::  levels     ! number of model levels
INTEGER, INTENT(IN) ::  halo_x     ! size of halo in "i" direction
INTEGER, INTENT(IN) ::  halo_y     ! size of halo in "j" direction
INTEGER, INTENT(IN) ::  field_type ! Defines the grid interpolation type
                                   ! of the input FIELD (u,v or w)
LOGICAL, INTENT(IN) :: l_vector    ! TRUE:  Data is a horizontal vector
                                   !        component
                                   ! FALSE: Data is a scalar

REAL, INTENT(INOUT) :: field(1-halo_x:row_length+halo_x,          &
    1-halo_y:rows+halo_y,levels)   ! Field to have its halos updated

INTEGER :: errorstatus

! Local variables
LOGICAL :: do_east    ! perform east swap
LOGICAL :: do_west    ! perform west swap
LOGICAL :: do_south   ! perform south swap
LOGICAL :: do_north   ! perform north swap


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SWAP_BOUNDS'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

errorstatus = 0

IF (halo_x == 0 .AND. halo_y == 0 .AND. rows == 1 &
                                  .AND. row_length == 1) THEN
    !Choose not to output warning as the SCM output becomes massive
    !errorstatus = -1
    ! CALL ereport('SWAP_BOUNDS',errorstatus,       &
    !            ' not initialised as this is a SCM.')
  HE_Initialised = .TRUE.
END IF


IF (.NOT. HE_Initialised) THEN
  errorstatus = 99
  CALL ereport('SWAP_BOUNDS',errorstatus,' not initialised')
END IF


CALL swap_bounds_hub( field, row_length, rows, levels,      &
                      halo_x, halo_y, field_type, l_vector)       


IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE swap_bounds
#endif
