! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! This module provides methods for filling external halos (those
!  around the edge of model domain) with sensible (copy of interior
!  points) numbers.  This is useful for LAM models when SWAPBOUNDS
!  doesn't fill these halos as they are normally filled with LBC
!  data. However, not all fields have LBC data applied to them, and
!  such fields need to have sensible numbers put in these external
!  halo regions.

MODULE fill_external_halos_mod

USE yomhook,      ONLY: lhook, dr_hook
USE parkind1,     ONLY: jprb, jpim
USE um_types,     ONLY: real32, real64
USE um_parvars,   ONLY: neighbour
USE um_parparams, ONLY: nodomain, pnorth, peast, psouth, pwest

IMPLICIT NONE

PRIVATE

PUBLIC :: fill_external_halos
INTERFACE fill_external_halos
MODULE PROCEDURE fill_external_halos_3D_dp,    & 
                 fill_external_halos_3D_sp,    &
                 fill_external_halos_3D_int,   &
                 fill_external_halos_3D_log,   & 
                 fill_external_halos_2D_log 

END INTERFACE


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FILL_EXTERNAL_HALOS_MOD'

CONTAINS


! Routine for filling in external halos using 3D double precision values
SUBROUTINE fill_external_halos_3D_dp(                             &
  field,row_length,rows,levels,                                   &
  halo_x, halo_y)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real64

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILL_EXTERNAL_HALOS_3D_DP'

#include "fill_external_halos.h"

RETURN
END SUBROUTINE fill_external_halos_3D_dp


! Routine for filling in external halos using 3D single precision values
SUBROUTINE fill_external_halos_3D_sp(                             &
  field,row_length,rows,levels,                                   &
  halo_x, halo_y)

IMPLICIT NONE

INTEGER, PARAMETER :: field_kind = real32

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILL_EXTERNAL_HALOS_3D_SP'

#include "fill_external_halos.h"

RETURN
END SUBROUTINE fill_external_halos_3D_sp


! Routine for filling in external halos using 3D integer
SUBROUTINE fill_external_halos_3D_int(                            &
  field,row_length,rows,levels,                                   &
  halo_x, halo_y)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILL_EXTERNAL_HALOS_3D_INT'

#define   FILL_EXTERNAL_HALOS_INTEGER 1
#include "fill_external_halos.h"
#undef    FILL_EXTERNAL_HALOS_INTEGER

RETURN
END SUBROUTINE fill_external_halos_3D_int


! Routine for filling in external halos using 3D logical
SUBROUTINE fill_external_halos_3D_log(                            &
  field,row_length,rows,levels,                                   &
  halo_x, halo_y)

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILL_EXTERNAL_HALOS_3D_LOG'

#define   FILL_EXTERNAL_HALOS_LOGICAL 1
#include "fill_external_halos.h"
#undef    FILL_EXTERNAL_HALOS_LOGICAL

RETURN
END SUBROUTINE fill_external_halos_3D_log


! Routine for filling in external halos using 2D logical
SUBROUTINE fill_external_halos_2D_log(                            &
  field,row_length,rows,levels,                                   &
  halo_x, halo_y)

USE ereport_mod,   ONLY: ereport

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length ! IN: number of points on a row
                                  !     (not including halos)
INTEGER, INTENT(IN) :: rows       ! IN: number of rows in a theta field
                                  !     (not including halos)
INTEGER, INTENT(IN) :: levels     ! IN: number of model levels
INTEGER, INTENT(IN) :: halo_x     ! IN: size of halo in "i" direction
INTEGER, INTENT(IN) :: halo_y     ! IN: size of halo in "j" direction

LOGICAL,               INTENT(INOUT) :: field(1-halo_x:row_length+halo_x,  &
                                              1-halo_y:rows+halo_y)
                                  ! IN/OUT : Field to have its halos updated


! Local variables

INTEGER  ::  i,j,k            ! loop indicies
INTEGER  ::  j_start,j_stop   ! loop indicies
INTEGER  ::  errorstatus

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILL_EXTERNAL_HALOS_2D_LOG'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (levels /= 1) THEN
  errorstatus = 1
  CALL ereport(RoutineName,errorstatus, ' levels/=1 in a 2D field')
END IF

CALL fill_external_halos_3D_log(field,row_length,rows,levels,halo_x,halo_y)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE fill_external_halos_2D_log


END MODULE fill_external_halos_mod
