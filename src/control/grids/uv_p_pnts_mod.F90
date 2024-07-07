! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Grids

MODULE uv_p_pnts_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UV_P_PNTS_MOD'

CONTAINS

SUBROUTINE uv_p_pnts_halo                                                      &
  ( u_single, v_single, cos_theta_longitude, sin_theta_longitude               &
  , global_row_length, proc_row_group, u_p, v_p )

! Purpose:
! Take a single level field of u and v and get u and v on the p grid.
! For use with diagnostic 3.230, set_seasalt and dms_flux
! and 3.436 Wind gust
!
! Code Description:
!   Language: Fortran 95.
!   This code is written to UM programming standards v10.3.

USE model_domain_mod,      ONLY: model_type, mt_global
USE atm_fields_bounds_mod, ONLY: udims, vdims, udims_s, vdims_s                &
                               , pdims, pdims_s, tdims

USE um_parvars,   ONLY: at_extremity
USE um_parparams, ONLY: pnorth, peast, psouth, pwest
USE u_to_p_mod,   ONLY: u_to_p
USE v_to_p_mod,   ONLY: v_to_p

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.
REAL, INTENT(IN) :: u_single ( udims_s%i_start:udims_s%i_end                   &
                             , udims_s%j_start:udims_s%j_end )

REAL, INTENT(IN) :: v_single ( vdims_s%i_start:vdims_s%i_end                   &
                             , vdims_s%j_start:vdims_s%j_end )

REAL, INTENT(IN) :: cos_theta_longitude ( tdims%i_start:tdims%i_end            &
                                        , tdims%j_start:tdims%j_end )
REAL, INTENT(IN) :: sin_theta_longitude ( tdims%i_start:tdims%i_end            &
                                        , tdims%j_start:tdims%j_end )

INTEGER,INTENT(IN) :: global_row_length  ! Number of points on a global row
INTEGER,INTENT(IN) :: proc_row_group

REAL, INTENT(OUT)  :: u_p ( pdims%i_start:pdims%i_end                          &
                          , pdims%j_start:pdims%j_end )
REAL, INTENT(OUT)  :: v_p ( pdims%i_start:pdims%i_end                          &
                          , pdims%j_start:pdims%j_end )

! The following are declared as arrays rather than scalars to maintain
! consistency of argument type, since they are passed to a general
! purpose subroutine in which they are array arguments
REAL :: mag_vector_np(1)
REAL :: mag_vector_sp(1)
REAL :: dir_vector_np(1)
REAL :: dir_vector_sp(1)

INTEGER :: i, j, i_field, y, z     ! Counters

! Temporarily create these as other subroutines use them
! can be removed when all endgame changes are made
INTEGER :: offx, offy, row_length, n_rows

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UV_P_PNTS_HALO'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Temporarily create these as other subroutines use them
! can be removed when all endgame changes are made
offx = pdims_s%i_end - pdims%i_end
offy = pdims_s%j_end - pdims%j_end

row_length = pdims%i_end
n_rows     = vdims%j_len

! Set allocated arrays to zero to make sure any data left in the memory
! allocated to them is purged
u_p(:,:) = 0.0
v_p(:,:) = 0.0

! Interpolate u and v to p grid.
CALL u_to_p                                                                    &
  ( u_single, udims_s%i_start, udims_s%i_end, udims_s%j_start, udims_s%j_end   &
  , pdims%i_start, pdims%i_end, pdims%j_start, pdims%j_end, 1                  &
  , at_extremity,u_p )


CALL v_to_p                                                                    &
  ( v_single, vdims_s%i_start, vdims_s%i_end, vdims_s%j_start, vdims_s%j_end   &
  , pdims%i_start, pdims%i_end, pdims%j_start, pdims%j_end, 1                  &
  , at_extremity,v_p )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE uv_p_pnts_halo



SUBROUTINE uv_p_pnts                                                           &
  ( u_single, v_single, cos_theta_longitude, sin_theta_longitude               &
  , global_row_length, proc_row_group, u_p, v_p )

! Purpose:
! Take a single level field of u and v and get u and v on the p grid.
! For use with diagnostic 3.230, set_seasalt and dms_flux
! and 3.436 Wind gust

! Code Description:
!   Language:Fortran 95.
!   This code is written to UM programming standards v10.3.

USE swap_bounds_2d_mv_mod, ONLY: swap_bounds_2d_mv
USE atm_fields_bounds_mod, ONLY:   udims, vdims, udims_s, vdims_s,     &
                                   pdims, pdims_s, tdims
USE swapable_field_mod, ONLY:     swapable_field_pointer_type
USE Field_Types, ONLY:             fld_type_u, fld_type_v

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.
REAL, INTENT(IN) :: u_single ( udims%i_start:udims%i_end                       &
                             , udims%j_start:udims%j_end )

REAL, INTENT(IN) :: v_single ( vdims%i_start:vdims%i_end                       &
                             , vdims%j_start:vdims%j_end )

REAL, INTENT(IN) :: cos_theta_longitude ( tdims%i_start:tdims%i_end            &
                                        , tdims%j_start:tdims%j_end )

REAL, INTENT(IN) :: sin_theta_longitude ( tdims%i_start:tdims%i_end            &
                                        , tdims%j_start:tdims%j_end )

INTEGER, INTENT(IN) :: global_row_length ! Number of points on a global row
INTEGER, INTENT(IN) :: proc_row_group

REAL, INTENT(OUT) :: u_p ( pdims%i_start:pdims%i_end                           &
                         , pdims%j_start:pdims%j_end )

REAL, INTENT(OUT) :: v_p ( pdims%i_start:pdims%i_end                           &
                         , pdims%j_start:pdims%j_end )

! Local alloctable arrays
REAL, ALLOCATABLE, TARGET :: u_halo(:,:)
REAL, ALLOCATABLE, TARGET :: v_halo(:,:)

! The following are declared as arrays rather than scalars to maintain
! consistency of argument type, since they are passed to a general
! purpose subroutine in which they are array arguments
REAL :: mag_vector_np(1)
REAL :: mag_vector_sp(1)
REAL :: dir_vector_np(1)
REAL :: dir_vector_sp(1)

INTEGER :: i, j, i_field, y, z     ! Counters

! Temporarily create these as other subroutines use them
! can be removed when all endgame changes are made
INTEGER :: offx, offy, row_length, n_rows

TYPE(swapable_field_pointer_type) :: fields_to_swap(2)   ! For multivar
                                                         ! swap_bounds
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'UV_P_PNTS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Temporarily create these as other subroutines use them
! can be removed when all endgame changes are made
offx = pdims_s%i_end-pdims%i_end
offy = pdims_s%j_end-pdims%j_end

row_length = pdims%i_end
n_rows     = vdims%j_len

! Set allocated arrays to zero to make sure any data left in the memory
! allocated to them is purged
u_p(:,:) = 0.0
v_p(:,:) = 0.0

! Allocate arrays to hold u and v on u/v grid with halos
ALLOCATE ( u_halo( udims_s%i_start:udims_s%i_end                               &
                 , udims_s%j_start:udims_s%j_end) )

ALLOCATE ( v_halo( vdims_s%i_start:vdims_s%i_end                               &
                 , vdims_s%j_start:vdims_s%j_end) )


! Set allocated arrays to zero to make sure any data left in the memory
! allocated to them is purged
u_halo(:,:) = 0.0
v_halo(:,:) = 0.0


! Set non halo points to input values of u and v
DO y=udims%j_start, udims%j_end
  DO z=udims%i_start, udims%i_end
    u_halo(z,y) = u_single(z,y)
  END DO
END DO

DO y=vdims%j_start, vdims%j_end
  DO z=vdims%i_start, vdims%i_end
    v_halo(z,y) = v_single(z,y)
  END DO
END DO

! Update halos for u_halo and v_halo
i_field = 0
i_field = i_field + 1
fields_to_swap(i_field) % field_2d    => u_halo(:,:)
fields_to_swap(i_field) % field_type  =  fld_type_u
fields_to_swap(i_field) % levels      =  1
fields_to_swap(i_field) % rows        =  udims%j_end
fields_to_swap(i_field) % vector      =  .TRUE.

i_field = i_field + 1
fields_to_swap(i_field) % field_2d    => v_halo(:,:)
fields_to_swap(i_field) % field_type  =  fld_type_v
fields_to_swap(i_field) % levels      =  1
fields_to_swap(i_field) % rows        =  n_rows
fields_to_swap(i_field) % vector      =  .TRUE.

CALL swap_bounds_2d_mv( fields_to_swap, i_field, pdims%i_end, offx, offy )


CALL uv_p_pnts_halo( u_halo, v_halo, cos_theta_longitude, sin_theta_longitude  &
                   , global_row_length, proc_row_group, u_p, v_p )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE uv_p_pnts

END MODULE uv_p_pnts_mod
