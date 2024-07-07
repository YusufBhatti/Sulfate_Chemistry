! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE therm_advec_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='THERM_ADVEC_MOD'

CONTAINS

SUBROUTINE therm_advec(NumPlevs, u_p, v_p, temp_plevs, thermal_advec)

! Description: Routine to calculate thermal advection
!
! Method: thermal_advec = -U . grad(temp) ; computed at theta points -
!                                     use average of u, v for each cell
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE atm_fields_bounds_mod
USE atm_fields_real_mod, ONLY: exner_theta_levels, theta, p_theta_levels

USE planet_constants_mod, ONLY: planet_radius
USE cderived_mod, ONLY: delta_lambda, delta_phi
USE field_types, ONLY: fld_type_p, fld_type_u, fld_type_v
USE nlsizes_namelist_mod, ONLY: model_levels, row_length, rows, n_rows, &
                                global_row_length, bl_levels
USE um_parvars, ONLY: offx, offy
USE trignometric_Mod, ONLY: cos_theta_latitude
USE mpp_conf_mod,  ONLY: swap_field_is_vector, swap_field_is_scalar

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumPLevs  ! Number of P levels for output

! u-wind on B-grid p levs
REAL, INTENT(IN) :: u_p(udims%i_start:udims%i_end, &
                         vdims%j_start:vdims%j_end, NumPLevs) 

! v-wind on B-grid p levs
REAL, INTENT(IN) :: v_p(udims%i_start:udims%i_end, &
                         vdims%j_start:vdims%j_end, NumPLevs) 

! Temperature p levs
REAL, INTENT(IN) :: temp_plevs(tdims%i_start:tdims%i_end, &
                               tdims%j_start:tdims%j_end, NumPLevs) 

! Thermal advection
REAL, INTENT(OUT) :: thermal_advec(tdims%i_start:tdims%i_end, &
                                   tdims%j_start:tdims%j_end, NumPLevs) 

! Local Constants:
REAL :: HdelX  ! 1/2*W-E grid spacing -  only valid for regular grid
REAL :: HdelY  ! 1/2*S-N grid spacing -  only valid for regular grid

! Local variables:

INTEGER :: i, j, k

!   Temp with halo
REAL,ALLOCATABLE :: t_halo(:,:)

! u_p, v_p with haloes
REAL,ALLOCATABLE :: u_p_halo(:,:)
REAL,ALLOCATABLE :: v_p_halo(:,:)

!   Temp gradients - longitudinal and latitudinal
REAL,ALLOCATABLE :: dtdx (:,:)
REAL,ALLOCATABLE :: dtdy (:,:)

REAL :: ave_up
REAL :: ave_vp


CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'THERM_ADVEC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE(dtdx(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
ALLOCATE(dtdy(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
ALLOCATE(t_halo(tdims_s%i_start:tdims_s%i_end, tdims_s%j_start:tdims_s%j_end))
ALLOCATE(u_p_halo(udims_s%i_start:udims_s%i_end, vdims_s%j_start:vdims_s%j_end))
ALLOCATE(v_p_halo(udims_s%i_start:udims_s%i_end, vdims_s%j_start:vdims_s%j_end))

thermal_advec(:,:,:) = 0.0
dtdx(:,:) = 0.0
dtdy(:,:) = 0.0
t_halo  (:,:) = 0.0
u_p_halo(:,:) = 0.0
v_p_halo(:,:) = 0.0

! Inverse latitude arc interval
HdelY = 0.5/(Planet_Radius * delta_phi)

DO  k = 1, NumPLevs

  ! Populate t_halo
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start , tdims%i_end
      t_halo(i,j) = temp_plevs(i,j,k)
    END DO
  END DO

  CALL swap_bounds(t_halo, tdims%i_len, tdims%j_len, 1, offx, offy, &
                   fld_type_p, swap_field_is_scalar)

  ! Populate u_p_halo, v_p_halo
  DO j = vdims%j_start, vdims%j_end
    DO i = udims%i_start , udims%i_end
      u_p_halo(i,j) = u_p(i,j,k)
      v_p_halo(i,j) = v_p(i,j,k)
    END DO
  END DO

  CALL swap_bounds(u_p_halo, udims%i_len, vdims%j_len, 1, offx, offy, &
                   fld_type_u,swap_field_is_vector)

  CALL swap_bounds(v_p_halo, udims%i_len, vdims%j_len, 1, offx, offy, &
                   fld_type_v,swap_field_is_vector)

  ! Compute dtdx at theta points
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start , tdims%i_end
      ! Inverse longitude arc interval
      HdelX = 0.5/(Planet_Radius * delta_lambda * cos_theta_latitude(i,j))
      dtdx(i,j) = (t_halo(i+1,j) -t_halo(i-1,j)) * HdelX
    END DO ! i
  END DO ! j

  ! Compute dtdy at theta points
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start , tdims%i_end
      dtdy(i,j) = (t_halo(i,j+1) -t_halo(i,j-1)) * HdelY
    END DO ! i
  END DO ! j

  ave_up = 0.0
  ave_vp = 0.0

  ! Compute thermal advection at theta points; haloes required to prevent 
  ! out of bounds
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start , tdims%i_end
      !Compute average u_p, v_p for the cell
      ave_up = 0.25*(u_p_halo(i-1,j-1) + u_p_halo(i,j-1) &
                                       + u_p_halo(i-1,j) + u_p_halo(i,j))
      ave_vp = 0.25*(v_p_halo(i-1,j-1) + v_p_halo(i,j-1) &
                                       + v_p_halo(i-1,j) + v_p_halo(i,j))
      ! Thermal advection
      thermal_advec(i,j,k) = -ave_up*dtdx(i,j) -ave_vp*dtdy(i,j)

    END DO ! i
  END DO ! j

END DO ! k

DEALLOCATE(dtdx)
DEALLOCATE(dtdy)
DEALLOCATE(t_halo)
DEALLOCATE(u_p_halo)
DEALLOCATE(v_p_halo)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE therm_advec

END MODULE therm_advec_mod
