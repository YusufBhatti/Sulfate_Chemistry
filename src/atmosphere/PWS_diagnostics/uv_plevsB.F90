! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE uv_plevsB_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UV_PLEVSB_MOD'

CONTAINS
!
!  Routine to calculate u, v on p levels, B grid

SUBROUTINE uv_plevsB(NumPLevs, PLevs_ws, UFields, VFields, UFieldsP, VFieldsP)

! Description: Routine to calculate u, v on p levels, B grid
!
! Method: vertical interpolation followed by horiz interp, 
!         using standard utility routines
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE atm_fields_bounds_mod
USE atm_fields_mod, ONLY: exner_rho_levels

USE uc_to_ub_mod    , ONLY: uc_to_ub
USE vc_to_vb_mod    , ONLY: vc_to_vb
USE vert_interp2_mod, ONLY: vert_interp2
USE p_to_u_mod      , ONLY: p_to_u
USE p_to_v_mod      , ONLY: p_to_v

USE planet_constants_mod, ONLY: kappa, p_zero
USE UM_ParVars          , ONLY: offx, offy
USE interpor_mod        , ONLY: interp_order_linear
USE nlsizes_namelist_mod, ONLY: model_levels, row_length, rows, n_rows, &
                                global_row_length

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumPLevs  ! Number of P levels for output

REAL, INTENT(IN) :: plevs_ws(NumPLevs) 

! u-wind on C-grid model levs
REAL, INTENT(IN) :: UFields(udims_s%i_start:udims_s%i_end, &
                            udims_s%j_start:udims_s%j_end, &
                            udims_s%k_start:udims_s%k_end) 

! v-wind on C-grid model levs
REAL, INTENT(IN) :: VFields(vdims_s%i_start:vdims_s%i_end, &
                            vdims_s%j_start:vdims_s%j_end, &
                            vdims_s%k_start:vdims_s%k_end) 

! u-wind on B-grid p levs
REAL, INTENT(OUT) :: UFieldsP(udims%i_start:udims%i_end, &
                              vdims%j_start:vdims%j_end, NumPLevs) 

! v-wind on B-grid p levs
REAL, INTENT(OUT) :: VFieldsP(udims%i_start:udims%i_end, &
                              vdims%j_start:vdims%j_end, NumPLevs) 

! Local variables:
INTEGER :: k

REAL :: pressure_pa ! Pressure in Pascal
REAL :: pressure_ex ! Exner pressure

! u-wind on C-grid 
REAL :: work_u(udims%i_start:udims%i_end, udims%j_start:udims%j_end, NumPLevs)
! v-wind on C-grid
REAL :: work_v(vdims%i_start:vdims%i_end, vdims%j_start:vdims%j_end, NumPLevs)

! exner at u, v points
REAL :: exner_at_u(udims%i_start:udims%i_end,  &
                   udims%j_start:udims%j_end,  &
                   udims%k_start:udims%k_end)
REAL :: exner_at_v(vdims%i_start:vdims%i_end,  &
                   vdims%j_start:vdims%j_end,  &
                   vdims%k_start:vdims%k_end)

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'UV_PLEVSB'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!  Calculate exner at u points, store in exner_at_u
CALL p_to_u(exner_rho_levels,                              &
            pdims_s%i_start,pdims_s%i_end,                 &
            pdims_s%j_start,pdims_s%j_end,                 &
            udims%i_start,udims%i_end,                     &
            udims%j_start,udims%j_end,                     &
            udims%k_start,udims%k_end,                     &
            exner_at_u)
                    
!  Calculate exner at v points, store in exner_at_v
CALL p_to_v(exner_rho_levels,                              &
            pdims_s%i_start,pdims_s%i_end,                 &
            pdims_s%j_start,pdims_s%j_end,                 &
            vdims%i_start,vdims%i_end,                     &
            vdims%j_start,vdims%j_end,                     &
            vdims%k_start,vdims%k_end,                     &
            exner_at_v)


         
DO  k = 1, NumPLevs

  pressure_pa = PLevs_ws(k)*100.0 ! Convert from hPa to Pa
  pressure_ex = ( pressure_pa /p_zero )**kappa
  
  ! Interpolate u and v to the required p levels
  CALL vert_interp2 (ufields, row_length, rows, model_levels,           &
                     pressure_ex, offx, offy, 0, 0,                     &
                     exner_at_u, interp_order_linear, work_u(:,:,k) )
  CALL vert_interp2 (vfields, row_length, n_rows, model_levels,         &
                     pressure_ex, offx, offy, 0, 0,                     &
                     exner_at_v, interp_order_linear, work_v(:,:,k) )
END DO

! Perform simple horizontal interpolation from 'C' to 'B' grid
CALL  uC_to_uB(work_u, row_length, rows, n_rows, NumPLevs,              &
                                        offx, offy, UFieldsP)  

CALL  vC_to_vB(work_v, rows, row_length, n_rows, NumPLevs,              &
               offx, offy, global_row_length,                           &
               VFieldsP, UFieldsP )


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE uv_plevsB

END MODULE uv_plevsB_mod

