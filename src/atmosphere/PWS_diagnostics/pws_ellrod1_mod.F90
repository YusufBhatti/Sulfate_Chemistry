! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Calculates Ellrod_ti1 turbulence index, stash section 20
!
! Programming standard; Unified Model Documentation Paper No. 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics

MODULE pws_ellrod1_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='PWS_ELLROD1_MOD'

CONTAINS

SUBROUTINE pws_ellrod1(u, v, p_theta_levels, p)

! Description:
!   Routine to calculate fields of Ellrod's first CAT indicator TI1 
!  (described in Ellrod and Knapp 1992), using model level u and v and a
!   finite difference method to obtain vertical derivatives. 
!   Units of TI1 are s^-2. 
!   Ref: Ellrod and Knapp (1992), 'An Objective Clear Air Turbulence 
!        Forecasting Technique: Verification and Operational Use', 
!        Weather Forecast., 7, pp 150-165.


USE yomhook,                    ONLY: lhook, dr_hook
USE parkind1,                   ONLY: jprb, jpim

USE pws_diags_mod,              ONLY: ellrodt1_turb_press_levels,            &
                                      ellrodt1_turb_press,                   &
                                      pws_ellrodt1_turb

USE u_to_p_mod,            ONLY: u_to_p
USE v_to_p_mod,            ONLY: v_to_p

USE uc_to_ub_mod, ONLY:  uc_to_ub
USE vc_to_vb_mod, ONLY:  vc_to_vb
USE pc_to_pb_mod, ONLY: pc_to_pb

USE field_types, ONLY: fld_type_u, fld_type_v
USE diff_mod, ONLY: diffx, diffy

USE atm_fields_bounds_mod, ONLY: udims, udims_s, vdims, vdims_s,             &
                                 pdims, pdims_s, tdims_s
USE level_heights_mod,     ONLY: r_rho_levels
USE model_domain_mod,           ONLY: l_regular
USE ereport_mod,                ONLY: ereport
USE missing_data_mod,           ONLY: rmdi, imdi
USE nlsizes_namelist_mod,       ONLY: model_levels,                          &
                                      global_row_length
USE pws_vert_wind_shear_sq_mod, ONLY: pws_vert_wind_shear_sq

IMPLICIT NONE

REAL, INTENT(IN)  :: u(udims_s%i_start:udims_s%i_end,                        &
                       udims_s%j_start:udims_s%j_end,                        &
                       udims_s%k_start:udims_s%k_end)

REAL, INTENT(IN)  :: v(vdims_s%i_start:vdims_s%i_end,                        &
                       vdims_s%j_start:vdims_s%j_end,                        &
                       vdims_s%k_start:vdims_s%k_end)


REAL, INTENT(IN)  :: p_theta_levels(tdims_s%i_start:tdims_s%i_end,           &
                                    tdims_s%j_start:tdims_s%j_end,           &
                                    tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(IN)  :: p(pdims_s%i_start:pdims_s%i_end,                        &
                      pdims_s%j_start:pdims_s%j_end,                         &
                      pdims_s%k_start:pdims_s%k_end)


! Local variables

REAL :: local_ub(udims%i_start:udims%i_end,                           &
                 vdims%j_start:vdims%j_end,                           &
                 udims%k_start:udims%k_end)   ! uwind on b grid

REAL :: local_vb(udims%i_start:udims%i_end,                           &
                 vdims%j_start:vdims%j_end,                           &
                 vdims%k_start:vdims%k_end)   ! vwind on b grid

REAL :: local_pb(udims%i_start:udims%i_end,                           &
                 vdims%j_start:vdims%j_end,                           &
                 pdims%k_start:pdims%k_end)   ! press on b grid

! horizontal diffs on model levels

REAL :: dudx(udims%i_start:udims%i_end,                               &
             vdims%j_start:vdims%j_end,                               &
                  pdims%k_start:pdims%k_end)  ! dudx on b grid

REAL :: dvdx(udims%i_start:udims%i_end,                               &
             vdims%j_start:vdims%j_end,                               &
                  pdims%k_start:pdims%k_end)  ! dvdx on b grid

REAL :: dudy(udims%i_start:udims%i_end,                               &
             vdims%j_start:vdims%j_end,                               &
                  pdims%k_start:pdims%k_end)  ! dudy on b grid

REAL :: dvdy(udims%i_start:udims%i_end,                               &
             vdims%j_start:vdims%j_end,                               &
                  pdims%k_start:pdims%k_end)  ! dvdy on b grid

! local horizontal diffs intepolated to requested pressure levels

REAL :: dudx_std

REAL :: dudy_std

REAL :: dvdx_std

REAL :: dvdy_std

! local variables required to calculate wind shear.

REAL :: r_rho_levels_uv(udims%i_start:udims%i_end,                    &
                          vdims%j_start:vdims%j_end,                  &
                          pdims%k_start:pdims%k_end)


REAL :: dudz(udims%i_start:udims%i_end,                               &
                          vdims%j_start:vdims%j_end,                  &
                          pdims%k_start:pdims%k_end)

REAL :: dvdz(udims%i_start:udims%i_end,                               &
                          vdims%j_start:vdims%j_end,                  &
                          pdims%k_start:pdims%k_end)

REAL :: shear_sq(udims%i_start:udims%i_end,                           &
                          vdims%j_start:vdims%j_end,                  &
                          pdims%k_start:pdims%k_end)

REAL :: shear_std

! Model level corresponding to pressure level we're interested in
INTEGER :: level

! Fraction of the way we are through the layer
REAL :: alpha

REAL :: st_def, sh_def, def, shearv

INTEGER :: i,j,k,l ! Loop counters
INTEGER :: ErrorStatus

CHARACTER (LEN=*), PARAMETER :: RoutineName='PWS_ELLROD1'

! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. l_regular) THEN
  ErrorStatus=100
  CALL EReport( RoutineName, ErrorStatus,                                    &
                "Cannot calculate Ellrod1 turbulence - " //                  &
                "var res and/or rotated grids are not supported" )
END IF

! interpolate u and v winds onto b grid equivalents.
CALL  uC_to_uB(                                                     &
      u(udims%i_start:udims%i_end,                                  &
        udims%j_start:udims%j_end,                                  &
        udims%k_start:udims%k_end),                                 &
  udims%i_len,udims%j_len,vdims%j_len,udims%k_len,                  &
  udims_s%halo_i,udims_s%halo_j,                                    &
  local_ub)


CALL  vC_to_vB(                                                   &
     v(vdims%i_start:vdims%i_end,                                 &
       vdims%j_start:vdims%j_end,                                 &
       vdims%k_start:vdims%k_end),                                &
pdims%j_len,pdims%i_len,vdims%j_len,vdims%k_len,                  &
vdims_s%halo_i,vdims_s%halo_j,                                    &
global_row_length,local_vb,local_ub)


! Calculate horizontal derivatives of wind components u and v

! For each model level, find horizontal derivatives
DO k = 1, model_levels


  ! Calculate latitudinal derivatives
  CALL DiffX( local_ub(:,:,k), dUdX(:,:,k), fld_type_u, ErrorStatus)  
  CALL DiffX( local_vb(:,:,k), dVdX(:,:,k), fld_type_v, ErrorStatus)  


  ! longitudinal derivatives
  CALL DiffY( local_ub(:,:,k), dUdY(:,:,k), fld_type_u, ErrorStatus)  
  CALL DiffY( local_vb(:,:,k), dVdY(:,:,k), fld_type_v, ErrorStatus)  



END DO

! get pressure and heights of rho levels on b grids

CALL  pc_to_pb(                                                     &
       r_rho_levels(pdims%i_start:pdims%i_end,                      &
         pdims%j_start:pdims%j_end,                                 &
         pdims%k_start:pdims%k_end),                                &
  pdims%i_len,pdims%j_len,vdims%j_len,pdims%k_len,                  &
  pdims_s%halo_i,pdims_s%halo_j,                                    &
  r_rho_levels_uv)

CALL  pc_to_pb(                                                     &
         p(pdims%i_start:pdims%i_end,                               &
         pdims%j_start:pdims%j_end,                                 &
         pdims%k_start:pdims%k_end),                                &
  pdims%i_len,pdims%j_len,vdims%j_len,pdims%k_len,                  &
  pdims_s%halo_i,pdims_s%halo_j,                                    &
  local_pb)

! calculate wind shear
DO k = 1, model_levels-1
  DO j = vdims%j_start,vdims%j_end
    DO i = udims%i_start,udims%i_end
      dudz(i,j,k) = (local_ub(i,j,k+1) - local_ub(i,j,k)) /         &
             (r_rho_levels_uv(i,j,k+1) - r_rho_levels_uv(i,j,k))
      dvdz(i,j,k) = (local_vb(i,j,k+1) - local_vb(i,j,k)) /         &
             (r_rho_levels_uv(i,j,k+1) - r_rho_levels_uv(i,j,k))
      shear_sq(i,j,k)=dudz(i,j,k)**2 + dvdz(i,j,k)**2
    END DO
  END DO
END DO

! For each standard (pressure) level gridpoint, determine the model level below 
! and above it using a divide and conquer algorithm. Then calculate du/dx, 
! du/dy, dv/dx and dv/dy on the standard level by linear interpolation of the 
! corresponding derivatives on the model level above and below the point. 

! Initialise
pws_ellrodt1_turb(:,:,:) = rmdi

! Interpolate wind shear squared onto requested pressure levels
DO l = 1, ellrodt1_turb_press_levels

  DO j = vdims%j_start, vdims%j_end
    DO i = udims%i_start, udims%i_end
      level = imdi
vert: DO k = pdims%k_start, pdims%k_end - 1
        IF ( (ellrodt1_turb_press(l) >  local_pb(i,j,k)  .AND.            &
             ellrodt1_turb_press(l) <= local_pb(i,j,k+1) ) .OR.           &
             (ellrodt1_turb_press(l) <  local_pb(i,j,k)  .AND.            &
             ellrodt1_turb_press(l) >= local_pb(i,j,k+1) ) ) THEN
              level = k
              alpha = (ellrodt1_turb_press(l) - local_pb(i,j,level))      &
                  /  (local_pb(i,j,level+1) - local_pb(i,j,level))
              EXIT vert
        END IF
      END DO vert

      IF ( (level /= imdi) .AND.                 &
           (dudx(i,j,level) /= rmdi) .AND.       & 
           (dudy(i,j,level) /= rmdi) .AND.       & 
           (dvdx(i,j,level) /= rmdi) .AND.       &
           (dvdy(i,j,level) /= rmdi) .AND.       &
           (dudx(i,j,level+1) /= rmdi) .AND.     & 
           (dudy(i,j,level+1) /= rmdi) .AND.     & 
           (dvdx(i,j,level+1) /= rmdi) .AND.     & 
           (dvdy(i,j,level+1) /= rmdi) )  THEN 

        dudx_std =                                                         & 
                dudx(i,j,level) + alpha *                                  &
               (dudx(i,j,level+1)  - dudx(i,j,level))
        dvdx_std =                                                         & 
                dvdx(i,j,level) + alpha *                                  &
               (dvdx(i,j,level+1)  - dvdx(i,j,level))
        dudy_std =                                                         & 
                dudy(i,j,level) + alpha *                                  &
               (dudy(i,j,level+1)  - dudy(i,j,level))
        dvdy_std =                                                         & 
                dvdy(i,j,level) + alpha *                                  &
               (dvdy(i,j,level+1)  - dvdy(i,j,level))
        shear_std =                                                        & 
                shear_sq(i,j,level) + alpha *                              &
               (shear_sq(i,j,level+1)  - shear_sq(i,j,level))

        St_Def = dUdX_std - dVdY_std 
        Sh_Def = dVdX_std + dUdY_std 
        Def = SQRT( (St_Def**2) + (Sh_Def**2) )
        ShearV = SQRT (shear_sq(i,j,level))
        pws_ellrodt1_turb(i,j,l) = ShearV * Def

      ELSE   ! leave ellrod as rmdi
        CYCLE
 
      END IF

    END DO
  END DO

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE pws_ellrod1


END MODULE pws_ellrod1_mod
