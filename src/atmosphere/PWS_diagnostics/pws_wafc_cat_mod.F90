! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine pws_wafc_cat_mod ---------------------------------------
!
!   Purpose: Calculates PWS WAFC Cat potential, stash section 20
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics
MODULE pws_wafc_cat_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_WAFC_CAT_MOD'

CONTAINS

SUBROUTINE pws_wafc_cat()

! Description
! Calculate WAFC CAT potential turbulence predictor
!

! Method:
!  1) Calculate shear turbulence predictor (TI1)
!     i) Calculate u and v winds on requested pressure levels (B grid winds)
!     ii) Calculate model level u and v winds on B grid
!     iii) Calculate pressure on rho levels at uv pts. B grid.
!     iv) calculate Ellrod's first CAT indicator
!  2) Calculate mountain wave part up to T+24
!     After that just duplicate the shear turbulence fields
!     i) need u and v comps of GWD stress on b grid.
!     ii) Calculate pressure on theta levels at uv pts. B grid.
!     iii) for pressure levels <=300hPa  calc MW predictor.
!  3) Combine CAT and MW fields by taking the max of the fields at each grid pt
!     for each level where both CAT and MW turbulence are requested 
!     and is <=300hPa

USE atm_fields_bounds_mod, ONLY: udims, vdims, pdims, udims_s, vdims_s, &
                                 pdims_s, tdims, tdims_s
USE level_heights_mod,     ONLY: r_rho_levels
USE missing_data_mod,      ONLY: rmdi
USE nlsizes_namelist_mod,  ONLY: row_length, rows, model_levels, &
                                 global_row_length

USE model_time_mod, ONLY: forecast_hrs

USE atm_fields_mod

USE pws_diags_mod, ONLY: wafc_cat_press_levs, wafc_cat_press,        &
                          pws_gwd_stress_lev_u,                      &
                          pws_gwd_stress_lev_v, pws_wafc_cat_diag

USE field_types, ONLY: fld_type_p, fld_type_u, fld_type_v

USE pws_wafc_ellrod1_mod, ONLY: pws_wafc_ellrod1
USE pws_mtn_stress_pref_mod,  ONLY: pws_mtn_stress_pref

USE uc_to_ub_mod, ONLY:  uc_to_ub
USE vc_to_vb_mod, ONLY:  vc_to_vb
USE pc_to_pb_mod, ONLY:  pc_to_pb

USE diff_mod, ONLY: diffp, diffx, diffy

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE uv_plevsB_mod, ONLY: uv_plevsB

IMPLICIT NONE


! Local variables
REAL  :: uwind_plev(udims%i_start:udims%i_end,                        &
                    vdims%j_start:vdims%j_end, wafc_cat_press_levs)

REAL  :: vwind_plev(udims%i_start:udims%i_end,                        &
                    vdims%j_start:vdims%j_end, wafc_cat_press_levs)

REAL :: local_ub(udims%i_start:udims%i_end,                           &
                 vdims%j_start:vdims%j_end,                           &
                 udims%k_start:udims%k_end)

REAL :: local_vb(udims%i_start:udims%i_end,                           &
                 vdims%j_start:vdims%j_end,                           &
                 vdims%k_start:vdims%k_end)

REAL :: local_pb(udims%i_start:udims%i_end,                           &
                 vdims%j_start:vdims%j_end,                           &
                 pdims%k_start:pdims%k_end)

REAL :: local_r_rho_levels_b(udims%i_start:udims%i_end,               &
                 vdims%j_start:vdims%j_end,                           &
                 pdims%k_start:pdims%k_end)

REAL :: local_p_theta_b(udims%i_start:udims%i_end,                    &
                 vdims%j_start:vdims%j_end,                           &
                 pdims%k_start:pdims%k_end)


REAL :: dzdp(udims%i_start:udims%i_end,                               &
             vdims%j_start:vdims%j_end)

REAL :: dudp(udims%i_start:udims%i_end,                               &
             vdims%j_start:vdims%j_end)

REAL :: dvdp(udims%i_start:udims%i_end,                               &
             vdims%j_start:vdims%j_end)

REAL :: dwdz(udims%i_start:udims%i_end,                               &
             vdims%j_start:vdims%j_end)

REAL :: dudx(udims%i_start:udims%i_end,                               &
             vdims%j_start:vdims%j_end)

REAL :: dvdx(udims%i_start:udims%i_end,                               &
             vdims%j_start:vdims%j_end)

REAL :: dudy(udims%i_start:udims%i_end,                               &
             vdims%j_start:vdims%j_end)

REAL :: dvdy(udims%i_start:udims%i_end,                               &
             vdims%j_start:vdims%j_end)

REAL  :: catpred(udims%i_start:udims%i_end,                           &
                 vdims%j_start:vdims%j_end, wafc_cat_press_levs)

REAL  :: gwdub(udims%i_start:udims%i_end,                             &
               vdims%j_start:vdims%j_end, pdims%k_len)

REAL  :: gwdvb(udims%i_start:udims%i_end,                             &
               vdims%j_start:vdims%j_end, pdims%k_len)

REAL  :: local_ptheta_b(udims%i_start:udims%i_end,                    &
                        vdims%j_start:vdims%j_end, pdims%k_len)

REAL  :: mwtpred(udims%i_start:udims%i_end,                           &
                 vdims%j_start:vdims%j_end, wafc_cat_press_levs)


INTEGER :: i,j,k, icode

CHARACTER (LEN=*), PARAMETER :: RoutineName='PWS_WAFC_CAT'
! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


icode=0
mwtpred(:,:,:) = rmdi ! initialise for all levels.

! ----------------------
!  1.i calculate winds u and v on requested pressure levels on b grid
! ----------------------

! Interp u, v on model levels C grid to p levels B grid
CALL uv_plevsB(wafc_cat_press_levs, wafc_cat_press,      &
               u, v, uwind_plev, vwind_plev)

! ----------------------
!  1.ii Calculate model level u and v winds on B grid
! ----------------------

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


! ----------------------
! 1. iii Calculate pressure on rho levels at uv pts. B grid.
! ----------------------

CALL  pc_to_pb(                                                   &
       r_rho_levels(1:row_length,                                 &
                    1:rows,                                       &
                    pdims%k_start:pdims%k_end),                   &
  pdims%i_len,pdims%j_len,vdims%j_len,pdims%k_len,                &
  pdims_s%halo_i,pdims_s%halo_j,                                  &
  local_r_rho_levels_b)

CALL  pc_to_pb(                                                   &
       p(pdims%i_start:pdims%i_end,                               &
         pdims%j_start:pdims%j_end,                               &
         pdims%k_start:pdims%k_end),                              &
  pdims%i_len,pdims%j_len,vdims%j_len,pdims%k_len,                &
  pdims_s%halo_i,pdims_s%halo_j,                                  &
  local_pb)

! ----------------------
! 1.iv) calculate Ellrod's first CAT indicator
! ----------------------


!-----------------------------------------------------------------------
! Calc. vertical derivatives of z, u & v wrt pressure using cubic spline
! Note that spline values must be monotonically increasing.
! Then calculate vertical wind shear dWdz

DO k=1,wafc_cat_press_levs

  ! next multiply by 100 to convert to pascals as stash uses hPa.
  wafc_cat_press(k) = wafc_cat_press(k)*100.0

  CALL DiffP(pdims%k_len, wafc_cat_press(k),                    &
             udims%i_start,udims%i_end,                         &
             vdims%j_start,vdims%j_end,                         &
             local_r_rho_levels_b, local_pb, dZdP)

  CALL DiffP(udims%k_len, wafc_cat_press(k),                    &
             udims%i_start,udims%i_end,                         &
             vdims%j_start,vdims%j_end,                         &
             local_ub, local_pb, dUdP)

  CALL DiffP(vdims%k_len, wafc_cat_press(k),                    &
             udims%i_start,udims%i_end,                         &
             vdims%j_start,vdims%j_end,                         &
             local_vb, local_pb, dVdP)

  dWdZ(:,:) = SQRT(dUdP(:,:)**2 + dVdP(:,:)**2)/ ABS(dZdP(:,:))


  !-----------------------------------------------------------------------
  ! Calculate horizontal derivatives
  ! latitudinal derivatives
  CALL DiffX( uwind_plev(:,:,k), dUdX, fld_type_u, icode )
  CALL DiffX( vwind_plev(:,:,k), dVdX, fld_type_v, icode)
  ! longitudinal derivatives
  CALL DiffY( uwind_plev(:,:,k), dUdY, fld_type_u, icode )
  CALL DiffY( vwind_plev(:,:,k), dVdY, fld_type_v, icode )

  !-----------------------------------------------------------------------
  ! Calculate CAT predictor (Ellrod's index TI1)

  CALL pws_wafc_ellrod1(uwind_plev(:,:,k),          &
                    vwind_plev(:,:,k),              &
                    dUdX,dUdY,dVdX,dVdY,dWdZ,       &
                    catpred(:,:,k) )

END DO


! ----------------------
! 2) Calculate mountain wave part up to T+24
!    i) need u and v comps of GWD stress on B grid.
! ----------------------

IF (forecast_hrs <= 24) THEN

  CALL  uC_to_uB(pws_gwd_stress_lev_u(:,:,1:),row_length,rows,vdims%j_len, &
                 pdims%k_len,pdims_s%halo_i,pdims_s%halo_j,                &
                 gwduB(:,:,:))


  CALL  vC_to_vB(pws_gwd_stress_lev_v(:,:,1:),rows,row_length,vdims%j_len, &
            pdims%k_len,pdims_s%halo_i,pdims_s%halo_j,                     &
            global_row_length,gwdvB(:,:,:),gwduB(:,:,:))

  ! ----------------------
  !   ii) Calculate pressure on theta levels at uv pts. B grid. 
  ! ----------------------

  CALL  pc_to_pb(                                                      &
         p_theta_levels(tdims%i_start:tdims%i_end,                     &
                        tdims%j_start:tdims%j_end,                     &
                        1:tdims%k_end),                                &
    tdims%i_len,tdims%j_len,vdims%j_len,pdims%k_len,                   &
    tdims_s%halo_i,tdims_s%halo_j,                                     &
    local_ptheta_b)

  ! ----------------------
  !     iii) for pressure levels <=300hPa  calc MW predictor
  ! ----------------------

  DO k=1,wafc_cat_press_levs
  
    IF (wafc_cat_press(k) <= 30000.0 ) THEN
      CALL pws_mtn_stress_pref(wafc_cat_press(k),                     &
                         uwind_plev(:,:,k), vwind_plev(:,:,k),        & 
                         gwduB, gwdvB, local_ptheta_b, mwtpred(:,:,k))
    END IF
  END DO

END IF
! -----------------------------------------------------------------------------
!  3) Combine CAT and MW fields by taking the max of the fields at each grid pt
!     for each level where both CAT and MW turbulence are requested 
!     and is <=300hPa
! -----------------------------------------------------------------------------

DO k=1,wafc_cat_press_levs
  DO j = vdims%j_start,vdims%j_end
    DO i = udims%i_start,udims%i_end
      IF (catpred(i,j,k) /= rmdi ) THEN
        pws_wafc_cat_diag(i,j,k) = MAX( mwtpred(i,j,k), catpred(i,j,k) )
      END IF
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE pws_wafc_cat

!--------------------------------------------------------------------------

END MODULE pws_wafc_cat_mod
