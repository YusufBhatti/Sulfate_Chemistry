! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Control routine to call Leonard term code

MODULE leonard_term_ctl_mod

IMPLICIT NONE

!
! Description:
!   Control routine to call code to calculate fluxes and increments due
!   to the Leonard terms
!
! Method:
!   Performs separate subroutine calls to calculate the Leonard term
!   vertical flux and increment for each conserved variable:
!     u      (zonal wind)
!     v      (meridional wind)
!     w      (vertical wind)
!     thetal (theta obtained if all condensate evaporates)
!     qw     (total water content q + qcl + qcf)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Diffusion and filtering
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName =                    &
                                        'LEONARD_TERM_CTL_MOD'


CONTAINS


SUBROUTINE leonard_term_ctl( etadot, u, v, w, theta, q, qcl, qcf,       &
                           rho_wet_rsq, exner_theta_levels,             &
                           exner_rho_levels,                            &
                           dtrdz_charney_grid, dtrdz_u, dtrdz_v,        &
                           u_inc_leonard, v_inc_leonard, w_inc_leonard, &
                           thetal_inc_leonard, qw_inc_leonard,          &
                           stashwork3 )

USE atm_fields_bounds_mod, ONLY: udims, udims_s, vdims, vdims_s,        &
                                 wdims, wdims_s,                        &
                                 tdims, tdims_s, tdims_l,               &
                                 pdims, pdims_s
USE nlsizes_namelist_mod,  ONLY: bl_levels
USE level_heights_mod,     ONLY: r_theta_levels, r_rho_levels
USE metric_terms_mod,      ONLY: deta_xi3_theta
USE timestep_mod,          ONLY: timestep
USE planet_constants_mod,  ONLY: lcrcp, lsrcp

USE turb_diff_mod,         ONLY: leonard_kl

USE leonard_term_vert_th_mod, ONLY: leonard_term_vert_th
USE leonard_term_vert_u_mod,  ONLY: leonard_term_vert_u
USE leonard_term_vert_v_mod,  ONLY: leonard_term_vert_v
USE p_to_t_mod,            ONLY: p_to_t
USE ereport_mod,           ONLY: ereport
USE errormessagelength_mod,ONLY: errormessagelength


USE stash_array_mod,       ONLY: len_stlist, stindex, stlist,           &
                                 num_stash_levels, stash_levels, si, sf
USE submodel_mod,          ONLY: atmos_im
USE um_parvars,            ONLY: at_extremity
USE um_stashcode_mod,      ONLY: stashcode_bl_sec

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

! Advective vertical velocity, i.e. w relative to the terrain-following
! grid.  This is the w which corresponds to vertical fluxes between
! grid boxes, so we use the horizonotal gradients of this to compute
! the Leonard term fluxes.
REAL, INTENT(IN) :: etadot( wdims_s%i_start:wdims_s%i_end,              &
                            wdims_s%j_start:wdims_s%j_end,              &
                            wdims_s%k_start:wdims_s%k_end )

! Start-of-timestep model winds, corresponding to momentum
REAL, INTENT(IN) :: u ( udims_s%i_start:udims_s%i_end,                  &
                        udims_s%j_start:udims_s%j_end,                  &
                        udims_s%k_start:udims_s%k_end )
REAL, INTENT(IN) :: v ( vdims_s%i_start:vdims_s%i_end,                  &
                        vdims_s%j_start:vdims_s%j_end,                  &
                        vdims_s%k_start:vdims_s%k_end )
REAL, INTENT(IN) :: w ( wdims_s%i_start:wdims_s%i_end,                  &
                        wdims_s%j_start:wdims_s%j_end,                  &
                        wdims_s%k_start:wdims_s%k_end )

! Start-of-timestep theta and moisture variables
REAL, INTENT(IN) :: theta ( tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(IN) :: q     ( tdims_l%i_start:tdims_l%i_end,              &
                            tdims_l%j_start:tdims_l%j_end,              &
                            tdims_l%k_start:tdims_l%k_end)
REAL, INTENT(IN) :: qcl   ( tdims_l%i_start:tdims_l%i_end,              &
                            tdims_l%j_start:tdims_l%j_end,              &
                            tdims_l%k_start:tdims_l%k_end)
REAL, INTENT(IN) :: qcf   ( tdims_l%i_start:tdims_l%i_end,              &
                            tdims_l%j_start:tdims_l%j_end,              &
                            tdims_l%k_start:tdims_l%k_end)

! Density, exner pressure, etc
REAL, INTENT(IN) :: rho_wet_rsq        ( pdims_s%i_start:pdims_s%i_end, &
                                         pdims_s%j_start:pdims_s%j_end, &
                                         pdims_s%k_start:pdims_s%k_end )
REAL, INTENT(IN) :: exner_theta_levels ( tdims_s%i_start:tdims_s%i_end, &
                                         tdims_s%j_start:tdims_s%j_end, &
                                         tdims_s%k_start:tdims_s%k_end )
REAL, INTENT(IN) :: exner_rho_levels   ( pdims_s%i_start:pdims_s%i_end, &
                                         pdims_s%j_start:pdims_s%j_end, &
                                         pdims_s%k_start:pdims_s%k_end + 1)

! Arrays of timestep / ( dz * rho * r^2 ), precalculated in the BL scheme
REAL, INTENT(IN) :: dtrdz_charney_grid ( pdims%i_start:pdims%i_end,     &
                                         pdims%j_start:pdims%j_end,     &
                                         bl_levels )
REAL, INTENT(IN) :: dtrdz_u            ( udims%i_start:udims%i_end,     &
                                         udims%j_start:udims%j_end,     &
                                         bl_levels )
REAL, INTENT(IN) :: dtrdz_v            ( vdims%i_start:vdims%i_end,     &
                                         vdims%j_start:vdims%j_end,     &
                                         bl_levels )

! Leonard term increments to conserved variables:
!                    Zonal wind
REAL, INTENT(OUT) :: u_inc_leonard      ( udims%i_start:udims%i_end,    &
                                          udims%j_start:udims%j_end,    &
                                          bl_levels )
!                    Meridional wind
REAL, INTENT(OUT) :: v_inc_leonard      ( vdims%i_start:vdims%i_end,    &
                                          vdims%j_start:vdims%j_end,    &
                                          bl_levels )
!                    Vertical wind
REAL, INTENT(OUT) :: w_inc_leonard      ( wdims%i_start:wdims%i_end,    &
                                          wdims%j_start:wdims%j_end,    &
                                          bl_levels )
!                    Liquid + ice water potential temperature
REAL, INTENT(OUT) :: thetal_inc_leonard ( tdims%i_start:tdims%i_end,    &
                                          tdims%j_start:tdims%j_end,    &
                                          bl_levels )
!                    Total water content q + qcl + qcf
REAL, INTENT(OUT) :: qw_inc_leonard     ( tdims%i_start:tdims%i_end,    &
                                          tdims%j_start:tdims%j_end,    &
                                          bl_levels )

! STASH workspace
REAL, INTENT(INOUT) :: stashwork3(*)


! Local arrays for storing derived fields:

! etadot, converted to ms-1, and with only small-size halos.
REAL :: etadot_mps     ( wdims_s%i_start:wdims_s%i_end,                 &
                         wdims_s%j_start:wdims_s%j_end,                 &
                         bl_levels )
! Conserved thermodynamic variables  (theta_l, then qw)
REAL :: field_n        ( tdims_s%i_start:tdims_s%i_end,                 &
                         tdims_s%j_start:tdims_s%j_end,                 &
                         bl_levels )
! Wet (total) density * r^2 on theta-levels, with halos
REAL :: rho_wet_rsq_tq ( tdims_s%i_start:tdims_s%i_end,                 &
                         tdims_s%j_start:tdims_s%j_end,                 &
                         bl_levels )
! 3-D field of Leonard term parameter on rho-levels, accounting for
! stability limiting which varies depending on etadot
REAL :: kl             ( tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end,                     &
                         1:bl_levels )


! Loop counters
INTEGER :: i, j, k

! STASH item numbers for the increment and flux diagnostics to be output
! inside calls to leonard_term_vert_th for different fields
INTEGER :: item_inc
INTEGER :: item_flx

! Stuff for STASH calls
INTEGER :: icode, item
INTEGER, PARAMETER :: im_index = 1
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LEONARD_TERM_CTL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! Interpolate rho_wet_rsq onto theta-levels.
! While there is already a copy of rho_wet interpolated onto theta-levels
! in atmos_physics2 (rho_wet_tq), that doesn't have halos, and we need
! halos for interpolation onto u, v points
CALL p_to_t( tdims%i_len, tdims%j_len,                                  &
             tdims_l%halo_i, tdims_l%halo_j,                            &
             tdims_s%halo_i, tdims_s%halo_j, bl_levels,                 &
             r_theta_levels(:,:,0:bl_levels+1),                         &
             r_rho_levels(:,:,1:bl_levels+1),                           &
             rho_wet_rsq(:,:,1:bl_levels+1),                            &
             rho_wet_rsq_tq  )
! Note: r_theta_levels and r_rho_levels are set in setcona_4a with
! large halo size (tdims_l%halo), while the fields we want to interp have
! the small halo size (tdims_s%halo)


! Convert etadot to ms-1
DO k = 1, bl_levels
  DO j = wdims_s%j_start, wdims_s%j_end
    DO i = wdims_s%i_start, wdims_s%i_end
      etadot_mps(i,j,k) = etadot(i,j,k) * deta_xi3_theta(i,j,k)
    END DO
  END DO
END DO


! ------------------------------------------------------------------
! Calculate kl at rho-points, accounting for stability limit
! ------------------------------------------------------------------

! Set Kl to zero at level 1 as it is only used between 
! levels 2:bl_levels
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    kl(i,j,1) = 0.0
  END DO
END DO

DO k = 2, bl_levels
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end

      ! Leonard term parameter is the min of the input leonard_kl
      ! and the max stable value   6 dz / ( dt dw )

      ! For dw we use the maximum horizontal finite difference that
      ! contributes to the flux at each rho-point.

      kl(i,j,k) = MIN( leonard_kl,                                      &
          6.0 * ( r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1) )     &
              / ( timestep * MAX(                                       &
                    ! Differences on theta-level above:
                    ABS( etadot_mps(i+1,j,k)   - etadot_mps(i,j,k)   ), &
                    ABS( etadot_mps(i-1,j,k)   - etadot_mps(i,j,k)   ), &
                    ABS( etadot_mps(i,j+1,k)   - etadot_mps(i,j,k)   ), &
                    ABS( etadot_mps(i,j-1,k)   - etadot_mps(i,j,k)   ), &
                    ! Differences on theta-level below:
                    ABS( etadot_mps(i+1,j,k-1) - etadot_mps(i,j,k-1) ), &
                    ABS( etadot_mps(i-1,j,k-1) - etadot_mps(i,j,k-1) ), &
                    ABS( etadot_mps(i,j+1,k-1) - etadot_mps(i,j,k-1) ), &
                    ABS( etadot_mps(i,j-1,k-1) - etadot_mps(i,j,k-1) )  &
                                )                                       &
                )    )
    END DO
  END DO
END DO


! ------------------------------------------------------------------
! Flux of zonal wind u
! ------------------------------------------------------------------

CALL leonard_term_vert_u( u, etadot_mps, rho_wet_rsq_tq,                &
                          dtrdz_u, u_inc_leonard, stashwork3 )

! ------------------------------------------------------------------
! Flux of meridional wind v
! ------------------------------------------------------------------

CALL leonard_term_vert_v( v, etadot_mps, rho_wet_rsq_tq,                &
                          dtrdz_v, v_inc_leonard, stashwork3 )

! ------------------------------------------------------------------
! Flux of vertical wind w
! ------------------------------------------------------------------

! Set item numbers for STASH calls inside the generic subroutine call
item_inc = 199
item_flx = 555
CALL leonard_term_vert_th( w(:,:,1:bl_levels), etadot_mps, kl,          &
                           rho_wet_rsq, dtrdz_charney_grid,             &
                           w_inc_leonard, exner_theta_levels,           &
                           exner_rho_levels,                            &
                           item_inc, item_flx, stashwork3 )

! ------------------------------------------------------------------
! Flux of liquid + ice water potential temperature
! ------------------------------------------------------------------

! Convert thermodynamic variables to theta_l for flux calculation
DO k = 1, bl_levels
  DO j = tdims_s%j_start, tdims_s%j_end
    DO i = tdims_s%i_start, tdims_s%i_end
      ! Store start-of-timestep theta_l in field_n
      field_n(i,j,k) = theta(i,j,k)                                     &
                     - ( lcrcp * qcl(i,j,k) + lsrcp * qcf(i,j,k) )      &
                       / exner_theta_levels(i,j,k)
    END DO
  END DO
END DO

! Set item numbers for STASH calls inside the generic subroutine call
item_inc = 197
item_flx = 556
CALL leonard_term_vert_th( field_n, etadot_mps, kl,                     &
                           rho_wet_rsq, dtrdz_charney_grid,             &
                           thetal_inc_leonard, exner_theta_levels,      &
                           exner_rho_levels,                            &
                           item_inc, item_flx, stashwork3 )

! ------------------------------------------------------------------
! Flux of total water content
! ------------------------------------------------------------------

! Convert moisture variables to total water
DO k = 1, bl_levels
  DO j = tdims_s%j_start, tdims_s%j_end
    DO i = tdims_s%i_start, tdims_s%i_end
      ! Store total water in field_n
      field_n(i,j,k) = q(i,j,k) + qcl(i,j,k) + qcf(i,j,k)
    END DO
  END DO
END DO

! Set item numbers for STASH calls inside the generic subroutine call
item_inc = 198
item_flx = 557
CALL leonard_term_vert_th( field_n, etadot_mps, kl,                     &
                           rho_wet_rsq,  dtrdz_charney_grid,            &
                           qw_inc_leonard, exner_theta_levels,          &
                           exner_rho_levels,                            &
                           item_inc, item_flx, stashwork3 )


! ------------------------------------------------------------------
! Output diagnostic of kl to STASH
! ------------------------------------------------------------------

icode = 0
item = 552
IF (icode <= 0 .AND. sf(item,stashcode_bl_sec)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d( stashwork3(si(item,stashcode_bl_sec,im_index)),     &
                    kl,                                                 &
                    tdims%i_end, tdims%j_end, bl_levels,                &
                    0, 0, 0, 0, at_extremity,                           &
                    stlist(1,stindex(1,item,stashcode_bl_sec,im_index)),&
                    len_stlist, stash_levels, num_stash_levels+1,       &
                    atmos_im, stashcode_bl_sec, item,                   &
                    icode, cmessage )
  IF (icode  >   0) THEN
    cmessage=": error in copydiag_3d(item)"//TRIM(cmessage)
    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE leonard_term_ctl

END MODULE leonard_term_ctl_mod
