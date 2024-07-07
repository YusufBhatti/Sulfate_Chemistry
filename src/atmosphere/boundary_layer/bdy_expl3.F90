! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: Calculate explicit fluxes taux and tauy on their
!           native (u or v) grids

!  Programming standard : UMDP 3

!  Documentation: UMDP 24.

!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
MODULE bdy_expl3_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'BDY_EXPL3_MOD'
CONTAINS

SUBROUTINE bdy_expl3 (                                                  &
! IN grid related variables
 bl_levels, l_calc_at_p,                                                &
 udimsi, vdimsi, udimsi_s, vdimsi_s, pdimsi,                            &
 udimso, vdimso,                                                        &
! IN SCM diags
 nSCMDpkgs,L_SCMDiags,                                                  &
! IN variables used in flux calculations
 u, v, u_0, v_0, rhokm_u_land, rhokm_v_land, flandfac_u, flandfac_v,    &
 rhokm_u_ssi, rhokm_v_ssi, fseafac_u, fseafac_v, flandg_u, flandg_v,    &
 zhnl, rdz_u, rdz_v, rhokm_u, rhokm_v, taux_fd_u, tauy_fd_v,            &
 rhogamu_u, rhogamv_v, f_ngstress_u, f_ngstress_v,                      &
! OUT explicit momentum fluxes
 taux_land, tauy_land, taux_ssi, tauy_ssi, taux, tauy                   &
 )

USE atm_fields_bounds_mod, ONLY: array_dims
USE bl_option_mod, ONLY: i_bl_vn, i_bl_vn_1a, on

USE jules_sea_seaice_mod, ONLY: l_ctile, buddy_sea
USE jules_surface_mod, ONLY: formdrag, explicit_stress
USE s_scmop_mod,       ONLY: default_streams,                           &
                             t_avg, d_bl, scmdiag_bl
USE scmoutput_mod,     ONLY: scmoutput
USE model_domain_mod,  ONLY: model_type, mt_single_column
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE ex_flux_uv_mod, ONLY: ex_flux_uv
USE mym_ex_flux_uv_mod, ONLY: mym_ex_flux_uv
IMPLICIT NONE

!  Inputs :-
INTEGER, INTENT(IN) ::                                                  &
 bl_levels
                 ! IN Max. no. of "boundary" levels
LOGICAL, INTENT(IN) :: l_calc_at_p
                 ! Flag for separate call to compute fluxes on p-grid,
                 ! if they are needed by the convection scheme
TYPE (array_dims), INTENT(IN) :: udimsi, vdimsi, udimsi_s, vdimsi_s, pdimsi
                 ! Dimensions of the inputs.  If calling
                 ! this routine on the p-grid, pdims or pdims_s
                 ! will be passed in in-place of all of these.
TYPE (array_dims), INTENT(IN) :: udimso, vdimso
                 ! Dimensions of the outputs.

! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN) ::                                                  &
 nSCMDpkgs             ! No of SCM diagnostics packages

LOGICAL, INTENT(IN) ::                                                  &
 L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

REAL, INTENT(IN) ::                                                     &
 u(udimsi_s%i_start:udimsi_s%i_end,udimsi_s%j_start:udimsi_s%j_end,     &
   bl_levels),                                                          &
 v(vdimsi_s%i_start:vdimsi_s%i_end,vdimsi_s%j_start:vdimsi_s%j_end,     &
   bl_levels),                                                          &
                 ! horizontal winds
 u_0(udimsi%i_start:udimsi%i_end,udimsi%j_start:udimsi%j_end),          &
 v_0(vdimsi%i_start:vdimsi%i_end,vdimsi%j_start:vdimsi%j_end),          &
                 ! surface currents
 rhokm_u_land(udimsi%i_start:udimsi%i_end,udimsi%j_start:udimsi%j_end), &
 rhokm_u_ssi(udimsi%i_start:udimsi%i_end,udimsi%j_start:udimsi%j_end),  &
 rhokm_v_land(vdimsi%i_start:vdimsi%i_end,vdimsi%j_start:vdimsi%j_end), &
 rhokm_v_ssi(vdimsi%i_start:vdimsi%i_end,vdimsi%j_start:vdimsi%j_end),  &
 rhokm_u(udimsi%i_start:udimsi%i_end,udimsi%j_start:udimsi%j_end,       &
         bl_levels),                                                    &
 rhokm_v(vdimsi%i_start:vdimsi%i_end,vdimsi%j_start:vdimsi%j_end,       &
         bl_levels),                                                    &
                 ! rho * Km terms on u/v grids
 taux_fd_u(udimsi%i_start:udimsi%i_end,udimsi%j_start:udimsi%j_end,     &
           bl_levels),                                                  &
 tauy_fd_v(vdimsi%i_start:vdimsi%i_end,vdimsi%j_start:vdimsi%j_end,     &
           bl_levels),                                                  &
                 ! explicit form drag
 flandfac_u(udimsi%i_start:udimsi%i_end,udimsi%j_start:udimsi%j_end),   &
 flandfac_v(vdimsi%i_start:vdimsi%i_end,vdimsi%j_start:vdimsi%j_end),   &
 fseafac_u(udimsi%i_start:udimsi%i_end,udimsi%j_start:udimsi%j_end),    &
 fseafac_v(vdimsi%i_start:vdimsi%i_end,vdimsi%j_start:vdimsi%j_end),    &
                 ! land and sea scaling factors for coastal tiling
 flandg_u(udimsi%i_start:udimsi%i_end,udimsi%j_start:udimsi%j_end),     &
 flandg_v(vdimsi%i_start:vdimsi%i_end,vdimsi%j_start:vdimsi%j_end),     &
                 ! land fraction on u/v grids
 rdz_u(udimsi%i_start:udimsi%i_end,udimsi%j_start:udimsi%j_end,         &
       2:bl_levels),                                                    &
 rdz_v(vdimsi%i_start:vdimsi%i_end,vdimsi%j_start:vdimsi%j_end,         &
       2:bl_levels),                                                    &
                 ! 1 / distance between levels
! Counter gradient stress terms for 1A Bl scheme
   rhogamu_u(udimsi%i_start:udimsi%i_end,udimsi%j_start:udimsi%j_end,   &
                 2:bl_levels),                                          &
   rhogamv_v(vdimsi%i_start:vdimsi%i_end,vdimsi%j_start:vdimsi%j_end,   &
                 2:bl_levels),                                          &
! Counter gradient stress terms for other BL schemes
   f_ngstress_u(udimsi%i_start:udimsi%i_end,udimsi%j_start:udimsi%j_end,&
                 2:bl_levels),                                          &
   f_ngstress_v(vdimsi%i_start:vdimsi%i_end,vdimsi%j_start:vdimsi%j_end,&
                 2:bl_levels),                                          &
   zhnl(pdimsi%i_start:pdimsi%i_end,pdimsi%j_start:pdimsi%j_end)
                   ! non-local boundary layer depth

! Outputs :-
REAL, INTENT(OUT) ::                                                    &
 taux_land(udimso%i_start:udimso%i_end,udimso%j_start:udimso%j_end),    &
                 ! Taux over land part of grid box.
 taux_ssi(udimso%i_start:udimso%i_end,udimso%j_start:udimso%j_end),     &
                 ! Taux over sea part of grid box.
 tauy_land(vdimso%i_start:vdimso%i_end,vdimso%j_start:vdimso%j_end),    &
                 ! Tauy over land part of grid box.
 tauy_ssi(vdimso%i_start:vdimso%i_end,vdimso%j_start:vdimso%j_end),     &
                 ! Tauy over sea part of grid box.
 taux(udimso%i_start:udimso%i_end,udimso%j_start:udimso%j_end,          &
      bl_levels),                                                       &
 tauy(vdimso%i_start:vdimso%i_end,vdimso%j_start:vdimso%j_end,          &
      bl_levels)
                 ! explicit momentum fluxes
!-----------------------------------------------------------------------
! local variables :-

CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'BDY_EXPL3'

REAL ::                                                                 &
 tau_grad_u(udimso%i_start:udimso%i_end,udimso%j_start:udimso%j_end,    &
            bl_levels),                                                 &
 tau_non_grad_u(udimso%i_start:udimso%i_end,udimso%j_start:udimso%j_end,&
            bl_levels),                                                 &
 tau_grad_v(vdimso%i_start:vdimso%i_end,vdimso%j_start:vdimso%j_end,    &
            bl_levels),                                                 &
 tau_non_grad_v(vdimso%i_start:vdimso%i_end,vdimso%j_start:vdimso%j_end,&
            bl_levels)
REAL, ALLOCATABLE :: zhnl_uv (:,:) ! zhnl on u or v points for ex_flux_uv call
INTEGER ::                                                              &
 i, j, offset

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-----------------------------------------------------------------------
IF ( model_type == mt_single_column .OR. l_calc_at_p ) THEN
  offset=0
ELSE
  offset=1
END IF

! surface fluxes
IF (l_ctile .AND. buddy_sea == on) THEN

  DO j = udimso%j_start, udimso%j_end
    DO i = udimso%i_start, udimso%i_end
      taux_land(i,j) = rhokm_u_land(i,j)* u(i,j,1) * flandfac_u(i,j)
      taux_ssi(i,j)  = rhokm_u_ssi(i,j) * ( u(i,j,1) - u_0(i,j) )       &
                                        * fseafac_u(i,j)
      taux(i,j,1) = flandg_u(i,j)*taux_land(i,j)                        &
                    + (1.0-flandg_u(i,j))*taux_ssi(i,j)
    END DO
  END DO

  DO j = vdimso%j_start, vdimso%j_end
    DO i = vdimso%i_start, vdimso%i_end
      tauy_land(i,j) = rhokm_v_land(i,j)* v(i,j,1) * flandfac_v(i,j)
      tauy_ssi(i,j)  = rhokm_v_ssi(i,j) * ( v(i,j,1) - v_0(i,j) )       &
                                        * fseafac_v(i,j)
      tauy(i,j,1) = flandg_v(i,j)*tauy_land(i,j)                        &
                    + (1.0-flandg_v(i,j))*tauy_ssi(i,j)
    END DO
  END DO

ELSE   ! Standard code

  DO j = udimso%j_start, udimso%j_end
    DO i = udimso%i_start, udimso%i_end
      taux_land(i,j) = rhokm_u_land(i,j) * u(i,j,1)
      taux_ssi(i,j) = rhokm_u_ssi(i,j) * ( u(i,j,1) - u_0(i,j) )

      taux(i,j,1) = flandg_u(i,j)*taux_land(i,j)                        &
                    + (1.0-flandg_u(i,j))*taux_ssi(i,j)
    END DO
  END DO

  DO j = vdimso%j_start, vdimso%j_end
    DO i = vdimso%i_start, vdimso%i_end
      tauy_land(i,j) = rhokm_v_land(i,j) * v(i,j,1)
      tauy_ssi(i,j) = rhokm_v_ssi(i,j) * ( v(i,j,1) - v_0(i,j) )

      tauy(i,j,1) = flandg_v(i,j)*tauy_land(i,j)                        &
                    + (1.0-flandg_v(i,j))*tauy_ssi(i,j)
    END DO
  END DO

END IF

! above surface fluxes
IF (i_bl_vn == i_bl_vn_1a ) THEN

  !-----------------------------------------------------------------------
  ! Calculation of explicit fluxes of U and V, 1A BL Scheme method
  !-----------------------------------------------------------------------
  CALL mym_ex_flux_uv(                                                  &
        udimsi, udimsi_s, udimso, bl_levels,                            &
        rdz_u, rhokm_u, rhogamu_u, u, taux_fd_u,                        &
        taux, tau_grad_u, tau_non_grad_u)
  CALL mym_ex_flux_uv(                                                  &
        vdimsi, vdimsi_s, vdimso, bl_levels,                            &
        rdz_v, rhokm_v, rhogamv_v, v, tauy_fd_v,                        &
        tauy, tau_grad_v, tau_non_grad_v)

ELSE
  !-----------------------------------------------------------------------
  ! 5.6 Calculation of explicit fluxes of U and V, other BL schemes method
  !-----------------------------------------------------------------------

  ! Copy ZHNL into an array defined on u-points
  ! In principle, ZHNL should be interpolated, but sensitivity is expected
  ! to be small so the adjacent p-point is used to avoid message passing
  ! This makes endgame formulation consistent with the current formulation
  ALLOCATE(zhnl_uv(udimso%i_start:udimso%i_end,udimso%j_start:udimso%j_end))
  DO j = udimso%j_start, udimso%j_end
    DO i = udimso%i_start, udimso%i_end
      zhnl_uv(i,j) = zhnl(i+offset,j)
    END DO
  END DO

  CALL ex_flux_uv (                                                     &
                         ! For U
    udimsi, udimsi_s, udimso, bl_levels,                                &
    u,zhnl_uv,rdz_u,rhokm_u,f_ngstress_u,taux_fd_u,taux,                &
    tau_grad_u,tau_non_grad_u                                           &
    )
  DEALLOCATE (zhnl_uv)

  ! Copy ZHNL into an array defined on v-points
  ! In principle, ZHNL should be interpolated, but sensitivity is expected
  ! to be small so the adjacent p-point is used to avoid message passing
  ! This makes endgame formulation consistent with the current formulation
  ALLOCATE(zhnl_uv(vdimso%i_start:vdimso%i_end,vdimso%j_start:vdimso%j_end))
  IF (model_type /= mt_single_column) THEN
    DO i = vdimso%i_start, vdimso%i_end
      DO j = vdimso%j_start, vdimso%j_end-1
        zhnl_uv(i,j) = zhnl(i,j+1)
      END DO
      zhnl_uv(i,vdimso%j_end) = zhnl(i,pdimsi%j_end)
    END DO
  ELSE
    DO i = vdimso%i_start, vdimso%i_end
      DO j = vdimso%j_start, vdimso%j_end
        zhnl_uv(i,j) = zhnl(i,j)
      END DO
    END DO
  END IF ! endgame

  CALL ex_flux_uv (                                                     &
                         ! For V
    vdimsi, vdimsi_s, vdimso, bl_levels,                                &
    v,zhnl_uv,rdz_v,rhokm_v,f_ngstress_v,tauy_fd_v,tauy,                &
    tau_grad_v,tau_non_grad_v                                           &
    )
  DEALLOCATE (zhnl_uv)

END IF

! Add additional orographic stress to surface stress over land

IF (formdrag  ==  explicit_stress) THEN

  IF ( l_calc_at_p ) THEN
    ! If this is an extra call to calculate stresses on the p-grid, just
    ! add the explicit orographic surface stresses to the main stress arrays,
    ! the same as is done on other levels in ex_flux_uv

    DO j = udimso%j_start, udimso%j_end
      DO i = udimso%i_start, udimso%i_end
        taux(i,j,1) = taux(i,j,1) + taux_fd_u(i,j,1)
      END DO
    END DO

    DO j = vdimso%j_start, vdimso%j_end
      DO i = vdimso%i_start, vdimso%i_end
        tauy(i,j,1) = tauy(i,j,1) + tauy_fd_v(i,j,1)
      END DO
    END DO

  ELSE
    ! If this is the main call to calculate stresses on their native u,v grids,
    ! add the explicit orographic surface stresses to separate arrays for land
    ! and sea; the main stress arrays get updated using these in ni_imp_ctl.

    DO j = udimso%j_start, udimso%j_end
      DO i = udimso%i_start, udimso%i_end
        IF (flandg_u(i,j) >  0.0) THEN
          taux_land(i,j) = taux_land(i,j) + taux_fd_u(i,j,1)
        END IF
        IF (flandg_u(i,j) <  1.0) THEN
          taux_ssi(i,j) = taux_ssi(i,j) + taux_fd_u(i,j,1)
        END IF
      END DO
    END DO

    DO j = vdimso%j_start, vdimso%j_end
      DO i = vdimso%i_start, vdimso%i_end
        IF (flandg_v(i,j) >  0.0) THEN
          tauy_land(i,j) = tauy_land(i,j) + tauy_fd_v(i,j,1)
        END IF
        IF (flandg_v(i,j) <  1.0) THEN
          tauy_ssi(i,j) = tauy_ssi(i,j) + tauy_fd_v(i,j,1)
        END IF
      END DO
    END DO

  END IF  ! ( l_calc_at_p )

END IF  ! (formdrag  ==  explicit_stress)

IF ( l_scmdiags(scmdiag_bl) .AND.                                       &
     model_type == mt_single_column .AND.                               &
     ( .NOT. l_calc_at_p ) ) THEN
     ! Avoid outputting SCM diags twice when doing
     ! an extra call to compute stresses on p-grid.
  CALL scmoutput(tau_grad_u,'taux_grad',                                &
       'Gradient part of u stress','kg/m/s2',                           &
       t_avg, d_bl, default_streams, '',routinename)

  CALL scmoutput(tau_non_grad_u,'taux_nongrad',                         &
       'Non-gradient part of u stress','kg/m/s2',                       &
       t_avg, d_bl, default_streams, '',routinename)

  CALL scmoutput(tau_grad_v,'tauy_grad',                                &
       'Gradient part of v stress','kg/m/s2',                           &
       t_avg, d_bl, default_streams, '',routinename)

  CALL scmoutput(tau_non_grad_v,'tauy_nongrad',                         &
       'Non-gradient part of v stress','kg/m/s2',                       &
       t_avg, d_bl, default_streams, '',routinename)
END IF ! scmdiag_bl / model_type


!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE bdy_expl3
END MODULE bdy_expl3_mod
