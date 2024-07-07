! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Contains various chunks of code from atm_step - the purpose of each
! section is indicated at the head of the section
! ENDGAME version

! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

SUBROUTINE atm_step_alloc_4A( &
cf_star, cfl_star, cff_star, q_star, qcl_star, qcf_star, &
frac_control, r_u, r_v, r_w, errorstatus, flag)

USE jules_surface_types_mod
USE planet_constants_mod, ONLY: cp

USE Atm_Step_local
USE atm_fields_bounds_mod, ONLY: pdims_s, tdims, tdims_s, tdims_l,     &
                                 udims, udims_s, vdims, vdims_s,       &
                                 wdims, wdims_s

USE jules_snow_mod, ONLY: nsmax
USE timestep_mod
USE level_heights_mod
USE bl_option_mod, ONLY: l_quick_ap2, i_bl_vn, i_bl_vn_1a
USE mym_option_mod, ONLY: bdy_tke, mymodel3
USE water_constants_mod, ONLY: lc
USE cv_run_mod, ONLY: l_conv_hist,                                     &
    l_conv_prog_group_1, l_conv_prog_group_2, l_conv_prog_group_3,     &
    l_conv_prog_precip
USE jules_sea_seaice_mod, ONLY: nice, nice_use, l_ctile
USE jules_surface_mod, ONLY: ISrfExCnvGust, IP_SrfExWithCnv
USE turb_diff_ctl_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE UM_ParVars
USE Control_Max_Sizes
USE Ereport_mod
USE lbc_mod

USE dynamics_input_mod, ONLY:  NumCycles,L_new_tdisc

USE gen_phys_inputs_mod,  ONLY: l_mr_physics

USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain

USE cloud_inputs_mod, ONLY: i_cld_area, i_cld_vn
USE pc2_constants_mod, ONLY: acf_cusack, i_cld_pc2
USE free_tracers_inputs_mod, ONLY: a_tracer_last, a_tracer_first
USE nlsizes_namelist_mod, ONLY:                                          &
    land_field, model_levels, n_cca_lev, n_rows, ntiles,                 &
    river_row_length, river_rows, row_length, rows, sm_levels,           &
    theta_off_size, tr_levels, tr_ukca, tr_vars

USE atm_fields_mod, ONLY: q, qcl, qcf, qcf2, qrain, qgraup, &
    cf_bulk, cf_liquid, cf_frozen, theta, exner_theta_levels, p_theta_levels, &
    frac_con1, frac_con2, frac_con3, frac_con4, frac_con5, frac_con6, &
    frac_con7, frac_con8, frac_con9, u, v, w

USE pc2_rhtl_mod, ONLY: pc2_rhtl
IMPLICIT NONE

! Subroutine arguments
REAL :: q_star  (tdims%i_start:tdims%i_end,                            &
                 tdims%j_start:tdims%j_end,                            &
                 tdims%k_start:tdims%k_end)
REAL :: qcl_star(tdims%i_start:tdims%i_end,                            &
                 tdims%j_start:tdims%j_end,                            &
                 tdims%k_start:tdims%k_end)
REAL :: qcf_star(tdims%i_start:tdims%i_end,                            &
                 tdims%j_start:tdims%j_end,                            &
                 tdims%k_start:tdims%k_end)

REAL :: cf_star (tdims%i_start:tdims%i_end,                            &
                 tdims%j_start:tdims%j_end,                            &
                 tdims%k_start:tdims%k_end)
REAL :: cfl_star(tdims%i_start:tdims%i_end,                            &
                 tdims%j_start:tdims%j_end,                            &
                 tdims%k_start:tdims%k_end)
REAL :: cff_star(tdims%i_start:tdims%i_end,                            &
                 tdims%j_start:tdims%j_end,                            &
                 tdims%k_start:tdims%k_end)

REAL :: frac_control(land_field,ntype)   !Forcing for land surface (3C)

REAL, TARGET :: R_u(udims_s%i_start:udims_s%i_end,                     &
                    udims_s%j_start:udims_s%j_end,                     &
                    udims_s%k_start:udims_s%k_end)
REAL, TARGET :: R_v(vdims_s%i_start:vdims_s%i_end,                     &
                    vdims_s%j_start:vdims_s%j_end,                     &
                    vdims_s%k_start:vdims_s%k_end)
REAL, TARGET :: R_w(wdims%i_start:wdims%i_end,                     &
                    wdims%j_start:wdims%j_end,                     &
                    wdims%k_start:wdims%k_end)

! Local variables

INTEGER :: errorstatus
CHARACTER(LEN=8) :: flag

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ATM_STEP_ALLOC_4A'

INTEGER :: alloc_stat

! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (flag == 'alloc_mr' ) THEN

  ! store physics changes for use in SL_moist_conserve

  ! only the pc2 section here is needed


  IF ( NumCycles > 1 .AND. i_cld_vn == i_cld_pc2 .AND. ErrorStatus == 0 ) THEN

    ! When cycling and PC2 is used cf_star etc need to be reset (at the
    ! beginning of each new cycle) to the value they had when they
    ! exited Physics1(). The following arrays hold these values.

    ALLOCATE ( cf_phys1 (tdims_s%i_start:tdims_s%i_end,              &
                         tdims_s%j_start:tdims_s%j_end,              &
                         tdims_s%k_start:tdims_s%k_end) )
    ALLOCATE ( cfl_phys1(tdims_s%i_start:tdims_s%i_end,              &
                         tdims_s%j_start:tdims_s%j_end,              &
                         tdims_s%k_start:tdims_s%k_end) )
    ALLOCATE ( cff_phys1(tdims_s%i_start:tdims_s%i_end,              &
                         tdims_s%j_start:tdims_s%j_end,              &
                         tdims_s%k_start:tdims_s%k_end) )
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          cf_phys1(i,j,k)  = cf_star(i,j,k)
          cfl_phys1(i,j,k) = cfl_star(i,j,k)
          cff_phys1(i,j,k) = cff_star(i,j,k)
        END DO
      END DO
    END DO
  END IF 

  ! ----------------------------------------

ELSE IF (flag == 'MixRatio') THEN

  ! only the pc2 bit here is needed

  IF (i_cld_vn == i_cld_pc2) THEN
    !  Calculate the intital Relative humidity wrt TL (rhts)
    ALLOCATE(rhts(row_length,rows,model_levels))
    CALL pc2_rhtl(halo_i, halo_j, offx, offy,                            &
                  row_length,rows,                                       &
                  rhts,l_mr_physics)

    !
    IF (i_cld_area == acf_cusack) THEN
      ! Cusack interpolation needs RH on large_levels from start of timestep,
      ! so take TL, qT and p to form RH after interpolation

      ALLOCATE(tlts(tdims%i_start:tdims%i_end,       &
                    tdims%j_start:tdims%j_end,       &
                                1:tdims%k_end))
      ALLOCATE(qtts(tdims%i_start:tdims%i_end,       &
                    tdims%j_start:tdims%j_end,       &
                                1:tdims%k_end))
      ALLOCATE(ptts(tdims%i_start:tdims%i_end,       &
                    tdims%j_start:tdims%j_end,       &
                                1:tdims%k_end))

      DO k = 1, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            tlts(i,j,k) = theta(i,j,k)                              &
                                * exner_theta_levels(i,j,k)         &
                                 - qcl(i,j,k) * lc / cp
            ptts(i,j,k) = p_theta_levels(i,j,k)
            qtts(i,j,k) = q(i,j,k) + qcl(i,j,k)
          END DO
        END DO
      END DO
    ELSE
      ALLOCATE(tlts(1,1,1))
      ALLOCATE(qtts(1,1,1))
      ALLOCATE(ptts(1,1,1))
    END IF !i_cld_area

  END IF ! i_cld_pc2

  ! NB: the star variables and R_u and R_v have not been set in the
  !     halo region yet.

  DO l=1,land_field
    frac_control(l,1)=frac_con1(l)
    frac_control(l,2)=frac_con2(l)
    frac_control(l,3)=frac_con3(l)
    frac_control(l,4)=frac_con4(l)
    frac_control(l,5)=frac_con5(l)
    frac_control(l,6)=frac_con6(l)
    frac_control(l,7)=frac_con7(l)
    frac_control(l,8)=frac_con8(l)
    frac_control(l,9)=frac_con9(l)
  END DO

  ! --------------------------------------------------


ELSE IF (flag == 'newtdisc' ) THEN

  ! Allocations required for new time discretisation

  IF ( l_new_tdisc ) THEN

    ALLOCATE(wetrho_r_sq_np1  (pdims_s%i_start:pdims_s%i_end,       &
                            pdims_s%j_start:pdims_s%j_end,          &
                            pdims_s%k_start:pdims_s%k_end),         &
                            stat=alloc_stat)
    IF (alloc_stat /= 0) CALL Ereport("ATM_STEP_ALLOC_4A",          &
                       ErrorStatus, "Unable to allocate.")

    ALLOCATE(dryrho_np1  (pdims_s%i_start:pdims_s%i_end,            &
                       pdims_s%j_start:pdims_s%j_end,               &
                       pdims_s%k_start:pdims_s%k_end),              &
                       stat=alloc_stat)
    IF (alloc_stat /= 0) CALL Ereport("ATM_STEP_ALLOC_4A",          &
                       ErrorStatus, "Unable to allocate.")

    ALLOCATE(u_np1    (udims_s%i_start:udims_s%i_end,               &
                       udims_s%j_start:udims_s%j_end,               &
                       udims_s%k_start:udims_s%k_end),              &
                       stat=alloc_stat)
    IF (alloc_stat /= 0) CALL Ereport("ATM_STEP_ALLOC_4A",          &
                       ErrorStatus, "Unable to allocate.")

    ALLOCATE(v_np1    (vdims_s%i_start:vdims_s%i_end,               &
                       vdims_s%j_start:vdims_s%j_end,               &
                       vdims_s%k_start:vdims_s%k_end),              &
                       stat=alloc_stat)
    IF (alloc_stat /= 0) CALL Ereport("ATM_STEP_ALLOC_4A",          &
                       ErrorStatus, "Unable to allocate.")

    ALLOCATE(w_np1    (wdims_s%i_start:wdims_s%i_end,               &
                       wdims_s%j_start:wdims_s%j_end,               &
                       wdims_s%k_start:wdims_s%k_end) ,             &
                       stat=alloc_stat)
    IF (alloc_stat /= 0) CALL Ereport("ATM_STEP_ALLOC_4A",          &
                       ErrorStatus, "Unable to allocate.")

    ! ENDGAME-only code
    ALLOCATE(thetav_np1(tdims_s%i_start:tdims_s%i_end,             &
                        tdims_s%j_start:tdims_s%j_end,             &
                        tdims_s%k_start:tdims_s%k_end),            &
                        stat=alloc_stat)
    IF (alloc_stat /= 0) CALL Ereport("ATM_STEP_ALLOC_4A",         &
                      ErrorStatus, "Unable to allocate.")

    ALLOCATE(etadot_np1(wdims_s%i_start:wdims_s%i_end,             &
                        wdims_s%j_start:wdims_s%j_end,             &
                        wdims_s%k_start:wdims_s%k_end),            &
                        stat=alloc_stat)
    IF (alloc_stat /= 0) CALL Ereport("ATM_STEP_ALLOC_4A",         &
                      ErrorStatus, "Unable to allocate.")


    ALLOCATE (m_v_np1(tdims_s%i_start:tdims_s%i_end,               &
                      tdims_s%j_start:tdims_s%j_end,               &
                      tdims_s%k_start:tdims_s%k_end),              &
                      stat=alloc_stat )
    IF (alloc_stat /= 0) CALL Ereport("ATM_STEP_ALLOC_4A",         &
                      ErrorStatus, "Unable to allocate.")


    ALLOCATE (m_cl_np1(tdims_s%i_start:tdims_s%i_end,              &
                       tdims_s%j_start:tdims_s%j_end,              &
                       tdims_s%k_start:tdims_s%k_end),             &
                       stat=alloc_stat )
    IF (alloc_stat /= 0) CALL Ereport("ATM_STEP_ALLOC_4A",         &
                      ErrorStatus, "Unable to allocate.")


    ALLOCATE (m_cf_np1(tdims_s%i_start:tdims_s%i_end,              &
                       tdims_s%j_start:tdims_s%j_end,              &
                       tdims_s%k_start:tdims_s%k_end),             &
                       stat=alloc_stat )
    IF (alloc_stat /= 0) CALL Ereport("ATM_STEP_ALLOC_4A",         &
                      ErrorStatus, "Unable to allocate.")



    ALLOCATE( m_cf2_np1(tdims_s%i_start:tdims_s%i_end,             &
                        tdims_s%j_start:tdims_s%j_end,             &
                        tdims_s%k_start:tdims_s%k_end),            &
                        stat=alloc_stat )
    IF (alloc_stat /= 0) CALL Ereport("ATM_STEP_ALLOC_4A",         &
                      ErrorStatus, "Unable to allocate.")


    ALLOCATE( m_r_np1(tdims_s%i_start:tdims_s%i_end,               &
                      tdims_s%j_start:tdims_s%j_end,               &
                      tdims_s%k_start:tdims_s%k_end) ,             &
                      stat=alloc_stat)
    IF (alloc_stat /= 0) CALL Ereport("ATM_STEP_ALLOC_4A",         &
                      ErrorStatus, "Unable to allocate.")



    ALLOCATE( m_gr_np1(tdims_s%i_start:tdims_s%i_end,              &
                       tdims_s%j_start:tdims_s%j_end,              &
                       tdims_s%k_start:tdims_s%k_end),             &
                       stat=alloc_stat )
    IF (alloc_stat /= 0) CALL Ereport("ATM_STEP_ALLOC_4A",         &
                      ErrorStatus, "Unable to allocate.")



    ALLOCATE( exner_surf_np1(pdims_s%i_start:pdims_s%i_end,        &
                             pdims_s%j_start:pdims_s%j_end),       &
                             stat=alloc_stat)
    IF (alloc_stat /= 0) CALL Ereport("ATM_STEP_ALLOC_4A",         &
                      ErrorStatus, "Unable to allocate.")


    ALLOCATE( exner_np1(pdims_s%i_start:pdims_s%i_end,             &
                        pdims_s%j_start:pdims_s%j_end,             &
                        pdims_s%k_start:pdims_s%k_end+1),          &
                        stat=alloc_stat)
    IF (alloc_stat /= 0) CALL Ereport("ATM_STEP_ALLOC_4A",         &
                         ErrorStatus, "Unable to allocate.")

    ! End of ENDGAME-only code


  ELSE

    ALLOCATE( u_np1(1,1,1) )
    ALLOCATE( v_np1(1,1,1) )
    ALLOCATE( w_np1(1,1,1) )
    ALLOCATE( dryrho_np1 (1,1,1) )

  END IF  ! L_new_tdisc

  ! ---------------------------------------------------

ELSE IF (flag == 'lbc_updt' ) THEN

  ! Code from atm_step required after updating primary fields with lam lbc data


  IF ( NumCycles > 1 ) THEN

    DEALLOCATE( bulk_cld_frac_phys1 )
    DEALLOCATE( bulk_cld_liq_phys1 )
    DEALLOCATE( bulk_cld_fr_phys1 )
    DEALLOCATE( area_cld_frac_phys1 )
    IF (i_bl_vn == i_bl_vn_1a) THEN
      IF (bdy_tke == mymodel3) THEN
        DEALLOCATE (cov_trb_phys1 )
        DEALLOCATE (qsq_trb_phys1 )
        DEALLOCATE (tsq_trb_phys1 )
      END IF
      DEALLOCATE( e_trb_phys1 )
    END IF
    DEALLOCATE( ti_phys1 )
    IF (.NOT. l_quick_ap2) DEALLOCATE( ti_gb_phys1 )
    DEALLOCATE( zh_phys1 )
    DEALLOCATE( z0msea_phys1 )
    DEALLOCATE( cca_phys1 )
    DEALLOCATE( ccb_phys1 )
    DEALLOCATE( cct_phys1 )

    IF (ISrfExCnvGust == IP_SrfExWithCnv) THEN
      DEALLOCATE( ddmfx_phys1 )
    END IF

    IF (l_conv_hist) THEN
      DEALLOCATE( past_conv_ht_phys1 )
      DEALLOCATE( past_precip_phys1 )
      DEALLOCATE( deep_flag_phys1 )
    END IF

    IF ( l_conv_prog_group_1 ) DEALLOCATE( conv_prog_1_phys1 )
    IF ( l_conv_prog_group_2 ) DEALLOCATE( conv_prog_2_phys1 )
    IF ( l_conv_prog_group_3 ) DEALLOCATE( conv_prog_3_phys1 )
    IF ( l_conv_prog_precip )  DEALLOCATE( conv_prog_precip_phys1 )

    IF ( L_ctile ) THEN
      DEALLOCATE( t_land_ctile_phys1 )
      DEALLOCATE( t_sea_ctile_phys1 )
      DEALLOCATE( t_sice_ctile_phys1 )
    END IF

    DEALLOCATE( t_surf_phys1 )
    DEALLOCATE( t_sf_tile_phys1 )
    DEALLOCATE( snow_tile_phys1 )
    DEALLOCATE( dolr_phys1 )

  END IF

  !deallocate extra microphysics variables
  DEALLOCATE (qcf2_star)
  DEALLOCATE (qrain_star)
  DEALLOCATE (qgraup_star)

  IF ( NumCycles >1 .AND. i_cld_vn == i_cld_pc2 ) THEN

    ! If cycling and i_cld_pc2
    ! is used deallocate cf_phys1, cfl_phys1, cff_phys1
    DEALLOCATE (cf_phys1)
    DEALLOCATE (cfl_phys1)
    DEALLOCATE (cff_phys1)

  END IF ! 


  ! ---------------------------------------------------


ELSE IF (flag == 'allocPC2' ) THEN

  ! Allocations for PC2

  ALLOCATE(t_inc_pres  (tdims%i_start:tdims%i_end,               &
                        tdims%j_start:tdims%j_end,               &
                                    1:tdims%k_end))
  ALLOCATE(q_inc_pres  (tdims%i_start:tdims%i_end,               &
                        tdims%j_start:tdims%j_end,               &
                                    1:tdims%k_end))
  ALLOCATE(qcl_inc_pres(tdims%i_start:tdims%i_end,               &
                        tdims%j_start:tdims%j_end,               &
                                    1:tdims%k_end))
  ALLOCATE(qcf_inc_pres(tdims%i_start:tdims%i_end,               &
                        tdims%j_start:tdims%j_end,               &
                                    1:tdims%k_end))
  ALLOCATE(cf_inc_pres (tdims%i_start:tdims%i_end,               &
                        tdims%j_start:tdims%j_end,               &
                                    1:tdims%k_end))
  ALLOCATE(cfl_inc_pres(tdims%i_start:tdims%i_end,               &
                        tdims%j_start:tdims%j_end,               &
                                    1:tdims%k_end))
  ALLOCATE(cff_inc_pres(tdims%i_start:tdims%i_end,               &
                        tdims%j_start:tdims%j_end,               &
                                    1:tdims%k_end))
  ALLOCATE(t_dini      (tdims%i_start:tdims%i_end,               &
                        tdims%j_start:tdims%j_end,               &
                                    1:tdims%k_end))
  ALLOCATE(q_dini      (tdims%i_start:tdims%i_end,               &
                        tdims%j_start:tdims%j_end,               &
                                    1:tdims%k_end ))
  ALLOCATE(qcl_dini    (tdims%i_start:tdims%i_end,               &
                        tdims%j_start:tdims%j_end,               &
                                    1:tdims%k_end ))
  ALLOCATE(qcf_dini    (tdims%i_start:tdims%i_end,               &
                        tdims%j_start:tdims%j_end,               &
                                    1:tdims%k_end ))
  ALLOCATE(cf_dini     (tdims%i_start:tdims%i_end,               &
                        tdims%j_start:tdims%j_end,               &
                                    1:tdims%k_end))
  ALLOCATE(cfl_dini    (tdims%i_start:tdims%i_end,               &
                        tdims%j_start:tdims%j_end,               &
                                    1:tdims%k_end))
  ALLOCATE(cff_dini    (tdims%i_start:tdims%i_end,               &
                        tdims%j_start:tdims%j_end,               &
                                    1:tdims%k_end))

  ! ---------------------------------------------------

ELSE IF ( flag == 'dealsmag' ) THEN

  ! Deallocations for Smagorinsky

  DEALLOCATE (delta_smag)
  DEALLOCATE (max_diff)
  DEALLOCATE (shear)
  DEALLOCATE (rneutml_sq)
  DEALLOCATE (visc_m)
  DEALLOCATE (visc_h)

  ! ---------------------------------------------------

ELSE IF (flag == 'dllocPC2' ) THEN

  ! Deallocations for PC2

  DEALLOCATE(cff_dini)
  DEALLOCATE(cfl_dini)
  DEALLOCATE(cf_dini)
  DEALLOCATE(qcf_dini)
  DEALLOCATE(qcl_dini)
  DEALLOCATE(q_dini)
  DEALLOCATE(t_dini)
  DEALLOCATE(cff_inc_pres)
  DEALLOCATE(cfl_inc_pres)
  DEALLOCATE(cf_inc_pres)
  DEALLOCATE(qcf_inc_pres)
  DEALLOCATE(qcl_inc_pres)
  DEALLOCATE(q_inc_pres)
  DEALLOCATE(t_inc_pres)

END IF ! flag
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE Atm_Step_alloc_4A
