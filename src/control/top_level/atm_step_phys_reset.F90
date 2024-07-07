! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Contains various chunks of code from atm_step - the purpose of each
! section is indicated at the head of the section
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Subroutine Interface:
SUBROUTINE atm_step_phys_reset( &
q_star, qcl_star, qcf_star, r_v, r_u, dolr, theta_star, gs1, &
cf_star, cfl_star, cff_star, tstar_tile, snodep_tile, flag)

USE atm_step_local
USE atm_fields_bounds_mod, ONLY : tdims, tdims_s, udims, udims_s,      &
                                  vdims, vdims_s

USE jules_sea_seaice_mod, ONLY: nice, nice_use, l_ctile
USE jules_surface_mod, ONLY: ISrfExCnvGust, IP_SrfExWithCnv
USE jules_surface_types_mod, ONLY: ntype, npft
USE cv_run_mod, ONLY: l_conv_hist,                                     &
    l_conv_prog_group_1, l_conv_prog_group_2, l_conv_prog_group_3,     &
    l_conv_prog_precip
USE jules_snow_mod, ONLY: nsmax
USE bl_option_mod, ONLY: l_quick_ap2, i_bl_vn, i_bl_vn_1a
USE mym_option_mod, ONLY: bdy_tke, mymodel3
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE UM_ParVars
USE Control_Max_Sizes
USE lbc_mod

USE dynamics_input_mod, ONLY:  NumCycles,                  &
    alpha_1,alpha_2,alpha_3,alpha_4,                       &
    alpha_1_2, alpha_2_2, alpha_3_2, alpha_4_2

USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain
USE cloud_inputs_mod, ONLY: i_cld_vn
USE pc2_constants_mod, ONLY: i_cld_pc2
USE free_tracers_inputs_mod, ONLY: a_tracer_last, a_tracer_first

USE nlsizes_namelist_mod, ONLY:                                          &
    land_field, model_levels, n_cca_lev, n_rows, ntiles,                 &
    river_row_length, river_rows, row_length, rows, sm_levels,           &
    theta_off_size, tr_levels, tr_ukca, tr_vars

USE atm_fields_mod, ONLY: cf_bulk, cf_frozen, cf_liquid, cf_area, cca, ccb,  &
    cct, ti, z0, zh, gs, ddmfx, deep_flag, past_precip, past_conv_ht, tstar, &
    tstar_land, tstar_sea, ti_cat, tstar_sice, tstar_sice_cat, conv_prog_1,  &
    conv_prog_2, conv_prog_3, conv_prog_precip, e_trb, tsq_trb, qsq_trb,     &
    cov_trb
! NOTE: tstar_tile, snodep_tile passed as args for now, 
!       rather than USEd from atm_fields_mod

IMPLICIT NONE


! Subroutine arguments

REAL, TARGET :: r_u(udims_s%i_start:udims_s%i_end,                     &
                    udims_s%j_start:udims_s%j_end,                     &
                    udims_s%k_start:udims_s%k_end)
REAL, TARGET :: r_v(vdims_s%i_start:vdims_s%i_end,                     &
                    vdims_s%j_start:vdims_s%j_end,                     &
                    vdims_s%k_start:vdims_s%k_end)

REAL :: theta_star(tdims%i_start:tdims%i_end,                          &
                   tdims%j_start:tdims%j_end,                          &
                   tdims%k_start:tdims%k_end)
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

REAL :: tstar_tile(land_field,ntiles)
REAL :: snodep_tile(land_field,ntiles)

REAL :: dOLR(row_length,rows)
REAL :: gs1(land_field)

! Local variables

CHARACTER(LEN=8) :: flag

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ATM_STEP_PHYS_RESET'

! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (flag == 'phy1rest') THEN

  ! Restore phys1 variables to be used as predictors.
  IF ( cycleno > 1 ) THEN

    ! reset weights after the first cycle
    alpha1 = alpha_1_2
    alpha2 = alpha_2_2
    alpha3 = alpha_3_2
    alpha4 = alpha_4_2

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)          &
!$OMP  SHARED(tdims,cf_bulk,cf_liquid,cf_frozen,cf_area,bulk_cld_frac_phys1, &
!$OMP  bulk_cld_liq_phys1,bulk_cld_fr_phys1,area_cld_frac_phys1, &
!$OMP  cf_star,cfl_star,cff_star,cf_phys1,cfl_phys1,cff_phys1, &
!$OMP  n_cca_lev,cca,cca_phys1,nice,p_ti,ti_phys1,i_cld_vn, &
!$OMP  e_trb, e_trb_phys1, tsq_trb, tsq_trb_phys1, qsq_trb, qsq_trb_phys1, &
!$OMP  cov_trb, cov_trb_phys1, i_bl_vn, bdy_tke)

!$OMP DO SCHEDULE(STATIC)
    DO k= 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        !CDIR NODEP
        DO i = tdims%i_start, tdims%i_end
          cf_bulk(i,j,k)   = bulk_cld_frac_phys1(i,j,k)
          cf_liquid(i,j,k) = bulk_cld_liq_phys1(i,j,k)
          cf_frozen(i,j,k) = bulk_cld_fr_phys1(i,j,k)
          cf_area(i,j,k)   = area_cld_frac_phys1(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

    IF (i_bl_vn == i_bl_vn_1a) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k= tdims%k_start, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          !CDIR NODEP
          DO i = tdims%i_start, tdims%i_end
            e_trb(i,j,k)   = e_trb_phys1(i,j,k)
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
      IF (bdy_tke == mymodel3) THEN
!$OMP DO SCHEDULE(STATIC)
        DO k= tdims%k_start, tdims%k_end
          DO j = tdims%j_start, tdims%j_end
            !CDIR NODEP
            DO i = tdims%i_start, tdims%i_end
              tsq_trb(i,j,k) = tsq_trb_phys1(i,j,k)
              qsq_trb(i,j,k) = qsq_trb_phys1(i,j,k)
              cov_trb(i,j,k) = cov_trb_phys1(i,j,k)
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT
      END IF
    END IF

    IF ( i_cld_vn == i_cld_pc2 ) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k=tdims%k_start, tdims%k_end
        DO j=tdims%j_start, tdims%j_end
          DO i=tdims%i_start, tdims%i_end
            cf_star(i,j,k)  = cf_phys1(i,j,k)
            cfl_star(i,j,k) = cfl_phys1(i,j,k)
            cff_star(i,j,k) = cff_phys1(i,j,k)
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF


!$OMP DO SCHEDULE(STATIC)
    DO k=1, n_cca_lev
      DO j=tdims%j_start, tdims%j_end
        !CDIR NODEP
        DO i=tdims%i_start, tdims%i_end
          cca(i,j,k) = cca_phys1(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
    DO k=1, nice
      DO j=tdims%j_start, tdims%j_end
        !CDIR NODEP
        DO i=tdims%i_start, tdims%i_end
          p_ti(i,j,k) = ti_phys1(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    DO j=tdims%j_start, tdims%j_end
      !CDIR NODEP
      DO i=tdims%i_start, tdims%i_end
        ccb(i,j) = ccb_phys1(i,j)
        cct(i,j) = cct_phys1(i,j)
      END DO
    END DO

    IF (.NOT. l_quick_ap2) THEN
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          ti(i,j,1) = ti_gb_phys1(i,j)
        END DO
      END DO
    END IF

    DO j=tdims%j_start, tdims%j_end
      !CDIR NODEP
      DO i=tdims%i_start, tdims%i_end
        z0(i,j)  = z0msea_phys1(i,j)
        zh(i,j)  = zh_phys1(i,j)
      END DO
    END DO

    DO i=1, land_field
      gs(i) = gs1(i)
    END DO

    IF (ISrfExCnvGust == IP_SrfExWithCnv) THEN
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          ddmfx(i,j) = ddmfx_phys1(i,j)
        END DO
      END DO
    END IF

    IF (l_conv_hist) THEN
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          deep_flag(i,j) = deep_flag_phys1(i,j)
          past_precip(i,j) = past_precip_phys1(i,j)
          past_conv_ht(i,j) = past_conv_ht_phys1(i,j)
        END DO
      END DO
    END IF

    IF ( l_conv_prog_group_1 ) THEN
      DO k=tdims%k_start, tdims%k_end
        DO j=tdims%j_start, tdims%j_end
          DO i=tdims%i_start, tdims%i_end
            conv_prog_1(i,j,k) = conv_prog_1_phys1(i,j,k)
          END DO
        END DO
      END DO
    END IF
    IF ( l_conv_prog_group_2 ) THEN
      DO k=tdims%k_start, tdims%k_end
        DO j=tdims%j_start, tdims%j_end
          DO i=tdims%i_start, tdims%i_end
            conv_prog_2(i,j,k) = conv_prog_2_phys1(i,j,k)
          END DO
        END DO
      END DO
    END IF
    IF ( l_conv_prog_group_3 ) THEN
      DO k=tdims%k_start, tdims%k_end
        DO j=tdims%j_start, tdims%j_end
          DO i=tdims%i_start, tdims%i_end
            conv_prog_3(i,j,k) = conv_prog_3_phys1(i,j,k)
          END DO
        END DO
      END DO
    END IF
    IF ( l_conv_prog_precip ) THEN
      DO k=tdims%k_start, tdims%k_end
        DO j=tdims%j_start, tdims%j_end
          DO i=tdims%i_start, tdims%i_end
            conv_prog_precip(i,j,k) = conv_prog_precip_phys1(i,j,k)
          END DO
        END DO
      END DO
    END IF


    DO j=tdims%j_start, tdims%j_end
      !CDIR NODEP
      DO i=tdims%i_start, tdims%i_end
        tstar(i,j) = t_surf_phys1(i,j)
        dOLR(i,j) = dolr_phys1(i,j)
      END DO
    END DO

    IF ( l_ctile ) THEN
      DO j=tdims%j_start, tdims%j_end
        !CDIR NODEP
        DO i=tdims%i_start, tdims%i_end
          tstar_land(i,j) = t_land_ctile_phys1(i,j)
          tstar_sea(i,j) = t_sea_ctile_phys1(i,j)
        END DO
      END DO

      DO k=1,nice_use
        DO j=tdims%j_start, tdims%j_end
          !CDIR NODEP
          DO i=tdims%i_start, tdims%i_end
            p_tstar_sice(i,j,k) = t_sice_ctile_phys1(i,j,k)
          END DO
        END DO
      END DO
    END IF

    DO j = 1, ntiles
      !CDIR NODEP
      DO i = 1, land_field
        tstar_tile (i,j) = t_sf_tile_phys1(i,j)
        snodep_tile(i,j) = snow_tile_phys1(i,j)
      END DO
    END DO

  ELSE

    ! namelist user defined alphas at first cycle
    alpha1 = alpha_1
    alpha2 = alpha_2
    alpha3 = alpha_3
    alpha4 = alpha_4

  END IF ! Cycle No

  ! ------------------------------------------------------------------

ELSE IF (flag == 'cyclrset') THEN

  IF ( numcycles > 1 ) THEN

    ! When cycling the following variables need to be reset (at the
    ! beginning of each new cycle) to the value they had when they
    ! exited Physics1().

    ALLOCATE( area_cld_frac_phys1(tdims%i_start:tdims%i_end,         &
                                  tdims%j_start:tdims%j_end,         &
                                  tdims%k_start:tdims%k_end) )
    ALLOCATE( bulk_cld_frac_phys1(tdims%i_start:tdims%i_end,         &
                                  tdims%j_start:tdims%j_end,         &
                                  tdims%k_start:tdims%k_end) )
    ALLOCATE( bulk_cld_liq_phys1 (tdims%i_start:tdims%i_end,         &
                                  tdims%j_start:tdims%j_end,         &
                                  tdims%k_start:tdims%k_end) )
    ALLOCATE( bulk_cld_fr_phys1  (tdims%i_start:tdims%i_end,         &
                                  tdims%j_start:tdims%j_end,         &
                                  tdims%k_start:tdims%k_end) )

    IF (i_bl_vn == i_bl_vn_1a) THEN
      ALLOCATE( e_trb_phys1(tdims%i_start:tdims%i_end,                  &
                            tdims%j_start:tdims%j_end,                  &
                            tdims%k_start:tdims%k_end) )
      IF (bdy_tke == mymodel3) THEN
        ALLOCATE( tsq_trb_phys1(tdims%i_start:tdims%i_end,              &
                                tdims%j_start:tdims%j_end,              &
                                tdims%k_start:tdims%k_end) )
        ALLOCATE( qsq_trb_phys1(tdims%i_start:tdims%i_end,              &
                                tdims%j_start:tdims%j_end,              &
                                tdims%k_start:tdims%k_end) )
        ALLOCATE( cov_trb_phys1(tdims%i_start:tdims%i_end,              &
                                tdims%j_start:tdims%j_end,              &
                                tdims%k_start:tdims%k_end) )
      END IF
    END IF

    ALLOCATE( ti_phys1(row_length, rows, nice) )
    IF (.NOT. l_quick_ap2) ALLOCATE( ti_gb_phys1(row_length, rows) )
    ALLOCATE( zh_phys1(row_length, rows) )
    ALLOCATE( z0msea_phys1(row_length, rows) )
    ALLOCATE( cca_phys1 (row_length, rows, n_cca_lev) )
    ALLOCATE( ccb_phys1 (row_length, rows ) )
    ALLOCATE( cct_phys1 (row_length, rows ) )

    IF (ISrfExCnvGust == IP_SrfExWithCnv) THEN
      ! downdraft mass flux at cloud base
      ! time level n-1 output from convection scheme required by time level n
      ! boundary layer scheme, hence needs resetting each cycle
      ALLOCATE( ddmfx_phys1(row_length, rows) )
    END IF

    IF (l_conv_hist) THEN
      ! convective history diagnostics
      ! time level n-1 output from convection scheme required by time level n
      ! convection scheme, hence need resetting each cycle
      ALLOCATE( deep_flag_phys1(row_length,rows) )
      ALLOCATE( past_precip_phys1(row_length,rows) )
      ALLOCATE( past_conv_ht_phys1(row_length,rows) )
    END IF

    IF ( l_conv_prog_group_1 ) THEN
      ALLOCATE( conv_prog_1_phys1(tdims%i_start:tdims%i_end,        &
                                  tdims%j_start:tdims%j_end,        &
                                  tdims%k_start:tdims%k_end) )
    END IF
    IF ( l_conv_prog_group_2 ) THEN
      ALLOCATE( conv_prog_2_phys1(tdims%i_start:tdims%i_end,        &
                                  tdims%j_start:tdims%j_end,        &
                                  tdims%k_start:tdims%k_end) )
    END IF
    IF ( l_conv_prog_group_3 ) THEN
      ALLOCATE( conv_prog_3_phys1(tdims%i_start:tdims%i_end,        &
                                  tdims%j_start:tdims%j_end,        &
                                  tdims%k_start:tdims%k_end) )
    END IF
    IF ( l_conv_prog_precip ) THEN
      ALLOCATE( conv_prog_precip_phys1(tdims%i_start:tdims%i_end,   &
                                       tdims%j_start:tdims%j_end,   &
                                       tdims%k_start:tdims%k_end) )
    END IF

    IF ( l_ctile ) THEN
      ALLOCATE( t_land_ctile_phys1(row_length,rows) )
      ALLOCATE( t_sea_ctile_phys1(row_length,rows) )
      ALLOCATE( t_sice_ctile_phys1(row_length,rows,nice_use) )
    END IF

    ALLOCATE( t_surf_phys1(row_length, rows) )
    ALLOCATE( t_sf_tile_phys1(land_field,ntiles) )
    ALLOCATE( snow_tile_phys1(land_field,ntiles) )
    ALLOCATE( dolr_phys1(row_length,rows) )

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)          &
!$OMP  SHARED(tdims,cf_bulk,cf_liquid,cf_frozen,cf_area,bulk_cld_frac_phys1, &
!$OMP  bulk_cld_liq_phys1,bulk_cld_fr_phys1,area_cld_frac_phys1, &
!$OMP  cf_star,cfl_star,cff_star,cf_phys1,cfl_phys1,cff_phys1, &
!$OMP  n_cca_lev,cca,cca_phys1,nice,p_ti,ti_phys1,ti,ti_cat, &
!$OMP  e_trb, e_trb_phys1, tsq_trb, tsq_trb_phys1, qsq_trb, qsq_trb_phys1, &
!$OMP  cov_trb, cov_trb_phys1, i_bl_vn, bdy_tke)

!$OMP DO SCHEDULE(STATIC)
    DO k= 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        !CDIR NODEP
        DO i = tdims%i_start, tdims%i_end
          bulk_cld_frac_phys1(i,j,k) = cf_bulk(i,j,k)
          bulk_cld_liq_phys1(i,j,k)  = cf_liquid(i,j,k)
          bulk_cld_fr_phys1(i,j,k)   = cf_frozen(i,j,k)
          area_cld_frac_phys1(i,j,k) = cf_area(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

    IF (i_bl_vn == i_bl_vn_1a) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k= tdims%k_start, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          !CDIR NODEP
          DO i = tdims%i_start, tdims%i_end
            e_trb_phys1(i,j,k) = e_trb(i,j,k)
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
      IF (bdy_tke == mymodel3) THEN
!$OMP DO SCHEDULE(STATIC)
        DO k= tdims%k_start, tdims%k_end
          DO j = tdims%j_start, tdims%j_end
            !CDIR NODEP
            DO i = tdims%i_start, tdims%i_end
              tsq_trb_phys1(i,j,k) = tsq_trb(i,j,k)
              qsq_trb_phys1(i,j,k) = qsq_trb(i,j,k)
              cov_trb_phys1(i,j,k) = cov_trb(i,j,k)
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT
      END IF
    END IF

!$OMP DO SCHEDULE(STATIC)
    DO k=1, n_cca_lev
      DO j=tdims%j_start, tdims%j_end
        !CDIR NODEP
        DO i=tdims%i_start, tdims%i_end
          cca_phys1(i,j,k) = cca(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP SINGLE
    ! The following code is intended to keep surface
    ! vars at the same timelevel at each cycle.
    IF (nice == 1) THEN
      ! Set sea ice category scheme D1 pointers as catogories are not present
      ! if nice = 1.
      p_ti=> ti
    ELSE
      p_ti=> ti_cat
    END IF
!$OMP END SINGLE
! Implicit barrier as we use p_ti below

!$OMP DO SCHEDULE(STATIC)
    DO k=1, nice
      DO j=tdims%j_start, tdims%j_end
        !CDIR NODEP
        DO i=tdims%i_start, tdims%i_end
          ti_phys1(i,j,k) = p_ti(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    DO j=tdims%j_start, tdims%j_end
      !CDIR NODEP
      DO i=tdims%i_start, tdims%i_end
        ccb_phys1(i,j)=ccb(i,j)
        cct_phys1(i,j)=cct(i,j)
      END DO
    END DO

    IF (.NOT. l_quick_ap2) THEN
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          ti_gb_phys1(i,j) = ti(i,j,1)
        END DO
      END DO
    END IF

    DO j=tdims%j_start, tdims%j_end
      !CDIR NODEP
      DO i=tdims%i_start, tdims%i_end
        z0msea_phys1(i,j) = z0(i,j)
        zh_phys1(i,j) = zh(i,j)
      END DO
    END DO

    IF (ISrfExCnvGust == IP_SrfExWithCnv) THEN
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          ddmfx_phys1(i,j) = ddmfx(i,j)
        END DO
      END DO
    END IF

    IF (l_conv_hist) THEN
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          deep_flag_phys1(i,j) = deep_flag(i,j)
          past_precip_phys1(i,j) = past_precip(i,j)
          past_conv_ht_phys1(i,j) = past_conv_ht(i,j)
        END DO
      END DO
    END IF

    IF (l_conv_prog_group_1) THEN
      DO k = tdims%k_start, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
!CDIR NODEP
          DO i = tdims%i_start, tdims%i_end
            conv_prog_1_phys1(i,j,k) = conv_prog_1(i,j,k)
          END DO
        END DO
      END DO
    END IF
    IF (l_conv_prog_group_2) THEN
      DO k = tdims%k_start, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
!CDIR NODEP
          DO i = tdims%i_start, tdims%i_end
            conv_prog_2_phys1(i,j,k) = conv_prog_2(i,j,k)
          END DO
        END DO
      END DO
    END IF
    IF (l_conv_prog_group_3) THEN
      DO k = tdims%k_start, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
!CDIR NODEP
          DO i = tdims%i_start, tdims%i_end
            conv_prog_3_phys1(i,j,k) = conv_prog_3(i,j,k)
          END DO
        END DO
      END DO
    END IF
    IF (l_conv_prog_precip) THEN
      DO k = tdims%k_start, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
!CDIR NODEP
          DO i = tdims%i_start, tdims%i_end
            conv_prog_precip_phys1(i,j,k) = conv_prog_precip(i,j,k)
          END DO
        END DO
      END DO
    END IF

    DO j=tdims%j_start, tdims%j_end
      !CDIR NODEP
      DO i=tdims%i_start, tdims%i_end
        t_surf_phys1(i,j) = tstar(i,j)
        dolr_phys1(i,j)   = dOLR(i,j)
      END DO
    END DO

    DO i=1, land_field
      gs1(i) = gs(i)
    END DO

    IF (nice_use == 1) THEN
      ! Set sea ice category scheme D1 pointers as categories are not present
      ! if nice_use = 1.
      p_tstar_sice => tstar_sice
    ELSE
      p_tstar_sice => tstar_sice_cat
    END IF

    IF ( l_ctile ) THEN
      DO j=tdims%j_start, tdims%j_end
        !CDIR NODEP
        DO i=tdims%i_start, tdims%i_end
          t_land_ctile_phys1(i,j) = tstar_land(i,j)
          t_sea_ctile_phys1(i,j) = tstar_sea(i,j)
        END DO
      END DO
      DO k=1,nice_use
        DO j=tdims%j_start, tdims%j_end
          !CDIR NODEP
          DO i=tdims%i_start, tdims%i_end
            t_sice_ctile_phys1(i,j,k) = p_tstar_sice(i,j,k)
          END DO
        END DO
      END DO

    END IF

    DO j = 1, ntiles
      !CDIR NODEP
      DO i = 1, land_field
        t_sf_tile_phys1(i,j) = tstar_tile(i,j)
        snow_tile_phys1(i,j) = snodep_tile(i,j)
      END DO
    END DO

  END IF  ! Num_Cycles > 1

END IF  ! Flag
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE atm_step_phys_reset
