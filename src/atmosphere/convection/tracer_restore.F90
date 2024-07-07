! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Restore tracers to their original arrays from tot_tracer

MODULE tracer_restore_mod

IMPLICIT NONE
!---------------------------------------------------------------------------
! Description: Copies all required tracers from tot_tracer back to separate
!              arrays

! method:

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection

! code description:
!   language: fortran 90
!   This code is written to umdp3 programming standards version 10.4
!---------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TRACER_RESTORE_MOD'

CONTAINS

SUBROUTINE  tracer_restore(row_length, rows, ntra_fld, ntra_lev   &
, tr_vars, tr_ukca, tot_tracer                                    &
, aerosol                                                         &
, dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6     &
, so2, so4_aitken, so4_accu, so4_diss                             &
, dms, nh3, soot_new, soot_agd, soot_cld, bmass_new, bmass_agd    &
, bmass_cld, ocff_new, ocff_agd, ocff_cld, nitr_acc, nitr_diss    &
, co2,  free_tracers, tracer_ukca                                 &
, ozone_tracer )


!$ USE omp_lib

! Classic aerosol flags
USE run_aerosol_mod, ONLY: l_sulpc_so2, l_soot, l_ocff, l_biomass,  &
                           l_sulpc_nh3, l_nitrate, l_sulpc_dms

USE dust_parameters_mod, ONLY:  l_dust, l_twobin_dust

USE ukca_option_mod,  ONLY: l_ukca_mode
USE ukca_d1_defs,     ONLY: l_ukca_plume_diags
USE ukca_scavenging_diags_mod, ONLY: ukca_plume_scav_final

USE carbon_options_mod, ONLY: l_co2_interactive     

USE rad_input_mod, ONLY: l_use_cariolle

USE cv_run_mod,  ONLY:  l_murk_conv,                                &
                        i_convection_vn,                            &
                        i_convection_vn_6a

USE atm_fields_bounds_mod, ONLY: tdims_s
USE nlsizes_namelist_mod,   ONLY: model_levels

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

!----------------------------------------------------------------------
! Subroutine arguments
INTEGER, INTENT(IN) ::       &
  row_length                 & ! Row length
 ,rows                       & ! Number of rows
 ,ntra_fld                   & ! Total number of tracers
 ,ntra_lev                   & ! Number of tracer levels either 1 or model
                               ! levels
 ,tr_vars                    & ! number of free tracers
 ,tr_ukca                      ! number of ukca tracers

REAL,INTENT(IN) ::           &
  tot_tracer(row_length, rows, ntra_lev, ntra_fld)

! Tracers
REAL,INTENT(INOUT) ::                         &
  aerosol(tdims_s%i_start:tdims_s%i_end,      & ! aerosol tracer
          tdims_s%j_start:tdims_s%j_end,      &
          tdims_s%k_start:tdims_s%k_end)      &
, dust_div1(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div1
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, dust_div2(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div2
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, dust_div3(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div3
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, dust_div4(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div4
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, dust_div5(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div5
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, dust_div6(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div6
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)
                                                !DIV1,DIV2,...,DIV6 DUST
REAL,INTENT(INOUT) ::                         &
  so2(tdims_s%i_start:tdims_s%i_end,          & ! SO2
      tdims_s%j_start:tdims_s%j_end,          &
      tdims_s%k_start:tdims_s%k_end)          &
, so4_aitken(tdims_s%i_start:tdims_s%i_end,   & ! SO4 aitken mode
             tdims_s%j_start:tdims_s%j_end,   &
             tdims_s%k_start:tdims_s%k_end)   &
, so4_accu(  tdims_s%i_start:tdims_s%i_end,   & ! SO4 accumulation mode
             tdims_s%j_start:tdims_s%j_end,   &
             tdims_s%k_start:tdims_s%k_end)   &
, so4_diss(tdims_s%i_start:tdims_s%i_end,     & ! SO4 dissipation mode
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, dms(tdims_s%i_start:tdims_s%i_end,          & ! DMS
      tdims_s%j_start:tdims_s%j_end,          &
      tdims_s%k_start:tdims_s%k_end)          &
, nh3(tdims_s%i_start:tdims_s%i_end,          & ! NH3
      tdims_s%j_start:tdims_s%j_end,          &
      tdims_s%k_start:tdims_s%k_end)          &
, soot_new(tdims_s%i_start:tdims_s%i_end,     & ! New soot
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, soot_agd(tdims_s%i_start:tdims_s%i_end,     & ! Aged soot
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, soot_cld(tdims_s%i_start:tdims_s%i_end,     & ! soot in cloud
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, bmass_new(tdims_s%i_start:tdims_s%i_end,    & ! New biomass
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, bmass_agd(tdims_s%i_start:tdims_s%i_end,    & ! Aged biomass
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, bmass_cld(tdims_s%i_start:tdims_s%i_end,    & ! Biomass in cloud
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)

REAL,INTENT(INOUT) ::                         &
  ocff_new(tdims_s%i_start:tdims_s%i_end,     & ! New ocff
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, ocff_agd(tdims_s%i_start:tdims_s%i_end,     & ! Aged ocff
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, ocff_cld(tdims_s%i_start:tdims_s%i_end,     & ! Cloud ocff
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, nitr_acc(tdims_s%i_start:tdims_s%i_end,     & ! nitrate accumulation mode
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, nitr_diss(tdims_s%i_start:tdims_s%i_end,    & ! nitrate dissipation mode
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, co2(tdims_s%i_start:tdims_s%i_end,          & ! CO2
      tdims_s%j_start:tdims_s%j_end,          &
      tdims_s%k_start:tdims_s%k_end)

REAL,INTENT(INOUT) ::                                                  &
  free_tracers(tdims_s%i_start:tdims_s%i_end,                          &
               tdims_s%j_start:tdims_s%j_end,                          &
               tdims_s%k_start:tdims_s%k_end, tr_vars)                 &
, tracer_ukca(tdims_s%i_start:tdims_s%i_end,                           &
              tdims_s%j_start:tdims_s%j_end,                           &
              tdims_s%k_start:tdims_s%k_end, tr_ukca)                  &
, ozone_tracer( tdims_s%i_start:tdims_s%i_end,          & ! Ozone tracer
                tdims_s%j_start:tdims_s%j_end,          &
                tdims_s%k_start:tdims_s%k_end)

!----------------------------------------------------------------------
! Local variables
INTEGER ::                &
  i, k,                   & ! loop counters
  ntra_tmp                  ! sum of tracers copied so far

INTEGER ::                &
  errorstatus                ! Error reporting variable

CHARACTER (LEN=errormessagelength) ::                             &
  cmessage

CHARACTER(LEN=*), PARAMETER :: RoutineName='TRACER_RESTORE'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------
! Copy from tot_tracer

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k, i)
!$OMP SINGLE
ntra_tmp = 0
!$OMP END SINGLE

IF (l_dust) THEN
!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    dust_div1(1:row_length, 1:rows, k) =            &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO

!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    dust_div2(1:row_length, 1:rows, k) =            &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO
  IF (.NOT. l_twobin_dust) THEN

!$OMP SINGLE
    ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
    DO k=1, model_levels
      dust_div3(1:row_length, 1:rows, k) =            &
              tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
    END DO
!$OMP END DO

!$OMP SINGLE
    ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
    DO k=1, model_levels
      dust_div4(1:row_length, 1:rows, k) =            &
               tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
    END DO
!$OMP END DO

!$OMP SINGLE
    ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
    DO k=1, model_levels
      dust_div5(1:row_length, 1:rows, k) =            &
               tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
    END DO
!$OMP END DO

!$OMP SINGLE
    ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
    DO k=1, model_levels
      dust_div6(1:row_length, 1:rows, k) =            &
               tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
    END DO
!$OMP END DO
  END IF
END IF    ! l_dust

IF (l_sulpc_so2) THEN
!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    so2       (1:row_length, 1:rows, k) =           &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO

!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
        so4_aitken(1:row_length, 1:rows, k) =           &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO

!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    so4_accu  (1:row_length, 1:rows, k) =           &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO

!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    so4_diss  (1:row_length, 1:rows, k) =           &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO

  IF (l_sulpc_dms) THEN

!$OMP SINGLE
    ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
    DO k=1, model_levels
      dms       (1:row_length, 1:rows, k) =         &
               tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
    END DO
!$OMP END DO
  END IF

  IF (l_sulpc_nh3) THEN

!$OMP SINGLE
    ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
    DO k=1, model_levels
      nh3       (1:row_length, 1:rows, k) =         &
               tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
    END DO
!$OMP END DO
  END IF
END IF

IF (l_soot) THEN

!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    soot_new  (1:row_length, 1:rows, k) =           &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO

!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    soot_agd  (1:row_length, 1:rows, k) =           &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO

!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    soot_cld  (1:row_length, 1:rows, k) =           &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO

END IF

IF (l_biomass) THEN
!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    bmass_new  (1:row_length, 1:rows, k) =          &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO

!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    bmass_agd  (1:row_length, 1:rows, k) =          &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO

!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    bmass_cld  (1:row_length, 1:rows, k) =          &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO
END IF

IF (l_co2_interactive) THEN

!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    co2       (1:row_length, 1:rows, k) =           &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO
END IF

IF (l_ocff) THEN

!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    ocff_new  (1:row_length, 1:rows, k) =           &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO

!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    ocff_agd  (1:row_length, 1:rows, k) =           &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO

!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    ocff_cld  (1:row_length, 1:rows, k) =           &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO
END IF

IF (l_nitrate) THEN

!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    nitr_acc  (1:row_length, 1:rows, k) =           &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO

!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    nitr_diss  (1:row_length, 1:rows, k) =           &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO

END IF

IF (l_use_cariolle) THEN
!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    ozone_tracer (1:row_length, 1:rows, k) =        &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO
END IF

IF (l_murk_conv) THEN
!$OMP SINGLE
  ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
  DO k=1, model_levels
    aerosol   (1:row_length, 1:rows, k) =           &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
  END DO
!$OMP END DO
END IF

IF (tr_vars > 0 ) THEN

  i = 0

!$OMP DO SCHEDULE(STATIC)
  DO k= ntra_tmp+1 , ntra_tmp + tr_vars
    i = k - ntra_tmp
    free_tracers(1:row_length, 1:rows, 1:model_levels, i)=          &
             tot_tracer  (1:row_length, 1:rows, 1:model_levels, k)
  END DO
!$OMP END DO

!$OMP SINGLE
  ntra_tmp = ntra_tmp + tr_vars
!$OMP END SINGLE
END IF

IF (tr_ukca > 0 ) THEN

  i = 0

!$OMP DO SCHEDULE(STATIC)
  DO k= ntra_tmp+1, ntra_tmp + tr_ukca
    i = k - ntra_tmp
    tracer_ukca(1:row_length, 1:rows, 1:model_levels, i)=          &
             tot_tracer  (1:row_length, 1:rows, 1:model_levels, k)
  END DO
!$OMP END DO

!$OMP SINGLE
  ntra_tmp = ntra_tmp + tr_ukca
!$OMP END SINGLE

END IF

!$OMP END PARALLEL

! Call routine to capture the change of tracers due to plume scavenging
! for removal diagnostics
IF (tr_ukca > 0 .AND. l_ukca_mode .AND. l_ukca_plume_diags .AND.   &
    i_convection_vn == i_convection_vn_6a) THEN

  ! Set final value of tracers
  CALL ukca_plume_scav_final(tr_ukca,tracer_ukca)

END IF

!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!----------------------------------------------------------------------
RETURN
END SUBROUTINE tracer_restore

END MODULE tracer_restore_mod
