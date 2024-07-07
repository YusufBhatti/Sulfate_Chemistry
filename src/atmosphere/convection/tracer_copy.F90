! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Copy the tracers to array tot_tracer

MODULE tracer_copy_mod

IMPLICIT NONE
!---------------------------------------------------------------------------
! Description: Copies all required tracers into tot_tracer

! method:

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection

! code description:
!   language: fortran 90
!   This code is written to umdp3 programming standards version 10.4
!---------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TRACER_COPY_MOD'

CONTAINS

SUBROUTINE  tracer_copy(row_length, rows, ntra_fld, ntra_lev      &
, tr_vars, tr_ukca,  aerosol                                      &
, dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6     &
, so2, so4_aitken, so4_accu, so4_diss                             &
, dms, nh3, soot_new, soot_agd, soot_cld, bmass_new, bmass_agd    &
, bmass_cld, ocff_new, ocff_agd, ocff_cld, nitr_acc, nitr_diss    &
, co2,  free_tracers, tracer_ukca                                 &
, ozone_tracer                                                    &
, tot_tracer )

USE ukca_option_mod,  ONLY: l_ukca_mode
USE ukca_d1_defs,     ONLY: l_ukca_plume_diags
USE ukca_scavenging_diags_mod, ONLY: ukca_plume_scav_initial
USE ukca_scavenging_mod, ONLY:  tracer_info

! Classic aerosol flags
USE run_aerosol_mod, ONLY: l_sulpc_so2, l_soot, l_ocff, l_biomass,  &
                           l_sulpc_nh3, l_nitrate, l_sulpc_dms

USE dust_parameters_mod, ONLY:  l_dust, l_twobin_dust

USE carbon_options_mod, ONLY: l_co2_interactive     

USE rad_input_mod, ONLY: l_use_cariolle
USE nlsizes_namelist_mod,   ONLY: model_levels

USE cv_run_mod,  ONLY:  l_murk_conv,                                &
                        i_convection_vn,                            &
                        i_convection_vn_6a

USE atm_fields_bounds_mod, ONLY: tdims_s

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

! Tracers
REAL,INTENT(IN) ::                         &
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
REAL,INTENT(IN) ::                         &
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

REAL,INTENT(IN) ::                         &
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

REAL,INTENT(IN) ::                                                     &
  free_tracers(tdims_s%i_start:tdims_s%i_end,                          &
               tdims_s%j_start:tdims_s%j_end,                          &
               tdims_s%k_start:tdims_s%k_end, tr_vars)                 &
, tracer_ukca(tdims_s%i_start:tdims_s%i_end,                           &
              tdims_s%j_start:tdims_s%j_end,                           &
              tdims_s%k_start:tdims_s%k_end, tr_ukca)                  &
, ozone_tracer( tdims_s%i_start:tdims_s%i_end,          & ! Ozone tracer
                tdims_s%j_start:tdims_s%j_end,          &
                tdims_s%k_start:tdims_s%k_end)

REAL,INTENT(INOUT) ::                                                     &
  tot_tracer(row_length, rows, ntra_lev, ntra_fld)
!----------------------------------------------------------------------
! Local variables
INTEGER ::                &
  ntra_tmp                  ! sum of tracers copied so far

INTEGER ::                &
  errorstatus                ! Error reporting variable

CHARACTER (LEN=errormessagelength) ::                             &
  cmessage

CHARACTER(LEN=*), PARAMETER :: RoutineName='TRACER_COPY'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------

! Call routine for plume scavenging diagnostics
IF (tr_ukca > 0 .AND. l_ukca_mode .AND. l_ukca_plume_diags .AND.  &
    i_convection_vn == i_convection_vn_6a) THEN

  ! Set initial value of tracers
  CALL ukca_plume_scav_initial(tr_ukca, tracer_ukca)

END IF


! copy arrays into tot_tracer

IF ( l_sulpc_so2        .OR.  l_soot          .OR.  &
     l_dust             .OR.  l_biomass       .OR.  &
     l_co2_interactive  .OR.  l_use_cariolle  .OR.  &
     l_ocff             .OR.  l_nitrate       .OR.  &
     l_murk_conv ) THEN

  ntra_tmp = 0

  IF (l_dust) THEN
    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       dust_div1(1:row_length, 1:rows, 1:model_levels)

    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       dust_div2(1:row_length, 1:rows, 1:model_levels)

    IF (.NOT. l_twobin_dust) THEN
      ntra_tmp = ntra_tmp + 1
      tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
         dust_div3(1:row_length, 1:rows, 1:model_levels)

      ntra_tmp = ntra_tmp + 1
      tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
         dust_div4(1:row_length, 1:rows, 1:model_levels)

      ntra_tmp = ntra_tmp + 1
      tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
         dust_div5(1:row_length, 1:rows, 1:model_levels)

      ntra_tmp = ntra_tmp + 1
      tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
         dust_div6(1:row_length, 1:rows, 1:model_levels)
    END IF
  END IF

  IF (l_sulpc_so2) THEN
    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       so2(1:row_length, 1:rows, 1:model_levels)

    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       so4_aitken(1:row_length, 1:rows, 1:model_levels)

    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       so4_accu(1:row_length, 1:rows, 1:model_levels)

    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       so4_diss(1:row_length, 1:rows, 1:model_levels)

    IF (l_sulpc_dms) THEN
      ntra_tmp = ntra_tmp + 1
      tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)=&
         dms(1:row_length, 1:rows, 1:model_levels)
    END IF

    IF (l_sulpc_nh3) THEN
      ntra_tmp = ntra_tmp + 1
      tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)=&
         nh3(1:row_length, 1:rows, 1:model_levels)
    END IF
  END IF   ! l_sulpc_so2

  IF (l_soot) THEN
    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       soot_new(1:row_length, 1:rows, 1:model_levels)

    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       soot_agd(1:row_length, 1:rows, 1:model_levels)

    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       soot_cld(1:row_length, 1:rows, 1:model_levels)
  END IF

  IF (l_biomass) THEN
    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       bmass_new(1:row_length, 1:rows, 1:model_levels)

    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       bmass_agd(1:row_length, 1:rows, 1:model_levels)

    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       bmass_cld(1:row_length, 1:rows, 1:model_levels)
  END IF

  IF (l_co2_interactive) THEN
    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       co2(1:row_length, 1:rows, 1:model_levels)
  END IF

  IF (l_ocff) THEN
    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       ocff_new(1:row_length, 1:rows, 1:model_levels)

    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       ocff_agd(1:row_length, 1:rows, 1:model_levels)

    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       ocff_cld(1:row_length, 1:rows, 1:model_levels)
  END IF

  IF (l_nitrate) THEN
    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       nitr_acc(1:row_length, 1:rows, 1:model_levels)

    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
       nitr_diss(1:row_length, 1:rows, 1:model_levels)
  END IF

  IF (l_use_cariolle) THEN
    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
          ozone_tracer(1:row_length, 1:rows, 1:model_levels)
  END IF

  IF (l_murk_conv) THEN
    ntra_tmp = ntra_tmp + 1
    tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
           aerosol(1:row_length, 1:rows, 1:model_levels)
  END IF

  IF ( tr_vars > 0 ) THEN
    tot_tracer(1:row_length, 1:rows, 1:tdims_s%k_end,&
               ntra_tmp+1:ntra_tmp+tr_vars) =&
       free_tracers(1:row_length, 1:rows, 1:tdims_s%k_end, 1:tr_vars)

     ntra_tmp = ntra_tmp + tr_vars
  END IF

  IF ( tr_ukca > 0 ) THEN
    tracer_info%i_ukca_first = ntra_tmp+1
    tot_tracer(1:row_length, 1:rows, 1:tdims_s%k_end,&
               ntra_tmp+1:ntra_tmp+tr_ukca) =&
         tracer_ukca(1:row_length, 1:rows, 1:tdims_s%k_end, 1:tr_ukca)

    ntra_tmp = ntra_tmp + tr_ukca
    tracer_info%i_ukca_last = ntra_tmp
  END IF

ELSE ! no soot, co2, sulphur, dust, biomass, OCFF, nitrate
     ! or murk
  ntra_tmp = 0
  IF ( tr_vars > 0 ) THEN
    tot_tracer  (1:row_length, 1:rows, 1:tdims_s%k_end, 1:tr_vars) = &
         free_tracers(1:row_length, 1:rows, 1:tdims_s%k_end, 1:tr_vars)
    ntra_tmp = tr_vars
  END IF

  IF ( tr_ukca > 0 ) THEN
    tracer_info%i_ukca_first = ntra_tmp+1
    tot_tracer(1:row_length, 1:rows, 1:tdims_s%k_end,&
    ntra_tmp+1:ntra_tmp+tr_ukca) = &
       tracer_ukca(1:row_length, 1:rows, 1:tdims_s%k_end,1:tr_ukca)
    ntra_tmp = ntra_tmp + tr_ukca
    tracer_info%i_ukca_last = ntra_tmp
  END IF

  IF ( tr_vars + tr_ukca == 0) THEN ! make sure the value is set.
    tot_tracer(:,:,:,:) = 0.0
  END IF
END IF !soot or co2, sulphur, dust, biomass, cariolle o3, ocff, nitrate or murk
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!----------------------------------------------------------------------
RETURN
END SUBROUTINE tracer_copy

END MODULE tracer_copy_mod
