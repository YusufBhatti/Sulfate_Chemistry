! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! All the calls to scavenging of classic aerosols and murk

MODULE all_scav_calls_mod

IMPLICIT NONE
!---------------------------------------------------------------------------
! Description: All the calls to scavenging of various tracers by convection

! method: Most of the calls are for CLASSIC aerosols but there is also a call
!         for murk.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection

! code description:
!   language: fortran 90
!   This code is written to umdp3 programming standards version 10.4
!---------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ALL_SCAV_CALLS_MOD'

CONTAINS

SUBROUTINE  all_scav_calls(row_length, rows, ccb, cct,               &
    conv_rain, conv_snow, rho, q_star, qcl_star, qcf_star,           &
    aerosol,                                                         &
    dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6,     &
    so2, so4_aitken, so4_accu, so4_diss,                             &
    dms, nh3, soot_agd, bmass_agd, ocff_agd, nitr_acc, nitr_diss  )

! Classic aerosol flags
USE run_aerosol_mod, ONLY: l_sulpc_so2, l_soot, l_ocff, l_biomass,  &
                           l_sulpc_nh3, l_nitrate
USE dust_parameters_mod, ONLY:  l_dust, krain_dust, ksnow_dust,          &
  l_twobin_dust

USE murk_inputs_mod, ONLY: l_murk_source

USE c_bm_con_mod, ONLY: krain_agedbmass, ksnow_agedbmass
USE c_ocff_con_mod, ONLY: krain_agedocff, ksnow_agedocff
USE c_st_con_mod, ONLY: krain_agedsoot, ksnow_agedsoot
USE c_nitrate_con_mod, ONLY: krain_nitracc, ksnow_nitracc, krain_nitrdiss, &
                             ksnow_nitrdiss

USE atm_fields_bounds_mod, ONLY:                                          &
   tdims_s, pdims_s, tdims

USE timestep_mod, ONLY: timestep

! Convection output diagnostics
USE cv_diagnostic_array_mod, ONLY:                                      &
        conscav_dust ,conscav_so4ait ,conscav_so4acc ,conscav_so4dis    &
       ,conscav_agedsoot ,conscav_agedbmass ,conscav_agedocff           &
       ,conscav_nitracc ,conscav_nitrdiss ,conwash_so2 ,conwash_nh3

! Subroutines
USE ncnwsh2_mod, ONLY: ncnwsh2
USE scnscv2_mod, ONLY: scnscv2
USE scnwsh2_mod, ONLY: scnwsh2

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

USE con_scav_mod, ONLY: con_scav
IMPLICIT NONE

!----------------------------------------------------------------------
! Subroutine arguments
INTEGER, INTENT(IN) ::       &
  row_length                 & ! Row length
, rows                         ! Number of rows

INTEGER, INTENT(IN) ::     &!
  ccb    (row_length,rows) &! Cnv Cld Base of highest layer on gridpoint
, cct    (row_length,rows)  ! Cnv Cld Top  of highest layer on gridpoint

REAL, INTENT(IN) ::             &
  conv_rain(row_length, rows)   &  ! Surface convective rainfall rate (kg/m2/s)
, conv_snow(row_length, rows)      ! Surface convective snowfall rate (kg/m2/s)

REAL, INTENT(IN) ::                            &
  rho(pdims_s%i_start:pdims_s%i_end,           & ! density *r*r (kg/m)
      pdims_s%j_start:pdims_s%j_end,           & !
      pdims_s%k_start:pdims_s%k_end)           & !
, q_star(tdims%i_start:tdims%i_end,            & ! water vapour (kg/kg)
         tdims%j_start:tdims%j_end,            &
                     1:tdims%k_end)            &
, qcl_star(tdims%i_start:tdims%i_end,          & ! cloud liquid water (kg/kg)
           tdims%j_start:tdims%j_end,          &
                       1:tdims%k_end)          &
, qcf_star(tdims%i_start:tdims%i_end,          & ! cloud ice water (kg/kg)
           tdims%j_start:tdims%j_end,          &
                       1:tdims%k_end)

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
, soot_agd(tdims_s%i_start:tdims_s%i_end,     & ! Aged soot
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, bmass_agd(tdims_s%i_start:tdims_s%i_end,    & ! Aged biomass
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)

REAL,INTENT(INOUT) ::                         &
  ocff_agd(tdims_s%i_start:tdims_s%i_end,     & ! Aged ocff
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, nitr_acc(tdims_s%i_start:tdims_s%i_end,     & ! nitrate accumulation mode
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, nitr_diss(tdims_s%i_start:tdims_s%i_end,    & ! nitrate dissipation mode
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)

!----------------------------------------------------------------------
! Local variables
LOGICAL :: l_scav_below_cloud   !control for scavenging levels

INTEGER ::                &
  ntra_tmp                  ! sum of tracers copied so far

INTEGER ::                &
  errorstatus                ! Error reporting variable

! Parameters for S Cycle tracer scavenging

REAL, PARAMETER ::                                                &
    krain_so4ait=0.3e-4                                           &
,   ksnow_so4ait=0.3e-4                                           &
,   krain_so4acc=0.3e-4                                           &
,   ksnow_so4acc=0.3e-4                                           &
,   krain_so4dis=0.3e-4                                           &
,   ksnow_so4dis=0.3e-4


CHARACTER (LEN=errormessagelength) ::                             &
  cmessage

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALL_SCAV_CALLS'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------
! Scavenging for aerosol murk (Not part of section 17: aerosol scheme)
IF (l_murk_source) THEN
  CALL con_scav( timestep                                           &
  ,    ccb, cct, conv_rain, conv_snow                               &
  ,    aerosol(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,&
                           1:tdims%k_end) )
END IF

! Section 17 : Aerosol scheme tracers
! Scavenging of other tracers should go here....
! Scavenge Mineral Dust tracers

IF (l_dust) THEN

  l_scav_below_cloud = .TRUE.

  CALL scnscv2( timestep,                                         &
           rho, q_star, qcl_star, qcf_star, dust_div1,            &
           Ccb, cct, conv_rain, conv_snow,                        &
           l_scav_below_cloud, krain_dust(1), ksnow_dust(1),      &
           Conscav_dust(1,1,1) )

  CALL scnscv2( timestep,                                         &
           rho, q_star, qcl_star, qcf_star, dust_div2,            &
           Ccb, cct, conv_rain, conv_snow,                        &
           l_scav_below_cloud, krain_dust(2), ksnow_dust(2),      &
           Conscav_dust(1,1,2) )

  ! Only do these when using six-bin dust
  IF (.NOT. l_twobin_dust) THEN
    CALL scnscv2( timestep,                                       &
             rho, q_star, qcl_star, qcf_star, dust_div3,          &
             Ccb, cct, conv_rain, conv_snow,                      &
             l_scav_below_cloud, krain_dust(3), ksnow_dust(3),    &
             Conscav_dust(1,1,3) )

    CALL scnscv2( timestep,                                       &
             rho, q_star, qcl_star, qcf_star, dust_div4,          &
             Ccb, cct, conv_rain, conv_snow,                      &
             l_scav_below_cloud, krain_dust(4), ksnow_dust(4),    &
             Conscav_dust(1,1,4) )

    CALL scnscv2( timestep,                                       &
             rho, q_star, qcl_star, qcf_star, dust_div5,          &
             Ccb, cct, conv_rain, conv_snow,                      &
             l_scav_below_cloud, krain_dust(5), ksnow_dust(5),    &
             Conscav_dust(1,1,5) )

    CALL scnscv2( timestep,                                       &
             rho, q_star, qcl_star, qcf_star, dust_div6,          &
             Ccb, cct, conv_rain, conv_snow,                      &
             l_scav_below_cloud, krain_dust(6), ksnow_dust(6),    &
             Conscav_dust(1,1,6) )
  END IF
END IF !L_DUST

! Scavenge Sulphur Cycle tracers

IF (l_sulpc_so2) THEN

  l_scav_below_cloud = .FALSE.

  ! Scavenge SO2
  CALL scnwsh2( timestep,                                         &
           rho, q_star, qcl_star, qcf_star, so2,                  &
           ccb, cct,                                              &
           conv_rain, conwash_so2 )

  ! Scavenge NH3 if present

  IF (l_sulpc_nh3) THEN

    ! Scavenge NH3
    CALL ncnwsh2( timestep,                                       &
           rho, q_star, qcl_star, qcf_star, nh3,                  &
           ccb, cct,                                              &
           conv_rain, conwash_nh3 )

  END IF          !End L_sulpc_nh3 condition

  ! Scavenge SO4_AIT
  CALL scnscv2( timestep,                                         &
           rho, q_star, qcl_star, qcf_star, so4_aitken,           &
           ccb, cct, conv_rain, conv_snow,                        &
           l_scav_below_cloud, krain_so4ait, ksnow_so4ait,        &
           conscav_so4ait )

  ! Scavenge SO4_ACC
  CALL scnscv2( timestep,                                         &
           rho, q_star, qcl_star, qcf_star, so4_accu,             &
           ccb, cct, conv_rain, conv_snow,                        &
           l_scav_below_cloud, krain_so4acc, ksnow_so4acc,        &
           conscav_so4acc )

  ! Scavenge SO4_DIS
  CALL scnscv2( timestep,                                         &
           rho, q_star, qcl_star, qcf_star, so4_diss,             &
           ccb, cct, conv_rain, conv_snow,                        &
           l_scav_below_cloud, krain_so4dis, ksnow_so4dis,        &
           conscav_so4dis )

END IF              !End L_sulpc_so2 condition


IF (l_soot) THEN  ! If soot modelling is included

  l_scav_below_cloud = .FALSE.

  ! Scavenge aged soot
  CALL scnscv2( timestep,                                         &
           rho, q_star, qcl_star, qcf_star, soot_agd,             &
           ccb, cct, conv_rain, conv_snow,                        &
           l_scav_below_cloud, krain_agedsoot, ksnow_agedsoot,    &
           conscav_agedsoot )

END IF  ! l_soot

IF (l_biomass) THEN  ! If biomass aerosol is included

  l_scav_below_cloud = .FALSE.

  ! Scavenge aged biomass aerosol
  CALL scnscv2( timestep,                                         &
           rho, q_star, qcl_star, qcf_star, bmass_agd,            &
           ccb, cct, conv_rain, conv_snow,                        &
           l_scav_below_cloud, krain_agedbmass, ksnow_agedbmass,  &
           conscav_agedbmass )

END IF  ! l_biomass

IF (l_ocff) THEN ! If fossil-fuel OC aerosol is included

  l_scav_below_cloud = .FALSE.

  ! Scavenge aged fossil-fuel organic carbon aerosol
  CALL scnscv2( timestep,                                         &
           rho, q_star, qcl_star, qcf_star, ocff_agd,             &
           ccb, cct, conv_rain, conv_snow,                        &
           l_scav_below_cloud, krain_agedocff, ksnow_agedocff,    &
           conscav_agedocff )

END IF  ! l_ocff

IF (l_nitrate) THEN ! If ammonium nitrate is included

  l_scav_below_cloud = .FALSE.

  ! Scavenge accumulation-mode ammonium nitrate
  CALL scnscv2( timestep,                                         &
           rho, q_star, qcl_star, qcf_star, nitr_acc,             &
           ccb, cct, conv_rain, conv_snow,                        &
           l_scav_below_cloud, krain_nitracc, ksnow_nitracc,      &
           conscav_nitracc )

  ! Scavenge dissolved ammonium nitrate
  CALL scnscv2( timestep,                                         &
           rho, q_star, qcl_star, qcf_star, nitr_diss,            &
           ccb, cct, conv_rain, conv_snow,                        &
           l_scav_below_cloud, krain_nitrdiss, ksnow_nitrdiss,    &
           conscav_nitrdiss )

END IF  ! l_nitrate
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!----------------------------------------------------------------------
RETURN
END SUBROUTINE all_scav_calls

END MODULE all_scav_calls_mod
