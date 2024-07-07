! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE populate_ancil_requests_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='POPULATE_ANCIL_REQUESTS_MOD'

!  Subroutine populate_ancil_requests

! Description:
!    Takes ancil requests from items namelist and sort them into the
!    ancil_request type in the order expected by the code and as previously 
!    mandated by the ancilmaster.  Also add each ancil file encountered  
!    to the array of ancil_files.
!
! ancil_order below is an array of ordered stash codes as was defined for the
! first 27 items within the retired ancilmaster, plus a few extra related to
! surface fields, to ensure surf temp and sst are calculated in determined order
!
! Method:
!    Simple sort based upon a known order
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Ancillaries
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.

CONTAINS

SUBROUTINE populate_ancil_requests(stashcode,         &
                                   section,           &
                                   item,              &
                                   period,            &
                                   interval,          &
                                   anc_filename,      &
                                   l_ancil_requested)

USE ancil_mod, ONLY: ancil_requests,  ancil_files,     &
    num_ancil_requests, num_ancil_files, add_ancil_request

USE filenamelength_mod, ONLY: filenamelength

USE um_stashcode_mod, ONLY:                                    &
    stashcode_lsm,            stashcode_orog,                     &
    stashcode_orog_var,       stashcode_orog_gdxx,                &
    stashcode_orog_gdxy,      stashcode_orog_gdyy,                &
    stashcode_ozone,          stashcode_mean_snow,                &
    stashcode_soil_temp,      stashcode_vol_smc_wilt,             &
    stashcode_vol_smc_cri,    stashcode_vol_smc_sat,              &
    stashcode_Ksat,           stashcode_thermal_capacity,         &
    stashcode_thermal_conduct,stashcode_veg_frac,                 &
    stashcode_z0,             stashcode_icefrac,                  &
    stashcode_tstar,          stashcode_icethick,                 &
    stashcode_surf_z_curr,    stashcode_surf_m_curr,              &
    stashcode_runoff_coast_out,                                   &
    stashcode_Ti_Mean,        stashcode_Ti_Sig,                   &
    stashcode_soil_suction,   stashcode_soil_moist,               &
    stashcode_SO2_emiss,      stashcode_dimethyl_sul_emiss,       &
    stashcode_total_aero_emiss,                                   &
    stashcode_total_aero,     stashcode_sil_orog_rough,           &
    stashcode_hlf_pk_to_trf_ht,                                   &
    stashcode_clapp_hb,       stashcode_frac_surf_type,           &
    stashcode_LAI,            stashcode_canopy_height,            &
    stashcode_can_conduct,    stashcode_disturb_frac_veg,         &
    stashcode_snw_free_alb_bs,stashcode_soil_carbon_cont,         &
    stashcode_z0m_soil,                                           &
    stashcode_land_frac,      stashcode_dust_parent_clay,         &
    stashcode_dust_parent_silt, stashcode_dust_parent_sand,       &
    stashcode_dust_soil_mf1,  stashcode_dust_soil_mf2,            &
    stashcode_dust_soil_mf3,  stashcode_dust_soil_mf4,            &
    stashcode_dust_soil_mf5,  stashcode_dust_soil_mf6,            &
    stashcode_biom_surf_em,   stashcode_biom_elev_em,             &
    stashcode_biom_elev_em_h1, stashcode_biom_elev_em_h2,         &
    stashcode_riv_storage,    stashcode_riv_sequence,             &
    stashcode_riv_direction,  stashcode_orog_x_grad,              &
    stashcode_orog_y_grad,    stashcode_unfilt_orog,              &
    stashcode_snowdep_grd_tile, stashcode_urbhgt,                 &
    stashcode_urbhwr,         stashcode_urbwrr,                   &
    stashcode_flake_depth

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments

! Array of stashcodes read from items namelist
INTEGER, INTENT(IN) :: stashcode(:)   
! Array of section numbers read from items namelist
INTEGER, INTENT(IN) :: section(:)   
! Array of item numbers read from items namelist
INTEGER, INTENT(IN) :: item(:)   
! Array of updating period read in from items namelist
INTEGER, INTENT(IN) :: period(:)
! Array of updating interval read in from items namelist
INTEGER, INTENT(IN) :: interval(:)
! Array of filenames read from items namelist
CHARACTER (LEN=filenamelength), INTENT(IN) :: anc_filename(:)
! Arrays of source values read from items namelist
LOGICAL, INTENT(IN) :: l_ancil_requested(:)   


! Local variables
INTEGER :: num_items
INTEGER, PARAMETER :: num_ordered_stash=67
! ordered stash codes
INTEGER            :: ancil_stash_order (num_ordered_stash) =   &   
         (/    stashcode_lsm,            stashcode_orog,                     &
               stashcode_orog_var,       stashcode_orog_gdxx,                &
               stashcode_orog_gdxy,      stashcode_orog_gdyy,                &
               stashcode_ozone,          stashcode_mean_snow,                &
               stashcode_soil_temp,      stashcode_vol_smc_wilt,             &
               stashcode_vol_smc_cri,    stashcode_vol_smc_sat,              &
               stashcode_Ksat,           stashcode_thermal_capacity,         &
               stashcode_thermal_conduct,stashcode_veg_frac,                 &
               stashcode_z0,             stashcode_icefrac,                  &
               stashcode_tstar,          stashcode_icethick,                 &
               stashcode_surf_z_curr,    stashcode_surf_m_curr,              &
               stashcode_runoff_coast_out,                                   &
               stashcode_Ti_Mean,        stashcode_Ti_Sig,                   &
               stashcode_soil_suction,   stashcode_soil_moist,               &
               stashcode_SO2_emiss,      stashcode_dimethyl_sul_emiss,       &
               stashcode_total_aero_emiss,                                   &
               stashcode_total_aero,     stashcode_sil_orog_rough,           &
               stashcode_hlf_pk_to_trf_ht,                                   &
               stashcode_clapp_hb,       stashcode_frac_surf_type,           &
               stashcode_LAI,            stashcode_canopy_height,            &
               stashcode_can_conduct,    stashcode_disturb_frac_veg,         &
               stashcode_snw_free_alb_bs,stashcode_soil_carbon_cont,         &
               stashcode_z0m_soil,                                           &
               stashcode_land_frac,      stashcode_dust_parent_clay,         &
               stashcode_dust_parent_silt, stashcode_dust_parent_sand,       &
               stashcode_dust_soil_mf1,  stashcode_dust_soil_mf2,            &
               stashcode_dust_soil_mf3,  stashcode_dust_soil_mf4,            &
               stashcode_dust_soil_mf5,  stashcode_dust_soil_mf6,            &
               stashcode_biom_surf_em,   stashcode_biom_elev_em,             &
               stashcode_biom_elev_em_h1, stashcode_biom_elev_em_h2,         &
               stashcode_riv_storage,    stashcode_riv_sequence,             &
               stashcode_riv_direction,  stashcode_orog_x_grad,              &
               stashcode_orog_y_grad,    stashcode_unfilt_orog,              &
               stashcode_urbhgt,         stashcode_urbhwr,                   &
               stashcode_urbwrr,         stashcode_snowdep_grd_tile,         &
               stashcode_flake_depth /)

INTEGER            :: i,j

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'POPULATE_ANCIL_REQUESTS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

num_items = SIZE(stashcode)

! Loop over list of ordered stash codes
DO i = 1, num_ordered_stash
  ! Check if any items requests match and if they are due to be configured
  ! from ancillary
  DO j = 1, num_items
    IF ( ancil_stash_order(i) == stashcode(j) .AND. &
        l_ancil_requested(j)) THEN
      CALL add_ancil_request(stashcode(j), section(j), item(j), period(j), &
           interval(j), anc_filename(j))
    END IF ! Stashcode check
  END DO ! End loop over items requests
END DO ! End loop over ordered stashcodes

! Append any remaining fields which were not in the ordered list of
! stashcodes
DO j = 1, num_items
  ! Only proceed if the items namelist stashcode is not already
  ! present in the ancil_request type
  IF ( (.NOT. ANY(ancil_requests%stashcode == stashcode(j))) .AND. &
       (l_ancil_requested(j))) THEN
    CALL add_ancil_request(stashcode(j), section(j), item(j), period(j), &
         interval(j), anc_filename(j))
  END IF ! Stashcode check
END DO ! End loop over items requests
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE populate_ancil_requests

END MODULE populate_ancil_requests_mod


