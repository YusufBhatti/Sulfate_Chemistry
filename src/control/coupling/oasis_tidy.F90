#if defined(MCT)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis_tidy()
  !
  ! Description: Tidy up at the end of an OASIS-based coupled run.
  !              Expected to be applicable to all OASIS versions
  !              including OASIS3-MCT.
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: Coupling
  !
  !--------------------------------------------------------------------

USE oasis_atm_data_mod, ONLY: max_transients, transient_o2a,        & 
  transient_o2a_count, transient_a2o, transient_a2o_count,          &
  taux, tauy, solar2d, blue2d, evap2d, longwave2d, sensible2d,      &
  pond_frac_n, pond_depth_n, heatflux, sublim, rainls, snowls,      &
  rainconv, snowconv, totalrain, totalsnow, riverout, calving,      &
  w10, topmeltn, fcondtopn, mslp, icenorth, icesouth, pCO2_out,     &
  dust_dep_out,                                                     &
  ocn_hice, ocn_hicen, ocn_freeze, ocn_freezen, ocn_freezen_first,  &
  ice_fract_cat_future, ocn_snowthick, ocn_snowthickn, ocn_sst,     & 
  ocn_sstfrz, ocn_sst_orig, ocn_u, ocn_v, ocn_icetn, ocn_icet_gb,   & 
  ocn_icekn, tstar_sicen, tstar_local, tstar_ssi, tstar_sice_local, &
  tstar_sice_cat_local, tstar_land_local, fland_loc, dms_conc_in,   &
  CO2flux_in, chloro_conc_in

USE um_types, ONLY: integer32

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(KIND=integer32) :: tc  ! transient field counter 

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='OASIS_TIDY'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Nullify coupling field pointers prior to deallocation of targets.
DO tc = 1, transient_o2a_count
  IF (ASSOCIATED(transient_o2a(tc)%field)) NULLIFY (transient_o2a(tc)%field)
END DO ! Over tc

DO tc = 1, transient_a2o_count
  IF (ASSOCIATED(transient_a2o(tc)%field)) NULLIFY (transient_a2o(tc)%field)
END DO ! Over tc

! Deallocate arrays used in the processing of
! outgoing coupling data
IF (ALLOCATED(taux)) DEALLOCATE(taux)
IF (ALLOCATED(tauy)) DEALLOCATE(tauy)
IF (ALLOCATED(solar2d)) DEALLOCATE(solar2d)
IF (ALLOCATED(blue2d)) DEALLOCATE(blue2d)
IF (ALLOCATED(evap2d)) DEALLOCATE(evap2d)
IF (ALLOCATED(longwave2d)) DEALLOCATE(longwave2d)
IF (ALLOCATED(sensible2d)) DEALLOCATE(sensible2d)
IF (ALLOCATED(pond_frac_n)) DEALLOCATE(pond_frac_n)
IF (ALLOCATED(pond_depth_n)) DEALLOCATE(pond_depth_n)
IF (ALLOCATED(heatflux)) DEALLOCATE(heatflux)
IF (ALLOCATED(sublim)) DEALLOCATE(sublim)
IF (ALLOCATED(rainls)) DEALLOCATE(rainls)
IF (ALLOCATED(snowls)) DEALLOCATE(snowls)
IF (ALLOCATED(rainconv)) DEALLOCATE(rainconv)
IF (ALLOCATED(snowconv)) DEALLOCATE(snowconv)
IF (ALLOCATED(totalrain)) DEALLOCATE(totalrain)
IF (ALLOCATED(totalsnow)) DEALLOCATE(totalsnow)
IF (ALLOCATED(riverout)) DEALLOCATE(riverout)
IF (ALLOCATED(calving)) DEALLOCATE(calving)
IF (ALLOCATED(w10)) DEALLOCATE(w10)
IF (ALLOCATED(topmeltn)) DEALLOCATE(topmeltn)
IF (ALLOCATED(fcondtopn)) DEALLOCATE(fcondtopn)
IF (ALLOCATED(mslp)) DEALLOCATE(mslp)
IF (ALLOCATED(icenorth)) DEALLOCATE(icenorth)
IF (ALLOCATED(icesouth)) DEALLOCATE(icesouth)
IF (ALLOCATED(pCO2_out)) DEALLOCATE(pCO2_out)
IF (ALLOCATED(dust_dep_out)) DEALLOCATE(dust_dep_out)


! Deallocate arrays used in the processing of
! incoming coupling data
IF (ALLOCATED(ocn_hice)) DEALLOCATE(ocn_hice)
IF (ALLOCATED(ocn_hicen)) DEALLOCATE(ocn_hicen)
IF (ALLOCATED(ocn_freeze)) DEALLOCATE(ocn_freeze)
IF (ALLOCATED(ocn_freezen)) DEALLOCATE(ocn_freezen)
IF (ALLOCATED(ocn_freezen_first)) DEALLOCATE(ocn_freezen_first)
IF (ALLOCATED(ice_fract_cat_future)) DEALLOCATE(ice_fract_cat_future)
IF (ALLOCATED(ocn_snowthick)) DEALLOCATE(ocn_snowthick)
IF (ALLOCATED(ocn_snowthickn)) DEALLOCATE(ocn_snowthickn)
IF (ALLOCATED(ocn_sst)) DEALLOCATE(ocn_sst)
IF (ALLOCATED(ocn_sstfrz)) DEALLOCATE(ocn_sstfrz)
IF (ALLOCATED(ocn_sst_orig)) DEALLOCATE(ocn_sst_orig)
IF (ALLOCATED(ocn_u)) DEALLOCATE(ocn_u)
IF (ALLOCATED(ocn_v)) DEALLOCATE(ocn_v)
IF (ALLOCATED(ocn_icetn)) DEALLOCATE(ocn_icetn)
IF (ALLOCATED(ocn_icet_gb)) DEALLOCATE(ocn_icet_gb)
IF (ALLOCATED(ocn_icekn)) DEALLOCATE(ocn_icekn)
IF (ALLOCATED(tstar_sicen)) DEALLOCATE(tstar_sicen)
IF (ALLOCATED(tstar_local)) DEALLOCATE(tstar_local)
IF (ALLOCATED(tstar_ssi)) DEALLOCATE(tstar_ssi)
IF (ALLOCATED(tstar_sice_local)) DEALLOCATE(tstar_sice_local)
IF (ALLOCATED(tstar_sice_cat_local)) DEALLOCATE(tstar_sice_cat_local)
IF (ALLOCATED(tstar_land_local)) DEALLOCATE(tstar_land_local)
IF (ALLOCATED(fland_loc)) DEALLOCATE(fland_loc)
IF (ALLOCATED(dms_conc_in)) DEALLOCATE(dms_conc_in)
IF (ALLOCATED(CO2flux_in)) DEALLOCATE(CO2flux_in)
IF (ALLOCATED(chloro_conc_in)) DEALLOCATE(chloro_conc_in)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN


END SUBROUTINE oasis_tidy
#endif
