#if defined(MCT)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis_updatecpl(cmessage)
  !
  ! Description:
  ! Update coupling field data into prognostic fields
  ! to ensure restartability. We do this simply by copying
  ! data from diagnostic (stash) areas at the appropriate timestep.
  ! Applicable to OASIS3-MCT and probably to all other couplers.
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: Coupling
  !
  !-------------------------------------------------------------------

USE oasis_atm_data_mod

USE atm_fields_bounds_mod, ONLY: pdims, wdims, tdims,                   &
                                 udims, vdims, wdims_s,                 &
                                 o3dims2

USE jules_snow_mod, ONLY: nsmax

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE jules_sea_seaice_mod, ONLY: nice, nice_use, l_ctile

USE UM_ParVars
USE Control_Max_Sizes
USE Decomp_DB
USE lbc_mod
USE coupling_control_mod,  ONLY: oasis_couple_freq_ao
USE free_tracers_inputs_mod, ONLY: a_tracer_last, a_tracer_first
USE nlsizes_namelist_mod, ONLY:                                         &
    land_field, len_tot, model_levels, n_cca_lev, n_obj_d1_max, n_rows, &
    ntiles, river_row_length, river_rows, row_length,                   &
    rows, sm_levels, st_levels, theta_off_size, tpps_ozone_levels,      &
    tr_lbc_ukca, tr_lbc_vars, tr_levels, tr_ukca, tr_vars
USE cv_run_mod, ONLY: l_param_conv
USE d1_array_mod, ONLY: d1
USE atm_fields_mod, ONLY: &
    c_solar, c_blue, c_longwave, c_taux, c_tauy, c_w10, c_sensible, c_sublim, &
    c_evap, c_fcondtopn, c_topmeltn, c_lsrain, c_lssnow, c_cvrain, c_cvsnow,  &
    c_riverout, c_calving, c_mslp, c_tstar_sicen, c_surf_CO2, c_dust_dep

USE carbon_options_mod, ONLY: l_co2_interactive

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


INTEGER :: i              ! Loop counter
INTEGER :: j              ! Loop counter
INTEGER :: k              ! Ice category counter
INTEGER :: oasis_error, ft
INTEGER :: icode          ! Error return code (=0 is OK)
INTEGER :: info           ! Return code from MPP
INTEGER :: gather_pe      ! Processor for gathering
CHARACTER(LEN=errormessagelength) :: cmessage ! OUT - Error return message

INTEGER :: pointer_top_d  ! Pointers to D1 for stash data
INTEGER :: pointer_fcond_d  ! multi-cat top/fcondtop respectively
INTEGER :: pointer_sublim_d ! multi-cat sublimation
INTEGER :: pointer_tstar_sicen_d ! sea-ice surface temp on cats

INTEGER :: index_j

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='OASIS_UPDATECPL'

! Our row lengths and column lengths must be the local domain sizes
! for each PE.

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

IF (l_taux) THEN
  ! Set up fields from the U grid
  DO j=1,oasis_jmt_u
    DO i=1,oasis_imt
      c_taux(i,j) =                                               &
       d1(ja_taux+i-1+((j-1)*oasis_imt))
    END DO
  END DO
END IF

IF (l_tauy) THEN
  ! Set up fields from the V grid
  DO j=1,oasis_jmt_v
    DO i=1,oasis_imt
      c_tauy(i,j) =                                               &
       d1(ja_tauy+i-1+((j-1)*oasis_imt))
    END DO
  END DO
END IF

IF (l_bluerad) THEN
  ! Set up blue radiation field.
  DO j=1,oasis_jmt
    DO i=1,oasis_imt
      c_blue(i,j) =                                               &
       d1(ja_blue+i-1+((j-1)*oasis_imt))
    END DO
  END DO
END IF

IF (l_heatflux) THEN
  ! Set up the varias fields specific to heatflux
  DO j=1,oasis_jmt
    DO i=1,oasis_imt
      c_solar(i,j) =                                              &
       d1(ja_solar+i-1+((j-1)*oasis_imt))

      c_longwave(i,j) =                                           &
       d1(ja_longwave+i-1+((j-1)*oasis_imt))

      c_sensible(i,j) =                                           &
       d1(ja_sensible+i-1+((j-1)*oasis_imt))
    END DO
  END DO
END IF

IF (l_evap2d) THEN
  DO j=1,oasis_jmt
    DO i=1,oasis_imt
      ! PME components
      c_evap(i,j) =                                               &
        d1(ja_evap+i-1+((j-1)*oasis_imt))
    END DO
  END DO
END IF

IF (l_train) THEN
  DO j=1,oasis_jmt
    DO i=1,oasis_imt
      c_lsrain(i,j) =                                             &
       d1(ja_lsrain+i-1+((j-1)*oasis_imt))
      IF (l_param_conv) c_cvrain(i,j) =                           &
       d1(ja_cvrain+i-1+((j-1)*oasis_imt))
    END DO
  END DO
END IF

IF (l_tsnow) THEN
  DO j=1,oasis_jmt
    DO i=1,oasis_imt
      c_lssnow(i,j) =                                             &
        d1(ja_lssnow+i-1+((j-1)*oasis_imt))
      IF (l_param_conv) c_cvsnow(i,j) =                           &
        d1(ja_cvsnow+i-1+((j-1)*oasis_imt))
    END DO
  END DO
END IF


IF (l_runoff) THEN
  DO j=1,oasis_jmt
    DO i=1,oasis_imt
      ! River runoff
      c_riverout(i,j) =                                           &
        d1(ja_riverout+i-1+((j-1)*oasis_imt))
    END DO
  END DO
END IF

IF (l_w10) THEN
  DO j=1,oasis_jmt
    DO i=1,oasis_imt
      ! 10m wind
      c_w10(i,j) =                                            &
        d1(ja_w10+i-1+((j-1)*oasis_imt))
    END DO
  END DO
END IF

IF (l_sublim(1)) THEN
  ! Sublimation (either single or multi-category)
  pointer_sublim_d=ja_sublim

  DO k=1,nice_use
    DO j=1,oasis_jmt
      DO i=1,oasis_imt
        c_sublim(i,j,k)=d1(pointer_sublim_d)
        pointer_sublim_d=pointer_sublim_d+1
      END DO
    END DO
  END DO

  IF (.NOT. l_ctile) THEN
    ! In the non-coastal tiling case we convert sublimation
    ! from a total to a rate. Note that if coastal tiling is
    ! not used (l_ctile false), field 3231 (sublimation total) is
    ! set by oasis_inita2o to sublim and needs to be converted to
    ! a rate here. If coastal tiling is used (l_ctile true), field
    ! 3353 (sublimation rate) is set by oasis_inita2o so no extra
    ! processing is required.
    ! (See UM7.0 or earlier swap_a2o for historical perspective
    ! on the way these fields are handled).
    DO k=1,nice_use
      DO j=1,oasis_jmt
        DO i=1,oasis_imt
          c_sublim(i,j,k)=c_sublim(i,j,k)/                          &
           REAL((oasis_couple_freq_ao(1)*3600)+oasis_couple_freq_ao(2)*60)
        END DO
      END DO
    END DO
  END IF
END IF

IF (l_topmeltn(1)) THEN
  ! Topmelt(multi-category)
  pointer_top_d=ja_topmeltn

  DO k=1,nice
    DO j=1,oasis_jmt
      DO i=1,oasis_imt
        c_topmeltn(i,j,k)=d1(pointer_top_d)
        pointer_top_d=pointer_top_d+1
      END DO
    END DO
  END DO
END IF

IF (l_fcondtopn(1)) THEN
  ! Fcondtop (multi-category)
  pointer_fcond_d=ja_fcondtopn
  DO k=1,nice
    DO j=1,oasis_jmt
      DO i=1,oasis_imt
        c_fcondtopn(i,j,k)=d1(pointer_fcond_d)
        pointer_fcond_d=pointer_fcond_d+1
      END DO
    END DO
  END DO
END IF

IF (l_tstar_sicen(1)) THEN
! Sea ice surface skin temp on categories (K)
  pointer_tstar_sicen_d=ja_tstar_sicen
  DO k=1,nice
    DO j=1,oasis_jmt
      DO i=1,oasis_imt
        c_tstar_sicen(i,j,k)=d1(pointer_tstar_sicen_d)
        pointer_tstar_sicen_d=pointer_tstar_sicen_d+1
      END DO
    END DO
  END DO
END IF

IF (l_mslp) THEN
  DO j=1,oasis_jmt
    DO i=1,oasis_imt
      ! Mean Sea Level Pressure
      c_mslp(i,j) =                                           &
        d1(ja_mslp+i-1+((j-1)*oasis_imt))
    END DO
  END DO
END IF

IF ( (l_pCO2) .AND. (l_co2_interactive) ) THEN
  ! MEDUSA -> UKCA coupling: pC02
  ! pCO2 will be calculated from the surface CO2 in OASIS3_PUT

  DO j=1,oasis_jmt
    index_j=ja_CO2+offx-1+(j-1+offy)*(oasis_imt+2*offx)
    DO i=1,oasis_imt
      c_surf_CO2(i,j)=d1(index_j+i)
    END DO
  END DO

END IF

IF (l_dust_dep) THEN
  ! MEDUSA -> UKCA coupling: total dust deposition
  DO j=1,oasis_jmt
    DO i=1,oasis_imt
      c_dust_dep(i,j)=d1(ja_dust_dep+i-1+((j-1)*oasis_imt))
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN


END SUBROUTINE oasis_updatecpl
#endif
