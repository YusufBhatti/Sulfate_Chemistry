#if defined(MCT)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis3_puta2o()

USE atm_fields_bounds_mod

USE water_constants_mod, ONLY: lc, lf
USE oasis_atm_data_mod
USE jules_sea_seaice_mod, ONLY: nice, nice_use
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE Control_Max_Sizes
USE Decomp_DB
USE lbc_mod
USE um_types
USE coupling_control_mod,  ONLY: l_oasis_icecalve
USE oasis_timers, ONLY: l_oasis_timers, putstarts, putends
USE nlsizes_namelist_mod, ONLY: ntiles

USE cv_run_mod, ONLY: l_param_conv
USE conversions_mod, ONLY: zerodegc   

USE atm_fields_mod, ONLY: &
    c_solar, c_blue, c_longwave, c_taux, c_tauy, c_w10, c_sensible, c_sublim, &
    c_evap, c_fcondtopn, c_topmeltn, c_lsrain, c_lssnow, c_cvrain, c_cvsnow,  &
    c_riverout, c_calving, c_mslp, c_tstar_sicen, c_surf_CO2, c_dust_dep,     &
    frac_land, snodep_tile

USE rad_input_mod, ONLY: co2_mmr
USE carbon_options_mod, ONLY: l_co2_interactive
USE ccarbon, ONLY: m_air, m_CO2
USE land_soil_dimensions_mod, ONLY: land_points, land_ice_points,       &
     land_index, land_ice_index

IMPLICIT NONE


!
! Description:
! Put data from atmosphere to be got by NEMO (or other Ocean model) 
! and/or CICE (or other sea-ice model.)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
!
!==================================================================



!     Subroutine arguments
LOGICAL (KIND=logical32) :: l_second_order ! Indicates if 2nd order
                          ! conservative is employed with OASIS3-MCT
                          ! in which case we need to calculate
                          ! and pass field gradients to the put call.

!     Local variables
! Some of the following constants really should
! come from somewhere central rather than being defined here
REAL :: latentHeatOfCond=lc
REAL :: latentHeatOfFus=lf
REAL :: greenland_volume
REAL :: antarctica_volume

INTEGER (KIND=integer64) :: i              ! Loop counter
INTEGER (KIND=integer64) :: j              ! Loop counter
INTEGER (KIND=integer64) :: k              ! Ice category counter
INTEGER (KIND=integer64) :: tc              ! Loop counter
INTEGER (KIND=integer32) :: oasis_info     ! OASIS 32 bit return code
INTEGER (KIND=integer64) :: oasis_jmt_this ! OASIS jmt size by grid cell type
INTEGER (KIND=integer64) :: ft             ! Field type
INTEGER (KIND=integer64) :: s_cols, s_rows ! Send buffer sizes

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='OASIS3_PUTA2O'

! Constant used for converting surface CO2 (in MMR) to ocean
! partial pressure of CO2 (in micro-atmospheres)
REAL, PARAMETER :: CO2conv_a2o = m_air * 1.0E6 / m_CO2

! Perform necessary processing as described, for example,
! in 'Operations 1' on  HadGEM3 coupling diagram. (See HadGEM3
! coupling documentation)

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (l_oasis_timers) CALL putstarts

IF (put_step_ao) THEN
  ! We set up coupling data based on prognostic contents,
  ! Prepare fields for putting to the coupler
  ! by copying to our temporary arrays.

  IF (l_bluerad) THEN
    DO j=1,oasis_jmt
      DO i=1,oasis_imt
        ! Blue band radiation
        blue2d(i,j)=c_blue(i,j)
      END DO
    END DO
  END IF

  IF (l_evap2d) THEN
    DO j=1,oasis_jmt
      DO i=1,oasis_imt
        ! Evap
        evap2d(i,j)=c_evap(i,j)
      END DO
    END DO
  END IF

  IF (l_heatflux) THEN
    ! Note here the dependence of heatflux on evap2d and solar2d
    DO j=1,oasis_jmt
      DO i=1,oasis_imt
        ! Copy various heat flux components
        solar2d(i,j)=c_solar(i,j)
        longwave2d(i,j)=c_longwave(i,j)
        sensible2d(i,j)=c_sensible(i,j)

        IF (fland_loc(i,j) == 1.0) THEN
          heatflux(i,j)=0.0
        ELSE
          heatflux(i,j)=longwave2d(i,j)                       &
              -(sensible2d(i,j)+latentHeatOfCond*evap2d(i,j))
          IF (l_bluerad) THEN
            heatflux(i,j)=solar2d(i,j)+heatflux(i,j)
          END IF
        END IF
      END DO
    END DO
  END IF

  IF (l_w10) THEN
    DO j=1,oasis_jmt
      DO i=1,oasis_imt
        ! 10m wind
        w10(i,j) = c_w10(i,j)
      END DO
    END DO
  END IF

  IF (l_train) THEN
    DO j=1,oasis_jmt
      DO i=1,oasis_imt
        rainls(i,j) = c_lsrain(i,j)
        IF (l_param_conv) THEN
          rainconv(i,j) = c_cvrain(i,j)
        ELSE
          rainconv(i,j) = 0.0
        END IF
        IF (fland_loc(i,j) == 1.0) THEN
          totalrain(i,j)=0.0
        ELSE
          totalrain(i,j)=rainls(i,j)+rainconv(i,j)
        END IF
      END DO
    END DO
  END IF

  IF (l_tsnow) THEN
    DO j=1,oasis_jmt
      DO i=1,oasis_imt
        snowls(i,j) = c_lssnow(i,j)
        IF (l_param_conv) THEN
          snowconv(i,j) = c_cvsnow(i,j)
        ELSE
          snowconv(i,j) = 0.0
        END IF
        IF (fland_loc(i,j) == 1.0) THEN
          totalsnow(i,j)=0.0
        ELSE
          totalsnow(i,j)=snowls(i,j)+snowconv(i,j)
        END IF
      END DO
    END DO
  END IF

  IF (l_runoff) THEN
    IF (l_oasis_icecalve) THEN
      DO j=1,oasis_jmt
        DO i=1,oasis_imt
          ! River runoff
          riverout(i,j) = c_riverout(i,j)
          calving(i,j)= c_calving(i,j)
          ! Riverout needs scaling to take take account of
          ! coastal tiling points
          ! Iceberg calving field also needs scaling
          IF (fland_loc(i,j) /= 1.0) THEN
            riverout(i,j)=riverout(i,j)/(1.0-fland_loc(i,j))
            calving(i,j)=calving(i,j)/(1.0-fland_loc(i,j))
            ! Add calving into runoff for the moment
            riverout(i,j)=riverout(i,j)+calving(i,j)
            ! Also add latent heat flux term
            heatflux(i,j)=heatflux(i,j)-latentHeatOfFus*calving(i,j)
          END IF
        END DO
      END DO
    ELSE
      DO j=1,oasis_jmt
        DO i=1,oasis_imt
          ! River runoff
          riverout(i,j) = c_riverout(i,j)
          ! Riverout needs scaling to take take account of
          ! coastal tiling points
          IF (fland_loc(i,j) /= 1.0) THEN
            riverout(i,j)=riverout(i,j)/(1.0-fland_loc(i,j))
          END IF
        END DO
      END DO
    END IF
  END IF
  
  
  
  ! ----------------------Sea ice flux fields
  
  ! In the 'time-travelling ice' coupling framework, these fluxes are passed
  ! as local values, produced by dividing through 'future' ice concentration
  ! which has just been read in from oasis3_geto2a.  This results in a more
  ! physically realistic apportionment of energy between grid cells, and
  ! because the ice concentration is a future value, energy is conserved.

  IF (l_sublim(1)) THEN
    ! Sublimation (single or multi-category)
    DO k=1,nice_use
      DO j=1,oasis_jmt
        DO i=1,oasis_imt
          sublim(i,j,k) = c_sublim(i,j,k)
          IF (ice_fract_cat_future(i,j,k) /= 0.0) THEN
            sublim(i,j,k) = sublim(i,j,k) / ice_fract_cat_future(i,j,k)
          END IF
        END DO
      END DO
    END DO
  END IF

  IF (l_topmeltn(1)) THEN
    ! Topmelt (category)
    DO k=1,nice
      DO j=1,oasis_jmt
        DO i=1,oasis_imt
          topmeltn(i,j,k) = c_topmeltn(i,j,k)
          IF (ice_fract_cat_future(i,j,k) /= 0.0) THEN
            topmeltn(i,j,k) = topmeltn(i,j,k) / ice_fract_cat_future(i,j,k)
          END IF
        END DO
      END DO
    END DO
  END IF

  IF (l_fcondtopn(1)) THEN
    ! Fcondtop (category)
    DO k=1,nice
      DO j=1,oasis_jmt
        DO i=1,oasis_imt
          fcondtopn(i,j,k) = c_fcondtopn(i,j,k)
          IF (ice_fract_cat_future(i,j,k) /= 0.0) THEN
            fcondtopn(i,j,k) = fcondtopn(i,j,k) / ice_fract_cat_future(i,j,k)
          END IF
        END DO
      END DO
    END DO
  END IF

  IF (l_tstar_sicen(1)) THEN
    ! Sea ice surface skin temp on categories (degC)
    DO k=1,nice
      DO j=1,oasis_jmt
        DO i=1,oasis_imt
          tstar_sicen(i,j,k) = c_tstar_sicen(i,j,k) - zerodegc
        END DO
      END DO
    END DO
  END IF

  IF (l_taux) THEN
    ! Set up fields from the U grid
    DO j=1,oasis_jmt_u
      DO i=1,oasis_imt
        taux(i,j)=c_taux(i,j)
      END DO
    END DO
  END IF

  IF (l_tauy) THEN
    ! Set up fields from the V grid
    DO j=1,oasis_jmt_v
      DO i=1,oasis_imt
        tauy(i,j)=c_tauy(i,j)
      END DO
    END DO
  END IF

  IF (l_mslp) THEN
    ! Set up mslp field
    DO j=1,oasis_jmt
      DO i=1,oasis_imt
        mslp(i,j)=c_mslp(i,j)
      END DO
    END DO
  END IF
  
  IF (l_icesheet) THEN
    ! Pass ice sheet mass over Greenland and Antarctica through the
    ! coupler as two single floats copied onto two 2D fields.
  
    ! DEPENDS ON: ice_sheet_mass
    CALL ice_sheet_mass( land_ice_index,land_ice_points, ntiles,            &
    land_points, land_index, frac_land, snodep_tile,                        &
    greenland_volume, antarctica_volume                )
  
    ! Set up icesheet fields
    DO j=1,oasis_jmt
       DO i=1,oasis_imt
         icenorth(i,j)=greenland_volume
         icesouth(i,j)=antarctica_volume
       END DO
    END DO
  END IF
  
  IF (l_pCO2) THEN
    ! Set up the partial pressure of CO2 for coupling
    ! Determine if we want interactive CO2
    IF (l_co2_interactive) THEN

      ! pCO2 can be determined from the surface CO2
      DO j=1,oasis_jmt
        DO i=1,oasis_imt
          pCO2_out(i,j)=CO2conv_a2o * c_surf_CO2(i,j)
        END DO
      END DO

    ELSE

      ! pCO2 is a constant everywhere, because CO2 is
      pCO2_out(:,:)=CO2conv_a2o * co2_mmr

    END IF
  
  END IF

  IF (l_dust_dep) THEN
    ! Set up the total dust deposition for coupling
    DO j=1,oasis_jmt
      DO i=1,oasis_imt
        dust_dep_out(i,j)=c_dust_dep(i,j)
      END DO
    END DO
  END IF
  
  ! In principle we should be able to PUT fields in any order, without 
  ! needing to adhere to a strict sequence known to both the sending
  ! and receiving component. However, our ability to do that may be dependent
  ! on platform/MPI implementation (essentially size of available MPI buffers)
  ! and in some cases (as seen on the IBM Power7) we're forced to observe
  ! a strict call sequence (which runs contrary to the FLUME principle that 
  ! source and target should not have to presume anything about each other).
  !===============================================================

  DO tc = 1, transient_a2o_count

    IF (transient_a2o(tc)%indx > 0) THEN

      l_second_order=.FALSE.

      ! This is an expected outgoing transient, what grid type is it on
      IF (transient_a2o(tc)%grid == "V") THEN
        ft = fld_type_v
        oasis_jmt_this = oasis_jmt_v
      ELSE
        IF (transient_a2o(tc)%grid == "U") THEN
          ft = fld_type_u
          oasis_jmt_this = oasis_jmt_u
        ELSE
          ft = fld_type_p
          oasis_jmt_this = oasis_jmt

          IF (transient_a2o(tc)%order == 2) THEN
            l_second_order=.TRUE.
          END IF

        END IF
      END IF

      s_cols = lasize(1,ft,halo_type_no_halo)
      s_rows = lasize(2,ft,halo_type_no_halo)

      ! DEPENDS ON: oasis3_put
      CALL oasis3_put(transient_a2o(tc)%field                            &
                      ,transient_a2o(tc)%oasis_id                        &
                      ,oasis_imt,oasis_jmt_this                          &
                      ,transient_a2o(tc)%oasis_info                      &
                      ,ft                                                &
                      ,s_cols,s_rows,put_step_ao,l_second_order          &
                      )

      ! We don't actually check the return code here because a
      ! successful put does not imply a successful get. So in general
      ! it seems more efficient to leave any error handling to the
      ! recieving component rather than potentially perform 
      ! rafts of semi-redundant error code checks here.  

    END IF

  END DO  ! For each potential transient

END IF  ! put_step_ao=true

IF (l_oasis_timers) CALL putends
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE oasis3_puta2o
#endif
