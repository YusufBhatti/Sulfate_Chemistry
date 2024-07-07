! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis_point_translist(nice,nice_use)
!
! Set up pointers to transients involved in coupling.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
!
!===========================================================

USE oasis_atm_data_mod
USE coupling_control_mod, ONLY: l_oasis_obgc, l_oasis_icecalve
USE jules_sea_seaice_mod, ONLY: l_saldep_freeze, l_sice_multilayers
USE ukca_option_mod,      ONLY: l_ukca_prim_moc
USE jules_radiation_mod,  ONLY: l_sea_alb_var_chl
USE ancil_mod,            ONLY: num_ancil_requests, ancil_requests
USE um_stashcode_mod,     ONLY: stashcode_ocnsrf_chlorophyll,  &
                                stashcode_dms_conc_sea
USE ereport_mod, ONLY: umMessage, ereport
USE umPrintMgr, ONLY: newline

IMPLICIT NONE

INTEGER, INTENT(IN) :: nice     ! Dims for multi-cat/layer ice fields
INTEGER, INTENT(IN) :: nice_use ! Dims for multi-cat/layer ice fields

! Local variables
INTEGER ::    icode             ! Error return code
INTEGER ::    iunmatched        ! Counter for unmatched coupling fields 

! Local loop indices
INTEGER :: tc
INTEGER :: i, n

CHARACTER(LEN=*), PARAMETER :: RoutineName='OASIS_POINT_TRANSLIST'

!---------------------------------------------------
! ALLOCATE space for incoming arrays and associated
! fields required in processing them. 
!---------------------------------------------------

!--------------------------------------------------
! Many fields may be allocated conditionally depending on which
! coupling fields are active. 
!
! Some fields, e.g. solar radiation, may be required under 
! a range of conditions, hence we need conditional allocation
! statements for those. Similarly multi-category fields are 
! allocated on detection of the presence of any category element
! of a particular field with a test to ensure that no attempt is
! made to repeatedly allocate them.  
! Other fields have a 1:1 relationship with the transient 
! definition and are simply allocated based on whether that coupling 
! field code is detected or not. 
!--------------------------------------------------

IF (l_oasis_icecalve) THEN 
  ALLOCATE(calving(oasis_imt,oasis_jmt))
  calving(:,:)=0.0
END IF
IF (l_sice_multilayers) THEN
  ALLOCATE(ocn_icet_gb(oasis_imt,oasis_jmt))
  ocn_icet_gb(:,:)=0.0
END IF
IF (l_saldep_freeze) THEN
  ALLOCATE(ocn_sstfrz(oasis_imt,oasis_jmt)) 
  ocn_sstfrz(:,:)=0.0
END IF

!-----------------------------------------------------------
! Now set up pointers to make processing of incoming
! and outgoing arrays easier later on. Array field_id must
! strictly match namcouple field ID number and be unique.
!-----------------------------------------------------------
iunmatched = 0 

DO tc = 1, transient_a2o_count

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Check potential outgoing atmos -> ocean fields
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (transient_a2o(tc)%field_id == vind_heatflux) THEN
    l_heatflux = .TRUE.

    IF (.NOT.ALLOCATED(solar2d)) THEN
      ALLOCATE(solar2d(oasis_imt,oasis_jmt))
      solar2d(:,:)=0.0
    END IF

    IF (.NOT.ALLOCATED(evap2d)) THEN
      ALLOCATE(evap2d(oasis_imt,oasis_jmt))
      evap2d(:,:)=0.0
    END IF

    ALLOCATE(longwave2d(oasis_imt,oasis_jmt))
    longwave2d(:,:)=0.0 

    ALLOCATE(sensible2d(oasis_imt,oasis_jmt))
    sensible2d(:,:)=0.0
 
    ALLOCATE(heatflux(oasis_imt,oasis_jmt))
    heatflux(:,:)=0.0

    transient_a2o(tc)%field => heatflux(:,:)
  END IF

  IF (transient_a2o(tc)%field_id == vind_bluerad) THEN
    l_bluerad = .TRUE.
    ALLOCATE(blue2d(oasis_imt,oasis_jmt))
    blue2d(:,:)=0.0
    transient_a2o(tc)%field => blue2d(:,:)
  END IF

  IF (transient_a2o(tc)%field_id == vind_solar) THEN
    l_solar = .TRUE.
    IF (.NOT.ALLOCATED(solar2d)) THEN
      ALLOCATE(solar2d(oasis_imt,oasis_jmt))
      solar2d(:,:)=0.0
    END IF
    transient_a2o(tc)%field => solar2d(:,:)
  END IF

  IF (transient_a2o(tc)%field_id == vind_runoff) THEN
    l_runoff = .TRUE.
    ALLOCATE(riverout(oasis_imt,oasis_jmt))
    riverout(:,:)=0.0
    transient_a2o(tc)%field => riverout(:,:)
  END IF

  IF (transient_a2o(tc)%field_id == vind_w10) THEN
    l_w10 = .TRUE.
    ALLOCATE(w10(oasis_imt,oasis_jmt))
    w10(:,:)=0.0
    transient_a2o(tc)%field => w10(:,:)
  END IF

  IF (transient_a2o(tc)%field_id == vind_train) THEN       
    l_train = .TRUE.
    ALLOCATE(rainls(oasis_imt,oasis_jmt))
    rainls(:,:)=0.0
    ALLOCATE(rainconv(oasis_imt,oasis_jmt))
    rainconv(:,:)=0.0
    ALLOCATE(totalrain(oasis_imt,oasis_jmt))
    totalrain(:,:)=0.0
    transient_a2o(tc)%field => totalrain(:,:)
  END IF

  IF (transient_a2o(tc)%field_id == vind_tsnow) THEN       
    l_tsnow = .TRUE.
    ALLOCATE(snowls(oasis_imt,oasis_jmt))
    snowls(:,:)=0.0
    ALLOCATE(snowconv(oasis_imt,oasis_jmt))
    snowconv(:,:)=0.0
    ALLOCATE(totalsnow(oasis_imt,oasis_jmt))
    totalsnow(:,:)=0.0
    transient_a2o(tc)%field => totalsnow(:,:)
  END IF

  IF (transient_a2o(tc)%field_id == vind_evap2d) THEN       
    l_evap2d = .TRUE.
    IF (.NOT.ALLOCATED(evap2d)) THEN
      ALLOCATE(evap2d(oasis_imt,oasis_jmt))
      evap2d(:,:)=0.0
    END IF
    transient_a2o(tc)%field => evap2d(:,:)
  END IF

  DO n = 1 , nice
    IF (transient_a2o(tc)%field_id == vind_topmeltn(n)) THEN
      l_topmeltn(n) = .TRUE.
      IF (.NOT.ALLOCATED(topmeltn)) THEN
        ALLOCATE(topmeltn(oasis_imt,oasis_jmt,nice))
        topmeltn(:,:,:)=0.0
      END IF
      transient_a2o(tc)%field => topmeltn(:,:,n)
    END IF

    IF (transient_a2o(tc)%field_id == vind_fcondtopn(n)) THEN
      l_fcondtopn(n) = .TRUE.
      IF (.NOT.ALLOCATED(fcondtopn)) THEN
        ALLOCATE(fcondtopn(oasis_imt,oasis_jmt,nice))
        fcondtopn(:,:,:)=0.0
      END IF
      transient_a2o(tc)%field => fcondtopn(:,:,n)
    END IF

    IF (transient_a2o(tc)%field_id == vind_tstar_sicen(n)) THEN
      l_tstar_sicen(n) = .TRUE.
      IF (.NOT.ALLOCATED(tstar_sicen)) THEN
        ALLOCATE(tstar_sicen(oasis_imt,oasis_jmt,nice))
        tstar_sicen(:,:,:)=0.0
      END IF
      transient_a2o(tc)%field => tstar_sicen(:,:,n)
    END IF

  END DO  ! Over n

  DO n = 1 , nice_use
    IF (transient_a2o(tc)%field_id == vind_sublim(n)) THEN
      l_sublim(n) = .TRUE.
      IF (.NOT.ALLOCATED(sublim)) THEN
        ALLOCATE(sublim(oasis_imt,oasis_jmt,nice_use))
        sublim(:,:,:)=0.0
      END IF
      transient_a2o(tc)%field => sublim(:,:,n)
    END IF
  END DO  ! Over n

  IF (transient_a2o(tc)%field_id == vind_taux) THEN
    l_taux = .TRUE.
    ALLOCATE(taux(oasis_imt,oasis_jmt_u))
    taux(:,:)=0.0
    transient_a2o(tc)%field => taux(:,:)
  END IF

  IF (transient_a2o(tc)%field_id == vind_tauy) THEN
    l_tauy = .TRUE.
    ALLOCATE(tauy(oasis_imt,oasis_jmt_v))
    tauy(:,:)=0.0
    transient_a2o(tc)%field => tauy(:,:)
  END IF

  IF (transient_a2o(tc)%field_id == vind_mslp) THEN
    l_mslp = .TRUE.
    ALLOCATE(mslp(oasis_imt,oasis_jmt))
    mslp(:,:)=0.0 
    transient_a2o(tc)%field => mslp(:,:)
  END IF
 
  IF (transient_a2o(tc)%field_id == vind_icesheet_n) THEN
    l_icesheet = .TRUE.
    ALLOCATE(icenorth(oasis_imt,oasis_jmt)) 
    icenorth(:,:)=0.0
    transient_a2o(tc)%field => icenorth(:,:)
  END IF

  IF (transient_a2o(tc)%field_id == vind_icesheet_s) THEN
    l_icesheet = .TRUE.
    ALLOCATE(icesouth(oasis_imt,oasis_jmt)) 
    icesouth(:,:)=0.0
    transient_a2o(tc)%field => icesouth(:,:)
  END IF
  
  IF (l_oasis_obgc) THEN
    ! Enables coupling to ocean biogeochemistry

    IF (transient_a2o(tc)%field_id == vind_pCO2) THEN
      l_pCO2 = .TRUE.
      ALLOCATE(pCO2_out(oasis_imt,oasis_jmt))
      pCO2_out(:,:)=0.0
      transient_a2o(tc)%field => pCO2_out(:,:)
    END IF

    IF (transient_a2o(tc)%field_id == vind_dust_dep) THEN
      l_dust_dep = .TRUE.
      ALLOCATE(dust_dep_out(oasis_imt,oasis_jmt))
      dust_dep_out(:,:)=0.0
      transient_a2o(tc)%field => dust_dep_out(:,:)
    END IF

  END IF

  ! Check to see that this coupling field has been 
  ! associated with a source UM field
  IF (.NOT.(ASSOCIATED(transient_a2o(tc)%field))) THEN
    iunmatched = iunmatched + 1
    icode = -1   ! Non terminal error code.
    WRITE(umMessage,'(A,I0)')                             &
      "Unmatched A=>O coupling field ID: ", transient_a2o(tc)%field_id
    CALL EReport(RoutineName, icode, umMessage)
  END IF


END DO ! A2O transient loop


DO tc = 1, transient_o2a_count

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Check potential incoming ocean -> atmos fields
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (transient_o2a(tc)%field_id == vind_ocn_sst) THEN
    ! Ocean SST will be required in all "normal" atmos-ocean
    ! coupled models
    l_ocn_sst = .TRUE.
    ALLOCATE(ocn_sst(oasis_imt,oasis_jmt))
    ocn_sst(:,:)=0.0
    ALLOCATE(ocn_sst_orig(oasis_imt,oasis_jmt))
    ocn_sst_orig(:,:)=0.0

    ! Aggregate ice values are always needed even if no ice coupling
    ! is employed. 
    ALLOCATE(ocn_hice(oasis_imt,oasis_jmt))
    ocn_hice(:,:)=0.0
    ALLOCATE(ocn_snowthick(oasis_imt,oasis_jmt))
    ocn_snowthick(:,:)=0.0
    ALLOCATE(ocn_freeze(oasis_imt,oasis_jmt))
    ocn_freeze(:,:)=0.0

    ! Tstar fields required in all cases with SST active. 
    ALLOCATE(tstar_local(oasis_imt,oasis_jmt))
    tstar_local(:,:)=0.0
    ALLOCATE(tstar_ssi(oasis_imt,oasis_jmt))
    tstar_ssi(:,:)=0.0
    ALLOCATE(tstar_sice_local(oasis_imt,oasis_jmt))
    tstar_sice_local(:,:)=0.0
    ALLOCATE(tstar_sice_cat_local(oasis_imt,oasis_jmt,nice_use))
    tstar_sice_cat_local(:,:,:)=0.0
    ALLOCATE(tstar_land_local(oasis_imt,oasis_jmt))
    tstar_land_local(:,:)=0.0

    transient_o2a(tc)%field => ocn_sst(:,:)
  END IF

  IF (transient_o2a(tc)%field_id == vind_ocn_sstfrz) THEN
    ! ocn_sstfrz may only be used as a coupling field if l_saldep_freeze
    ! is true, in which case the target array has been allocated above. 
    l_ocn_sstfrz = .TRUE.
    transient_o2a(tc)%field => ocn_sstfrz(:,:)
  END IF

  DO n = 1 , nice
    IF (transient_o2a(tc)%field_id == vind_ocn_freezen(n)) THEN
      l_ocn_freezen(n) = .TRUE.
      IF (.NOT.ALLOCATED(ocn_freezen)) THEN
        ALLOCATE(ocn_freezen(oasis_imt,oasis_jmt,nice))
        ocn_freezen(:,:,:)=0.0
      END IF

      ! ocn_freezen_first is always required in some capacity if 
      ! ocn_freezen is true but may not necessarily be a full 
      ! separate coupling field hence we need to ensure it is 
      ! allocated here too. 
      IF (.NOT.ALLOCATED(ocn_freezen_first)) THEN
        ALLOCATE(ocn_freezen_first(oasis_imt,oasis_jmt,nice))
        ocn_freezen_first(:,:,:)=0.0
      END IF

      ! This field is read from ocn_freezen_first, but is only required for
      ! dividing through fluxes in oasis_puta2o, so is only defined in
      ! oasis_atm_data_mod.
      IF (.NOT.ALLOCATED(ice_fract_cat_future)) THEN
        ALLOCATE(ice_fract_cat_future(oasis_imt,oasis_jmt,nice))
        ice_fract_cat_future(:,:,:)=0.0
      END IF

      transient_o2a(tc)%field => ocn_freezen(:,:,n)
    END IF

    IF (transient_o2a(tc)%field_id == vind_ocn_freezen_first(n)) THEN
      l_ocn_freezen_first(n) = .TRUE.
      IF (.NOT.ALLOCATED(ocn_freezen_first)) THEN
        ALLOCATE(ocn_freezen_first(oasis_imt,oasis_jmt,nice))
        ocn_freezen_first(:,:,:)=0.0
      END IF
      transient_o2a(tc)%field => ocn_freezen_first(:,:,n)
    END IF

    IF (transient_o2a(tc)%field_id == vind_ocn_snowthickn(n)) THEN
      l_ocn_snowthickn(n) = .TRUE.
      IF (.NOT.ALLOCATED(ocn_snowthickn)) THEN
        ALLOCATE(ocn_snowthickn(oasis_imt,oasis_jmt,nice))
        ocn_snowthickn(:,:,:)=0.0
      END IF
      transient_o2a(tc)%field => ocn_snowthickn(:,:,n)
      transient_o2a(tc)%ice_min = .TRUE.
    END IF

    IF (transient_o2a(tc)%field_id == vind_ocn_hicen(n)) THEN
      l_ocn_hicen(n) = .TRUE.
      IF (.NOT.ALLOCATED(ocn_hicen)) THEN
        ALLOCATE(ocn_hicen(oasis_imt,oasis_jmt,nice))
        ocn_hicen(:,:,:)=0.0
      END IF
      transient_o2a(tc)%field => ocn_hicen(:,:,n)
      transient_o2a(tc)%ice_min = .TRUE.
    END IF

    IF (transient_o2a(tc)%field_id == vind_ocn_icetn(n)) THEN
      l_ocn_icetn(n) = .TRUE.
      IF (.NOT.ALLOCATED(ocn_icetn)) THEN
        ALLOCATE(ocn_icetn(oasis_imt,oasis_jmt,nice))
        ocn_icetn(:,:,:)=0.0
      END IF
      transient_o2a(tc)%field => ocn_icetn(:,:,n)
      transient_o2a(tc)%ice_min = .TRUE.
    END IF

    IF (transient_o2a(tc)%field_id == vind_ocn_icekn(n)) THEN
      l_ocn_icekn(n) = .TRUE.
      IF (.NOT.ALLOCATED(ocn_icekn)) THEN
        ALLOCATE(ocn_icekn(oasis_imt,oasis_jmt,nice))
        ocn_icekn(:,:,:)=0.0
      END IF
      transient_o2a(tc)%field => ocn_icekn(:,:,n)
      transient_o2a(tc)%ice_min = .TRUE.
    END IF

    IF (transient_o2a(tc)%field_id == vind_pond_frac_n(n)) THEN
      l_pond_frac_n = .TRUE.
      IF (.NOT.ALLOCATED(pond_frac_n)) THEN
        ALLOCATE(pond_frac_n(oasis_imt,oasis_jmt,nice))
        pond_frac_n(:,:,:)=0.0
      END IF
      transient_o2a(tc)%field => pond_frac_n(:,:,n)
    END IF

    IF (transient_o2a(tc)%field_id == vind_pond_depth_n(n)) THEN
      l_pond_depth_n = .TRUE.
      IF (.NOT.ALLOCATED(pond_depth_n)) THEN
        ALLOCATE(pond_depth_n(oasis_imt,oasis_jmt,nice))
        pond_depth_n(:,:,:)=0.0
      END IF
      transient_o2a(tc)%field => pond_depth_n(:,:,n)
    END IF

  END DO 

  IF (transient_o2a(tc)%field_id == vind_ocn_u) THEN    
    ! Velocity components will appear in all normal atmos-ocean coupled models
    l_ocn_u = .TRUE.
    ALLOCATE(ocn_u(oasis_imt,oasis_jmt_u))
    ocn_u(:,:)=0.0
    transient_o2a(tc)%field => ocn_u(:,:)
  END IF

  IF (transient_o2a(tc)%field_id == vind_ocn_v) THEN
    ! Velocity components will appear in all normal atmos-ocean coupled models
    l_ocn_v = .TRUE.
    ALLOCATE(ocn_v(oasis_imt,oasis_jmt_v))
    ocn_v(:,:)=0.0
    transient_o2a(tc)%field => ocn_v(:,:)
  END IF

  IF (transient_o2a(tc)%field_id == vind_ocn_freeze) THEN
    l_ocn_freeze = .TRUE.
    IF (.NOT.ALLOCATED(ocn_freeze)) THEN
      ALLOCATE(ocn_freeze(oasis_imt,oasis_jmt))
      ocn_freeze(:,:)=0.0
    END IF
    transient_o2a(tc)%field => ocn_freeze(:,:)
  END IF

  IF (l_oasis_obgc) THEN
    ! Enables coupling to ocean biogeochemistry

    IF (transient_o2a(tc)%field_id == vind_dms_conc) THEN
      ! Check DMS conc. is not updated from ancillaries
      DO i = 1, num_ancil_requests
        IF ((ancil_requests(i)%stashcode == stashcode_dms_conc_sea).AND.  &
            (ancil_requests(i)%interval > 0)) THEN
          icode=1
          WRITE(umMessage,'(A)')                             &
            "Both DMS conc. coupling and DMS conc."//newline//&
            "ancillary updating have been requested."
          CALL EReport(RoutineName, icode, umMessage)
        END IF   
      END DO
 
      l_dms_conc = .TRUE.
      ALLOCATE(dms_conc_in(oasis_imt,oasis_jmt))
      dms_conc_in(:,:)=0.0
      transient_o2a(tc)%field => dms_conc_in(:,:)
    END IF

    IF (transient_o2a(tc)%field_id == vind_CO2flux) THEN
      l_CO2flux = .TRUE.
      ALLOCATE(CO2flux_in(oasis_imt,oasis_jmt))
      CO2flux_in(:,:)=0.0
      transient_o2a(tc)%field => CO2flux_in(:,:)
    END IF

    IF (transient_o2a(tc)%field_id == vind_chloro_conc) THEN

      ! Chlorophyll coupling field is present.
      IF (.NOT.(l_ukca_prim_moc.OR.l_sea_alb_var_chl)) THEN
        ! If we don't have the relevent options active then 
        ! we can't couple to this field.  
        icode=1
        WRITE(umMessage,'(A)')                             &
          "Chlorophyll coupling field requested"//newline//&
          "But neither  l_ukca_prim_moc or l_sea_alb_var_chl are active."
        CALL EReport(RoutineName, icode, umMessage)
      END IF 

      ! Check chlorophyll is not updated from ancillaries
      DO i = 1, num_ancil_requests
        IF ((ancil_requests(i)%stashcode == stashcode_ocnsrf_chlorophyll).AND. &
            (ancil_requests(i)%interval > 0)) THEN
          icode=1
          WRITE(umMessage,'(A)')                             &
            "Both chlorophyll coupling and chlorophyll"//newline//&
            "ancillary updating have been requested."
          CALL EReport(RoutineName, icode, umMessage)
        END IF   
      END DO

      l_chloro_conc = .TRUE.
      ALLOCATE(chloro_conc_in(oasis_imt,oasis_jmt))
      chloro_conc_in(:,:)=0.0
      transient_o2a(tc)%field => chloro_conc_in(:,:)
    END IF

  END IF

  ! Check to see that this coupling field has been 
  ! associated with a target UM field
  IF (.NOT.(ASSOCIATED(transient_o2a(tc)%field))) THEN
    iunmatched = iunmatched + 1
    icode = -1   ! Non terminal error code.
    WRITE(umMessage,'(A,I0)')                             &
      "Unmatched O=>A coupling field ID: ", transient_o2a(tc)%field_id
    CALL EReport(RoutineName, icode, umMessage)
  END IF

END DO ! Over tc

IF (iunmatched > 0) THEN
  icode = iunmatched ! Terminal error code.
  WRITE(umMessage,'(I0,A)')                             &
    iunmatched, " coupling fields are unmatched in the UM. "//newline//&
    "See UM output for full details."
  CALL EReport(RoutineName, icode, umMessage)
END IF



END SUBROUTINE oasis_point_translist
