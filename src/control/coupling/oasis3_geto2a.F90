#if defined(MCT)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis3_geto2a()

USE atm_fields_bounds_mod

! OASIS3-MCT is the only supported coupler 
USE mod_prism, ONLY: PRISM_NotDef

USE water_constants_mod, ONLY: tfs,tm
USE oasis_atm_data_mod, ONLY: get_step_oa, l_ocn_sst, l_dms_conc, l_CO2flux,  &
    l_ocn_freezen, l_chloro_conc,                                             &
    ocn_sst, ocn_sst_orig, ocn_sstfrz, ocn_freeze, ocn_freezen, ocn_hice,     &
    ocn_hicen, ocn_snowthick, ocn_snowthickn, tstar_local, tstar_ssi,         &
    tstar_sice_local, tstar_sice_cat_local, tstar_land_local,                 &
    ice_fract_cat_future, ocn_freezen_first, pond_frac_n, pond_depth_n,       & 
    ocn_icetn, ocn_icekn, ocn_icet_gb, fland_loc, ocn_u,ocn_v,                &
    vind_ocn_sst, vind_ocn_freezen, vind_ocn_freezen_first, vind_ocn_hicen,   &
    vind_ocn_snowthickn, vind_ocn_u, vind_ocn_v, vind_dms_conc, vind_co2flux, &
    vind_chloro_conc,                                                         &
    oasis_imt, oasis_jmt, oasis_jmt_u, oasis_jmt_v,                           &
    transient_o2a, transient_o2a_count,                                       &
    dms_conc_in, CO2flux_in, chloro_conc_in,                                  &
    TI_min, TI_max, aicenmin, hi_min
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars 
USE ereport_mod, ONLY: ereport
USE Control_Max_Sizes
USE cderived_mod, ONLY: delta_lambda
USE Decomp_DB

USE lbc_mod
USE UM_types
USE c_kappai, ONLY: kappai,kappai_snow,rhosnow
USE gen_phys_inputs_mod, ONLY: l_leads_temp_prog

USE jules_sea_seaice_mod, ONLY: nice, nice_use, l_sice_multilayers, &
                                 l_ctile, l_sice_meltponds_cice,     &
                                 l_saldep_freeze
USE oasis_timers, ONLY: l_oasis_timers, getstarts, getends
USE oasis_operations_mod, ONLY: oasis_get_ok
USE atm_fields_bounds_mod, ONLY: udims, vdims

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY: global_row_length, row_length, &
    rows, n_rows

USE atm_fields_mod, ONLY: &
    tstar, tstar_land, tstar_sea, tstar_sice, tstar_sice_cat, &
    ti, ti_cat, ice_k_cat, ice_fract_cat, ice_thick_cat,      &
    snodep_sea, snodep_sea_cat,ice_fraction, ice_thickness,   &
    pond_frac_cat, pond_depth_cat,dms_conc, CO2flux, sstfrz,  &
    chloro_sea, u_sea, v_sea

USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE

!
! Description:
! Receive incoming coupling data from the ocean (typically NEMO) and/or 
! sea-ice (typically CICE) to be used to drive the atmosphere.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
!
!=======================================================================

!     Local variables
INTEGER (KIND=integer64) :: i, j, k, n, tc ! Loop counters
INTEGER (KIND=integer64) :: r_cols, r_rows ! Receive buffer dimensions

INTEGER (KIND=integer64) :: ft             ! Field type index
INTEGER (KIND=integer32) :: icode          ! UM error code

! OASIS3-MCT 32 bit return codes:
INTEGER (KIND=integer32) :: infosst, infou, infov
INTEGER (KIND=integer32) :: infofreezen, infofreezen_first, infosnowthickn
INTEGER (KIND=integer32) :: infohicen, infoicetn, infoicekn
INTEGER (KIND=integer32) :: info_dms_conc, info_CO2flux, info_chloro_conc

INTEGER (KIND=integer64) :: oasis_jmt_this   ! J dim dependent on field type

CHARACTER(LEN=errormessagelength) :: Cmessage
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='OASIS3_GETO2A'

! Initialise potential incoming fields
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (l_oasis_timers) CALL getstarts

! If this is a timestep when we are genuinely expecting to
! recieve incoming fields, then check each potential field.
IF (get_step_oa) THEN

  IF (l_dms_conc) THEN
    dms_conc_in(:,:) = 0.0
  END IF
  IF (l_CO2flux) THEN
    CO2flux_in(:,:) = 0.0
  END IF
  IF (l_chloro_conc) THEN
    chloro_conc_in(:,:) = 0.0
  END IF

  IF (l_ocn_sst) THEN
    ! Save UM surface temperature values in order to overlay 
    ! dummy values received over land from the coupler with true values. 
    IF (l_ctile) THEN
      DO j = 1, oasis_jmt
        DO i = 1, oasis_imt
          ocn_sst(i,j)=tstar_sea(i,j)
          ocn_sst_orig(i,j)=tstar_sea(i,j)
        END DO
      END DO
    ELSE
      DO j = 1, oasis_jmt
        DO i = 1, oasis_imt
          ocn_sst(i,j)=tstar(i,j)
          ocn_sst_orig(i,j)=tstar(i,j)
        END DO
      END DO
    END IF

    ! Aggregate ice values are always needed with SST even if no 
    ! ice coupling is employed. 
    ocn_freeze(:,:) = 0.0
    ocn_hice(:,:) = 0.0
    ocn_snowthick(:,:) = 0.0

    ! Similarly, tstar fields are also always needed with SST 
    tstar_local(:,:) = 0.0
    tstar_ssi(:,:) = 0.0

    DO j = 1, oasis_jmt
      DO i = 1, oasis_imt
        tstar_sice_local(i,j)=tstar_sice(i,j,1)
        tstar_land_local(i,j)=tstar_land(i,j)
      END DO
    END DO

    IF (nice_use  >   1) THEN
      DO k = 1, nice_use
        DO j = 1, oasis_jmt
          DO i = 1, oasis_imt
            tstar_sice_cat_local(i,j,k)=tstar_sice_cat(i,j,k)
          END DO
        END DO
      END DO
    ELSE    ! TSTAR_SICE_CAT doesn't exist
      DO j = 1, oasis_jmt
        DO i = 1, oasis_imt
          tstar_sice_cat_local(i,j,1)=tstar_sice(i,j,1)
        END DO
      END DO
    END IF
  END IF  ! if l_ocn_sst = TRUE

  IF (l_ocn_freezen(1)) THEN
    ! In all configurations involving seaice coupling, 
    ! freezen, snowthickn and hicen must all be present
    ! at the same time. Hence we can test for the presence 
    ! of any one of these fields when initialising all of them.

    ! Initalise future ice conc. (To be read from ocn_freezen_first).
    ice_fract_cat_future(:,:,:) = 0.0

    DO k=1,nice
      DO j = 1, oasis_jmt
        DO i = 1, oasis_imt
          ocn_freezen(i,j,k)=ice_fract_cat(i,j,k)
          ocn_freezen_first(i,j,k) = ice_fract_cat_future(i,j,k)
          ocn_hicen(i,j,k)=ice_thick_cat(i,j,k)
          ocn_snowthickn(i,j,k)=snodep_sea_cat(i,j,k)
        END DO
      END DO
    END DO
  END IF

  IF (l_sice_meltponds_cice) THEN
    DO k=1,nice
      DO j = 1, oasis_jmt
        DO i = 1, oasis_imt
          pond_frac_n(i,j,k)=pond_frac_cat(i,j,k)
          pond_depth_n(i,j,k)=pond_depth_cat(i,j,k)
        END DO
      END DO
    END DO
  END IF

  IF (l_saldep_freeze) THEN
    DO j = 1, oasis_jmt
      DO i = 1, oasis_imt
        ocn_sstfrz(i,j) = sstfrz(i,j)
      END DO
    END DO
  END IF

  IF (l_sice_multilayers) THEN
    ocn_icet_gb(:,:) = 0.0
    DO k=1,nice
      DO j = 1, oasis_jmt
        DO i = 1, oasis_imt
          ocn_icetn(i,j,k)=ti_cat(i,j,k)
          ocn_icekn(i,j,k)=ice_k_cat(i,j,k)
        END DO
      END DO
    END DO
  END IF

  ! Initialise coupler return codes, since some fields may not be involved 
  ! in coupling. We use PRISM_NotDef which is defined by OASIS3-MCT vn2+ 
  infosst = PRISM_NotDef
  infou = PRISM_NotDef
  infov = PRISM_NotDef
  infofreezen = PRISM_NotDef
  infosnowthickn = PRISM_NotDef
  infohicen  = PRISM_NotDef
  infoicetn = PRISM_NotDef
  infoicekn = PRISM_NotDef
  info_dms_conc = PRISM_NotDef
  info_CO2flux = PRISM_NotDef
  info_chloro_conc = PRISM_NotDef

  DO tc = 1, transient_o2a_count

    ! It is usually important to check fields in the 
    ! right sequence, at least on the platforms tested to date.  
    ! For a "few" fields, we can get away with a random order, but
    ! usually for the numbers of fields involved in realistic models
    ! the MPI buffers appear unable to cope.
    IF (transient_o2a(tc)%indx > 0) THEN

        ! This is an expexted incoming transient, what grid type is it on
      IF (transient_o2a(tc)%grid == "V") THEN
        ft = fld_type_v
        oasis_jmt_this = oasis_jmt_v
      ELSE
        IF (transient_o2a(tc)%grid == "U") THEN
          ft = fld_type_u
          oasis_jmt_this = oasis_jmt_u
        ELSE
          ft = fld_type_p
          oasis_jmt_this = oasis_jmt
        END IF
      END IF

      r_cols = oasis_imt        !lasize(1,ft,halo_type_no_halo)
      r_rows = oasis_jmt_this   !lasize(2,ft,halo_type_no_halo)

      ! DEPENDS ON: oasis3_get
      CALL oasis3_get(transient_o2a(tc)%field                           &
                     ,transient_o2a(tc)%oasis_id                        &
                     ,oasis_imt,oasis_jmt_this                          &
                     ,transient_o2a(tc)%oasis_info,ft                   &
                     ,r_cols,r_rows,get_step_oa                         &
                     )

      IF (oasis_get_ok(transient_o2a(tc)%oasis_info)) THEN

        ! Check if we need some standard processing on ice fields
        ! to impose meaningful minima.
        IF (transient_o2a(tc)%ice_min) THEN

          DO j = 1, oasis_jmt
            DO i = 1, oasis_imt
              transient_o2a(tc)%field(i,j) =      &
                         MAX(0.0,transient_o2a(tc)%field(i,j))
            END DO
          END DO

        END IF

        IF (transient_o2a(tc)%field_id == vind_ocn_sst) THEN 
          infosst = transient_o2a(tc)%oasis_info
        ELSE IF (transient_o2a(tc)%field_id == vind_ocn_freezen(1)) THEN 
          infofreezen = transient_o2a(tc)%oasis_info
        ELSE IF (transient_o2a(tc)%field_id == vind_ocn_freezen_first(1)) THEN 
          infofreezen_first = transient_o2a(tc)%oasis_info
        ELSE IF (transient_o2a(tc)%field_id == vind_ocn_hicen(1)) THEN 
          infohicen = transient_o2a(tc)%oasis_info
        ELSE IF (transient_o2a(tc)%field_id == vind_ocn_snowthickn(1)) THEN 
          infosnowthickn = transient_o2a(tc)%oasis_info
        ELSE IF (transient_o2a(tc)%field_id == vind_ocn_u) THEN 
          infou = transient_o2a(tc)%oasis_info
        ELSE IF (transient_o2a(tc)%field_id == vind_ocn_v) THEN 
          infov = transient_o2a(tc)%oasis_info
        ELSE IF (transient_o2a(tc)%field_id == vind_dms_conc) THEN 
          info_dms_conc = transient_o2a(tc)%oasis_info
        ELSE IF (transient_o2a(tc)%field_id == vind_CO2flux) THEN 
          info_CO2flux = transient_o2a(tc)%oasis_info
        ELSE IF (transient_o2a(tc)%field_id == vind_chloro_conc) THEN 
          info_chloro_conc = transient_o2a(tc)%oasis_info
        END IF
   

      END IF ! Info indicates data received

    END IF ! Transient in use

  END DO ! Over tc


  ! Now we have received all the incoming ocean and/or seaice coupling
  ! fields which we are going to receive on this TS.
  ! We may need to do some manipulation depending on what we've been sent.


  ! Perform necessary processing as described, for example,
  ! in 'Operations 5' on  HadGEM3 coupling diagram. (See HadGEM3
  ! coupling documentation)

  ! All processing depends on having received ocean SST data.
  ! If this hasn't happened, then there's nothing to be done.
  ! Hence, once we've established that we've received sst then we 
  ! begin checking other fields. 
  ! This has been coded to avoid the assumption that all sources are 
  ! identical - e.g. we cater for some fields potentially being from 
  ! model components and some fields from ancillary or restart files 
  ! and some fields using EXPOUT in the namcouple file all at the 
  ! same timestep.
  IF (oasis_get_ok(infosst)) THEN

    ! Our first stage of processing deals with incoming seaice data
    ! and thus depends on us having received the necessary fields.
    IF (oasis_get_ok(infofreezen) .AND.          &
        oasis_get_ok(infohicen) .AND.            &
        oasis_get_ok(infosnowthickn)) THEN

      DO k=1,nice
        DO j = 1, oasis_jmt
          DO i = 1, oasis_imt

            IF (l_sice_multilayers) THEN
              ! Ensure ice temperature (in K) and effective conductivity are
              ! always positive.

              ! Ensure consistency between location of ice and ice temperature
              ! and conductivity fields

              IF (ocn_freezen(i,j,k) >  0.0 .AND.                  &
                  (ocn_icetn(i,j,k) == 0.0 .OR. ocn_icekn(i,j,k) == 0.0)) THEN

                CALL umPrint( 'PROBLEM IN OASIS3_GETO2A',src='oasis3_geto2a')
                CALL umPrint( 'Missing ice T or K data',src='oasis3_geto2a')
                WRITE(umMessage,'(A,1X,I4,1X,I4,1X,I4)') 'i,j,k=',i,j,k

                CALL umPrint(umMessage,src='oasis3_geto2a')

                WRITE(umMessage,'(A,1X,E18.10,1X,E18.10,1X,E18.10)')     &
                  'ocn_freezen,ocn_icetn,ocn_icekn=',                    &
                  ocn_freezen(i,j,k),ocn_icetn(i,j,k),ocn_icekn(i,j,k)

                CALL umPrint(umMessage,src='oasis3_geto2a')

                !!!! Remove ice for now
                ocn_freezen(i,j,k) = 0.0 ! other fields will be zeroed below
              END IF

            END IF ! multilayers is true

            ! Here we ensure that regardless of what has happened to the ice
            ! fraction and GBM ice depth during the coupling,
            ! the minimum local ice depth is 1cm (as in CICE). This is done by
            ! squashing the ice to a smaller area whilst conserving volume (so
            ! the GBM ice depth is unchanged).
            !
            ! Note that freezen is unchanged unless the local ice depth (i.e.
            ! hicen/freezen is less than 1.0e-02, in which case freezen is
            ! decreased so that hicen/freezen=1.0e-02

            ! If using multilayers....
            ! **CHANGE the way min local ice depth is enforced and increase 
            ! min ice depth to hi_min, set at top.  Get rid of ice below that 
            ! value - do not squash
            ! **INCREASE min ice concentration
            !
            ! Do not apply these limits to first-order ice concentration,
            ! which is used to produce an energy-conserving field of local 
            ! fluxes going to the ocean

            IF (l_sice_multilayers) THEN
              hi_min = 0.2
              IF (ocn_freezen(i,j,k) > 0.0) THEN
                IF (ocn_hicen(i,j,k)/ocn_freezen(i,j,k) < hi_min) THEN
                  ocn_freezen(i,j,k) = 0.0  ! other fields will be zeroed next
                END IF
              END IF
              aicenmin=1.0e-02
            ELSE  ! old method
              hi_min = 0.01
              ocn_freezen(i,j,k)=MIN(ocn_freezen(i,j,k),            &
               (ocn_hicen(i,j,k)/hi_min))
              aicenmin=2.0e-04
            END IF

            ! Also impose a minimum ice fraction
            IF (ocn_freezen(i,j,k) < aicenmin) THEN

              ocn_freezen(i,j,k)=0.0
              ocn_hicen(i,j,k)=0.0
              ocn_snowthickn(i,j,k)=0.0
              IF (l_sice_multilayers) THEN
                ocn_icetn(i,j,k)=0.0
                ocn_icekn(i,j,k)=0.0
              END IF
              IF (l_sice_meltponds_cice) THEN
                pond_frac_n(i,j,k)=0.0
                pond_depth_n(i,j,k)=0.0
              END IF 
            END IF

            ocn_freeze(i,j)=ocn_freeze(i,j)+ocn_freezen(i,j,k)

          END DO  ! Over I
        END DO  ! Over J
      END DO  ! Over K

      DO k=1,nice
        DO j = 1, oasis_jmt
          DO i = 1, oasis_imt

            ! Also reduce category ice fractions if aggregate > 1
            ! We don't do this to first-order ice concentration
            ! because it might mess up conservation of energy
            IF (ocn_freeze(i,j) > 1.0) THEN  
              ocn_freezen(i,j,k)=ocn_freezen(i,j,k)/ocn_freeze(i,j)
            END IF

            ocn_hicen(i,j,k)=ocn_hicen(i,j,k)+                    &
               ocn_snowthickn(i,j,k)*(kappai/kappai_snow)
            ocn_hice(i,j)=ocn_hice(i,j)+ocn_hicen(i,j,k)
            ocn_snowthick(i,j)=ocn_snowthick(i,j)+                &
               ocn_snowthickn(i,j,k)

            IF (ocn_freezen(i,j,k) >  0.0) THEN
              ocn_hicen(i,j,k)=ocn_hicen(i,j,k)/                 &
                 ocn_freezen(i,j,k)
              ocn_snowthickn(i,j,k)=ocn_snowthickn(i,j,k)*       &
                 rhosnow/ocn_freezen(i,j,k)
            END IF

            IF (l_sice_multilayers) THEN
              ! Ice T and K fields are passed as GBM (i.e. multiplied by
              ! category ice concentration).  Convert them back to local
              ! calues here.
              IF (ocn_freezen(i,j,k) >  0.0) THEN
                ocn_icetn(i,j,k)=ocn_icetn(i,j,k)/ocn_freezen(i,j,k)
                ocn_icekn(i,j,k)=ocn_icekn(i,j,k)/ocn_freezen(i,j,k)
                ! ELSE not needed as fields will have been zeroed above
              END IF
            END IF

          END DO    ! Over I
        END DO    ! Over J
      END DO   ! Over K

      DO j = 1, oasis_jmt
        DO i = 1, oasis_imt
          IF (ocn_freeze(i,j) >  0.0) THEN
            ocn_freeze(i,j)=MIN(ocn_freeze(i,j),1.0)
            ocn_snowthick(i,j)=ocn_snowthick(i,j)*                &
             rhosnow/ocn_freeze(i,j)
            ocn_hice(i,j)=ocn_hice(i,j)/ocn_freeze(i,j)
          END IF
        END DO
      END DO

      IF (l_sice_multilayers) THEN

        ! Calculate GBM icet and icek - ignoring open water
        ! Don't need to do this if nice=1 as the pointer for TI_GB and TI will
        ! be indentical and likewise for ICE_K and ICE_K_CAT
        IF (nice > 1) THEN
 
          ! Set ocn_icetn to min or max where it exceeds these bounds.
          DO j = 1, oasis_jmt
            DO i = 1, oasis_imt
              DO k = 1,nice
                IF (ocn_icetn(i,j,k) > TI_max) THEN
                  ocn_icetn(i,j,k) = TI_max
                END IF
                IF ((ocn_icetn(i,j,k) < TI_min) .AND. &
                    (ocn_icetn(i,j,k) /= 0.0)) THEN
                   ocn_icetn(i,j,k) = TI_min
                END IF
               END DO
             END DO
          END DO
          
          DO j = 1, oasis_jmt
            DO i = 1, oasis_imt
              IF (ocn_freeze(i,j) >  0.0) THEN
                ocn_icet_gb(i,j) = 0.0
                DO k = 1,nice
                  ocn_icet_gb(i,j) = ocn_icet_gb(i,j) +                   &
                            (ocn_freezen(i,j,k)/ocn_freeze(i,j) ) *     &
                             ocn_icetn(i,j,k)
                END DO
              ELSE
                ocn_icet_gb(i,j)=ocn_icetn(i,j,1) ! zero (open ocean) or land
              END IF
            END DO
          END DO
        END IF
      END IF

      ! Process surface temperature fields here (as used to be done
      ! in SWAPA20 pre UM version 7.0).

      IF (nice_use  >   1) THEN
        DO j = 1, oasis_jmt
          DO i = 1, oasis_imt
            IF (fland_loc(i,j) <  1.0 .AND. ocn_freeze(i,j) > 0.0 ) THEN
              tstar_sice_local(i,j) = 0.0
              DO k = 1, nice_use
                tstar_sice_cat_local(i,j,k) =                           &
                              MIN(tstar_sice_cat_local(i,j,k),tm)
                tstar_sice_local(i,j) = tstar_sice_local(i,j) +         &
                              (ocn_freezen(i,j,k)/ocn_freeze(i,j) ) *   &
                              tstar_sice_cat_local(i,j,k)
              END DO
            END IF
          END DO
        END DO
      ELSE
        DO j = 1, oasis_jmt
          DO i = 1, oasis_imt
            IF (fland_loc(i,j) <  1.0 .AND. ocn_freeze(i,j) > 0.0 ) THEN
              tstar_sice_cat_local(i,j,1) =                             &
                              MIN(tstar_sice_cat_local(i,j,1),tm)
              tstar_sice_local(i,j) = tstar_sice_cat_local(i,j,1)
            END IF
          END DO
        END DO
      END IF

      IF (nice_use  >   1) THEN
        DO k = 1, nice_use
          DO j = 1, oasis_jmt
            DO i = 1, oasis_imt
              tstar_sice_cat(i,j,k)=tstar_sice_cat_local(i,j,k)
            END DO
          END DO
        END DO
      END IF

      ! Correct meltpond fields for coupling errors that can result in
      ! meltpond fractions being marginally < 0 or marginally > 1:
      IF (l_sice_meltponds_cice) THEN
        DO k=1,nice
          DO j = 1, oasis_jmt
            DO i = 1, oasis_imt
              IF (pond_frac_n(i,j,k) < 0.0) THEN
                pond_frac_n(i,j,k) = 0.0
                pond_depth_n(i,j,k) = 0.0
              END IF
              IF (pond_frac_n(i,j,k) > 1.0) THEN
                pond_frac_n(i,j,k) = 1.0
              END IF
              IF (pond_depth_n(i,j,k) < 0.0) THEN
                pond_depth_n(i,j,k) = 0.0
              END IF
            END DO
          END DO
        END DO
      END IF

      IF (l_saldep_freeze) THEN
        DO j = 1, oasis_jmt
          DO i = 1, oasis_imt
            sstfrz(i,j) = ocn_sstfrz(i,j)
          END DO
        END DO
      END IF

      DO k = 1,nice
        DO j = 1, oasis_jmt
          DO i = 1, oasis_imt
            ice_fract_cat(i,j,k)=ocn_freezen(i,j,k)
            
            ! Read in 'future' ice concentration (regridded using 1st-order 
            ! remapping) which will be used to divide through sea ice 
            ! fluxes in oasis3_puta2o
            ice_fract_cat_future(i,j,k) = ocn_freezen_first(i,j,k)
            
            ice_thick_cat(i,j,k)=ocn_hicen(i,j,k)
            snodep_sea_cat(i,j,k)=ocn_snowthickn(i,j,k)
            IF (l_sice_meltponds_cice) THEN
              pond_frac_cat(i,j,k)=pond_frac_n(i,j,k)
              pond_depth_cat(i,j,k)=pond_depth_n(i,j,k)
            END IF
          END DO
        END DO
      END DO

      IF (l_sice_multilayers) THEN
        DO k = 1,nice
          DO j = 1, oasis_jmt
            DO i = 1, oasis_imt
              ti_cat(i,j,k)=ocn_icetn(i,j,k)
              ice_k_cat(i,j,k)=ocn_icekn(i,j,k)
            END DO
          END DO
        END DO

        DO j = 1, oasis_jmt
          DO i = 1, oasis_imt
            ti(i,j,1)=ocn_icet_gb(i,j)
          END DO
        END DO
      END IF

    ELSE 

      IF (infofreezen /= PRISM_NotDef .AND. &
          (.NOT.oasis_get_ok(infofreezen))) THEN
        WRITE(Cmessage, "(A,I0)") "Unexpected infofreezen value:",      & 
                                     infofreezen
        icode = 500 + ABS(infofreezen)
        CALL Ereport("OASIS3_geto2a", icode, Cmessage)
      END IF

      IF (infohicen /= PRISM_NotDef .AND. &
          (.NOT.oasis_get_ok(infohicen))) THEN
        WRITE(Cmessage, "(A,I0)") "Unexpected infohicen value:",        & 
                                     infohicen
        icode = 500 + ABS(infohicen)
        CALL Ereport("OASIS3_geto2a", icode, Cmessage)
      END IF

      IF (infosnowthickn /= PRISM_NotDef .AND. &
          (.NOT.oasis_get_ok(infosnowthickn))) THEN
        WRITE(Cmessage, "(A,I0)") "Unexpected infosnowthickn value:",   &
                                     infosnowthickn
        icode = 500 + ABS(infosnowthickn)
        CALL Ereport("OASIS3_geto2a", icode, Cmessage)
      END IF

    END IF ! Our multi cat seaice return codes indicated things to do. 

    !==================================================================
    ! Here we deal with general SST processing.
    ! These operations may involve contributions arising from
    ! seaice, but are coded such that if seaice coupling
    ! is not involved, then the necessary terms (e.g. ocn_freeze)
    ! will already have been initialised to zero and so may appear
    ! in the following calculations without having any impact.
    !==================================================================

    ! The incoming SST is expected in K NOT degrees C.
    ! Any necessary transforms will be done in the PSMILE
    ! controlled via the coupler or in the sending
    ! model component. (i.e. it is not this component's
    ! responsibility to know anything about the sending
    ! component.)
    DO j = 1, oasis_jmt
      DO i = 1, oasis_imt

        ! Special form of masking for T points - if
        ! value returned from coupling is approx abs zero
        ! we assume this is an unchanged point masked during
        ! the coupling process and we reinstate
        ! our original value.
        IF (ocn_sst(i,j) <  1.0) THEN
          ocn_sst(i,j)=ocn_sst_orig(i,j)
        END IF
      END DO
    END DO

    DO j = 1, oasis_jmt
      DO i = 1, oasis_imt
        tstar_ssi(i,j)=ocn_sst(i,j)

        IF (fland_loc(i,j) <  1.0) THEN

          IF (l_leads_temp_prog.eqv..FALSE.) THEN
            IF (ocn_freeze(i,j) >  0.0) ocn_sst(i,j)=tfs
          END IF

          tstar_ssi(i,j)=ocn_freeze(i,j)*tstar_sice_local(i,j)+  &
             (1.0-ocn_freeze(i,j))*ocn_sst(i,j)
        END IF

        tstar_local(i,j)=fland_loc(i,j)*tstar_land_local(i,j)+    &
           (1.0-fland_loc(i,j))*tstar_ssi(i,j)
      END DO
    END DO

    ! Having got our new values we need to move them
    ! to the appropriate arrays where they'll be employed next
    ! time the appropriate calculation is called.
    DO j = 1, oasis_jmt
      DO i = 1, oasis_imt

        IF (l_ctile) THEN
          tstar_sea(i,j)=ocn_sst(i,j)
        END IF

        tstar(i,j)=tstar_local(i,j)

        tstar_sice(i,j,1)=tstar_sice_local(i,j)
        ice_fraction(i,j,1)=ocn_freeze(i,j)
        ice_thickness(i,j,1)=ocn_hice(i,j)
        snodep_sea(i,j,1)=ocn_snowthick(i,j)
      END DO
    END DO

  ELSE 

    ! If we didn't get a valid return code for SST
    IF (infosst /= PRISM_NotDef) THEN
      WRITE(Cmessage, "(A,I0)") "Unexpected infosst value:", infosst
      icode = 500 + ABS(infosst)
      CALL Ereport("OASIS3_geto2a", icode, Cmessage)
    END IF

  END IF ! SST coupling field was received.


  ! Update our surface currents (U,V)
  IF (oasis_get_ok(infou) .AND.                            &
      oasis_get_ok(infov)) THEN

    ! Now we have to make sure that things are
    ! consistent on the North polar row. (We shouldn't care
    ! about the South for ocean and seaice fields because there's 
    ! no coupling over Antarctica for these) To do this, under 
    ! EG, we discard our incoming N polar V values and calculate new
    ! ones based on our "top row" of U values. We need to do this 
    ! because remapped values from the polar singularity 
    ! are rarely uniform and often unphysical especially when 
    ! mapping from a tripolar grid.
    ! Of course we only need to do this for a global model
    ! (or at least a model which covers either pole!)

    IF (model_type == mt_global) THEN

      ! DEPENDS ON: Correct_Polar_UV
      CALL Correct_Polar_UV(                 &
           delta_lambda, global_row_length,  &
           ocn_u,ocn_v,                      &
           oasis_jmt_u,oasis_jmt_v)

    END IF


    ! Incoming currents expected in m/s
    ! We need to index u_sea and v_sea using udims and vdims rather than 1:imt
    ! and 1:jmt to avoid the bodge of changing array start/end indexing by 
    ! passing them through the argument list as was previously done. 
    u_sea(udims%i_start:udims%i_end,udims%j_start:udims%j_end) = &
      ocn_u(1:oasis_imt,1:oasis_jmt_u)
    v_sea(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end) = &
      ocn_v(1:oasis_imt,1:oasis_jmt_v)

  ELSE

    IF (infou /= PRISM_NotDef .AND. (.NOT.oasis_get_ok(infou))) THEN
      WRITE(Cmessage, "(A,I0)") "Unexpected infou value: ", infou
      icode = 500 + ABS(infou)
      CALL Ereport("OASIS3_geto2a", icode, Cmessage) 
    END IF

    IF (infov /= PRISM_NotDef .AND. (.NOT.oasis_get_ok(infov))) THEN
      WRITE(Cmessage, "(A,I0)") "Unexpected infov value: ", infov
      icode = 500 + ABS(infov)
      CALL Ereport("OASIS3_geto2a", icode, Cmessage) 
    END IF

  END IF ! We received incoming u and v

  IF (l_dms_conc) THEN

    ! Check to see if we received DMS concentration via the coupler 
    IF (oasis_get_ok(info_dms_conc)) THEN

      ! DMS concentration from MEDUSA
      DO j = 1, oasis_jmt
        DO i = 1, oasis_imt
          dms_conc(I,J)=dms_conc_in(I,J)
        END DO
      END DO

    ELSE
      IF (info_dms_conc /= PRISM_NotDef) THEN
        ! If we didn't get a valid return code for DMS_CONC
        WRITE(Cmessage, "(A,I0)") "Unexpected info_dms_conc value:", &
                                                         info_dms_conc
        icode = 500 + ABS(info_dms_conc)
        CALL Ereport("OASIS3_geto2a", icode, Cmessage)
      END IF
    END IF
   
  END IF

  IF (l_CO2flux) THEN

    ! Check to see if we received CO2 flux via the coupler 
    IF (oasis_get_ok(info_CO2flux)) THEN

      ! CO2 flux from MEDUSA
      DO j = 1, oasis_jmt
        DO i = 1, oasis_imt
          CO2flux(I,J)=CO2flux_in(I,J)
        END DO
      END DO

    ELSE
      IF (info_CO2flux /= PRISM_NotDef) THEN
        ! If we didn't get a valid return code for CO2flux
        WRITE(Cmessage, "(A,I0)") "Unexpected info_CO2flux value:", &
                                                         info_CO2flux
        icode = 500 + ABS(info_CO2flux)
        CALL Ereport("OASIS3_geto2a", icode, Cmessage)
      END IF
    END IF

  END IF

  IF (l_chloro_conc) THEN


    ! Check to see if we received chlorophyll via the coupler 
    IF (oasis_get_ok(info_chloro_conc)) THEN

      ! Chlorophyll from MEDUSA
      DO j = 1, oasis_jmt
        DO i = 1, oasis_imt
          chloro_sea(I,J)=chloro_conc_in(I,J)
        END DO
      END DO

    ELSE
      IF (info_chloro_conc /= PRISM_NotDef) THEN
        ! If we didn't get a valid return code for chlorophyll
        WRITE(Cmessage, "(A,I0)") "Unexpected info_chloro_conc value:", &
                                                         info_chloro_conc
        icode = 500 + ABS(info_chloro_conc)
        CALL Ereport("OASIS3_geto2a", icode, Cmessage)
      END IF
    END IF
   
  END IF

END IF ! If get_step_oa = true

IF (l_oasis_timers) CALL getends
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE oasis3_geto2a
#endif
