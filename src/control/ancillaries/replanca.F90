! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!   Subroutine replanca in Module replanca_mod --------------------
!
!   Purpose:  Updates ancillary fields as specified ancil_requests array.
!     Tests whether update is required for each field, allowing for
!     dependencies between fields. Uses LOOKUP_ANCILA array to find data for
!     appropriate time, reads a record and checks for current data
!     type. Reads second record if time interpolation required. Updates
!     the field. Under lcal360 the 360 day rather than the Gregorian calender
!     is used.
!     -------------------------------------------------------------
!
!     Code Owner: Please refer to the UM file CodeOwners.txt
!     This file belongs in section: Ancillaries

MODULE replanca_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='REPLANCA_MOD'

CONTAINS

SUBROUTINE replanca(i_year,i_month,i_day,i_hour,                              &
                    i_minute,i_second,i_day_number,                           &
                    ancil_reftime,offset_steps,                               &
                    p_field,p_rows,u_field,v_field,d1,land,                   &
                    a_step,land_field,steps_per_hr,                           &
                    fland_ctile,                                              &
                    tstar_land_ctile,tstar_sea_ctile,                         &
                    tstar_sice_ctile,                                         &
                    ice_fraction,tstar,tstar_anom,                            &
                    sm_levels,dz_soil,smc_updated,                            &
                    ns_space,first_lat,                                       &
                    len1_lookup,len_fixhd,len_inthd,                          &
                    len_realhd,len_d1,                                        &
                    nlookups,                                                 &
                    icode,cmessage,lcal360)        ! Intent Out

USE water_constants_mod, ONLY: tfs, tm
USE mask_compression, ONLY: compress_to_mask

USE ancilcta_namelist_mod, ONLY:                                               &
  l_sstanom, l_amipii_ice_processing, use_lookup_dates_anc_time_interp

USE jules_sea_seaice_mod, ONLY: l_ctile
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE lookup_addresses
USE t_int_mod,   ONLY: t_int
USE t_int_c_mod, ONLY: t_int_c
USE gen_phys_inputs_mod,  ONLY: l_leads_temp_prog
USE ukca_option_mod, ONLY: l_ukca
USE planet_constants_mod, ONLY: l_planet_grey_surface
USE conversions_mod, ONLY: isec_per_day

USE um_stashcode_mod ! Stash codes
USE ancil_mod,      ONLY: num_ancil_requests, ancil_files, ancil_requests, &
                          find_ancil_file_by_stash, find_ancil_req_by_stash
USE ancil_headers_mod, ONLY: fixhd_ancila, inthd_ancila, lookup_ancila,  &
                             lookup_start_ancila
USE ppxlook_mod, ONLY: exppxi

USE cancila_mod, ONLY: nlookup, lookup_step, levels, d1_anciladd, steps, &
                       update

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE missing_data_mod, ONLY: rmdi, imdi

USE ereport_mod, ONLY: ereport

USE errormessagelength_mod, ONLY: errormessagelength

USE readflds_mod, ONLY: readflds

USE items_nml_mod, ONLY: period_years, period_months

IMPLICIT NONE


LOGICAL :: lcal360

INTEGER ::                                                                     &
       i_year,                                                                 &
                          ! Current Model Time
       i_month,                                                                &
                          !   "      "     "
       i_day,                                                                  &
                          !   "      "     "
       i_hour,                                                                 &
                          !   "      "     "
       i_minute,                                                               &
                          !   "      "     "
       i_second,                                                               &
                          !   "      "     "
       i_day_number,                                                           &
       ancil_reftime(6),                                                       &
                          ! Reference time for ancillary updating
       offset_steps,                                                           &
                          ! Offset in timesteps of ref. from basis


       a_step,land_field,steps_per_hr,                                         &


       p_field,                                                                &
                          ! Size of horizontal fields
       p_rows,                                                                 &
                          !
       u_field,                                                                &
                          !   "  "      "         "
       v_field,                                                                &
                          !   "  "      "         "
       nlookups,                                                               &
                          ! Number of lookup tables
       len_d1             ! Size of primary data array

INTEGER ::                                                                     &
       len1_lookup,                                                            &
                          ! First dimension of lookup table
       len_fixhd,                                                              &
                          ! Length of headers in data sets
       len_inthd,                                                              &
                          !
       len_realhd,                                                             &
                          !
       lookup(len1_lookup,nlookups),                                           &
                          ! Data set lookup tables
       sm_levels          ! number of soil levels

REAL ::                                                                        &
       d1(len_d1),                                                             &
                          ! INOUT  Primary data array used to hold
                          !        all fields except TSTAR and
                          !        ICE_FRACTION
       ice_fraction(p_field),                                                  &
                          ! INOUT  Ice frac of sea part of grid
                          !        box, updated if requested
       fland_ctile(land_field),                                                &
                          ! IN     Fractional land on land pts.
       fland_g(p_field),                                                       &
                          ! WORK   Frac land over all points.
       tstar(p_field),                                                         &
                          ! INOUT  TSTAR:updated if requested
       tstar_land_ctile(p_field),                                              &
                          ! INOUT  as above, but for land.
       tstar_sea_ctile(p_field),                                               &
                          ! INOUT  as above, but for open sea.
       tstar_sice_ctile(p_field),                                              &
                          ! INOUT  as above, but for sea-ice.
       tstar_anom(p_field),                                                    &
                          ! INOUT  SST anomaly,formed in recon;
                          !        added if requested in model run
       ns_space,                                                               &
                          ! NS latitude spacing
       first_lat,                                                              &
                          ! Latitude of first gridpoint
       dz_soil(sm_levels) ! OUT soil thicknesses

LOGICAL ::                                                                     &
       land(p_field),                                                          &
                          ! WORK LAND mask
       sea(p_field),                                                           &
                          ! WORK SEA mask
       ltstar_sice,                                                            &
                          ! IN TRUE if TSTAR_SICE has been read in
                          !         from input dump.
                          !         If FALSE set to TSTAR_SEA.
       smc_updated        ! OUT T if smc updated

INTEGER ::                                                                     &
       icode,                                                                  &
                          ! Return code
       iounit             !OUT I/O unit passed out in RECON mode

CHARACTER(LEN=errormessagelength) ::                                           &
       cmessage           ! Error message
! ----------------------------------------------------------------

!     Local real arrays

REAL ::                                                                        &
       ancil1(p_field),                                                        &
                          ! Buffers to hold values of ancillary
                          ! data for time interpolation.
       ancil2(p_field),                                                        &
                          !
       ancil_data(p_field),                                                    &
                          ! Field of ancillary data held prior
                          ! to selective updating.
       snow_change(p_field),                                                   &
                          ! Fractional time of change of
                          ! snow cover
       ice_extent(p_field,2),                                                  &
                          ! Fractional time of change
                          ! of ice cover
       pres_value(p_field),                                                    &
                          ! Prescribed value of data when
                          ! controlling field is zero.
       no_ice_extent(p_field),                                                 &
                          ! Indicator for no sea ice
                          ! =0 if ice cover
       tstar_land(p_field),                                                    &
                          !Temporary store for land surface temp.
       tstar_sea(p_field),                                                     &
                          !as above, but for open sea.
       tstar_sice(p_field),                                                    &
                          !as above, but for sea-ice.
       tstar_ssi(p_field)
                          !as above, but for sea mean.

!     Local variables

INTEGER ::                                                                     &
       i,                                                                      &
                          !
       i1,                                                                     &
                          !
       i2,                                                                     &
                          !
       i3,                                                                     &
       id,                                                                     &
                          !
       im,                                                                     &
                          !
       iy,                                                                     &
                          !
       k,                                                                      &
       l,                                                                      &
                          ! Land index
       file_num               !

INTEGER :: irec          ! Loop counter
INTEGER :: ds_index      ! Identifier for data sets / ancil files
INTEGER :: ds_index1      ! Identifier for data sets / ancil files
INTEGER :: nftin          ! FTN number for ancillary file
INTEGER :: stash_item     ! stash item number for ancillary field

INTEGER ::                                                                     &
       interval,                                                               &
                          ! Interval between data times
       step,                                                                   &
                          ! Number of data times skipped.
       months,                                                                 &
                          ! Used in calculation of position
                          ! of data required.
       hours,                                                                  &
                          !
       period,                                                                 &
                          ! Period of periodic data
       start_month,                                                            &
                          !
       level,                                                                  &
                          !
       ancil_ref_days,                                                         &
                          ! Ancil.reference time in whole days
       ancil_ref_secs,                                                         &
                          ! Ancil.reference time in extra seconds
       day,sec,                                                                &
                          ! Times relative to reference time
       day1,sec1,                                                              &
                          ! Times relative to reference time
       incr_sec,                                                               &
                          ! Increment in sec
       len1,                                                                   &
       iend,                                                                   &
       ii,row_length,j
INTEGER ::                                                                     &
       i_year1,                                                                &
                          ! Copy of Current Model Time year
       i_month1,                                                               &
                          !   "      "     "          month
       i_day1,                                                                 &
                          !   "      "     "          day
       i_hour1            !   "      "     "          hour

INTEGER ::                                                                     &
       update_months      ! update frequency (months) if Gregorian
LOGICAL ::                                                                     &
       lgreg_monthly      ! True for Gregorian monthly updating

! *IF -DEF,CAL360

INTEGER ::                                                                     &
       i_year_basis,                                                           &
                          ! Basis Model Time
       i_month_basis,                                                          &
                          !   "     "     "
       i_day_basis,                                                            &
                          !   "     "     "
       i_hour_basis,                                                           &
                          !   "     "     "
       i_minute_basis,                                                         &
                          !   "     "     "
       i_second_basis,                                                         &
                          !   "     "     "
       i_day_number_basis

! *ENDIF


INTEGER ::                                                                     &
       i_year_ref,                                                             &
                          ! Reference Time
       i_month_ref,                                                            &
                          !    "       "
       i_day_ref,                                                              &
                          !    "       "
       i_hour_ref,                                                             &
                          !    "       "
       i_minute_ref,                                                           &
                          !    "       "
       i_second_ref
                          !    "       "


LOGICAL ::                                                                     &
       linterpolate,                                                           &
                          ! Indicates whether time
                          ! interpolation needed.
       lt_int_c,                                                               &
                          ! Indicates use of controlled time
                          ! interpolation
       lmismatch,                                                              &
                          ! Used in header checks
       lice_fraction,                                                          &
                          !
       lsnow_depth,                                                            &
                          !
       single_time,                                                            &
                          ! Indicates that only one time is
                          ! available in data set
       periodic,                                                               &
                          ! Data set is periodic
       regular,                                                                &
                          ! Interval between data times in
                          ! dataset is regular in model timesteps.
       lice_depth

REAL ::                                                                        &
       zero,                                                                   &
                          !
       time1,                                                                  &
                          ! Times if data used in time interpolation
       time2,                                                                  &
                          !
       targ_time,                                                              &
                          !Target time for time interpolation
       lat_p                                                                  
                          ! latitude of point

INTEGER :: ancil_req_num_icefrac  ! address of sea ice fraction in ancil request
                                  ! array
INTEGER :: ancil_req_num_surftemp ! address of surface temp in ancil request 
                                  ! array

INTEGER :: ftn_icefrac       ! fortran unit for sea ice fraction

INTEGER :: file_num_icefrac  ! address for sea ice fraction in 
                             ! lookup_start_ancila
INTEGER :: file_num_surftemp ! address for sea ice fraction in 
                             ! lookup_start_ancila

INTEGER :: lb_code
INTEGER :: bot_level

REAL, PARAMETER :: dummy_real = 1.0 ! a REAL number- used only by TRANSFER to 
                                    ! interrogate the REAL representation of 
                                    ! lookup_ancila


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='REPLANCA'

!   1.  Initialisation for atmosphere
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode=0
iounit=0
smc_updated=.FALSE.
update_months=0
incr_sec = 0

! These are printed out and possibly not always required.
hours    = imdi
interval = imdi
period   = imdi

!     Set up surface temperatures:

IF (l_ctile) THEN
  DO i=1,p_field
    tstar_land(i)=tstar_land_ctile(i)
    tstar_sea(i)=tstar_sea_ctile(i)
    tstar_sice(i)=tstar_sice_ctile(i)
    IF (ice_fraction(i) <= 0.0) THEN
      tstar_ssi(i)=tstar_sea(i)
    ELSE
      tstar_ssi(i)=ice_fraction(i)*tstar_sice(i)                               &
        +(1.0-ice_fraction(i))*tstar_sea(i)
    END IF
  END DO
ELSE
  DO i=1,p_field
    tstar_land(i)=tstar(i)
    tstar_ssi(i)=tstar(i)
  END DO
END IF


!     Initialise ANCIL1/2. Includes Halos for MPP runs.
DO i=1,p_field
  ancil1(i)=0.0
  ancil2(i)=0.0
END DO
!   1.1 Set logical UPDATE for each ancillary field independently

DO irec=1,num_ancil_requests

  update(irec)=.FALSE.
  IF (steps(irec) /= 0) THEN
    !         UPDATE(IREC)=MOD(A_STEP,STEPS(IREC)) == 0
    update(irec)=(MOD(a_step+offset_steps,steps(irec)) == 0                  &
                   .OR. a_step == 0)                                         &
                    .AND. ancil_requests(irec)%period >  0                   &
                     .AND. d1_anciladd(irec) >  1
  END IF

  !   1.05 Copy ancillary updating reference time to local variables
  i_year_ref   = ancil_reftime(1)
  i_month_ref  = ancil_reftime(2)
  i_day_ref    = ancil_reftime(3)
  i_hour_ref   = ancil_reftime(4)
  i_minute_ref = ancil_reftime(5)
  i_second_ref = ancil_reftime(6)
  !        and convert to reference days & secs
  ! DEPENDS ON: time2sec
  CALL time2sec(i_year_ref,i_month_ref,i_day_ref,                              &
                i_hour_ref,i_minute_ref,i_second_ref,                          &
                0,0,ancil_ref_days,ancil_ref_secs,lcal360)

  IF (.NOT. lcal360) THEN

    !   1.11 Set logical UPDATE for Gregorian calender updates at monthly
    !        or yearly intervals. NB STEPS value set to 1 day in INANCILA

    IF (ancil_requests(irec)%period == period_years .OR. & 
         ancil_requests(irec)%period == period_months) THEN

      months=i_month+i_year*12-(i_month_ref+i_year_ref*12)

      update_months= ancil_requests(irec)%interval*       &
        ((3-ancil_requests(irec)%period)/2 *12 +          & 
         1-(3-ancil_requests(irec)%period)/2)

      update(irec)=MOD(months,update_months) == 0 .AND. i_day == 1

    END IF

  END IF !  (.NOT.LCAL360)

END DO

!  1.2 Allow for dependencies between fields

! Sea surface temperature must be updated when sea ice is updated
i1 = 0
i2 = 0
DO irec=1,num_ancil_requests

  IF (ancil_requests(irec)%stashcode == stashcode_icefrac)  i1 = irec
  IF (ancil_requests(irec)%stashcode == stashcode_surftemp) i2 = irec

  IF (i1 > 0 .AND. i2 > 0) THEN
    update(i2) = update(i1) .OR. update(i2)
    ancil_req_num_icefrac = i1  ! for later usage
    EXIT
  END IF

END DO

! Both surface current components must be updated together
i1 = 0
i2 = 0
DO irec=1,num_ancil_requests

  IF (ancil_requests(irec)%stashcode == stashcode_surf_z_curr) i1 = irec
  IF (ancil_requests(irec)%stashcode == stashcode_surf_m_curr) i2 = irec

  IF (i1 > 0 .AND. i2 > 0) THEN
    update(i1) = update(i1) .OR. update(i2)
    update(i2) = update(i1)
    EXIT
  END IF

END DO

!  Select method of time interpolation for SST. The interpolation
!  allows for sea ice if ice data is available at the same times
!  as the temperature data. Otherwise linear interpolation is used.
lt_int_c=.TRUE.

ancil_req_num_surftemp = find_ancil_req_by_stash(stashcode_surftemp)
! Check if ancil request exists for surface_temp. Return value of IMDI
! indicates no request found.
IF (ancil_req_num_surftemp /= imdi) THEN

  ! As ancil request exists, is a surface_temp update required this timestep?
  IF (update(ancil_req_num_surftemp)) THEN
    file_num_surftemp = find_ancil_file_by_stash(stashcode_surftemp)
    file_num_icefrac  = find_ancil_file_by_stash(stashcode_icefrac)
    ftn_icefrac   = ancil_files(file_num_icefrac) % unit_num
    ! Turn off time controlled interpolation if ice fraction ancil file time 
    ! indicator suggests a single time field
    IF (fixhd_ancila(10,file_num_icefrac) == 0) lt_int_c=.FALSE.
  
    IF (lt_int_c) THEN
      DO i=21,41
        IF (fixhd_ancila(i,file_num_icefrac) /= &
             fixhd_ancila(i,file_num_surftemp)) THEN
          lt_int_c=.FALSE.
          WRITE(umMessage,'(A,A,A)')' WARNING:controlled time interp for SST', &
        ' not available: Mismatch in SST and SEA-ICE ancillary data',          &
        ' times in FIXED HEADER'
          CALL umPrint(umMessage,src='replanca')
          WRITE(umMessage,'(A,I4,A,I0)')' position=',i,' SEA-ICE=', &
               fixhd_ancila(i,file_num_icefrac)
          CALL umPrint(umMessage,src='replanca')
          WRITE(umMessage,'(A,I4,A,I0)')' position=',i,' SST=',     &
               fixhd_ancila(i,file_num_surftemp)
          CALL umPrint(umMessage,src='replanca')
        END IF
      END DO
    END IF
  END IF
END IF


! Read in fractional land field
! Set up global fractional land field
IF (l_ctile) THEN
  l=0
  DO i=1,p_field
    fland_g(i)=0.0
    IF (land(i)) THEN
      l=l+1
      fland_g(i)=fland_ctile(l)
    END IF
  END DO
ELSE
  DO i=1,p_field
    IF (land(i)) THEN
      fland_g(i)=1.0
    ELSE
      fland_g(i)=0.0
    END IF
  END DO
END IF

DO i=1,p_field
  sea(i)=.FALSE.
  IF (fland_g(i) <  1.0)sea(i)=.TRUE.
END DO

!  Loop over ancillary fields(atmosphere)

DO irec = 1, num_ancil_requests

  stash_item = ancil_requests(irec)%stashcode
  lice_depth = stash_item == stashcode_icethick  ! required for LAMIPII

  IF (.NOT. update(irec)) THEN
    CYCLE
  END IF
  nftin = ancil_files(ancil_requests(irec)%ancil_file_num)%unit_num
  file_num  = ancil_requests(irec)%ancil_file_num

  IF (lice_depth .AND. l_amipii_ice_processing) THEN

    ! Uses ice fraction set earlier in field loop.
    ! WARNING this will fail if the order of ancillary fields is ever
    !         changed so that ice-depth precedes ice fraction
    ! Note : For complete sea ice cover
    !        Arctic ice depth    = 2m
    !        Antarctic ice depth = 1m
    ! For ice concentrations less than 1. ice depth is 1 or 2 times conc.
    ! This results in similar values to those from runs using ancillary
    ! files containing ice depths set to 1 or 2m.

    row_length=p_field/p_rows
    DO i=1,p_rows
      ! work out latitude in radians
      lat_p=first_lat+ns_space*(i+datastart(2)-offy-1)
      DO j=1,row_length
        ii=j+(i-1)*row_length
        ancil_data(ii)=0.0
        IF (ice_fraction(ii) >  0.0) THEN
          IF (lat_p >  0.0) THEN   ! Arctic ice depth
            ancil_data(ii)=2.0*ice_fraction(ii)
          ELSE                     ! Antarctic ice depth
            ancil_data(ii)=1.0*ice_fraction(ii)
          END IF
        END IF
      END DO
    END DO
    !     Sea ice thickness
    !        Update over all sea points (all sea ice points are the only
    !        ones strictly required, but this cannot be determined easily)

    DO i=1,p_field
      IF (sea(i)) THEN
        d1(d1_anciladd(irec)+i-1)=ancil_data(i)
      END IF
    END DO
    CYCLE ! Finished with ice_thickness field, now cycle main loop to next field
  END IF
  !   Update required for field

  WRITE(umMessage,*)'REPLANCA: UPDATE REQUIRED FOR FIELD',irec
  CALL umPrint(umMessage,src='replanca')

  IF ( fixhd_ancila(10,file_num) < 0 .OR. fixhd_ancila(10,file_num)  >  2) THEN
    icode = 700 + irec
    cmessage = 'REPLANCA: error in fixed header(10) of ancillary file'
    GO TO 9999
  END IF

  !     Check whether more than one data time available in data set

  single_time = fixhd_ancila(10,file_num) == 0

  !     Set default values for time interpolation

  linterpolate=.TRUE.
  IF (single_time) THEN
    linterpolate=.FALSE.
  END IF

  IF (stash_item == stashcode_soil_temp) linterpolate=.FALSE.
  IF (stash_item >= stashcode_vol_smc_wilt .AND.                          &
      stash_item <= stashcode_thermal_conduct) linterpolate=.FALSE.
  IF (stash_item == stashcode_iceberg_calving) linterpolate=.FALSE.

  !   2.1 Find position of input record

  !   Default settings of search parameters if only one time present

  IF (single_time) THEN
    step=0
  ELSE

    lgreg_monthly=.FALSE.

    IF (.NOT. lcal360) THEN
      IF (ancil_requests(irec)%period == period_years .OR.  &
           ancil_requests(irec)%period == period_months) THEN
        lgreg_monthly=.TRUE.
        update_months= ancil_requests(irec)%interval*        &
        ((3-ancil_requests(irec)%period)/2 *12 +             &
        1-(3-ancil_requests(irec)%period)/2)
      END IF
    END IF

    periodic = fixhd_ancila(10,file_num) == 2
    regular = .TRUE.

    IF (.NOT. lcal360) THEN
      regular = fixhd_ancila(35,file_num) == 0 .AND. &
           fixhd_ancila(36,file_num) == 0
      ! i.e. data at intervals of days/hours & non-periodic
      IF (periodic) regular = regular .AND. fixhd_ancila(37,file_num) == 0
      ! i.e. data at intervals of hours & periodic
    END IF

    ! Error checking on time information.

    IF ( fixhd_ancila(35,file_num)  <   0 .OR.                                 &
         fixhd_ancila(36,file_num)  <   0 .OR.                                 &
         fixhd_ancila(36,file_num)  >  12 .OR.                                 &
         regular .AND. ( fixhd_ancila(37,file_num) <  0 .OR.                   &
         fixhd_ancila(37,file_num)  > 31  .OR.                                 &
         fixhd_ancila(38,file_num)  <   0 .OR.                                 &
         fixhd_ancila(38,file_num) >  24 ) ) THEN
      ! FIXHD_ANCILA(39-40) are not used by REPLANCA.
      ! FIXHD_ANCILA(35-37) have already been used if not CAL360.
      icode = 700 + irec
      WRITE(cmessage,*) 'REPLANCA: Error in validity time interval',       &
                        ' given in ancil file (FIXHD_ANCILA(35-38))'
      GO TO 9999
    END IF

    IF ( fixhd_ancila(21,file_num)  <   0 .AND. .NOT. periodic                 &
         .OR. .NOT. ( regular .AND. periodic ) .AND.                       &
      !If it is REGULAR & PERIODIC more detailed check is applied below
     ( fixhd_ancila(22,file_num)  <   0  .OR.                                  &
       fixhd_ancila(22,file_num)  >   12 .OR.                                  &
       fixhd_ancila(23,file_num)  <   0  .OR.                                  &
       fixhd_ancila(23,file_num)  >   31 .OR.                                  &
       fixhd_ancila(24,file_num)  <   0  .OR.                                  &
       fixhd_ancila(24,file_num)  >   24 .OR.                                  &
       fixhd_ancila(25,file_num)  <   0  .OR.                                  &
       fixhd_ancila(25,file_num)  >   60 .OR.                                  &
       fixhd_ancila(26,file_num)  <   0  .OR.                                  &
       fixhd_ancila(26,file_num)  >   60 ) ) THEN
      icode = 700 + irec
      WRITE(cmessage,*) 'REPLANCA: Error in first validity time given',    &
                    ' in ancillary file (FIXHD_ANCILA(21-26))'
      GO TO 9999
    END IF

    period=-1
    IF (.NOT. periodic) THEN

      !  If data taken from full time series of input data.

      ! DEPENDS ON: time2sec
      CALL time2sec(i_year,i_month,i_day,i_hour,i_minute,i_second          &
                   ,ancil_ref_days,ancil_ref_secs,day,sec ,lcal360)

      !  Adjust time to middle of updating interval

      IF (.NOT. lgreg_monthly) THEN
        sec=sec+steps(irec)*1800/steps_per_hr

        ! If start-up, adjust for offset of reference time from initial time,
        ! & update with values for half a period before first standard 
        ! update.
        IF (a_step == 0) THEN
          day1 = day
          sec1 = sec
          incr_sec=-3600*MOD(offset_steps,steps(irec))/steps_per_hr
          ! DEPENDS ON: time_df
          CALL time_df(day1,sec1,0,incr_sec,day,sec)
        END IF

      ELSE
        im=MOD(i_month+update_months-1,12) + 1
        iy=i_year+(i_month+update_months-1)/12
        ! DEPENDS ON: time2sec
        CALL time2sec(iy,im,i_day,i_hour ,i_minute,i_second                &
                     ,ancil_ref_days,ancil_ref_secs,day1,sec1,lcal360)
        IF (MOD(day+day1,2) == 0) THEN
          day=(day+day1)/2
          sec=(sec+sec1)/2
        ELSE
          day=(day+day1-1)/2
          sec=(sec+sec1+isec_per_day)/2
        END IF
        ! If start-up, adjust for offset of reference time from
        ! initial time, & update with values for half a period before
        ! first standard update.
        IF (a_step == 0) THEN
          day1 = day
          sec1 = sec
          incr_sec=-3600*MOD(offset_steps,steps(irec))/steps_per_hr
          ! DEPENDS ON: time_df
          CALL time_df(day1,sec1,0,incr_sec,day,sec)
        END IF
      END IF


      IF (regular) THEN
        !  2.1.1  Standard cases:360 day calender;
        !         or Gregorian calendar with
        !         interval between data times in days or hours
        !         updating interval may be regular in model timesteps,
        !         or (LGREG_MONTHLY=T) irregular in model timesteps,

        hours=sec/3600+day*24
        !  FInd time(in hours) of first ancillary data on file
        ! DEPENDS ON: time2sec
        CALL time2sec(fixhd_ancila(21,file_num),fixhd_ancila(22,file_num),    &
                      fixhd_ancila(23,file_num),fixhd_ancila(24,file_num),    &
                      fixhd_ancila(25,file_num),fixhd_ancila(26,file_num),    &
                      ancil_ref_days,ancil_ref_secs,day,sec, lcal360)
        hours=hours-sec/3600-day*24

        IF (hours <  0) THEN
          icode=400+irec
          cmessage='REPLANCA: Current time precedes start time of data'
          GO TO 9999
        END IF

        !  Find interval(in hours) between ancillary data on file
        interval=fixhd_ancila(35,file_num)*8640+fixhd_ancila(36,file_num)*720+ &
                 fixhd_ancila(37,file_num)*24+fixhd_ancila(38,file_num)

        ! Do not interpolate in time if data time exactly matches model time

        IF (MOD(hours,interval) == 0) THEN
          linterpolate=.FALSE.
        END IF

        step=hours/interval
        targ_time=REAL(hours)
        time1=step*interval
        time2=(step+1)*interval

      ELSE

        ! 2.1.2 Gregorian calender;ancillary data interval is in months or
        !       years,which is irregular in model timesteps.
        ! original code is inaccurate for this section - corrected code
        ! under use_lookup_dates_anc_time_interp makes use of dates in lookup 
        ! headers. For a real calendar year the mid-point of each month is 
        ! different in terms of its hour and day. The old inaccurate method 
        ! assumes the hour and day are taken from the fixhd_ancila values. 
        ! These are only usually correct for the first month on the ancillary
        ! file.

        !  Adjust YMD time to middle of updating interval

        i_year1=i_year
        i_month1=i_month
        i_day1=i_day
        i_hour1=i_hour
        ! DEPENDS ON: sec2time
        CALL sec2time(day,sec,ancil_ref_days,ancil_ref_secs,               &
                      i_year,i_month,i_day,                                &
                      i_hour,i_minute,i_second,i_day_number, lcal360)

        !  Find interval(in months) between ancillary data on file
        interval=fixhd_ancila(35,file_num)*12+fixhd_ancila(36,file_num)
        months=i_year*12+i_month
        start_month=fixhd_ancila(21,file_num)*12+fixhd_ancila(22,file_num)
        months=months-start_month
        !  Check for time within month
        ! corrected code uses pp lookup header
        IF (use_lookup_dates_anc_time_interp) THEN   
          step=months/interval
          i2=nlookup(irec)+lookup_step(irec)*step
          i1=i2+lookup_start_ancila(file_num,1)-1
          ! Check against day and hour of actual lookup header
          !   not first field
          IF ((i_day*24+i_hour) <                                          &
               (lookup_ancila(3,i1)*24+lookup_ancila(4,i1))) THEN
            months=months-1
          END IF
        ELSE              ! old less accurate code uses FIXHD_ANCILA
          IF ((i_day*24+i_hour) <                                          &
               (fixhd_ancila(23,file_num)*24+fixhd_ancila(24,file_num))) THEN
            months=months-1
          END IF
        END IF ! use_lookup_dates_anc_time_interp

        IF (months <  0) THEN
          icode=400+irec
          cmessage='REPLANCA: Current time precedes start time of data'
          GO TO 9999
        END IF

        !  Adjust YMD time back to start of updating interval

        i_year=i_year1
        i_month=i_month1
        i_day=i_day1
        i_hour=i_hour1

        step=months/interval

        IF (use_lookup_dates_anc_time_interp) THEN       ! corrected code
          targ_time=REAL(sec)/3600+REAL(day*24)
          ! correct calculation of dates uses lookup table dates
          !   not fixhd_ancila date
          i2=nlookup(irec)+lookup_step(irec)*step
          i1=i2+lookup_start_ancila(file_num,1)-1
          i_year1=lookup_ancila(1,i1)
          i_month1=lookup_ancila(2,i1)
          i_day1=lookup_ancila(3,i1)
          i_hour1=lookup_ancila(4,i1)
          ! DEPENDS ON: time2sec
          CALL time2sec(i_year1,i_month1,i_day1,i_hour1,                     &
                        fixhd_ancila(25,file_num),fixhd_ancila(26,file_num), &
                        ancil_ref_days,ancil_ref_secs,day,sec,               &
                        lcal360)
          time1=REAL(sec)/3600+REAL(day*24)
          ! I1+LOOKUP_STEP(irec) correct pointer to next field as
          ! multi-level fields are possible.
          ! First check next lookup is within lookups
          IF (i1+lookup_step(irec) < nlookups) THEN
            ! Second check next lookup is the same field.
            IF (lookup_ancila(item_code,i1) ==                             &
              lookup_ancila(item_code,i1+lookup_step(irec))) THEN
              i_year1=lookup_ancila(1,i1+lookup_step(irec))
              i_month1=lookup_ancila(2,i1+lookup_step(irec))
              i_day1=lookup_ancila(3,i1+lookup_step(irec))
              i_hour1=lookup_ancila(4,i1+lookup_step(irec))
            ELSE
              icode = 500
              WRITE(cmessage,'(A,I5,A,I5)')                                &
              'REPLANCA: error finding next lookup entry for field:',      &
                 irec, ' stashcode ', ancil_requests(irec)%stashcode
              GO TO 9999
            END IF
          ELSE
            icode = 500
            WRITE(cmessage,'(A,I5,A,I5)')                                  &
            'REPLANCA: error due to lookup out of range for field:',       &
               irec, ' stashcode ', ancil_requests(irec)%stashcode
            GO TO 9999
          END IF

          ! DEPENDS ON: time2sec
          CALL time2sec(i_year1,i_month1,i_day1,i_hour1,                      &
                        fixhd_ancila(25,file_num),fixhd_ancila(26,file_num),  &
                        ancil_ref_days,ancil_ref_secs,day,sec, lcal360)
          time2=REAL(sec)/3600+REAL(day*24)

        ELSE  ! use_lookup_dates_anc_time_interp test - old inaccurate code 
              ! using FIXHD_ANCILA
          ! NB INTERVAL may be > 1 month
          months=step*interval
          ! Calculate data times for time interpolation
          targ_time=REAL(sec)/3600+REAL(day*24)
          im=MOD(fixhd_ancila(22,file_num)+months-1,12)+1
          iy=fixhd_ancila(21,file_num)+(months+fixhd_ancila(22,file_num)-1)/12
          ! DEPENDS ON: time2sec
          CALL time2sec(iy,im,fixhd_ancila(23,file_num),                     &
                        fixhd_ancila(24,file_num),                           &
                        fixhd_ancila(25,file_num),fixhd_ancila(26,file_num), &
                        ancil_ref_days,ancil_ref_secs,day,sec,               &
                        lcal360)
          time1=REAL(sec)/3600+REAL(day*24)
          im=MOD(fixhd_ancila(22,file_num)+months+interval-1,12)+1
          iy=fixhd_ancila(21,file_num)+                                      &
               (months+interval+fixhd_ancila(22,file_num)-1)/12
          ! DEPENDS ON: time2sec
          CALL time2sec(iy,im,fixhd_ancila(23,file_num),                     &
                        fixhd_ancila(24,file_num),                           &
                        fixhd_ancila(25,file_num),fixhd_ancila(26,file_num), &
                        ancil_ref_days,ancil_ref_secs,day,sec, lcal360)
          time2=REAL(sec)/3600+REAL(day*24)
        END IF     ! end use_lookup_dates_anc_time_interp test

        ! Do not interpolate in time if data time
        !   exactly matches model time
        IF (targ_time == time1) THEN
          linterpolate=.FALSE.
        END IF

      END IF ! End of REGULAR/not REGULAR

    ELSE  ! PERIODIC data

      !  2.2   If data is taken from ancillary periodic data.

      ! DEPENDS ON: time2sec
      CALL time2sec(i_year,i_month,i_day,i_hour,i_minute,i_second,         &
                    ancil_ref_days,ancil_ref_secs,day,sec, lcal360)

      !  Adjust time to middle of updating interval

      IF (.NOT. lgreg_monthly) THEN
        sec=sec+steps(irec)*1800/steps_per_hr

        !  If start-up, adjust for offset of reference time
        !  from initial time, & update with values for
        !  half a period before first standard update.
        IF (a_step == 0) THEN
          day1 = day
          sec1 = sec
          incr_sec=-3600*MOD(offset_steps,steps(irec))/steps_per_hr
          ! DEPENDS ON: time_df
          CALL time_df(day1,sec1,0,incr_sec,day,sec)
        END IF

      ELSE
        im=MOD(i_month+update_months-1,12) + 1
        iy=i_year+(i_month+update_months-1)/12
        ! DEPENDS ON: time2sec
        CALL time2sec(iy,im,i_day,i_hour ,i_minute,i_second               &
                     ,ancil_ref_days,ancil_ref_secs,day1,sec1,lcal360)
        IF (MOD(day+day1,2) == 0) THEN
          day=(day+day1)/2
          sec=(sec+sec1)/2
        ELSE
          day=(day+day1-1)/2
          sec=(sec+sec1+isec_per_day)/2
        END IF

        !  If start-up, adjust for offset of reference time
        !  from initial time, & update with values for
        !  half a period before first standard update.
        IF (a_step == 0) THEN
          day1 = day
          sec1 = sec
          incr_sec=-3600*MOD(offset_steps,steps(irec))/steps_per_hr
          ! DEPENDS ON: time_df
          CALL time_df(day1,sec1,0,incr_sec,day,sec)
        END IF
      END IF


      !  Adjust YMD time to middle of updating interval

      i_year1=i_year
      i_month1=i_month
      i_day1=i_day
      i_hour1=i_hour
      ! DEPENDS ON: sec2time
      CALL sec2time(day,sec,ancil_ref_days,ancil_ref_secs,                &
                    i_year,i_month,i_day,                                 &
                    i_hour,i_minute,i_second,i_day_number,lcal360)

      IF (regular) THEN
        ! 2.2.1 Standard cases:
        !  1) 360 day calender, with allowed periods of
        !     1 day, 1 month or 1 year;
        !  2) Gregorian calender with update in hours,
        !     and period of data 1 day.
        !  For both updating interval and number of
        !  data times to be skipped in data set calculated in hours.

        hours=sec/3600+day*24
        interval=fixhd_ancila(35,file_num)*8640+fixhd_ancila(36,file_num)*720+ &
                 fixhd_ancila(37,file_num)*24+fixhd_ancila(38,file_num)

        period=inthd_ancila(3,file_num)*interval

        !    Do not allow non-standard periods
        IF (lcal360) THEN
          IF (period /= 8640 .AND. period /= 720 .AND. period /= 24) THEN
            icode=600+irec
            cmessage='REPLANCA: Non-standard period for periodic data'
            GO TO 9999
          END IF
        ELSE
          IF (period /= 24) THEN
            icode=600+irec
            cmessage='REPLANCA: Non-standard period for periodic data'
            GO TO 9999
          END IF
        END IF

        IF (period == 24) THEN
          ! Ancillary data interval in hour(s), period is 1 day

          iy=i_year
          im=i_month
          id=i_day
          IF (i_hour <  fixhd_ancila(24,file_num)) hours=hours+24

        ELSE IF (period == 720) THEN
          ! Ancillary data interval in day(s) or hours , period is 1 month

          iy=i_year
          im=i_month
          id=fixhd_ancila(23,file_num)
          IF ((i_day*24+i_hour) <  (fixhd_ancila(23,file_num)*24+              &
                                    fixhd_ancila(24,file_num))) THEN
             hours=hours+720
          END IF

        ELSE IF (period == 8640) THEN
          ! Ancillary data interval in month(s) or days or hours,
          !   period is 1 year

          iy=i_year
          im=fixhd_ancila(22,file_num)
          id=fixhd_ancila(23,file_num)
          IF ((i_month*720+i_day*24+i_hour) <                              &
             (fixhd_ancila(22,file_num)*720+fixhd_ancila(23,file_num)*24+  &
                                        fixhd_ancila(24,file_num))) THEN
             hours=hours+8640
          END IF
        END IF

        ! DEPENDS ON: time2sec
        CALL time2sec(iy,im,id,fixhd_ancila(24,file_num),                      &
                      fixhd_ancila(25,file_num),fixhd_ancila(26,file_num),     &
                      ancil_ref_days,ancil_ref_secs,day,sec,                   &
                      lcal360)
        hours=hours-sec/3600-day*24

        ! Do not interpolate in time if data time exactly matches model time

        IF (MOD(hours,interval) == 0) THEN
          linterpolate=.FALSE.
        END IF
        step=hours/interval
        targ_time=REAL(hours)
        time1=step*interval
        time2=(step+1)*interval

      ELSE  ! non regular case

        !  2.2.2 Gregorian calender,and data interval is in months,
        !        period is 1 year
        !        Updating interval and number of data times to be skipped
        !        calculated in months.

        targ_time=REAL(sec)/3600+REAL(day*24)
        interval=fixhd_ancila(36,file_num)+fixhd_ancila(35,file_num)*12
        period=inthd_ancila(3,file_num)*interval
        IF (period /= 12) THEN
          icode=600+irec
          cmessage='REPLANCA: Non-standard period for periodic data'
          GO TO 9999
        END IF
        !  Difference between date now (month) &
        !    first date ancil file (month)
        months=i_month-fixhd_ancila(22,file_num)

        IF (use_lookup_dates_anc_time_interp) THEN
          ! correct code to use lookup header dates
          ! Correctly use day and hour from lookup header not fixhd_ancila 
          ! which contains values for first field on ancillary file only.
          step=months/interval
          i2=nlookup(irec)+lookup_step(irec)*step
          i1=i2+lookup_start_ancila(file_num,1)-1
          !  Check for time within month - using ppheader information
          IF ((i_day*24+i_hour) <  (lookup_ancila(3,i1)*24+                &
                                    lookup_ancila(4,i1))) THEN
            months=months-1
          END IF
          IF (months <  0) THEN
            months=months+12
          END IF
          ! recalculate STEP
          step=months/interval
          ! NB INTERVAL may be > 1 month
          months=step*interval
          iy=i_year
          im=MOD(fixhd_ancila(22,file_num)+months-1,12)+1
          IF (im >  i_month) iy=iy-1
          i2=nlookup(irec)+lookup_step(irec)*step
          i1=i2+lookup_start_ancila(file_num,1)-1
          ! DEPENDS ON: time2sec
          CALL time2sec(iy,im,lookup_ancila(3,i1),lookup_ancila(4,i1),         &
                        fixhd_ancila(25,file_num),fixhd_ancila(26,file_num),   &
                        ancil_ref_days,ancil_ref_secs,day,sec,lcal360)
          time1=REAL(sec)/3600+REAL(day*24)
          !  Calculate  TIME2 for second ancillary data time
          !  set IY correctly for time interpolation calculations
          iy=i_year
          im=MOD(fixhd_ancila(22,file_num)+months+interval-1,12)+1
          IF (im <  i_month) iy=iy+1
          i1=(im-1)/interval
          i2=nlookup(irec)+lookup_step(irec)*i1
          i1=i2+lookup_start_ancila(file_num,1)-1
          ! DEPENDS ON: time2sec
          CALL time2sec(iy,im,lookup_ancila(3,i1),lookup_ancila(4,i1),        &
                        fixhd_ancila(25,file_num),fixhd_ancila(26,file_num),  &
                        ancil_ref_days,ancil_ref_secs,day,sec,lcal360)
          time2=REAL(sec)/3600+REAL(day*24)

        ELSE   ! original code inaccurate use of FIXHD_ANCILA dates
          !  Check for time within month
          IF ((i_day*24+i_hour) < (fixhd_ancila(23,file_num)*24+               &
                                   fixhd_ancila(24,file_num))) THEN
            months=months-1
          END IF
          IF (months <  0) THEN
            months=months+12
          END IF

          step=months/interval
          ! NB INTERVAL may be > 1 month
          months=step*interval
          !  Calculate TIME1 for first ancillary data time
          !  set IY correctly for time interpolation calculations
          iy=i_year
          im=MOD(fixhd_ancila(22,file_num)+months-1,12)+1
          IF (im >  i_month) iy=iy-1
          ! DEPENDS ON: time2sec
          CALL time2sec(iy,im,fixhd_ancila(23,file_num),                     &
                        fixhd_ancila(24,file_num),                           &
                        fixhd_ancila(25,file_num),fixhd_ancila(26,file_num), &
                        ancil_ref_days,ancil_ref_secs,day,sec,lcal360)
          time1=REAL(sec)/3600+REAL(day*24)
          !  Calculate  TIME2 for second ancillary data time
          !    set IY correctly for time interpolation calculations
          iy=i_year
          im=MOD(fixhd_ancila(22,file_num)+months+interval-1,12)+1
          IF (im <  i_month) iy=iy+1
          ! DEPENDS ON: time2sec
          CALL time2sec(iy,im,fixhd_ancila(23,file_num),                     &
                        fixhd_ancila(24,file_num),                           &
                        fixhd_ancila(25,file_num),fixhd_ancila(26,file_num), &
                        ancil_ref_days,ancil_ref_secs,day,sec, lcal360)
          time2=REAL(sec)/3600+REAL(day*24)
        END IF  ! end use_lookup_dates_anc_time_interp test

        ! Do not interpolate in time if data time exactly matches model time

        IF (targ_time == time1) THEN
          linterpolate=.FALSE.
        END IF

      END IF  ! regular/non-regular


      !  Adjust YMD time back to start of updating interval

      i_year=i_year1
      i_month=i_month1
      i_day=i_day1
      i_hour=i_hour1


    END IF  ! non-periodic/periodic

    IF (linterpolate) THEN
      WRITE(umMessage,*)                                                  &
      ' REPLANCA - time interpolation for ancillary field, stashcode ',   &
                                          ancil_requests(irec)%stashcode
      CALL umPrint(umMessage,src='replanca')
      WRITE(umMessage,*)' targ_time,time1,time2 ',targ_time,time1,time2
      CALL umPrint(umMessage,src='replanca')
      WRITE(umMessage,*)' hours,int,period ',hours,interval,period
      CALL umPrint(umMessage,src='replanca')
    END IF

  END IF ! singletime/non-singletime

  !  2.3   Check STASH Code

  i2=nlookup(irec)+lookup_step(irec)*step

  i1=lookup_ancila(item_code, i2+lookup_start_ancila(file_num,1)-1)

  lmismatch=.FALSE.

  WRITE(umMessage,*)' Information used ',                                  &
    'in checking ancillary data set:',                                     &
    ' position of lookup table in dataset:',i2
  CALL umPrint(umMessage,src='replanca')

  WRITE(umMessage,*)' Position of first ',                                 &
    'lookup table referring to ',                                          &
    'data type ',nlookup(irec)
  CALL umPrint(umMessage,src='replanca')

  WRITE(umMessage,*)' Interval between lookup ',                           &
    'tables referring to data ',                                           &
    'type ', lookup_step(irec),' Number of steps', step
  CALL umPrint(umMessage,src='replanca')

  WRITE(umMessage,*)' STASH code in dataset ',i1,                          &
    '  STASH code requested ', ancil_requests(irec)%stashcode
  CALL umPrint(umMessage,src='replanca')

  WRITE(umMessage,*)'''start'' position of lookup tables for dataset ',    &
    file_num, 'in overall lookup array ' ,lookup_start_ancila(file_num,1)
  CALL umPrint(umMessage,src='replanca')

  IF (i1 /= ancil_requests(irec)%stashcode) THEN
    WRITE(umMessage,*) i1, ancil_requests(irec)%stashcode, irec
    CALL umPrint(umMessage,src='replanca')
    lmismatch=.TRUE.
  END IF

  !  Error exit if checks fail

  IF (lmismatch) THEN
    icode=200+irec
    cmessage='REPLANCA: PP HEADERS ON ANCILLARY FILE DO NOT MATCH'
    GO TO 9999
  END IF

  IF (linterpolate .AND. .NOT. single_time) THEN
    !  Check time interpolation factors
    IF (targ_time <  time1 .OR. targ_time >  time2) THEN
      WRITE(umMessage,*)' Information used in interpolation/replacement:'
      CALL umPrint(umMessage,src='replanca')
      WRITE(umMessage,*)' Time of first data=', time1
      CALL umPrint(umMessage,src='replanca')
      WRITE(umMessage,*)' Validity Time for update=', targ_time
      CALL umPrint(umMessage,src='replanca')
      WRITE(umMessage,*)' Time of second data=', time2
      CALL umPrint(umMessage,src='replanca')

      icode=500+irec
      cmessage='REPLANCA: TIME INTERPOLATION ERROR'
      GO TO 9999
    END IF
  END IF

  !  3   Loop over levels of ancillary data for field I
  !  Reset pointer for dataset

  !  Includes loop over X and Y components of surface currents

  lice_fraction= stash_item == stashcode_icefrac     ! ancil 27, stash 31
  lsnow_depth=   stash_item == stashcode_snow_amount ! ancil 9, stash 23
  lice_depth=    stash_item == stashcode_icethick    ! ancil 29, stash 32

  ! Find whether we are on theta level 0
  bot_level = -1
  lb_code = -1
  icode   = 0
  lb_code = exppxi(1, ancil_requests(irec)%stashcode/1000,             &
                   MOD(ancil_requests(irec)%stashcode,1000),9, icode,cmessage)

  IF (lb_code > 0) THEN
    ! DEPENDS ON: levcod
    CALL levcod(lb_code,bot_level,icode,cmessage)
  END IF

  IF (icode /= 0) THEN
    GO TO 9999
  END IF

  IF (bot_level /= 0) THEN
    bot_level = 1
  END IF

  DO level = 1, levels(irec)

    !  Do not go through loop for ice edge or snow edge
    IF ((lice_fraction .OR. lsnow_depth) .AND. level == 2) THEN
      CYCLE
    END IF

    !  3.1 Read data for single level of ancillary field.

    IF (.NOT. lice_fraction) THEN

      ! AMIPII case ice depth field not read from ancillary file
      IF (.NOT. (lice_depth .AND. l_amipii_ice_processing)) THEN

        CALL readflds(nftin,                                                   &
                      1,                                                       &
                      i2,                                                      &
                      lookup_ancila(:,lookup_start_ancila(file_num,1):),       &
                      ancil1,                                                  &
                      fixhd_ancila(:,file_num),                                &
                      1,                                                       &
                      icode,                                                   &
                      cmessage)

      END IF
      IF (icode /= 0) THEN
        icode=irec+100
        iounit=nftin
        cmessage='REPLANCA :I/O ERROR '
        GO TO 9999
      END IF

    ELSE

      ! If ice-fraction,read fractional time field as well
      !   UNLESS IT IS A SINGLE TIME FIELD
      ! If snow-depth,read fractional time field as well only if time
      !   interpolation required.

      IF (.NOT. single_time .AND. .NOT. l_amipii_ice_processing) THEN
        IF (lookup_ancila(item_code,i2+lookup_start_ancila(file_num,1))    &
                                    == stashcode_ice_edge_ancil) THEN

          CALL readflds(nftin,                                                 &
                        2,                                                     &
                        i2,                                                    &
                        lookup_ancila(:,lookup_start_ancila(file_num,1):),     &
                        ice_extent,                                            &
                        fixhd_ancila(:,file_num),                              &
                        1,                                                     &
                        icode,                                                 &
                        cmessage)

          IF (icode /= 0) THEN
            icode=irec+100
            iounit=nftin
            cmessage='REPLANCA :I/O ERROR '
            GO TO 9999
          END IF

        ELSE
          icode=irec+100
          iounit=nftin
          cmessage='REPLANCA :ICE CHANGE DATA MISSING'
          GO TO 9999
        END IF
      ELSE  ! single time or l_amipii_ice_processing - ie no time change field

        CALL readflds(nftin,                                                   &
                      1,                                                       &
                      i2,                                                      &
                      lookup_ancila(:,lookup_start_ancila(file_num,1):),       &
                      ice_extent,                                              &
                      fixhd_ancila(:,file_num),                                &
                      1,                                                       &
                      icode,                                                   &
                      cmessage)

        IF (icode /= 0) THEN
          icode=irec+100
          iounit=nftin
          cmessage='REPLANCA: I/O ERROR'
          GO TO 9999
        END IF
      END IF
    END IF

    IF (lsnow_depth .AND. linterpolate) THEN
      IF (lookup_ancila(item_code,i2+lookup_start_ancila(file_num,1))     &
                                  == stashcode_snow_edge) THEN

        CALL readflds(nftin,                                                   &
                      1,                                                       &
                      i2+1,                                                    &
                      lookup_ancila(:,lookup_start_ancila(file_num,1):),       &
                      snow_change,                                             &
                      fixhd_ancila(:,file_num),                                &
                      1,                                                       &
                      icode,                                                   &
                      cmessage)

        IF (icode /= 0) THEN
          icode=irec+100
          iounit=nftin
          cmessage='REPLANCA :I/O ERROR '
          GO TO 9999
        END IF

      ELSE
        icode=irec+100
        iounit=nftin
        cmessage='REPLANCA :SNOW CHANGE DATA MISSING'
        GO TO 9999
      END IF
    END IF

    ! If sea surface temperature or other ice fields, read ice fraction
    !  and fractional time field if not already present and if required
    !  by time interpolation.  Similar if ice depth needed.

    !  ancref=28  stash=24   Sea-surface Temperature
    !  ancref=29  stash=32   Sea-ice Thickness
    !  ancref=27  stash=31   Sea-ice Fraction
    !             stash=38   Ice edge ancillary

    IF ( ancil_requests(irec)%stashcode == stashcode_icethick .OR.     &
        (ancil_requests(irec)%stashcode == stashcode_surftemp          &
         .AND. lt_int_c)) THEN

      ! ancil_req_num_icefrac = address of sea ice frac in ancil request array
      ! ftn_icefrac = fortran unit for sea ice fraction
      !   both previously determined

      IF (.NOT. update(ancil_req_num_icefrac)) THEN
        i3 = nlookup(ancil_req_num_icefrac)                &
              + lookup_step(ancil_req_num_icefrac)*step    &
              + lookup_start_ancila(file_num_icefrac,1)

        IF ( lookup_ancila(item_code,i3) == stashcode_ice_edge_ancil) THEN

          CALL readflds(ftn_icefrac,                                           &
                        2,                                                     &
                        nlookup(ancil_req_num_icefrac)                         &
                         + lookup_step(ancil_req_num_icefrac)*step,            &
                        lookup_ancila(:,lookup_start_ancila(file_num_icefrac,  &
                                                            1):),              &
                        ice_extent,                                            &
                        fixhd_ancila(:,file_num_icefrac),                      &
                        1,                                                     &
                        icode,                                                 &
                        cmessage)

          IF (icode /= 0) THEN
            icode=irec+100
            iounit=nftin
            cmessage='REPLANCA :I/O ERROR '
            GO TO 9999
          END IF
          IF ( TRANSFER(lookup_ancila(bmdi,i3-1),                          &
                        dummy_real)  /=  rmdi ) THEN
            icode = 700 + irec
            WRITE(cmessage,*)&
            'REPLANCA: RMDI in lookup of ancillary file of times',         &
                    ' of sea-ice change not standard'
            GO TO 9999
          END IF

        ELSE
          icode=irec+100
          iounit=nftin
          cmessage='REPLANCA :ICE FIELD DATA MISSING'
          GO TO 9999
        END IF
      END IF
    END IF

    !  3.3 If time interpolation required, read second record

    IF (linterpolate) THEN

      i1=i2+ lookup_step(irec)
      IF (i1 <= fixhd_ancila(152,file_num)) THEN

        ! If the two fields differ we have a problem with the ancillary.
        IF ( lookup_ancila(item_code,lookup_start_ancila(file_num,1)+i1-1) /=  &
             lookup_ancila(item_code,                                          &
                           lookup_start_ancila(file_num,1)+i2-1) ) THEN
          icode=irec+100
          cmessage='REPLANCA: start and end fields are different.'
          GO TO 9999
        END IF

        ! AMIP II and ice depth don't read in ice depth field
        IF (.NOT. (l_amipii_ice_processing .AND. lice_depth)) THEN

          CALL readflds(nftin,                                                 &
                        1,                                                     &
                        i1,                                                    &
                        lookup_ancila(:,lookup_start_ancila(file_num,1):),     &
                        ancil2,                                                &
                        fixhd_ancila(:,file_num),                              &
                        1,                                                     &
                        icode,                                                 &
                        cmessage)
        END IF
        IF (icode /= 0) THEN
          icode=irec+300
          iounit=nftin
          cmessage='REPLANCA :I/O ERROR '
          GO TO 9999
        END IF

      ELSE !end of data on file

        !   If end of data has been reached go back to the start, if data is
        !   periodic. Otherwise cancel time interpolation

        IF (periodic) THEN

          i1 = nlookup(irec) + level - 1
          ! If the two fields differ we have a problem with the ancillary.
          IF ( lookup_ancila(item_code,                                    &
                             lookup_start_ancila(file_num,1)+i1-1) /=      &
               lookup_ancila(item_code,                                    &
                             lookup_start_ancila(file_num,1)+i2-1) ) THEN
            icode=irec+100
            cmessage='REPLANCA: start and end fields are different.'
            GO TO 9999
          END IF

          CALL readflds(nftin,                                                 &
                        1,                                                     &
                        i1,                                                    &
                        lookup_ancila(:,lookup_start_ancila(file_num,1):),     &
                        ancil2,                                                &
                        fixhd_ancila(:,file_num),                              &
                        1,                                                     &
                        icode,                                                 &
                        cmessage)

          IF (icode /= 0) THEN
            icode=irec+300
            iounit=nftin
            cmessage='REPLANCA :I/O ERROR '
            GO TO 9999
          END IF
        ELSE
          ! We switch off time interpolation if we have reached
          ! the end of the file.
          CALL umPrint(                                                    &
            'REPLANCA: Reached end of ancillary, '//                       &
            'switched off time interpolation.',src='replanca-replanca1a')
          linterpolate=.FALSE.
        END IF
      END IF! End of position on file test

      icode=0
    END IF ! End LINTERPOLATE

    !  3.4 Perform time interpolation

    IF (linterpolate) THEN

      zero=0.0

      !  Select appropriate time interpolation for each field
      !  Snowdepth: set equal to zero if no snow cover

      IF (lsnow_depth) THEN
        DO i=1,p_field
          pres_value(i)=zero
        END DO

        ! For the call to T_INT_C, need to know BMDI is OK for SNOW_CHANGE
        !  which was read in from position I2+1.
        IF ( TRANSFER(lookup_ancila(bmdi,lookup_start_ancila(file_num,1)+i2),  &
                                    dummy_real)  /=  rmdi ) THEN
          icode = 700 + irec
          WRITE(cmessage,*)                                                &
           'REPLANCA: RMDI in lookup of ancil file of times',              &
                    ' of snow change non-standard '
          GO TO 9999
        END IF

        CALL t_int_c (ancil1,time1,ancil2,time2,ancil_data,                &
             targ_time,p_field,snow_change,ancil1,pres_value)

        ! Ice fraction: ice depth set equal to zero if no ice
      ELSE IF (ancil_requests(irec)%stashcode == stashcode_icefrac .OR.    &
               ancil_requests(irec)%stashcode == stashcode_icethick) THEN

        IF (ancil_requests(irec)%stashcode == stashcode_icefrac) THEN
          ! For the call to T_INT_C, need to know BMDI
          ! is OK for ICE_EXTENT(1,2) which was read in from position I1+1
          IF (.NOT. l_amipii_ice_processing) THEN
            IF ( TRANSFER(lookup_ancila(bmdi,                                 &
                                        lookup_start_ancila(file_num,1)+i1),  &
                                        dummy_real)  /=  rmdi ) THEN
              icode = 700 + irec
              WRITE(cmessage,*)                                            &
               'REPLANCA: RMDI in lookup of ancil file of times',          &
                     ' of sea-ice chge non-standard'
              GO TO 9999
            END IF
          END IF

          IF (l_amipii_ice_processing) THEN
            ! linear uncontrolled time interpolation
            CALL t_int (ice_extent,time1,ancil2,time2,ancil_data,          &
                        targ_time,p_field)

            ! For AMIP II strictly ice concentrations should range between
            ! 0.0 and 1.0 but because of assumptions on T* in the boundary
            ! layer and radiation schemes ice conc's are restricted to
            ! 0.3 to 1.0. This will allow SSTs in areas of less than 30%
            ! ice to be used rather than TFS=-1.8C.

            DO i=1,p_field
              IF (ancil_data(i) <  0.3) ancil_data(i)=0.0
              IF (ancil_data(i) >  1.0) ancil_data(i)=1.0
            END DO

          ELSE       ! non AMIPII_ice_processing option
            DO i=1,p_field
              pres_value(i)=0
            END DO

            CALL t_int_c (ice_extent,time1,ancil2,time2,ancil_data,        &
                          targ_time,p_field,ice_extent(1,2),ice_extent,    &
                          pres_value)

          END IF     ! end AMIPII_ice_processing test

        ELSE IF (ancil_requests(irec)%stashcode == stashcode_icethick) THEN

          DO i=1,p_field
            pres_value(i)=0
          END DO

          CALL t_int_c (ancil1,time1,ancil2,time2,ancil_data,              &
                        targ_time,p_field,ice_extent(1,2),ice_extent,      &
                        pres_value)

        END IF

        ! Sea surface temperature, set equal to TFS if ice present

      ELSE IF (ancil_requests(irec)%stashcode == stashcode_surftemp &
                                                  .AND. lt_int_c) THEN
        IF (l_amipii_ice_processing) THEN

          CALL t_int (ancil1,time1,ancil2,time2,                          &
                                             ancil_data,targ_time,p_field)
          ! remove any T below TFS
          DO i=1,p_field
            IF (ancil_data(i) <  tfs)  ancil_data(i)=tfs
          END DO

        ELSE     ! non AMIPII_ice_processing option

          IF (.NOT. l_leads_temp_prog) THEN
            DO i=1,p_field
              pres_value(i)=tfs

              ! Set no_ice_extent indicator for controlled SST interpolation
              IF (ice_extent(i,1) == 0) THEN
                no_ice_extent(i)=1.0
              ELSE
                no_ice_extent(i)=0.0
              END IF
            END DO

            CALL t_int_c (ancil1,time1,ancil2,time2,ancil_data,            &
                          targ_time,p_field,ice_extent(1,2),no_ice_extent, &
                          pres_value)
          ELSE
            CALL t_int (ancil1,time1,ancil2,time2,ancil_data,              &
                                                         targ_time,p_field)
          END IF

        END IF   ! end AMIPII_ice_processing test
        ! Otherwise linear interpolation in time, unless mdi
        !   present at either time.

      ELSE

        ! Time interpolation checks the data against the standard mdi
        !   - check that the field is labelled as using the same one.
        !  (It is to ensure the correct I1 here that I3 is used above.)

        IF ( TRANSFER(lookup_ancila(bmdi,lookup_start_ancila(file_num,1)+i1-1),&
                                    dummy_real)  /=  rmdi .OR.             &
             TRANSFER(lookup_ancila(bmdi,lookup_start_ancila(file_num,1)+i2-1),&
                                    dummy_real)  /=  rmdi ) THEN
          WRITE(umMessage, *) 'LOOKUPS:',                                  &
          TRANSFER(lookup_ancila(bmdi,lookup_start_ancila(file_num,1)+i1-1),   &
                                    dummy_real),                           &
          TRANSFER(lookup_ancila(bmdi,lookup_start_ancila(file_num,1)+i2-1),   &
                                    dummy_real)
          CALL umPrint(umMessage,src='replanca')
          icode = 700 + irec
          cmessage = 'REPLANCA: MDI in lookup of ancil file is non-standard'
          GO TO 9999
        END IF

        len1=p_field
        !   Ozone, test for zonal mean or full field
        IF (ancil_requests(irec)%stashcode == stashcode_ozone) THEN
          IF (lookup_ancila(lbnpt,                                         &
                            lookup_start_ancila(file_num,1)+i2-1) == 1) THEN
            len1=p_rows
          END IF
          !   Cariolle ozone, test for zonal mean or full field.
          !   Currently same test as for conventional ozone.
        ELSE IF (ancil_requests(irec)%stashcode >= stashcode_ozone_tracer .AND.&
                 ancil_requests(irec)%stashcode <= stashcode_o3_colo3) THEN
          IF (lookup_ancila(lbnpt,                                         &
                            lookup_start_ancila(file_num,1)+i2-1) == 1) THEN
            len1=p_rows
          END IF
        END IF

        CALL t_int(ancil1,time1,ancil2,time2,ancil_data,targ_time,len1)

      END IF ! End Lsnow_depth

      ! If no interpolation, copy data into final array

    ELSE ! no interpolation
      IF (lice_fraction) THEN
        IF (l_amipii_ice_processing) THEN
          DO i=1,p_field

            ancil_data(i)=ice_extent(i,1)

            ! For AMIP II strictly ice concentrations should range between
            ! 0.0 and 1.0 but because of assumptions on T* in the boundary
            ! layer and radiation schemes ice conc's are restricted to
            ! 0.3 to 1.0. This will allow SSTs in areas of less than 30%
            ! ice to be used rather than TFS=-1.8C.

            IF (ancil_data(i) <  0.3) ancil_data(i)=0.0
            IF (ancil_data(i) >  1.0) ancil_data(i)=1.0

          END DO
        ELSE           ! non AMIP II option
          DO i=1,p_field
            ancil_data(i)=ice_extent(i,1)
          END DO
        END IF           ! end of AMIPII_ice_processing test
      ELSE IF (l_amipii_ice_processing .AND. &
               ancil_requests(irec)%stashcode == stashcode_surftemp) THEN
        DO i=1,p_field
          ancil_data(i)=ancil1(i)
          IF (ancil_data(i) <  tfs) ancil_data(i)=tfs
        END DO
      ELSE
        DO i=1,p_field
          ancil_data(i)=ancil1(i)
        END DO
      END IF
    END IF !End interpolate/no interpolate

    !  3.5 Updating action for each field at each level
    !      Fields replaced except that Sea Surface Temperature may be
    !      incremented. Take apropriate action for each field.

    IF (ancil_requests(irec)%stashcode == stashcode_lsm .OR.                &
       ancil_requests(irec)%stashcode == stashcode_orog .OR.                &
       ancil_requests(irec)%stashcode == stashcode_ozone  .OR.              &
       ancil_requests(irec)%stashcode == stashcode_SO2_emiss .OR.           &
       ancil_requests(irec)%stashcode == stashcode_dimethyl_sul_emiss .OR.  &
       ancil_requests(irec)%stashcode == stashcode_total_aero_emiss .OR.    &
       ancil_requests(irec)%stashcode == stashcode_total_aero .OR.          &
      (ancil_requests(irec)%stashcode >= stashcode_user_anc_sing1 .AND.     &
       ancil_requests(irec)%stashcode <= stashcode_user_anc_sing20 .AND.    &
       l_ukca) .OR.                                                         &
      (ancil_requests(irec)%stashcode >= stashcode_ammonia_gas_emiss .AND.  &
       ancil_requests(irec)%stashcode <= stashcode_soot_hi_lev) .OR.        &
      (ancil_requests(irec)%stashcode >= stashcode_3d_nat_so2_em .AND.      &
       ancil_requests(irec)%stashcode <= stashcode_hi_SO2_emiss_emiss) .OR. &
       ancil_requests(irec)%stashcode == stashcode_CO2_surf_emiss .OR.      &
      (ancil_requests(irec)%stashcode >= stashcode_user_anc_mult1 .AND.     &
       ancil_requests(irec)%stashcode <= stashcode_user_anc_mult20) .OR.    &
      (ancil_requests(irec)%stashcode >= stashcode_dust_parent_clay .AND.   &
       ancil_requests(irec)%stashcode <= stashcode_dust_soil_mf6) .OR.      &
       ancil_requests(irec)%stashcode == stashcode_biom_surf_em  .OR.       &
       ancil_requests(irec)%stashcode == stashcode_biom_elev_em  .OR.       &
       ancil_requests(irec)%stashcode == stashcode_biom_elev_em_h1 .OR.     &
       ancil_requests(irec)%stashcode == stashcode_biom_elev_em_h2 .OR.     &
       ancil_requests(irec)%stashcode == stashcode_dms_conc_sea  .OR.       &
       ancil_requests(irec)%stashcode == stashcode_3d_ho2_conc   .OR.       &
      (ancil_requests(irec)%stashcode >= stashcode_clim_biogenic_aero .AND. &
       ancil_requests(irec)%stashcode <= stashcode_clim_delta_aero) .OR.    &
      (ancil_requests(irec)%stashcode >= stashcode_ozone_tracer .AND.       &
       ancil_requests(irec)%stashcode <= stashcode_o3_colo3) .OR.           &
       ancil_requests(irec)%stashcode == stashcode_ocff_surf_emiss .OR.     &
       ancil_requests(irec)%stashcode == stashcode_ocff_hilev_emiss .OR.    &
       ancil_requests(irec)%stashcode == stashcode_unfilt_orog .OR.         &
      (ancil_requests(irec)%stashcode >= stashcode_urbhgt .AND.             &
       ancil_requests(irec)%stashcode <= stashcode_urbwrr) ) THEN

      !  3.5.0 Updates at all points

      len1=p_field
      !   Ozone, test for zonal mean or full field
      IF (ancil_requests(irec)%stashcode == stashcode_ozone) THEN
        IF (lookup_ancila(lbnpt,lookup_start_ancila(file_num,1)+i2-1) == 1) THEN
          len1=p_rows
        END IF

        !   Cariolle ozone, test for zonal mean or full field.
        !   Currently same test as for conventional ozone.
      ELSE IF (ancil_requests(irec)%stashcode >= stashcode_ozone_tracer .AND.  &
               ancil_requests(irec)%stashcode <= stashcode_o3_colo3) THEN
        IF (lookup_ancila(lbnpt,lookup_start_ancila(file_num,1)+i2-1) == 1) THEN
          len1=p_rows
        END IF
      END IF

      DO i=1,len1
        d1(d1_anciladd(irec)+i-1+(level-bot_level)*len1)=ancil_data(i)
      END DO

      ! Copy level 1 to level 0
      IF (bot_level == 0 .AND. level == 1) THEN
        DO i=1,len1
          d1(d1_anciladd(irec)+i-1)=ancil_data(i)
        END DO
      END IF

      !  3.5.1 Updates over all land points

    ELSE IF ((ancil_requests(irec)%stashcode >= stashcode_orog_var .AND.   & ! (34)
         ancil_requests(irec)%stashcode <= stashcode_orog_gdyy)            & ! (37)
       .OR.  ancil_requests(irec)%stashcode == stashcode_soil_temp         & ! (20)
       .OR.  ancil_requests(irec)%stashcode == stashcode_snow_amount       & ! (23)
       .OR.  ancil_requests(irec)%stashcode == stashcode_vol_smc_wilt      & ! (40)
       .OR.  ancil_requests(irec)%stashcode == stashcode_vol_smc_cri       & ! (41)
       .OR.  ancil_requests(irec)%stashcode == stashcode_vol_smc_sat       & ! (43)
       .OR.  ancil_requests(irec)%stashcode == stashcode_Ksat              & ! (44)
       .OR.  ancil_requests(irec)%stashcode == stashcode_thermal_capacity  & ! (46)
       .OR.  ancil_requests(irec)%stashcode == stashcode_thermal_conduct   & ! (47)
       .OR.  ancil_requests(irec)%stashcode == stashcode_rough_length      & ! (26)
       .OR.  ancil_requests(irec)%stashcode == stashcode_veg_frac          & ! (50)
       .OR.  ancil_requests(irec)%stashcode == stashcode_runoff_coast_out  & ! (93)
       .OR.  ancil_requests(irec)%stashcode == stashcode_soil_suction      & ! (48)
       .OR.  ancil_requests(irec)%stashcode == stashcode_soil_moist        & !  (9)
       .OR. (ancil_requests(irec)%stashcode >= stashcode_user_anc_sing1 .AND.  &
             ancil_requests(irec)%stashcode <= stashcode_user_anc_sing20 .AND. &
             .NOT. l_ukca)                                     & !(301-320)
       .OR.  ancil_requests(irec)%stashcode == stashcode_sil_orog_rough    & ! (17)
       .OR.  ancil_requests(irec)%stashcode == stashcode_hlf_pk_to_trf_ht  & ! (18)
       .OR.  ancil_requests(irec)%stashcode == stashcode_orog_x_grad       & !  (5)
       .OR.  ancil_requests(irec)%stashcode == stashcode_orog_y_grad       & !  (6)
       .OR.  ancil_requests(irec)%stashcode == stashcode_clapp_hb          & !(207)
       .OR.  ancil_requests(irec)%stashcode == stashcode_can_conduct       & !(213)
       .OR.  ancil_requests(irec)%stashcode == stashcode_frac_surf_type    & !(216)
       .OR.  ancil_requests(irec)%stashcode == stashcode_lai               & !(217)
       .OR.  ancil_requests(irec)%stashcode == stashcode_canopy_height     & !(218)
       .OR.  ancil_requests(irec)%stashcode == stashcode_disturb_frac_veg  & !(219)
       .OR.  ancil_requests(irec)%stashcode == stashcode_disturb_frac_veg_prev & !(286)
       .OR.  ancil_requests(irec)%stashcode == stashcode_wood_prod_fast     & !(287)
       .OR.  ancil_requests(irec)%stashcode == stashcode_wood_prod_med      & !(288)
       .OR.  ancil_requests(irec)%stashcode == stashcode_wood_prod_slow     & !(289)
       .OR.  ancil_requests(irec)%stashcode == stashcode_nitrogen_deposition & !(447)
       .OR.  ancil_requests(irec)%stashcode == stashcode_crop_frac         & !(448)
       .OR.  ancil_requests(irec)%stashcode == stashcode_pasture_frac      & !(458)
       .OR.  ancil_requests(irec)%stashcode == stashcode_snw_free_alb_bs   & !(220)
       .OR.  ancil_requests(irec)%stashcode == stashcode_soil_carbon_cont  & !(223)
       .OR.  ancil_requests(irec)%stashcode == stashcode_top_ind_mean      & !(274)
       .OR.  ancil_requests(irec)%stashcode == stashcode_top_ind_stddev    & !(275)
       .OR.  ancil_requests(irec)%stashcode == stashcode_rgrain            & !(231)
       .OR.  ancil_requests(irec)%stashcode == stashcode_snow_tile         & !(240)
       .OR.  ancil_requests(irec)%stashcode == stashcode_snow_grnd         & !(242)
       .OR.  ancil_requests(irec)%stashcode == stashcode_snowdep_grd_tile  & !(376)
       .OR.  ancil_requests(irec)%stashcode == stashcode_snowpack_bk_dens  & !(377)
       .OR.  ancil_requests(irec)%stashcode == stashcode_nsnow_layrs_tiles & !(380)
       .OR.  ancil_requests(irec)%stashcode == stashcode_snow_laythk_tiles & !(381)
       .OR.  ancil_requests(irec)%stashcode == stashcode_snow_ice_tile     & !(382)
       .OR.  ancil_requests(irec)%stashcode == stashcode_snow_liq_tile     & !(383)
       .OR.  ancil_requests(irec)%stashcode == stashcode_snow_T_tile       & !(384)
       .OR.  ancil_requests(irec)%stashcode == stashcode_snow_laydns_tiles & !(385)
       .OR.  ancil_requests(irec)%stashcode == stashcode_snow_grnsiz_tiles & !(386)
       .OR.  ancil_requests(irec)%stashcode == stashcode_surf_sw_alb       & !(243)
       .OR.  ancil_requests(irec)%stashcode == stashcode_surf_vis_alb      & !(244)
       .OR.  ancil_requests(irec)%stashcode == stashcode_surf_nir_alb      & !(245)
       .OR.  ancil_requests(irec)%stashcode == stashcode_z0m_soil          & ! (97)
       .OR.  ancil_requests(irec)%stashcode == stashcode_flake_depth) THEN   !(291)

      !  If not reconfiguration, set snowdepth values at all land points
      !  Reset TSTAR to TM if snow cover present

      IF (lsnow_depth) THEN
        DO i=1,p_field
          IF (land(i)) THEN
            d1(d1_anciladd(irec)+i-1)=ancil_data(i)
            IF (tstar_land(i) >  tm .AND. ancil_data(i) >  0.0) THEN
              tstar_land(i)=tm
            END IF
          END IF
        END DO

        !  Set all other fields , which are stored at land points only

      ELSE
        !  If field is a single time soil moisture, only update if
        !  model time matches data time and then deactivate to prevent
        !  any further updating

        IF (ancil_requests(irec)%stashcode == stashcode_soil_moist .AND.     &
                fixhd_ancila(10,file_num) == 0 ) THEN

          IF (lookup_ancila(lbyr,                                            &
                            lookup_start_ancila(file_num,1)) == i_year .AND. &
             lookup_ancila(lbmon,                                            &
                           lookup_start_ancila(file_num,1)) == i_month .AND. &
             lookup_ancila(lbdat,                                            &
                           lookup_start_ancila(file_num,1)) == i_day .AND.   &
             lookup_ancila(lbhr,                                             &
                           lookup_start_ancila(file_num,1)) == i_hour) THEN

            WRITE(umMessage,*) 'Updating soil moisture at ',                 &
                   i_year,i_month,i_day,i_hour,                              &
            ' for level ',TRANSFER(                                          &
             lookup_ancila(blev,lookup_start_ancila(file_num,1)+level-1),    &
                                         dummy_real)
            CALL umPrint(umMessage,src='replanca')
            dz_soil(level)=TRANSFER(                                         &
              lookup_ancila(blev,lookup_start_ancila(file_num,1)+level-1),   &
                                         dummy_real)

            CALL compress_to_mask(ancil_data,d1(d1_anciladd(irec)+           &
                                   (level-1)*land_field),land,p_field,i)
            ! Switch off to prevent further attempts to update
            ancil_requests(irec)%period=0
            steps(irec)=0
            ! Set flag to indicate that soil moisture has been updated
            smc_updated=.TRUE.
          ELSE
            WRITE(umMessage,*) 'Update of soil moisture skipped'
            CALL umPrint(umMessage,src='replanca')
          END IF

        ELSE
          ! other fields
          CALL compress_to_mask(ancil_data,d1(d1_anciladd(irec)+  &
                               (level-1)*land_field),land,p_field,i)
        END IF

      END IF

      ! Iceberg calving for the OASIS coupler:
    ELSE IF (ancil_requests(irec)%stashcode == stashcode_iceberg_calving) THEN
      DO i=1,u_field
        d1(d1_anciladd(irec)+i-1)=ancil_data(i)
      END DO

      !  3.5.2 Ice fraction
    ELSE IF (ancil_requests(irec)%stashcode == stashcode_icefrac) THEN
      DO i=1,p_field
        ice_fraction(i)=0.0
        IF (sea(i)) THEN
          ice_fraction(i)=ancil_data(i)
        END IF
      END DO

      ! Reduce TSTAR to TFS where ice fraction greater than zero
      ! Required at present because radiation and boundary layer codes
      ! assume T* is TFS and ignore any value set in TSTAR.

      IF (.NOT. l_leads_temp_prog) THEN
        DO i=1,p_field
          IF (ice_fraction(i) >  0.0) THEN
            tstar_ssi(i)=MIN(tstar_ssi(i),tfs)
          END IF
        END DO
      END IF

      !  3.5.3 Sea surface temperatures for atmosphere, adds anomaly field
      !        (if defined) to updated ancillary where no sea-ice present,
      !        otherwise updates field directly from ancillary

    ELSE IF (ancil_requests(irec)%stashcode == stashcode_surftemp) THEN

      DO i=1,p_field
        IF (l_ctile .OR. ice_fraction(i) <= 0.0) THEN
          IF (sea(i)) THEN
            IF (l_sstanom) THEN
              tstar_sea(i)=ancil_data(i)+tstar_anom(i)
            ELSE
              tstar_sea(i)=ancil_data(i)
            END IF
            IF (ice_fraction(i) <= 0.0) tstar_ssi(i)=tstar_sea(i)
          END IF
        END IF
      END DO

      ! Sea-point fields, where no special action is required:
    ELSE IF (ancil_requests(irec)%stashcode == stashcode_icethick) THEN
      ! Sea-ice Thickness
      DO i=1,p_field
        IF (sea(i)) THEN
          d1(d1_anciladd(irec)+i-1)=ancil_data(i)
        END IF
      END DO

      ! Ocean Nr. Surface Chlorophyll content
    ELSE IF (ancil_requests(irec)%stashcode == &
         stashcode_ocnsrf_chlorophyll) THEN
      DO i=1,p_field
        d1(d1_anciladd(irec)+i-1)=ancil_data(i)
      END DO

      !  3.5.5 Surface currents

    ELSE IF (ancil_requests(irec)%stashcode == stashcode_surf_z_curr) THEN
      DO i=1,u_field
        d1(d1_anciladd(irec)+i-1)=ancil_data(i)
      END DO

    ELSE IF (ancil_requests(irec)%stashcode == stashcode_surf_m_curr) THEN
      DO i=1,v_field
        d1(d1_anciladd(irec)+i-1)=ancil_data(i)
      END DO

    ELSE IF (ancil_requests(irec)%stashcode >= stashcode_riv_sequence .AND. &
             ancil_requests(irec)%stashcode <= stashcode_riv_storage) THEN
      ! River Routing - not yet supported in UM
      icode=750
      cmessage='REPLANCA: ERROR Trying to use Riv Route Ancils'
      GO TO 9999

    ELSE
      ! if we get to here, then an update request has been put in but the 
      ! code above does not know how to deal with it. This should be an
      ! immediate fatal error:
      WRITE (cmessage,'(A,I6,A)')                                  &
           ' REPLANCA: ERROR - stash code ',               &
           ancil_requests(irec)%stashcode,' omitted in field,stashcode tests'
      icode = 30
      CALL ereport( routinename, icode, cmessage )

    END IF !End tests on stashcodes

    i2=i2+1

  END DO !  End loop over levels
END DO ! End main loop over ancil fields

IF (l_planet_grey_surface) THEN
  ! Do not adjust tstar if using an idealised grey surface.
ELSE IF (l_ctile) THEN
  DO i=1,p_field
    IF (sea(i) .AND. ice_fraction(i) >  0.0) THEN
      IF (l_leads_temp_prog .OR. l_amipii_ice_processing) THEN

        tstar_ssi(i)=ice_fraction(i)*tstar_sice(i)                             &
          +(1.0-ice_fraction(i))*tstar_sea(i)

      ELSE

        tstar_sea(i)=tfs
        tstar_sice(i)=(tstar_ssi(i)                                            &
          -(1.0-ice_fraction(i))*tstar_sea(i))/ice_fraction(i)

      END IF
    END IF

    tstar(i)=fland_g(i)*tstar_land(i)                                          &
      +(1.0-fland_g(i))*tstar_ssi(i)
  END DO
ELSE
  DO i=1,p_field
    IF (land(i)) THEN
      tstar(i)=tstar_land(i)
    ELSE
      tstar(i)=tstar_ssi(i)
    END IF
  END DO
END IF

!     Set up surface temperatures:
IF (l_ctile) THEN
  DO i=1,p_field
    tstar_land_ctile(i)=tstar_land(i)
    tstar_sea_ctile(i)=tstar_sea(i)
    ! The use of TSTAR_SICE appears to cause problems in
    ! some configurations (e.g. seasonal). Possibly because
    ! of the use of inconsistent ancillary fields.
    ! Hence this is commented out but retained for reference.
    ! [See also equivalent change in replanca-rcf_replanca.F90]
    !tstar_sice_ctile(i)=tstar_sice(i)
  END DO
END IF

9999 CONTINUE

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE replanca

END MODULE replanca_mod
