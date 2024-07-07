! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Read in Ancillary Fields
!
!  Subroutine replanca_rcf_replanca  - Read in Ancillary Fields
!
! Description:
!   Read in Ancillary Fields
!
! Method:
!   This routine reads in the required ancillary fields. Any time
!   interpolation is done as required. Both 360 and 365 day calendars
!   are catered for. The ANCILmaster files are used in the
!   reconfiguration.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!

SUBROUTINE replanca_rcf_replanca(i_year,i_month,i_day,i_hour,    &
                    i_minute,i_second,                           &
                    p_field,p_rows,u_field,v_field,rr_field,     &
                    land_field,                                  &
                    d1,land,                                     &
                    ice_fraction,tstar,fland,                    &
                    tstar_land_ctile,tstar_sea_ctile,            &
                    tstar_sice_ctile,                            &
                    tstar_anom,                                  &
                    ns_space,first_lat,                          &
                    len1_lookup,len_fixhd,len_inthd,             &
                    len_d1,fixhd,inthd,                          &
                    lookup,rlookup,lookup_start,                 &
                    nlookups,                                    &
                    icode,cmessage)

USE ancil_mod, ONLY:                                             &
    AncF_UnitNo,                                                 &
    ancil_add,                                                   &
    ancil_files,                                                 &
    ancil_requests,                                              &
    Levels,                                                      &
    Lookup_Step,                                                 &
    NLookup,                                                     &
    num_ancil_files,                                             &
    num_ancil_requests,                                          &
    find_ancil_req_by_stash,                                     &
    find_ancil_file_by_stash

USE decomp_params, ONLY:                                          &
    Decomp_rcf_output

USE gen_phys_inputs_mod, ONLY:  l_leads_temp_prog
USE nlstcall_mod, ONLY: lcal360

USE ancilcta_namelist_mod, ONLY:                                 &
    l_amipii_ice_processing,                                      &
    use_lookup_dates_anc_time_interp,                                    &
    l_sstanom

USE UM_ParCore, ONLY:                                            &
    mype

USE umPrintMgr, ONLY:                                            &
    umPrint,                        &
    umMessage,                      &
    PrintStatus, PrStatus_Diag

USE ukca_option_mod, ONLY:                                       &
    l_ukca

USE jules_sea_seaice_mod, ONLY: l_ctile

USE jules_snow_mod, ONLY: nsmax

USE io

USE water_constants_mod, ONLY: tfs,tm

USE lookup_addresses

USE Rcf_Exppx_Mod, ONLY:      &
    Rcf_Exppx

USE Rcf_Ppx_Info_Mod, ONLY:    &
    STM_Record_Type

USE Rcf_Grid_Type_Mod, ONLY: &
    Output_Grid

USE Rcf_Level_Code_Mod, ONLY: &
    Rcf_Level_Code


USE t_int_mod, ONLY: t_int
USE t_int_c_mod, ONLY: t_int_c

USE um_stashcode_mod, ONLY:                                    &
    stashcode_land_frac, stashcode_tstar, stashcode_icefrac,      &
    stashcode_icethick, stashcode_soil_temp, stashcode_vol_smc_wilt, &
    stashcode_vol_smc_cri, stashcode_vol_smc_sat, stashcode_Ksat, &
    stashcode_thermal_capacity,stashcode_thermal_conduct,         &
    stashcode_mean_snow, stashcode_riv_sequence,                  &
    stashcode_riv_direction, stashcode_riv_storage,               &
    stashcode_ozone,stashcode_ozone_tracer,stashcode_o3_colo3,    &
    stashcode_lsm, stashcode_orog,                                &
    stashcode_SO2_emiss,stashcode_dimethyl_sul_emiss,             &
    stashcode_total_aero_emiss, stashcode_total_aero,             &
    stashcode_soot_hi_lev,stashcode_ammonia_gas_emiss,            &
    stashcode_3d_nat_so2_em,stashcode_hi_SO2_emiss_emiss,         &
    stashcode_CO2_surf_emiss, stashcode_soil_massfrac6,           &
    stashcode_dust_parent_clay,stashcode_biom_surf_em,            &
    stashcode_biom_elev_em,stashcode_dms_conc_sw,                 &
    stashcode_biom_elev_em_h1, stashcode_biom_elev_em_h2,         &
    stashcode_clim_biogenic_aero,stashcode_clim_delta_aero,       &
    stashcode_ocff_surf_emiss, stashcode_ocff_hilev_emiss,        &
    stashcode_unfilt_orog,stashcode_ice_edge_inancil,             &
    stashcode_veg_frac,stashcode_z0,stashcode_runoff_coast_out,   &
    stashcode_soil_suction,stashcode_soil_moist,                  &
    stashcode_sil_orog_rough,stashcode_hlf_pk_to_trf_ht,          &
    stashcode_orog_x_grad,stashcode_orog_y_grad,                  &
    stashcode_clapp_hb,stashcode_frac_surf_type,                  &
    stashcode_soil_carbon_cont, stashcode_Ti_Mean,                &
    stashcode_Ti_Sig,             stashcode_rgrain,               &
    stashcode_tsurf_elev_surft,                                   &
    stashcode_snow_tile,          stashcode_snow_grnd,            &
    stashcode_snowdep_grd_tile,   stashcode_snowpack_bk_dens,     &
    stashcode_nsnow_layrs_tiles,  stashcode_snow_laythk_tiles,    &
    stashcode_snow_ice_tile,      stashcode_snow_liq_tile,        &
    stashcode_snow_T_tile,        stashcode_snow_laydns_tiles,    &
    stashcode_snow_grnsiz_tiles,                                  &
    stashcode_surf_sw_alb,stashcode_surf_nir_alb,                 &
    stashcode_iceberg_calving,stashcode_chlorophyll,              &
    stashcode_surf_z_curr,stashcode_surf_m_curr,                  &
    stashcode_urbhgt, stashcode_urbhwr, stashcode_urbwrr,         &
    stashcode_z0m_soil, stashcode_flake_depth,                    &
    stashcode_nitrogen_deposition, stashcode_soilnitro_dpm,       &
    stashcode_soilnitro_rpm, stashcode_soilnitro_bio,             &
    stashcode_soilnitro_hum, stashcode_soilcarb_rpm,              &
    stashcode_soilcarb_dpm, stashcode_soilcarb_bio,               &
    stashcode_soilcarb_hum,stashcode_soil_inorgnit,               &
    stashcode_crop_frac, stashcode_pasture_frac

USE Ereport_Mod, ONLY:                                           &
    Ereport

USE missing_data_mod, ONLY: imdi,rmdi

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER :: i_year          ! Curent Model Time
INTEGER :: i_month         !   "      "     "
INTEGER :: i_day           !   "      "     "
INTEGER :: i_hour          !   "      "     "
INTEGER :: i_minute        !   "      "     "
INTEGER :: i_second        !   "      "     "

INTEGER :: p_field         ! Size of horizontal fields
INTEGER :: p_rows          !
INTEGER :: u_field         !   "  "      "         "
INTEGER :: v_field         !   "  "      "         "
INTEGER :: rr_field        !   "  "      "         "
INTEGER :: Land_Field      !   "  "      "         "

INTEGER :: nlookups        ! Number of lookup tables
INTEGER :: Len_FixHd       ! Length of Fixed Header
INTEGER :: Len_IntHd       ! Length of Integer Header
INTEGER :: Len1_Lookup     ! First dimension of Lookup
INTEGER :: len_d1          ! Size of primary data array

INTEGER :: fixhd(len_fixhd,num_ancil_files)    ! Anc fixed headers
INTEGER :: inthd(len_inthd,num_ancil_files)    ! Anc Integer Headers
INTEGER :: lookup(len1_lookup,nlookups)        ! Anc Lookup Tables
INTEGER :: lookup_start(num_ancil_files)       ! Start of lookup tables for
                                            ! anc files opened.
INTEGER :: icode                        ! Return Code
INTEGER :: land_frac_pos                ! Land fraction
!                                      ! ancil position


REAL :: NS_Space              ! NS latitude spacing
REAL :: First_Lat             ! Latitude of first gridpoint
REAL :: d1(Len_D1)            ! INOUT Data array to hold fields
                           !       except TStar and Ice Fraction
REAL :: Ice_Fraction(P_Field) ! INOUT Ice Fraction, updated if
                           !       requested
REAL :: TStar (P_Field)       ! INOUT T Star, updated if requested
REAL :: TStar_Land_Ctile (P_Field)
!                            ! INOUT T*_land, updated if requested
REAL :: TStar_Sea_Ctile (P_Field)
!                            ! INOUT T*_sea, updated if requested
REAL :: TStar_Sice_Ctile (P_Field)
!                            ! INOUT T*_sice, updated if requested
REAL :: TStar_Anom(P_Field)   ! INOUT SST Anomaly, formed in Recon.
                           !       Added in model run if requested

REAL :: RLookup(Len1_Lookup,NLookups) ! REAL copy of Lookup Table

LOGICAL :: Land (P_Field)     ! Land mask

CHARACTER (LEN=errormessagelength) :: CMessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'REPLANCA_RCF_REPLANCA'

! Local Variables

! Buffers to hold values of ancillary data for time interpolation.
! Field of ancillary data held prior to selective updating.
REAL, ALLOCATABLE :: ancil1(:)
REAL, ALLOCATABLE :: ancil2(:)
REAL, ALLOCATABLE :: ancil_data(:)

REAL :: snow_change(p_field)   ! Fractional time of change of
                            ! snow cover
REAL :: ice_extent(p_field,2)  ! Fractional time of change
                            ! of ice cover
REAL :: pres_value(p_field)    ! Prescribed value of data when
                            ! controlling field is zero.
REAL :: no_ice_extent(p_field) ! Indicator for no sea ice
                            ! =0 if ice cover

REAL :: TStar_Land (P_Field)  ! T*_land, updated if requested
REAL :: TStar_Sea (P_Field)   ! T*_sea, updated if requested
REAL :: TStar_Sice (P_Field)  ! T*_sice, updated if requested
REAL :: TStar_Ssi (P_Field)   ! T*_ssi, updated if requested
REAL :: Flandg (P_Field)      ! Land fraction in gridbox
REAL :: Fland (Land_Field)    ! Land fraction in gridbox
!                            ! (on land-only points).
INTEGER :: i,j,l              ! Loop indices

INTEGER :: i1,i2,i3,ii       ! Array indices
INTEGER :: id                !
INTEGER :: im                !
INTEGER :: iy                !
INTEGER :: field             ! Current Ancil Ref Number.
INTEGER :: file_num          ! Anc File number
INTEGER :: irec              ! Anc request number

INTEGER :: interval          ! Interval between data times
INTEGER :: step              ! Number of data times skipped.
INTEGER :: months            ! )Used in calculation of position
INTEGER :: hours             ! )of data required
INTEGER :: period            ! Period of periodic data
INTEGER :: start_month       !
INTEGER :: level             ! Loop index for levels
INTEGER :: ancil_ref_days    ! Ancil.reference time in whole days
INTEGER :: ancil_ref_secs    ! Ancil.reference time in extra seconds
INTEGER :: day,sec           ! Times relative to reference time
INTEGER :: LEN
INTEGER :: row_length
INTEGER :: i_year1          ! Copy of Curent Model Time year
INTEGER :: i_month1         !   "      "     "          month
INTEGER :: i_day1           !   "      "     "          day
INTEGER :: i_hour1          !   "      "     "          hour

INTEGER :: anc_file_open    !  Stores which ancillary file is open.

INTEGER :: fileancil_sif  ! ancillary file number for sea-ice-fraction
INTEGER :: fileancil_sst  ! ancillary file number for SST

INTEGER, PARAMETER     :: filenameprovided=1
INTEGER :: field_size       !  Stores field size allowing use of
                         !  P_FIELD or RR_FIELD in calls.

LOGICAL :: linterpolate      ! Indicates whether time
                          ! interpolation needed.
LOGICAL :: lt_int_c          ! Indicates use of controlled time
                          ! interpolation
LOGICAL :: lmismatch         ! Used in header checks
LOGICAL :: single_time       ! Indicates that only one time is
                          ! available in data set
LOGICAL :: periodic          ! Data set is periodic
LOGICAL :: regular           ! Interval between data times in
                          ! dataset is regular in model timesteps.
LOGICAL :: lice_depth        ! T : Field is Ice Depth
LOGICAL :: lice_fraction     ! T : Field is Ice Fraction
LOGICAL :: lsnow_depth       ! T : Field is Snow depth

LOGICAL :: Sea (P_Field)      ! Sea mask

REAL :: zero         !
REAL :: time1        !  Times if data used in time interpolation
REAL :: time2        !
REAL :: time         !  Target time for time interpolation
REAL :: lat_p        !  Latitude of point


! STASHmaster entry
TYPE (STM_Record_Type), POINTER :: STM_Rec
INTEGER :: bot_level

! Dr Hook 
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! 0.  Set up local arrays for ANCIL data to required size.
field_size=MAX(p_field,rr_field)

! Set positions of seaicefrac and sst ancil files to imdi to allow checks to be
! made on actually finding them.
fileancil_sif = imdi
fileancil_sst = imdi

IF (PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
  WRITE(umMessage,'(A,I12)') "Rcf_replanca_rcf_replanca: chosen FIELD_SIZE=",&
        field_size
  CALL umPrint(umMessage,src='replanca-rcf_replanca')
END IF

ALLOCATE ( ancil1(field_size) )
ALLOCATE ( ancil2(field_size) )
ALLOCATE ( ancil_data(field_size) )

! 1.  Initialisation for atmosphere

icode=0
cmessage=' '
anc_file_open = 0

! Initialise ANCIL1/2. Includes Halos for MPP runs.
ancil1(:)=0.0
ancil2(:)=0.0

! Read in fractional land field first:
! use ANCIL_DATA as temporary storage

IF (l_ctile) THEN

  land_frac_pos = find_ancil_file_by_stash(stashcode_land_frac)
  ! find_ancil_file_by_stash return imdi if request not found
  IF (land_frac_pos > 0) THEN !land_frac configure requested
    CALL file_open (AncF_UnitNo,                                          &
                    ancil_files(land_frac_pos)%filename,                  &
                    LEN_TRIM(ancil_files(land_frac_pos)%filename),0,      &
                    filenameprovided,icode)

    ! DEPENDS ON: rcf_readflds
    CALL Rcf_ReadFlds                                                     &
                (AncF_UnitNo,1,                                           &
                 nlookup(find_ancil_req_by_stash(stashcode_land_frac)),   &
                 lookup(1,lookup_start(land_frac_pos)),                   &
                 len1_lookup,ancil_data,p_field,fixhd(1,land_frac_pos),   &
                 icode,cmessage)

    CALL file_close (AncF_UnitNo,                                         &
                     ancil_files(land_frac_pos)%filename,                 &
                     LEN_TRIM(ancil_files(land_frac_pos)%filename),       &
                     filenameprovided,icode)

    CALL umPrint( 'READ IN LAND FRACTION',src='replanca-rcf_replanca')

    DO i=1,p_field
      d1(ancil_add(land_frac_pos)+i-1)=ancil_data(i)
    END DO

    DO i=1,p_field
      flandg(i)=0.0
      IF (land(i)) flandg(i)=ancil_data(i)
    END DO

  ELSE             ! Land frac already read in from input dump

    l=0
    DO i=1,p_field
      flandg(i)=0.0
      IF (land(i)) THEN
        l=l+1
        flandg(i)=fland(l)
      END IF
    END DO
  END IF

  DO i=1,p_field
    ! If land or sea fraction is less than machine tolerance print warning
    IF (land(i) .AND. flandg(i) <  EPSILON(1.0)) THEN
      CALL umPrint('*****************WARNING********************', &
          src='replanca-rcf_replanca')
      CALL umPrint('LAND FRACTION IS LESS THAN MACHINE TOLERANCE', &
          src='replanca-rcf_replanca')
    END IF
    IF (.NOT. land(i) .AND. 1.0-flandg(i) <  EPSILON(1.0)) THEN
      CALL umPrint('*****************WARNING********************', &
          src='replanca-rcf_replanca')
      CALL umPrint('SEA FRACTION IS LESS THAN MACHINE TOLERANCE', &
          src='replanca-rcf_replanca')
    END IF

    IF (flandg(i) <= 0.0 .AND. land(i)) THEN
      CALL umPrint('*ERROR* a) LAND FRAC & LAND MASK INCONSISTENT', &
          src='replanca-rcf_replanca')
      icode = 800
      cmessage='replanca_rcf_replanca:ERR:LAND FRAC & MASK ARE ' &
      //'INCONSISTENT'
    END IF
    IF (flandg(i) >  0.0 .AND. .NOT. land(i)) THEN
      CALL umPrint('*ERROR* b) LAND FRAC & LAND MASK INCONSISTENT' &
          ,src='replanca-rcf_replanca')
      icode = 801
      cmessage='replanca_rcf_replanca:ERR:LAND FRAC & MASK ARE ' &
      //'INCONSISTENT'
    END IF

  END DO

ELSE                     ! Not coastal tiling:
  DO i=1,p_field
    IF (land(i)) THEN
      flandg(i)=1.0
    ELSE
      flandg(i)=0.0
    END IF
  END DO
END IF                   ! End of coastal tiling loop


DO i=1,p_field
  IF (flandg(i) <  1.0) THEN
    sea(i)=.TRUE.
  ELSE
    sea(i)=.FALSE.
  END IF
END DO


! Set up surface temperatures:
IF (l_ctile) THEN
  DO i=1,p_field
    tstar_land(i)=tstar_land_ctile(i)
    tstar_sea(i)=tstar_sea_ctile(i)
    tstar_sice(i)=tstar_sice_ctile(i)
    IF (ice_fraction(i) <= 0.0) THEN
      tstar_ssi(i)=tstar_sea(i)
    ELSE
      tstar_ssi(i)=ice_fraction(i)*tstar_sice(i)                &
             +(1.0-ice_fraction(i))*tstar_sea(i)
    END IF
  END DO
ELSE
  DO i=1,p_field
    tstar_land(i)=tstar(i)
    tstar_ssi(i)=tstar(i)
  END DO
END IF


! set up fileancil links for requried fields.
fileancil_sst = find_ancil_file_by_stash(stashcode_tstar)
fileancil_sif = find_ancil_file_by_stash(stashcode_icefrac)

! Initialise for valid time interpolation in reconfiguration mode.
ancil_ref_days = 0
ancil_ref_secs = 0

!  1.2 Allow for dependencies between fields
! Sea surface temperature must be updated when sea ice is updated

! Both surface current components must be updated together
! this was addressed and set up in inancila


!  Select method of time interpolation for SST. The interpolation
!  allows for sea ice if ice data is available at the same times
!  as the temperature data. Otherwise linear interpolation is used.

lt_int_c=.TRUE.

! Check if tstar has been requested
IF (find_ancil_req_by_stash(stashcode_tstar) /= imdi) THEN
  IF ( fileancil_sif == imdi ) THEN
    lt_int_c=.FALSE.
  ELSE
    ! FixHD item 10 is ancil time series type, 0 is single time.
    IF (fixhd(10,fileancil_sif) == 0) lt_int_c=.FALSE.
  END IF

  IF (lt_int_c) THEN
    DO i=21,41
      IF (fixhd(i,fileancil_sif) /=                          &
          fixhd(i,fileancil_sst)) THEN
        lt_int_c=.FALSE.
        WRITE(umMessage,'(3A)')' WARNING: controlled time interpolation ',&
        'for SST not available: Mismatch in SST and SEA-ICE ',    &
        'ancillary data times in FIXED HEADER'
        CALL umPrint(umMessage,src='replanca-rcf_replanca')
        WRITE(umMessage,'(A,I4,A,I6)')' position=',i,                     &
                    ' SEA-ICE=',fixhd(i,fileancil_sif)
        CALL umPrint(umMessage,src='replanca-rcf_replanca')
        WRITE(umMessage,'(A,I4,A,I6)')' position=',i,                     &
                    ' SST    =',fixhd(i,fileancil_sst)
        CALL umPrint(umMessage,src='replanca-rcf_replanca')
      END IF
    END DO
  END IF
END IF


!  Loop over ancillary fields(atmosphere)

DO irec=1,num_ancil_requests

  field =    ancil_requests(irec)%stashcode
  lice_depth=field == stashcode_icethick  ! required for l_amipii_ice_processing

  file_num = find_ancil_file_by_stash(field)

  IF ( file_num == imdi) THEN
    WRITE (cmessage,'(A,I4)')                                  &
           'Unable to identify ancillary to use for field', field
    icode = 30
    CALL Ereport( RoutineName, icode, Cmessage )
  END IF

  IF (file_num  /=  anc_file_open) THEN ! File is not open

    ! Only close file if anc_file_open is greater than 0. (e.g. first open).
    IF (anc_file_open > 0) THEN

      CALL file_close (AncF_UnitNo,                                   &
                       ancil_files(anc_file_open)%filename,           &
                       LEN_TRIM(ancil_files(anc_file_open)%filename), &
                       filenameprovided,icode)
    END IF

    CALL file_open (AncF_UnitNo,                                &
                    ancil_files(file_num)%filename,             &
                    LEN_TRIM(ancil_files(file_num)%filename),   &
                    0,filenameprovided,icode)

    IF (icode /= 0) THEN
      CALL umPrint( ' problem with file_open.',src='replanca-rcf_replanca')
      WRITE(umMessage,'(A,I5)') ' icode ',icode
      CALL umPrint(umMessage,src='replanca-rcf_replanca')
      WRITE(umMessage,'(2A)') ' cmessage ',cmessage
      CALL umPrint(umMessage,src='replanca-rcf_replanca')
      GO TO 9999
    END IF

    !       Record ancillary file opened
    anc_file_open = file_num

  END IF ! If file not already opened

  IF (lice_depth .AND. l_amipii_ice_processing) THEN

    ! Uses ice fraction set earlier in field loop.
    ! WARNING this will fail if the order of ancillary fields is ever
    ! changed so that ice-depth preceeds ice fraction
    ! Note : For complete sea ice cover
    !        Arctic ice depth    = 2m
    !        Antarctic ice depth = 1m
    ! For ice concentrations less than 1. ice depth is 1 or 2 times conc.
    ! This results in similar values to those from runs using ancillary
    ! files containing ice depths set to 1 or 2m.

    row_length=p_field/p_rows
    DO i=1,p_rows
      ! work out latitude in radians
      lat_p=first_lat-ns_space*(i-1)
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
    !       Sea ice thickness
    !       Update over all sea points (all sea ice points are the only
    !       ones strictly required, but this cannot be determined easily)

    DO i=1,p_field
      IF (sea(i)) THEN
        d1(ancil_add(irec)+i-1)=ancil_data(i)
      END IF
    END DO
    CYCLE ! Finished with ice_thickness field, now cycle main loop to next field
  END IF
  !     Update required for field

  IF ( mype == 0 ) THEN
    WRITE(umMessage,'(A,I6,2A)')'replanca_rcf_replanca: UPDATE REQUIRED'&
    //' FOR FIELD',field,' : ', ancil_requests(irec)%stash_name
    CALL umPrint(umMessage,src='replanca-rcf_replanca')
  END IF
  IF ( fixhd(10,file_num) < 0 .OR. fixhd(10,file_num) > 2 ) THEN
    WRITE(umMessage,'(A,I5)') ' file ',file_num
    CALL umPrint(umMessage,src='replanca-rcf_replanca')
    WRITE(umMessage,'(A,I5)') ' FIXHD(10,file) ',fixhd(10,file_num)
    CALL umPrint(umMessage,src='replanca-rcf_replanca')
    icode = 700 + file_num
    cmessage = 'replanca_rcf_replanca: Error in fixed '// &
    'header(10) of ancillary file'
    GO TO 9999
  END IF

  !     Check whether more than one data time available in data set

  single_time=fixhd(10,file_num) == 0

  !     Set default values for time interpolation

  linterpolate=.TRUE.

  IF (single_time) THEN
    linterpolate=.FALSE.
  END IF

  IF (field == stashcode_soil_temp            .OR.   &
      field == stashcode_vol_smc_wilt         .OR.   &
      field == stashcode_vol_smc_cri          .OR.   &
      field == stashcode_vol_smc_sat          .OR.   &
      field == stashcode_Ksat                 .OR.   &
      field == stashcode_thermal_capacity     .OR.   &
      field == stashcode_thermal_conduct) THEN
    linterpolate=.FALSE.
  END IF

  !  2.1 Find position of input record

  !     Default settings of search parameters if only one time present

  IF (single_time) THEN
    step=0
  ELSE

    periodic=fixhd(10,file_num) == 2
    regular=.TRUE.

    IF (.NOT. lcal360) THEN
      regular=fixhd(35,file_num) == 0 .AND. fixhd(36,file_num) == 0
      ! i.e. data at intervals of days/hours & non-periodic
      IF (periodic) regular=regular .AND. fixhd(37,file_num) == 0
      ! i.e. data at intervals of hours & periodic
    END IF

    !         Error checking on time information.

    IF ( fixhd(35,file_num) < 0 .OR.                            &
         fixhd(36,file_num) < 0 .OR. fixhd(36,file_num) > 12 .OR.   &
         regular .AND. (fixhd(37,file_num) < 0 .OR. fixhd(37,file_num) > 31 &
        .OR. fixhd(38,file_num) < 0 .OR. fixhd(38,file_num) > 24) ) THEN
      !           FIXHD(39-40) are not used by replanca_rcf_replanca.
      !           FIXHD(35-37) have already been used if not CAL360.
      icode = 700 + field
      cmessage = 'replanca_rcf_replanca: Error in validity '&
      //'time interval given in ancillary file'
      IF (lhook) CALL dr_hook(RoutineName, zhook_out, zhook_handle)
      RETURN
    END IF

    IF ( fixhd(21,file_num) < 0 .AND. .NOT. periodic              &
        .OR. .NOT. ( regular .AND. periodic ) .AND.                &
      ! If it is REGULAR & PERIODIC more detailed check is applied below
        ( fixhd(22,file_num) < 0 .OR. fixhd(22,file_num) > 12 .OR.     &
          fixhd(23,file_num) < 0 .OR. fixhd(23,file_num) > 31 .OR.     &
          fixhd(24,file_num) < 0 .OR. fixhd(24,file_num) > 24 .OR.     &
          fixhd(25,file_num) < 0 .OR. fixhd(25,file_num) > 60 .OR.     &
          fixhd(26,file_num) < 0 .OR. fixhd(26,file_num) > 60 ) ) THEN
      icode = 700 + field
      cmessage = 'replanca_rcf_replanca: Error in first '// &
      'validity time given in ancillary file'
      IF (lhook) CALL dr_hook(RoutineName, zhook_out, zhook_handle)
      RETURN
    END IF

    IF (.NOT. periodic) THEN

      !         If data taken from full time series of input data.

      ! DEPENDS ON: time2sec
      CALL time2sec(i_year,i_month,i_day,i_hour                &
                   ,i_minute,i_second                          &
                   ,ancil_ref_days,ancil_ref_secs,day,sec      &
                   ,lcal360)

      IF (regular) THEN

        !  2.1.1  Standard cases:360 day calender;
        !  2.1.1  or Gregorian calendar with
        !         interval between data times in days or hours
        !         updating interval may be regular in model timesteps,
        !         or (LGREG_MONTHLY=T) irregular in model timesteps,

        hours=sec/3600+day*24
        !  FInd time(in hours) of first ancillary data on file
        ! DEPENDS ON: time2sec
        CALL time2sec(fixhd(21,file_num),fixhd(22,file_num),          &
               fixhd(23,file_num),fixhd(24,file_num),                 &
               fixhd(25,file_num),fixhd(26,file_num),                 &
               ancil_ref_days,ancil_ref_secs,day,sec,         &
               lcal360)
        hours=hours-sec/3600-day*24

        IF (hours <  0) THEN
          icode=400+field
          cmessage='replanca_rcf_replanca: Current time '&
          //'precedes start time of data'
          IF (lhook) CALL dr_hook(RoutineName, zhook_out, zhook_handle)
          RETURN
        END IF

        !  FInd interval(in hours) between ancillary data on file
        interval=fixhd(35,file_num)*8640+fixhd(36,file_num)*720+      &
                 fixhd(37,file_num)*24+fixhd(38,file_num)

        ! Do not interpolate in time if data time exactly matches model time

        IF (MOD(hours,interval) == 0) THEN
          linterpolate=.FALSE.
        END IF

        step=hours/interval
        time=REAL(hours)
        time1=step*interval
        time2=(step+1)*interval

      ELSE

        !  2.1.2 Gregorian calender;ancillary data interval is in months or
        !        years,which is irregular in model timesteps.
        !  original code is inaccurate for this section - corrected code under
        !  use_lookup_dates_anc_time_interp makes use of dates in lookup headers
        !  For a real calendar year the mid-point of each month is different
        !  in terms of its hour and day. The old inaccurate method assumes
        !  the hour and day are taken from the fixhd values. These are only
        !  usually correct for the first month on the ancillary file.

        !  FInd interval(in months) between ancillary data on file
        interval=fixhd(35,file_num)*12+fixhd(36,file_num)
        months=i_year*12+i_month
        start_month=fixhd(21,file_num)*12+fixhd(22,file_num)
        months=months-start_month
        ! Check for time within month
        ! corrected code uses pp lookup header
        IF (use_lookup_dates_anc_time_interp) THEN   
          step=months/interval
          i2=nlookup(irec)+lookup_step(irec)*step
          i1=i2+lookup_start(file_num)-1
          ! Check against day and hour of actual lookup header not first field
          IF ((i_day*24+i_hour) <                             &
              (lookup(3,i1)*24+lookup(4,i1))) THEN
            months=months-1
          END IF
        ELSE           ! old less accurate code uses FIXHD
          IF ((i_day*24+i_hour) <                             &
              (fixhd(23,file_num)*24+fixhd(24,file_num))) THEN
            months=months-1
          END IF
        END IF ! use_lookup_dates_anc_time_interp

        IF (months <  0) THEN
          icode=400+field
          cmessage='replanca_rcf_replanca: Current time '&
          //'precedes start time of data'
          IF (lhook) CALL dr_hook(RoutineName, zhook_out, zhook_handle)
          RETURN
        END IF


        step=months/interval

        IF (use_lookup_dates_anc_time_interp) THEN       ! corrected code
          time=REAL(sec)/3600+REAL(day*24)
          ! correct calculation of dates uses lookup table dates not fixhd date
          i2=nlookup(Irec)+lookup_step(irec)*step
          i1=i2+lookup_start(file_num)-1
          i_year1=lookup(1,i1)
          i_month1=lookup(2,i1)
          i_day1=lookup(3,i1)
          i_hour1=lookup(4,i1)
          ! DEPENDS ON: time2sec
          CALL time2sec(i_year1,i_month1,i_day1,i_hour1,      &
              fixhd(25,file_num),fixhd(26,file_num),                  &
              ancil_ref_days,ancil_ref_secs,day,sec,          &
              lcal360)
          time1=REAL(sec)/3600+REAL(day*24)
          ! I1+LOOKUP_STEP(irec) correct pointer to next field as multi-level 
          ! fields are possible.
          IF (i1+lookup_step(irec) < nlookups) THEN
            ! Second check next lookup is the same field.
            IF (lookup(item_code,i1) ==                       &
                lookup(item_code,i1+lookup_step(irec))) THEN
              i_year1=lookup(1,i1+lookup_step(irec))
              i_month1=lookup(2,i1+lookup_step(irec))
              i_day1=lookup(3,i1+lookup_step(irec))
              i_hour1=lookup(4,i1+lookup_step(irec))
            ELSE
              icode = 500
              WRITE(cmessage,'(A,I6)') 'REPLANCA: error finding'&
              //' next lookup entry for field:', field
              GO TO 9999
            END IF
          ELSE
            icode = 500
            WRITE(cmessage,'(A,I6)') 'REPLANCA: error due to '&
            //'lookup out of range for field:', field
            GO TO 9999
          END IF

          ! DEPENDS ON: time2sec
          CALL time2sec(i_year1,i_month1,i_day1,i_hour1,      &
              fixhd(25,file_num),fixhd(26,file_num),                  &
              ancil_ref_days,ancil_ref_secs,day,sec,          &
              lcal360)
          time2=REAL(sec)/3600+REAL(day*24)

        ELSE   ! use_lookup_dates_anc_time_interp test - old inaccurate 
               ! code using FIXHD
          ! NB INTERVAL may be > 1 month
          months=step*interval
          ! Calculate data times for time interpolation
          time=REAL(sec)/3600+REAL(day*24)
          im=MOD(fixhd(22,file_num)+months-1,12)+1
          iy=fixhd(21,file_num)+(months+fixhd(22,file_num)-1)/12
          ! DEPENDS ON: time2sec
          CALL time2sec(iy,im,fixhd(23,file_num),fixhd(24,file_num),  &
              fixhd(25,file_num),fixhd(26,file_num),                  &
              ancil_ref_days,ancil_ref_secs,day,sec,          &
              lcal360)
          time1=REAL(sec)/3600+REAL(day*24)
          im=MOD(fixhd(22,file_num)+months+interval-1,12)+1
          iy=fixhd(21,file_num)+(months+interval+fixhd(22,file_num)-1)/12
          ! DEPENDS ON: time2sec
          CALL time2sec(iy,im,fixhd(23,file_num),fixhd(24,file_num),  &
              fixhd(25,file_num),fixhd(26,file_num),                  &
              ancil_ref_days,ancil_ref_secs,day,sec,          &
              lcal360)
          time2=REAL(sec)/3600+REAL(day*24)
        END IF     ! end use_lookup_dates_anc_time_interp test

        ! Do not interpolate in time if data time exactly matches model time

        IF (time == time1) THEN
          linterpolate=.FALSE.
        END IF

      END IF ! End of REGULAR/not REGULAR

    ELSE  ! PERIODIC data

      !  2.2   If data is taken from ancillary periodic data.

      ! DEPENDS ON: time2sec
      CALL time2sec(i_year,i_month,i_day,i_hour,             &
                    i_minute,i_second,                       &
                    ancil_ref_days,ancil_ref_secs,day,sec,   &
                    lcal360)

      IF (regular) THEN
        !  2.2.1 Standard cases:1) 360 day calender, with allowed periods of
        !        1 day, 1 month or 1 year;
        !
        !        2) Gregorian calender with update in hours,and period of
        !        data 1 day.
        !
        !        For both updating interval and number of
        !        data times to be skipped in data set calculated in hours.

        hours=sec/3600+day*24
        interval=fixhd(35,file_num)*8640+fixhd(36,file_num)*720+      &
                 fixhd(37,file_num)*24+fixhd(38,file_num)

        period=inthd(3,file_num)*interval

        !    Do not allow non-standard periods
        IF (lcal360) THEN
          IF (period/=8640 .AND. period/=720 .AND. period/=24) THEN
            icode=600+field
            cmessage='replanca_rcf_replanca: Non-standard '&
            //'period for periodic data'
            IF (lhook) CALL dr_hook(RoutineName, zhook_out, zhook_handle)
            RETURN
          END IF
        ELSE
          IF (period /= 24) THEN
            icode=600+field
            cmessage='replanca_rcf_replanca: Non-standard '&
            //'period for periodic data'
            IF (lhook) CALL dr_hook(RoutineName, zhook_out, zhook_handle)
            RETURN
          END IF
        END IF
        IF (period == 24) THEN
          ! Ancillary data interval in hour(s), period is 1 day

          iy=i_year
          im=i_month
          id=i_day
          IF (i_hour < fixhd(24,file_num)) hours=hours+24

        ELSE IF (period == 720) THEN
          ! Ancillary data interval in day(s) or hours , period is 1 month

          iy=i_year
          im=i_month
          id=fixhd(23,file_num)
          IF ((i_day*24+i_hour) <                            &
              (fixhd(23,file_num)*24+fixhd(24,file_num))) THEN
            hours=hours+720
          END IF

        ELSE IF (period == 8640) THEN
          ! Ancillary data interval in month(s)or days or hours, period is 1 
          ! year
          iy=i_year
          im=fixhd(22,file_num)
          id=fixhd(23,file_num)
          IF ((i_month*720+i_day*24+i_hour)<(fixhd(22,file_num)*720 &
                        +fixhd(23,file_num)*24+fixhd(24,file_num))) THEN
            hours=hours+8640
          END IF

        END IF

        ! DEPENDS ON: time2sec
        CALL time2sec(iy,im,id,fixhd(24,file_num),                &
                 fixhd(25,file_num),fixhd(26,file_num),               &
                 ancil_ref_days,ancil_ref_secs,day,sec,       &
                 lcal360)
        hours=hours-sec/3600-day*24

        ! Do not interpolate in time if data time exactly matches model time

        IF (MOD(hours,interval) == 0) THEN
          linterpolate=.FALSE.
        END IF
        step=hours/interval
        time=REAL(hours)
        time1=step*interval
        time2=(step+1)*interval

      ELSE  ! non regular case

        !  2.2.2 Gregorian calender,and data interval is in months,
        !        period is 1 year
        !        Updating interval and number of data times to be skipped
        !        calculated in months.

        time=REAL(sec)/3600+REAL(day*24)
        interval=fixhd(36,file_num)+fixhd(35,file_num)*12
        period=inthd(3,file_num)*interval
        IF (period /= 12) THEN
          icode=600+field
          cmessage='replanca_rcf_replanca: Non-standard '// &
          'period for periodic data'
          IF (lhook) CALL dr_hook(RoutineName, zhook_out, zhook_handle)
          RETURN
        END IF
        !  Difference between date now (month) & first date ancil file (month)
        months=i_month-fixhd(22,file_num)

        IF (use_lookup_dates_anc_time_interp) THEN 
          ! Correctly use day and hour from lookup header not fixhd which
          ! contains values for first field on ancillary file only.
          step=months/interval
          i2=nlookup(irec)+lookup_step(irec)*step
          i1=i2+lookup_start(file_num)-1
          !  Check for time within month - using ppheader information
          IF ((i_day*24+i_hour)<(lookup(3,i1)*24+lookup(4,i1))) THEN
            months=months-1
          END IF
          IF (months <  0) THEN
            months=months+12
          END IF
          ! recalculate STEP
          step=months/interval
          ! NB INTERVAL qmay be > 1 month
          months=step*interval
          iy=i_year
          im=MOD(fixhd(22,file_num)+months-1,12)+1
          IF (im >  i_month) iy=iy-1
          i2=nlookup(irec)+lookup_step(irec)*step
          i1=i2+lookup_start(file_num)-1
          ! DEPENDS ON: time2sec
          CALL time2sec(iy,im,lookup(3,i1),lookup(4,i1),      &
              fixhd(25,file_num),fixhd(26,file_num),                  &
              ancil_ref_days,ancil_ref_secs,day,sec,lcal360)
          time1=REAL(sec)/3600+REAL(day*24)
          !  Calculate  TIME2 for second ancillary data time
          !  set IY correctly for time interpolation calculations
          iy=i_year
          im=MOD(fixhd(22,file_num)+months+interval-1,12)+1
          IF (im <  i_month) iy=iy+1
          i1=(im-1)/interval
          i2=nlookup(irec)+lookup_step(irec)*i1
          i1=i2+lookup_start(file_num)-1
          ! DEPENDS ON: time2sec
          CALL time2sec(iy,im,lookup(3,i1),lookup(4,i1),      &
              fixhd(25,file_num),fixhd(26,file_num),                  &
              ancil_ref_days,ancil_ref_secs,day,sec,lcal360)
          time2=REAL(sec)/3600+REAL(day*24)

        ELSE   ! original code inaccurate use of FIXHD dates
          !  Check for time within month
          IF ((i_day*24+i_hour) <                             &
              (fixhd(23,file_num)*24+fixhd(24,file_num))) THEN
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
          im=MOD(fixhd(22,file_num)+months-1,12)+1
          IF (im >  i_month) iy=iy-1
          ! DEPENDS ON: time2sec
          CALL time2sec(iy,im,fixhd(23,file_num),fixhd(24,file_num),  &
              fixhd(25,file_num),fixhd(26,file_num),                  &
              ancil_ref_days,ancil_ref_secs,day,sec,          &
              lcal360)
          time1=REAL(sec)/3600+REAL(day*24)
          !  Calculate  TIME2 for second ancillary data time
          !  set IY correctly for time interpolation calculations
          iy=i_year
          im=MOD(fixhd(22,file_num)+months+interval-1,12)+1
          IF (im <  i_month) iy=iy+1
          ! DEPENDS ON: time2sec
          CALL time2sec(iy,im,fixhd(23,file_num),fixhd(24,file_num),  &
              fixhd(25,file_num),fixhd(26,file_num),                  &
              ancil_ref_days,ancil_ref_secs,day,sec,          &
              lcal360)
          time2=REAL(sec)/3600+REAL(day*24)
        END IF  ! end use_lookup_dates_anc_time_interp test

        ! Do not interpolate in time if data time exactly matches model time

        IF (time == time1) THEN
          linterpolate=.FALSE.
        END IF

      END IF  ! regular/non-regular

    END IF  ! non-periodic/periodic

  END IF ! singletime/non-singletime

  !  2.3   Check STASH Code

  i2=nlookup(irec)+lookup_step(irec)*step

  i1=lookup(item_code,i2+lookup_start(file_num)-1)

  lmismatch=.FALSE.

  IF (PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
    WRITE(umMessage,'(A,I12)') ' Information used in checking ancillary '&
    //'data set: position of lookup table in dataset:',i2
    CALL umPrint(umMessage,src='replanca-rcf_replanca')
    WRITE(umMessage,'(A,I12)')' Position of first lookup table referring'&
    //' to data type ',nlookup(irec)
    CALL umPrint(umMessage,src='replanca-rcf_replanca')
    WRITE(umMessage,'(2(A,I8))') &
        ' Interval between lookup tables referring'&
    //' to data type ', lookup_step(irec),' Number of steps', step
    CALL umPrint(umMessage,src='replanca-rcf_replanca')
    WRITE(umMessage,'(2(A,I8))') ' STASH code in dataset ',i1,          &
          '  STASH code requested ',ancil_requests(irec) % stashcode
    CALL umPrint(umMessage,src='replanca-rcf_replanca')
    WRITE(umMessage,'(A,I12)') '''Start'' position of lookup tables for '&
    //'dataset in overall lookup array ' ,lookup_start(file_num)
    CALL umPrint(umMessage,src='replanca-rcf_replanca')
  END IF

  IF (i1 /= ancil_requests(irec) % stashcode) THEN
    WRITE(umMessage,'(A,3I8)') &
    ' I1,ancil_requests(irec) % stashcode', &
    i1,ancil_requests(irec)%stashcode,irec
    CALL umPrint(umMessage,src='replanca-rcf_replanca')
    lmismatch=.TRUE.
  END IF

  !  Error exit if checks fail

  IF (lmismatch) THEN
    icode=200+field
    cmessage='replanca_rcf_replanca: PP HEADERS ON ANCILLARY '//&
    'FILE DO NOT MATCH'
    IF (lhook) CALL dr_hook(RoutineName, zhook_out, zhook_handle)
    RETURN
  END IF

  IF (linterpolate .AND. .NOT. single_time) THEN
    !  Check time interpolation factors
    IF (time < time1 .OR. time > time2) THEN
      CALL umPrint(' Information used in interpolation/replacement:', &
          src='replanca-rcf_replanca')
      WRITE(umMessage,'(A,I8)')' Time of first data=', time1
      CALL umPrint(umMessage,src='replanca-rcf_replanca')
      WRITE(umMessage,'(A,I8)')' Validity Time for update=', time
      CALL umPrint(umMessage,src='replanca-rcf_replanca')
      WRITE(umMessage,'(A,I8)')' Time of second data=', time2
      CALL umPrint(umMessage,src='replanca-rcf_replanca')

      icode=500+field
      cmessage='replanca_rcf_replanca: TIME INTERPOLATION ERROR'
      IF (lhook) CALL dr_hook(RoutineName, zhook_out, zhook_handle)
      RETURN
    END IF
  END IF

  !  3   Loop over levels of ancillary data for field I
  !  Reset pointer for dataset

  !  Includes loop over X and Y components of surface currents

  lice_fraction=field == stashcode_icefrac
  lsnow_depth=field == stashcode_mean_snow
  lice_depth=field == stashcode_icethick

  ! Find whether we need to handle the extra level due to ancillary not 
  ! supporting theta level 0.
  STM_Rec => Rcf_Exppx( 1, ancil_requests(irec) % section, &
                           ancil_requests(irec) % item )

  IF ( STM_Rec % lb_code > 0) THEN
    CALL Rcf_Level_Code( STM_Rec % lb_code, bot_level, Output_Grid )
  ELSE
    bot_level = 1
  END IF

  ! Force bottom level to be 1 (if not 0)
  IF (bot_level /= 0) THEN
    bot_level = 1
  END IF

  ! Only want to loop upto levels in ancillary.
  DO level=1,levels(irec)

    !  Do not go through loop for ice edge or snow edge

    IF (.NOT. ((lice_fraction .OR. lsnow_depth) .AND. level == 2)) THEN

      !       Check to see if field is one of the River Routing ones.
      IF ( field  ==  stashcode_riv_storage .OR.                      &
          field  ==  stashcode_riv_sequence .OR.                      &
          field  ==  stashcode_riv_direction ) THEN

        IF (PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
          CALL umPrint( 'Rcf_replanca_rcf_replanca: '&
              //'Resetting Ancil field size to RR_FIELD')
        END IF

        field_size=rr_field

      ELSE ! If not River Routing, then assume field size is PFIELD

        field_size=p_field

      END IF

      !  3.1 Read data for single level of ancillary field.

      IF (.NOT. lice_fraction) THEN
        ! AMIPII case ice depth field not read from ancillary file
        IF (.NOT. (lice_depth .AND. l_amipii_ice_processing)) THEN

          ! DEPENDS ON: rcf_readflds
          CALL Rcf_ReadFlds                                   &
              (AncF_UnitNo,1,i2,lookup(1,lookup_start(file_num)), &
              len1_lookup,ancil1,field_size,fixhd(1,file_num),    &
              icode,cmessage)

        END IF

        IF (icode > 0) THEN
          icode=field+100
          cmessage='replanca_rcf_replanca :I/O ERROR '
          GO TO 9999
        END IF

      ELSE

        !  If ice-fraction,read fractional time field as well
        !        UNLESS IT IS A SINGLE TIME FIELD
        !  If snow-depth,read fractional time field as well only if time
        !  interpolation required.

        IF (.NOT. single_time .AND. .NOT. l_amipii_ice_processing) THEN
          IF (lookup(item_code,i2+lookup_start(file_num)) == 38) THEN

            ! DEPENDS ON: rcf_readflds
            CALL Rcf_ReadFlds                                   &
                (AncF_UnitNo,2,i2,lookup(1,lookup_start(file_num)), &
                len1_lookup,ice_extent,field_size,fixhd(1,file_num),&
                icode,cmessage)

            IF (icode > 0) THEN
              icode=field+100
              cmessage='replanca_rcf_replanca :I/O ERROR '
              GO TO 9999
            END IF

          ELSE
            icode=field+100
            cmessage='replanca_rcf_replanca :ICE CHANGE DATA MISSING'
            GO TO 9999
          END IF
        ELSE ! single time or l_amipii_ice_processing - ie no time change field

          ! DEPENDS ON: rcf_readflds
          CALL Rcf_ReadFlds                                   &
              (AncF_UnitNo,1,i2,lookup(1,lookup_start(file_num)), &
              len1_lookup,ice_extent,field_size,fixhd(1,file_num),&
              icode,cmessage)

          IF (icode >  0) THEN
            icode=field+100
            cmessage='replanca_rcf_replanca :I/O ERROR '
            GO TO 9999
          END IF
        END IF
      END IF

      IF (lsnow_depth .AND. linterpolate) THEN
        IF (lookup(item_code,i2+lookup_start(file_num)) == 27) THEN

          ! DEPENDS ON: rcf_readflds
          CALL Rcf_ReadFlds                                      &
              (AncF_UnitNo,1,i2+1,                               &
              lookup(1,lookup_start(file_num)),                      &
              len1_lookup,snow_change,field_size,fixhd(1,file_num),  &
              icode,cmessage)

          IF (icode >  0) THEN
            icode=field+100
            cmessage='replanca_rcf_replanca :I/O ERROR '
            GO TO 9999
          END IF

        ELSE
          icode=field+100
          cmessage='replanca_rcf_replanca :SNOW CHANGE DATA MISSING'
          GO TO 9999
        END IF
      END IF

      !  If sea surface temperature or other ice fields, read ice fraction
      !  and fractional time field if not already pressent and if required
      !  by time interpolation.

      IF (field == stashcode_icethick .OR.           &
         (field == stashcode_tstar .AND. lt_int_c) ) THEN
        IF (.NOT. ANY( ancil_requests%stashcode == stashcode_icefrac)) THEN
          i3 = nlookup(fileancil_sif) + lookup_step(fileancil_sif)*step  &
              + lookup_start(fileancil_sif)

          IF ( lookup(item_code,i3)  ==  38 ) THEN

            ! DEPENDS ON: rcf_readflds
            CALL Rcf_ReadFlds                                 &
                (AncF_UnitNo,2,                               &
                nlookup(fileancil_sif)+lookup_step(fileancil_sif)*step,  &
                lookup(1,lookup_start(fileancil_sif)),   &
                len1_lookup,ice_extent,                       &
                field_size,fixhd(1,fileancil_sif),       &
                icode,cmessage)

            IF (icode /= 0) THEN
              icode=field+100
              cmessage='replanca_rcf_replanca :I/O ERROR '
              GO TO 9999
            END IF
            IF ( rlookup(bmdi,i3-1)  /=  rmdi ) THEN
              icode = 700 + field
              cmessage = 'replanca_rcf_replanca: nonstandard lookup'&
                  //' RMDI in ancil file sea-ice chge times'
              GO TO 9999
            END IF

          ELSE
            icode=field+100
            cmessage='replanca_rcf_replanca :ICE FIELD DATA MISSING'
            GO TO 9999
          END IF
        END IF
      END IF

      !  3.3 If time interpolation required, read second record

      IF (linterpolate) THEN

        i1=i2+ lookup_step(irec)
        IF (i1 <= fixhd(152,file_num)) THEN

          IF ( lookup(item_code,lookup_start(file_num)+i1-1) /= &
               lookup(item_code,lookup_start(file_num)+i2-1) ) THEN
            icode=field+100
            cmessage='replanca_rcf_replanca: start and end fields' // &
                     ' are different.'
            GO TO 9999
          END IF

          ! AMIP II and ice depth don't read in ice depth field
          IF (.NOT. (l_amipii_ice_processing .AND. lice_depth)) THEN

            ! DEPENDS ON: rcf_readflds
            CALL Rcf_ReadFlds                                 &
                (AncF_UnitNo,1,i1,                            &
                lookup(1,lookup_start(file_num)),                 &
                len1_lookup,ancil2,field_size,fixhd(1,file_num),  &
                icode,cmessage)

          END IF

          IF (icode /= 0) THEN
            icode=field+100
            cmessage='replanca_rcf_replanca :I/O ERROR '
            GO TO 9999
          END IF

        ELSE !end of data on file

          !   If end of data has been reached go back to the start.If data is
          !   periodic.
          !   Otherwise cancel time interpolation

          IF (periodic) THEN

            i1 = nlookup(irec) + level - 1
            IF ( lookup(item_code,lookup_start(file_num)+i1-1) /= &
                 lookup(item_code,lookup_start(file_num)+i2-1) ) THEN
              icode=field+100
              cmessage='replanca_rcf_replanca: start and end fields'// &
                       ' are different.'
              GO TO 9999
            END IF

            ! DEPENDS ON: rcf_readflds
            CALL Rcf_ReadFlds                               &
                (AncF_UnitNo,1,i1,                          &
                lookup(1,lookup_start(file_num)),               &
                len1_lookup,ancil2,field_size,fixhd(1,file_num),&
                icode,cmessage)

            IF (icode /= 0) THEN
              icode=field+100
              cmessage='replanca_rcf_replanca :I/O ERROR '
              GO TO 9999
            END IF
          ELSE
            WRITE(umMessage,'(A,A)') &
                'REPLANCA: Reached end of ancillary ',   &
                'switched off time interpolation.'
            CALL umPrint(umMessage,src='replanca-rcf_replanca')
            linterpolate=.FALSE.
          END IF
        END IF ! End of position on file test

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
          IF ( rlookup(bmdi,lookup_start(file_num)+i2) /= rmdi ) THEN
            icode = 700 + field
            cmessage = 'replanca_rcf_replanca: nonstandard lookup '&
                //'RMDI in ancil file snow chge times'
            GO TO 9999
          END IF

          CALL t_int_c (ancil1,time1,ancil2,time2,ancil_data,      &
              time,field_size,snow_change,ancil1,pres_value)

          ! Ice fraction: ice depth set equal to zero if no ice

        ELSE IF (field == stashcode_icefrac .OR.      &
                  field == stashcode_icethick) THEN
          IF (field == stashcode_icefrac) THEN
            ! For the call to T_INT_C, need to know BMDI is OK for
            ! ICE_EXTENT(1,2) which was read in from position I1+1
            IF (.NOT. l_amipii_ice_processing) THEN
              IF (rlookup(bmdi,lookup_start(file_num)+i1) /= rmdi) THEN
                icode = 700 + field
                cmessage = 'replanca_rcf_replanca: nonstandard lookup'&
                    //' RMDI in ancil file sea-ice chge times'
                GO TO 9999
              END IF
            END IF

            IF (l_amipii_ice_processing) THEN
              ! linear uncontrolled time interpolation
              CALL t_int (ice_extent,time1,ancil2,time2,ancil_data, &
                  time,field_size)

              ! For AMIP II strictly ice concentrations should range between
              ! 0.0 and 1.0 but because of assumptions on T* made by the
              ! boundary layer and radiation schemes ice concentrations are 
              ! restricted to 0.3 to 1.0. This will allow SSTs in areas of less 
              ! than 30% ice to be used rather than TFS=-1.8C.

              DO i=1,field_size
                IF (ancil_data(i) <  0.3) ancil_data(i)=0.0
                IF (ancil_data(i) >  1.0) ancil_data(i)=1.0
              END DO

            ELSE       ! non AMIPII option
              DO i=1,field_size
                pres_value(i)=0
              END DO

              CALL t_int_c (ice_extent,time1,ancil2,time2,ancil_data, &
                  time,field_size,ice_extent(1,2),ice_extent,pres_value)

            END IF     ! end AMIPII test

          ELSE IF (field == stashcode_icethick) THEN

            DO i=1,field_size
              pres_value(i)=0
            END DO

            CALL t_int_c (ancil1,time1,ancil2,time2,ancil_data,       &
                time,field_size,ice_extent(1,2),ice_extent,pres_value)

          END IF

          ! Sea surface temperature, set equal to TFS if ice present

        ELSE IF (field == stashcode_tstar .AND. lt_int_c) THEN
          IF (l_amipii_ice_processing) THEN

            CALL t_int (ancil1,time1,ancil2,time2,ancil_data,         &
                time,field_size)
            ! remove any T below TFS
            DO i=1,field_size
              IF (ancil_data(i) <  tfs)  ancil_data(i)=tfs
            END DO

          ELSE     ! non AMIPII option

            DO i=1,field_size
              pres_value(i)=tfs

              ! Set no_ice_extent indicator for controlled SST interpolation
              IF (ice_extent(i,1) == 0) THEN
                no_ice_extent(i)=1.0
              ELSE
                no_ice_extent(i)=0.0
              END IF
            END DO

            CALL t_int_c (ancil1,time1,ancil2,time2,ancil_data,       &
                time,field_size,ice_extent(1,2),no_ice_extent,pres_value)

          END IF   ! end AMIPII test
          ! Otherwise linear interpolation in time, unless missing
          ! data indicator present at either time.

        ELSE

          ! Time interpolation checks the data against the standard missing data
          ! indicator - check that the field is labelled as using the same 
          ! one. (It is to have the right I1 here that I3 is used above.)
          IF ( rlookup(bmdi,lookup_start(file_num)+i1-1) /= rmdi .OR.    &
              rlookup(bmdi,lookup_start(file_num)+i2-1) /= rmdi ) THEN
            WRITE(umMessage,'(A,2F12.1)') 'LOOKUPS:',                 &
                rlookup(bmdi,lookup_start(file_num)+i1-1),                &
                rlookup(bmdi,lookup_start(file_num)+i2-1)
            CALL umPrint(umMessage,src='replanca-rcf_replanca')
            icode = 700 + field
            cmessage = 'replanca_rcf_replanca: Nonstandard MDI in '&
                //'lookup of ancil file'
            GO TO 9999
          END IF

          LEN=field_size
          !   Ozone, test for zonal mean or full field
          IF (field == stashcode_ozone) THEN
            IF (lookup(lbnpt,lookup_start(file_num)+i2-1) == 1) THEN
              LEN=p_rows
            END IF
            !   Tropopause-based ozone, test for zonal mean or full field.
            !   Currently the same test as for conventional ozone.
            !       ELSE IF (FIELD == 110) THEN
            !         IF (LOOKUP(LBNPT,LOOKUP_START(file_num)+I2-1) == 1) THEN
            !           LEN=P_ROWS
            !         END IF
            !   Cariolle ozone, test for zonal mean or full field.
            !   Currently same test as for conventional ozone.
          ELSE IF (field >= stashcode_ozone_tracer .AND.    &
                   field <= stashcode_o3_colo3) THEN
            IF (lookup(lbnpt,lookup_start(file_num)+i2-1) == 1) THEN
              LEN=p_rows
            END IF
          END IF

          CALL t_int(ancil1,time1,ancil2,time2,ancil_data,            &
              time,LEN)

        END IF ! End Lsnow_depth

        ! If no interpolation, copy data into final array

      ELSE ! no interpolation
        IF (lice_fraction) THEN
          IF (l_amipii_ice_processing) THEN
            DO i=1,field_size

              ancil_data(i)=ice_extent(i,1)

              ! For AMIP II strictly ice concentrations should range between
              ! 0.0 and 1.0 but because of assumptions on T* made by the
              ! boundary layer and radiation schemes ice concentrations are 
              ! restricted to 0.3 to 1.0. This will allow SSTs in areas of less 
              ! than 30% ice to be used rather than TFS=-1.8C.

              IF (ancil_data(i) <  0.3) ancil_data(i)=0.0
              IF (ancil_data(i) >  1.0) ancil_data(i)=1.0

            END DO
          ELSE           ! non AMIP II option
            DO i=1,field_size
              ancil_data(i)=ice_extent(i,1)
            END DO
          END IF           ! end of AMIPII test
        ELSE IF (l_amipii_ice_processing .AND. field == stashcode_tstar) THEN
          DO i=1,field_size
            ancil_data(i)=ancil1(i)
            IF (ancil_data(i) <  tfs) ancil_data(i)=tfs
          END DO
        ELSE
          DO i=1,field_size
            ancil_data(i)=ancil1(i)
          END DO

        END IF
      END IF !End interpolate/no interpolate

      !  3.5 Updating action for each field at each level
      !      Fields replaced except that Sea Surface Temperature may be
      !      incremented. Take appropriate action for each field.

      IF (field == stashcode_lsm .OR. field == stashcode_orog .OR.    &
          field == stashcode_ozone .OR. field == stashcode_SO2_emiss  &
          .OR. field == stashcode_dimethyl_sul_emiss                  &
        !  .OR.FIELD == 41.OR.FIELD == 42.OR.FIELD == 43              &
          .OR. field == stashcode_total_aero_emiss                    &
          .OR. field ==   stashcode_total_aero                        &
                                   ! multi-level murk
          .OR. (field >= 301 .AND. field <= 320 .AND. l_ukca)         &
                                   ! single-level user ancillaries
          .OR. (field >= stashcode_ammonia_gas_emiss .AND.            &
               field <= stashcode_soot_hi_lev)                        &
                                   !NH3,soot aerosol emissions
          .OR. (field >= stashcode_3d_nat_so2_em .AND.                &
               field <= stashcode_hi_SO2_emiss_emiss)                 &
                                   !Sulphur cycle
          .OR. field == stashcode_CO2_surf_emiss                      &
                                   !CO2 EMISSIONS
        !  .OR.FIELD == 82                                            &
                                   !HADCM2 sulphate aerosol
          .OR. (field >= 321 .AND. field <= 340)                      &
                                   !multi-level user ancillaries
          .OR. (field >= stashcode_dust_parent_clay .AND.             &
               field <= stashcode_soil_massfrac6)                     &
                                   !mineral dust fields
          .OR. (field >= stashcode_biom_surf_em .AND.                 &
               field <= stashcode_biom_elev_em)                       &
                                   !Biomass emissions
          .OR. field == stashcode_dms_conc_sw                         &
                                   !Seawater DMS concentration
          .OR. (field >= stashcode_biom_elev_em_h1 .AND.              &
                field <= stashcode_biom_elev_em_h2)                   &
                                   !Injection heights for biomass emiss
          .OR. (field >= stashcode_riv_sequence .AND.                 &
               field <=stashcode_riv_storage)                         &
                                   !River routing
          .OR. (field >= stashcode_clim_biogenic_aero .AND.           &
               field <= stashcode_clim_delta_aero)                    &
                                   !Aerosol climatologies
          .OR. (field >= stashcode_ozone_tracer .AND.                 &
               field <= stashcode_o3_colo3)                           &
                                   !Cariolle ozone ancillaries
          .OR. (field >=stashcode_ocff_surf_emiss  .AND.              &
               field <= stashcode_ocff_hilev_emiss)                   &
                                   !OCFF emissions
          .OR. field ==  stashcode_unfilt_orog                        &
                                                 !Unfiltered orography
          .OR. (field >= stashcode_urbhgt .AND.                        &
               field <= stashcode_urbwrr)                             &
                                   !Urban morphology
                        ) THEN

        !  3.5.0 Updates at all points

        LEN=field_size
        !   Ozone, test for zonal mean or full field
        IF (field == stashcode_ozone) THEN
          IF (lookup(lbnpt,lookup_start(file_num)+i2-1) == 1) THEN
            LEN=p_rows
          END IF
          !   Tropopause-based ozone, test for zonal mean or full field.
          !   Currently the same test as for conventional ozone.
          !     ELSE IF (FIELD == 110) THEN
          !       IF (LOOKUP(LBNPT,LOOKUP_START(file_num)+I2-1) == 1) THEN
          !         LEN=P_ROWS
          !       END IF
          !   Cariolle ozone, test for zonal mean or full field.
          !   Currently same test as for conventional ozone.
        ELSE IF (field >= stashcode_ozone_tracer .AND.   &
                 field <= stashcode_o3_colo3) THEN
          IF (lookup(lbnpt,lookup_start(file_num)+i2-1) == 1) THEN
            LEN=p_rows
          END IF
        END IF

        DO i=1,LEN
          d1(ancil_add(irec)+i-1+(level-bot_level)*LEN)=ancil_data(i)
        END DO

        ! Copy theta level 1 to theta level 0
        IF (level == 1 .AND. bot_level == 0) THEN
          DO i=1,LEN
            d1(ancil_add(irec)+i-1)=ancil_data(i)
          END DO
        END IF

        !  3.5.1 Updates over all land points

      ELSE IF ((field > stashcode_orog  .AND.                       &
                field < stashcode_ice_edge_inancil )                &
          .OR. (field == stashcode_mean_snow .OR.                   &
                field== stashcode_soil_temp   .OR.                  &
                field == stashcode_vol_smc_wilt .OR.                &
                field == stashcode_vol_smc_cri)                     &
          .OR. (field >= stashcode_vol_smc_sat .AND.                &
                field <=  stashcode_thermal_conduct )               &
          .OR. (field == stashcode_veg_frac   .OR.                  &
                field== stashcode_z0)                               &
          .OR. (field == stashcode_runoff_coast_out) .OR.           &
               (field == stashcode_soil_suction     .OR.            &
                field == stashcode_soil_moist)                      &
          .OR. (field >= 301 .AND. field <= 320 .AND. .NOT. l_ukca) &
                                ! single level user ancillaries
          .OR. (field == stashcode_sil_orog_rough  .OR.             &
                field == stashcode_hlf_pk_to_trf_ht)                &
                                ! Orographic roughness
          .OR. (field == stashcode_orog_x_grad    .OR.              &
                field == stashcode_orog_y_grad )                    &
                                ! Orographic X & Y gradients
          .OR. (field == stashcode_clapp_hb)                        &
                                ! MOSES-I
          .OR. (field >= stashcode_frac_surf_type .AND.             &
                field <= stashcode_soil_carbon_cont)                &
                                ! MOSES-II
          .OR. (field >= stashcode_Ti_Mean .AND.                    &
                field <= stashcode_Ti_Sig)                          &
                                ! LSH Topographic index fields
          .OR. (field == stashcode_land_frac)                       &
                                !COASTAL TILING

                                ! Pool Soil Carbon 
          .OR. (field == stashcode_soilcarb_dpm)                 &
          .OR. (field == stashcode_soilcarb_rpm)                 &
          .OR. (field == stashcode_soilcarb_bio)                 &
          .OR. (field == stashcode_soilcarb_hum)                 &
                                ! Pool Soil Nitrogen 
          .OR. (field == stashcode_soilnitro_dpm)                &
          .OR. (field == stashcode_soilnitro_rpm)                &
          .OR. (field == stashcode_soilnitro_bio)                &
          .OR. (field == stashcode_soilnitro_hum)                &
                                ! Nitrogen Depositon          
          .OR. (field == stashcode_nitrogen_deposition)          &
                                ! Inorganic Nitrogen Pool
          .OR. (field == stashcode_soil_inorgnit)                &
                                ! Landuse Fracs
          .OR. (field == stashcode_crop_frac)                    &
          .OR. (field == stashcode_pasture_frac)                 &
                                
          .OR. (field == stashcode_rgrain)                       &
          .OR. (field == stashcode_snow_tile)                    &
          .OR. (field == stashcode_snow_grnd)                    &
          .OR. (field == stashcode_snowdep_grd_tile)             &
          .OR. (field == stashcode_snowpack_bk_dens)             &
          .OR. (field == stashcode_nsnow_layrs_tiles)            &
          .OR. (field == stashcode_snow_laythk_tiles)            &
          .OR. (field == stashcode_snow_ice_tile)                &
          .OR. (field == stashcode_snow_liq_tile)                &
          .OR. (field == stashcode_snow_T_tile)                  &
          .OR. (field == stashcode_snow_laydns_tiles)            &
          .OR. (field == stashcode_snow_grnsiz_tiles)            &
          .OR. (field == stashcode_tsurf_elev_surft)             &
                                !Snow Fields
          .OR. (field >= stashcode_surf_sw_alb  .AND.               &
                field <= stashcode_surf_nir_alb)                    &
                                ! Land surface albedos obs/clim
          .OR. (field == stashcode_z0m_soil)                        &
          .OR. (field == stashcode_flake_depth)                     &
         ) THEN
        !  Set default value of Z0 over sea
        IF (field == stashcode_z0) THEN
          DO i=1,p_field
            IF (sea(i)) THEN
              d1(ancil_add(irec)+i-1+(level-1)*p_field)=10.0e-4
            END IF
          END DO
        END IF

        DO i=1,p_field
          IF (land(i)) THEN
            d1(ancil_add(irec)+i-1+(level-1)*p_field)=ancil_data(i)
          END IF
        END DO

        !   Reset TSTAR to TM if snow cover present

        IF (lsnow_depth) THEN
          DO i=1,p_field
            IF (land(i) .AND. ancil_data(i) >  0.0) THEN
              IF (tstar_land(i) >  tm) tstar_land(i)=tm
            END IF
          END DO
        END IF

        ! Iceberg calving for the OASIS coupler:
      ELSE IF (field == stashcode_iceberg_calving) THEN
        DO i=1,u_field
          d1(ancil_add(irec)+i-1)=ancil_data(i)
        END DO

        !  3.5.2 Ice fraction

      ELSE IF (field == stashcode_icefrac) THEN

        DO i=1,p_field
          ice_fraction(i)=0.0
          IF (sea(i)) THEN
            ice_fraction(i)=ancil_data(i)
          END IF
        END DO

        !  Reduce TSTAR to TFS where ice fraction greater than zero
        ! Required at present because radiation and boundary layer codes
        ! assume T* is TFS and ignore any value set in TSTAR.

        IF (.NOT. l_leads_temp_prog) THEN
          DO i=1,p_field
            IF (ice_fraction(i) >  0.0) THEN
              tstar_ssi(i)=MIN(tstar_ssi(i),tfs)
            END IF
          END DO
        END IF

        !  3.5.3 Sea surface temperatures for atmosphere, allow fields to be
        !        incremented rather than replaced

      ELSE IF (field == stashcode_tstar) THEN

        IF (l_sstanom) THEN
          DO i=1,p_field
            tstar_anom(i)=0.0
          END DO
        END IF

        DO i=1,p_field
          !           Calculate SST anomalies over all open sea points,
          !           but ignore anomalies over grid points with sea-ice during
          !           the forecast stage when updating SSTs
          IF (.NOT. land(i)) THEN
            IF (l_sstanom) THEN
              tstar_anom(i)=tstar_ssi(i)-ancil_data(i)
            ELSE
              tstar_sea(i)=ancil_data(i)
              IF (ice_fraction(i) <= 0.0) tstar_ssi(i)=tstar_sea(i)
            END IF
          END IF
        END DO

        !  3.5.4  Sea-point fields, where no special action is required:
        ! Sea-ice Thickness
      ELSE IF (field == stashcode_icethick) THEN
        DO i=1,p_field
          IF (sea(i)) THEN
            d1(ancil_add(irec)+i-1)=ancil_data(i)
          END IF
        END DO

        ! Ocean Nr. Surface Chlorophyll content
      ELSE IF (field == stashcode_chlorophyll) THEN
        DO i=1,p_field
          d1(ancil_add(irec)+i-1)=ancil_data(i)
        END DO

        !  3.5.5 Surface currents

      ELSE IF (field == stashcode_surf_z_curr) THEN
        DO i=1,u_field
          d1(ancil_add(irec)+i-1)=ancil_data(i)
        END DO

      ELSE IF (field == stashcode_surf_m_curr) THEN
        DO i=1,v_field
          d1(ancil_add(irec)+i-1)=ancil_data(i)
        END DO

      ELSE
        ! if we get to here, then an update request has been put in but the 
        ! code above does not know how to deal with it. This should be an
        ! immediate fatal error:
         WRITE (cmessage,'(A,I6,A)')                                  &
             ' replanca_rcf_replanca: ERROR - FIELD ',               &
             field,' omitted in field,stashcode tests'
        icode = 30
        CALL Ereport( routinename, icode, cmessage )

      END IF !End tests on FIELD numbers

      i2=i2+1

    END IF
  END DO  !  End loop over levels
END DO !  End main loop over ancil fields

IF (l_ctile) THEN
  DO i=1,p_field
    IF (sea(i) .AND. ice_fraction(i) >  0.0) THEN

      IF (l_leads_temp_prog .OR. l_amipii_ice_processing) THEN
        tstar_sice(i)=MIN(tstar_sice(i),tfs)
        tstar_ssi(i)=ice_fraction(i)*tstar_sice(i)                &
          +(1.0-ice_fraction(i))*tstar_sea(i)
      ELSE
        tstar_sea(i)=tfs
        tstar_sice(i)=(tstar_ssi(i)                               &
          -(1.0-ice_fraction(i))*tstar_sea(i))/ice_fraction(i)
      END IF

    END IF
    !
    tstar(i)=flandg(i)*tstar_land(i)                              &
      +(1.0-flandg(i))*tstar_ssi(i)
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

! Set up surface temperatures:
IF (l_ctile) THEN
  DO i=1,p_field
    tstar_land_ctile(i)=tstar_land(i)
    tstar_sea_ctile(i)=tstar_sea(i)
    ! Ensure consistency with equivalent code in
    ! replanca_rcf_replanca-rpanca1a.F90. Also helps to avoid
    ! crazy values of TSTAR_SICE in reconfigured dumps.
    ! TSTAR_SICE_CTILE(I)=TSTAR_SICE(I)
  END DO
END IF

! Deallocate the temporary storage arrays used for ancil data.
DEALLOCATE ( ancil1 )
DEALLOCATE ( ancil2 )
DEALLOCATE ( ancil_data )

! Only close file if anc_file_open is greater than 0.
IF (anc_file_open > 0) THEN
  CALL file_close (AncF_UnitNo,                                   &
                   ancil_files(anc_file_open)%filename,           &
                   LEN_TRIM(ancil_files(anc_file_open)%filename), &
                   filenameprovided,icode)

  anc_file_open = 0
END IF

9999  CONTINUE
IF (lhook) CALL dr_hook(RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE replanca_rcf_replanca
