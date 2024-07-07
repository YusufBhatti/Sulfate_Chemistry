#if defined(MCT)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis_inita2o(cmessage)
!
! Description: Initialise pointers to STASH diagnostic arrays
!              which are used when coupling with external ocean
!              and ice models. Any such arrays thus effectively become
!              prognostic variables in the context of coupled models. 
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
!
!--------------------------------------------------------------------

! Variable definitions for use with OASIS3-MCT
USE OASIS_atm_data_mod

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE atm_fields_bounds_mod
USE atm_fields_mod, ONLY: frac_land
USE mask_compression, ONLY: expand_from_mask
USE UM_ParVars
USE UM_ParParams, ONLY: halo_type_no_halo
USE atm_land_sea_mask, ONLY:  atmos_landmask_local
USE jules_sea_seaice_mod, ONLY: &
    nice, nice_use, l_sice_multilayers, l_ctile
USE submodel_mod, ONLY: atmos_im
USE oasis_timers, ONLY: l_oasis_timers, initstarts, initends

USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, len1_lookup,                            &
    len_dumphist, len_fixhd, len_tot, model_levels, mpp_len1_lookup,   &
    n_cca_lev, n_obj_d1_max, sm_levels, st_levels,                     &
    tpps_ozone_levels, tr_lbc_ukca, tr_lbc_vars, tr_ukca, tr_vars
USE cv_run_mod, ONLY: l_param_conv
USE mpp_conf_mod, ONLY: swap_field_is_scalar

USE carbon_options_mod, ONLY: l_co2_interactive
USE missing_data_mod, ONLY: imdi, rmdi  
USE errormessagelength_mod, ONLY: errormessagelength
USE umprintmgr, ONLY: newline
USE ereport_mod, ONLY: ereport
IMPLICIT NONE

!
CHARACTER(LEN=errormessagelength) :: cmessage    ! OUT - Error return message
! ----------------------------------------------------------------------
!
!  Local variables

INTEGER ::    icode       ! Error return code

INTEGER :: icode_warn     ! Local warning flag
INTEGER :: icode_abort    ! Local missing field counter/abort flag
INTEGER :: process_code   ! Processing code
INTEGER :: freq_code      ! Frequency code
INTEGER :: start_step,end_step,period ! Start, end and period step
INTEGER :: gridpt_code,weight_code ! Gridpt and weighting codes
INTEGER :: bottom_level,top_level ! Bottom and top input level
INTEGER :: grid_n,grid_s,grid_w,grid_e ! Grid corner definitions
INTEGER :: stashmacro_tag ! STASHmacro tag number
INTEGER :: field_code     ! Temporary STASH field code holder
INTEGER :: sec_code       ! Temporary STASH section code holder

INTEGER :: nlandpt        ! Dummy argument for expand_from_mask
                          ! output field.
INTEGER :: i, j, n, tc    ! Local loop indices

REAL, ALLOCATABLE :: fland(:,:)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='OASIS_INITA2O'

!----------------------------------------------------------------------
!     Set grid definition information (undefined as search is on
!     STASHmacro tag number)
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (l_oasis_timers) CALL initstarts

process_code=imdi
freq_code   =imdi
start_step  =imdi
end_step    =imdi
period      =imdi
gridpt_code =imdi
weight_code =imdi
bottom_level=imdi
top_level   =imdi
grid_n      =imdi
grid_s      =imdi
grid_e      =imdi
grid_w      =imdi

! Initialise flags for outgoing coupling fields
l_heatflux = .FALSE.
l_bluerad = .FALSE.
l_solar = .FALSE.
l_runoff = .FALSE.
l_icesheet = .FALSE.
l_w10 = .FALSE.
l_train = .FALSE.
l_tsnow = .FALSE.
l_evap2d = .FALSE.
l_taux = .FALSE.
l_tauy = .FALSE.
l_sublim(:) = .FALSE.
l_topmeltn(:) = .FALSE.
l_fcondtopn(:) = .FALSE.
l_mslp = .FALSE. 
l_pond_frac_n(:) = .FALSE.
l_pond_depth_n(:) = .FALSE.
l_tstar_sicen(:) = .FALSE.
l_pCO2 = .FALSE.
l_dust_dep = .FALSE.

! Initialise flags for incoming coupling fields
l_ocn_sst = .FALSE.
l_ocn_sstfrz = .FALSE.
l_ocn_u = .FALSE.
l_ocn_v = .FALSE.
l_ocn_freeze = .FALSE.
l_ocn_icetn(:) = .FALSE.
l_ocn_icekn(:) = .FALSE.
l_ocn_freezen(:) = .FALSE.
l_ocn_freezen_first(:) = .FALSE.
l_ocn_snowthickn(:) = .FALSE.
l_ocn_hicen(:) = .FALSE.
l_dms_conc = .FALSE.
l_CO2flux = .FALSE.
l_chloro_conc = .FALSE.

! Set up dimensions for use with incoming and outgoing
! transient arrays.
oasis_imt=lasize(1,fld_type_p,halo_type_no_halo)
oasis_jmt=lasize(2,fld_type_p,halo_type_no_halo)
oasis_jmt_u=lasize(2,fld_type_u,halo_type_no_halo)
oasis_jmt_v=lasize(2,fld_type_v,halo_type_no_halo)


! Set up fractional land points - we need a single halo for these
! because of POTENTIAL need to calculate 2nd order gradients
! of transient fields.
ALLOCATE(fland_loc(0:oasis_imt+1,0:oasis_jmt+1))

ALLOCATE(fland(oasis_imt,oasis_jmt))

CALL expand_from_mask(fland,         &
 frac_land,                                                       &
 atmos_landmask_local,                                            &
 lasize(1,fld_type_p,halo_type_no_halo)*                          &
 lasize(2,fld_type_p,halo_type_no_halo), nlandpt)

! Initialise our extended land fraction field.       
fland_loc(:,:) = 0.0
! Ensure land fraction is zero on any missing data points.
DO j=1,oasis_jmt
  DO i=1,oasis_imt
    IF (fland(I,J) /= RMDI) THEN
      fland_loc(I,J)=fland(i,j)
    END IF
  END DO
END DO

! It's important that top and bottom global domain halo rows outside
! the true model domain are treated as land in gradient calculations.
fland_loc(0:oasis_imt+1,0) = 1.0
fland_loc(0:oasis_imt+1,oasis_jmt+1) = 1.0

! Ensure our land sea mask halos are up to date.
CALL Swap_Bounds(fland_loc, oasis_imt, oasis_jmt, 1,        &
      1, 1, fld_type_p,swap_field_is_scalar)

! Process our transient list to identify those fields which need
! any special treatment
! Nothing currently required for ENDGame.

! Set up pointers for the transients we're using in this run
! DEPENDS ON: OASIS_point_translist
CALL OASIS_point_translist(nice,nice_use)

!----------------------------------------------------------------------
! Get address for each field from its STASH section/item code
! and STASHmacro tag if a diagnostic, or from its primary pointer
! if prognostic or ancillary field
! Atmosphere -> Ocean (tag=10)
!
stashmacro_tag=10

! Error handling is arranged such that an attempt is made to set up
! pointers for all fields regardless of whether any are found to be
! missing. For any missing fields, a suitable message is printed, 
! but the code does not abort until all required fields have been 
! examined. Thus, a complete list of missing fields is provided 
! without the user having to repeat the: "add missing field, attempt 
! run, add missing field, attempt run" cycle. 

icode_abort = 0
icode = 0

IF (l_taux) THEN

  sec_code = 3
  ! Get the surface X-component of windstress (U grid)
  ! The actual field we need depends on whether coastal tiling 
  ! is active or not.
  IF (l_ctile) THEN
    field_code = 392
  ELSE
    field_code = 219
  END IF

  ! DEPENDS ON: findptr
  CALL findptr(atmos_im,sec_code,field_code,                         &
    process_code,freq_code,start_step,end_step,period,               &
    gridpt_code,weight_code,                                         &
    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,              &
    stashmacro_tag,imdi,ja_taux,                                     &
    icode,cmessage)

  IF (icode /= 0 .OR. ja_taux <= 0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')           &
          "Coupling field not enabled - taux",               & 
                      sec_code, field_code
    icode_warn = -1
    CALL EReport( RoutineName, icode_warn, cmessage)
    icode_abort = icode_abort + 1 
  END IF

END IF ! taux requested


IF (l_tauy) THEN

  sec_code = 3
  ! Get the surface Y-component of windstress (V grid)
  ! The actual field we need depends on whether coastal tiling 
  ! is active or not.
  IF (l_ctile) THEN
    field_code = 394
  ELSE
    field_code = 220
  END IF

  ! DEPENDS ON: findptr
  CALL findptr(                                                      &
    atmos_im,sec_code,field_code,                                    &
    process_code,freq_code,start_step,end_step,period,               &
    gridpt_code,weight_code,                                         &
    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,              &
    stashmacro_tag,imdi,ja_tauy,                                     &
    icode,cmessage)

  IF (icode /= 0 .OR. ja_tauy <= 0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')         &
          "Coupling field not enabled - tauy",             & 
                      sec_code, field_code  
    icode_warn = -1    
    CALL EReport( RoutineName, icode_warn, cmessage)
    icode_abort = icode_abort + 1 
  END IF

END IF ! tauy requested


IF (l_heatflux) THEN
  ! Set up pointers to various components specifically invloved in
  ! the heatflux field.

  ! Net integrated downward solar on atmos grid
  sec_code = 1
  field_code = 203

  ! DEPENDS ON: findptr
  CALL findptr(                                                    &
    atmos_im,sec_code,field_code,                                  &
    process_code,freq_code,start_step,end_step,period,             &
    gridpt_code,weight_code,                                       &
    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,            &
    stashmacro_tag,imdi,ja_solar,                                  &
    icode,cmessage)

  IF (icode /= 0 .OR. ja_solar <= 0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')        &
          "Coupling field not enabled - Net solar",       & 
                      sec_code, field_code
    icode_warn = -1    
    CALL EReport( RoutineName, icode_warn, cmessage)
    icode_abort = icode_abort + 1 
  END IF

  ! Net downward longwave on atmos grid
  sec_code = 2
  field_code = 203

  ! DEPENDS ON: findptr
  CALL findptr(                                                    &
    atmos_im,sec_code,field_code,                                  &
    process_code,freq_code,start_step,end_step,period,             &
    gridpt_code,weight_code,                                       &
    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,            &
    stashmacro_tag,imdi,ja_longwave,                               &
    icode,cmessage)

  IF (icode /= 0 .OR. ja_longwave <= 0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')         &
          "Coupling field not enabled - Longwave",         & 
                     sec_code, field_code
    icode_warn=-1
    CALL EReport( RoutineName, icode_warn, cmessage)
    icode_abort = icode_abort + 1 
  END IF

  ! Sensible heat on atmos grid, area mean over open sea
  sec_code=3
  field_code=228

  ! DEPENDS ON: findptr
  CALL findptr(                                                   &
    atmos_im,sec_code,field_code,                                 &
    process_code,freq_code,start_step,end_step,period,            &
    gridpt_code,weight_code,                                      &
    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,           &
    stashmacro_tag,imdi,ja_sensible,                              &
    icode,cmessage)

  IF (icode /= 0 .OR. ja_sensible <= 0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')        &
          "Coupling field not enabled - Sensible heat",   & 
                      sec_code, field_code
    icode_warn=-1
    CALL EReport( RoutineName, icode_warn, cmessage)
    icode_abort = icode_abort + 1
  END IF

END IF   ! l_heatflux=true

IF (l_bluerad) THEN

  ! Net downward blueband solar on atmos grid
  ! Actual blue-band field used depends on coastal tiling switch.
  sec_code=1
  IF (l_ctile) THEN
    field_code = 260
  ELSE
    field_code = 204
  END IF

  ! DEPENDS ON: findptr
  CALL findptr(                                                  &
    atmos_im,sec_code,field_code,                                &
    process_code,freq_code,start_step,end_step,period,           &
    gridpt_code,weight_code,                                     &
    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,          &
    stashmacro_tag,imdi,ja_blue,                                 &
    icode,cmessage)

  IF (icode /= 0 .OR. ja_blue <= 0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')       &
          "Coupling field not enabled - Blue solar atm", & 
                      sec_code, field_code
    icode_warn=-1
    CALL EReport( RoutineName, icode_warn, cmessage)
    icode_abort = icode_abort + 1 
  END IF

END IF   ! l_bluerad=true


IF (l_evap2d) THEN

  ! Surface evaporation over sea weighted by fractional leads
  sec_code = 3
  field_code = 232

  ! DEPENDS ON: findptr
  CALL findptr(atmos_im,sec_code,field_code,                       &
    process_code,freq_code,start_step,end_step,period,             &
    gridpt_code,weight_code,                                       &
    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,            &
    stashmacro_tag,imdi,ja_evap,                                   &
    icode,cmessage)

  IF (icode /= 0 .OR. ja_evap==0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')       &
          "Coupling field not enabled - Evap over sea",  & 
                       sec_code, field_code
    icode_warn=-1
    CALL EReport( RoutineName, icode_warn, cmessage)
    icode_abort = icode_abort + 1 
  END IF

END IF


IF (l_tsnow) THEN

  ! Large-scale snowfall rate on atmos grid
  sec_code = 4
  field_code = 204

  ! DEPENDS ON: findptr
  CALL findptr(atmos_im, sec_code,field_code,            &
    process_code,freq_code,start_step,end_step,period,   &
    gridpt_code,weight_code,                             &
    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,  &
    stashmacro_tag,imdi,ja_lssnow,                       &
    icode,cmessage)

  IF (icode /= 0 .OR. ja_lssnow <= 0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')       &
          "Coupling field not enabled - LS Snow",        & 
                      sec_code, field_code
    icode_warn=-1
    CALL EReport( RoutineName, icode_warn, cmessage)
    icode_abort = icode_abort + 1 
  END IF

  IF (l_param_conv) THEN

    ! Convective snowfall rate on atmos grid
    sec_code = 5
    field_code = 206

    ! DEPENDS ON: findptr
    CALL findptr(atmos_im, sec_code,field_code,            &
      process_code,freq_code,start_step,end_step,period,   &
      gridpt_code,weight_code,                             &
      bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,  &
      stashmacro_tag,imdi,ja_cvsnow,                       &
      icode,cmessage)

    IF (icode /= 0 .OR. ja_cvsnow <= 0) THEN
      WRITE(cmessage,'(A," STASH code:",I0,",",I0)')       &
          "Coupling field not enabled - Conv Snow",        & 
                      sec_code, field_code
      icode_warn=-1
      CALL EReport( RoutineName, icode_warn, cmessage)
      icode_abort = icode_abort + 1 
    END IF

  END IF

END IF  ! Total snow


IF (l_train) THEN

  ! Large-scale rainfall rate on atmos grid
  sec_code = 4
  field_code = 203

  ! DEPENDS ON: findptr
  CALL findptr(atmos_im, sec_code,field_code,            &
    process_code,freq_code,start_step,end_step,period,   &
    gridpt_code,weight_code,                             &
    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,  &
    stashmacro_tag,imdi,ja_lsrain,                       &
    icode,cmessage)

  IF (icode /= 0 .OR. ja_lsrain <= 0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')       &
          "Coupling field not enabled - LS Rain",        & 
                      sec_code, field_code
    icode_warn=-1
    CALL EReport( RoutineName, icode_warn, cmessage)
    icode_abort = icode_abort + 1 
  END IF

  IF (l_param_conv) THEN

    ! Convective rainfall rate on atmos grid
    sec_code = 5
    field_code = 205
 
   ! DEPENDS ON: findptr
    CALL findptr(atmos_im, sec_code, field_code,           &
      process_code,freq_code,start_step,end_step,period,   &
      gridpt_code,weight_code,                             &
      bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,  &
      stashmacro_tag,imdi,ja_cvrain,                       &
      icode,cmessage)

    IF (icode /= 0 .OR. ja_cvrain <= 0) THEN
      WRITE(cmessage,'(A," STASH code:",I0,",",I0)')       &
            "Coupling field not enabled - Conv Rain",      & 
                        sec_code, field_code
       icode_warn=-1
      CALL EReport( RoutineName, icode_warn, cmessage)
      icode_abort = icode_abort + 1 
    END IF
  END IF

END IF   ! Total rain

IF (l_runoff) THEN

  ! River oputflow on atmos grid
  sec_code = 26
  field_code = 004

  ! DEPENDS ON: findptr
  CALL findptr(atmos_im, sec_code,field_code,            &
    process_code,freq_code,start_step,end_step,period,   &
    gridpt_code,weight_code,                             &
    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,  &
    stashmacro_tag,imdi,ja_riverout,                     &
    icode,cmessage)

  IF (icode /= 0 .OR. ja_riverout <= 0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')       &
          "Coupling field not enabled - River Outflow",  & 
                      sec_code, field_code
    icode_warn=-1
    CALL EReport( RoutineName, icode_warn, cmessage)
    icode_abort = icode_abort + 1 
  END IF

END IF


IF (l_w10) THEN

  ! 10m wind on atmos grid
  sec_code = 3
  field_code = 230

  ! DEPENDS ON: findptr
  CALL findptr(atmos_im, sec_code,field_code,            &
    process_code,freq_code,start_step,end_step,period,   &
    gridpt_code,weight_code,                             &
    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,  &
    stashmacro_tag,imdi,ja_w10,                          &
    icode,cmessage)

  IF (icode /= 0 .OR. ja_w10 <= 0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')       &
          "Coupling field not enabled - 10m wind ATMOS", & 
                      sec_code, field_code
    icode_warn=-1
    CALL EReport( RoutineName, icode_warn, cmessage)
    icode_abort = icode_abort + 1 
  END IF

END IF


IF (l_sublim(1)) THEN
  ! Sublimation rate
  ! Note that field 3353 is a rate(kg/m2/s) and 3231 is a total (kg/ts).
  ! Sublimation usually needs to be passed to ocean models (eg. NEMO)
  ! as rate. So if coastal tiling is not used, the sublimation total
  ! field is converted to a rate in oasis_updatecpl.
  ! (See UM7.0 or earlier swap_a2o for historical perspective).

  sec_code = 3
  IF (l_ctile) THEN
    IF (nice_use  ==  1) THEN
      field_code = 353 ! Single category field
    ELSE
      field_code = 509 ! Multi-category field
    END IF
  ELSE
    field_code = 231
  END IF

  ! DEPENDS ON: findptr
  CALL findptr(atmos_im, sec_code,field_code,            &
    process_code,freq_code,start_step,end_step,period,   &
    gridpt_code,weight_code,                             &
    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,  &
    stashmacro_tag,imdi,ja_sublim,                       &
    icode,cmessage)

  IF (icode /= 0 .OR. ja_sublim <= 0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')       &
          "Coupling field not enabled - Sublim",         &
                      sec_code, field_code
    icode_warn=-1
    CALL EReport( RoutineName, icode_warn, cmessage)
    icode_abort = icode_abort + 1 
  END IF

END IF


IF (l_fcondtopn(1)) THEN

  ! Multi category FCONDTOP
  ! NOTE: We're assuming that only multi category fields are being used!
  IF (l_sice_multilayers) THEN

    ! Multilayers coupling: fcondtop = the sea ice surface downwards
    ! heat flux (surf_ht_flux_sice in UM)
    sec_code = 3
    field_code = 510

    ! DEPENDS ON: findptr
    CALL findptr(atmos_im, sec_code,field_code,           &
      process_code,freq_code,start_step,end_step,period,  &
      gridpt_code,weight_code,                            &
      bottom_level,top_level,grid_n,grid_s,grid_w,grid_e, &
      stashmacro_tag,imdi,ja_fcondtopn,                   &
      icode,cmessage)

  ELSE

    ! 'Zero layer' coupling: fcondtop = the sea ice conductive flux
    !  through the ice (sea_ice_htf in UM, or botmelt in old
    !  UM ocean language)
    sec_code = 3
    field_code = 256

    ! DEPENDS ON: findptr
    CALL findptr(atmos_im, sec_code, field_code,          &
      process_code,freq_code,start_step,end_step,period,  &
      gridpt_code,weight_code,                            &
      bottom_level,top_level,grid_n,grid_s,grid_w,grid_e, &
      stashmacro_tag,imdi,ja_fcondtopn,                   &
      icode,cmessage)

  END IF

  IF (icode /= 0 .OR. ja_fcondtopn <= 0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')        &
          "Coupling field not enabled - FCONDTOPN",       &
                      sec_code, field_code
    icode_warn=-1
    CALL EReport( RoutineName, icode_warn, cmessage)
    icode_abort = icode_abort + 1 
  END IF

END IF


IF (l_topmeltn(1)) THEN

  ! Multi category TOPMELT
  ! NOTE: We're assuming that only multi category fields are being used!
  sec_code = 3
  field_code = 257

  ! DEPENDS ON: findptr
  CALL findptr(atmos_im, sec_code, field_code,           &
    process_code,freq_code,start_step,end_step,period,   &
    gridpt_code,weight_code,                             &
    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,  &
    stashmacro_tag,imdi,ja_topmeltn,                     &
    icode,cmessage)

  IF (icode /= 0 .OR. ja_topmeltn <= 0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')       &
          "Coupling field not enabled - TOPMELTN",       &
                      sec_code, field_code
    icode_warn=-1
    CALL EReport( RoutineName, icode_warn, cmessage)
    icode_abort = icode_abort + 1 
  END IF

END IF


IF (l_tstar_sicen(1)) THEN

  ! Sea ice surface skin temp on categories (K)
  ! NOTE: We're assuming that only multi category fields are being used!
  sec_code = 0
  field_code = 441

  ! DEPENDS ON: findptr
  CALL findptr(atmos_im, sec_code,field_code,            &
    process_code,freq_code,start_step,end_step,period,   &
    gridpt_code,weight_code,                             &
    bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,  &
    stashmacro_tag,imdi,ja_tstar_sicen,                  &
    icode,cmessage)

  IF (icode /= 0 .OR. ja_tstar_sicen <= 0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')       &
          "Coupling field not enabled - TSTAR_SICEN",    & 
                     sec_code, field_code
    icode_warn=-1
    CALL EReport( RoutineName, icode_warn, cmessage)
    icode_abort = icode_abort + 1 
  END IF

END IF


IF (l_mslp) THEN

  ! Get the Mean Sea Level Pressure
  sec_code = 16
  field_code = 222

  ! DEPENDS ON: findptr
  CALL findptr(atmos_im,sec_code,field_code,                       &
  process_code,freq_code,start_step,end_step,period,               &
  gridpt_code,weight_code,                                         &
  bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,              &
  stashmacro_tag,imdi,ja_mslp,                                     &
  icode,cmessage)

  IF (icode /= 0 .OR. ja_mslp <= 0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')                 &
          "Coupling field not enabled - mslp",                     &
                     sec_code, field_code
    icode_warn=-1
    CALL EReport( RoutineName, icode_warn, cmessage) 
    icode_abort = icode_abort + 1 
  END IF

END IF ! mslp requested

IF ( (l_pCO2) .AND. (l_co2_interactive) ) THEN

  ! Set up partial pressure of CO2
  sec_code = 0
  field_code = 252

  ! MEDUSA -> UKCA coupling: pC02 
  ! DEPENDS ON: findptr
  CALL findptr(atmos_im,sec_code,field_code,                        &
               process_code,freq_code,start_step,end_step,period,   &
               gridpt_code,weight_code,                             &
               bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,  &
               stashmacro_tag,imdi,ja_CO2,                          &
               icode,cmessage)

  IF (icode /= 0 .OR. ja_CO2 <= 0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')                  &
          "Coupling field not enabled - pCO2",                      &
                     sec_code, field_code
    icode_warn=-1
    CALL EReport( RoutineName, icode_warn, cmessage) 
    icode_abort = icode_abort + 1 
  END IF

END IF ! pCO2 required

IF (l_dust_dep) THEN

  ! Set up partial total dust deposition
  sec_code = 3
  field_code = 440

  ! MEDUSA -> UKCA coupling: total dust deposition flux
  ! DEPENDS ON: findptr
  CALL findptr(atmos_im,sec_code,field_code,                        &
               process_code,freq_code,start_step,end_step,period,   &
               gridpt_code,weight_code,                             &
               bottom_level,top_level,grid_n,grid_s,grid_w,grid_e,  &
               stashmacro_tag,imdi,ja_dust_dep,                     &
               icode,cmessage)

  IF (icode /= 0 .OR. ja_dust_dep == 0) THEN
    WRITE(cmessage,'(A," STASH code:",I0,",",I0)')                  &
          "Coupling field not enabled - dust deposition",           &
                     sec_code, field_code
    icode_warn=-1
    CALL EReport( RoutineName, icode_warn, cmessage) 
    icode_abort = icode_abort + 1 
  END IF

END IF  ! Dust required

IF (icode_abort /= 0) THEN
  ! If any of the necessary fields are missing we need to issue
  ! an appropriate error and shut down.
  ! We need a positive error code for that.
  icode=1
  WRITE(cmessage,'(A,I0,A)') newline, icode_abort,               &
        " item(s) missing from STASH coupling fields."//newline//&
        "Please refer to standard output for full details."
  CALL EReport(RoutineName, icode, cmessage)
END IF

IF (l_oasis_timers) CALL initends
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE oasis_inita2o
#endif
