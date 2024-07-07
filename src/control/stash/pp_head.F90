! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE PPHEAD------------------------------------------
!
!    Creates a 64 word PP header from the the following:-
!    1)  PP_XREF (PP cross-reference array record for this sect/item)
!    2)  FIXED length header
!    3)  INTEGER constants array
!    4)  REAL constants array
!    5)  Some input arguments
!
!    Programming standard: U M DOC  Paper 3 vn8.3
!
!    INTERFACE and ARGUMENTS:------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH

SUBROUTINE pp_head(                                               &
    im_ident,fixhd,inthd,realhd,                                  &
    len_inthd,len_realhd,ie,IS,gr,                                &
    lfullfield,real_level,pseudo_level,                           &
    samples,start,start_or_verif_time,end_or_data_time,pp_len,    &
    extraw,pp_int_head,pp_real_head,n_cols_out,num_words,         &
    len_buf_words,n_rows_out,nrow_in,srow_in,wcol_in,ecol_in,     &
    lbproc_comp,                                                  &
    sample_prd,fcst_prd,comp_accrcy,packing_type,                 &
    st_grid,iwa,zseak_rho,ck_rho,zseak_theta,ck_theta,            &
    soil_thickness,                                               & 
    model_levels,levindex,rotate,elf,                             &
    icode,cmessage)

USE model_file, ONLY: &
    mf_data_missing
USE yomhook, ONLY:  &
    lhook,          &
    dr_hook
USE parkind1, ONLY: &
    jprb,           &
    jpim
USE IAU_mod, ONLY: &
    l_iau,          &
    IAU_StartMin,   &
    IAU_EndMin
USE UM_ParVars
USE Control_Max_Sizes

USE pws_icing_pot_mod, ONLY: thickness
USE pws_dust_diag_mod, ONLY: bottom_of_layer, top_of_layer
USE c_model_id_mod, ONLY: model_id
USE lookup_addresses
USE submodel_mod, ONLY: atmos_im
USE model_id_mod, ONLY: itab
USE missing_data_mod, ONLY: imdi, rmdi

USE cppxref_mod, ONLY:                                            &
    ppx_field_code,ppx_lbvc_code,ppx_lv_code,                     &
    ppx_meto8_fieldcode,ppx_meto8_levelcode,ppx_meto8_surf,       &
    ppx_data_type,ppx_full_level,ppx_half_level,                  &   
    ppx_lb_code, ppx_lt_code, len_lbproc_c

USE nlstcall_mod, ONLY: model_analysis_mins
USE ppxlook_mod, ONLY: exppxi
USE stparam_mod, ONLY: st_riv_grid, st_uv_grid, st_cv_grid,    &
    st_zu_grid, st_mu_grid, st_cu_grid, st_zt_grid, st_scalar, &
    st_mt_grid, block_size, zonal_mean_base, field_mean_base,  &
    global_mean_base, merid_mean_base

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE stochastic_physics_run_mod, ONLY: stph_nens
USE model_time_mod, ONLY: &
    iau_dtresetstep, stepim
USE errormessagelength_mod, ONLY: errormessagelength

USE um_version_mod, ONLY: um_version_int

USE nlsizes_namelist_mod, ONLY: len_fixhd

USE packing_codes_mod, ONLY: PC_No_Packing, PC_RunLength_Packing,              &
    PC_WGDOS_Packing

IMPLICIT NONE

CHARACTER(LEN=errormessagelength) :: cmessage !OUT OUT MESSAGE FROM ROUTINE

INTEGER, INTENT(IN) ::      &
    start_or_verif_time(7), &! verif time/start time for means etc
    end_or_data_time(7),    &! data time/end time for means etc
    samples,                &! no of samples in period (timeseries)
    im_ident,               &! internal model identifier
    pp_len,                 &! length of the lookup table
    len_inthd,              &! length of the integer constants
    len_realhd,             &! length of the real constants
    fixhd(len_fixhd),       &! array of fixed constants
    inthd(len_inthd),       &! array of integer constants
    st_grid,                &! stash horizontal grid type
    model_levels,           &! no of model levels
    levindex,               &! level index
    n_rows_out,             &! pphoriz_out=n_rows_out*n_cols_out+extra
    n_cols_out,             &! pphoriz_out=n_cols_out*n_rows_out+extra
    nrow_in,srow_in,        &! the most nrthrly/southerly row.
    wcol_in,ecol_in,        &! the most westerly/easterly column
    pseudo_level,           &! output pp pseudo-level
    comp_accrcy,            &! packing accuracy in power of 2
    packing_type,           &! 0 = no packing, 1 = wgdos
    num_words,              &! number of 64 bit words to hold data
    extraw,                 &! number of extra-data words
    len_buf_words,          &! number of 64 bit words (rounded to 512)
    iwa,                    &! start word address.
    ie,                     &! item number
    IS,                     &! section number
    gr,                     &! grid point code
    lbproc_comp(len_lbproc_c)! subcomponents(0/1) to make up lbproc

LOGICAL, INTENT(IN) ::      &
    start,                  &! flag to control update for verif/start time
    lfullfield               ! TRUE if output field on full horiz domain

REAL, INTENT(IN) ::              &
    fcst_prd,                    &! forecast period
    realhd(len_realhd),          &! real header
    real_level,                  &! output pp level(real)
    sample_prd,                  &! sampling period in hours for time mean
    zseak_rho    (model_levels), &! vert coeff zsea on rho levels
    ck_rho       (model_levels), &! vert coeff ck   on rho levels
    zseak_theta(0:model_levels), &! vert coeff zsea on theta levels
    ck_theta   (0:model_levels), &! vert coeff ck   on theta levels   
    soil_thickness(model_levels)  ! soil thickness  on soil levels

INTEGER, INTENT(OUT) ::          &
    icode,                       &! return code from the routine
    pp_int_head(pp_len)           ! integer lookup table
REAL, INTENT(OUT) ::             &
    pp_real_head(pp_len)          ! Real Lookup table
!
! ---------------------------------------------------------------------
!    DEFINE LOCAL VARIABLES

CHARACTER(LEN=*) :: RoutineName
PARAMETER (   RoutineName='PP_HEAD')
INTEGER ::                      &
    pp_lbfc,                    &!  m08 level code
    pp_lbtyp,                   &!  m08 field type code
    pp_lblev,                   &!  m08 field level code
    pp_iproj,                   &!  m08 projection number
    pp_lbvc,                    &!  vertical coord type
    ii,                         &!  local counter
    int_level,                  &!  integer value of level
    k,                          &!  local counter
    ia,ib,ic,                   &!  component codes to make up lbtim
    mean_code,                  &!  spatial averaging code derived from gr
    lvcode,                     &!  lv code
    exptcode,                   &!  integer coded experiment name   
    j,                          &!  local counter 
    lev_b,                      &!  bottom deep soil level   
    lev_t                        !  top deep soil level

REAL :: part_sum                 ! partial sum for deep soil levels 
REAL :: result                   ! final sum for deep soil levels 
REAL :: bu_lev                   ! boundary upper level 
REAL :: br_lev                   ! boundary lower level       


! Local scalars:

LOGICAL ::                       &
    elf,                         &
    rotate

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!
!     Construct PP header
!
!  Timestamps ----------------------------------------------------------
!
!
!  Set up time info dependent on start flag.
!  For all but time series start will be TRUE so all time information
!  will be set up from FIXHD in effect, but for time series start
!  will be set up by TEMPORAL and passed in, so that dump headers are
!  set correctly for such fields.
!  Note: end_or_data_time will be updated from current model time in
!        FIXHD(28-34) for time means/accumulations etc.
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (start) THEN    ! start timestep so update start time
  pp_int_head(lbyr)=start_or_verif_time(1)
  pp_int_head(lbmon)=start_or_verif_time(2)
  pp_int_head(lbdat)=start_or_verif_time(3)
  pp_int_head(lbhr)=start_or_verif_time(4)
  pp_int_head(lbmin)=start_or_verif_time(5)
  pp_int_head(lbsec)=start_or_verif_time(6)
  !        PP_INT_HEAD(LBDAY)=start_or_verif_time(7)
  ! If the IAU scheme was used on the previous timestep to add a complete
  ! increment at the nominal analysis time, reset verification minute to
  ! zero for fields not in sections 0, 15 or 16. Required so we can obtain
  ! "analysis" fields from the physics modules.
  ! UM6.5 - model_analysis_hrs replaced by model_analysis_mins
  IF (l_iau .AND. im_ident == atmos_im) THEN
    IF ( IAU_StartMin == IAU_EndMin              .AND.            &
         IAU_StartMin == model_analysis_mins     .AND.            &
         stepim(atmos_im) == iau_dtresetstep + 1     .AND.        &
         IS /= 0                                 .AND.            &
         IS /= 15                                .AND.            &
         IS /= 16 ) THEN
      pp_int_head(lbmin) = 0
    END IF
  END IF
END IF
pp_int_head(lbyrd)=end_or_data_time(1)
pp_int_head(lbmond)=end_or_data_time(2)
pp_int_head(lbdatd)=end_or_data_time(3)
pp_int_head(lbhrd)=end_or_data_time(4)
pp_int_head(lbmind)=end_or_data_time(5)
pp_int_head(lbsecd)=end_or_data_time(6)
!      PP_INT_HEAD(LBDAYD)=end_or_data_time(7)
!
!  Secondary time information ------------------------------------------
!
! LBTIM is 100*IA+10*IB+IC - this encodes the time processing type
!
ia=INT(sample_prd)           ! Sampling period in whole hours
IF (sample_prd == 0.0) THEN   ! NB: may be a fraction of an hour
  ib=1                       ! Forecast field
ELSE
  IF (ia == 0) THEN
    ia=1                     ! 0 < sample_prd < 1 counts as 1 hour
  END IF
  ib=2                       ! Time mean or accumulation
END IF
ic=fixhd(8)                  ! Calendar (1: Gregorian, 2: 360 day)
!
pp_int_head(lbtim)=100*ia+10*ib+ic
pp_int_head(lbft)=fcst_prd
!
!  Data length ---------------------------------------------------------
!
pp_int_head(lblrec)=num_words
!
!  Grid code (determined from dump fixed-length header) ----------------
!
IF (samples == 0) THEN
  !       Field is not a timeseries
  ! plane Cartesian grids
  IF (fixhd(11) == 9) THEN
    pp_int_head(lbcode)=9
  ELSE ! UM default lat-lon grid indicated by fixhd(11)=1
    IF (fixhd(4) <  100) THEN
      pp_int_head(lbcode)=1  
    ELSE
      pp_int_head(lbcode)=101
    END IF
  END IF
ELSE
  !       Field is a timeseries
  pp_int_head(lbcode)=31300
  IF (fixhd(8) == 1) THEN
    !         Calendar --  1: Gregorian
    pp_int_head(lbcode)=pp_int_head(lbcode)+20
  ELSE IF (fixhd(8) == 2) THEN
    !         Calendar -- 360 day (Model Calendar)
    pp_int_head(lbcode)=pp_int_head(lbcode)+23
  ELSE
    !         Unknown calendar. Fail.
    icode=2
    cmessage='PPHEAD: unknown calender type in fixhd(8)'
  END IF
END IF
!
!  Hemispheric subregion indicator -------------------------------------
!
IF (samples >  0 .OR. .NOT. lfullfield) THEN
  !  Field is a timeseries/trajectory or subdomain of the full model area
  pp_int_head(lbhem)=3
ELSE IF (fixhd(4) <  100) THEN
  !  Otherwise, use the value for the full model area encoded in the dump
  pp_int_head(lbhem)=fixhd(4)
ELSE
  pp_int_head(lbhem)=fixhd(4)-100
END IF
!
!  Field dimensions (rows x cols) --------------------------------------
!
pp_int_head(lbrow)=n_rows_out
pp_int_head(lbnpt)=n_cols_out
!
!  'Extra data' length (now accomodates timeseries sampling data) ------
!
pp_int_head(lbext)=extraw
!
!  Packing method indicator (new definition introduced at vn2.8)--------
IF (packing_type == PC_WGDOS_Packing) THEN   ! WGDOS packing
  pp_int_head(lbpack)= PC_WGDOS_Packing
ELSE IF (packing_type == PC_RunLength_Packing) THEN! Run length encoding
  pp_int_head(lbpack)= PC_RunLength_Packing
ELSE IF (packing_type == PC_No_Packing) THEN! No packing
  pp_int_head(lbpack)= PC_No_Packing
ELSE IF (packing_type == mf_data_missing) THEN! No idea
  pp_int_head(lbpack)=mf_data_missing ! Not set code
ELSE
  icode=1
  cmessage='PPHEAD  Packing type undefined'
  pp_int_head(lbpack)= PC_No_Packing
END IF
!
!  PP header release no ------------------------------------------------
!
pp_int_head(lbrel)=3   ! for seconds in LBSEC(D), not day_number
!      PP_INT_HEAD(LBREL)=2
!
!  Primary fieldcode (some hardwiring for ELF winds) -------------------
!  Secondary fieldcode not used currently
!
pp_lbfc=exppxi(im_ident, IS, ie, ppx_field_code,                  &
               icode, cmessage)
IF (elf .AND. .NOT. rotate) THEN  ! ELF winds are in x,y direction
  IF (pp_lbfc == 56) pp_lbfc=48
  IF (pp_lbfc == 57) pp_lbfc=49
END IF
pp_int_head(lbfc)=pp_lbfc
pp_int_head(lbcfc)=0
!
!  Processing code (encodes several things in one field) ---------------
!
pp_int_head(lbproc)=0

DO ii=len_lbproc_c,1,-1
  pp_int_head(lbproc)=pp_int_head(lbproc)*2+lbproc_comp(ii)
END DO
!
!  Vertical coordinate type --------------------------------------------
!  Vertical coordinate type for reference level not coded
!
pp_lbvc=exppxi(im_ident, IS, ie, ppx_lbvc_code,                   &
               icode, cmessage)
pp_int_head(lbvc)=pp_lbvc
pp_int_head(lbrvc)=0

! [Note that, although most diagnostics are defined over int_level=
! (1:model_levels), and hence int_level=LevIndex, some variables are
! defined over int_level=(0:model_levels). For this special case,
! LevIndex is (1:model_levels+1), since LevIndex is defined starting at
! 1, and int_level != LevIndex.]
int_level = real_level+0.00001  ! ensure no rounding problems
!
!  Experiment number was coded from environment variable RUNID for
!  non-operational; set to RUN_INDIC_OP for operational use.
!  Under Rose itab is used for all model runs and expt/runid is retired.

IF (itab /= imdi) THEN
  pp_int_head(lbexp)= itab
ELSE
  ! If itab is imdi an ereport warning is printed by check_configid() 
  ! routine after configid namelist is first read in
  pp_int_head(lbexp)= 0
END IF
!
!  Direct access dataset start address and no of records ---------------
!
pp_int_head(lbegin)=iwa
pp_int_head(lbnrec)=len_buf_words
!
!  Operational fieldsfile projection no, fieldtype + level codes -------
!  These are hardwired according to model resolution
!
IF (inthd(6) == 192) THEN
  pp_iproj=802
ELSE IF (inthd(6) == 288) THEN
  pp_iproj=800
ELSE IF (inthd(6) == 96) THEN
  pp_iproj=870
ELSE IF (inthd(6) == 432) THEN
  pp_iproj=800
ELSE
  pp_iproj=900
END IF
pp_lbtyp=exppxi(im_ident, IS, ie, ppx_meto8_fieldcode,            &
               icode, cmessage)
lvcode=exppxi(im_ident, IS, ie, ppx_lv_code,                      &
             icode, cmessage)
IF (real_level == -1.0) THEN
  pp_lblev=exppxi(im_ident, IS, ie, ppx_meto8_levelcode,          &
               icode, cmessage)   ! levelcode 9999 or 8888
ELSE
  IF (im_ident  ==  atmos_im) THEN
    IF (lvcode == ppx_half_level .AND. int_level == 0) THEN
      !             This is a surface level: reset lblev
      pp_lblev=ppx_meto8_surf
    ELSE
      pp_lblev=int_level
    END IF
  ELSE
    pp_lblev=int_level
  END IF
END IF
pp_int_head(lbproj)=pp_iproj
pp_int_head(lbtyp)=pp_lbtyp
pp_int_head(lblev)=pp_lblev
!
!  Reserved slots for future expansion ---------------------------------
!
pp_int_head(lbrsvd1)=0
pp_int_head(lbrsvd2)=0
pp_int_head(lbrsvd3)=0
pp_int_head(lbrsvd4)=0
!
!  Add ensemble member to field header (if defined) --------------------
!
IF (stph_nens /= imdi) pp_int_head(lbrsvd4)=stph_nens
!
! Generate model version_id
pp_int_head(lbsrce)= um_version_int*10000 + model_id
!
! Data type - extract from PPXREF
pp_int_head(data_type)=exppxi(im_ident, IS, ie, ppx_data_type,    &
                              icode, cmessage)
!
!  Address within dump or PP file --------------------------------------
!
IF (iwa==mf_data_missing) THEN
  ! If we don't know the file location
  ! we also cannot know the data offset
  pp_int_head(naddr)=mf_data_missing
ELSE
  ! The data addr is the addr in the file
  ! from the data start location, except for fieldsfiles
  ! where is is stated to be the same as LBEGIN.
  pp_int_head(naddr)=iwa
END IF
!
!  LBUSER3 is not currently used (ie set to 0).
!
pp_int_head(lbuser3)=0
!
!  STASH section/item code ---------------------------------------------
!
pp_int_head(item_code)=IS*1000+ie
!
!  STASH pseudo-level (for fields which have pseudo-levels defined) ----
!
pp_int_head(lbplev)=pseudo_level
!
!  Spare for user's use ------------------------------------------------
!
pp_int_head(lbuser6)=0
pp_int_head(model_code) = im_ident
!
!  Reserved for future PP package use ----------------------------------
!
pp_real_head(brsvd3)=0.0
pp_real_head(brsvd4)=0.0
pp_real_head(bdatum)=0.0
pp_real_head(bacc)=comp_accrcy ! packing accuracy stored as real
!
!  Vertical grid description -------------------------------------------
!  Level and reference level
!
lev_b = exppxi(im_ident,is,ie,ppx_lb_code,icode,cmessage) 
lev_t = exppxi(im_ident,is,ie,ppx_lt_code,icode,cmessage) 

IF (pp_lbvc >= 126 .AND. pp_lbvc <= 139) THEN ! Special codes
  !                                                  (surf botttom,
  !                                                   top all zero)
  pp_real_head(blev)=0.0
  pp_real_head(bhlev)=0.0
  pp_real_head(brlev)=0.0
  pp_real_head(bhrlev)=0.0
  pp_real_head(bulev)=0.0
  pp_real_head(bhulev)=0.0
  IF (is == 20 .AND. ie == 59) THEN
    ! Set lbproc, blev and brlev for 2000-5000ft dust concentration 20/059
    ! This has LBVC=129 which is the surface; technically this is untrue but
    ! most section 20 diagnostics have this issue
    pp_real_head(brlev) = bottom_of_layer ! 2000 feet in metres
    pp_real_head(blev) = top_of_layer     ! 5000 feet in metres
    pp_int_head(lbproc) = pp_int_head(lbproc) + 2048
  END IF

ELSE IF (pp_lbvc == 9 .OR. pp_lbvc == 65) THEN ! Hybrid/ETA levels
  lvcode=exppxi(im_ident, IS, ie, ppx_lv_code,                    &
    icode, cmessage)

  ! From vn5.2:
  ! height of model level k above mean sea level is
  !       z(i,j,k) = zsea(k) + C(k)*zorog(i,j)
  ! bulev,bhulev      zsea,C of upper layer boundary
  ! blev ,bhlev       zsea,C of level
  ! brlev,bhrlev      zsea,C of lower level boundary
  ! The level here can refer to either a theta or rho level, with
  ! layer boundaries defined by surrounding rho or theta levels.
  !
  ! [Note assumption that the top of the atmosphere model,
  !  ie eta_theta_levels(model_levels) = 1.0, is equidistant from a
  !  bounding p level above (not represented explicitly) and the p
  !  level below (at eta_rho_levels(model_levels) ).]

  IF (lvcode == ppx_half_level) THEN ! theta level (& w)

    IF (int_level >= model_levels) THEN ! top level
      pp_real_head(bulev) =  zseak_theta(model_levels) * 2.0      &
                                   - zseak_rho(model_levels)
      pp_real_head(bhulev)=     Ck_theta(model_levels) * 2.0      &
                                   -    Ck_rho(model_levels)
    ELSE
      pp_real_head(bulev) = zseak_rho(int_level+1)
      pp_real_head(bhulev)=    Ck_rho(int_level+1)
    END IF                             ! top level

    pp_real_head(blev) = zseak_theta(int_level)
    pp_real_head(bhlev)=    Ck_theta(int_level)

    ! Note that the lowest theta level (=1) has a lower layer boundary at
    ! the surface, as set explicitly in the model interface to physics.
    IF (int_level <= 1) THEN            ! bottom level
      pp_real_head(brlev) = 0.0     ! zsea at/below surface
      pp_real_head(bhrlev)= 1.0     ! C    at/below surface
    ELSE
      pp_real_head(brlev) = zseak_rho(int_level)
      pp_real_head(bhrlev)=    Ck_rho(int_level)
    END IF                              ! bottom level

  ELSE IF (lvcode == ppx_full_level) THEN ! rho level (& u,v,p)

    IF (int_level >  model_levels) THEN ! p above top level
      pp_real_head(bulev) = zseak_theta(model_levels) * 2.0       &
                                  - zseak_rho(model_levels)
      pp_real_head(bhulev)=    Ck_theta(model_levels) * 2.0       &
                                   -   Ck_rho(model_levels)
      pp_real_head(blev) = pp_real_head(bulev)
      pp_real_head(bhlev)= pp_real_head(bhulev)
    ELSE
      pp_real_head(bulev) = zseak_theta(int_level)
      pp_real_head(bhulev)=    Ck_theta(int_level)
      pp_real_head(blev)  = zseak_rho(int_level)
      pp_real_head(bhlev) =    Ck_rho(int_level)
    END IF                              ! p above top level

    IF (int_level <= 0) THEN            ! bottom level
      pp_real_head(brlev) = 0.0    ! zsea at/below surface
      pp_real_head(bhrlev)= 1.0    ! C    at/below surface
    ELSE
      pp_real_head(brlev) = zseak_theta(int_level-1)
      pp_real_head(bhrlev)=    Ck_theta(int_level-1)
    END IF                              ! bottom level

  ELSE                ! Illegal lvcode
    icode=1
    cmessage=' Inconsistent vertical coordinate codes in '//&
        'STASHmaster for this output field: LevelT indicates '//&
        'that this variable is on model levels, but LBVC used '//&
        'for pp header label does not.'
    WRITE(umMessage,'(A,A,A,4I6)') RoutineName,cmessage,  &
            ' section,item,LevelT,pp_lbvc=',IS,ie,lvcode,pp_lbvc
    CALL umPrint(umMessage,src='pp_head')

  END IF               ! Test on lvcode
ELSE IF(pp_lbvc == 6 .AND. lev_b == 8 .AND. lev_t == 9) THEN 
! Deep soil levels. BLEV should be set to the mid-point of the layer  
  result = 0.0   
  part_sum = 0.0 
  bu_lev = 0.0 
  br_lev = 0.0 
  IF (int_level == 1) THEN ! first level   
    result = 0.5 * soil_thickness(int_level) 
    br_lev = soil_thickness(int_level) 
  ELSE   
    DO j = 1, int_level - 1   
      part_sum = part_sum + soil_thickness(j)  
      bu_lev = bu_lev + soil_thickness(j) 
      br_lev = bu_lev + soil_thickness(j + 1) 
    END DO   
    result = part_sum + 0.5 * soil_thickness(int_level) 
  END IF   
  pp_real_head(blev)= result 
  pp_real_head(bulev)= bu_lev 
  pp_real_head(brlev)= br_lev  
  pp_real_head(bhlev)=0.0   
  pp_real_head(bhrlev)=0.0       
  pp_real_head(bhulev)=0.0 
 
ELSE IF(pp_lbvc == 8 .AND. IS == 20) THEN
! Special height levels for PWS diagnostics section, e.g. Thickness
  IF (ie == 1 .OR. ie == 2) THEN
    pp_real_head(brlev)= 1000.0  ! 1000-nnn mb Thickness
    pp_int_head(lbproc)=pp_int_head(lbproc) + 1024
  ELSE
    pp_real_head(brlev)= 0.0
  END IF
  IF (ie == 1) THEN
    pp_real_head(blev)= 500.0  ! 1000-500 mb Thickness
  ELSE IF (ie == 2) THEN
    pp_real_head(blev)= 850.0  ! 1000-850 mb Thickness
  ELSE 
    pp_real_head(blev)= real_level 
  END IF
  ! Set BRLEV for WAFC Icing Potential. Should be the upper boundary (closest
  ! to the surface) of the pressure level in hPa, and the real_level is 
  ! central in the region, which has width of "thickness".
  ! Hence this becomes:
  ! real_level (hPa) + thickness (in Pa) / 100.0 (convert to hPa) / 2.0 (half
  ! thickness of layer)
  IF (ie == 43) THEN
    pp_real_head(brlev) = real_level + thickness / 200.0
  END IF
  pp_real_head(bhlev)=0.0
  pp_real_head(bhrlev)=0.0
  pp_real_head(bulev)=0.0
  pp_real_head(bhulev)=0.0
ELSE IF(is == 20 .AND. ie == 12) THEN
! Convective cloud depth has BLEV=0.0
  pp_real_head(blev)=0.0
  pp_real_head(bhlev)=0.0 
  pp_real_head(brlev)=0.0  
  pp_real_head(bhrlev)=0.0 
  pp_real_head(bulev)=0.0  
  pp_real_head(bhulev)=0.0 
ELSE IF(is == 20 .AND. ie == 3) THEN
! Wind at 10 metres
  pp_real_head(blev)=10.0
  pp_real_head(bhlev)=0.0 
  pp_real_head(brlev)=0.0  
  pp_real_head(bhrlev)=0.0 
  pp_real_head(bulev)=0.0  
  pp_real_head(bhulev)=0.0 
ELSE 
  pp_real_head(blev)=real_level 
  pp_real_head(bhlev)=0.0 
  pp_real_head(brlev)=0.0  ! The boundary levels 
  pp_real_head(bhrlev)=0.0 ! are not known 
  pp_real_head(bulev)=0.0  ! for pressure 
  pp_real_head(bhulev)=0.0 ! levels. 
END IF 
! 
!  Horizontal grid description -----------------------------------------
!  Position of pole (from dump fixed-length header)
!  Grid orientation (hardwired 0.0)
!  Origin and spacing of grid (depends on output grid type)
!
pp_real_head(bplat)=realhd(5)
pp_real_head(bplon)=realhd(6)
pp_real_head(bgor)=0.0
IF (samples >  0) THEN   ! Indicates a timeseries/trajectory
  pp_real_head(bzx)=0.0
  pp_real_head(bdx)=0.0
  pp_real_head(bzy)=0.0
  pp_real_head(bdy)=0.0
ELSE
  IF (st_grid == st_riv_grid) THEN
    pp_real_head(bdy) = 180.0/pp_int_head(lbrow)
    pp_real_head(bzy) = realhd(3) - pp_real_head(bdy)*0.5
    pp_real_head(bdx) = 360.0/pp_int_head(lbnpt)
    pp_real_head(bzx) = realhd(4) - pp_real_head(bdx)*0.5
  ELSE

    IF (st_grid == st_uv_grid .OR. st_grid == st_cv_grid .OR.        &
        st_grid == st_zu_grid .OR. st_grid == st_mu_grid) THEN
      pp_real_head(bzy)=realhd(3)-realhd(2)    ! V pts
    ELSE
      pp_real_head(bzy)=realhd(3)-realhd(2)/2.0   ! Zeroth Lat BZY
    END IF

    IF (st_grid == st_uv_grid .OR. st_grid == st_cu_grid .OR.        &
        st_grid == st_zu_grid .OR. st_grid == st_mu_grid) THEN
      pp_real_head(bzx)=realhd(4)-realhd(1)   ! U points
    ELSE
      pp_real_head(bzx)=realhd(4)-realhd(1)/2.0   ! Zeroth Long BZX
    END IF
    
    pp_real_head(bdx)=realhd(1) ! Long intvl BDX
    pp_real_head(bdy)=realhd(2) ! Lat intvl BDY
  END IF
  !
  ! Add on offset for fields not starting from the origin (sub-areas)
  !
        ! If this is variable horizontal grid information
        ! then the horizontal start points are NOT equally spaced.
        ! The existing code is meaningless in this case so
        ! we must get our actual start point some other way!
  pp_real_head(bzy)=pp_real_head(bzy)                           &
      +(srow_in-1)*pp_real_head(bdy)
  pp_real_head(bzx)=pp_real_head(bzx)                           &
      +(wcol_in-1)*pp_real_head(bdx)
  IF (pp_real_head(bzx) >= 360.0) THEN
    pp_real_head(bzx)=pp_real_head(bzx)-360.0
  END IF

  !
  ! If horizontal averaging has been applied to the output field,
  ! set BDX and/or BDY to the full (sub)domain extent which was processed.
  ! If the input field was intrinsically non-2D (eg. zonal), assume that
  ! the collapsed dimension(s) covered the full model domain.
  !
  mean_code=(gr/block_size)*block_size
  IF (st_grid == st_zt_grid .OR. st_grid == st_zu_grid            &
      .OR. st_grid == st_scalar) THEN
    pp_real_head(bdx)=REAL(inthd(6))*pp_real_head(bdx)
  ELSE IF (mean_code == zonal_mean_base .OR.                       &
      mean_code == field_mean_base .OR.                           &
      mean_code == global_mean_base) THEN
    pp_real_head(bdx)=ABS(REAL(ecol_in-wcol_in))*pp_real_head(bdx)
  END IF
  !
  IF (st_grid == st_mt_grid .OR. st_grid == st_mu_grid            &
      .OR. st_grid == st_scalar) THEN
    pp_real_head(bdy)=REAL(inthd(7))*pp_real_head(bdy)
  ELSE IF (mean_code == merid_mean_base .OR.                       &
      mean_code == field_mean_base .OR.                           &
      mean_code == global_mean_base) THEN
    pp_real_head(bdy)=ABS(REAL(nrow_in-srow_in))*pp_real_head(bdy)

  END IF
END IF
IF ( (fixhd(117) > 1 .AND. fixhd(122) > 1)) THEN
  pp_real_head(bdx) = rmdi
  pp_real_head(bdy) = rmdi
  pp_real_head(bzx) = rmdi
  pp_real_head(bzy) = rmdi
END IF
!
! Missing data indicator (from PARAMETER) ------------------------------
! MKS scaling factor (unity as model uses SI units throughout)
!
pp_real_head(bmdi)=rmdi
pp_real_head(bmks)=1.0
!

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pp_head
