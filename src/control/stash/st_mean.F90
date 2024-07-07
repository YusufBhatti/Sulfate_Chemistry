! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Routine: ST_MEAN --------------------------------------------------
!
!    Purpose: Extracts derived diagnostics from climate mean data in the
!             D1 array, using the special mean diagnostic sections 21-24
!             to define the required diagnostics.
!             Designed to be called from within the means subroutine.
!             This routine is closely modelled on routines ST_DIAG1 and
!             ST_DIAG2, but here the functionality is merged into a
!             single routine, although using the same underlying PCRs
!             DYN_DIAG and PHY_DIAG.
!
!
!    Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!    External documentation:
!      UM Doc Paper C0 - The top-level control system
!
!
!    Interface and arguments: ------------------------------------------
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: stash

SUBROUTINE st_mean (                                              &
             sn,intmean,                                          &
             icode,cmessage)

USE level_heights_mod
USE trignometric_mod, ONLY: sec_v_latitude, tan_v_latitude,      &
                             sec_theta_latitude
USE cderived_mod, ONLY: delta_lambda, delta_phi
USE dyn_coriolis_mod, ONLY: f3_at_v

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE atm_fields_bounds_mod, ONLY: tdims
USE atm_fields_mod, ONLY: u, v, w, rho, theta, q, qcl, qcf,       &
                          exner_rho_levels, p, p_theta_levels,    &
                          exner_theta_levels, pstar

USE dump_headers_mod, ONLY: rh_deltaEW, rh_deltaNS, rh_baselat,        &
                            rh_baselong, rh_rotlat, rh_rotlong,        &
                            a_fixhd, a_realhd
USE UM_ParVars
USE Control_Max_Sizes
USE phy_diag_mod, ONLY: phy_diag
USE calc_pmsl_inputs_mod, ONLY: npmsl_height, l_pmsl_sor
USE submodel_mod, ONLY: atmos_sm, atmos_im
USE stash_array_mod, ONLY:                                             &
    stindex, stlist, num_stash_levels, stash_levels, si, sf
USE nlsizes_namelist_mod, ONLY: bl_levels, global_row_length,          &
                                global_rows, mpp_len1_lookup, n_rows,  &
                                row_length, rows, theta_field_size,    &
                                tr_vars, u_field_size, v_field_size

USE model_time_mod, ONLY: &
    forecast_hrs
USE errormessagelength_mod, ONLY: errormessagelength

USE dyn_diag_mod, ONLY: dyn_diag
IMPLICIT NONE



INTEGER ::                                                        &
   sn                                                             &
             ! IN  - Section number for mean diagnostics
  ,intmean                                                        &
             ! IN  - Size of STASHWORK array

  ,icode     ! OUT - Return code from routine

CHARACTER(LEN=errormessagelength) :: cmessage 
             ! OUT - Return message if failure occurred
!
! ----------------------------------------------------------------------
!
!  Dynamically allocated workspace for STASH processing
!
REAL :: stashwork(intmean)

! Dummy length for _diag routines
INTEGER :: Dummy_len
PARAMETER (Dummy_len=1)

!  Pressures for DYN_DIAG
REAL ::                                                           &
        ucomp_press(num_stash_levels)                             &
       ,vcomp_press(num_stash_levels)                             &
       ,cat_prob_press(num_stash_levels)                          &
       ,pv_theta(num_stash_levels)                                &
                                    ! requested theta levels
                                    ! for pv.
       ,pv_press(num_stash_levels)                                &
                                    ! requested p levels
       ,theta_on_pv(num_stash_levels)                             &
                                      ! requested pv levels
       ,t_press(num_stash_levels)                                 &
       ,w_press(num_stash_levels)                                 &
       ,q_press(num_stash_levels)                                 &
       ,heavy_press(num_stash_levels)                             &
       ,z_press(num_stash_levels)                                 &
       ,Dummy_levels(Dummy_len)

INTEGER :: Dummy_levels_int(Dummy_len)

!  Dummy array for DYN_DIAG
REAL ::                                                           &
       pstar_old(theta_field_size)

!  Pressures for PHY_DIAG
REAL ::                                                           &
        t_p_press(num_stash_levels)                               &
       ,hts_press(num_stash_levels)                               &
       ,q_p_ice_press(num_stash_levels)                           &
       ,q_p_w_press(num_stash_levels)                             &
       ,wbpt_press(num_stash_levels)                              &
       ,th_adv_press(num_stash_levels)                            &
       ,dummy_tr_press(tr_vars+1,num_stash_levels)
!  Trig functions
REAL ::                                                           &
        nmost_lat,                                                &
        wmost_long,                                               &
        ew_space,                                                 &
        ns_space,                                                 &
        phi_pole,                                                 &
        lambda_pole,                                              &
        lat_step_inverse,                                         &
        long_step_inverse
!
!  Local variables
!
INTEGER ::                                                        &
        i,k,                                                      &
        isl,                                                      &
        ni,                                                       &
        level,                                                    &
        ucomp_p_levs,                                             &
        vcomp_p_levs,                                             &
        wcomp_p_levs,                                             &
        cat_prob_levs,                                            &
        pv_theta_levs,                                            &
        pv_press_levs,                                            &
        theta_on_pv_levs,                                         &
        dummy_tr_levs(tr_vars+1),                                 &
        tr_theta_field_size,                                      &
                             ! Dummy size for tracer fields
        hts_levs,                                                 &
                      !  no. pressure levels for heights
        t_p_levs,                                                 &
                      !  no. pressure levels for temperature
        q_p_ice_levs,                                             &
                      !  no. pressure levels for humdity wrt ice
        q_p_w_levs,                                               &
                      !  no. pressure levels for humdity wrt water
        im_index      !  Internal Model Index for stash arrays
INTEGER ::                                                        &
        pt201,pt202,pt203,pt204,pt205,pt206,pt207,pt208,pt209,    &
  pt210,pt211,pt212,pt213,pt214,pt215,pt216,pt217,pt218,pt219,    &
  pt220,pt221,pt222,pt223,pt224,pt225,pt226,pt227,pt228,pt229,    &
  pt260,pt261,pt262,pt263,pt264,pt265,pt266,                      &
  pt256                                                           &
  ,pt270,pt271

LOGICAL :: rotate_uv,rotate_max_uv                                   &
        ,sf_tracer(tr_vars+1)

LOGICAL, PARAMETER :: off = .FALSE.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ST_MEAN'
!
! NB: mean P_EXNER has been calculated in MEANCTL1
!
!     Set to atmosphere internal model
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
im_index = 1

! DYN_DIAG and PHY_DIAG do not initialise MPP halos
!  DIR$ CACHE_BYPASS STASHWORK
DO i=1,intmean
  stashwork(i)=0.0
END DO

! ----------------------------------------------------------------------
!  2. Set up levels and pointer information for call to DYN_DIAG
!     NOTE: Item nos differ from section 15.
!
isl=stindex(1,201,sn,im_index)
IF (isl >  0) THEN
  IF (stlist(10,isl) <  0) THEN
    IF (stlist(11,isl) == 2) THEN
      ni=-stlist(10,isl)
      ucomp_p_levs=stash_levels(1,ni)
      DO k=1,ucomp_p_levs
        ucomp_press(k)=stash_levels(k+1,ni)/1000.0
      END DO
    ELSE
      icode=1
      cmessage='ST_MEAN : STASH_LEVELS not pressure for U_COMP'
      GO TO 9999
    END IF
  ELSE
    cmessage='ST_MEAN : STASH_LEVELS not a LEVEL list for U_COMP'
    icode=1
    GO TO 9999
  END IF
ELSE
  ucomp_p_levs=1
END IF

isl=stindex(1,202,sn,im_index)
IF (isl >  0) THEN
  IF (stlist(10,isl) <  0) THEN
    IF (stlist(11,isl) == 2) THEN
      ni=-stlist(10,isl)
      vcomp_p_levs=stash_levels(1,ni)
      DO k=1,vcomp_p_levs
        vcomp_press(k)=stash_levels(k+1,ni)/1000.0
      END DO
    ELSE
      icode=1
      cmessage='ST_MEAN : STASH_LEVELS not pressure for V_COMP'
      GO TO 9999
    END IF
  ELSE
    cmessage='ST_MEAN : STASH_LEVELS not a LEVEL list for V_COMP'
    icode=1
    GO TO 9999
  END IF
ELSE
  vcomp_p_levs=1
END IF

pt201=si(201,sn,im_index)
pt202=si(202,sn,im_index)
pt205=si(205,sn,im_index)
pt206=si(206,sn,im_index)
pt207=si(207,sn,im_index)
pt208=si(208,sn,im_index)
pt209=si(209,sn,im_index)
pt260=si(260,sn,im_index)
pt261=si(261,sn,im_index)
pt262=si(262,sn,im_index)
pt263=si(263,sn,im_index)
pt264=si(264,sn,im_index)
pt265=si(265,sn,im_index)
pt266=si(266,sn,im_index)
pt270=si(270,sn,im_index)
pt271=si(271,sn,im_index)

! ----------------------------------------------------------------------
!  3. Call DYN_DIAG for section SN diagnostics
!
! Set flags to control wind rotation according to whether ELF grid
IF (a_fixhd(4) == 3 .OR. a_fixhd(4) == 103) THEN  ! ELF Grid
  rotate_uv=.TRUE.
  rotate_max_uv=.TRUE.
ELSE
  rotate_uv=.FALSE.
  rotate_max_uv=.FALSE.
END IF
nmost_lat=a_realhd(3)
wmost_long=a_realhd(4)
ns_space=a_realhd(2)
ew_space=a_realhd(1)
phi_pole=a_realhd(5)
lambda_pole=a_realhd(6)

!---------------------------------------------------------------------
! Beware! STASH item no.s for meaned Sections are not all identical
!         with item no.s for individual diagnostic sections (15,16) !!
!         section 15  here
!         201         201
!         202         202
!---------------------------------------------------------------------


CALL Dyn_diag(                                                    &
! Primary data: in
        exner_rho_levels                                                &
       ,rho,u,v,w                                                       &
       ,exner_theta_levels                                              &
       ,theta                                                           &
! Grid sizes and definition: in
      ,rows,n_rows,row_length,bl_levels                                 &
      ,global_rows,global_row_length                                    &
      ,theta_field_size,u_field_size,v_field_size                       &
      ,eta_theta_levels,eta_rho_levels                                  &
! Grid coordinates: in
      ,delta_lambda,delta_phi                                           &
      ,a_realhd(rh_deltaEW),a_realhd(rh_deltaNS)                        &
      ,a_realhd(rh_baselat),a_realhd(rh_baselong)                       &
      ,a_realhd(rh_rotlat),a_realhd(rh_rotlong)                         &
! Pre-calculated grid associated arrays: in
      ,r_at_u,r_at_v                                                    &
      , r_theta_levels, r_rho_levels, sec_v_latitude                    &
      , tan_v_latitude, sec_theta_latitude, f3_at_v                     &
! Time information: in
      ,forecast_hrs                                                     &
! Theta levels for output arrays
      , dummy_levels                                                    &
! Pressure levels for output arrays: in
      ,Dummy_levels,Dummy_levels                                        &
      ,ucomp_press,vcomp_press,Dummy_levels                             &
      ,Dummy_levels,Dummy_levels                                        &
      ,Dummy_levels,Dummy_levels,Dummy_levels                           &
      ,Dummy_levels,Dummy_levels,Dummy_levels                           &
! Model levels    for output arrays: in
      ,Dummy_levels_int,Dummy_levels_int                                &
      ,Dummy_levels                                                     &
      ,Dummy_levels_int,Dummy_levels_int                                &
! Flags to request each diagnostic output field: in
! wind related diagnostics
      ,off,off                                                          &
      ,sf(201,sn),sf(202,sn)                                            &
      ,off,off,off                                                      &
      ,off,off                                                          &    
      ,off,off                                                          &      
! PV related diagnostics
      ,off,off,off                                                      &
      ,off,off,off                                                      &
! test fields
       ,off,off,off,off                                                 &
! flux diagnostics
      ,sf(260,sn),sf(261,sn),sf(262,sn),sf(263,sn)                      &
      ,sf(264,sn),sf(265,sn),sf(266,sn)                                 &
! height and height level diagnostics
      ,off,off                                                          &  
      ,off,off,off                                                      &
      ,off,off,off                                                      &
! other diagnostics
      ,sf(270,sn),sf(271,sn)                                            &
! Flags for wind rotation (lam grid): in
      ,off                                                              &
! Diagnostics lengths: in
      ,Dummy_len,Dummy_len                                              &
      ,Dummy_len,Dummy_len                                              &
      ,ucomp_p_levs,vcomp_p_levs,Dummy_len                              &
      ,Dummy_len,Dummy_len                                              &
      ,Dummy_len,Dummy_len                                              &
      ,Dummy_len,Dummy_len,Dummy_len,Dummy_len                          &
      ,Dummy_len,Dummy_len,Dummy_len,Dummy_len                          &
! Diagnostic arrays: out
! wind related diagnostics
      ,Stashwork,Stashwork                                              &
      ,Stashwork(pt201),Stashwork(pt202)                                &
      ,Stashwork,Stashwork,Stashwork                                    &
      ,Stashwork,Stashwork                                              &
      ,Stashwork,Stashwork                                              &
! PV related diagnostics
      ,Stashwork,Stashwork,Stashwork                                    &
      ,Stashwork,Stashwork,Stashwork                                    &
! test fields
      ,Stashwork,Stashwork,Stashwork,Stashwork                          &
! flux diagnostics
      ,Stashwork(pt260),Stashwork(pt261),Stashwork(pt262)               &
      ,Stashwork(pt263),Stashwork(pt264),Stashwork(pt265)               &
      ,Stashwork(pt266)                                                 &
! height and height level diagnostics
      ,Stashwork,Stashwork                                              &
      ,Stashwork,Stashwork,Stashwork                                    &
      ,Stashwork,Stashwork,Stashwork                                    &
! other diagnostics
      ,Stashwork(pt270),Stashwork(pt271)                                &
      )
! ----------------------------------------------------------------------
!  5. Set up levels and pointer information for call to PHY_DIAG
!     NOTE: Item nos differ from section 16.
!
isl=stindex(1,212,sn,im_index)
IF (isl >  0) THEN
  IF (stlist(10,isl) <  0) THEN
    IF (stlist(11,isl) == 2) THEN
      ni=-stlist(10,isl)
      t_p_levs=stash_levels(1,ni)
      DO k=1,t_p_levs
        t_p_press(k)=stash_levels(k+1,ni)/1000.0
      END DO
    ELSE
      icode=1
      cmessage='ST_MEAN : STASH_LEVELS not pressure for T_P  '
      GO TO 9999
    END IF
  ELSE
    cmessage='ST_MEAN : STASH_LEVELS not a LEVEL list for T_P  '
    icode=1
    GO TO 9999
  END IF
ELSE
  t_p_levs=1
END IF

isl=stindex(1,211,sn,im_index)
IF (isl >  0) THEN
  IF (stlist(10,isl) <  0) THEN
    IF (stlist(11,isl) == 2) THEN
      ni=-stlist(10,isl)
      hts_levs=stash_levels(1,ni)
      DO k=1,hts_levs
        hts_press(k)=stash_levels(k+1,ni)/1000.0
      END DO
    ELSE
      icode=1
      cmessage='ST_MEAN : STASH_LEVELS not pressure for HTS  '
      GO TO 9999
    END IF
  ELSE
    cmessage='ST_MEAN : STASH_LEVELS not a LEVEL list for HTS  '
    icode=1
    GO TO 9999
  END IF
ELSE
  hts_levs=1
END IF

isl=stindex(1,213,sn,im_index)
IF (isl >  0) THEN
  IF (stlist(10,isl) <  0) THEN
    IF (stlist(11,isl) == 2) THEN
      ni=-stlist(10,isl)
      q_p_ice_levs=stash_levels(1,ni)
      DO k=1,q_p_ice_levs
        q_p_ice_press(k)=stash_levels(k+1,ni)/1000.0
      END DO
    ELSE
      icode=1
      cmessage='ST_MEAN : STASH_LEVELS not pressure for QPIce'
      GO TO 9999
    END IF
  ELSE
    cmessage='ST_MEAN : STASH_LEVELS not a LEVEL list for QPIce'
    icode=1
    GO TO 9999
  END IF
ELSE
  q_p_ice_levs=1
END IF
! ------------ Extract required pressures for R H wrt water ------------
isl=stindex(1,256,sn,im_index)
IF (isl >  0) THEN
  IF (stlist(10,isl) <  0) THEN
    IF (stlist(11,isl) == 2) THEN
      ni=-stlist(10,isl)
      q_p_w_levs=stash_levels(1,ni)
      DO k=1,q_p_w_levs
        q_p_w_press(k)=stash_levels(k+1,ni)/1000.0
      END DO
    ELSE
      icode=1
      cmessage='ST_MEAN : STASH_LEVELS not pressure for Q_PWater '
      GO TO 9999
    END IF
  ELSE
    cmessage='ST_MEAN : STASH_LEVELS not a LEVEL list for Q_PWat'
    icode=1
    GO TO 9999
  END IF
ELSE
  q_p_w_levs=1
END IF

DO k=1,tr_vars+1
  sf_tracer(k)=.FALSE.   ! Set tracer logical indicators to false
END DO

tr_theta_field_size=1  ! Set size for tracer arrays to 1
!
IF (sf(211,sn)) THEN
  sf(210,sn)=.TRUE. !making sure model half heights switched
END IF               !on if heights on pressure surface is reqd

pt211=si(211,sn,im_index)
pt212=si(212,sn,im_index)
pt213=si(213,sn,im_index)
pt214=si(214,sn,im_index)
pt215=si(215,sn,im_index)
pt216=si(216,sn,im_index)
pt217=si(217,sn,im_index)
pt218=si(218,sn,im_index)
pt219=si(219,sn,im_index)
pt220=si(220,sn,im_index)
pt224=si(224,sn,im_index)
pt256=si(256,sn,im_index)

! ----------------------------------------------------------------------
!  6. Call PHY_DIAG for section SN diagnostics
!

!---------------------------------------------------------------------
! Beware! STASH item no.s for meaned Sections are not identical with
!         item no.s for individual diagnostic sections (15,16) !!
!         section 16  here
!         202         211
!         203         212
!         204         213
!         222         224
!---------------------------------------------------------------------

CALL Phy_diag(                                                    &
! Primary data: in
       pstar,p,rho,u,v,w                                                &
      ,theta,q                                                          &
      ,qcl,qcf                                                          &
      ,p_theta_levels                                                   &
      ,exner_rho_levels,exner_theta_levels                              &
! Grid sizes and definition: in
      ,rows,n_rows,row_length,bl_levels                                 &
      ,theta_field_size,u_field_size,v_field_size                       &
      ,global_row_length, global_rows                                   &
      ,delta_lambda,delta_phi,sec_theta_latitude                        &
! Control information: in
      ,npmsl_height,l_pmsl_sor                                          &
! Pressure levels for output arrays: in
      ,hts_press,t_p_press                                              &
      ,q_p_ice_press,q_p_w_press,Dummy_levels                           &
! Flags to request each diagnostic output field: in
      ,off                                                              &
      ,off                                                              &
      ,sf(211,sn),sf(212,sn),sf(213,sn),off                             &
      ,off,off                                                          &
      ,sf(224,sn)                                                       &
      ,off,sf(256,sn),sf(1,20),sf(2,20),sf(28,20)                       &
! Diagnostics lengths: in
      ,hts_levs,t_p_levs,q_p_ice_levs,q_p_w_levs                        &
      ,Dummy_len                                                        &
! Diagnostic arrays: out
      ,Stashwork                                                        &
      ,Stashwork                                                        &
      ,stashwork(pt211),stashwork(pt212),stashwork(pt213)               &
      ,Stashwork                                                        &
      ,Stashwork,Stashwork                                              &
      ,stashwork(pt224)                                                 &
      ,Stashwork,stashwork(pt256)                                       &
      )

! ----------------------------------------------------------------------
!  7. Call STASH to perform processing for merged DYN_DIAG/PHY_DIAG
!     diagnostics for mean section SN
!

! DEPENDS ON: stash
CALL stash(atmos_sm,atmos_im,sn,stashwork,                              &
                                 icode,cmessage)

9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE st_mean
