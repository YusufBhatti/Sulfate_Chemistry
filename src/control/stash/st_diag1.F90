! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine ST_DIAG1 ----------------------------------------------
!
!   Purpose: Calculates STASH output diagnostics from 'dynamical'
!      fields, ie the wind fields.
!      Called at timestep 0 and at the end of timesteps.
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
!   Documentation : Unified Model Documentation Paper No P0
!               and Unified Model documentation paper No C4
!
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: stash

SUBROUTINE st_diag1( int15,                                       &
  t_incr_diagnostic,q_incr_diagnostic,qcl_incr_diagnostic,        &
                    icode,cmessage)

USE level_heights_mod
USE trignometric_mod, ONLY: sec_v_latitude, tan_v_latitude,      &
                             sec_theta_latitude
USE cderived_mod, ONLY: delta_lambda, delta_phi
USE dyn_coriolis_mod, ONLY: f3_at_v
USE jules_snow_mod, ONLY: nsmax
USE jules_sea_seaice_mod, ONLY: nice, nice_use

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE Control_Max_Sizes
USE lbc_mod
USE dump_headers_mod, ONLY: rh_deltaEW, rh_deltaNS, rh_baselat,        &
                            rh_baselong, rh_rotlat, rh_rotlong,        &
                            a_realhd
USE submodel_mod, ONLY: atmos_sm, atmos_im
USE stash_array_mod, ONLY:                                                     &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf
USE cppxref_mod, ONLY: ppx_rotate_code, ppx_elf_rotated
USE ppxlook_mod, ONLY: exppxi

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE free_tracers_inputs_mod, ONLY: a_tracer_last, a_tracer_first
USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, bl_levels, global_row_length,           &
    global_rows, land_field, len1_lookup,                              &
    len_dumphist, len_fixhd, len_tot, model_levels, mpp_len1_lookup,   &
    n_cca_lev, n_obj_d1_max, n_rows, ntiles,                           &
    river_row_length, river_rows, row_length, rows, sm_levels,         &
    theta_field_size, theta_off_size, tr_levels, tr_ukca, tr_vars,     &
    u_field_size, v_field_size

USE model_time_mod, ONLY: &
    forecast_hrs

USE atm_fields_mod, ONLY: exner_rho_levels, rho, u, v, w,  &
                          exner_theta_levels, theta

USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_lam, mt_cyclic_lam, mt_bi_cyclic_lam

USE dyn_diag_mod, ONLY: dyn_diag
IMPLICIT NONE


INTEGER ::                                                        &
        int15,                                                    &
                          ! Dummy for STASH_MAXLEN(15)
        icode              ! Out return code : 0 Normal exit
!                               !                 : >0 Error exit

CHARACTER(LEN=errormessagelength) ::                              &
        cmessage          ! Out error message if ICODE > 0

REAL, INTENT(IN) ::                                  &
  t_incr_diagnostic(row_length,rows,model_levels)    & ! dT qt_bal_cld
 ,q_incr_diagnostic(row_length,rows,model_levels)    & ! dq qt_bal_cld
 ,qcl_incr_diagnostic(row_length,rows,model_levels)    ! dqcl qt_bal_cld


CHARACTER(LEN=*) :: RoutineName
PARAMETER ( RoutineName = 'ST_DIAG1')


!  Dynamically allocated workspace for stash processing
INTEGER ::                                                        &
        ucomB_model(num_stash_levels)                             &
                                           ! 'B' grid model levs
       ,vcomB_model(num_stash_levels)                             &
                                           ! 'B' grid model levs
       ,Htheta_model(num_stash_levels)                            &
       ,Hrho_model(num_stash_levels)


REAL ::                                                           &
        stashwork(int15)                                          &
       ,ucomB_press(num_stash_levels)                             &
                                           ! 'B' grid press levs
       ,vcomB_press(num_stash_levels)                             &
                                           ! 'B' grid press levs
       ,ucomp_press(num_stash_levels)                             &
       ,vcomp_press(num_stash_levels)                             &
       ,wcomp_press(num_stash_levels)                             &
       ,testd_press(num_stash_levels)                             &
       ,testd_model(num_stash_levels)                             &
       ,pv_theta(num_stash_levels)                                &
       ,p_height(num_stash_levels)                                &
       ,theta_height(num_stash_levels)                            &
       ,rho_height(num_stash_levels)                              &
       ,w_height(num_stash_levels)                                &
       ,u_height(num_stash_levels)                                &
       ,v_height(num_stash_levels)                                &
       ,pv_press(num_stash_levels)

! Local variables

INTEGER ::                                                        &
        i,                                                        &
        ni,                                                       &
        k,                                                        &
        isl,                                                      &
        bl,                                                       &
        tl,                                                       &
        level,                                                    &
        irotu,irotv,                                              &
        ucomB_m_levs,                                             &
        vcomB_m_levs,                                             &
        ucomB_p_levs,                                             &
        vcomB_p_levs,                                             &
        ucomp_p_levs,                                             &
        Htheta_m_levs,Hrho_m_levs,                                &
        p_h_levs,theta_h_levs,rho_h_levs,                         &
        w_h_levs,u_h_levs,v_h_levs,                               &
        vcomp_p_levs,                                             &
        wcomp_p_levs,                                             &
        pv_theta_levs, pv_press_levs,                             &
        testd_p_levs,testd_m_levs                                 &
       ,im_ident                                                  &
       ,im_index                                                  &
                      !  Internal Model Index for Stash Arrays
       ,item                                                      &
                      !  STASH item no.
       ,sect          !  STASH section

PARAMETER( sect = 15 ) ! for 'dynamical' related diagnostics

INTEGER ::                                                        &
        pt002,pt003,                                              &
        pt201,pt202,pt203,pt204,pt205,pt206,pt207,pt208,pt209,    &
        pt210,pt211,pt212,pt213,pt214,pt215,pt216,pt217,pt218,    &
        pt219,pt220,pt221,pt222,pt223,pt224,pt225,pt226,pt227,    &
        pt228,pt229,pt230,                                        &
        pt231,pt232,pt233,pt234,pt235,pt236,pt237,                &
        pt238,pt239,pt240,pt241                                   &
       ,pt242,pt243,pt244,pt245,pt246                             &
       ,pt260,pt261,pt262,pt263,pt264,pt265,pt266                 &
       ,pt101,pt102,pt108,pt119                                   &
       ,pt127,pt142,pt143,pt144                                   &
       ,pt270,pt271

LOGICAL ::                                                        &
! Flags for wind rotation (lam grid)
       rot_uvcomB_press

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!  Internal Structure:

!     Set to atmosphere internal model
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
im_ident = 1
im_index = 1
icode = 0
cmessage=''


!   Section 15  Dynamics diagnostics
!
!   Local workspace definitions
!  ---------------------------------------------------------------------
!       call DYN_DIAG to calculate dynamical diagnostics and
!       call STASH to process output
!  ---------------------------------------------------------------------
!
!  This section of code contains numbers used to
!  check on the type of level which determines
!  the interpolation,
!  the codes are as follows (all for stashlist entry 11)
!  1 -- model levels
!  2 -- Pressure Levels
!  3 -- Height Levels
!  4 -- Theta Levels
!  5 -- Potential Vorticity Levels

! -------------------Extract Reqd Model levs for u on B grid (15,002)--

item = 2    ! u component on B grid, model levels
isl=stindex(1,item,sect,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  ucomB_m_levs=stash_levels(1,ni)
  DO k =1,ucomB_m_levs
    ucomB_model(k)=stash_levels(k+1,ni)  ! integer levels
  END DO
ELSE
  ucomB_m_levs=0
END IF

! -------------------Extract Reqd Model levs for v on B grid (15,003)--

item = 3    ! v component on B grid, model levels
isl=stindex(1,item,sect,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  vcomB_m_levs=stash_levels(1,ni)
  DO k =1,ucomB_m_levs
    vcomB_model(k)=stash_levels(k+1,ni)  ! integer levels
  END DO
ELSE
  vcomB_m_levs=0
END IF

! -------------------Extract Reqd Pressures for u_comB_p-------------

isl=stindex(1,201,15,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  ucomB_p_levs=stash_levels(1,ni)
  DO k =1,ucomB_p_levs
    ucomB_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  ucomB_p_levs=1
END IF

! -------------------Extract Reqd Pressures for v_comB_p-------------

isl=stindex(1,202,15,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  vcomB_p_levs=stash_levels(1,ni)
  DO k =1,vcomB_p_levs
    vcomB_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  vcomB_p_levs=1
END IF
! -------------------Extract Reqd Pressures for U_COMP_P-------------

isl=stindex(1,243,15,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  ucomp_p_levs=stash_levels(1,ni)
  DO k =1,ucomp_p_levs
    ucomp_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  ucomp_p_levs=1
END IF

! -------------------Extract Reqd Pressures for V_COMP_P-------------

isl=stindex(1,244,15,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  vcomp_p_levs=stash_levels(1,ni)
  DO k =1,vcomp_p_levs
    vcomp_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  vcomp_p_levs=1
END IF

! -------------------Extract Reqd Pressures for W_COMP_P-------------

isl=stindex(1,242,15,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  wcomp_p_levs=stash_levels(1,ni)
  DO k =1,wcomp_p_levs
    wcomp_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  wcomp_p_levs=1
END IF
! ----------Extract required thetas for Potn_vort on theta ----

isl=stindex(1,214,15,im_index)
IF (isl >  0) THEN
  IF (stlist(10,isl) <  0) THEN
    IF (stlist(11,isl) == 4) THEN
      ni = -stlist(10,isl)
      pv_theta_levs = stash_levels(1,ni)
      DO k = 1,pv_theta_levs
        pv_theta(k) = stash_levels(k+1,ni)/1000.0
        ! ***** levels are stored as integers so divide by a thousand **
      END DO
    ELSE
      cmessage =                                                  &
       ' ST_DIAG1 level not theta for pv_theta: exit routine'
      WRITE(umMessage,*) RoutineName,cmessage
      CALL umPrint(umMessage,src='st_diag1')
      icode = 1
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
  ELSE
    cmessage =                                                    &
   ' ST_DIAG1 level not a LEVELS list for PV_Theta: exit routine'
    WRITE(umMessage,*) RoutineName,cmessage
    CALL umPrint(umMessage,src='st_diag1')
    icode = 1
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
ELSE
  pv_theta_levs = 1
END IF

! ----------Extract required pressures for Potn_vort on press ----

isl=stindex(1,229,15,im_index)
IF (isl >  0) THEN
  IF (stlist(10,isl) <  0) THEN
    IF (stlist(11,isl) == 2) THEN
      ni = -stlist(10,isl)
      pv_press_levs = stash_levels(1,ni)
      DO k = 1,pv_press_levs
        pv_press(k) = stash_levels(k+1,ni)/1000.0
        ! ***** levels are stored as integers so divide by a thousand **
      END DO
    ELSE
      cmessage =                                                  &
       ' ST_DIAG1 level not pressure for pv_press: exit routine'
      WRITE(umMessage,*) RoutineName,cmessage
      CALL umPrint(umMessage,src='st_diag1')
      icode = 1
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
  ELSE
    cmessage =                                                    &
   ' ST_DIAG1 level not a LEVELS list for PV_press exit routine'
    WRITE(umMessage,*) RoutineName,cmessage
    CALL umPrint(umMessage,src='st_diag1')
    icode = 1
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
ELSE
  pv_press_levs = 1
END IF

! ----------Extract required PVs for Theta on pv -----------------

isl=stindex(1,230,15,im_index)
IF (isl >  0) THEN
  cmessage =                                                     &
' PV theta_on_pv (230,15) not enabled: exit diagnostic routine '
  WRITE(umMessage,*) RoutineName,cmessage
  CALL umPrint(umMessage,src='st_diag1')
  icode = -1
END IF

! ----------Extract required pressures for CAT_PROB_SINGLE--------------

isl=stindex(1,205,15,im_index)
IF (isl >  0) THEN
  cmessage =                                                     &
'CAT_PROB_SINGLE (205,15) not enabled: exit diagnostic routine '
  WRITE(umMessage,*) RoutineName,cmessage
  CALL umPrint(umMessage,src='st_diag1')
  icode = -1
END IF


! -------------------Extract Reqd Pressures for Test Diagnostic 233--

isl=stindex(1,233,15,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  testd_p_levs=stash_levels(1,ni)
  DO k =1,testd_p_levs
    testd_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  testd_p_levs=1
END IF

! -------------------Extract Reqd Model levs for Test Diagnostic 234--

isl=stindex(1,234,15,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  testd_m_levs=stash_levels(1,ni)
  DO k =1,testd_m_levs
    testd_model(k)=stash_levels(k+1,ni)  ! Converts to real
  END DO
ELSE
  testd_m_levs=1
END IF

!------ Extract Reqd model levels for Height on theta levels --------

isl = stindex(1,101,15,im_index)
IF ( isl  >   0 ) THEN
  ni = -stlist(10,isl)
  Htheta_m_levs = stash_levels(1,ni)
  DO k = 1, Htheta_m_levs
    Htheta_model(k) = stash_levels(k+1,ni)  ! integer levels
  END DO
ELSE
  Htheta_m_levs = 0
END IF

!------ Extract Reqd model levels for Height on rho levels --------

isl = stindex(1,102,15,im_index)
IF ( isl  >   0 ) THEN
  ni = -stlist(10,isl)
  Hrho_m_levs = stash_levels(1,ni)
  DO k = 1, Hrho_m_levs
    Hrho_model(k) = stash_levels(k+1,ni)  ! integer levels
  END DO
ELSE
  Hrho_m_levs = 0
END IF

!-------------------Extract Reqd Heights for p -------------

IF ( sf(108,15) ) THEN
  isl=stindex(1,108,15,im_index)
  IF ( isl  >   0 ) THEN
    ni       = -stlist(10,isl)
    p_h_levs = stash_levels(1,ni)
    DO k=1,p_h_levs
      p_height(k) = stash_levels(k+1,ni)/1000.0
    END DO
  ELSE
    p_h_levs = 1
  END IF
END IF

!-------------------Extract Reqd Heights for theta -------------

isl=stindex(1,119,15,im_index)
IF ( isl  >   0 ) THEN
  ni           = -stlist(10,isl)
  theta_h_levs = stash_levels(1,ni)
  DO k=1,theta_h_levs
    theta_height(k) = stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  theta_h_levs = 1
END IF


!-------------------Extract Reqd Heights for rho -------------

isl=stindex(1,127,15,im_index)
IF ( isl  >   0 ) THEN
  ni         = -stlist(10,isl)
  rho_h_levs = stash_levels(1,ni)
  DO k=1,rho_h_levs
    rho_height(k) = stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  rho_h_levs = 1
END IF

!-------------------Extract Reqd Heights for w -------------

isl=stindex(1,142,15,im_index)
IF ( isl  >   0 ) THEN
  ni       = -stlist(10,isl)
  w_h_levs = stash_levels(1,ni)
  DO k=1,w_h_levs
    w_height(k) = stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  w_h_levs = 1
END IF

!-------------------Extract Reqd Heights for u -------------

isl=stindex(1,143,15,im_index)
IF ( isl  >   0 ) THEN
  ni         = -stlist(10,isl)
  u_h_levs = stash_levels(1,ni)
  DO k=1,u_h_levs
    u_height(k) = stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  u_h_levs = 1
END IF

!-------------------Extract Reqd Heights for v -------------

isl=stindex(1,144,15,im_index)
IF ( isl  >   0 ) THEN
  ni         = -stlist(10,isl)
  v_h_levs = stash_levels(1,ni)
  DO k=1,v_h_levs
    v_height(k) = stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  v_h_levs = 1
END IF


IF (icode == 0) THEN        ! Check error code

  !----------------------------------------------------------------------
  ! Set up flags for rotating winds for a subset of lam diagnostics. This
  ! converts u,v of native equatorial lat-long lam grid back to u,v with
  ! respect to standard lat-long grid. Note that rotation of winds
  ! is not performed generically by STASH, but code needs to be
  ! introduced specifically for each (u,v) pair. The rotation flag in
  ! STASHmaster needs to be consistent and hence the flag is checked
  ! before performing rotation in dyn_diag. [Hence also, rotation can be
  ! suppressed by a STASHmaster change.]
  ! Note no explicit error trapping here.

  rot_uvcomB_press = .FALSE.               ! default

  IF (model_type == mt_lam .OR.                                   &
      model_type == mt_cyclic_lam .OR.                            &
      model_type == mt_bi_cyclic_lam) THEN

    ! winds (B grid) on pressure levels =201/202
    IF (sf(201,sect) .AND. sf(202,sect)) THEN

      irotu=exppxi(im_ident,sect,201,ppx_rotate_code,              &
                   icode,cmessage)
      irotv=exppxi(im_ident,sect,202,ppx_rotate_code,              &
                   icode,cmessage)

      IF ( irotu == ppx_elf_rotated .AND.                           &
          irotv == ppx_elf_rotated ) THEN     
        rot_uvcomB_press = .TRUE.
      END IF

    END IF     ! 201/202

  END IF ! lam only
  ! ------------------Set up Pointers for STASHWORK -------------------

  pt002=si(  2,sect,im_index)
  pt003=si(  3,sect,im_index)
  pt201=si(201,15,im_index)
  pt202=si(202,15,im_index)
  pt203=si(203,15,im_index)
  pt204=si(204,15,im_index)
  pt205=si(205,15,im_index)
  pt206=si(206,15,im_index)
  pt207=si(207,15,im_index)
  pt208=si(208,15,im_index)
  pt209=si(209,15,im_index)
  pt210=si(210,15,im_index)
  pt211=si(211,15,im_index)
  pt212=si(212,15,im_index)
  pt213=si(213,15,im_index)
  pt214=si(214,15,im_index)
  pt215=si(215,15,im_index)
  pt216=si(216,15,im_index)
  pt217=si(217,15,im_index)
  pt218=si(218,15,im_index)
  pt219=si(219,15,im_index)
  pt220=si(220,15,im_index)
  pt221=si(221,15,im_index)
  pt222=si(222,15,im_index)
  pt223=si(223,15,im_index)
  pt224=si(224,15,im_index)
  pt225=si(225,15,im_index)
  pt226=si(226,15,im_index)
  pt227=si(227,15,im_index)
  pt228=si(228,15,im_index)
  pt229=si(229,15,im_index)
  pt230=si(230,15,im_index)
  pt231=si(231,15,im_index)
  pt232=si(232,15,im_index)
  pt233=si(233,15,im_index)
  pt234=si(234,15,im_index)
  pt235=si(235,15,im_index)
  pt236=si(236,15,im_index)
  pt237=si(237,15,im_index)
  pt238=si(238,15,im_index)
  pt239=si(239,15,im_index)
  pt240=si(240,15,im_index)
  pt241=si(241,15,im_index)
  pt242=si(242,15,im_index)
  pt243=si(243,sect,im_index)
  pt244=si(244,sect,im_index)
  pt245=si(245,sect,im_index)
  pt246=si(246,sect,im_index)
  pt260=si(260,sect,im_index)
  pt261=si(261,sect,im_index)
  pt262=si(262,sect,im_index)
  pt263=si(263,sect,im_index)
  pt264=si(264,sect,im_index)
  pt265=si(265,sect,im_index)
  pt266=si(266,sect,im_index)
  pt101 = si(101,sect,im_index)
  pt102 = si(102,sect,im_index)
  pt108 = si(108,sect,im_index)
  pt119 = si(119,sect,im_index)
  pt127 = si(127,sect,im_index)
  pt142 = si(142,sect,im_index)
  pt143 = si(143,sect,im_index)
  pt144 = si(144,sect,im_index)
  pt270=si(270,sect,im_index)
  pt271=si(271,sect,im_index)

  ! Initialise STASHWORK array because DYN_DIAG does not initialise halos
  !  DIR$ CACHE_BYPASS STASHWORK
  DO i=1,int15
    stashwork(i)=0.0
  END DO

  CALL Dyn_diag(                                                    &
  ! Primary data: in
          exner_rho_levels                                                &
         ,rho,u,v,w                                                       &
         ,exner_theta_levels                                              &
        ,theta                                                            &
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
        ,pv_theta                                                         &
  ! Pressure levels for output arrays: in
        ,ucomB_press,vcomB_press                                          &
        ,ucomp_press,vcomp_press,wcomp_press                              &
        ,testd_press                                                      &
        ,p_height,theta_height,rho_height                                 &
        ,w_height,u_height,v_height                                       &
        ,pv_press                                                         &
  ! Model levels    for output arrays: in
        ,ucomB_model,vcomB_model                                          &
        ,testd_model                                                      &
        ,Htheta_model,Hrho_model                                          &
  ! Flags to request each diagnostic output field: in
  ! wind related diagnostics
        ,sf(  2,sect),sf(  3,sect)                                        &
        ,sf(201,sect),sf(202,sect)                                        &
        ,sf(243,sect),sf(244,sect),sf(242,sect)                           &
        ,sf(212,sect),sf(213,sect)                                        &
        ,sf(245,sect),sf(246,sect)                                        &
  ! PV related diagnostics
        ,sf(214,sect),sf(215,sect),sf(216,sect)                           &
        ,sf(217,sect),sf(218,sect),sf(229,sect)                           &
  ! test fields
        ,sf(231,sect),sf(232,sect),sf(233,sect),sf(234,sect)              &
  ! flux diagnostics
        ,sf(260,sect),sf(261,sect),sf(262,sect),sf(263,sect)              &
        ,sf(264,sect),sf(265,sect),sf(266,sect)                           &
  ! height and height level diagnostics
        ,sf(101,sect),sf(102,sect)                                        &
        ,sf(108,sect),sf(119,sect),sf(127,sect)                           &
        ,sf(143,sect),sf(144,sect),sf(142,sect)                           &
  ! other diagnostics
        ,sf(270,sect),sf(271,sect)                                        &
  ! Flags for wind rotation (lam grid): in
        ,rot_uvcomB_press                                                 &
  ! Diagnostics lengths: in
        ,ucomB_m_levs,vcomB_m_levs                                        &
        ,ucomB_p_levs,vcomB_p_levs                                        &
        ,ucomp_p_levs,vcomp_p_levs,wcomp_p_levs                           &
        ,pv_theta_levs, pv_press_levs                                     &
        ,testd_p_levs,testd_m_levs                                        &
        ,Htheta_m_levs,Hrho_m_levs                                        &
        ,p_h_levs,theta_h_levs,rho_h_levs,w_h_levs,u_h_levs,v_h_levs      &
  ! Diagnostic arrays: out
  ! wind related diagnostics
        ,Stashwork(pt002),Stashwork(pt003)                                &
        ,Stashwork(pt201),Stashwork(pt202)                                &
        ,Stashwork(pt243),Stashwork(pt244),Stashwork(pt242)               &
        ,Stashwork(pt212),Stashwork(pt213)                                &
        ,Stashwork(pt245),Stashwork(pt246)                                &
  ! PV related diagnostics
        ,Stashwork(pt214),Stashwork(pt215),Stashwork(pt216)               &
        ,Stashwork(pt217),Stashwork(pt218),Stashwork(pt229)               &
  ! test fields
        ,Stashwork(pt231),Stashwork(pt232)                                &
        ,Stashwork(pt233),Stashwork(pt234)                                &
  ! flux diagnostics
        ,Stashwork(pt260),Stashwork(pt261),Stashwork(pt262)               &
        ,Stashwork(pt263),Stashwork(pt264),Stashwork(pt265)               &
        ,Stashwork(pt266)                                                 &
  ! height and height level diagnostics
        ,Stashwork(pt101),Stashwork(pt102)                                &
        ,Stashwork(pt108),Stashwork(pt119),Stashwork(pt127)               &
        ,Stashwork(pt143),Stashwork(pt144),Stashwork(pt142)               &
  ! other diagnostics
        ,Stashwork(pt270),Stashwork(pt271)                                &
        )

  ! Increments from qt_bal_cld call

  item = 181    !  (181,15)  dT increment
  IF ( sf(item,sect) ) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),             &
          t_incr_diagnostic,                                        &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 182    !  (182,15)  dq increment
  IF ( sf(item,sect) ) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),             &
          q_incr_diagnostic,                                        &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)


  item = 183    !  (183,15)  dqcl increment
  IF ( sf(item,sect) ) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),             &
          qcl_incr_diagnostic,                                      &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  IF (icode /= 0) THEN
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  ! DEPENDS ON: stash
  CALL stash(atmos_sm,atmos_im,15,stashwork,                            &
             icode,cmessage)

END IF   ! Check on error code ICODE

! Error or warning exit
IF (icode /= 0) THEN

  CALL Ereport(RoutineName,icode,cmessage)
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE st_diag1
