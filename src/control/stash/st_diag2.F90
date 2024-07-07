! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!   Purpose : To provide the interface for PHY_DIAG
!
!
!   Programming standard; Unified Model Documentation Paper No. 3
!                         version no. 1, dated 15/01/90
!
!   Logical components covered : D4
!
!   System task : P0
!
!   Documentation : Unified Model Documentation Paper No P0
!
!
!   Arguments
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: stash

SUBROUTINE st_diag2( int16,                                       & 
                      icode,cmessage)

USE atm_fields_bounds_mod

USE trignometric_mod, ONLY: sec_theta_latitude
USE cderived_mod, ONLY: delta_lambda, delta_phi
USE jules_snow_mod, ONLY: nsmax
USE jules_sea_seaice_mod, ONLY: nice, nice_use
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE Control_Max_Sizes
USE lbc_mod
USE phy_diag_mod, ONLY: phy_diag
USE calc_pmsl_inputs_mod, ONLY: npmsl_height, l_pmsl_sor
USE dump_headers_mod, ONLY: a_realhd
USE submodel_mod, ONLY: atmos_sm, atmos_im
USE stash_array_mod, ONLY:                                             &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf
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

USE atm_step_local, ONLY: t_inc_pres, q_inc_pres, qcl_inc_pres,        &
                          qcf_inc_pres, cf_inc_pres, cfl_inc_pres,     &
                          cff_inc_pres, t_dini, q_dini, qcl_dini,      &
                          qcf_dini, cf_dini, cfl_dini, cff_dini

USE atm_fields_mod, ONLY: exner_rho_levels, exner_theta_levels, &
    p_theta_levels, p, rho, u, v, w, theta, q, qcl, qcf, pstar

USE cloud_inputs_mod, ONLY: i_cld_vn
USE pc2_constants_mod, ONLY: i_cld_pc2

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER ::                                                        &
        int16                                                     &
                          ! STASHWORK size = STASH_MAXLEN(16)
       ,icode             ! Out return code : 0 Normal exit
                          !                 :>0 Error exit

CHARACTER(LEN=errormessagelength) ::                              &
        cmessage          ! Out error message if ICODE > 0

!  Locally dynamically allocated work area

REAL ::                                                           &
        stashwork(int16),                                         &
        t_p_press(num_stash_levels),                              &
        hts_press(num_stash_levels),                              &
        RHice_press(num_stash_levels),                            &
        RHwat_press(num_stash_levels),                            &
        wbpt_press(num_stash_levels),                             &
        th_adv_press(num_stash_levels),                           &
        press_levs(num_stash_levels),                             &
        tr_press(tr_vars+1,num_stash_levels),                     &
        T_m(theta_field_size,model_levels),                       &
                                        ! T on model levels
        hts_theta(theta_field_size,model_levels),                 &
                                        ! hts on theta levels
        hts_rho(theta_field_size,model_levels),                   &
                                        ! hts on rho levels
        qc(theta_field_size,model_levels),                        &
                                        ! cloud water content
        qT(theta_field_size,model_levels)
                                        ! total specific humidity

! Local variables

INTEGER ::                                                        &
        i,j,                                                      &
        ni,                                                       &
        k,                                                        &
        isl,                                                      &
        bl,                                                       &
        tl,                                                       &
        level,                                                    &
        last_point,                                               &
        first_point                                               &
       ,itr                                                       &
       ,stash_tr_first                                            &
       ,stash_tr_last                                             &
       ,im_ident                                                  &
                      !  Internal Model Identifier
       ,im_index                                                  &
                      !  Internal Model Index for Stash arrays
      ,item                                                       &
                      ! STASH item
      ,sect           ! STASH section
PARAMETER( sect = 16 ) ! for 'physics' end of timestep diagnostics

INTEGER ::                                                        &
        t_p_levs,                                                 &
        hts_levs, h2_p_levs,                                      &
        RHice_levs,RHwat_levs,                                    &
        wbpt_levs,                                                &
        th_adv_p_levs                                             &
       ,tr_press_levs(tr_vars+1)                                  &
       ,tr_theta_field_size   ! size for DA of tracer fields

INTEGER ::                                                        &
        h2_ind(num_stash_levels)

INTEGER ::                                                        &
        pt_tracer(tr_vars+1)
INTEGER ::                                                        &
        pt201,pt202,pt203,pt204,pt205,pt206,pt207,pt208,pt209,    &
  pt210,pt211,pt212,pt213,pt214,pt215,pt216,pt217,pt218,pt219,    &
  pt220,pt221,pt222,pt223,pt224,pt225,pt255,pt256

REAL ::                                                           &
        ew_space,                                                 &
        ns_space

LOGICAL ::                                                        &
        sf_tracer(tr_vars+1)

! 3d work array for calculating separate +/- PC2 increments
! only to be allocated if it is needed.
REAL, ALLOCATABLE ::  work3d(:,:,:)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ST_DIAG2'

!  Initialisation
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
first_point = 1
last_point  = theta_field_size

!  Internal Structure:

!     Set to atmosphere internal model
im_ident = 1
im_index = 1

! Allocate work3 array if calculating +/- increments for cfl,cff,qcl,qcf.
IF (sf(140,sect) .OR. sf(141,sect) .OR.                           &
    sf(142,sect) .OR. sf(143,sect) .OR.                           &
    sf(146,sect) .OR. sf(147,sect) .OR.                           &
    sf(148,sect) .OR. sf(149,sect) .OR.                           &
    sf(150,sect) .OR. sf(151,sect) .OR.                           &
    sf(156,sect) .OR. sf(157,sect) ) THEN
  ALLOCATE ( work3d(row_length,rows,model_levels) )
END IF

! ----- Calculate additional diagnostic quantities------------------

! ------------ Extract required pressures for T_P ----------------------
ns_space=a_realhd(2)
ew_space=a_realhd(1)

isl=stindex(1,203,16,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  t_p_levs=stash_levels(1,ni)
  DO k =1,t_p_levs
    t_p_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  t_p_levs=1
END IF

! ------------ Extract required pressures for Heights ------------------

isl=stindex(1,202,16,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  hts_levs=stash_levels(1,ni)
  DO k =1,hts_levs
    hts_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  hts_levs=1
END IF


! ------------ Extract required pressures for R H wrt ice --------------

isl=stindex(1,204,16,im_index)
IF (isl >  0) THEN
  IF (stlist(10,isl) <  0) THEN
    IF (stlist(11,isl) == 2) THEN
      ni=-stlist(10,isl)
      RHice_levs=stash_levels(1,ni)
      DO k =1,RHIce_levs
        RHice_press(k)=stash_levels(k+1,ni)/1000.0
      END DO
    ELSE
      cmessage='ST_DIAG2 : Level not pressure for RHice'
      icode=1
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
  ELSE
    cmessage='ST_DIAG2 : Level not a levels list for RHice'
    icode=1
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
ELSE
  RHice_Levs=1
END IF

! ------------ Extract required pressures for R H wrt water ------------
isl=stindex(1,256,16,im_index)
IF (isl >  0) THEN
  IF (stlist(10,isl) <  0) THEN
    IF (stlist(11,isl) == 2) THEN
      ni=-stlist(10,isl)
      RHwat_levs=stash_levels(1,ni)
      DO k =1,RHwat_levs
        RHwat_press(k)=stash_levels(k+1,ni)/1000.0
      END DO
    ELSE
      cmessage='ST_DIAG2 : Level not pressure for RHwat'
      icode=1
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
  ELSE
    cmessage='ST_DIAG2 : Level not a levels list for RHwat'
    icode=1
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
ELSE
  RHwat_levs=1
END IF
! ------------ Extract required pressures for Wet bulb pot temp---------

isl=stindex(1,205,16,im_index)
IF (isl >  0) THEN
  IF (stlist(10,isl) <  0) THEN
    IF (stlist(11,isl) == 2) THEN
      ni=-stlist(10,isl)
      wbpt_levs=stash_levels(1,ni)
      DO k =1,wbpt_levs
        wbpt_press(k)=stash_levels(k+1,ni)/1000.0
      END DO
    ELSE
      cmessage='ST_DIAG2 : Level not pressure for WBPT'
      icode=1
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
  ELSE
    cmessage='ST_DIAG2 : Level not a levels list for WBPT'
    icode=1
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
ELSE
  wbpt_levs=1
END IF
! ------------ Extract required pressures for Thermal advection --------

isl=stindex(1,219,16,im_index)
IF (isl >  0) THEN
  IF (stlist(10,isl) <  0) THEN
    IF (stlist(11,isl) == 2) THEN
      ni=-stlist(10,isl)
      th_adv_p_levs=stash_levels(1,ni)
      DO k =1,th_adv_p_levs
        th_adv_press(k)=stash_levels(k+1,ni)/1000.0
      END DO
    ELSE
      cmessage='ST_DIAG2 : Level not pressure for THADV'
      icode=1
      IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
  ELSE
    cmessage='ST_DIAG2 : Level not a levels list for THADV'
    icode=1
    IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
ELSE
  th_adv_p_levs=1
END IF

! ------------ Extract required pressures for Tracers ------------------

stash_tr_first=226
stash_tr_last=254
!     ITR is a count of tracers found to be using this diagnostic
itr=0
IF (tr_vars >  0) THEN

  !     Initialize PT_TRACER and SF_TRACER arrays
  DO i=1,tr_vars
    pt_tracer(i)=1
    sf_tracer(i)=.FALSE.
  END DO

  DO j=stash_tr_first,stash_tr_last
    isl=stindex(1,j,16,im_index)
    IF (isl >  0) THEN
      IF (stlist(10,isl) <  0) THEN
        IF (stlist(11,isl) == 2) THEN
          itr=itr+1
          ni=-stlist(10,isl)
          tr_press_levs(itr)=stash_levels(1,ni)
          DO k=1,tr_press_levs(itr)
            tr_press(itr,k)=stash_levels(k+1,ni)/1000.0
          END DO
          pt_tracer(itr)=si(j,16,im_index)
          sf_tracer(itr)=sf(j,16)
        ELSE
          cmessage='ST_DIAG2 : Level not pressure for Tracers'
          icode=1
          IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
          RETURN
        END IF
      ELSE
        cmessage='ST_DIAG2 : Level not a levels list for Tracers'
        icode=1
        IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
        RETURN
      END IF
    END IF
  END DO
END IF

!     Set last (or only) values in tracer pointer arrays
pt_tracer(tr_vars+1)=1
sf_tracer(tr_vars+1)=.FALSE.

!     Set size of TR_P_FIELD_DA depending on whether any tracers
!     are using this diagnostic. (Used for dynamic allocation of
!     tracer arrays in phy_diag).
IF (itr >  0) THEN
  tr_theta_field_size=theta_field_size
ELSE
  tr_theta_field_size=1

END IF

! ----- check height available for calculation of height**2 -----------


IF (sf(224,16)) THEN
  IF (.NOT. sf(202,16)) THEN
    cmessage='ST_DIAG2 : ERROR h**2 requires H at same timestep'
    icode=1
    GO TO 9999
  ELSE
    isl=stindex(1,224,16,im_index)
    IF (isl >  0) THEN
      ni=-stlist(10,isl)
      h2_p_levs=stash_levels(1,ni)
      DO k =1,h2_p_levs
        press_levs(k)=stash_levels(k+1,ni)/1000.0
        DO i=1,hts_levs
          IF (press_levs(k) == hts_press(i)) THEN
            h2_ind(k)=i
          END IF
        END DO
      END DO
    ELSE
      h2_p_levs=1
    END IF
  END IF
ELSE
  h2_p_levs=1
END IF


pt201=si(201,16,im_index)
pt202=si(202,16,im_index)
pt203=si(203,16,im_index)
pt204=si(204,16,im_index)
pt205=si(205,16,im_index)
pt206=si(206,16,im_index)
pt207=si(207,16,im_index)
pt208=si(208,16,im_index)
pt209=si(209,16,im_index)
pt210=si(210,16,im_index)
pt211=si(211,16,im_index)
pt212=si(212,16,im_index)
pt213=si(213,16,im_index)
pt214=si(214,16,im_index)
pt215=si(215,16,im_index)
pt216=si(216,16,im_index)
pt217=si(217,16,im_index)
pt218=si(218,16,im_index)
pt219=si(219,16,im_index)
pt220=si(220,16,im_index)
pt221=si(221,16,im_index)
pt222=si(222,16,im_index)
pt223=si(223,16,im_index)
pt224=si(224,16,im_index)
pt225=si(225,16,im_index)
pt255=si(255,16,im_index)
pt256=si(256,16,im_index)

! Initialise STASHWORK as PHY_DIAG avoids MPP halo calculations
!  DIR$ CACHE_BYPASS STASHWORK
DO i=1,int16
  stashwork(i)=0.0
END DO

CALL Phy_diag(                                                    &
! Primary data: in
       pstar,p,rho,u,v,w                                                &
      ,theta,q,qcl,qcf                                                  &
      ,p_theta_levels                                                   &
      ,exner_rho_levels,exner_theta_levels                              &
! Grid sizes and definition: in
      ,rows,n_rows,row_length,bl_levels                                 &
      ,theta_field_size,u_field_size,v_field_size                       &
      ,global_row_length, global_rows                                   &
      ,delta_lambda,delta_phi,sec_theta_latitude                        &
! Control information: in
      ,npmsl_height,L_pmsl_sor                                          &
! Pressure levels for output arrays: in
      ,hts_press,t_p_press                                              &
      ,RHice_press,RHwat_press,wbpt_press                               &
! Flags to request each diagnostic output field: in
      ,sf(  4,16),sf(201,16),sf(202,16),sf(203,16),sf(204,16)           &
      ,sf(205,16),sf(206,16),sf(207,16),sf(222,16),sf(255,16)           &
      ,sf(256,16),sf(1,20),sf(2,20),sf(28,20)                           &
! Diagnostics lengths: in
      ,hts_levs,t_p_levs,RHice_levs,RHwat_levs,wbpt_levs                &
! Diagnostic arrays: out
      ,T_m                                                              &
      ,hts_theta,stashwork(pt202),stashwork(pt203)                      &
      ,stashwork(pt204),stashwork(pt205),qc,qT,stashwork(pt222)         &
      ,hts_rho,stashwork(pt256)                                         &
      )



item = 4      !  (4,16) is T on model levels
IF ( sf(item,sect) ) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3D(stashwork(si(item,sect,im_index)),T_m,         &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        im_ident,sect,item,                                       &
        icode,cmessage)

END IF   ! SF(  item,sect)

! These physics increments will only have been set up and calculated if
! the PC2 cloud scheme is run
IF (i_cld_vn == i_cld_pc2) THEN

  ! PC2 cloud scheme diagnostics for condensation due to vertical
  ! advection

  item = 181      !  (181,16)  Temperature change from condensation
  IF ( sf(item,sect) ) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),t_inc_pres,  &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 182      !  (182,16)  Vapour change from condensation
  IF ( sf(item,sect) ) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),q_inc_pres,  &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 183      !  (183,16)  Liquid change from condensation
  IF ( sf(item,sect) ) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),qcl_inc_pres,&
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 150      !  (150,16)  Liquid change condensation positive
  IF ( sf(item,sect) ) THEN
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work3d(i,j,k) = MAX(0.0,qcl_inc_pres(i,j,k))
        END DO
      END DO
    END DO
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),work3d,      &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 151      !  (151,16)  Liquid change condensation negative
  IF ( sf(item,sect) ) THEN
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work3d(i,j,k) = MIN(0.0,qcl_inc_pres(i,j,k))
        END DO
      END DO
    END DO
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),work3d,      &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 184      !  (184,16)  Ice change from condensation
  IF ( sf(item,sect) ) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),qcf_inc_pres,&
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 192      !  (192,16)  Total cloud change from condensation
  IF ( sf(item,sect) ) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),cf_inc_pres, &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 193      !  (193,16)  Liquid cloud change from condensation
  IF ( sf(item,sect) ) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),cfl_inc_pres,&
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 156      !  (156,16)  Liquid cloud change cond'n positive
  IF ( sf(item,sect) ) THEN
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work3d(i,j,k) = MAX(0.0,cfl_inc_pres(i,j,k))
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),work3d,      &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 157      !  (157,16)  Liquid cloud change cond'n positive
  IF ( sf(item,sect) ) THEN
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work3d(i,j,k) = MIN(0.0,cfl_inc_pres(i,j,k))
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),work3d,      &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 194      !  (194,16)  Ice cloud change from condensation
  IF ( sf(item,sect) ) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),cff_inc_pres,&
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  ! PC2 cloud scheme diagnostics for condensation due to initialisation

  item = 161      !  (161,16)  Temperature change from initialisation
  IF ( sf(item,sect) ) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),t_dini,      &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 162      !  (162,16)  Vapour change from initialisation
  IF ( sf(item,sect) ) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),q_dini,      &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 163      !  (163,16)  Liquid change from initialisation
  IF ( sf(item,sect) ) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),qcl_dini,    &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 140      !  (140,16)  Liquid change initialisation positive
  IF ( sf(item,sect) ) THEN
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work3d(i,j,k) = MAX(0.0,qcl_dini(i,j,k))
        END DO
      END DO
    END DO
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),work3d,      &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 141      !  (141,16)  Liquid change initialisation negative
  IF ( sf(item,sect) ) THEN
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work3d(i,j,k) = MIN(0.0,qcl_dini(i,j,k))
        END DO
      END DO
    END DO
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),work3d,      &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 164      !  (164,16)  Ice change from initialisation
  IF ( sf(item,sect) ) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),qcf_dini,    &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 142      !  (142,16)  Ice change initialisation positive
  IF ( sf(item,sect) ) THEN
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work3d(i,j,k) = MAX(0.0,qcf_dini(i,j,k))
        END DO
      END DO
    END DO
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),work3d,      &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 143      !  (143,16)  Ice change initialisation negative
  IF ( sf(item,sect) ) THEN
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work3d(i,j,k) = MIN(0.0,qcf_dini(i,j,k))
        END DO
      END DO
    END DO
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),work3d,      &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 172      !  (172,16)  Total cloud change from initialisation
  IF ( sf(item,sect) ) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),cf_dini,     &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 173      ! (173,16)  Liquid cloud change from initialisation
  IF ( sf(item,sect) ) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),cfl_dini,    &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 146      !  (146,16)  Liquid cloud change init'n positive
  IF ( sf(item,sect) ) THEN
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work3d(i,j,k) = MAX(0.0,cfl_dini(i,j,k))
        END DO
      END DO
    END DO
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),work3d,      &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 147      !  (147,16)  Liquid cloud change init'n negative
  IF ( sf(item,sect) ) THEN
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work3d(i,j,k) = MIN(0.0,cfl_dini(i,j,k))
        END DO
      END DO
    END DO
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),work3d,      &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 174      !  (174,16)  Ice cloud change from initialisation
  IF ( sf(item,sect) ) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),cff_dini,    &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 148     !  (148,16)  Ice cloud change init'n positive
  IF ( sf(item,sect) ) THEN
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work3d(i,j,k) = MAX(0.0,cff_dini(i,j,k))
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),work3d,      &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

  item = 149     !  (149,16)  Ice cloud change init'n negative
  IF ( sf(item,sect) ) THEN
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work3d(i,j,k) = MIN(0.0,cff_dini(i,j,k))
        END DO
      END DO
    END DO

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3D(stashwork(si(item,sect,im_index)),work3d,      &
          row_length,rows,model_levels,0,0,0,0, at_extremity,       &
          stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
          stash_levels,num_stash_levels+1,                          &
          im_ident,sect,item,                                       &
          icode,cmessage)

  END IF   ! SF(  item,sect)

END IF ! i_cld_vn == i_cld_pc2

item = 201    !  (201,16) is geopotential height on theta levels

IF ( sf(item,sect) ) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3D(stashwork(si(item,sect,im_index)),hts_theta,   &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        im_ident,sect,item,                                       &
        icode,cmessage)

END IF   ! SF(  item,sect)

item = 206    !  (206,16) Cloud water content (qc) on q levels

IF ( sf(item,sect) ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3D(stashwork(si(item,sect,im_index)),qc,          &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        im_ident,sect,item,                                       &
        icode,cmessage)

END IF   ! SF(  item,sect)

item = 207    !  (207,16) Total specific humidity (qT) on q levels

IF ( sf(item,sect) ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3D(stashwork(si(item,sect,im_index)),qT,          &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        im_ident,sect,item,                                       &
        icode,cmessage)

END IF   ! SF(  item,sect)

item = 255    !  (255,16) is geopotential height on rho levels

IF ( sf(item,sect) ) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3D(stashwork(si(item,sect,im_index)),hts_rho,     &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
        stash_levels,num_stash_levels+1,                          &
        im_ident,sect,item,                                       &
        icode,cmessage)

END IF   ! SF(  item,sect)

IF (sf(140,sect) .OR. sf(141,sect) .OR.                           &
    sf(142,sect) .OR. sf(143,sect) .OR.                           &
    sf(146,sect) .OR. sf(147,sect) .OR.                           &
    sf(148,sect) .OR. sf(149,sect) .OR.                           &
    sf(150,sect) .OR. sf(151,sect) .OR.                           &
    sf(156,sect) .OR. sf(157,sect) ) THEN
  DEALLOCATE ( work3d )
END IF


! DEPENDS ON: stash
CALL stash(atmos_sm,atmos_im,16,stashwork,                        &
           icode,cmessage)

9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE st_diag2
