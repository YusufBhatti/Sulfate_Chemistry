! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  End of Timestep climate diagnostics control routine
!
! Subroutine Interface:
SUBROUTINE st_diag3( stashwork,int30,                             &
    energy_corr_now, sin_v_latitude,                              &
    wet_to_dry_n, co2_mmr,                                        &
    icode,cmessage)

USE level_heights_mod
USE trignometric_mod, ONLY: sin_theta_latitude, cos_theta_latitude,&
                            sin_theta_longitude, cos_theta_longitude
USE cderived_mod, ONLY: delta_lambda, delta_phi
USE rad_mask_trop_mod, ONLY: max_trop_level,min_trop_level
USE jules_snow_mod, ONLY: nsmax
USE jules_sea_seaice_mod, ONLY: nice, nice_use

USE atm_fields_bounds_mod, ONLY:                                  &
  pdims, tdims, tdims_s, udims, vdims, wdims

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport
USE UM_ParVars
USE Control_Max_Sizes
USE lbc_mod
USE check_prod_levs_mod, ONLY: check_prod_levs
USE eot_diag_mod, ONLY: eot_diag
USE cppxref_mod, ONLY: ppxref_codelen, ppxref_charlen
USE dump_headers_mod, ONLY: a_realhd, rh_tot_mass_init
USE submodel_mod, ONLY: atmos_sm, atmos_im
USE stash_array_mod, ONLY:                                                     &
    nitems, stindex, stlist, num_stash_levels, stash_levels, si, sf
USE free_tracers_inputs_mod, ONLY: a_tracer_last, a_tracer_first

USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, bl_levels, global_row_length,           &
    land_field, len1_lookup, len_dumphist,                             &
    len_fixhd, len_tot, model_levels, mpp_len1_lookup, n_cca_lev,      &
    n_obj_d1_max, n_rows, ntiles, river_row_length,                    &
    river_rows, row_length, rows, sm_levels, theta_field_size,         &
    theta_off_size, tr_levels, tr_ukca, tr_vars, u_field_size,         &
    v_field_size

USE atm_fields_mod, ONLY: exner_rho_levels, exner_theta_levels, &
    p_theta_levels, p, rho, u, v, w, theta, q, qcl, qcf, pstar, co2, dryrho

USE eot_increments_mod, ONLY:                                                 &
    eot_inc_t, eot_inc_rho, eot_inc_u, eot_inc_v, eot_inc_w,                  &
    eot_inc_q, eot_inc_qcl, eot_inc_qcf, eot_inc_qrain, eot_inc_qgraup,       &
    eot_inc_qcf2, eot_inc_m_v, eot_inc_m_cl, eot_inc_m_cf, eot_inc_m_r,       &
    eot_inc_m_gr, eot_inc_m_cf2, eot_inc_cf, eot_inc_cfl, eot_inc_cff,        &
    eot_inc_theta, eot_inc_thetav, eot_inc_dryrho

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


! Description:
! Calculates end of timestep diagnostics (held in STASH section 30)
! of which there are three types:
! 1) Model level fields and products
! 2) Pressure level fields and products ie U V T UU VV TT UT VT
! 3) Single level fields eg mountain torque and column integrations
!    eq atmospheric mass
!
! Method:
!   This routine sets controls the EoT diagnostics, and sets up
!   array sizes etc for a call to EOT_DIAG which calculates the
!   fields. STASHWORK is filled with the required data, and
!   passed up to atm_step2, ready for a call to stash.
!   The stash items are split into centuries, as follows:
!
!   000s   Standard Model level fields and products
!   100s   Other Model level fields
!   200s   Standard Pressure level fields and products
!   300s   Other Pressure level fields
!   400s   Surface fields and Column integrals
!
!  Standard diagnostics are assigned a number, which are:
!  1 U wind, 2 V wind, 3 W wind, 4 Temperature, 5 Q, 6 RH, 7 Z, 8 Omega
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: stash
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
!


INTEGER, PARAMETER ::                                             &
  npress_diags=8                                                  &
                           ! No. of diags for products on pressure levels
, nmodel_diags=7           ! No. of diags on modellevels

!   Scalar  arguments with intent(in):

!   Array  arguments with intent(in):
REAL :: energy_corr_now                                              &
, wet_to_dry_n(tdims_s%i_start:tdims_s%i_end,                        &
               tdims_s%j_start:tdims_s%j_end,                        &
               tdims_s%k_start:tdims_s%k_end)


REAL :: dry_atmosphere_mass

!   Scalar arguments with intent(InOut):
INTEGER ::                                                        &
        int30,                                                    &
                          ! Dummy for STASH_MAXLEN(30)
        icode             ! Out return code : 0 Normal exit
                          !                 : >0 Error exit

CHARACTER(LEN=errormessagelength) ::                              &
        cmessage          ! Out error message if ICODE > 0

!   Array  arguments with intent(InOut):
REAL ::                                                           &
  stashwork(*)

!   Scalar arguments with intent(out):
INTEGER ::                                                        &
  u_m_levs,                                                       &
  v_m_levs,                                                       &
  w_m_levs,                                                       &
  t_m_levs,                                                       &
  q_m_levs,                                                       &
  z_m_levs,                                                       &
  ke_m_levs,                                                      &
  om_m_levs,                                                      &
  wbig_m_levs,                                                    &
  wbig2_m_levs,                                                   &
  dry_mass_m_levs,                                                &
  Rh_m_levs,                                                      &
  u_mm_levs,                                                      &
  v_mm_levs,                                                      &
  w_mm_levs,                                                      &
  t_mm_levs,                                                      &
  q_mm_levs,                                                      &
  z_mm_levs,                                                      &
  ke_mm_levs,                                                     &
  teot_m_levs,                                                    &
  u_inc_levs,                                                     &
  v_inc_levs,                                                     &
  w_inc_levs,                                                     &
  t_inc_levs,                                                     &
  q_inc_levs,                                                     &
  qcl_inc_levs,                                                   &
  qcf_inc_levs,                                                   &
  qrain_inc_levs,                                                 &
  qgraup_inc_levs,                                                &
  qqcf2_inc_levs,                                                 &
  m_v_inc_levs,                                                   &
  m_cl_inc_levs,                                                  &
  m_cf_inc_levs,                                                  &
  m_r_inc_levs,                                                   &
  m_gr_inc_levs,                                                  &
  m_cf2_inc_levs,                                                 &
  cf_inc_levs,                                                    &
  cfl_inc_levs,                                                   &
  cff_inc_levs,                                                   &
  rho_inc_levs,                                                   &
  theta_inc_levs,                                                 &
  thetav_inc_levs,                                                &
  dryrho_inc_levs,                                                &
  u2_inc_levs,                                                    &
  v2_inc_levs,                                                    &
  w2_inc_levs,                                                    &
  t2_inc_levs,                                                    &
  q2_inc_levs,                                                    &
  qcl2_inc_levs,                                                  &
  qcf2_inc_levs,                                                  &
  rho2_inc_levs,                                                  &
  u_p_levs,                                                       &
  v_p_levs,                                                       &
  w_p_levs,                                                       &
  t_p_levs,                                                       &
  q_p_levs,                                                       &
  rh_p_levs,                                                      &
  z_p_levs,                                                       &
  om_p_levs,                                                      &
  heavy_p_levs,                                                   &
  w_pt_levs,                                                      &
  t_pt_levs,                                                      &
  q_pt_levs,                                                      &
  rh_pt_levs,                                                     &
  z_pt_levs,                                                      &
  om_pt_levs,                                                     &
  heavy_pt_levs,                                                  &
  tv_p_levs,                                                      &
  tvom_p_levs,                                                    &
  tem_p_levs,                                                     &
  vstarbar_p_levs,                                                &
  wstarbar_p_levs,                                                &
  fy_p_levs,                                                      &
  fz_p_levs,                                                      &
  divf_p_levs,                                                    &
  heatflux_p_levs,                                                &
  momflux_p_levs,                                                 &
  field1_p_levs,                                                  &
  field2_p_levs,                                                  &
  tc_track_p_levs,                                                &
  n_levels                                                        &
  ,im_index               !  Internal Model Index for Stash Arrays

INTEGER ::                                                        &
  pt001,pt002,pt003,pt004,pt005,pt006,pt007,pt008,                &
  pt101,pt102,pt103,pt104,pt105,pt106,pt107,                      &
  pt111,pt112,pt113,pt114,pt115,                                  &
  pt171,pt172,pt173,pt174,pt175,pt176,pt177,pt178,                &
  pt181,pt182,pt183,pt184,pt185,pt186,pt187,pt188,pt189,pt190,    &
  pt191,pt192,pt193,pt194,pt195,pt196,pt197,pt198,pt199,pt200,    &
  pt201,pt202,pt203,pt204,pt205,pt206,pt207,pt208,                &
  pt293,pt294,pt295,pt296,pt297,pt298,                            &
  pt301,pt302,pt303,pt304,                                        &
  pt310,pt311,pt312,pt313,pt314,pt315,pt316,                      &
  pt401,pt402,pt403,pt404,pt405,pt406,pt407,pt408,pt409,pt410,    &
  pt411,pt412,pt413,pt414,pt415,pt416,pt417,pt418,pt419,pt420,    &
  pt421,pt422,pt423,pt424,pt425,pt426,pt427,pt428,pt429,pt430,    &
  pt431,pt432,pt433,pt434,pt435,pt436,pt437,pt438,pt439           &
 ,pt440,pt441,pt442,pt451,pt452,pt453,pt454,                      &
  pt455,pt456,pt457,pt458,pt459,pt460,pt461,pt462,pt463,pt464,    &
  pt465,pt466,pt467,pt901,pt902,pt903
!   Array  arguments with intent(out):
REAL ::                                                           &
  u_press(num_stash_levels)                                       &
  ,v_press(num_stash_levels)                                      &
  ,w_press(num_stash_levels)                                      &
  ,t_press(num_stash_levels)                                      &
  ,q_press(num_stash_levels)                                      &
  ,rh_press(num_stash_levels)                                     &
  ,z_press(num_stash_levels)                                      &
                             ! Z pressure levels
  ,om_press(num_stash_levels)                                     &
  ,heavy_press(num_stash_levels)                                  &
  ,w_presst(num_stash_levels)                                     &
  ,t_presst(num_stash_levels)                                     &
  ,q_presst(num_stash_levels)                                     &
  ,rh_presst(num_stash_levels)                                    &
  ,z_presst(num_stash_levels)                                     &
  ,om_presst(num_stash_levels)                                    &
  ,heavy_presst(num_stash_levels)                                 &
  ,press_levs(num_stash_levels)                                   &
  ,tem_press(num_stash_levels)                                    &
  ,vstarbar_press(num_stash_levels)                               &
  ,wstarbar_press(num_stash_levels)                               &
  ,fy_press(num_stash_levels)                                     &
  ,fz_press(num_stash_levels)                                     &
  ,divf_press(num_stash_levels)                                   &
  ,heatflux_press(num_stash_levels)                               &
  ,momflux_press(num_stash_levels)                                &
  ,tc_track_press(num_stash_levels)
LOGICAL :: u_m_list(model_levels)                                    &
  ,v_m_list(model_levels)                                         &
  ,w_m_list(model_levels)                                         &
  ,t_m_list(model_levels)                                         &
  ,q_m_list(model_levels)                                         &
  ,z_m_list(model_levels)                                         &
  ,ke_m_list(model_levels)                                        &
  ,om_m_list(model_levels)                                        &
  ,wbig_m_list(model_levels)                                      &
  ,wbig2_m_list(model_levels)                                     &
  ,RH_m_list(model_levels)                                        &
  ,dry_mass_m_list(model_levels)                                  &
  ,u_mm_list(model_levels)                                        &
  ,v_mm_list(model_levels)                                        &
  ,w_mm_list(model_levels)                                        &
  ,t_mm_list(model_levels)                                        &
  ,q_mm_list(model_levels)                                        &
  ,z_mm_list(model_levels)                                        &
  ,ke_mm_list(model_levels)                                       &
  ,teot_m_list(model_levels)                                      &
  ,u_inc_list(model_levels)                                       &
  ,v_inc_list(model_levels)                                       &
  ,w_inc_list(model_levels)                                       &
  ,t_inc_list(model_levels)                                       &
  ,q_inc_list(model_levels)                                       &
  ,qcl_inc_list(model_levels)                                     &
  ,qcf_inc_list(model_levels)                                     &
  ,qrain_inc_list(model_levels)                                   &
  ,qgraup_inc_list(model_levels)                                  &
  ,qqcf2_inc_list(model_levels)                                   &
  ,m_v_inc_list(model_levels)                                     &
  ,m_cl_inc_list(model_levels)                                    &
  ,m_cf_inc_list(model_levels)                                    &
  ,m_r_inc_list(model_levels)                                     &
  ,m_gr_inc_list(model_levels)                                    &
  ,m_cf2_inc_list(model_levels)                                   &
  ,cf_inc_list(model_levels)                                      &
  ,cfl_inc_list(model_levels)                                     &
  ,cff_inc_list(model_levels)                                     &
  ,rho_inc_list(model_levels)                                     &
  ,theta_inc_list(model_levels)                                   &
  ,thetav_inc_list(model_levels)                                  &
  ,dryrho_inc_list(model_levels)                                  &
  ,u2_inc_list(model_levels)                                      &
  ,v2_inc_list(model_levels)                                      &
  ,w2_inc_list(model_levels)                                      &
  ,t2_inc_list(model_levels)                                      &
  ,q2_inc_list(model_levels)                                      &
  ,qcl2_inc_list(model_levels)                                    &
  ,qcf2_inc_list(model_levels)                                    &
  ,rho2_inc_list(model_levels)

REAL ::    co2_mmr

LOGICAL :: prod_m_list(model_levels,nmodel_diags,nmodel_diags)

INTEGER ::                                                        &
  prod_IND(num_stash_levels*2,npress_diags,npress_diags),         &
  prod_p_levs(npress_diags,npress_diags)

! Local parameters:

! Local scalars:
INTEGER ::                                                        &
  i,                                                              &
  ni,                                                             &
  k,                                                              &
  isl,                                                            &
  bl,                                                             &
  tl,                                                             &
  level,                                                          &
  code30,nii,nij,j

CHARACTER(LEN=2) :: cfield1,cfield2

!     Check TEM diagnostic levels match requested levels (True = OK)
      LOGICAL :: l_match

!     ErrorStatus - Error flag (0 = OK)
      INTEGER::    ErrorStatus

!     Error message
      CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ST_DIAG3'

! Local dynamic arrays:
REAL ::                                                           &
  field1_press(num_stash_levels)                                  &
  ,field2_press(num_stash_levels)

REAL:: sin_v_latitude (row_length, n_rows)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of header

!     Set to atmosphere internal model
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
im_index = 1

!--------------------Extract Reqd Model Levels for U-------------
isl=stindex(1,001,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    u_m_list,stash_levels,num_stash_levels+1,icode,cmessage)
  u_m_levs=0
  DO i=1,model_levels
    IF (u_m_list(i)) u_m_levs=u_m_levs+1
  END DO
ELSE
  DO i=1,model_levels
    u_m_list(i)=.FALSE.
  END DO
  u_m_levs=1
END IF
!--------------------Extract Reqd Model Levels for V-------------
isl=stindex(1,002,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    v_m_list,stash_levels,num_stash_levels+1,icode,cmessage)
  v_m_levs=0
  DO i=1,model_levels
    IF (v_m_list(i)) v_m_levs=v_m_levs+1
  END DO
ELSE
  DO i=1,model_levels
    v_m_list(i)=.FALSE.
  END DO
  v_m_levs=1
END IF

!--------------------Extract Reqd Model Levels for W-------------
isl=stindex(1,003,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    w_m_list,stash_levels,num_stash_levels+1,icode,cmessage)
  w_m_levs=0
  DO i=1,model_levels
    IF (w_m_list(i)) w_m_levs=w_m_levs+1
  END DO
ELSE
  DO i=1,model_levels
    w_m_list(i)=.FALSE.
  END DO
  w_m_levs=1
END IF

!--------------------Extract Reqd Model Levels for T-------------
isl=stindex(1,004,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    t_m_list,stash_levels,num_stash_levels+1,icode,cmessage)
  t_m_levs=0
  DO i=1,model_levels
    IF (t_m_list(i)) t_m_levs=t_m_levs+1
  END DO
ELSE
  DO i=1,model_levels
    t_m_list(i)=.FALSE.
  END DO
  t_m_levs=1
END IF

!--------------------Extract Reqd Model Levels for Q-------------
isl=stindex(1,005,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    q_m_list,stash_levels,num_stash_levels+1,icode,cmessage)
  q_m_levs=0
  DO i=1,model_levels
    IF (q_m_list(i)) q_m_levs=q_m_levs+1
  END DO
ELSE
  DO i=1,model_levels
    q_m_list(i)=.FALSE.
  END DO
  q_m_levs=1
END IF

!--------------------Extract Reqd Model Levels for Z-------------
isl=stindex(1,006,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    z_m_list,stash_levels,num_stash_levels+1,icode,cmessage)
  z_m_levs=0
  DO i=1,model_levels
    IF (z_m_list(i)) z_m_levs=z_m_levs+1
  END DO
ELSE
  DO i=1,model_levels
    z_m_list(i)=.FALSE.
  END DO
  z_m_levs=1
END IF
!--------------------Extract Reqd Model Levels for KE-------------
isl=stindex(1,007,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    ke_m_list,stash_levels,num_stash_levels+1,icode,cmessage)
  ke_m_levs=0
  DO i=1,model_levels
    IF (ke_m_list(i)) ke_m_levs=ke_m_levs+1
  END DO
ELSE
  DO i=1,model_levels
    ke_m_list(i)=.FALSE.
  END DO
  ke_m_levs=1
END IF
!--------------------Extract Reqd Model Levels for Omega-------------
isl=stindex(1,008,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    om_m_list,stash_levels,num_stash_levels+1,icode,cmessage)
  om_m_levs=0
  DO i=1,model_levels
    IF (om_m_list(i)) om_m_levs=om_m_levs+1
  END DO
ELSE
  DO i=1,model_levels
    om_m_list(i)=.FALSE.
  END DO
  om_m_levs=1
END IF
!--------------Extract Reqd Model Levels for Dry mass weighting----
isl=stindex(1,115,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
  dry_mass_M_LIST,stash_levels,num_stash_levels+1,icode,cmessage)
  dry_mass_m_levs=0
  DO i=1,model_levels
    IF (dry_mass_m_list(i)) dry_mass_m_levs=dry_mass_m_levs+1
  END DO
ELSE
  DO i=1,model_levels
    dry_mass_m_list(i)=.FALSE.
  END DO
  dry_mass_m_levs=1
END IF

!-Extract Reqd model levels for products and store in prod_m_list

DO i=1,nmodel_diags
  DO j=1,nmodel_diags
    code30=i*10+j
    IF (sf(code30,30)) THEN
      isl=stindex(1,code30,30,im_index)
      ! DEPENDS ON: set_levels_list
      CALL set_levels_list(model_levels,nitems,stlist(1,isl),     &
        prod_m_list(1,i,j),                                       &
        stash_levels,num_stash_levels+1,icode,cmessage)
    ELSE
      DO k=1,model_levels
        prod_M_LIST(k,i,j)=.FALSE.
      END DO
    END IF
  END DO
END DO


!--------------------Extract Reqd Model Levels for wbig----------
isl=stindex(1,112,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    Wbig_M_LIST,stash_levels,num_stash_levels+1,icode,cmessage)
  Wbig_M_LEVS=0
  DO i=1,model_levels
    IF (Wbig_M_LIST(i)) Wbig_M_LEVS=Wbig_M_LEVS+1
  END DO
ELSE
  Wbig_M_LEVS=1
END IF

!--------------------Extract Reqd Model Levels for wbig2----------
isl=stindex(1,114,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    Wbig2_M_LIST,stash_levels,num_stash_levels+1,icode,cmessage)
  Wbig2_M_LEVS=0
  DO i=1,model_levels
    IF (Wbig2_M_LIST(i)) Wbig2_M_LEVS=Wbig2_M_LEVS+1
  END DO
ELSE
  DO i=1,model_levels
    Wbig2_M_LIST(i)=.FALSE.
  END DO
  Wbig2_M_LEVS=1
END IF
!--------------------Extract Reqd Model Levels for RH ----------
isl=stindex(1,113,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    rh_m_list,stash_levels,num_stash_levels+1,icode,cmessage)
  rh_m_levs=0
  DO i=1,model_levels
    IF (rh_m_list(i)) rh_m_levs=rh_m_levs+1
  END DO
ELSE
  rh_m_levs=1
END IF

!--------------------Extract Reqd Model Levels for U mass weighted
isl=stindex(1,101,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    u_mm_list,stash_levels,num_stash_levels+1,icode,cmessage)
  u_mm_levs=0
  DO i=1,model_levels
    IF (u_mm_list(i)) u_mm_levs=u_mm_levs+1
  END DO
ELSE
  DO i=1,model_levels
    u_mm_list(i)=.FALSE.
  END DO
  u_mm_levs=1
END IF
!--------------------Extract Reqd Model Levels for V mass weighted
isl=stindex(1,102,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    v_mm_list,stash_levels,num_stash_levels+1,icode,cmessage)
  v_mm_levs=0
  DO i=1,model_levels
    IF (v_mm_list(i)) v_mm_levs=v_mm_levs+1
  END DO
ELSE
  DO i=1,model_levels
    v_mm_list(i)=.FALSE.
  END DO
  v_mm_levs=1
END IF

!--------------------Extract Reqd Model Levels for W mass weighted
isl=stindex(1,103,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    w_mm_list,stash_levels,num_stash_levels+1,icode,cmessage)
  w_mm_levs=0
  DO i=1,model_levels
    IF (w_mm_list(i)) w_mm_levs=w_mm_levs+1
  END DO
ELSE
  DO i=1,model_levels
    w_mm_list(i)=.FALSE.
  END DO
  w_mm_levs=1
END IF

!--------------------Extract Reqd Model Levels for T mass weighted
isl=stindex(1,104,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    t_mm_list,stash_levels,num_stash_levels+1,icode,cmessage)
  t_mm_levs=0
  DO i=1,model_levels
    IF (t_mm_list(i)) t_mm_levs=t_mm_levs+1
  END DO
ELSE
  DO i=1,model_levels
    t_mm_list(i)=.FALSE.
  END DO
  t_mm_levs=1
END IF

!--------------------Extract Reqd Model Levels for Q mass weighted
isl=stindex(1,105,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    q_mm_list,stash_levels,num_stash_levels+1,icode,cmessage)
  q_mm_levs=0
  DO i=1,model_levels
    IF (q_mm_list(i)) q_mm_levs=q_mm_levs+1
  END DO
ELSE
  DO i=1,model_levels
    q_mm_list(i)=.FALSE.
  END DO
  q_mm_levs=1
END IF

!--------------------Extract Reqd Model Levels for Z mass weighted
isl=stindex(1,106,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    z_mm_list,stash_levels,num_stash_levels+1,icode,cmessage)
  z_mm_levs=0
  DO i=1,model_levels
    IF (z_mm_list(i)) z_mm_levs=z_mm_levs+1
  END DO
ELSE
  DO i=1,model_levels
    z_mm_list(i)=.FALSE.
  END DO
  z_mm_levs=1
END IF

!--------------------Extract Reqd Model Levels for KE mass weighted
isl=stindex(1,107,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    ke_mm_list,stash_levels,num_stash_levels+1,icode,cmessage)
  ke_mm_levs=0
  DO i=1,model_levels
    IF (ke_mm_list(i)) ke_mm_levs=ke_mm_levs+1
  END DO
ELSE
  DO i=1,model_levels
    ke_mm_list(i)=.FALSE.
  END DO
  ke_mm_levs=1
END IF


!--------------------Extract Reqd Model Levels for T at EOT
isl=stindex(1,111,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    teot_m_list,stash_levels,num_stash_levels+1,icode,cmessage)
  teot_m_levs=0
  DO i=1,model_levels
    IF (teot_m_list(i)) teot_m_levs=teot_m_levs+1
  END DO
ELSE
  DO i=1,model_levels
    teot_m_list(i)=.FALSE.
  END DO
  teot_m_levs=1
END IF

!--------------------Extract Reqd Model Levels for T increment
isl=stindex(1,181,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    t_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  t_inc_levs=0
  DO i=1,model_levels
    IF (t_inc_list(i)) t_inc_levs=t_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    t_inc_list(i)=.FALSE.
  END DO
  t_inc_levs=1
END IF

!--------------------Extract Reqd Model Levels for Q increment
isl=stindex(1,182,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    q_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  q_inc_levs=0
  DO i=1,model_levels
    IF (q_inc_list(i)) q_inc_levs=q_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    q_inc_list(i)=.FALSE.
  END DO
  q_inc_levs=1
END IF

!--------------------Extract Reqd Model Levels for QCL increment
isl=stindex(1,183,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    qcl_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  qcl_inc_levs=0
  DO i=1,model_levels
    IF (qcl_inc_list(i)) qcl_inc_levs=qcl_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    qcl_inc_list(i)=.FALSE.
  END DO
  qcl_inc_levs=1
END IF

!--------------------Extract Reqd Model Levels for QCF increment
isl=stindex(1,184,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    qcf_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  qcf_inc_levs=0
  DO i=1,model_levels
    IF (qcf_inc_list(i)) qcf_inc_levs=qcf_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    qcf_inc_list(i)=.FALSE.
  END DO
  qcf_inc_levs=1
END IF

!--------------------Extract Reqd Model Levels for U increment
isl=stindex(1,185,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    u_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  u_inc_levs=0
  DO i=1,model_levels
    IF (u_inc_list(i)) u_inc_levs=u_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    u_inc_list(i)=.FALSE.
  END DO
  u_inc_levs=1
END IF

!--------------------Extract Reqd Model Levels for V increment
isl=stindex(1,186,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    v_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  v_inc_levs=0
  DO i=1,model_levels
    IF (v_inc_list(i)) v_inc_levs=v_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    v_inc_list(i)=.FALSE.
  END DO
  v_inc_levs=1
END IF

!--------------------Extract Reqd Model Levels for W increment
isl=stindex(1,187,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    w_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  w_inc_levs=0
  DO i=1,model_levels
    IF (w_inc_list(i)) w_inc_levs=w_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    w_inc_list(i)=.FALSE.
  END DO
  w_inc_levs=1
END IF

!--------------------Extract Reqd Model Levels for RHO increment
isl=stindex(1,188,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    rho_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  rho_inc_levs=0
  DO i=1,model_levels
    IF (rho_inc_list(i)) rho_inc_levs=rho_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    rho_inc_list(i)=.FALSE.
  END DO
  rho_inc_levs=1
END IF

!--------------------Extract Reqd Model Levels for THETA increment
isl=stindex(1,901,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    theta_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  theta_inc_levs=0
  DO i=1,model_levels
    IF (theta_inc_list(i)) theta_inc_levs=theta_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    theta_inc_list(i)=.FALSE.
  END DO
  theta_inc_levs=1
END IF

!--------------------Extract Reqd Model Levels for THETAV increment
isl=stindex(1,902,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    thetav_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  thetav_inc_levs=0
  DO i=1,model_levels
    IF (thetav_inc_list(i)) thetav_inc_levs=thetav_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    thetav_inc_list(i)=.FALSE.
  END DO
  thetav_inc_levs=1
END IF

!--------------------Extract Reqd Model Levels for DRYRHO increment
isl=stindex(1,903,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    dryrho_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  dryrho_inc_levs=0
  DO i=1,model_levels
    IF (dryrho_inc_list(i)) dryrho_inc_levs=dryrho_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    dryrho_inc_list(i)=.FALSE.
  END DO
  dryrho_inc_levs=1
END IF

!--------------------Extract Reqd Model Levels for qrain increment
isl=stindex(1,189,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
         qrain_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  qrain_inc_levs = 0
  DO i=1,model_levels
    IF (qrain_inc_list(i)) qrain_inc_levs=qrain_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    qrain_inc_list(i)=.FALSE.
  END DO
  qrain_inc_levs = 1
END IF

!--------------------Extract Reqd Model Levels for qgraup increment
isl=stindex(1,190,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
         qgraup_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  qgraup_inc_levs = 0
  DO i=1,model_levels
    IF (qgraup_inc_list(i)) qgraup_inc_levs=qgraup_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    qgraup_inc_list(i)=.FALSE.
  END DO
  qgraup_inc_levs = 1
END IF

!--------------------Extract Reqd Model Levels for qcf2 increment
isl=stindex(1,191,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
         qqcf2_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  qqcf2_inc_levs = 0
  DO i=1,model_levels
    IF (qqcf2_inc_list(i)) qqcf2_inc_levs=qqcf2_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    qqcf2_inc_list(i)=.FALSE.
  END DO
  qqcf2_inc_levs = 1
END IF

!--------------------Extract Reqd Model Levels for cf increment
isl=stindex(1,192,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
         cf_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  cf_inc_levs = 0
  DO i=1,model_levels
    IF (cf_inc_list(i)) cf_inc_levs=cf_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    cf_inc_list(i)=.FALSE.
  END DO
  cf_inc_levs = 1
END IF

!--------------------Extract Reqd Model Levels for cfl increment
isl=stindex(1,193,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
         cfl_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  cfl_inc_levs = 0
  DO i=1,model_levels
    IF (cfl_inc_list(i)) cfl_inc_levs=cfl_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    cfl_inc_list(i)=.FALSE.
  END DO
  cfl_inc_levs = 1
END IF

!--------------------Extract Reqd Model Levels for cff increment
isl=stindex(1,194,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
         cff_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  cff_inc_levs = 0
  DO i=1,model_levels
    IF (cff_inc_list(i)) cff_inc_levs=cff_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    cff_inc_list(i)=.FALSE.
  END DO
  cff_inc_levs = 1
END IF

!--------------------Extract Reqd Model Levels for m_v increment
isl=stindex(1,195,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
         m_v_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  m_v_inc_levs = 0
  DO i=1,model_levels
    IF (m_v_inc_list(i)) m_v_inc_levs=m_v_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    m_v_inc_list(i)=.FALSE.
  END DO
  m_v_inc_levs = 1
END IF

!--------------------Extract Reqd Model Levels for m_cl increment
isl=stindex(1,196,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
         m_cl_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  m_cl_inc_levs = 0
  DO i=1,model_levels
    IF (m_cl_inc_list(i)) m_cl_inc_levs=m_cl_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    m_cl_inc_list(i)=.FALSE.
  END DO
  m_cl_inc_levs = 1
END IF

!--------------------Extract Reqd Model Levels for m_cf increment
isl=stindex(1,197,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
         m_cf_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  m_cf_inc_levs = 0
  DO i=1,model_levels
    IF (m_cf_inc_list(i)) m_cf_inc_levs=m_cf_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    m_cf_inc_list(i)=.FALSE.
  END DO
  m_cf_inc_levs = 1
END IF

!--------------------Extract Reqd Model Levels for m_r increment
isl=stindex(1,198,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
         m_r_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  m_r_inc_levs = 0
  DO i=1,model_levels
    IF (m_r_inc_list(i)) m_r_inc_levs=m_r_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    m_r_inc_list(i)=.FALSE.
  END DO
  m_r_inc_levs = 1
END IF

!--------------------Extract Reqd Model Levels for m_gr increment
isl=stindex(1,199,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
         m_gr_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  m_gr_inc_levs = 0
  DO i=1,model_levels
    IF (m_gr_inc_list(i)) m_gr_inc_levs=m_gr_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    m_gr_inc_list(i)=.FALSE.
  END DO
  m_gr_inc_levs = 1
END IF

!--------------------Extract Reqd Model Levels for m_cf2 increment
isl=stindex(1,200,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
         m_cf2_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  m_cf2_inc_levs = 0
  DO i=1,model_levels
    IF (m_cf2_inc_list(i)) m_cf2_inc_levs=m_cf2_inc_levs+1
  END DO
ELSE
  DO i=1,model_levels
    m_cf2_inc_list(i)=.FALSE.
  END DO
  m_cf2_inc_levs = 1
END IF

!--------------------Extract Reqd Model Levels for T increment **2
isl=stindex(1,171,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    t2_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  t2_inc_levs=0
  DO i=1,model_levels
    IF (t2_inc_list(i)) t2_inc_levs=t2_inc_levs+1
  END DO
ELSE
  t2_inc_levs=1
END IF

!--------------------Extract Reqd Model Levels for Q increment **2
isl=stindex(1,172,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    q2_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  q2_inc_levs=0
  DO i=1,model_levels
    IF (q2_inc_list(i)) q2_inc_levs=q2_inc_levs+1
  END DO
ELSE
  q2_inc_levs=1
END IF

!--------------------Extract Reqd Model Levels for QCL increment **2
isl=stindex(1,173,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    qcl2_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  qcl2_inc_levs=0
  DO i=1,model_levels
    IF (qcl2_inc_list(i)) qcl2_inc_levs=qcl2_inc_levs+1
  END DO
ELSE
  qcl2_inc_levs=1
END IF

!--------------------Extract Reqd Model Levels for QCF increment **2
isl=stindex(1,174,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    qcf2_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  qcf2_inc_levs=0
  DO i=1,model_levels
    IF (qcf2_inc_list(i)) qcf2_inc_levs=qcf2_inc_levs+1
  END DO
ELSE
  qcf2_inc_levs=1
END IF

!--------------------Extract Reqd Model Levels for U increment **2
isl=stindex(1,175,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    u2_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  u2_inc_levs=0
  DO i=1,model_levels
    IF (u2_inc_list(i)) u2_inc_levs=u2_inc_levs+1
  END DO
ELSE
  u2_inc_levs=1
END IF

!--------------------Extract Reqd Model Levels for V increment **2
isl=stindex(1,176,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    v2_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  v2_inc_levs=0
  DO i=1,model_levels
    IF (v2_inc_list(i)) v2_inc_levs=v2_inc_levs+1
  END DO
ELSE
  v2_inc_levs=1
END IF
!--------------------Extract Reqd Model Levels for W increment **2
isl=stindex(1,177,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    w2_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  w2_inc_levs=0
  DO i=1,model_levels
    IF (w2_inc_list(i)) w2_inc_levs=w2_inc_levs+1
  END DO
ELSE
  w2_inc_levs=1
END IF

!--------------------Extract Reqd Model Levels for RHO increment **2
isl=stindex(1,178,30,im_index)
IF (isl >  0) THEN
  ! DEPENDS ON: set_levels_list
  CALL set_levels_list(model_levels,nitems,stlist(1,isl),         &
    rho2_inc_list,stash_levels,num_stash_levels+1,icode,cmessage)
  rho2_inc_levs=0
  DO i=1,model_levels
    IF (rho2_inc_list(i)) rho2_inc_levs=rho2_inc_levs+1
  END DO
ELSE
  rho2_inc_levs=1
END IF

!--------------------Extract Reqd Pressures for U_COMP_P-------------

isl=stindex(1,201,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  u_p_levs=stash_levels(1,ni)
  DO k =1,u_p_levs
    u_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  u_p_levs=1
END IF

!--------------------Extract Reqd Pressures for V_COMP_P-------------

isl=stindex(1,202,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  v_p_levs=stash_levels(1,ni)
  DO k =1,v_p_levs
    v_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  v_p_levs=1
END IF

!--------------------Extract Reqd Pressures for W_COMP_P-------------

isl=stindex(1,203,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  w_p_levs=stash_levels(1,ni)
  DO k =1,w_p_levs
    w_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  w_p_levs=1
END IF

!--------------------Extract Reqd Pressures for T on wind grid-------

isl=stindex(1,204,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  t_p_levs=stash_levels(1,ni)
  DO k =1,t_p_levs
    t_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  t_p_levs=1
END IF


!--------------------Extract Reqd Pressures for q on wind grid-------

isl=stindex(1,205,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  q_p_levs=stash_levels(1,ni)
  DO k =1,q_p_levs
    q_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  q_p_levs=1
END IF

!--------------------Extract Reqd Pressures for RH-----

isl=stindex(1,206,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  rh_p_levs=stash_levels(1,ni)
  DO k =1,rh_p_levs
    rh_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  rh_p_levs=1
END IF

!--------------------Extract Reqd Pressures for geopotential height-----

isl=stindex(1,207,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  z_p_levs=stash_levels(1,ni)
  DO k =1,z_p_levs
    z_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  z_p_levs=1
END IF

!--------------------Extract Reqd Pressures for Omega---

isl=stindex(1,208,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  om_p_levs=stash_levels(1,ni)
  DO k =1,om_p_levs
    om_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  om_p_levs=1
END IF

!--------------------Extract Reqd Pressures for W_COMP_P on T points--

isl=stindex(1,293,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  w_pt_levs=stash_levels(1,ni)
  DO k =1,w_pt_levs
    w_presst(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  w_pt_levs=1
END IF

!--------------------Extract Reqd Pressures for T on T points-------

isl=stindex(1,294,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  t_pt_levs=stash_levels(1,ni)
  DO k =1,t_pt_levs
    t_presst(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  t_pt_levs=1
END IF


!--------------------Extract Reqd Pressures for q on T points-------

isl=stindex(1,295,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  q_pt_levs=stash_levels(1,ni)
  DO k =1,q_pt_levs
    q_presst(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  q_pt_levs=1
END IF

!--------------------Extract Reqd Pressures for RH on T points-----

isl=stindex(1,296,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  rh_pt_levs=stash_levels(1,ni)
  DO k =1,rh_pt_levs
    rh_presst(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  rh_pt_levs=1
END IF

!-----Extract Reqd Pressures for geopotential height on T points-----

isl=stindex(1,297,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  z_pt_levs=stash_levels(1,ni)
  DO k =1,z_pt_levs
    z_presst(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  z_pt_levs=1
END IF

!--------------------Extract Reqd Pressures for Omega on T points---

isl=stindex(1,298,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  om_pt_levs=stash_levels(1,ni)
  DO k =1,om_pt_levs
    om_presst(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  om_pt_levs=1
END IF

!--------------------Extract Reqd Pressures for Vstarbar---

isl=stindex(1,310,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  vstarbar_p_levs=stash_levels(1,ni)
  DO k =1,vstarbar_p_levs
    vstarbar_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  vstarbar_p_levs=1
END IF

!-----Temporary solution define TEM Pressures = Vstarbar---
      TEM_P_LEVS = VSTARBAR_P_LEVS
      IF(TEM_P_LEVS > 1) THEN
        DO K =1,TEM_P_LEVS
          TEM_PRESS(K) = VSTARBAR_PRESS(K)
        END DO
      END IF

!-----Can do something cleverer with CHECK_PROD_LEVS style-
!     code to mask output diagnostics onto TEM_PRESS levels
!     Meanwhile simply check that levels are exact match---

!--------------------Extract Reqd Pressures for Wstarbar---

isl=stindex(1,311,30,im_index)
IF (isl >  0 .AND. TEM_P_LEVS > 1) THEN
  ni=-stlist(10,isl)
  wstarbar_p_levs=stash_levels(1,ni)
  l_match = .TRUE.
  DO k =1,wstarbar_p_levs
    wstarbar_press(k)=stash_levels(k+1,ni)/1000.0
    IF (wstarbar_press(k) /= tem_press(k)) THEN
      l_match = .FALSE.

      ErrorStatus =  1        ! Error – Fail
      Cmessage='Wstarbar not on TEM pressure levels - '           &
      //'requested pressure levels do not match, diagnostic aborted'

      CALL Ereport(Routinename,ErrorStatus,Cmessage)

      EXIT
    END IF
  END DO
ELSE
  wstarbar_p_levs=1
END IF

!----Extract Reqd Pressures for EP Flux (Phi component)---

isl=stindex(1,312,30,im_index)
IF (isl >  0 .AND. TEM_P_LEVS > 1) THEN
  ni=-stlist(10,isl)
  fy_p_levs=stash_levels(1,ni)
  l_match = .TRUE.
  DO K =1,fy_p_levs
    fy_press(k)=stash_levels(k+1,ni)/1000.0
    IF (fy_press(k) /= tem_press(k)) THEN
      l_match = .FALSE.

      ErrorStatus =  1        ! Error – Fail
      Cmessage='Fy not on TEM pressure levels - '                 &
      //'requested pressure levels do not match, diagnostic aborted'

      CALL Ereport(Routinename,ErrorStatus,Cmessage)

      EXIT
    END IF
  END DO
ELSE
  fy_p_levs=1
END IF

!----Extract Reqd Pressures for EP Flux (z component)---

isl=stindex(1,313,30,im_index)
IF (isl >  0 .AND. TEM_P_LEVS > 1) THEN
  ni=-stlist(10,isl)
  fz_p_levs=stash_levels(1,ni)
  l_match = .TRUE.
  DO K =1,fz_p_levs
    fz_press(k)=stash_levels(k+1,ni)/1000.0
    IF (fz_press(k) /= tem_press(k)) THEN
      l_match = .FALSE.

      ErrorStatus =  1        ! Error – Fail
      Cmessage='Fz not on TEM pressure levels - '                 &
      //'requested pressure levels do not match, diagnostic aborted'

      CALL Ereport(Routinename,ErrorStatus,Cmessage)

      EXIT
    END IF
  END DO
ELSE
  fz_p_levs=1
END IF

!-----Extract Reqd Pressures for Divergence of EP Flux---

isl=stindex(1,314,30,im_index)
IF (isl >  0 .AND. TEM_P_LEVS > 1) THEN
  ni=-stlist(10,isl)
  divf_p_levs=stash_levels(1,ni)
  l_match = .TRUE.
  DO K =1,divf_p_levs
    divf_press(k)=stash_levels(k+1,ni)/1000.0
    IF (divf_press(k) /= tem_press(k)) THEN
      l_match = .FALSE.

      ErrorStatus =  1        ! Error – Fail
      Cmessage='Div.F not on TEM pressure levels - '              &
      //'requested pressure levels do not match, diagnostic aborted'

      CALL Ereport(Routinename,ErrorStatus,Cmessage)

      EXIT
    END IF
  END DO
ELSE
  divf_p_levs=1
END IF

!-----Extract Reqd Pressures for Meridional heat flux---

isl=stindex(1,315,30,im_index)
IF (isl >  0 .AND. TEM_P_LEVS > 1) THEN
  ni=-stlist(10,isl)
  heatflux_p_levs=stash_levels(1,ni)
  l_match = .TRUE.

  DO k =1,heatflux_p_levs
    heatflux_press(k)=stash_levels(k+1,ni)/1000.0
    IF (heatflux_press(k) /= tem_press(k)) THEN
      l_match = .FALSE.

      ErrorStatus =  1        ! Error – Fail
      Cmessage='Heat Flux not on TEM pressure levels - '          &
      //'requested pressure levels do not match, diagnostic aborted'

      CALL Ereport(Routinename,ErrorStatus,Cmessage)

      EXIT
    END IF
  END DO
ELSE
  heatflux_p_levs=1
END IF


!-----Extract Reqd Pressures for Meridional momentum flux---

isl=stindex(1,316,30,im_index)
IF (isl >  0 .AND. TEM_P_LEVS > 1) THEN
  ni=-stlist(10,isl)
  momflux_p_levs=stash_levels(1,ni)
  l_match = .TRUE.
  DO k =1,momflux_p_levs
    momflux_press(k)=stash_levels(k+1,ni)/1000.0
    IF (momflux_press(k) /= tem_press(k)) THEN
      l_match = .FALSE.

      ErrorStatus =  1        ! Error – Fail
      Cmessage='Momentum Flux not on TEM pressure levels - '      &
      //'requested pressure levels do not match, diagnostic aborted'

      CALL Ereport(Routinename,ErrorStatus,Cmessage)

      EXIT
    END IF
  END DO
ELSE
  momflux_p_levs=1
END IF

!--------------------Extract Reqd Pressures for TRACK-------------

isl=stindex(1,457,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  tc_track_p_levs=stash_levels(1,ni)
  DO k =1,tc_track_p_levs
    tc_track_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  tc_track_p_levs=1
END IF

!--------------------Extract Reqd Pressures for products---


DO i=1,npress_diags
  DO j=1,npress_diags
    IF (i  ==  1) THEN
      cfield1='u'
      field1_press=u_press
      field1_p_levs=u_p_levs
    ELSE IF (i  ==  2) THEN
      cfield1='v'
      field1_press=v_press
      field1_p_levs=v_p_levs
    ELSE IF (i  ==  3) THEN
      cfield1='w'
      field1_press=w_press
      field1_p_levs=w_p_levs
    ELSE IF (i  ==  4) THEN
      cfield1='t'
      field1_press=t_press
      field1_p_levs=t_p_levs
    ELSE IF (i  ==  5) THEN
      cfield1='q'
      field1_press=q_press
      field1_p_levs=q_p_levs
    ELSE IF (i  ==  6) THEN
      cfield1='rh'
      field1_press=rh_press
      field1_p_levs=rh_p_levs
    ELSE IF (i  ==  7) THEN
      cfield1='z'
      field1_press=z_press
      field1_p_levs=z_p_levs
    ELSE IF (i  ==  8) THEN
      cfield1='om'
      field1_press=om_press
      field1_p_levs=om_p_levs
    END IF

    IF (j  ==  1) THEN
      cfield2='u'
      field2_press=u_press
      field2_p_levs=u_p_levs
    ELSE IF (j  ==  2) THEN
      cfield2='v'
      field2_press=v_press
      field2_p_levs=v_p_levs
    ELSE IF (j  ==  3) THEN
      cfield2='w'
      field2_press=w_press
      field2_p_levs=w_p_levs
    ELSE IF (j  ==  4) THEN
      cfield2='t'
      field2_press=t_press
      field2_p_levs=t_p_levs
    ELSE IF (j  ==  5) THEN
      cfield2='q'
      field2_press=q_press
      field2_p_levs=q_p_levs
    ELSE IF (j  ==  6) THEN
      cfield2='rh'
      field2_press=rh_press
      field2_p_levs=rh_p_levs
    ELSE IF (j  ==  7) THEN
      cfield2='z'
      field2_press=z_press
      field2_p_levs=z_p_levs
    ELSE IF (j  ==  8) THEN
      cfield2='om'
      field2_press=om_press
      field2_p_levs=om_p_levs
    END IF
    code30=200+i*10+j
    IF (sf(code30,30)) THEN
      IF ((.NOT. sf(200+i,30)) .OR. (.NOT. sf(200+j,30))) THEN
        cmessage=' '
        cmessage='ST_DIAG3 : '//cfield1//'*'//cfield2//           &
          ' error '//cfield1//' and '//cfield2//                  &
          ' must be requested'
        icode=1
        GO TO 9999
      ELSE
        ni=-stlist(10,stindex(1,code30,30,im_index))
        nii=-stlist(10,stindex(1,200+i,30,im_index))
        nij=-stlist(10,stindex(1,200+j,30,im_index))

        CALL check_prod_levs(                                     &
          num_stash_levels,stash_levels,ni,                       &
          field1_p_levs,field1_press,                             &
          field2_p_levs,field2_press,                             &
          prod_P_LEVS(i,j),prod_IND(1,i,j))
      END IF
    ELSE
      prod_P_LEVS(i,j)=1
    END IF
  END DO
END DO

!     Do checks for other pressure fields

IF (sf(302,30) .AND. .NOT.                                          &
  (sf(245,30) .AND. sf(204,30) .AND. sf(205,30))) THEN
  cmessage='ST_DIAG3 : Virtual temperature'//                     &
          ' error T*Q must also be requested'
  icode=1
  GO TO 9999
END IF

IF (sf(303,30) .AND. .NOT.                                          &
  (sf(245,30) .AND. sf(204,30) .AND.                                 &
  sf(205,30) .AND. sf(208,30))) THEN
  cmessage='ST_DIAG3 : Virtual temperature*omega'//               &
          ' error T*Q and omega must also be requested'
  icode=1
  GO TO 9999
END IF



!--------------------Extract Reqd Pressures for Heavyside function---

isl=stindex(1,301,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  heavy_p_levs=stash_levels(1,ni)
  DO k =1,heavy_p_levs
    heavy_press(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  heavy_p_levs=1
END IF

isl=stindex(1,304,30,im_index)
IF (isl >  0) THEN
  ni=-stlist(10,isl)
  heavy_pt_levs=stash_levels(1,ni)
  DO k =1,heavy_pt_levs
    heavy_presst(k)=stash_levels(k+1,ni)/1000.0
  END DO
ELSE
  heavy_pt_levs=1
END IF

!Lset up level information for tv
IF (sf(302,30)) THEN
  IF (sf(204,30) .AND. sf(205,30)) THEN
    tv_p_levs=prod_p_levs(4,5)
  ELSE
    tv_p_levs=1
    cmessage='ST_DIAG3 : T and Q must be selected to calculate '//&
      'virtual temperature'
    icode=1
    GO TO 9999
  END IF
ELSE
  tv_p_levs=1
END IF

!Lset up level information for tv*om
IF (sf(303,30)) THEN
  IF (sf(204,30) .AND. sf(205,30) .AND. sf(208,30)) THEN
    tvom_p_levs=prod_p_levs(4,5)
  ELSE
    tvom_p_levs=1
    cmessage='ST_DIAG3 : T,Q and Omega'//                         &
      'must be selected to calculate '//                          &
      'virtual temperature'
    icode=1
    GO TO 9999
  END IF
ELSE
  tvom_p_levs=1
END IF



!-------------------Set up Pointers for STASHWORK -------------------
pt001=si(001,30,im_index)
pt002=si(002,30,im_index)
pt003=si(003,30,im_index)
pt004=si(004,30,im_index)
pt005=si(005,30,im_index)
pt006=si(006,30,im_index)
pt007=si(007,30,im_index)
pt008=si(008,30,im_index)
pt101=si(101,30,im_index)
pt102=si(102,30,im_index)
pt103=si(103,30,im_index)
pt104=si(104,30,im_index)
pt105=si(105,30,im_index)
pt106=si(106,30,im_index)
pt107=si(107,30,im_index)
pt111=si(111,30,im_index)
pt112=si(112,30,im_index)
pt113=si(113,30,im_index)
pt114=si(114,30,im_index)
pt115=si(115,30,im_index)
pt171=si(171,30,im_index)
pt172=si(172,30,im_index)
pt173=si(173,30,im_index)
pt174=si(174,30,im_index)
pt175=si(175,30,im_index)
pt176=si(176,30,im_index)
pt177=si(177,30,im_index)
pt178=si(178,30,im_index)
pt181=si(181,30,im_index)
pt182=si(182,30,im_index)
pt183=si(183,30,im_index)
pt184=si(184,30,im_index)
pt185=si(185,30,im_index)
pt186=si(186,30,im_index)
pt187=si(187,30,im_index)
pt188=si(188,30,im_index)
pt189=si(189,30,im_index)
pt190=si(190,30,im_index)
pt191=si(191,30,im_index)
pt192=si(192,30,im_index)
pt193=si(193,30,im_index)
pt194=si(194,30,im_index)
pt195=si(195,30,im_index)
pt196=si(196,30,im_index)
pt197=si(197,30,im_index)
pt198=si(198,30,im_index)
pt199=si(199,30,im_index)
pt200=si(200,30,im_index)
pt201=si(201,30,im_index)
pt202=si(202,30,im_index)
pt203=si(203,30,im_index)
pt204=si(204,30,im_index)
pt205=si(205,30,im_index)
pt206=si(206,30,im_index)
pt207=si(207,30,im_index)
pt208=si(208,30,im_index)
pt293=SI(293,30,im_index)
pt294=SI(294,30,im_index)
pt295=SI(295,30,im_index)
pt296=SI(296,30,im_index)
pt297=SI(297,30,im_index)
pt298=SI(298,30,im_index)
pt301=si(301,30,im_index)
pt302=si(302,30,im_index)
pt303=si(303,30,im_index)
pt304=SI(304,30,im_index)
pt310=si(310,30,im_index)
pt311=si(311,30,im_index)
pt312=si(312,30,im_index)
pt313=si(313,30,im_index)
pt314=si(314,30,im_index)
pt315=si(315,30,im_index)
pt316=si(316,30,im_index)
pt401=si(401,30,im_index)
pt402=si(402,30,im_index)
pt403=si(403,30,im_index)
pt404=si(404,30,im_index)
pt405=si(405,30,im_index)
pt406=si(406,30,im_index)
pt407=si(407,30,im_index)
pt408=si(408,30,im_index)
pt409=si(409,30,im_index)
pt410=si(410,30,im_index)
pt411=si(411,30,im_index)
pt412=si(412,30,im_index)
pt413=si(413,30,im_index)
pt414=si(414,30,im_index)
pt415=si(415,30,im_index)
pt416=si(416,30,im_index)
pt417=si(417,30,im_index)
pt418=si(418,30,im_index)
pt419=si(419,30,im_index)
pt420=si(420,30,im_index)
pt421=si(421,30,im_index)
pt422=si(422,30,im_index)
pt423=si(423,30,im_index)
pt424=si(424,30,im_index)
pt425=si(425,30,im_index)
pt426=si(426,30,im_index)
pt427=si(427,30,im_index)
pt428=si(428,30,im_index)
pt429=si(429,30,im_index)
pt430=si(430,30,im_index)
pt431=si(431,30,im_index)
pt432=si(432,30,im_index)
pt433=si(433,30,im_index)
pt434=si(434,30,im_index)
pt435=si(435,30,im_index)
pt436=si(436,30,im_index)
pt437=si(437,30,im_index)
pt438=si(438,30,im_index)
pt439=si(439,30,im_index)
pt440=si(440,30,im_index)
pt441=si(441,30,im_index)
pt442=si(442,30,im_index)
pt451=si(451,30,im_index)
pt452=si(452,30,im_index)
pt453=si(453,30,im_index)
pt454=si(454,30,im_index)
! FV-TRACK output for cyclone tracking
pt455=SI(455,30,im_index)
pt456=SI(456,30,im_index)
pt457=SI(457,30,im_index)
pt458=SI(458,30,im_index)
pt459=SI(459,30,im_index)
pt460=SI(460,30,im_index)
! Column total diagnostics
pt461=SI(461,30,im_index)
pt462=SI(462,30,im_index)
pt463=SI(463,30,im_index)
! Scalar diagnostic
pt464=si(464,30,im_index)
pt465=si(465,30,im_index)
pt466=si(466,30,im_index)
pt467=si(467,30,im_index)
! EG Prognostic Increments
pt901=si(901,30,im_index)
pt902=si(902,30,im_index)
pt903=si(903,30,im_index)

! Initialise STASHWORK array because EOT_DIAG does not initialise halos
!  DIR$ CACHE_BYPASS STASHWORK
DO i=1,int30
  stashwork(i)=0.0
END DO

dry_atmosphere_mass = a_realhd(rh_tot_mass_init)
CALL eot_diag(                                                          &
! Primary data: in
        pstar,p,rho,u,v,w                                               &
        ,theta,q,qcl,qcf,co2,dryrho                                     &
        ,p_theta_levels ,exner_rho_levels                               &
        ,exner_theta_levels,energy_corr_now                             &
        ,eot_inc_u, eot_inc_v, eot_inc_w, eot_inc_t                     &
        ,eot_inc_q, eot_inc_qcl, eot_inc_qcf, eot_inc_rho               &
        ,eot_inc_qrain, eot_inc_qgraup, eot_inc_qcf2                    &
        ,eot_inc_cf, eot_inc_cfl, eot_inc_cff                           &
        ,eot_inc_m_v, eot_inc_m_cl, eot_inc_m_cf                        &
        ,eot_inc_m_r, eot_inc_m_gr, eot_inc_m_cf2                       &
        ,eot_inc_theta, eot_inc_thetav, eot_inc_dryrho                  &
        ,wet_to_dry_n                                                   &
        ,co2_mmr, dry_atmosphere_mass                                   &
! Grid sizes and definition: in
        ,rows,n_rows,row_length,bl_levels                               &
        ,theta_field_size,u_field_size,v_field_size                     &
        ,sin_theta_latitude,cos_theta_latitude                          &
        ,sin_v_latitude                                                 &
        ,sin_theta_longitude,cos_theta_longitude                        &
        ,delta_lambda ,delta_phi,global_row_length                      &
! Pressure levels for output arrays: in
        ,u_press,v_press,w_press,t_press,q_press,rh_press               &
        ,z_press,om_press,w_presst,t_presst,q_presst                    &
        ,rh_presst,z_presst,om_presst,heavy_press,heavy_presst          &
        ,tem_press,tc_track_press                                       &
! Flags to request each diagnostic output field: in
        ,sf(001,30),sf(002,30),sf(003,30),sf(004,30),sf(005,30)         &
        ,sf(006,30),sf(007,30),sf(008,30)                               &
        ,sf(101,30),sf(102,30),sf(103,30),sf(104,30),sf(105,30)         &
        ,sf(106,30),sf(107,30)                                          &
        ,sf(111,30),sf(112,30),sf(114,30),sf(113,30),sf(115,30)         &
        ,sf(181,30),sf(182,30),sf(183,30),sf(184,30),sf(185,30)         &
        ,sf(186,30),sf(187,30),sf(188,30),sf(189,30),sf(190,30)         &
        ,sf(191,30),sf(192,30),sf(193,30),sf(194,30),sf(195,30)         &
        ,sf(196,30),sf(197,30),sf(198,30),sf(199,30),sf(200,30)         &
        ,sf(171,30),sf(172,30),sf(173,30),sf(174,30),sf(175,30)         &
        ,sf(176,30),sf(177,30),sf(178,30)                               &
        ,sf(201,30),sf(202,30),sf(203,30)                               &
        ,sf(204,30),sf(205,30),sf(206,30),sf(207,30),sf(208,30)         &
        ,sf(293,30),sf(294,30),sf(295,30),sf(296,30),sf(297,30)         &
        ,sf(298,30),sf(301,30),sf(302,30),sf(303,30),sf(304,30)         &
        ,sf(310,30),sf(311,30),sf(312,30),sf(313,30),sf(314,30)         &
        ,sf(315,30),sf(316,30)                                          &
        ,sf(401,30),sf(402,30),sf(403,30),sf(404,30),sf(405,30)         &
        ,sf(406,30),sf(407,30),sf(408,30),sf(409,30),sf(410,30)         &
        ,sf(411,30),sf(412,30),sf(413,30),sf(414,30),sf(415,30)         &
        ,sf(416,30),sf(417,30),sf(418,30),sf(419,30),sf(420,30)         &
        ,sf(421,30),sf(422,30),sf(423,30),sf(424,30),sf(425,30)         &
        ,sf(426,30),sf(427,30),sf(428,30),sf(429,30),sf(430,30)         &
        ,sf(431,30),sf(432,30),sf(433,30),sf(434,30),sf(435,30)         &
        ,sf(436,30),sf(437,30),sf(438,30),sf(439,30)                    &
        ,sf(440,30),sf(441,30)                                          &
        ,sf(451,30),sf(452,30),sf(453,30),sf(454,30)                    &
        ,sf(442,30),sf(455,30),sf(457,30),sf(459,30)                    &
        ,sf(461,30),sf(462,30),sf(463,30),sf(464,30),sf(465,30)         &
        ,sf(466,30),sf(467,30),sf(901,30),sf(902,30),sf(903,30)         &
! Diagnostics lengths: in
        ,u_m_levs,v_m_levs,w_m_levs,t_m_levs,q_m_levs,z_m_levs          &
        ,ke_m_levs,om_m_levs,u_mm_levs,v_mm_levs,w_mm_levs,t_mm_levs    &
        ,q_mm_levs,z_mm_levs,ke_mm_levs,teot_m_levs,wbig_m_levs         &
        ,wbig2_m_levs,RH_m_levs,dry_mass_m_levs                         &
        ,t_inc_levs,q_inc_levs,qcl_inc_levs,qcf_inc_levs                &
        ,qrain_inc_levs,qgraup_inc_levs,qqcf2_inc_levs                  &
        ,cf_inc_levs,cfl_inc_levs,cff_inc_levs                          &
        ,m_v_inc_levs,m_cl_inc_levs,m_cf_inc_levs,m_r_inc_levs          &
        ,m_gr_inc_levs,m_cf2_inc_levs                                   &
        ,u_inc_levs,v_inc_levs,w_inc_levs,rho_inc_levs                  &
        ,theta_inc_levs,thetav_inc_levs,dryrho_inc_levs                 &
        ,t2_inc_levs,q2_inc_levs,qcl2_inc_levs,qcf2_inc_levs            &
        ,u2_inc_levs,v2_inc_levs,w2_inc_levs,rho2_inc_levs              &
        ,u_p_levs,v_p_levs,w_p_levs, t_p_levs,q_p_levs,rh_p_levs        &
        ,z_p_levs,om_p_levs,w_pt_levs,t_pt_levs,q_pt_levs,rh_pt_levs    &
        ,z_pt_levs,om_pt_levs,heavy_p_levs,heavy_pt_levs                &
        ,tv_p_levs,tvom_p_levs                                          &
        ,tem_p_levs                                                     &
        ,tc_track_p_levs, prod_p_levs                                   &
                     ! pressure levels required for each product
        ,npress_diags,nmodel_diags                                      &
! Tropopause index bounds
        ,min_trop_level,max_trop_level                                  &
! Model levels indicies: in
        ,u_m_list,v_m_list,w_m_list,t_m_list,q_m_list,z_m_list,ke_m_list&
        ,om_m_list,u_mm_list,v_mm_list,w_mm_list,t_mm_list,q_mm_list    &
        ,z_mm_list,ke_mm_list,teot_m_list,wbig_m_list,wbig2_m_list      &
        ,RH_m_list,dry_mass_m_list                                      &
        ,t_inc_list,q_inc_list,qcl_inc_list,qcf_inc_list                &
        ,qrain_inc_list,qgraup_inc_list,qqcf2_inc_list                  &
        ,cf_inc_list,cfl_inc_list,cff_inc_list                          &
        ,m_v_inc_list,m_cl_inc_list,m_cf_inc_list,m_r_inc_list          &
        ,m_gr_inc_list,m_cf2_inc_list                                   &
        ,u_inc_list,v_inc_list,w_inc_list,rho_inc_list                  &
        ,theta_inc_list,thetav_inc_list,dryrho_inc_list                 &
        ,t2_inc_list,q2_inc_list,qcl2_inc_list,qcf2_inc_list            &
        ,u2_inc_list,v2_inc_list,w2_inc_list,rho2_inc_list              &
        ,prod_m_list                                                    &
                     ! index of model levels required for product
! Product levels: in
        ,prod_ind                                                       &
                  ! index of pressure levels required for product
! Diagnostic arrays: out
        ,si(1,30,im_index) ,Stashwork                                   &
        ,Stashwork(pt001),Stashwork(pt002),Stashwork(pt003)             &
        ,Stashwork(pt004),Stashwork(pt005),Stashwork(pt006)             &
        ,Stashwork(pt007),Stashwork(pt008)                              &
        ,Stashwork(pt101),Stashwork(pt102),Stashwork(pt103)             &
        ,Stashwork(pt104),Stashwork(pt105),Stashwork(pt106)             &
        ,Stashwork(pt107),Stashwork(pt111),Stashwork(pt112)             &
        ,Stashwork(pt114),Stashwork(pt113),Stashwork(pt115)             &
        ,Stashwork(pt181),Stashwork(pt182),Stashwork(pt183)             &
        ,Stashwork(pt184),Stashwork(pt185),Stashwork(pt186)             &
        ,Stashwork(pt187),Stashwork(pt188),Stashwork(pt189)             &
        ,Stashwork(pt190),Stashwork(pt191),Stashwork(pt192)             &
        ,Stashwork(pt193),Stashwork(pt194),Stashwork(pt195)             &
        ,Stashwork(pt196),Stashwork(pt197),Stashwork(pt198)             &
        ,Stashwork(pt199),Stashwork(pt200)                              &
        ,Stashwork(pt171),Stashwork(pt172),Stashwork(pt173)             &
        ,Stashwork(pt174),Stashwork(pt175),Stashwork(pt176)             &
        ,Stashwork(pt177),Stashwork(pt178)                              &
        ,Stashwork(pt201),Stashwork(pt202),Stashwork(pt203)             &
        ,Stashwork(pt204),Stashwork(pt205),Stashwork(pt206)             &
        ,Stashwork(pt207),Stashwork(pt208)                              &
        ,Stashwork(pt293),Stashwork(pt294),Stashwork(pt295)             &
        ,Stashwork(pt296),Stashwork(pt297),Stashwork(pt298)             &
        ,Stashwork(pt301),Stashwork(pt302),Stashwork(pt303)             &
        ,Stashwork(pt304),Stashwork(pt310),Stashwork(pt311)             &
        ,Stashwork(pt312),Stashwork(pt313),Stashwork(pt314)             &
        ,Stashwork(pt315),Stashwork(pt316)                              &
        ,Stashwork(pt401),Stashwork(pt402),Stashwork(pt403)             &
        ,Stashwork(pt404),Stashwork(pt405),Stashwork(pt406)             &
        ,Stashwork(pt407),Stashwork(pt408),Stashwork(pt409)             &
        ,Stashwork(pt410),Stashwork(pt411),Stashwork(pt412)             &
        ,Stashwork(pt413),Stashwork(pt414),Stashwork(pt415)             &
        ,Stashwork(pt416),Stashwork(pt417),Stashwork(pt418)             &
        ,Stashwork(pt419),stashwork(pt420),Stashwork(pt421)             &
        ,Stashwork(pt422),stashwork(pt423),Stashwork(pt424)             &
        ,Stashwork(pt425),stashwork(pt426),Stashwork(pt427)             &
        ,Stashwork(pt428),stashwork(pt429),Stashwork(pt430)             &
        ,Stashwork(pt431),stashwork(pt432),Stashwork(pt433)             &
        ,Stashwork(pt434),stashwork(pt435),Stashwork(pt436)             &
        ,Stashwork(pt437),stashwork(pt438),Stashwork(pt439)             &
        ,Stashwork(pt440),Stashwork(pt441)                              &
        ,Stashwork(pt451),Stashwork(pt452),Stashwork(pt453)             &
        ,Stashwork(pt454),Stashwork(pt442)                              &
! FV-TRACK STASH diagnostics for cyclone tracking
        ,Stashwork(pt455),Stashwork(pt456),Stashwork(pt457)             &
        ,Stashwork(pt458),Stashwork(pt459),Stashwork(pt460)             &
        ,Stashwork(pt461),Stashwork(pt462),Stashwork(pt463)             &
        ,Stashwork(pt464),Stashwork(pt465),Stashwork(pt466)             &
        ,Stashwork(pt467)                                               &
! EG Prognostic Increments
        ,Stashwork(pt901),Stashwork(pt902),Stashwork(pt903))

IF (icode /= 0) THEN
  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! DEPENDS ON: stash
CALL stash(atmos_sm,atmos_im,30,stashwork,                            &
           icode,cmessage)

9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE st_diag3
