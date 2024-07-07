! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: gravity wave drag

MODULE diagnostics_gwd_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'DIAGNOSTICS_GWD_MOD'
CONTAINS

SUBROUTINE diagnostics_gwd(                                       &
                       row_length, rows, model_levels             &
,                      n_rows, u_rows, v_rows                     &
,                      off_x, off_y                               &
,                      at_extremity                               &
,                      u, v, R_u, R_v, T_inc                      &
,                      u_incr_diagnostic,v_incr_diagnostic        &
,                      t_incr_diagnostic                          &
,                      exner_theta_levels                         &
,                      stress_ud     ,  stress_vd                 &
,                      stress_ud_satn,  stress_vd_satn            &
,                      stress_ud_wake,  stress_vd_wake            &
,                      du_dt_satn    ,  dv_dt_satn                &
,                      du_dt_wake   ,  dv_dt_wake                 &
,                      u_s_d, v_s_d, nsq_s_d                      &
,                      num_lim_d, num_fac_d                       &
,                      fr_d, bld_d ,bldt_d                        &
, tausx_d, tausy_d, taus_scale_d                                  &
, sd_orog_land, land_sea_mask, land_points                        &
, orog_grad_xx_land, orog_grad_xy_land, orog_grad_yy_land         &
, orog_slope_d, orog_anis_d, orog_dir_d                           &
, gwspec_eflux, gwspec_sflux, gwspec_wflux                        &
, gwspec_nflux, gwspec_ewacc, gwspec_nsacc,                       &
     STASHwork                                                    &
     )


! Purpose:
!          Calculates diagnostics and outputs them.
!
! Method:
!   Item by item copying of diagnostic arrays made available from gwd
! scheme into the corresponding STASHwork (master array for output)
! diagnostic space. All sizes and address pointers have been
! calculated by the STASH initialisation routines previously. Each
! diagnostic is protected by a STASHflag logical switch which is
! re-set every timestep.
!   Note that diagnostics for du_dt,dv_dt are available on all
! model rho levels. Stress diagnostics are now available on all model
! theta_levels - running from level 0 at the ground to
! level=model_levels at the model lid.
! New single level diagnostics (6214-6222) are output at theta points
! in the horizontal, i.e. they are output at the points on which the
! GWD scheme operates. 6233 and 6234 are the same.
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:                                   &
    udims, vdims, wdims, tdims, pdims,                             &
    udims_s, vdims_s, wdims_s, tdims_s, pdims_s
! The following are for PWS diagnostics which depend on GWD diagnostics
USE pws_diags_mod, ONLY: pws_gwd_stress_lev_u, pws_gwd_stress_lev_v, &
                         flag_wafc_cat, flag_mtn_wave_turb
USE mask_compression, ONLY: expand_from_mask
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY:                                                     &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!
!     Fixed starting theta-level (e.g. for STASH diagnostic arrays)
INTEGER, PARAMETER :: tkfix1start    = 1

! Arguments with Intent IN. ie: Input variables.

LOGICAL ::                                                        &
  at_extremity(4)  ! Indicates if this processor is at north,
                   ! south, east or west of the processor grid

INTEGER ::                                                        &
  row_length                                                      &
                   ! number of points on a row
, rows                                                            &
                   ! number of rows in a theta field
, n_rows                                                          &
                   ! number of rows in a v field
, model_levels                                                    &
                   ! number of model levels
, number_format                                                   &
                   ! switch controlling number format diagnostics
                   ! are written out in. See PP_WRITE for details.
, land_points               ! Number of land points

INTEGER ::                                                        &
  off_x                                                           &
                      !IN. size of small halo in x direction
, off_y               !IN. size of small halo in y direction

INTEGER ::             &
   u_rows              & ! rows on u grid
,  v_rows                ! rows on v grid

INTEGER ::                                                        &
  INDEX

REAL ::                                                           &
  exner  ! Exner pressure surface value


! Primary Arrays used in all models
REAL       ::                           &!intent(inout):
 r_u(udims_s%i_start:udims_s%i_end,    & !u wind increment
     udims_s%j_start:udims_s%j_end,    &
     udims_s%k_start:udims_s%k_end)    &
,r_v(vdims_s%i_start:vdims_s%i_end,    & !v wind increment
     vdims_s%j_start:vdims_s%j_end,    &
     vdims_s%k_start:vdims_s%k_end)    &
,T_inc(tdims%i_start:tdims%i_end,      & !temperature increment
       tdims%j_start:tdims%j_end,      &
       1:tdims%k_end)

REAL ::                              &!block for intent(inout)
u(udims_s%i_start:udims_s%i_end,    & ! primary u field (ms**-1)
  udims_s%j_start:udims_s%j_end,    &
  udims_s%k_start:udims_s%k_end)    &
,v(vdims_s%i_start:vdims_s%i_end,    & ! primary v field (ms**-1)
  vdims_s%j_start:vdims_s%j_end,    &
  vdims_s%k_start:vdims_s%k_end)

REAL ::                                            &
    exner_theta_levels                            &
                (tdims_s%i_start:tdims_s%i_end,   & ! Exner on theta level
                  tdims_s%j_start:tdims_s%j_end,  & !
                  tdims_s%k_start:tdims_s%k_end)

REAL ::                                                           &
  stress_ud      (row_length,  rows,0:model_levels)               &
, stress_vd      (row_length,n_rows,0:model_levels)               &
, stress_ud_satn (row_length,  rows,0:model_levels)               &
, stress_vd_satn (row_length,n_rows,0:model_levels)               &
, stress_ud_wake (row_length,  rows,0:model_levels)               &
, stress_vd_wake (row_length,n_rows,0:model_levels)               &
, du_dt_satn     (row_length,  rows,  model_levels)               &
, dv_dt_satn     (row_length,n_rows,  model_levels)               &
, du_dt_wake     (row_length, rows,   model_levels)               &
, dv_dt_wake     (row_length,n_rows,  model_levels)               &
, gwspec_eflux(udims%i_start:udims%i_end,                         &
               udims%j_start:udims%j_end,                         &
                 tkfix1start:tdims%k_end)  & ! IN Fp in E azimuth
, gwspec_sflux(vdims%i_start:vdims%i_end,                         &
               vdims%j_start:vdims%j_end,                         &
                 tkfix1start:tdims%k_end)  & ! IN Fp in S azimuth
, gwspec_wflux(udims%i_start:udims%i_end,                         &
               udims%j_start:udims%j_end,                         &
                 tkfix1start:tdims%k_end)  & ! IN Fp in W azimuth
, gwspec_nflux(vdims%i_start:vdims%i_end,                         &
               vdims%j_start:vdims%j_end,                         &
                 tkfix1start:tdims%k_end)  & ! IN Fp in N azimuth
, gwspec_ewacc(udims%i_start:udims%i_end,                         &
               udims%j_start:udims%j_end,                         &
               udims%k_start:udims%k_end)  & ! IN Accel of U wind
, gwspec_nsacc(vdims%i_start:vdims%i_end,                         &
               vdims%j_start:vdims%j_end,                         &
               vdims%k_start:vdims%k_end)  & ! IN Accel of V wind
! Note: space is only allocated for _incr_diagnostic arrays in calling
!       routine if coresponding STASHflags are set .true.
      , u_incr_diagnostic(udims%i_start:udims%i_end,                    &
                          udims%j_start:udims%j_end,                    &
                          udims%k_start:udims%k_end)                    &
      , v_incr_diagnostic(vdims%i_start:vdims%i_end,                    &
                          vdims%j_start:vdims%j_end,                    &
                          vdims%k_start:vdims%k_end)                    &
      , t_incr_diagnostic(tdims%i_start:tdims%i_end,                    &
                          tdims%j_start:tdims%j_end,                    &
                          tkfix1start:tdims%k_end)                      &
      , num_lim_d(row_length, rows)                                     &
      , num_fac_d(row_length, rows)                                     &
      , u_s_d  (row_length, rows)                                       &
      , v_s_d  (row_length, rows)                                       &
      , nsq_s_d  (row_length, rows)                                     &
      , fr_d   (row_length, rows)                                       &
      , bld_d  (row_length, rows)                                       &
      , bldt_d (row_length, rows)                                       &
      , tausx_d (row_length, u_rows)                                    &
      , tausy_d (row_length, v_rows)                                    &
      , taus_scale_d (row_length, rows)                                 &
      , orog_slope_d  (row_length,rows)                                 &
      , orog_anis_d (row_length,rows)                                   &
      , orog_dir_d  (row_length,rows)                                   &
      , work_array(row_length, rows)                                    &
      , sd_orog_land(land_points)                                       &
                                  ! orographic std dev on land points
      , orog_grad_xx_land(land_points)                                  &
                                       ! sigma_xx on land points
      , orog_grad_xy_land(land_points)                                  &
                                       ! sigma_xy on land points
      , orog_grad_yy_land(land_points) ! sigma_yy on land points

LOGICAL, INTENT(IN) ::                                            &
  land_sea_mask (row_length, rows)  ! land_sea_mask


! Diagnostic variables
REAL ::                                                           &
 STASHwork(*)     ! STASH workspace

! Local variables
INTEGER ::                                                        &
  i, j, k, l                                                      &
                ! loop indices
,    icode                                                        &
                          ! Return code  =0 Normal exit  >1 Error
 ,item                                                            &
                  ! STASH item
 ,sect            ! STASH section
PARAMETER( sect = 6 ) ! for gravity wave drag

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*) :: RoutineName
PARAMETER ( RoutineName='DIAGNOSTICS_GWD')

INTEGER ::                                                        &
  im_index        ! internal model index

INTEGER ::                                                        &
  gwspec_eflux_p_levels,                                          &
  stress_ud_p_levels,                                             &
  du_dt_satn_p_levels,                                            &
  gwspec_wflux_p_levels,                                          &
  gwspec_ewacc_p_levels

REAL ::                                                           &
  interp_data_3(row_length*rows*model_levels)                     &
, interp_data_3_u(row_length*u_rows*model_levels)                 &
, interp_data_3_v(row_length*v_rows*model_levels)                 &
, sd_orog     (row_length,rows)                                   &
                                 ! orog std dev full field
, orog_grad_xx(row_length,rows)                                   &
                                 !     sigma_xx full field
, orog_grad_xy(row_length,rows)                                   &
                                 !     sigma_xy full field
, orog_grad_yy(row_length,rows)  !     sigma_yy full field


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ND pp codes temporarily kept around :
!
!      PP_code_r                =   1
!      PP_code_p                =   8
!      PP_code_eta              =  10   !Note: currently no PP code
!!                             for eta, so set to 10 which is sigma.
!      PP_code_T                =  16
!      PP_code_msl              = 128
!      PP_code_surf             = 129
!      PP_code_du_dt_satn       =  68
!      PP_code_dv_dt_satn       =  69
!      PP_code_gwd_stress_u     = 161
!      PP_code_gwd_stress_v     = 162
!      PP_code_u_tend           = 975
!      PP_code_v_tend           = 976
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
icode = 0 ! Initialise error status
im_index = 1

! ----------------------------------------------------------------------
! Section 1.  Diagnostic Calculation and output.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 1.1 Diagnostics for GWD
! ----------------------------------------------------------------------

!-----------------------------------------------------------------------
! u,v increment diagnostics= modified - previous
!-----------------------------------------------------------------------

item = 185  ! u increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,udims,u_incr_diagnostic,R_u)
  DO k=1,model_levels
    DO j=udims%j_start,udims%j_end
      DO i=udims%i_start,udims%i_end
        u_incr_diagnostic(i,j,k) = R_u(i,j,k) -                   &
                                      u_incr_diagnostic(i,j,k)
      END DO !i
    END DO  !j
  END DO   !k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       u_incr_diagnostic,                                        &
       row_length,u_rows,model_levels,0,0,0,0, at_extremity,     &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 185)"//cmessage
  END IF

END IF  !  sf(item,sect)

item = 186  ! v increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,vdims,v_incr_diagnostic,R_v)
  DO k=1,model_levels
    DO j=vdims%j_start,vdims%j_end
      DO i=vdims%i_start,vdims%i_end
        v_incr_diagnostic(i,j,k) = R_v(i,j,k) -                     &
                                        v_incr_diagnostic(i,j,k)
      END DO !i
    END DO  !j
  END DO   !k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       v_incr_diagnostic,                                        &
       row_length,v_rows,model_levels,0,0,0,0, at_extremity,     &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 186)"//cmessage
  END IF

END IF  !  sf(item,sect)

!-----------------------------------------------------------------------
! temperature increment diagnostic= modified - previous
!-----------------------------------------------------------------------

item = 181  ! T increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,tdims,t_incr_diagnostic,T_inc)
  DO k=tkfix1start,model_levels
    DO j=tdims%j_start,tdims%j_end
      DO i=tdims%i_start,tdims%i_end
        t_incr_diagnostic(i,j,k) = T_inc(i,j,k) -                 &
                                      t_incr_diagnostic(i,j,k)
      END DO !i
    END DO  !j
  END DO   !k
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
       t_incr_diagnostic,                                        &
       row_length,rows,model_levels,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,sect,item,                                       &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 181)"//cmessage
  END IF

END IF !  sf(item,sect)

! ----------------------------------------------------------------------
!  U component of stress
! ----------------------------------------------------------------------
! Item 201 stress_ud

IF (sf(201,6)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(201,6,im_index)),stress_ud,      &
       row_length,u_rows,model_levels+1,0,0,0,0,                 &
       at_extremity,                                             &
       stlist(1,stindex(1,201,6,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,6,201,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(stress_ud)"
    GO TO 9999
  END IF


END IF

IF (flag_wafc_cat .OR. flag_mtn_wave_turb) THEN
  pws_gwd_stress_lev_u(:,:,:) = stress_ud(:,:,:)
END IF

IF (sf(241,6)) THEN
  stress_ud_p_levels = stash_levels(1,                           &
     -stlist(10,stindex(1,241,6,im_index)))

  ! Interpolate stress_ud on model levels/u points on the c grid onto
  ! pressure levels specified in the namelist domain profile for stash 6115
  ! and uv positions (B grid). Copy interpolated diagnostic into STASHwork
  ! for stash processing.
  ! DEPENDS ON: interp_2_press_at_uv_pos_b_grid
  CALL Interp_2_Press_At_UV_Pos_B_Grid(                          &
       im_index, 241,6, stress_ud,                               &
       row_length, rows, n_rows, stress_ud_p_levels,             &
       model_levels, off_x,off_y,                                &
       exner_theta_levels(:,:,1:tdims%k_end),                    &
                                   STASHwork(si(241,6,im_index)))
END IF



! ----------------------------------------------------------------------
!  V component of stress
! ----------------------------------------------------------------------
! Item 202 stress_vd

IF (sf(202,6)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(202,6,im_index)),stress_vd,      &
       row_length,v_rows,model_levels+1,0,0,0,0,                 &
       at_extremity,                                             &
       stlist(1,stindex(1,202,6,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,6,202,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(stress_vd)"
    GO TO 9999
  END IF

END IF

IF (flag_wafc_cat .OR. flag_mtn_wave_turb) THEN
  pws_gwd_stress_lev_v(:,:,:) = stress_vd(:,:,:)
END IF

! ----------------------------------------------------------------------
!  Orographic standard deviation
! ----------------------------------------------------------------------
! Item 203 Orographic standard deviation

IF (sf(203,6)) THEN

  CALL expand_from_mask(sd_orog,sd_orog_land, land_sea_mask,        &
      row_length*rows,land_points)

  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(203,6,im_index)),                   &
       sd_orog,                                                  &
       row_length,rows,0,0,0,0,at_extremity,                     &
       atmos_im,6,203,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="Error in copydiag(sd_orog)"
    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
!  Orographic sigma_xx field
! ----------------------------------------------------------------------
! Item 204 Orographic sigma_xx field

IF (sf(204,6)) THEN


  CALL expand_from_mask(orog_grad_xx,orog_grad_xx_land,          &
               land_sea_mask,row_length*rows,land_points)

  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(204,6,im_index)),                   &
       orog_grad_xx,                                             &
       row_length,rows,0,0,0,0,at_extremity,                     &
       atmos_im,6,204,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="Error in copydiag(orog_grad_xx)"
    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
!  Orographic sigma_xy field
! ----------------------------------------------------------------------
! Item 205 Orographic sigma_xy field

IF (sf(205,6)) THEN

  CALL expand_from_mask(orog_grad_xy,orog_grad_xy_land,             &
      land_sea_mask,row_length*rows,land_points)

  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(205,6,im_index)),                   &
       orog_grad_xy,                                             &
       row_length,rows,0,0,0,0,at_extremity,                     &
       atmos_im,6,205,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="Error in copydiag(orog_grad_xy)"
    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
!  Orographic sigma_yy field
! ----------------------------------------------------------------------
! Item 206 Orographic sigma_yy field

IF (sf(206,6)) THEN

  CALL expand_from_mask(orog_grad_yy,orog_grad_yy_land,             &
      land_sea_mask,row_length*rows,land_points)

  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(206,6,im_index)),                   &
       orog_grad_yy,                                             &
       row_length,rows,0,0,0,0,at_extremity,                     &
       atmos_im,6,206,                                           &
       icode,cmessage)

  IF (icode  >   0) THEN
    cmessage="Error in copydiag(orog_grad_yy)"
    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
!  U component of satn stress
! ----------------------------------------------------------------------
! Item 223 stress_ud_satn

IF (sf(223,6)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(223,6,im_index)),stress_ud_satn, &
       row_length,u_rows,model_levels+1,0,0,0,0,                 &
       at_extremity,                                             &
       stlist(1,stindex(1,223,6,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,6,223,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(stress_ud_satn)"
    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
!  V component of satn stress
! ----------------------------------------------------------------------
! Item 224 stress_vd_satn

IF (sf(224,6)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(224,6,im_index)),stress_vd_satn, &
       row_length,v_rows,model_levels+1,0,0,0,0,                 &
       at_extremity,                                             &
       stlist(1,stindex(1,224,6,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,6,224,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(stress_vd_satn)"
    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
!  U component  of wake stress
! ----------------------------------------------------------------------
! Item 227 stress_ud_wake

IF (sf(227,6)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(227,6,im_index)),stress_ud_wake, &
       row_length,u_rows,model_levels+1,0,0,0,0,                 &
       at_extremity,                                             &
       stlist(1,stindex(1,227,6,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,6,227,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(stress_ud_wake)"
    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
!  V component of wake stress
! ----------------------------------------------------------------------
! Item 228 stress_vd_wake

IF (sf(228,6)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(228,6,im_index)),stress_vd_wake, &
       row_length,v_rows,model_levels+1,0,0,0,0,                 &
       at_extremity,                                             &
       stlist(1,stindex(1,228,6,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,6,228,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(stress_vd_wake)"
    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
! du/dt saturation component
! ----------------------------------------------------------------------
! Item 207 du_dt_satn

IF (sf(207,6)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(207,6,im_index)),du_dt_satn,     &
       row_length,u_rows,model_levels,0,0,0,0,                   &
       at_extremity,                                             &
       stlist(1,stindex(1,207,6,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,6,207,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(du_dt_satn)"
    GO TO 9999
  END IF

END IF


IF (sf(247,6)) THEN
  du_dt_satn_p_levels = stash_levels(1,                          &
     -stlist(10,stindex(1,247,6,im_index)))

  ! Interpolate du_dt_satn on model levels/u points on the c grid onto
  ! pressure levels specified in the namelist domain profile for stash 6115
  ! and uv positions (B grid). Copy interpolated diagnostic into STASHwork
  ! for stash processing.
  ! DEPENDS ON: interp_2_press_at_uv_pos_b_grid
  CALL Interp_2_Press_At_UV_Pos_B_Grid(                          &
       im_index, 247,6, du_dt_satn,                              &
       row_length, rows, n_rows, du_dt_satn_p_levels,            &
       model_levels, off_x,off_y,                                &
       exner_theta_levels(:,:,1:tdims%k_end),                    &
                                   STASHwork(si(247,6,im_index)))
END IF



! ----------------------------------------------------------------------
!  dv/dt saturation component
! ----------------------------------------------------------------------
! Item 208 dv_dt_satn

IF (sf(208,6)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(208,6,im_index)),dv_dt_satn,     &
       row_length,v_rows,model_levels,0,0,0,0,                   &
       at_extremity,                                             &
       stlist(1,stindex(1,208,6,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,6,208,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(dv_dt_satn)"
    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
!  du/dt wake component
! ----------------------------------------------------------------------
! Item 231 du_dt_wake

IF (sf(231,6)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(231,6,im_index)),du_dt_wake,     &
       row_length,u_rows,model_levels,0,0,0,0,                   &
       at_extremity,                                             &
       stlist(1,stindex(1,231,6,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,6,231,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(du_dt_wake)"
    GO TO 9999
  END IF

END IF


! ----------------------------------------------------------------------
!  dv/dt wake component
! ----------------------------------------------------------------------
! Item 232 dv_dt_wake

IF (sf(232,6)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(232,6,im_index)),dv_dt_wake,     &
       row_length,v_rows,model_levels,0,0,0,0,                   &
       at_extremity,                                             &
       stlist(1,stindex(1,232,6,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,6,232,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(dv_dt_wake)"
    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
!  u latest.
! ----------------------------------------------------------------------
! Item 002 u

IF (sf(002,6)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i,l)                                                        &
!$OMP SHARED(model_levels,u_rows,row_length,interp_data_3_u,u,R_u)
  DO k = 1, model_levels
    DO j = 1, u_rows
      DO i = 1, row_length
        l = i + (j-1)*row_length + (k-1)*u_rows*row_length
        interp_data_3_u(l) = u(i,j,k) + R_u(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(002,6,im_index)),interp_data_3_u, &
       row_length,u_rows,model_levels,0,0,0,0,                    &
       at_extremity,                                              &
       stlist(1,stindex(1,002,6,im_index)),len_stlist,            &
       stash_levels,num_stash_levels+1,                           &
       atmos_im,6,002,                                            &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(u+R_u)"
    GO TO 9999
  END IF

END IF



! ----------------------------------------------------------------------
!  v latest
! ----------------------------------------------------------------------
! Item 003 v

IF (sf(003,6)) THEN

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i,l)                                                        &
!$OMP SHARED(model_levels,v_rows,row_length,interp_data_3_v,v,R_v)
  DO k = 1, model_levels
    DO j = 1, v_rows
      DO i = 1, row_length
        l = i + (j-1)*row_length + (k-1)*v_rows*row_length
        interp_data_3_v(l) = v(i,j,k) + R_v(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(003,6,im_index)),interp_data_3_v, &
       row_length,v_rows,model_levels,0,0,0,0,                    &
       at_extremity,                                              &
       stlist(1,stindex(1,003,6,im_index)),len_stlist,            &
       stash_levels,num_stash_levels+1,                           &
       atmos_im,6,003,                                            &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(v+R_v)"
    GO TO 9999
  END IF

END IF


! ----------------------------------------------------------------------
!  U_S diag
! ----------------------------------------------------------------------
! Item 214 U_S_D

IF (sf(214,6)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(214,6,im_index)),                   &
       u_s_d,                                                    &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,6,214,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag(u_s)"
    GO TO 9999
  END IF

END IF
! ----------------------------------------------------------------------
!  V_S diag
! ----------------------------------------------------------------------
! Item 215 V_S_d

IF (sf(215,6)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(215,6,im_index)),                   &
       v_s_d,                                                    &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,6,215,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag(v_s)"
    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
!  NSQ_S diag
! ----------------------------------------------------------------------
! Item 216 nsq_s_d

IF (sf(216,6)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(216,6,im_index)),                   &
       nsq_s_d,                                                  &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,6,216,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag(nsq_s)"
    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
!  Froude number diag
! ----------------------------------------------------------------------
! Item 217 FR_d

IF (sf(217,6)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(217,6,im_index)),                   &
       fr_d,                                                     &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,6,217,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag(Fr)"
    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
!  Blocked Layer Depth diag
! ----------------------------------------------------------------------
! Item 218 BLD_d

IF (sf(218,6)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(218,6,im_index)),                   &
       bld_d,                                                    &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,6,218,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag(BLD)"
    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
!  % of time of blocked layer
! ----------------------------------------------------------------------
! Item 222 BLDT_D

IF (sf(222,6)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(222,6,im_index)),                   &
       bldt_d,                                                   &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,6,222,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag(BLDT)"

    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
!  % of time numerical limiter invoked (4A only)
! ----------------------------------------------------------------------
! Item 233 NUM_LIM_D

IF (sf(233,6)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(233,6,im_index)),                   &
       num_lim_d,                                                &
   row_length,rows,0,0,0,0, at_extremity,                        &
       atmos_im,6,233,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag(num_lim)"

    GO TO 9999
  END IF

END IF



! ----------------------------------------------------------------------
!  % redn. of flow-blocking stress after numerical limiter invoked (4A only)
! ----------------------------------------------------------------------
! Item 234 NUM_FAC_D

IF (sf(234,6)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(234,6,im_index)),                   &
       num_fac_d,                                                &
   row_length,rows,0,0,0,0, at_extremity,                        &
       atmos_im,6,234,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag(num_fac)"

    GO TO 9999
  END IF

END IF


! ----------------------------------------------------------------------
!  Orographic X-component of the total surface stress
!  Same as 6201, level 1 output, but here the pp headers are correctly
!  labelled
! ----------------------------------------------------------------------
! Item 235 TAUSX_D

IF (sf(235,6)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(235,6,im_index)),                   &
       tausx_d,                                                  &
       row_length,u_rows,0,0,0,0, at_extremity,                  &
       atmos_im,6,235,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag(tausx)"

    GO TO 9999
  END IF

END IF


! ----------------------------------------------------------------------
!  Orographic Y-component of the total surface stress
!  Same as 6202, level 1 output, but here the pp headers are correctly
!  labelled
! ----------------------------------------------------------------------
! Item 236 TAUSY_D

IF (sf(236,6)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(236,6,im_index)),                   &
       tausy_d,                                                  &
       row_length,v_rows,0,0,0,0, at_extremity,                  &
       atmos_im,6,236,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag(tausy)"

    GO TO 9999
  END IF

END IF


! ----------------------------------------------------------------------
!  Orographic surface stress Froude number dependent scaling factor (4A only)
! ----------------------------------------------------------------------
! Item 237 TAUS_SCALE_D

IF (sf(237,6)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(237,6,im_index)),                   &
       taus_scale_d,                                             &
   row_length,rows,0,0,0,0, at_extremity,                        &
       atmos_im,6,237,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag(taus_scale)"

    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
!  Slope of sub-grid scale orography (5A scheme only)
! ----------------------------------------------------------------------
! Item 248 orog_slope_d

IF (sf(248,6)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(248,6,im_index)),                   &
       orog_slope_d,                                             &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,6,248,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag(OROG_SLOPE_D)"

    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
!  Anisotropy of sub-grid scale orography (5A scheme only)
! ----------------------------------------------------------------------
! Item 249 orog_anis_d

IF (sf(249,6)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(249,6,im_index)),                   &
       orog_anis_d,                                              &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,6,249,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag(OROG_ANIS_D)"

    GO TO 9999
  END IF

END IF

! ----------------------------------------------------------------------
!  Orientation of sub-grid scale orography (5A scheme only)
! ----------------------------------------------------------------------
! Item 250 orog_dir_d

IF (sf(250,6)) THEN
  ! DEPENDS ON: copydiag
  CALL copydiag(STASHwork(si(250,6,im_index)),                   &
       orog_dir_d,                                               &
       row_length,rows,0,0,0,0, at_extremity,                    &
       atmos_im,6,250,                                           &
       icode,cmessage)


  IF (icode  >   0) THEN
    cmessage="Error in copydiag(OROG_DIR_D)"

    GO TO 9999
  END IF

END IF

!---------------------------------------------------------------------
! Section 2.
! SPECTRAL (NON-OROGRAPHIC) GRAVITY WAVE SCHEME DIAGNOSTICS
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!  USSP Eastward flux component
! --------------------------------------------------------------------
! Item 101 GWSPEC_EFLUX(i,j,k)
IF (sf(101,6)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(101,6,im_index)),gwspec_eflux,   &
       row_length,rows,model_levels,0,0,0,0,                     &
       at_extremity,                                             &
       stlist(1,stindex(1,101,6,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,6,101,                                           &
       icode,cmessage)
  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(GWSPEC_EFLUX)"
    GO TO 9999
  END IF
END IF

IF (sf(111,6)) THEN
  gwspec_eflux_p_levels = stash_levels(1,                        &
     -stlist(10,stindex(1,111,6,im_index)))

  ! Interpolate GWSPEC_EFLUX on model levels/u points on the c grid onto
  ! pressure levels specified in the namelist domain profile for stash 6115
  ! and uv positions (B grid). Copy interpolated diagnostic into STASHwork
  ! for stash processing.
  ! DEPENDS ON: interp_2_press_at_uv_pos_b_grid
  CALL Interp_2_Press_At_UV_Pos_B_Grid(                          &
       im_index, 111,6, gwspec_eflux,                            &
       row_length, rows, n_rows, gwspec_eflux_p_levels,          &
       model_levels, off_x,off_y,                                &
       exner_theta_levels(:,:,1:tdims%k_end),                    &
                                   STASHwork(si(111,6,im_index)))
END IF
! --------------------------------------------------------------------
!  USSP Southward flux component
! --------------------------------------------------------------------
! Item 102 GWSPEC_SFLUX(i,j,k)
IF (sf(102,6)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(102,6,im_index)),gwspec_sflux,   &
       row_length,n_rows,model_levels,0,0,0,0,                   &
       at_extremity,                                             &
       stlist(1,stindex(1,102,6,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,6,102,                                           &
       icode,cmessage)
  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(GWSPEC_SFLUX)"
    GO TO 9999
  END IF
END IF
! --------------------------------------------------------------------
!  USSP Westward flux component
! --------------------------------------------------------------------
! Item 103 GWSPEC_WFLUX(i,j,k)
IF (sf(103,6)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(103,6,im_index)),gwspec_wflux,   &
       row_length,rows,model_levels,0,0,0,0,                     &
       at_extremity,                                             &
       stlist(1,stindex(1,103,6,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,6,103,                                           &
       icode,cmessage)
  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(GWSPEC_WFLUX)"
    GO TO 9999
  END IF
END IF

IF (sf(113,6)) THEN
  gwspec_wflux_p_levels = stash_levels(1,                        &
     -stlist(10,stindex(1,113,6,im_index)))

  ! Interpolate GWSPEC_WFLUX on model levels/u points on the c grid onto
  ! pressure levels specified in the namelist domain profile for stash 6115
  ! and uv positions (B grid). Copy interpolated diagnostic into STASHwork
  ! for stash processing.
  ! DEPENDS ON: interp_2_press_at_uv_pos_b_grid
  CALL Interp_2_Press_At_UV_Pos_B_Grid(                          &
       im_index, 113,6, gwspec_wflux,                            &
       row_length, rows, n_rows, gwspec_wflux_p_levels,          &
       model_levels, off_x,off_y,                                &
       exner_theta_levels(:,:,1:tdims%k_end),                    &
                                   STASHwork(si(113,6,im_index)))
END IF



! --------------------------------------------------------------------
!  USSP Northward flux component
! --------------------------------------------------------------------
! Item 104 GWSPEC_NFLUX(i,j,k)
IF (sf(104,6)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(104,6,im_index)),gwspec_nflux,   &
       row_length,n_rows,model_levels,0,0,0,0,                   &
       at_extremity,                                             &
       stlist(1,stindex(1,104,6,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,6,104,                                           &
       icode,cmessage)
  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(GWSPEC_NFLUX)"
    GO TO 9999
  END IF
END IF
! --------------------------------------------------------------------
!  USSP EW acceleration component
! --------------------------------------------------------------------
! Item 105 GWSPEC_EWACC(i,j,k)
IF (sf(105,6)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(105,6,im_index)),gwspec_ewacc,   &
       row_length,rows,model_levels,0,0,0,0,                     &
       at_extremity,                                             &
       stlist(1,stindex(1,105,6,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,6,105,                                           &
       icode,cmessage)
  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(GWSPEC_EWACC)"
    GO TO 9999
  END IF
END IF

IF (sf(115,6)) THEN
  gwspec_ewacc_p_levels = stash_levels(1,                        &
     -stlist(10,stindex(1,115,6,im_index)))

  ! Interpolate GWSPEC_EWACC on model levels/u points on the c grid onto
  ! pressure levels specified in the namelist domain profile for stash 6115
  ! and uv positions (B grid). Copy interpolated diagnostic into STASHwork
  ! for stash processing.
  ! DEPENDS ON: interp_2_press_at_uv_pos_b_grid
  CALL Interp_2_Press_At_UV_Pos_B_Grid(                          &
       im_index, 115,6, gwspec_ewacc,                            &
       row_length, rows, n_rows, gwspec_ewacc_p_levels,          &
       model_levels, off_x,off_y,                                &
       exner_theta_levels(:,:,1:tdims%k_end),                    &
                                   STASHwork(si(115,6,im_index)))
END IF



! --------------------------------------------------------------------
!  USSP NS acceleration component
! --------------------------------------------------------------------
! Item 106 GWSPEC_NSACC(i,j,k)
IF (sf(106,6)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d(STASHwork(si(106,6,im_index)),gwspec_nsacc,   &
       row_length,n_rows,model_levels,0,0,0,0,                   &
       at_extremity,                                             &
       stlist(1,stindex(1,106,6,im_index)),len_stlist,           &
       stash_levels,num_stash_levels+1,                          &
       atmos_im,6,106,                                           &
       icode,cmessage)
  IF (icode  >   0) THEN
    cmessage="Error in copydiag_3d(GWSPEC_NSACC)"
    GO TO 9999
  END IF
END IF

9999 CONTINUE                  !  single point error handling.
IF (icode /= 0) THEN

  CALL Ereport(RoutineName,icode,Cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE diagnostics_gwd
END MODULE diagnostics_gwd_mod
