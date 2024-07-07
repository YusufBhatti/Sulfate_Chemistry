! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! SCM convection diagnostics output

MODULE scm_diag_conv_mod

IMPLICIT NONE
!---------------------------------------------------------------------------
! Description: Calls to scmoutput to output convective diagnostics from the SCM

! method:

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection

! code description:
!   language: fortran 90
!   This code is written to umdp3 programming standards version 10.4
!---------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SCM_DIAG_CONV_MOD'

CONTAINS

SUBROUTINE scm_diag_conv( n_cca_lev, rows, row_length,                    &
                          ntml, ntpar,                                    &
                          ccb, cct, lcbase, lctop, conv_type,             &
                          nscmdpkgs, l_scmdiags,                          &
                          cumulus, l_shallow, l_congestus, l_mid_level,   &
                          wstar, wthvs, cclwp, cca_2d, lcca, deep_flag,   &
                          past_precip, past_conv_ht,                      &
                          dubydt_conv_u, dvbydt_conv_v,                   &
                          theta_star, exner_theta_levels,                 &
                          cca, ccw  )

! Convective diagnostic output arrays 
USE cv_diagnostic_array_mod, ONLY:                                      &
        precip_deep, precip_shall ,precip_mid ,precip_cong ,cape_out    &
       ,deep_ind ,shallow_ind ,congestus_ind ,congestus_ind2 ,mid_ind   &
       ,freeze_diag ,kterm_diag                                         &
       ,wstar_up_diag ,wstar_dn_diag ,mb1_diag ,mb2_diag ,wqt_cb        &
       ,wthetal_cb ,wqt_inv ,wthetal_inv ,sh_top ,sh_base ,cg_top       &
       ,cg_base ,cg_term ,cape_ts_diag ,ind_cape_reduced_diag           &
       ,uw_dp, vw_dp ,uw_shall ,vw_shall ,uw_mid ,vw_mid ,up_flux       &
       ,dwn_flux ,entrain_up ,detrain_up ,entrain_dwn ,detrain_dwn      &
       ,up_flux_half                                                    &
       ,T_incr_diag_conv ,q_incr_diag_conv                              &
       ,qcl_incr_diag_conv ,qcf_incr_diag_conv                          &
       ,cf_liquid_incr_diag_conv ,cf_frozen_incr_diag_conv              &
       ,bulk_cf_incr_diag_conv                                          &
       ,mf_deep ,mf_congest ,mf_shall ,mf_midlev                        &
       ,dt_deep ,dt_congest ,dt_shall ,dt_midlev                        &
       ,dq_deep ,dq_congest ,dq_shall ,dq_midlev                        &
       ,du_deep ,du_congest ,du_shall ,du_midlev                        &
       ,dv_deep ,dv_congest ,dv_shall ,dv_midlev                        &
       ,wqt_flux_sh ,wql_flux_sh ,wthetal_flux_sh ,wthetav_flux_sh      &
       ,conv_rain_3d ,conv_snow_3d                                      &
       ,qcl_incr_inhom_diag ,qcf_incr_inhom_diag                        &
       ,bulk_cf_incr_inhom_diag ,cf_liquid_incr_inhom_diag              &
       ,cf_frozen_incr_inhom_diag, deep_cfl_limited, mid_cfl_limited    &
       ,w2p,ep1, ep2, ep3, dt_dd, dq_dd, area_ud, area_dd

USE cv_run_mod,  ONLY:                                                  &   
    i_convection_vn, i_convection_vn_5a, i_convection_vn_6a,            &
    l_mom, iconv_congestus, l_3d_cca


USE cv_stash_flg_mod, ONLY: flg_w_eqn

USE timestep_mod, ONLY: timestep

! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:                                        &
    udims, vdims, wdims, tdims, pdims, tdims_s,                         &
    ScmRowLen, ScmRow

USE s_scmop_mod, ONLY: default_streams,                                &
                       t_inst, t_avg, d_sl, d_wet, d_all, d_point,     &
                       scmdiag_conv, scmdiag_incs
USE scmoutput_mod, ONLY: scmoutput


USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

!----------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN)  ::  &
  n_cca_lev              & ! Number of levels for conv cloud
                           ! amount: 1 for 2D, nlevs for 3D.
 ,rows                   & ! Number of rows
 ,row_length               ! Row length

INTEGER, INTENT(IN) ::          &
  ntml(row_length, rows)        & ! Top level of surface mixed layer
 ,ntpar(row_length, rows)       & ! Top level of initial parcel ascent
 ,ccb    (row_length,rows)      & ! Cnv Cld Base of highest layer on gridpoint
 ,cct    (row_length,rows)      & ! Cnv Cld Top  of highest layer on gridpoint
 ,lcbase (row_length,rows)      & ! Cnv Cld Base of lowest  layer on gridpoint
 ,lctop  (row_length,rows)        ! Cnv Cld Top  of lowest  layer on gridpoint


INTEGER, INTENT(IN) :: conv_type(row_length, rows)
              ! Integer index describing convective type:
              !    0=no convection
              !    1=non-precipitating shallow
              !    2=drizzling shallow
              !    3=warm congestus
              !    ...
              !    8=deep convection

! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN)  ::  &
  nscmdpkgs                ! No of SCM diagnostics packages

LOGICAL, INTENT(IN)  ::  &
  l_scmdiags(nscmdpkgs)    ! Logicals for SCM diagnostics packages

LOGICAL, INTENT(IN) ::             &
  cumulus(row_length, rows)        & ! Logical switch from boundary
 ,l_shallow(row_length, rows)      & ! Logical switch for shallow Cu
 ,l_congestus(row_length, rows)    & ! Logical switch for congestus
 ,l_mid_level(row_length,rows)       ! true if mid level convection

REAL,INTENT(IN) ::               &
  wstar(row_length, rows)        & ! Mixed layer convective velocity scale
 ,wthvs(row_length, rows)        & ! Surface flux of THV
 ,cclwp(row_length,rows)         & ! Cloud Condensed water path (kg/m^2)
 ,cca_2d(row_length,rows)        & ! Cnv.Cld Amount (2d) with no anvil
 ,lcca(row_length,rows)          & ! Conv. Cloud Amount of lowest Convective
                                   ! Cloud Layer on gridpoint. With no anvil
                                   ! scheme modification (0-1)
 ,deep_flag(row_length,rows)     & ! Value between 0.0 and 1.0
                                   ! 1 if deep convection last timestep
 ,past_precip(row_length,rows)   & ! Rate of convective precip (kg/m2/s)
                                   ! decayed over time.
 ,past_conv_ht(row_length,rows)    ! Height of previous convection (m)

! Convective momentum transport tendencies on natives grids
REAL, INTENT(IN) :: dubydt_conv_u ( udims%i_start:udims%i_end,      &
                                    udims%j_start:udims%j_end,      &
                                    udims%k_start:udims%k_end )
REAL, INTENT(IN) :: dvbydt_conv_v ( vdims%i_start:vdims%i_end,      &
                                    vdims%j_start:vdims%j_end,      &
                                    vdims%k_start:vdims%k_end )
! (note that in the SCM these are just copies of the fields on the p-grid).

REAL,INTENT(IN) ::                                   &
  theta_star ( tdims%i_start: tdims%i_end,           & ! Theta after convection
               tdims%j_start: tdims%j_end,           & !
               tdims%k_start: tdims%k_end )          & !
, exner_theta_levels(tdims_s%i_start:tdims_s%i_end,  & ! Exner on theta level
                     tdims_s%j_start:tdims_s%j_end,  & !
                     tdims_s%k_start:tdims_s%k_end)  & !
, cca(row_length,rows,n_cca_lev)                     & ! Cnv.Cld Amount (0-1)
, ccw(tdims%i_end,tdims%j_end,tdims%k_end)             ! Cnv.Cld Water (kg/kg)

!---------------------------------------------------------------------------
! local variables.
!---------------------------------------------------------------------------

! loop counters
INTEGER ::                    &
  i, j, k, ii, jj, iScm, jScm

! Work arrays SCM diagnostic output
REAL :: TmpScm2d(ScmRowLen,ScmRow)
REAL :: TmpScm3d(ScmRowLen,ScmRow,tdims%k_end)

CHARACTER(LEN=*), PARAMETER :: RoutineName='SCM_DIAG_CONV'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------



IF (l_scmdiags(scmdiag_conv) .OR.                                    &
        l_scmdiags(scmdiag_incs)) THEN

  IF (flg_w_eqn) THEN
    !       Stash 5,196
    CALL scmoutput(w2p,'w2p',                                            &
             'w^2 Simpson&Wiggert','(m/s)^2',                            &
             t_avg,d_all,default_streams,'',routinename)

    DO k=1, tdims%k_end
      DO j=tdims%j_start, tdims%j_end
        jScm = j - tdims%j_start + 1
        DO i=tdims%i_start, tdims%i_end
          iScm = i - tdims%i_start + 1
          IF (w2p(i,j,k) > 0.0) THEN
            TmpScm3d(iScm,jScm,k) = SQRT(w2p(i,j,k))
          ELSE
            TmpScm3d(iScm,jScm,k) = 0.0
          END IF
        END DO ! i
      END DO ! j
    END DO ! k
      
    !     Stash 5,197
    CALL scmoutput(TmpScm3d,'wp',                        &
         'w Simpson&Wiggert','(m/s)',                    &
         t_avg,d_all,default_streams,'',routinename)
  END IF
  !-------------------------

  !     Stash 5,163
  CALL scmoutput(qcl_incr_inhom_diag,'dqcl_cinh',          &
       'QCL inc: conv inhom','kg/kg',                      &
       t_avg,d_wet,default_streams,'',routinename)

  !   Stash 5,164
  CALL scmoutput(qcf_incr_inhom_diag,'dqcf_cinh',          &
       'QCF inc: conv inhom','kg/kg',                      &
       t_avg,d_wet,default_streams,'',routinename)

  !   Stash 5,172
  CALL scmoutput(bulk_cf_incr_inhom_diag,'dbcf_cinh',      &
       'Bulk cloud frac inc: conv inhom','Fraction',       &
       t_avg,d_wet,default_streams,'',routinename)

  !   Stash 5,173
  CALL scmoutput(cf_liquid_incr_inhom_diag,'dcfl_cinh',    &
       'Liquid cloud frac inc: conv inhom','Fraction',     &
       t_avg,d_wet,default_streams,'',routinename)

  !   Stash 5,174
  CALL scmoutput(cf_frozen_incr_inhom_diag,'dcff_cinh',    &
       'Frozen cloud frac inc: conv inhom','Fraction',     &
       t_avg,d_wet,default_streams,'',routinename)

  !   Stash 5,181
  CALL scmoutput(t_incr_diag_conv,'dt_conv',               &
       'T increment convection','K',                       &
       t_avg,d_all,default_streams,'',routinename)

  !   Stash 5,182
  CALL scmoutput(q_incr_diag_conv,'dq_conv',               &
       'q increment convection','kg/kg',                   &
       t_avg,d_wet,default_streams,'',routinename)

  !   Stash 5,183
  CALL scmoutput(qcl_incr_diag_conv,'dqcl_conv',          &
       'QCL increment convection','kg/kg',                &
       t_avg,d_wet,default_streams,'',routinename)

  !   Stash 5,184
  CALL scmoutput(qcf_incr_diag_conv,'dqcf_conv',          &
       'QCF increment convection','kg/kg',                &
       t_avg,d_wet,default_streams,'',routinename)

  !  Stash 5,185
  DO k = 1, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        ! Calc u increment from tendency; u-grid same size as T grid in SCM
        TmpScm3d(i,j,k) = dubydt_conv_u(i,j,k) * timestep
      END DO ! i
    END DO ! j
  END DO ! k
  CALL scmoutput(TmpScm3d,'du_conv',     &
       'U increment convection','m/s',           &
       t_avg,d_all,default_streams,'',routinename)

  !  Stash 5,186
  DO k = 1, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        ! Calc u increment from tendency; u-grid same size as T grid in SCM
        TmpScm3d(i,j,k) = dvbydt_conv_v(i,j,k) * timestep
      END DO ! i
    END DO ! j
  END DO ! k
  CALL scmoutput(TmpScm3d,'dv_conv',     &
       'V increment convection','m/s',           &
       t_avg,d_all,default_streams,'',routinename)

  !   Stash 5,192
  CALL scmoutput(bulk_cf_incr_diag_conv,'dbcf_conv',     &
   'Bulk cloud frac increment convection','Fraction',    &
   t_avg,d_wet,default_streams,'',routinename)

  CALL scmoutput(cf_liquid_incr_diag_conv,'dcfl_conv',   &
   'Liquid cloud frac increment convection','Fraction',  &
   t_avg,d_wet,default_streams,'',routinename)

  !  Stash 5,194
  CALL scmoutput(cf_frozen_incr_diag_conv,'dcff_conv',   &
   'Frozen cloud frac increment convection','Fraction',  &
   t_avg,d_wet,default_streams,'',routinename)

  ! Note: for some reason the convective theta increment is
  ! output near the end of atmos_physics2, instead of here.
  ! It would make much more sense to output it here.
  ! However, this would change the order of output of the
  ! SCM diags, so that the KGO checks in rose stem fail.
  ! For now (during development) we want to be able to tell
  ! if we've actually changed the SCM model-evolution, and so
  ! need to preserve the order of scm diag output.  So leave
  ! the scmoutput call for theta_inc in atmos_physics2 for now.
!  CALL SCMoutput(theta_incr_diag_conv,'dth_conv',        &
!       'Convective increment - theta','K',               &
!       t_avg,d_all,default_streams,'',RoutineName)


END IF ! scmdiag_conv .OR. scmdiag_incs

!-----------------------------------------------------------------------
!     SCM Convection Diagnostics Package
!-----------------------------------------------------------------------
IF (l_scmdiags(scmdiag_conv)) THEN

  !       Stash 5,209
  DO k = 1, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        TmpScm3d(i,j,k) = theta_star(i,j,k) * exner_theta_levels(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k
  CALL scmoutput(TmpScm3d,'T_afterconv',                                 &
           'Temperature after convection','K',                           &
           t_avg,d_all,default_streams,'',routinename)

  !   Stash 5,273
  TmpScm2d=REAL(ntml)
  CALL scmoutput(TmpScm2d,'NTML',            &
       'Top level of surface mixed layer','Model level',     &
       t_inst,d_point,default_streams,'',routinename)

  !   Stash 5,274
  TmpScm2d = REAL(ntpar)
  CALL scmoutput(TmpScm2d,'NTPAR',               &
       'Top level of initial parcel ascent','Model level',       &
       t_inst,d_point,default_streams,'',routinename)

  !   Stash 5,275
  CALL scmoutput(freeze_diag,'freeze_lev',       &
       'Freezing level','Model level',           &
       t_inst,d_point,default_streams,'',routinename)

  IF (i_convection_vn == i_convection_vn_6a ) THEN

    DO j=1,rows
      DO i=1,row_length
        TmpScm2d(i,j)=REAL(conv_type(i,j))
      END DO
    END DO

  ELSE

    DO j=1,rows
      DO i=1,row_length
        IF (cumulus(i,j)) THEN
          TmpScm2d(i,j)=1.0
        ELSE
          TmpScm2d(i,j)=0.0
        END IF
      END DO ! i
    END DO ! j

  END IF


  CALL scmoutput(TmpScm2d,'ind_cumulus',                     &
       'Indicator for cumulus convection','Indicator',       &
       t_inst,d_point,default_streams,'',routinename)

  !   Stash 5,267
  CALL scmoutput(deep_cfl_limited,'cfl_limited_deep',        &
       'Indicator for CFL limited deep','Indicator',         &
       t_inst,d_point,default_streams,'',routinename)

  !   Stash 5,268
  CALL scmoutput(mid_cfl_limited,'cfl_limited_mid',      &
       'Indicator for CFL limited mid','Indicator',      &
       t_inst,d_point,default_streams,'',routinename)

  !   Stash 5,269
  CALL scmoutput(deep_ind,'ind_deep',                    &
       'Indicator for deep convection','Indicator',      &
       t_inst,d_point,default_streams,'',routinename)

  !   Stash 5,270
  DO j=1,rows
    DO i=1,row_length
      IF (l_shallow(i,j)) THEN
        TmpScm2d(i,j)=1.0
      ELSE
        TmpScm2d(i,j)=0.0
      END IF
    END DO ! i
  END DO ! j

  CALL scmoutput(TmpScm2d,'ind_shallow',                     &
       'Indicator for shallow convection','Indicator',       &
       t_inst,d_point,default_streams,'',routinename)

  !   Stash 5,272
  DO j=1,rows
    DO i=1,row_length
      IF (l_mid_level(i,j)) THEN
        TmpScm2d(i,j)=1.0
      ELSE
        TmpScm2d(i,j)=0.0
      END IF
    END DO ! i
  END DO ! j

  CALL scmoutput(TmpScm2d,'ind_midconv',                     &
       'Indicator for mid-level convection','Indicator',     &
       t_inst,d_point,default_streams,'',routinename)

  IF (iconv_congestus == 1) THEN
    !    Stash 5,310
    DO j=1,rows
      DO i=1,row_length
        IF (l_congestus(i,j)) THEN
          TmpScm2d(i,j)=1.0
        ELSE
          TmpScm2d(i,j)=0.0
        END IF
      END DO ! i
    END DO ! j

    CALL scmoutput(TmpScm2d,'ind_congest',                   &
     'Indicator for congestus convection','Indicator',       &
     t_inst,d_point,default_streams,'',routinename)
  END IF

  !   Stash 5,276
  TmpScm2d=kterm_diag
  CALL scmoutput(TmpScm2d,'kterm_deep',                      &
       'Deep convection termination level','Model level',    &
       t_inst,d_point,default_streams,'',routinename)

  TmpScm2d = REAL(cct)
  CALL scmoutput(TmpScm2d,'cct',                                     &
       'Convective Cloud Top Level of Highest Convective Layer','-', &
       t_avg,d_sl,default_streams,'',routinename)

  TmpScm2d = REAL(ccb)
  CALL scmoutput(TmpScm2d,'ccb',                                     &
       'Convective Cloud Base Level of Highest Convective Layer','-',&
       t_avg,d_sl,default_streams,'',routinename)

  TmpScm2d = REAL(lcbase)
  CALL scmoutput(TmpScm2d,'lcbase',                                  &
       'Convective Cloud Base Level of Lowest Convective Layer','-', &
       t_avg,d_sl,default_streams,'',routinename)
    
  TmpScm2d = REAL(lctop)
  CALL scmoutput(TmpScm2d,'lctop',                                   &
       'Convective Cloud Top Level of Lowest Convective Layer','-',  &
       t_avg,d_sl,default_streams,'',routinename)

  CALL scmoutput(lcca,'lcca',                                &
       'Convective Cloud Amount at base of Lowest '//        &
       'Convective Cloud Layer (no anvil)','-',              &
       t_avg,d_sl,default_streams,'',routinename)

  !   Stash 5,277
  CALL scmoutput(precip_deep,'precip_deep',              &
       'Deep convective precipitation','kg/m2/s',        &
       t_avg,d_sl,default_streams,'',routinename)

  !   Stash 5,278
  CALL scmoutput(precip_shall,'precip_shall',                &
       'Shallow convective precipitation','kg/m2/s',         &
       t_avg,d_sl,default_streams,'',routinename)

  !   Stash 5,279
  CALL scmoutput(precip_mid,'precip_mid',                    &
       'Mid level convective precipitation','kg/m2/s',       &
       t_avg,d_sl,default_streams,'',routinename)

  IF (iconv_congestus == 1) THEN
    !    Stash 5,280
    CALL scmoutput(precip_cong,'precip_cong',            &
     'Congestus convective precipitation','kg/m2/s',     &
     t_avg,d_sl,default_streams,'',routinename)
  END IF

  !   Stash 5,250
  CALL scmoutput(up_flux,'up_massflux',          &
       'updraught mass flux','Pa/s',             &
       t_avg,d_all,default_streams,'',routinename)

  !   Stash 5,251
  CALL scmoutput(dwn_flux,'down_massflux',           &
       'downdraught mass flux','Pa/s',               &
       t_avg,d_all,default_streams,'',routinename)

  !   Stash 5,252
  CALL scmoutput(entrain_up,'entrain_up',            &
       'updraught entrainment rate','Pa/s',          &
       t_avg,d_all,default_streams,'',routinename)

  !   Stash 5,254
  CALL scmoutput(entrain_dwn,'entrain_dw',           &
       'downdraught entrainment rate','Pa/s',        &
       t_avg,d_all,default_streams,'',routinename)

  !   Stash 5,253
  CALL scmoutput(detrain_up,'detrain_up',            &
       'updraught detrainment rate','Pa/s',          &
       t_avg,d_all,default_streams,'',routinename)

  !   Stash 5,255
  CALL scmoutput(detrain_dwn,'detrain_dw',           &
       'downdraught detrainment rate','Pa/s',        &
       t_avg,d_all,default_streams,'',routinename)

  IF (l_mom) THEN

    !   Stash 5,258
    CALL scmoutput(uw_dp,'uw_dp',                            &
     'x_comp of stress from deep convection','kg/m/s2',      &
     t_avg,d_all,default_streams,'',routinename)

    !   Stash 5,259
    CALL scmoutput(vw_dp,'vw_dp',                            &
     'y_comp of stress from deep convection','kg/m/s2',      &
     t_avg,d_all,default_streams,'',routinename)

    !   Stash 5,260
    CALL scmoutput(uw_shall,'uw_shall',                      &
     'x_comp of stress from shallow convection','kg/m/s2',   &
     t_avg,d_all,default_streams,'',routinename)

    !   Stash 5,261
    CALL scmoutput(vw_shall,'vw_shall',                      &
     'y_comp of stress from shallow convection','kg/m/s2',   &
     t_avg,d_all,default_streams,'',routinename)

    !   Stash 5,249
    CALL scmoutput(up_flux_half,'up_halfmassflux',   &
     'updraught mass flux on half levs','Pa/s',      &
     t_avg,d_all,default_streams,'',routinename)

  END IF

  !   Stash 5,227
  CALL scmoutput(conv_rain_3d,'conv_rain_3d',        &
       'Convective Rainfall Flux','kg/m2/s',         &
       t_avg,d_wet,default_streams,'',routinename)

  !   Stash 5,228
  CALL scmoutput(conv_snow_3d,'conv_snow_3d',        &
       'Convective Snowfall Flux','kg/m2/s',         &
       t_avg,d_wet,default_streams,'',routinename)


  ! Extra added for testing

  CALL scmoutput(wstar,'wstar',                      &
       'Sub cloud convective velocity scale','m/s',  &
       t_avg,d_sl,default_streams,'',routinename)

  CALL scmoutput(wthvs,'wthvs',                      &
       'wthetav flux at surface','K m/s',            &
       t_avg,d_sl,default_streams,'',routinename)

  TmpScm3d(:,:,:) = 0.0
  IF (l_3d_cca) THEN
    DO k=1, n_cca_lev
      DO j=1, rows
        DO i=1, row_length
          TmpScm3d(i,j,k) = cca(i,j,k)
        END DO
      END DO
    END DO
  ELSE
    DO j=1, rows
      DO i=1, row_length
        IF (ccb(i,j) /= 0) THEN
          DO k=ccb(i,j), cct(i,j)
            TmpScm3d(i,j,k) = cca(i,j,1)
          END DO
        END IF
      END DO
    END DO
  END IF

  ! Stash 5,212
  ! Convective Cloud Amount
  CALL scmoutput(TmpScm3d,'cca',                     &
       'Convective Cloud Amount','Fraction',         &
       t_avg,d_all,default_streams,'',routinename)

  ! Stash 5,213
  ! Convective Cloud Water
  CALL scmoutput(ccw,'ccw',                          &
       'Convective Cloud Water','kg/kg',             &
       t_avg,d_wet,default_streams,'',routinename)

  ! Condensed water path convective cloud
  CALL scmoutput(cclwp,'cclwp',                      &
       'Condensed Cloud Water Path','Kg/m2',         &
       t_avg,d_sl,default_streams,'',routinename)

  ! Unsure as to why cca_2d was limited to 0.5, this is here because
  ! this is the way the full UM outputs the diagnostic. A reason as
  ! to why a 0.5 cap should be found.
  DO j=1, rows
    DO i=1, row_length
      TmpScm2d(i,j) = MIN(0.5, cca_2d(i,j))
    END DO
  END DO

  CALL scmoutput(TmpScm2d,'cca_2d',                  &
       '2d Convective Cloud Amount','Fraction',      &
       t_avg,d_sl,default_streams,'',routinename)

  ! History diagnostics

  CALL scmoutput(deep_flag,'deep_flag',              &
       'History of deep convection',' ',             &
      t_avg,d_sl,default_streams,'',routinename)

  CALL scmoutput(past_precip,'past_precip',          &
       'History of convective precip','kg/m2/s',     &
       t_avg,d_sl,default_streams,'',routinename)

  CALL scmoutput(past_conv_ht,'past_conv_ht',        &
       'History of convective depth','m',            &
       t_avg,d_sl,default_streams,'',routinename)

  ! DD info
  CALL scmoutput(dt_dd,'dt_conv_dd',             &
       'dT/dt from conv DD','Ks-1',              &
       t_avg,d_all,default_streams,'',routinename)

  CALL scmoutput(dq_dd,'dq_conv_dd',             &
       'dq/dt from conv DD','kg/kg/s',           &
       t_avg,d_all,default_streams,'',routinename)

  CALL scmoutput(area_ud,'area_ud',                       &
       'Fractional area of updraught', 'dimensionless',   &
       t_avg,d_all,default_streams,'',routinename)

  CALL scmoutput(area_dd,'area_dd',                       &
       'Fractional area of downdraught', 'dimensionless', &
       t_avg,d_all,default_streams,'',routinename)

END IF ! scmdiag_conv

!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!----------------------------------------------------------------------
RETURN
END SUBROUTINE scm_diag_conv

END MODULE scm_diag_conv_mod
