! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Set flags for use in convection code

MODULE cv_stash_flg_mod

! Description:
!   Module containing stash flags used in top level convection routine.
!
! Method:
!   Declares all flags and contains a subroutine used to set the flags.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.
!
! Declarations:

USE science_fixes_mod,  ONLY: l_fix_conv_diags_var

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
SAVE

!=====================================================================
! Logical flags used to control output/setting of various convection
! diagnostics
!=====================================================================

LOGICAL ::           &
  l_apply_diag         ! TRUE if sub-step for calculating diagnostics

LOGICAL ::           &
  flg_up_flx         & ! stash flag for updraught mass flux
 ,flg_up_flx_half    & ! stash flag for updraught mass flux
 ,flg_dwn_flx        & ! stash flag for downdraught mass flux
 ,flg_entr_up        & ! stash flag for updraught entrainment
 ,flg_entr_dwn       & ! stash flag for downdraught entrainmnent
 ,flg_detr_up        & ! stash flag for updraught detrainment
 ,flg_detr_dwn       & ! stash flag for downdraught detrainmnent
 ,flg_conv_rain_3d   & ! stash flag for 3d conv rainfall
 ,flg_conv_snow_3d     ! stash flag for 3d conv snowfall

! CMT diagnostics

LOGICAL ::           &
  flg_uw_dp          & ! stash flag for deep x-comp stress
 ,flg_vw_dp          & ! stash flag for deep y-comp stress
 ,flg_uw_shall       & ! stash flag for shallow x stress
 ,flg_vw_shall       & ! stash flag for shallow y stress
 ,flg_uw_mid         & ! stash flag for mid-level x stress
 ,flg_vw_mid           ! stash flag for mid-level y stress
! Increment diagnostics:

LOGICAL ::           &
  l_qcl_incr_cinh    & ! liquid cloud condensate qCL (inhom)
 ,l_qcf_incr_cinh    & ! frozen cloud condensate qCF (inhom)
 ,l_cfl_incr_cinh    & ! liquid cloud amount cf_liquid (inhom)
 ,l_cff_incr_cinh    & ! frozen cloud amount cf_frozen (inhom)
 ,l_bcf_incr_cinh    & ! total (bulk) cloud amount bulk_cf (inhom)

 ,l_theta_incr_conv  & ! potential temperature increment
 ,l_T_incr_conv      & ! temperature
 ,l_q_incr_conv      & ! humidity
 ,l_qcl_incr_conv    & ! liquid cloud condensate qCL
 ,l_qcf_incr_conv    & ! frozen cloud condensate qCF
 ,l_cfl_incr_conv    & ! liquid cloud amount cf_liquid
 ,l_cff_incr_conv    & ! frozen cloud amount cf_frozen
 ,l_bcf_incr_conv    & ! total (bulk) cloud amount bulk_cf
 ,l_T_conv_only      & ! temperature from convection only
 ,l_q_conv_only        ! humidity from convection only

! DD increments

LOGICAL ::           &
  flg_dt_dd          & ! stash flag for dT/dt from DD
 ,flg_dq_dd          & ! stash flag for dq/dt from DD
 ,flg_du_dd          & ! stash flag for du/dt from DD
 ,flg_dv_dd            ! stash flag for dv/dt from DD

! Fractional areas
LOGICAL ::           &
  flg_area_ud        & ! stash flag for updraught fractional area
 ,flg_area_dd          ! stash flag for downdraught fractional area


! 5A turbulence based schemes

LOGICAL ::           &
  flg_wqt_flux       & ! stash flag for w'qt'
 ,flg_wql_flux       & ! stash flag for w'ql'
 ,flg_wthetal_flux   & ! stash flag for w'thetal'
 ,flg_wthetav_flux   & ! stash flag for w'thetav'

 ,flg_deep_tops      & ! stash flag for deep tops

 ,flg_mf_deep        & ! stash flag for deep mass flux
 ,flg_mf_congest     & ! stash flag for congestus mass flux
 ,flg_mf_shall       & ! stash flag for shallow mass flux
 ,flg_mf_midlev      & ! stash flag for mid_level mass flux

 ,flg_dt_deep        & ! stash flag for deep dT
 ,flg_dt_congest     & ! stash flag for congestus dT
 ,flg_dt_shall       & ! stash flag for shallow dT
 ,flg_dt_midlev      & ! stash flag for mid_level dT

 ,flg_dq_deep        & ! stash flag for deep dq
 ,flg_dq_congest     & ! stash flag for congestus dq
 ,flg_dq_shall       & ! stash flag for shallow dq
 ,flg_dq_midlev      & ! stash flag for mid_level dq

 ,flg_du_deep        & ! stash flag for deep du
 ,flg_du_congest     & ! stash flag for congestus du
 ,flg_du_shall       & ! stash flag for shallow du
 ,flg_du_midlev      & ! stash flag for mid_level du

 ,flg_dv_deep        & ! stash flag for deep dv
 ,flg_dv_congest     & ! stash flag for congestus dv
 ,flg_dv_shall       & ! stash flag for shallow dv
 ,flg_dv_midlev        ! stash flag for mid_level dv

! Flag for Simpson & Wiggert w-equation diagnostic
LOGICAL ::           &
  flg_w_eqn            ! stash flag for w_eqn


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CV_STASH_FLG_MOD'

CONTAINS

!=========================================================================
! Full model version - sets flags according to user requirements
! i.e. stash settings held in array sf
!=========================================================================
SUBROUTINE set_convection_output_flags()

USE jules_surface_mod,  ONLY: ISrfExCnvGust, IP_SrfExWithCnv
USE model_domain_mod,   ONLY: model_type, mt_single_column
USE stash_array_mod,    ONLY: sf
USE cosp_input_mod,     ONLY: l_cosp

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_CONVECTION_OUTPUT_FLAGS'

!---------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (model_type /= mt_single_column) THEN

  flg_conv_rain_3d = (sf(227,5) .AND. l_apply_diag) .OR. l_cosp
  flg_conv_snow_3d = (sf(228,5) .AND. l_apply_diag) .OR. l_cosp

  flg_up_flx = sf(250,5) .AND. l_apply_diag

  flg_up_flx_half = ( sf(249,5) .OR. sf(246,5) ) .AND. l_apply_diag

  flg_dwn_flx     = ( sf(251,5) .OR. (ISrfExCnvGust == IP_SrfExWithCnv) )      &
                    .AND. l_apply_diag

  flg_entr_up  = sf(252,5) .AND. l_apply_diag
  flg_detr_up  = sf(253,5) .AND. l_apply_diag
  flg_entr_dwn = sf(254,5) .AND. l_apply_diag
  flg_detr_dwn = sf(255,5) .AND. l_apply_diag

  flg_uw_dp    = sf(258,5) .AND. l_apply_diag
  flg_vw_dp    = sf(259,5) .AND. l_apply_diag
  flg_uw_shall = sf(260,5) .AND. l_apply_diag
  flg_vw_shall = sf(261,5) .AND. l_apply_diag
  flg_uw_mid   = sf(263,5) .AND. l_apply_diag
  flg_vw_mid   = sf(264,5) .AND. l_apply_diag

  ! PC2 increments

  l_qcl_incr_cinh = ( sf(163,5) .OR. sf(140,5) .OR. sf(141,5) .OR.             &
                      sf(150,5) .OR. sf(151,5) ).AND. l_apply_diag
  l_qcf_incr_cinh = (sf(164,5) .OR. (sf(142,5) .OR. sf(143,5)))                &
                    .AND. l_apply_diag
  l_bcf_incr_cinh = sf(172,5) .AND. l_apply_diag
  l_cfl_incr_cinh = ( sf(173,5) .OR. sf(146,5) .OR. sf(147,5) .OR.             &
                      sf(156,5) .OR. sf(157,5) .OR. sf(158,5) )                &
                    .AND. l_apply_diag
  l_cff_incr_cinh = ( sf(174,5) .OR. sf(148,5) .OR. sf(149,5) )                &
                    .AND. l_apply_diag
  l_bcf_incr_conv = sf(192,5) .AND. l_apply_diag
  l_cfl_incr_conv = (sf(193,5) .OR. (sf(156,5) .OR. sf(157,5) .OR. sf(158,5))) &
                    .AND. l_apply_diag
  l_cff_incr_conv = sf(194,5) .AND. l_apply_diag

  ! Convection increments

  l_theta_incr_conv = .FALSE.  ! Theta increment only output by the SCM
  l_T_incr_conv   = ( sf(181,5) .OR. sf(187,5) .OR. sf(161,5) )                &
                    .AND. l_apply_diag
  l_q_incr_conv   = ( sf(182,5) .OR. sf(188,5) .OR. sf(162,5) )                &
                    .AND. l_apply_diag
  l_qcl_incr_conv = ( sf(183,5) .OR. sf(150,5) .OR. sf(151,5) .OR. sf(152,5) ) &
                    .AND. l_apply_diag
  l_qcf_incr_conv = sf(184,5) .AND. l_apply_diag

  l_T_conv_only   = ((l_fix_conv_diags_var .AND. sf(187,5))                    &
                    .OR. sf(161,5)) .AND. l_apply_diag
  l_q_conv_only   = ((l_fix_conv_diags_var .AND. sf(188,5))                    &
                    .OR. sf(162,5)) .AND. l_apply_diag


  ! 5A convection code

  flg_wqt_flux     = (sf(290,5) .OR. sf(304,5) .OR. sf(306,5) .OR.             &
                      sf(429,5) .OR. sf(430,5) .OR. sf(431,5) .OR. sf(432,5))  &
                     .AND. l_apply_diag

  flg_wql_flux     = sf(291,5) .AND. l_apply_diag

  flg_wthetal_flux = (sf(292,5) .OR. sf(305,5) .OR. sf(307,5) .OR.             &
                      sf(425,5) .OR. sf(426,5) .OR. sf(427,5) .OR. sf(428,5))  &
                     .AND. l_apply_diag

  flg_wthetav_flux = sf(293,5) .AND. l_apply_diag
  flg_deep_tops    = sf(319,5) .AND. l_apply_diag
  flg_mf_deep      = sf(320,5) .AND. l_apply_diag
  flg_mf_congest   = sf(321,5) .AND. l_apply_diag

  flg_mf_shall     = (sf(322,5) .OR. sf(417,5) .OR. sf(418,5) .OR.             &
                      sf(419,5) .OR. sf(420,5)) .AND. l_apply_diag
  flg_mf_midlev    = sf(323,5) .AND. l_apply_diag
  flg_dt_deep      = sf(324,5) .AND. l_apply_diag
  flg_dt_congest   = sf(325,5) .AND. l_apply_diag
  flg_dt_shall     = (sf(326,5) .OR. sf(409,5) .OR. sf(410,5) .OR.             &
                      sf(411,5) .OR. sf(412,5)) .AND. l_apply_diag
  flg_dt_midlev    = sf(327,5) .AND. l_apply_diag
  flg_dq_deep      = sf(328,5) .AND. l_apply_diag
  flg_dq_congest   = sf(329,5) .AND. l_apply_diag
  flg_dq_shall     = (sf(330,5) .OR. sf(413,5) .OR. sf(414,5) .OR.             &
                      sf(415,5) .OR. sf(416,5)) .AND. l_apply_diag
  flg_dq_midlev    = sf(331,5) .AND. l_apply_diag
  flg_du_deep      = sf(332,5) .AND. l_apply_diag
  flg_du_congest   = sf(333,5) .AND. l_apply_diag
  flg_du_shall     = sf(334,5) .AND. l_apply_diag
  flg_du_midlev    = sf(335,5) .AND. l_apply_diag
  flg_dv_deep      = sf(336,5) .AND. l_apply_diag
  flg_dv_congest   = sf(337,5) .AND. l_apply_diag
  flg_dv_shall     = sf(338,5) .AND. l_apply_diag
  flg_dv_midlev    = sf(339,5) .AND. l_apply_diag

  ! This is also needed for the calculation of an updraught area (5,229)
  flg_w_eqn        = (sf(196,5) .OR. sf(197,5) .OR. sf(229,5))                 &
                     .AND. l_apply_diag

  flg_dt_dd        = sf(198,5) .AND. l_apply_diag
  flg_dq_dd        = sf(199,5) .AND. l_apply_diag
  flg_du_dd        = sf(175,5) .AND. l_apply_diag
  flg_dv_dd        = sf(176,5) .AND. l_apply_diag

  flg_area_ud      = sf(229,5) .AND. l_apply_diag
  flg_area_dd      = sf(230,5) .AND. l_apply_diag

ELSE

  ! Setting of flags defined in module

  l_apply_diag     = .TRUE.

  flg_conv_rain_3d = .TRUE.
  flg_conv_snow_3d = .TRUE.

  flg_up_flx_half  = .TRUE.
  flg_up_flx       = .TRUE.
  flg_dwn_flx      = .TRUE.
  flg_entr_up      = .TRUE.
  flg_detr_up      = .TRUE.
  flg_entr_dwn     = .TRUE.
  flg_detr_dwn     = .TRUE.

  flg_uw_dp        = .TRUE.
  flg_vw_dp        = .TRUE.
  flg_uw_shall     = .TRUE.
  flg_vw_shall     = .TRUE.
  flg_uw_mid       = .TRUE.
  flg_vw_mid       = .TRUE.

  l_qcl_incr_cinh  = .TRUE.
  l_qcf_incr_cinh  = .TRUE.
  l_bcf_incr_cinh  = .TRUE.
  l_cfl_incr_cinh  = .TRUE.
  l_cff_incr_cinh  = .TRUE.
  l_bcf_incr_conv  = .TRUE.
  l_cfl_incr_conv  = .TRUE.
  l_cff_incr_conv  = .TRUE.

  l_theta_incr_conv= .TRUE.
  l_T_incr_conv    = .TRUE.
  l_q_incr_conv    = .TRUE.
  l_qcl_incr_conv  = .TRUE.
  l_qcf_incr_conv  = .TRUE.

  l_t_conv_only = .FALSE.         ! Not output by the SCM for some reason.
  l_q_conv_only = .FALSE.

  flg_dt_dd = .TRUE.
  flg_dq_dd = .TRUE.

  ! 5A & 6A convection schemes
  flg_du_dd    = .TRUE.
  flg_dv_dd    = .TRUE.
  flg_area_ud  = .TRUE.
  flg_area_dd  = .TRUE.
  flg_wqt_flux = .TRUE.
  flg_wql_flux = .TRUE.

  flg_wthetal_flux = .TRUE.
  flg_wthetav_flux = .TRUE.
  flg_deep_tops    = .TRUE.

  flg_mf_deep    = .TRUE.
  flg_mf_congest = .TRUE.
  flg_mf_shall   = .TRUE.
  flg_mf_midlev  = .TRUE.
  flg_dt_deep    = .TRUE.
  flg_dt_congest = .TRUE.
  flg_dt_shall   = .TRUE.
  flg_dt_midlev  = .TRUE.
  flg_dq_deep    = .TRUE.
  flg_dq_congest = .TRUE.
  flg_dq_shall   = .TRUE.
  flg_dq_midlev  = .TRUE.
  flg_du_deep    = .TRUE.
  flg_du_congest = .TRUE.
  flg_du_shall   = .TRUE.
  flg_du_midlev  = .TRUE.
  flg_dv_deep    = .TRUE.
  flg_dv_congest = .TRUE.
  flg_dv_shall   = .TRUE.
  flg_dv_midlev  = .TRUE.
  flg_w_eqn      = .TRUE.

END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE set_convection_output_flags

END MODULE cv_stash_flg_mod
