! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_r_s_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_R_S_MOD'

CONTAINS
SUBROUTINE eg_r_s(                                                    &
           theta_star,q_star,qcl_star,qcf_star,qcf2_star,             &
           qrain_star, qgraup_star , theta_star_n ,m_star_n,          &
           mcl_star_n, mcf_star_n, mcf2_star_n, mrain_star_n,         &
           mgraup_star_n, l_physics, l_skeb2)

USE parkind1,                  ONLY: jpim, jprb       !DrHook
USE yomhook,                   ONLY: lhook, dr_hook   !DrHook
USE planet_constants_mod,      ONLY: epsln => repsilon
USE um_parcore,           ONLY: mype, nproc
USE eg_q_to_mix_mod
USE nlsizes_namelist_mod, ONLY: global_row_length, global_rows, model_levels
USE atm_fields_bounds_mod
USE field_types
USE gen_phys_inputs_mod,  ONLY: l_mr_physics
USE mphys_inputs_mod,     ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain
USE um_is_nan_mod,        ONLY: um_has_nan
USE check_nan_inf_mod, ONLY: fail_if_nan_inf

USE umPrintMgr

USE fields_rhs_mod,       ONLY: r_u_p2, r_v_p2, r_w_p2 , r_u_p2_n,    &
                                r_v_p2_n , r_w_p2_n, s_u,s_v,s_w,     &
                                s_thetav,s_m_v,s_m_cl,s_m_cf,         &
                                s_m_cf2,s_m_r,s_m_gr,                 &
                                r_u_skeb,r_v_skeb

IMPLICIT NONE
!
! Description:
!
!    Computes the fast increment S from the increment R and the
!     * state as returned from atmos_physics2

!    theta* is converted into theta_vd. The increment is
!    then computed by subtracting the value before the
!    call to atmos_physics2 saved in R_X_SP
!
!    At the end a simple implementation of boundary layer drag
!    is optionally added.
!
! Method:
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments


LOGICAL :: l_physics

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_R_S'


!increments from atmos_physics2:
REAL ::      theta_star(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::          m_star(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::        mcl_star(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::        mcf_star(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::       mcf2_star(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::      mrain_star(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::     mgraup_star(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::          q_star(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::        qcl_star(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::        qcf_star(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::       qcf2_star(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::      qrain_star(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::     qgraup_star(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

! previously stored increments
! (before atmos physics 2 aka F1SP predictor:
REAL ::        m_star_n(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::      mcl_star_n(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::      mcf_star_n(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::     mcf2_star_n(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::    mrain_star_n(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::   mgraup_star_n(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::    theta_star_n(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)


REAL    :: inv_epsilon
INTEGER :: i,j,k
INTEGER :: ierr

! local diagnostics
REAL :: max_s_theta,  min_s_theta,   av_s_theta,   has_nan_s_theta
REAL :: max_s_m_v,    min_s_m_v,     av_s_m_v,     has_nan_s_m_v
REAL :: max_s_m_cl,   min_s_m_cl,    av_s_m_cl,    has_nan_s_m_cl
REAL :: max_s_m_cf,   min_s_m_cf,    av_s_m_cf,    has_nan_s_m_cf
REAL :: max_s_m_cf2,  min_s_m_cf2,   av_s_m_cf2,   has_nan_s_m_cf2
REAL :: max_s_m_rain, min_s_m_rain,  av_s_m_rain,  has_nan_s_m_rain
REAL :: max_s_m_graup,min_s_m_graup, av_s_m_graup, has_nan_s_m_graup
REAL :: min_s_u,      max_s_u,       av_s_u,       has_nan_s_u
REAL :: min_s_v,      max_s_v,       av_s_v,       has_nan_s_v
REAL :: min_s_w,      max_s_w,       av_s_w,       has_nan_s_w

REAL :: fieldsize
LOGICAL :: l_skeb2

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

inv_epsilon = 1.0/epsln

IF (.NOT. l_physics) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF (l_mr_physics) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(tdims,m_star,mcl_star,mcf_star,q_star,qcl_star,   &
!$OMP& qcf_star,l_mcr_qcf2,l_mcr_qrain,l_mcr_qgraup,mcf2_star,  &
!$OMP& qcf2_star,mrain_star,qrain_star,mgraup_star,qgraup_star)
  DO k=tdims%k_start,tdims%k_end
    DO j=tdims%j_start,tdims%j_end
      DO i=tdims%i_start,tdims%i_end
        m_star  (i,j,k) = q_star  (i,j,k)
        mcl_star(i,j,k) = qcl_star(i,j,k)
        mcf_star(i,j,k) = qcf_star(i,j,k)
        IF (l_mcr_qcf2  ) mcf2_star  (i,j,k) = qcf2_star  (i,j,k)
        IF (l_mcr_qrain ) mrain_star  (i,j,k) = qrain_star  (i,j,k)
        IF (l_mcr_qgraup) mgraup_star(i,j,k) = qgraup_star(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

ELSE
  CALL eg_q_to_mix                                                    &
                  (tdims,tdims,                                       &
                   q_star, qcl_star, qcf_star,                        &
                   qcf2_star, qrain_star, qgraup_star,                &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,             &
                   m_star, mcl_star, mcf_star,                        &
                   mcf2_star, mrain_star, mgraup_star, swap_in=.FALSE.)
END IF

!$OMP PARALLEL  PRIVATE(i,j,k) SHARED(udims, s_u,r_u_p2,              &
!$OMP&            r_u_p2_n,vdims,s_v,r_v_p2,r_v_p2_n,pdims,s_w,       &
!$OMP&            r_w_p2,r_w_p2_n,r_u_skeb,r_v_skeb,tdims_s,          &
!$OMP&            theta_star,m_star,inv_epsilon,s_m_v,s_m_cl,mcl_star,&
!$OMP&            mcl_star_n,s_m_cf,mcf_star,mcf_star_n,s_thetav,     &
!$OMP&            tdims,mrain_star,mrain_star_n,mgraup_star,          &
!$OMP&            l_mcr_qcf2,l_mcr_qrain, l_mcr_qgraup,s_m_r,s_m_gr,  &
!$OMP&            s_m_cf2,l_skeb2,theta_star_n,m_star_n,              &
!$OMP&            mgraup_star_n,mcf2_star,mcf2_star_n) DEFAULT(NONE)

IF (l_skeb2) THEN
!$OMP DO SCHEDULE (STATIC)
  DO k=udims%k_start,udims%k_end
    DO j=udims%j_start,udims%j_end
      DO i=udims%i_start,udims%i_end
        s_u(i,j,k) = r_u_p2(i,j,k) - r_u_p2_n(i,j,k) + r_u_skeb(i,j,k)
      END DO
    END DO

    DO j=vdims%j_start,vdims%j_end
      DO i=vdims%i_start,vdims%i_end
        s_v(i,j,k) = r_v_p2(i,j,k) - r_v_p2_n(i,j,k) + r_v_skeb(i,j,k)
      END DO
    END DO

    !   zeroth level for s_w is done below
    DO j=pdims%j_start,pdims%j_end
      DO i=pdims%i_start,pdims%i_end
        s_w(i,j,k) = r_w_p2(i,j,k) - r_w_p2_n(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
ELSE
!$OMP DO SCHEDULE (STATIC)
  DO k=udims%k_start,udims%k_end
    DO j=udims%j_start,udims%j_end
      DO i=udims%i_start,udims%i_end
        s_u(i,j,k) = r_u_p2(i,j,k) - r_u_p2_n(i,j,k)
      END DO
    END DO

    DO j=vdims%j_start,vdims%j_end
      DO i=vdims%i_start,vdims%i_end
        s_v(i,j,k) = r_v_p2(i,j,k) - r_v_p2_n(i,j,k)
      END DO
    END DO

    !   zeroth level for s_w is done below
    DO j=pdims%j_start,pdims%j_end
      DO i=pdims%i_start,pdims%i_end
        s_w(i,j,k) = r_w_p2(i,j,k) - r_w_p2_n(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

!$OMP DO SCHEDULE (STATIC)
DO k = tdims%k_start, tdims%k_end
   DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        s_thetav(i,j,k) = theta_star(i,j,k)                       &
                             *(1.0 + m_star(i,j,k)*inv_epsilon)   &
                         -theta_star_n(i,j,k)                     &
                             *(1.0 + m_star_n(i,j,k)*inv_epsilon)

        s_m_v(i,j,k)  = m_star(i,j,k) - m_star_n(i,j,k)

        s_m_cl(i,j,k) = mcl_star(i,j,k) - mcl_star_n(i,j,k)

        s_m_cf(i,j,k) = mcf_star(i,j,k) - mcf_star_n(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

IF (l_mcr_qcf2 ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims%k_start,tdims%k_end
    DO j=pdims%j_start,pdims%j_end
      DO i=pdims%i_start,pdims%i_end
        s_m_cf2(i,j,k) = mcf2_star(i,j,k) - mcf2_star_n(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (l_mcr_qrain ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims%k_start,tdims%k_end
    DO j=pdims%j_start,pdims%j_end
      DO i=pdims%i_start,pdims%i_end
        s_m_r(i,j,k) = mrain_star(i,j,k) - mrain_star_n(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (l_mcr_qgraup) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k=tdims%k_start,tdims%k_end
    DO j=pdims%j_start,pdims%j_end
      DO i=pdims%i_start,pdims%i_end
        s_m_gr(i,j,k) = mgraup_star(i,j,k) - mgraup_star_n(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF


!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    s_w(i,j,0) = r_w_p2(i,j,0) - r_w_p2_n(i,j,0)
  END DO
END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


IF (printstatus >=  prstatus_diag) THEN

  fieldsize       = global_row_length * global_rows * model_levels
  max_s_u         = MAXVAL(s_u(udims%i_start:udims%i_end,             &
                               udims%j_start:udims%j_end,             &
                               udims%k_start:udims%k_end))
  min_s_u         = MINVAL(s_u(udims%i_start:udims%i_end,             &
                               udims%j_start:udims%j_end,             &
                               udims%k_start:udims%k_end))
  av_s_u          =    SUM(s_u(udims%i_start:udims%i_end,             &
                               udims%j_start:udims%j_end,             &
                               udims%k_start:udims%k_end))/fieldsize
  has_nan_s_u = 0.0
  IF(um_has_nan(s_u(udims%i_start:udims%i_end,                        &
                    udims%j_start:udims%j_end,                        &
                    udims%k_start:udims%k_end))) has_nan_s_u = 1.0

  fieldsize       = global_row_length * (global_rows+1) * model_levels
  max_s_v         = MAXVAL(s_v(vdims%i_start:vdims%i_end,             &
                               vdims%j_start:vdims%j_end,             &
                               vdims%k_start:vdims%k_end))
  min_s_v         = MINVAL(s_v(vdims%i_start:vdims%i_end,             &
                               vdims%j_start:vdims%j_end,             &
                               vdims%k_start:vdims%k_end))
  av_s_v          =    SUM(s_v(vdims%i_start:vdims%i_end,             &
                               vdims%j_start:vdims%j_end,             &
                               vdims%k_start:vdims%k_end))/fieldsize
  has_nan_s_v = 0.0
  IF(um_has_nan(s_v(vdims%i_start:vdims%i_end,                        &
                    vdims%j_start:vdims%j_end,                        &
                    vdims%k_start:vdims%k_end))) has_nan_s_v = 1.0

  fieldsize       = global_row_length * global_rows *(model_levels+1)
  max_s_w         = MAXVAL(s_w(wdims%i_start:wdims%i_end,             &
                               wdims%j_start:wdims%j_end,             &
                               wdims%k_start:wdims%k_end))
  min_s_w         = MINVAL(s_w(wdims%i_start:wdims%i_end,             &
                               wdims%j_start:wdims%j_end,             &
                               wdims%k_start:wdims%k_end))
  av_s_w          =    SUM(s_w(wdims%i_start:wdims%i_end,             &
                               wdims%j_start:wdims%j_end,             &
                               wdims%k_start:wdims%k_end))/fieldsize

  has_nan_s_w = 0.0
  IF(um_has_nan(s_w(wdims%i_start:wdims%i_end,                        &
                    wdims%j_start:wdims%j_end,                        &
                    wdims%k_start:wdims%k_end))) has_nan_s_w = 1.0

  fieldsize       = global_row_length * global_rows *(model_levels+1)
  max_s_theta     = MAXVAL(s_thetav(tdims%i_start:tdims%i_end,        &
                                    tdims%j_start:tdims%j_end,        &
                                    tdims%k_start:tdims%k_end))
  min_s_theta     = MINVAL(s_thetav(tdims%i_start:tdims%i_end,        &
                                    tdims%j_start:tdims%j_end,        &
                                    tdims%k_start:tdims%k_end))
  av_s_theta      =    SUM(s_thetav(tdims%i_start:tdims%i_end,        &
                                    tdims%j_start:tdims%j_end,        &
                                    tdims%k_start:tdims%k_end))/fieldsize
  has_nan_s_theta = 0.0
  IF(um_has_nan(s_thetav(tdims%i_start:tdims%i_end,                   &
                         tdims%j_start:tdims%j_end,                   &
                         tdims%k_start:tdims%k_end))) has_nan_s_theta = 1.0


  max_s_m_v       = MAXVAL(s_m_v(tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
  min_s_m_v       = MINVAL(s_m_v(tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
  av_s_m_v        =    SUM(s_m_v(tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))/fieldsize
  has_nan_s_m_v = 0.0
  IF(um_has_nan(s_m_v(tdims%i_start:tdims%i_end,                      &
                      tdims%j_start:tdims%j_end,                      &
                      tdims%k_start:tdims%k_end))) has_nan_s_m_v = 1.0

  max_s_m_cl      = MAXVAL(s_m_cl(tdims%i_start:tdims%i_end,          &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
  min_s_m_cl      = MINVAL(s_m_cl(tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end))
  av_s_m_cl       =    SUM(s_m_cl(tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end))/fieldsize
  has_nan_s_m_cl = 0.0
  IF(um_has_nan(s_m_cl(tdims%i_start:tdims%i_end,                     &
                       tdims%j_start:tdims%j_end,                     &
                       tdims%k_start:tdims%k_end))) has_nan_s_m_cl = 1.0


  max_s_m_cf      = MAXVAL(s_m_cf(tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end))
  min_s_m_cf      = MINVAL(s_m_cf(tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end))
  av_s_m_cf       =    SUM(s_m_cf(tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end))/fieldsize
  has_nan_s_m_cf = 0.0
  IF(um_has_nan(s_m_cf(tdims%i_start:tdims%i_end,                     &
                       tdims%j_start:tdims%j_end,                     &
                       tdims%k_start:tdims%k_end))) has_nan_s_m_cf = 1.0


  IF (l_mcr_qcf2) THEN
    max_s_m_cf2   = MAXVAL(s_m_cf2(tdims%i_start:tdims%i_end,         &
                                   tdims%j_start:tdims%j_end,         &
                                   tdims%k_start:tdims%k_end))
    min_s_m_cf2   = MINVAL(s_m_cf2(tdims%i_start:tdims%i_end,         &
                                   tdims%j_start:tdims%j_end,         &
                                   tdims%k_start:tdims%k_end))
    av_s_m_cf2    =    SUM(s_m_cf2(tdims%i_start:tdims%i_end,         &
                                   tdims%j_start:tdims%j_end,         &
                                   tdims%k_start:tdims%k_end))/fieldsize
    has_nan_s_m_cf2 = 0.0
    IF(um_has_nan(s_m_cf2(tdims%i_start:tdims%i_end,                  &
                          tdims%j_start:tdims%j_end,                  &
                          tdims%k_start:tdims%k_end))) has_nan_s_m_cf2 = 1.0

  END IF

  IF (l_mcr_qrain) THEN
    max_s_m_rain  = MAXVAL(s_m_r(tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
    min_s_m_rain  = MINVAL(s_m_r(tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
    av_s_m_rain   =    SUM(s_m_r(tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))/fieldsize
    has_nan_s_m_rain = 0.0
    IF(um_has_nan(s_m_r(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end))) has_nan_s_m_rain = 1.0

  END IF

  IF (l_mcr_qgraup) THEN
    max_s_m_graup = MAXVAL(s_m_gr(tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end))
    min_s_m_graup = MINVAL(s_m_gr(tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end))
    av_s_m_graup  =    SUM(s_m_gr(tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end))/fieldsize
    has_nan_s_m_graup = 0.0
    IF(um_has_nan(s_m_gr(tdims%i_start:tdims%i_end,                   &
                         tdims%j_start:tdims%j_end,                   &
                         tdims%k_start:tdims%k_end))) has_nan_s_m_graup = 1.0

  END IF

  CALL gc_rmax(1,nproc,ierr,max_s_u)
  CALL gc_rmin(1,nproc,ierr,min_s_u)
  CALL gc_rmax(1,nproc,ierr,has_nan_s_u)
  CALL gc_rsum(1,nproc,ierr, av_s_u)
  CALL gc_rmax(1,nproc,ierr,max_s_v)
  CALL gc_rmin(1,nproc,ierr,min_s_v)
  CALL gc_rmax(1,nproc,ierr,has_nan_s_v)
  CALL gc_rsum(1,nproc,ierr, av_s_v)
  CALL gc_rmax(1,nproc,ierr,max_s_w)
  CALL gc_rmin(1,nproc,ierr,min_s_w)
  CALL gc_rmax(1,nproc,ierr,has_nan_s_w)
  CALL gc_rsum(1,nproc,ierr, av_s_w)
  CALL gc_rmax(1,nproc,ierr,max_s_theta)
  CALL gc_rmin(1,nproc,ierr,min_s_theta)
  CALL gc_rmax(1,nproc,ierr,has_nan_s_theta)
  CALL gc_rsum(1,nproc,ierr, av_s_theta)
  CALL gc_rmax(1,nproc,ierr,max_s_m_v)
  CALL gc_rmin(1,nproc,ierr,min_s_m_v)
  CALL gc_rmax(1,nproc,ierr,has_nan_s_m_v)
  CALL gc_rsum(1,nproc,ierr, av_s_m_v)
  CALL gc_rmax(1,nproc,ierr,max_s_m_cl)
  CALL gc_rmin(1,nproc,ierr,min_s_m_cl)
  CALL gc_rmax(1,nproc,ierr,has_nan_s_m_cl)
  CALL gc_rsum(1,nproc,ierr, av_s_m_cl)
  CALL gc_rmax(1,nproc,ierr,max_s_m_cf)
  CALL gc_rmin(1,nproc,ierr,min_s_m_cf)
  CALL gc_rmax(1,nproc,ierr,has_nan_s_m_cf)
  CALL gc_rsum(1,nproc,ierr, av_s_m_cf)
  IF (l_mcr_qcf2) THEN
    CALL gc_rmax(1,nproc,ierr,max_s_m_cf2)
    CALL gc_rmin(1,nproc,ierr,min_s_m_cf2)
    CALL gc_rmax(1,nproc,ierr,has_nan_s_m_cf2)
    CALL gc_rsum(1,nproc,ierr, av_s_m_cf2)
  END IF
  IF (l_mcr_qrain) THEN
    CALL gc_rmax(1,nproc,ierr,max_s_m_rain)
    CALL gc_rmin(1,nproc,ierr,min_s_m_rain)
    CALL gc_rmax(1,nproc,ierr,has_nan_s_m_rain)
    CALL gc_rsum(1,nproc,ierr, av_s_m_rain)
  END IF
  IF (l_mcr_qgraup) THEN
    CALL gc_rmax(1,nproc,ierr,max_s_m_graup)
    CALL gc_rmin(1,nproc,ierr,min_s_m_graup)
    CALL gc_rmax(1,nproc,ierr,has_nan_s_m_graup)
    CALL gc_rsum(1,nproc,ierr, av_s_m_graup)
  END IF

  IF (mype == 0) THEN
    CALL umPrint('********************************************'//&
        '************************************************',src='eg_R_S')
    CALL umPrint('Fast physics sources for ENDGame from atmos_'//&
        'physics2:',src='eg_R_S')
    CALL umPrint('             min                      max   '//&
        '                   average (non-bit reproduc'//&
        'ing  (1=has NaN 0= no NaN)',src='eg_R_S')
    WRITE(umMessage,'(A,4E32.16)') 's_u      :',min_s_u    ,max_s_u            &
                                      , av_s_u    ,has_nan_s_u
    CALL umPrint(umMessage,src='eg_R_S')
    WRITE(umMessage,FMT='(A,4E32.16)') 's_v      :',min_s_v    ,max_s_v        &
                                      , av_s_v    ,has_nan_s_v
    CALL umPrint(umMessage,src='eg_R_S')
    WRITE(umMessage,FMT='(A,4E32.16)') 's_w      :',min_s_w    ,max_s_w        &
                                      , av_s_w    ,has_nan_s_w
    CALL umPrint(umMessage,src='eg_R_S')
    WRITE(umMessage,FMT='(A,4E32.16)') 's_thetav :',min_s_theta,max_s_theta    &
                                      , av_s_theta    ,has_nan_s_theta
    CALL umPrint(umMessage,src='eg_R_S')
    WRITE(umMessage,FMT='(A,4E32.16)') 's_m_v    :',min_s_m_v  ,max_s_m_v      &
                                      , av_s_m_v    ,has_nan_s_m_v
    CALL umPrint(umMessage,src='eg_R_S')
    WRITE(umMessage,FMT='(A,4E32.16)') 's_m_cl   :',min_s_m_cl ,max_s_m_cl     &
                                      , av_s_m_cl    ,has_nan_s_m_cl
    CALL umPrint(umMessage,src='eg_R_S')
    WRITE(umMessage,FMT='(A,4E32.16)') 's_m_cf   :',min_s_m_cf ,max_s_m_cf     &
                                      , av_s_m_cf    ,has_nan_s_m_cf
    CALL umPrint(umMessage,src='eg_R_S')
    IF (l_mcr_qcf2) THEN  
      WRITE(umMessage,'(A,4E32.16)') 's_m_cf2  :',min_s_m_cf2    &
                      ,max_s_m_cf2    , av_s_m_cf2    ,has_nan_s_m_cf2
      CALL umPrint(umMessage,src='eg_R_S')
    END IF
    IF (l_mcr_qgraup) THEN
      WRITE(umMessage,'(A,4E32.16)') 's_m_graup:',min_s_m_graup  &
          ,max_s_m_graup  , av_s_m_graup    ,has_nan_s_m_graup
      CALL umPrint(umMessage,src='eg_R_S')
    END IF
    IF (l_mcr_qrain) THEN
      WRITE(umMessage,'(A,4E32.16)') 's_m_rain :',min_s_m_rain   &
          ,max_s_m_rain   , av_s_m_rain    ,has_nan_s_m_rain
      CALL umPrint(umMessage,src='eg_R_S')
    END IF

    CALL umPrint('***********************************************'//&
        '*********************************************',src='eg_R_S')
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_r_s
END MODULE eg_r_s_mod
