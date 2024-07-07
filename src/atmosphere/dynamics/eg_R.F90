! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_r_mod
USE umPrintMgr
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_R_MOD'

CONTAINS
SUBROUTINE eg_r(theta_star, m_star, mcl_star,mcf_star,mgraup_star,    &
                mrain_star, mcf2_star ,m_v, theta )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE planet_constants_mod, ONLY: epsln=>repsilon
USE um_parcore,       ONLY: mype, nproc
USE nlsizes_namelist_mod, ONLY: global_row_length, global_rows, model_levels
USE atm_fields_bounds_mod
USE field_types
USE um_is_nan_mod,    ONLY: um_has_nan
USE check_nan_inf_mod, ONLY: fail_if_nan_inf

USE fields_rhs_mod,  ONLY:  r_m_v ,r_m_cl  ,r_m_cf  ,r_m_gr,       &
                            r_m_r  ,r_m_cf2,r_theta

USE mphys_inputs_mod,      ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain


IMPLICIT NONE
!
! Description: compute the slow physics source terms for the dynamics
!
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
REAL ::      theta_star(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::           theta(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::          m_star(tdims%i_start:tdims%i_end,                    &
                        tdims%j_start:tdims%j_end,                    &
                        tdims%k_start:tdims%k_end)

REAL ::             m_v(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

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

REAL :: inv_epsln

INTEGER :: ierr
INTEGER :: i  ! Looper
INTEGER :: j  ! Looper
INTEGER :: k  ! Looper

! local diagnostics
REAL :: max_r_theta  ,min_r_theta,   av_r_theta,   has_nan_r_theta
REAL :: max_r_m_v    ,min_r_m_v,     av_r_m_v,     has_nan_r_m_v
REAL :: max_r_m_cl   ,min_r_m_cl,    av_r_m_cl,    has_nan_r_m_cl
REAL :: max_r_m_cf   ,min_r_m_cf,    av_r_m_cf,    has_nan_r_m_cf
REAL :: max_r_m_cf2  ,min_r_m_cf2,   av_r_m_cf2,   has_nan_r_m_cf2
REAL :: max_r_m_rain ,min_r_m_rain,  av_r_m_rain,  has_nan_r_m_rain
REAL :: max_r_m_graup,min_r_m_graup, av_r_m_graup, has_nan_r_m_graup

REAL :: fieldsize

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_R'

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

inv_epsln = 1.0/epsln

!
! atmos_physics1 returns increments.
! So no conversion is required (apart from dry -> virtual dry)
! Note: calculations done in the halos to avoid halo exchange of result.
!

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)        &
!$OMP& SHARED(tdims, r_theta, theta_star,                               &
!$OMP&        m_star, m_v, inv_epsln, theta, r_m_v, r_m_cl, mcl_star,   &
!$OMP&        r_m_cf, mcf_star, l_mcr_qcf2, r_m_cf2, mcf2_star,         &
!$OMP&        l_mcr_qrain, r_m_r, mrain_star, l_mcr_qgraup, r_m_gr,     &
!$OMP&        mgraup_star)
DO k = tdims%k_start, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      r_theta(i,j,k) = theta_star(i,j,k) * (1.0 + (m_star(i,j,k) +     &
                       m_v(i,j,k) ) *inv_epsln) + theta(i,j,k) *      &
                       m_star(i,j,k)*inv_epsln

      r_m_v (i,j,k) = m_star  (i,j,k)
      r_m_cl(i,j,k) = mcl_star(i,j,k)
      r_m_cf(i,j,k) = mcf_star(i,j,k)

      IF (l_mcr_qcf2)   r_m_cf2(i,j,k) = mcf2_star  (i,j,k)
      IF (l_mcr_qrain)  r_m_r  (i,j,k) = mrain_star (i,j,k)
      IF (l_mcr_qgraup) r_m_gr (i,j,k) = mgraup_star(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO


! Everything below here is diagnostic...
IF (printstatus >=  prstatus_diag) THEN

  fieldsize       = global_row_length * global_rows *(model_levels+1)

  max_r_theta     = MAXVAL(r_theta(tdims%i_start:tdims%i_end,        &
                                   tdims%j_start:tdims%j_end,        &
                                   tdims%k_start:tdims%k_end))
  min_r_theta     = MINVAL(r_theta(tdims%i_start:tdims%i_end,        &
                                   tdims%j_start:tdims%j_end,        &
                                   tdims%k_start:tdims%k_end))
  av_r_theta      =    SUM(r_theta(tdims%i_start:tdims%i_end,        &
                                   tdims%j_start:tdims%j_end,        &
                                   tdims%k_start:tdims%k_end))/fieldsize
  has_nan_r_theta = 0.0
  IF(um_has_nan(r_theta(tdims%i_start:tdims%i_end,                   &
                        tdims%j_start:tdims%j_end,                   &
                        tdims%k_start:tdims%k_end))) has_nan_r_theta = 1.0

  max_r_m_v       = MAXVAL(r_m_v(tdims%i_start:tdims%i_end,          &
                                 tdims%j_start:tdims%j_end,          &
                                 tdims%k_start:tdims%k_end))
  min_r_m_v       = MINVAL(r_m_v(tdims%i_start:tdims%i_end,          &
                                 tdims%j_start:tdims%j_end,          &
                                 tdims%k_start:tdims%k_end))
  av_r_m_v        =    SUM(r_m_v(tdims%i_start:tdims%i_end,          &
                                 tdims%j_start:tdims%j_end,          &
                                 tdims%k_start:tdims%k_end))/fieldsize
  has_nan_r_m_v = 0.0
  IF(um_has_nan(r_m_v(tdims%i_start:tdims%i_end,                     &
                      tdims%j_start:tdims%j_end,                     &
                      tdims%k_start:tdims%k_end))) has_nan_r_m_v = 1.0

  max_r_m_cl      = MAXVAL(r_m_cl(tdims%i_start:tdims%i_end,         &
                                  tdims%j_start:tdims%j_end,         &
                                  tdims%k_start:tdims%k_end))
  min_r_m_cl      = MINVAL(r_m_cl(tdims%i_start:tdims%i_end,         &
                                  tdims%j_start:tdims%j_end,         &
                                  tdims%k_start:tdims%k_end))
  av_r_m_cl       =    SUM(r_m_cl(tdims%i_start:tdims%i_end,         &
                                  tdims%j_start:tdims%j_end,         &
                                  tdims%k_start:tdims%k_end))/fieldsize
  has_nan_r_m_cl = 0.0
  IF(um_has_nan(r_m_cl(tdims%i_start:tdims%i_end,                    &
                       tdims%j_start:tdims%j_end,                    &
                       tdims%k_start:tdims%k_end))) has_nan_r_m_cl = 1.0

  max_r_m_cf      = MAXVAL(r_m_cf(tdims%i_start:tdims%i_end,         &
                                  tdims%j_start:tdims%j_end,         &
                                  tdims%k_start:tdims%k_end))
  min_r_m_cf      = MINVAL(r_m_cf(tdims%i_start:tdims%i_end,         &
                                  tdims%j_start:tdims%j_end,         &
                                  tdims%k_start:tdims%k_end))
  av_r_m_cf       =    SUM(r_m_cf(tdims%i_start:tdims%i_end,         &
                                  tdims%j_start:tdims%j_end,         &
                                  tdims%k_start:tdims%k_end))/fieldsize
  has_nan_r_m_cf = 0.0
  IF(um_has_nan(r_m_cf(tdims%i_start:tdims%i_end,                    &
                       tdims%j_start:tdims%j_end,                    &
                       tdims%k_start:tdims%k_end))) has_nan_r_m_cf = 1.0

  IF (l_mcr_qcf2) THEN
    max_r_m_cf2   = MAXVAL(r_m_cf2(tdims%i_start:tdims%i_end,        &
                                   tdims%j_start:tdims%j_end,        &
                                   tdims%k_start:tdims%k_end))
    min_r_m_cf2   = MINVAL(r_m_cf2(tdims%i_start:tdims%i_end,        &
                                   tdims%j_start:tdims%j_end,        &
                                   tdims%k_start:tdims%k_end))
    av_r_m_cf2    =    SUM(r_m_cf2(tdims%i_start:tdims%i_end,        &
                                   tdims%j_start:tdims%j_end,        &
                                   tdims%k_start:tdims%k_end))/fieldsize
    has_nan_r_m_cf2 = 0.0
    IF(um_has_nan(r_m_cf2(tdims%i_start:tdims%i_end,                 &
                          tdims%j_start:tdims%j_end,                 &
                          tdims%k_start:tdims%k_end))) has_nan_r_m_cf2 = 1.0

  END IF
  IF (l_mcr_qrain) THEN
    max_r_m_rain  = MAXVAL(r_m_r(tdims%i_start:tdims%i_end,          &
                                 tdims%j_start:tdims%j_end,          &
                                 tdims%k_start:tdims%k_end))
    min_r_m_rain  = MINVAL(r_m_r(tdims%i_start:tdims%i_end,          &
                                 tdims%j_start:tdims%j_end,          &
                                 tdims%k_start:tdims%k_end))
    av_r_m_rain   =    SUM(r_m_r(tdims%i_start:tdims%i_end,          &
                                 tdims%j_start:tdims%j_end,          &
                                 tdims%k_start:tdims%k_end))/fieldsize
    has_nan_r_m_rain = 0.0
    IF(um_has_nan(r_m_r(tdims%i_start:tdims%i_end,                   &
                        tdims%j_start:tdims%j_end,                   &
                        tdims%k_start:tdims%k_end))) has_nan_r_m_rain = 1.0

  END IF
  IF (l_mcr_qgraup) THEN
    max_r_m_graup = MAXVAL(r_m_gr(tdims%i_start:tdims%i_end,         &
                                  tdims%j_start:tdims%j_end,         &
                                  tdims%k_start:tdims%k_end))
    min_r_m_graup = MINVAL(r_m_gr(tdims%i_start:tdims%i_end,         &
                                  tdims%j_start:tdims%j_end,         &
                                  tdims%k_start:tdims%k_end))
    av_r_m_graup  =    SUM(r_m_gr(tdims%i_start:tdims%i_end,         &
                                  tdims%j_start:tdims%j_end,         &
                                  tdims%k_start:tdims%k_end))/fieldsize
    has_nan_r_m_graup = 0.0
    IF(um_has_nan(r_m_gr(tdims%i_start:tdims%i_end,                  &
                         tdims%j_start:tdims%j_end,                  &
                         tdims%k_start:tdims%k_end))) has_nan_r_m_graup = 1.0

  END IF

  CALL gc_rmax(1,nproc,ierr,max_r_theta)
  CALL gc_rmin(1,nproc,ierr,min_r_theta)
  CALL gc_rmax(1,nproc,ierr,has_nan_r_theta)
  CALL gc_rsum(1,nproc,ierr, av_r_theta)
  CALL gc_rmax(1,nproc,ierr,max_r_m_v)
  CALL gc_rmin(1,nproc,ierr,min_r_m_v)
  CALL gc_rmax(1,nproc,ierr,has_nan_r_m_v)
  CALL gc_rsum(1,nproc,ierr, av_r_m_v)
  CALL gc_rmax(1,nproc,ierr,max_r_m_cl)
  CALL gc_rmin(1,nproc,ierr,min_r_m_cl)
  CALL gc_rmax(1,nproc,ierr,has_nan_r_m_cl)
  CALL gc_rsum(1,nproc,ierr, av_r_m_cl)
  CALL gc_rmax(1,nproc,ierr,max_r_m_cf)
  CALL gc_rmin(1,nproc,ierr,min_r_m_cf)
  CALL gc_rmax(1,nproc,ierr,has_nan_r_m_cf)
  CALL gc_rsum(1,nproc,ierr, av_r_m_cf)
  IF (l_mcr_qcf2) THEN
    CALL gc_rmax(1,nproc,ierr,max_r_m_cf2)
    CALL gc_rmin(1,nproc,ierr,min_r_m_cf2)
    CALL gc_rmax(1,nproc,ierr,has_nan_r_m_cf2)
    CALL gc_rsum(1,nproc,ierr, av_r_m_cf2)
  END IF
  IF (l_mcr_qrain) THEN
    CALL gc_rmax(1,nproc,ierr,max_r_m_rain)
    CALL gc_rmin(1,nproc,ierr,min_r_m_rain)
    CALL gc_rmax(1,nproc,ierr,has_nan_r_m_rain)
    CALL gc_rsum(1,nproc,ierr, av_r_m_rain)
  END IF
  IF (l_mcr_qgraup) THEN
    CALL gc_rmax(1,nproc,ierr,max_r_m_graup)
    CALL gc_rmin(1,nproc,ierr,min_r_m_graup)
    CALL gc_rmax(1,nproc,ierr,has_nan_r_m_graup)
    CALL gc_rsum(1,nproc,ierr, av_r_m_graup)
  END IF

  IF (mype == 0) THEN
    WRITE(umMessage,FMT='(A,4E32.16)') 'r_thetav :',min_r_theta,max_r_theta,   &
                                            av_r_theta,has_nan_r_theta
    CALL umPrint(umMessage,src='eg_R')
    WRITE(umMessage,FMT='(A,4E32.16)') 'r_m_v    :',min_r_m_v  ,max_r_m_v,     &
                                              av_r_m_v ,has_nan_r_m_v
    CALL umPrint(umMessage,src='eg_R')
    WRITE(umMessage,FMT='(A,4E32.16)') 'r_m_cl   :',min_r_m_cl ,max_r_m_cl,    &
                                              av_r_m_cl,has_nan_r_m_cl
    CALL umPrint(umMessage,src='eg_R')
    WRITE(umMessage,FMT='(A,4E32.16)') 'r_m_cf   :',min_r_m_cf ,max_r_m_cf,    &
                                               av_r_m_cf,has_nan_r_m_cf
    CALL umPrint(umMessage,src='eg_R')
    IF (l_mcr_qcf2) THEN
      WRITE(umMessage,FMT='(A,4E32.16)') 'r_m_cf2  :',min_r_m_cf2,   &
          max_r_m_cf2,   &
          av_r_m_cf2,has_nan_r_m_cf2
      CALL umPrint(umMessage,src='eg_R')
    END IF
    IF (l_mcr_qgraup) THEN
      WRITE(umMessage,FMT='(A,4E32.16)') 'r_m_graup:',min_r_m_graup, &
          max_r_m_graup ,&
          av_r_m_graup,has_nan_r_m_graup
      CALL umPrint(umMessage,src='eg_R')
    END IF
    IF (l_mcr_qrain) THEN
      WRITE(umMessage,FMT='(A,4E32.16)') 'r_m_rain :',min_r_m_rain , &
          max_r_m_rain , &
          av_r_m_rain,has_nan_r_m_rain
      CALL umPrint(umMessage,src='eg_R')
    END IF

    CALL umPrint( '============================================'// &
        '========================================',src='eg_R')
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_r
END MODULE eg_r_mod
