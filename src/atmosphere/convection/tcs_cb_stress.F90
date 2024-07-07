! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
MODULE tcs_cb_stress_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'TCS_CB_STRESS_MOD'
CONTAINS

SUBROUTINE tcs_cb_stress (conv_type, n_npts                              &
                          ,timestep                                      &
                          ,uw0, vw0, mb, zlcl                            &
                          ,du_start, dv_start, dz_cb, fng_nlclp1         &
                          ,uw_cb, vw_cb )

! Description:
!   Calculate cloud base uw and vw stresses using turbulence methods.
!

!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection

! Modules

USE conv_type_defs       ! convective type definitions
USE tcs_cmt_params_sh    ! shallow CMT parameters
USE tcs_cmt_params_cg    ! congetsus CMT parameters
USE tcs_cmt_params_dp    ! deep CMT parameters

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
IMPLICIT NONE


! Arguments:
INTEGER, INTENT(IN) :: &
  conv_type            & ! Indicator for type of convection
 ,n_npts                 ! number of convecting columns

REAL, INTENT(IN) ::    &
  timestep               ! Convection timestep (s)

REAL, INTENT(IN) ::    &
  uw0(n_npts)          & ! surface uw stress (N/m2)
 ,vw0(n_npts)          & ! surface vw stress (N/m2)
 ,mb(n_npts)           & ! Cloud base mass flux (Pa/s)
 ,zlcl(n_npts)         & ! Exact height of LCL (m) ? (or model level height)
 ,du_start(n_npts)     & ! Change in U across cloud base at start of step
 ,dv_start(n_npts)     & ! Change in V across cloud base at start of step
 ,dz_cb(n_npts)        & ! depth of layer for cloud base jump (m)
 ,fng_nlclp1(n_npts)     ! Value of non-gradient function at nlcl+1

REAL, INTENT(OUT) ::   &
  uw_cb(n_npts)        & ! uw at cloud base (N/m2)
 ,vw_cb(n_npts)          ! vw at cloud base (N/m2)

!-------------------------------------------------------------------------------
! Local declarations:

INTEGER ::  i       ! loop counter

REAL  ::          &
  gamma_cmt       & ! Constant for cloud base CMT
 ,delta_cmt       & ! Constant for cloud base CMT
 ,alpha           & ! Factor for in solution of du/dv
 ,beta_u          & ! Factor for in solution of du
 ,beta_v          & ! Factor for in solution of dv
 ,betau_r_alpha   & ! beta_u/alpha
 ,betav_r_alpha   & ! beta_v/alpha
 ,alpha_dt        & ! dt * alpha
 ,fng_scale         !

REAL  ::          &
  du_avg(n_npts)  & ! Average value of jump du through timestep
 ,dv_avg(n_npts)    ! Average value of jump dv through timestep

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TCS_CB_STRESS'



!-------------------------------------------------------------------------------
! Based in input flag decide on CMT parameters
!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
SELECT CASE (conv_type)       ! type of convection

CASE (shallow_conv)         ! Number 1
  gamma_cmt = gamma_cmt_shall
  delta_cmt = delta_cmt_shall

CASE (congestus_conv)       ! number 2
  gamma_cmt = gamma_cmt_shall
  delta_cmt = delta_cmt_shall

CASE (deep_conv)            ! number 3
  gamma_cmt = gamma_cmt_deep
  delta_cmt = delta_cmt_deep

  !     Case (mid_conv)             ! may be used in future ?
  !       gamma_cmt = gamma_cmt_mid
  !       delta_cmt = delta_cmt_mid

CASE DEFAULT
  WRITE(umMessage,*) ' tcs_cb_stress unallowed convection type ',conv_type
  CALL umPrint(umMessage,src='tcs_cb_stress')

END SELECT                    ! type of convection

!-------------------------------------------------------------------------------
! Calculate jumps in vorticity components across cloud base. Uses an implicit
! technique assuming du and dv vary as exp(-t/tau) through a timestep.
! An explicit treatment can lead to problems.
!
!   d(du)/dt = -alpha*du + beta
!
! Solution  du=du(t_start)*exp(-alpha*dt)+(beta/alpha)(1.-exp(-alpha*dt))
!
! Requires average solution over timestep
!
!    du_avg = integral (du(t))dt)/timestep
!
! du_avg = [(beta/alpha)-du(t_start)](exp(-alpha*dt)-1)/(dt*alpha)+beta/alpha
!-------------------------------------------------------------------------------
! Note this code does not exactly resemble the documentation for the old
! turbulence based CMT (see Section 5.5 of UM doc paper 27).
! This code uses a different definition of gamma and delta.
!  gamma_cmt = 1/delta    in UM doc
!  delta_cmt = (1 -gamma/delta)  in UM doc
!-------------------------------------------------------------------------------

DO i=1,n_npts

  fng_scale = (fng_nlclp1(i) -1.0)*(zlcl(i)/dz_cb(i))

  alpha = gamma_cmt*mb(i)*(1.0-fng_scale)/dz_cb(i)

  beta_u  = (delta_cmt*(1.0 - fng_scale) - 1.0)*uw0(i)/zlcl(i)
  beta_v  = (delta_cmt*(1.0 - fng_scale) - 1.0)*vw0(i)/zlcl(i)

  betau_r_alpha = beta_u/alpha
  betav_r_alpha = beta_v/alpha
  alpha_dt = alpha*timestep

  du_avg(i) = (betau_r_alpha - du_start(i))*(EXP(-alpha_dt)-1.0)/alpha_dt &
                + betau_r_alpha
  dv_avg(i) = (betav_r_alpha - dv_start(i))*(EXP(-alpha_dt)-1.0)/alpha_dt &
                 + betav_r_alpha

  ! cloud base stress components

  uw_cb(i) = delta_cmt*uw0(i) - gamma_cmt*mb(i)*zlcl(i)*du_avg(i)/dz_cb(i)
  vw_cb(i) = delta_cmt*vw0(i) - gamma_cmt*mb(i)*zlcl(i)*dv_avg(i)/dz_cb(i)

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

!-------------------------------------------------------------------------------

END SUBROUTINE tcs_cb_stress
END MODULE tcs_cb_stress_mod
