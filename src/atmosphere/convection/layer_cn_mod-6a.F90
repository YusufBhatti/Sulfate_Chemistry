! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Calculates layer dependent constants

MODULE layer_cn_6a_mod

IMPLICIT NONE

! mutually exclusive convection indicator types.
INTEGER, PARAMETER ::  deep = 1          ! indicator all points are deep
INTEGER, PARAMETER ::  congestus = 2     ! indicator all points are congestus
INTEGER, PARAMETER ::  shallow = 3       ! indicator all points are shallow
INTEGER, PARAMETER ::  midlevel = 4      ! indicator all points are mid
!
! Description:
!   Calculates the following layer dependent constants:
!   * pressure (k, K+1/2 and k+1)
!   * layer thickness (k, k+1/2, k+1)
!   * entrainment coefficients (k+1/4, k+3/4)
!   * detrainment coefficients (k)
!
! Method:
!      See Unified Model documentation paper 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LAYER_CN_6A_MOD'

CONTAINS

SUBROUTINE layer_cn_6a(k, npnts, nlev,                                  &
                       mdet_on,                                         &
                       ntml, ntpar, start_lev,                          &
                       exner_layer_boundaries, exner_layer_centres,     &
                       p_layer_boundaries, p_layer_centres,             &
                       z_rho,                                           &
                       conv_prog_precip,                                &
                       recip_pstar, entrain_coef, rhum,                 &
                       zk, zkp12, zkp1,                                 &
                       thek, qek, thekp1, qekp1,                        &
                       thpk, qpk,                                       &
                       wsc_o_mb, qsat_lcl, w_max,                       &
                       conv_indicator,                                  &
                       bconv,                                           &
                       ! Out
                       pk, pkp1, exk, exkp1,                            &
                       delpk, delpkp12, delpkp1,                        &
                       delp_uv_k, delp_uv_kp1,                          &
                       ekp14, ekp34, amdetk                             &
                       )

USE planet_constants_mod, ONLY: cp, g
USE cv_run_mod, ONLY:                                                   &
    ent_fac_dp, ent_fac_md, amdet_fac, ent_opt_dp, ent_opt_md,          &
    ent_dp_power, ent_md_power, mdet_opt_dp, icvdiag, cldbase_opt_sh,   &
    w_cape_limit,                                                       &
    prog_ent_grad, prog_ent_int, prog_ent_max, prog_ent_min 
USE cv_param_mod, ONLY: ae2, refdepth_dp, refqsat, sh_wstar_closure,    &
    sh_grey_closure, entcoef
USE cv_dependent_switch_mod, ONLY: l_var_entrain, l_new_det,            &
    l_const_ent

USE water_constants_mod, ONLY: lc, lf
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! Variables which are input
!----------------------------------------------------------------------
! Vector lengths and loop counters

INTEGER,INTENT(IN) :: k             ! present model layer
INTEGER,INTENT(IN) :: npnts         ! Number of points
INTEGER,INTENT(IN) :: nlev          ! Number of model levels for calculations
INTEGER,INTENT(IN) :: mdet_on       ! flag for adaptive mixing detrainment
                                    !  on = 1, off = 0
INTEGER,INTENT(IN) :: ntml(npnts)   ! Number of levels in the surface-based
                                    ! turbulently mixed layer
INTEGER,INTENT(IN) :: ntpar(npnts)  ! Top of initial parcel ascent
INTEGER,INTENT(IN) :: start_lev(npnts)  ! Initiation level

! Field on model levels
REAL,INTENT(IN) :: exner_layer_boundaries(npnts,0:nlev) ! Exner ratio at layer
                                                        ! boundary starting
                                                        ! at level k-1/2
REAL,INTENT(IN) :: exner_layer_centres(npnts,0:nlev)    ! Exner function
                                                        ! at layer centre
REAL,INTENT(IN) :: p_layer_centres(npnts,0:nlev)        ! Pressure
                                                        ! at layer centre (Pa)
REAL,INTENT(IN) :: p_layer_boundaries(npnts,0:nlev)     ! Pressure
                                                        ! at layer boundary (Pa)
REAL,INTENT(IN) :: z_rho(npnts,nlev)                    ! height of rho levels
                                                        ! (m)
REAL,INTENT(IN) :: conv_prog_precip(npnts,nlev) ! Surface precipitation based
                                                ! 3d convective prognostic in
                                                ! kg/m2/s 
! Fields on a single level
REAL,INTENT(IN) :: recip_pstar(npnts) ! Reciprocal of pstar array (1/Pa)
REAL,INTENT(IN) :: entrain_coef(npnts)! entrainment coefficients
REAL,INTENT(IN) :: rhum(npnts)        ! Relative humidity at level K
REAL,INTENT(IN) :: zk(npnts)          ! height on k
REAL,INTENT(IN) :: zkp12(npnts)       ! height on k+1/2
REAL,INTENT(IN) :: zkp1(npnts)        ! height on k+1
REAL,INTENT(IN) :: thek(npnts)        ! theta for environment on k
REAL,INTENT(IN) :: qek(npnts)         ! q for environment on k
REAL,INTENT(IN) :: thekp1(npnts)      ! theta for environment on k+1
REAL,INTENT(IN) :: qekp1(npnts)       ! q for environment on k+1
REAL,INTENT(IN) :: thpk(npnts)        ! parcel theta on k
REAL,INTENT(IN) :: qpk(npnts)         ! parcel q on k
REAL,INTENT(IN) :: wsc_o_mb(npnts)    ! Convective velocity scale divided
                                      ! by cloud base mass flux mb
REAL,INTENT(IN) :: qsat_lcl(npnts)    ! qsat at the LCL
REAL,INTENT(IN) :: w_max(npnts)       ! maximum large-scale w in column

INTEGER,INTENT(IN) :: conv_indicator  ! type of convection for layer calls.
LOGICAL,INTENT(IN) :: bconv(npnts)    ! Mask for points at which
                                      ! convection is occurring

!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------
REAL,INTENT(OUT) :: pk(npnts)         ! pressure at mid-point of layer k (Pa)
REAL,INTENT(OUT) :: pkp1(npnts)       ! pressure at mid-point of layer k+1 (Pa)
REAL,INTENT(OUT) :: exk(npnts)        ! Exner ratio at mid-point of layer k
REAL,INTENT(OUT) :: exkp1(npnts)      ! Exner ratio at mid-point of layer k+1
REAL,INTENT(OUT) :: delpk(npnts)      ! pressure difference across layer k (Pa)
REAL,INTENT(OUT) :: delpkp12(npnts)   ! pressure diff. across layer k+1/2 (Pa)
REAL,INTENT(OUT) :: delpkp1(npnts)    ! pressure diff. across layer k+1 (Pa)
REAL,INTENT(OUT) :: delp_uv_k(npnts)  ! pressure difference across UV
                                      ! layer k (Pa)
REAL,INTENT(OUT) :: delp_uv_kp1(npnts)! pressure difference across UV
                                      ! layer k+1 (Pa)

REAL,INTENT(OUT) :: ekp14(npnts)    ! Entrainment coefficient at level k+1/4
                                    ! multiplied by appropriate layer thickness
REAL,INTENT(OUT) :: ekp34(npnts)    ! Entrainment coefficient at level k+3/4
                                    ! multiplied by appropriate layer thickness
REAL,INTENT(OUT) :: amdetk(npnts)   ! Mixing detrainment coefficient at level k
                                    ! multiplied by appropriate layer thickness

!----------------------------------------------------------------------
! Local variables
!----------------------------------------------------------------------

INTEGER :: i              ! loop counter

REAL :: aekp14            ! Used in calculation of entrainment rate
REAL :: aekp34            ! Used in calculation of entrainment rate
REAL :: pkp12(npnts)      ! Pressure at upper boundary of layer k (Pa)
REAL :: pntml(npnts)      ! Pressure at upper boundary of layer ntml (Pa)
REAL :: delp_cld(npnts)   ! thickness of cloud layer (Pa)
REAL :: delz_cld(npnts)   ! thickness of cloud layer (m)
REAL :: delpkp14(npnts)   ! thickness of layer for k+1/4 (Pa)
REAL :: delpkp34(npnts)   ! thickness of layer for k+3/4 (Pa)
REAL :: amdet_func(npnts) ! detrainment function of RH
REAL :: tmp_entrain_coef(npnts)! entrainment coefficients
REAL :: r_mdet_depth      ! Reciprocal of the scale depth used in variable
                          ! mixing detrainment calculation (m^-1)
REAL :: mdet_fac          ! constant used in mixing detrainment calculation
                          ! (unitless)
REAL :: factor
REAL :: factor1
REAL :: factor2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LAYER_CN_6A'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Initialise entrainment and detrainment rates to zero and hence they
! will default to 0.0 if not explicitly set.
!----------------------------------------------------------------------
DO i=1,npnts
  ekp14(i)  = 0.0
  ekp34(i)  = 0.0
  amdetk(i) = 0.0
END DO

!----------------------------------------------------------------------
! Set constant ae used in calculation of entrainment and detrainment
! rates depending upon level.
!----------------------------------------------------------------------
IF (conv_indicator == deep) THEN
  !Deep only
  aekp14 = ent_fac_dp*ae2
  aekp34 = ent_fac_dp*ae2
ELSE
  !Used only for mid-level
  aekp14 = ent_fac_md*ae2
  aekp34 = ent_fac_md*ae2
END IF

!---------------------------------------------------------------------
! Calculate pressurea and pressure thicknesses
! NB pressure differences are calculated at k+1 - k to give positive
! thicknesses.
!---------------------------------------------------------------------
DO i=1,npnts
  pk(i)           = p_layer_centres(i,k)
  pkp12(i)        = p_layer_boundaries(i,k)
  pkp1(i)         = p_layer_centres(i,k+1)
  exk(i)          = exner_layer_centres(i,k)
  exkp1(i)        = exner_layer_centres(i,k+1)
  delpk(i)        = p_layer_boundaries(i,k-1)   - p_layer_boundaries(i,k)
  delpkp1(i)      = p_layer_boundaries(i,k)     - p_layer_boundaries(i,k+1)
  delpkp12(i)     = pk(i)                       - pkp1(i)
  delpkp14(i)     = pk(i)                       - p_layer_boundaries(i,k)
  delpkp34(i)     = p_layer_boundaries(i,k)     - pkp1(i)
  delp_uv_k(i)    = p_layer_centres(i,k-1)      - p_layer_centres(i,k)
  delp_uv_kp1(i)  = p_layer_centres(i,k)        - p_layer_centres(i,k+1)
END DO

IF (conv_indicator == shallow) THEN
  ! Cloud thickness  - only used for shallow convection
  DO i=1,npnts
    delp_cld(i)  = p_layer_boundaries(i,ntml(i)) -                      &
                   p_layer_boundaries(i,ntpar(i))
    pntml(i)     = p_layer_boundaries(i,ntml(i))
  END DO
END IF

IF ((conv_indicator == deep) .AND. (mdet_opt_dp == 2 .OR. ent_opt_dp == 4  &
              .OR.  ent_opt_dp == 5)) THEN
  ! Cloud thickness  - only used for deep convection with particular options.
  DO i=1,npnts
    delz_cld(i)  = z_rho(i,ntpar(i)+1)-z_rho(i,ntml(i)+1)
  END DO
END IF

! ---------------------------------------------------------------------
! Calculate entrainment coefficients multiplied by approppriate
! layer thickness.
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Shallow convection entrainment coefficients
! ---------------------------------------------------------------------
! Exponential decrease of entrainment rate with height above NTPAR
! is stopped to ensure massflux continues to decrease
! significantly with height (characteristic of shallow)

IF (conv_indicator == shallow) THEN
  IF (cldbase_opt_sh == sh_wstar_closure) THEN
    ! Original entrainment rates
    DO i=1,npnts
      IF (k >= ntml(i)) THEN
        ekp14(i) = wsc_o_mb(i) * delpkp14(i) * 0.03 * EXP( -1.0*MIN( 1.0, &
                  (pntml(i)-pk(i)   )/delp_cld(i) ) ) / delp_cld(i)
        ekp34(i) = wsc_o_mb(i) * delpkp34(i) * 0.03 * EXP( -1.0*MIN( 1.0, &
                  (pntml(i)-pkp12(i))/delp_cld(i) ) ) / delp_cld(i)
      END IF    ! level
      IF (k == ntml(i)) ekp14(i) = 0.0
    END DO
  ELSE IF (cldbase_opt_sh == sh_grey_closure) THEN
    ! Increased entrainment for grey zone
    DO i=1,npnts
      IF (k >= ntml(i)) THEN
        ekp14(i) = wsc_o_mb(i) * delpkp14(i) * 0.06 * EXP( -1.5*MIN( 1.0, &
                  (pntml(i)-pk(i)   )/delp_cld(i) ) ) / delp_cld(i)
        ekp34(i) = wsc_o_mb(i) * delpkp34(i) * 0.06 * EXP( -1.5*MIN( 1.0, &
                  (pntml(i)-pkp12(i))/delp_cld(i) ) ) / delp_cld(i)
      END IF    ! level
      IF (k == ntml(i)) ekp14(i) = 0.0
    END DO
  END IF
  ! ---------------------------------------------------------------------
  ! Congestus convection entrainment coefficients
  ! ---------------------------------------------------------------------
  ! 0.5/z rates
ELSE IF (conv_indicator == congestus) THEN
  DO i=1,npnts
    IF (k > ntml(i)) THEN
      ! Is there any reason why k=ntml isn't set?
      ekp34(i) = 0.5*2.0 * (zkp1(i)-zkp12(i))/(zkp12(i)+zkp1(i))
      ekp14(i) = 0.5*2.0 * (zkp12(i) - zk(i))/(zkp12(i) + zk(i))
    END IF    ! type of convection and level
  END DO

  ! ---------------------------------------------------------------------
  ! Deep convection entrainment coefficients
  ! ---------------------------------------------------------------------
ELSE IF (conv_indicator == deep) THEN

  SELECT CASE (ent_opt_dp)
  CASE (0)     ! original Ap/(p*)^2  style entrainment
    DO i=1,npnts
      ekp14(i) = entcoef * aekp14 * pk(i)    * delpkp14(i) *            &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = entcoef * aekp34 * pkp12(i) * delpkp34(i) *            &
                 recip_pstar(i) * recip_pstar(i)
    END DO
  CASE (1)     !  n/z style entrainment  where n = ent_fac
    DO i=1,npnts
      ekp14(i) = ent_fac_dp *2.0 * (zkp12(i) - zk(i))/(zkp12(i) + zk(i))
      ekp34(i) = ent_fac_dp *2.0 * (zkp1(i)-zkp12(i))/(zkp12(i)+zkp1(i))
    END DO

  CASE (2)     ! Higher entrainment near surface (various versions have
              ! been used during GA3.0 plus testing) current code has
              ! version known as New 3.
              ! factor * original Ap/(p*)^2
    DO i=1,npnts
      IF (pk(i) > 50000.0) THEN
        factor1 = 1.0 + 1.25*(1.0-(100000.0 - pk(i))/50000.0)
      ELSE
        factor1 = 1.0
      END IF
      IF (pkp12(i) > 50000.0) THEN
        factor2 = 1.0 + 1.25*(1.0-(100000.0 - pkp12(i))/50000.0)
      ELSE
        factor2 = 1.0
      END IF
      ekp14(i) = factor1 * entcoef * aekp14 * pk(i)    * delpkp14(i) *  &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = factor2 * entcoef * aekp34 * pkp12(i) * delpkp34(i) *  &
                 recip_pstar(i) * recip_pstar(i)
    END DO

  CASE (3)     ! factor * (A/p*)*((p/p*)^m)  style entrainment
    DO i=1,npnts
      ekp14(i) = entcoef * aekp14 * delpkp14(i) * recip_pstar(i) *      &
                (pk(i)    * recip_pstar(i))**ent_dp_power
      ekp34(i) = entcoef * aekp34 * delpkp34(i) * recip_pstar(i) *      &
                (pkp12(i) * recip_pstar(i))**ent_dp_power
    END DO

  CASE (4)     ! variable n/z style entrainment
    DO i=1,npnts
      tmp_entrain_coef(i) = refdepth_dp/delz_cld(i)
    END DO
    DO i=1,npnts
      ekp34(i) = tmp_entrain_coef(i) *2.0 * (zkp1(i)-zkp12(i))          &
                               / (zkp12(i)+zkp1(i))
      ekp14(i) = tmp_entrain_coef(i) *2.0 * (zkp12(i) - zk(i))          &
                               / (zkp12(i) + zk(i))
    END DO

  CASE (5)     ! variable p/(p*)^2  style entrainment
    DO i=1,npnts
      tmp_entrain_coef(i) = refdepth_dp/delz_cld(i)
    END DO
    DO i=1,npnts
      ekp14(i) = tmp_entrain_coef(i) *                                  &
                 entcoef * aekp14 * pk(i)    * delpkp14(i) *            &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = tmp_entrain_coef(i) *                                  &
                 entcoef * aekp34 * pkp12(i) * delpkp34(i) *            &
                 recip_pstar(i) * recip_pstar(i)
    END DO

  CASE(6)     ! variable p/(p*)^2 entrainment dependent on precipitation
              ! based 3d convective prognostic at the initiation level.
    DO i=1,npnts 
      IF (w_max(i) < w_cape_limit) THEN
        IF (conv_prog_precip(i,start_lev(i)) > 1.e-10) THEN
          tmp_entrain_coef(i) = prog_ent_grad *                         &
                                LOG10(conv_prog_precip(i,start_lev(i)) *&
                                refqsat/qsat_lcl(i)) + prog_ent_int
          tmp_entrain_coef(i) = MIN(                                    &
                                MAX(tmp_entrain_coef(i), prog_ent_min), &
                                prog_ent_max)
        ELSE
          tmp_entrain_coef(i) = prog_ent_max
        END IF
      ELSE
        ! If w_cape_limit is exceeded turn off entrainment scaling
        tmp_entrain_coef(i) = 1.0
      END IF
    END DO
    DO i=1,npnts 
      ekp14(i) = tmp_entrain_coef(i) *                                  &
                 entcoef * aekp14 * pk(i)    * delpkp14(i) *            &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = tmp_entrain_coef(i) *                                  &
                 entcoef * aekp34 * pkp12(i) * delpkp34(i) *            &
                 recip_pstar(i) * recip_pstar(i)
    END DO

  CASE(7)     ! variable p/(p*)^2 entrainment dependent on precipitation
              ! based 3d convective prognostic at the current level.
    DO i=1,npnts 
      IF (w_max(i) < w_cape_limit) THEN
        IF (conv_prog_precip(i,k) > 1.e-10) THEN
          tmp_entrain_coef(i) = prog_ent_grad *                         &
                                LOG10(conv_prog_precip(i,k) *           &
                                refqsat/qsat_lcl(i)) + prog_ent_int
          tmp_entrain_coef(i) = MIN(                                    &
                                MAX(tmp_entrain_coef(i), prog_ent_min), &
                                prog_ent_max)
        ELSE
          tmp_entrain_coef(i) = prog_ent_max
        END IF
      ELSE
        ! If w_cape_limit is exceeded turn off entrainment scaling
        tmp_entrain_coef(i) = 1.0
      END IF
    END DO
    DO i=1,npnts 
      ekp14(i) = tmp_entrain_coef(i) *                                  &
                 entcoef * aekp14 * pk(i)    * delpkp14(i) *            &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = tmp_entrain_coef(i) *                                  &
                 entcoef * aekp34 * pkp12(i) * delpkp34(i) *            &
                 recip_pstar(i) * recip_pstar(i)
    END DO
  END SELECT  ! test on ent_opt_dp

  ! ---------------------------------------------------------------------
  ! Mid-Level convection entrainment coefficients
  ! ---------------------------------------------------------------------
ELSE

  SELECT CASE (ent_opt_md)
  CASE (0)     ! orignal Ap/(p*)^2  style entrainment

    DO i=1,npnts
      ekp14(i) = entcoef * aekp14 * pk(i)    * delpkp14(i) *            &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = entcoef * aekp34 * pkp12(i) * delpkp34(i) *            &
                 recip_pstar(i) * recip_pstar(i)
    END DO

  CASE (1)     !  n/z style entrainment  where n = ent_fac

    DO i=1,npnts
      ekp34(i) = ent_fac_md *2.0 * (zkp1(i)-zkp12(i))/(zkp12(i)+zkp1(i))
      ekp14(i) = ent_fac_md *2.0 * (zkp12(i) - zk(i))/(zkp12(i) + zk(i))
    END DO

  CASE (2)     ! New 3 profile - higher near surface (option not used
              ! during GA3.0 plus testing).
    DO i=1,npnts
      IF (pk(i) > 50000.0) THEN
        factor1 = 1.0 + 1.25*(1.0-(100000.0 - pk(i))/50000.0)
      ELSE
        factor1 = 1.0
      END IF
      IF (pkp12(i) > 50000.0) THEN
        factor2 = 1.0 + 1.25*(1.0-(100000.0 - pkp12(i))/50000.0)
      ELSE
        factor2 = 1.0
      END IF
      ekp14(i) = factor1 * entcoef * aekp14 * pk(i)    * delpkp14(i) *  &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = factor2 * entcoef * aekp34 * pkp12(i) * delpkp34(i) *  &
                 recip_pstar(i) * recip_pstar(i)
    END DO

  CASE (3)     ! factor * (A/p*)*((p/p*)^m)  style entrainment
    DO i=1,npnts
      ekp14(i) = entcoef * aekp14 * delpkp14(i) * recip_pstar(i) *      &
                (pk(i)    * recip_pstar(i))**ent_md_power
      ekp34(i) = entcoef * aekp34 * delpkp34(i) * recip_pstar(i) *      &
                (pkp12(i) * recip_pstar(i))**ent_md_power
    END DO

  CASE(6)     ! variable p/(p*)^2 entrainment dependent on precipitation
              ! based 3d convective prognostic at the initiation level.
    DO i=1,npnts 
      IF (w_max(i) < w_cape_limit) THEN
        IF (conv_prog_precip(i,start_lev(i)) > 1.e-10) THEN
          tmp_entrain_coef(i) = prog_ent_grad *                         &
                                LOG10(conv_prog_precip(i,start_lev(i)) *&
                                refqsat/qsat_lcl(i)) + prog_ent_int
          tmp_entrain_coef(i) = MIN(                                    &
                                MAX(tmp_entrain_coef(i), prog_ent_min), &
                                prog_ent_max)
        ELSE
          tmp_entrain_coef(i) = prog_ent_max
        END IF
      ELSE
        ! If w_cape_limit is exceeded turn off entrainment scaling
        tmp_entrain_coef(i) = 1.0
      END IF
    END DO
    DO i=1,npnts 
      ekp14(i) = tmp_entrain_coef(i) *                                  &
                 entcoef * aekp14 * pk(i)    * delpkp14(i) *            &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = tmp_entrain_coef(i) *                                  &
                 entcoef * aekp34 * pkp12(i) * delpkp34(i) *            &
                 recip_pstar(i) * recip_pstar(i)
    END DO

  CASE(7)     ! variable p/(p*)^2 entrainment dependent on precipitation
              ! based 3d convective prognostic at the current level.
    DO i=1,npnts 
      IF (w_max(i) < w_cape_limit) THEN
        IF (conv_prog_precip(i,k) > 1.e-10) THEN
          tmp_entrain_coef(i) = prog_ent_grad *                         &
                                LOG10(conv_prog_precip(i,k) *           &
                                refqsat/qsat_lcl(i)) + prog_ent_int
          tmp_entrain_coef(i) = MIN(                                    &
                                MAX(tmp_entrain_coef(i), prog_ent_min), &
                                prog_ent_max)
        ELSE
          tmp_entrain_coef(i) = prog_ent_max
        END IF
      ELSE
        ! If w_cape_limit is exceeded turn off entrainment scaling
        tmp_entrain_coef(i) = 1.0
      END IF
    END DO
    DO i=1,npnts 
      ekp14(i) = tmp_entrain_coef(i) *                                  &
                 entcoef * aekp14 * pk(i)    * delpkp14(i) *            &
                 recip_pstar(i) * recip_pstar(i)
      ekp34(i) = tmp_entrain_coef(i) *                                  &
                 entcoef * aekp34 * pkp12(i) * delpkp34(i) *            &
                 recip_pstar(i) * recip_pstar(i)
    END DO

  END SELECT  ! test on ent_opt

END IF        ! type of convection


! ---------------------------------------------------------------------
! If variable entrainment then overwrite previously calculated
! entrainment values when there is a valid (positive) entrain_coef
! entrain_coef is calculated in the diagnosis (conv_diag_comp).
! ---------------------------------------------------------------------
IF ( (icvdiag == 4) .OR. (icvdiag == 5) ) THEN
  IF ( (conv_indicator == shallow) .OR.     &
       (conv_indicator == congestus) .OR.   &
       (conv_indicator == deep) ) THEN
    IF (l_const_ent) THEN         ! no height dependence
      DO i=1,npnts
        IF (entrain_coef(i)  > 0.0) THEN
          ekp14(i) = entrain_coef(i) * (zkp12(i) - zk(i))
          ekp34(i) = entrain_coef(i) * (zkp1(i)-zkp12(i))
        END IF
        IF ((conv_indicator == shallow) .OR.  &
            (conv_indicator == congestus )) THEN
          IF (k == ntml(i)) ekp14(i) = 0.0
        END IF
      END DO
    ELSE      ! rates vary with height
      DO i=1,npnts
        IF (entrain_coef(i)  > 0.0) THEN

          factor = 1.0

          ekp34(i) = entrain_coef(i) *2.0 * (zkp1(i)-zkp12(i))          &
                                  *factor / (zkp12(i)+zkp1(i))
          ekp14(i) = entrain_coef(i) *2.0 * (zkp12(i) - zk(i))          &
                                  *factor / (zkp12(i) + zk(i))
        END IF    ! test on entrain_coef
        IF ((conv_indicator == shallow) .OR.                           &
            (conv_indicator == congestus )) THEN
          IF (k == ntml(i)) ekp14(i) = 0.0
        END IF
      END DO
    END IF    ! test on l_const_ent
  END IF      ! conv type
END IF


! ---------------------------------------------------------------------
! Calculate mixing detrainment coefficient multiplied by appropriate
! layer thickness.
! ---------------------------------------------------------------------
! NB the mixing detrainment is always zero for the first convecting
! level.

! ---------------------------------------------------------------------
! Shallow convection mixing detrainment coefficients
! ---------------------------------------------------------------------
IF (conv_indicator == shallow) THEN

  IF (l_new_det) THEN    ! Use new relationship for detrainment
    DO i=1,npnts
      IF ( k > ntml(i) ) THEN
        IF (entrain_coef(i) > 0.0) THEN    ! alter detrainment
          ! Trying 2* entrainment rates for shallow
          amdetk(i) = (1.0 + 1.0) * (ekp14(i) + ekp34(i))
          !    But set to 1 if RH is greater than 100%
          IF (rhum(i)  >   1.0) THEN
            amdetk(i) = 1.0 * (ekp14(i) + ekp34(i))
          END IF
        END IF
      END IF
    END DO    ! npnts
  ELSE IF (cldbase_opt_sh == sh_grey_closure) THEN
    ! Increased entrainment for grey zone
    DO i=1,npnts
      IF (k  > ntml(i) ) THEN
        amdetk(i) = 1.5*(ekp14(i) + ekp34(i))
      END IF
    END DO    ! npnts
  ELSE  !Original detrainment rates
    DO i=1,npnts
      IF ( k > ntml(i) ) THEN
        IF (rhum(i)  <=  0.85) THEN
          amdetk(i) = (1.0 + 0.3) * (ekp14(i) + ekp34(i))
        ELSE IF (rhum(i)  >   1.0) THEN
          amdetk(i) =  1.0        * (ekp14(i) + ekp34(i))
        ELSE  ! 0.85 <  rhum <= 1.0
          amdetk(i) = (1.0 + (0.3/0.15) * (1.0-rhum(i)))                &
                                  * (ekp14(i) + ekp34(i))
        END IF
      END IF
    END DO    ! npnts
  END IF


  ! ---------------------------------------------------------------------
  ! Congestus convection mixing detrainment coefficients
  ! ---------------------------------------------------------------------
  ! adaptive mixing detrainment used in all cases
ELSE IF (conv_indicator == congestus) THEN
  DO i=1,npnts
    IF ( k > ntml(i) .AND. rhum(i) <= 1.0 ) THEN
      amdetk(i) = (1.0-rhum(i)) * (ekp14(i) + ekp34(i))
    END IF
  END DO !npnts

  ! ---------------------------------------------------------------------
  ! Deep convection mixing detrainment coefficients
  ! ---------------------------------------------------------------------
ELSE IF (conv_indicator == deep) THEN
  SELECT CASE (mdet_opt_dp)
  CASE (0)     ! Original (1.-1./(ent_fac_dp*ae2))
    DO i=1,npnts
      IF ( k > ntml(i)) THEN
        amdetk(i) = (ekp14(i) + ekp34(i)) * (1.0-1.0/aekp34)
      END IF
    END DO !npnts
  CASE (1)     ! Adaptive mixing detrainment
    DO i=1,npnts
      IF ( k > ntml(i) .AND. rhum(i) <= 1.0 ) THEN
        amdetk(i) = amdet_fac*(ekp14(i) + ekp34(i))*(1-rhum(i))
      END IF
    END DO !npnts
  CASE (2)     ! Variable mixing detrainment
    r_mdet_depth = 1.0/(refdepth_dp - 2000.0)
    mdet_fac     = 5.0/3.0
    DO i=1,npnts
      IF ( k > ntml(i) ) THEN
        factor    = MAX(mdet_fac - r_mdet_depth * delz_cld(i), 1.0/3.0)
        amdetk(i) = amdet_fac*(ekp14(i) + ekp34(i)) * factor
      END IF
    END DO !npnts
  END SELECT  ! mdet_dpt_dp

  ! ---------------------------------------------------------------------
  ! Mid-level convection mixing detrainment coefficients
  ! ---------------------------------------------------------------------
ELSE
  DO i=1,npnts
    IF (bconv(i)) THEN
      amdetk(i) = (ekp14(i) + ekp34(i)) * (1.0-1.0/aekp34)
    END IF
  END DO !npnts
END IF  ! type of convection scheme

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE layer_cn_6a

END MODULE layer_cn_6a_mod
