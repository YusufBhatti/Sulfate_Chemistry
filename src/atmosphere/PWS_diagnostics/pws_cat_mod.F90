! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine pws_cat_mod ---------------------------------------
!
!   Purpose: Calculates PWS CAT diags, stash section 20
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics
MODULE pws_cat_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_CAT_MOD'

CONTAINS

SUBROUTINE pws_cat(u,v,p, u_press, v_press, press_lev,k)

! Description
! Calculate the cat diagnostics, 16,17 and 18.
! based upon wind shears, using Dutton forumla

USE atm_fields_bounds_mod

USE nlsizes_namelist_mod,  ONLY: global_row_length
USE level_heights_mod,     ONLY: r_theta_levels,r_rho_levels
USE planet_constants_mod,  ONLY: kappa, p_zero, planet_radius
USE missing_data_mod,      ONLY: imdi, rmdi
USE nlsizes_namelist_mod,  ONLY: row_length, rows, model_levels

USE pws_diags_mod, ONLY: pws_cat_turb, flag_cat_turb

USE uc_to_ub_mod, ONLY: uc_to_ub
USE vc_to_vb_mod, ONLY: vc_to_vb
USE pc_to_pb_mod, ONLY: pc_to_pb

USE diff_mod, ONLY: diffp, diffx, diffy

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE um_parvars
USE Field_Types,  ONLY: fld_type_u, fld_type_v

IMPLICIT NONE

REAL, INTENT(IN)  :: u(udims_s%i_start:udims_s%i_end,                 &
                       udims_s%j_start:udims_s%j_end,                 &
                       udims_s%k_start:udims_s%k_end)


REAL, INTENT(IN) :: v(vdims_s%i_start:vdims_s%i_end,                  &
                      vdims_s%j_start:vdims_s%j_end,                  &
                      vdims_s%k_start:vdims_s%k_end)


REAL, INTENT(IN) :: p(pdims_s%i_start:pdims_s%i_end,                  &
                    pdims_s%j_start:pdims_s%j_end,                    &
                    pdims_s%k_start:pdims_s%k_end+1)

REAL, INTENT(IN) :: u_press(udims%i_start:udims%i_end,                &
                            vdims%j_start:vdims%j_end)

REAL, INTENT(IN) :: v_press(udims%i_start:udims%i_end,                &
                            vdims%j_start:vdims%j_end)

REAL, INTENT(IN) :: press_lev

INTEGER, INTENT(IN) :: k

! Local variables

REAL :: local_ub(udims%i_start:udims%i_end,                           &
                 vdims%j_start:vdims%j_end,                           &
                 udims%k_start:udims%k_end)

REAL :: local_vb(udims%i_start:udims%i_end,                           &
                 vdims%j_start:vdims%j_end,                           &
                 vdims%k_start:vdims%k_end)

REAL :: local_r_rho_levels_b(udims%i_start:udims%i_end,               &
                 vdims%j_start:vdims%j_end,                           &
                 pdims%k_start:pdims%k_end)

REAL :: local_pb(udims%i_start:udims%i_end,                           &
                 vdims%j_start:vdims%j_end,                           &
                 pdims%k_start:pdims%k_end)

REAL :: dzdp(udims%i_start:udims%i_end,                               &
                 vdims%j_start:vdims%j_end)

REAL :: dudp(udims%i_start:udims%i_end,                               &
                 vdims%j_start:vdims%j_end)

REAL :: dvdp(udims%i_start:udims%i_end,                               &
                 vdims%j_start:vdims%j_end)

REAL :: dwdz(udims%i_start:udims%i_end,                               &
                 vdims%j_start:vdims%j_end)

REAL :: dudx(udims%i_start:udims%i_end,                               &
                 vdims%j_start:vdims%j_end)

REAL :: dvdx(udims%i_start:udims%i_end,                               &
                 vdims%j_start:vdims%j_end)

REAL :: dudy(udims%i_start:udims%i_end,                               &
                 vdims%j_start:vdims%j_end)

REAL :: dvdy(udims%i_start:udims%i_end,                               &
                 vdims%j_start:vdims%j_end)

REAL :: catpred(udims%i_start:udims%i_end,                            &
                 vdims%j_start:vdims%j_end)

INTEGER :: icode

CHARACTER (LEN=*), PARAMETER :: RoutineName='PWS_CAT'
! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! -----------------------------------

icode=0

!interpolate input fields to B grid versions.
  CALL  uC_to_uB(                                                   &
       u(udims%i_start:udims%i_end,                                 &
         udims%j_start:udims%j_end,                                 &
         udims%k_start:udims%k_end),                                &
  udims%i_len,udims%j_len,vdims%j_len,udims%k_len,                  &
  udims_s%halo_i,udims_s%halo_j,                                    &
  local_ub)

  CALL  vC_to_vB(                                                   &
       v(vdims%i_start:vdims%i_end,                                 &
         vdims%j_start:vdims%j_end,                                 &
         vdims%k_start:vdims%k_end),                                &
  pdims%j_len,pdims%i_len,vdims%j_len,vdims%k_len,                  &
  vdims_s%halo_i,vdims_s%halo_j,                                    &
  global_row_length,local_vb,local_ub)

  CALL  pc_to_pb(                                                   &
         r_rho_levels(1:row_length,                                 &
                        1:rows,                                     &
                        pdims%k_start:pdims%k_end),                 &
  pdims%i_len,pdims%j_len,vdims%j_len,pdims%k_len,                  &
  pdims_s%halo_i,pdims_s%halo_j,                                    &
  local_r_rho_levels_b)

  CALL  pc_to_pb(                                                   &
         p(pdims%i_start:pdims%i_end,                               &
         pdims%j_start:pdims%j_end,                                 &
         pdims%k_start:pdims%k_end),                                &
  pdims%i_len,pdims%j_len,vdims%j_len,pdims%k_len,                  &
  pdims_s%halo_i,pdims_s%halo_j,                                    &
  local_pb)


!---- VERTICAL DERIVATIVES ----
! Find vertical derivatives of z, u & v wrt pressure using cubic spline
! Note that spline values must be monotonically increasing.

 CALL DiffP(pdims%k_len, press_lev,                    &
            udims%i_start,udims%i_end,                 &
            vdims%j_start,vdims%j_end,                 &
            local_r_rho_levels_b, local_pb, dZdP)

 CALL DiffP(udims%k_len, press_lev,                    &
            udims%i_start,udims%i_end,                 &
            vdims%j_start,vdims%j_end,                 &
            local_ub, local_pb, dUdP)

 CALL DiffP(vdims%k_len, press_lev,                    &
            udims%i_start,udims%i_end,                 &
            vdims%j_start,vdims%j_end,                 &
            local_vb, local_pb, dVdP)

dWdZ(:,:) = SQRT(dUdP(:,:)**2 + dVdP(:,:)**2) &
                                   / ABS(dZdP(:,:))

! calculate the horizontal derivs of U and V wrt to long and lat.

 CALL DiffX( u_press, dUdX, fld_type_u, icode )

 CALL DiffX( v_press, dVdX, fld_type_v, icode)

 CALL DiffY( u_press, dUdY, fld_type_u, icode )

 CALL DiffY( v_press, dVdY, fld_type_v, icode )

! calculate cat predictor using Dutton.

 CALL DuttonCAT(u_press, v_press,             &
                dUdX, dUdY, dVdX, dVdY, dWdZ, &
                catpred, icode )

pws_cat_turb(:,:,k) = catpred (:,:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE pws_cat

! ----------------------------------------

SUBROUTINE DuttonCAT( UStdLev,     &  ! in 1,2
                      VStdLev,     &  ! in
                      dUdX,        &  ! in
                      dUdY,        &  ! in
                      dVdX,        &  ! in
                      dVdY,        &  ! in
                      dWdZ,        &  ! in
                      CATPred,     &  ! inout
                      ErrorStatus )   ! inout

! Description:
!   Routine to convert values of horizontal and vertical wind shear into
!   a CAT predictor between 0 and 7.5
!
! Method:
!   Calculate HWS and VWS from derivatives of wind components.  Use
!   Dutton's empirical formula to combine HWS and VWS (see Dutton, 1980:
!   Probability forecasts of clear air turbulence based on numerical
!   output: Met Mag v109 pp293-310).  Interpolate to values between 0
!   and 7.5.
!
! IMPORTANT NOTE:
!   The interpolation knots and Dutton's formula are empirically
!   derived, so adjustments are necessary with significant model changes
!   to maintain consistent output on AvSigWX charts.
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE atm_fields_bounds_mod
USE missing_data_mod, ONLY: rmdi
USE trignometric_mod, ONLY: sin_v_latitude
USE conversions_mod, ONLY: recip_pi_over_180

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine Arguments:


REAL, INTENT(IN) :: UStdLev (udims%i_start:udims%i_end,                   &
                             vdims%j_start:vdims%j_end)   ! U & V wind comps on
REAL, INTENT(IN) :: VStdLev (udims%i_start:udims%i_end,                   &
                             vdims%j_start:vdims%j_end)   !   std level
REAL, INTENT(IN) :: dUdX(udims%i_start:udims%i_end,                   &
                         vdims%j_start:vdims%j_end),     &
                    dUdY(udims%i_start:udims%i_end,                   &
                         vdims%j_start:vdims%j_end)    ! Horiz. derivs of u

REAL, INTENT(IN) :: dVdX(udims%i_start:udims%i_end,                   &
                         vdims%j_start:vdims%j_end),     &
                    dVdY(udims%i_start:udims%i_end,                   &
                         vdims%j_start:vdims%j_end)    ! Horiz. derivs of v

REAL, INTENT(IN) :: dWdZ(udims%i_start:udims%i_end,                   &
                         vdims%j_start:vdims%j_end)    ! Vert. deriv of wspeed

REAL, INTENT(OUT) :: CATPred(udims%i_start:udims%i_end,                   &
                               vdims%j_start:vdims%j_end)    ! CAT diagnostic

INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "DUTTONCAT"
INTEGER,PARAMETER :: knots=3     ! No knots used to interp. from Dutton

! Local Variables:
INTEGER :: i, j                  ! Loop counters
INTEGER :: kk, kref, itest       ! For finding relevant knots for interp
REAL :: alat                     ! Latitude
REAL :: alpha                    ! Interpolation coefficient
REAL :: vmag                     ! Wind speed
REAL :: CATProb = rmdi           ! Final CAT predictor
REAL :: empire                   ! Dutton's empirical CAT predictor
REAL :: ShearH                   ! Horizontal wind shear (HWS)
REAL :: ShearV                   ! Vertical wind shear   (VWS)

REAL :: eval(3)   = (/5.0, 7.5 , 65.0/) ! Values of Dutton's
                                        ! empirical indicator.
REAL :: catval(3) = (/0.0, 1.75,  7.5/) ! Corresponding values of
                                        ! final CAT predictor
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header --------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO j = vdims%j_start,vdims%j_end
  !--------------------------------------- calculate latitude
  alat = ASIN(sin_v_latitude(1,j)) * recip_pi_over_180

  DO i = udims%i_start,udims%i_end

    IF ( (dUdX(i,j) == rmdi) .OR.  &
         (dUdY(i,j) == rmdi) .OR.  &
         (dVdX(i,j) == rmdi) .OR.  &
         (dVdY(i,j) == rmdi) ) THEN
      CATPred(i,j) = rmdi
    ELSE
      !  Requires ShearH in units of (m/s)/(100km),ShearV in units of
      !  (m/s)/km.  SI units are (m/s)/m, so need to take
      !  100000.*ShearH and 1000.*ShearV.
      ShearV = 1000.0 * dWdZ(i,j)
      vmag  = UStdLev(i,j)**2 + VStdLev(i,j)**2
      IF ( vmag > 0.0 ) THEN
        ShearH = (100000.0/vmag) * &
         (  UStdLev(i,j) * VStdLev(i,j) * dUdX(i,j) &
          - UStdLev(i,j) * UStdLev(i,j) * dUdY(i,j) &
          + VStdLev(i,j) * VStdLev(i,j) * dVdX(i,j) &
          - UStdLev(i,j) * VStdLev(i,j) * dVdY(i,j) )
      ELSE
        ShearH = 0.0
      END IF
      ! vector form for shearH implies a change of sign if in S.Hemi
      IF ( alat < 0.0 ) THEN
        ShearH = -ShearH
      END IF

      empire = 10.5+1.5*(ShearH+ShearV)
      CATProb = catval(knots)
      itest = 0
      ! kref must be initialised here for vectorisation
      kref=0

      IF ( empire < eval(1) ) THEN
        CATProb = catval(1)
        itest = 2
      END IF
      !      DO kk = 2,knots
      !        IF ( (empire < eval(kk)) .AND.  &
      !             (itest == 0) ) THEN
      !          kref  = kk
      !          itest = 1
      !        END IF
      !      END DO
      ! the loop is manually unrolled for vectorisation, knowing that knots=3
      !
      IF (( empire < eval(2)) .AND. ( itest==0)) THEN
        kref=2
        itest=1
      END IF

      IF (( empire < eval(3)) .AND. ( itest==0)) THEN
        kref=3
        itest=1
      END IF

      IF ( itest == 1 ) THEN
        alpha   = (empire-eval(kref-1))/(eval(kref)-eval(kref-1))
        CATProb = alpha*catval(kref) + (1.0-alpha)*catval(kref-1)
      END IF

      CATPred(i,j) = CATProb
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE DuttonCAT

END MODULE pws_cat_mod
