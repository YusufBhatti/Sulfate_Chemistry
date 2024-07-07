! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE calc_3d_cca_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CALC_3D_CCA_MOD'
CONTAINS

!  Subroutine CALC_3D_CCA: Calculates a conv. cld amt on theta model levels.
!
!  Subroutine Interface:

SUBROUTINE calc_3d_cca                                                        &
  ( np_field, npnts, nlev, n_cca_lev, nbl, cld_base, cld_top, p_lyr_bnds      &
  , frz_lev, cca_2d, cca_3d, z_theta, z_rho )

USE cv_param_mod, ONLY:                                                     &
    deep_dp, anv_pressure, anv_height, anv_model_levels,                    &
    anv_limited_pressure_depth

USE cv_run_mod, ONLY:                                                       &
    l_cloud_deep, tower_factor, anvil_factor, anv_opt, l_ccrad

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE
!
! Description:
! ------------
! Calculates a 3D convective cloud amount (i.e. on theta model levels) from
! the 2D convective cloud amount array according to parameters specified in
! the gui/namelist and the position of cloud base, cloud top and freezing level.
!
! Method:
! -------
! The 2D convective cloud amount is expanded into the vertical by applying it
! between the cloud base and top with the additional constraints that:
!
!         (ia)  If the cloud base is in the boundary layer (Original)
!         (ib)  If the cloud base is below the freezing level (CCRad)
!         (ii)  Cloud top is above the freezing level and
!         (iii) The cloud is more than 500mb deep
!
! Then the cloud below the freezing level will be multiplied by TOWER_FACTOR,
! and the cloud above the freezing level will be linearly
! (model level/height/pressure(default)) increased to cloud top where it will
! be equal to the 2D fraction * ANVIL_FACTOR.
!
! NOTE: ***** The above method needs to be rewritten if these mods are********
!       ***** Implemented*****************************************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!

!-----------------------------------------------------------------------------
!   Scalar arguments with intent(in):
!-----------------------------------------------------------------------------

INTEGER, INTENT(IN)  :: npnts        ! Number of points
INTEGER, INTENT(IN)  :: np_field     ! Full data length
INTEGER, INTENT(IN)  :: nlev         ! Number of levels
INTEGER, INTENT(IN)  :: n_cca_lev    ! Number of cca levels
INTEGER, INTENT(IN)  :: nbl          ! Number of Boundary layer levels

!-----------------------------------------------------------------------------
!   Array  arguments with intent(in):
!-----------------------------------------------------------------------------

REAL,    INTENT(IN)  :: z_theta    (np_field,   nlev) ! z (th layer centres)
REAL,    INTENT(IN)  :: z_rho      (np_field,   nlev) ! z (rh level  bounds)
REAL,    INTENT(IN)  :: p_lyr_bnds (np_field, 0:nlev)
                                                ! Pressure on layer
                                                ! boundaries (rho levels-1)

REAL,    INTENT(IN)  :: cca_2d   (npnts) ! 2D convective cloud amount
INTEGER, INTENT(IN)  :: cld_top  (npnts) ! Conv. cloud top  (theta level)
INTEGER, INTENT(IN)  :: cld_base (npnts) ! Conv. cloud base (theta level)
INTEGER, INTENT(IN)  :: frz_lev  (npnts) ! Freezing level   (theta level)




!-----------------------------------------------------------------------------
!   Array  arguments with intent(out):
!-----------------------------------------------------------------------------

REAL,    INTENT(OUT) :: cca_3d(np_field, n_cca_lev)
                                         ! Convective cloud amount on
                                         ! model levels (theta levels)


! Local variables:
! -----------------
INTEGER             :: i, k             ! Loop counters


INTEGER :: anv_lev     ! Base level of 'anvil'
REAL    :: anv_dep     ! Anvil depth in model levels
REAL    :: anv_p_dep   ! Anvil depth in pressure
REAL    :: anv_z_dep   ! Anvil depth in metres
REAL    :: anv_p_base  ! Anvil base pressure (rho-level)

REAL    :: p_cld_base  ! Pressure at lowest  cloud layer BOUNDARY
REAL    :: p_cld_top   ! Pressure at highest cloud layer BOUNDARY



LOGICAL :: cbct_crit (npnts) ! .TRUE. if cloud base/top are sensible
LOGICAL :: dep_crit  (npnts) ! .TRUE. if depth criteria met
LOGICAL :: base_crit (npnts) ! .TRUE. if cloud base criteria met
LOGICAL :: anv_on    (npnts) ! .TRUE. if all anvil criteria met
INTEGER :: tp_of_lp  (npnts) ! Index of cloud top, required so that CCRad
                             ! correction can be reverted

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_3D_CCA'

!-----------------------------------------------------------------------------
! Code Statements
!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!---------------------------------------------------------------------------
! 1.0 Initialise local arrays
!---------------------------------------------------------------------------

DO i=1, npnts
  cbct_crit(i)  = .FALSE.
  dep_crit(i)   = .FALSE.
  base_crit(i)  = .FALSE.
  anv_on(i)     = .FALSE.
  tp_of_lp(i)   = 0
END DO



!---------------------------------------------------------------------------
! 2.0 Option changes/bugfixes specified by CCRad
!---------------------------------------------------------------------------
IF (l_ccrad) THEN

  DO i=1, npnts
    ! Check for sensible Cloud base/top
    cbct_crit(i)  = ((cld_top(i)  >= cld_base(i)) .AND.                     &
                     (cld_top(i)  /= 0)           .AND.                     &
                     (cld_base(i) /= 0))

    ! Check for cloud base/top above/below freezing level
    base_crit(i)  = ((cld_base(i) < frz_lev(i))   .AND.                     &
                     (cld_top(i)  > frz_lev(i)))

    ! Change index of cloud top (Bug fix)
    tp_of_lp(i)   = cld_top(i)
  END DO

ELSE
  ! Original test criteria
  DO i=1, npnts
    cbct_crit(i) = .TRUE.
    base_crit(i) = ((cld_base(i) < nbl)           .AND.                     &
                    (cld_top(i)  > frz_lev(i)))
    tp_of_lp(i)  = cld_top(i)-1
  END DO

END IF


! Locate grid points where depth criteria is satisfied
IF (l_cloud_deep) THEN
  DO i=1, npnts
    p_cld_base  = p_lyr_bnds(i,MAX(cld_base(i)-1,0))
    p_cld_top   = p_lyr_bnds(i,cld_top(i))
    dep_crit(i) = (p_cld_base - p_cld_top) >= deep_dp
  END DO
ELSE
  DO i=1,npnts
    dep_crit(i) = .TRUE.
  END DO
END IF



!---------------------------------------------------------------------------
! 3.0 Locate grid points where all anvil criteria are satisfied
!---------------------------------------------------------------------------
DO i=1, npnts
  IF (base_crit(i) .AND. dep_crit(i)) THEN
    anv_on(i) = .TRUE.
  END IF
END DO


!---------------------------------------------------------------------------
! 4.0 Apply CCA Profiles
!---------------------------------------------------------------------------
DO i=1, npnts
  IF ( cbct_crit(i)       .AND.                                             &
      (cca_2d(i) > 0.0) ) THEN

    IF (anv_on(i)) THEN

      !---------------------------------------------------------------------
      ! 4.1a Cloud satisfies anvil criteria: Apply Anvil
      !---------------------------------------------------------------------
      SELECT CASE(anv_opt)
      CASE (anv_height)
        !-----------------------------------------------------------------
        ! CCA increases with height from freezing level to cloud-top
        !-----------------------------------------------------------------
        anv_lev   = MAX(cld_base(i), frz_lev(i))
        anv_z_dep = z_rho(i,cld_top(i)) - z_rho(i,anv_lev-1)

        DO k=anv_lev, cld_top(i)

          cca_3d(i,k) =                                                   &
                  (anvil_factor - tower_factor)*cca_2d(i)                 &
                * (z_theta(i,k) - z_rho(i,anv_lev-1)) / anv_z_dep         &
                + (cca_2d(i) * tower_factor)

          IF (cca_3d(i,k) >= 1.0) THEN
            cca_3d(i,k) = 0.99
          END IF

        END DO


      CASE (anv_model_levels)
        !-----------------------------------------------------------------
        ! CCA increases with model level from freezing level to cloud-top:
        ! (original code)
        !-----------------------------------------------------------------
        anv_lev = MAX(cld_base(i), frz_lev(i))
        anv_dep = cld_top(i) - anv_lev

        DO k=anv_lev, tp_of_lp(i)

          cca_3d(i,k) =                                                   &
                  (anvil_factor - tower_factor) * cca_2d(i)               &
                * (k - anv_lev + 1) / anv_dep                             &
                + (cca_2d(i) * tower_factor)

          IF (cca_3d(i,k)  >=  1.0) THEN
            cca_3d(i,k) = 0.99
          END IF

        END DO


      CASE (anv_pressure)
        !-----------------------------------------------------------------
        ! CCA increases with pressure from freezing level to cloud-top
        !-----------------------------------------------------------------
        anv_lev   = MAX(cld_base(i), frz_lev(i))
        anv_p_dep = p_lyr_bnds(i,anv_lev-1) - p_lyr_bnds(i,cld_top(i))

        DO k=anv_lev, tp_of_lp(i)

          cca_3d(i,k) =                                                   &
                (anvil_factor - tower_factor)*cca_2d(i)                   &
              * (p_lyr_bnds(i,anv_lev-1) -  p_lyr_bnds(i,k)) / anv_p_dep  &
              + (cca_2d(i) * tower_factor)

          IF (cca_3d(i,k) >= 1.0) THEN
            cca_3d(i,k) = 0.99
          END IF

        END DO

      CASE (anv_limited_pressure_depth)
        !-----------------------------------------------------------------
        ! Pressure based, but limit anvil depth to 5000.0 pa
        !-----------------------------------------------------------------
        anv_p_base = p_lyr_bnds(i,cld_top(i)) + 5000.0

        ! Ensure that Anvil is at least 2 levels deep, so only loop to
        ! layer below cloud top so that it will be 2 levels deep even
        ! if k is set at top of loop
        DO k=1, tp_of_lp(i)-1
          IF (p_lyr_bnds(i,k) > anv_p_base) THEN
            anv_lev = k
          END IF
        END DO

        anv_p_dep = p_lyr_bnds(i,anv_lev-1) - p_lyr_bnds(i,cld_top(i))

        DO k=anv_lev, tp_of_lp(i)

          cca_3d(i,k) =                                                   &
                (anvil_factor - tower_factor)*cca_2d(i)                   &
              * (p_lyr_bnds(i,anv_lev-1) -  p_lyr_bnds(i,k)) / anv_p_dep  &
              + (cca_2d(i) * tower_factor)

          IF (cca_3d(i,k) >= 1.0) THEN
            cca_3d(i,k) = 0.99
          END IF

        END DO

      END SELECT

      !---------------------------------------------------------------------
      ! 4.1b Cloud satisfies anvil criteria: Apply Tower below anvil base
      !---------------------------------------------------------------------
      DO k=cld_base(i), anv_lev-1
        cca_3d(i,k) = tower_factor * cca_2d(i)
      END DO

    ELSE

      ! Anvil criteria not met
      DO k=cld_base(i), tp_of_lp(i)
        cca_3d(i,k) = cca_2d(i)
      END DO

    END IF      ! End test on anvil criteria

    !-----------------------------------------------------------------------
    ! Finally check there is no cloud below cloud base or above cloud top:
    ! (original code)
    !-----------------------------------------------------------------------
    IF (.NOT. l_ccrad) THEN
      DO k=1, (cld_base(i)-1)
        cca_3d(i,k) = 0.0
      END DO

      IF (cld_top(i) > 0) THEN
        DO k=cld_top(i), n_cca_lev
          cca_3d(i,k) = 0.0
        END DO
      END IF

    END IF      ! l_ccrad
  END IF      ! cca_2d > 0 and sensible ccb/cct
END DO      ! loop over npnts


!
!=============================================================================
!  End of anvil calculation
!=============================================================================
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE calc_3d_cca
END MODULE calc_3d_cca_mod
