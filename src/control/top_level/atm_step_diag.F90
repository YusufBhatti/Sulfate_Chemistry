! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Contains various chunks of code from atm_step - mainly related to
! calculation of diagnostics - the purpose of each section is indicated
! at the head of the section

! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

SUBROUTINE atm_step_diag(flag)

USE atm_step_local
USE atm_fields_bounds_mod, ONLY: tdims

USE jules_sea_seaice_mod, ONLY: nice, nice_use

USE rad_input_mod, ONLY: l_radiation

USE g_wave_input_mod, ONLY: l_gwd, l_use_ussp
USE bl_option_mod, ONLY: i_bl_vn, i_bl_vn_0

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE UM_ParVars
USE Control_Max_Sizes
USE lbc_mod
USE dynamics_testing_mod, ONLY: L_Physics
USE pc2_constants_mod, ONLY: rhcpt_horiz_var, rhcpt_tke_based
USE jules_hydrology_mod,  ONLY: l_hydrology

USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY: stash_maxlen, sf

USE gen_phys_inputs_mod, ONLY: l_use_methox
USE mphys_inputs_mod, ONLY: l_rain
USE electric_inputs_mod, ONLY: l_use_electric
USE cloud_inputs_mod, ONLY: i_rhcpt
USE eng_corr_inputs_mod, ONLY: l_emcorr
USE cv_run_mod, ONLY: l_param_conv  ! Convection scheme switch
USE nlsizes_namelist_mod, ONLY: model_levels, row_length, rows

USE atm_fields_mod, ONLY: q, qcl, exner_theta_levels, theta,              &
    ti, ice_fraction, ice_thickness, ti_cat, ice_fract_cat,               &
    ice_thick_cat, tstar_sice, snodep_sea, tstar_sice_cat, snodep_sea_cat

IMPLICIT NONE

INTEGER :: flag

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ATM_STEP_DIAG'

! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (flag == 1) THEN

  ! Allocate diagnostic space for stash
  IF (l_radiation) THEN
    ALLOCATE (stashwork1(stash_maxlen(1,atmos_im)))
    ALLOCATE (stashwork2(stash_maxlen(2,atmos_im)))
  ELSE
    ALLOCATE (stashwork1(1))
    ALLOCATE (stashwork2(1))
  END IF
  IF (l_rain .OR. l_use_methox) THEN
    ALLOCATE (stashwork4(stash_maxlen(4,atmos_im)))
  ELSE
    ALLOCATE (stashwork4(1))
  END IF
  IF (l_use_electric) THEN
    ALLOCATE (stashwork21(stash_maxlen(21,atmos_im)))
  ELSE
    ALLOCATE (stashwork21(1))
  END IF
  IF (l_gwd .OR. l_use_ussp) THEN
    ALLOCATE (stashwork6(stash_maxlen(6,atmos_im)))
  ELSE
    ALLOCATE (stashwork6(1))
  END IF
  IF (l_emcorr) THEN   ! only if energy correction required
    ALLOCATE (stashwork14(stash_maxlen(14,atmos_im)))
  ELSE
    ALLOCATE (stashwork14(1))
  END IF   ! l_emcorr

  IF (nice  ==  1) THEN
    ! Set sea ice category scheme D1 pointers as categories are not present
    ! if nice = 1.
    p_ti        => ti                    ! used in atmos_physics2
    p_ice_fract => ice_fraction          ! used in atmos_physics2
    p_ice_thick => ice_thickness         ! used in atmos_physics2
  ELSE
    p_ti        => ti_cat
    p_ice_fract => ice_fract_cat
    p_ice_thick => ice_thick_cat
  END IF

  IF (nice_use  ==  1) THEN
    ! Set sea ice category scheme D1 pointers as categories are not fully used
    ! if nice_use = 1.
    p_tstar_sice => tstar_sice            ! used in atmos_physics1&2
    p_snodep_sice => snodep_sea           ! used in atmos_physics1
    p_ice_thick_rad => ice_thickness      ! used in atmos_physics1
    p_ice_fract_rad => ice_fraction       ! used in atmos_physics1
  ELSE
    p_tstar_sice => tstar_sice_cat
    p_snodep_sice => snodep_sea_cat
    p_ice_thick_rad => ice_thick_cat
    p_ice_fract_rad => ice_fract_cat
  END IF

  ! ---------------------------------------------------------


ELSE IF (flag == 4) THEN

  ! Apply diagnostics only at last cycle but need to allocate for all cycles.
  ! Allocate diagnostic space for stash
  IF (.NOT. ALLOCATED(stashwork3)) THEN
    IF ( i_bl_vn /= i_bl_vn_0 ) THEN
      ALLOCATE (stashwork3(stash_maxlen(3,atmos_im)))
    ELSE
      ALLOCATE (stashwork3(1))
    END IF ! i_bl_vn
  END IF

  IF (.NOT. ALLOCATED(stashwork5)) THEN
    IF (l_param_conv) THEN
      ALLOCATE (stashwork5(stash_maxlen(5,atmos_im)))
    ELSE
      ALLOCATE (stashwork5(1))
    END IF   ! on l_param_conv
  END IF

  IF (.NOT. ALLOCATED(stashwork8)) THEN
    IF (l_hydrology) THEN
      ALLOCATE (stashwork8(stash_maxlen(8,atmos_im)))
    ELSE
      ALLOCATE (stashwork8(1))
    END IF
  END IF

  IF (.NOT. ALLOCATED(stashwork9)) THEN
    ALLOCATE (stashwork9(stash_maxlen(9,atmos_im)))
  END IF
  IF (.NOT. ALLOCATED(stashwork19)) THEN
    ALLOCATE (stashwork19(stash_maxlen(19,atmos_im)))
  END IF
  IF (.NOT. ALLOCATED(stashwork26)) THEN
    ALLOCATE (stashwork26(stash_maxlen(26,atmos_im)))
  END IF

  IF (i_rhcpt == rhcpt_horiz_var .OR. i_rhcpt == rhcpt_tke_based) THEN
    ! Dimension diagnostic 3D RHcrit array
    rhc_row_length = row_length
    rhc_rows = rows
  ELSE
    ! RHcrit will be a 1D parametrized array input from user interface
    rhc_row_length = 1
    rhc_rows = 1
  END IF    ! i_rhcpt

  IF ( cycleno == 1 .AND. l_physics ) THEN
    ALLOCATE ( rhcpt(rhc_row_length, rhc_rows, model_levels) )
  END IF


  ! ---------------------------------------------------------


ELSE IF (flag == 8) THEN

  IF (sf(181,15) ) THEN
    ALLOCATE ( t_incr_diagnostic(tdims%i_start:tdims%i_end,        &
                                 tdims%j_start:tdims%j_end,        &
                                             1:tdims%k_end) )
    DO k=            1, tdims%k_end
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          t_incr_diagnostic(i,j,k) = theta(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( t_incr_diagnostic(1,1,1) )
  END IF

  IF (sf(182,15) ) THEN
    ALLOCATE ( q_incr_diagnostic(tdims%i_start:tdims%i_end,        &
                                 tdims%j_start:tdims%j_end,        &
                                             1:tdims%k_end) )
    DO k=            1, tdims%k_end
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          q_incr_diagnostic(i,j,k) = q(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( q_incr_diagnostic(1,1,1) )
  END IF

  IF (sf(183,15) ) THEN
    ALLOCATE ( qcl_incr_diagnostic(tdims%i_start:tdims%i_end,        &
                                   tdims%j_start:tdims%j_end,        &
                                               1:tdims%k_end) )
    DO k=            1, tdims%k_end
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          qcl_incr_diagnostic(i,j,k) = qcl(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( qcl_incr_diagnostic(1,1,1) )
  END IF


  ! ---------------------------------------------------------


ELSE IF (flag == 9) THEN

   ! calculate changes to T , q etc
  IF (sf(181,15) ) THEN
    DO k=            1, tdims%k_end
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          t_incr_diagnostic(i,j,k)= (theta(i,j,k)-             &
                        t_incr_diagnostic(i,j,k)) * exner_theta_levels(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  END IF

  IF (sf(182,15) ) THEN
    DO k=            1, tdims%k_end
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          q_incr_diagnostic(i,j,k) = q(i,j,k) - q_incr_diagnostic(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  END IF

  IF (sf(183,15) ) THEN
    DO k=            1, tdims%k_end
      DO j=tdims%j_start, tdims%j_end
        DO i=tdims%i_start, tdims%i_end
          qcl_incr_diagnostic(i,j,k) = qcl(i,j,k) - qcl_incr_diagnostic(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
  END IF


  ! ---------------------------------------------------------

  ! Allocate stash space for Nudging
ELSE IF (flag == 39 ) THEN
  ALLOCATE ( stashwork39(STASH_maxlen(39,atmos_im)) )
END IF
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE atm_step_diag
