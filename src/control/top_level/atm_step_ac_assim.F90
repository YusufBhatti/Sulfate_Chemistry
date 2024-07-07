! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Contains atm_step code for AC assimilation
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Subroutine Interface:
SUBROUTINE atm_step_ac_assim( &
obs_flag_len, obs_len, obs_flag, obs,                 &
q_star, qcl_star, qcf_star, theta_star, ntml, cumulus,&
errorstatus)

USE submodel_mod, ONLY: atmos_sm, atmos_im
USE stash_array_mod, ONLY: stash_maxlen
USE atm_step_local
USE atm_fields_bounds_mod, ONLY : pdims, tdims, tdims_s

USE level_heights_mod
USE trignometric_mod
USE jules_snow_mod, ONLY: nsmax

USE dyn_var_res_Mod, ONLY: glambda_p, phi_p

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE Control_Max_Sizes
USE lbc_mod
USE ac_ctl_mod, ONLY: ac_ctl

USE acp_namel_mod,  ONLY: l_ac

USE nlstcall_mod, ONLY: lassimilation

USE free_tracers_inputs_mod, ONLY: a_tracer_last, a_tracer_first
USE jules_sea_seaice_mod, ONLY: nice, nice_use
USE jules_surface_types_mod, ONLY: ntype, npft
USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, land_field, len1_lookup,                &
    len_dumphist, len_fixhd, len_tot, model_levels,                    &
    mpp_len1_lookup, n_cca_lev, n_obj_d1_max, n_rows,                  &
    ntiles, river_row_length, river_rows, row_length, rows, sm_levels, &
    theta_field_size, theta_off_size, tr_levels, tr_ukca, tr_vars

USE atm_fields_mod, ONLY: p, p_theta_levels, exner_theta_levels, &
    cf_area, cf_bulk, cf_liquid, cf_frozen, pstar

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


! Subroutine arguments

INTEGER :: obs_flag_len
INTEGER :: obs_len
INTEGER :: obs_flag(obs_flag_len)
REAL    :: obs(obs_len)

REAL :: theta_star(tdims%i_start:tdims%i_end,            &
                   tdims%j_start:tdims%j_end,            &
                   tdims%k_start:tdims%k_end)

REAL :: q_star    (tdims%i_start:tdims%i_end,            &
                   tdims%j_start:tdims%j_end,            &
                   tdims%k_start:tdims%k_end)
REAL :: qcl_star  (tdims%i_start:tdims%i_end,            &
                   tdims%j_start:tdims%j_end,            &
                   tdims%k_start:tdims%k_end)
REAL :: qcf_star  (tdims%i_start:tdims%i_end,            &
                   tdims%j_start:tdims%j_end,            &
                   tdims%k_start:tdims%k_end)

INTEGER :: ntml   (pdims%i_start:pdims%i_end, pdims%j_start:pdims%j_end)
LOGICAL :: cumulus(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)

INTEGER :: errorstatus

! Local variables

CHARACTER(LEN=errormessagelength) :: cmessage  ! Error message if return code >0
CHARACTER(LEN=*) :: RoutineName
PARAMETER (RoutineName='ATM_STEP_AC_ASSIM')
LOGICAL, PARAMETER :: l_mr_acctl    = .FALSE. ! Use mixing ratios for ac_ctl

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (l_ac .AND. lassimilation .AND. errorstatus == 0) THEN
  ! DO AC assimilation

  ! copy non halo values of _star moisture variables into WORK arrays
  ! for passing to AC_CTL where they can be interpreted more
  ! conveniently as 2-d arrays (and halos not required)
  ALLOCATE ( stashwork18(stash_maxlen(18,atmos_im)) )
  ALLOCATE ( work_q  (tdims%i_start:tdims%i_end,        &
                      tdims%j_start:tdims%j_end,        &
                                  1:tdims%k_end) )
  ALLOCATE ( work_qcl(tdims%i_start:tdims%i_end,        &
                      tdims%j_start:tdims%j_end,        &
                                  1:tdims%k_end) )
  ALLOCATE ( work_qcf(tdims%i_start:tdims%i_end,        &
                      tdims%j_start:tdims%j_end,        &
                                  1:tdims%k_end) )

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)           &
!$OMP SHARED(tdims, work_q, q_star, work_qcl, qcl_star, work_qcf, qcf_star)
  DO k =            1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        work_q(i,j,k)   = q_star(i,j,k)
        work_qcl(i,j,k) = qcl_star(i,j,k)
        work_qcf(i,j,k) = qcf_star(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO


  CALL ac_ctl(stash_maxlen(18,atmos_im),theta_field_size,               &
                        theta_star(:,:,1:), work_q, work_qcl, work_qcf, &
                        obs_flag,obs,obs_flag_len,obs_len,              &
                        p, p_theta_levels, exner_theta_levels,          &
                        cos_theta_latitude,                             &
                        cf_area, cf_bulk, cf_liquid, cf_frozen,         &
                        pstar,                            &
                        ntml,                            &
                        cumulus,                           &
                        stashwork18, glambda_p(lambda_start), phi_p,    &
                        l_mr_acctl,errorstatus,cmessage)


  ! moved from out of AC_CTL
  ! DEPENDS ON: stash
  CALL stash(atmos_sm, atmos_im, 18, stashwork18,                  &
                          errorstatus,cmessage)

  ! Check error condition
  IF (errorstatus >  0) THEN

    CALL ereport(RoutineName,errorstatus,cmessage)
  END IF

  ! copy back new values of _star moisture variables
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)           &
!$OMP SHARED(tdims, work_q, q_star, work_qcl, qcl_star, work_qcf, qcf_star)
  DO k =              1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        q_star  (i,j,k) = work_q(i,j,k)
        qcl_star(i,j,k) = work_qcl(i,j,k)
        qcf_star(i,j,k) = work_qcf(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  DEALLOCATE (work_q)
  DEALLOCATE (work_qcl)
  DEALLOCATE (work_qcf)
  DEALLOCATE (stashwork18)

END IF      ! lassimilation and L_AC

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE atm_step_ac_assim


