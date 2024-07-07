! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  --------------- SUBROUTINE IN_ACCTL --------------------------------
!
!   Purpose : Calls the initialisation program for the data assimilation
!   assimilation section
!
!   Programming standard; U M Documentation Paper No. 3
!
!----------------------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation

SUBROUTINE in_acctl(                                              &
                  icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE dyn_var_res_Mod, ONLY:  phi_p,glambda_p

USE atm_fields_bounds_mod

USE UM_ParParams
USE Control_Max_Sizes
USE ac_init_mod, ONLY: ac_init
USE acp_namel_mod, ONLY: l_ac
USE submodel_mod, ONLY: atmos_im
USE nlstcall_mod, ONLY: model_basis_time
USE dump_headers_mod, ONLY: a_realhd, a_levdepc
USE nlsizes_namelist_mod, ONLY:                                          &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,      &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,       &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst,   &
    a_len_inthd, a_len_realhd, bl_levels, len1_lookup,                   &
    len_dumphist, len_fixhd, model_levels, mpp_len1_lookup,              &
    n_cca_lev, n_rows, row_length, rows, sm_levels, st_levels,           &
    tpps_ozone_levels, tr_lbc_ukca, tr_lbc_vars, tr_levels, tr_ukca,     &
    tr_vars

USE model_time_mod, ONLY: &
    secs_per_stepim

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


INTEGER  ::     icode              ! OUT: Error return code
CHARACTER(LEN=errormessagelength) :: cmessage     ! OUT: Error return message


INTEGER :: no_lambda_p

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IN_ACCTL'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

IF (l_ac) THEN

  no_lambda_p=SIZE(glambda_p)
  no_lambda_p=MAX(no_lambda_p,1)

  CALL ac_init (bl_levels, tr_levels,                               &
                rows, n_rows, row_length,                           &
                secs_per_stepim(atmos_im),                          &
                model_basis_time(1), model_basis_time(2),           &
                model_basis_time(3), model_basis_time(4),           &
                model_basis_time(5),                                &
                a_realhd(1), a_realhd(2), a_realhd(3),              &
                a_realhd(4), a_realhd(5), a_realhd(6),              &
                a_levdepc(1), a_levdepc(1),                         &
                glambda_p(1),phi_p,no_lambda_p,                     &
                icode,cmessage)

END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE in_acctl

