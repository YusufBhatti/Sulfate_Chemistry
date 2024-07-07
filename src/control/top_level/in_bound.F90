! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine IN_BOUND
!
!   Purpose : Takes as input,the code defining whether updates of
!    boundary data are required.  The physical files required are
!    identified, and the headers lookup tables are read into common
!    blocks.  Reads the update intervals from the boundary datasets.
!    Where the update interval is in months or years, the check will be
!    made daily.
!
!   Programming standard; Unified Model Documentation Paper No. 3
!   version no. 1, dated 15/01/90
!
!   Logical components covered : C720
!
!   System task : C7
!
!   Documentation : Unified Model Documentation Paper No C7
!
!
!    Arguments
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

SUBROUTINE in_bound (                                             &
   a_len1_levdepcda,a_len2_levdepcda,                             &
   a_len1_rowdepcda,a_len2_rowdepcda,                             &
   a_len1_coldepcda,a_len2_coldepcda,                             &
           icode,cmessage)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE atm_fields_bounds_mod
USE UM_ParVars
USE Control_Max_Sizes
USE lbc_mod

USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, len1_lbc_comp_lookup, len1_lookup,      &
    len_dumphist, len_fixhd, model_levels,                             &
    mpp_len1_lookup, n_cca_lev, sm_levels, st_levels,                  &
    tpps_ozone_levels, tr_lbc_ukca, tr_lbc_vars, tr_ukca, tr_vars

USE errormessagelength_mod, ONLY: errormessagelength
USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE


INTEGER ::                                                        &
         a_len1_levdepcda,                                        &
                             ! IN : copy of A_LEN1_LEVDEPC
         a_len2_levdepcda,                                        &
                             ! IN : copy of A_LEN2_LEVDEPC
         a_len1_rowdepcda,                                        &
                             ! IN : copy of A_LEN1_ROWDEPC
         a_len2_rowdepcda,                                        &
                             ! IN : copy of A_LEN2_ROWDEPC
         a_len1_coldepcda,                                        &
                             ! IN : copy of A_LEN1_COLDEPC
         a_len2_coldepcda,                                        &
                             ! IN : copy of A_LEN2_COLDEPC
        icode            ! Return code = 0 Normal Exit
!                              !    "     "  > 0 Error Exit

CHARACTER(LEN=errormessagelength) :: cmessage  ! Error message if ICODE > 0

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IN_BOUND'
!

!       Internal Structure

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

IF (model_type /= mt_global) THEN

  ! Call inbound for the atmosphere

  ! DEPENDS ON: inbounda
  CALL inbounda(                                                    &
     a_len1_levdepcda,a_len2_levdepcda,                             &
     a_len1_rowdepcda,a_len2_rowdepcda,                             &
     a_len1_coldepcda,a_len2_coldepcda)

END IF

!   4   End of routine

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE in_bound
