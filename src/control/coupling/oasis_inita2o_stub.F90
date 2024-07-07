#if !defined(MCT)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
SUBROUTINE oasis_inita2o(cmessage)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE atm_fields_bounds_mod

USE ereport_mod, ONLY: ereport
USE UM_ParParams
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, len1_lookup,                            &
    len_dumphist, len_fixhd, len_tot, model_levels, mpp_len1_lookup,   &
    n_cca_lev, n_obj_d1_max, sm_levels, st_levels, tpps_ozone_levels,  &
    tr_lbc_ukca, tr_lbc_vars, tr_ukca, tr_vars

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!
! Description: Stub routine for oasis_inita2o.
!              Run time controls should mean this routine is never
!              actually called.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
!
!---------------------------------------------------------------------
!
CHARACTER(LEN=errormessagelength) :: cmessage    ! OUT - Error return message
! ----------------------------------------------------------------------
CHARACTER(LEN=errormessagelength)   :: Message
CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'OASIS_INITA2O'
INTEGER             :: ErrorStat        ! Return code:
                                        !   0 = Normal exit
                                        ! +ve = Fatal Error
                                        ! -ve = Warning

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

ErrorStat   = 1
Message     = 'OASIS3-MCT Routines unavailable - see output.'

WRITE(umMessage,'(A)') '**ERROR**: oasis_inita2o unavailable.'
CALL umPrint(umMessage,src='oasis_inita2o_stub')
WRITE(umMessage,'(A)') 'Check MCT cpp key is set'
CALL umPrint(umMessage,src='oasis_inita2o_stub')


CALL ereport(RoutineName, ErrorStat, Message)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
END SUBROUTINE oasis_inita2o
#endif
