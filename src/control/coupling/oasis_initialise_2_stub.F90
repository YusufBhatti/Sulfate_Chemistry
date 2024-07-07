#if !defined(MCT)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine: oasis_initialise_2
!
!  Purpose: Second part of initialisation for OASIS3-MCT (dummy)
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Coupling
!  -------------------------------------------------------------------

SUBROUTINE oasis_initialise_2 (cmessage)

USE atm_fields_mod
USE atm_fields_bounds_Mod

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

! Declarations:
USE io
USE um_types
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE Control_Max_Sizes
USE Decomp_DB
! Make sure um_sleep is picked up from portio.
USE io_dependencies
USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, len1_lookup,                            &
    len_dumphist, len_fixhd, len_tot, model_levels, mpp_len1_lookup,   &
    n_cca_lev, n_obj_d1_max, sm_levels, st_levels, tpps_ozone_levels,  &
    tr_lbc_ukca, tr_lbc_vars, tr_ukca, tr_vars

USE errormessagelength_mod, ONLY: errormessagelength
USE umprintmgr, ONLY: umPrint, umMessage

IMPLICIT NONE


!
! ----------------------------------------------------------------------
!
!

! Subroutine arguments:
CHARACTER(LEN=errormessagelength) :: cmessage    ! OUT - Error return message

! Local variables:

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'OASIS_INITIALISE_2'
CHARACTER(LEN=errormessagelength)   :: Message

INTEGER :: ErrorStat

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

ErrorStat   = 1
Message     = 'OASIS3-MCT Routines unavailable - see output.'

WRITE(umMessage,'(A)')     &
      '**ERROR**: oasis_initialise_2 called but is unavailable.'
CALL umPrint(umMessage,src='oasis_initialise_2_stub')
WRITE(umMessage,'(A)') 'Check MCT cpp key is set'
CALL umPrint(umMessage,src='oasis_initialise_2_stub')

CALL ereport(RoutineName, ErrorStat, Message)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE oasis_initialise_2
#endif
