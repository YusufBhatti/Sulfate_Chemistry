#if !defined(MCT)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis_updatecpl(cmessage)

USE atm_fields_bounds_mod

USE jules_snow_mod, ONLY: nsmax
USE jules_sea_seaice_mod, ONLY: nice, nice_use

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE Control_Max_Sizes
USE Decomp_DB
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE free_tracers_inputs_mod, ONLY: a_tracer_last, a_tracer_first
USE nlsizes_namelist_mod, ONLY:                                         &
    land_field, len_tot, model_levels, n_cca_lev, n_obj_d1_max, n_rows, &
    ntiles, river_row_length, river_rows, row_length,                   &
    rows, sm_levels, st_levels, theta_off_size, tpps_ozone_levels,      &
    tr_lbc_ukca, tr_lbc_vars, tr_levels, tr_ukca, tr_vars

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


!
!
! Description:
! Stub routine for oasis_updatecpl to allow compilation in
! stand-alone atmosphere jobs. Run time logic should
! mean this version never actually gets called.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
!
!=====================================================================

CHARACTER(LEN=errormessagelength) :: cmessage ! OUT - Error return message

CHARACTER(LEN=errormessagelength)   :: Message
CHARACTER(LEN=*), PARAMETER   :: RoutineName = 'OASIS_UPDATECPL'
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

CALL umPrint( '**ERROR**: oasis_updatecpl is unavailable.', &
    src='oasis_updatecpl_stub')
CALL umPrint( 'Check MCT cpp key is set',src='oasis_updatecpl_stub')


CALL ereport(RoutineName, ErrorStat, Message)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

END SUBROUTINE oasis_updatecpl
#endif
