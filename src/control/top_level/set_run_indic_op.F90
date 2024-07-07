! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Subroutine: SET_RUN_INDIC_OP-----------------------------------
!
!    Purpose: Interface routine which formerly set RUN_INDIC_OP in dump
!             header. Called by INITIAL.
!             Under Rose itab is used directly from model_id_mod
!
!    Programming standard: UM Doc Paper 3
!
!    Logical components covered: C0
!
!    Project task: C0
!
!    -------------------------------------------------------------------
!    Interface and arguments: ------------------------------------------
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Top Level

SUBROUTINE set_run_indic_op(                                      &
              icode,cmessage)
!
! ----------------------------------------------------------------------

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams

USE model_id_mod, ONLY: itab
USE missing_data_mod, ONLY: imdi
USE dump_headers_mod, ONLY: a_fixhd
USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, len1_lookup,                            &
    len_dumphist, len_fixhd, mpp_len1_lookup
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!
!  Arguments
!
INTEGER :: icode             ! Work - Internal return code
CHARACTER(LEN=errormessagelength) :: cmessage    ! Work - Internal error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_RUN_INDIC_OP'
!
!  Local variables: none
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!  Under Rose itab becomes the default for all model runs 
!  and run_indic_op is retired.

IF (itab /= imdi) a_fixhd(6) = itab
!
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE set_run_indic_op
