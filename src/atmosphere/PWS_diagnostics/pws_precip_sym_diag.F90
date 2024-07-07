! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE pws_precip_sym_diag_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_PRECIP_SYM_DIAG_MOD'

CONTAINS
!
!   Subroutine pws_precip_sym_diag ---------------------------------------
!
!   Purpose: Calculates precip symbol for PWS, stash section 20
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
!   Documentation : Unified Model Documentation Paper No ??
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics

SUBROUTINE pws_precip_sym_diag(STASHwork,                &
                         icode,cmessage)


USE atm_fields_bounds_mod
USE atm_fields_mod

USE jules_snow_mod, ONLY: nsmax
USE jules_sea_seaice_mod, ONLY: nice, nice_use

USE vert_interp_mdi_mod, ONLY: vert_interp_mdi

USE um_stashcode_mod, ONLY: stashcode_pws_sec, stashcode_pws_precip_sym
USE pws_diags_mod,    ONLY:                                       &
    pws_precip_sym_ls_rain,pws_precip_sym_ls_snow,                &
    pws_precip_sym_conv_rain,pws_precip_sym_conv_snow,            &
    pws_precip_sym_t1p5m, pws_precip_sym

USE nlsizes_namelist_mod
USE lbc_mod
USE um_parvars
USE interpor_mod, ONLY: interp_order_linear
USE free_tracers_inputs_mod, ONLY: a_tracer_last, a_tracer_first

USE conversions_mod, ONLY: zerodegc

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim


IMPLICIT NONE


! Local variables


CHARACTER (LEN=*), PARAMETER :: RoutineName = 'PWS_PRECIP_SYM_DIAG'

! S/R arguments
REAL, INTENT(INOUT) :: stashwork(*)        ! STASH workspace
INTEGER :: icode ! Return code : 0 Normal exit, >0 Error exit
CHARACTER(LEN=errormessagelength) :: cmessage  ! Error message if ICODE > 0

REAL   ::  PrecipDyn(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
REAL   ::  PrecipCon(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
REAL   ::  PrecipSnow(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)
REAL   ::  PrecipTot(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)

! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set precip codes
PrecipDyn  = pws_precip_sym_ls_rain + pws_precip_sym_ls_snow
PrecipCon  = pws_precip_sym_conv_rain + pws_precip_sym_conv_snow

PrecipSnow = pws_precip_sym_ls_snow + pws_precip_sym_conv_snow
PrecipTot  = PrecipDyn + PrecipCon

WHERE ( PrecipTot > 0.0 )
  pws_precip_sym = PrecipTot + 1000.0
END WHERE
WHERE ( PrecipSnow > 0.0 )
  pws_precip_sym  = PrecipTot + 2000.0
END WHERE
WHERE ( (pws_precip_sym_conv_snow > 0.0) .AND. &
       ( pws_precip_sym_t1p5m > 2.0+zerodegc) )
  pws_precip_sym  = PrecipTot + 3000.0
END WHERE

WHERE ( (PrecipTot > 0.0) .AND. (PrecipDyn > PrecipCon) )
  pws_precip_sym  = pws_precip_sym + 10000.0
END WHERE
WHERE ( (PrecipTot > 0.0) .AND. (PrecipDyn <= PrecipCon) )
  pws_precip_sym  = pws_precip_sym + 20000.0
END WHERE

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE pws_precip_sym_diag

END MODULE pws_precip_sym_diag_mod
