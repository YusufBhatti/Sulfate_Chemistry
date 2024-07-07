! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE pws_thickness_diag_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_THICKNESS_DIAG_MOD'

CONTAINS
!
!   Subroutine pws_thickness_diag ---------------------------------------
!
!   Purpose: Calculates thickness diagnostics for PWS, stash section 20
!
!   Programming standard; Unified Model Documentation Paper No. 3
!
!   Documentation : Unified Model Documentation Paper No ??
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics

SUBROUTINE pws_thickness_diag(STASHnumber,                &
                         icode,cmessage)


! needed by typ_atm_fields.h:
USE atm_fields_bounds_mod
USE jules_snow_mod, ONLY: nsmax
USE jules_sea_seaice_mod, ONLY: nice, nice_use

USE vert_interp_mdi_mod, ONLY: vert_interp_mdi

USE um_stashcode_mod, ONLY: stashcode_pws_sec,                    &
    stashcode_pws_thickness500, stashcode_pws_thickness850
USE pws_diags_mod,    ONLY:                                       &
    pws_thickness, pws_geopht_1000, pws_geopht_850,           &
    pws_geopht_500

USE nlsizes_namelist_mod
USE lbc_mod
USE um_parvars
USE interpor_mod, ONLY: interp_order_linear
USE free_tracers_inputs_mod, ONLY: a_tracer_last, a_tracer_first

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim


IMPLICIT NONE


! Local variables

INTEGER, PARAMETER :: STASH_thkn_500 = 20001
INTEGER, PARAMETER :: STASH_thkn_850 = 20002
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'PWS_THICKNESS_DIAG'

! S/R arguments
INTEGER, INTENT(IN) :: STASHnumber        ! STASH request
INTEGER :: icode ! Return code : 0 Normal exit, >0 Error exit
CHARACTER(LEN=errormessagelength) :: cmessage  ! Error message if ICODE > 0

! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate requested thickness 
SELECT CASE(STASHnumber)

  CASE(STASH_thkn_500)   ! 1000-500mb Thickness
pws_thickness = -pws_geopht_1000 + pws_geopht_500

  CASE(STASH_thkn_850)   ! 1000-850mb Thickness
pws_thickness = -pws_geopht_1000 + pws_geopht_850

  CASE DEFAULT
   icode = -STASHnumber        ! Warning
   cmessage='Incorrect STASH request'
   CALL ereport(RoutineName,icode,cmessage)

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE pws_thickness_diag

END MODULE pws_thickness_diag_mod
