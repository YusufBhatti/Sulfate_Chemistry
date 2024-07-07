! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Calculate the total number of tracers

MODULE tracer_total_var_mod

IMPLICIT NONE
!---------------------------------------------------------------------------
! Description: Called from convection to work out the total number of tracer
! variables 

! method:

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection

! code description:
!   language: fortran 90
!   This code is written to umdp3 programming standards version 10.4
!---------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TRACER_TOTAL_VAR_MOD'

CONTAINS

SUBROUTINE  tracer_total_var(tr_vars, tr_ukca,              &
                             ntra_fld, ntra_lev )

! Classic aerosol flags
USE run_aerosol_mod, ONLY: l_sulpc_so2, l_soot, l_ocff, l_biomass,  &
                           l_sulpc_nh3, l_nitrate, l_sulpc_dms

USE dust_parameters_mod, ONLY:  l_dust, l_twobin_dust

USE carbon_options_mod, ONLY: l_co2_interactive     

USE rad_input_mod, ONLY: l_use_cariolle

USE cv_run_mod,  ONLY:  l_murk_conv

USE atm_fields_bounds_mod, ONLY: tdims

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

IMPLICIT NONE

!----------------------------------------------------------------------
! Subroutine arguments
INTEGER, INTENT(IN)  ::      &
  tr_vars                    & ! number of free tracers
 ,tr_ukca                      ! number of ukca tracers

INTEGER, INTENT(INOUT) ::    &
  ntra_fld                   & ! Total number of tracers
 ,ntra_lev                     ! Number of tracer levels either 1 or model
                               ! levels

!----------------------------------------------------------------------
! Local variables

INTEGER ::                &
  errorstatus                ! Error reporting variable

CHARACTER (LEN=errormessagelength) ::                             &
  cmessage

CHARACTER(LEN=*), PARAMETER :: RoutineName='TRACER_TOTAL_VAR'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------
! Set up array dimensions for total tracer array (free + sulphur cycle
! tracers) so that convective transport of all tracers is done

IF ( l_sulpc_so2        .OR.  l_soot          .OR.  &
     l_dust             .OR.  l_biomass       .OR.  &
     l_co2_interactive  .OR.  l_use_cariolle  .OR.  &
     l_ocff             .OR.  l_nitrate       .OR.  &
     l_murk_conv ) THEN

  ntra_lev = tdims%k_end
  ntra_fld = 0        !Initialise to zero

  IF (l_dust) THEN
    IF (l_twobin_dust) THEN
      ntra_fld = ntra_fld + 2    !ADD 2 DUST SIZE CLASSES
    ELSE
      ntra_fld = ntra_fld + 6    !ADD 6 DUST SIZE CLASSES
    END IF
  END IF

      IF (l_sulpc_so2) THEN
        ntra_fld = ntra_fld + 4    !Add SO2 + 3 SO4 modes
        IF (l_sulpc_nh3) THEN
          ntra_fld = ntra_fld + 1  !Add NH3 field
        END IF
        IF (l_sulpc_dms) THEN
          ntra_fld = ntra_fld + 1  !Add DMS field
        END IF
      END IF

      IF (l_soot) THEN
        ntra_fld = ntra_fld + 3    !Add 3 modes of soot
      END IF

      IF (l_biomass) THEN
        ntra_fld = ntra_fld + 3    !Add 3 modes of biomass aerosol
      END IF

      IF (l_co2_interactive) THEN
        ntra_fld = ntra_fld + 1    !Add CO2 field
      END IF

      IF (l_ocff) THEN
        ntra_fld = ntra_fld + 3    !Add 3 modes of fossil-fuel OC
      END IF

      IF (l_nitrate) THEN
        ntra_fld = ntra_fld + 2    !Add 2 modes of ammonium nitrate
      END IF

      IF (l_use_cariolle) THEN
        ntra_fld = ntra_fld + 1    !Add Cariolle Ozone tracer field
      END IF

      IF (l_murk_conv) THEN
        ntra_fld = ntra_fld + 1    !ADD MURK aerosol
      END IF

      IF (tr_vars > 0) THEN
        ntra_fld = ntra_fld + tr_vars
      END IF

      IF (tr_ukca > 0) THEN
        ntra_fld = ntra_fld + tr_ukca
      END IF

ELSE
  IF (tr_vars == 0 .AND. tr_ukca == 0) THEN
    ntra_fld = 1             ! can't have 0 sized arrays
    ntra_lev = 1
  ELSE
    ntra_fld = tr_vars + tr_ukca
    ntra_lev = tdims%k_end
  END IF
END IF


!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!----------------------------------------------------------------------
RETURN
END SUBROUTINE tracer_total_var

END MODULE tracer_total_var_mod
