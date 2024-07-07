! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Reads the run_free_tracers namelist
!
MODULE rcf_readnl_runfreetracers_mod

IMPLICIT NONE

! Description:
!  rcf_readnl_runfreetracers reads in the run_free_tracers namelist.
!
! Method:
!  Data read in via unit number passed as argument.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!  Language: FORTRAN 95.
!  This runcode is written to UMDP3 v8.6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_READNL_RUNFREETRACERS_MOD'

CONTAINS

SUBROUTINE rcf_readnl_runfreetracers (nft)

USE free_tracers_inputs_mod,  ONLY: run_free_tracers, i_free_tracer,        &
                                    print_nlist_run_free_tracers,           &
                                    a_max_trvars, i_free_tracer_lbc,        &
                                    read_nml_run_free_tracers
USE umPrintMgr,               ONLY: PrintStatus, PrStatus_Oper
USE UM_parcore,               ONLY: mype
USE run_aerosol_mod,          ONLY: l_soot, l_biomass, l_sulpc_so2, l_ocff, &
                                    l_nitrate
USE dust_parameters_mod,      ONLY: i_dust
USE ukca_option_mod,          ONLY: l_ukca
USE nlsizes_namelist_mod,     ONLY: model_levels, tr_levels, tr_vars
USE ukca_option_mod,          ONLY: tc_lbc_ukca
USE ukca_tracer_stash,        ONLY: a_max_ukcavars
USE missing_data_mod,         ONLY: imdi
USE rcf_grid_type_mod,        ONLY: output_grid

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)  :: nft  ! unit number

! Local variables
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_READNL_RUNFREETRACERS'

INTEGER :: i
INTEGER :: tr_lbc_vars
INTEGER :: tr_lbc_ukca
LOGICAL :: l_aero
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! Read the namelist
CALL read_nml_run_free_tracers(nft)

IF ( PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  CALL print_nlist_run_free_tracers()
END IF

! If atmospheric tracers (section 33) are being used
! then count the number of tracers actually being used
tr_vars=0
DO i=1,a_max_trvars
  IF (i_free_tracer(i) /= 0)tr_vars = tr_vars + 1
END DO

! Find the number of tracers being used by counting the non-zero elements
! of the arrays that list them
tr_lbc_vars=0
tr_lbc_ukca=0
DO i=1,a_max_trvars
  IF (i_free_tracer_lbc(i) /= 0) tr_lbc_vars = tr_lbc_vars + 1
END DO
DO i=1,a_max_ukcavars
  IF (tc_lbc_ukca(i) /= 0) tr_lbc_ukca = tr_lbc_ukca + 1
END DO

! If no tracers have been selected and if aerosols (section 17) aren't being
! used and if UKCA (section 34) isn't being used then reset the number of
! tracer levels to 1 - otherwise, use model_levels

l_aero = (l_soot .OR. l_biomass .OR.                                 &
         (i_dust /= imdi .AND. i_dust /=0 ) .OR.                     &
          l_sulpc_so2 .OR. l_ocff .OR. l_nitrate)

IF (tr_vars==0 .AND. tr_lbc_vars==0 .AND. tr_lbc_ukca==0              &
    .AND. .NOT. l_aero .AND. .NOT. l_ukca) THEN
  tr_levels=1
ELSE
  tr_levels=model_levels
END IF
output_grid % tr_levels         = tr_levels

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_readnl_runfreetracers
END MODULE rcf_readnl_runfreetracers_mod
