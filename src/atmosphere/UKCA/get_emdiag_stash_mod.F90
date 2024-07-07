! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Get STASH item number (in section 50) for a given diagnostics
!
! Method:
!   Get the name of the diagnostic as input argument and use a
!   CASE SELECT statement to return the corresponding STASH
!   item number in section 50.
!
! Part of the UKCA model, a community model supported by the
! Met Office and NCAS, with components provided initially
! by The University of Cambridge, University of Leeds and
! The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE get_emdiag_stash_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GET_EMDIAG_STASH_MOD'

CONTAINS

FUNCTION get_emdiag_stash (diag_name)

USE ukca_constants
USE ukca_emiss_mode_mod, ONLY: aero_ems_species
USE ereport_mod,     ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE parkind1,        ONLY: jpim,  jprb     ! DrHook
USE yomhook,         ONLY: lhook, dr_hook  ! DrHook

IMPLICIT NONE

! Function arguments
CHARACTER (LEN=10), INTENT(IN) :: diag_name ! name of the diagnostic
                                            ! (e.g. emission field)
! Function return value
REAL :: get_emdiag_stash

! Local variables
INTEGER                            :: ierr
CHARACTER (LEN=errormessagelength) :: cmessage

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0  ! DrHook tracing entry
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1  ! DrHook tracing exit
REAL    (KIND=jprb)            :: zhook_handle   ! DrHook tracing

CHARACTER(LEN=*), PARAMETER :: RoutineName='GET_EMDIAG_STASH'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

ierr = 0

! Return STASH item in Section 50 for a given UKCA diagnostic.
! The item numbers below should agree with STASHmaster_A.
SELECT CASE ( diag_name )
CASE ('NO        ')
  get_emdiag_stash = 156

CASE ('CH4       ')
  get_emdiag_stash = 157

CASE ('CO        ')
  get_emdiag_stash = 158

CASE ('HCHO      ')
  get_emdiag_stash = 159

CASE ('C2H6      ')
  get_emdiag_stash = 160

CASE ('C3H8      ')
  get_emdiag_stash = 161

CASE ('Me2CO     ')
  get_emdiag_stash = 162

CASE ('MeCHO     ')
  get_emdiag_stash = 163

CASE ('C5H8      ')
  get_emdiag_stash = 164

CASE ('C4H10     ')
  get_emdiag_stash = 165

CASE ('C2H4      ')
  get_emdiag_stash = 166

CASE ('C3H6      ')
  get_emdiag_stash = 167

CASE ('TOLUENE   ')
  get_emdiag_stash = 168

CASE ('oXYLENE   ')
  get_emdiag_stash = 169

CASE ('CH3OH     ')
  get_emdiag_stash = 170

CASE ('H2        ')
  get_emdiag_stash = 171

CASE ('NO_aircrft')
  get_emdiag_stash = 172

CASE ('Monoterp  ')
  get_emdiag_stash = 211

CASE ('NVOC      ')
  get_emdiag_stash = 212

CASE ('NH3       ')
  get_emdiag_stash = 213

CASE ('DMS       ')
  get_emdiag_stash = 214

CASE ('SO2_low   ')
  get_emdiag_stash = 215

CASE ('SO2_high  ')
  get_emdiag_stash = 216

CASE ('SO2_nat   ')
  get_emdiag_stash = 217

CASE DEFAULT
  ! Report error unless this is an aerosol emission (BC_fossil:, etc)
  ! No diagnostics for these source emissions because the resulting
  ! mode-by-mode emissions are available in section 38.
  IF (.NOT. ANY(aero_ems_species == diag_name))  ierr = 1
END SELECT

! Report error if diagnostics not found
IF (ierr == 1) THEN
  cmessage = 'Unexpected UKCA emiss diagnostic: ' // TRIM (diag_name)
  CALL ereport ('GET_EMDIAG_STASH', ierr, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)

RETURN

END FUNCTION get_emdiag_stash

END MODULE get_emdiag_stash_mod
