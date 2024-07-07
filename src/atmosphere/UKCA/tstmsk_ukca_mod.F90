! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Checks option code for section 34 (UKCA prognostics)

MODULE tstmsk_ukca_mod

IMPLICIT NONE
PRIVATE
PUBLIC :: tstmsk_ukca

! Description:
!   This subroutine checks if a UKCA prognostic variable (section 34)
!   is turned on or not. For variables which are on, it sets the
!   logical lmask to true.
!
! Method:
!
!   If UKCA is off, then the item is off and code exits.
!   If UKCA is on, then the option codes which are inputs to this
!   subroutine are tested.
!
!   If n30 is 0 then the item is controlled by the choice of chemistry
!   scheme in use.
!   Each option code then corresponds to a specific chemistry scheme
!   and which option code is checked depends on the choice of chemistry
!   scheme in use. For example if StratTrop chemitry is on, then option
!   code n5 is used. If n5=0 then this item is off, if n5=1 it is on
!   and if n5=2 it is on only when the optional aerosol chemistry
!   extension is in use.
!
!   If an unknown chemistry scheme is on, this subroutine calls ereport
!   with a fatal error.
!
!   If n30 is 1 then the item is controlled by the choice of aerosol
!   scheme in use.
!   If n30=1 and l_ukca_mode is false then the code exits, leaving
!   lmask false
!   If l_ukca_mode is true, then the option code tested depends on
!   the value of i_mode_setup
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: Fortran 
!   This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TSTMSK_UKCA_MOD'

CONTAINS

SUBROUTINE tstmsk_ukca(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,      &
  n13,n14,n15,n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,  &
  n29,n30,lmask)

USE ukca_option_mod,       ONLY: l_ukca, i_ukca_chem,               &
    l_ukca_mode, i_mode_setup, l_ukca_ageair
USE ukca_chem_schemes_mod, ONLY:                                    &
    i_ukca_chem_trop, i_ukca_chem_raq, i_ukca_chem_tropisop,        &
    i_ukca_chem_strattrop, i_ukca_chem_strat, i_ukca_chem_offline,  &
    i_ukca_chem_offline_be, i_ukca_chem_off
USE ereport_mod,           ONLY: ereport
USE parkind1,              ONLY: jpim, jprb
USE yomhook,               ONLY: lhook, dr_hook

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Control is by means of 30 option codes n1 to n30 for each
! STASHmaster entry

! 30 option codes are inputs to the routine
INTEGER, INTENT(IN) :: n1,  n2, n3, n4, n5, n6, n7, n8, n9,n10
INTEGER, INTENT(IN) :: n11,n12,n13,n14,n15,n16,n17,n18,n19,n20
INTEGER, INTENT(IN) :: n21,n22,n23,n24,n25,n26,n27,n28,n29,n30

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Local variables
INTEGER                       :: errcode   ! error code for ereport
INTEGER                       :: ntotal    ! sum of option codes

CHARACTER(LEN=*), PARAMETER  :: RoutineName='TSTMSK_UKCA'
CHARACTER(LEN=errormessagelength)           :: cmessage  ! used for ereport

LOGICAL, INTENT(INOUT) :: lmask

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! If all option codes are zero then we set lmask to TRUE
! for backward compatibility with other tstmsk sections

ntotal= n1+ n2+ n3+ n4+ n5+ n6+ n7+ n8+ n9+n10+ &
       n11+n12+n13+n14+n15+n16+n17+n18+n19+n20+ &
       n21+n22+n23+n24+n25+n26+n27+n28+n29+n30

IF (ntotal == 0) THEN
  lmask=.TRUE.
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! If UKCA is off, all UKCA related items are turned off so return

IF (.NOT. l_ukca) THEN
  lmask=.FALSE.
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

SELECT CASE(n30)

CASE (0)

  ! If n30=0 this is controlled by the UKCA chemistry scheme
  ! n1 = Age-of-air tracer only
  ! n2 = BE Tropospheric
  ! n3 = BE RAQ
  ! n4 = NR TropIsop
  ! n5 = NR StratTrop
  ! n6 = NR Strat
  ! n7 = NR Offline
  ! n8 = BE Offline
  SELECT CASE (i_ukca_chem)
  CASE (i_ukca_chem_off)  ! No Chem, only Age-air available
    ! Check option code 1
    CALL check_chem(n1,lmask)
  CASE (i_ukca_chem_trop) ! Standard Tropospheric (BE)
    ! Check option code 2
    CALL check_chem(n2, lmask)
  CASE (i_ukca_chem_raq) ! Regional Air Quality
    ! Check option code 3
    CALL check_chem(n3, lmask)
  CASE (i_ukca_chem_tropisop) ! Tropospheric + isoprene
    ! Check option code 4
    CALL check_chem(n4, lmask)
  CASE (i_ukca_chem_strattrop) ! Statosphere - troposphere chemistry
    ! Check option code 5
    CALL check_chem(n5, lmask)
  CASE (i_ukca_chem_strat)  ! Stratospheric chemistry
    ! Check option code 6
    CALL check_chem(n6, lmask)
  CASE (i_ukca_chem_offline)  ! Offline oxidants chemistry
    ! Check option code 7
    CALL check_chem(n7, lmask)
  CASE (i_ukca_chem_offline_be)  ! BE Offline oxidants chemistry
    ! Check option code 8
    CALL check_chem(n8, lmask)
  CASE DEFAULT
    IF ( .NOT. l_ukca_ageair ) THEN
       ! If not a known code and age-of-air not selected, 
       ! this is an error, call ereport
      errcode = 100+ABS(i_ukca_chem)
      cmessage='Unknown chemistry scheme in '//routinename
      CALL ereport(routinename,errcode,cmessage)
    END IF
  END SELECT

! n1 = age of air - available to all schemes or on its own
!    Only check if item not already selected by other options,
!    otherwise there is a risk of resetting
  IF ( (.NOT. lmask) .AND. l_ukca_ageair ) THEN
     SELECT CASE (n1)
     CASE (0)
       lmask = .FALSE.
     CASE (1)
       lmask = .TRUE.
     CASE DEFAULT
       ! This is an error, call ereport
        errcode = 1
        cmessage='Unknown n1 option code in '//routinename
        CALL ereport(routinename,errcode,cmessage)
     END SELECT
  END IF

CASE (1)

  ! If n30=1 it is controlled by the value of i_mode_setup
  ! If GLOMAP-mode is off, then off so return with lmask
  ! FALSE
  IF (.NOT. l_ukca_mode) THEN
    lmask = .FALSE.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  ! we now call check_glomap with the correct option code.
  SELECT CASE (i_mode_setup)
  CASE (1)                          ! SO4 and Sea-salt in 4 soluble modes
    CALL check_glomap(n1, lmask)
  CASE (2)                          ! SO4, BC, OC, and Sea-salt
                                    ! in 4 soluble and 1 insoluble modes
    CALL check_glomap(n2, lmask)
  CASE (6)                          ! Dust in two insoluble modes
    CALL check_glomap(n6, lmask)
  CASE (8)                          ! SO4, BC, OC, Sea-salt, and dust
                                    ! in 4 soluble and 3 insoluble modes
    CALL check_glomap(n8, lmask)
  CASE (10)                         ! SO4, BC, OC, Sea-salt, NO3 and NH4
                                    ! in 4 soluble and 1 insoluble modes
    CALL check_glomap(n10, lmask)
  CASE DEFAULT
       ! This is an error as there is no code for this aerosol
       ! scheme, call ereport
    errcode = 100+ABS(i_ukca_chem)
    cmessage='Unknown value of i_mode_setup in '//routinename
    CALL ereport(routinename,errcode,cmessage)
  END SELECT

CASE DEFAULT
     ! This is an error, call ereport
  errcode = 200
  cmessage='Unknown value of n30 in tstmsk_ukca'
  CALL ereport('tstmsk_ukca',errcode,cmessage)
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE tstmsk_ukca

! -----------------------------------------------------------------------------
! The same logic is used for all chemistry schemes
! but check a different option code 
! -----------------------------------------------------------------------------
SUBROUTINE check_chem(n, lmask)

USE ukca_option_mod,        ONLY: l_ukca_chem_aero, l_ukca,           &
                                  i_ukca_photol, l_ukca_h2o_feedback
USE ukca_photo_scheme_mod,  ONLY: i_ukca_nophot
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE parkind1,               ONLY: jpim, jprb
USE yomhook,                ONLY: lhook, dr_hook

IMPLICIT NONE

! option code being tested
INTEGER,INTENT(IN) :: n
! logical value returned by this subroutine
LOGICAL,INTENT(INOUT) :: lmask

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Local variables
INTEGER                           :: errcode   ! error code for ereport
CHARACTER(LEN=errormessagelength) :: cmessage  ! used for ereport
CHARACTER(LEN=*), PARAMETER       :: RoutineName='CHECK_CHEM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE (n)

! If the option code is zero it is always off
CASE (0)
  lmask = .FALSE.
! If the option code is one it is always on
CASE (1)
  lmask = .TRUE.
! only on if aerosol chemistry on
CASE (2)
  lmask = l_ukca_chem_aero 
CASE (3)
! Was used for 'new' emission system (NetCDF)
! Temporarily deactivated as there is only
! one emission system now 
  errcode = 3
  cmessage='Invalid STASHmaster Option code= 3 encountered'
  CALL ereport(routinename,errcode,cmessage)
CASE (4)
! only if photolysis is on
  lmask = (i_ukca_photol /= i_ukca_nophot) 
CASE (5)
! only if water vapour feedback is on
  lmask = l_ukca_h2o_feedback
CASE DEFAULT
    ! This is an error, call ereport
  errcode = 6
  cmessage='Unknown option code in '//routinename
  CALL ereport(routinename,errcode,cmessage)
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE check_chem

! -----------------------------------------------------------------------------
! The same logic is used for all values of i_mode_setup
! but check a different option code 
! -----------------------------------------------------------------------------
SUBROUTINE check_glomap(n, lmask)

USE ukca_option_mod,        ONLY: l_ukca_trophet, l_ukca_arg_act, &
                                  l_ukca_aie1, l_ukca_aie2, l_ukca_radaer
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE parkind1,               ONLY: jpim, jprb
USE yomhook,                ONLY: lhook, dr_hook

IMPLICIT NONE

! option code being tested
INTEGER,INTENT(IN) :: n
! logical value returned by this subroutine
LOGICAL,INTENT(INOUT) :: lmask

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Local variables
INTEGER                           :: errcode   ! error code for ereport
CHARACTER(LEN=errormessagelength) :: cmessage  ! used for ereport
CHARACTER(LEN=*), PARAMETER  :: RoutineName='CHECK_GLOMAP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check option code n
SELECT CASE (n)
CASE (0)
  ! If n is zero, this is always off
  lmask = .FALSE.
CASE (1)
  ! If n is 1 then it is always on
  lmask = .TRUE.
CASE (2)
  ! If n is 2 then it is on if the tropospheric heterogeneous
  ! chemistry is on
  lmask = l_ukca_trophet
CASE (3)
  ! If n is 3 then it is on if ACTIVATE is on
  lmask = l_ukca_arg_act
CASE (4)
  ! If n is 4 then it is on if aie1 or aie2 is on or
  ! ACTIVATE is on
  lmask = l_ukca_aie1 .OR. l_ukca_aie2 .OR. l_ukca_arg_act
CASE (5)
  ! If n is 5 then it is on if RADAER is on
  lmask = l_ukca_radaer
CASE DEFAULT
    ! This is an error, call ereport
  errcode = 901
  cmessage='Unknown option code in '//routinename
  CALL ereport(routinename,errcode,cmessage)
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE check_glomap

END MODULE tstmsk_ukca_mod
