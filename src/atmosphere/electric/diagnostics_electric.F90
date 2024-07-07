! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE diagnostics_electric_mod

! Purpose: Calculates diagnostics for the electric scheme

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Electric

USE atm_fields_bounds_mod,  ONLY: tdims
USE electric_constants_mod, ONLY: i,j,k
USE ereport_mod,            ONLY: ereport
USE submodel_mod,           ONLY: atmos_im
USE stash_array_mod,        ONLY: si, sf
USE errormessagelength_mod, ONLY: errormessagelength

! Dr Hook modules
USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim

IMPLICIT NONE


! Method:
!  Required level lists and logical switches are determined by the
! calling routine from STASH requests and STASHflags.
! Intercepted arrays - calculated unconditionally - and diagnostic
! arrays - dependent on STASHflags - are input from the large scale
! precipitation routines called previously. Each diagnostic is simply
! copied into the STASHwork array to be passed on to STASH for
! output processing, except where indicated.

!  Diagnostics currently available: (in order calculated)

! STASH items (all section 21)
! ------------------------------------------------------------
! 100: lightning flash rate (s-1)
! 101: Storm location flag
! 102: Graupel Water Path (GWP; kg m-2)
! 103: Total Ice Water Path  (TIWP; kg m-2)
! 104: Total number of flashes
! 105: Flash rate from graupel flux (McCaul Scheme)
! 106: Flash rate from storm ice (McCaul Scheme)
!-------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DIAGNOSTICS_ELECTRIC_MOD'

CONTAINS

SUBROUTINE diagnostics_electric( flash, storm_field, gwp, tiwp, num_flashes,  &
                                 fr1_mc, fr2_mc, stashwork,                   &
                                 at_extremity )

IMPLICIT NONE

REAL, INTENT(IN) :: flash( tdims%i_start : tdims%i_end,                       &
                           tdims%j_start : tdims%j_end )

REAL, INTENT(IN) :: gwp(   tdims%i_start : tdims%i_end,                       &
                           tdims%j_start : tdims%j_end )

REAL, INTENT(IN) :: tiwp(  tdims%i_start : tdims%i_end,                       &
                           tdims%j_start : tdims%j_end )

LOGICAL, INTENT(IN) :: storm_field (tdims%i_start : tdims%i_end,              &
                                  tdims%j_start : tdims%j_end )

LOGICAL, INTENT(IN) :: at_extremity(4)
! Indicates if this processor is at north, south, east or west of the
! processor grid)

REAL, INTENT(IN) :: num_flashes( tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end )

REAL, INTENT(IN) :: fr1_mc( tdims%i_start : tdims%i_end,                      &
                            tdims%j_start : tdims%j_end )

REAL, INTENT(IN) :: fr2_mc( tdims%i_start : tdims%i_end,                      &
                            tdims%j_start : tdims%j_end )


REAL, INTENT(INOUT) :: stashwork(*)
! STASH workspace for electric section

!============================================================================
! Local variables
!============================================================================

INTEGER :: icode
INTEGER :: item

INTEGER, PARAMETER :: sect = 21

CHARACTER(LEN=errormessagelength) :: cmessage

CHARACTER(*), PARAMETER :: RoutineName='DIAGNOSTICS_ELECTRIC'

INTEGER :: im_index        ! internal model index

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!==================================================================
! Start the subroutine
!==================================================================

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode    = 0 ! Initialise error status
im_index = 1

!==================================================================
!  Output diagnostics
!  Copy diagnostic information to STASHwork for STASH processing
!==================================================================

!--------------------------------------------------------------------------
! Item 100: Flash rate
!--------------------------------------------------------------------------

item = 100

! Copy flash rate  to STASHwork

IF (icode <= 0 .AND. sf(item, sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag( stashwork(si(item, sect, im_index)), flash,                 &
                 tdims%i_end, tdims%j_end, 0, 0, 0, 0, at_extremity,         &
                 atmos_im, sect, item, icode, cmessage               )

  IF (icode >  0) THEN
    cmessage=":error in copydiag(flash)"
  END IF

END IF

!--------------------------------------------------------------------------
! Item 101: Storm location flag
!--------------------------------------------------------------------------

item = 101

IF (icode <= 0 .AND. sf(item, sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag( stashwork(si(item, sect, im_index)), storm_field,           &
                 tdims%i_end, tdims%j_end, 0, 0, 0, 0, at_extremity,         &
                 atmos_im, sect, item, icode, cmessage               )

  IF (icode >  0) THEN
    cmessage=":error in copydiag(storm_field)"
  END IF

END IF

!--------------------------------------------------------------------------
! Item 102: Graupel Water Path (GWP)
!--------------------------------------------------------------------------

item = 102

! Copy GWP  to STASHwork

IF (icode <= 0 .AND. sf(item, sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag( stashwork(si(item, sect, im_index)), gwp,                   &
                 tdims%i_end, tdims%j_end, 0, 0, 0, 0, at_extremity,         &
                 atmos_im, sect, item, icode, cmessage               )

  IF (icode >  0) THEN
    cmessage=":error in copydiag(gwp)"
  END IF

END IF

!--------------------------------------------------------------------------
! Item 103: Total Ice Water Path (GWP)
!--------------------------------------------------------------------------

item = 103

! Copy TIWP  to STASHwork

IF (icode <= 0 .AND. sf(item, sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag( stashwork(si(item, sect, im_index)), tiwp,                  &
                 tdims%i_end, tdims%j_end, 0, 0, 0, 0, at_extremity,         &
                 atmos_im, sect, item, icode, cmessage               )

  IF (icode >  0) THEN
    cmessage=":error in copydiag(tiwp)"
  END IF

END IF

!--------------------------------------------------------------------------
! Item 104: NUMBER OF FLASHES
!--------------------------------------------------------------------------

item = 104

! Copy number of flashes to STASHwork

IF (icode <= 0 .AND. sf(item, sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag( stashwork(si(item, sect, im_index)), num_flashes,           &
                 tdims%i_end, tdims%j_end, 0, 0, 0, 0, at_extremity,         &
                 atmos_im, sect, item, icode, cmessage               )

  IF (icode >  0) THEN
    cmessage=":error in copydiag(num_flashes)"
  END IF

END IF

!--------------------------------------------------------------------------
! Item 105: Flash rate from flux (McCaul scheme)
!--------------------------------------------------------------------------

item = 105

! Copy flash rate  to STASHwork

IF (icode <= 0 .AND. sf(item, sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag( stashwork(si(item, sect, im_index)), fr1_mc,                &
                 tdims%i_end, tdims%j_end, 0, 0, 0, 0, at_extremity,         &
                 atmos_im, sect, item, icode, cmessage               )

  IF (icode >  0) THEN
    cmessage=":error in copydiag(flash)"
  END IF

END IF

!--------------------------------------------------------------------------
! Item 106: Flash rate from ice water path (McCaul scheme)
!--------------------------------------------------------------------------

item = 106

! Copy flash rate  to STASHwork

IF (icode <= 0 .AND. sf(item, sect)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag( stashwork(si(item, sect, im_index)), fr2_mc,                &
                 tdims%i_end, tdims%j_end, 0, 0, 0, 0, at_extremity,         &
                 atmos_im, sect, item, icode, cmessage               )

  IF (icode >  0) THEN
    cmessage=":error in copydiag(flash)"
  END IF

END IF

! Single point exception handling
IF (icode /= 0) THEN

  CALL ereport(RoutineName,icode,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE diagnostics_electric

END MODULE diagnostics_electric_mod
