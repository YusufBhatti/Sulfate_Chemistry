! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Purpose: Species labelled 'CF': ASAD will treat the species as a constant
!     but will call this routine so that the user may set the
!     values differently at each gridpoint for example. Currently
!     used for setting the water vapour and CO2 concentrations, as well
!     as the offline oxidants.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_FYINIT
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!     Interface
!     On entry, the following will be set:
!              species - character name of species to set.
!                        Will be the same as listed in chch.d file
!              klen    - length of array, y.
!
!     On exit, the following must be set:
!              y       - Array of points to set for the species.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE asad_inicnt_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_INICNT_MOD'

CONTAINS

SUBROUTINE asad_inicnt( species, y, klen, nlev )

USE asad_mod,              ONLY: wp, co2, tnd
USE ukca_option_mod,       ONLY: l_ukca_offline, l_ukca_offline_be
USE ukca_constants,        ONLY: c_oh, c_o3, c_no3, c_ho2
USE ukca_chem_offline,     ONLY: o3_offline, oh_offline,        &
                                 no3_offline, ho2_offline
USE atm_fields_bounds_mod, ONLY: array_dims, tdims
USE carbon_options_mod,    ONLY: l_co2_interactive
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: klen      ! No of spatial points
INTEGER, INTENT(IN) :: nlev      ! Model level

CHARACTER (LEN=10), INTENT(IN)  :: species  ! Species char strng

REAL, INTENT(OUT)   :: y(klen)   ! Species concentration this may
                                 ! be in volumetric mixing ratio
                                 ! units (H2O) or as mass mixing
                                 ! ratio (offline oxidants)

!       Local variables

INTEGER :: errcode                ! Variable passed to ereport
INTEGER :: row_length             ! row_length for theta field
INTEGER :: rows                   ! rows for theta field

CHARACTER (LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_INICNT'


!     1.  Copy water, CO2, and offline oxidants (if required) into ASAD array.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

row_length = tdims%i_end
rows = tdims%j_end

IF (species(1:4) == 'CO2 ' .AND. l_co2_interactive) THEN
!  The CO2 field is set to the UM prognostic if it is available and
!  is used in the chemical scheme.
  y(:) = co2(:)
ELSE IF ( species(1:4) == 'H2O ' ) THEN
  ! Note that wp is in units of volumetric mixing ratio
  y(:) = wp(:)
ELSE IF (l_ukca_offline .OR. l_ukca_offline_be) THEN
  ! These species are converted from mass mixing ratio to vmr
  IF (species(1:4) == 'OH  ') THEN
    y(:) = RESHAPE(oh_offline(1:row_length,1:rows,nlev),(/klen/))
    y(:) = y(:)/c_oh
  ELSE IF (species(1:4) == 'O3  ') THEN
    y(:) = RESHAPE(o3_offline(1:row_length,1:rows,nlev),(/klen/))
    y(:) = y(:)/c_o3
  ELSE IF (species(1:4) == 'NO3 ') THEN
    y(:) = RESHAPE(no3_offline(1:row_length,1:rows,nlev),(/klen/))
    y(:) = y(:)/c_no3
  ELSE IF (species(1:4) == 'HO2 ') THEN
    y(:) = RESHAPE(ho2_offline(1:row_length,1:rows,nlev),(/klen/))
    y(:) = y(:)/c_ho2
  ELSE
    errcode = 125
    cmessage = ' Species '//species//' is not treated by this routine'
    CALL ereport('ASAD_INICNT',errcode,cmessage)
  END IF
ELSE
  errcode=124
  cmessage= ' Species '//species//' not treated by this routine'
  CALL ereport('ASAD_INICNT',errcode,cmessage)
END IF

! Convert to molecules/cm^3 from vmr
y(:) = y(:)*tnd(:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE asad_inicnt
END MODULE asad_inicnt_mod
