! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE ukca_emiss_diags_mode_mod

USE parkind1, ONLY: jprb,  jpim
USE yomhook,  ONLY: lhook, dr_hook

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY:                                                     &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf

IMPLICIT NONE
PRIVATE
PUBLIC :: ukca_emiss_diags_mode

!  Description:
!    Produce emission diagnostics for GLOMAP-mode emissions under 
!    the new UKCA emission system (based on NetCDF input emissions).
!
!  Method: 
!    Information about all available diagnostics is contained in the
!    mode_diag array, which links each STASH item code to mode and component.
!    This array is filled by ukca_emiss_diags_mode_init on the first timestep.
!    The routine ukca_emiss_diags_mode loops over each diagnostic, finds
!    corresponding mass emissions and copies the diags elements of this into a
!    local array em_diags. em_diags is then copied to stash in copydiag_3d.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

TYPE :: type_mode_diag_struct
  INTEGER :: item  ! stash item number
  INTEGER :: mode  ! diag's mode
  INTEGER :: component  ! diag's component
END TYPE type_mode_diag_struct

INTEGER, PARAMETER :: num_mode_diags = 13
TYPE(type_mode_diag_struct) :: mode_diag(num_mode_diags)

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_EMISS_DIAGS_MODE_MOD'

CONTAINS


SUBROUTINE ukca_emiss_diags_mode_init

USE ukca_mode_setup,  ONLY: cp_su, cp_bc, cp_oc, cp_cl, cp_du

IMPLICIT NONE

INTEGER :: icount
CHARACTER (LEN=errormessagelength) :: cmessage

icount = 0
mode_diag(:)%item = -1
mode_diag(:)%mode = -1
mode_diag(:)%component = -1

! SO4 to aitken-sol
icount = icount + 1
mode_diag(icount)%item = 201
mode_diag(icount)%mode = 2
mode_diag(icount)%component = cp_su

! SO4 to accum-sol
icount = icount + 1
mode_diag(icount)%item = 202
mode_diag(icount)%mode = 3
mode_diag(icount)%component = cp_su

! SO4 to coarse-sol
icount = icount + 1
mode_diag(icount)%item = 203
mode_diag(icount)%mode = 4
mode_diag(icount)%component = cp_su

! sea-salt to accum-sol
icount = icount + 1
mode_diag(icount)%item = 204
mode_diag(icount)%mode = 3
mode_diag(icount)%component = cp_cl

! sea-salt to coarse-sol
icount = icount + 1
mode_diag(icount)%item = 205
mode_diag(icount)%mode = 4
mode_diag(icount)%component = cp_cl

! black carbon to aitken-sol
icount = icount + 1
mode_diag(icount)%item = 206
mode_diag(icount)%mode = 2
mode_diag(icount)%component = cp_bc

! black carbon to aitken-ins
icount = icount + 1
mode_diag(icount)%item = 207
mode_diag(icount)%mode = 5
mode_diag(icount)%component = cp_bc

! organic carbon to aitken-sol
icount = icount + 1
mode_diag(icount)%item = 208
mode_diag(icount)%mode = 2
mode_diag(icount)%component = cp_oc

! organic carbon to aitken-ins
icount = icount + 1
mode_diag(icount)%item = 209
mode_diag(icount)%mode = 5
mode_diag(icount)%component = cp_oc

! dust to accum-sol
icount = icount + 1
mode_diag(icount)%item = 210
mode_diag(icount)%mode = 3
mode_diag(icount)%component = cp_du

! dust to accum-ins
icount = icount + 1
mode_diag(icount)%item = 211
mode_diag(icount)%mode = 6
mode_diag(icount)%component = cp_du

! dust to coarse-sol
icount = icount + 1
mode_diag(icount)%item = 212
mode_diag(icount)%mode = 4
mode_diag(icount)%component = cp_du

! dust to coarse-ins
icount = icount + 1
mode_diag(icount)%item = 213
mode_diag(icount)%mode = 7
mode_diag(icount)%component = cp_du

! check that number of diagnostics matches the length of the array
IF (icount /= num_mode_diags) THEN
  cmessage = "Number of mode diagnostics wrong: icount /= num_mode_diags"
  CALL ereport ('UKCA_EMISS_DIAGS_MODE_INIT', icount, cmessage)
END IF

END SUBROUTINE ukca_emiss_diags_mode_init


SUBROUTINE ukca_emiss_diags_mode (row_length, rows, model_levels, area, &
    len_stashwork, stashwork)

USE ukca_mode_setup,    ONLY: mm

USE ukca_emiss_mod,     ONLY: emissions, num_em_flds

USE um_parvars,         ONLY: at_extremity
USE um_stashcode_mod,   ONLY: stashcode_glomap_sec
USE yomhook,            ONLY: lhook, dr_hook
USE parkind1,           ONLY: jprb, jpim

IMPLICIT NONE


! Subroutine arguments
INTEGER, INTENT(IN)    :: row_length        ! Model dimensions
INTEGER, INTENT(IN)    :: rows
INTEGER, INTENT(IN)    :: model_levels
REAL, INTENT(IN) :: area(row_length, rows, model_levels)  ! Area of grid cell

INTEGER, INTENT(IN)    :: len_stashwork     ! Length of diagnostics array
REAL,    INTENT(INOUT) :: stashwork (len_stashwork) ! Diagnostics array

! Local variables
INTEGER    :: k, l, n      ! counters / indices
INTEGER    :: section      ! stash section
INTEGER    :: item         ! stash item
INTEGER    :: icode = 0    ! local error status
INTEGER    :: im_index     ! internal model index
INTEGER    :: imode
INTEGER    :: icp
INTEGER    :: ilev
INTEGER    :: err_code

REAL :: em_diags (row_length, rows, model_levels)

LOGICAL, SAVE :: lfirst = .TRUE.  ! Indicator of first call to this routine
LOGICAL       :: lfound = .TRUE.  ! Has corresponding emission been found

CHARACTER(LEN=errormessagelength) :: cmessage  ! Error return message

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EMISS_DIAGS_MODE'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (lfirst) THEN
  CALL ukca_emiss_diags_mode_init
  lfirst = .FALSE.
END IF

im_index = 1
section = stashcode_glomap_sec

! ------------------------------------------------------------------------
! For each diagnostic add contributions from each source emission, then
! copy to corresponding stash item.
DO k = 1, num_mode_diags
  imode = mode_diag(k)%mode
  icp = mode_diag(k)%component
  item = mode_diag(k)%item

  IF (sf(item, section)) THEN
    ! Loop through all emissions to store the diagnostics for this mode/cpt
    em_diags(:,:,:) = 0.0
    lfound = .FALSE.
    DO l = 1, num_em_flds

      IF ((emissions(l)%tracer_name == 'mode_emiss') .AND. &
          (emissions(l)%moment == 3) .AND. &
          (emissions(l)%mode == imode) .AND. &
          (emissions(l)%component == icp)) THEN

        lfound = .TRUE.
        IF (emissions(l)%three_dim) THEN
          em_diags(:,:,:) = em_diags(:,:,:) + emissions(l)%diags(:,:,:)
        ELSE
          ! For 2D emissions, diags stored in level 1 but vertical scaling
          ! stored and we use it to re-project onto levels. Column total of 
          ! vert_scaling_3d is always equal to 1.
          DO ilev = 1, model_levels
            em_diags(:,:,ilev) = em_diags(:,:,ilev) +   &
              emissions(l)%diags(:,:,1) * emissions(l)%vert_scaling_3d(:,:,ilev)
          END DO
        END IF  ! three_dim
      END IF  ! if this emission matches the diagnostic

    END DO  ! loop over emissions

    ! Raise a warning (and skip to next diag) if no emissions for this diag
    IF (.NOT. lfound) THEN
      cmessage = 'No mode emissions corresponding to section 38 ' // &
                 'diagnostic request'
      err_code = -item
      CALL ereport ('UKCA_EMISS_DIAGS_MODE', err_code, cmessage)
      CYCLE
    END IF

    ! Convert from kg/m2/s to mol/gridbox/s  
    em_diags(:,:,:) = em_diags(:,:,:) * area(:,:,:) / mm(icp)

    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d (stashwork (si(item,section,im_index)),           &
             em_diags (:,:,:),                                         &
             row_length, rows,model_levels, 0,0,0,0, at_extremity,     &
             stlist(1,stindex(1,item,section,im_index)), len_stlist,   &
             stash_levels, num_stash_levels+1,                         &
             atmos_im, section, item, icode, cmessage)

    IF (icode >  0) THEN
      err_code = section*1000+item
      CALL ereport ('UKCA_EMISS_DIAGS_MODE', err_code, cmessage)
    END IF

  END IF  ! if diagnostic requested
END DO  ! loop over diags
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_emiss_diags_mode

END MODULE ukca_emiss_diags_mode_mod
