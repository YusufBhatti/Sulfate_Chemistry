! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE ukca_emiss_diags_mod

USE parkind1,     ONLY: jprb,  jpim
USE yomhook,      ONLY: lhook, dr_hook
USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY:                                                     &
    len_stlist, stindex, stlist, num_stash_levels, stash_levels, si, sf

IMPLICIT NONE

!  Description:
!    Produce emission diagnostics for the new UKCA emission system
!    (based on NetCDF input emissions).
!
!  Method:
!    Go through the list of emitted species in the chemistry scheme
!    (em_chem_spec). Look for all emission diagnostics present for
!    each of these species in emissions(:)%diags, add them up and
!    store them as a total emission diagnostic in the allocatable
!    array 'em_diags' (note that this code allows the use of several
!    emission fields for a given tracer). Call UPDATE_EMDIAGSTRUCT
!    to copy  'em_diags' to the relevant pointers of the 'emdiags'
!    structure. Finally, go diagnostic by diagnostic, if selected
!    then copy to STASHwork via calls to COPYDIAG (for 2D diags)
!    and COPYDIAG_3D (for 3D diags).
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


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_EMISS_DIAGS_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE ukca_emiss_diags (row_length, rows, model_levels,        &
    len_stashwork, stashwork)

USE ukca_d1_defs,         ONLY: UKCA_diag_sect, em_chem_spec,       &
                                n_chem_emissions, n_3d_emissions
USE ukca_option_mod,              ONLY: l_ukca_primbcoc
USE ukca_emiss_mod,               ONLY: emissions, num_em_flds
USE ukca_update_emdiagstruct_mod, ONLY: update_emdiagstruct, emdiags
USE get_emdiag_stash_mod,         ONLY: get_emdiag_stash

USE ereport_mod,        ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE um_parvars,         ONLY: at_extremity
USE yomhook,            ONLY: lhook, dr_hook
USE parkind1,           ONLY: jprb, jpim

IMPLICIT NONE


! Subroutine arguments
INTEGER, INTENT(IN)    :: row_length        ! Model dimensions
INTEGER, INTENT(IN)    :: rows
INTEGER, INTENT(IN)    :: model_levels
INTEGER, INTENT(IN)    :: len_stashwork     ! Length of diagnostics array

REAL,    INTENT(INOUT) :: stashwork (len_stashwork) ! Diagnostics array

! Local variables
INTEGER    :: k, l, n      ! counters / indices
INTEGER    :: section      ! stash section
INTEGER    :: item         ! stash item
INTEGER    :: icode = 0    ! local error status
INTEGER    :: im_index     ! internal model index
INTEGER    :: errcode      ! Error code for ereport
REAL, ALLOCATABLE :: em_diags (:,:,:)

CHARACTER(LEN=errormessagelength) :: cmessage ! Error return message

INTEGER (KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER :: zhook_out = 1
REAL    (KIND=jprb)            :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EMISS_DIAGS'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ------------------------------------------------------------------------
! Loop through all emissions to store the diagnostics
DO k = 1, n_chem_emissions + n_3d_emissions

  ! Do initial loop over num_em_flds to get correct size of em_diags
  ! array that is required. This is needed here as BC & OC biomass
  ! netCDF fields are either 2D or 3D, so we cannot assume size of
  ! the array based on the name of the emitted species alone.
  size_loop: DO l = 1, num_em_flds    ! fields in emissions structure
     IF (em_chem_spec (k) == emissions(l)%tracer_name) THEN
        ALLOCATE(em_diags(row_length, rows,                          &
                 SIZE(emissions(l)%diags(:,:,:),DIM=3)))
        EXIT size_loop
     END IF
  END DO size_loop
  
  ! Set total diagnostic to zero every time we look for a new tracer
  em_diags (:,:,:) = 0.0

  ! Get index by matching emiss_chem_spec with emissions%tracer_name
  n = -1
  DO l = 1, num_em_flds    ! fields in emissions structure
    IF (em_chem_spec (k) == emissions(l)%tracer_name) THEN
      n = l
      ! When a field is found in the emissions structure then add
      ! the corresponding diagnostic to 'em_diags'

      ! First check sizes of arrays are correct to trap any incorrect
      ! array allocation. This should be trapped above, but this only
      ! checks on the first instance of %tracer_name in the emissions
      ! array, so be careful and check here again.
      IF (SIZE(em_diags(:,:,:),DIM=3) /=                              &
          SIZE(emissions(n)%diags(:,:,:),DIM=3)) THEN
         icode    = n
         cmessage = "Incorrect number of levels defined for " //      &
                    TRIM (em_chem_spec (k)) // " emissions diagnostic"

         CALL ereport ('UKCA_EMISS_DIAGS', icode, cmessage)
      END IF
      
      em_diags (:,:,:) = em_diags (:,:,:) + emissions(n)%diags(:,:,:)
    END IF
  END DO

  ! Report error if no emission field found in emissions structure
  ! Exception is BC/OC which can be missing if l_ukca_primbcoc=false
  IF (n == -1) THEN
    IF (.NOT. (                                                    &
        ((em_chem_spec(k)(1:2) == 'BC') .AND. .NOT. l_ukca_primbcoc) .OR.  &
        ((em_chem_spec(k)(1:2) == 'OC') .AND. .NOT. l_ukca_primbcoc) ) ) THEN
      icode    = -n
      cmessage = "Emission field " // TRIM (em_chem_spec (k)) //         &
                 " not found in emissions structure"

      CALL ereport ('UKCA_EMISS_DIAGS', icode, cmessage)
    END IF
  END IF

  ! Store 'em_diags' in the pointer of the emdiags structure
  ! that corresponds to the current diagnostic.
  CALL update_emdiagstruct (                                 &
         row_length, rows, model_levels,             &
         em_diags, em_chem_spec (k))

  ! Deallocate em_diags to be able to use it for a new diagnostic
  DEALLOCATE (em_diags)

END DO    ! 1, n_chem_emissions + n_3d_emissions


!---------------------------------------------------------------------
! Check that item numbers in UKCA section have been selected
! via STASH. If so then call copydiag or copydiag_3d.
im_index = 1
section  = UKCA_diag_sect

!---------------------------------------------------------------------
! Sec 50, item 156: NO surface emissions
item = get_emdiag_stash ('NO        ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_no (:,:),                                &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 157: CH4 surface emissions
item = get_emdiag_stash ('CH4       ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_ch4 (:,:),                               &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 158: CO surface emissions
item = get_emdiag_stash ('CO        ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_co (:,:),                                &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 159: HCHO surface emissions
item = get_emdiag_stash ('HCHO      ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_hcho (:,:),                              &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 160: C2H6 surface emissions
item = get_emdiag_stash ('C2H6      ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_c2h6 (:,:),                              &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 161: C3H8 surface emissions
item = get_emdiag_stash ('C3H8      ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_c3h8 (:,:),                              &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 162: Me2CO surface emissions
item = get_emdiag_stash ('Me2CO     ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_me2co (:,:),                             &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 163: MeCHO surface emissions
item = get_emdiag_stash ('MeCHO     ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_mecho (:,:),                             &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 164: C5H8 surface emissions
item = get_emdiag_stash ('C5H8      ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_c5h8 (:,:),                              &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 165: C4H10 surface emissions
item = get_emdiag_stash ('C4H10     ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_c4h10 (:,:),                             &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 166: C2H4 surface emissions
item = get_emdiag_stash ('C2H4      ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_c2h4 (:,:),                              &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 167: C3H6 surface emissions
item = get_emdiag_stash ('C3H6      ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_c3h6 (:,:),                              &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 168: TOLUENE surface emissions
item = get_emdiag_stash ('TOLUENE   ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_tol (:,:),                               &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 169: oXYLENE surface emissions
item = get_emdiag_stash ('oXYLENE   ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_oxyl (:,:),                              &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 170: CH3OH surface emissions
item = get_emdiag_stash ('CH3OH     ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_ch3oh (:,:),                             &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 171: H2 surface emissions
item = get_emdiag_stash ('H2        ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_h2 (:,:),                                &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 172: NO_aircrft 3D emissions
item = get_emdiag_stash ('NO_aircrft')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork (si(item,section,im_index)),           &
           emdiags%em_no_air (:,:,:),                                &
           row_length, rows,model_levels, 0,0,0,0, at_extremity,     &
           stlist(1,stindex(1,item,section,im_index)), len_stlist,   &
           stash_levels, num_stash_levels+1,                         &
           atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 211: Monoterp surface emissions
item = get_emdiag_stash ('Monoterp  ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_montrp (:,:),                            &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 212: NVOC surface emissions
item = get_emdiag_stash ('NVOC      ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_nvoc (:,:),                              &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 213: NH3 surface emissions
item = get_emdiag_stash ('NH3       ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_nh3 (:,:),                               &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 214: DMS surface emissions
item = get_emdiag_stash ('DMS       ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_dms (:,:),                               &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 215: SO2 surface emissions
item = get_emdiag_stash ('SO2_low   ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_so2low (:,:),                            &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 216: SO2 high level emissions
item = get_emdiag_stash ('SO2_high  ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag
  CALL copydiag (stashwork (si(item,section,im_index)),              &
                 emdiags%em_so2hi (:,:),                             &
                 row_length, rows, 0,0,0,0, at_extremity,            &
                 atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

!---------------------------------------------------------------------
! Sec 50, item 217: SO2 natural 3D emissions
item = get_emdiag_stash ('SO2_nat   ')
IF (sf(item, section)) THEN

  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork (si(item,section,im_index)),           &
           emdiags%em_so2nat (:,:,:),                                &
           row_length, rows,model_levels, 0,0,0,0, at_extremity,     &
           stlist(1,stindex(1,item,section,im_index)), len_stlist,   &
           stash_levels, num_stash_levels+1,                         &
           atmos_im, section, item, icode, cmessage)

  IF (icode >  0) THEN
    errcode = section*1000+item
    CALL ereport ('UKCA_EMISS_DIAGS', errcode, cmessage)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ukca_emiss_diags

END MODULE ukca_emiss_diags_mod
