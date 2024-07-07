! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description: Ouput UKCA diagnostics on pressure levels
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: UKCA
!
!  Code description:
!    Language: Fortran
!
! ---------------------------------------------------------------------
MODULE ukca_plev_diags_mod
! ----------------------------------------------------------------------
! Processing Diagnostics on pressure levels: section 34 to section 51
! and section 50 to section 52
! ----------------------------------------------------------------------

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_PLEV_DIAGS_MOD'

CONTAINS

SUBROUTINE ukca_plev_diags(l_ukca_chem, l_ukca_chem_plev, l_ukca_asad_plev, &
    row_length, rows, model_levels, all_ntp, exner_theta_levels,            &
    tracer_ukca_um, STASHwork50)

USE atm_fields_bounds_mod,  ONLY: tdims_s
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE nlsizes_namelist_mod,   ONLY: tr_ukca
! Things from the STASH system we require 
USE stash_array_mod, ONLY: si, sf, stindex, stlist, stash_levels, &
    stash_maxlen, num_stash_levels
USE submodel_mod, ONLY: atmos_sm, atmos_im

USE ukca_calc_plev_diag_mod, ONLY: ukca_calc_plev_diag
! Get section for UKCA prognostics and  diagnostics
USE ukca_d1_defs, ONLY: ukca_sect, ukca_diag_sect, ukca_s34_plev, &
    ukca_s50_plev
USE ukca_ntp_mod, ONLY: ntp_type, stash2ntpindex
USE ukca_option_mod, ONLY: tr_ukca_a
USE ukca_tracer_stash,    ONLY: a_max_ukcavars
USE umPrintMgr, ONLY: umMessage, umPrint
USE version_mod, ONLY: nitemp

! Dr Hook
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Input arguments
! Non transported prognostics
TYPE(ntp_type), INTENT(IN) ::  all_ntp(:)
! x, y and z sizes of arrays
INTEGER, INTENT(IN) ::  row_length, rows, model_levels
REAL, INTENT(IN) ::  exner_theta_levels                                 &
                           (tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            1:tdims_s%k_end)
REAL, INTENT(IN) ::  tracer_ukca_um                             &
              (tdims_s%i_start:tdims_s%i_end,                   &
               tdims_s%j_start:tdims_s%j_end,                   &
               tdims_s%k_start:tdims_s%k_end,1:tr_ukca)
REAL, INTENT(IN) ::  STASHwork50(:)

! Logicals
LOGICAL, INTENT(IN) ::  l_ukca_chem ! Is UKCA chemistry on?
LOGICAL, INTENT(IN) ::  l_ukca_chem_plev  ! Are section 34 on p levels on?
LOGICAL, INTENT(IN) ::  l_ukca_asad_plev  ! Are section 50 on p levels on?

! Local variables 
! Allocatable arrays to hold STASHwork data
! STASHwork51 and STASHwork52
REAL, ALLOCATABLE :: STASHwork51(:)
REAL, ALLOCATABLE :: STASHwork52(:)

! To hold data to write to STASH
REAL, ALLOCATABLE :: diag_temp(:,:,:)   

! Error return message
CHARACTER(LEN=errormessagelength)   :: cmessage

INTEGER    :: sect_in, sect_out ! section -IN/OUT
INTEGER    :: itm_in, itm_out   ! item -IN/OUT
INTEGER    :: n_plevs           ! number of output pressure levels
INTEGER    :: pt_pdiag, pt_diag ! pointers to diag in STASHwork
INTEGER    :: isl, ni           ! Stash indices
INTEGER    :: ntpi              ! Index of NTP
INTEGER    :: errcode, icode    ! error codes 

! Section and item number of Heavyside function
INTEGER, PARAMETER :: itm_hvside = 999
INTEGER, PARAMETER :: sec_hvside = 51

! pointers to Heaviside in STASHwork
INTEGER    :: pt_hvside

INTEGER    :: im_index          ! internal model index
INTEGER    :: um_tr_index       ! index in um tracer array

LOGICAL, SAVE :: l_chk_hvside = .FALSE. ! Check if Heavyside diag is ON

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_PLEV_DIAGS'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!- End of Header ---------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set up allocatable arrays needed for STASH. Need section 51 
! allocatables to hold heaviside diagnostic for section 52 ones 
IF (l_ukca_chem_plev .AND. .NOT. l_ukca_asad_plev) THEN
  ALLOCATE (STASHwork51(STASH_maxlen(51,atmos_im)))
  STASHwork51=0.0
END IF

IF (l_ukca_asad_plev) THEN
  ALLOCATE (STASHwork51(STASH_maxlen(51,atmos_im)))
  STASHwork51=0.0
  ALLOCATE (STASHwork52(STASH_maxlen(52,atmos_im)))
  STASHwork52=0.0
END IF

im_index = 1
! ----------------------------------------------------------------------
! Interpolate section 34 on model levels, sending results to
! the section 51
! ----------------------------------------------------------------------
IF (l_ukca_chem_plev) THEN
  sect_in = ukca_sect
  sect_out = ukca_s34_plev  ! section 51
  l_chk_hvside = .TRUE.

  DO itm_out = 1, nitemp        ! Loop over all items in section 51
   
    ! If there is a STASH request for this item we need to do
    ! vertical interpolation and add results to stashwork
    IF (sf(itm_out,sect_out) .AND. itm_out /= itm_hvside) THEN

      ! The item number used to get the data is the same as the 
      ! item number being output (in a different section)
      itm_in = itm_out          ! the same item number

! Heavyside diag required
      IF (l_chk_hvside .AND. .NOT. sf(itm_hvside,sec_hvside)) THEN
        errcode = sec_hvside*1000+itm_hvside
        CALL ereport(ModuleName//':'//RoutineName,errcode, &
              'Heavyside diagnostic not requested with P level diags')
      END IF
      
      ! Pointer to Heavyside diag in STASHwork
      pt_hvside = si(itm_hvside,sec_hvside,im_index)

      ! Pointer to item in STASHwork
      pt_pdiag = si(itm_out,sect_out,im_index)

      ! Stash indices
      isl = stindex(1,itm_out,sect_out,im_index) 
      ni = -stlist(10,isl)
      
      ! set number of output pressure levels
      n_plevs = stash_levels(1,ni)

      ! take tracers out of array directly
      IF (itm_out <= a_max_ukcavars) THEN
    
          ! is this tracer on?
          IF (tr_ukca_a(itm_out)) THEN
              um_tr_index = COUNT(tr_ukca_a(1:itm_out))
              ALLOCATE(diag_temp(1:row_length,1:rows,1:model_levels))
              diag_temp = tracer_ukca_um(1:row_length,1:rows,1:model_levels, &
                                         um_tr_index)
          END IF
        
      ELSE
          ! This will be a non transported prognostic
          
          ! Look up the index. This routine has an ereport
          ! if it fails to find the index
          ntpi = stash2ntpindex(all_ntp, sect_in, itm_in)
          ALLOCATE(diag_temp(1:row_length,1:rows,1:model_levels))
          diag_temp = all_ntp(ntpi)%data_3d(1:row_length,1:rows,1:model_levels)

      END IF 
      
      ! If the above block found data then the diag_temp array 
      ! is allocated. Use, then deallocate.
      IF (ALLOCATED(diag_temp) ) THEN
          CALL ukca_calc_plev_diag( num_stash_levels,                &
            stash_levels(:,ni), row_length, rows,                    &
            model_levels, n_plevs,                                   &
            diag_temp(1:row_length,1:rows,1:model_levels),           &
            exner_theta_levels(1:row_length,1:rows,1:model_levels),  &
            STASHwork51(pt_pdiag), STASHwork51(pt_hvside))

          DEALLOCATE(diag_temp)

      ! Otherwise we failed to find the data. Report an error
      ELSE

        WRITE(umMessage,'(A43,I7,A7,I7)')                      &
          'Item for P_LEVS not found IN: ',                    &
           sect_in*1000+itm_in,' OUT: ',sect_out*1000+itm_out
        CALL umPrint(umMessage,src=ModuleName//':'//RoutineName)
        cmessage = 'Item for P_LEVS not found'
        errcode = sect_out*1000+itm_out
        CALL ereport(ModuleName//':'//RoutineName, &
                     errcode, cmessage)
      END IF

      ! we only need to check the Heaviside diagnostic once
      l_chk_hvside = .FALSE.

    END IF  ! sf(itm_out) was requested
  END DO  ! loop for all items in sec 51

! Write results if only section 34 was required, otherwise STASHwork51
! will be written together with STASHwork52
  IF (.NOT. l_ukca_asad_plev) THEN

! DEPENDS ON: stash
    ! intialise values of cmessage and icode to prevent false errors
    cmessage=' '
    icode = 0
    CALL stash(atmos_sm,atmos_im,sect_out,STASHwork51,                &
      icode,cmessage)

    IF (icode >  0) THEN
      cmessage=" Error writing Plev diags to STASH:"//cmessage
      errcode=sect_out
      CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
    END IF
  END IF
END IF ! l_ukca_chem

! ----------------------------------------------------------------------
! 5.3.2 Interpolate section 50 on model levels, sending results to
! the section 52 on pressure levels
! ----------------------------------------------------------------------
IF (l_ukca_chem .AND. l_ukca_asad_plev) THEN

  sect_in = UKCA_diag_sect
  sect_out = UKCA_s50_plev  ! section 52
  l_chk_hvside = .TRUE.

  DO itm_out = 1, nitemp        ! Loop over all items in section 50

    IF (sf(itm_out,sect_out)) THEN
      itm_in = itm_out          ! the same item number

! Equivalent diag required
      IF ( itm_in /= itm_hvside .AND. .NOT. sf(itm_in,sect_in)) THEN
        errcode = sect_in*1000+itm_in
        CALL ereport(ModuleName//':'//RoutineName,errcode,         &
           'UKCA model_lev diagnostic required when P_lev version requested')
      END IF

! Heavyside diag required
      IF (l_chk_hvside .AND. .NOT. sf(itm_hvside,sec_hvside)) THEN
        errcode = sec_hvside*1000+itm_hvside
        CALL ereport(ModuleName//':'//RoutineName,errcode,  &
        'Heavyside diagnostic not requested with P level diags')
      END IF
      pt_hvside = si(itm_hvside,sec_hvside,im_index)

      n_plevs = 0
      pt_pdiag = si(itm_out,sect_out,im_index)   ! Pointer to item
                                                 ! in STASHwork

      isl = stindex(1,itm_out,sect_out,im_index) ! Stash indices
      ni = -stlist(10,isl)
      n_plevs = stash_levels(1,ni)        ! num output press levels
      pt_diag = si(itm_in,sect_in,im_index)
      ALLOCATE(diag_temp(row_length, rows,model_levels))
      diag_temp = RESHAPE(STASHwork50(pt_diag:pt_diag +          &
                         (row_length*rows*model_levels)-1),      &
                         (/ row_length, rows,model_levels /))

      CALL ukca_calc_plev_diag( num_stash_levels,                &
        stash_levels(:,ni), row_length, rows,                    &
        model_levels, n_plevs,                                   &
        diag_temp(1:row_length,1:rows,1:model_levels),           &
        exner_theta_levels(1:row_length,1:rows,1:model_levels),  &
        STASHwork52(pt_pdiag), STASHwork51(pt_hvside))

      IF ( ALLOCATED(diag_temp) ) DEALLOCATE(diag_temp)
      l_chk_hvside = .FALSE.

    END IF  ! sf(itm_out) was requested
  END DO  ! loop for all items in sec 50

  ! intialise values of cmessage and icode to prevent false errors
  cmessage=' '
  icode = 0
! DEPENDS ON: stash
  CALL stash(atmos_sm,atmos_im,sect_out,STASHwork52,             &
   icode,cmessage)

  IF (icode >  0) THEN
    cmessage=" Error writing Plev diags to STASH:"//cmessage
    errcode=sect_out
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF

! DEPENDS ON: stash
  CALL stash(atmos_sm,atmos_im,ukca_s34_plev,STASHwork51,        &
    icode,cmessage)

  IF (icode >  0) THEN
    cmessage=" Error writing Plev diags to STASH:"//cmessage
    errcode=ukca_s34_plev
    CALL ereport(ModuleName//':'//RoutineName,errcode,cmessage)
  END IF

END IF

! Deallocate STASHwork arrays
IF (ALLOCATED(STASHwork51)) DEALLOCATE(STASHwork51)
IF (ALLOCATED(STASHwork52)) DEALLOCATE(STASHwork52)
  
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ukca_plev_diags
END MODULE ukca_plev_diags_mod
