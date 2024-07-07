! *****************************COPYRIGHT*******************************
!
! (c) [University of Oxford] [2011]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
!  Description:  Module to define type ukca_cdnc_struct, the structure
!                used by UKCA_CDNC, with intialisation and read routines
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components provided by The University of Cambridge,
!  University of Leeds, University of Oxford, and The Met Office.
!  See: www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
MODULE ukca_cdnc_mod

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr
USE errormessagelength_mod, ONLY: errormessagelength
USE submodel_mod, ONLY: submodel_for_sm, atmos_im
USE d1_array_mod, ONLY: d1_section, d1_item, d1_address, d1_length,        &
                        d1_no_levels, d1_halo_type, d1, d1_addr, no_obj_d1
IMPLICIT NONE
PUBLIC
SAVE

! Main structure holding all the variables needed
! to retrieve Cloud Droplet Number Concentration
! from D1, as calculated by UKCA_ACTIVATE

TYPE ukca_cdnc_struct

  ! Stash code, D1 address, number of levels and total length
  ! of the D1 fields corresponding to CDNC
  INTEGER :: stashcode_cdnc
  INTEGER :: d1_address_cdnc
  INTEGER :: d1_nlevs_cdnc
  INTEGER :: d1_length_cdnc

  ! Stash code, D1 address, number of levels and total length
  ! of the D1 fields corresponding to CDNC3
  INTEGER :: stashcode_cdnc3
  INTEGER :: d1_address_cdnc3
  INTEGER :: d1_nlevs_cdnc3
  INTEGER :: d1_length_cdnc3

  ! Cloud droplet number concentration (m^-3)
  REAL, ALLOCATABLE :: cdnc(:, :, :)
  ! <Cloud droplet number concentration^-1/3>
  REAL, ALLOCATABLE :: cdnc3(:, :, :)

END TYPE ukca_cdnc_struct

! Dimensions of CDNC arrays - will be either (row_length, rows,
! model_levels) or (1,1,1) depending on whether aerosol indirect
! effects are activated or not
INTEGER :: cdnc_dim1
INTEGER :: cdnc_dim2
INTEGER :: cdnc_dim3

! Structure for UKCA CDNC/radiation/ppn interaction
TYPE (ukca_cdnc_struct) :: ukca_cdnc

! The following prognostics are expected in section 34.
INTEGER, PARAMETER :: stashc_cdnc  = 968
INTEGER, PARAMETER :: stashc_cdnc3 = 967

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_CDNC_MOD'

CONTAINS

! ######################################################################

SUBROUTINE ukca_cdnc_init(ierr,cmessage)

IMPLICIT NONE

INTEGER, INTENT(INOUT)            :: ierr      ! Error indicator 
                                               !(0 is OK, >0 error)
CHARACTER(LEN=errormessagelength), INTENT(INOUT) :: cmessage  ! Error message

! Local variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CDNC_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise STASH code for CDNC prognostics


ierr = 0
ukca_cdnc%stashcode_cdnc = stashc_cdnc
ukca_cdnc%stashcode_cdnc3 = stashc_cdnc3
IF (PrintStatus >= PrStatus_Diag) THEN
  WRITE(umMessage,'(A,I6)') 'ukca_cdnc%stashcode_cdnc: ',           &
                  ukca_cdnc%stashcode_cdnc
  CALL umPrint(umMessage,src='ukca_cdnc_mod')
  WRITE(umMessage,'(A,I6)') 'ukca_cdnc%stashcode_cdnc3: ',          &
                          ukca_cdnc%stashcode_cdnc3
  CALL umPrint(umMessage,src='ukca_cdnc_mod')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_cdnc_init

! ######################################################################

SUBROUTINE ukca_cdnc_get(                                         &
                           first_call,                            &
                           ierr,                                  &
                           cmessage)

USE UM_ParVars
USE UM_ParParams, ONLY: halo_type_no_halo
USE atm_fields_bounds_mod,  ONLY: tdims_s
USE ukca_d1_defs,  ONLY: ukca_sect

USE nlsizes_namelist_mod, ONLY: &
    len_tot, n_obj_d1_max, row_length, rows, tr_levels

IMPLICIT NONE

INTEGER, INTENT(OUT)            :: ierr        ! Error indicator 
                                               !(0 is OK, >0 error)
CHARACTER(LEN=errormessagelength), INTENT(OUT) :: cmessage    ! Error message
LOGICAL, INTENT(IN)             :: first_call  ! Indicates first call:
                                      !  complete setup has to be done

! Local variables
INTEGER :: i
INTEGER :: j
INTEGER :: i_obj

! Atmosphere submodel index
INTEGER :: m_atm_modl

! Logical for checking whether all prognostics have been found
LOGICAL :: l_missing

! Expected size
INTEGER :: buffer_size

! Local variable for levels, to handle the fact that
!  levels = tr_levels in NewDyn but tr_levels+1 in ENDGame runs
INTEGER :: tr_levs

! Temporary buffer
REAL, ALLOCATABLE :: temp_buf(:,:,:)

! Variables for tagged prognostics in D1.
INTEGER :: section       ! stash section
INTEGER :: item          ! stash item
INTEGER :: levs          ! No of levels
INTEGER :: LEN           ! length of field
INTEGER :: addr          ! address in D1
INTEGER :: halotyp       ! halo type

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CDNC_GET'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ierr = 0
m_atm_modl = submodel_for_sm(atmos_im)

! Expected number of levels for required fields in D1
tr_levs  = tdims_s%k_len
buffer_size = row_length * rows * tr_levs

! For first call, check that all needed prognostics are available
! at the expected dimensions. Retain their D1 address for later calls.

IF (first_call) THEN

  ! Initialise d1 address in information structure to "not found" (-1)
  ukca_cdnc%d1_address_cdnc = -1
  ukca_cdnc%d1_address_cdnc3 = -1

  ! Go through D1 and retain some information about the
  ! prognostics we need.
  DO i=1,no_obj_D1(m_atm_modl)
    section = d1_addr(D1_section,i,m_atm_modl)
    IF (section /= ukca_sect) CYCLE
    item    = d1_addr(D1_item,i,m_atm_modl)
    levs    = d1_addr(d1_no_levels,i,m_atm_modl)
    LEN     = d1_addr(d1_length,i,m_atm_modl)
    addr    = d1_addr(d1_address,i,m_atm_modl)
    halotyp = d1_addr(d1_halo_type,i,m_atm_modl)


    IF (item == ukca_cdnc%stashcode_cdnc) THEN
      IF (halotyp /= halo_type_no_halo) THEN
        ierr = 700
        cmessage = 'ukca_cdnc_get: non-zero halo for CDNC'
        IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, &
                                zhook_handle)
        RETURN
      END IF
      ukca_cdnc%d1_address_cdnc = addr
      ukca_cdnc%d1_nlevs_cdnc   = levs
      ukca_cdnc%d1_length_cdnc  = LEN
      IF (PrintStatus >= PrStatus_Diag) THEN
        WRITE(umMessage,'(A,I10)')'UKCA_CDNC: address in D1 ',addr
        CALL umPrint(umMessage,src='ukca_cdnc_mod')
      END IF
    END IF

    IF (item == ukca_cdnc%stashcode_cdnc3) THEN
      IF (halotyp /= halo_type_no_halo) THEN
        ierr = 701
        cmessage = 'ukca_cdnc_get: non-zero halo for CDNC3'
        IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, &
                                zhook_handle)
        RETURN
      END IF
      ukca_cdnc%d1_address_cdnc3 = addr
      ukca_cdnc%d1_nlevs_cdnc3   = levs
      ukca_cdnc%d1_length_cdnc3  = LEN
      IF (PrintStatus >= PrStatus_Diag) THEN
        WRITE(umMessage,'(A,I10)')'UKCA_CDNC3: address in D1 ',addr
        CALL umPrint(umMessage,src='ukca_cdnc_mod')
      END IF
    END IF

    IF (ukca_cdnc%d1_address_cdnc > 0 .AND.                       &
        ukca_cdnc%d1_address_cdnc3 > 0) EXIT

  END DO ! i_obj (D1 items)

  ! Check that all prognostics have been found.
  l_missing = .FALSE.

  IF (ukca_cdnc%d1_address_cdnc == -1) THEN
    WRITE(umMessage,'(A,I6,A)')                                 &
     'Prognostic ', ukca_cdnc%stashcode_cdnc,'not found in D1.'
    CALL umPrint(umMessage,src='ukca_cdnc_mod')
    l_missing = .TRUE.
  END IF

  IF (ukca_cdnc%d1_address_cdnc3 == -1) THEN
    WRITE(umMessage,'(A,I6,A)')                                &
     'Prognostic ', ukca_cdnc%stashcode_cdnc3,'not found in D1.'
    CALL umPrint(umMessage,src='ukca_cdnc_mod')
    l_missing = .TRUE.
  END IF

  IF (l_missing) THEN
    ierr = 702
    cmessage =                                                   &
     'ukca_cdnc_get: Prog(s) needed for UKCA are missing from D1'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  ! Check number of levels and total size of arrays in D1
  IF (ukca_cdnc%d1_nlevs_cdnc /= tr_levs) THEN
    ierr = 705
    WRITE(umMessage,'(A,I4,A,I4)')'Expecting ', tr_levs,         &
      ' levels, got ', ukca_cdnc%d1_nlevs_cdnc
    CALL umPrint(umMessage,src='ukca_cdnc_mod')
    cmessage =                                                   &
       'ukca_cdnc_get: Unexpected number of levels in D1 prog.'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  IF (ukca_cdnc%d1_length_cdnc /= buffer_size) THEN
    ierr = 706
    WRITE(umMessage, '(A,I10,A,I10)')'Expecting ',buffer_size,    &
      ' elements, got ',ukca_cdnc%d1_length_cdnc
    CALL umPrint(umMessage,src='ukca_cdnc_mod')
    cmessage =                                                    &
      'ukca_cdnc_get: Unexpected total size of D1 prog.'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  IF (ukca_cdnc%d1_nlevs_cdnc3 /= tr_levs) THEN
    ierr = 705
    WRITE(umMessage,'(A,I4,A,I4)')'Expecting ', tr_levs,         &
        ' levels, got ', ukca_cdnc%d1_nlevs_cdnc3
    CALL umPrint(umMessage,src='ukca_cdnc_mod')
    cmessage =                                                    &
       'ukca_cdnc_get: Unexpected number of levels in D1 prog.'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  IF (ukca_cdnc%d1_length_cdnc3 /= buffer_size) THEN
    ierr = 706
    WRITE(umMessage,'(A,I10,A,I10)')'Expecting ',buffer_size,     &
     ' elements, got ', ukca_cdnc%d1_length_cdnc3
    CALL umPrint(umMessage,src='ukca_cdnc_mod')
    cmessage =                                                    &
      'ukca_cdnc_get: Unexpected total size of D1 prog.'
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

END IF ! first_call

ALLOCATE(temp_buf(row_length,rows,tdims_s%k_start:tdims_s%k_end))
temp_buf(:,:,:) = 0.0

! Cloud droplet number concentration.

! Read in full vector from D1 (will be tr_levs=0:tr_levels for ENDGame), 
! but always use 1:tr_levels in any configuration
temp_buf(:,:,:)      =  RESHAPE(d1(ukca_cdnc%d1_address_cdnc :    &
                                ukca_cdnc%d1_address_cdnc +       &
                                ukca_cdnc%d1_length_cdnc - 1),    &
                               (/row_length, rows, tr_levs/))
ukca_cdnc%cdnc(:,:,:) = temp_buf(:,:,1:tr_levels)

! <Cloud droplet number concentration^-1/3>
temp_buf(:,:,:) = 0.0
temp_buf(:,:,:)      = RESHAPE(d1(ukca_cdnc%d1_address_cdnc3 :    &
                                 ukca_cdnc%d1_address_cdnc3 +     &
                                 ukca_cdnc%d1_length_cdnc3 - 1),  &
                               (/row_length, rows, tr_levs/))
ukca_cdnc%cdnc3(:,:,:) = temp_buf(:,:,1:tr_levels)

DEALLOCATE(temp_buf)

IF (PrintStatus >= PrStatus_Diag) THEN
  WRITE(umMessage,'(A)') 'CDNC_GET:'
  CALL umPrint(umMessage,src='ukca_cdnc_mod')
  WRITE(umMessage,'(A,I10)') 'ukca_cdnc%d1_address_cdnc',         &
                       ukca_cdnc%d1_address_cdnc
  CALL umPrint(umMessage,src='ukca_cdnc_mod')
  WRITE(umMessage,'(A,I10)') 'ukca_cdnc%d1_address_cdnc3',        &
              ukca_cdnc%d1_address_cdnc3
  CALL umPrint(umMessage,src='ukca_cdnc_mod')
  IF (ALLOCATED(ukca_cdnc%cdnc)) THEN  
    WRITE(umMessage,'(A)')'ukca_cdnc%cdnc is allocated'
    CALL umPrint(umMessage,src='ukca_cdnc_mod')
  END IF
  WRITE(umMessage,'(A,I10)') 'cdnc: ',SIZE(ukca_cdnc%cdnc)
  CALL umPrint(umMessage,src='ukca_cdnc_mod')
  WRITE(umMessage,'(A,2E12.4)')'cdnc: ',MAXVAL(ukca_cdnc%cdnc),   &
                      MINVAL(ukca_cdnc%cdnc)
  CALL umPrint(umMessage,src='ukca_cdnc_mod')
  WRITE(umMessage,'(A,I10)') 'cdnc3: ',SIZE(ukca_cdnc%cdnc3)
  CALL umPrint(umMessage,src='ukca_cdnc_mod')
  WRITE(umMessage,'(A,2E12.4)')'cdnc3: ',MAXVAL(ukca_cdnc%cdnc3), &
                       MINVAL(ukca_cdnc%cdnc3)
  CALL umPrint(umMessage,src='ukca_cdnc_mod')
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_cdnc_get

END MODULE ukca_cdnc_mod
