! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Shared code to open land sea mask ancillary file and setup header
!
MODULE Rcf_Open_LSM_Ancil_mod

IMPLICIT NONE

! Description:
!   Opens the land sea mask ancillary specified in the items namelist.
!
! Method:
!   Determine filename, open and check headers
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_OPEN_LSM_ANCIL_MOD'

CONTAINS

! Subroutine interface:
SUBROUTINE Rcf_Open_LSM_Ancil( Hdr_Anc, Output_Grid, lsm_filename )

USE ancil_mod, ONLY:              &
    ancil_files,                  &
    find_ancil_file_by_stash

USE Ereport_Mod, ONLY:            &
    Ereport

USE errormessagelength_mod, ONLY: &
    errormessagelength

USE file_manager, ONLY:           &
    assign_file_unit              

USE filenamelength_mod, ONLY:     &
    filenamelength

USE io, ONLY:                     &
    file_open

USE lookup_addresses, ONLY:       &
    lbcode, lbproc

USE Rcf_Grid_Type_Mod, ONLY:      &
    grid_type

USE Rcf_Items_Mod, ONLY:          &
    Num_items

USE Rcf_ReadUMhdr_mod, ONLY:      &
    Rcf_ReadUMhdr

USE Rcf_UMhead_Mod, ONLY:         &
    um_header_type

USE UM_ParCore, ONLY:             &
    mype

USE um_stashcode_mod, ONLY:       &
    stashcode_lsm

USE umPrintMgr, ONLY:             &
    umPrint,                      &
    umMessage,                    &
    PrintStatus,                  &
    PrStatus_Normal

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (grid_type), TARGET,       INTENT(IN)  :: Output_Grid
TYPE (um_header_type),TARGET,   INTENT(OUT) :: hdr_anc ! Header for ancillary
CHARACTER (LEN=filenamelength), INTENT(OUT) :: lsm_filename

! Local Data
INTEGER                      :: i           ! Looper
INTEGER                      :: errorstatus
INTEGER                      :: lsm_file_number
INTEGER, PARAMETER           :: filenameprovided=1
INTEGER, PARAMETER           :: lbcode_reg_lat_lon     = 1  
INTEGER, PARAMETER           :: lbcode_reg_lat_lon_rot = 101
                                            ! Valid values of lbcode
INTEGER, PARAMETER           :: lbproc_no_processing   = 0
                                            ! Instantaneous fields
CHARACTER (LEN=*), PARAMETER ::                                   &
                   RoutineName = 'RCF_OPEN_LSM_ANCIL'
CHARACTER (LEN=errormessagelength) :: cmessage
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------
!   First, determine the name of the LSM ancillary
!----------------------------------------------------------------

lsm_file_number = find_ancil_file_by_stash(stashcode_lsm)
lsm_filename = ancil_files(lsm_file_number)%filename

!---------------------------------------------------------------
!   Then open the ancillary file
!---------------------------------------------------------------

  IF ( PrintStatus >= PrStatus_Normal .AND. mype == 0 ) THEN
    WRITE(umMessage,'(A)')                                        &
         'Reading in Ancillary Land/Sea Mask'
    CALL umPrint(umMessage,src=RoutineName)
  END IF

  CALL assign_file_unit(lsm_filename, Hdr_Anc % UnitNum, handler="portio")

  CALL File_Open( Hdr_Anc % UnitNum, lsm_filename,                &
       LEN_TRIM(lsm_filename),0,filenameprovided,errorstatus)
  IF ( errorstatus /= 0 ) THEN
    cmessage = 'Error Opening Land/Sea Mask Ancillary'
    CALL Ereport( RoutineName, errorstatus, cmessage )
  END IF
  
  CALL Rcf_ReadUMhdr( Hdr_Anc )

!---------------------------------------------------------------
!   Check header is consistent with expected resolution.
!---------------------------------------------------------------

  IF (hdr_anc % intc(6) /= Output_Grid % glob_p_row_length) THEN
    errorstatus = 41
    cmessage = 'Incorrect E-W resolution for Land/Sea Mask Ancillary'
    CALL Ereport( RoutineName, errorstatus, cmessage )
  END IF

  IF (hdr_anc % intc(7) /= Output_Grid % glob_p_rows) THEN
    errorstatus = 42
    cmessage = 'Incorrect N-S resolution for Land/Sea Mask Ancillary'
    CALL Ereport( RoutineName, errorstatus, cmessage )
  END IF

!---------------------------------------------------------------
!   Setup the ancillary field before reading in from file.
!
!   Note that some variable resolution ancillary files have the 
!   second digit of lbcode set to 1. The current interpretation 
!   of UMDP F3 is that this is incorrect and these fields would 
!   fail the "is_valid_field()" rcf_setup_field. For this 
!   reason, we need to adjust the headers in that case.
!
!   Similarly, some land sea masks contain non-zero values of
!   lbproc (i.e. non-instantaneous fields). 
!   Again, in that case, reset these to zero.
!
!---------------------------------------------------------------
  
  IF ( ( Hdr_Anc % Lookup(lbcode, 1)                              &
                   == lbcode_reg_lat_lon + 10 ) .OR.              &
       ( Hdr_Anc % Lookup(lbcode, 1)                              &
                   == lbcode_reg_lat_lon_rot + 10 ) ) THEN
    
    errorstatus = -100
    WRITE(cmessage,'(A,I5,A,I5)')                                 &
         'Resetting lbcode in land mask ancil from ',             &
         Hdr_Anc % Lookup(lbcode, 1), ' to ',                     &
         Hdr_Anc % Lookup(lbcode, 1) - 10
    CALL Ereport( RoutineName, errorstatus, cmessage )
    Hdr_Anc % Lookup(lbcode, 1) = Hdr_Anc % Lookup(lbcode, 1) - 10

  END IF


  IF ( Hdr_Anc % Lookup(lbproc, 1)                                &
                 /= lbproc_no_processing ) THEN

    errorstatus = -200
    WRITE(cmessage,'(A,I5,A,I5)')                                 &
         'Resetting lbproc in land mask ancil from ',             &
         Hdr_Anc % Lookup(lbproc, 1),' to ',lbproc_no_processing
    CALL Ereport( RoutineName, errorstatus, cmessage )
    Hdr_Anc % Lookup(lbproc, 1) = lbproc_no_processing

  END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Open_LSM_Ancil
END MODULE Rcf_Open_LSM_Ancil_mod
