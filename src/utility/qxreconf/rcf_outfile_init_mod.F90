! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE Rcf_Outfile_Init_Mod

!  Subroutine Rcf_Files_Init - initialisation of output file
!
! Description:
!   Open output file and setup output header default sizes.
!
! Method:
!   Output dump opened with File_Open. Output header default sizes 
!   setup from input sizes.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_OUTFILE_INIT_MOD'

CONTAINS

SUBROUTINE Rcf_Outfile_Init( hdr_in, hdr_out )

USE Ereport_Mod, ONLY: Ereport

USE errormessagelength_mod, ONLY: errormessagelength

USE file_manager, ONLY: assign_file_unit

USE filenamelength_mod, ONLY: filenamelength

USE io

USE io_constants, ONLY: ioOpenReadWrite

USE missing_data_mod, ONLY: imdi

USE model_file, ONLY: model_file_open

USE nlcfiles_namelist_mod, ONLY: astart

USE nlsizes_namelist_mod, ONLY:                              &
    a_len_inthd,        a_len_realhd,      a_len2_levdepc,   &
    a_len2_rowdepc,     a_len2_coldepc,    a_len2_flddepc,   &
    a_len_extcnst

USE Rcf_Address_Mod, ONLY: Rcf_Address

USE Rcf_Grid_Type_Mod, ONLY: Output_Grid

USE Rcf_NRecon_Mod, ONLY: &
    DumpProgLevs,      &
    PrimDataLen

USE rcf_nlist_recon_technical_mod, ONLY:                             &
    interp_all_fields,                                               &
    len_dumphist,                                                    &
    select_output_fields,                                            &
    tstmsk_to_decide,                                                &
    defined_by_namelist

USE Rcf_UMhead_Mod, ONLY: um_header_type

USE Submodel_Mod, ONLY: &
    Internal_Model_List, &
    N_Internal_Model

USE UM_ParCore, ONLY: mype

USE umPrintMgr, ONLY: &
    umPrint,          &
    umMessage,        &
    PrintStatus,      &
    PrStatus_Normal

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (Um_Header_type), INTENT(IN)       :: hdr_in
TYPE (Um_Header_type), INTENT(INOUT)    :: hdr_out

! Local variables
CHARACTER (LEN=*), PARAMETER            :: RoutineName='RCF_OUTFILE_INIT'
CHARACTER (LEN=errormessagelength)      :: Cmessage
CHARACTER (LEN=filenamelength)          :: DumpName
INTEGER                                 :: ErrorStatus
INTEGER                                 :: err
INTEGER                                 :: i
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise errorstatus
errorstatus = 0

!-----------------------------------------------------------------
! Set addressing information for output dump.
!-----------------------------------------------------------------
! Define submodel and section/version configuration
! DEPENDS ON: setmodl
CALL setmodl(ErrorStatus,cmessage)
IF (ErrorStatus  /=  0) THEN
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! Output dump addressing:
 CALL rcf_address(hdr_in % lookup)

! Question: Do we need to zero IMDIs in the FLH?

!-------------------------------------------------------------------
! Setup Output Dump header sizes (from Input) and open file
!-------------------------------------------------------------------

! Ensure output dump location exists.
! Using dir_create would be preferable, but we don't know whether
! this is only the basename of the file or the full path to it.
IF (mype == 0 ) THEN
  CALL Shell('mkdir -p $(dirname '//TRIM(astart)//')', 20+LEN_TRIM(astart))
END IF

CALL assign_file_unit(astart, hdr_out % UnitNum, handler="portio")
CALL Model_File_Open ( hdr_out % UnitNum, astart, read_write=ioOpenReadWrite, &
                       error=err)

IF ( err /= 0 ) THEN
  Cmessage    = 'Failed to Open Output Dump'
  ErrorStatus = 20
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN

  CALL umPrint(  '',src='rcf_files_init_mod')
  WRITE(umMessage,'(2A)') 'Output dump : ',TRIM(astart)
  CALL umPrint(umMessage,src='rcf_files_init_mod')


END IF

! Initially set some values as in input dump (ensure sizes are >= 0 for
! allocating)
hdr_out % LenIntC      = MAX(hdr_in % LenIntC,0)
hdr_out % LenRealC     = MAX(hdr_in % LenRealC,0)
hdr_out % Len2LevDepC  = MAX(hdr_in % Len2LevDepC,0)
hdr_out % Len2RowDepC  = MAX(hdr_in % Len2RowDepC,0)
hdr_out % Len1Lookup   = MAX(hdr_in % Len1Lookup,0)
hdr_out % LenCompFldI1 = MAX(hdr_in % LenCompFldI1,0)
hdr_out % LenCompFldI2 = MAX(hdr_in % LenCompFldI2,0)
hdr_out % LenCompFldI3 = MAX(hdr_in % LenCompFldI3,0)
hdr_out % Len2ColDepC  = MAX(hdr_in % Len2ColDepC,0)
hdr_out % Len1FldsOfC  = MAX(hdr_in % Len1FldsOfC,0)
hdr_out % Len2FldsOfC  = MAX(hdr_in % Len2FldsOfC,0)
hdr_out % LenExtraC    = MAX(hdr_in % LenExtraC,0)
hdr_out % LenHistFile  = MAX(hdr_in % LenHistFile,0)

! Other values may be overwritten from Namelist RECON
IF ( a_len_inthd /= imdi ) THEN
  hdr_out % LenIntC = a_len_inthd
END IF
IF ( a_len_realhd /= imdi ) THEN
  hdr_out % LenRealC = a_len_realhd
END IF
IF ( a_len2_levdepc /= imdi ) THEN
  hdr_out % Len2LevDepC = a_len2_levdepc
END IF
IF ( a_len2_rowdepc /= imdi ) THEN
  hdr_out % Len2RowDepC = a_len2_rowdepc
END IF
IF ( a_len2_coldepc /= imdi ) THEN
  hdr_out % Len2ColDepC = a_len2_coldepc
END IF
IF ( a_len2_flddepc /= imdi ) THEN
  hdr_out % Len2FldsOfC = a_len2_flddepc
END IF
IF ( a_len_extcnst /= imdi ) THEN
  hdr_out % LenExtraC = a_len_extcnst
END IF
hdr_out % LenHistFile = len_dumphist

! Check that sizes of integer constants are big enough for vn5.0
!  dumps
IF (hdr_out % LenIntC < 46) THEN
  ErrorStatus = 30
  Cmessage = 'Length of Integer Constants needs to be at least '//&
               '46 for vn5.0 and higher dumps'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! Other values can be deduced from Output_Grid
! Row length is always same between P and U - even for LAMs where we need to
! keep similar structure to global.
hdr_out % Len1ColDepC  = Output_Grid % glob_p_row_length
! Rows are different between Endgame and ND due to V rows being greater than P.
hdr_out % Len1RowDepC  = MAX(Output_Grid % glob_p_rows, &
                             Output_Grid % glob_v_rows)

hdr_out % Len1LevDepC  = Output_Grid % model_levels + 1

! Others have to be deduced from STASH generated addressing.
hdr_out % Len2Lookup   = 0
hdr_out % LenData      = 0
DO i = 1, N_Internal_Model
  hdr_out % Len2Lookup = hdr_out % Len2Lookup +                    &
                         DumpProgLevs( Internal_Model_List( i ) )
  hdr_out % LenData    = hdr_out % LenData +                       &
                         PrimDataLen( Internal_Model_List( i ) )
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Outfile_Init
END MODULE Rcf_Outfile_Init_Mod
