! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Sets the output dump fixed header.

MODULE Rcf_Setup_FixHd_Mod

IMPLICIT NONE

!  Subroutine Rcf_Setup_FixHd - sets the output dump fixed header.
!
! Description:
!   Sets up the fixed header for the output dump based on namelist
!   information.
!
! Method:
!   UMDP F3 defines the fixed header.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SETUP_FIXHD_MOD'

CONTAINS

SUBROUTINE Rcf_Setup_FixHd( Hdr_In, Hdr_Out )

USE Rcf_Grid_Type_Mod, ONLY: &
    Output_Grid

USE Rcf_HeadAddress_Mod     ! Huge amounts of this...

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_UMhead_Mod, ONLY: &
    Um_header_type,    &
    LenFixHd

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Oper

USE rcf_nlist_recon_technical_mod, ONLY: &
    reset_data_time

USE nlstcall_mod, ONLY: &
    lcal360

USE model_domain_mod, ONLY: &
    model_type,             &
    mt_bi_cyclic_lam,       &
    mt_cyclic_lam,          &
    l_cartesian,            &
    output_grid_stagger

USE lam_config_inputs_mod, ONLY: &
    delta_lon, delta_lat, frstlata, polelata

USE Rcf_headers_Mod, ONLY: &
    FixHd

USE io_configuration_mod, ONLY: &
    io_field_padding

USE missing_data_mod, ONLY: imdi

USE errormessagelength_mod, ONLY: errormessagelength

USE um_version_mod, ONLY: um_version_int

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (Um_header_type), TARGET      :: Hdr_In
TYPE (Um_header_type), TARGET      :: Hdr_Out

! Local parameters:
INTEGER, PARAMETER :: no_version = -9999 ! Value returned by function
                                         ! get_umversion if environment
                                         ! variable VN not set
! Local Variables
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_SETUP_FIXHD'
CHARACTER (LEN=errormessagelength) :: Cmessage
INTEGER, POINTER             :: FixHd_In( : )
INTEGER, POINTER             :: FixHd_Out( : )
INTEGER                      :: llll
INTEGER                      :: i
INTEGER                      :: ipos
INTEGER                      :: ErrorStatus
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

EXTERNAL Date_Time

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Setup pointer to fixhd arrays
FixHd_In  => Hdr_In % FixHd
FixHd_Out => Hdr_Out % FixHd

!-----------------------------------------------------------------
! First Setup all parts of fixed header to IMDI
!-----------------------------------------------------------------
FixHd_Out( : ) = imdi

!-----------------------------------------------------------------
! Data Set Format Version Number
!-----------------------------------------------------------------
FixHd_Out( FH_Version ) = FH_Version_Value
!-----------------------------------------------------------------
! Times
!-----------------------------------------------------------------

! Current Time
CALL Date_Time(FixHd_Out( FH_CTYear  ), FixHd_Out( FH_CTMonth ), &
         FixHd_Out( FH_CTDay   ), FixHd_Out( FH_CTHour  ), &
         FixHd_Out( FH_CTMinute), FixHd_Out( FH_CTSecond) )

! Validity Time - copied from input.
FixHd_Out( FH_VTYear : FH_VTDayNo ) = FixHd_In( FH_VTYear : FH_VTDayNo )

! Initial Time - either copied from input or from validity (if RESET)
IF (reset_data_time) THEN
  FixHd_Out( FH_DTYear : FH_DTDayNo ) =                       &
                   FixHd_Out( FH_VTYear : FH_VTDayNo )
  IF (PrintStatus >= PrStatus_Oper) THEN
    WRITE(umMessage,'('' ANALYSIS TIME RESET '')')
    CALL umPrint(umMessage,src='rcf_setup_fixhd_mod')
    WRITE(umMessage,'('' Analysis:'',7I9)') (FixHd_Out(i), i = FH_DTYear, &
                                                  FH_DTDayNo)
    CALL umPrint(umMessage,src='rcf_setup_fixhd_mod')
  END IF
ELSE
  FixHd_Out( FH_DTYear : FH_DTDayNo ) =                       &
                   FixHd_In( FH_DTYear : FH_DTDayNo )
END IF

!-------------------------------------------------------------------
! Other Fixhd stuff - non-addressing
!-------------------------------------------------------------------

! Atmos or Ocean?
FixHd_Out( FH_SubModel ) = FixHd_In( FH_Submodel )

! Vertical Co-ord type
FixHd_Out( FH_VertCoord ) = FH_VertCoord_CP

! Grid Staggering 
FixHd_Out( FH_GridStagger ) = output_grid_stagger

! Horizontal grid type:
IF (l_cartesian) THEN
  ! A LAM is assumed. The FLH header cannot distinguish between (bi-)cyclic
  ! and non-cyclic Cartesian domains, so for this we must refer to the
  ! input namelist. We also reuse the concept of 'wrapping LAMs' to denote
  ! a cyclic domain, although the header does not currently distinguish
  ! between (E-W) cyclic and bi-cyclic domains.
  IF (model_type == mt_cyclic_lam .OR. model_type == mt_bi_cyclic_lam) THEN
    FixHd_Out( FH_HorizGrid ) = FH_HorizGrid_LamWrap
  ELSE
    FixHd_Out( FH_HorizGrid ) = FH_HorizGrid_LamNoWrap
  END IF

ELSE

  FixHd_Out( FH_HorizGrid ) = Rcf_Setup_FixHd_HorizGrid()

END IF

! Projection number:
IF (l_cartesian) THEN
  FixHd_Out( FH_CoordSystem ) = Fh_CoordSystem_Cartesian
ELSE
  FixHd_Out( FH_CoordSystem ) = FH_CoordSystem_RegLonLat
END IF

! Instantaneous Data
FixHd_Out( FH_Dataset      ) = FH_Dataset_InstDump
FixHd_Out( FH_RunId        ) = FixHd_In( FH_RunId        )
FixHd_Out( FH_ExptNo       ) = FixHd_In( FH_ExptNo       )

! Calendar - 360 or 365 day
IF (LCal360) THEN
  FixHd_Out( FH_CalendarType ) = 2     ! 360
ELSE
  FixHd_Out( FH_CalendarType ) = 1     ! 365
END IF

FixHd_Out( FH_AncilDataId ) = FixHd_In( FH_AncilDataId )

! Dump model version number
FixHd_Out( FH_ModelVersion ) = um_version_int

!------------------------------------------------------------------
! Addressing information
!------------------------------------------------------------------
! Integer constants
ipos= LenFixHd+1
FixHd_Out( FH_IntCStart )=ipos
FixHd_Out( FH_IntCSize  )=hdr_out % LenIntC

! Real constants
ipos=ipos+FixHd_Out( FH_IntCSize )
FixHd_Out( FH_RealCStart )=ipos
FixHd_Out( FH_RealCSize  )=hdr_out % LenRealC

! Level dependent constants
ipos=ipos+FixHd_Out( FH_RealCSize )
llll=hdr_out % Len1LevDepC*hdr_out % Len2LevDepC
IF (llll == 0) THEN
  FixHd_Out( FH_LevDepCStart )=imdi
  FixHd_Out( FH_LevDepCSize1 )=imdi
  FixHd_Out( FH_LevDepCSize2 )=imdi
ELSE
  FixHd_Out( FH_LevDepCStart )=ipos
  FixHd_Out( FH_LevDepCSize1 )=hdr_out % Len1LevDepC
  FixHd_Out( FH_LevDepCSize2 )=hdr_out % Len2LevDepC
END IF

! Row dependent constants
ipos=ipos+llll
llll=hdr_out % Len2RowDepC*hdr_out % Len1RowDepC
IF (llll == 0) THEN
  FixHd_Out( FH_RowDepCStart )=imdi
  FixHd_Out( FH_RowDepCSize1 )=imdi
  FixHd_Out( FH_RowDepCSize2 )=imdi
ELSE
  FixHd_Out( FH_RowDepCStart )=ipos
  FixHd_Out( FH_RowDepCSize1 )=hdr_out % Len1RowDepC
  FixHd_Out( FH_RowDepCSize2 )=hdr_out % Len2RowDepC
END IF

! Column dependent constants
ipos=ipos+llll
llll=hdr_out % Len1ColDepC*hdr_out % Len2ColDepC
IF (llll == 0) THEN
  FixHd_Out( FH_ColDepCStart )=imdi
  FixHd_Out( FH_ColDepCSize1 )=imdi
  FixHd_Out( FH_ColDepCSize2 )=imdi
ELSE
  FixHd_Out( FH_ColDepCStart )=ipos
  FixHd_Out( FH_ColDepCSize1 )=hdr_out % Len1ColDepC
  FixHd_Out( FH_ColDepCSize2 )=hdr_out % Len2ColDepC
END IF

! Fields of constants
ipos=ipos+llll
llll=hdr_out % Len1FldsOfC*hdr_out % Len2FldsOfC
IF (llll == 0) THEN
  FixHd_Out( FH_FldsOfCStart )=imdi
  FixHd_Out( FH_FldsOfCSize1 )=imdi
  FixHd_Out( FH_FldsOfCSize2 )=imdi
ELSE
  FixHd_Out( FH_FldsOfCStart )=ipos
  FixHd_Out( FH_FldsOfCSize1 )=hdr_out % Len1FldsOfC
  FixHd_Out( FH_FldsOfCSize2 )=hdr_out % Len2FldsOfC
END IF

! Extra constants
ipos=ipos+llll
llll=hdr_out % LenExtraC
IF (llll == 0) THEN
  FixHd_Out( FH_ExtraCStart )=imdi
  FixHd_Out( FH_ExtraCSize )=imdi
ELSE
  FixHd_Out( FH_ExtraCStart )=ipos
  FixHd_Out( FH_ExtraCSize )=hdr_out % LenExtraC
END IF

! Temp history record
ipos=ipos+llll
llll=hdr_out % LenHistFile
IF (llll == 0) THEN
  FixHd_Out( FH_HistStart )=imdi
  FixHd_Out( FH_HistSize )=imdi
ELSE
  FixHd_Out( FH_HistStart )=ipos
  FixHd_Out( FH_HistSize )=hdr_out % LenHistFile
END IF

! Compress index 1
ipos=ipos+llll
llll=hdr_out % LenCompFldI1
IF (llll == 0) THEN
  FixHd_Out( FH_CompFldI1Start )=imdi
  FixHd_Out( FH_CompFldI1Size )=imdi
ELSE
  FixHd_Out( FH_CompFldI1Start )=ipos
  FixHd_Out( FH_CompFldI1Size )=hdr_out % LenCompFldI1
END IF
! Compress index 2
ipos=ipos+llll
llll=hdr_out % LenCompFldI2
IF (llll == 0) THEN
  FixHd_Out( FH_CompFldI2Start )=imdi
  FixHd_Out( FH_CompFldI2Size )=imdi
ELSE
  FixHd_Out( FH_CompFldI2Start )=ipos
  FixHd_Out( FH_CompFldI2Size )=hdr_out % LenCompFldI2
END IF
! Compress index 3
ipos=ipos+llll
llll=hdr_out % LenCompFldI3
IF (llll == 0) THEN
  FixHd_Out( FH_CompFldI3Start )=imdi
  FixHd_Out( FH_CompFldI3Size )=imdi
ELSE
  FixHd_Out( FH_CompFldI3Start )=ipos
  FixHd_Out( FH_CompFldI3Size )=hdr_out % LenCompFldI3
END IF

! Lookup
ipos=ipos+llll
llll=hdr_out % Len1Lookup*hdr_out % Len2Lookup
IF (llll == 0) THEN
  FixHd_Out( FH_LookupStart )=imdi
  FixHd_Out( FH_LookupSize1 )=imdi
  FixHd_Out( FH_LookupSize2 )=imdi
ELSE
  FixHd_Out( FH_LookupStart )=ipos
  FixHd_Out( FH_LookupSize1 )=hdr_out % Len1Lookup
  FixHd_Out( FH_LookupSize2 )=hdr_out % Len2Lookup
  !      For dumps, add no of prognostic fields to FixHd_Out
  !      In Reconfiguration hdr_out % Len2Lookup = No of prognostic fields
  IF (FixHd_Out( FH_DataSet ) == FH_DataSet_InstDump ) THEN
    FixHd_Out( FH_NumProgFields ) = hdr_out % Len2Lookup
  END IF
END IF

! Model data
ipos=ipos+llll
llll=hdr_out % LenData
IF (llll == 0) THEN
  FixHd_Out( FH_DataStart )=imdi
  FixHd_Out( FH_DataSize )=imdi
ELSE
          ! make sure the data starts on a sector bndry
  fixhd_Out( FH_DataStart )=    &
       ((ipos+io_field_padding-1)/io_field_padding)*io_field_padding+1
  FixHd_Out( FH_DataSize )=hdr_out % LenData
END IF

!-----------------------------------------------------------------
! Finally, do we need to overwrite with namelist/module values?
!-----------------------------------------------------------------
DO i = 1, LenFixHd
  IF ( FixHd(i) /= imdi ) THEN
    IF ( PrintStatus >= PrStatus_Oper ) THEN
      WRITE(umMessage,*) 'FixHd(',i,') has been reset from ', FixHd_out(i), &
                  ' to ', FixHd(i)
      CALL umPrint(umMessage,src='rcf_setup_fixhd_mod')
    END IF

    FixHd_Out(i) = FixHd(i)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

CONTAINS

INTEGER FUNCTION Rcf_Setup_FixHd_HorizGrid( ) RESULT(horiz_grid_type)

IMPLICIT NONE

CHARACTER (LEN=*), PARAMETER      :: routinename='RCF_SETUP_FIXHD_HORIZGRID'
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( Output_Grid % Global ) THEN
  horiz_grid_type= FH_HorizGrid_Global
ELSE IF ( (frstlata > 89.99 ) .AND.                                     &
         (delta_lon *  Output_Grid % glob_p_row_length > 359.99) .AND. &
         (delta_lat * (Output_Grid % glob_p_rows-1) > 89.99) .AND.     &
         (delta_lat * (Output_Grid % glob_p_rows-1) < 90.01)) THEN

  ! Warn user that this option has not been used for several years
  ! and output should be checked carefully.
  Cmessage = 'N. hemisphere-only grid - check carefully'
  ErrorStatus = -10
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  ! These are read in from lam_config and reset here.
  delta_lon = 360.0/REAL(Output_Grid % glob_p_row_length)
  delta_lat = 90.0/REAL(Output_Grid % glob_p_rows-1)
  horiz_grid_type=FH_HorizGrid_NH

ELSE IF ((frstlata < -89.99) .AND.                                      &
        (delta_lon *  Output_Grid % glob_p_row_length > 359.99) .AND.  &
        (delta_lat * (Output_Grid % glob_p_rows-1) > 89.99) .AND.      &
        (delta_lat * (Output_Grid % glob_p_rows-1) < 90.01)) THEN


  ! Warn user that this option has not been used for several years
  ! and output should be checked carefully.

  Cmessage = 'S. hemisphere-only grid - check carefully'
  ErrorStatus = -20
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  ! These are read in from lam_config and reset here.
  delta_lon = 360.0/REAL(Output_Grid % glob_p_row_length)
  delta_lat = 90.0/REAL(Output_Grid % glob_p_rows-1)
  horiz_grid_type= FH_HorizGrid_SH


  ! The Fixhd has no means of identifying whether a grid is a cyclic
  ! or bicyclic lam. A wrapping LAM could be E-W N-S or both and and code
  ! treats wrapping as seperate to cyclic, though the headers do not.

  ! It is not clear whether this potential functionality has any use outside
  ! the Cartesian case, which is coded for separately.


ELSE IF (Output_Grid % glob_p_row_length * delta_lon > 359.99) THEN
  IF (polelata > 89.99) THEN
    horiz_grid_type= FH_HorizGrid_LamWrap
  ELSE
    horiz_grid_type= FH_HorizGrid_LamWrapEq
  END IF
ELSE
  IF (polelata > 89.99) THEN
    horiz_grid_type= FH_HorizGrid_LamNoWrap
  ELSE
    horiz_grid_type= FH_HorizGrid_LamNoWrapEq
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION Rcf_Setup_FixHd_HorizGrid

END SUBROUTINE Rcf_Setup_FixHd
END MODULE Rcf_Setup_FixHd_Mod
