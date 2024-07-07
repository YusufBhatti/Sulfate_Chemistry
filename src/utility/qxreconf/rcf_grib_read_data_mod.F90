! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Read, from file, and Decode 1 field's worth of GRIB data

MODULE Rcf_Grib_Read_Data_Mod

IMPLICIT NONE

! SUBROUTINE Rcf_Grib_Read_Data : Read raw info from file
!
! Description:
!   This routine is used to handle the reading of raw GRIB encoded data
!   in. This data is then passed to DECODE to be decoded.
!
! Method:
!   Recieve position markers and Storage arrays (or Types) from
!   calling routine.
!   Use SETPOS8 to set position in file to start of GRIB record.
!   Read a block of Raw data.
!   Call DECODE to decode raw Data.
!   Store the decoded data
!   Use the decoded length to calculate the next start point and move on
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GRIB_READ_DATA_MOD'

CONTAINS
SUBROUTINE Rcf_Grib_Read_Data(Unit_Num ,Current,DATA,Len_Data,        &
                              Pos_in_File)

USE Rcf_GRIB_Block_Params_Mod, ONLY: &
    Grib_Record,     &       ! derived type for storing a GRIB header
    LenArrayMax,     &       ! Max size of buffer for data
    p_B4Undef

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Min,           &       ! =1 Minimum output
    PrStatus_Normal,        &       ! =2 Short informative output
    PrStatus_Oper,          &       ! =3 Full informative output
    PrStatus_Diag                   ! =4 Extra Diagnostic output

USE UM_ParCore, ONLY: &
    mype

USE EReport_Mod, ONLY:     &
    EReport                         ! Subroutine to generate error msgs

USE io

USE missing_data_mod, ONLY: rmdi 

USE errormessagelength_mod, ONLY: errormessagelength

USE ukmo_grib_mod, ONLY: ukmo_decode

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Subroutine arguments

!< Scalar arguments with intent(in):>
INTEGER, INTENT(IN)              :: Unit_Num  ! unit number file is on
INTEGER, INTENT(IN)              :: Len_Data  ! Length of data to read
INTEGER, INTENT(IN)              :: Pos_in_File ! position to start
                                                ! reading data from

!< Array  arguments with intent(InOut):>
TYPE (Grib_Record),POINTER       :: Current        !\ Pointer to current
                                                   !/ grib record
!< Array  arguments with intent(out):>
REAL, INTENT(OUT)                :: DATA(Len_Data) ! Array deocded field
                                                   ! is returned in.

! Local constants
! This was usedfor debugging. It allowed all non debugging messages to
! be turned off.
INTEGER , PARAMETER              :: PrStatus_Debug = 0 ! debug

! Local variables

CHARACTER (LEN=*), PARAMETER     :: RoutineName='RCF_GRIB_READ_DATA'

CHARACTER (LEN=errormessagelength)     :: Cmessage      ! used for EReport
INTEGER                          :: ErrorStatus   ! used for EReport
INTEGER                          :: len_read_data ! used for buffin
INTEGER                          :: Error         ! Error Status
REAL                             :: RError         ! Error Status
INTEGER                          :: act_io_len

!=======================================================================
!  Variables used for calls to 'DECODE'
!=======================================================================

INTEGER , PARAMETER    :: LenVertCord = 3000    ! Len of Vert_Coords
INTEGER , PARAMETER    :: LenBitmap   = LenArrayMax ! Len of bitmap used
INTEGER , PARAMETER    :: LenQuasi    =    1    ! Len of array Quasi

!Variables used when calling DECODE
INTEGER                :: Bitmap(LenBitmap)    ! Array for bitmap of
                                               ! non full field data
INTEGER                :: Quasi(LenQuasi)      ! desc of Quasi-reg grid
INTEGER                :: I_Record(Len_Data)   ! Int array holding
                                               ! Character array of
                                               !encoded GRIB data
REAL                   :: VertCoords(LenVertCord) ! Vertical coord params
CHARACTER(LEN=1)       :: C_Record(Len_Data) ! Encoded data read into
                                               ! this using Buffin
! Lets size this to 2 8-byte INTEGERs to support GRIB1 and GRIB2.
INTEGER                :: section0(2)
INTEGER                :: grib_version

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!=======================================================================
!  Read binary form of record into CRecord work array
!=======================================================================

C_Record(:) = " "                           ! Make sure buffer is clear
act_io_len  = 0

! Set Position in file to be read
CALL SetPos8( Unit_Num, Pos_in_File, Error)

IF ( Error /= 0 ) THEN
  WRITE(Cmessage,'(A,I8)') 'Failed trying to SetPos to ',        &
                               Pos_in_File
  ErrorStatus = 30
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

! Read the section 0 block of data into a 8 byte integer.
CALL Buffin (Unit_Num, section0,                 &
             2, act_io_len, RError)

IF ( RError > 0.0 ) THEN
  Cmessage    = 'Failed trying to read section 0 from GRIB file'
  ErrorStatus = 40
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

len_read_data = 0
grib_version = IBITS(section0(1),0,8)
IF (grib_version == 1) THEN
  ! Within section 0 of a GRIB file 5:7 bytes contains the length.
  len_read_data = IBITS(section0(1), 8, 24)
ELSE IF (grib_version == 2) THEN
  ! Within section 0 of a GRIB2 file 9:16 bytes contains the length.
  len_read_data = section0(2)
ELSE
  WRITE(umMessage,'(A,I0)') "GRIB version found is ", grib_version
  CALL umPrint(umMessage,src='rcf_grib_read_data_mod')
  Cmessage    = 'Unknown GRIB version whilst extracting length of message.'
  ErrorStatus = 41
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

IF (len_read_data > len_data) THEN
  WRITE(umMessage,'(A,I0)') "GRIB message maximum allowed length ", len_data
  CALL umPrint(umMessage,src='rcf_grib_read_data_mod')
  WRITE(umMessage,'(A,I0)') "GRIB message length read in ", len_read_data
  CALL umPrint(umMessage,src='rcf_grib_read_data_mod')
  Cmessage    = 'GRIB message is bigger than expected maximum buffer size'
  ErrorStatus = 41
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

! Set Position in file to be read again to reread data.
CALL SetPos8( Unit_Num, Pos_in_File, Error)

! Read the block of data into CRecord
CALL Buffin (Unit_Num, C_Record(:),                 &
             len_read_data, 1, act_io_len, RError)

IF ( RError > 0.0 ) THEN
  Cmessage    = 'Failed trying to read record from GRIB file'
  ErrorStatus = 42
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

! Report on progress
IF ( PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
  WRITE(umMessage,*) "Whilst trying to read full record tried for ",       &
               Len_read_Data, "and got", act_io_len + 8
  CALL umPrint(umMessage,src='rcf_grib_read_data_mod')
END IF


!=======================================================================
!  Prepare for and Call DECODE routine
!=======================================================================
! Initialise some of the arrays used in the call to DECODE
I_Record(:)           = 0
Current % Num_Fp      = 0
Current % Num_vert    = 0
Current % Num_Bitmap  = 0
Current % Num_Quasi   = 0

! Despite not being used DECODE expects this variable to be 0
Current % Block_4(p_B4Undef) = 0

! Transfer the character array that the encoded data is in into
! an Integer array which DECODE expects

I_Record(:Len_read_Data) = TRANSFER(C_Record(:Len_read_Data),                 &
                                    I_Record(1), Len_read_Data)
CALL ukmo_decode(DATA, Len_Data, Current % Num_Fp,               &
            VertCoords, LenVertCord, Current % Num_vert,         &
            Bitmap, LenBitmap, Current % Num_Bitmap,             &
            Quasi, LenQuasi, Current % Num_Quasi,                &
            Current % Block_0,                                   &
            Current % Block_1,                                   &
            Current % Block_2,                                   &
            Current % Block_4,                                   &
            Current % Block_R,                                   &
            Current % Computed_Keys,                             &
            I_Record, act_io_len)

! If vertical coordinates exist, store in grib record
IF (Current % Num_vert /= 0) THEN
  IF ( mype == 0 ) THEN
    WRITE(umMessage,*) "Storing vertical coordinates"
    CALL umPrint(umMessage,src='rcf_grib_read_data_mod')
  END IF
  ALLOCATE(Current % VertCoords(Current % Num_vert))
  Current % VertCoords=VertCoords(1:Current % Num_vert)
END IF

! Set missing data values properly in bitmaps
IF (Current % Num_Bitmap /= 0) THEN
  IF ( mype == 0 ) THEN
    WRITE(umMessage,*) "Setting MDI properly for bitmap fields"
    CALL umPrint(umMessage,src='rcf_grib_read_data_mod')
  END IF
  WHERE (DATA == Current % Block_R(2)) DATA=rmdi
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_Grib_Read_Data
END MODULE Rcf_Grib_Read_Data_Mod
