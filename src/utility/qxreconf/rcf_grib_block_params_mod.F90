! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Define the 'Grib_record' type and the parameters for the 'blocks'

MODULE rcf_GRIB_Block_Params_Mod

! Description: Parameters defining the elements of the block data read
!              from GRIB files.
!              Plus the derived types in which the data is stored
!              within reconfiguration
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!

!Pull in some parameters from other modules

USE Rcf_GRIB_Lookups_Mod, ONLY:    &
  LenDesc

USE ukmo_grib_mod, ONLY: grib_api_computed_keys

IMPLICIT NONE

!=======================================================================
!  Block Sizes
!=======================================================================

INTEGER , PARAMETER    :: Len_Block_0 =    4   !\ The Blocks used
INTEGER , PARAMETER    :: Len_Block_1 = 1000   !/ to store the
INTEGER , PARAMETER    :: Len_Block_2 =   20   !\ GRIB header
INTEGER , PARAMETER    :: Len_Block_3 =    2   !/ Info
INTEGER , PARAMETER    :: Len_Block_4 =    2   !\                      .
INTEGER , PARAMETER    :: Len_Block_R =   20   !/                      .

!=======================================================================
!  Variables Usedfor calls to 'DECODE'
!=======================================================================

INTEGER , PARAMETER    :: LenArrayMax =15000000 ! Max Array size used
                                                ! for work arrays and
                                                ! such in DECODE

REAL                   :: FpData(LenArrayMax)   ! Array decoded field is
                                                ! returned in.
LOGICAL                :: LgData(LenArrayMax)   ! Logical Array for
                                                ! Logical fields (LSM)

!=======================================================================
!  The Derived Type used to hold the field header Information
!=======================================================================

TYPE Grib_Record

  INTEGER               :: Start_pos            !\ Start pos in file
  INTEGER               :: StashCode            !  Stash Code for field
  INTEGER               :: Data_Type            !  Type of data

  INTEGER               :: Num_Fp               !\ Number of elements
  INTEGER               :: Num_Vert             !/ returned by DECODE
  INTEGER               :: Num_Bitmap           !\ in the respective
  INTEGER               :: Num_Quasi            !/ arrays.

  INTEGER               :: Block_0(Len_Block_0) !\ The blocks holding
  INTEGER               :: Block_1(Len_Block_1) !/ the header info
  INTEGER               :: Block_2(Len_Block_2) !\ returned by DECODE
  INTEGER               :: Block_3(Len_Block_3) !/
  INTEGER               :: Block_4(Len_Block_4) !\

  REAL                  :: Block_R(Len_Block_R)

  TYPE(grib_api_computed_keys) :: computed_keys ! Requested computed
                                                ! keys

  REAL,POINTER          :: VertCoords(:)        ! Points to dynamically
                                                ! allocated array

  CHARACTER(LEN=lenDesc) :: Desc                ! Description of param

  TYPE(Grib_Record), POINTER     :: Next        !\ Pointers to the next
  TYPE(Grib_Record), POINTER     :: Prev        !/ and prev in list

END TYPE Grib_Record

INTEGER , PARAMETER    :: Grb_Data_Real =  1   !/
INTEGER , PARAMETER    :: Grb_Data_Int  =  2   !/
INTEGER , PARAMETER    :: Grb_Data_Log  =  3   !/
!=======================================================================
!  The Derived Type to /mark/ the lists of fields
!=======================================================================
! note- if the list has only one member then begin and end will point
! to the same place.

TYPE List_Marker
  TYPE(Grib_Record),POINTER      :: Begin
  TYPE(Grib_Record),POINTER      :: END
  INTEGER                        :: LstCount
END TYPE List_Marker

!=======================================================================
!  Block 0
!=======================================================================

INTEGER, PARAMETER  :: p_Ed_no       =  1  ! GRIB edition number
INTEGER, PARAMETER  :: p_Tbl_Vers_No =  2  ! Tables version no.
INTEGER, PARAMETER  :: p_Mes_Len     =  3  ! Total message (record) len
INTEGER, PARAMETER  :: p_B0_undef    =  4  ! currently undefined

!=======================================================================
!  Block 1
!=======================================================================

INTEGER, PARAMETER  :: p_Orig_cntr   =  1  ! Originating center
INTEGER, PARAMETER  :: p_Gen_Proc_ID =  2  ! Generating process ID
INTEGER, PARAMETER  :: p_Grd_ID      =  3  ! Grid ID no.
INTEGER, PARAMETER  :: p_Blck_ID     =  4  ! Block ID flags (table 1)
INTEGER, PARAMETER  :: p_Param_ID    =  5  ! Parameter ID (Table 2)
INTEGER, PARAMETER  :: p_Lvl_Type    =  6  ! Type of Level ID (table 3)
INTEGER, PARAMETER  :: p_Lvl_Desc_1  =  7  ! 1st level desription param
INTEGER, PARAMETER  :: p_Lvl_Desc_2  =  8  ! 2nd level desription param
                                           ! either are 0 if not needed
                                           ! otherwise specify hieght,
                                           ! pressure etc
INTEGER, PARAMETER  :: p_Ref_Year    =  9  ! Year of century, ref time
INTEGER, PARAMETER  :: p_Ref_Month   = 10  ! Month of ref time
INTEGER, PARAMETER  :: p_Ref_Day     = 11  ! Day of Reference Time
INTEGER, PARAMETER  :: p_Ref_Hour    = 12  !
INTEGER, PARAMETER  :: p_Ref_Min     = 13  !
INTEGER, PARAMETER  :: p_Time_Unit   = 14  ! Indicator of time unit
                                           ! (Table 4)
INTEGER, PARAMETER  :: p_Time_Int_1  = 15  ! 1st time interval if reqd
INTEGER, PARAMETER  :: p_Time_Int_2  = 16  ! 2nd time interval (P2)
INTEGER, PARAMETER  :: p_Time_Range  = 17  ! Time range indicator
                                           ! Table 5
INTEGER, PARAMETER  :: p_No_in_Avg   = 18  ! No. included in Average
                                           ! if reqd for time range
INTEGER, PARAMETER  :: p_Ref_Cent    = 19  ! Century of ref time
INTEGER, PARAMETER  :: p_Dec_Scale   = 20  ! Decimal Scale factor
INTEGER, PARAMETER  :: p_Len_Blck_1  = 21  ! Length of block 1 in octets

! elements 22-33 are reserved for future use

INTEGER, PARAMETER  :: p_Sub_C_ID    = 34  ! Sub center ID no.

! octets 36 onwards contain originating center octets.
! 1 octet is returned in each array element. The number of originating
! center octets is Length of block 1 - 40
! If center = 74 and Sub center Id is 2 then 36-79 hold superstash pp
! data in 32 bit words.

!=======================================================================
!  Block 2
!=======================================================================

! The 1st 3 elements of block 2 are always the same.
! The rest vary depending on grid type.
! At present, only lat/long, Gaussian, mercator grid params
! are specified here. i.e polar stereographic are not mentioned

INTEGER, PARAMETER  :: p_V_Coord_Ps =  1  ! no. Vertical coord params
INTEGER, PARAMETER  :: p_PV_PL      =  2  ! PV or PL if level is hybrid
                                          ! or quasi regular
INTEGER, PARAMETER  :: p_Data_Rep   =  3  ! Data representation type

! Assuming Lat/long, Gaussian, Mercator or Rotated Lat/long

INTEGER, PARAMETER  :: p_Pnts_Prll  =  4  ! No. of points along parall
INTEGER, PARAMETER  :: p_Pnts_Merid =  5  ! No. of points along merid
INTEGER, PARAMETER  :: p_LatGrdPnt1 =  6  ! Lat of 1st Grid point
INTEGER, PARAMETER  :: p_LonGrdPnt1 =  7  ! Long of 1st Grid point
INTEGER, PARAMETER  :: p_Res_Flg    =  8  ! resolution and cmpnt flags
                                          ! as defined in Table 7
INTEGER, PARAMETER  :: p_LatExtrmPt =  9  ! Lat of extreme point
INTEGER, PARAMETER  :: p_LonExtrmPt = 10  ! Long of extreme point
INTEGER, PARAMETER  :: p_IncrPrll   = 11  ! Increment along Parrallel
INTEGER, PARAMETER  :: p_IncrMerid  = 12  ! Increment along Meridian
INTEGER, PARAMETER  :: p_ScanModeF  = 13  ! Scanning mode flags as
                                          ! defined in Table 8
INTEGER, PARAMETER  :: p_S_Pole_Lat = 14  ! Lat of S Pole -Iff rotated
INTEGER, PARAMETER  :: p_S_Pole_Lon = 15  ! Long of S Pole -Iff rotated

! Elements 16 to 20 are undefined

!=======================================================================
!  Block 3
!=======================================================================

! Block 3 is currently not held in the GRIB_Record derived data type
! and is not defined in/by DECODE

!=======================================================================
!  Block 4
!=======================================================================

INTEGER, PARAMETER  :: p_PackFlag  =  1  ! Packing Flag
INTEGER, PARAMETER  :: p_B4Undef   =  2  ! Input value, unchanged

!=======================================================================
!  Block R
!=======================================================================

REAL,    PARAMETER  :: p_RotAngl   =  1  ! Angle of rotation for rotated
                                         ! Lat long grid
REAL,    PARAMETER  :: p_GMDI      =  2  ! Missing data Indicator used
                                         ! for bitmap grids

! Elements 3 to 20 are currently undefined

!=======================================================================
!  Miscellaneous values set up as parameters
!=======================================================================

! Values to do with dump header addressing
INTEGER, PARAMETER  :: p_Len_IntC    = 46 ! normally from start dump
INTEGER, PARAMETER  :: p_Len_RealC   = 38 ! ditto -or- Recon namelist
INTEGER, PARAMETER  :: p_Len2_LevDepC = 8 ! ditto -or- Recon namelist
INTEGER, PARAMETER  :: p_Len1_RowDepC = 0 ! ditto -or- Recon namelist
INTEGER, PARAMETER  :: p_Len2_RowDepC = 0 ! ditto -or- Recon namelist
INTEGER, PARAMETER  :: p_Len1_ColDepC = 0 ! ditto -or- Recon namelist
INTEGER, PARAMETER  :: p_Len2_ColDepC = 0 ! ditto -or- Recon namelist
INTEGER, PARAMETER  :: p_Len1_FldsofC = 0 ! ditto -or- Recon namelist
INTEGER, PARAMETER  :: p_Len2_FldsofC = 0 ! ditto -or- Recon namelist
INTEGER, PARAMETER  :: p_Len_ExtraC   = 0 ! ditto -or- Recon namelist
INTEGER, PARAMETER  :: p_Len_HistFile = 0 ! ditto -or- Recon namelist
INTEGER, PARAMETER  :: p_Len_CompFldI1= 0 ! ditto -or- Recon namelist
INTEGER, PARAMETER  :: p_Len_CompFldI2= 0 ! ditto -or- Recon namelist
INTEGER, PARAMETER  :: p_Len_CompFldI3= 0 ! ditto -or- Recon namelist
INTEGER, PARAMETER  :: p_Len1_Lookup = 64 ! ditto -or- Recon namelist

!=======================================================================
!  Values Used in GRIB Table 2 (Parameter IDs)
!=======================================================================

!Those used by ECMWF (the non-standard ones)
INTEGER, PARAMETER  :: EID_Geopotential   = 129
INTEGER, PARAMETER  :: EID_Temperature    = 130
INTEGER, PARAMETER  :: EID_U_Wind         = 131
INTEGER, PARAMETER  :: EID_V_Wind         = 132
INTEGER, PARAMETER  :: EID_Spec_Humidity  = 133
INTEGER, PARAMETER  :: EID_Surf_Press     = 134
INTEGER, PARAMETER  :: EID_W_Wind         = 135
INTEGER, PARAMETER  :: EID_Surf_Temp      = 139
INTEGER, PARAMETER  :: EID_Soil_Moisture  = 140
INTEGER, PARAMETER  :: EID_Log_Surf_Press = 152

!=======================================================================
!  Values Used in GRIB Table 3 (Type and Value of level)
!=======================================================================

INTEGER, PARAMETER  :: Tb3_Surface  = 1   ! Ground or Water Surface
INTEGER, PARAMETER  :: Tb3_Pressure = 100 ! Isobaric Level type
INTEGER, PARAMETER  :: Tb3_Mid_IsoB = 101 ! layer between Isobaric
INTEGER, PARAMETER  :: Tb3_Mean_Sea = 102 ! Mean Sea Level
INTEGER, PARAMETER  :: Tb3_Hybrid   = 109 ! Hybrid Level
INTEGER, PARAMETER  :: Tb3_GrndLay  = 112 ! Below ground layer

!=======================================================================
!  Done
!=======================================================================

END MODULE rcf_GRIB_Block_Params_Mod
