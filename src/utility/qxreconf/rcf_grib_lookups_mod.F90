! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  List parameters and cross reference info for GRIB data

MODULE rcf_GRIB_lookups_Mod

! Description: This module stores the tables (in the form of derived
!              types) and parameters used to cross reference STASH IDs
!              with those used by other centres.
!
!              To add data for another centre just add another column
!              and use 'Is_not_def' for any ID already in the table
!              that any centre has no equivalent for.
!
!              Note also that for any new pressure level fields
!              rcf_grib_control_mod should be set up to indicate
!              whether the grib pressure level fields are ordered
!              away from surface (ECMWF).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!=======================================================================
!  The 'Lists' - used for storage of GRIB header data
!=======================================================================
! These parameters define
! 1) How many linked lists are used to store the fields being read
!    (1 list per variable - multiple 'fields' or 'levels' per list)
! 2) Which list each variable is stored in.
! 3) The numerical order of the lists defines the order in which they
!    Are written to the intermediary UM formatted dump
!
!
! Use all of RCF_STASHCODES for stash magic numbers
USE um_stashcode_mod

! Use all of RCF_ECMWFCODES for ECMWF magic numbers
USE Rcf_ECMWFcodes_Mod
USE Rcf_UKMOcodes_Mod

IMPLICIT NONE
INTEGER , PARAMETER :: grib_max_fields    = 45 ! no. of lists

INTEGER , PARAMETER              :: grib_LandMask_field    =  1
INTEGER , PARAMETER              :: grib_Orog_field        =  2
INTEGER , PARAMETER              :: grib_Surf_Temp_field   =  3
INTEGER , PARAMETER              :: grib_Surf_Pres_field   =  4
INTEGER , PARAMETER              :: grib_U_field           =  5
INTEGER , PARAMETER              :: grib_V_field           =  6
INTEGER , PARAMETER              :: grib_W_field           =  7
INTEGER , PARAMETER              :: grib_Temp_field        =  8
INTEGER , PARAMETER              :: grib_Q_field           =  9
INTEGER , PARAMETER              :: grib_Rel_Hum_field     = 10
INTEGER , PARAMETER              :: grib_Exner_field       = 11
INTEGER , PARAMETER              :: grib_Soil_Temp_field   = 12
INTEGER , PARAMETER              :: grib_Soil_Moist_field  = 13
INTEGER , PARAMETER              :: grib_ozone_field       = 14
INTEGER , PARAMETER              :: grib_NOX_field         = 15
INTEGER , PARAMETER              :: grib_CH4_field         = 16
INTEGER , PARAMETER              :: grib_CO_field          = 17
INTEGER , PARAMETER              :: grib_HCHO_field        = 18
INTEGER , PARAMETER              :: grib_GO3_field         = 19
INTEGER , PARAMETER              :: grib_NO2_field         = 20
INTEGER , PARAMETER              :: grib_NO_field          = 21
INTEGER , PARAMETER              :: grib_SO2_field         = 22
INTEGER , PARAMETER              :: grib_HNO3_field        = 23
INTEGER , PARAMETER              :: grib_PAN_field         = 24
INTEGER , PARAMETER              :: grib_C2H6_field        = 25
INTEGER , PARAMETER              :: grib_C3H8_field        = 26
INTEGER , PARAMETER              :: grib_OMFRSH_field      = 27
INTEGER , PARAMETER              :: grib_OMAGD_field       = 28
INTEGER , PARAMETER              :: grib_BCFRSH_field      = 29
INTEGER , PARAMETER              :: grib_BCAGD_field       = 30
INTEGER , PARAMETER              :: grib_UMDUST1_field     = 31
INTEGER , PARAMETER              :: grib_UMDUST2_field     = 32
INTEGER , PARAMETER              :: grib_UMDUST3_field     = 33
INTEGER , PARAMETER              :: grib_UMDUST4_field     = 34
INTEGER , PARAMETER              :: grib_UMDUST5_field     = 35
INTEGER , PARAMETER              :: grib_UMDUST6_field     = 36
INTEGER , PARAMETER              :: grib_UMSO4AITK_field   = 37
INTEGER , PARAMETER              :: grib_UMSO4ACCU_field   = 38
INTEGER , PARAMETER              :: grib_Snow_field        = 39
INTEGER , PARAMETER              :: grib_SeaIce_field      = 40
INTEGER , PARAMETER              :: grib_QCL_field         = 41
INTEGER , PARAMETER              :: grib_QCF_field         = 42
INTEGER , PARAMETER              :: grib_CC_field          = 43
INTEGER , PARAMETER              :: grib_QRAIN_field       = 44
INTEGER , PARAMETER              :: grib_QSNOW_field       = 45

INTEGER , PARAMETER              :: grib_Misc_field        =  0
! The misc field is used to store header data not stored in the
! other lists

!=======================================================================
!  Parameters : The Originating centre no.s
!=======================================================================
! These parameters correspond the the values of the Originating centres
! supplied by the GRIB record - Block 1 - Octet 1. They are used within
! rcf_grib_assign.F90 to select which column in table A (below) to use
! to find the variables STASH code from its given GRIB parameter code
! (Table 2 or Block 1 - Octet 5)

INTEGER, PARAMETER  :: GrbOrigECMWF  = 98  ! ECMWF
INTEGER, PARAMETER  :: GrbOrigUKMO   = 74  ! Met Office
INTEGER, PARAMETER  :: GrbOrigIntrnl = -1  ! Internal data (created by
                                           ! GRIB code not read)
!=======================================================================
!  Parameters : The Table no.s
!=======================================================================
! These are set so that entries can be distinguished by table number.

INTEGER, PARAMETER :: GrbTblECMWFstd  = 128 ! Standard table
INTEGER, PARAMETER :: GrbTblECMWFgems = 210 ! GEMS Chemistry table
INTEGER, PARAMETER :: GrbTblUKMOgems  = 141 ! UKMO variables from GEMS

!=======================================================================
!  Table A: The Parameter ID cross reference table
!  Part 1 : Defining the table
!=======================================================================

INTEGER, PARAMETER  :: p_Max_Cols    =  4  ! Number of columns in table
                                           ! 1 per Centre
                                           ! (inc UK Met O stashcodes)

INTEGER, PARAMETER  :: p_Max_Rows    = 53  ! Current no. of paramters
                                           ! listed. No real relation to
                                           ! no. of lists

INTEGER, PARAMETER  :: Is_not_def    = -1  ! Used when no valid ID is
                                           ! available at the centre in
                                           ! question.

INTEGER , PARAMETER :: lenDesc       = 20  ! length of description

TYPE Cross_Ref_Table
  CHARACTER(LEN=lenDesc) :: DescText           ! Textual desc of param
  INTEGER           :: CrossRefIDs(p_Max_Cols) ! The columns of ID nos
  INTEGER           :: Table_No(p_Max_Cols)    ! Grib table number for field
  INTEGER           :: List_No                 ! List Number for dynamic lists
END TYPE Cross_Ref_Table

TYPE (Cross_Ref_Table)           :: Param_ID_CrossRef(p_Max_Rows)

!=======================================================================
!  Parameters : The ID columns
!=======================================================================
! These parameters define which column (in Table A) to use for a given
! set of variable ID no.s - e.g. column 1 is used to hold STASH no.s
! and column 2 for ECMWF's version of Table 2 parameter IDs

INTEGER, PARAMETER  :: p_STASH_IDCol  =  1  ! Column used to store STASH
                                            ! ID no's
INTEGER, PARAMETER  :: p_ECMWF_IDCol  =  2  ! Column used to store ECMWF
                                            ! ID no's
INTEGER, PARAMETER  :: p_UKMO_IDCol   =  3  ! Column used to store UKMO
                                            ! ID no's in grib files

! The last column is used to store locally defined ID no.s for data
! which will be generated internally by the GRIB code itself and not
! read in from the GRIB data.
INTEGER, PARAMETER  :: p_Intrnl_IDCol =  p_Max_Cols

!=======================================================================
!  Table A: The Parameter ID cross reference table
!  Part 2 : Populating the Table
!=======================================================================
! To add a definition for a new parameter, just append another block of
! data definitions.
! To add a new centre to get data from, you must add the ID code to
! each (all) of the CrossRefIDs.

! Example
!Data  Param_ID_CrossRef(1) % CrossRefIDs / -STASH-,-ECMWF-,-UKMO-,-INTERNAL- /
!Data  Param_ID_CrossRef(1) % DescText    /'Description here    '/
!Data  Param_ID_CrossRef(1) % List_No     /list no. parameter/
!Data  Param_ID_CrossRef(1) % Table_No    / -STASH-,-ECMWF-,-UKMO-,-INTERNAL- /

DATA  Param_ID_CrossRef( 1) % CrossRefIDs / stashcode_u, &
                                            ECMWFcode_u, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef( 1) % DescText    /'U Wind              '/
DATA  Param_ID_CrossRef( 1) % List_No     /grib_U_field/
DATA  Param_ID_CrossRef( 1) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef( 2) % CrossRefIDs / stashcode_v, &
                                            ECMWFcode_v, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef( 2) % DescText    /'V Wind              '/
DATA  Param_ID_CrossRef( 2) % List_No     /grib_V_field/
DATA  Param_ID_CrossRef( 2) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef( 3) % CrossRefIDs / stashcode_t,     &
                                            ECMWFcode_T,     &
                                            Is_not_def,      &
                                            Is_not_def/
DATA  Param_ID_CrossRef( 3) % DescText    /'Temperature         '/
DATA  Param_ID_CrossRef( 3) % List_No     /grib_Temp_field/
DATA  Param_ID_CrossRef( 3) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef( 4) % CrossRefIDs / stashcode_q,        &
                                            ECMWFcode_spec_hum, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef( 4) % DescText    /'Specific Humidity   '/
DATA  Param_ID_CrossRef( 4) % List_No     /grib_Q_field/
DATA  Param_ID_CrossRef( 4) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/


DATA  Param_ID_CrossRef( 5) % CrossRefIDs / stashcode_soil_temp, &
                                            ECMWFcode_soil_temp_1, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef( 5) % DescText    /'Soil Temp Lvl 1     '/
DATA  Param_ID_CrossRef( 5) % List_No     /grib_Soil_Temp_field/
DATA  Param_ID_CrossRef( 5) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef( 6) % CrossRefIDs / stashcode_soil_temp, &
                                            ECMWFcode_soil_temp_2, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef( 6) % DescText    /'Soil Temp Lvl 2     '/
DATA  Param_ID_CrossRef( 6) % List_No     /grib_Soil_Temp_field/
DATA  Param_ID_CrossRef( 6) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef( 7) % CrossRefIDs / stashcode_soil_temp, &
                                            ECMWFcode_soil_temp_3, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef( 7) % DescText    /'Soil Temp Lvl 3     '/
DATA  Param_ID_CrossRef( 7) % List_No     /grib_Soil_Temp_field/
DATA  Param_ID_CrossRef( 7) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef( 8) % CrossRefIDs / stashcode_soil_temp, &
                                            ECMWFcode_soil_temp_4, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef( 8) % DescText    /'Soil Temp Lvl 4     '/
DATA  Param_ID_CrossRef( 8) % List_No     /grib_Soil_Temp_field/
DATA  Param_ID_CrossRef( 8) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef( 9) % CrossRefIDs / stashcode_tstar,     &
                                            ECMWFcode_skin_temp, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef( 9) % DescText    /'Skin Temperature    '/
DATA  Param_ID_CrossRef( 9) % List_No     /grib_Surf_Temp_field/
DATA  Param_ID_CrossRef( 9) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(10) % CrossRefIDs / stashcode_tstar, &
                                            Is_not_def,      &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(10) % DescText    /'Surface Temperature '/
DATA  Param_ID_CrossRef(10) % List_No     /grib_Surf_Temp_field/
DATA  Param_ID_CrossRef(10) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/


DATA  Param_ID_CrossRef(11) % CrossRefIDs / stashcode_lsm, &
                                            ECMWFcode_lsm, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(11) % DescText    /'Land/Sea Mask       '/
DATA  Param_ID_CrossRef(11) % List_No     /grib_LandMask_field/
DATA  Param_ID_CrossRef(11) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(12) % CrossRefIDs / stashcode_orog,   &
                                            ECMWFcode_geopot, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(12) % DescText    /'Geopotential        '/
DATA  Param_ID_CrossRef(12) % List_No     /grib_Orog_field/
DATA  Param_ID_CrossRef(12) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(13) % CrossRefIDs / stashcode_w, &
                                            ECMWFcode_w, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(13) % DescText    /'W Wind              '/
DATA  Param_ID_CrossRef(13) % List_No     /grib_W_field/
DATA  Param_ID_CrossRef(13) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(14) % CrossRefIDs / stashcode_pstar, &
                                            ECMWFcode_pstar, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(14) % DescText    /'Surface Pressure    '/
DATA  Param_ID_CrossRef(14) % List_No     /grib_Surf_Pres_field/
DATA  Param_ID_CrossRef(14) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(15) % CrossRefIDs / stashcode_q, &
                                            Is_not_def,  &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(15) % DescText    /'Relative Humidity   '/
DATA  Param_ID_CrossRef(15) % List_No     /grib_Rel_Hum_field/
DATA  Param_ID_CrossRef(15) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(16) % CrossRefIDs / stashcode_exner, &
                                            Is_not_def,  &
                                            Is_not_def,  &
                                            1/
DATA  Param_ID_CrossRef(16) % DescText    /'Exner on P levels   '/
DATA  Param_ID_CrossRef(16) % List_No     /grib_Exner_field/
DATA  Param_ID_CrossRef(16) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,  &
                                           is_not_def/

DATA  Param_ID_CrossRef(17) % CrossRefIDs / stashcode_pstar, &
                                            ECMWFcode_log_p, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(17) % DescText    /'Log Surface Pressure'/
DATA  Param_ID_CrossRef(17) % List_No     /grib_Surf_Pres_field/
DATA  Param_ID_CrossRef(17) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(18) % CrossRefIDs / stashcode_soil_moist, &
                                            ECMWFcode_soil_moist_1, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(18) % DescText    /'Vol Soil Moist Lvl 1'/
DATA  Param_ID_CrossRef(18) % List_No     /grib_Soil_Moist_field/
DATA  Param_ID_CrossRef(18) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(19) % CrossRefIDs / stashcode_soil_moist, &
                                            ECMWFcode_soil_moist_2, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(19) % DescText    /'Vol Soil Moist Lvl 2'/
DATA  Param_ID_CrossRef(19) % List_No     /grib_Soil_Moist_field/
DATA  Param_ID_CrossRef(19) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(20) % CrossRefIDs / stashcode_soil_moist, &
                                            ECMWFcode_soil_moist_3, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(20) % DescText    /'Vol Soil Moist Lvl 3'/
DATA  Param_ID_CrossRef(20) % List_No     /grib_Soil_Moist_field/
DATA  Param_ID_CrossRef(20) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(21) % CrossRefIDs / stashcode_soil_moist, &
                                            ECMWFcode_soil_moist_4, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(21) % DescText    /'Vol Soil Moist Lvl 4'/
DATA  Param_ID_CrossRef(21) % List_No     /grib_Soil_Moist_field/
DATA  Param_ID_CrossRef(21) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(22) % CrossRefIDs / stashcode_ozone, &
                                            ECMWFcode_ozone, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(22) % DescText    /'Ozone               '/
DATA  Param_ID_CrossRef(22) % List_No     /grib_Ozone_field/
DATA  Param_ID_CrossRef(22) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(23) % CrossRefIDs / stashcode_NO2, &
                                            ECMWFcode_NOX, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(23) % DescText    /'Nitrogen Oxides     '/
DATA  Param_ID_CrossRef(23) % List_No     /grib_NOX_field/
DATA  Param_ID_CrossRef(23) % Table_No    /is_not_def,      &
                                           grbTblECMWFgems, &
                                           Is_not_def,      &
                                           is_not_def/

DATA  Param_ID_CrossRef(24) % CrossRefIDs / stashcode_CH4, &
                                            ECMWFcode_CH4, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(24) % DescText    /'Methane             '/
DATA  Param_ID_CrossRef(24) % List_No     /grib_CH4_field/
DATA  Param_ID_CrossRef(24) % Table_No    /is_not_def,      &
                                           grbTblECMWFgems, &
                                           Is_not_def,      &
                                           is_not_def/

DATA  Param_ID_CrossRef(25) % CrossRefIDs / stashcode_CO, &
                                            ECMWFcode_CO, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(25) % DescText    /'Carbon monoxide     '/
DATA  Param_ID_CrossRef(25) % List_No     /grib_CO_field/
DATA  Param_ID_CrossRef(25) % Table_No    /is_not_def,      &
                                           grbTblECMWFgems, &
                                           Is_not_def,      &
                                           is_not_def/

DATA  Param_ID_CrossRef(26) % CrossRefIDs / stashcode_HCHO, &
                                            ECMWFcode_HCHO, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(26) % DescText    /'Formaldehyde        '/
DATA  Param_ID_CrossRef(26) % List_No     /grib_HCHO_field/
DATA  Param_ID_CrossRef(26) % Table_No    /is_not_def,      &
                                           grbTblECMWFgems, &
                                           Is_not_def,      &
                                           is_not_def/

DATA  Param_ID_CrossRef(27) % CrossRefIDs / stashcode_O3,  &
                                            ECMWFcode_GO3, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(27) % DescText    /'GEMS ozone          '/
DATA  Param_ID_CrossRef(27) % List_No     /grib_GO3_field/
DATA  Param_ID_CrossRef(27) % Table_No    /is_not_def,      &
                                           grbTblECMWFgems, &
                                           Is_not_def,      &
                                           is_not_def/

DATA  Param_ID_CrossRef(28) % CrossRefIDs / stashcode_NO2, &
                                            ECMWFcode_NO2, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(28) % DescText    /'Nitrogen dioxide    '/
DATA  Param_ID_CrossRef(28) % List_No     /grib_NO2_field/
DATA  Param_ID_CrossRef(28) % Table_No    /is_not_def,      &
                                           grbTblECMWFgems, &
                                           Is_not_def,      &
                                           is_not_def/

DATA  Param_ID_CrossRef(29) % CrossRefIDs / stashcode_NO,  &
                                            Is_not_def,    &
                                            UKMOcode_NO,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(29) % DescText    /'Nitrogen monoxide   '/
DATA  Param_ID_CrossRef(29) % List_No     /grib_NO_field/
DATA  Param_ID_CrossRef(29) % Table_No    /is_not_def,      &
                                           Is_not_def,      &
                                           GrbTblUKMOgems,  &
                                           is_not_def/

DATA  Param_ID_CrossRef(30) % CrossRefIDs / stashcode_SO2, &
                                            Is_not_def,    &
                                            UKMOcode_SO2,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(30) % DescText    /'Sulphur dioxide     '/
DATA  Param_ID_CrossRef(30) % List_No     /grib_SO2_field/
DATA  Param_ID_CrossRef(30) % Table_No    /is_not_def,      &
                                           Is_not_def,      &
                                           GrbTblUKMOgems,  &
                                           is_not_def/

DATA  Param_ID_CrossRef(31) % CrossRefIDs / stashcode_HNO3, &
                                            Is_not_def,     &
                                            UKMOcode_HNO3, &
                                            Is_not_def/
DATA  Param_ID_CrossRef(31) % DescText    /'Nitric acid         '/
DATA  Param_ID_CrossRef(31) % List_No     /grib_HNO3_field/
DATA  Param_ID_CrossRef(31) % Table_No    /is_not_def,      &
                                           Is_not_def,      &
                                           GrbTblUKMOgems , &
                                           is_not_def/

DATA  Param_ID_CrossRef(32) % CrossRefIDs / stashcode_PAN, &
                                            Is_not_def,    &
                                            UKMOcode_PAN, &
                                            Is_not_def/
DATA  Param_ID_CrossRef(32) % DescText    /'PAN                 '/
DATA  Param_ID_CrossRef(32) % List_No     /grib_PAN_field/
DATA  Param_ID_CrossRef(32) % Table_No    /is_not_def,      &
                                           Is_not_def,      &
                                           GrbTblUKMOgems,  &
                                           is_not_def/

DATA  Param_ID_CrossRef(33) % CrossRefIDs / stashcode_C2H6, &
                                            Is_not_def,     &
                                            UKMOcode_C2H6, &
                                            Is_not_def/
DATA  Param_ID_CrossRef(33) % DescText    /'Ethane              '/
DATA  Param_ID_CrossRef(33) % List_No     /grib_C2H6_field/
DATA  Param_ID_CrossRef(33) % Table_No    /is_not_def,      &
                                           Is_not_def,      &
                                           GrbTblUKMOgems,  &
                                           is_not_def/

DATA  Param_ID_CrossRef(34) % CrossRefIDs / stashcode_C3H8, &
                                            Is_not_def,     &
                                            UKMOcode_C3H8, &
                                            Is_not_def/
DATA  Param_ID_CrossRef(34) % DescText    /'Propane             '/
DATA  Param_ID_CrossRef(34) % List_No     /grib_C3H8_field/
DATA  Param_ID_CrossRef(34) % Table_No    /is_not_def,      &
                                           Is_not_def,      &
                                           GrbTblUKMOgems,  &
                                           is_not_def/

DATA  Param_ID_CrossRef(35) % CrossRefIDs / stashcode_mmr_ocff_fr, &
                                            ECMWFcode_OMFRSH, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(35) % DescText    /'Fresh/hydrophobic OC'/
DATA  Param_ID_CrossRef(35) % List_No     /grib_OMFRSH_field/
DATA  Param_ID_CrossRef(35) % Table_No    /is_not_def,      &
                                           grbTblECMWFgems, &
                                           Is_not_def,      &
                                           is_not_def/

DATA  Param_ID_CrossRef(36) % CrossRefIDs / stashcode_mmr_ocff_ag, &
                                            ECMWFcode_OMAGD, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(36) % DescText    /'Aged/hydrophillic OC'/
DATA  Param_ID_CrossRef(36) % List_No     /grib_OMAGD_field/
DATA  Param_ID_CrossRef(36) % Table_No    /is_not_def,      &
                                           grbTblECMWFgems, &
                                           Is_not_def,      &
                                           is_not_def/

DATA  Param_ID_CrossRef(37) % CrossRefIDs / stashcode_mmr_bc_fr, &
                                            ECMWFcode_BCFRSH, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(37) % DescText    /'Fresh/hydrophobic BC'/
DATA  Param_ID_CrossRef(37) % List_No     /grib_BCFRSH_field/
DATA  Param_ID_CrossRef(37) % Table_No    /is_not_def,      &
                                           grbTblECMWFgems, &
                                           Is_not_def,      &
                                           is_not_def/

DATA  Param_ID_CrossRef(38) % CrossRefIDs / stashcode_mmr_bc_ag, &
                                            ECMWFcode_BCAGD, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(38) % DescText    /'Aged/hydrophillic BC'/
DATA  Param_ID_CrossRef(38) % List_No     /grib_BCAGD_field/
DATA  Param_ID_CrossRef(38) % Table_No    /is_not_def,      &
                                           grbTblECMWFgems, &
                                           Is_not_def,      &
                                           is_not_def/

DATA  Param_ID_CrossRef(39) % CrossRefIDs / stashcode_dust1_mmr, &
                                            Is_not_def,          &
                                            UKMOcode_UMDUST1,    &
                                            Is_not_def/
DATA  Param_ID_CrossRef(39) % DescText    /'Dust div 1          '/
DATA  Param_ID_CrossRef(39) % List_No     /grib_UMDUST1_field/
DATA  Param_ID_CrossRef(39) % Table_No    /is_not_def,      &
                                           Is_not_def,      &
                                           GrbTblUKMOgems,  &
                                           is_not_def/

DATA  Param_ID_CrossRef(40) % CrossRefIDs / stashcode_dust2_mmr, &
                                            Is_not_def,          &
                                            UKMOcode_UMDUST2,    &
                                            Is_not_def/
DATA  Param_ID_CrossRef(40) % DescText    /'Dust div 2          '/
DATA  Param_ID_CrossRef(40) % List_No     /grib_UMDUST2_field/
DATA  Param_ID_CrossRef(40) % Table_No    /is_not_def,           &
                                           Is_not_def,           &
                                           GrbTblUKMOgems,       &
                                           is_not_def/

DATA  Param_ID_CrossRef(41) % CrossRefIDs / stashcode_dust3_mmr, &
                                            Is_not_def,          &
                                            UKMOcode_UMDUST3,    &
                                            Is_not_def/
DATA  Param_ID_CrossRef(41) % DescText    /'Dust div 3          '/
DATA  Param_ID_CrossRef(41) % List_No     /grib_UMDUST3_field/
DATA  Param_ID_CrossRef(41) % Table_No    /is_not_def,           &
                                           Is_not_def,           &
                                           GrbTblUKMOgems,       &
                                           is_not_def/

DATA  Param_ID_CrossRef(42) % CrossRefIDs / stashcode_dust4_mmr, &
                                            Is_not_def,          &
                                            UKMOcode_UMDUST4,    &
                                            Is_not_def/
DATA  Param_ID_CrossRef(42) % DescText    /'Dust div 4          '/
DATA  Param_ID_CrossRef(42) % List_No     /grib_UMDUST4_field/
DATA  Param_ID_CrossRef(42) % Table_No    /is_not_def,           &
                                           Is_not_def,           &
                                           GrbTblUKMOgems,       &
                                           is_not_def/

DATA  Param_ID_CrossRef(43) % CrossRefIDs / stashcode_dust5_mmr, &
                                            Is_not_def,          &
                                            UKMOcode_UMDUST5,    &
                                            Is_not_def/
DATA  Param_ID_CrossRef(43) % DescText    /'Dust div 5          '/
DATA  Param_ID_CrossRef(43) % List_No     /grib_UMDUST5_field/
DATA  Param_ID_CrossRef(43) % Table_No    /is_not_def,           &
                                           Is_not_def,           &
                                           GrbTblUKMOgems,       &
                                           is_not_def/

DATA  Param_ID_CrossRef(44) % CrossRefIDs / stashcode_dust6_mmr, &
                                            Is_not_def,          &
                                            UKMOcode_UMDUST6,   &
                                            Is_not_def/
DATA  Param_ID_CrossRef(44) % DescText    /'Dust div 6          '/
DATA  Param_ID_CrossRef(44) % List_No     /grib_UMDUST6_field/
DATA  Param_ID_CrossRef(44) % Table_No    /is_not_def,           &
                                           Is_not_def,           &
                                           GrbTblUKMOgems,       &
                                           is_not_def/

DATA  Param_ID_CrossRef(45) % CrossRefIDs / stashcode_mmr_so4_aitken, &
                                            Is_not_def,               &
                                            UKMOcode_UMSO4AITK,       &
                                            Is_not_def/
DATA  Param_ID_CrossRef(45) % DescText    /'SO4 aitken mode     '/
DATA  Param_ID_CrossRef(45) % List_No     /grib_UMSO4AITK_field/
DATA  Param_ID_CrossRef(45) % Table_No    /is_not_def,                &
                                           Is_not_def,                &
                                           GrbTblUKMOgems,            &
                                           is_not_def/

DATA  Param_ID_CrossRef(46) % CrossRefIDs / stashcode_mmr_so4_accum, &
                                            Is_not_def,              &
                                            UKMOcode_UMSO4ACCU,     &
                                            Is_not_def/
DATA  Param_ID_CrossRef(46) % DescText    /'SO4 accum. mode     '/
DATA  Param_ID_CrossRef(46) % List_No     /grib_UMSO4ACCU_field/
DATA  Param_ID_CrossRef(46) % Table_No    /is_not_def,              &
                                           Is_not_def,              &
                                           GrbTblUKMOgems,          &
                                           is_not_def/

DATA  Param_ID_CrossRef(47) % CrossRefIDs / stashcode_mean_snow,  &
                                            ECMWFcode_snow_depth, &
                                            Is_not_def,           &
                                            is_not_def/
DATA  Param_ID_CrossRef(47) % DescText    /'Snow Depth          '/
DATA  Param_ID_CrossRef(47) % List_No     /grib_Snow_field/
DATA  Param_ID_CrossRef(47) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(48) % CrossRefIDs / stashcode_icefrac, &
                                            ECMWFcode_icefrac, &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(48) % DescText    /'Seaice Fraction     '/
DATA  Param_ID_CrossRef(48) % List_No     /grib_SeaIce_field/
DATA  Param_ID_CrossRef(48) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(49) % CrossRefIDs / stashcode_qcl,        &
                                            ECMWFcode_qcl,        &
                                            Is_not_def,           &
                                            Is_not_def/
DATA  Param_ID_CrossRef(49) % DescText    /'QCL                 '/
DATA  Param_ID_CrossRef(49) % List_No     /grib_QCL_field/
DATA  Param_ID_CrossRef(49) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(50) % CrossRefIDs / stashcode_qcf,        &
                                            ECMWFcode_qcf,        &
                                            Is_not_def,           &
                                            Is_not_def/
DATA  Param_ID_CrossRef(50) % DescText    /'QCF                 '/
DATA  Param_ID_CrossRef(50) % List_No     /grib_QCF_field/
DATA  Param_ID_CrossRef(50) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(51) % CrossRefIDs / stashcode_area_cf, &
                                            ECMWFcode_cc,      &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(51) % DescText    /'Cloud Cover         '/
DATA  Param_ID_CrossRef(51) % List_No     /grib_CC_field/
DATA  Param_ID_CrossRef(51) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(52) % CrossRefIDs / stashcode_qrain,   &
                                            ECMWFcode_qrain,   &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(52) % DescText    /'QRAIN               '/
DATA  Param_ID_CrossRef(52) % List_No     /grib_QRAIN_field/
DATA  Param_ID_CrossRef(52) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

DATA  Param_ID_CrossRef(53) % CrossRefIDs / stashcode_qgraup,  &
                                            ECMWFcode_qsnow,   &
                                            Is_not_def,  &
                                            Is_not_def/
DATA  Param_ID_CrossRef(53) % DescText    /'QSNOW               '/
DATA  Param_ID_CrossRef(53) % List_No     /grib_QSNOW_field/
DATA  Param_ID_CrossRef(53) % Table_No    /is_not_def,     &
                                           grbTblECMWFstd, &
                                           Is_not_def,     &
                                           is_not_def/

END MODULE rcf_GRIB_lookups_Mod
