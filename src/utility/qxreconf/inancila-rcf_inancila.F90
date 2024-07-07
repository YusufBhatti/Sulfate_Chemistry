! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Read in Headers from Ancillary Files
!
!  Subroutine inancila_rcf_inancila  - Read in Headers from Ancillary Files
!
! Description:
!   Read in the headers & lookuptables from ancillary files.
!
! Method:
!    For ancillary files that are required, the files are opened
!    and the headers and look-up tables read in.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
MODULE inancila_rcf_inancila_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INANCILA_RCF_INANCILA_MOD'

CONTAINS

SUBROUTINE inancila_rcf_inancila(len_fixhd,         &
                                 len_inthd,         &
                                 len_realhd,        &
                                 len1_levdepc,      &
                                 len2_levdepc,      &
                                 fixhd,             &
                                 inthd,             &
                                 realhd,            &
                                 lookup,            &
                                 a_realhd,          &
                                 a_levdepc,         &
                                 nlookups,          &
                                 lookup_start,      &
                                 len1_lookup,       &
                                 row_length,        &
                                 loc_row_length,    &
                                 p_rows,            &
                                 loc_p_rows,        &
                                 u_rows,            &
                                 r_row_length,      &
                                 r_rows,            &
                                 loc_r_row_length,  &
                                 loc_r_rows,        &
                                 p_levels,          &
                                 tr_levels,         &
                                 st_levels,         &
                                 sm_levels,         &
                                 ozone_levels,      &
                                 tpps_ozone_levels, &
                                 Ancil_Add,         &
                                 icode,cmessage)

USE ozone_inputs_mod, ONLY:                                  &
     zon_av_ozone, zon_av_tpps_ozone

USE um_stashcode_mod, ONLY:                                       &
    stashcode_prog_sec,           stashcode_soil_moist,           &
    stashcode_ozone,              stashcode_soil_temp,            &
    stashcode_frac_surf_type,     stashcode_snowdep_grd_tile,     &
    stashcode_snow_tile,          stashcode_snow_grnd,            &
    stashcode_tsurf_elev_surft,                                   &
    stashcode_rgrain,             stashcode_snowpack_bk_dens,     &
    stashcode_nsnow_layrs_tiles,  stashcode_snow_laythk_tiles,    &
    stashcode_snow_ice_tile,      stashcode_snow_liq_tile,        &
    stashcode_snow_t_tile,        stashcode_snow_laydns_tiles,    &
    stashcode_snow_grnsiz_tiles,                                  &
    stashcode_TppsOzone,          stashcode_LAI,                  &
    stashcode_canopy_height,      stashcode_total_aero_emiss,     &
    stashcode_total_aero,         stashcode_3d_nat_so2_em,        &
    stashcode_3d_ozone_mixrat,    stashcode_clim_biogenic_aero,   &
    stashcode_clim_delta_aero,    stashcode_ozone_tracer,         &
    stashcode_o3_colo3,           stashcode_riv_sequence,         &
    stashcode_riv_storage,        stashcode_lsm,                  &
    stashcode_mean_snow,          stashcode_icefrac


USE Ancil_mod, ONLY:                                             &
    AncF_UnitNo,                  levels,                         &
    nlookup,                      lookup_step,                    &
    ancil_requests,                num_ancil_files,                &
    ancil_files,                   num_ancil_requests

USE UM_ParCore, ONLY:                                            &
    mype

USE umPrintMgr

USE nlsizes_namelist_mod, ONLY:                                  &
    ntiles

USE nlstcall_mod, ONLY:                                          &
    lcal360

USE ancilcta_namelist_mod, ONLY:                                 &
    l_sstanom

USE near_equal_real_mod, ONLY:                                   &
    near_equal_real

USE Ereport_Mod, ONLY:                                           &
    Ereport

USE Rcf_HeadAddress_Mod, ONLY:                                   &
    SoilDepths,                                                   &
    FH_Dataset,                                                   &
    FH_Dataset_Ancil

USE filenamelength_mod, ONLY:  &
    filenamelength

USE jules_surface_types_mod
USE jules_surface_mod, ONLY: l_aggregate
USE jules_snow_mod, ONLY: nsmax
USE land_tile_ids,  ONLY: surface_type_ids, ml_snow_type_ids
USE max_dimensions, ONLY: ntype_max, snow_layers_max
USE clmchfcg_scenario_mod, ONLY: nsulpat
USE io
USE lookup_addresses


USE Rcf_Exppx_Mod, ONLY:      &
    Rcf_Exppx

USE Rcf_Ppx_Info_Mod, ONLY:    &
    STM_Record_Type

USE Rcf_Grid_Type_Mod, ONLY: &
    Output_Grid

USE Rcf_Level_Code_Mod, ONLY: &
    Rcf_Level_Code

USE missing_data_mod, ONLY: imdi

USE submodel_mod, ONLY: atmos_im

USE errormessagelength_mod, ONLY: errormessagelength

USE ancil_check_mod, ONLY:          &
    ancil_check_horizontal_grid,    &
    ancil_check_grid_stagger

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER :: Len_FixHd     ! Length of fixed header   ) in
INTEGER :: Len_IntHd     ! Length of Integer header ) ancillary
INTEGER :: Len_RealHd    ! Length of Real header    ) files
INTEGER :: Len1_LevDepC  ! ) First and second dimension of model
INTEGER :: Len2_LevDepC  ! ) level dependent constants array
INTEGER :: Len1_Lookup   ! First dimension for lookup tables.

INTEGER :: NLookups       ! Total no of lookup entries in
                          ! ancillary files to be read in.

!     Arrays for ancillary file headers
INTEGER :: FixHd(Len_FixHd,num_ancil_files)  ! Fixed H
INTEGER :: IntHd(Len_IntHd,num_ancil_files)  ! Integer H
REAL    :: RealHd(Len_RealHd,num_ancil_files)! Real H
INTEGER :: Lookup(Len1_Lookup,NLookups)      ! Lookup Table

!     Arrays for dump headers
REAL   :: A_RealHd(Len_RealHd)   ! Real Header
REAL   :: A_LevDepC(Len1_LevDepC,Len2_LevDepC)
                          ! Level dependent constants array

INTEGER :: Row_Length     ! Global row length
INTEGER :: Loc_Row_Length ! Local  row length

INTEGER :: P_Rows         ! Global no of p rows
INTEGER :: Loc_P_Rows     ! Local  no of p rows
INTEGER :: U_Rows         ! Global no of u rows
INTEGER :: R_Rows         ! Global no of r rows
INTEGER :: R_Row_Length   ! Global length of r rows
INTEGER :: loc_R_Rows     ! Local no of r rows
INTEGER :: loc_R_Row_Length ! Local length of r rows
INTEGER :: P_Levels       ! No of model levels
INTEGER :: TR_Levels      ! No of Tracer levels
INTEGER :: ST_Levels      ! No of Soil Temperature levels
INTEGER :: SM_Levels      ! No of Soil Moisture levels
INTEGER :: Ozone_Levels   ! No of Ozone levels
INTEGER ::tpps_ozone_levels
!                  No of ozone levels in tropopause-based ozone dataset

INTEGER :: Lookup_Start(num_ancil_files)
                          ! Pointer to first lookup entry for
                          ! each ancillary file
INTEGER :: Ancil_Add(num_ancil_requests)
                          ! Addresses in work space for anc data

INTEGER  :: Icode         ! Return Code

CHARACTER (LEN=errormessagelength) :: CMessage  !  Error Message if ICode > 0
CHARACTER (LEN=20) :: umFormat  !  Format code for umMessage

! Local variables
INTEGER, PARAMETER     :: filenameprovided=1
INTEGER :: i,j,j1,k,l,i_stash  ! Loop indices
INTEGER :: Irec          ! Loop over record
INTEGER :: Section
INTEGER :: Item
INTEGER :: stashcode
INTEGER :: Lookups, tLookups
INTEGER :: Start_Block
INTEGER :: Ozone_Row_Length
INTEGER :: tpps_ozone_row_length !! row length for tpps_ozone
INTEGER, PARAMETER :: Dummy = 1
INTEGER :: anc_file_no
INTEGER :: check_ps_levels
INTEGER :: expected_lbplev(ntype_max*snow_layers_max)

REAL :: LevDepC( (P_levels+1) * 4 )
REAL :: ColDepC(row_length+1),RowDepC(p_rows+1)

LOGICAL :: l_vert_mismatch
LOGICAL :: Check_Fail    ! flag for more complex checks.
LOGICAL :: l_soilm_ancil = .FALSE.  ! soil moisture from ancilary

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'INANCILA_RCF_INANCILA'
CHARACTER (LEN=filenamelength) ::  AncFileName  ! Anc file name

! STASHmaster entry
TYPE (STM_Record_Type), POINTER :: STM_Rec
INTEGER :: bot_level
INTEGER :: um_levels

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! CL Internal Structure

icode=0
cmessage=' '
!  1.  Initialisation for atmosphere model

!     Set no of levels for ancillary fields.
DO irec=1,num_ancil_requests

  SELECT CASE (ancil_requests(irec) % stashcode)
  CASE (stashcode_ozone)                    ! Ozone
    levels (irec) = ozone_levels
  CASE (stashcode_soil_temp)                ! Deep Soil Temps
    levels (irec) = st_levels
  CASE (stashcode_soil_moist)               ! Multi-layer hydrology
    levels (irec) = sm_levels
  CASE (stashcode_total_aero_emiss,    &
        stashcode_total_aero,          &
        stashcode_3d_nat_so2_em : stashcode_3d_ozone_mixrat, &
        321 : 340,  &    ! user ancil stashcodes
        stashcode_clim_biogenic_aero : stashcode_clim_delta_aero, &
        stashcode_ozone_tracer :stashcode_o3_colo3)
    levels (irec) = p_levels
    !          CASE (xxx)                ! not in ancilmaster
    !            levels (irec) = nsulpat
  CASE (stashcode_frac_surf_type)
    levels (irec) = ntype
  CASE (stashcode_LAI, stashcode_canopy_height)
    levels (irec) = npft
  CASE (stashcode_TppsOzone)       ! TppsOzone not in ancilmaster
    levels(irec) = tpps_ozone_levels
  CASE (stashcode_snowdep_grd_tile,              &
        stashcode_snowpack_bk_dens,              &
        stashcode_snow_tile,                     &
        stashcode_snow_grnd,                     &
        stashcode_rgrain,                        &
        stashcode_tsurf_elev_surft,              &
        stashcode_nsnow_layrs_tiles)
    levels(irec) = ntiles
  CASE (stashcode_snow_laythk_tiles,             &
        stashcode_snow_ice_tile,                 &
        stashcode_snow_liq_tile,                 &
        stashcode_snow_T_tile,                   &
        stashcode_snow_laydns_tiles,             &
        stashcode_snow_grnsiz_tiles)
    levels(irec) = ntiles*nsmax
  CASE DEFAULT                     ! Single level
    levels (irec) = 1
  END SELECT

END DO

! CL 1.4 Read headers

lookups=0

DO i=1,num_ancil_files

  ! C  Initialise LOOKUP_START (=0 implies file I not required)
  lookup_start(i)=0

  ! Open Ancillary file if required and read in headers & lookup table

  IF (mype == 0) THEN

    CALL umPrint('',src='inancila-rcf_inancila.F90')
    WRITE(umMessage,'(A,I5,A,A)')'Ancillary File ',i,' : ',&
          TRIM(ancil_files(i)%filename)
    CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')

  END IF

  ! Read headers for physical files required

  !      Open the ancillary file

  CALL File_Open (AncF_UnitNo, ancil_files(i)%filename,        &
    LEN_TRIM(ancil_files(i)%filename) ,0,filenameprovided,icode)

  IF (icode /= 0) THEN
    WRITE (Cmessage, '(A,A)') 'Problem opening Ancillary file ', &
       TRIM(ancil_files(i)%filename)
    CALL Ereport ( RoutineName, icode, Cmessage )
  END IF

  ! 1.4.1 Buffer in fixed length header record


  CALL Setpos (AncF_UnitNo, 0, icode)

  IF (icode /= 0) THEN
    WRITE (Cmessage, '(A,A)') 'Problem in SETPOS for Ancillary file ', &
        TRIM(ancil_files(i)%filename)
    CALL Ereport ( RoutineName, icode, Cmessage )
  END IF! Check icode

  !       Read in fixed header to get array dimensions
  ! DEPENDS ON: read_flh
  CALL read_flh(AncF_UnitNo,fixhd(1,i),len_fixhd,icode,cmessage)
  IF (icode >  0) THEN
    WRITE(umMessage,'(A,A)') ' Error in reading fixed header for file ',&
        TRIM(ancil_files(i)%filename)
    CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
    GO TO 9999   !  Return
  END IF

  ! C       Check for negative dimensions
  IF (fixhd(101,i) <= 0) fixhd(101,i)=1
  IF (fixhd(106,i) <= 0) fixhd(106,i)=1
  IF (fixhd(111,i) <= 0) fixhd(111,i)=1
  IF (fixhd(112,i) <= 0) fixhd(112,i)=1
  IF (fixhd(151,i) <= 0) fixhd(151,i)=1
  IF (fixhd(152,i) <= 0) fixhd(152,i)=1
  IF (fixhd(161,i) <= 0) fixhd(161,i)=1

  ! Check for valid Ancil file.
  IF (fixhd(FH_Dataset,i) /= FH_Dataset_Ancil) THEN
    WRITE(CMessage,'(A,A)') 'Invalid fixed header for ancillary file ',&
        TRIM(ancil_files(i)%filename)

    icode = 10
    CALL Ereport ( RoutineName, icode, Cmessage )
  END IF

  CALL ancil_check_grid_stagger(model_grid_stagger = output_grid%grid_stagger,&
                                ancil_grid_stagger = fixhd(9,i),              &
                                ancil_file = ancil_files(i))

  ! C Set start position of boundary fields for file
  lookup_start(i)=lookups+1

  IF (lookups+fixhd(152,i) >  nlookups) THEN
    WRITE(umMessage,'(A,A)')                                           &
                'No room in LOOKUP table for Ancillary File ',       &
                TRIM(ancil_files(i)%filename)
    CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
    cmessage='inancila_rcf_inancila: Insufficient space for LOOKUP headers'
    icode=14
    GO TO 9999   !  Return
  END IF


  CALL setpos(AncF_UnitNo, 0, icode)

  IF (icode /= 0) THEN
    WRITE (Cmessage, '(3A)')                         &
      'Problem in SETPOS for Ancillary file ',             &
    TRIM(ancil_files(i)%filename),' before READHEAD.'
    CALL Ereport ( RoutineName, icode, Cmessage )
  END IF! Check icode


  ! DEPENDS ON: readhead
  IF (fixhd(116,i) > 1) THEN
    CALL readhead(AncF_UnitNo,                                      &
                fixhd(1,i),len_fixhd,                             &
                inthd(1,i),fixhd(101,i),                          &
                realhd(1,i),fixhd(106,i),                         &
                levdepc,fixhd(111,i),fixhd(112,i),                &
                rowdepc,fixhd(116,i),fixhd(117,i),                &
                coldepc,fixhd(121,i),fixhd(122,i),                &
                dummy,dummy,dummy,                                &
                dummy,dummy,                                      &
                dummy,dummy,                                      &
                dummy,dummy,                                      &
                dummy,dummy,                                      &
                dummy,dummy,                                      &
                lookup(1,lookups+1),fixhd(151,i),fixhd(152,i),    &
                fixhd(161,i),                                     &
                start_block,icode,cmessage)
  ELSE
    CALL readhead(AncF_UnitNo,                                      &
                fixhd(1,i),len_fixhd,                             &
                inthd(1,i),fixhd(101,i),                          &
                realhd(1,i),fixhd(106,i),                         &
                levdepc,fixhd(111,i),fixhd(112,i),                &
                dummy,dummy,dummy,                                &
                dummy,dummy,dummy,                                &
                dummy,dummy,dummy,                                &
                dummy,dummy,                                      &
                dummy,dummy,                                      &
                dummy,dummy,                                      &
                dummy,dummy,                                      &
                dummy,dummy,                                      &
                lookup(1,lookups+1),fixhd(151,i),fixhd(152,i),    &
                fixhd(161,i),                                     &
                start_block,icode,cmessage)
  END IF

  IF (icode >  0) THEN
    WRITE(umMessage,'(A,A)') 'ERROR in READHEAD for Ancillary File ',  &
        TRIM(ancil_files(i)%filename)
    CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
    GO TO 9999   !   Return
  END IF

  !      Close the ancillary file

  CALL File_Close (AncF_UnitNo,ancil_files(i)%filename ,         &
    LEN_TRIM(ancil_files(i)%filename) ,filenameprovided,icode)

  IF (icode /= 0) THEN
    WRITE (Cmessage, '(A,A)') 'Problem closing Ancillary file ', &
                   TRIM(ancil_files(i)%filename  )
    CALL Ereport ( RoutineName, icode, Cmessage )
  END IF

  ! Flexible tiles: Check tile configuration in ancillaries with tiles
  IF ( .NOT. l_aggregate ) THEN
    DO i_stash = 1, ancil_files(i) % num_stash
      ! Need to check all stashcodes contained in the file
      stashcode = ancil_files(i) % stashcodes(i_stash)
      section   = stashcode / 1000
      item      = stashcode - (section * 1000)
      stm_rec => rcf_exppx( atmos_im, section, item )

      ! Check for tiles and determine pseudo level type
      SELECT CASE ( stm_rec % pl_code )
      CASE (7)
        check_ps_levels = ntype
      CASE (8)
        check_ps_levels = npft
      CASE (9)
        check_ps_levels = ntiles
      CASE (11)
        check_ps_levels = nsmax * ntiles
      CASE DEFAULT
        check_ps_levels = imdi
      END SELECT

      !   (*_type_ids set in rcf_read_namelists)
      IF (stm_rec % pl_code == 11) THEN
        expected_lbplev(1:check_ps_levels) = ml_snow_type_ids(1:check_ps_levels)
      ELSE
        expected_lbplev(1:check_ps_levels) = surface_type_ids(1:check_ps_levels)
      END IF

      ! Only check if check_ps_levels has been set
      IF ( check_ps_levels > 0 ) THEN

        ! Find first occurrence of stash item in lookup table,
        !  label previous item as tlookups so the set of pseudo-levels
        !  associated with the stashcode are in index positions tlookup+1 to
        !  tlookup+check_ps_levels
        tlookups = -1
        DO j=lookups+1, lookups+fixhd(152,i)
          IF (lookup(item_code,j) == stashcode) THEN
            tlookups = j-1
            EXIT
          END IF
        END DO
        IF (tlookups < 0) THEN
          WRITE(cmessage,'(A,I0,A)') "Expected surface field stash ", item,   &
             " was not found in ancillary: " // newline //                    &
             ancil_files(i)%filename
          icode = 98000 + item
          CALL ereport ( routinename, icode, cmessage )
        END IF

        ! Assume all remaining fields with given stashcode are in pseudo-level
        ! order after the first one.
        DO j = 1, check_ps_levels
          IF ( expected_lbplev(j) /= lookup(lbplev, tlookups+j) ) THEN
            WRITE(umMessage,'(A,A)')                                   &
               'Land surface configuration does not match ',           &
               TRIM( ancil_files(i)%filename )
            CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
            WRITE(umFormat,'(A,I0,A)') '(A23,',check_ps_levels,'(1X,I6))'
            WRITE(umMessage,umFormat) 'Ancillary configuration',       &
               lookup(lbplev,tlookups+1:tlookups+check_ps_levels)
            CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
            WRITE(umMessage,umFormat) 'Namelist  configuration',       &
               expected_lbplev(1:check_ps_levels)
            CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
            cmessage = 'Land surface configuration does not match ancillary'
            icode = 99000 + item
            CALL ereport ( routinename, icode, cmessage )
          END IF
        END DO

      END IF
    END DO
  END IF

  !   Check calendar indicator is correct if this is set in the ancillary
  IF (fixhd(8,i) /=  imdi) THEN
    IF ((     lcal360 .AND. fixhd(8,i) /= 2) .OR.                  &
        (.NOT. lcal360 .AND. fixhd(8,i) /= 1) ) THEN
      icode = 100 + i
      cmessage = 'inancila_rcf_inancila : '                        &
          // 'Wrong calendar set in Ancillary File'
      IF (mype == 0) THEN
        CALL umPrint('  ******** Error in inancila_rcf_inancila ********',&
            src='inancila-rcf_inancila.F90')
        WRITE(umMessage,'(A,A)')                                        &
            '  Wrong calendar setting in Ancillary File ',               &
            TRIM(ancil_files(i)%filename )
        CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
        IF (lcal360) THEN
          CALL umPrint(                                          &
             '  Model run is set up for 360 day calendar.',     &
             src='inancila-rcf_inancila.F90')
          CALL umPrint(                                          &
             '  Ancillary File is for 365 day calendar.',       &
             src='inancila-rcf_inancila.F90')
        ELSE
          CALL umPrint(                                          &
             '  Model run is set up for 365 day calendar.',     &
             src='inancila-rcf_inancila.F90')
          CALL umPrint(                                          &
             '  Ancillary File is for 360 day calendar.',       &
             src='inancila-rcf_inancila.F90')
        END IF
        CALL umPrint( '  Rerun with correct ancillary file.',    &
           src='inancila-rcf_inancila.F90')
      END IF
      GO TO 9999   !  Return
    END IF
  ELSE
    IF (mype == 0) THEN
      CALL umPrint(' Unspecified calendar type in ancillary file',&
         src='inancila-rcf_inancila.F90')
    END IF
  END IF


  ! Check validity of Ancil Header array sizes
  IF (fixhd(101,i) /=  len_inthd) THEN
    icode = 17
    umMessage = 'inancila_rcf_inancila: Invalid Integer Header length'
    CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
    WRITE(umMessage,'(2A)') 'Error reading ancil file : ', &
                              TRIM(ancil_files(i)%filename)
    CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
    WRITE(umMessage,'(2(A,1X,I3))') ' FLH(101)  : ',fixhd(101,i), &
                                    ', len_inthd : ',len_inthd
    CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
    cmessage = 'Integer Header length in file does not match expected value'
    GO TO 9999
  END IF
  IF (fixhd(106,i) /=  len_realhd) THEN
    icode = 18
    umMessage = 'inancila_rcf_inancila: Invalid Real Header length'
    CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
    WRITE(umMessage,'(2A)') 'Error reading ancil file :- ', &
                              TRIM(ancil_files(i)%filename)
    CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
    WRITE(umMessage,'(2(A,1X,I3))') ' FLH(106)   : ',fixhd(106,i), &
                                    ', len_realhd : ',len_realhd
    CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
    cmessage = 'Real Header length in file does not match expected value'
    GO TO 9999
  END IF


  ! 1.4.2 Check horizontal grid
  !
  ! Value of integer `lookups` keeps track of the index (in the master 
  ! lookup array) of the last lookup of the previous ancil file to be read.  
  ! Pass in slice of the master lookup array which corresponds to the lookups 
  ! for the current ancil file.
  CALL ancil_check_horizontal_grid(                                 &
           lookup(:, lookups+1:lookups+fixhd(152,i)),               &
           ancil_files(i), len1_lookup,                             &
           output_grid%glob_p_rows,  output_grid%glob_p_row_length, &
           output_grid%glob_u_rows,  output_grid%glob_u_row_length, &
           output_grid%glob_v_rows,  output_grid%glob_v_row_length, &
           output_grid%glob_r_rows,  output_grid%glob_r_row_length)

  ! CL 1.4.3 Buffer in real constants

  IF (fixhd(105,i) > 0 .AND.                                     &
        fixhd(117,i) < 1 .AND. fixhd(122,i) < 1) THEN

    ! C Check validity of real header and print out information

    ! Only check if ancillary and model is on C grid (with P at poles)
    IF (fixhd(9,i) == 3 .AND. output_grid % grid_stagger == 3) THEN
      DO j=1,6
        IF (near_equal_real(realhd(j,i),a_realhd(j))) THEN
           DO i_stash = 1, ancil_files(i)%num_stash
             stashcode = ancil_files(i)%stashcodes(i_stash)
             IF (stashcode /=  stashcode_ozone         &
                  .OR. (j /= 1 .AND. j /= 4)) THEN
               WRITE(umMessage,'(A)') ' Inconsistency in Real Headers.'
               CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
               WRITE(umMessage,'(A)') ' Real Header Values.'
               CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
               DO k=1,6
                 WRITE(umMessage,*) k,' Anc File ',realhd(k,i),     &
                      ' Model Dump ',a_realhd(k)
                 CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
               END DO
               icode=8
               cmessage='inancila_rcf_inancila: REAL header Error.'
               GO TO 9999
             END IF
           END DO ! End loop over stashcodes in file
        END IF
      END DO
    END IF

  END IF

  ! CL 1.4.4 Buffer in level dependent constants if required
  ! C        Not retained in model after initial check

  IF (fixhd(110,i) >  0) THEN

    ! The following stash items may have multi-level data
    !  Ozone, tpps_ozone, soil temperature  should contain multi-level data.
    !  Soil moisture, snow depth, fractional snow time and soil moisture in
    !  layers may possibly also have multi level data.
    !  Aerosols, murkiness, user ancils may also have multi level data.
   DO i_stash = 1, ancil_files(i)%num_stash
     ! Need to check all stashcodes contained in the file
     ! so loop over all stash in file
     stashcode = ancil_files(i)%stashcodes(i_stash)
    IF ( stashcode == stashcode_ozone  .OR.   &
         stashcode == stashcode_total_aero_emiss  .OR.  &
         stashcode == stashcode_total_aero  .OR.        &
         (stashcode > 320   .AND.             &
          stashcode < 341)  .OR.               &
                            !  User Ancillary fields
         (stashcode > 479   .AND.             &
         stashcode < 488) ) THEN
                           !  Cariolle Ozone

      ! Check that ancillary file is set up for correct vertical levels

      IF (fixhd(111,i)-1 /= p_levels) THEN
        icode=110
        WRITE (cmessage, '(A, I4, A, I4)')                          &
        'Ancillary File set up for wrong no of model levels. Anc ', &
        fixhd(111,i)-1, ' Model ',p_levels
        CALL Ereport ( RoutineName, icode, Cmessage )
      END IF

      l_vert_mismatch = .FALSE.

      ! Check eta_theta and eta_rho

      DO j=1,p_levels+1
        IF (near_equal_real( levdepc(j), a_levdepc(j,1) )) THEN
          l_vert_mismatch = .TRUE.
          EXIT
        END IF
      END DO

      DO j=1,p_levels
        IF (near_equal_real( levdepc(fixhd(111,i)+j), a_levdepc(j,2) )) THEN
          l_vert_mismatch = .TRUE.
          EXIT
        END IF
      END DO

      ! Abort if there is a mis-match

      IF (l_vert_mismatch) THEN
        WRITE(umMessage,'(A)')                                     &
          'Mismatch in vertical levels between model and Ancillary File.'
        CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
        WRITE(umMessage,'(A)')                                     &
           'Anc File : ',TRIM(ancil_files(i)%filename)
        CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
        WRITE(umMessage,'(A)') 'Eta_Theta - Model'
        CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
        WRITE(umMessage,'(5F10.7)') (a_levdepc(k,1),k=1,p_levels+1)
        CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
        WRITE(umMessage,'(A)') 'Eta_Theta - Anc File'
        CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
        WRITE(umMessage,'(5F10.7)') (levdepc(k),k=1,p_levels+1)
        CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
        WRITE(umMessage,'(A)') 'Eta_Rho   - Model'
        CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
        WRITE(umMessage,'(5F10.7)') (a_levdepc(k,2),k=1,p_levels)
        CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
        WRITE(umMessage,'(A)') 'Eta_Rho   - Anc File'
        CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
        WRITE(umMessage,'(5F10.7)') (levdepc(p_levels+1+k),k=1,p_levels)
        CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
        icode=11
        WRITE (cmessage, '(2A)') 'Mismatch in LEVDEPC ',          &
        'between model and Ancillary File.'
        CALL Ereport ( RoutineName, icode, Cmessage )
      END IF

      !! tropopause-based ozone   NOT in ancil master
      !           ELSE IF (anc_file(I) % anc_file_number  ==  25) THEN
                   !! no checks to run

                 !   Soil Moisture File
    ELSE IF (stashcode == stashcode_mean_snow .OR.    &
              stashcode == stashcode_soil_moist .OR.  &
              stashcode == stashcode_snowdep_grd_tile) THEN
      ! Check Soil Moisture levels

      IF (PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
        CALL umPrint('',src='inancila-rcf_inancila.F90')
        WRITE(umMessage,'(A,I3)') 'SoilDepths = ',SoilDepths
        CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
        WRITE(umMessage,'(A,I3)') 'SM_Levels  = ',sm_levels
        CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
        DO j=1,sm_levels
          WRITE(umMessage,*) 'model ',a_levdepc(j,SoilDepths),       &
          ' anc ',levdepc(fixhd(111,i)*3+j)
          CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
        END DO
      END IF

      ! The smow ancillary file can be used for snow but not soil moisture
      ! so only check soil moisture if it is required.
      DO k=1,num_ancil_requests
        IF ( ancil_requests(k)%stashcode == stashcode_soil_moist ) THEN
          l_soilm_ancil=.TRUE.
        END IF
      END DO

      DO j=1,sm_levels
        IF (near_equal_real(levdepc(fixhd(111,i)*3+j),            &
                     a_levdepc(j,SoilDepths))                     &
                     .AND. l_soilm_ancil) THEN
          icode=12
          cmessage='inancila_rcf_inancila: error in LEVDEPC.'
          GO TO 9999
        END IF
      END DO

      !   Deep Soil Temperature File
    ELSE IF (stashcode == stashcode_soil_temp) THEN
      IF (PrintStatus >= PrStatus_Diag .AND. mype == 0) THEN
        CALL umPrint('',src='inancila-rcf_inancila.F90')
        WRITE(umMessage,'(A,I3)') 'SoilDepths = ',SoilDepths
        CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
        WRITE(umMessage,'(A,I3)') 'st_levels  = ',st_levels
        CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
        DO j=1,st_levels
          WRITE(umMessage,*) 'model ',a_levdepc(j,SoilDepths),         &
          ' anc ',levdepc(fixhd(111,i)*3+j)
          CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
        END DO
      END IF

      DO j=1,st_levels
        IF (near_equal_real(levdepc(fixhd(111,i)*3+j),                 &
                    a_levdepc(j,SoilDepths))) THEN
          icode=12
          cmessage='inancila_rcf_inancila: error in LEVDEPC.'
          GO TO 9999
        END IF
      END DO

    END IF  !  If I
  END DO

  END IF  !  If Fixhd(110,I) > 0

  !        Add on no of lookup entries
  lookups=lookups+fixhd(152,i)

END DO    ! Loop over ancillary files (I)

! CL 1.5 Set positions in main data blocks

item=1
DO i=1,num_ancil_requests
  Ancil_add(i)=item

  ! Find whether we need to handle the extra level due to ancillary not
  ! supporting theta level 0.
  STM_Rec => Rcf_Exppx( 1, ancil_requests(i) % section, &
       ancil_requests(i) % item  )

  IF ( STM_Rec % lb_code > 0) THEN
    CALL Rcf_Level_Code( STM_Rec % lb_code, bot_level, Output_Grid )
  ELSE
    bot_level = 1
  END IF

  ! Force bottom level to be 1 (if not 0)
  IF (bot_level /= 0) THEN
    bot_level = 1
  END IF
  ! um_levels is the number of levels in the UM (not in the ancillary)
  um_levels = levels(i)-bot_level+1

  SELECT CASE ( ancil_requests(i) % stashcode )
  CASE (stashcode_ozone)        ! Ozone
    IF ( zon_av_ozone ) THEN
      item = item + loc_p_rows * um_levels
    ELSE
      item = item + loc_row_length * loc_p_rows * um_levels
    END IF

  CASE (stashcode_TppsOzone)       ! tropopause-based ozone
    IF ( zon_av_tpps_ozone ) THEN
      item = item + loc_p_rows * um_levels
    ELSE
      item = item + loc_row_length * loc_p_rows * um_levels
    END IF

  CASE (stashcode_riv_sequence:stashcode_riv_storage) ! River Routing
    item = item + loc_R_Row_Length * loc_R_Rows * um_levels

  CASE DEFAULT ! All 'normal' ancil fields.
    item = item + loc_row_length * loc_p_rows * um_levels

  END SELECT
END DO

! CL 1.6 Set positions of data
DO i=1,num_ancil_requests
  anc_file_no = ancil_requests(i)%ancil_file_num

  nlookup(i) =0
  lookup_step(i)=0

  ! C If LOOKUP_START=0 for anc_file_no, no fields required.
  IF (lookup_start(anc_file_no) >  0) THEN

    IF (  PrintStatus > PrStatus_Diag .AND. mype == 0 ) THEN
      WRITE(umMessage,'(A,I5)') ' lookup_start > 0 for stashcode ',  &
           ancil_requests(i)%stashcode
      CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
    END IF

    DO j=lookup_start(anc_file_no),lookups
      ! Find first occurrence of stash item in lookup table
      IF (lookup(item_code,j) == ancil_requests(i) % stashcode) THEN
        nlookup(i)=j-lookup_start(anc_file_no)+1
        EXIT
      END IF

    END DO

    ! C Find second occurence of data to set LOOKUP_STEP

    lookup_step(i)=0

    IF (j <  lookups) THEN

      DO j1=j+levels(i),lookups
        IF (lookup(item_code,j1) == ancil_requests(i) % stashcode) THEN
          lookup_step(i)=                                         &
               j1-nlookup(i)-lookup_start(anc_file_no)+1
          EXIT
        END IF
      END DO

    END IF
    ! Check Ancillary files are consistent since we assume:
    ! for month in month_list:
    !   for STASH in STASH_list:
    !     for level in level_psuedolevel_list:
    ! First check if its periodic since we can only check these at
    ! this point and whether it has multiple fields defined.
    IF (fixhd(10,anc_file_no) == 2 .AND. lookup_step(i) > 0) THEN
      ! Loop over the fields which should only differ by time.
      DO j = lookup_start(anc_file_no) + nlookup(i)-1,               &
           lookup_start(anc_file_no) + fixhd(152,anc_file_no) - 1, &
           lookup_step(i)
        ! If the STASH is different we have a problem.
        IF (lookup(item_code,j) /= ancil_requests(i) % stashcode ) THEN
          icode = 53
          WRITE (cmessage,'(A,A)')                                  &
               'Incorrect structure for ancillary file ',               &
               TRIM(ancil_files(anc_file_no)%filename)
          CALL ereport ( routinename, icode, cmessage )
          EXIT
        END IF
      END DO
    END IF
  END IF

END DO

! SET LEVELS=2 FOR ICE FRACTION AND SNOW DEPTH, TO INDICATE PRESCENCE
! fractional time fields

DO i=1,num_ancil_requests
  IF (ancil_requests(i)%stashcode == stashcode_mean_snow ) THEN
    levels(i)=2      !   Snow Depth
  END IF
  IF (ancil_requests(i)%stashcode == stashcode_icefrac ) THEN
    levels(i)=2           !   Ice Fraction
  END IF
END DO


IF (  PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN

  CALL umPrint('',src='inancila-rcf_inancila.F90')
  CALL umPrint(' Summary from inancila_rcf_inancila ', &
       src='inancila-rcf_inancila.F90')

  CALL umPrint(' LOOKUP_START ',src='inancila-rcf_inancila.F90')
  DO i=1,num_ancil_files
    WRITE(umMessage,'(A,I7)')'lookup_start ',lookup_start(i)
    CALL umPrint(umMessage,src='inancila-rcf_inancila.F90')
  END DO

  CALL umPrint(' Ancil Ref Numbers ',src='inancila-rcf_inancila.F90')
  DO i=1,num_ancil_files
    CALL umPrint(TRIM(str(i)), &
         src='inancila-rcf_inancila.F90')
  END DO

  CALL umPrint(' LOOKUP_STEP ',src='inancila-rcf_inancila.F90')
  DO i=1,SIZE(lookup_step)
    CALL umPrint(TRIM(str(lookup_step(i))),src='inancila-rcf_inancila.F90')
  END DO

  CALL umPrint(' NLOOKUP ',src='inancila-rcf_inancila.F90')
  DO i=1,SIZE(nlookup)
    CALL umPrint(TRIM(str(nlookup(i))),src='inancila-rcf_inancila.F90')
  END DO

  CALL umPrint(' LEVELS ',src='inancila-rcf_inancila.F90')
  DO i=1,SIZE(lookup_step)
    CALL umPrint(TRIM(str(levels(i))),src='inancila-rcf_inancila.F90')
  END DO

  CALL umPrint(' ANCIL_ADD ',src='inancila-rcf_inancila.F90')
  DO i=1,SIZE(lookup_step)
    CALL umPrint(TRIM(str(ancil_add(i))),src='inancila-rcf_inancila.F90')
  END DO
  CALL umPrint('',src='inancila-rcf_inancila.F90')

END IF

9999 CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE inancila_rcf_inancila

END MODULE inancila_rcf_inancila_mod
