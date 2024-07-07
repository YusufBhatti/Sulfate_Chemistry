! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Ancillary processing for the Atmosphere model

MODULE Rcf_Ancil_Atmos_Mod

!  Subroutine Rcf_Ancil_Atmos  - Ancillary processing for Atmosphere
!
! Description:
!    Controls all ancillary processing for the atmosphere model
!
! Method:
!    1. Reads in the ancilmaster records to get information on the
!       ancillary fields and files.
!    2. Determines workspace required for ancillary lookups and data.
!    3. Calls inancila_rcf_inancila to read in ancillary lookups.
!    4. Calls replanca_rcf_replanca to read in ancillary data.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_ANCIL_ATMOS_MOD'

CONTAINS

SUBROUTINE Rcf_Ancil_Atmos ( Hdr_In, Hdr_Out,                           &
                             Fields_In, Field_Count_In,                 &
                             Fields_Out, Field_Count_Out, data_source )

USE Ancil_mod, ONLY: &
    nlookup,             lookup_step,         &
    levels,              ancil_add,           &
    ancil_requests,       num_ancil_requests,  &
    num_ancil_files

USE inancila_rcf_inancila_mod, ONLY: inancila_rcf_inancila

USE calc_nlookups_Mod, ONLY: &
    calc_nlookups

USE Rcf_calc_len_ancil_Mod, ONLY: &
    Rcf_calc_len_ancil

USE rcf_interpolate_Mod, ONLY: &
    rcf_interpolate

USE Rcf_Items_Mod, ONLY:       &
    Num_Items

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    newline,                &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Diag

USE submodel_mod, ONLY: atmos_im

USE Rcf_Grid_Type_Mod, ONLY:   &
    Input_Grid,             &
    Output_Grid

USE Rcf_HeadAddress_Mod, ONLY: &
    fh_vtyear,              &
    fh_vtmonth,             &
    fh_vtday,               &
    fh_vthour,              &
    fh_vtminute,            &
    fh_vtsecond,            &
    rc_latspacing,          &
    rc_firstlat

USE Rcf_UMhead_Mod, ONLY:  &
    Um_header_type,         &
    LenFixHd

USE Rcf_Lsm_Mod, ONLY:     &
    glob_lsm_out,           &
    local_lsm_out,          &
    glob_land_out,          &
    local_land_out

USE decomp_params, ONLY:    &
    Decomp_rcf_input,       &
    Decomp_rcf_output

USE Rcf_Alloc_Field_mod, ONLY:  &
    Rcf_Alloc_Field,             &
    Rcf_Dealloc_Field

USE rcf_read_field_mod, ONLY: &
    Rcf_Read_Field

USE rcf_write_field_mod, ONLY: &
    Rcf_Write_Field

USE Rcf_Locate_mod, ONLY: &
    Rcf_Locate

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE UM_ParVars, ONLY:      &
    current_decomp_type,    &
    g_datastart,            &
    change_decomposition

USE UM_ParCore, ONLY: &
    mype

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_Set_Interp_Flags_Mod, ONLY: &
    interp_h_only,                   &
    interp_copy

USE Rcf_Interp_Weights_Mod, ONLY: &
    h_int_active

USE um_stashcode_mod, ONLY: &
    stashcode_tstar,           &
    stashcode_land_frac,       &
    stashcode_tstar_land,      &
    stashcode_tstar_sea,       &
    stashcode_tstar_sice,      &
    stashcode_icefrac,         &
    stashcode_tstar_anom,      &
    stashcode_prog_sec,        &
    stashcode_surf_z_curr,     &
    stashcode_surf_m_curr,     &
    stashcode_lsm

USE Rcf_Exppx_Mod, ONLY: &
    Rcf_Exppx

USE jules_sea_seaice_mod, ONLY: l_ctile

USE nlsizes_namelist_mod, ONLY: &
    tpps_ozone_levels

USE ancilcta_namelist_mod, ONLY:&
    l_sstanom

USE mask_compression, ONLY: compress_to_mask

USE cppxref_mod, ONLY: ppx_atm_compressed

USE rcf_data_source_mod, ONLY: &
    data_source_type,          &
    already_processed

USE items_nml_mod, ONLY:       &
    input_dump,                &
    ancillary_file,            &
    set_to_mdi

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments

TYPE (um_header_type),      INTENT(IN) :: hdr_in
TYPE (um_header_type),      INTENT(IN) :: hdr_out
TYPE (field_type), POINTER             :: fields_in (:)
TYPE (field_type), POINTER             :: fields_out (:)
INTEGER ,                   INTENT(IN) :: field_count_in
INTEGER ,                   INTENT(IN) :: Field_count_out
TYPE (data_source_type), POINTER       :: data_source(:)


! ----------------------------------------------------------------
! Arrays to store headers from ancillary files
! (Could go into a module & add USE to inancila_rcf_inancila & replanca_rcf_replanca ?)

! Length of fixed header is defined in umhead

INTEGER, PARAMETER :: LenInthd_anc  = 15
INTEGER, PARAMETER :: LenRealhd_anc = 6

! The following arrays are 2-dimensional, the second dimension
! being the number of atmosphere ancillary files.

INTEGER, ALLOCATABLE :: fixhd_ancil(:,:)
INTEGER, ALLOCATABLE :: inthd_ancil(:,:)
REAL   , ALLOCATABLE :: realhd_ancil(:,:)
INTEGER, ALLOCATABLE :: lookup_ancil(:,:)

! ----------------------------------------------------------------
! Variables & Work arrays for ancillary processing

INTEGER                            :: nlookups
INTEGER                            :: len_ancil
INTEGER                            :: ipt
REAL,    ALLOCATABLE :: ancil_data(:)
REAL,    ALLOCATABLE :: ancil_land_temp(:)
INTEGER, ALLOCATABLE :: lookup_start(:)

! ----------------------------------------------------------------
! Local Variables

INTEGER :: pos_IceFrac_In       ! Position of Ice Fraction in Input dump
INTEGER :: pos_IceFrac_Out      ! Position of Ice Fraction in Output dump
INTEGER :: pos_Tstar_In         ! Position of T Star in Input dump
INTEGER :: pos_Tstar_Out        ! Position of T Star in Output dump
INTEGER :: pos_Land_Frac_In     ! Position of land frac in Input dump
INTEGER :: pos_Land_Frac_Out    ! Position of land frac in Output dump
INTEGER :: pos_Tstar_Land_In    ! Position of T Star_Land in Input dump
INTEGER :: pos_Tstar_Land_Out   ! Position of T Star_Land in Output
INTEGER :: pos_Tstar_Sea_In     ! Position of T Star_Sea in Input
INTEGER :: pos_Tstar_Sea_Out    ! Position of T Star_Sea in Output
INTEGER :: pos_Tstar_Sice_In    ! Position of T Star_Sice in Input
INTEGER :: pos_Tstar_Sice_Out   ! Position of T Star_Sice in Output
INTEGER :: pos_Tstar_Anom       ! Position of T Star Anomaly in Input dump
INTEGER :: pos


TYPE (field_type), POINTER   :: IceFrac_In
TYPE (field_type), POINTER   :: TStar_In
TYPE (field_type), POINTER   :: Land_Frac_In
TYPE (field_type), POINTER   :: TStar_Land_In
TYPE (field_type), POINTER   :: TStar_Sea_In
TYPE (field_type), POINTER   :: TStar_Sice_In
TYPE (field_type), POINTER   :: IceFrac_Out
TYPE (field_type), POINTER   :: TStar_Out
TYPE (field_type), POINTER   :: Land_Frac_Out
TYPE (field_type), POINTER   :: TStar_Land_Out
TYPE (field_type), POINTER   :: TStar_Sea_Out
TYPE (field_type), POINTER   :: TStar_Sice_Out
TYPE (field_type), POINTER   :: TStar_Anom

TYPE (field_type)            :: dummy
TYPE (field_type)            :: dummy_anc
TYPE (field_type), TARGET    :: dummy_IceFr
TYPE (field_type), TARGET    :: dummy_TStar
TYPE (field_type), TARGET    :: dummy_Land_Frac
TYPE (field_type), TARGET    :: dummy_TStar_Land
TYPE (field_type), TARGET    :: dummy_TStar_Sea
TYPE (field_type), TARGET    :: dummy_TStar_Sice
TYPE (field_type), TARGET    :: dummy_TStar_Anom

INTEGER      :: i,j,k          !  Loop indices
INTEGER      :: i_full
INTEGER      :: i_land
INTEGER      :: irec
INTEGER      :: ianc_Mask
INTEGER      :: ianc_IceFrac
INTEGER      :: ianc_TStar
INTEGER      :: ianc_u_curr
INTEGER      :: ianc_v_curr
INTEGER      :: ErrorStatus
INTEGER      :: orig_decomp

LOGICAL      :: multi_pe_in       ! flag to confirm 'gather' req'd

CHARACTER (LEN=errormessagelength)           :: Cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_ANCIL_ATMOS'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! ----------------------------------------------------------------

! ===============================================
! Set up list of ancillary fields to be read in
! This may  erroreport if some dependencises are missing.

DO i=1,num_ancil_requests
  ! set up stash names for each item requuested.
  dummy_anc % stashmaster => rcf_exppx (atmos_im,  &
                                       ancil_requests(i)%section, &
                                       ancil_requests(i)%item )

  ancil_requests(i)%stash_name = dummy_anc % stashmaster % NAME

  ! Check if both surface currents are selected
  IF (ancil_requests(i)%stashcode == stashcode_surf_z_curr .AND.          &
       ALL(ancil_requests%stashcode /= stashcode_surf_m_curr ) ) THEN
    WRITE (cmessage,'(A)')                                                &
      'Ancil configure request for Stashcode 28 (SURFACE ZONAL CURRENT '  &
      // newline //                                                       &
      'AFTER TIMESTEP) but no request made for Stashcode 29 (SURFACE'     &
      // newline //                                                       &
      'MERID CURRENT AFTER TIMESTEP). If one field is configured from'    &
      // newline //                                                       &
      'ancil the othershould also be configured from ancil.'              &
      // newline // newline //                                            &
      'Please use items namelist to configure stashcode 29 from ancil.'
    ErrorStatus = 30
    CALL Ereport ( RoutineName, ErrorStatus, cmessage)
  ELSE IF (ancil_requests(i)%stashcode == stashcode_surf_m_curr .AND.     &
      ALL(ancil_requests%stashcode /= stashcode_surf_z_curr ) ) THEN
    WRITE (cmessage,'(A)')                                                &
      'Ancil configure request for Stashcode 29 (SURFACE MERID CURRENT '  &
      // newline //                                                       &
      'AFTER TIMESTEP) but no request made for Stashcode 28 (SURFACE'     &
      // newline //                                                       &
      'ZONAL CURRENT AFTER TIMESTEP). If one field is configured from'    &
      // newline //                                                       &
      'ancil the other field should also be configured from ancil.'       &
      // newline // newline //                                            &
      'Please use items namelist to configure stashcode 28 from ancil.'
    ErrorStatus = 32
    CALL Ereport ( RoutineName, ErrorStatus, cmessage)
  END IF

  ! Sea surface temperature must be updated when sea ice is updated
  ! If sst ancil not already provide to be read in then namelist input has no
  ! means of identifying the require ancillary.

  IF (ancil_requests(i)%stashcode == stashcode_icefrac .AND.     &
      ALL(ancil_requests%stashcode /= stashcode_tstar ) ) THEN
    WRITE (cmessage,'(A)')                                  &
      'User needs to supply and configure SST ancillary file.'
    ErrorStatus = 35
    CALL Ereport ( RoutineName, ErrorStatus, cmessage)
  END IF

  ! Sea surface temperature anomaly switches on climatological sst
  IF ( l_sstanom .AND. ALL(ancil_requests%stashcode /= stashcode_tstar ) ) THEN
    WRITE (cmessage,'(A)')                                  &
      'User needs to supply and configure SST ancillary file.'
    ErrorStatus = 36
    CALL Ereport ( RoutineName, ErrorStatus, cmessage)
  END IF

END DO

IF ( PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
  IF (num_ancil_requests > 0 ) THEN
    WRITE(umMessage,*)
    CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
    WRITE(umMessage,'(I5,A)') num_ancil_requests , &
                         ' Ancillary fields to be read in:'
    CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
    DO i=1,num_ancil_requests
      WRITE(umMessage,'(I5,A)') ancil_requests(i)%stashcode,     &
                           ancil_requests(i)%stash_name
      CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
    END DO
  END IF
END IF


! Determine the number of lookup entries to be read in from ancillary files
CALL calc_nlookups (nlookups)


IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
  WRITE(umMessage,'(A,I0)') ' rcf_ancil_atmos : nlookups ',nlookups
  CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
END IF

! Determine the workspace (len_ancil) required for the ancillaries.
CALL rcf_calc_len_ancil (Output_Grid % loc_p_field,           &
                         Output_Grid % loc_r_field,           &
                         Output_Grid % loc_p_rows,  len_ancil )


IF (PrintStatus >= PrStatus_Normal .AND. mype == 0 ) THEN
  WRITE(umMessage,'(A,I0)') ' rcf_ancil_atmos : len_ancil ',len_ancil
  CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
END IF

! ===============================================

ALLOCATE ( fixhd_ancil(LenFixHd,num_ancil_files) )
ALLOCATE ( inthd_ancil(LenInthd_anc,num_ancil_files) )
ALLOCATE ( realhd_ancil(LenRealhd_anc,num_ancil_files) )

fixhd_ancil(:,:) = 0
inthd_ancil(:,:) = 0
realhd_ancil (:,:) = 0.0

! ===============================================

ALLOCATE ( lookup_ancil(Hdr_Out % Len1LookUp, nlookups) )
ALLOCATE ( nlookup(num_ancil_requests) )
ALLOCATE ( lookup_step(num_ancil_requests) )
ALLOCATE ( levels(num_ancil_requests) )
ALLOCATE ( ancil_add(num_ancil_requests) )
ALLOCATE ( lookup_start(num_ancil_files) )

ancil_add(:) = 0

CALL inancila_rcf_inancila(LenFixHd,                                    &
     LenInthd_anc,                                         &
     LenRealhd_anc,                                        &
     Hdr_Out % Len1LevDepC,                                &
     Hdr_Out % Len2LevDepC,                                &
     fixhd_ancil,inthd_ancil,realhd_ancil,lookup_ancil,    &
     Hdr_Out % RealC,                                      &
     Hdr_Out % LevDepC,                                    &
     nlookups,                                             &
     lookup_start,                                         &
     Hdr_Out % Len1LookUp,                                 &
     Output_Grid % glob_p_row_length,                      &
     Output_Grid % loc_p_row_length,                       &
     Output_Grid % glob_p_rows,                            &
     Output_Grid % loc_p_rows,                             &
     Output_Grid % glob_u_rows,                            &
     Output_Grid % glob_r_row_length,                      &
     Output_Grid % glob_r_rows,                            &
     Output_Grid % loc_r_row_length,                       &
     Output_Grid % loc_r_rows,                             &
     Output_Grid % model_levels,                           &
     Output_Grid % tr_levels,                              &
     Output_Grid % st_levels,                              &
     Output_Grid % sm_levels,                              &
     Output_Grid % ozone_levels,                           &
     tpps_ozone_levels,                                    &
     ancil_add,                                            &
     ErrorStatus,cmessage)

IF (ErrorStatus /= 0) THEN
  WRITE(umMessage,*) ' Error in inancila_rcf_inancila'
  CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
  WRITE(umMessage,*) ' CMESSAGE ',cmessage
  CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
  WRITE(umMessage,*) ' ErrorStatus ',ErrorStatus
  CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
  CALL Ereport ( RoutineName, ErrorStatus, cmessage)
END IF

! ===============================================
! Locate Ice fraction in Output dump

CALL Rcf_Locate ( stashcode_prog_sec, stashcode_icefrac,             &
                  Fields_Out, Field_Count_Out, Pos_IceFrac_Out,      &
                  zero_ok_arg=.TRUE. )

IF ( Pos_IceFrac_Out /= 0 ) THEN  ! Ice Fraction in output dump

  IceFrac_Out => Fields_Out (Pos_IceFrac_Out)
  CALL Rcf_Alloc_Field ( IceFrac_Out )

  ! Locate Ice fraction in input dump

  CALL Rcf_Locate ( stashcode_prog_sec, stashcode_icefrac,           &
                    Fields_In, Field_Count_In, Pos_IceFrac_In,       &
                    zero_ok_arg=.TRUE. )

  IF ( Pos_IceFrac_In /= 0) THEN  ! Ice Fraction in input dump

    IceFrac_In => Fields_In (Pos_IceFrac_In)
    CALL Rcf_Alloc_Field ( IceFrac_In )

    ! Read in Ice Fraction

    CALL Rcf_Read_Field ( IceFrac_in, Hdr_In, Decomp_rcf_input )

    ! Set interpolation
    IF (h_int_active) THEN
      IceFrac_In % interp = interp_h_only
    ELSE
      IceFrac_In % interp = interp_copy
    END IF

    CALL rcf_interpolate( IceFrac_In, IceFrac_Out, Input_Grid,   &
                          Output_Grid, dummy, dummy )

    CALL Rcf_Dealloc_Field ( IceFrac_In )

  ELSE ! Ice Fraction not in input dump

    !  Ice_Fraction is not in the input dump. It is assumed here
    !  that an ITEMS namelist with SOURCE=2 or similar exists for
    !  the Ice Fraction Stash Code. The rcf will abort in Create_Dump
    !  if the Output Ice Fraction field cannot be initialised.

  END IF

ELSE ! Ice Fraction not in output dump

  ! replanca_rcf_replanca expects an Ice Fraction field of length p_field
  ! even if not required - set up and allocate.

  Dummy_IceFr % level_size = Output_Grid % loc_p_field
  Dummy_IceFr % levels     = 1
  Dummy_IceFr % stashmaster => Rcf_Exppx( 1, 0, stashcode_icefrac )

  IceFrac_Out => Dummy_IceFr
  CALL Rcf_Alloc_Field (IceFrac_out)

END IF

! ===============================================

! Locate TStar in output dump.

CALL Rcf_Locate ( stashcode_prog_sec, stashcode_tstar,               &
                  Fields_Out, Field_Count_Out, Pos_TStar_Out,        &
                  zero_ok_arg=.TRUE. )

IF ( Pos_TStar_Out /= 0 ) THEN  ! TStar in output dump

  TStar_Out => Fields_Out (Pos_TStar_Out)
  CALL Rcf_Alloc_Field ( TStar_Out )

  CALL Rcf_Locate ( stashcode_prog_sec, stashcode_tstar,             &
                    Fields_In, Field_Count_In, Pos_TStar_In,         &
                    zero_ok_arg=.TRUE. )

  IF ( Pos_TStar_In /= 0) THEN  ! TStar in input dump

    TStar_In => Fields_In (Pos_TStar_In)
    CALL Rcf_Alloc_Field ( TStar_In )

    CALL Rcf_Read_Field ( TStar_In, Hdr_In, Decomp_rcf_input )

    ! Set interpolation
    IF (h_int_active) THEN
      TStar_In % interp = interp_h_only
    ELSE
      TStar_In % interp = interp_copy
    END IF

    CALL rcf_interpolate( TStar_In, TStar_Out, Input_Grid,     &
                          Output_Grid, dummy, dummy )

    CALL Rcf_Dealloc_Field ( TStar_In )

  ELSE ! TStar not in input dump

    !  TStar is not in the input dump. It is assumed here
    !  that an ITEMS namelist with SOURCE=2 or similar exists for
    !  the TStar Stash Code. The rcf will abort in Create_Dump
    !  if the Output TStar field cannot be initialised.

  END IF

ELSE ! TStar not in output dump

  ! replanca_rcf_replanca expects a TStar field of length p_field
  ! even if not required - set up and allocate.

  Dummy_TStar % level_size = Output_Grid % loc_p_field
  Dummy_TStar % levels     = 1
  Dummy_TStar % stashmaster => Rcf_Exppx( 1, 0, stashcode_tstar )

  TStar_Out => Dummy_TStar
  CALL Rcf_Alloc_Field (TStar_out)

END IF

! ===============================================
! Locate Land_Frac in output dump.

CALL Rcf_Locate ( stashcode_prog_sec, stashcode_land_frac,           &
                  Fields_Out, Field_Count_Out, Pos_Land_Frac_Out,    &
                  zero_ok_arg=.TRUE. )

IF ( Pos_Land_Frac_Out /= 0 ) THEN  ! Land_Frac in output dump

  Land_Frac_Out => Fields_Out (Pos_Land_Frac_Out)
  CALL Rcf_Alloc_Field ( Land_Frac_Out )

  ! Check not reading land fraction from ancillary file:

  IF (.NOT. ANY( ancil_requests%stashcode == stashcode_land_frac)) THEN

    CALL Rcf_Locate ( stashcode_prog_sec, stashcode_land_frac,       &
                      Fields_In, Field_Count_In, Pos_Land_Frac_In,   &
                      zero_ok_arg=.TRUE. )

    IF ( Pos_Land_Frac_In /= 0) THEN  ! Land_Frac in input dump

      Land_Frac_In => Fields_In (Pos_Land_Frac_In)
      CALL Rcf_Alloc_Field ( Land_Frac_In )

      CALL Rcf_Read_Field ( Land_Frac_In, Hdr_In, Decomp_rcf_input )

      ! Dont allow horizonal interpolation for Land_Frac:
      IF (h_int_active) THEN
        ErrorStatus = 20
        WRITE (CMessage, '(2A)') 'Horizontal interpolation of ' &
        ,'land fraction not allowed: Need ancillary'
        CALL Ereport ( RoutineName, ErrorStatus, cmessage)
      ELSE
        Land_Frac_In % interp = interp_copy
      END IF

      CALL rcf_interpolate( Land_Frac_In, Land_Frac_Out, Input_Grid, &
                          Output_Grid, dummy, dummy )


      CALL Rcf_Dealloc_Field ( Land_Frac_In )

    ELSE ! Land frac not in input dump

      ErrorStatus = 22
      WRITE (CMessage, '(2A)') 'Land fraction is not in the input dump ' &
      ,'=> an ancillary must be provided.'
      CALL Ereport ( RoutineName, ErrorStatus, cmessage)

    END IF

  END IF ! Land_Frac read from ancillary

ELSE ! Land_Frac not in output dump

  WRITE(umMessage,*)'Land Frac is not in output dump => setting to dummy'
  CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')

  ! replanca_rcf_replanca expects a Land_Frac field of length p_field
  ! even if not required - set up and allocate.

  Dummy_Land_Frac % level_size = Output_Grid % loc_p_field
  Dummy_Land_Frac % levels     = 1
  Dummy_Land_Frac % stashmaster => Rcf_Exppx(1,0,stashcode_land_frac)

  Land_Frac_Out => Dummy_Land_Frac
  CALL Rcf_Alloc_Field (Land_Frac_out)

END IF

! ===============================================
! Locate Tstar_Land in output dump.

CALL Rcf_Locate ( stashcode_prog_sec, stashcode_tstar_land,          &
                  Fields_Out, Field_Count_Out, Pos_Tstar_Land_Out,   &
                  zero_ok_arg=.TRUE. )

IF ( Pos_Tstar_Land_Out /= 0 ) THEN  ! Tstar_Land in output dump

  Tstar_Land_Out => Fields_Out (Pos_Tstar_Land_Out)
  CALL Rcf_Alloc_Field ( Tstar_Land_Out )

  CALL Rcf_Locate ( stashcode_prog_sec, stashcode_tstar_land,        &
                    Fields_In, Field_Count_In, Pos_Tstar_Land_In,    &
                    zero_ok_arg=.TRUE. )

  IF ( Pos_Tstar_Land_In /= 0) THEN  ! Tstar_Land in input dump

    Tstar_Land_In => Fields_In (Pos_Tstar_Land_In)
    CALL Rcf_Alloc_Field ( Tstar_Land_In )

    CALL Rcf_Read_Field ( Tstar_Land_In, Hdr_In, Decomp_rcf_input )

    ! Set interpolation
    IF (h_int_active) THEN
      Tstar_Land_In % interp = interp_h_only
    ELSE
      Tstar_Land_In % interp = interp_copy
    END IF

    CALL rcf_interpolate( Tstar_Land_In, Tstar_Land_Out, Input_Grid, &
                          Output_Grid, dummy, dummy )

    CALL Rcf_Dealloc_Field ( Tstar_Land_In )

  ELSE ! Tstar_Land not in input dump

    !  Tstar_Land is not in the input dump. It is assumed here
    !  that an ITEMS namelist with SOURCE=2 or similar exists for
    !  the Tstar_Land Stash Code. The rcf will abort in Create_Dump
    !  if the Output Tstar_Land field cannot be initialised.

    WRITE(umMessage,*)'Land T* is not in input dump => setting to gbm T*'
    CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
    Tstar_Land_out % DATA (:,:) = Tstar_out % DATA(:,:)

  END IF

ELSE ! Tstar_Land not in output dump

  WRITE(umMessage,*)'Land T* is not in output dump => setting to dummy'
  CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')

  ! replanca_rcf_replanca expects a Tstar_Land field of length p_field
  ! even if not required - set up and allocate.

  Dummy_Tstar_Land % level_size = Output_Grid % loc_p_field
  Dummy_Tstar_Land % levels     = 1
  Dummy_Tstar_Land % stashmaster => Rcf_Exppx(1,0,stashcode_tstar_land)

  Tstar_Land_Out => Dummy_Tstar_Land
  CALL Rcf_Alloc_Field (Tstar_Land_out)

END IF

! ===============================================
! Locate Tstar_Sea in output dump.

CALL Rcf_Locate ( stashcode_prog_sec, stashcode_tstar_sea,           &
                  Fields_Out, Field_Count_Out, Pos_Tstar_Sea_Out,    &
                  zero_ok_arg=.TRUE. )

IF ( Pos_Tstar_Sea_Out /= 0 ) THEN  ! Tstar_Sea in output dump

  Tstar_Sea_Out => Fields_Out (Pos_Tstar_Sea_Out)
  CALL Rcf_Alloc_Field ( Tstar_Sea_Out )

  CALL Rcf_Locate ( stashcode_prog_sec, stashcode_tstar_sea,         &
                    Fields_In, Field_Count_In, Pos_Tstar_Sea_In,     &
                    zero_ok_arg=.TRUE. )

  IF ( Pos_Tstar_Sea_In /= 0) THEN  ! Tstar_Sea in input dump

    Tstar_Sea_In => Fields_In (Pos_Tstar_Sea_In)
    CALL Rcf_Alloc_Field ( Tstar_Sea_In )

    CALL Rcf_Read_Field ( Tstar_Sea_In, Hdr_In, Decomp_rcf_input )

    ! Set interpolation
    IF (h_int_active) THEN
      Tstar_Sea_In % interp = interp_h_only
    ELSE
      Tstar_Sea_In % interp = interp_copy
    END IF

    CALL rcf_interpolate( Tstar_Sea_In, Tstar_Sea_Out, Input_Grid, &
                          Output_Grid, dummy, dummy )

    CALL Rcf_Dealloc_Field ( Tstar_Sea_In )

  ELSE ! Tstar_Sea not in input dump

    !  Tstar_Sea is not in the input dump. It is assumed here
    !  that an ITEMS namelist with SOURCE=2 or similar exists for
    !  the Tstar_Sea Stash Code. The rcf will abort in Create_Dump
    !  if the Output Tstar_Sea field cannot be initialised.

    WRITE(umMessage,*)'Open sea T* is not in input dump => setting to gbm T*'
    CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
    Tstar_Sea_out % DATA (:,:) = Tstar_out % DATA(:,:)

  END IF

ELSE ! Tstar_Sea not in output dump

  ! replanca_rcf_replanca expects a Tstar_Sea field of length p_field
  ! even if not required - set up and allocate.

  Dummy_Tstar_Sea % level_size = Output_Grid % loc_p_field
  Dummy_Tstar_Sea % levels     = 1
  Dummy_Tstar_Sea % stashmaster => Rcf_Exppx(1,0,stashcode_tstar_sea)

  Tstar_Sea_Out => Dummy_Tstar_Sea
  CALL Rcf_Alloc_Field (Tstar_Sea_out)

END IF

! ===============================================
! Locate Tstar_Sice in output dump.

CALL Rcf_Locate ( stashcode_prog_sec, stashcode_tstar_sice,          &
                  Fields_Out, Field_Count_Out, Pos_Tstar_Sice_Out,   &
                  zero_ok_arg=.TRUE. )

IF ( Pos_Tstar_Sice_Out /= 0 ) THEN  ! Tstar_Sice in output dump

  Tstar_Sice_Out => Fields_Out (Pos_Tstar_Sice_Out)
  CALL Rcf_Alloc_Field ( Tstar_Sice_Out )

  CALL Rcf_Locate ( stashcode_prog_sec, stashcode_tstar_sice,        &
                    Fields_In, Field_Count_In, Pos_Tstar_Sice_In,    &
                    zero_ok_arg=.TRUE. )

  IF ( Pos_Tstar_Sice_In /= 0) THEN  ! Tstar_Sice in input dump

    Tstar_Sice_In => Fields_In (Pos_Tstar_Sice_In)
    CALL Rcf_Alloc_Field ( Tstar_Sice_In )

    CALL Rcf_Read_Field ( Tstar_Sice_In, Hdr_In, Decomp_rcf_input )

    ! Set interpolation
    IF (h_int_active) THEN
      Tstar_Sice_In % interp = interp_h_only
    ELSE
      Tstar_Sice_In % interp = interp_copy
    END IF

    CALL rcf_interpolate( Tstar_Sice_In, Tstar_Sice_Out, Input_Grid, &
                          Output_Grid, dummy, dummy )

    CALL Rcf_Dealloc_Field ( Tstar_Sice_In )

  ELSE ! Tstar_Sice not in input dump

    !  Tstar_Sice is not in the input dump. It is assumed here
    !  that an ITEMS namelist with SOURCE=2 or similar exists for
    !  the Tstar_Sice Stash Code. The rcf will abort in Create_Dump
    !  if the Output Tstar_Sice field cannot be initialised.

    WRITE(umMessage,*)'Sea-ice T* is not in input dump => setting to gbm T*'
    CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
    Tstar_Sice_out % DATA (:,:) = Tstar_out % DATA(:,:)

  END IF

ELSE ! Tstar_Sice not in output dump

  ! replanca_rcf_replanca expects a Tstar_Sice field of length p_field
  ! even if not required - set up and allocate.

  Dummy_Tstar_Sice % level_size = Output_Grid % loc_p_field
  Dummy_Tstar_Sice % levels     = 1
  Dummy_Tstar_Sice % stashmaster => Rcf_Exppx(1,0,stashcode_tstar_sice)

  Tstar_Sice_Out => Dummy_Tstar_Sice
  CALL Rcf_Alloc_Field (Tstar_Sice_out)

END IF
! ===============================================
! Locate TStar in output dump.

CALL Rcf_Locate ( stashcode_prog_sec, stashcode_tstar_anom,             &
                  Fields_Out, Field_Count_Out, Pos_TStar_Anom,          &
                  zero_ok_arg=.TRUE. )

IF ( Pos_TStar_Anom /= 0 ) THEN  ! TStar Anomaly in output dump


  IF ( PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
    WRITE(umMessage,*) 'Doing SST anomaly.'
    CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
    WRITE(umMessage,*) 'Pos_Tstar_Anom=',Pos_Tstar_Anom
    CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
  END IF
  Tstar_Anom => Fields_Out (Pos_Tstar_Anom)
  CALL Rcf_Alloc_Field ( Tstar_Anom )

ELSE ! TStar Anomaly not in output dump

  ! replanca_rcf_replanca expects a TStar Anomaly field of length p_field
  ! even if not required - set up and allocate.

  IF ( PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
    WRITE(umMessage,*)'TStar Anomaly is not in output dump => setting to dummy'
    CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
  END IF
  Dummy_TStar_Anom % level_size = Output_Grid % loc_p_field
  Dummy_TStar_Anom % levels     = 1
  Dummy_TStar_Anom % stashmaster =>                                     &
                             Rcf_Exppx( 1, 0, stashcode_tstar_anom )

  TStar_Anom => Dummy_TStar_Anom
  CALL Rcf_Alloc_Field (TStar_Anom)

END IF

! ===============================================
ALLOCATE ( ancil_data(len_ancil) )
ancil_data(:) = 0.0

! DEPENDS ON: replanca_rcf_replanca
CALL replanca_rcf_replanca (                              &
     Hdr_Out % FixHd (FH_VTYear),                         &
     Hdr_Out % FixHd (FH_VTMonth),                        &
     Hdr_Out % FixHd (FH_VTDay),                          &
     Hdr_Out % FixHd (FH_VTHour),                         &
     Hdr_Out % FixHd (FH_VTMinute),                       &
     Hdr_Out % FixHd (FH_VTSecond),                       &
     Output_Grid % loc_p_field,                           &
     Output_Grid % loc_p_rows,                            &
     Output_Grid % loc_u_field,                           &
     Output_Grid % loc_v_field,                           &
     Output_Grid % loc_r_field,                           &
     Output_Grid % loc_land_field,                        &
     ancil_data,                                          &
     local_lsm_out,                                       &
     IceFrac_Out % DATA,                                  &
     TStar_Out % DATA,                                    &
     Land_Frac_Out % DATA,                                &
     TStar_Land_Out % DATA,                               &
     TStar_Sea_Out % DATA,                                &
     TStar_Sice_Out % DATA,                               &
     TStar_Anom % DATA,                                   &
     Hdr_Out % RealC (RC_LatSpacing),                     &
     Hdr_Out % RealC (RC_FirstLat),                       &
     Hdr_Out % Len1LookUp, LenFixHd,                      &
     LenInthd_anc,                                        &
     len_ancil,fixhd_ancil,                               &
     inthd_ancil,lookup_ancil,lookup_ancil,               &
     lookup_start,nlookups,                               &
     ErrorStatus,cmessage)

IF (ErrorStatus /= 0) THEN
  WRITE(umMessage,*) ' Error in replanca_rcf_replanca'
  CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
  WRITE(umMessage,*) ' CMESSAGE ',cmessage
  CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
  WRITE(umMessage,*) ' ErrorStatus ',ErrorStatus
  CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
  CALL Ereport ( RoutineName, ErrorStatus, cmessage)
END IF

! ===============================================
! Write ancillary fields into output dump

DO i=1,num_ancil_requests
  IF (ancil_requests(i)%stashcode == stashcode_lsm) THEN
    ianc_Mask = i
  END IF
  IF (ancil_requests(i)%stashcode == stashcode_icefrac) THEN
    ianc_IceFrac = i
  END IF
  IF (ancil_requests(i)%stashcode == stashcode_tstar) THEN
    ianc_TStar = i
  END IF
END DO


! Setup output decomposition for writing.
orig_decomp = current_decomp_type
IF (orig_decomp /= decomp_rcf_output) THEN
  CALL Change_Decomposition( decomp_rcf_output )
END IF

DO k=1,num_ancil_requests

  IF (ancil_requests(k)%stashcode == stashcode_lsm  .OR.  &
      ancil_requests(k)%stashcode == stashcode_icefrac  .OR.  &
      ancil_requests(k)%stashcode == stashcode_tstar) THEN

    !     The ancillary fields : Land-Sea mask, Ice Fraction and TStar
    !     are not written to the output dump in this loop.

    CYCLE

  ELSE ! write this ancillary to the output dump
    !

    ! Locate position of ancillary field in output dump
    IF (ancil_add(k) /= 0) THEN
      CALL Rcf_Locate ( ancil_requests(k)%section,             &
                        ancil_requests(k)%item,                &
                        Fields_Out, Field_Count_Out, Pos)
    END IF

    ! For Land fields, compress field to land points first

    IF (Fields_out(Pos) % stashmaster % grid_type ==                &
                                        ppx_atm_compressed) THEN
      ! Allocate array used in compress to mask.
      ALLOCATE(ancil_land_temp(Fields_Out(pos) % level_size))

      DO j=1, Fields_out(pos) % levels

        i_full = ancil_add(k)+(j-1) * Output_Grid % loc_p_field
        i_land = ancil_add(k)+(j-1) * Fields_Out(pos) % level_size

        ! Use temporary array since this sets INTENT(OUT)
        CALL compress_to_mask (                &
             ancil_data( i_full ),             &
             ancil_land_temp,                  &
             local_lsm_out,                    &
             Output_Grid % loc_p_field,        &
             local_land_out )

        DO i = 1, Fields_Out(pos) % level_size
          ancil_data(i_land + i - 1) = ancil_land_temp(i)
        END DO

      END DO  !  Loop over j
      DEALLOCATE(ancil_land_temp)

    END IF  !  If land-only field

    ! ===============================================

    ! Write the ancillary fields out to the dump
    multi_pe_in=.TRUE.
    ! DEPENDS ON: rcf_writflds
    CALL Rcf_writflds (Hdr_Out % UnitNum,                        &
                       Fields_out(pos) % levels,                 &
                       Fields_Out(pos) % dump_pos,               &
                       Hdr_Out % Lookup,                         &
                       Hdr_Out % Len1LookUp,                     &
                       ancil_data(ancil_add(k)),                 &
                       Fields_out(pos) % level_size,             &
                       Hdr_Out % FixHd,                          &
                       ErrorStatus,cmessage,multi_pe_in)

    IF (ErrorStatus /= 0) THEN
      WRITE(umMessage,*) ' Error in WRITFLDS for Ancillary Field ',k
      CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
      WRITE(umMessage,*) ' CMESSAGE ',cmessage
      CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
      WRITE(umMessage,*) ' ErrorStatus ',ErrorStatus
      CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
      CALL Ereport ( RoutineName, ErrorStatus, cmessage)
    END IF

  END IF  ! if K etc

END DO

! Change back to original if different.
IF (current_decomp_type /= orig_decomp) THEN
  CALL Change_Decomposition( orig_decomp )
END IF

DEALLOCATE ( lookup_ancil )
DEALLOCATE ( nlookup )
DEALLOCATE ( lookup_step )
DEALLOCATE ( levels )
DEALLOCATE ( ancil_add )

DEALLOCATE ( fixhd_ancil )
DEALLOCATE ( inthd_ancil )
DEALLOCATE ( realhd_ancil )

DEALLOCATE ( ancil_data )
DEALLOCATE ( lookup_start )

! ===============================================
! We want to mark fields which have been processed here if its
! possible they will be processed again later - the check used to be performed
! in rcf_create_dump.  For some reason tstar and icefrac are being
! treated as already processed when labelled as set_to_mdi even though we want
! to set them to the values calculated here. Lets just keep the same logic here
! as was found in rcf_create_dump.

!Write TStar to the output dump

IF ( Pos_TStar_Out /= 0 ) THEN

  IF (mype == 0) THEN
    IF (PrintStatus >= PrStatus_Diag ) THEN
      WRITE(umMessage,*) 'Writing TStar to output dump.'
      CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
    END IF
  END IF

  CALL Rcf_Write_Field ( TStar_Out, Hdr_Out, Decomp_rcf_output )
  IF (data_source(pos_tstar_out) % source == input_dump .OR. &
      data_source(pos_tstar_out) % source == set_to_mdi) THEN
    data_source(pos_tstar_out) % source = already_processed
  END IF

END IF

CALL Rcf_Dealloc_Field ( TStar_Out )

! ===============================================

CALL Rcf_Dealloc_Field ( Land_Frac_Out )

! ===============================================

! Write TStar_land to the output dump

IF ( Pos_TStar_Land_Out /= 0 ) THEN

  IF (mype == 0) THEN
    IF (PrintStatus >= PrStatus_Diag ) THEN
      WRITE(umMessage,*) 'Writing TStar_Land to output dump.'
      CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
    END IF
  END IF

  CALL Rcf_Write_Field ( TStar_Land_Out, Hdr_Out, Decomp_rcf_output )
  IF (data_source(pos_tstar_land_out) % source == input_dump) THEN
    data_source(pos_tstar_land_out) % source = already_processed
  END IF

END IF

CALL Rcf_Dealloc_Field ( TStar_Land_Out )

! ===============================================

! Write TStar_sea to the output dump

IF ( Pos_TStar_Sea_Out /= 0 ) THEN

  IF (mype == 0) THEN
    IF (PrintStatus >= PrStatus_Diag ) THEN
      WRITE(umMessage,*) 'Writing TStar_Sea to output dump.'
      CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
    END IF
  END IF

  CALL Rcf_Write_Field ( TStar_Sea_Out, Hdr_Out, Decomp_rcf_output )
  IF (data_source(pos_tstar_sea_out) % source == input_dump) THEN
    data_source(pos_tstar_sea_out) % source = already_processed
  END IF

END IF

CALL Rcf_Dealloc_Field ( TStar_Sea_Out )

! ===============================================

! Write TStar_sice to the output dump

IF ( Pos_TStar_Sice_Out /= 0 ) THEN

  IF (mype == 0) THEN
    IF (PrintStatus >= PrStatus_Diag ) THEN
      WRITE(umMessage,*) 'Writing TStar_Sice to output dump.'
      CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
    END IF
  END IF

  CALL Rcf_Write_Field ( TStar_Sice_Out, Hdr_Out, Decomp_rcf_output )
  IF (data_source(pos_tstar_sice_out) % source == input_dump) THEN
    data_source(pos_tstar_sice_out) % source = already_processed
  END IF

END IF

CALL Rcf_Dealloc_Field ( TStar_Sice_Out )

! ===============================================
! ===============================================

! Write Ice Fraction to the dump

IF ( Pos_IceFrac_Out /= 0 ) THEN

  IF (mype == 0) THEN
    IF (PrintStatus >= PrStatus_Diag ) THEN
      WRITE(umMessage,*) 'Writing Ice Fraction to output dump.'
      CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
    END IF
  END IF

  CALL Rcf_Write_Field ( IceFrac_Out, Hdr_Out, Decomp_rcf_output )
  IF (data_source(pos_icefrac_out) % source == input_dump .OR. &
      data_source(pos_icefrac_out) % source == set_to_mdi) THEN
    data_source(pos_icefrac_out) % source = already_processed
  END IF

END IF

CALL Rcf_Dealloc_Field ( IceFrac_Out )

! ===============================================

! Write T* Anomaly to the dump

IF ( Pos_TStar_Anom /= 0 ) THEN

  IF (mype == 0) THEN
    IF (PrintStatus >= PrStatus_Diag ) THEN
      WRITE(umMessage,*) 'Writing TStar Anomaly to output dump.'
      CALL umPrint(umMessage,src='rcf_ancil_atmos_mod')
    END IF
  END IF

  CALL Rcf_Write_Field ( TStar_Anom, Hdr_Out, Decomp_rcf_output )
  IF (data_source(pos_tstar_anom) % source == input_dump) THEN
    data_source(pos_tstar_anom) % source = already_processed
  END IF

END IF

CALL Rcf_Dealloc_Field ( TStar_Anom )

! ===============================================


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Ancil_Atmos
END MODULE Rcf_Ancil_Atmos_Mod
