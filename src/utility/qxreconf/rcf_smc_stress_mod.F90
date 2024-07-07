! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE Rcf_smc_stress_Mod

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE 

PRIVATE :: lhook, dr_hook, jpim, jprb

! Description:
!   This module calculate soil moisture stress fields,
!   and then interpolate smc_stress into the new domain.
!   This smc_stress will be converted into the soil moisture
!   in the post process routine.
!
! Method:
!   Read input of smc, smc at wilting and critical point.
!   Calculate smc stress using these inputs and then interpolate
!   smc_stress to the new domain. The smc in the new domain will be
!   then calculated from smc_stress in the post process routine.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

!-----------------------------------------------------------------------

! Declare a dummy (unused) field_type object that we point the soil 
! properties to if these are missing or aren't read in. This is never 
! allocated, but allows us to use rcf_dealloc_field on all the soil 
! properties at the end of the main subroutine. We declare this at the 
! module level so we don't need to pass it through the subroutine calls.

TYPE( field_type ), TARGET, SAVE :: dummy_field 

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SMC_STRESS_MOD'

CONTAINS

SUBROUTINE Rcf_smc_stress(smc_field, fields_in, field_count_in, hdr_in,     &
                          fields_out, field_count_out, hdr_out, data_source)
 
! Subroutine Rcf_smc_stress - calculate the soil moisture stress
 
USE umPrintMgr, ONLY:       &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal,        &
    newline

USE UM_ParCore, ONLY:       &
    nproc,                  &
    mype

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Grid_Type_Mod, ONLY: &
    Input_Grid

USE decomp_params, ONLY: &
    decomp_rcf_input,    &
    decomp_rcf_output

USE rcf_nlist_recon_science_mod, ONLY: &
    use_smc_stress

USE Rcf_Interp_Weights_Mod, ONLY: &
    h_int_active

USE Rcf_V_Int_Ctl_Mod, ONLY: &
    v_int_active_soil

USE Rcf_Data_Source_Mod, ONLY: &
    data_source_type

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Alloc_Field_Mod, ONLY:&
    Rcf_DeAlloc_Field

USE um_stashcode_mod, ONLY: &
    stashcode_vol_smc_wilt, &
    stashcode_vol_smc_cri,  &
    stashcode_vol_smc_sat,  &
    stashcode_prog_sec

USE water_constants_mod, ONLY: rho_water

USE missing_data_mod, ONLY: rmdi

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Arguments
TYPE( field_type), INTENT(INOUT)     :: smc_field

TYPE( field_type ), POINTER          :: fields_in(:)
TYPE( field_type ), POINTER          :: fields_out(:)
TYPE (data_source_type), POINTER     :: data_source(:)

TYPE( field_type ), POINTER          :: smc_wlt
TYPE( field_type ), POINTER          :: smc_cri
TYPE( field_type ), POINTER          :: smc_sat
TYPE( field_type ), POINTER          :: smc_wlt_anc
TYPE( field_type ), POINTER          :: smc_cri_anc
TYPE( field_type ), POINTER          :: smc_sat_anc

TYPE( um_header_type ), INTENT(IN)   :: hdr_in
TYPE( um_header_type ), INTENT(IN)   :: hdr_out
INTEGER, INTENT(IN)                  :: field_count_in
INTEGER, INTENT(IN)                  :: field_count_out

! Local variables
INTEGER                              :: i            ! Loop variables
INTEGER                              :: j
INTEGER                              :: k
INTEGER                              :: decomp_old   ! old decomposition
INTEGER                              :: stat         ! gcom status

INTEGER                              :: smc_wlt_srce ! source of soil
INTEGER                              :: smc_cri_srce ! properties for 
INTEGER                              :: smc_sat_srce ! target files
INTEGER                              :: dummy_srce   ! and dummy value 
                                                     ! for input files
INTEGER                              :: int_use_smc_stress
INTEGER                              :: info

CHARACTER (LEN=errormessagelength)  :: cmessage
INTEGER                             :: errorstatus
CHARACTER (LEN=*), PARAMETER        :: RoutineName = 'RCF_SMC_STRESS'
LOGICAL                             :: l_fields_missing

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  WRITE(umMessage,'(A)') 'Processing soil moisture (stashcode 9) '
  CALL umPrint(umMessage,src=routinename)
END IF

!------------------------------------------------------------------
! Find and read input/output soil properties
! If any of these are missing these routines will reset 
! use_smc_stress to .false.
!------------------------------------------------------------------

! initialise flag for missing fields

l_fields_missing=.FALSE.
! input wilting point (if present): smc_wlt
CALL read_field_if_present(stashcode_prog_sec, stashcode_vol_smc_wilt,     &
     fields_in, field_count_in, Hdr_In, decomp_rcf_input, smc_wlt,         &
     data_source, dummy_srce, l_fields_missing)
! input critical point (if present): smc_cri
CALL read_field_if_present(stashcode_prog_sec, stashcode_vol_smc_cri,      &
     fields_in, field_count_in, Hdr_In, decomp_rcf_input, smc_cri,         &
     data_source, dummy_srce, l_fields_missing)
! input saturation point (if present): smc_sat
CALL read_field_if_present(stashcode_prog_sec, stashcode_vol_smc_sat,      &
     fields_in, field_count_in, Hdr_In, decomp_rcf_input, smc_sat,         &
     data_source, dummy_srce, l_fields_missing)
! output wilting point (if present and from ancil): smc_wlt_anc
CALL read_field_if_present(stashcode_prog_sec, stashcode_vol_smc_wilt,     &
     fields_out, field_count_out, Hdr_Out, decomp_rcf_output, smc_wlt_anc, &
     data_source, smc_wlt_srce, l_fields_missing, l_out_fld=.TRUE.)
! output critical point (if present and from ancil): smc_cri_anc
CALL read_field_if_present(stashcode_prog_sec, stashcode_vol_smc_cri,      &
     fields_out, field_count_out, Hdr_Out, decomp_rcf_output, smc_cri_anc, &
     data_source, smc_cri_srce, l_fields_missing, l_out_fld=.TRUE.)
! output saturation point (if present and from ancil): smc_sat_anc
CALL read_field_if_present(stashcode_prog_sec, stashcode_vol_smc_sat,      &
     fields_out, field_count_out, Hdr_Out, decomp_rcf_output, smc_sat_anc, &
     data_source, smc_sat_srce, l_fields_missing, l_out_fld=.TRUE.)

IF (l_fields_missing) THEN

  ! Print a warning if we have switched use_smc_stress to false
  IF ( use_smc_stress .AND. PrintStatus >= PrStatus_Normal ) THEN
    WRITE(cmessage,'(A)') &
       'Input and output soil properties are not all present. '// newline//&
       ' Not converting soil moisture to soil stress and '     // newline//&
       ' resetting namelist variable use_smc_stress to .false.'
    errorstatus = -10
    CALL ereport(routinename, errorstatus, cmessage)
  END IF

  use_smc_stress=.FALSE.
END IF

!------------------------------------------------------------------
! Before we perform main calculations, check whether we really need
! to do this. We will skip these calculations and reset 
! use_smc_stress to .false. if all of the following are true:
!  (i) Horizontal interpolation is off: h_int_active=.false.
!  AND:
!  (ii) Vertical soil levels unchanged: v_int_active_soil=.false.
!  AND
!  (iii) Input and output soil properties are identical.
!------------------------------------------------------------------

!------------------------------------------------------------------
! Only perform these checks at all if use_smc_stress is still true 
! i.e. if all required soil properties are in input/output files
!------------------------------------------------------------------
IF (use_smc_stress) THEN

!  Check conditions (i) and (ii)
  IF ( (.NOT. h_int_active) .AND. (.NOT. v_int_active_soil) ) THEN

    !------------------------------------------------------------------
    ! Check condition (iii): Because we check each soil property in turn, 
    ! we start by setting int_use_smc_stress=0. If any points in input 
    ! file and ancillary differ, or the field source is neither the input 
    ! file nor an ancillary, then check_input_output_differ will set
    ! int_use_smc_stress to 1.
    ! If int_use_smc_stress remains 0, we reset use_smc_stress to .FALSE.
    ! and warn if necessary.
    !------------------------------------------------------------------

    int_use_smc_stress=0
    CALL check_input_output_differ(smc_wlt, smc_wlt_anc,   &
         smc_wlt_srce, int_use_smc_stress)
    CALL check_input_output_differ(smc_cri, smc_cri_anc,   &
         smc_cri_srce, int_use_smc_stress)
    CALL check_input_output_differ(smc_sat, smc_sat_anc,   &
         smc_sat_srce, int_use_smc_stress)

    ! Check for differences on any processor, and ensure all processors
    ! agree
    CALL gc_imax(1, nproc, info, int_use_smc_stress)

    IF (int_use_smc_stress == 0) THEN
      ! Then the input field and ancillary field are identical

      ! Indicate that the main soil calculation is no longer needed
      IF ( use_smc_stress .AND. PrintStatus >= PrStatus_Normal ) THEN
        WRITE(cmessage,'(A)') &
          'Input and output soil properties are identical. '    // newline//&
          ' Not converting soil moisture to soil stress and '   // newline//&
          ' resetting namelist variable use_smc_stress to .false.'
        errorstatus = -10
        CALL ereport(routinename, errorstatus, cmessage)
      END IF

      use_smc_stress = .FALSE.
    END IF

  END IF ! (.NOT. h_int_active) .AND. (.NOT. v_int_active_soil)

END IF ! use_smc_stres

!------------------------------------------------------------------
! Main calculations: only performed if use_smc_stress is still true
!------------------------------------------------------------------

IF (use_smc_stress) THEN

! Convert smc to VOL SMC:
  DO j = 1, Input_Grid % sm_levels
    smc_field % DATA(:,j) = smc_field % DATA(:,j)/ &
         ( Input_Grid % soil_depths(j)*rho_water )
  END DO

! Calculate smc_stress and use the new data for smc:
  DO j = 1, smc_field % levels

    DO k=1, smc_field % level_size

      IF (smc_cri%DATA(k,1)-smc_wlt%DATA(k,1) == 0 .OR.  &
          smc_cri%DATA(k,1)== rmdi .OR.                  &
          smc_wlt%DATA(k,1)== rmdi .OR.                  &
          smc_field % DATA(k,j) == rmdi) THEN
        smc_field % DATA(k,j) = rmdi
      ELSE
        smc_field % DATA(k,j)= &
             (smc_field % DATA(k,j)-smc_wlt % DATA(k,1))/&
             (smc_cri % DATA(k,1)-smc_wlt % DATA(k,1))
      END IF
    END DO
  END DO
  IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
    WRITE(umMessage,'(A)') &
         'Input Soil moisture has been converted to soil stress.'
    CALL umPrint(umMessage,src=routinename)
  END IF
END IF

! Deallocate the input/output fields
CALL Rcf_DeAlloc_Field( smc_sat_anc )
CALL Rcf_DeAlloc_Field( smc_cri_anc )
CALL Rcf_DeAlloc_Field( smc_wlt_anc )
CALL Rcf_DeAlloc_Field( smc_sat )
CALL Rcf_DeAlloc_Field( smc_cri )
CALL Rcf_DeAlloc_Field( smc_wlt )
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_smc_stress

SUBROUTINE read_field_if_present(stashcode_sec, stashcode_item,&
     fields, field_count, hdr, decomp, field, data_source,     &
     source, l_fields_missing, l_out_fld)

! Subroutine read_field_if_present
!
! Description:
!   This subroutine is used to check for the presence of the soil 
!   properties in the input file or the output file and open/read
!   the fields if present. 
!
! Method:
!   Note that if the fields are missing, we flag this with a logical
!

USE Rcf_UMhead_Mod, ONLY:      &
    um_header_type

USE Rcf_Locate_Mod, ONLY:      &
    Rcf_Locate

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Data_Source_Mod, ONLY: &
    data_source_type

USE items_nml_mod, ONLY: &
    Ancillary_File

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)                :: stashcode_sec
INTEGER, INTENT(IN)                :: stashcode_item
INTEGER, INTENT(IN)                :: field_count
INTEGER, INTENT(IN)                :: decomp
INTEGER, INTENT(OUT)               :: source

TYPE( field_type ), POINTER        :: fields(:)
TYPE( field_type ), POINTER        :: field
TYPE( um_header_type ), INTENT(IN) :: hdr
TYPE (data_source_type), POINTER   :: data_source(:)

LOGICAL, OPTIONAL, INTENT(IN)      :: l_out_fld
LOGICAL, INTENT(INOUT)             :: l_fields_missing
                                      ! Can be passed in as true
                                      ! if previous calls flagged
                                      ! missing fields
! Local variables
INTEGER                            :: pos
LOGICAL                            :: l_output_field
LOGICAL                            :: l_read_field

CHARACTER (LEN=*), PARAMETER :: RoutineName='READ_FIELD_IF_PRESENT'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set variable l_output_field
IF ( PRESENT( l_out_fld ) ) THEN
  l_output_field=l_out_fld
ELSE
  l_output_field=.FALSE.
END IF

! First attempt to locate the required field
CALL Rcf_Locate(stashcode_sec, stashcode_item,        &
                fields, field_count, pos, zero_ok_arg=.TRUE.)

! If present, allocate and read the field
IF ( pos /= 0 ) THEN

  IF (l_output_field) THEN
!   For output fields, we only need to read these if ancillaries
    source =  Data_Source( pos ) % Source
    IF ( source == Ancillary_File ) THEN
      l_read_field=.TRUE.
    ELSE
      l_read_field=.FALSE.
    END IF
  ELSE
!   Return a dummy value of source for fields in the 
!   input dump as this is not needed in rcf_smc_stress
    source = -99
!   Always read the input fields (if present)
    l_read_field=.TRUE.
  END IF

! Allocate and read fields
  IF (l_read_field) THEN
    field => fields(pos)
    CALL Rcf_Alloc_Field( field )
    CALL Rcf_Read_Field(  field, hdr, decomp )
  ELSE
    field => dummy_field
  END IF
ELSE
! If not present, reset l_fields_missing to .true.
  l_fields_missing=.TRUE.
  field => dummy_field
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE read_field_if_present

SUBROUTINE check_input_output_differ(soil_prop_in, soil_prop_anc,   &
                                     soil_prop_srce, int_use_smc_stress)

! Subroutine check_input_output_differ
!
! Description:
!   This subroutine is used to check whether the soil properties
!   differ between the input and output dumps. If they do, 
!   use_smc_stress is reset to .true.
!
! Method:
!   Assume that these are not identical unless all soil properties 
!   are either copied from the input file or using ancillary data 
!   identical to that in the input file.
!   Note this routine is only called when the input and output grid
!   are the same.
!

USE items_nml_mod, ONLY: &
    Input_Dump, Ancillary_File

IMPLICIT NONE

! Arguments

TYPE( field_type ), POINTER        :: soil_prop_in   ! Input file soil property
TYPE( field_type ), POINTER        :: soil_prop_anc  !  Ancillary soil property
INTEGER, INTENT(IN)                :: soil_prop_srce 
INTEGER, INTENT(INOUT)             :: int_use_smc_stress

CHARACTER (LEN=*), PARAMETER :: RoutineName='CHECK_INPUT_OUTPUT_DIFFER'

! Local variables
INTEGER                            :: n_soil_points  ! Loop variable
REAL, PARAMETER                    :: tol=1.0e-11

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( soil_prop_srce == Ancillary_File ) THEN
  
  n_soil_points = soil_prop_in % level_size

  IF ( ANY ( ABS(soil_prop_in % DATA(1:n_soil_points,1) -         &
                 soil_prop_anc % DATA(1:n_soil_points,1)) > tol ) ) THEN
    ! Ancillary data differs from input data.
    int_use_smc_stress = 1
  END IF
      
ELSE IF ( soil_prop_srce /= Input_Dump ) THEN
  ! Source for soil prop not input file or ancil.
  int_use_smc_stress = 1
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_input_output_differ

END MODULE Rcf_smc_stress_Mod
