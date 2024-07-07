! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Reads the ITEMS namelists

MODULE rcf_readnl_items_mod

IMPLICIT NONE

!  Subroutine Rcf_Readnl_Items - Reads the ITEMS namelists
!
! Description:
!   Reads the ITEMS namelists to control the source of required fields.
!
! Method:
!   Read initially into fixed size array and transferred to
!   to correctly sized dynamically allocated array later.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_READNL_ITEMS_MOD'

CONTAINS

SUBROUTINE rcf_readnl_items (unit_number)

USE ancil_mod,                     ONLY: &
    max_items,                           &
    stash_num_max

USE errormessagelength_mod,        ONLY: &
    errormessagelength

USE filenamelength_mod,            ONLY: &
    filenamelength

USE items_nml_mod,                 ONLY: &
    netcdf_varname_len,                  &
    read_nml_items,                      &
    ancillary_file,                      &
    input_dump,                          &
    set_to_const,                        &
    netcdf_file

USE populate_ancil_requests_mod, ONLY: &
    populate_ancil_requests

USE rcf_interp_weights_mod,        ONLY: &
    h_int_method,                        &
    nearest_neighbour

USE rcf_items_mod,                 ONLY: &
    Source_Array,                        &
    Sctn_Array,                          & 
    Item_Array,                          &
    Area_Array,                          &
    Upas_Array,                          &
    Upaa_Array,                          &
    Uprc_Array,                          &
    Upaf_Array,                          &
    Upnv_Array,                          &
    Num_Items

USE rcf_lsm_mod,                   ONLY: &
    lsm_fixed_value,                     &
    lsm_source

USE rcf_nlist_recon_technical_mod, ONLY: &
    select_output_fields,                &
    interp_all_fields

USE um_stashcode_mod,              ONLY: &
    stashcode_lsm,                       &
    stashcode_prog_sec

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER      ::  unit_number

! Local variables
INTEGER      ::  i ! Loop counter
! temporary arrays
INTEGER                        :: stashcode_temp ( max_items )
INTEGER                        :: section_temp ( max_items )
INTEGER                        :: item_temp ( max_items )
INTEGER                        :: source_temp ( max_items )
INTEGER                        :: domain_temp ( max_items )
INTEGER                        :: user_prog_section_temp ( max_items )
INTEGER                        :: user_prog_item_temp ( max_items )
INTEGER                        :: period_temp ( max_items )
INTEGER                        :: interval_temp ( max_items )
REAL                           :: real_constant_temp ( max_items )
CHARACTER (LEN=filenamelength) :: items_filename_temp ( max_items )
LOGICAL                        :: update_anc_temp ( max_items )
LOGICAL                        :: anc_configure_temp ( max_items ) = .FALSE.
CHARACTER (LEN=netcdf_varname_len) :: netcdf_varname_temp ( max_items )

! Arguments
CHARACTER(LEN=errormessagelength) :: reason_for_items_warning

CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_READNL_ITEMS'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

num_items = 0 
CALL read_nml_items(unit_number,         num_items,              &
                    stashcode_temp,      section_temp,           &
                    item_temp,           source_temp,            &    
                    domain_temp,         user_prog_section_temp, &   
                    user_prog_item_temp, period_temp,            &
                    interval_temp,       real_constant_temp,     &
                    update_anc_temp,     items_filename_temp,    &
                    netcdf_varname_temp)

IF ( num_items > 0 ) THEN
  ! Items namelists should be discarded in the following two cases
  IF (select_output_fields == interp_all_fields) THEN
    WRITE(reason_for_items_warning, '(A,I0,A)') &
         "select_output_fields = ", select_output_fields,  &
         " (interp all fields)"
    CALL rcf_print_items_discarded_warning()
    ! Reset num_items to 0 as all items namelist are discarded
    num_items = 0
  END IF
  IF ( h_int_method == nearest_neighbour) THEN
    WRITE(reason_for_items_warning, '(A,I0,A)') &
         "h_int_method = ", h_int_method, " (nearest neighbour)"
    CALL rcf_print_items_discarded_warning()
    ! Reset num_items to 0 as all items namelist are discarded
    num_items = 0
  END IF
END IF

DO i = 1, num_items
  ! Check where the land-sea mask is being read from.
  ! This is used later to set the size of land_field: whether read from
  ! ancillary, copied from the input dump (if no horizontal interpolation),
  ! set to a fixed value (e.g. 0 for aquaplanets) or something else.
  IF ( (section_temp (i) == stashcode_prog_sec) .AND.          &
       (item_temp (i) == stashcode_lsm) ) THEN
    lsm_source = source_temp(i)
    IF (source_temp(i) == set_to_const) THEN
      lsm_fixed_value = real_constant_temp(i)
    END IF
  END IF
  IF (source_temp(i) == ancillary_file) anc_configure_temp(i) = .TRUE.
END DO

!-----------------------------------------------------------------
! Allocate proper space for items list and copy data
!-----------------------------------------------------------------
IF ( Num_Items > 0 ) THEN
  ! Take temporary arrays read from items namelists and pass into source
  ! arrays and ancil request types
  ALLOCATE( sctn_array  ( num_items ) )
  ALLOCATE( item_array  ( num_items ) )
  ALLOCATE( source_array( num_items ) )
  ALLOCATE( area_array  ( num_items ) )
  ALLOCATE( uprc_array  ( num_items ) )
  ALLOCATE( upaf_array  ( num_items ) )
  ALLOCATE( upnv_array  ( num_items ) )
  ALLOCATE( upaa_array  ( num_items ) )
  ALLOCATE( upas_array  ( num_items ) )

  sctn_array  ( 1 : num_items ) = section_temp( 1: num_items )
  item_array  ( 1 : num_items ) = item_temp( 1: num_items )
  source_array( 1 : num_items ) = source_temp( 1: num_items )
  area_array  ( 1 : num_items ) = domain_temp( 1: num_items )
  uprc_array  ( 1 : num_items ) = real_constant_temp( 1: num_items )
  upaf_array  ( 1 : num_items ) = items_filename_temp( 1: num_items )
  upnv_array  ( 1 : num_items ) = netcdf_varname_temp( 1: num_items )
  upaa_array  ( 1 : num_items ) = user_prog_item_temp( 1: num_items )
  upas_array  ( 1 : num_items ) = user_prog_section_temp( 1: num_items )

  CALL populate_ancil_requests( stashcode_temp ( 1:num_items ),       &
       section_temp ( 1:num_items ),         &
       item_temp ( 1:num_items ),            &
       period_temp ( 1:num_items ),          &
       interval_temp ( 1:num_items ),        &
       items_filename_temp ( 1:num_items ),  &
       anc_configure_temp( 1:num_items) )
END IF !  if num_items > 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

CONTAINS

SUBROUTINE rcf_print_items_discarded_warning()
! Internal subroutine to avoid code duplication, 
! prints warning message if the items namelists are discarded
USE ereport_mod, ONLY: &
    ereport
USE umPrintMgr,  ONLY: &
    newline
USE UM_ParCore,  ONLY: &
    mype

IMPLICIT NONE
! Local variables
CHARACTER(LEN=errormessagelength+512) :: cmessage_long
INTEGER :: errorstatus

WRITE(cmessage_long, '(A)')                                                 &
 "All settings specified in ITEMS namelists will be discarded due "         &
 // newline // "to the following setting:"                                  &
 // newline //                                                              &
 reason_for_items_warning                                                   &
 // newline //                                                              &
 "Default source of all fields will be set to INPUT DUMP." //newline//      &
 "Note that depending on science setup and which fields are available"      &
 // newline //                                                              &
 "in the dump that the source of the field may be altered by the logic"     &
 // newline //                                                              &
 "in rcf_reset_data_source. A table of all fields and their corresponding"  &
 // newline //                                                              &
 "source values can be found in the main standard output when running"      &
 // newline //                                                              &
 "with PrintStatus >= PrStatus_Oper"                                  
errorstatus = -10
CALL ereport( routinename, errorstatus, cmessage_long )

END SUBROUTINE rcf_print_items_discarded_warning

END SUBROUTINE  Rcf_readnl_items

END MODULE rcf_readnl_items_mod
