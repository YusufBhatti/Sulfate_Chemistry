! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Ancillary processing for the Atmosphere model

MODULE rcf_netcdf_atmos_mod

!  Subroutine rcf_netcdf_atmos  - NetCDF processing for Atmosphere
!
! Description:
!    Controls all NetCDF processing for the atmosphere model
!
! Method:
!    1. Open NetCDF and obtain global attributes
!    2. Loop over fields and match with stashcode attribute
!         or netcdf_varname from namelist:item
!    3. Obtain pos using rcf_locate
!    4. Time interpolation
!    5. Write field to dump
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!    This file belongs in section: Reconfiguration_NetCDF
!
! Code Description:
!    Language: FORTRAN 90
!    This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_NETCDF_ATMOS_MOD'

CONTAINS

SUBROUTINE rcf_netcdf_atmos ( hdr_out, fields_out, field_count_out )

! USE statements here

USE ereport_mod,              ONLY: ereport
USE errormessagelength_mod,   ONLY: errormessagelength
USE decomp_params,            ONLY: decomp_rcf_output

USE glomap_clim_get_netcdffile_rec_mod, ONLY: &
    glomap_clim_get_netcdffile_rec

USE glomap_clim_netcdf_io_mod,          ONLY: &
    netcdf_fclose,                            &
    netcdf_fopen,                             &
    netcdf_get_data,                          &
    netcdf_get_file_info,                     &
    netcdf_get_var_att,                       &
    netcdf_get_var_info,                      &
    netcdf_var_check_dims

USE glomap_clim_netcdf_parameter_mod, ONLY: &
    ncdf_file_unit_no,                      &
    number_ncdf_files,                      &
    max_dims_per_var,                       &
    max_ncdf_files,                         &
    stash_name_len,                         &
    ncdf_list

USE items_nml_mod,            ONLY: &
    netcdf_varname_len,             &
    netcdf_unset

USE missing_data_mod,         ONLY: imdi

USE parkind1,                 ONLY: &
    jpim,                           &
    jprb

USE rcf_field_type_mod,       ONLY: field_type
USE rcf_grid_type_mod,        ONLY: output_grid

USE rcf_headaddress_mod,      ONLY: &
    fh_vtyear,                      &
    fh_vtmonth,                     &
    fh_vtday,                       &
    fh_vthour,                      &
    fh_vtminute,                    &
    fh_vtsecond

USE rcf_items_mod,            ONLY: &
    upaf_array

USE rcf_locate_mod,           ONLY: rcf_locate

USE rcf_umhead_mod,           ONLY: um_header_type

USE t_int_mod,                ONLY: t_int

USE UM_ParVars,               ONLY: &
    current_decomp_type,            &
    change_decomposition

USE umPrintMgr,               ONLY: &
    umPrint,                        &
    umMessage,                      &
    PrintStatus,                    &
    PrStatus_Normal,                &
    PrStatus_Diag

USE yomhook,                  ONLY: &
    lhook,                          &
    dr_hook

IMPLICIT NONE

! Arguments

TYPE (um_header_type),      INTENT(IN) :: hdr_out
TYPE (field_type), POINTER             :: fields_out (:)
INTEGER ,                   INTENT(IN) :: field_count_out

! ----------------------------------------------------------------
! Local Variables

! GLOBAL ATTRIBUTE 'update_type' read in from NetCDF file !
! 0 Single time
! 1 Time series
! 2 Periodic time series
INTEGER            :: ncdf_update_type (max_ncdf_files)
CHARACTER (LEN=1)  :: char_ncdf_update_type

INTEGER  :: pos
INTEGER  :: varid           ! Variable loop index
INTEGER  :: l               ! Dummy loop index
INTEGER  :: i,j,k,z         ! Loop indices
INTEGER  :: ndims           ! Nr of dims in this variable
INTEGER  :: IOStatus        ! Return code from 'Inquire'
INTEGER  :: ErrorStatus
INTEGER  :: orig_decomp

INTEGER  :: nDimensions     ! no. dimensions of dataset
INTEGER  :: nVariables      ! no. variables of dataset
INTEGER  :: nAttributes     ! no. global attributes of dataset
INTEGER  :: unlimitedDimId  ! returned ID of the unlimited dimension
                            ! unlimitedDimId returns -1 if there
                            ! are no unlimited dimensions

INTEGER  :: stashcode_att_array(1) ! VARIABLE ATTRIBUTE from file

INTEGER  :: xdimsize(max_dims_per_var) ! Size of dims in this variable

INTEGER  :: ireco(2)          ! NetCDF time record indices: prev & next
INTEGER  :: creco(2,6)        ! Time records: 1st dim is for prev & next
                              ! 2nd dim is for yr, mon, day, hr, min, sec

LOGICAL  :: l_exact_match     ! T if no time interpolation required
LOGICAL  :: l_exist           ! T if exists.
LOGICAL  :: l_first           ! T if reading first netcdf field
LOGICAL  :: lfoundstash   
LOGICAL  :: l_noreco
LOGICAL  :: multi_pe_in       ! flag to confirm 'gather' req'd

REAL     :: ftreco(2)         ! Time of prev/next time records (fract
                              ! hours since 0:00 h on 1st Jan)

REAL     :: fhournow          ! Current time (fractional hours since
                              ! 0:00 h on 1st Jan)

REAL, ALLOCATABLE  :: tmpdat1    (:,:)
REAL, ALLOCATABLE  :: tmpdat2    (:,:)
REAL, ALLOCATABLE  :: fielddata1 (:,:)
REAL, ALLOCATABLE  :: fielddata2 (:,:)

CHARACTER (LEN=netcdf_varname_len)   :: variable_name ! Name of NetCDF variable
CHARACTER (LEN=errormessagelength)   :: cmessage
CHARACTER (LEN=*), PARAMETER         :: RoutineName='RCF_NETCDF_ATMOS'

INTEGER (KIND=jpim), PARAMETER       :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER       :: zhook_out = 1
REAL    (KIND=jprb)                  :: zhook_handle

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! ===============================================
! Setup output decomposition for writing.
orig_decomp = current_decomp_type
IF (orig_decomp /= decomp_rcf_output) THEN
  CALL Change_Decomposition( decomp_rcf_output )
END IF

! ===============================================
! Set up list of NetCDF fields to be read in
! erroreport if some dependencises are missing.

l_first = .TRUE.
DO z=1, number_ncdf_files
  !---------------------------
  !    Open the netcdf file
  !---------------------------
  CALL netcdf_fopen (ncdf_list(z)%ncdf_file,ncdf_file_unit_no)
  
  !    Get Global attribute 'update_type'
  !    We adopt the same conventions as for ancillary files:
  !    0 Single time , 1 Time series , 2 Periodic time series
  l_exist = .FALSE.
  CALL netcdf_get_var_att (ncdf_file_unit_no, 1, 'GLOBAL', 'update_type',     &
       char_ncdf_update_type, l_exist)
  IF ( .NOT. l_exist) THEN
    ErrorStatus = 54003
    cmessage = 'Global attribute update_type not present in file ' //         &
                TRIM(ncdf_list(z)%ncdf_file)
    CALL ereport( RoutineName, ErrorStatus, cmessage )
  END IF

  READ (char_ncdf_update_type, '(I1)') ncdf_update_type(z)
  !  Check the update indicator is valid
  IF ( (ncdf_update_type(z) /= 0)  .AND.                                      &
       (ncdf_update_type(z) /= 1)  .AND.                                      &
       (ncdf_update_type(z) /= 2) ) THEN
    ErrorStatus = 54500 + ncdf_update_type(z)
    cmessage = 'Single time (0) Time series (1) and Periodic (2) are ' //     &
               'only options allowed. File: ' // TRIM(ncdf_list(z)%ncdf_file)
    CALL ereport( RoutineName, ErrorStatus, cmessage )
  END IF
  
  !---------------------------
  !    Determine the number of variables in NetCDf file 
  !---------------------------
  CALL netcdf_get_file_info ( ncdf_file_unit_no, nDimensions, nVariables,     &
                              nAttributes, unlimitedDimId) 
  
  ! Loop over all fields within file ncdf_list(z)%ncdf_file
  DO i=1,ncdf_list(z)%fields_in_file
    !---------------------------
    !    Loop over number of Variables in NetCDF file
    !---------------------------
    !  Go through all variables in netcdf file to find the variable
    !  which matches required stashcode or else has matching varname
    lfoundstash = .FALSE.
    varid_nVariables:DO varid = 1, nVariables
      ! Set variable_name to ' ' so that NETCDF_GET_VAR_INFO returns
      ! the variable name associated with the variable ID.
      variable_name = ' '
      l = varid       ! cannot pass loop counter as arg declared INOUT
      CALL netcdf_get_var_info (ncdf_file_unit_no, l, variable_name,          &
                                ndims=ndims)
      
      ! Skip 1D variables
      IF ( ndims < 2 ) CYCLE varid_nVariables ! Cycle loop over nVariables
      
      ! variable_name not allowed to be named "unset" in NetCDF file
      IF ( variable_name == netcdf_unset ) THEN
        ErrorStatus = 54999
        cmessage = 'Varname not allowed to be named "unset" in File: ' //     &
              TRIM( upaf_array(i) )
        CALL ereport( RoutineName, ErrorStatus, cmessage )
      END IF
      
      ! Some variables (e.g dimension related) will not contain the
      ! 'stashcode' attribute
      l_exist = .FALSE.
      CALL netcdf_get_var_att (ncdf_file_unit_no, 1, variable_name,           &
                              'stashcode', stashcode_att_array, l_exist)
      
      IF (l_exist) THEN
        
        IF ( ncdf_list(z)%stashcode(i) == stashcode_att_array(1) ) THEN
          lfoundstash = .TRUE.
          WRITE(umMessage,'(A,A)') 'Found variable: ', TRIM(variable_name)
          CALL umPrint(umMessage,src=RoutineName)
          EXIT varid_nVariables ! Exit loop over nVariables
        ELSE
          ! Try next variable
          CYCLE varid_nVariables ! Cycle loop over nVariables
        END IF
        
      ELSE
        
        ! Variable does not contain 'stashcode' attribute
        ! Try matching with ncdf_list(z)%varname(i) from namelist
        IF ( TRIM( variable_name ) == TRIM( ncdf_list(z)%varname(i) ) ) THEN
          WRITE(umMessage,'(A,A,A)') 'Found variable : ',                     &
                                      TRIM(variable_name),                    &
                                     'from "netcdf_varname" in the namelist.'
          CALL umPrint(umMessage,src=RoutineName)
          lfoundstash = .TRUE.
          EXIT  varid_nVariables ! Exit loop over nVariables
        ELSE
          ! Try next variable
          CYCLE varid_nVariables ! Cycle loop over nVariables
        END IF
        
      END IF ! (l_exist) 'stashcode'
    
    END DO varid_nVariables ! varid = 1, nVariables
    
    IF (.NOT. lfoundstash) THEN
      ! Not able to obtain information for this stashcode
      ErrorStatus = ncdf_list(z)%stashcode(i)
      cmessage = 'Stashcode not present in file ' //                          &
                  TRIM(ncdf_list(z)%ncdf_file)
      CALL ereport( RoutineName, ErrorStatus, cmessage )
    END IF
    
    !---------------------------
    !    Determine pos using rcf_locate
    !---------------------------
    CALL Rcf_Locate ( ncdf_list(z)%section(i), ncdf_list(z)%item(i),          &
                      fields_out, field_count_out, pos)
    
    ! ===============================================
    ! Time interploation stuff
    ! Returns:
    !  Time record indices,
    !  Readable time of records (last/next and yr,mon,day,hr,min,sec)
    !  Current model time in fract hours since base time (0:00 h on 1st Jan)
    !  Prev/next time in hours since the same base time
        
    CALL glomap_clim_get_netcdffile_rec ( output_grid%glob_p_row_length,      &
                                          output_grid%glob_p_rows,            &
                                          ncdf_file_unit_no,                  &
                                          i, ncdf_update_type(z),             &
                                          (/ hdr_out%fixhd(fh_vtyear),        &
                                          hdr_out%fixhd(fh_vtmonth),          &
                                          hdr_out%fixhd(fh_vtday),            &
                                          hdr_out%fixhd(fh_vthour),           &
                                          hdr_out%fixhd(fh_vtminute),         &
                                          hdr_out%fixhd(fh_vtsecond) /),      &
                                          l_first, ireco, creco, fhournow,    &
                                          ftreco, l_exact_match )
    
    ! Check that matching time records are found
    IF ( ANY(ireco(:) == -1) ) THEN  
      WRITE (umMessage,'(A,6(1x,I4),1x,F15.5)') 'At time:',                   &
                                               hdr_out%fixhd(fh_vtyear),      &
                                               hdr_out%fixhd(fh_vtmonth),     &
                                               hdr_out%fixhd(fh_vtday),       &
                                               hdr_out%fixhd(fh_vthour),      &
                                               hdr_out%fixhd(fh_vtminute),    &
                                               hdr_out%fixhd(fh_vtsecond),    &
                                               fhournow
      CALL umPrint(umMessage,src=RoutineName)
      CALL netcdf_fclose (ncdf_file_unit_no)
      cmessage = 'Required time not present in file ' // ncdf_list(z)%ncdf_file
      ErrorStatus = 54900 + i
      CALL ereport( RoutineName, ErrorStatus, cmessage )
    END IF
    
    ! ===============================================
    ! Check if variable dimensions match 
    CALL netcdf_var_check_dims (output_grid%glob_p_row_length,                &
                                output_grid%glob_p_rows,                      &
                                output_grid%model_levels, ncdf_file_unit_no,  &
                                variable_name, xdimsize)
    
    ALLOCATE ( fielddata1 ( ( output_grid%loc_p_row_length *                  &
                              output_grid%loc_p_rows ), xdimsize(3) ) )
    
    ! ===============================================
    ! Time interpolation
    
    l_noreco = .FALSE.  ! initially assume that records are found
    IF ( l_exact_match ) THEN
      !     No interpolation if exact time match. Get only values
      !     for register ireco(1).
      CALL netcdf_get_data (output_grid%glob_p_row_length,                    &
                            output_grid%glob_p_rows,                          &
                            output_grid%loc_p_row_length,                     &
                            output_grid%loc_p_rows,                           &
                            xdimsize(3), ncdf_file_unit_no,                   &
                            decomp_rcf_output,                                &
                            ireco(1), variable_name, fielddata1, l_noreco)
      
      IF (l_noreco) THEN  
        CALL netcdf_fclose (ncdf_file_unit_no)
        cmessage = 'Error reading data for variable: ' // TRIM(variable_name)
        ErrorStatus = 54800 + i
        CALL ereport( RoutineName, ErrorStatus, cmessage )
      END IF
    ELSE
      ALLOCATE ( tmpdat1 ( ( output_grid%loc_p_row_length *                   &
                             output_grid%loc_p_rows ), xdimsize(3) ) )
      
      ALLOCATE ( tmpdat2 ( ( output_grid%loc_p_row_length *                   &
                             output_grid%loc_p_rows ), xdimsize(3) ) )
      
      ! Read in data for previous matching time. Stop if record not found.
      CALL netcdf_get_data (output_grid%glob_p_row_length,                    &
                            output_grid%glob_p_rows,                          &
                            output_grid%loc_p_row_length,                     &
                            output_grid%loc_p_rows,                           &
                            xdimsize(3), ncdf_file_unit_no,                   &
                            decomp_rcf_output,                                &
                            ireco(1), variable_name, tmpdat1, l_noreco)
      
      IF (l_noreco) THEN  
        CALL netcdf_fclose (ncdf_file_unit_no)
        cmessage = 'Error reading data for variable: ' // TRIM(variable_name)
        ErrorStatus = 54700 + i
        CALL ereport( RoutineName, ErrorStatus, cmessage )
      END IF
      
      ! Read in data for next matching time. Stop if record not found.
      CALL netcdf_get_data (output_grid%glob_p_row_length,                    &
                            output_grid%glob_p_rows,                          &
                            output_grid%loc_p_row_length,                     &
                            output_grid%loc_p_rows,                           &
                            xdimsize(3), ncdf_file_unit_no,                   &
                            decomp_rcf_output,                                &
                            ireco(2), variable_name, tmpdat2, l_noreco)
      
      IF (l_noreco) THEN  
        CALL netcdf_fclose (ncdf_file_unit_no)
        cmessage = 'Error reading data for variable: ' // TRIM(variable_name)
        ErrorStatus = 54600 + i
        CALL ereport( RoutineName, ErrorStatus, cmessage )
      END IF
      
      ! Time interpolate
      CALL t_int ( tmpdat1, ftreco(1), tmpdat2, ftreco(2), fielddata1,        &
                   fhournow, output_grid%loc_p_row_length *                   &
                   output_grid%loc_p_rows * xdimsize(3) )
      
      DEALLOCATE(tmpdat2)
      DEALLOCATE(tmpdat1)
      
    END IF ! l_exact_match
    
    ! ===============================================
    ! Match field size with that expected by d1
    !  end game has a zeroth level to confuse things
        
    ALLOCATE ( fielddata2 (                                                   &
                    ( output_grid%loc_p_row_length *output_grid%loc_p_rows ), &
                      fields_out(pos)%levels ) )
    
    IF ( fields_out(pos)%levels == output_grid%model_levels + 1 ) THEN
      fielddata2(:,1) = 0.0
      DO k = 1, output_grid%model_levels
        fielddata2(:,k+1) = fielddata1(:,k)
      END DO
    ELSE
      fielddata2 = fielddata1
    END IF
    
    ! ===============================================
    ! Write the ancillary fields out to the dump
    multi_pe_in=.TRUE.
    ! DEPENDS ON: rcf_writflds
    CALL Rcf_writflds (hdr_out%unitnum,                                       &
                       fields_out(pos)%levels,                                &
                       fields_out(pos)%dump_pos,                              &
                       hdr_out%lookup,                                        &
                       hdr_out%len1lookup,                                    &
                       fielddata2,                                            &
                       fields_out(pos)%level_size,                            &
                       hdr_out%fixhd,                                         &
                       ErrorStatus,cmessage,multi_pe_in)
    
    DEALLOCATE(fielddata2)
    DEALLOCATE(fielddata1)
    
    IF (ErrorStatus /= 0) THEN
      WRITE(umMessage,'(A,I5)') ' Error in RCF_WRITFLDS for NetCDF field ', i
      CALL umPrint(umMessage,src=RoutineName)
      
      CALL ereport ( RoutineName, ErrorStatus, cmessage)
    END IF
    
    l_first = .FALSE.
  END DO
  
  !---------------------------
  !    Close the netcdf file
  !---------------------------
  CALL netcdf_fclose (ncdf_file_unit_no)
END DO

! Change back to original if different.
IF (current_decomp_type /= orig_decomp) THEN
  CALL Change_Decomposition( orig_decomp )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
END SUBROUTINE rcf_netcdf_atmos

END MODULE rcf_netcdf_atmos_mod
