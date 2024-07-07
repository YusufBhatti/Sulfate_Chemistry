! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_netcdf_init_mod

!  Subroutine rcf_netcdf_init - Identify climatology fields to read from netcdf
!
! Description:
!    Controls all NetCDF processing for the atmosphere model
!
! Method:
!    1. Generate list of NetCDF files and fields within these files
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!    This file belongs in section: Reconfiguration_NetCDF
!
! Code Description:
!    Language: FORTRAN 90
!    This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_NETCDF_INIT_MOD'

CONTAINS

SUBROUTINE rcf_netcdf_init ( )

USE ereport_mod,                      ONLY: &
    ereport

USE errormessagelength_mod,           ONLY: &
    errormessagelength

USE glomap_clim_netcdf_parameter_mod, ONLY: &
    number_ncdf_files,                      &
    ncdf_list

USE items_nml_mod,                    ONLY: &
    netcdf_file

USE parkind1,                         ONLY: &
    jpim,                                   &
    jprb

USE rcf_exppx_mod,                    ONLY: &
    rcf_exppx

USE rcf_field_type_mod,               ONLY: &
    field_type

USE rcf_items_mod,                    ONLY: &
    source_array,                           &
    sctn_array,                             &
    item_array,                             &
    upaf_array,                             &
    upnv_array,                             &
    num_items

USE submodel_mod,                     ONLY: &
    atmos_im

USE UM_ParCore,                       ONLY: &
    mype

USE um_stashcode_mod,                 ONLY: &
    stashcode_glomap_clim_sec

USE umPrintMgr,                       ONLY: &
    umPrint,                                &
    umMessage,                              &
    PrintStatus,                            &
    PrStatus_Normal

USE yomhook,                          ONLY: &
    lhook,                                  &
    dr_hook

IMPLICIT NONE

! Local Variables

INTEGER  :: cf                ! Current file
INTEGER  :: i,j               ! Loop indices
INTEGER  :: ErrorStatus       ! ErrorStatus

LOGICAL  :: l_item_processed  ! F if newly encountered file

LOGICAL  :: l_first = .TRUE.  ! Only print the message once

CHARACTER (LEN=errormessagelength)   :: cmessage
CHARACTER (LEN=*), PARAMETER         :: RoutineName='RCF_NETCDF_INIT'

TYPE (field_type)                    :: dummy_ncdf

INTEGER (KIND=jpim), PARAMETER       :: zhook_in  = 0
INTEGER (KIND=jpim), PARAMETER       :: zhook_out = 1
REAL    (KIND=jprb)                  :: zhook_handle

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! ==========================================
! Put information from ITEMS into ncdf_list
! ==========================================

DO i=1, num_items
  IF ( source_array(i) /= netcdf_file ) CYCLE ! Only source==10
  
  IF ( l_first ) THEN
    IF ( PrintStatus >= PrStatus_Normal .AND. mype == 0 ) THEN
      WRITE(umMessage,'(A,A)') RoutineName, ' - NetCDF files to be opened : '
      CALL umPrint(umMessage,src=RoutineName)
    END IF
    
    l_first = .FALSE.
  END IF
    
  l_item_processed = .FALSE. ! Reset flag
  
  !---------------------------------------------
  ! Loop over netcdf files already encountered
  !   will be skipped if number_ncdf_files==0
  !---------------------------------------------
  no_ncdf_file:DO j=1, number_ncdf_files
    IF ( TRIM(upaf_array(i)) == ncdf_list(j)%ncdf_file ) THEN
      ! If filename the same then add new stash to type
      ncdf_list(j)%fields_in_file = ncdf_list(j)%fields_in_file + 1
      cf                          = j
      l_item_processed            = .TRUE.
      EXIT no_ncdf_file
    END IF
  END DO no_ncdf_file
  
  !-------------------------------------------------------------------------
  ! If items request was not processed it must be a newly encountered file
  !-------------------------------------------------------------------------
  IF ( .NOT. l_item_processed) THEN
    
    !---------------------
    ! Populate ncdf_file
    !---------------------
    number_ncdf_files            = number_ncdf_files + 1
    cf                           = number_ncdf_files
    ncdf_list(cf)%fields_in_file = 1
    ncdf_list(cf)%ncdf_file      = TRIM(upaf_array(i))
    
    !-------------------------------------------------
    ! Print out which NetCDF files need to be opened
    !-------------------------------------------------
    IF ( PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
      WRITE (umMessage,'(A,I4,1X,A)') 'File No ', number_ncdf_files,          &
                                       TRIM( upaf_array(i) )
      CALL umPrint(umMessage,src=RoutineName)
    END IF
    
  END IF ! ( .NOT. l_item_processed)
  
  !------------------------------------------------------------
  ! Populate varname, stashcode, section, item and stash_name
  !------------------------------------------------------------
  ncdf_list(cf)%varname    ( ncdf_list(cf)%fields_in_file ) = upnv_array(i)
  ncdf_list(cf)%stashcode  ( ncdf_list(cf)%fields_in_file ) = ( sctn_array(i) &
                                                     * 1000 ) + item_array(i)
  ncdf_list(cf)%section    ( ncdf_list(cf)%fields_in_file ) = sctn_array(i)
  ncdf_list(cf)%item       ( ncdf_list(cf)%fields_in_file ) = item_array(i)
  
  !---------------------------------------------
  ! set up stash names for each item requested
  !---------------------------------------------
  dummy_ncdf%stashmaster => rcf_exppx (atmos_im, sctn_array(i), item_array(i))
  
  ncdf_list(cf)%stash_name ( ncdf_list(cf)%fields_in_file ) =                 &
                                                   dummy_ncdf%stashmaster%NAME
  
  ! Section 54 is a simple case
  !   Further testing/development is required for other prognostic sections
  IF ( ncdf_list(cf)%section(ncdf_list(cf)%fields_in_file) /=                 &
       stashcode_glomap_clim_sec ) THEN
    
    ErrorStatus = 54000 + ncdf_list(cf)%section(ncdf_list(cf)%fields_in_file)
    cmessage = 'Only Section 54 has been tested for reading ' //              &
               'NetCDF data to the reconfiguration dump'
    CALL ereport( RoutineName, ErrorStatus, cmessage )
    
  END IF
END DO ! Loop over num_items

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
END SUBROUTINE rcf_netcdf_init

END MODULE rcf_netcdf_init_mod
