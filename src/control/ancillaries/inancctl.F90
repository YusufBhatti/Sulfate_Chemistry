! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine inancctl in Module inancctl_mod
!
!   Programming standard; UMDP3
!
!   Purpose : Takes as input,the code defining the frequency of update
!             of ancillary fields as set by the user interface.
!             Converts them into a list of numbers of timesteps after
!             which each field must be updated, and calculates the
!             frequency with which this list must be interogated.
!             Where the update interval is in months or years,
!             the check will be carried out each day. The physical
!             files required are also determined by input code,
!             and the headers and lookup tables are read in.
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Ancillaries

MODULE inancctl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INANCCTL_MOD'

CONTAINS

SUBROUTINE inancctl(                                              &
                  icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE atm_fields_bounds_mod
USE UM_ParVars
USE decomp_params, ONLY: decomp_standard_atmos
USE Control_Max_Sizes
USE Decomp_DB
USE lookup_addresses
USE nlstgen_mod,  ONLY: steps_per_periodim, secs_per_periodim
USE dump_headers_mod, ONLY: a_fixhd, a_realhd, a_levdepc
USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY: nitems, si
USE nlstcall_mod, ONLY: model_basis_time, ancil_reftime, lcal360
USE umPrintMgr, ONLY: umPrint, umMessage

USE ancilcta_namelist_mod, ONLY: nancil_lookupsa

USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, len1_lookup,                            &
    len_dumphist, len_fixhd, model_levels, mpp_len1_lookup, n_cca_lev, &
    ntiles, ozone_levels, pp_len_inthd,                                &
    pp_len_realhd, sm_levels, st_levels, tpps_ozone_levels,            &
    tr_lbc_ukca, tr_lbc_vars, tr_levels, tr_ukca, tr_vars

USE ancil_mod, ONLY: max_items

USE model_time_mod, ONLY: ancillary_stepsim

USE ancil_headers_mod, ONLY: allocate_ancil_headers
USE inancila_mod, ONLY: inancila
USE populate_ancil_requests_mod, ONLY: populate_ancil_requests
USE items_nml_mod, ONLY: read_nml_items, netcdf_varname_len
USE ancilcta_namelist_mod, ONLY: read_nml_ancilcta

USE errormessagelength_mod, ONLY: errormessagelength
USE filenamelength_mod,     ONLY: filenamelength

USE file_manager, ONLY: get_file_unit_by_id

IMPLICIT NONE


!   Arguments
!

INTEGER :: icode      ! Out return code :0 Nor al Exit; :>0 Error

CHARACTER(LEN=errormessagelength) :: cmessage   ! Out error message if ICODE >0


! Local variables
INTEGER :: im_index       ! internal model index for STASH arrays
INTEGER :: i
INTEGER :: steps_per_hr   ! steps per hour for atmos sub_model
INTEGER :: decomp_type   ! domain decomposition type
INTEGER :: shared_unit

INTEGER  :: num_items
INTEGER  :: stashcode_temp ( max_items )
INTEGER  :: section_temp ( max_items )
INTEGER  :: item_temp ( max_items )
INTEGER  :: source_temp ( max_items )
INTEGER  :: domain_temp ( max_items )
INTEGER  :: user_prog_section_temp ( max_items )
INTEGER  :: user_prog_item_temp ( max_items )
INTEGER  :: period_temp( max_items )
INTEGER  :: interval_temp( max_items )
REAL     :: real_constant_temp ( max_items )
CHARACTER(LEN=filenamelength) :: items_filename_temp( max_items )
CHARACTER (LEN=netcdf_varname_len) :: netcdf_varname_temp ( max_items )
LOGICAL  :: update_anc_temp ( max_items )

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INANCCTL'

!  initialise reference time for time interpolation of ancillaries.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF ( ancil_reftime(1) == 0 .AND. ancil_reftime(2) == 0            &
  .AND. ancil_reftime(3) == 0 .AND. ancil_reftime(4) == 0         &
  .AND. ancil_reftime(5) == 0 .AND. ancil_reftime(5) == 0 ) THEN
  WRITE(umMessage,*)' ANCIL_REFTIME set same as MODEL_BASIS_TIME'
  CALL umPrint(umMessage,src='inancctl')
  ancil_reftime(1) = model_basis_time(1)
  ancil_reftime(2) = model_basis_time(2)
  ancil_reftime(3) = model_basis_time(3)
  ancil_reftime(4) = model_basis_time(4)
  ancil_reftime(5) = model_basis_time(5)
  ancil_reftime(6) = model_basis_time(6)
  WRITE(umMessage,*)' ANCIL_REFTIME = MODEL_BASIS_TIME = ',ancil_reftime
  CALL umPrint(umMessage,src='inancctl')
ELSE
  WRITE(umMessage,'(A,I10)')' ancil_reftime set by User Interface'
  CALL umPrint(umMessage,src='inancctl')
END IF
WRITE(umMessage,'(A,I4)') 'ancil_reftime(1) = ', ancil_reftime(1)
CALL umPrint(umMessage,src='inancctl')
WRITE(umMessage,'(A,I4)') 'ancil_reftime(2) = ', ancil_reftime(2)
CALL umPrint(umMessage,src='inancctl')
WRITE(umMessage,'(A,I4)') 'ancil_reftime(3) = ', ancil_reftime(3)
CALL umPrint(umMessage,src='inancctl')
WRITE(umMessage,'(A,I4)') 'ancil_reftime(4) = ', ancil_reftime(4)
CALL umPrint(umMessage,src='inancctl')
WRITE(umMessage,'(A,I4)') 'ancil_reftime(5) = ', ancil_reftime(5)
CALL umPrint(umMessage,src='inancctl')
WRITE(umMessage,'(A,I4)') 'ancil_reftime(6) = ', ancil_reftime(6)
CALL umPrint(umMessage,src='inancctl')

!  Set up internal model identifier and STASH index
im_index = 1
! Check that current decomposition is correct for ancillaries
! for this sub-model
decomp_type=decomp_standard_atmos
IF (current_decomp_type  /=  decomp_type) THEN
  CALL change_decomposition(decomp_type,icode)
  IF (icode  >   0) THEN
    WRITE(umMessage,*) 'INANCCTL : Error'
    CALL umPrint(umMessage,src='inancctl')
    WRITE(umMessage,*) 'Call to CHANGE_DECOMPOSITION failed with ',      &
               'decomposition type ',decomp_type
    CALL umPrint(umMessage,src='inancctl')
    cmessage='INANCCTL;Unsupported decomposition for MPP code'
    GO TO 9999
  END IF
END IF

steps_per_hr = 3600*steps_per_periodim(atmos_im)/                 &
                      secs_per_periodim(atmos_im)

!   Read in control information from namelist
shared_unit = get_file_unit_by_id("shared", handler="fortran")
CALL read_nml_ancilcta(shared_unit)
! Read in items namelist
CALL read_nml_items(shared_unit,           num_items,              &
                    stashcode_temp,        section_temp,           &
                    item_temp,             source_temp,            &    
                    domain_temp,           user_prog_section_temp, &   
                    user_prog_item_temp,   period_temp,            &
                    interval_temp,         real_constant_temp,     &
                    update_anc_temp,       items_filename_temp,    &
                    netcdf_varname_temp)

! Convert items namelist arrays into ancil request and ancil file objects
CALL populate_ancil_requests(stashcode_temp(1:num_items),      &
                             section_temp(1:num_items),        &
                             item_temp(1:num_items),           &
                             period_temp(1:num_items),         &
                             interval_temp(1:num_items),       &
                             items_filename_temp(1:num_items), &
                             update_anc_temp(1:num_items))

! Allocate arrays for ancil headers
CALL allocate_ancil_headers()

CALL inancila (len_fixhd,pp_len_inthd,pp_len_realhd,a_len1_levdepc,    &
               a_len2_levdepc, a_fixhd,a_realhd,a_levdepc,             &
               nancil_lookupsa,len1_lookup,           &
                     model_levels,tr_levels,                                 &
                     st_levels,sm_levels,ntiles,                             &
                     ozone_levels,tpps_ozone_levels,                         &
!                    all ancillaries assumed to be in section 0
                     si(1,0,im_index),                                       &
                                           ! si for atmos_im, sect 0
                     nitems,                                                 &
                     ancillary_stepsim(atmos_im),steps_per_hr,               &
                     icode,cmessage,lcal360)

IF (icode >  0) THEN
  WRITE(umMessage,*) 'INANCCTL: Error return from INANCILA ',icode
  CALL umPrint(umMessage,src='inancctl')
  WRITE(umMessage,*) cmessage
  CALL umPrint(umMessage,src='inancctl')
  GO TO 9999            ! Jump to end of routine
END IF

9999 CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE inancctl

END MODULE inancctl_mod
