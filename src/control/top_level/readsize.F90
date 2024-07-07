! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Subroutine READSIZE --------------------------------------------
!
!      Purpose:
!    This routine reads the length required for the dimensioning
!    of all main data arrays that have dimensions set via questions
!    in the USER INTERFACE.
!    It is called by the shell program UM_SHELL and passes the
!    data via modules to that routine. Data is then passed by
!    argument lists to allow for proper dynamic allocation in the
!    portable model.
!    The unit on which the start dump is opened is stored for use
!    elsewhere.
!
!   Programming standard:
!
!   ------------------------------------------------------------------
!    Interface and arguments: ------------------------------------------
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Top Level

SUBROUTINE readsize(dump_unit)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
! JULES
USE Check_IOStat_mod
USE ereport_mod, ONLY: ereport
USE Atmos_Max_Sizes
USE UM_ParParams
USE lbc_mod
USE io, ONLY: buffin
USE io_constants, ONLY: ioOpenReadOnly, ioNameProvided, ioDefault, ioAllLocal
USE model_file
USE UM_ParVars
USE UM_ParCore, ONLY: nproc
USE calc_ntiles_mod, ONLY: calc_ntiles
USE jules_surface_mod, ONLY: l_aggregate
USE jules_sea_seaice_mod, ONLY: nice, nice_use, l_sice_multilayers, l_ctile
USE jules_surface_types_mod, ONLY: npft, nnvg
USE land_tile_ids, ONLY: set_tile_id_arrays
USE mym_option_mod, ONLY: tke_levels, shcu_levels
USE bl_option_mod, ONLY: i_bl_vn, i_bl_vn_1a
USE missing_data_mod, ONLY: rmdi
USE rad_input_mod, ONLY: cusack_aero_hgt, aero_bl_levels
USE cloud_inputs_mod, ONLY: rhcrit
! module for aerosol emissions options
USE run_aerosol_mod, ONLY: run_aerosol_check

USE nlcfiles_namelist_mod, ONLY: astart
USE history, ONLY: history_file_exists, checkpoint_dump_im, blank_file_name
USE filenamelength_mod, ONLY: filenamelength
USE submodel_mod, ONLY: atmos_im

USE file_manager, ONLY: assign_file_unit, release_file_unit

USE nlsizes_namelist_mod, ONLY:                                          &
    a_len2_coldepc, a_len1_flddepc, a_len2_flddepc, a_len2_levdepc,      &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst,   &
    a_len_inthd, a_len_realhd, bl_levels, cloud_levels,                  &
    global_row_length, global_rows, land_field, model_levels, ntiles,    &
    ozone_levels, pp_len_inthd, pp_len_realhd, river_row_length,         &
    river_rows, sm_levels, st_levels, tpps_ozone_levels, tr_levels,      &
    tr_vars

USE errormessagelength_mod, ONLY: errormessagelength

USE io_configuration_mod, ONLY: io_alltoall_readflds

USE umPrintMgr, ONLY: umPrint, PrStatus_Oper

IMPLICIT NONE

!  Local variables

INTEGER, INTENT(OUT) :: dump_unit ! Unit number for ASTART
INTEGER                       :: errorstatus
CHARACTER (LEN=*), PARAMETER  :: RoutineName='READSIZE'
CHARACTER (LEN=errormessagelength)            :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


INTEGER, POINTER ::   fixhdr(:)
INTEGER :: len_io
INTEGER :: ioLocality
REAL ::   a
INTEGER, ALLOCATABLE  :: inthdr(:)
CHARACTER(LEN=filenamelength) :: dump_filename


!  1.1 Read size information from namelists

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!     Initialise TPPS_OZONE_LEVELS from OZONE_LEVELS
!     Prevents failures when using 'preset=nan'
tpps_ozone_levels = 0

! Read restart file header to get the value of the size of a number of arrays
! The name of the file to use is  based on whether it is a crun or an nrun. For
! an nrun, history_file_exists is false, for a crun history_file_exists is
! true. For nrun use the astart file and for a crun use the file specified
! by checkpoint_dump_im instead.
IF (history_file_exists .AND. &
    TRIM(checkpoint_dump_im(atmos_im)) /= blank_file_name) THEN
  dump_filename = checkpoint_dump_im(atmos_im)
ELSE
  dump_filename = astart
END IF

CALL assign_file_unit(dump_filename, dump_unit, handler="portio")

! In order to use the all-to-all method, each rank must be able to independently
! read different parts of the dump. This requires the "ioAllLocal" setting for
! ioLocality
IF (io_alltoall_readflds .AND. nproc>1 ) THEN
  CALL umPrint('Opening Start Dump using the all-to-all method',               &
               level=PrStatus_Oper,src=RoutineName)
  ioLocality = ioAllLocal
ELSE
  CALL umPrint('Opening Start Dump using the standard method',                 &
               level=PrStatus_Oper,src=RoutineName)
  ioLocality = ioDefault
END IF

CALL model_file_open(dump_unit, dump_filename,                  &
                     read_write=ioOpenReadOnly,                 &
                     name_in_environ=ioNameProvided,            &
                     ioLocality=ioLocality)

CALL loadHeader(dump_unit,FixedHeader)
fixhdr=>attachHeader(dump_unit,FixedHeader)

ALLOCATE (inthdr(fixhdr(101)))
IF (fixhdr(100) >  0) THEN
  CALL buffin(dump_unit,inthdr,fixhdr(101))
END IF

global_ROW_LENGTH=inthdr(6)
global_ROWS=inthdr(7)
model_levels=inthdr(8)
! inthdr(21) contains the value used for "mdi"
IF (inthdr(25) == inthdr(21)) THEN
  land_field=0
ELSE
  land_field=inthdr(25)
END IF
IF (inthdr(11) == inthdr(21)) THEN
  cloud_levels=model_levels
ELSE
  cloud_levels=inthdr(11)
END IF
IF (inthdr(10) == inthdr(21)) THEN
  st_levels=0
ELSE
  st_levels=inthdr(10)
END IF
IF (inthdr(28) == inthdr(21)) THEN
  sm_levels=0
ELSE
  sm_levels=inthdr(28)
END IF
IF (inthdr(13) == inthdr(21)) THEN
  bl_levels=0
ELSE
  bl_levels=inthdr(13)
END IF
IF (inthdr(26) == inthdr(21)) THEN
  ozone_levels=0
ELSE
  ozone_levels=inthdr(26)
END IF
IF (inthdr(12) == inthdr(21)) THEN
  tr_levels=0
ELSE
  tr_levels=inthdr(12)
END IF
IF (inthdr(14) == inthdr(21)) THEN
  tr_vars=0
ELSE
  tr_vars=inthdr(14)
END IF
IF (inthdr(19) == inthdr(21)) THEN
  river_row_length=0
ELSE
  river_row_length=inthdr(19)
END IF
IF (inthdr(20) == inthdr(21)) THEN
  river_rows=0
ELSE
  river_rows=inthdr(20)
END IF

a_len_inthd=fixhdr(101)
a_len_realhd=fixhdr(106)
pp_len_inthd=fixhdr(101)
pp_len_realhd=fixhdr(106)
IF (fixhdr(112) == inthdr(21)) THEN
  a_len2_levdepc=0
ELSE
  a_len2_levdepc=fixhdr(112)
END IF
IF (fixhdr(117) == inthdr(21)) THEN
  a_len2_rowdepc=0
ELSE
  a_len2_rowdepc=fixhdr(117)
END IF
IF (fixhdr(122) == inthdr(21)) THEN
  a_len2_coldepc=0
ELSE
  a_len2_coldepc=fixhdr(122)
END IF

! Read fields of constants sizes from header
IF (fixhdr(126) == inthdr(21)) THEN
  a_len1_flddepc=0
ELSE
  a_len1_flddepc=fixhdr(126)
END IF
IF (fixhdr(127) == inthdr(21)) THEN
  a_len2_flddepc=0
ELSE
  a_len2_flddepc=fixhdr(127)
END IF

IF (fixhdr(131) == inthdr(21)) THEN
  a_len_extcnst=0
ELSE
  a_len_extcnst=fixhdr(131)
END IF
IF (fixhdr(141) == inthdr(21)) THEN
  a_len_cfi1=0
ELSE
  a_len_cfi1=fixhdr(141)
END IF
IF (fixhdr(143) == inthdr(21)) THEN
  a_len_cfi2=0
ELSE
  a_len_cfi2=fixhdr(143)
END IF
IF (fixhdr(145) == inthdr(21)) THEN
  a_len_cfi3=0
ELSE
  a_len_cfi3=fixhdr(145)
END IF

! Check Wet levels in dump is same as model_levels
IF (inthdr(9) /=  model_levels) THEN
  errorstatus = 35
  WRITE (Cmessage, '(A, I3, A, I3)')                                          &
          'No. of wet levels in dump not equal to no. of model levels'        &
          , inthdr(9) , " != ", model_levels
  CALL ereport( routinename, errorstatus, cmessage )
END IF

DEALLOCATE(inthdr)

! Check that dynamically provided sizes are not larger than their
! statically defined maximums.
IF (global_row_length > row_length_max) THEN
  errorstatus = 10
  cmessage = 'Row length is larger than maximum defined in'       &
           //' Atmos_Max_Sizes'

  CALL ereport( routinename, errorstatus, cmessage )
END IF

IF (global_rows > rows_max) THEN
  errorstatus = 20
  cmessage = 'Number of rows is larger than maximum defined in'   &
           //' Atmos_Max_Sizes'

  CALL ereport( routinename, errorstatus, cmessage )
END IF

IF (model_levels > model_levels_max) THEN
  errorstatus = 30
  cmessage = 'Number of levels is larger than maximum defined in' &
           //' Atmos_Max_Sizes'

  CALL ereport( routinename, errorstatus, cmessage )
END IF

CALL calc_ntiles(l_aggregate,npft,nnvg,ntiles)
CALL set_tile_id_arrays ( )

! In preparation for NICE_USE being set in namelist
IF (nice_use /= 1 .AND. nice_use /= nice ) THEN
  errorstatus = 40
  cmessage = 'Problem with sea ice category settings' &
           //'NICE_USE must either equal 1 or NICE'

  CALL ereport( routinename, errorstatus, cmessage )
END IF

! Checks for sea ice settings - in preparation for these logicals being
! set in namelist. Perform this check in readsize() after nice has been
! read in by namelist and l_ctile and l_sice_multilayers have been read
! in from JULES_SEA_SEAICE in readlsta()

IF (.NOT. l_ctile .AND. nice_use > 1) THEN
  errorstatus = 10
  cmessage = 'Problem with sea ice settings:'                  &
       // 'Either set l_ctile=T or nice_use=1'
  CALL ereport( routinename, errorstatus, cmessage )
END IF
IF (l_sice_multilayers .AND. nice_use /= nice) THEN
  errorstatus = 20
  cmessage = 'Problem with sea ice settings:'                   &
       //'Either set l_sice_multilayers=F or nice_use=nice'
  CALL ereport( routinename, errorstatus, cmessage )
END IF

IF ( i_bl_vn == i_bl_vn_1a ) THEN
  ! A negative value for "tke_levels" means it should default to bl_levels.
  IF (tke_levels < 0 .OR. tke_levels > bl_levels) THEN
    tke_levels = bl_levels
    cmessage = 'WARNING: the value of tke_levels has been '       &
             //'reset to bl_levels'
    errorstatus = -1
    CALL ereport(RoutineName, errorstatus, cmessage)
  END IF

  ! A negative value for "shcu_levels" means it should default to tke_levels.
  IF (shcu_levels < 0) THEN
    shcu_levels = tke_levels
    cmessage = 'WARNING: the value of shcu_levels has been '      &
             //'reset to tke_levels'
    errorstatus = -1
    CALL ereport(RoutineName, errorstatus, cmessage)
  END IF
END IF

IF (ANY(rhcrit(1:cloud_levels) == rmdi)) THEN
  WRITE (cmessage,'(A62)') 'Unset values of RHcrit detected. '      &
                         //'Please check size in namelist.'
  errorstatus = 98
  CALL ereport(RoutineName, errorstatus, cmessage)
END IF

! Set radiation aerosol switch that depends on bl_levels
IF (cusack_aero_hgt==2) aero_bl_levels = bl_levels

! Check run_aerosol - needs model_levels
CALL run_aerosol_check

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE readsize
