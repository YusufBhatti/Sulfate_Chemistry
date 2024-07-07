! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_nlist_recon_technical_mod

! Description:
!   Defines variables for the RECON_TECHNICAL namelist along with
!   the associated subroutines for reading, printing and checking
!   the namelist variables.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v6 programming standards.
!
! Summary of contents:
!  * Module variables and namelist defintion
!  * print_nlist_recon_technical()
!  * rcf_assign_vars_recon_technical()
!  * rcf_check_recon_technical()
!  * read_nml_recon_technical()

USE missing_data_mod, ONLY: imdi

USE filenamelength_mod, ONLY: filenamelength

USE umPrintMgr, ONLY:       &
    umPrint,                &
    newline

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

PRIVATE

! Make namelist handling routines public.
PUBLIC :: read_nml_recon_technical,                                          &
          rcf_check_recon_technical,                                         &
          rcf_assign_vars_recon_technical,                                   &
          print_nlist_recon_technical

! Make namelist items public.
PUBLIC ::                                                                    &
          dump_pack,                                                         &
          input_dump_type,                                                   &
          select_output_fields,                                              &
          output_field_list,                                                 &
          reset_data_time,                                                   &
          l_trans,                                                           &
          var_recon,                                                         &
          l_validity_lookup_u,                                               &
          l_rcf_init_flexi,                                                  &
          l_basic_interp,                                                    &
          ainitial,                                                          &
          transp

! Make useful parameters public.
PUBLIC ::                                                                    &
          max_num_fields,                                                    &
          um_input_dump,                                                     &
          grib_input_dump,                                                   &
          grib2ff_input_dump,                                                &
          tstmsk_to_decide,                                                  &
          interp_all_fields,                                                 &
          defined_by_namelist,   &
          len_dumphist
          
CHARACTER (LEN=*), PARAMETER, PRIVATE :: &
  ModuleName = "RCF_NLIST_RECON_TECHNICAL_MOD"

! Maximum size of array of output field stashcodes
INTEGER, PARAMETER :: max_num_fields = 50

!----------------------------------------------------
! Items set from RECON_TECHNICAL namelist
!----------------------------------------------------

! Dumping packing option
! 1 => STASHmaster controlled packing for diagnostic and primary fields.
! 2 => Unpacked primary fields. STASHmaster-packed diagnostics.
! 3 => Unpacked primary and diagnostic fields
INTEGER       :: dump_pack = imdi

! File type of input dump
INTEGER       :: input_dump_type = imdi
! 1=> Standard UM input dump
! 2=> ECMWF GRIB format
! 3=> Output from grib2ff program

! Reset data time to verification time
LOGICAL       :: reset_data_time = .FALSE.  

! Transplating of data
LOGICAL       :: l_trans = .FALSE.

! Run the reconfiguration in the VAR suite
LOGICAL       :: var_recon = .FALSE.

! Set the validity time of all lookup headers to be
! the validity time of the input U field.  This may
! be neccessary for VAR reconfigurations where one 
! cannot guarantee that input fields are instantaneous
LOGICAL       :: l_validity_lookup_u

! Flexi tiles: Perform general initialisation of tiled land surface fields
LOGICAL       :: l_rcf_init_flexi = .FALSE.

! Select how the reconfiguration decides which fields
! to process from the input dump
INTEGER       :: select_output_fields = imdi
! 0 => allow tstmsk routine to decide
! 3 => interp all fields found in input dump
! 4 => use array of namelists read in from namelist to 
!      determine which fields are present in output dump      

! Input dump filename
CHARACTER(LEN=filenamelength) :: ainitial = 'unset'

! Transplanted data input filename
CHARACTER(LEN=filenamelength) :: transp   = 'unset'

! Perform basic interpolation withoutout pre/post interpolation
! processing and checks
LOGICAL       :: l_basic_interp = .FALSE.

! List of stashcodes to include in the output dump
INTEGER       :: output_field_list(max_num_fields) = imdi

NAMELIST /recon_technical/                             &
 dump_pack, input_dump_type, reset_data_time, l_trans, &
 var_recon,select_output_fields, ainitial,             &
 transp, l_validity_lookup_u, l_rcf_init_flexi,        &
 l_basic_interp, output_field_list

!----------------------------------------------------
! Items set in module
!----------------------------------------------------

! No. of convective levels
INTEGER       :: conv_levels = 0

! Length of history File
INTEGER       :: len_dumphist = 0

! 1=> Standard UM input dump
INTEGER, PARAMETER :: um_input_dump = 1

! 2=> ECMWF GRIB format
INTEGER, PARAMETER :: grib_input_dump = 2

! 3=> Output from grib2ff program
INTEGER, PARAMETER :: grib2ff_input_dump = 3

! Parameters for the select_output_fields integer selector:
!
! Standard option of allowing reconfiguration to use
! the tstmsk routine to cross reference the science schemes
! and options chosen in the namelists with the version masks 
! and the option codes for each field from the STASHmaster.
INTEGER, PARAMETER :: tstmsk_to_decide = 0

! Interpolate all valid fields found in the input dump
INTEGER, PARAMETER :: interp_all_fields = 3

! Use array of stashcodes to determine which fields should
! be present in the output dump
INTEGER, PARAMETER :: defined_by_namelist = 4

CONTAINS

!-----------------------------------------------------------------------

SUBROUTINE print_nlist_recon_technical()

IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
INTEGER :: i
CHARACTER (LEN=*), PARAMETER  :: RoutineName='PRINT_NLIST_RECON_TECHNICAL'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL umPrint('Contents of namelist recon_technical', &
    src='rcf_nlist_recon_technical_mod')

WRITE(lineBuffer,'(A,I0)')' dump_pack = ',dump_pack
CALL umPrint(lineBuffer,src='rcf_nlist_recon_technical_mod')
WRITE(lineBuffer,'(A,I0)')' input_dump_type = ',input_dump_type
CALL umPrint(lineBuffer,src='rcf_nlist_recon_technical_mod')
WRITE(lineBuffer,'(A,I8)')' select_output_fields = ',select_output_fields
CALL umPrint(lineBuffer,src='rcf_nlist_recon_technical_mod')
WRITE(lineBuffer,'(A,L1)')' l_validity_lookup_u = ',l_validity_lookup_u
CALL umPrint(lineBuffer,src='rcf_nlist_recon_technical_mod')
WRITE(lineBuffer,'(A,L1)')' reset_data_time = ',reset_data_time
CALL umPrint(lineBuffer,src='rcf_nlist_recon_technical_mod')
WRITE(lineBuffer,'(A,L1)')' l_trans = ',l_trans
CALL umPrint(lineBuffer,src='rcf_nlist_recon_technical_mod')
WRITE(lineBuffer,'(A,L1)')' Var_Recon = ',var_recon
CALL umPrint(lineBuffer,src='rcf_nlist_recon_technical_mod')
WRITE(lineBuffer,'(A,A)')' ainitial = ',TRIM(ainitial)
CALL umPrint(lineBuffer,src='rcf_nlist_recon_technical_mod')
WRITE(lineBuffer,'(A,A)')' transp = ',TRIM(transp)
CALL umPrint(lineBuffer,src='rcf_nlist_recon_technical_mod')
WRITE(lineBuffer,'(A,L1)')' l_rcf_init_flexi = ',l_rcf_init_flexi
CALL umPrint(lineBuffer,src='rcf_nlist_recon_technical_mod')
WRITE(lineBuffer,'(A,L1)')' l_basic_interp = ',l_basic_interp
CALL umPrint(lineBuffer,src='rcf_nlist_recon_technical_mod')
IF (select_output_fields == defined_by_namelist) THEN
  DO i = 1, max_num_fields
    WRITE(lineBuffer,'(A,I0,A,I0)')' select_output_fields(',i,') = ', &
         output_field_list(i)
    CALL umPrint(lineBuffer,src='rcf_nlist_recon_technical_mod')
  END DO
END IF
CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='rcf_nlist_recon_technical_mod')
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE print_nlist_recon_technical

!-----------------------------------------------------------------------

SUBROUTINE rcf_assign_vars_recon_technical()
USE Rcf_Grid_Type_Mod, ONLY: Output_Grid
IMPLICIT NONE

! Set variables as required
Output_Grid % conv_levels   = conv_levels

END SUBROUTINE rcf_assign_vars_recon_technical

!-----------------------------------------------------------------------

SUBROUTINE rcf_check_recon_technical()

USE chk_opts_mod,     ONLY: &
    chk_var,                &
    def_src

USE Ereport_Mod, ONLY: &
    Ereport

USE nlcfiles_namelist_mod, ONLY: astart

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Local variables
CHARACTER (LEN=errormessagelength) :: chk_string
CHARACTER (LEN=errormessagelength) :: Cmessage
INTEGER :: ErrorStatus
LOGICAL :: pass_log(3) = .FALSE.
LOGICAL :: found_non_field = .FALSE.
INTEGER :: i

CHARACTER (LEN=*),  PARAMETER :: RoutineName='RCF_CHECK_RECON_TECHNICAL'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
INTEGER,            PARAMETER :: max_poss_STASH_ID = 99999
INTEGER,            PARAMETER :: min_poss_STASH_ID = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = ModuleName//':'//RoutineName

IF (select_output_fields == interp_all_fields    .OR.                        &
    select_output_fields == defined_by_namelist) THEN
  IF (input_dump_type == grib_input_dump) THEN
    ErrorStatus = 11
    WRITE (cmessage, '(A)')'Interpolating ALL fields found in input ' //      &
         'dump or ' // newline // 'selecting subset of fields using ' //      &
         'namelist not supported with GRIB  '// newline // 'input dump.'
    CALL ereport(TRIM(ModuleName//':'//RoutineName), ErrorStatus, cmessage)
  END IF

  IF (input_dump_type == grib2ff_input_dump) THEN
    ErrorStatus = 22
    WRITE (cmessage, '(A)')'Interpolating ALL fields found in input ' //      &
         'dump or ' // newline // 'selecting subset of fields using ' //      &
         'namelist not tested using GRIB2FF '// newline // 'output ' //       &
         'as input dump.'
    CALL ereport(TRIM(ModuleName//':'//RoutineName), ErrorStatus, cmessage)
  END IF
END IF

IF ( TRIM(astart) == TRIM(ainitial) ) THEN
  ErrorStatus = 33
  WRITE (cmessage, '(A)') "Input and output dump are the same file." //&
       newline // &
       "Input Dump (ainitial): " // TRIM(ainitial) // newline // &
       "Output Dump (astart): " // TRIM(astart)
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF ( select_output_fields == tstmsk_to_decide .AND. &
     l_basic_interp ) THEN
  ErrorStatus = 44
  WRITE(cmessage, '(A,I0,A)') "If using select_output_fields = ", &
       tstmsk_to_decide, newline // &
       "then l_basic_interp MUST be false. Please see metadata help for"// & 
       newline // &
       "namelist item select_output_fields"
  CALL ereport(TRIM(ModuleName//':'//RoutineName), ErrorStatus, cmessage)
END IF

WRITE(chk_string, '(A,I0,A,I0,A,I0,A)')'[',grib2ff_input_dump,',', & 
     grib_input_dump,',',um_input_dump,']'
CALL chk_var(input_dump_type, "input_dump_type", chk_string, &
             report_pass=pass_log(1) )
CALL chk_var(dump_pack, "dump_pack", "[1,2,3]", report_pass=pass_log(2) )
WRITE(chk_string, '(A,I0,A,I0,A,I0,A)')'[',tstmsk_to_decide,',', & 
     interp_all_fields,',',defined_by_namelist,']'
CALL chk_var(select_output_fields, "select_output_fields", &
             chk_string, report_pass=pass_log(3) )

IF (ANY(.NOT. pass_log)) THEN
  ErrorStatus = 55
  WRITE(cmessage, '(I0,A)') COUNT(.NOT. pass_log),                       &
       ' errors have been found in the recon_technical namelist'         &
       // newline //                                                     &
       'Please see output for details and then correct the namelist'
  CALL ereport(TRIM(ModuleName//':'//RoutineName), ErrorStatus, cmessage)
END IF

! Check output_field_list has only potentially valid STASH numbers with the
! rest of the array padded consistently with imdi
DO i=1, max_num_fields
  IF (output_field_list(i) > max_poss_STASH_ID .OR. &
      output_field_list(i) < min_poss_STASH_ID )    THEN
    IF ( output_field_list(i) == imdi ) THEN
      found_non_field = .TRUE.
    ELSE
      ErrorStatus = 45
      WRITE(cmessage, '(4(A,I0))') "Found STASH no. outside possible " //     &
           "range [", min_poss_STASH_ID, ":", max_poss_STASH_ID, "] "         &
           // newline // "in list of output fields to put in dump."           &
           // newline //                                                      &
           "output_field_list(", i ,") has value: ", output_field_list(i)
      CALL ereport(TRIM(ModuleName//':'//RoutineName), ErrorStatus, cmessage)
    END IF
  ELSE
    IF ( found_non_field ) THEN
      ErrorStatus = 47
      WRITE(cmessage, '(2(A,I0))') "Found possible stashcode after an " //    &
           "imdi in list of output fields to put in dump." // newline //      &
           "Please check list provided for an imdi." // newline //            &
           "Otherwise there has been an error reading the list." // newline //&
           "output_field_list(", i ,") has value: ", output_field_list(i)
      CALL ereport(TRIM(ModuleName//':'//RoutineName), ErrorStatus, cmessage)
    END IF
  END IF
END DO

def_src = ''

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_check_recon_technical

!-----------------------------------------------------------------------

SUBROUTINE read_nml_recon_technical(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat
USE errormessagelength_mod, ONLY: errormessagelength
USE setup_namelist, ONLY: setup_nml_type
USE filenamelength_mod, ONLY: filenamelength

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode

CHARACTER(LEN=errormessagelength) :: iomessage

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 3 + max_num_fields
INTEGER, PARAMETER :: n_log = 6
INTEGER, PARAMETER :: n_chars = 2 * filenamelength

TYPE my_namelist
  SEQUENCE
  INTEGER :: dump_pack
  INTEGER :: input_dump_type
  INTEGER :: select_output_fields
  INTEGER :: output_field_list(max_num_fields)
  LOGICAL :: reset_data_time
  LOGICAL :: l_trans
  LOGICAL :: var_recon
  LOGICAL :: l_validity_lookup_u
  LOGICAL :: l_rcf_init_flexi
  LOGICAL :: l_basic_interp
  CHARACTER(LEN=filenamelength)  :: ainitial
  CHARACTER(LEN=filenamelength)  :: transp
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

CHARACTER (LEN=*), PARAMETER  :: RoutineName='READ_NML_RECON_TECHNICAL'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                    n_log_in=n_log, n_chars_in=n_chars)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=recon_technical, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist recon_technical", iomessage)

  my_nml % dump_pack                 = dump_pack
  my_nml % input_dump_type           = input_dump_type
  my_nml % select_output_fields      = select_output_fields
  my_nml % output_field_list         = output_field_list
  my_nml % reset_data_time           = reset_data_time
  my_nml % l_trans                   = l_trans
  my_nml % var_recon                 = var_recon
  my_nml % l_validity_lookup_u       = l_validity_lookup_u
  my_nml % l_rcf_init_flexi          = l_rcf_init_flexi
  my_nml % l_basic_interp            = l_basic_interp
  my_nml % ainitial                  = ainitial
  my_nml % transp                    = transp

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  dump_pack                 = my_nml % dump_pack
  input_dump_type           = my_nml % input_dump_type
  select_output_fields      = my_nml % select_output_fields
  output_field_list         = my_nml % output_field_list
  reset_data_time           = my_nml % reset_data_time
  l_trans                   = my_nml % l_trans
  var_recon                 = my_nml % var_recon
  l_validity_lookup_u       = my_nml % l_validity_lookup_u
  l_rcf_init_flexi          = my_nml % l_rcf_init_flexi
  l_basic_interp            = my_nml % l_basic_interp
  ainitial                  = my_nml % ainitial
  transp                    = my_nml % transp

END IF

CALL mpl_type_free(mpl_nml_type,icode)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_nml_recon_technical


END MODULE rcf_nlist_recon_technical_mod
