! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Description
! When a stream is due to be reinitialised, the old file is closed,
! the post-processing system is informed if requested, and a new file
! is initialised.
!
!  Method : 1. Files to be processed on each call are controlled by
!              the file unit reinitialisation switch set by ppctl_init 
!              and at regular intervals thereafter.
!           2. Completed files are closed, and if post-processing is
!              selected a zero-length indicator .arch file is created
!           3. Where streams are to be reinitialised, a new file is
!              created with an appropriate name
!
! Code Description:
!   Language: FORTRAN 90 + common extensions.
!   This code is written to UMDP3 v6 programming standards.


SUBROUTINE ppctl_reinit (        &
     i_ao,icode,cmessage )

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE filenamelength_mod, ONLY:        &
    filenamelength
USE model_file, ONLY:                &
    model_file_open,                  &
    model_file_close
USE io_constants, ONLY:               &
    ioFileTypeMPIIO,                  &
    ioOpenReadWrite,                  &
    ioNoDelete
USE io_configuration_mod, ONLY: l_postp
USE um_ParVars
USE Control_Max_Sizes
USE dump_headers_mod, ONLY: a_fixhd, a_inthd, a_realhd, a_levdepc,     &
                            a_rowdepc, a_coldepc
USE submodel_mod, ONLY: atmos_im
USE iau_mod,     ONLY: l_iau
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE filename_generation_mod, ONLY: get_filename
USE file_manager, ONLY: &
    init_file_loop, um_file_type, assign_file_unit, release_file_unit
USE nlstcall_mod,         ONLY: model_analysis_mins, l_fastrun                 

USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, intf_len2_coldepc, intf_len2_levdepc,   &
    intf_len2_rowdepc, intf_lookupsa, len1_lookup,                     &
    len_dumphist, len_fixhd, max_intf_model_levels,                    &
    max_lbcrow_length, max_lbcrows, mpp_len1_lookup, n_intf_a,         &
    pp_len_inthd, pp_len_realhd

USE model_time_mod, ONLY: &
    secs_per_stepim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER          , INTENT(IN)    :: i_ao     ! Atmosphere/Ocean indicator
INTEGER          , INTENT(OUT)    :: icode    ! Return code from routine
CHARACTER(LEN=errormessagelength), INTENT(OUT)  :: cmessage 
                                             ! Return message if failure
                                                ! occurred
!
!----------------------------------------------------------------------
INTEGER :: name_offset
INTEGER :: arch_unit

CHARACTER(LEN=filenamelength) :: oldppfile ! Previous PPfile on unit
CHARACTER(LEN=filenamelength) :: archfile  ! Archiving file
TYPE(um_file_type), POINTER :: pp_file ! Pointer to file object

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PPCTL_REINIT'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Catch calls specifying wrong submodel id
IF (i_ao /= atmos_im) THEN
  CALL umPrint('PPCTL_REINIT: UM only supports Atmos model', &
      src='ppctl_reinit')
END IF

! If the run is using an analysis calculate the offset to use in any
! relative filenames
IF (l_fastrun .OR. l_iau) THEN
  name_offset = - (model_analysis_mins * 60)/secs_per_stepim(atmos_im)
ELSE
  name_offset = 0
END IF

! ----------------------------------------------------------------------
!  1.1 Loop over all valid FORTRAN units and select those to be
!      initialised this timestep as set by SETTSCTL
!
NULLIFY(pp_file)
pp_file => init_file_loop(handler="portio")
DO WHILE (ASSOCIATED(pp_file))

  IF (pp_file % meta % initialise) THEN
    !
    !  1.2   Save the name of the existing file before regenerating it
    !        
    !
    oldppfile = pp_file % filename

    ! Generate the archiving flag file
    IF (l_postp) THEN
      IF (INDEX(oldppfile, "Reserved unit") == 0) THEN
        archfile = TRIM(oldppfile(INDEX(oldppfile,'/',.TRUE.)+1:))//'.arch'
        CALL assign_file_unit(archfile, arch_unit, HANDLER="fortran")
        OPEN(UNIT=arch_unit, FILE=archfile)
        CLOSE(UNIT=arch_unit)
        CALL release_file_unit(arch_unit, HANDLER="fortran")
      END IF
    END IF

    !  1.3 Construct filename.
    !
    CALL get_filename(pp_file % meta % filename_base,             &
                      pp_file % filename,                         &
                      relative_offset_steps = name_offset,        &
                      reinit_steps = pp_file % meta % init_steps)
    !
    !  1.4 If a previous file existed it will be open (and hence, "managed")
    !      so make sure we close it before opening the new file
    !

    IF (pp_file % pp_meta % managed) THEN
      CALL model_file_close(pp_file % UNIT, oldppfile, &
          delete=ioNoDelete,error=icode)
    ELSE
      WRITE(umMessage,'(A,I4,A,A)')                                          &
          'PPCTL_REINIT: Warning tried to close unopen file on unit ',       &
          pp_file % UNIT,' name=',oldppfile
      CALL umPrint(umMessage,src='ppctl_reinit')

    END IF


    WRITE(umMessage,'(A,A,A,I3)')'PPCTL_REINIT: Opening new file ',          &
       TRIM(pp_file % filename),' on unit ', pp_file % UNIT
    CALL umPrint(umMessage,src='ppctl_reinit')

    CALL Model_File_Open(pp_file % UNIT, pp_file % filename,                 &
        read_write=ioOpenReadWrite, error=icode,allowRemap=.TRUE.,           &
        fileType=ioFileTypeMPIIO)

    IF (icode /= 0) THEN
      cmessage='PPCTL_REINIT   : Error opening new PPfile'
      GO TO 9999   !  Return
    END IF

    !
    !  1.6 Initialise the direct access LOOKUP headers if OUTPUT file
    !
    !  (a) PP files
    IF (pp_file % pp_meta % init_file_type == 'p' .OR. &
        pp_file % pp_meta % init_file_type == 'c') THEN

      IF (pp_file % meta % is_output_file) THEN
        WRITE(umMessage,'(A,I3)')                                     &
            'PPCTL_REINIT: Initialising new file on unit ',           &
            pp_file % UNIT
        CALL umPrint(umMessage,src='ppctl_reinit')
        ! DEPENDS ON: init_pp
        CALL init_pp(pp_file % UNIT, pp_file % pp_meta % init_file_type,     &
           len1_lookup, pp_file % pp_meta % reserved_headers,                &
           a_fixhd,a_inthd,a_realhd,a_levdepc,                               &
           a_rowdepc,a_coldepc,                                              &
           a_len_inthd,a_len_realhd,a_len1_levdepc,                          &
           a_len2_levdepc,a_len1_rowdepc,a_len2_rowdepc,                     &
           a_len1_coldepc,a_len2_coldepc,                                    &
           pp_len_inthd, pp_len_realhd,                                      &
           icode,cmessage)
      END IF
    END IF

    ! Reset counter of written fields and activate flag to show that this 
    ! file has been initialised
    pp_file % pp_meta % last_written_field = 0
    pp_file % meta % partially_written = .TRUE.

  END IF ! file requires reinitialisation

  ! Increment file loop pointer for next iteration
  pp_file => pp_file % next

END DO 

9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ppctl_reinit

