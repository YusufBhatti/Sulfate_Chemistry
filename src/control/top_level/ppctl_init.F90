! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Description:
!   Initialises the setup of fieldsfiles at the beginning of an NRUN or CRUN
!
! Method:
!           1. Opens all active input and output files on initial
!              call to routine.
!           2. Files to be processed on each call are controlled by
!              the file unit reinitialisation switch set on step 0 and 
!              at regular intervals thereafter.
!           3. All active PP are
!              intialised with user-specified number of LOOKUP headers
!              (fields).
!
! Code Description:
!   Language: FORTRAN 90 + common extensions.
!   This code is written to UMDP3 v6 programming standards.

SUBROUTINE ppctl_init(            &
   i_ao,icode,cmessage )

USE yomhook,                 ONLY: lhook, dr_hook
USE parkind1,                ONLY: jprb, jpim
USE model_file,              ONLY: model_file_open
USE io_constants,            ONLY: ioopenreadwrite, iofiletypempiio
USE submodel_mod,            ONLY: atmos_im
USE dump_headers_mod,        ONLY: a_fixhd, a_inthd, a_realhd,         &
                                   a_levdepc, a_rowdepc, a_coldepc
USE iau_mod,                 ONLY: l_iau
USE umPrintMgr,              ONLY: umPrint, umMessage
USE file_manager,            ONLY: init_file_loop, um_file_type
USE filename_generation_mod, ONLY: get_filename
USE errormessagelength_mod,  ONLY: errormessagelength
USE nlstcall_mod,            ONLY: model_analysis_mins, l_fastrun   
USE model_time_mod,          ONLY: secs_per_stepim

USE nlsizes_namelist_mod,    ONLY:                                     &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, intf_len2_coldepc, intf_len2_levdepc,   &
    intf_len2_rowdepc, intf_lookupsa, len1_lookup,                     &
    len_dumphist, len_fixhd, max_intf_model_levels,                    &
    max_lbcrow_length, max_lbcrows, mpp_len1_lookup, n_intf_a,         &
    pp_len_inthd, pp_len_realhd

IMPLICIT NONE


INTEGER, INTENT(IN)  :: i_ao     ! Atmosphere/Ocean indicator
INTEGER, INTENT(OUT) :: icode    ! Return code from routine
CHARACTER(LEN=errormessagelength), &
         INTENT(OUT) :: cmessage ! Return message if failure occurred

TYPE(um_file_type), POINTER :: pp_file ! Pointer to file objects
CHARACTER(LEN=*), PARAMETER :: RoutineName = "PPCTL_INIT"

INTEGER :: name_offset

! Dr-Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

IF (i_ao /= atmos_im) THEN
  CALL umPrint('PPCTL_INIT: UM only supports Atmos model',src=routinename)
END IF

! If the run is using an analysis calculate the offset to use in any
! relative filenames
IF (l_fastrun .OR. l_iau) THEN
  name_offset = - (model_analysis_mins * 60)/secs_per_stepim(atmos_im)
ELSE
  name_offset = 0
END IF

NULLIFY(pp_file)
pp_file => init_file_loop(handler="portio")
DO WHILE (ASSOCIATED(pp_file))

  IF (pp_file % meta % partially_written) THEN

    ! A CRUN will start with active files that exist on disk. This
    ! call reads the file's lookup and header to initialise the lookup
    ! and header arrays.

    ! DEPENDS ON: init_pp_crun
    CALL init_pp_crun(pp_file % UNIT, pp_file % filename, &
       len1_lookup,                                       &
       pp_file % pp_meta % reserved_headers,              &
       pp_file % pp_meta % init_file_type)

  ELSE

    ! Otherwise check whether or not the file needs to be initialised at the 
    ! start of the run
    IF (pp_file % meta % initialise) THEN

      IF (pp_file % meta % init_steps /= 0) THEN
        ! For re-initialising output streams a suitable filename must be 
        ! generated (if not the fixed filename is already set)
        CALL get_filename(pp_file % meta % filename_base,             &
                          pp_file % filename,                         &
                          relative_offset_steps = name_offset,        &
                          reinit_steps = pp_file % meta % init_steps)
      END IF

      ! Open the file and initialise the PP headers
      WRITE(umMessage,'(A,A,A,I3)')'PPCTL_INIT: Opening new file ',   &
          TRIM(pp_file % filename),' on unit ', pp_file % UNIT
      CALL umPrint(umMessage,src=routinename)

      CALL model_file_open(pp_file % UNIT, pp_file % filename,        &
          read_write=ioOpenReadWrite, fileType=ioFileTypeMPIIO)

      ! DEPENDS ON: init_pp
      CALL init_pp(pp_file % UNIT, pp_file % pp_meta % init_file_type,    &
          len1_lookup, pp_file % pp_meta % reserved_headers,              &
          a_fixhd,a_inthd,a_realhd,a_levdepc,                             &
          a_rowdepc,a_coldepc,                                            &
          a_len_inthd,a_len_realhd,a_len1_levdepc,                        &
          a_len2_levdepc,a_len1_rowdepc,a_len2_rowdepc,                   &
          a_len1_coldepc,a_len2_coldepc,pp_len_inthd,                     &
          pp_len_realhd,icode,cmessage)
      IF (icode >  0) GO TO 9999

      ! Set the field counter - this is a new file so should be zero, and
      ! the flag which indicates the file is active
      pp_file % pp_meta % last_written_field = 0
      pp_file % meta % partially_written = .TRUE.

    END IF ! File needs to be initialised
  END IF ! File is partially written
    
  ! Increment the pointer to the next file
  pp_file => pp_file % next
END DO

9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ppctl_init

