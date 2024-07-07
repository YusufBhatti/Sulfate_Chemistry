! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

SUBROUTINE ppctl_init_climate_means(     &
   i_ao,meanlev,icode,cmessage )

! Description:
!   Initialises a climate mean PP file.
!
! Method:
!   One file is initialised per call, based on the mean level and
!   the submodel.
!
!  External documentation: UM documentation papers:
!      C0 - The top-level control system;
!      C4 - Storage handling and diagnostic system
!      C5 - Control of means calculations
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 standards.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

USE yomhook,            ONLY: lhook, dr_hook
USE parkind1,           ONLY: jprb, jpim
USE filenamelength_mod, ONLY: filenamelength
USE model_file,         ONLY: model_file_open
USE io_constants, ONLY: ioOpenReadWrite
USE um_ParVars
USE Control_Max_Sizes
USE nlstcall_mod, ONLY: lcal360
USE nlstgen_mod, ONLY: &
    pp_len2_meanim, mean_filename_base, dumpfreqim, meanfreqim
USE dump_headers_mod, ONLY: a_fixhd, a_inthd, a_realhd, a_levdepc,     &
                            a_rowdepc, a_coldepc
USE submodel_mod, ONLY: atmos_im
USE filename_generation_mod, ONLY: get_filename
USE file_manager, ONLY: get_file_by_id, um_file_type
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, len1_lookup,                            &
    len_dumphist, len_fixhd, mpp_len1_lookup, pp_len_inthd, pp_len_realhd

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER           ,INTENT(IN)  :: i_ao      ! Atmosphere/Ocean indicator
INTEGER           ,INTENT(IN)  :: meanlev   ! Mean level indicator
INTEGER           ,INTENT(OUT) :: icode     ! Return code from routine
CHARACTER(LEN=errormessagelength) ,INTENT(OUT) :: cmessage  
                                            ! Return message on failure
TYPE(um_file_type), POINTER    :: um_file   ! Mean file object
!
! ----------------------------------------------------------------------

INTEGER :: i_mean            ! Mean period counter
INTEGER :: mean_steps_offset ! Offset for naming of files

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PPCTL_INIT_CLIMATE_MEANS'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

IF (i_ao /= atmos_im) THEN
  CALL umPrint('PPCTL_INIT_CLIMATE_MEANS: UM only supports Atmos model', &
      src='ppctl_init_climate_means')
  icode = 1
  GO TO 9999
END IF

!  1.1  Construct PPfile name from model information using defined
!       naming convention (OUTPUT files) and save this filename to
!       the file object (reserved earlier in readcntl)

! Note: Since these mean files are produced at the end of the period they 
!       represent, the filename generated must be offset by whatever the mean
!       period is (otherwise they would be labelled with the name of the next
!       meaning period instead)
IF (meanlev == 0) THEN
  mean_steps_offset = 0
ELSE
  mean_steps_offset = dumpfreqim(atmos_im)
  DO i_mean = 1, meanlev
    mean_steps_offset = mean_steps_offset*meanfreqim(i_mean, atmos_im)
  END DO
END IF

NULLIFY(um_file)
um_file => get_file_by_id("aomean", handler="portio")

CALL get_filename(TRIM(mean_filename_base(meanlev)),          &
                  um_file % filename,                         &
                  absolute_offset_steps = -mean_steps_offset, &
                  meanlev = meanlev,                          &
                  reinit_steps = um_file % meta % init_steps)

!  1.2 Open named file

WRITE(umMessage,'(A,A,A,I3)')'PPCTL_INIT_CLIMATE_MEANS: Opening new file ',  &
    um_file % filename,' on unit ',um_file % UNIT
CALL umPrint(umMessage,src='ppctl_init_climate_means')

CALL model_file_open(um_file % UNIT, um_file % filename, &
                     read_write=ioOpenReadWrite,error=icode)
IF (icode /= 0) THEN
  cmessage='PPCTL_INIT_CLIMATE_MEANS: Error opening new PPfile'
  GO TO 9999   !  Return
END IF

!  1.3 Initialise the direct access LOOKUP headers (OUTPUT file)
WRITE(umMessage,'(A,I3)')                                     &
   'PPCTL_INIT_CLIMATE_MEANS: Initialising new file on unit ',um_file % UNIT
CALL umPrint(umMessage,src='ppctl_init_climate_means')
! DEPENDS ON: init_pp
CALL init_pp(um_file % UNIT, um_file % pp_meta % init_file_type, &
   len1_lookup,pp_len2_meanim(meanlev,atmos_im),                 &
   a_fixhd,a_inthd,a_realhd,a_levdepc,                           &
   a_rowdepc,a_coldepc,                                          &
   a_len_inthd,a_len_realhd,a_len1_levdepc,                      &
   a_len2_levdepc,a_len1_rowdepc,a_len2_rowdepc,                 &
   a_len1_coldepc,a_len2_coldepc,                                &
   pp_len_inthd, pp_len_realhd,                                  &
   icode,cmessage)
IF (icode >  0) GO TO 9999
um_file % pp_meta % last_written_field = 0

9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ppctl_init_climate_means

