! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  finalise the reconfiguration - ties up loose ends

MODULE Rcf_Finalise_Mod

!   Subroutine Rcf_Finalise - final tidyup tasks
!
! Description:
! This module performs tidy up tasks including closing files and
! calling Timer for output.
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_FINALISE_MOD'

CONTAINS

SUBROUTINE Rcf_Finalise( hdr_in, hdr_out )

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE file_manager, ONLY: &
    release_file_unit

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_Grid_Type_Mod, ONLY: &
    Input_Grid

USE rcf_nlist_recon_technical_mod, ONLY: &
    ainitial,                            &
    grib_input_dump,                     &
    input_dump_type

USE nlstcall_mod, ONLY:   &
    LTimer

USE io, ONLY: file_close
USE io_constants, ONLY: ioNameProvided, ioNameInEnv, ioNoDelete
USE model_file, ONLY: model_file_close
USE nlcfiles_namelist_mod, ONLY: astart

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (um_header_type), INTENT(INOUT)  :: hdr_in
TYPE (um_header_type), INTENT(INOUT)  :: hdr_out

! Local variables
INTEGER                               :: err
INTEGER                               :: ErrorStatus
CHARACTER (LEN=*), PARAMETER          :: RoutineName='RCF_FINALISE'
CHARACTER (LEN=errormessagelength)    :: Cmessage
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-------------------------------------------------------------------
! Close open files
!-------------------------------------------------------------------
err = 0
! Output dump

CALL Model_File_Close( hdr_out % UnitNum, astart, &
                       name_in_environ=ioNameProvided, delete=ioNoDelete, &
                       error=err)

IF ( err /= 0 ) THEN
  Cmessage    = 'Failed to Close Output Dump'
  ErrorStatus = 10
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

CALL release_file_unit ( hdr_out % UnitNum, handler="portio" )

! Input dump
IF (Input_Grid % Rotated .AND.                                                &
     input_dump_type /= grib_input_dump ) THEN

  CALL File_Close( hdr_in % UnitNum, 'RECONTMP', 8, &
                   ioNameInEnv, ioNoDelete, err)
ELSE

  CALL File_Close( hdr_in % UnitNum, ainitial, &
                   name_in_environ=ioNameProvided, delete=ioNoDelete, error=err)
END IF

IF ( err /= 0 ) THEN
  Cmessage    = 'Failed to Close Input Dump'
  ErrorStatus = 20
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

CALL release_file_unit ( hdr_in % UnitNum, handler="portio" )

!-------------------------------------------------------------------
! Produce Timer output
!-------------------------------------------------------------------
IF (LTimer) THEN
  CALL Timer("Reconfigure",2)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Finalise
END MODULE Rcf_Finalise_Mod
