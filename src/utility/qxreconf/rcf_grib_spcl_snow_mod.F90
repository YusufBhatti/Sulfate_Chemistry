! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Convert snow amount from m to mm

MODULE rcf_grib_spcl_snow_mod
IMPLICIT NONE

! SUBROUTINE Rcf_Grib_Spcl_Snow

! Description:
!   This routine multiplies the values by 1000 to convert
!   snow amount from m to mm.
!   At present there is a double check to ensure the data
!   came from ECMWF as not all input data will be volumetric.
!
! Method:
!   Multiply by 1000.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.4 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GRIB_SPCL_SNOW_MOD'

CONTAINS
SUBROUTINE rcf_grib_spcl_snow(current,fielddata)

USE rcf_grib_block_params_mod, ONLY:                                          &
  grib_record,                                                                 &
  lenarraymax,                                                                 &
  p_orig_cntr

USE rcf_grib_lookups_mod, ONLY:                                               &
  grborigecmwf

USE umPrintMgr, ONLY:                                                         &
    umPrint,                                                                  &
    printstatus,                                                              &
    prstatus_normal,                                                          &
    prstatus_diag

USE um_parcore, ONLY:                                                         &
  mype

USE ereport_mod, ONLY:                                                        &
  ereport

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Subroutine arguments

!< Array  arguments with intent(InOut):>
TYPE (grib_record),POINTER          :: current
REAL, INTENT(INOUT)                 :: fielddata(lenarraymax)

! Local constants
CHARACTER (LEN=*), PARAMETER     :: RoutineName='RCF_GRIB_SPCL_SNOW'

! Local variables

CHARACTER (LEN=errormessagelength)    :: cmessage   ! used for EReport
INTEGER                          :: errorstatus   ! used for EReport
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!=======================================================================
!  Routine Code Start :
!=======================================================================

! double check it's an ECMWF field
IF (current % block_1(p_orig_cntr) == grborigecmwf) THEN

  IF ( printstatus >= prstatus_diag  ) THEN
    IF ( mype == 0 ) THEN
      CALL umPrint( "Converting snow amount from m of water to kg/m2", &
          src='rcf_grib_spcl_snow_mod')
    END IF
  END IF

  ! Multiply values by 1000
  fielddata(:) = 1000 * fielddata(:)

ELSE
  WRITE (cmessage,*) 'Data did not come from ECMWF'
  errorstatus = 10
  CALL ereport( routinename, errorstatus, cmessage )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE rcf_grib_spcl_snow
END MODULE rcf_grib_spcl_snow_mod
