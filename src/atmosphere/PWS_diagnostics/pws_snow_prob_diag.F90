! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE pws_snow_prob_diag_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_SNOW_PROB_DIAG_MOD'

CONTAINS
!
!  Routine to calculate snow probability diagnostic
!  according to Boyden Method based on 850,1000mb thicknesses

SUBROUTINE pws_snow_prob_diag(STASHnumber,     &  ! in
                              icode,cmessage)     ! inout

! Description:
!
! Method: Boyden CL, 1964: A comparison of snow predictors. Met.Mag,93,353-365
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP 003

! needed by typ_atm_fields.h:
USE atm_fields_bounds_mod
USE nlsizes_namelist_mod, ONLY: row_length, rows

USE um_stashcode_mod, ONLY: stashcode_pws_sec,                    &
    stashcode_pws_snow_prob
USE pws_diags_mod,    ONLY:                                       &
    pws_snow_prob, pws_geopht_1000, pws_geopht_850

USE missing_data_mod, ONLY: rmdi
USE um_parvars

USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN) :: STASHnumber        ! STASH request
INTEGER :: icode ! Return code : 0 Normal exit, >0 Error exit
CHARACTER(LEN=errormessagelength) :: cmessage  ! Error message if ICODE > 0

! Local Constants:
INTEGER, PARAMETER :: STASH_snow_prob = 20028
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'PWS_SNOW_PROB_DIAG'

! dr_hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header --------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate Snow Probability by Boyden Method
IF (STASHnumber == STASH_snow_prob) THEN

pws_snow_prob = 5220.0 +3.86666*pws_geopht_1000 -4.0*pws_geopht_850

WHERE ( pws_geopht_1000 >=  1.0e8 )
  pws_snow_prob = 0.0
END WHERE

! Limit to percentage probability
WHERE (  pws_snow_prob <=    0.0 )
  pws_snow_prob = 0.0
END WHERE
WHERE (  pws_snow_prob >=  100.0 )
  pws_snow_prob = 100.0
END WHERE

ELSE

icode = -STASHnumber        ! Warning
cmessage='Incorrect STASH request'
CALL ereport(RoutineName,icode,cmessage)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE pws_snow_prob_diag

END MODULE pws_snow_prob_diag_mod
