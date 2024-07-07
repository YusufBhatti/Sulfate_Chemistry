! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Determine STASH input length per vertical level for prog var

MODULE Rcf_Address_Length_Mod

!  Subroutine Rcf_Address_Length - determines field size
!
! Description:
!    Calculates size of field for output dump addressing (single level).
!
! Method:
!    Calculates field sizes based on Grid_Type from stashmaster
!    Generally single level, but old LBC needs level info.
!
!    Based on UM 4.5/5.0 code.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_ADDRESS_LENGTH_MOD'

CONTAINS

SUBROUTINE Rcf_Address_Length( igpl, halo_type, LEN )

USE rimtypes
USE lbc_mod

USE Rcf_Lsm_Mod, ONLY: &
    glob_land_out

USE ozone_inputs_mod, ONLY:  &
    zon_av_ozone

USE Rcf_Grid_Type_Mod, ONLY: &
    Output_Grid

USE UM_ParVars

USE nlsizes_namelist_mod, ONLY: &
    tr_vars

USE Rcf_Level_Code_Mod, ONLY: &
    Rcf_Level_Code

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Subroutine arguments:
INTEGER, INTENT(IN)   :: igpl      ! Grid code
INTEGER, INTENT(IN)   :: halo_type ! code from stashmaster
INTEGER, INTENT(OUT)  :: LEN       ! Length

! Local scalars:
INTEGER               :: ip   !  pressure levels
INTEGER               :: it   !  tracer levels
INTEGER               :: ix1  !  leftmost point
INTEGER               :: ix2  !  rightmost point
INTEGER               :: iy1  !  lower point
INTEGER               :: iy2  ! upper point

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_ADDRESS_LENGTH'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!- End of Header ---------------------------------------------------

! Determine row/column nos. for global domain on output grid
! DEPENDS ON: lltorc
CALL lltorc(igpl,90,-90,0,360,iy1,iy2,ix1,ix2, Output_Grid)


SELECT CASE (igpl)
CASE ( 21 )
  LEN= glob_land_out     ! Land compressed

CASE ( 22 )
  IF (zon_av_ozone) THEN   !  Zonal
    LEN = iy2-iy1+1
  ELSE                   !  Full fields
    LEN = (ix2-ix1+1)*(iy2-iy1+1)
  END IF

CASE (23)
  LEN = Output_Grid % glob_r_field

CASE ( 25 )              ! Old LBC
  CALL Rcf_Level_Code( 2, ip, Output_Grid)   ! pressure levels
  CALL Rcf_Level_Code(11, it, Output_Grid)   ! tracer levels
  LEN=( Output_Grid % glob_p_rows +                                &
        Output_Grid % glob_p_row_length                            &
        -2*rimwidtha(rima_type_norm))*2*rimwidtha(rima_type_norm)* &
      ( 1+5*ip+tr_vars*it)-2*ip*4*rimwidtha(rima_type_norm)


CASE DEFAULT
  LEN =(ix2-ix1+1)*(iy2-iy1+1)

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Address_Length

END MODULE Rcf_Address_Length_Mod
