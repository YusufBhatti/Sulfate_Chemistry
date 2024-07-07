! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Function Interface
INTEGER FUNCTION get_fld_type (grid_type_code)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
USE cppxref_mod, ONLY:  &
    ppx_atm_tall,       &
    ppx_atm_tland,      &
    ppx_atm_tsea,       &
    ppx_atm_tzonal,     &
    ppx_atm_tmerid,     &
    ppx_atm_compressed, &   
    ppx_atm_ozone,      &
    ppx_atm_lbc_theta,  &
    ppx_atm_lbc_orog,   &
    ppx_atm_cuall,      &
    ppx_atm_lbc_u,      &
    ppx_atm_cvall,      &
    ppx_atm_lbc_v,      &
    ppx_atm_uall,       &
    ppx_atm_uland,      & 
    ppx_atm_usea,       &
    ppx_atm_uzonal,     &
    ppx_atm_umerid,     &
    ppx_atm_scalar,     &
    ppx_atm_river

IMPLICIT NONE


! Description:
! Takes a STASH grid type code, and returns which type of
! grid this is - mass or wind grid.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

INTEGER, INTENT(IN)  ::  grid_type_code     ! IN : STASH grid type code

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GET_FLD_TYPE'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF ( (grid_type_code  ==  ppx_atm_tall) .OR.                      &
     (grid_type_code  ==  ppx_atm_tland) .OR.                     &
     (grid_type_code  ==  ppx_atm_tsea) .OR.                      &
     (grid_type_code  ==  ppx_atm_tzonal) .OR.                    &
     (grid_type_code  ==  ppx_atm_tmerid) .OR.                    &
     (grid_type_code  ==  ppx_atm_compressed) .OR.                &
     (grid_type_code  ==  ppx_atm_ozone) .OR.                     &
     (grid_type_code  ==  ppx_atm_lbc_theta)  .OR.                &
     (grid_type_code  ==  ppx_atm_scalar)  .OR.                   &
     (grid_type_code  ==  ppx_atm_lbc_orog) ) THEN
  get_fld_type=fld_type_p
ELSE IF                                                           &
   ( (grid_type_code  ==  ppx_atm_cuall) .OR.                     &
     (grid_type_code  ==  ppx_atm_lbc_u) ) THEN
  get_fld_type=fld_type_u
ELSE IF                                                           &
   ( (grid_type_code  ==  ppx_atm_cvall) .OR.                     &
     (grid_type_code  ==  ppx_atm_lbc_v)  .OR.                    &
     (grid_type_code  ==  ppx_atm_uall ) ) THEN
  get_fld_type=fld_type_v
ELSE IF                                                           &
   ( (grid_type_code  ==  ppx_atm_uall) .OR.                      &
     (grid_type_code  ==  ppx_atm_uland) .OR.                     &
     (grid_type_code  ==  ppx_atm_usea) .OR.                      &
     (grid_type_code  ==  ppx_atm_uzonal) .OR.                    &
     (grid_type_code  ==  ppx_atm_umerid) ) THEN
  get_fld_type=fld_type_v  ! This is actually the B U/V (velocity)
                           ! grid, but it has the same sizes as
                           ! the C V grid.
ELSE IF                                                           &
   ( (grid_type_code  ==  ppx_atm_river) ) THEN
  get_fld_type=fld_type_r  ! River routing grid

ELSE
  get_fld_type=fld_type_unknown
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END FUNCTION get_fld_type

