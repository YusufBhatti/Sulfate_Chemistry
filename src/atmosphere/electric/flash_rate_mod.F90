! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE flash_rate_mod

! Purpose: Calls subroutines to calculate thunderstorm flash rate
!          dependent on the method requested by the user.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Electric

USE atm_fields_bounds_mod,  ONLY: wdims_s, tdims
USE electric_inputs_mod,    ONLY: electric_method, em_gwp, em_mccaul
USE electric_constants_mod, ONLY: minus5

USE fr_gwp_mod,             ONLY: fr_gwp
USE fr_mccaul_mod,          ONLY: fr_mccaul
USE calc_wp_below_t_mod,    ONLY: calc_wp_below_t

! Dr Hook modules
USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FLASH_RATE_MOD'

CONTAINS

SUBROUTINE flash_rate( storm_field, qgraup, t, w, rhodz, gwp, tiwp, flash,     &
                       fr1_mc, fr2_mc )

IMPLICIT NONE

REAL, INTENT(IN) :: qgraup(     tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end,                   &
                                            1 : tdims%k_end      )

REAL, INTENT(IN) :: t(          tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end,                   &
                                            1 : tdims%k_end      )

REAL, INTENT(IN) :: w(        wdims_s%i_start : wdims_s%i_end,                 &
                              wdims_s%j_start : wdims_s%j_end,                 &
                              wdims_s%k_start : wdims_s%k_end    )

REAL, INTENT(IN) :: rhodz(      tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end,                   &
                                            1 : tdims%k_end      )

REAL, INTENT(IN) :: gwp(        tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end      )

REAL, INTENT(IN) :: tiwp(       tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end      )

REAL, INTENT(INOUT) :: flash(   tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end      )

REAL, INTENT(INOUT) :: fr1_mc(  tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end      )

REAL, INTENT(INOUT) :: fr2_mc(  tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end      )

LOGICAL, INTENT(IN) :: storm_field( tdims%i_start : tdims%i_end,               &
                                    tdims%j_start : tdims%j_end      )

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FLASH_RATE'

!-----------------------------------------------------------
! Local variables
!-----------------------------------------------------------
REAL :: gwp_m5(  tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end      )


!==================================================================
! Start the subroutine
!==================================================================

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!====================================================================
! Determine flash rate according to electric scheme
!====================================================================

SELECT CASE (electric_method)

CASE (em_gwp) ! GWP relation

  ! First generate gwp below -5 C threshold
  CALL calc_wp_below_t(qgraup, rhodz, t, minus5, gwp_m5 )

  ! Next call the flash rate
  CALL fr_gwp(storm_field, gwp_m5, flash)

CASE (em_mccaul) ! McCaul et al (2009), Weather and Forecasting

  CALL fr_mccaul( storm_field, t, w, qgraup, gwp, tiwp, flash, fr1_mc, fr2_mc)

CASE DEFAULT ! GWP relation
             ! (this is your failsafe option in case electric method
             !  is not defined)

  ! First generate gwp below -5 C threshold
  CALL calc_wp_below_t(qgraup, rhodz, t, minus5, gwp_m5 )

  ! Next call the flash rate
  CALL fr_gwp(storm_field, gwp_m5, flash)

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE flash_rate

END MODULE flash_rate_mod
