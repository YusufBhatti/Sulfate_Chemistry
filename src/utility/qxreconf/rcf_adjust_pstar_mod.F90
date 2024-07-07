! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Performs orographic adjustment of P*

MODULE Rcf_Adjust_Pstar_Mod

!  Subroutine Rcf_Adjust_Pstar    Adjusts P* when an ancillary
!                                 orography has been used with an
!                                 interpolated P*
!
! Description:
!   A similar method to that used in 4.5 is employed, with a
!   change to a height based vertical co-ordinate frame.
!
! Method:
!   Assuming a constant lapse rate, a surface temperature (Ts) is
!   calculated based on a temperature at a reference height (Tr)
!              Ts = Tr + lapse * (Zr - Z0_at_old_orography)
!
!   P* is then adjusted thus:-
!     P*' = P* [ Ts - lapse ( Zr - Z0_new_orography) ]  ** g/(R lapse)
!              [ ----------------------------------- ]
!              [                 Ts                  ]
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_ADJUST_PSTAR_MOD'

CONTAINS

SUBROUTINE Rcf_Adjust_Pstar( pstar, t, orog_in, orog_out, Input_Grid,  &
                             Output_Grid, r_theta_levels)

USE Rcf_Grid_Type_Mod, ONLY: &
    Grid_Type

USE Rcf_Field_Type_Mod, ONLY: &
    Field_Type

USE Rcf_Interp_Weights_Mod, ONLY: &
    h_int_active

USE Rcf_Set_Interp_Flags_Mod, ONLY: &
    interp_h_only,                   &
    interp_copy,                     &
    interp_no_op

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Interpolate_Mod, ONLY: &
    Rcf_Interpolate

USE Rcf_Field_Equals_Mod, ONLY: &
    Rcf_Field_Equals

USE planet_constants_mod, ONLY: g, planet_radius, r, lapse

USE atmos_model_working_constants_mod, ONLY: &
    upperheight

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( Grid_Type ), INTENT(IN)       :: Input_Grid
TYPE( Grid_Type ), INTENT(IN)       :: Output_Grid
TYPE( Field_Type ), INTENT(INOUT)   :: pstar
TYPE( Field_Type ), INTENT(INOUT)   :: orog_in
TYPE( Field_Type ), INTENT(IN)      :: orog_out
TYPE( Field_Type ), INTENT(IN)      :: t
REAL, INTENT(IN)                    :: r_theta_levels(                 &
                                           output_grid % loc_p_field,  &
                                       0 : output_grid % model_levels )

! Local Variables/Parameters
TYPE( Field_Type )    :: orog_interp         ! interpolated input
                                             ! orography
TYPE( Field_Type )    :: dummy               ! dummy field
REAL                  :: g_over_lapse_r
REAL                  :: Tr                  ! T at ref. level.
REAL                  :: Zr                  ! Height at ref. level.
REAL                  :: Z0i                 ! Interpolated orog.
REAL                  :: Z0o                 ! Output orog.
INTEGER               :: ref_lev             ! first theta level above
                                             ! upperheight
INTEGER               :: i                   ! Looper

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_ADJUST_PSTAR'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!--------------------------------------------------------------------
! Set up the orography fields to obtain interpolated orography
!--------------------------------------------------------------------
CALL Rcf_Field_Equals( orog_interp, orog_out )

IF (h_int_active) THEN
  orog_in % interp   = interp_h_only
ELSE
  orog_in % interp   = interp_copy
END IF

orog_interp % interp = interp_no_op

CALL Rcf_Alloc_Field( orog_interp )
CALL Rcf_Interpolate( orog_in, orog_interp, Input_Grid, Output_Grid, &
                      dummy, dummy )

! Similar calculation in calc_pmsl (could be merged into a common function?)
DO i = 1, output_grid % model_levels
  ! Upperheight is approximated to be where temperature is free from surface
  ! effects.  NB Model levels are terrain following near ground so height would be
  ! lower over land.
  ref_lev = i
  IF ( output_grid % z_top_of_model * output_grid % eta_theta_levels(i) &
       > upperheight ) EXIT
END DO

!--------------------------------------------------------------------
! Can now do the p* adjustment calculation as specified above
!--------------------------------------------------------------------
g_over_lapse_r = g/(lapse * r)
DO i = 1, pstar % level_size
  Tr  = t % DATA( i, ref_lev )
  Zr  = r_theta_levels( i, ref_lev ) - planet_radius
  Z0i = orog_interp % DATA( i, 1 )     ! interpolated orography
  Z0o = orog_out % DATA( i, 1 )        ! output orography

  pstar % DATA( i, 1 ) = pstar % DATA( i, 1 ) *                       &
                      ( ( Tr + lapse * (Zr - Z0o) ) /                 &
                        ( Tr + lapse * (Zr - Z0i) ) ) ** g_over_lapse_r
END DO


!--------------------------------------------------------------------
! Tidy up
!--------------------------------------------------------------------
CALL Rcf_Dealloc_Field( orog_interp )


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Adjust_Pstar

END MODULE Rcf_Adjust_Pstar_Mod
