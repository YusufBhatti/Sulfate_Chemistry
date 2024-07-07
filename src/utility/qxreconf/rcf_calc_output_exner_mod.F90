! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Rcf_Calc_Output_Exner_Mod - calculate output dump exner pressure

MODULE Rcf_Calc_Output_Exner_Mod
IMPLICIT NONE

!  Subroutine Rcf_Calc_Output_Exner
!
! Description:
!    This module calculates exner pressure for the output dump
!
! Method:
!    This code is derived from that used in the New Dynamics
!    reconfiguration program. First P* is calculated on the input
!    grid and then interpolated horizontally onto the output grid.
!    This is then used to calculate the 1st level of exner which
!    in turn is used to calculate the higher levels.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_CALC_OUTPUT_EXNER_MOD'

CONTAINS

SUBROUTINE Rcf_Calc_Output_Exner( fields_in, field_count_in, orog_in, &
                                  hdr_in, orog_out, t, q, Exner,      &
                                  orog_source, r_rho_levels,          &
                                  r_theta_levels )
USE UM_ParCore, ONLY: &
    mype

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,                &
    PrStatus_Normal

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_UMhead_Mod, ONLY: &
    UM_header_type

USE Rcf_Grid_Type_Mod, ONLY: &
    Input_Grid,          &
    Output_Grid

USE um_stashcode_mod, ONLY: &
    stashcode_pstar,           &
    stashcode_prog_sec

USE Rcf_Exppx_Mod, ONLY: &
    Rcf_Exppx

USE Rcf_Set_Interp_Flags_Mod, ONLY: &
    interp_h_only,                   &
    interp_copy,                     &
    interp_no_op

USE Submodel_Mod, ONLY:             &
    atmos_im

USE Rcf_Calc_P_Star_Mod, ONLY: &
    Rcf_Calc_P_Star

USE decomp_params, ONLY: &
    Decomp_rcf_input

USE Rcf_Interp_Weights_Mod, ONLY: &
    h_int_active

USE Rcf_Interpolate_Mod, ONLY: &
    Rcf_Interpolate

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Adjust_Pstar_Mod, ONLY: &
    Rcf_Adjust_Pstar

USE items_nml_mod, ONLY: &
    Ancillary_File

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE rcf_field_equals_mod, ONLY: &
    rcf_field_equals

USE planet_constants_mod, ONLY: g, cp, kappa, pref, repsilon

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)                :: field_count_in
TYPE( field_type ), POINTER        :: fields_in(:)
TYPE( field_type ), INTENT(INOUT)  :: orog_in
TYPE( field_type ), INTENT(IN)     :: orog_out
TYPE( field_type ), INTENT(IN)     :: t
TYPE( field_type ), INTENT(IN)     :: q
TYPE( field_type ), INTENT(INOUT)  :: Exner
TYPE( UM_header_type ), INTENT(IN) :: hdr_in
REAL, INTENT(IN)                   :: r_rho_levels(               &
                                      output_grid % loc_p_field,  &
                                      0 : output_grid % model_levels+1)
REAL, INTENT(IN)                   :: r_theta_levels(             &
                                      output_grid % loc_p_field,  &
                                      0 : output_grid % model_levels+1)
INTEGER, INTENT(IN)                :: orog_source

! Local variables
INTEGER                            :: i          ! Looper
INTEGER                            :: k          ! Looper
INTEGER                            :: pos        ! position in fields_in
INTEGER                            :: k_off      ! Offset due to 0th level in T
REAL                               :: temp       ! temporary value
REAL                               :: weight1    ! averaging weight
REAL                               :: weight2    ! averaging weight
REAL                               :: deltaz     ! vertical spacing
TYPE( field_type )                 :: pstar_in
TYPE( field_type )                 :: pstar_out
TYPE( field_type )                 :: dummy

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_CALC_OUTPUT_EXNER'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
  WRITE(umMessage,*) 'Calculating Exner'
  CALL umPrint(umMessage,src='rcf_calc_output_exner_mod')
END IF

!--------------------------------------------------------------------
! Setup p_star_in field - either read or calculate from input dump
!--------------------------------------------------------------------
! Check input dump for p_star
CALL Rcf_Locate( stashcode_prog_sec, stashcode_Pstar,                &
                 fields_in, field_count_in, pos, zero_ok_arg = .TRUE. )

! If present, read in, else calculate.
IF (pos /= 0) THEN
  Pstar_in  = fields_in(pos)
  CALL Rcf_Alloc_Field( Pstar_in )
  CALL Rcf_Read_Field( Pstar_in, Hdr_In, decomp_rcf_input )

ELSE

  CALL rcf_field_equals(pstar_in, input_grid)
  pstar_in % stashmaster     => Rcf_Exppx(atmos_im, 0, stashcode_pstar)

  CALL Rcf_Alloc_Field( pstar_in )
  CALL Rcf_Calc_P_Star( Input_Grid, fields_in, field_count_in, hdr_in, &
                        decomp_rcf_input, orog_in, pstar_in )

END IF

! Only need to do horizontal interpolation of P* if h_int_active
IF (h_int_active) THEN
  pstar_in % interp          = interp_h_only
ELSE
  pstar_in % interp          = interp_copy
END IF

!--------------------------------------------------------------------
! Setup p_star_out field - interpolated from p_star_in
!--------------------------------------------------------------------
CALL rcf_field_equals(pstar_out, output_grid)
pstar_out % stashmaster     => Rcf_Exppx( atmos_im, 0, stashcode_pstar )

CALL Rcf_Alloc_Field( pstar_out )
CALL Rcf_Interpolate( pstar_in, pstar_out, Input_Grid, Output_Grid, &
                      dummy, dummy)

!--------------------------------------------------------------------
! If required, adjust the p* field for orographic effects
!--------------------------------------------------------------------
IF ( orog_source == Ancillary_File ) THEN
  CALL Rcf_Adjust_Pstar( pstar_out, t, orog_in, orog_out, Input_Grid, &
                         Output_Grid, r_theta_levels )
END IF

!--------------------------------------------------------------------
! calculate exner
! equation is
!
! exner(k) (Cp T_v(k-1) + g weight1 (r(k)-r(k-1)) =
! exner(k-1) (Cp T_v(k-1) - (1-weight1) (r(k)-r(k-1))
!
! where the weights represent the interpolation coefficients
! required to calculate exner at a T point.
! Note: In the new integration scheme Theta below level 1 is assumed
!       to be equal to theta at level 1. The equation is thus
!
! exner(1) (Cp T_v(1) + g (r(k)-r(k-1)) =
! exner* (Cp T_v(1) )
!--------------------------------------------------------------------

IF (t % bottom_level == 0 .AND. q % bottom_level == 0) THEN
  k_off = 1
ELSE
  k_off = 0
END IF

!--------------------------------------------------------------------
! Rho level 1
!--------------------------------------------------------------------

k = 1
DO i = 1, pstar_out % level_size

  temp = t % DATA(i,k+k_off) * (1.0 + (1.0/repsilon -1.0) * q % DATA(i,k+k_off))
  exner % DATA(i,k) = (cp * temp *                                &
                      (pstar_out % DATA(i,1) / pref)**kappa ) /   &
                      (cp * temp + g * (r_rho_levels(i,k) -       &
                        r_theta_levels(i,k-1) ) )

END DO

!--------------------------------------------------------------------
! Rho levels 2 to model_levels + 1
!--------------------------------------------------------------------
DO k = 2, exner % levels
  IF ( k <= Output_Grid % model_levels + 1 ) THEN
    DO i = 1, exner % level_size
      IF (k > Output_Grid % model_levels) THEN

        ! extra pressure levels is same height above top theta as
        ! pressure below is below top theta - hence weights are 0.5
        weight1 = 0.5
        deltaz  = 2.0 * (r_theta_levels(i,k-1) - r_rho_levels(i,k-1))

      ELSE

        deltaz  = r_rho_levels(i,k) - r_rho_levels(i,k-1)
        weight1 = (r_rho_levels(i,k) - r_theta_levels(i,k-1) ) / deltaz

      END IF ! (K > model_levels)

      weight2 = 1.0 - weight1

      temp = t % DATA(i,k-1+k_off) * (1.0 + (1.0/repsilon -1.0)                 &
                             * q % DATA(i,k-1+k_off))

      exner % DATA(i,k) = (cp * temp - g * weight1 * deltaz) *         &
                          exner % DATA(i,k-1) /                        &
                          (cp * temp + g * weight2 * deltaz)

    END DO
  ELSE    ! k is > Output_Grid % model_levels + 1
    DO i = 1, exner % level_size

      ! extra pressure levels is same height above top theta as
      ! pressure below is below top theta - hence weights are 0.5
      weight1 = 0.5
      deltaz = 2.0 * (r_theta_levels(i,k-1) - r_rho_levels(i,k-1))

      weight2 = 1.0 - weight1

      temp = t % DATA(i,k-1+k_off)
      exner % DATA(i,k) = (cp * temp - g * weight1 * deltaz) *       &
                          exner % DATA(i,k-1) /                      &
                          (cp * temp + g * weight2 * deltaz)

    END DO
  END IF ! wet/dry layer
END DO ! k


!--------------------------------------------------------------------
! Clear up memory for the pstar fields
!--------------------------------------------------------------------
CALL Rcf_Dealloc_Field( pstar_in )
CALL Rcf_Dealloc_Field( pstar_out )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Calc_Output_Exner
END MODULE Rcf_Calc_Output_Exner_Mod
