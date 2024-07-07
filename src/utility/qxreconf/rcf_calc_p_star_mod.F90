! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Rcf_Calc_P_Star  - calculates P* from rho and exner

MODULE Rcf_Calc_P_Star_Mod

!  Subroutine Rcf_Calc_P_Star - calculates P* from rho and exner
!
! Description:
!    This subroutine calculates P* (surface pressure)
!
! Method:
!     Code derived New Dynamics code. Uses rho and exner and
!     calculates heights internally based on grid parameters and
!     orography.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_CALC_P_STAR_MOD'

CONTAINS

SUBROUTINE Rcf_Calc_P_Star( grid, fields, field_count, hdr, decomp,  &
                            orog, p_star )

USE um_stashcode_mod, ONLY: &
    stashcode_rho,             &
    stashcode_exner,           &
    stashcode_p,               &
    stashcode_prog_sec

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Exner_P_Convs_Mod, ONLY: &
    Rcf_Conv_Exner_P

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_Generate_Heights_Mod, ONLY: &
    height_gen_original,             &
    height_gen_smooth,               &
    height_gen_linear

USE Rcf_field_equals_mod, ONLY: &
    Rcf_Field_Equals

USE planet_constants_mod, ONLY: g, planet_radius

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)                :: field_count   ! number of fields
INTEGER, INTENT(IN)                :: decomp        ! decomposition
TYPE( grid_type ), INTENT(IN)      :: grid
TYPE( field_type ), POINTER        :: fields(:)
TYPE( field_type ), INTENT(IN)     :: orog          ! orography field
TYPE( UM_header_type ), INTENT(IN) :: hdr           ! dump header
TYPE( field_type ), INTENT(INOUT)  :: p_star        ! surface pressure

! Local variables
INTEGER                             :: i             ! looper
INTEGER                             :: pos           ! field position
REAL                                :: first_rho_level
TYPE( field_type )                  :: exner
TYPE( field_type ), POINTER         :: rho

INTEGER                             :: ErrorStatus
CHARACTER (LEN=*), PARAMETER        :: RoutineName = 'RCF_CALC_P_STAR'
CHARACTER (LEN=errormessagelength)  :: Cmessage

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-------------------------------------------------------------------
! Find rho in the fields array and read it in
!-------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_rho,                  &
                 fields, field_count, pos )
rho => fields(pos)
CALL Rcf_Alloc_Field( rho )
CALL Rcf_Read_Field( rho, hdr, decomp )

!------------------------------------------------------------------
! Find, read and convert exner to P - OR find and read P.
! Will be P in exner variable eventually however.
!------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_exner,                &
                 fields, field_count, pos, zero_ok_arg = .TRUE. )

IF (pos /= 0) THEN
  CALL Rcf_Field_Equals( exner, fields(pos) )
  CALL Rcf_Alloc_Field( exner )
  CALL Rcf_Read_Field( exner, hdr, decomp )
  CALL Rcf_Conv_Exner_P( exner )
ELSE        ! have P in dump
  CALL Rcf_Locate( stashcode_prog_sec, stashcode_p,                  &
                   fields, field_count, pos )
  CALL Rcf_Field_Equals( exner, fields(pos) )
  CALL Rcf_Alloc_Field( exner )
  CALL Rcf_Read_Field( exner, hdr, decomp )
END IF

!-------------------------------------------------------------------
! Calculate the first rho level height and thus P*
!-------------------------------------------------------------------
DO i = 1, orog % level_size
  IF ( grid % height_gen_method == height_gen_original .OR.            &
       grid % height_gen_method == height_gen_linear ) THEN

    first_rho_level = orog % DATA(i,1) + planet_radius +               &
                      grid % eta_rho_levels(1) * grid % z_top_of_model

  ELSE IF ( grid % height_gen_method == height_gen_smooth ) THEN

    first_rho_level = planet_radius + grid % eta_rho_levels(1) *       &
                      grid % z_top_of_model + orog % DATA(i,1) *       &
                      (1.0 - grid % eta_rho_levels(1) /                &
            grid % eta_rho_levels(grid % first_constant_r_rho_level))**2

  ELSE
    ErrorStatus = 10
    Cmessage = 'Height generation method unknown'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF


  p_star % DATA(i,1) = exner % DATA(i,1) + g * rho % DATA(i,1) *    &
                      (first_rho_level -                            &
                      (orog % DATA(i,1) + planet_radius)) /         &
                      ( first_rho_level * first_rho_level )
END DO

!-------------------------------------------------------------------
! Remove the no longer needed exner and rho fields
!-------------------------------------------------------------------
CALL Rcf_Dealloc_Field( exner )
CALL Rcf_Dealloc_Field( rho )
NULLIFY( rho )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Calc_P_Star
END MODULE Rcf_Calc_P_Star_Mod

