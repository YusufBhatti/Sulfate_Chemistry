! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_derv_dry_rho_mod

IMPLICIT NONE
!  Subroutine Rcf_Derv_Dry_Rho_Mod

! Description:
!   Derive dry density from exner and density using the equation of state

! Method:
!   An implementation of the code previously (at VN8.4) located
!   in atm_step_4A.F90 l.1109-l.1118

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_DERV_DRY_RHO_MOD'

CONTAINS

SUBROUTINE rcf_derv_dry_rho( fields_in, field_count_in,     &
                              fields_out, field_count_out,  &
                              hdr_in,hdr_out,               &
                              dry_rho )

USE rcf_locate_mod, ONLY: &
    rcf_locate

USE rcf_alloc_field_mod, ONLY: &
    rcf_alloc_field,            &
    rcf_dealloc_field

USE rcf_read_field_mod, ONLY: &
    rcf_read_field

USE um_stashcode_mod, ONLY: &
    stashcode_mv,  stashcode_theta, stashcode_prog_sec, &
    stashcode_orog, stashcode_thetavd, stashcode_exner

USE rcf_field_type_mod, ONLY: &
    field_type

USE rcf_umhead_mod, ONLY: &
    um_header_type

USE umPrintMgr

USE um_parcore, ONLY: &
    mype

USE decomp_params, ONLY: &
    decomp_rcf_output,                           &
    decomp_rcf_input

USE rcf_interp_weights_mod, ONLY:               &
    h_int_active

USE rcf_v_int_ctl_mod, ONLY:                    &
    v_int_active

USE rcf_set_interp_flags_mod, ONLY:             &
    interp_h_only,                               &
    interp_v_only,                               &
    interp_copy,                                 &
    interp_all,                                  &
    interp_no_op

USE rcf_grid_type_mod, ONLY:                    &
    grid_type,                                   &
    input_grid, output_grid

USE rcf_interpolate_mod, ONLY:                  &
    rcf_interpolate

USE planet_constants_mod, ONLY: p_zero,kappa,r

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER       :: fields_out(:)
TYPE( field_type ), POINTER       :: fields_in(:)
TYPE( um_header_type), INTENT(IN) :: hdr_out
TYPE( um_header_type), INTENT(IN) :: hdr_in
TYPE( field_type ), INTENT(INOUT), TARGET :: dry_rho
INTEGER, INTENT(IN)               :: field_count_out
INTEGER, INTENT(IN)               :: field_count_in

! Internal variables
TYPE( field_type ), POINTER       :: exner
TYPE( field_type ), POINTER       :: thetavd
TYPE( field_type ), POINTER       :: orog_in, orog_out
TYPE( field_type )                :: orog_interp
TYPE( field_type )                :: dummy

CHARACTER (LEN=*), PARAMETER      :: RoutineName='RCF_DERV_DRY_RHO'

INTEGER                           :: interp_option ! interpolation to perform


INTEGER                           :: pos   ! position in array
INTEGER                           :: i,j,k,field ! loop index

INTEGER                           :: theta_pos
INTEGER                           :: thetavd_pos
INTEGER                           :: mv_pos
INTEGER                           :: dummy_pos

REAL :: intw_w2rho

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. printstatus >= prstatus_normal ) THEN
  WRITE (umMessage,FMT='(a)') 'Deriving dry_rho '
  CALL umPrint(umMessage, src='rcf_derv_dry_rho')
END IF


CALL rcf_locate(stashcode_prog_sec, stashcode_thetavd,                       &
                  fields_out, field_count_out, pos)
thetavd => fields_out(pos)
CALL rcf_alloc_field( thetavd )
CALL rcf_read_field( thetavd, hdr_out, decomp_rcf_output )

CALL rcf_locate(stashcode_prog_sec, stashcode_exner,                         &
                  fields_out, field_count_out, pos)
exner => fields_out(pos)
CALL rcf_alloc_field( exner )
CALL rcf_read_field( exner, hdr_out, decomp_rcf_output )



DO i = 1, dry_rho % levels

  intw_w2rho   = ( output_grid % eta_rho_levels(i)                           &
                  -output_grid % eta_theta_levels(i-1) ) /                   &
                 ( output_grid % eta_theta_levels(i)                         &
                  -output_grid % eta_theta_levels(i-1) )

  dry_rho % DATA(:,i) = p_zero * exner % DATA(:,i) ** ((1.0 - kappa)/kappa)  &
                        / (r *(       intw_w2rho *thetavd % DATA(:,i+1) +    &
                               (1.0 - intw_w2rho)*thetavd % DATA(:,i)))

END DO

CALL rcf_dealloc_field( exner )
CALL rcf_dealloc_field( thetavd )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_derv_dry_rho
END MODULE rcf_derv_dry_rho_mod
