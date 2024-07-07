! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_derv_exner_surf_mod

IMPLICIT NONE
!  Subroutine Rcf_Derv_Exner_Surf_Mod

! Description:
!   Derive surface Exner field using hydrostatic extrapolation

! Method:
!   This follows the implementation of VN8.4. There pstar was computed in
!    src/control/top_level/setcona_4A.F90 l.823-828, using Calc_P_star
!   This p_star value was then converted to Exner in atm_step_4A.F90
!    l.996-1003
!   Both is implemented below.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_DERV_EXNER_SURF_MOD'

CONTAINS

SUBROUTINE rcf_derv_exner_surf( fields_in, field_count_in,    &
                                fields_out, field_count_out,  &
                                hdr_in,hdr_out,               &
                                exner_surf )

USE rcf_locate_mod, ONLY: &
    rcf_locate

USE rcf_alloc_field_mod, ONLY: &
    rcf_alloc_field,            &
    rcf_dealloc_field

USE rcf_read_field_mod, ONLY: &
    rcf_read_field

USE um_stashcode_mod, ONLY: &
    stashcode_prog_sec, &
    stashcode_pstar, stashcode_orog,&
    stashcode_exner, stashcode_rho

USE rcf_field_type_mod, ONLY: &
    field_type

USE rcf_umhead_mod, ONLY: &
    um_header_type

USE umPrintMgr

USE decomp_params, ONLY: &
    decomp_rcf_output

USE rcf_grid_type_mod, ONLY:                    &
    grid_type,                                   &
    output_grid

USE planet_constants_mod, ONLY: recip_kappa, kappa, p_zero, g

USE rcf_exppx_mod, ONLY: &
    rcf_exppx

USE submodel_mod, ONLY: &
    atmos_im

USE rcf_generate_heights_mod, ONLY: rcf_generate_heights

USE um_parcore, ONLY: &
    mype

USE rcf_field_equals_mod, ONLY:                 &
    rcf_field_equals

USE rcf_set_interp_flags_mod, ONLY:             &
    interp_no_op

USE cppxref_mod, ONLY: ppx_atm_tall

USE stparam_mod, ONLY:     &
    st_levels_model_theta, &
    st_levels_model_rho

USE missing_data_mod, ONLY: imdi

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
TYPE( field_type ), INTENT(INOUT), TARGET :: exner_surf
INTEGER, INTENT(IN)               :: field_count_out
INTEGER, INTENT(IN)               :: field_count_in

! Internal variables
TYPE( field_type )                :: pstar
TYPE( field_type ), POINTER       :: exner
TYPE( field_type ), POINTER       :: rho_out
TYPE( field_type ), POINTER       :: orog_out

CHARACTER (LEN=*), PARAMETER      :: RoutineName='RCF_DERV_EXNER_SURF'

INTEGER                           :: pos      ! position in array

REAL, ALLOCATABLE :: xi3_at_rho(:,:)
REAL, ALLOCATABLE :: xi3_at_theta(:,:)

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. printstatus >= prstatus_normal ) THEN
  WRITE (umMessage,FMT='(a)') 'Deriving exner_surf. NOTE: ' &
       // 'constant g approximation!'
  CALL umPrint(umMessage,src='rcf_derv_exner_surf')
END IF

!----------------------------------------------------------------------
! Setup Orographies
!----------------------------------------------------------------------

  ! orogoraphy in output dump, must exist...
CALL rcf_locate( stashcode_prog_sec, stashcode_orog,             &
               fields_out,field_count_out,pos)
orog_out => fields_out(pos)
CALL rcf_alloc_field( orog_out )
CALL rcf_read_field( orog_out, hdr_out, decomp_rcf_output )

!--------------------------------------------------------------------
! Setup p_star_in field - either read or calculate from input dump
!--------------------------------------------------------------------

  ! compute it as previously done from setcona
pstar % levels          = 1
pstar % rows            = output_grid % loc_p_rows
pstar % row_len         = output_grid % loc_p_row_length
pstar % level_size      = output_grid % loc_p_field
pstar % glob_rows       = output_grid % glob_p_rows
pstar % glob_row_len    = output_grid % glob_p_row_length
pstar % glob_level_size = output_grid % glob_p_field
pstar % dump_pos        = imdi
pstar % interp          = interp_no_op
pstar % stashmaster     => rcf_exppx(atmos_im, 0, stashcode_pstar)
CALL rcf_alloc_field( pstar )

CALL rcf_locate( stashcode_prog_sec, stashcode_exner,          &
                 fields_out,field_count_out,pos)
exner => fields_out(pos)
CALL rcf_alloc_field( exner )
CALL rcf_read_field(  exner, hdr_out, decomp_rcf_output )


CALL rcf_locate( stashcode_prog_sec, stashcode_rho,            &
                fields_out,field_count_out,pos)
rho_out => fields_out(pos)
CALL rcf_alloc_field( rho_out )
CALL rcf_read_field(  rho_out, hdr_out, decomp_rcf_output )

ALLOCATE (xi3_at_rho  (output_grid % loc_p_row_length *        &
                       output_grid % loc_p_rows,               &
                       0:output_grid % model_levels+1))
ALLOCATE (xi3_at_theta(output_grid % loc_p_row_length *        &
                       output_grid % loc_p_rows,               &
                       0:output_grid % model_levels+1))

CALL rcf_generate_heights( output_grid,orog_out,               &
                           ppx_atm_tall, st_levels_model_theta,&
                           xi3_at_theta , pstar % level_size )

CALL rcf_generate_heights( output_grid,orog_out,               &
                           ppx_atm_tall, st_levels_model_rho,  &
                           xi3_at_rho , pstar % level_size )

pstar%DATA(:,1) = p_zero*exner%DATA(:,1)**recip_kappa          &
                            + g * rho_out%DATA(:,1)            &
                                * (xi3_at_rho  (:,1) -         &
                                   xi3_at_theta(:,0) ) /       &
                                  (xi3_at_rho  (:,1) *         &
                                   xi3_at_rho  (:,1) )
DEALLOCATE (xi3_at_theta)
DEALLOCATE (xi3_at_rho)
CALL rcf_dealloc_field( rho_out )
CALL rcf_dealloc_field( exner )

exner_surf%DATA(:,:) = (pstar%DATA(:,:)/p_zero)**kappa

CALL rcf_dealloc_field( pstar )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_derv_exner_surf
END MODULE rcf_derv_exner_surf_mod
