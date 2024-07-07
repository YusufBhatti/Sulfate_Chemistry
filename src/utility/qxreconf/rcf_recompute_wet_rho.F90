#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE rcf_recompute_wet_rho_mod

IMPLICIT NONE

!  Subroutine Rcf_Recompute_Wet_Rho_Mod

! Description:
!   Derive wet density r**2 from dry density and mixing ratios

! Method:
!   An implementation of the code previously (at VN8.4) located
!   in atm_step_4A.F90 l.1138-l.1162

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_RECOMPUTE_WET_RHO_MOD'

CONTAINS

SUBROUTINE rcf_recompute_wet_rho( fields_out, field_count_out,  &
                                  hdr_out )

USE rcf_locate_mod, ONLY: &
    rcf_locate

USE rcf_alloc_field_mod, ONLY: &
    rcf_alloc_field,            &
    rcf_dealloc_field

USE rcf_read_field_mod, ONLY: &
    rcf_read_field

USE um_stashcode_mod, ONLY: &
    stashcode_mv,  stashcode_mcl,  stashcode_mcf,   &
    stashcode_mr,  stashcode_mgr,  stashcode_mcf2,  &
    stashcode_dry_rho, stashcode_prog_sec, &
    stashcode_orog, stashcode_rho

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

USE rcf_grid_type_mod, ONLY:                    &
    grid_type,                                   &
    input_grid, output_grid

USE planet_constants_mod, ONLY: p_zero,kappa,r

USE rcf_generate_heights_mod, ONLY: rcf_generate_heights

USE rcf_write_field_mod, ONLY: &
    rcf_write_field

USE cppxref_mod, ONLY: ppx_atm_tall,    &
                       ppx_atm_cvall,   &
                       ppx_atm_cuall

USE stparam_mod, ONLY:     &
    st_levels_model_rho

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER       :: fields_out(:)
TYPE( um_header_type), INTENT(IN) :: hdr_out
INTEGER, INTENT(IN)               :: field_count_out


! Internal variables
TYPE( field_type ), POINTER       :: mv
TYPE( field_type ), POINTER       :: mcl
TYPE( field_type ), POINTER       :: mcf
TYPE( field_type ), POINTER       :: mr
TYPE( field_type ), POINTER       :: mgr
TYPE( field_type ), POINTER       :: mcf2
TYPE( field_type ), POINTER       :: dry_rho
TYPE( field_type ), POINTER       :: rho
TYPE( field_type ), POINTER       :: orog_in, orog_out
TYPE( field_type )                :: orog_interp
TYPE( field_type )                :: dummy

CHARACTER (LEN=*), PARAMETER      :: RoutineName='RCF_RECOMPUTE_WET_RHO'

INTEGER                           :: pos   ! position in array
INTEGER                           :: i,j,k,field ! loop index

REAL, ALLOCATABLE :: xi3_at_rho(:,:)

REAL :: intw_w2rho, mixing_r

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. printstatus >= prstatus_normal ) THEN
  WRITE (umMessage,FMT='(a)') 're-deriving wet rho'
  CALL umPrint(umMessage,src='rcf_recompute_wet_rho')
END IF

! orogoraphy in output dump, must exist...
CALL rcf_locate( stashcode_prog_sec, stashcode_orog,                    &
               fields_out,field_count_out,pos)
orog_out => fields_out(pos)
CALL rcf_alloc_field( orog_out )
CALL rcf_read_field( orog_out, hdr_out, decomp_rcf_output )


CALL rcf_locate(stashcode_prog_sec, stashcode_mv,                       &
                  fields_out, field_count_out, pos)
mv => fields_out(pos)
CALL rcf_alloc_field( mv )
CALL rcf_read_field(  mv, hdr_out, decomp_rcf_output )


CALL rcf_locate(stashcode_prog_sec, stashcode_mcl,                      &
                  fields_out, field_count_out, pos)
mcl => fields_out(pos)
CALL rcf_alloc_field( mcl )
CALL rcf_read_field(  mcl, hdr_out, decomp_rcf_output )


CALL rcf_locate(stashcode_prog_sec, stashcode_mcf,                      &
                  fields_out, field_count_out, pos)
mcf => fields_out(pos)
CALL rcf_alloc_field( mcf )
CALL rcf_read_field(  mcf, hdr_out, decomp_rcf_output )


CALL rcf_locate(stashcode_prog_sec, stashcode_mr,                       &
                  fields_out, field_count_out, pos)
mr => fields_out(pos)
CALL rcf_alloc_field( mr )
CALL rcf_read_field(  mr, hdr_out, decomp_rcf_output )


CALL rcf_locate(stashcode_prog_sec, stashcode_mgr,                      &
                  fields_out, field_count_out, pos)
mgr => fields_out(pos)
CALL rcf_alloc_field( mgr )
CALL rcf_read_field(  mgr, hdr_out, decomp_rcf_output )


CALL rcf_locate(stashcode_prog_sec, stashcode_mcf2,                     &
                  fields_out, field_count_out, pos)
mcf2 => fields_out(pos)
CALL rcf_alloc_field( mcf2 )
CALL rcf_read_field(  mcf2, hdr_out, decomp_rcf_output )


CALL rcf_locate(stashcode_prog_sec, stashcode_dry_rho,                  &
                  fields_out, field_count_out, pos)
dry_rho => fields_out(pos)
CALL rcf_alloc_field( dry_rho )
CALL rcf_read_field(  dry_rho, hdr_out, decomp_rcf_output )


CALL rcf_locate(stashcode_prog_sec, stashcode_rho,                      &
                  fields_out, field_count_out, pos)
rho => fields_out(pos)
CALL rcf_alloc_field( rho )
CALL rcf_read_field(  rho, hdr_out, decomp_rcf_output )

ALLOCATE (xi3_at_rho  (output_grid % loc_p_row_length *                 &
                         output_grid % loc_p_rows,                      &
                         0:output_grid % model_levels+1))

CALL rcf_generate_heights( output_grid,orog_out,                        &
                             ppx_atm_tall, st_levels_model_rho,         &
                             xi3_at_rho , rho % level_size )

DO i = 1, dry_rho % levels

  intw_w2rho   = ( output_grid % eta_rho_levels(i)                      &
                  -output_grid % eta_theta_levels(i-1) ) /              &
                 ( output_grid % eta_theta_levels(i)                    &
                  -output_grid % eta_theta_levels(i-1) )

  DO j = 1, rho % level_size

    mixing_r =               intw_w2rho *(mv   % DATA(j,i+1)+          &
                                          mr   % DATA(j,i+1)+          &
                                          mgr  % DATA(j,i+1)+          &
                                          mcl  % DATA(j,i+1)+          &
                                          mcf  % DATA(j,i+1)+          &
                                          mcf2 % DATA(j,i+1)           &
                                         ) +                           &
                      (1.0 - intw_w2rho)*(mv   % DATA(j,i)+            &
                                          mr   % DATA(j,i)+            &
                                          mgr  % DATA(j,i)+            &
                                          mcl  % DATA(j,i)+            &
                                          mcf  % DATA(j,i)+            &
                                          mcf2 % DATA(j,i) )

    rho % DATA(j,i) = dry_rho % DATA(j,i) * (1.0+ mixing_r )

    rho % DATA(j,i) = rho % DATA(j,i) * (xi3_at_rho(j,i)*xi3_at_rho(j,i))

  END DO
END DO

DEALLOCATE (xi3_at_rho)

CALL rcf_write_field( rho, hdr_out, decomp_rcf_output )

CALL rcf_dealloc_field( rho )
CALL rcf_dealloc_field( dry_rho )
CALL rcf_dealloc_field( mcf2 )
CALL rcf_dealloc_field( mgr )
CALL rcf_dealloc_field( mr )
CALL rcf_dealloc_field( mcf )
CALL rcf_dealloc_field( mcl )
CALL rcf_dealloc_field( mv )
CALL rcf_dealloc_field( orog_out )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_recompute_wet_rho
END MODULE rcf_recompute_wet_rho_mod
#endif
