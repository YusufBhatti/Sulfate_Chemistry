! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Calculate diagnostic quantities from the initial atmosphere dump
!
! Subroutine Interface:

SUBROUTINE InitDiag(                                              &
 Dummy)

USE atm_fields_bounds_mod, ONLY : pdims, tdims, udims, vdims, wdims
USE atm_step_local, ONLY: land_pts_trif, npft_trif, dim_cs1, dim_cs2, &
                          t_inc_pres, q_inc_pres, qcl_inc_pres,       &
                          qcf_inc_pres, cf_inc_pres, cfl_inc_pres,    &
                          cff_inc_pres, t_dini, q_dini, qcl_dini,     &
                          qcf_dini, cf_dini, cfl_dini, cff_dini
USE eot_increments_mod, ONLY:                                         &
    eot_incs_init, eot_incs_calc, eot_incs_dealloc, eot_inc_rho
USE trignometric_mod, ONLY: sin_v_latitude
USE wet_to_dry_n_calc_mod
USE jules_surface_types_mod, ONLY: ntype, npft

USE jules_snow_mod, ONLY: nsmax

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE Control_Max_Sizes
USE lbc_mod
USE submodel_mod, ONLY: atmos_sm, atmos_im
USE stash_array_mod, ONLY: stash_maxlen, sf
USE free_tracers_inputs_mod, ONLY: a_tracer_last, a_tracer_first
USE jules_sea_seaice_mod, ONLY: nice, nice_use

USE dynamics_testing_mod, ONLY: l_physics
USE cloud_inputs_mod,  ONLY: i_cld_vn
USE pc2_constants_mod, ONLY: i_cld_off, i_cld_pc2

USE nlsizes_namelist_mod, ONLY:                                        &
    a_len1_coldepc, a_len1_flddepc, a_len1_levdepc, a_len1_rowdepc,    &
    a_len2_coldepc, a_len2_flddepc, a_len2_levdepc, a_len2_lookup,     &
    a_len2_rowdepc, a_len_cfi1, a_len_cfi2, a_len_cfi3, a_len_extcnst, &
    a_len_inthd, a_len_realhd, land_field, len1_lookup,                &
    len_dumphist, len_fixhd, len_tot, model_levels,                    &
    mpp_len1_lookup, n_cca_lev, n_obj_d1_max, n_rows,                  &
    ntiles, river_row_length, river_rows, row_length, rows, sm_levels, &
    theta_off_size, tr_levels, tr_ukca, tr_vars
USE atm_fields_mod, ONLY: u, v, q, qcl, qcf, qcf2, qrain, qgraup
USE field_types, ONLY: fld_type_u, fld_type_v
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod, ONLY: swap_field_is_vector
USE pws_diags_mod, ONLY: flag_upd_helicity_5k

USE d1_array_mod, ONLY: d1

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!
! Description:
!   InitDiag processes diagnostic requests from the initial atmosphere
!   dump, including both prognostic variables resident in D1 - UM
!   STASH section 0 - and fields derived from physics and dynamics
!   variables, as calculated in UM STASH sections 15 and 16.
!
! Method:
!   1. Call STASH to process diagnostics requests for sections 0,33,34
!   2. Process dynamics-derived diagnostics (section 15):
!    2a. Check validity of diagnostic requests (avoids re-checking
!        within every later call of St_diag1).
!    2b. Call St_diag1 as interface to: Dyn_diag calculations of
!        dynamics derived variables and STASH (section 15).
!   3. Process physics-derived diagnostics (section 16):
!    3a. Check validity of diagnostic requests (avoids re-checking
!        within every later call of St_diag2).
!    3b. Call St_diag2 as interface to: Phy_diag calculations of
!        physics derived variables and STASH (section 16).
!   4. Process climate diagnostics (section 30):
!    4a. Check validity of diagnostic requests (avoids re-checking
!        within every later call of St_diag3).
!    4b. Call St_diag3 as interface to: EOT_diag calculations of
!        climate variables and STASH (section 30).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!

! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER :: Dummy            ! Not used, needed to end arg list

!   Array  arguments with intent(in):

!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):

! Local parameters:
CHARACTER(LEN=*), PARAMETER :: RoutineName='INITDIAG'

! Local scalars:
!   ErrorStatus
INTEGER :: ErrorStatus          ! Error flag (0 = OK)
CHARACTER(LEN=errormessagelength) :: CMessage  ! Error message if return code >0

INTEGER ::    &
 im_index,    & !  Internal Model Index for stash arrays
 i,j,k          !  Loop indices

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Local dynamic arrays:
REAL, ALLOCATABLE :: t_incr_diagnostic(:,:,:)
REAL, ALLOCATABLE :: q_incr_diagnostic(:,:,:)
REAL, ALLOCATABLE :: qcl_incr_diagnostic(:,:,:)

! Local dynamic arrays for section 30
REAL :: energy_corr_now
REAL, ALLOCATABLE :: STASHwork30(:)
REAL :: co2_mmr
REAL :: tot_dry_mass_final

!- End of header

!   0. Initialisation

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
ErrorStatus = 0
Cmessage=''
im_index = 1

!----------------------------------------------------------------------
!   1. Call STASH to process diagnostics requests for section 0.

! DEPENDS ON: stash
CALL stash(atmos_sm,atmos_im,0,d1,                                 &
 ErrorStatus,Cmessage)

! -------------------------------------------------------------------
!   1a. Call STASH to process diagnostics requests for section 33.

IF ( sf(0,33) )                                                    &
! DEPENDS ON: stash
  CALL stash(atmos_sm,atmos_im,33,d1,                              &
          ErrorStatus,Cmessage)

! -------------------------------------------------------------------
!   1b. Call STASH to process diagnostics requests for section 34.

IF ( sf(0,34) )                                                    &
! DEPENDS ON: stash
  CALL stash(atmos_sm,atmos_im,34,d1,                              &
          ErrorStatus,Cmessage)

!----------------------------------------------------------------------
!   2. Process dynamics-derived diagnostics (section 15):
!    2a. Check validity of diagnostic requests (avoids re-checking
!        within every later call of St_diag1).
!    2b. Call St_diag1 as interface to: Dyn_diag calculations of
!        dynamics derived variables and STASH (section 15).

IF ( sf(0,15) .AND. ErrorStatus == 0) THEN

  ! ----------------------------------------------------------------------
  ! Add ability to get increments from qt_bal_cld call and output in section 15
  ! This is only valid from atm_step; here zeroed to avoid potential errors.
  ! ----------------------------------------------------------------------

  IF (sf(181,15) ) THEN
    ALLOCATE ( T_incr_diagnostic(row_length,rows,model_levels) )
    DO k=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          T_incr_diagnostic(i,j,k) = 0.0
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( T_incr_diagnostic(1,1,1) )
  END IF

  IF (sf(182,15) ) THEN
    ALLOCATE ( q_incr_diagnostic(row_length,rows,model_levels) )
    DO k=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          q_incr_diagnostic(i,j,k) = 0.0
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( q_incr_diagnostic(1,1,1) )
  END IF

  IF (sf(183,15) ) THEN
    ALLOCATE ( qcl_incr_diagnostic(row_length,rows,model_levels) )
    DO k=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          qcl_incr_diagnostic(i,j,k) = 0.0
        END DO ! i
      END DO ! j
    END DO ! k
  ELSE
    ALLOCATE ( qcl_incr_diagnostic(1,1,1) )
  END IF

  ! DEPENDS ON: st_diag1
  CALL St_diag1( stash_maxlen(15,im_index),                       &
    t_incr_diagnostic,q_incr_diagnostic,qcl_incr_diagnostic,       &
    ErrorStatus,Cmessage)

  DEALLOCATE ( t_incr_diagnostic)
  DEALLOCATE ( q_incr_diagnostic)
  DEALLOCATE ( qcl_incr_diagnostic)

END IF      ! Diagnostics required for this section

!----------------------------------------------------------------------
!   3. Process physics-derived diagnostics (section 16):
!    3a. Check validity of diagnostic requests (avoids re-checking
!        within every later call of St_diag2).
!    3b. Call St_diag2 as interface to: Phy_diag calculations of
!        physics derived variables and STASH (section 16).

IF ( sf(0,16) .AND. ErrorStatus == 0) THEN

  IF ( i_cld_vn == i_cld_pc2 ) THEN

    ALLOCATE(t_inc_pres  (row_length,rows,model_levels))
    ALLOCATE(q_inc_pres  (row_length,rows,model_levels))
    ALLOCATE(qcl_inc_pres(row_length,rows,model_levels))
    ALLOCATE(qcf_inc_pres(row_length,rows,model_levels))
    ALLOCATE(cf_inc_pres (row_length,rows,model_levels))
    ALLOCATE(cfl_inc_pres(row_length,rows,model_levels))
    ALLOCATE(cff_inc_pres(row_length,rows,model_levels))
    ALLOCATE(t_dini      (row_length,rows,model_levels))
    ALLOCATE(q_dini      (row_length,rows,model_levels))
    ALLOCATE(qcl_dini    (row_length,rows,model_levels))
    ALLOCATE(qcf_dini    (row_length,rows,model_levels))
    ALLOCATE(cf_dini     (row_length,rows,model_levels))
    ALLOCATE(cfl_dini    (row_length,rows,model_levels))
    ALLOCATE(cff_dini    (row_length,rows,model_levels))

  END IF

  ! DEPENDS ON: st_diag2
  CALL St_diag2( stash_maxlen(16,im_index),                       &
    ErrorStatus,Cmessage)

  IF ( i_cld_vn == i_cld_pc2 ) THEN

    DEALLOCATE(cff_dini)
    DEALLOCATE(cfl_dini)
    DEALLOCATE(cf_dini)
    DEALLOCATE(qcf_dini)
    DEALLOCATE(qcl_dini)
    DEALLOCATE(q_dini)
    DEALLOCATE(t_dini)
    DEALLOCATE(cff_inc_pres)
    DEALLOCATE(cfl_inc_pres)
    DEALLOCATE(cf_inc_pres)
    DEALLOCATE(qcf_inc_pres)
    DEALLOCATE(qcl_inc_pres)
    DEALLOCATE(q_inc_pres)
    DEALLOCATE(t_inc_pres)
   
  END IF

END IF      ! Diagnostics required for this section

!----------------------------------------------------------------------
!   4. Process climate diagnostics (section 30):
!    4a. Check validity of diagnostic requests (avoids re-checking
!        within every later call of St_diag3).
!    4b. Call St_diag3 as interface to: EOT_diag calculations of
!        climate variables and STASH (section 30).

IF ( (sf(0,30) .OR. flag_upd_helicity_5k) .AND. ErrorStatus == 0) THEN

  CALL eot_incs_init()
  CALL eot_incs_calc()

  energy_corr_now=0.0

  ! Initialise arrays as they will end up with rubbish in them otherwise
  eot_inc_rho(:,:,:)     = 0.0
  co2_mmr                = 0.0
  tot_dry_mass_final     = 0.0

  ! Calculate wet_to_dry_n for total column integrals
  CALL wet_to_dry_n_calc(q,qcl,qcf,qcf2,qrain,qgraup)

! Need full swap_bounds for u,v for zeroth timestep 
! to avoid processor edge effects in section 30 wind diagnostics.

  CALL swap_bounds(u,row_length,rows,model_levels,                &
                   offx,offy,fld_type_u,swap_field_is_vector)

  CALL swap_bounds(v,row_length,n_rows,model_levels,              &
                   offx,offy,fld_type_v,swap_field_is_vector)

  ! size of diagnostic space
  ALLOCATE (STASHwork30(STASH_maxlen(30,atmos_im)))

  ! DEPENDS ON: st_diag3
  CALL St_diag3(STASHwork30,stash_maxlen(30,im_index),             &
    energy_corr_now, sin_v_latitude,                               &
    wet_to_dry_n, co2_mmr,                                         &
    ErrorStatus,Cmessage)

  CALL eot_incs_dealloc()

  CALL destroy_wet_to_dry_n()

END IF      ! Diagnostics required for this section

!----------------------------------------------------------------------

! Check error condition
IF (ErrorStatus >  0) THEN

  CALL Ereport(RoutineName,ErrorStatus,Cmessage)
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE InitDiag
