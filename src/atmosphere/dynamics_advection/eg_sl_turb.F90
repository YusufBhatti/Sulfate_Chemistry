! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sl_turb_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_SL_TURB_MOD'

CONTAINS
SUBROUTINE eg_sl_turb(                                                  &
              turb_array_size,                                          &
              row_length, rows, halo_i, halo_j, datastart, g_i_pe,      &
              e_trb, tsq_trb, qsq_trb, cov_trb,                         &
              error_code )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE mym_option_mod, ONLY: high_order_scheme_adv_turb,                   &
     monotone_scheme_adv_turb, l_high_adv_turb, l_mono_adv_turb,        &
     bdy_tke, mymodel3, tke_levels
USE um_parcore,        ONLY: mype, nproc
USE nlsizes_namelist_mod,  ONLY: global_row_length
USE um_parvars,        ONLY: nproc_x,nproc_y, at_extremity,             &
                             gc_proc_row_group, gc_proc_col_group
USE level_heights_mod, ONLY: eta_theta_levels
USE atm_fields_bounds_mod, ONLY: tdims, tdims_l
USE eg_interpolation_eta_pmf_mod, ONLY: eg_interpolation_eta_pmf
USE departure_pts_mod, ONLY: depart_xi1_w, depart_xi2_w, depart_xi3_w
USE ereport_mod, ONLY: ereport
USE Field_Types, ONLY: fld_type_p, fld_type_w
USE mpp_conf_mod,ONLY: swap_field_is_scalar
USE dynamics_input_mod, ONLY: l_sl_bc_correction

IMPLICIT NONE
!
! Description:
!   Find departure point timelevel n for turbulent variables.
!
!
! Method: ENDGame formulation version 1.01,
!         section 7.3.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_SL_TURB'

! Model dimensions

INTEGER, INTENT(IN) :: row_length, rows

! MPP options

INTEGER, INTENT(IN) ::                                                  &
  halo_i,                                                               &
                     ! Size of halo in i.
  halo_j,                                                               &
                     ! Size of halo in j.
  datastart(3),                                                         &
                     ! First gridpoints held by this processor.
  g_i_pe(1-halo_i:global_row_length+halo_i)
                     ! processor on my processor-row
                     ! holding a given value in i direction

! Loop index bounds for arrays defined on p, u, v points respectively

INTEGER, INTENT(IN) :: turb_array_size

! Timelevel n arrival point quantities

REAL, INTENT (INOUT) :: e_trb  (tdims%i_start:tdims%i_end,              &
                                tdims%j_start:tdims%j_end,              &
                                tdims%k_start:tke_levels)
REAL, INTENT (INOUT) :: tsq_trb(tdims%i_start:tdims%i_end,              &
                                tdims%j_start:tdims%j_end,              &
                                tdims%k_start:tke_levels)
REAL, INTENT (INOUT) :: qsq_trb(tdims%i_start:tdims%i_end,              &
                                tdims%j_start:tdims%j_end,              &
                                tdims%k_start:tke_levels)
REAL, INTENT (INOUT) :: cov_trb(tdims%i_start:tdims%i_end,              &
                                tdims%j_start:tdims%j_end,              &
                                tdims%k_start:tke_levels)

INTEGER, INTENT(INOUT) :: error_code  ! Non-zero on exit if error detected.

! Local variables

INTEGER :: i,j,k, number_of_inputs
INTEGER :: k_int_linear ! Linear interpolation is used at departure
                        ! points in this layer and below.
                        ! (Optional argument for subroutine
                        !  eg_interpolation_eta.)

! tmp & dummy arrays

REAL :: super_array(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,    &
                    0:tke_levels-1, turb_array_size )

REAL :: super_array_out(1:row_length,1:rows,                            &
                        0:tke_levels-1, turb_array_size)

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO j=tdims%j_start, tdims%j_end
  DO i=tdims%i_start, tdims%i_end
    super_array(i,j,0,1)  = e_trb(i,j,1)
  END DO
END DO
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(tke_levels,tdims,super_array,e_trb)
DO k=1, tke_levels-1
  DO j=tdims%j_start, tdims%j_end
    DO i=tdims%i_start, tdims%i_end
      super_array(i,j,k,1)  = e_trb(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

IF (bdy_tke == mymodel3) THEN
  DO j=tdims%j_start, tdims%j_end
    DO i=tdims%i_start, tdims%i_end
      super_array(i,j,0,2) = tsq_trb(i,j,1)
      super_array(i,j,0,3) = qsq_trb(i,j,1)
      super_array(i,j,0,4) = cov_trb(i,j,1)
    END DO
  END DO
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(tke_levels,tdims,super_array,tsq_trb,qsq_trb,cov_trb)
  DO k=1, tke_levels-1
    DO j=tdims%j_start, tdims%j_end
      DO i=tdims%i_start, tdims%i_end
        super_array(i,j,k,2) = tsq_trb(i,j,k)
        super_array(i,j,k,3) = qsq_trb(i,j,k)
        super_array(i,j,k,4) = cov_trb(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! DEPENDS ON: swap_bounds
CALL Swap_Bounds(                                                       &
             super_array, row_length, rows,                             &
             turb_array_size*(tke_levels),                              &
             halo_i, halo_j, fld_type_p, swap_field_is_scalar)

number_of_inputs = turb_array_size

! Set layers over which linear interpolation is used
IF (l_sl_bc_correction) THEN
  k_int_linear=2
ELSE
  k_int_linear=1
END IF

CALL eg_interpolation_eta_pmf(                                          &
                     eta_theta_levels,fld_type_w,                       &
                     number_of_inputs,                                  &
                     row_length, rows, tke_levels,                      &
                     rows,                                              &
                     row_length, rows, tke_levels,                      &
                     high_order_scheme_adv_turb,                        &
                     monotone_scheme_adv_turb,                          &
                     l_high_adv_turb, l_mono_adv_turb,                  &
                     depart_xi3_w, depart_xi1_w,                        &
                     depart_xi2_w, mype, nproc, nproc_x, nproc_y,       &
                     halo_i, halo_j,                                    &
                     global_row_length, datastart, at_extremity,        &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,      &
                     0, 0, error_code,                                  &
                     super_array, super_array_out,                      &
                     k_int_linear_in=k_int_linear)

! Unpack super array
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(tke_levels,tdims,super_array_out,e_trb)
DO k=0, tke_levels-1
  DO j=tdims%j_start, tdims%j_end
    DO i=tdims%i_start, tdims%i_end
      e_trb(i,j,k)  = super_array_out(i,j,k,1)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

IF (bdy_tke == 3) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(tke_levels,tdims,super_array_out,tsq_trb,qsq_trb,cov_trb)
  DO k=0, tke_levels-1
    DO j=tdims%j_start, tdims%j_end
      DO i=tdims%i_start, tdims%i_end
        tsq_trb(i,j,k)  = super_array_out(i,j,k,2)
        qsq_trb(i,j,k)  = super_array_out(i,j,k,3)
        cov_trb(i,j,k)  = super_array_out(i,j,k,4)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_sl_turb
END MODULE eg_sl_turb_mod
