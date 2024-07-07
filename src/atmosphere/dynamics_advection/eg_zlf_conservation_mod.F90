! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_zlf_conservation_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                &
                  ModuleName='EG_ZLF_CONSERVATION_MOD'

CONTAINS
SUBROUTINE eg_zlf_conservation(tracers_in, tracers_out,                &
              tr_array_size, tr_zlf_size,                              &
              row_length, rows, n_rows, model_levels, halo_i, halo_j,  &
              offx, offy, datastart, g_i_pe, high_order_scheme,        &
              monotone_scheme,  l_high, l_mono,                        &
              zlf_conserv_option, zlf_cfl_top_level,                   &
              rho_n, rho_np1, error_code                               )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE um_parcore,            ONLY: mype, nproc
USE nlsizes_namelist_mod,  ONLY: global_row_length
USE um_parvars,        ONLY: nproc_x,nproc_y, at_extremity,           &
                             gc_proc_row_group, gc_proc_col_group

USE level_heights_mod, ONLY: eta_theta_levels, eta_rho_levels
USE atm_fields_bounds_mod
USE eg_interpolation_eta_pmf_mod, ONLY: eg_interpolation_eta_pmf
USE horiz_grid_mod
USE metric_terms_mod
USE departure_pts_mod
USE Field_Types
USE dynamics_input_mod, ONLY: l_sl_bc_correction, apply_ocf, apply_adas
USE enforce_mono_mz_mod
USE eg_mass_conserv_mod
USE umPrintMgr, ONLY: umMessage, umPrint
USE eg_zlf_mod, ONLY: eg_zlf_setup,  eg_zlf_zero_rim,            &
                      eg_zlf_overwrite_rim,                      &
                      zlf_zeros_rim_default, zlf_zeros_rim_step, &
                      zlf_linear_boundaries, zlf_np_overwrite,   &
                      zlf_np_linear, zlf_para_set 
USE atm_step_local, ONLY: first_atmstep_call
IMPLICIT NONE
!
! Description:
!    Apply the ZLF scheme (or LAM mass conservation) for the tracer array
!          tracers_in  -->  tracers_out
!
!    Note that conservation is applied only for part of the tracer array
!    which is tracers_all(:,:,:,1:tr_zlf_size) while the rest of
!    the array tracers_all(:,:,:,tr_zlf_size+1:tr_array_size) is
!    simply advected without conservation/ZLF.
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
!   This code is written to UM programming standards

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_ZLF_CONSERVATION'


INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels,    &
                       tr_array_size, tr_zlf_size,                &
                       zlf_conserv_option,                        &
                       halo_i, halo_j, offx, offy, datastart(3),  &
                       g_i_pe(1-halo_i:global_row_length+halo_i), &
                       high_order_scheme, monotone_scheme,        &
                       zlf_cfl_top_level
LOGICAL, INTENT(IN) :: l_high,  l_mono

REAL, INTENT (IN) ::   tracers_in(tdims_l%i_start:tdims_l%i_end,  &
                                  tdims_l%j_start:tdims_l%j_end,  &
                                  tdims_l%k_start:tdims_l%k_end,  &
                                  tr_array_size                   )
REAL, INTENT (OUT) :: tracers_out(1-offx:row_length+offx,         &
                                  1-offy:rows+offy,               &
                                  0:model_levels, tr_array_size   )

INTEGER, INTENT (INOUT) :: error_code  

REAL, INTENT(IN) ::   rho_n(pdims_s%i_start:pdims_s%i_end,        &
                            pdims_s%j_start:pdims_s%j_end,        &
                            pdims_s%k_start:pdims_s%k_end)
REAL, INTENT(IN) :: rho_np1(pdims_s%i_start:pdims_s%i_end,        &
                            pdims_s%j_start:pdims_s%j_end,        &
                            pdims_s%k_start:pdims_s%k_end)

! Local variables

INTEGER :: i,j,k, tr_nozlf_size, number_of_inputs 
INTEGER :: k_int_linear 
LOGICAL :: l_high_in  
LOGICAL :: l_mono_in 

! local arrays

REAL, ALLOCATABLE ::       tr_zlf(:,:,:,:)
REAL, ALLOCATABLE ::      tr_zlf2(:,:,:,:)
REAL, ALLOCATABLE ::     tr_nozlf(:,:,:,:)
REAL, ALLOCATABLE ::   tr_zlf_out(:,:,:,:)
REAL, ALLOCATABLE ::     tr_zlf_n(:,:,:,:)
REAL, ALLOCATABLE ::  tr_zlf_out2(:,:,:,:)
REAL, ALLOCATABLE :: tr_nozlf_out(:,:,:,:)
REAL, ALLOCATABLE :: mz_temp_qsmin(:)

LOGICAL :: l_conserv_smooth_lap 
INTEGER :: na_zlf,na_nozlf,offxs,offys,i_tr
INTEGER :: zlf_conserv_option_local


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

tr_nozlf_size = tr_array_size - tr_zlf_size

! split the full array into zlf array and no-zlf array
! Note that the tracers not treated with ZLF are stacked first:
!     i.e.,  tracers_in(:,:,:,1:tr_nozlf_size)
! and the tracers treated with ZLF are stacked second:
!    i.e.,  tracers_in(:,:,:,tr_nozlf_size+1:tr_array_size)
!

IF ( first_atmstep_call .AND. (.NOT. zlf_para_set) ) THEN
  CALL eg_zlf_setup(tdims%k_start, tdims%k_end)
END IF

zlf_conserv_option_local = zlf_conserv_option 
          
na_zlf   = MAX(1,tr_zlf_size         )
na_nozlf = MAX(1,tr_nozlf_size       )
offxs    = MAX(1-tdims_s%i_start,offx)
offys    = MAX(1-tdims_s%j_start,offy)
      
ALLOCATE(tr_nozlf(     tdims_l%i_start:tdims_l%i_end,           &
                       tdims_l%j_start:tdims_l%j_end,           &
                       tdims_l%k_start:tdims_l%k_end, na_nozlf )) 
ALLOCATE(tr_nozlf_out( tdims_s%i_start:tdims_s%i_end,           &
                       tdims_s%j_start:tdims_s%j_end,           &
                       tdims_s%k_start:tdims_s%k_end, na_nozlf ))   
                                                            
ALLOCATE( tr_zlf(      tdims_l%i_start:tdims_l%i_end,           &
                       tdims_l%j_start:tdims_l%j_end,           &
                       tdims_l%k_start:tdims_l%k_end, na_zlf   ))      
ALLOCATE(tr_zlf2(      tdims_l%i_start:tdims_l%i_end,           &
                       tdims_l%j_start:tdims_l%j_end,           &
                       tdims_l%k_start:tdims_l%k_end, na_zlf   ))
ALLOCATE(tr_zlf_out(   tdims_s%i_start:tdims_s%i_end,           &
                       tdims_s%j_start:tdims_s%j_end,           &
                       tdims_s%k_start:tdims_s%k_end, na_zlf   ))
ALLOCATE(tr_zlf_n(     tdims_s%i_start:tdims_s%i_end,           &
                       tdims_s%j_start:tdims_s%j_end,           &
                       tdims_s%k_start:tdims_s%k_end, na_zlf   ))
ALLOCATE(tr_zlf_out2(  tdims_s%i_start:tdims_s%i_end,           &
                       tdims_s%j_start:tdims_s%j_end,           &
                       tdims_s%k_start:tdims_s%k_end, na_zlf   ))     
ALLOCATE(mz_temp_qsmin(na_zlf))
  
 
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i_tr,k,j,i)                &
!$OMP SHARED(tr_nozlf_size, tdims_l, tr_nozlf, tracers_in, tr_zlf_size, &
!$OMP        tr_zlf)
IF ( tr_nozlf_size > 0 ) THEN 
  DO i_tr = 1, tr_nozlf_size
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims_l%k_start, tdims_l%k_end
      DO j = tdims_l%j_start, tdims_l%j_end
        DO i = tdims_l%i_start, tdims_l%i_end            
          tr_nozlf(i,j,k,i_tr) = tracers_in(i,j,k,i_tr)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
END IF
  
IF ( tr_zlf_size > 0 ) THEN
  DO i_tr = 1, tr_zlf_size
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims_l%k_start, tdims_l%k_end
      DO j = tdims_l%j_start, tdims_l%j_end
        DO i = tdims_l%i_start, tdims_l%i_end
         tr_zlf(i,j,k,i_tr) = tracers_in(i,j,k,tr_nozlf_size+i_tr)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO                           
END IF
!$OMP END PARALLEL

! Set layers over which linear interpolation is used
IF (l_sl_bc_correction) THEN
  k_int_linear=2
ELSE
  k_int_linear=1
END IF

!
!  Apply the ZLF scheme to impose mass conservative advection
!  This scheme involves the following steps:
!
!  (a) zeros the external halos and part of the rim zone
!  (b) advect the modified field
!  (c) apply the appropriate monotonicity scheme
!  (d) aplly the appropriate mass conservation scheme
!  (e) overwite the rim zone with the original field
! 
                 
IF ( tr_zlf_size > 0 ) THEN ! apply ZLF for "tr_zlf" tarcers
      
  IF ( zlf_zeros_rim_default ) THEN        
    ! tr_zlf2 = tr_zlf  
    ! this original array without the zeroing
    ! to be used for overwriting the rim
    ! after advecting the modified field
       
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i_tr, k, j, i)       &
!$OMP SHARED(tr_zlf_size, tdims_l, tr_zlf2, tr_zlf)
    DO i_tr = 1, tr_zlf_size
!$OMP DO SCHEDULE(STATIC)
      DO k = tdims_l%k_start, tdims_l%k_end
        DO j = tdims_l%j_start, tdims_l%j_end
          DO i = tdims_l%i_start, tdims_l%i_end
            tr_zlf2(i,j,k,i_tr) = tr_zlf(i,j,k,i_tr)
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO
!$OMP END PARALLEL
       
    CALL eg_zlf_zero_rim(tr_zlf,                      &
                         1-halo_i, row_length+halo_i, &
                         1-halo_j, rows+halo_j,       &
                         0, model_levels,             &
                         tr_zlf_size,                 &
                         zlf_cfl_top_level,           &
                         zlf_conserv_option_local     )
  END IF 
     
  ! This returns a high-order solution for zlf_ed-tracers
  
  CALL eg_interpolation_eta_pmf(                                  &
                 eta_theta_levels,fld_type_w,                     &
                 tr_zlf_size,                                     &
                 row_length, rows, model_levels+1,                &
                 rows,                                            &
                 row_length, rows, model_levels+1,                &
                 high_order_scheme, monotone_scheme,              &
                 l_high, l_mono, depart_xi3_w, depart_xi1_w,      &
                 depart_xi2_w, mype, nproc, nproc_x, nproc_y,     &
                 halo_i, halo_j,                                  &
                 global_row_length, datastart, at_extremity,      &
                 g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                 offxs, offys, error_code,                        &
                 tr_zlf, tr_zlf_out,                              &
                 k_int_linear_in=k_int_linear)
                        
  ! this to return a the low-order (linear) solution for zlf tracers
                      
  l_high_in=.FALSE.
  l_mono_in=.TRUE.
  CALL eg_interpolation_eta_pmf(                                   &
                  eta_theta_levels,fld_type_w,                     &
                  tr_zlf_size,                                     &
                  row_length, rows, model_levels+1,                &
                  rows,                                            &
                  row_length, rows, model_levels+1,                &
                  0, 1,                                            &
                  l_high_in,l_mono_in,                             &
                  depart_xi3_w, depart_xi1_w,                      &
                  depart_xi2_w, mype, nproc, nproc_x, nproc_y,     &
                  halo_i, halo_j,                                  &
                  global_row_length, datastart, at_extremity,      &
                  g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                  offxs, offys, error_code,                        &
                  tr_zlf, tr_zlf_out2,                             &
                  k_int_linear_in=k_int_linear)
                      
  IF ( zlf_zeros_rim_step .AND. zlf_linear_boundaries ) THEN  
      
    CALL eg_zlf_overwrite_rim(tr_zlf_out,                    &
                              tr_zlf_out2,                   &
                              1-offxs, row_length+offxs,     &
                              1-offys, rows+offys,           &
                              0, model_levels,               &
                              1-offxs, row_length+offxs,     &
                              1-offys, rows+offys,           &
                              tr_zlf_size,                   &
                              zlf_np_linear                  )
                                        
  END IF
      
   ! copy "tr_zlf" in smaller array to be used in the
   ! calculation of the mass at time step n
      
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i_tr, k, j, i)          &
!$OMP SHARED(na_zlf, tdims, tr_zlf_n, tr_zlf)
  DO i_tr = 1, na_zlf
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          tr_zlf_n(i,j,k,i_tr) = tr_zlf(i,j,k,i_tr)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
!$OMP END PARALLEL
     
  !  apply the conservation step
  !  zlf_conserv_option_local==1 apply OCF scheme
  !  zlf_conserv_option_local==2 apply ADAS scheme
      
  IF ( zlf_conserv_option_local==apply_ocf ) THEN   
        
    CALL enforce_conserv_OCF(tr_zlf_out,                &
                             tr_zlf_out2,               &
                             tr_zlf_n, rho_n, rho_np1,  &     
                             1-offxs, row_length+offxs, &
                             1-offys, rows+offys,       &
                             0, model_levels,           &
                             1, row_length,             &
                             1, rows,                   &
                             0, model_levels,           &
                             tr_zlf_size                )
      
 
  ELSE IF ( zlf_conserv_option_local==apply_adas ) THEN
         
    mz_temp_qsmin(:)     = 0.0                   

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i_tr, k, j, i)            &
!$OMP SHARED(na_zlf, tdims_s, tr_zlf_out2)
    DO i_tr = 1, na_zlf
!$OMP DO SCHEDULE(STATIC)
      DO k = tdims_s%k_start, tdims_s%k_end
        DO j = tdims_s%j_start, tdims_s%j_end
          DO i = tdims_s%i_start, tdims_s%i_end
            tr_zlf_out2(i,j,k,i_tr) = 0.0 ! this now contains zero sources
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO
!$OMP END PARALLEL

    L_conserv_smooth_lap = .TRUE.
         
    CALL eg_mass_conservation_fix(rho_n, rho_np1,       &
                                  tr_zlf_n,             &
                                  tr_zlf_out,           &
                                  tr_zlf_out2,          &
                                  mz_temp_qsmin,        &
                                  tr_zlf_size,          &
                                  L_conserv_smooth_lap  )
  END IF  
      
  ! overwrite the rim zone with the fields at time step n
       
  IF ( zlf_zeros_rim_step  ) THEN  
      
    CALL eg_zlf_overwrite_rim(tr_zlf_out,                  &
                              tr_zlf2,                     &
                              1-offxs, row_length+offxs,   &
                              1-offys, rows+offys,         &
                              0, model_levels,             &
                              1-halo_i, row_length+halo_i, &
                              1-halo_j, rows+halo_j,       &
                              tr_zlf_size,                 &
                              zlf_np_overwrite             )
  END IF
      
END IF  ! end zlf tracers
 
! This returns a high-order non-conservative solution for 
! non_zlf_ed-tracers
      
IF ( tr_nozlf_size > 0 ) THEN  
                    
  CALL eg_interpolation_eta_pmf(                                 &
                eta_theta_levels,fld_type_w,                     &
                tr_nozlf_size,                                   &
                row_length, rows, model_levels+1,                &
                rows,                                            &
                row_length, rows, model_levels+1,                &
                high_order_scheme, monotone_scheme,              &
                l_high, l_mono, depart_xi3_w, depart_xi1_w,      &
                depart_xi2_w, mype, nproc, nproc_x, nproc_y,     &
                halo_i, halo_j,                                  &
                global_row_length, datastart, at_extremity,      &
                g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                offxs, offys, error_code,                        &
                tr_nozlf, tr_nozlf_out,                          &
                k_int_linear_in=k_int_linear)
END IF
      
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i_tr, k, j, i)             &
!$OMP SHARED(tr_nozlf_size, tdims, tracers_out, tr_nozlf_out,   &
!$OMP        tr_zlf_size, tr_zlf_out)
IF ( tr_nozlf_size > 0 ) THEN 
  DO i_tr = 1, tr_nozlf_size
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          tracers_out(i,j,k,i_tr) = tr_nozlf_out(i,j,k,i_tr)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO      
END IF
IF ( tr_zlf_size > 0 ) THEN
  DO i_tr = 1, tr_zlf_size      
!$OMP DO SCHEDULE(STATIC)
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          tracers_out(i,j,k,tr_nozlf_size+i_tr) = tr_zlf_out(i,j,k,i_tr)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END DO
END IF
!$OMP END PARALLEL
      
DEALLOCATE( mz_temp_qsmin )
DEALLOCATE( tr_zlf_out2   )
DEALLOCATE( tr_zlf_n      )
DEALLOCATE( tr_zlf_out    )
DEALLOCATE( tr_zlf2       )
DEALLOCATE( tr_zlf        )
DEALLOCATE( tr_nozlf_out  )
DEALLOCATE( tr_nozlf      )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_zlf_conservation
END MODULE eg_zlf_conservation_mod
