! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine apply_w_forcing
!
MODULE apply_w_forcing_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='APPLY_W_FORCING_MOD'

CONTAINS

SUBROUTINE apply_w_forcing(fdims1,fdims2,fdims3,w_force,rho,field,field_inc)

! Purpose: Apply w forcing to input field
!
! Method: 
!    Use w profile to generate a large w forcing to the input field.
!    Calculation of this forcing is the same as used in the Met Office 
!    LEM/MONC cloud resolving model.

USE timestep_mod,             ONLY: timestep
USE atm_fields_bounds_mod,    ONLY: array_dims, wdims

USE horiz_grid_mod,           ONLY: intw_w2rho, intw_rho2w

USE level_heights_mod,        ONLY: z_ref_rho

USE nlsizes_namelist_mod,     ONLY: global_row_length, global_rows

USE global_2d_sums_mod,       ONLY: global_2d_sums

USE umPrintMgr,               ONLY: PrintStatus, PrStatus_Diag,            &
                                    umMessage, umPrint

USE yomhook,                  ONLY: lhook, dr_hook
USE parkind1,                 ONLY: jprb, jpim

IMPLICIT NONE

! Input
! ---------------------------------------------------------------------
TYPE(array_dims),      INTENT(IN) :: fdims1   ! array bounds for field data
TYPE(array_dims),      INTENT(IN) :: fdims2   ! for rho
TYPE(array_dims),      INTENT(IN) :: fdims3   ! output field

REAL, INTENT(IN)   ::                    &
  w_force(wdims%k_start:wdims%k_end)     & ! w profile (m/s) 
 ,rho(fdims2%i_start:fdims2%i_end,     & ! dry density (on uv levels kg/m3)
        fdims2%j_start:fdims2%j_end,   &
        fdims2%k_start:fdims2%k_end)   &
 ,field(fdims1%i_start:fdims1%i_end,   & ! input field theta or water vapour
        fdims1%j_start:fdims1%j_end,   & ! may have a halo
        fdims1%k_start:fdims1%k_end)

! Output
! ---------------------------------------------------------------------
! - increment to the field
REAL, INTENT(OUT)             ::          &
  field_inc(fdims3%i_start:fdims3%i_end,  & ! input field e.g. theta or 
            fdims3%j_start:fdims3%j_end,  & ! water vapour (no halos)
            fdims3%k_start:fdims3%k_end)

! Local variables
! ---------------------------------------------------------------------
INTEGER :: i, j, k

REAL    ::            & 
  flux_above          &  ! Flux from above for w forcing
 ,flux_below          &  ! Flux from below for w forcing
 ,tot_inc                ! Total increment for level

REAL    ::                                  &
  field_2D_ave(fdims1%k_start:fdims1%k_end)   & ! mean vertical profile of field
 ,rho_2D_ave(fdims2%k_start:fdims2%k_end)     & ! mean vertical profile of rho
 ,rho_th(fdims2%k_start:fdims2%k_end)         & ! profile of rho on theta levels
 ,w_uv(fdims2%k_start:fdims2%k_end)             ! profile of w on uv/rho levels

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER                       :: ErrorStatus

CHARACTER(LEN=*), PARAMETER   :: RoutineName='APPLY_W_FORCING'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
! ---------------------------------------------------------------------

! initialise output array to zero
! This may go from level 0 (surface) to model levels but we only want to 
! calculate actual values from level 1, so initialising to zero.
DO k = fdims3%k_start, fdims3%k_end
  DO j = fdims3%j_start, fdims3%j_end
    DO i = fdims3%i_start, fdims3%i_end
      field_inc(i,j,k) = 0.0
    END DO
  END DO
END DO

! Method being used like LEM/MONC applying a mean w forcing using a centred
! method and density weighting. 

! Need grid mean profiles of input field and rho
! Sum the field on each grid-level

CALL global_2d_sums(field, fdims1%i_len - 2 * fdims1%halo_i,            &
                           fdims1%j_len - 2 * fdims1%halo_j,            &
                           fdims1%halo_i, fdims1%halo_j,                &
                           fdims1%k_len, field_2d_ave)

CALL global_2d_sums(rho, fdims2%i_len - 2 * fdims2%halo_i,              &
                         fdims2%j_len - 2 * fdims2%halo_j,              &
                         fdims2%halo_i, fdims2%halo_j,                  &
                         fdims2%k_len, rho_2d_ave)
   
! All gridboxes same area so mean is sum divided by number of points
DO k = fdims1%k_start, fdims1%k_end
  field_2d_ave(k) = field_2d_ave(k) / (global_row_length*global_rows)
END DO
DO k = fdims2%k_start, fdims2%k_end
  rho_2d_ave(k) = rho_2d_ave(k) / (global_row_length*global_rows)
END DO

! Interpolate w to UV levels and rho to theta levels
! W held on levels 0 to model levels
! w_uv held on rho levels which will not have a zeroth level
! rho_th - may have a zeroth level but don't want to calculate a value
! for that level.
! Ensuring loop starts from level 1 and cannot go beyond dimensions for 
! rho_2d_ave (model_levels-1)

DO k = 1, fdims2%k_end-1

  w_uv(k) = intw_w2rho(k,1) * w_force(k) + intw_w2rho(k,2) * w_force(k-1)

  rho_th(k) = intw_rho2w(k,1) * rho_2d_ave(k+1) +                       &
                intw_rho2w(k,2) * rho_2d_ave(k)
     
END DO

! Set top theta level to top rho level density
! Hopefully w forcing near model top will be zero so this will not matter.

rho_th(fdims3%k_end) =  rho_2d_ave(fdims3%k_end)
k=fdims3%k_end
w_uv(k) = intw_w2rho(k,1) * w_force(k) + intw_w2rho(k,2) * w_force(k-1)

! Calculate increment for each level. Same for all gridpoints
DO k = 1, fdims3%k_end-1
  ! Exactly as LEM discretisation but UM advection scheme is semi-Lagranian
  ! so does not match LEM/MONC.

  flux_above =  -0.5*timestep*w_uv(k+1)* rho_2d_ave(k+1)*             &
                  ( field_2d_ave(k+1)  - field_2d_ave(k) )/           &
                  ( z_ref_rho(k+1) - z_ref_rho(k) )
  flux_below =  -0.5*timestep*w_uv(k)* rho_2d_ave(k)*                 &
                  ( field_2d_ave(k)  - field_2d_ave(k-1) )/           &
                  ( z_ref_rho(k+1) - z_ref_rho(k) ) 
  tot_inc = (flux_above + flux_below)/rho_th(k)
        
  DO j = fdims3%j_start, fdims3%j_end
    DO i = fdims3%i_start, fdims3%i_end
      field_inc(i,j,k) = tot_inc
    END DO
  END DO
END DO


! Optional print out for checking 
IF (PrintStatus == PrStatus_Diag) THEN

  WRITE(umMessage,'(A,A)')  ' Level  w force   w_uv    Increment, field_ave,', &
     'field, rho'
  CALL umPrint(umMessage,src=RoutineName)

  DO k=1,fdims3%k_end
    WRITE(umMessage,'(i3,6e14.6)')  k, w_force(k), w_uv(k),  &
               field_inc(1,1,k), field_2d_ave(k),field(1,1,k), rho_2d_ave(k)
    CALL umPrint(umMessage,src=RoutineName)
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE apply_w_forcing

END MODULE apply_w_forcing_mod
