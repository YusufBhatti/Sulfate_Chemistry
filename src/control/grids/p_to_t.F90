! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! SUBROUTINE P_TO_T
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Grids
MODULE p_to_t_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='P_TO_T_MOD'

CONTAINS

SUBROUTINE p_to_t                                                              &
  ( row_length, rows, halo_i, halo_j, halo_i_data, halo_j_data, levels         &
  , r_theta_levels, r_rho_levels, field_in, field_out )


! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v10.3 programming standards.
!   This routine is thread safe and can be called from an 
!   OMP parallel region

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE compute_chunk_size_mod, ONLY: compute_chunk_size

IMPLICIT NONE

! Arguments with intent in. ie: input variables.
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: levels
INTEGER, INTENT(IN) :: halo_i
INTEGER, INTENT(IN) :: halo_j
INTEGER, INTENT(IN) :: halo_i_data
INTEGER, INTENT(IN) :: halo_j_data

REAL, INTENT(IN)  :: r_theta_levels ( 1-halo_i:row_length+halo_i               &
                                    , 1-halo_j:rows+halo_j, 0:levels+1 )

REAL, INTENT(IN)  :: r_rho_levels   ( 1-halo_i:row_length+halo_i               &
                                    , 1-halo_j:rows+halo_j, levels+1 )

REAL, INTENT(IN)  :: field_in ( 1-halo_i_data:row_length+halo_i_data           &
                              , 1-halo_j_data:rows+halo_j_data, levels+1 )


! Arguments with intent out. ie: output variables.
REAL, INTENT(OUT) :: field_out ( 1-halo_i_data:row_length+halo_i_data          &
                               , 1-halo_j_data:rows+halo_j_data, levels )


! Local variables.

INTEGER :: i, j, k
INTEGER :: ompt_start, ompt_end

REAL :: weight_1
REAL :: weight_2
REAL :: weight_3

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'P_TO_T'

! ----------------------------------------------------------------------
! Interpolate field_in (on P grid) to T grid (field_out).
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! This procedure is now thread safe. i.e. it can be called from an
! OMP parallel region. Compute_chunk_size will determine the start and end
! iterations for each thread as if it were a static schedule.
! If called from a serial region, start and end are 1 and levels

ompt_start = 1
ompt_end   = levels
! Only call compute_chunk_size if compiling with OMP
! The procedure call is protected by the optional compile
! sentinel
!$ CALL compute_chunk_size(1,levels,ompt_start,ompt_end)

DO k=ompt_start, ompt_end
  DO j = 1-halo_j_data, rows+halo_j_data
    DO i = 1-halo_i_data, row_length+halo_i_data

      weight_1 = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
      weight_2 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
      weight_3 = r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k)
      field_out (i,j, k) =                                                     &
               weight_2/weight_1 * field_in (i,j,k+1)                          &
             + weight_3/weight_1 * field_in (i,j,k)
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE p_to_t

END MODULE p_to_t_mod
