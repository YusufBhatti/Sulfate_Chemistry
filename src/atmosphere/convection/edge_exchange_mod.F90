! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Purpose: To propagate across processor boundaries, via halo swapping.
!
!  Method: E-W and N-S halos and halo corners are copied into the
!  body of work arrays, to then be swapped into the haloes of
!  neighbouring processors by a standard 2D swap-bounds routine.
!  Input from swapping is then combined with non-halo data in the
!  I/O array(s).
!  This procedure *assumes* that the non-halo domain is at least
!  twice the size of the halo width, in each direction.
!
!  Programming standard : UMDP 3
!
!  Documentation: UMDP 025
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Convection
!---------------------------------------------------------------------
MODULE edge_exchange_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'EDGE_EXCHANGE_MOD'
CONTAINS

SUBROUTINE edge_exchange( n_halo_fields, l_vector, l_max_not_add, io_array )

USE atm_fields_bounds_mod,  ONLY : pdims, pdims_s
USE field_types,            ONLY : fld_type_p
USE swap_bounds_2d_mv_mod,  ONLY : swap_bounds_2d_mv
USE swapable_field_mod,     ONLY : swapable_field_pointer_type

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) :: n_halo_fields

LOGICAL, INTENT(IN) ::   l_vector(     n_halo_fields)         &
                       , l_max_not_add(n_halo_fields)

REAL, INTENT(INOUT) ::                                        &
                    io_array( pdims_s%i_start:pdims_s%i_end,  &
                              pdims_s%j_start:pdims_s%j_end,  &
                              n_halo_fields )

! Local variables
REAL, TARGET ::                                               &
                  work_array( pdims_s%i_start:pdims_s%i_end,  &
                              pdims_s%j_start:pdims_s%j_end,  &
                              3 * n_halo_fields )

! counter for swap_bounds
INTEGER :: i_field

TYPE(swapable_field_pointer_type) :: fields_to_swap(3 * n_halo_fields)
! for multivar swap_bounds

INTEGER :: i, j, i_hf, hf_pos
! Loop indexes

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EDGE_EXCHANGE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! zero the work array
!-----------------------
DO i_hf = 1, 3 * n_halo_fields
  DO j = pdims_s%j_start, pdims_s%j_end
    DO i = pdims_s%i_start, pdims_s%i_start
      work_array(i, j, i_hf) = 0.0
    END DO
  END DO
END DO

! copy the 'in' halos to the 'work' edges
!------------------------------------------
DO i_hf = 1, n_halo_fields

  hf_pos = 3 * (i_hf-1)

! E-W
!-----
  hf_pos = hf_pos + 1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims_s%i_start, pdims%i_start - 1
      work_array(i + pdims_s%halo_i,j,hf_pos) = io_array(i,j,i_hf)
    END DO
  END DO
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_end + 1, pdims_s%i_end
      work_array(i - pdims_s%halo_i,j,hf_pos) = io_array(i,j,i_hf)
    END DO
  END DO

! N-S
!-----
  hf_pos = hf_pos + 1
  DO j = pdims_s%j_start, pdims%j_start - 1
    DO i = pdims%i_start, pdims%i_end
      work_array(i,j + pdims_s%halo_j,hf_pos) = io_array(i,j,i_hf)
    END DO
  END DO
  DO j = pdims%j_end + 1, pdims_s%j_end
    DO i = pdims%i_start, pdims%i_end
      work_array(i,j - pdims_s%halo_j,hf_pos) = io_array(i,j,i_hf)
    END DO
  END DO

! halo corners
!---------------
  hf_pos = hf_pos + 1
  DO j = pdims_s%j_start, pdims%j_start - 1
    DO i = pdims_s%i_start, pdims%i_start - 1
      work_array(i + pdims_s%halo_i,    &
                 j + pdims_s%halo_j,hf_pos) = io_array(i,j,i_hf)
    END DO
  END DO
  DO j = pdims%j_end + 1, pdims_s%j_end
    DO i = pdims_s%i_start, pdims%i_start - 1
      work_array(i + pdims_s%halo_i,    &
                 j - pdims_s%halo_j,hf_pos) = io_array(i,j,i_hf)
    END DO
  END DO
  DO j = pdims_s%j_start, pdims%j_start - 1
    DO i = pdims%i_end + 1, pdims_s%i_end
      work_array(i - pdims_s%halo_i,    &
                 j + pdims_s%halo_j,hf_pos) = io_array(i,j,i_hf)
    END DO
  END DO
  DO j = pdims%j_end + 1, pdims_s%j_end
    DO i = pdims%i_end + 1, pdims_s%i_end
      work_array(i - pdims_s%halo_i,    &
                 j - pdims_s%halo_j,hf_pos) = io_array(i,j,i_hf)
    END DO
  END DO

END DO

! swap bounds in 2D
!--------------------
! fill up the fields array
!----------------------------
i_field = 0

DO i_hf = 1, n_halo_fields

  hf_pos = 3 * (i_hf-1)

  hf_pos = hf_pos + 1
  i_field = i_field + 1
  fields_to_swap(i_field) % field_2d    => work_array(:,:,hf_pos)
  fields_to_swap(i_field) % field_type  =  fld_type_p
  fields_to_swap(i_field) % levels      =  1
  fields_to_swap(i_field) % rows        =  pdims%j_end
  fields_to_swap(i_field) % vector      =  l_vector( i_hf )

  hf_pos = hf_pos + 1
  i_field = i_field + 1
  fields_to_swap(i_field) % field_2d    => work_array(:,:,hf_pos)
  fields_to_swap(i_field) % field_type  =  fld_type_p
  fields_to_swap(i_field) % levels      =  1
  fields_to_swap(i_field) % rows        =  pdims%j_end
  fields_to_swap(i_field) % vector      =  l_vector( i_hf )

  hf_pos = hf_pos + 1
  i_field = i_field + 1
  fields_to_swap(i_field) % field_2d    => work_array(:,:,hf_pos)
  fields_to_swap(i_field) % field_type  =  fld_type_p
  fields_to_swap(i_field) % levels      =  1
  fields_to_swap(i_field) % rows        =  pdims%j_end
  fields_to_swap(i_field) % vector      =  l_vector( i_hf )

END DO

! single call to swap bounds
!------------------------------
CALL swap_bounds_2d_mv( fields_to_swap, i_field,                      &
                        pdims%i_end, pdims_s%halo_i, pdims_s%halo_j)


! add/max the 'work' halos with the 'in' edges
!-----------------------------------------------
DO i_hf = 1, n_halo_fields

  hf_pos = 3 * (i_hf-1)

! E-W
!-----
  hf_pos = hf_pos + 1
  IF (l_max_not_add( i_hf )) THEN
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims_s%i_start, pdims%i_start - 1
        io_array(i + pdims_s%halo_i,j,i_hf) = &
        MAX( io_array(i + pdims_s%halo_i,j,i_hf), work_array(i,j,hf_pos) )
      END DO
    END DO
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_end + 1, pdims_s%i_end
        io_array(i - pdims_s%halo_i,j,i_hf) = &
        MAX( io_array(i - pdims_s%halo_i,j,i_hf), work_array(i,j,hf_pos) )
      END DO
    END DO
  ELSE
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims_s%i_start, pdims%i_start - 1
        io_array(i + pdims_s%halo_i,j,i_hf) = &
        io_array(i + pdims_s%halo_i,j,i_hf) + work_array(i,j,hf_pos)
      END DO
    END DO
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims%i_end + 1, pdims_s%i_end
        io_array(i - pdims_s%halo_i,j,i_hf) = &
        io_array(i - pdims_s%halo_i,j,i_hf) + work_array(i,j,hf_pos)
      END DO
    END DO
  END IF

! N-S
!-----
  hf_pos = hf_pos + 1
  IF (l_max_not_add( i_hf )) THEN
    DO j = pdims_s%j_start, pdims%j_start - 1
      DO i = pdims%i_start, pdims%i_end
        io_array(i,j + pdims_s%halo_j,i_hf) = &
        MAX( io_array(i,j + pdims_s%halo_j,i_hf), work_array(i,j,hf_pos) )
      END DO
    END DO
    DO j = pdims%j_end + 1, pdims_s%j_end
      DO i = pdims%i_start, pdims%i_end
        io_array(i,j - pdims_s%halo_j,i_hf) = &
        MAX( io_array(i,j - pdims_s%halo_j,i_hf), work_array(i,j,hf_pos) )
      END DO
    END DO
  ELSE
    DO j = pdims_s%j_start, pdims%j_start - 1
      DO i = pdims%i_start, pdims%i_end
        io_array(i,j + pdims_s%halo_j,i_hf) = &
        io_array(i,j + pdims_s%halo_j,i_hf) + work_array(i,j,hf_pos)
      END DO
    END DO
    DO j = pdims%j_end + 1, pdims_s%j_end
      DO i = pdims%i_start, pdims%i_end
        io_array(i,j - pdims_s%halo_j,i_hf) = &
        io_array(i,j - pdims_s%halo_j,i_hf) + work_array(i,j,hf_pos)
      END DO
    END DO
  END IF

! halo corners
!---------------
  hf_pos = hf_pos + 1
  IF (l_max_not_add( i_hf )) THEN
    DO j = pdims_s%j_start, pdims%j_start - 1
      DO i = pdims_s%i_start, pdims%i_start - 1
        io_array(i + pdims_s%halo_i,j + pdims_s%halo_j,i_hf) = &
        MAX( io_array(i + pdims_s%halo_i,j + pdims_s%halo_j,i_hf),  &
        work_array(i,j,hf_pos) )
      END DO
    END DO
    DO j = pdims%j_end + 1, pdims_s%j_end
      DO i = pdims_s%i_start, pdims%i_start - 1
        io_array(i + pdims_s%halo_i,j - pdims_s%halo_j,i_hf) = &
        MAX( io_array(i + pdims_s%halo_i,j - pdims_s%halo_j,i_hf),  &
        work_array(i,j,hf_pos) )
      END DO
    END DO
    DO j = pdims_s%j_start, pdims%j_start - 1
      DO i = pdims%i_end + 1, pdims_s%i_end
        io_array(i - pdims_s%halo_i,j + pdims_s%halo_j,i_hf) = &
        MAX( io_array(i - pdims_s%halo_i,j + pdims_s%halo_j,i_hf),  &
        work_array(i,j,hf_pos) )
      END DO
    END DO
    DO j = pdims%j_end + 1, pdims_s%j_end
      DO i = pdims%i_end + 1, pdims_s%i_end
        io_array(i - pdims_s%halo_i,j - pdims_s%halo_j,i_hf) = &
        MAX( io_array(i - pdims_s%halo_i,j - pdims_s%halo_j,i_hf),  &
        work_array(i,j,hf_pos) )
      END DO
    END DO
  ELSE
    DO j = pdims_s%j_start, pdims%j_start - 1
      DO i = pdims_s%i_start, pdims%i_start - 1
        io_array(i + pdims_s%halo_i,j + pdims_s%halo_j,i_hf) = &
        io_array(i + pdims_s%halo_i,j + pdims_s%halo_j,i_hf) + &
        work_array(i,j,hf_pos)
      END DO
    END DO
    DO j = pdims%j_end + 1, pdims_s%j_end
      DO i = pdims_s%i_start, pdims%i_start - 1
        io_array(i + pdims_s%halo_i,j - pdims_s%halo_j,i_hf) = &
        io_array(i + pdims_s%halo_i,j - pdims_s%halo_j,i_hf) + &
        work_array(i,j,hf_pos)
      END DO
    END DO
    DO j = pdims_s%j_start, pdims%j_start - 1
      DO i = pdims%i_end + 1, pdims_s%i_end
        io_array(i - pdims_s%halo_i,j + pdims_s%halo_j,i_hf) = &
        io_array(i - pdims_s%halo_i,j + pdims_s%halo_j,i_hf) + &
        work_array(i,j,hf_pos)
      END DO
    END DO
    DO j = pdims%j_end + 1, pdims_s%j_end
      DO i = pdims%i_end + 1, pdims_s%i_end
        io_array(i - pdims_s%halo_i,j - pdims_s%halo_j,i_hf) = &
        io_array(i - pdims_s%halo_i,j - pdims_s%halo_j,i_hf) + &
        work_array(i,j,hf_pos)
      END DO
    END DO
  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE edge_exchange
END MODULE edge_exchange_mod
