! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine P_TO_V  for calculating variables held at p points
! at v points. 
!
! This routine does interior points of array not halos,
! but requires halo information to be set.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Grids
!
! global code has E-W wrap around
! LAM code has most east u point set to zero
MODULE p_to_v_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'P_TO_V_MOD'

CONTAINS

SUBROUTINE p_to_v                                                              &
  ( array_on_p_points, ini_start,ini_end, inj_start,inj_end                    &
  , outi_start,outi_end, outj_start,outj_end, outk_start,outk_end              &
  , array_on_v_points )

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE compute_chunk_size_mod, ONLY: compute_chunk_size

USE model_domain_mod, ONLY: model_type, mt_single_column

IMPLICIT NONE

INTEGER, INTENT(IN) :: ini_start,ini_end
INTEGER, INTENT(IN) :: inj_start,inj_end
INTEGER, INTENT(IN) :: outi_start,outi_end
INTEGER, INTENT(IN) :: outj_start,outj_end
INTEGER, INTENT(IN) :: outk_start,outk_end

REAL, INTENT(IN)  :: array_on_p_points( ini_start:ini_end                      &
                                      , inj_start:inj_end                      &
                                      , outk_start:outk_end )

REAL, INTENT(OUT) :: array_on_v_points( outi_start:outi_end                    &
                                      , outj_start:outj_end                    &
                                      , outk_start:outk_end )

! Local variables
INTEGER :: i, j, k
INTEGER  ::   ompt_start   ! start of thread's work
INTEGER  ::   ompt_end     ! end of thread's work

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'P_TO_V'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! this procedure is now thread safe. i.e. it can be called from an
! OMP parallel region. Compute_chunk_size will determine the start and end
! iterations for each thread as if it were a static schedule.
! If called from a serial region, start and end are 1 and levels

ompt_start = outk_start
ompt_end = outk_end
! only call compute_chunk_size if compiling with OMP
! The procedure call is protected by the optional compile
! sentinel
!$ CALL compute_chunk_size(outk_start,outk_end,ompt_start,ompt_end)

SELECT CASE (model_type)

CASE (mt_single_column)

  DO k=ompt_start, ompt_end
    DO j=outj_start, outj_end
      DO i=outi_start, outi_end
        array_on_v_points(i,j,k) = array_on_p_points(i,j,k)
      END DO
    END DO
  END DO

CASE DEFAULT

  DO k=ompt_start, ompt_end
    DO j=outj_start, outj_end
      DO i=outi_start, outi_end
        array_on_v_points(i,j,k) = 0.5 * ( array_on_p_points(i,j,k)            &
                                         + array_on_p_points(i,j+1,k) )
      END DO
    END DO
  END DO

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE p_to_v

END MODULE p_to_v_mod
