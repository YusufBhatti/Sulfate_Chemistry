! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine U_TO_p for calculating variables held at u points
! at p points. 
!
! This routine does interior points of array not halos,
! but requires halo information to be set.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Grids
MODULE u_to_p_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'U_TO_P_MOD'

CONTAINS

SUBROUTINE u_to_p                                                              &
  ( array_on_u_points, ini_start, ini_end, inj_start, inj_end                  &
  , outi_start, outi_end, outj_start, outj_end, levels, at_extremity           &
  , array_on_p_points )

USE um_parparams,       ONLY: nodomain, pnorth, peast, psouth, pwest

USE model_domain_mod,   ONLY: model_type, mt_single_column
USE compute_chunk_size_mod, ONLY: compute_chunk_size

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) :: ini_start, ini_end
INTEGER, INTENT(IN) :: inj_start, inj_end
INTEGER, INTENT(IN) :: outi_start, outi_end
INTEGER, INTENT(IN) :: outj_start, outj_end
INTEGER, INTENT(IN) :: levels


LOGICAL, INTENT(IN) :: at_extremity(4) ! Indicates if this processor is at
                                       ! north, south, east or west of the
                                       ! processor grid

REAL, INTENT(IN)  :: array_on_u_points( ini_start:ini_end                      &
                                      , inj_start:inj_end                      &
                                      , levels )

REAL, INTENT(OUT) :: array_on_p_points( outi_start:outi_end                    &
                                      , outj_start:outj_end                    &
                                      , levels )

! local variables
INTEGER :: i, j, k
INTEGER :: j0, j1  
INTEGER :: ompt_start   ! start of thread's work
INTEGER :: ompt_end     ! end of thread's work

! Parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'U_TO_P'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! this procedure is now thread safe. i.e. it can be called from an
! OMP parallel region. Compute_chunk_size will determine the start and end
! iterations for each thread as if it were a static schedule.
! If called from a serial region, start and end are 1 and levels

ompt_start = 1
ompt_end = levels
! only call compute_chunk_size if compiling with OMP
! The procedure call is protected by the optional compile
! sentinel
!$ CALL compute_chunk_size(1,levels,ompt_start,ompt_end)

j0 = outj_start
j1 = outj_end


SELECT CASE (model_type)

CASE (mt_single_column)
  DO k=ompt_start, ompt_end
    DO j=j0, j1
      DO i=outi_start, outi_end
        array_on_p_points(i,j,k) = array_on_u_points(i,j,k)
      END DO
    END DO
  END DO

CASE DEFAULT
  DO k=ompt_start, ompt_end
    DO j=j0, j1
      DO i=outi_start, outi_end
        array_on_p_points(i,j,k) = 0.5 * ( array_on_u_points(i,j,k)            &
                                         + array_on_u_points(i-1,j,k) )
      END DO
    END DO
  END DO

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE u_to_p

END MODULE u_to_p_mod
