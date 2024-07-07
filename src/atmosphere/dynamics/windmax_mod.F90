! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE windmax_mod

USE um_parcore,     ONLY: mype, nproc
USE atm_fields_bounds_mod
USE missing_data_mod, ONLY: rmdi
USE dynamics_input_mod, ONLY: not_conserved

USE filenamelength_mod, ONLY: filenamelength
USE file_manager, ONLY: assign_file_unit, release_file_unit
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

INTEGER, PUBLIC :: eg_unit
CHARACTER(LEN=filenamelength), PUBLIC, PARAMETER :: windmax_file="windmax.dat"
PUBLIC :: print_windmax_1, print_windmax_2

PRIVATE


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='WINDMAX_MOD'

CONTAINS

SUBROUTINE print_windmax_1(u,v,w,conserve_dry_mass)

IMPLICIT NONE
!
! Description:
!
!
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

INTEGER, INTENT(IN) :: conserve_dry_mass ! Method for conserving dry mass

REAL, INTENT(IN) :: u(udims_s%i_start:udims_s%i_end,       &
   udims_s%j_start:udims_s%j_end,udims_s%k_start:udims_s%k_end)
REAL, INTENT(IN) :: v(vdims_s%i_start:vdims_s%i_end,       &
   vdims_s%j_start:vdims_s%j_end,vdims_s%k_start:vdims_s%k_end)
REAL, INTENT(IN) :: w(wdims_s%i_start:wdims_s%i_end,       &
   wdims_s%j_start:wdims_s%j_end,wdims_s%k_start:wdims_s%k_end)

REAL :: vels(3)
INTEGER :: istat,i,j,k

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_WINDMAX_1'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

vels=RMDI

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)  &
!$OMP& SHARED(vdims,u,v,w,udims,wdims) &
!$OMP& REDUCTION(MAX:vels)

!$OMP DO SCHEDULE(STATIC)
DO k=udims%k_start,udims%k_end
  DO j=udims%j_start,udims%j_end
    DO i=udims%i_start,udims%i_end
      vels(1)=MAX(vels(1),ABS(u(i,j,k)))
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k=vdims%k_start,vdims%k_end
  DO j=vdims%j_start,vdims%j_end
    DO i=vdims%i_start,vdims%i_end
      vels(2)=MAX(vels(2),ABS(v(i,j,k)))
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k=wdims%k_start,wdims%k_end
  DO j=wdims%j_start,wdims%j_end
    DO i=wdims%i_start,wdims%i_end
      vels(3)=MAX(vels(3),ABS(w(i,j,k)))
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

! the gcom call below calculates the maximum of each element of vels
! across all processors and sends the result to pe0 only
CALL gc_rmax_single_task(3, nproc, istat, vels, 0)

IF ( mype == 0 ) THEN
  IF (conserve_dry_mass == not_conserved) THEN
    WRITE(UNIT=eg_unit,FMT='(I7,3E25.16)') 0, vels(1), vels(2), vels(3)
  ELSE
    WRITE(UNIT=eg_unit,FMT='(I7,4E25.16)') 0, vels(1), vels(2), vels(3), 1.0
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

END SUBROUTINE


SUBROUTINE print_windmax_2(u,v,w,conserve_dry_mass,mass_fix_factor)

USE timestep_mod,   ONLY: timestep_number

IMPLICIT NONE

INTEGER, INTENT(IN) :: conserve_dry_mass

REAL, INTENT(IN) :: u(udims_s%i_start:udims_s%i_end,       &
   udims_s%j_start:udims_s%j_end,udims_s%k_start:udims_s%k_end)
REAL, INTENT(IN) :: v(vdims_s%i_start:vdims_s%i_end,       &
   vdims_s%j_start:vdims_s%j_end,vdims_s%k_start:vdims_s%k_end)
REAL, INTENT(IN) :: w(wdims_s%i_start:wdims_s%i_end,       &
   wdims_s%j_start:wdims_s%j_end,wdims_s%k_start:wdims_s%k_end)

REAL, INTENT(IN) :: mass_fix_factor

REAL :: vels(3)
INTEGER :: istat,i,j,k

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_WINDMAX_2'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

vels=RMDI

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k)  &
!$OMP& SHARED(vdims,u,v,w,udims,wdims)  &
!$OMP& REDUCTION(MAX:vels)

!$OMP DO SCHEDULE(STATIC)
DO k=udims%k_start,udims%k_end
  DO j=udims%j_start,udims%j_end
    DO i=udims%i_start,udims%i_end
      vels(1)=MAX(vels(1),ABS(u(i,j,k)))
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k=vdims%k_start,vdims%k_end
  DO j=vdims%j_start,vdims%j_end
    DO i=vdims%i_start,vdims%i_end
      vels(2)=MAX(vels(2),ABS(v(i,j,k)))
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k=wdims%k_start,wdims%k_end
  DO j=wdims%j_start,wdims%j_end
    DO i=wdims%i_start,wdims%i_end
      vels(3)=MAX(vels(3),ABS(w(i,j,k)))
    END DO
  END DO
END DO
!$OMP END DO 
!$OMP END PARALLEL

! the gcom call below calculates the maximum of each element of vels
! across all processors and sends the result to pe0 only
CALL gc_rmax_single_task(3, nproc, istat, vels, 0)

IF ( mype == 0 ) THEN
  IF (conserve_dry_mass == not_conserved) THEN

    WRITE(UNIT=eg_unit,FMT='(I7,3E25.16)') timestep_number,           &
                                           vels(1), vels(2), vels(3)
  ELSE
    WRITE(UNIT=eg_unit,FMT='(I7,4E25.16)') timestep_number,           &
                                           vels(1), vels(2), vels(3), &
                                           mass_fix_factor
  END IF
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)
END SUBROUTINE
END MODULE
