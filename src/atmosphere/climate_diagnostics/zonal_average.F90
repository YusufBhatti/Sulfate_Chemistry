! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE zonal_average_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ZONAL_AVERAGE_MOD'

CONTAINS


SUBROUTINE zonal_average(field, zonal_field,ini_istart, ini_iend,     &
                                            ini_jstart, ini_jend,     &
                         global_row_length,levs, gc_proc_row_group)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Climate Diagnostics

IMPLICIT NONE

INTEGER :: y, z, levs
INTEGER :: ini_istart, ini_iend
INTEGER :: ini_jstart, ini_jend

INTEGER :: info ! Error variable for RVECSUMR
INTEGER :: global_row_length, gc_proc_row_group

REAL    :: field(ini_istart:ini_iend,ini_jstart:ini_jend,levs)
REAL    :: zonal_field(ini_jstart:ini_jend,levs), lat_result(levs)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ZONAL_AVERAGE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO y=ini_jstart, ini_jend ! latitude loop

  CALL gcg_rvecsumr(ini_iend-ini_istart+1,                      &
                    ini_iend-ini_istart+1, 1, levs,             &
                    field(ini_istart:ini_iend, y, 1:levs),      &
                    gc_proc_row_group, info, lat_result)

  DO z=1,levs
    zonal_field(y,z) = lat_result(z) / global_row_length
  END DO

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE zonal_average

END MODULE zonal_average_mod
