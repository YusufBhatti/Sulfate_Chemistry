! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Gathers a field from many processors to one processor

! Subroutine Interface:
SUBROUTINE gather_field(local_field,    global_field,     &
                        local_row_len,  local_rows,       &
                        global_row_len, global_rows,      &
                        grid_type,      halo_type,        &
                        gather_pe,      proc_group )


USE mpp_conf_mod, ONLY: gcom_coll_limit
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore, ONLY: nproc
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

!
! Description:
! Interface to potentially 2 methods of taking a field that is decomposed
! over a group of processors and gathering the data so a single processor
! contains the entire global field.

! Method:
! For C96_1C GCOM and MPL versions are available. The choice is made
! on a user definable threshold that is given in the gui/namelist.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! Subroutine Arguments:
INTEGER :: local_row_len   ! length of rows in local part of field
INTEGER :: local_rows      ! number of rows in local part of field
INTEGER :: global_row_len  ! length of rows in global field
INTEGER :: global_rows     ! number of rows in global field
INTEGER :: grid_type       ! type (p,u or v) of grid
INTEGER :: halo_type       ! halo type (hence width) of grid
INTEGER :: gather_pe       ! processor to gather global field to
INTEGER :: proc_group      ! group id of processors involved here

REAL    :: local_field(local_row_len*local_rows)
                                         ! local part of field
REAL    :: global_field(global_row_len*global_rows)
                                         ! (on pe gather_pe) global field

! Local variables
INTEGER :: icode      ! for call to sync

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='GATHER_FIELD'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

#if defined(UTILIO) || defined(RECON_SERIAL)

! Only GCOM version available
! DEPENDS ON: gather_field_gcom
CALL gather_field_gcom(local_field,    global_field,       &
                       local_row_len,  local_rows,         &
                       global_row_len, global_rows,        &
                       grid_type,      halo_type,          &
                       gather_pe,      proc_group )


#else

! Use GCOM version if number of processors less than or equal to
! threshold
IF (nproc <= gcom_coll_limit) THEN

  ! DEPENDS ON: gather_field_gcom
  CALL gather_field_gcom(local_field,    global_field,     &
                         local_row_len,  local_rows,       &
                         global_row_len, global_rows,      &
                         grid_type,      halo_type,        &
                         gather_pe,      proc_group )


ELSE

  ! DEPENDS ON: gather_field_mpl
  CALL gather_field_mpl(local_field,    global_field,      &
                        local_row_len,  local_rows,        &
                        global_row_len, global_rows,       &
                        grid_type,      halo_type,         &
                        gather_pe,      proc_group )


END IF
#endif

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE gather_field
