! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE um_submodel_init_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UM_SUBMODEL_INIT_MOD'

CONTAINS
!
! A routine to initialise model for submodel and internal model coupling
!
! Subroutine Interface:
SUBROUTINE UM_Submodel_Init(ErrorStatus)

USE d1_array_mod, ONLY: alt_n_submodel_partition, alt_n_submodel_partition_max
USE submodel_mod, ONLY:                                                        &
    internal_model_list, atmos_sm, submodel_for_sm, n_internal_for_sm,         &
    submodel_for_im,                                                           &
    submodel_partition_index, submodel_partition_list,                         &
    n_submodel_partition_max, atmos_im
USE check_iostat_mod
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

IMPLICIT NONE
!
! Description:
!   UM_Submodel_Init initialises the model with information specifying
!   internal model and submodel partitions for the run,
!   We can hard-wire these as we only have atmos.
!   The code is kept to support stash proc arrays.
!
! Method:
!   The routine reads information from the user interface, providing
!   lists of internal models and their associated submodel data
!   partitions. This is required in both the reconfiguration and the
!   model as a prior step to calculating addressing in STASH_PROC.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered:
! System Task:
!
! Subroutine arguments

!   ErrorStatus
INTEGER ::   ErrorStatus          ! Error flag (0 = OK)

! Local parameters:

! Local scalars:
INTEGER ::                                                        &
 s                                                                &
                 ! submodel loop
,i                                                                &
                 ! internal model loop
,sm                                                               &
                 ! submodel identifier
,im                                                               &
                 ! internal model identifier
,sm_prev                                                          &
                 ! previous submodel identifier
,im_prev         ! previous internal model identifier

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Local dynamic arrays:

! Function & Subroutine calls: None

!- End of header

! IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Only Atmosphere Model catered for so no need to read in NSUBMODL.
internal_model_list(1) = atmos_im
submodel_for_im(1)     = atmos_sm

im = atmos_im
sm = atmos_sm
submodel_partition_list(1) = sm
submodel_for_sm(im) = 1
submodel_partition_index(1)=sm
n_internal_for_sm(sm)=1

! Need a copy of No of submodels for use by d1_array_mod.
alt_n_submodel_partition=1

IF (alt_n_submodel_partition_max /= n_submodel_partition_max) THEN
  WRITE(umMessage,*)'UM_Submodel_In: Mismatch in parameters '
  CALL umPrint(umMessage,src='um_submodel_init')
  WRITE(umMessage,*)'N_SUBMODEL_PARTITION_MAX and '
  CALL umPrint(umMessage,src='um_submodel_init')
  WRITE(umMessage,*)'ALT_N_SUBMODEL_PARTITION_MAX. '
  CALL umPrint(umMessage,src='um_submodel_init')
  WRITE(umMessage,*)'They should be identical '
  CALL umPrint(umMessage,src='um_submodel_init')
  ErrorStatus=1
END IF
! IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE UM_Submodel_Init
END MODULE um_submodel_init_mod
