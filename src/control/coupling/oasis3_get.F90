! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis3_get(data_64,var_id32,l_cols,l_rows,oasis_info,ft &
   ,r_cols,r_rows,get_step)

  !
  !  Description: Interface to OASIS3-MCT "get" operation.
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: Coupling
  !===========================================================

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
      PrintStatus,          &
      PrStatus_Diag

USE mod_prism, ONLY: prism_get_proto
USE um_parvars
USE um_types

USE oasis_atm_data_mod, ONLY: prism_nsec
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

REAL (KIND=real64), INTENT(OUT) :: data_64(l_cols,l_rows,1) ! The received
                                                            ! transient

INTEGER (KIND=integer64), INTENT(IN) :: l_cols, l_rows  ! Local domain size

INTEGER (KIND=integer32), INTENT(OUT) :: oasis_info ! OASIS return code

INTEGER (KIND=integer64), INTENT(IN) :: ft          ! Field type

INTEGER (KIND=integer64), INTENT(IN) :: r_cols, r_rows  ! Receive buffer size

LOGICAL (KIND=logical64), INTENT(IN) :: get_step ! Exchange step flag

INTEGER (KIND=integer32), INTENT(IN) :: var_id32

! Local variables
REAL (KIND=real64) :: data_recv(r_cols,r_rows) ! Incoming data

INTEGER (KIND=integer32) :: prism_nsec32

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='OASIS3_GET'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

oasis_info = 0

! We only actually do the gets and related processing on a genuine
! exchange step for performance reasons.
IF (get_step) THEN

  ! We ensure the incoming data array is initially zero in order
  ! to ensure OASIS generated statistics make sense (i.e. do not
  ! include values from uninitialised areas of memory.)
  data_recv(:,:) = 0.0

  prism_nsec32=prism_nsec

  ! GET the incoming coupling field
  CALL prism_get_proto(var_id32,prism_nsec32,data_recv,oasis_info)

  IF ( PrintStatus >= PrStatus_Diag ) THEN
    IF (oasis_info /= 0) THEN
      WRITE(umMessage,'(A,1X,I6)') "Return code from prism_get_proto",oasis_info
      CALL umPrint(umMessage,src='oasis3_get')
      WRITE(umMessage,'(A,1X,I6,1X,I8)') "Field: ", var_id32,prism_nsec
      CALL umPrint(umMessage,src='oasis3_get')
    END IF
  END IF

  ! Nothing to do here beyond copying the receive buffer
  ! straight to our local domain.
  data_64(1:l_cols,1:l_rows,1)=data_recv(1:r_cols,1:r_rows)

END IF  ! get_step = true

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE oasis3_get
