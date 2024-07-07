! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis3_put(data_64,var_id32,l_cols,l_rows,oasis_info,ft           &
     ,s_cols,s_rows,put_step_ao,l_second_order)
  !
  !  Description: Interface to OASIS3/OASIS3-MCT "put" operation.
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: Coupling
  !
  !===========================================================

USE umPrintMgr, ONLY:       &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Diag

USE mod_prism, ONLY: prism_put_proto

USE oasis_grad_calc_mod, ONLY: oasis_grad_calc

USE um_parvars
USE um_types

USE oasis_atm_data_mod, ONLY: prism_nsec
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

REAL (KIND=real64) :: data_64(l_cols,l_rows) ! The transient
REAL (KIND=real64) :: grad_ew_64(l_cols,l_rows) ! EW gradient transient
REAL (KIND=real64) :: grad_ns_64(l_cols,l_rows) ! NS gradient transient
                                                         ! to send
INTEGER (KIND=integer64), INTENT(IN) :: l_cols, l_rows     ! Loc. domain sizes

INTEGER (KIND=integer32), INTENT(OUT) :: oasis_info        ! OASIS return code

INTEGER (KIND=integer64), INTENT(IN) :: ft                 ! Field type

INTEGER (KIND=integer64), INTENT(IN) :: s_cols, s_rows     ! Send buffer sizes

LOGICAL (KIND=logical64), INTENT(IN) :: put_step_ao       ! Ocean coupling step

LOGICAL (KIND=logical32), INTENT(INOUT) :: l_second_order ! 2nd order coupling

INTEGER (KIND=integer32), INTENT(IN) :: var_id32

! Local variables
REAL (KIND=real64) :: data_send(s_cols,s_rows) ! Outgoing data buffer
REAL (KIND=real64) :: grad_ew_send(s_cols,s_rows) ! Outgoing data buffer
REAL (KIND=real64) :: grad_ns_send(s_cols,s_rows) ! Outgoing data buffer

INTEGER (KIND=integer32) :: prism_nsec32

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='OASIS3_PUT'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

IF (put_step_ao) THEN

  ! If this field is using 2nd order conservative regridding and
  ! OASIS3-MCT, then we need to calculate the gradients to pass to
  ! OASIS3-MCT

  IF (l_second_order) THEN

    CALL oasis_grad_calc(data_64,grad_ew_64,grad_ns_64)

  END IF

  ! Move the data to send direct to our send buffer.
  data_send(1:l_cols,1:l_rows)=data_64(1:s_cols,1:s_rows)

  prism_nsec32 = prism_nsec

  IF (l_second_order) THEN

    grad_ew_send(1:l_cols,1:l_rows)=grad_ew_64(1:s_cols,1:s_rows)
    grad_ns_send(1:l_cols,1:l_rows)=grad_ns_64(1:s_cols,1:s_rows)

    ! PUT the outgoing coupling field with gradients. Note that 
    ! it's important that the 2nd term relates to the gradient in the 
    ! latitude direction and the 3rd term to the gradient in the longitude 
    ! direction for consistency with SCRIP formatting of 2nd order weights 
    ! files. In principle we could use a 4th and 5th term for diagonal
    ! gradients, but in practice we never employ weights files which 
    ! cater for those terms. 
    CALL prism_put_proto(var_id32,prism_nsec32,data_send,  &
                                               oasis_info,   &
                                               grad_ns_send, &
                                               grad_ew_send)
  ELSE

    ! PUT the outgoing coupling field
    CALL prism_put_proto(var_id32,prism_nsec32,data_send,oasis_info)
  END IF

  IF ( PrintStatus >= PrStatus_Diag ) THEN
    IF (oasis_info /= 0) THEN
      WRITE(umMessage,'(A,1X,I6)') "Return code from prism_put_proto: ", &
                                      oasis_info
      CALL umPrint(umMessage,src='oasis3_put')
      WRITE(umMessage,'(A,1X,I6,1X,I8)') "Field: ", var_id32,prism_nsec
      CALL umPrint(umMessage,src='oasis3_put')
    END IF
  END IF


END IF  ! put_step_ao = true.

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE oasis3_put
