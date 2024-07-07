! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Decode the STASH level code

! *****************************************************************
! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
! Any Changes to this routine must be accompanied with equivalent
! changes to rcf_level_code_mod.F90 and rcf_address_mod.F90
! *****************************************************************

! Subroutine Interface:
SUBROUTINE levcod(ilin,ilout,ErrorStatus,cmessage)

USE rad_input_mod, ONLY: h_swbands,h_lwbands

USE model_domain_mod, ONLY: output_grid_stagger,                    &
                            FH_GridStagger_C

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY:                                     &
    bl_levels, cloud_levels, model_levels, ozone_levels, sm_levels, &
    st_levels, tr_levels
USE errormessagelength_mod, ONLY: errormessagelength

USE stash_model_mod, ONLY: stlevgwdrag, botvdifflev, topvdifflev

IMPLICIT NONE

! Description:
!   Sets ILOUT to an appropriate level size according to the value of IL
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 90
!    Written to UMDP3 programming standards version 8.3.
!

! Subroutine arguments:
!   Scalar arguments with intent(in):
INTEGER :: ilin    ! Model level code

!   Scalar arguments with intent(out):
INTEGER :: ilout   ! An actual level
CHARACTER(LEN=errormessagelength) :: cmessage

! Local scalars:
INTEGER :: i
INTEGER :: j

! ErrorStatus:
INTEGER :: ErrorStatus

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LEVCOD'

!- End of Header ---------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (ilin == 1) THEN
  ! First atmos level
  ilout=1
ELSE IF (ilin == 2) THEN
  ! Top atmos level
  ilout=model_levels
ELSE IF (ilin == 3) THEN
  ! Top wet level
  ilout=model_levels
ELSE IF (ilin == 4) THEN
  ilout=model_levels-1
ELSE IF (ilin == 5) THEN
  ! First boundary layer level
  ilout=MIN(1,bl_levels)
ELSE IF (ilin == 6) THEN
  ! Last boundary layer level
  ilout=bl_levels
ELSE IF (ilin == 7) THEN
  ilout=bl_levels+1
ELSE IF (ilin == 8) THEN
  ! First soil level
  ilout=MIN(1,st_levels)
ELSE IF (ilin == 9) THEN
  ! Last soil level
  ilout=st_levels
ELSE IF (ilin == 10) THEN
  ! First tracer level
  ilout=model_levels-tr_levels+1
ELSE IF (ilin == 11) THEN
  ! Last tracer level
  ilout=model_levels
ELSE IF (ilin == 12) THEN
  ilout=model_levels+1
ELSE IF (ilin == 13) THEN
  ! First gravity wave drag level
  ilout=StLevGWdrag
ELSE IF (ilin == 14) THEN
  ! Last gravity wave drag level
  ilout=model_levels
ELSE IF (ilin == 15) THEN
  ilout=BotVDiffLev
ELSE IF (ilin == 16) THEN
  ilout=TopVDiffLev-1
ELSE IF (ilin == 17) THEN
  ilout=TopVDiffLev
ELSE IF (ilin == 18) THEN
  ilout=bl_levels-1
ELSE IF (ilin == 19) THEN
  ilout=model_levels+1
ELSE IF (ilin == 20) THEN
  ilout=MIN(2,st_levels)
ELSE IF (ilin == 21) THEN
  ilout=1
ELSE IF (ilin == 22) THEN
  !       ILOUT=NLEVSO - ocean levels no longer supported
  ilout=1
ELSE IF (ilin == 23) THEN
  ilout=ozone_levels
ELSE IF (ilin == 24) THEN
  ilout=model_levels*h_swbands
ELSE IF (ilin == 25) THEN
  ilout=(model_levels+1)*h_swbands
ELSE IF (ilin == 26) THEN
  ilout=model_levels*h_swbands
ELSE IF (ilin == 27) THEN
  ilout=model_levels*h_lwbands
ELSE IF (ilin == 28) THEN
  ilout=(model_levels+1)*h_lwbands
ELSE IF (ilin == 29) THEN
  ilout=model_levels*h_lwbands
ELSE IF (ilin == 30) THEN
  ilout=2
ELSE IF (ilin == 32) THEN
  ilout=h_swbands
ELSE IF (ilin == 33) THEN
  ilout=h_lwbands
ELSE IF (ilin == 34) THEN
  ilout=sm_levels
ELSE IF (ilin == 35) THEN
  ilout=cloud_levels
ELSE IF (ilin == 38) THEN
  ! Surface level on Charney-Phillips theta grid
  ilout=0
ELSE IF (ilin == 39) THEN
  ! Number of ISCCP simulator levels
  ilout=7
ELSE IF (ilin == 40) THEN
  ! Fields which are different between ENDGAME and New Dynamics regarding the
  ! whether it starts at theta level 0 in ENDGAME but 1 in New Dynamics.
  IF (output_grid_stagger == FH_GridStagger_C) THEN
    ilout=1
  ELSE
    ilout=0
  END IF

ELSE
  WRITE(umMessage,'(A35,I6)')'LEVCOD: IMPOSSIBLE LEVEL CODE FOUND ',ilin
  CALL umPrint(umMessage,src='levcod')
  WRITE(cmessage,'(A35,I6)')'LEVCOD: IMPOSSIBLE LEVEL CODE FOUND ',ilin
  ErrorStatus=1
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE levcod
