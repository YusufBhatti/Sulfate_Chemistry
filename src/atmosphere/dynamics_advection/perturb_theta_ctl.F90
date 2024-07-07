! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine perturb_theta_ctl
SUBROUTINE perturb_theta_ctl(                                     &
                         row_length_in, rows_in, model_levels_in, &
                         global_row_length_in, global_rows_in,    &
                         at_extremity, offx_in, offy_in,          &
                         IntRand_Seed, l_datastart )

! Purpose:
!          Perturb both theta and thetavd at the bit level using 
!          a (common) set of random numbers
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE atm_fields_mod, ONLY: theta, thetav

IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.
INTEGER, INTENT(IN) :: ROW_LENGTH_in     ! No of points per local row
INTEGER, INTENT(IN) :: ROWS_in           ! No of local (theta) rows
INTEGER, INTENT(IN) :: global_ROW_LENGTH_in
INTEGER, INTENT(IN) :: global_ROWS_in
INTEGER, INTENT(IN) :: MODEL_LEVELS_in   ! No of model levels
INTEGER, INTENT(IN) :: Offx_in    ! standard halo size in East-West
INTEGER, INTENT(IN) :: Offy_in    ! standard halo size in North-South
LOGICAL, INTENT(IN) :: At_extremity(4)
INTEGER :: IntRand_Seed
INTEGER, INTENT(IN) :: l_datastart(2)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PERTURB_THETA_CTL'


!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
! DEPENDS ON: perturb_theta
CALL perturb_theta( theta, thetav,                                &
                         row_length_in, rows_in, model_levels_in, &
                         global_row_length_in, global_rows_in,    &
                         at_extremity, offx_in, offy_in,          &
                         IntRand_Seed, l_datastart )

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE perturb_theta_ctl
