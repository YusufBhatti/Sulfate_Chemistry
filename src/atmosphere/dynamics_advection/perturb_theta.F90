! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine perturb_theta
SUBROUTINE perturb_theta(theta, thetav,                           &
                         row_length, rows, model_levels,          &
                         global_row_length, global_rows,          &
                         at_extremity, offx, offy, IntRand_Seed,  &
                         l_datastart )

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
USE UM_ParParams
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE


! Arguments with Intent IN. ie: Input variables.
INTEGER, INTENT(IN) :: row_length     ! No of points per local row
INTEGER, INTENT(IN) :: rows           ! No of local (theta) rows
INTEGER, INTENT(IN) :: model_levels   ! No of model levels
INTEGER, INTENT(IN) :: Offx    ! standard halo size in East-West
INTEGER, INTENT(IN) :: Offy    ! standard halo size in North-South
LOGICAL, INTENT(IN) :: At_extremity(4)
INTEGER, INTENT(IN) :: IntRand_Seed   ! Seed for random numbers
INTEGER, INTENT(IN) :: l_datastart(2)
INTEGER, INTENT(IN) :: global_ROW_LENGTH
INTEGER, INTENT(IN) :: global_ROWs

REAL, INTENT (INOUT) ::                                           &
  theta(1-offx:row_length+offx, 1-offy:rows+offy,                 &
        0:model_levels)

REAL, INTENT (INOUT) ::                                           &
  thetav(1-offx:row_length+offx, 1-offy:rows+offy,                &
        0:model_levels)

! Local Variables
INTEGER :: i,j,k    ! loop variables
INTEGER :: seed_len ! Length of random number seed vector
INTEGER :: j_start,j_end
INTEGER :: gi,gj

INTEGER, ALLOCATABLE :: seed(:) ! Vector random number seed

REAL :: RandomNumbers(global_row_length,global_rows,model_levels)
REAL :: eps_machine
REAL :: factor
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PERTURB_THETA'



!----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
IF (model_type == mt_global .AND. at_extremity(Psouth)) THEN
  j_start=2
ELSE
  j_start=1
END IF
IF (model_type == mt_global .AND. at_extremity(PNorth)) THEN
  j_end=rows-1
ELSE
  j_end=rows
END IF

WRITE(umMessage,*)'PERTURBING THETA FIELD'
CALL umPrint(umMessage,src='perturb_theta')
eps_machine=EPSILON(1.0)
WRITE(umMessage,*)' MACHINE epsilon: ', eps_machine
CALL umPrint(umMessage,src='perturb_theta')
WRITE(umMessage,*)' IntRand_Seed:    ', IntRand_Seed
CALL umPrint(umMessage,src='perturb_theta')

! Initialise random number seed.
! (For some reason, the seeds n*4, n*4+1, n*4+2 and n*4+3 all
! give the same random numbers on the IBM. To get around this,
! we multiply the supplied seed by 4.)
CALL RANDOM_SEED (SIZE = seed_len)
ALLOCATE (seed(seed_len))
seed(:) = IntRand_seed * 4
CALL RANDOM_SEED (put = seed)
DEALLOCATE (seed)

! Get random numbers in the range 0-1:
CALL RANDOM_NUMBER (RandomNumbers)

DO k=1,model_levels
  DO j=j_start,j_end
    gj=j+l_datastart(2)-1
    DO i=1, row_length
      gi=i+l_datastart(1)-1
      factor = (2.0*RandomNumbers(gi,gj,k) - 1.0) * eps_machine
      theta(i,j,k)  = theta(i,j,k)  + theta(i,j,k)  * factor
      thetav(i,j,k) = thetav(i,j,k) + thetav(i,j,k) * factor  
    END DO
  END DO
END DO

! Copy level 1 back into ENDGame level 0
DO j=j_start,j_end
  gj=j+l_datastart(2)-1
  DO i=1, row_length
    gi=i+l_datastart(1)-1
    theta(i,j,0)  = theta(i,j,1)
    thetav(i,j,0) = thetav(i,j,1)
  END DO
END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE perturb_theta
