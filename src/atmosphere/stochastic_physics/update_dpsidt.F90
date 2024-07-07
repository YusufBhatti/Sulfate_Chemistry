! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stochastic Physics
MODULE update_dpsidt_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UPDATE_DPSIDT_MOD'

CONTAINS


SUBROUTINE update_dpsidt(alpha, dpsidtc, dpsidts, stph_n1, stph_n2,     &
                         nlim, gspect, zero, icode, info, firsttimestep)

! This routine updates the spherical harmonic coefficients using
! a first-order auto-regressive process. gspect is the wave-number
! dependent power spectrum generated at the start of the run in
! subroutine "backscatter_spectrum"

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE stochastic_physics_run_mod, ONLY:                                   &
    firsttimestep_true, firsttimestep_crun, stphseed,                   &
    stph_skeb2_data_present, stph_skeb2_data_check
USE nlstcall_nrun_as_crun_mod, ONLY: l_nrun_as_crun
USE stph_seed_mod, ONLY: stph_seed_sync
USE umprintmgr, ONLY: umprint, printstatus, prstatus_normal
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: stph_n1, stph_n2, nlim, icode, info, zero, firsttimestep
REAL, INTENT(IN) :: alpha, gspect(nlim)
REAL, INTENT(INOUT) :: dpsidtc(0:stph_n2,stph_n1:stph_n2)
REAL, INTENT(INOUT) :: dpsidts(0:stph_n2,stph_n1:stph_n2)

!     Local variables
REAL, SAVE :: oneminusalpha, squarerootalpha, sdx2
INTEGER :: m, n

! Random numbers on each processor
REAL :: rand_nums(0:stph_n2,stph_n1:stph_n2)
REAL :: rand_nums2(0:stph_n2,stph_n1:stph_n2)

! Local error reporting
INTEGER :: err
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_DPSIDT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! It is important that rand_nums and rand_nums2 have the same random
! number between all the PES. Thus we retrieve the seed values from
! PE zero, then broadcast the values to all the PEs which then set
! the new seed values. 
CALL stph_seed_sync(zero)

! Now all the PEs should be at the same point in the random sequence.
! Generate a separate set of random numbers for sin and cosine coeffs
CALL RANDOM_NUMBER(rand_nums)
CALL RANDOM_NUMBER(rand_nums2)

! First call in both NRUN and CRUN
IF (firsttimestep == firsttimestep_true .OR.                            &
    firsttimestep == firsttimestep_crun) THEN

  ! Set up AR1-process parameters and save
  oneminusalpha=1-alpha
  ! SQRT(alpha) ensures that the noise (random nums) decorrelates
  ! faster than the timestep
  squarerootalpha=SQRT(alpha)
  ! Twice standard deviation of random numbers
  ! These have a variance of (1/12)
  sdx2 = 2.0 * SQRT(0.083333333333333)

END IF

! If this is the first time step, and we are not reading the seed from
! a file, then we initialise dpsidtc and dpsidts with random values
! If it is not the first timestep, or it is the first timestep but we have 
! read in coefficient values from a file, then update the existing values 
! with random perturbations
IF ( firsttimestep == firsttimestep_true .AND.                         &
    l_nrun_as_crun .AND.                                               &
   (stph_skeb2_data_check == stph_skeb2_data_present) )  THEN
  cmessage = 'update_dpsidt: l_nrun_as_crun is TRUE:' //  &
    'Updating based on existing values'
  err = -1
  CALL ereport(RoutineName, err, cmessage)
END IF

IF ( firsttimestep == firsttimestep_true .AND.                         &
   .NOT. ((stphseed == 2 .OR. l_nrun_as_crun) .AND.                    &
   (stph_skeb2_data_check == stph_skeb2_data_present) ) ) THEN
  IF (printstatus >= prstatus_normal) THEN
    CALL umPrint('update_dpsidt: updating for first nrun', src='update_dpsidt')
  END IF

  !  Initialise with random value around expected mean
  DO n = stph_n1, stph_n2
    DO m = 0, n
      dpsidtc(m,n) = sdx2 * (rand_nums(m,n)-0.5)*gspect(n)
      dpsidts(m,n) = sdx2 * (rand_nums2(m,n)-0.5)*gspect(n)
    END DO
  END DO
ELSE
  !  Continues with AR1 update using existing values of dpsidtc/s
  !  (RANDOM-0.5) provides white noise with mean of 0 and var. of (1/12)
  IF (printstatus >= prstatus_normal) THEN
    CALL umPrint('update_dpsidt: updating based on existing values', &
                 src='update_dpsidt')
  END IF
  
  DO n = stph_n1, stph_n2
    DO m = 0, n
      dpsidtc(m,n) = oneminusalpha*dpsidtc(m,n) + squarerootalpha       &
                     *(rand_nums(m,n)-0.5)*gspect(n)
      dpsidts(m,n) = oneminusalpha*dpsidts(m,n) + squarerootalpha       &
                     *(rand_nums2(m,n)-0.5)*gspect(n)
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE update_dpsidt

END MODULE update_dpsidt_mod
