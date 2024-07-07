! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Description: This routine updates the spherical harmonic coefficients
!               using a first-order auto-regressive process.
!               pspect is the wavenumber-dependent power spectrum
!               generated at the start of the run in
!               subroutine "pattern_spectrum2"
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stochastic Physics
MODULE update_pattern_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UPDATE_PATTERN_MOD'

CONTAINS

SUBROUTINE update_pattern(alpha, coeffc, coeffs, nlim,                  &
                          pspect, zero, icode, info)


USE stochastic_physics_run_mod, ONLY: stph_n2

USE atm_step_local,      ONLY: first_atmstep_call

USE UM_ParVars
USE UM_ParCore, ONLY: mype, nproc
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY: umPrint, umMessage


IMPLICIT NONE

INTEGER, INTENT(IN) :: nlim, icode, info, zero
REAL, INTENT(IN) :: alpha, pspect(nlim)
REAL, INTENT(INOUT) :: coeffc(0:stph_n2,1:stph_n2)
REAL, INTENT(INOUT) :: coeffs(0:stph_n2,1:stph_n2)
!     Local variables
REAL, SAVE :: oneminusalpha, squarerootalpha
INTEGER :: m, n

! Random seeds 
INTEGER, ALLOCATABLE :: iranseed(:)
! Size of the random seed array
INTEGER :: seedcount      

! Random numbers on each processor
REAL :: rand_nums(0:stph_n2,1:stph_n2), rand_nums2(0:stph_n2,1:stph_n2)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_PATTERN'
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! It is important that rand_nums and rand_nums2 have the same random
! number between all the PES. Thus we retrieve the seed values from
! PE zero, then broadcast the values to all the PEs which then set
! the new seed values. 

! Get size of random seed and allocate array
CALL RANDOM_SEED(SIZE=seedcount)
IF (.NOT. ALLOCATED(iranseed)) ALLOCATE(iranseed(seedcount))

! Get the random seed values from processor zero and broadcast them to
! all the PEs.
IF ( mype == zero ) THEN
  CALL RANDOM_SEED(GET=iranseed(1:seedcount))
END IF

CALL gc_ibcast(3245, seedcount, zero, nproc, icode, iranseed)

! Set the random seed for all PEs apart from the zero one. 
! and deallocate the array
IF (mype /= 0) THEN
  CALL RANDOM_SEED(PUT=iranseed(1:seedcount))
END IF

IF (ALLOCATED(iranseed)) DEALLOCATE(iranseed)

! Now all the PEs should be at the same point in the random sequence.
! Generate a separate set of random numbers for sin and cosine coeffs
CALL RANDOM_NUMBER(rand_nums)
CALL RANDOM_NUMBER(rand_nums2)

IF (first_atmstep_call) THEN

  ! Set up AR1-process parameters and save
  oneminusalpha = 1-alpha
  ! SQRT(alpha) ensures that the noise (random nums) decorrelates
  ! faster than the timestep
  squarerootalpha = SQRT(alpha)

END IF

CALL umPrint('update_pattern: updating coeffc and coeffs', &
             src='update_pattern')

DO n = 1, stph_n2
  DO m = 0, n
    !  (rand_nums() - 0.5) provides white noise with mean of 0 and var of (1/12)
    coeffc(m,n) = oneminusalpha*coeffc(m,n) + squarerootalpha           &
                   *(rand_nums(m,n)-0.5)*pspect(n)

    coeffs(m,n) = oneminusalpha*coeffs(m,n) + squarerootalpha           &
                   *(rand_nums2(m,n)-0.5)*pspect(n)
  END DO
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE update_pattern

END MODULE update_pattern_mod
