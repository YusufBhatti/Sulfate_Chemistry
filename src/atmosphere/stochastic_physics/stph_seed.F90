! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stochastic Physics
MODULE stph_seed_mod

USE umPrintMgr
USE dump_headers_mod, ONLY: ih_sp_seed, a_inthd, ih_stph_n1,            &
                            ih_stph_n2, ih_stochastic_flag
USE UM_ParVars

IMPLICIT NONE

CONTAINS

SUBROUTINE stph_seed_copy_from_dump

USE UM_ParCore, ONLY: mype
USE stochastic_physics_run_mod, ONLY: stph_seed_present

IMPLICIT NONE
! local temporary arrays
CHARACTER(LEN=256) :: cmessage    ! OUT error message

! random seeds
INTEGER, ALLOCATABLE :: iranseed(:)
INTEGER :: i, tam

IF (mype == 0) THEN
  CALL umPrint('Reading iranseed from dump headers',src='stph_seed')

#if defined(IBM_XL_FORTRAN)
  ! IBM documentation suggests using generator=2 to improve random
  ! cycling period and calculation speed.
  ! http://publib.boulder.ibm.com/infocenter/comphelp/v8v101/index.jsp \
  ! ?topic=/com.ibm.xlf101a.doc/xlfcr/runpgms.htm
  CALL RANDOM_SEED(generator=2)
#endif
  !  Get size of random seed for read/write/allocate
  CALL RANDOM_SEED(SIZE=tam)
  IF (.NOT. ALLOCATED(iranseed)) ALLOCATE(iranseed(tam))

  iranseed(:) = 0.0
  
  IF (a_inthd(ih_stochastic_flag) >= stph_seed_present) THEN
    iranseed(1) = a_inthd(ih_sp_seed)
    
    WRITE(umMessage, '(A,I30)') 'Value of iranseed read from dump ', iranseed(1)
    CALL umPrint(umMessage,src='stph_seed')
  END IF
  
  CALL RANDOM_SEED(put=iranseed(1:tam))
    
  DEALLOCATE(iranseed)
END IF
  
END SUBROUTINE stph_seed_copy_from_dump

SUBROUTINE stph_seed_copy_to_dump

USE UM_ParCore, ONLY: mype
USE stochastic_physics_run_mod, ONLY: stph_n1, stph_n2, stph_header_flag

IMPLICIT NONE

! local temporary arrays
CHARACTER(LEN=256) :: cmessage    ! OUT error message

! random seeds
INTEGER, ALLOCATABLE :: iranseed(:)
INTEGER :: i, tam

IF (mype == 0) THEN
#if defined(IBM_XL_FORTRAN)
  ! IBM documentation suggests using generator=2 to improve random
  ! cycling period and calculation speed.
  ! http://publib.boulder.ibm.com/infocenter/comphelp/v8v101/index.jsp \
  ! ?topic=/com.ibm.xlf101a.doc/xlfcr/runpgms.htm
  CALL umPrint('Copying iranseed to dump headers',src='stph_seed')
  
  CALL RANDOM_SEED(generator=2)
#endif
  !  Get size of random seed for read/write/allocate
  CALL RANDOM_SEED(SIZE=tam)
  IF (.NOT. ALLOCATED(iranseed)) ALLOCATE(iranseed(tam))

  iranseed(:) = 0.0
  
  CALL RANDOM_SEED(get=iranseed(1:tam))

  ! put seed value into integer header array for outputting to dump
  a_inthd(ih_stochastic_flag) = stph_header_flag
  a_inthd(ih_stph_n1) = stph_n1
  a_inthd(ih_stph_n2) = stph_n2
  a_inthd(ih_sp_seed) = iranseed(1)
  
  WRITE(umMessage, '(A,I30)')                                           &
    'Value of iranseed copied to the dump headers', iranseed(1)
  CALL umPrint(umMessage,src='stph_seed')
  

  DEALLOCATE(iranseed)
END IF

END SUBROUTINE stph_seed_copy_to_dump

! Handles STOCHASTIC PHYSICS Random Seed:
! Synchronise current seed from a specific PE to all others
SUBROUTINE stph_seed_sync(pe)

USE UM_ParCore, ONLY: mype, nproc

IMPLICIT NONE 

! Which PE should be used as the master 
INTEGER, INTENT(IN) :: pe

! Size of the random seed array
INTEGER :: seedcount   
INTEGER :: icode

! Random seeds 
INTEGER, ALLOCATABLE :: iranseed(:)

! --------------------------------------------------------------------

IF (nproc > 1) THEN
  ! Get size of random seed  and allocate array
  CALL RANDOM_SEED(SIZE=seedcount)

  IF (.NOT. ALLOCATED(iranseed)) ALLOCATE(iranseed(seedcount))

  ! Get the random seed values from the master processor and broadcast them to
  ! all the PEs.
  IF ( mype == pe ) THEN
    CALL RANDOM_SEED(GET=iranseed(1:seedcount))
  END IF  

  ! Do the broadcast
  CALL gc_ibcast(3245, seedcount, pe, nproc, icode, iranseed)

  ! Set the random seed for all other PEs 
  IF (mype /= pe) THEN
    CALL RANDOM_SEED(PUT=iranseed(1:seedcount))
  END IF

  IF (ALLOCATED(iranseed)) DEALLOCATE(iranseed)
END IF  

END SUBROUTINE stph_seed_sync

SUBROUTINE stph_seed_gen_from_date(dt_in)

IMPLICIT NONE 

INTEGER, INTENT(IN), OPTIONAL :: dt_in(8) ! date/time info; if not provided
                                          ! current time will be used
INTEGER :: dt(8)
INTEGER :: tam      ! size of the random seed
INTEGER :: iarg     ! random seed input argument
INTEGER :: max_iarg ! maximum integer size of seed argument

! random seeds
INTEGER, ALLOCATABLE :: iranseed(:)
INTEGER, ALLOCATABLE :: prevseed(:)

! random numbers for calculating seed
REAL, ALLOCATABLE :: rnum(:)

#if defined(IBM_XL_FORTRAN)
  ! IBM documentation suggests using generator=2 to improve random
  ! cycling period and calculation speed.
  ! http://publib.boulder.ibm.com/infocenter/comphelp/v8v101/index.jsp \
  ! ?topic=/com.ibm.xlf101a.doc/xlfcr/runpgms.htm
  CALL RANDOM_SEED(generator=2)
#endif
CALL RANDOM_SEED(SIZE=tam)

IF (PRESENT(dt_in)) THEN
  dt = dt_in
ELSE
  CALL DATE_AND_TIME(values=dt)
END IF

IF (.NOT. ALLOCATED(prevseed)) ALLOCATE(prevseed(tam))
IF (.NOT. ALLOCATED(iranseed)) ALLOCATE(iranseed(tam))
IF (.NOT. ALLOCATED(rnum)) ALLOCATE(rnum(tam))

! Use full date in randomising seed to reduce chance of
!  recycling from one run to the next.
! The formula calculates the days since ~2000AD and adds
!  in time suitably inflated to fully change the seed.
! Only use last two digits of year to prevent numerical
!  overflow at some date in the future.
! A random number generated from this seed is used to
!  multiply the seed again
dt(1) = dt(1) - 100*INT(0.01*dt(1))
iarg = (dt(3) - 32075 +                                             &
         1461*(dt(1) + 4800 + (dt(2) - 14)/12)/4 +                  &
         367*(dt(2) - 2 - (dt(2)-14)/12*12)/12 -                    &
         3*((dt(1)+4900+(dt(2)-14)/12)/100)/4)*1000 +               &
         dt(8)**2.86 + dt(7)**3.79 + dt(5)**5.12 +                  &
         dt(6)**3.24
! Constrain iarg in a range to prevent numerical overflow
! 2**8 < iarg < SQRT(HUGE(INT)), since iarg is raised to n<2
max_iarg = FLOOR(SQRT(1.0*HUGE(iarg)))
iarg = MOD(iarg, max_iarg)
iarg = MAX(iarg, 256)
! Generate initial random number set to use as power
prevseed(:) = iarg
CALL RANDOM_SEED(put=prevseed(1:tam))
CALL RANDOM_NUMBER(rnum)

#if defined(IBM_XL_FORTRAN)
! Range of seed from ~ 37640 => 9.22E18 when IBM GENERATOR=2
iranseed(:)=iarg**(1.9 + 0.1*rnum(:))
#else
! Range of seed from 0 to 2**31 (32-bit Int)
iranseed(:)=iarg*rnum(:)
#endif
! Set final seed
CALL RANDOM_SEED(put=iranseed(1:tam))

IF (printstatus  >=  prstatus_normal) THEN
    CALL umPrint('STPH SETUP: Size of Random Seed',src='stph_setup')
    WRITE(umMessage,'(I10)') tam
  IF (PRESENT(dt_in)) THEN
    CALL umPrint( 'STPH SETUP: Random seed generated using model '    &
                // 'date')
  ELSE
    CALL umPrint('STPH SETUP: Random seed generated using computer '  &
                // 'date', src='stph_setup')
  END IF
  WRITE(umMessage,'(100I30)') iranseed
  CALL umPrint(umMessage,src='stph_setup')
END IF

DEALLOCATE(prevseed)
DEALLOCATE(iranseed)
DEALLOCATE(rnum)

END SUBROUTINE stph_seed_gen_from_date

END MODULE stph_seed_mod
