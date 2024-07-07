! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Purpose random number generation routine
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model
!
!    Programming standards :
!      Fortran 95 , UMDP3
!----------------------------------------------------------------
!
!     Pseudo-random real numbers, Gaussian distribution of
!     mean m, standard deviation s.
!
!     Knuth D E 1981 The Art of Computer Programming (Volume 2) Addison-Wesley
!
!     Park, SK and Miller K.W. 1988 Communications of the ACM, 
!     vol 31 pp 1192-1201
!
!     G. E. P. Box and Mervin E. Muller, 
!     A Note on the Generation of Random Normal Deviates, 
!     The Annals of Mathematical Statistics (1958), Vol. 29, No. 2 pp. 610–611

MODULE random_num_gen

USE random_num_var

IMPLICIT NONE

CONTAINS

SUBROUTINE random_initial(iseed)

!     Initialize the random number generator ran1_jc to give
!     repeatable sequence
!
!     This routine *needs* to be called first when using the
!     random num generator routines 

IMPLICIT NONE

! Argument
INTEGER :: iseed

! ----------------------
! Zero the arrays of the ran1_jc function
DO jr=1, tabs
  iv(jr) = 0
END DO
iy = 0

! Initialize the ran1_jc function with the user provided value:
drive = - INT(iseed)      ! drive must be negative when
                          ! initializing ran1_jc.
RETURN
END SUBROUTINE random_initial


! -----------------------------------------------------------------

SUBROUTINE random_state(drive_out, iv_out, iy_out)

!     Save the state of the random number generator

IMPLICIT NONE
! Arguments to extract from the memory of the random routines
INTEGER, INTENT(OUT) :: drive_out,iv_out(tabs),iy_out

! Copy the state values straight into the _out variables.
DO jr=1, tabs
  iv_out(jr) = iv(jr)
END DO

drive_out = drive
iy_out = iy

RETURN

END SUBROUTINE random_state

! ------------------------------------------------------------------

SUBROUTINE restore_random_state(drive_in,iv_in,iy_in)

!     Restore the state of the random number generator

IMPLICIT NONE

INTEGER, PARAMETER :: tabs = 32

! Arguments to put into the memory of the random routines
INTEGER, INTENT(IN) :: drive_in,iv_in(tabs),iy_in

! Set the values within the number generator
DO jr=1, tabs
  iv(jr) = iv_in(jr)
END DO

drive = drive_in
iy    = iy_in

RETURN

END SUBROUTINE restore_random_state

! ------------------------------------------------------------------

REAL FUNCTION random_func(m,s)

USE ran1_jc_mod, ONLY: ran1_jc

IMPLICIT NONE

! Arguments
REAL :: m,s         ! (mean, standard deviation) of the distribution.

! Local variables
LOGICAL, SAVE :: flag
REAL,    SAVE :: keeper
REAL          :: factor,square,value1,value2

DATA flag/.FALSE./

IF (flag) THEN 
  random_func = m+s*keeper       
  flag        = .FALSE.      
ELSE 
  DO                  
    value1 = 2.0*ran1_jc()-1.0                                    
    value2 = 2.0*ran1_jc()-1.0                 
    square = value1**2+value2**2                      
    IF (square < 1.0 .AND. square > 0.0) EXIT
  END DO
                                         
  factor      = SQRT(-2.0*LOG(square)/square) 
  keeper      = value1*factor
  random_func = m+s*(value2*factor)
  flag        = .TRUE.                
END IF

RETURN

END FUNCTION random_func

!-------------------------------------

END MODULE random_num_gen

