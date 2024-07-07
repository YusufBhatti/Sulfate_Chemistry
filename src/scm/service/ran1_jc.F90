! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model
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

!    Programming standards :
!      Fortran 95
!------------------------------------------------------------------
MODULE ran1_jc_mod

USE random_num_var

IMPLICIT NONE


CONTAINS

REAL FUNCTION ran1_jc()

IMPLICIT NONE

INTEGER, PARAMETER :: &
  a    = 16807        &
, m    = 2147483647   &
, q    = 127773       &
, r    = 2836         &
, ndiv = 1 + (m-1)/tabs

REAL, PARAMETER ::   &
  am     = 1.0/m     &
, eps    = 1.2e-7    &
, lfvlt1 = 1.0-eps

INTEGER :: k

IF (drive <= 0 .OR. iy == 0) THEN 
  drive = MAX(-drive,1)    

  DO jr=tabs+8, 1, -1      
    k = drive/q
    drive = a*(drive-k*q)-r*k

    IF (drive <  0) THEN
      drive = drive + m
    END IF

    IF (jr <= tabs) THEN  
      iv(jr) = drive
    END IF

  END DO
  iy = iv(1)
END IF

k = drive/q                 
drive = a*(drive-k*q)-r*k  
                              
IF (drive <  0) THEN
  drive = drive+m
END IF

jr  = 1 + iy/ndiv            
iy  = iv(jr)                  
                              
iv(jr) = drive
ran1_jc = MIN(am*iy,lfvlt1)                        

RETURN

END FUNCTION ran1_jc

END MODULE ran1_jc_mod
