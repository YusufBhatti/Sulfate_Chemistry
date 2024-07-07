! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE scmoutput_mod

IMPLICIT NONE

CONTAINS


! Stub version of Single-Column Model diagnostic output routine

SUBROUTINE scmoutput                                                          &
  ( x, sname, lname, units, timprof, domprof, streams, sname2                 &
  , calling_routine )


IMPLICIT NONE


! Description:
!   SCM output stub routine for non-scm model configurations
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model
!
! Code Description:
!   Language: Fortran77 with some bits of Fortran90

INTEGER, INTENT(IN) :: domprof         ! Domain profile for the diagnostic
REAL,    INTENT(IN) :: x(*)            ! Input variable
! NOTE: with scmoutput inside a module, the input array x needs to be
! declared with assumed size (*).  Declaring it with assumed shape (:), (:,:)
! or (:,:,:) causes compile errors in the full UM wherever different-ranked
! arrays are passed into scmoutput.
! In the actual SCM version of this routine, we don't run into this problem,
! because in that case scmoutput knows what shape all the arrays passed to it
! should have (so x is declared as an automatic array, for which fortran
! allows rank mis-matches in argument lists).

CHARACTER(LEN=*), INTENT(IN) :: sname  ! Diag Short name, (unique)
CHARACTER(LEN=*), INTENT(IN) :: lname  ! Diag Long name
CHARACTER(LEN=*), INTENT(IN) :: units  ! Diag Units
CHARACTER(LEN=*), INTENT(IN) :: sname2 ! Short name of another, previously
                                       ! defined diagnostic which will be
                                       ! used in the construction of this
                                       ! one according to the time profile

CHARACTER(LEN=*), INTENT(IN) :: calling_routine
                                       ! Routine that has called scmoutput

INTEGER, INTENT(IN) :: timprof         ! Diag time profile
INTEGER, INTENT(IN) :: streams         ! An encoded integer specifying
                                       ! which output streams the diagnostic
                                       ! is to go to

RETURN

END SUBROUTINE scmoutput


END MODULE scmoutput_mod
