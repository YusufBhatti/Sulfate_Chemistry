! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINES ADDINC  RELAXC ---------------------------------------
!
!    Purpose : Add analysis increments to model field.
!    Purpose : Scale analysis increments by relaxation coefficients.
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!
!
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: AC Assimilation
MODULE relaxc_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RELAXC_MOD'

CONTAINS

SUBROUTINE relaxc (incr,LEN,npts,nrows,relax,icode,cmessage)

!     SCALE AC INCREMENTS (INCR) BY RELAXATION COEFFICIENT (RELAX)
!     WHICH IS A FUNCTION OF ROW
!     DIMENSIONS GIVEN BY LEN,NPTS,NROWS


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Atmos_Max_Sizes
USE UM_ParParams
USE comobs_mod, ONLY: nobtypmx
USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER :: nrows,npts,LEN
REAL ::                                                           &
               incr(LEN),                                         &
               relax(nrows)
!
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER :: IS,jrow,jpt

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RELAXC'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (ltimer_ac) CALL timer('RELAXC  ',3)

!     MULTIPLY BY RELAXATION COEFFICIENT
IS=1
DO jrow=1,nrows
  DO jpt=0,npts-1
    incr(IS+jpt) = incr(IS+jpt)*relax(jrow)
  END DO
  IS=is+npts
END DO

IF (ltimer_ac) CALL timer ('RELAXC  ',4)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE relaxc
END MODULE relaxc_mod
