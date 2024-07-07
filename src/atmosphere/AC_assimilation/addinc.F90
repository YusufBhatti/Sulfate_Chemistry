! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINES ADDINC  RELAXC ---------------------------------------
!
!    Purpose : Add analysis increments to model field.
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!
!
!    ARGUMENTS
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE addinc_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ADDINC_MOD'

CONTAINS

SUBROUTINE addinc (field,incr,LEN,npts,nrows,navt,icode,cmessage)

!  ADD AC INCREMENTS (INCR) TO MODEL PROGNOSTIC VARIABLES (FIELD)
!  WHICH IS A FUNCTION OF ROW
!  DIMENSIONS GIVEN BY LEN,NPTS,NROWS

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Atmos_Max_Sizes
USE UM_ParParams
USE comobs_mod, ONLY: nobtypmx

USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength
USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE

INTEGER :: npts,nrows,LEN
REAL ::                                                           &
               field(LEN),                                        &
               incr(LEN)
!
INTEGER :: icode
INTEGER :: navt
CHARACTER(LEN=errormessagelength) :: cmessage

INTEGER :: IS,irow,ipt

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ADDINC'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (ltimer_ac) CALL timer('ADDINC  ',3)

!   ADD INCREMENTS TO MODEL FIELD
IF (model_type == mt_global) THEN
  !   FOR GLOBAL DOMAIN

  DO ipt=1,LEN
    field(ipt) = field(ipt)+incr(ipt)
  END DO

ELSE
  !   USE ONLY THOSE VALUES NOT NEAR EDGE OF LIMITED AREA

  DO irow=1,nrows
    IS=(irow-1)*npts
    DO ipt=1,npts
      field(IS+ipt) = field(IS+ipt)+incr(IS+ipt)
    END DO
  END DO
  !
END IF ! IF GLOBAL

! ensure RH increments don't generate -ve moisture points
IF (navt == 4) THEN
  DO ipt=1,LEN
    IF (field(ipt) <= 0.0)field(ipt)=0.0
  END DO
END IF

IF (ltimer_ac) CALL timer('ADDINC  ',4)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE addinc
END MODULE addinc_mod
