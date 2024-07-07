! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_check_sl_domain_mod
USE umPrintMgr
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_CHECK_SL_DOMAIN_MOD'

CONTAINS
SUBROUTINE eg_check_sl_domain(                                        &
                        depart_phi, depart_lambda,                    &
                        row_length, rows_depart, model_levels,        &
                        max_lambda, min_lambda,                       &
                        max_phi, min_phi,                             &
                        length, width)

USE conversions_mod, ONLY: pi
USE parkind1,    ONLY: jpim, jprb       !DrHook
USE yomhook,     ONLY: lhook, dr_hook   !DrHook
USE ereport_mod, ONLY: ereport

USE model_domain_mod, ONLY: model_type, mt_global, mt_lam, mt_cyclic_lam       &
                          , mt_bi_cyclic_lam

IMPLICIT NONE
!
! Description:
!          Checks data lies inside model domain.
!
!
! Method:
!
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_CHECK_SL_DOMAIN'

INTEGER  ::                                                           &
  row_length                                                          &
               ! Dimension of u_adv array in i direction.
, rows_depart                                                         &
               ! Dimension of depart arrays in j direction.
, model_levels
               ! Dimension of u_adv array in k direction.


REAL ::                                                               &
  max_lambda                                                          &
, max_phi                                                             &
, min_lambda                                                          &
, min_phi

! Arguments with Intent IN/OUT.

REAL ::                                                               &
                ! Departure point co-ordinates.
  depart_lambda (row_length, rows_depart, model_levels)               &
, depart_phi (row_length, rows_depart, model_levels)

! Make these optional for time being !
REAL, OPTIONAL :: Length, Width

! Local Variables.

! scalars

INTEGER :: i, j, k      ! Loop indices

INTEGER ::                                                            &
  cnt

INTEGER :: ierr
! External Routines: None

! Functions: None

! Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check optional arguments

IF (model_type == mt_cyclic_lam .OR.                                &
    model_type == mt_bi_cyclic_lam ) THEN

  ! check horizontal in the periodic in x LAM at u points.

  IF (.NOT. PRESENT(width) .OR. .NOT. PRESENT(length)) THEN
    ierr = 1
    CALL ereport("eg_check_sl_domain", ierr,                          &
        " missing length or width in call" )
  END IF
END IF

! ----------------------------------------------------------------------
! Section 1.   Check points lie inside model domain.
! ----------------------------------------------------------------------

IF (model_type == mt_global) THEN
  ! check horizontal
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,rows_depart,row_length,   &
!$OMP& depart_lambda,min_lambda,max_lambda)
  DO k = 1, model_levels
    DO j = 1, rows_depart
      DO i = 1, row_length

        IF (depart_lambda(i,j,k)  <   min_lambda ) THEN
          depart_lambda(i,j,k) = pi*2.0 + depart_lambda(i,j,k)
        ELSE IF (depart_lambda(i,j,k)  >=  max_lambda ) THEN
          depart_lambda(i,j,k) = depart_lambda(i,j,k) - pi*2.0
        END IF

      END DO
    END DO
  END DO
!$OMP END PARALLEL DO 

ELSE IF (model_type == mt_lam) THEN

  ! check horizontal in the LAM at u points.

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,rows_depart,row_length,depart_phi,   &
!$OMP& min_phi,max_phi,depart_lambda,min_lambda,max_lambda)
  DO k = 1, model_levels
    DO j = 1, rows_depart
      DO i = 1, row_length
        IF (depart_phi(i,j,k)  <   min_phi ) THEN
          depart_phi(i,j,k) = min_phi
        ELSE IF (depart_phi(i,j,k) >  max_phi) THEN
          depart_phi(i,j,k) = max_phi
        END IF
        IF (depart_lambda(i,j,k)  <   min_lambda ) THEN
          depart_lambda(i,j,k) = min_lambda
        ELSE IF (depart_lambda(i,j,k) >  max_lambda) THEN
          depart_lambda(i,j,k) = max_lambda
        END IF
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO 

ELSE IF (model_type == mt_cyclic_lam) THEN

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP& SHARED(model_levels,rows_depart,row_length,depart_phi,   &
!$OMP& min_phi,max_phi,depart_lambda,min_lambda,max_lambda,length)
  DO k = 1, model_levels
    DO j = 1, rows_depart
      DO i = 1, row_length
        IF (depart_phi(i,j,k)  <   min_phi ) THEN
          depart_phi(i,j,k) = min_phi
        ELSE IF (depart_phi(i,j,k) >  max_phi) THEN
          depart_phi(i,j,k) = max_phi
        END IF
        IF (depart_lambda(i,j,k)  <   min_lambda ) THEN
          depart_lambda(i,j,k) = depart_lambda(i,j,k) + length
        ELSE IF (depart_lambda(i,j,k) >= max_lambda) THEN
          depart_lambda(i,j,k) = depart_lambda(i,j,k) - length
        END IF
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO 

ELSE IF (model_type == mt_bi_cyclic_lam) THEN
  ! check horizontal in the periodic in both x and y LAM at u points
  cnt = 0

  DO k = 1, model_levels
    DO j = 1, rows_depart
      DO i = 1, row_length
        IF (depart_phi(i,j,k)  <   min_phi) THEN
          depart_phi(i,j,k) = depart_phi(i,j,k) + width
          cnt = cnt + 1
        ELSE IF (depart_phi(i,j,k)  >   max_phi) THEN
          depart_phi(i,j,k) = depart_phi(i,j,k) - width
          cnt = cnt + 1
        END IF
        IF (depart_lambda(i,j,k)  <   min_lambda) THEN
          depart_lambda(i,j,k) = depart_lambda(i,j,k) + length
        ELSE IF (depart_lambda(i,j,k)  >   max_lambda) THEN
          depart_lambda(i,j,k) = depart_lambda(i,j,k) - length
        END IF
      END DO
    END DO
  END DO
  IF ( cnt > 0 ) THEN
    WRITE(umMessage,*)'WARNING moving NS wrap point on ',cnt,' grid points'
    CALL umPrint(umMessage,src='eg_check_sl_domain')
  END IF
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_check_sl_domain
END MODULE eg_check_sl_domain_mod
