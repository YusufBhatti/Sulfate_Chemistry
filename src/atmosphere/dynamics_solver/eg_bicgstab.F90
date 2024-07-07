! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE eg_BiCGStab_mod

USE um_parvars,            ONLY: offx, offy
USE nlsizes_namelist_mod,  ONLY: global_row_length

USE global_2d_sums_mod,    ONLY: global_2d_sums

USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim
USE atm_fields_bounds_mod
USE eg_inner_prod_mod
USE gmres1_coef_mod
USE eg_calc_ax_mod
USE eg_precon_mod
USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod,            ONLY: ereport,newline
USE errorurl_mod,           ONLY: errorurl1, errorurl2
USE gcr_input_mod,          ONLY: ItMax => GCR_max_iterations
USE um_is_nan_mod,          ONLY: um_is_nan
USE um_types,               ONLY: Real32, Real64


IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_BICGSTAB_MOD'

CONTAINS

SUBROUTINE eg_BiCGStab(x,b,tol,pre_type, l_rel_tol,             &
           row_length, rows, n_rows, model_levels,              &
           sc_err_min, init_err, fin_err, no_its, l_inc)

IMPLICIT NONE

!
! Description: Postconditioned BiCGstab algorithm from Wesseling (2002)
!              to Solve Ax = b
!
! Method:
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

REAL,    PARAMETER   :: small = 1.0e-16

! Array dimensions

INTEGER, INTENT(IN)  :: row_length
INTEGER, INTENT(IN)  :: rows
INTEGER, INTENT(IN)  :: model_levels
INTEGER, INTENT(IN)  :: n_rows
INTEGER, INTENT(IN)  :: pre_type

! Data for parallel code and domain decomposition

LOGICAL, INTENT(IN)  :: l_rel_tol, l_inc

REAL,    INTENT(IN)  :: sc_err_min
REAL,    INTENT(OUT) :: init_err, fin_err
INTEGER, INTENT(OUT) :: no_its

REAL, INTENT(IN)    ::                                                   &
    b(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

REAL, INTENT(INOUT) ::                                                   &
    x(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

! Temp work space arrays

REAL                ::                                                   &
    r(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

! Original input value of x, which is saved in case the solver is restarted

REAL                ::                                                   &
    x_in(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

! Residual error and "Stabilizer" used in the iterations

REAL                ::                                                   &
    cr(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

REAL                ::                                                   &
    Ax(1-offx:row_length+offx,1-offy:rows+offy,model_levels),            &
    p(1-offx:row_length+offx,1-offy:rows+offy,model_levels),             &
    t(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

REAL                ::                                                   &
    v(1-offx:row_length+offx,1-offy:rows+offy,model_levels),             &
    s(1-offx:row_length+offx,1-offy:rows+offy,model_levels),             &
    cs(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

! Absolute tolerance
REAL                :: tol

! Temp variables
REAL                :: rho, alf
REAL                :: omg, bet
REAL                :: nrm

! Variables used to define the normalized error
REAL                :: err, sc_err, fin_tol
! Iteration count
INTEGER             :: it, i, j, k
INTEGER             :: icode

! Logicals to record solver convergence and restart if unconverged
LOGICAL             :: l_solver_converged
LOGICAL             :: l_forced_restart

CHARACTER(LEN=errormessagelength)  :: cmessage, errorurl

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_BICGSTAB'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

errorurl = newline//"This is a common point for the model to fail if it"//  &
           newline//"has ingested or developed NaNs or infinities"//        &
           newline//"elsewhere in the code."//                              &
           newline//TRIM(errorurl1)//newline//TRIM(errorurl2)

! Back-up x in case needed for a restart

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,x,x_in)
DO k = 1, model_levels
  DO j = 1, rows
    DO i = 1, row_length
      x_in(i,j,k)  = x(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! First attempt to run the solver
! Set logicals for one pass through using standard initialisation
l_forced_restart=.FALSE.
l_solver_converged=.FALSE.

! Define error with respect to the magnitude of the forcing
sc_err = SQRT(eg_inner_prod(b,b,pdims,pdims_s))

! don't let sc_err get too small
sc_err = MAX(sc_err, sc_err_min)

DO WHILE (.NOT. l_solver_converged)

  ! Calculate initial residual

  IF ( .NOT. l_inc ) THEN

    CALL eg_Calc_Ax(Ax,x,pdims,pdims_s)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,r,Ax,cr,b)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          r(i,j,k)  = b(i,j,k) - Ax(i,j,k)
          cr(i,j,k) = r(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  ELSE
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,cr,b)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          cr(i,j,k) = b(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Note that if the solver has previously not converged after ItMax
  ! iterations, we have one more attempt by changing the direction of cr
  ! As good a way as any without changing the magnitude too much is to 
  ! set cr=ABS(cr)

  IF (l_forced_restart) THEN

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,cr)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          cr(i,j,k) = ABS(cr(i,j,k))
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF
        
! Calculate initial error

  rho = eg_inner_prod(r,r,pdims,pdims_s)
  err = SQRT(rho)

  init_err = err/sc_err
  fin_tol  = tol
  IF ( l_rel_tol ) fin_tol = init_err*tol

  alf = 1.0
  omg = 1.0
  nrm = 1.0

  DO it = 1, itmax

    IF ( it > 1 ) rho = eg_inner_prod(r,cr,pdims,pdims_s)

    bet = ( rho/nrm )*( alf/omg )

    IF ( it == 1 ) THEN
      CALL eg_precon(p,r,pre_type)
    ELSE

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,t,r,v,bet,omg)
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            t(i,j,k)   = r(i,j,k) - bet*omg*v(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

      CALL eg_precon(s,t,pre_type)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,p,s,bet)
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            p(i,j,k) = s(i,j,k) + bet*p(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
    END IF

    CALL eg_Calc_Ax(v,p,pdims,pdims_s)

    nrm = eg_inner_prod(cr,v,pdims,pdims_s)

    alf = rho/nrm

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,s,r,alf,v)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          s(i,j,k)   = r(i,j,k) - alf*v(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    CALL eg_precon(cs,s,pre_type)

    CALL eg_calc_Ax(t,cs,pdims,pdims_s)

    omg = gmres1_coef(s,t,err,pdims,pdims_s)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,x,alf,    &
!$OMP             p,omg,cs,r,s,t)
    DO k = 1, model_levels
      x(:,:,k)   = x(:,:,k) + alf*p(:,:,k) + omg*cs(:,:,k)
      DO j = 1, rows
        DO i = 1, row_length
          r(i,j,k)   = s(i,j,k) - omg*t(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    nrm = rho

    ! Calculate residual 

    err = SQRT(err)/sc_err

    ! Check for "bad" values

    IF ( um_is_nan(err) ) THEN
      icode = 1
      WRITE(cmessage,'(A,I6,A)') 'NaNs in error term in BiCGstab after ', &
           it, ' iterations '//TRIM(errorurl)
      CALL ereport('EG_BICGSTAB',icode,cmessage)
    END IF

    IF ( ABS(omg) < small ) THEN
      icode = 11
      WRITE(cmessage,'(A,I6,A)') 'Convergence failure in BiCGstab after ', &
           it, ' iterations: omg is too small'//                           &
           newline//TRIM(errorurl1)//newline//TRIM(errorurl2)
      CALL ereport('EG_BICGSTAB',icode,cmessage)
    END IF

    ! Check residual for convergence

    IF ( L_inc ) THEN
      IF ( err < fin_tol .AND. it > 1 ) THEN
        ! Converged, no need to restart solver
        l_solver_converged=.TRUE.
        EXIT ! EXIT loop over it = 1, itmax
      END IF
    ELSE
      IF ( err < fin_tol ) THEN
        ! Converged, no need to restart solver
        l_solver_converged=.TRUE.
        EXIT ! EXIT loop over it = 1, itmax
      END IF
    END IF

  END DO

  ! Also trap convergence failure

  IF (.NOT. l_solver_converged) THEN
    ! If this is the first convergence failure, force a restart of the solver
    IF (.NOT. l_forced_restart) THEN
      icode = -2
      ! Print out a warning message
      WRITE(cmessage,'(A,I6,A)') 'Convergence failure in BiCGstab after ', &
           itmax, ' iterations: attempting to restart solver'
      CALL ereport('EG_BICGSTAB_MIXED_PREC',icode,cmessage)

      ! Reset x() to its initial value

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,x,x_in)
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            x(i,j,k)  = x_in(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

      ! Finally, set flag to restart solver
      l_forced_restart=.TRUE.
    ELSE
      ! If this is not the first convergence failure then fail
      icode = 2
      ! Include shortened error message
      ! (as it is unlikely that this point of failure is due to NaNs)
      cmessage='Convergence failure in BiCGstab after restart'// &
           TRIM(errorurl1)//newline//TRIM(errorurl2)
      CALL ereport('EG_BICGSTAB_MIXED_PREC',icode,cmessage)
    END IF
  END IF ! .NOT. l_solver_converged
END DO 

fin_err = err
no_its  = it

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_BiCGStab

!!!!!!!!!!!!!!!!!!!!!!! End of BiCGStab in double !!!

SUBROUTINE eg_BiCGStab_mixed_prec(x,b,tol,pre_type, l_rel_tol,          &
           row_length, rows, n_rows, model_levels,              &
           sc_err_min, init_err, fin_err, no_its, l_inc)

IMPLICIT NONE

!
! Description: Postconditioned BiCGstab algorithm from Wesseling (2002)
!              to Solve Ax = b
!              In single precision. Args passed in are double,
!              temporaries are single
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

REAL(KIND=real64),    PARAMETER   :: small = 1.0e-16

! Array dimensions

INTEGER, INTENT(IN)    :: row_length
INTEGER, INTENT(IN)    :: rows
INTEGER, INTENT(IN)    :: model_levels
INTEGER, INTENT(IN)    :: n_rows
INTEGER, INTENT(IN)    :: pre_type

! Data for parallel code and domain decomposition

LOGICAL, INTENT(IN)  :: l_rel_tol, l_inc

REAL(KIND=real64),    INTENT(IN)  :: sc_err_min
REAL(KIND=real64),    INTENT(OUT) :: init_err, fin_err
INTEGER, INTENT(OUT) :: no_its

REAL(KIND=real64), INTENT(IN)    ::                                      &
    b(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
REAL(KIND=real64), INTENT(INOUT) ::                                      &
    x(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

! Temp work space arrays in single precision

REAL(KIND=real64)                ::                                      &
    r(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

! Original input value of x, which is saved in case the solver is restarted

REAL(KIND=real64)                ::                                      &
    x_in(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

! Residual error and "Stabilizer" used in the iterations

REAL(KIND=real32)                ::                                      &
    cr(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
REAL(KIND=real32)                ::                                      &
    cs(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

REAL(KIND=real32) ::                                                     &
    Ax(1-offx:row_length+offx,1-offy:rows+offy,model_levels),            &
    p(1-offx:row_length+offx,1-offy:rows+offy,model_levels)


REAL(KIND=real32)                ::                                      &
    v(1-offx:row_length+offx,1-offy:rows+offy,model_levels),             &
    s(1-offx:row_length+offx,1-offy:rows+offy,model_levels),             &
    t(1-offx:row_length+offx,1-offy:rows+offy,model_levels)


! Absolute tolerance
REAL(KIND=real64)   :: tol

! Temp variables
REAL (KIND=real32)  :: rho, alf, bet
REAL (KIND=real32)  :: omg
REAL (KIND=real32)  :: nrm

! Variables used to define the normalized error
REAL(KIND=real64)   ::  fin_tol,  sc_err
REAL(KIND=real32)   :: err

! Iteration count
INTEGER             :: it, i, j, k
INTEGER             :: icode

! Logicals to record solver convergence and restart if unconverged
LOGICAL             :: l_solver_converged
LOGICAL             :: l_forced_restart

CHARACTER(LEN=errormessagelength)   :: cmessage, errorurl

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_BICGSTAB_MIXED_PREC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

errorurl = newline//"This is a common point for the model to fail if it"//  &
           newline//"has ingested or developed NaNs or infinities"//        &
           newline//"elsewhere in the code."//                              &
           newline//TRIM(errorurl1)//newline//TRIM(errorurl2)

! Back-up x in case needed for a restart

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,x,x_in)
DO k = 1, model_levels
  DO j = 1, rows
    DO i = 1, row_length
      x_in(i,j,k)  = x(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! First attempt to run the solver
! Set logicals for one pass through using standard initialisation
l_forced_restart=.FALSE.
l_solver_converged=.FALSE.

! Define error with respect to the magnitude of the forcing
sc_err = SQRT(eg_inner_prod(b,b,pdims,pdims_s))

! don't let sc_err get too small
sc_err = MAX(sc_err, sc_err_min)


DO WHILE (.NOT.l_solver_converged )

  ! Calculate initial residual
  
  IF ( .NOT. l_inc ) THEN
    
    CALL eg_Calc_Ax(Ax,x,pdims,pdims_s)
    
    
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,r,Ax,cr,b)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          r(i,j,k)  = b(i,j,k) - Ax(i,j,k)
          cr(i,j,k) = r(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  ELSE

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,cr,b)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          cr(i,j,k) = b(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

! Note that if the solver has previously not converged after ItMax
! iterations, we have one more attempt by changing the direction of cr
! As good a way as any without changing the magnitude too much is to
! set cr=ABS(cr)

  IF (l_forced_restart) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,cr)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          cr(i,j,k) = ABS(cr(i,j,k))
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Calculate initial error
  rho = eg_inner_prod(r,r,pdims,pdims_s)
  
  err = SQRT(rho)

  init_err = err/sc_err
  fin_tol  = tol
  IF ( l_rel_tol ) fin_tol = init_err*tol

  alf = 1.0
  omg = 1.0
  nrm = 1.0
  
  DO it = 1, itmax

    IF ( it > 1 ) rho = eg_inner_prod(r,cr,pdims,pdims_s)
    
    bet = ( rho/nrm )*( alf/omg )
    
    IF ( it == 1 ) THEN
      
      CALL eg_precon(p,r,pre_type)
    ELSE

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,t,r,v,bet,omg)
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            t(i,j,k)   = r(i,j,k) - bet*omg*v(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

      CALL eg_precon(s,t,pre_type)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,p,s,bet)
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            p(i,j,k) = s(i,j,k) + bet*p(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
    END IF


    CALL eg_Calc_Ax(v,p,pdims,pdims_s)
    nrm = eg_inner_prod(cr,v,pdims,pdims_s)
    alf = rho/nrm

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,s,r,alf,v)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          s(i,j,k)   = r(i,j,k) - alf*v(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
    CALL eg_precon(cs,s,pre_type)

    CALL eg_calc_Ax(t,cs,pdims,pdims_s)

    omg = gmres1_coef(s,t,err,pdims,pdims_s)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,x,alf,    &
!$OMP             p,omg,cs,r,s,t)
    DO k = 1, model_levels
      x(:,:,k)   = x(:,:,k) + alf*p(:,:,k) + omg*cs(:,:,k)
      DO j = 1, rows
        DO i = 1, row_length
          r(i,j,k)   = s(i,j,k) - omg*t(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
    nrm = rho

    ! Calculate residual

    err = SQRT(err)/sc_err

    ! Check for "bad" values

    IF ( um_is_nan(err) ) THEN
      icode = 1
      WRITE(cmessage,'(A,I6,A)') 'NaNs in error term in BiCGstab after ', &
           it, ' iterations '//TRIM(errorurl)
      CALL ereport('EG_BICGSTAB',icode,cmessage)
    END IF

    IF ( ABS(omg) < small ) THEN
      icode = 11
      WRITE(cmessage,'(A,I6,A)') 'Convergence failure in BiCGstab after ', &
           it, ' iterations: omg is too small'//                           &
           newline//TRIM(errorurl1)//newline//TRIM(errorurl2)
      CALL ereport('EG_BICGSTAB',icode,cmessage)
    END IF

    ! Check residual for convergence

    IF ( err < fin_tol ) THEN
      ! Converged, no need to restart solver
      l_solver_converged=.TRUE.
      EXIT ! EXIT loop over it = 1, itmax
    END IF

  END DO

  ! Also trap convergence failure
  IF (.NOT. l_solver_converged) THEN
    ! If this is the first convergence failure, force a restart of the solver
    IF (.NOT. l_forced_restart) THEN
      icode = -2
      ! Print out a warning message
      WRITE(cmessage,'(A,I6,A)') 'Convergence failure in BiCGstab after ', &
           itmax, ' iterations: attempting to restart solver'
      CALL ereport('EG_BICGSTAB_MIXED_PREC',icode,cmessage)
      ! Reset x() to its initial value
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k) &
!$OMP             SHARED(model_levels,rows,row_length,x,x_in)
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            x(i,j,k)  = x_in(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
      ! Finally, set flag to restart solver
      l_forced_restart=.TRUE.
    ELSE
      ! If this is not the first convergence failure then fail
      icode = 2
      ! Include shortened error message
      ! (as it is unlikely that this point of failure is due to NaNs)
      cmessage='Convergence failure in BiCGstab after restart'// &
           TRIM(errorurl1)//newline//TRIM(errorurl2)
      CALL ereport('EG_BICGSTAB_MIXED_PREC',icode,cmessage)
    END IF
  END IF ! .NOT. l_solver_converged
END DO

fin_err = err
no_its  = it

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_BiCGStab_mixed_prec


END MODULE eg_BiCGStab_mod
