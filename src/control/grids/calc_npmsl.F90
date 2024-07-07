! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine calc_npmsl
MODULE calc_npmsl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALC_NPMSL_MOD'

CONTAINS

SUBROUTINE calc_npmsl                                                          &
  ( pmsl_local, pstar_local, phi_star_local, thetas_local, tm_local            &
  , cos_p_latitude_local, delta_lambda, delta_phi                              &
  , row_length_local, rows_local, row_length, rows                             &
  , me, n_proc, n_procx, n_procy, off_x, off_y, neighbour, at_extremity        &
  , all_proc_group, npmsl_height )

!-------------------------------------------------------------

! Purpose:
!          Modifies Pressure at mean sea level as calculated
!          in Calc_PMSL to account for errors due to high
!          orography

! Method:
!          As described in section 4.5 of UMDP 80. R. Rawlins

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Grids

! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards. v10.3


USE planet_constants_mod, ONLY: g, planet_radius,                              &
                                r, cp, kappa, p_zero, recip_kappa
USE conversions_mod, ONLY: pi
USE yomhook,         ONLY: lhook, dr_hook
USE parkind1,        ONLY: jprb, jpim
USE UM_ParParams
USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_lam

IMPLICIT NONE

INTEGER, INTENT(IN) :: off_x
INTEGER, INTENT(IN) :: off_y
INTEGER, INTENT(IN) :: me              ! Processor number
INTEGER, INTENT(IN) :: n_proc
INTEGER, INTENT(IN) :: n_procx
INTEGER, INTENT(IN) :: n_procy
INTEGER, INTENT(IN) :: neighbour(4)    ! Array with the Ids of the four
                                       ! neighbours in the horizontal plane
INTEGER, INTENT(IN) :: all_proc_group  ! Group id for all processors
LOGICAL, INTENT(IN) :: at_extremity(4) ! Indicates if this processor is at
                                       ! north, south, east or west of the
                                       ! processor grid
REAL, INTENT(IN) :: delta_lambda
REAL, INTENT(IN) :: delta_phi
REAL, INTENT(IN) :: npmsl_height       ! Orographic height above which
                                       ! relaxation occurs

! --------------------------------------------------
! In the following variable definitions,
!   thetas = theta at Surface
!   thetam = theta at Mean Sea Level
! The same naming convention applies to t and exner
! --------------------------------------------------

! Define variables local to PE
INTEGER, INTENT(IN) :: row_length_local    ! Number of points on a row
INTEGER, INTENT(IN) :: rows_local          ! Number of rows of data

REAL, INTENT(IN)    :: thetas_local( 1-off_x:row_length_local+off_x            &
                                   , 1-off_y:rows_local+off_y )

REAL, INTENT(IN)    :: pstar_local          (row_length_local, rows_local)
REAL, INTENT(IN)    :: phi_star_local       (row_length_local, rows_local)
REAL, INTENT(IN)    :: tm_local             (row_length_local, rows_local)
REAL, INTENT(IN)    :: cos_p_latitude_local (row_length_local, rows_local)
REAL, INTENT(INOUT) :: pmsl_local           (row_length_local, rows_local)


! Define PE0's variables
INTEGER, INTENT(IN) :: row_length    ! Number of points on a row
INTEGER, INTENT(IN) :: rows          ! Number of rows of data

REAL :: pmsl     (row_length, rows)
REAL :: pstar    (row_length, rows)
REAL :: phi_star (row_length, rows)
REAL :: thetas   (row_length, rows)
REAL :: tm       (row_length, rows)
REAL :: cos_p_latitude (row_length, rows)

! Define local variables
REAL :: exnert_local(row_length_local, rows_local)
REAL :: exners  (row_length, rows)
REAL :: exnerm  (row_length, rows)
REAL :: exnert  (row_length, rows, 2)  ! New exner pressure
REAL :: avexner 
REAL :: orog    (row_length, rows)
REAL :: thetam  (row_length, rows)
REAL :: fug     
REAL :: fvg    
REAL :: f       (row_length, rows)
REAL :: sphlo1  (row_length, rows)
REAL :: sphlo2  (row_length, rows)

REAL :: sphla1         ! Latitundinal spherical term
REAL :: sphla2         ! Latitundinal spherical term squared
REAL :: orr            ! Over relaxation factor
REAL :: oorr           ! Over relaxation factor at previous iteration
REAL :: rho            ! Spherical radius of convergence
REAL :: n              ! Typical number of points above npmsl_height


INTEGER ::  niter      ! Number of iterations
INTEGER ::  i,j,kk,mm
INTEGER ::  ip1,im1,jp1,jm1


! Workspace usage
REAL :: aa  (row_length, rows)
REAL :: bb  (row_length, rows)
REAL :: pre (row_length, rows)
REAL :: fac

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_NPMSL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Constants
n     = 15
rho   = COS(pi/n)
niter = 20


! DEPENDS ON: gather_field
CALL gather_field                                                              &
  ( pstar_local, pstar                                                         &
  , row_length_local, rows_local                                               &
  , row_length, rows                                                           &
  , fld_type_p, halo_type_no_halo                                              &
  , 0, all_proc_group )

! DEPENDS ON: gather_field
CALL gather_field                                                              &
  ( pmsl_local, pmsl                                                           &
  , row_length_local, rows_local                                               &
  , row_length, rows                                                           &
  , fld_type_p, halo_type_no_halo                                              &
  , 0, all_proc_group )

! DEPENDS ON: gather_field
CALL gather_field                                                              &
  ( thetas_local, thetas                                                       &
  , row_length_local+2*off_x, rows_local+2*off_y                               &
  , row_length, rows                                                           &
  , fld_type_p, halo_type_single                                               &
  , 0, all_proc_group )

! DEPENDS ON: gather_field
CALL gather_field                                                              &
  ( tm_local, tm                                                               &
  , row_length_local, rows_local                                               &
  , row_length, rows                                                           &
  , fld_type_p, halo_type_no_halo                                              &
  , 0, all_proc_group )

! DEPENDS ON: gather_field
CALL gather_field                                                              &
  ( phi_star_local, phi_star                                                   &
  , row_length_local, rows_local                                               &
  , row_length, rows                                                           &
  , fld_type_p, halo_type_no_halo                                              &
  , 0, all_proc_group )

! DEPENDS ON: gather_field
CALL gather_field                                                              &
  ( cos_p_latitude_local, cos_p_latitude                                       &
  , row_length_local, rows_local                                               &
  , row_length, rows                                                           &
  , fld_type_p, halo_type_no_halo                                              &
  , 0, all_proc_group )

! Calculate new pmsl on PE0
IF (me == 0) THEN

  sphla1 = 1.0/(delta_phi*planet_radius)
  sphla2 = (delta_phi*planet_radius)**2
! OpenMP around this loop changes pmsl values for two tests when using cray
! compiler version 2.2.1 as it changes the order of calculations for sphlo2.
  DO j=1, rows
    DO i=1, row_length
      f(i,j)      = 0.0
      exners(i,j) = (pstar(i,j)/p_zero)**kappa
      exnerm(i,j) = (pmsl(i,j)/p_zero)**kappa
      orog(i,j)   = phi_star(i,j)/g
      thetam(i,j) = tm(i,j)/exnerm(i,j)
      ! The check below only applies very close to the poles and the values
      ! of the two variables set to arbritrary numbers don't get used in the
      ! correction step
      IF (ABS(cos_p_latitude(i,j)) < 1.0e-8) THEN
        sphlo1(i,j) = 1.0e8
        sphlo2(i,j) = 0.0
      ELSE
        sphlo1(i,j) = 1.0/(delta_lambda*planet_radius*cos_p_latitude(i,j))
        sphlo2(i,j) = (delta_lambda*planet_radius*cos_p_latitude(i,j))**2
      END IF

      pre(i,j)      = (sphlo2(i,j)*sphla2)/((2*sphlo2(i,j))+(2*sphla2))
      exnert(i,j,1) = exnerm(i,j)
      exnert(i,j,2) = exnerm(i,j)
    END DO
  END DO

  ! Calculate geostrophic wind
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,ip1,im1,fac,jp1,jm1,fug,fvg)          &
!$OMP SHARED(rows,row_length,model_type,sphlo1,phi_star,cp,thetas,exners,      &
!$OMP sphla1,aa,bb,thetam,f)
!$OMP DO SCHEDULE(STATIC)
  DO j=1, rows
    DO i=1, row_length
      IF (i == 1) THEN
        ip1 = i+1
        IF (model_type == mt_lam) THEN
          im1 = i
          fac = 1.0
        ELSE
          im1 = row_length
          fac = 0.5
        END IF

      ELSE IF (i == row_length) THEN
        IF (model_type == mt_lam) THEN
          ip1 = i
          fac = 1.0
        ELSE
          ip1 = 1
          fac = 0.5
        END IF
        im1 = i-1

      ELSE
        ip1 = i+1
        im1 = i-1
        fac = 0.5
      END IF

      fvg = fac*sphlo1(i,j)*((phi_star(ip1,j)-phi_star(im1,j))            &
               + cp*thetas(i,j)*(exners(ip1,j)-exners(im1,j)))

      IF (j == 1) THEN
        jp1 = j+1
        jm1 = j
        fac = 1.0
      ELSE IF (j == rows) THEN
        jp1 = j
        jm1 = j-1
        fac = 1.0
      ELSE
        jp1 = j+1
        jm1 = j-1
        fac = 0.5
      END IF

      fug = -fac*sphla1*((phi_star(i,jp1)-phi_star(i,jm1))                &
               + cp*thetas(i,j)*(exners(i,jp1)-exners(i,jm1)))


  ! Calculate RHS of equation F

      aa(i,j) = fvg/thetam(i,j)
      bb(i,j) = fug/thetam(i,j)
    END DO
  END DO
!$OMP END DO

  IF (model_type == mt_lam) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j=2, rows-1
      DO i=2, row_length-1
        f(i,j) = (0.5/cp)*sphlo1(i,j)*(aa(i+1,j)-aa(i-1,j))                    &
               - (0.5/cp)*sphla1*(bb(i,j+1)-bb(i,j-1))
      END DO
    END DO
!$OMP END DO
  ELSE
!$OMP DO SCHEDULE(STATIC)
    DO j=2, rows-1
      DO i=1, row_length
        IF (i == 1) THEN
          ip1 = i+1
          im1 = row_length
        ELSE IF (i == row_length) THEN
          ip1 = 1
          im1 = i-1
        ELSE
          ip1 = i+1
          im1 = i-1
        END IF
        f(i,j) = (0.5/cp)*sphlo1(i,j)*(aa(ip1,j)-aa(im1,j))                    &
               - (0.5/cp)*sphla1*(bb(i,j+1)-bb(i,j-1))
      END DO
    END DO
!$OMP END DO

  END IF
!$OMP END PARALLEL

  ! Relaxation to find exner twiddle
  DO kk=1, niter
    IF (kk == 1) THEN
      orr = 1.0
    ELSE IF (kk == 2) THEN
      orr = 1.0/(1.0-((rho**2)/2.0))
    ELSE
      orr = 1.0/(1.0-((oorr*rho**2)/4.0))
    END IF
    IF (model_type == mt_lam) THEN
      DO j=3, rows-2
        DO i=3, row_length-2
          IF (orog(i,j) > npmsl_height) THEN
            avexner       = pre(i,j)*((1.0/sphla2*(exnert(i+1,j,1)             &
                          + exnert(i-1,j,2)))                                  &
                          + (1.0/sphlo2(i,j)*(exnert(i,j+1,1)                  &
                          + exnert(i,j-1,2)))-f(i,j))

            exnert(i,j,2) = exnert(i,j,1)+orr*(avexner - exnert(i,j,1))
          END IF
        END DO
      END DO
    ELSE
      DO j=3, rows-2
        DO i=1, row_length
          IF (orog(i,j) >  npmsl_height) THEN
            IF (i == 1) THEN
              ip1 = i+1
              im1 = row_length
              mm  = 1
            ELSE IF (i == row_length) THEN
              ip1 = 1
              im1 = i-1
              mm  = 2
            ELSE
              ip1 = i+1
              im1 = i-1
              mm  = 2
            END IF
            avexner       = pre(i,j)*((1.0/sphla2*(exnert(ip1,j,1)             &
                          + exnert(im1,j,mm)))                                 &
                          + (1.0/sphlo2(i,j)*(exnert(i,j+1,1)                  &
                          + exnert(i,j-1,2)))-f(i,j))

            exnert(i,j,2) = exnert(i,j,1)+orr*(avexner - exnert(i,j,1))
          END IF
        END DO
      END DO
    END IF

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j)                  &
!$OMP SHARED(rows,row_length,exnert)
    DO j=1, rows
      DO i=1, row_length
        exnert(i,j,1) = exnert(i,j,2)
      END DO
    END DO
!$OMP END PARALLEL DO
    oorr = orr

  END DO

END IF  ! me == 0


! Distribute from exnert on PE0 to exnert_local
! DEPENDS ON: scatter_field
CALL scatter_field                                                             &
  ( exnert_local, exnert                                                       &
  , row_length_local, rows_local                                               &
  , row_length, rows                                                           &
  , fld_type_p, halo_type_no_halo                                              &
  , 0, all_proc_group )


! Calculate new pmsl
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j) SHARED(rows_local,                &
!$OMP pmsl_local,p_zero,exnert_local,recip_kappa,row_length_local)
DO j=1, rows_local
  DO i=1, row_length_local
    pmsl_local(i,j) = p_zero*(exnert_local(i,j)**(recip_kappa))
  END DO
END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE calc_npmsl

END MODULE calc_npmsl_mod
