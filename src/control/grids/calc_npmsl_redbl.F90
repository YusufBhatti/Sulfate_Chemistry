! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine Calc_NPMSL_Redbl
MODULE calc_npmsl_redbl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CALC_NPMSL_REDBL_MOD'

CONTAINS

SUBROUTINE calc_npmsl_redbl                                                    &
  ( pmsl, pstar, phi_star_local, thetas, tm, cos_p_latitude, delta_lambda      &
  , delta_phi, row_length, rows, global_row_length, global_rows                &
  , me, n_proc, n_procx, n_procy, off_x, off_y, halo_i, halo_j, neighbour      &
  , at_extremity, all_proc_group, npmsl_height )

!-------------------------------------------------------------


! Purpose:
!          Modifies Pressure at mean sea level as calculated
!          in Calc_PMSL to account for errors due to high
!          orography
!
! Method:
!          Extrapolates as in UMDP 80 but revised smoothing algorithm.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Grids
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards. v10.3

USE swap_bounds_2d_mv_mod, ONLY: swap_bounds_2d_mv
USE mpp_conf_mod,          ONLY: swap_field_is_scalar
USE planet_constants_mod, ONLY: g, planet_radius, &
                                cp, kappa, p_zero, recip_kappa
USE swapable_field_mod,  ONLY: swapable_field_pointer_type
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE UM_ParParams

USE model_domain_mod, ONLY: model_type, mt_lam

IMPLICIT NONE

INTEGER, INTENT(IN) :: off_x
INTEGER, INTENT(IN) :: off_y
INTEGER, INTENT(IN) :: halo_i
INTEGER, INTENT(IN) :: halo_j
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
!  thetas = theta at Surface
!  thetam = theta at Mean Sea Level
! The same naming convention applies to t and exner
! --------------------------------------------------

!  Define variables local to PE


INTEGER, INTENT(IN) :: row_length        ! Number of points on a row
INTEGER, INTENT(IN) :: rows              ! Number of rows of data
INTEGER, INTENT(IN) :: global_row_length ! Number of points on a row
INTEGER, INTENT(IN) :: global_rows       ! Number of rows of data

!REAL         :: field_global(global_row_length, global_rows)
REAL, INTENT(IN)    :: pstar (row_length, rows)
REAL, INTENT(IN)    :: tm    (row_length, rows)

REAL, INTENT(IN)    :: phi_star_local (row_length, rows)
REAL, INTENT(IN)    :: cos_p_latitude (row_length, rows)

REAL, INTENT(IN)    :: thetas( 1-off_x:row_length+off_x                        &
                             , 1-off_y:rows+off_y)

REAL, INTENT(INOUT) :: pmsl  (row_length, rows)

! Define local variables
REAL :: orog   
REAL :: thetam (row_length, rows)
REAL :: fug   
REAL :: fvg  
REAL :: sphlo1 (row_length, rows)
REAL :: exnerm (row_length, rows)

REAL         :: exnert   (1-off_x:row_length+off_x, 1-off_y:rows+off_y)
REAL, TARGET :: phi_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y)
REAL, TARGET :: exners   (1-off_x:row_length+off_x, 1-off_y:rows+off_y)

REAL, TARGET, ALLOCATABLE :: avexner(:,:)
REAL, TARGET, ALLOCATABLE :: sphlo2(:,:)
REAL,         ALLOCATABLE :: f(:,:)
REAL,         ALLOCATABLE :: exnert_big(:,:)

REAL ::  sphla1            ! Latitundinal spherical term
REAL ::  sphla2            ! Latitundinal spherical term squared

INTEGER :: niter           ! Number of iterations
INTEGER :: small           ! Extended halo size for sweeping
INTEGER :: loop
INTEGER :: counter
INTEGER :: index_n
INTEGER :: index_s
INTEGER :: index_w
INTEGER :: index_e


INTEGER :: i, j, kk, mm
INTEGER :: ip1, im1, jp1, jm1
INTEGER :: j_start, j_stop, i_start, i_stop

LOGICAL, ALLOCATABLE :: check(:,:)


! Workspace usage

REAL, TARGET :: aa (1-off_x:row_length+off_x, 1-off_y:rows+off_y)
REAL, TARGET :: bb (1-off_x:row_length+off_x, 1-off_y:rows+off_y)
REAL, TARGET, ALLOCATABLE :: pre(:,:)
REAL ::  fac

! Variables for multivariate swapbounds
INTEGER :: i_field
TYPE(swapable_field_pointer_type) :: fields_to_swap(3)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_NPMSL_REDBL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Constants
! Set extended sweeping halo - make sure choice is valid
small = MIN(halo_i,halo_j)

ALLOCATE ( check   (1-small:row_length+small, 1-small:rows+small) )
ALLOCATE ( pre     (1-small:row_length+small, 1-small:rows+small) )
ALLOCATE ( avexner (1-small:row_length+small, 1-small:rows+small) )
ALLOCATE ( sphlo2  (1-small:row_length+small, 1-small:rows+small) )

niter  = global_row_length
sphla1 = 1.0/(delta_phi*planet_radius)
sphla2 = (delta_phi*planet_radius)**2

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,orog)       &
!$OMP SHARED(exnerm)                                                     &
!$OMP SHARED(rows,row_length,pstar,pmsl,phi_star_local,tm,exners,p_zero, &
!$OMP phi_star,thetam,sphlo1,sphlo2,pre,exnert,avexner,delta_lambda,g,   &
!$OMP planet_radius,cos_p_latitude,check,kappa,npmsl_height,sphla2)
DO j=1, rows
  DO i=1, row_length

    exners(i,j) = (pstar(i,j) / p_zero)**kappa
    exnerm(i,j) = (pmsl(i,j)  / p_zero)**kappa
    orog        = phi_star_local(i,j) / g
    check(i,j)  = .FALSE.

    IF ( orog > npmsl_height ) check(i,j) = .TRUE.

    phi_star(i,j) = phi_star_local(i,j)
    thetam(i,j)   = tm(i,j) / exnerm(i,j)

    ! The check below only applies very close to the poles and the values
    ! of the two variables set to arbitrary numbers don't get used in the
    ! correction step
    IF ( ABS(cos_p_latitude(i,j)) < 1.0e-8 ) THEN
      sphlo1(i,j) = 1.0e8
      sphlo2(i,j) = 0.0
    ELSE
      sphlo1(i,j) = 1.0/(delta_lambda*planet_radius*cos_p_latitude(i,j))
      sphlo2(i,j) = (delta_lambda*planet_radius*cos_p_latitude(i,j))**2
    END IF
    pre(i,j)     = (sphlo2(i,j)*sphla2) / ((2*sphlo2(i,j))+(2*sphla2))
    exnert(i,j)  = exnerm(i,j)
    avexner(i,j) = exnerm(i,j)
  END DO
END DO
!$OMP END PARALLEL DO

i_field = 1
fields_to_swap(i_field) % field_2d    => avexner(:,:)
fields_to_swap(i_field) % field_type  =  fld_type_p
fields_to_swap(i_field) % levels      =  1
fields_to_swap(i_field) % rows        =  rows
fields_to_swap(i_field) % vector      =  .FALSE.
i_field = i_field + 1
fields_to_swap(i_field) % field_2d    => pre(:,:)
fields_to_swap(i_field) % field_type  =  fld_type_p
fields_to_swap(i_field) % levels      =  1
fields_to_swap(i_field) % rows        =  rows
fields_to_swap(i_field) % vector      =  .FALSE.
i_field = i_field + 1
fields_to_swap(i_field) % field_2d    => sphlo2(:,:)
fields_to_swap(i_field) % field_type  =  fld_type_p
fields_to_swap(i_field) % levels      =  1
fields_to_swap(i_field) % rows        =  rows
fields_to_swap(i_field) % vector      =  .FALSE.

CALL swap_bounds_2d_mv(fields_to_swap, i_field, row_length,  &
                       small, small )

! DEPENDS ON: swap_bounds
CALL swap_bounds (check, row_length, rows, 1, small, small,  &
                  fld_type_p, swap_field_is_scalar)

i_field = 1
fields_to_swap(i_field) % field_2d    => phi_star(:,:)
fields_to_swap(i_field) % field_type  =  fld_type_p
fields_to_swap(i_field) % levels      =  1
fields_to_swap(i_field) % rows        =  rows
fields_to_swap(i_field) % vector      =  .FALSE.
i_field = i_field + 1
fields_to_swap(i_field) % field_2d    => exners(:,:)
fields_to_swap(i_field) % field_type  =  fld_type_p
fields_to_swap(i_field) % levels      =  1
fields_to_swap(i_field) % rows        =  rows
fields_to_swap(i_field) % vector      =  .FALSE.

CALL swap_bounds_2d_mv(fields_to_swap, i_field, row_length,  &
                       off_x,off_y)

! Calculate geostrophic wind
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(j,i,ip1,im1,jp1,jm1,fac,fug,fvg) &
!$OMP SHARED(rows,row_length,at_extremity,model_type,sphlo1,phi_star,    &
!$OMP cp,thetas,exners,sphla1,aa,bb,thetam) SCHEDULE(STATIC)
DO j=1, rows
  DO i=1, row_length
    IF (i == 1 .AND. at_extremity(pwest)) THEN
      ip1 = i+1
      IF (model_type == mt_lam) THEN
        im1 = i
        fac = 1.0
      ELSE
        im1 = i-1
        fac = 0.5
      END IF

    ELSE IF (i == row_length .AND. at_extremity(peast)) THEN
      IF (model_type == mt_lam) THEN
        ip1 = i
        fac = 1.0
      ELSE
        ip1 = i+1
        fac = 0.5
      END IF
      im1 = i-1

    ELSE
      ip1 = i+1
      im1 = i-1
      fac = 0.5
    END IF

    fvg = fac*sphlo1(i,j)*((phi_star(ip1,j)-phi_star(im1,j))              &
             + cp*thetas(i,j)*(exners(ip1,j)-exners(im1,j)))

    IF (j == 1 .AND. at_extremity(psouth)) THEN
      jp1 = j+1
      jm1 = j
      fac = 1.0
    ELSE IF (j == rows .AND. at_extremity(pnorth)) THEN
      jp1 = j
      jm1 = j-1
      fac = 1.0
    ELSE
      jp1 = j+1
      jm1 = j-1
      fac = 0.5
    END IF

    fug = -fac*sphla1*((phi_star(i,jp1)-phi_star(i,jm1))                  &
             + cp*thetas(i,j)*(exners(i,jp1)-exners(i,jm1)))

! Calculate RHS of equation F

    aa(i,j) = fvg / thetam(i,j)
    bb(i,j) = fug / thetam(i,j)

  END DO
END DO
!$OMP END PARALLEL DO

i_field = 1
fields_to_swap(i_field) % field_2d    => aa(:,:)
fields_to_swap(i_field) % field_type  =  fld_type_p
fields_to_swap(i_field) % levels      =  1
fields_to_swap(i_field) % rows        =  rows
fields_to_swap(i_field) % vector      =  .FALSE.
i_field = i_field + 1
fields_to_swap(i_field) % field_2d    => bb(:,:)
fields_to_swap(i_field) % field_type  =  fld_type_p
fields_to_swap(i_field) % levels      =  1
fields_to_swap(i_field) % rows        =  rows
fields_to_swap(i_field) % vector      =  .FALSE.

CALL swap_bounds_2d_mv(fields_to_swap, i_field, row_length,  &
                       off_x,off_y)

IF (at_extremity(PSouth)) THEN
  j_start = 2
ELSE
  j_start = 1
END IF
IF (at_extremity(PNorth)) THEN
  j_stop = rows-1
ELSE
  j_stop = rows
END IF

i_start = 1
i_stop  = row_length
IF (model_type == mt_lam .AND. at_extremity(Pwest)) THEN
  i_start = 2
END IF
IF (model_type == mt_lam .AND. at_extremity(Peast)) THEN
  i_stop = row_length-1
END IF

ALLOCATE ( f(1-small:row_length+small,1-small:rows+small) )

IF (model_type == mt_lam) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(j,i)                  &
!$OMP SHARED(j_stop,i_stop,i_start,j_start,sphlo1,sphla1,aa,bb,f,cp)
  DO j=j_start, j_stop
    DO i=i_start, i_stop
      f(i,j) = (0.5/cp)*sphlo1(i,j)*(aa(i+1,j)-aa(i-1,j))                      &
             - (0.5/cp)*sphla1*(bb(i,j+1)-bb(i,j-1))
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(j,i,ip1,im1)          &
!$OMP SHARED(j_stop,i_stop,i_start,j_start,sphlo1,sphla1,aa,bb,f,cp)
  DO j=j_start, j_stop
    DO i=i_start, i_stop
      ip1 = i+1
      im1 = i-1
      f(i,j) = (0.5/cp)*sphlo1(i,j)*(aa(ip1,j)-aa(im1,j))                      &
             - (0.5/cp)*sphla1*(bb(i,j+1)-bb(i,j-1))
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! DEPENDS ON: swap_bounds
CALL swap_bounds (f, row_length, rows, 1, small, small,                        &
                  fld_type_p, swap_field_is_scalar)

IF (at_extremity(PSouth)) THEN
  j_start = 3
ELSE
  j_start = 1
END IF
IF (at_extremity(PNorth)) THEN
  j_stop = rows-2
ELSE
  j_stop = rows
END IF

i_start = 1
i_stop  = row_length
IF (model_type == mt_lam .AND. at_extremity(Pwest)) i_start = 3
IF (model_type == mt_lam .AND. at_extremity(Peast)) i_stop = row_length-2

! Relaxation to find exner twiddle
ALLOCATE ( exnert_big(1-small:row_length+small, 1-small:rows+small) )

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(j,i)                  &
!$OMP SHARED(rows,row_length,exnert_big,exnert)
DO j=1, rows
  DO i=1, row_length
    exnert_big(i,j) = exnert(i,j)
  END DO
END DO
!$OMP END PARALLEL DO

counter = 0

outside: DO loop=1, ((niter/small) +1)

  ! DEPENDS ON: swap_bounds
  CALL swap_bounds ( exnert_big, row_length, rows, 1, small                    &
                   , small, fld_type_p, swap_field_is_scalar )

  ! Loop over the sweeping area
  DO kk=1, small

    ! Reset indices each time we step in
    IF (at_extremity(PSouth)) THEN
      index_s = j_start
    ELSE
      index_s = j_start - small + kk
    END IF

    IF (at_extremity(PNorth)) THEN
      index_n = j_stop
    ELSE
      index_n = j_stop + small - kk
    END IF

    index_w = i_start - small + kk
    index_e = i_stop  + small - kk

    IF (model_type == mt_lam .AND. at_extremity(Pwest)) index_w = i_start
    IF (model_type == mt_lam .AND. at_extremity(Peast)) index_e = i_stop

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j) SHARED(index_e,index_w,      &
!$OMP index_n,index_s,check,pre,sphla2,exnert_big,sphlo2,f,avexner)
!$OMP DO SCHEDULE(STATIC)
    DO j=index_s, index_n
      DO i=index_w, index_e
        IF (check(i,j)) THEN
          avexner(i,j) = pre(i,j)*((1.0/sphla2*(exnert_big(i+1,j)              &
                       + exnert_big(i-1,j)))                                   &
                       + (1.0/sphlo2(i,j)*(exnert_big(i,j+1)                   &
                       + exnert_big(i,j-1)))-f(i,j))
        END IF
      END DO
    END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
    DO j=index_s, index_n
      DO i=index_w, index_e
        IF (check(i,j)) THEN
          exnert_big(i,j) = avexner(i,j)
        END IF
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    counter = counter + 1
    IF (counter == niter) EXIT outside
  END DO  ! kk=1,small
END DO outside ! loop


! Calculate new pmsl
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(j,i)                 &
!$OMP SHARED(rows,row_length,pmsl,exnert_big,p_zero,recip_kappa)
DO j=1, rows
  DO i=1, row_length
    pmsl(i,j) = p_zero * (exnert_big(i,j)**(recip_kappa))
  END DO
END DO
!$OMP END PARALLEL DO

DEALLOCATE(exnert_big)
DEALLOCATE(f)
DEALLOCATE(sphlo2)
DEALLOCATE(avexner)
DEALLOCATE(pre)
DEALLOCATE(check)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE calc_npmsl_redbl

END MODULE calc_npmsl_redbl_mod
