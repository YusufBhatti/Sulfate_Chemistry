! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Set eta_theta_levels and eta_rho_levels using the recon_idealised namelist

MODULE rcf_ideal_set_eta_levels_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_IDEAL_SET_ETA_LEVELS_MOD'

CONTAINS

! Subroutine rcf_ideal_set_eta_levels
!
! Description:
!   Set/overwrite eta_theta_levels and eta_rho_levels according to the
!   recon_idealised namelist choices.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.

SUBROUTINE rcf_ideal_set_eta_levels (model_levels, eta_theta_levels,   &
                                     eta_rho_levels)

USE umPrintMgr,   ONLY: &
    newline,            &
    PrintStatus,        &
    PrStatus_Normal,    &
    PrStatus_Oper,      &
    PrStatus_Diag,      &
    umMessage,          &
    umPrint

USE UM_ParCore,   ONLY: &
    mype

USE rcf_nlist_recon_idealised_mod, ONLY: &
    big_factor,                          &
    big_layers,                          &
    first_theta_height,                  &
    grid_number,                         &
    height_domain,                       &
    thin_theta_height,                   &
    transit_layers,                      &
    vert_grid_ratio
    
USE rcf_ideal_vgrid_mod, ONLY: &
    vert_regular,              &
    vert_quadratic_theta,      &
    vert_bi_quadratic,         &
    vert_quadratic_uvtheta,    &
    vert_schar,                &
    vert_dwd,                  &
    vert_stretch_plus_regular, &
    vert_quad_stretch_thin,    &
    vert_geometric_theta,      &
    vert_dump,                 &
    vert_idl_um_grid

USE ereport_mod, ONLY: &
    ereport

USE yomhook,   ONLY: &
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

INTEGER, INTENT(IN) :: model_levels    ! Output grid model_levels
REAL, POINTER :: eta_theta_levels(:)   ! Output grid eta_theta_levels
REAL, POINTER :: eta_rho_levels(:)     ! Output grid eta_rho_levels

! Local variables
INTEGER :: k
INTEGER :: icode            ! Error code
INTEGER :: regular_layers   ! No of regular (not stretched) levels
INTEGER :: crit_level       ! Level at which magnification is below some value
REAL    :: mag_factor       ! Magnification factor for transition region from 
                            ! regular levels to big levels
REAL    :: mag1             ! Temp. value used with mag_value1 and mag_factor1
REAL    :: mag3             ! Temp. value used with mag_value3 and mag_factor3
REAL    :: delta_eta        ! Difference between levels
REAL    :: eta_ratio        ! Ratio between successive eta levels
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'RCF_IDEAL_SET_ETA_LEVELS'

! Local parameters for idealised configurations
REAL, PARAMETER :: mag_value1 = 1.04     ! value for 38 levels
REAL, PARAMETER :: mag_value2 = 1.1001   ! value for 38 levels
REAL, PARAMETER :: mag_value3 = 1.1      ! value for 38 levels
REAL, PARAMETER :: top_mag_value = 1.5   ! value for 38 levels
REAL, PARAMETER :: mag_factor1 = 0.01    ! value for 38 levels
REAL, PARAMETER :: mag_factor2 = 0.02    ! value for 38 levels
REAL, PARAMETER :: mag_factor3 = 0.02    ! value for 38 levels

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!---------------------------------
! Set eta_theta and eta_rho levels
!---------------------------------

SELECT CASE(grid_number)

CASE(vert_regular)
! Regular grid

  delta_eta = 1.0/model_levels
  eta_theta_levels(0) = 0.0
  DO k = 1, model_levels
    eta_theta_levels(k) = k * delta_eta
    eta_rho_levels(k)   = 0.5 * ( eta_theta_levels(k-1) + eta_theta_levels(k) )
  END DO

  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
    CALL umPrint('Regular vertical grid',src=RoutineName)
    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE(umMessage,'(2(A,E16.8))')                                      &
     ' eta_theta_levels(model_levels): ', eta_theta_levels(model_levels),  &
     ', height_domain (metres): ', height_domain
      CALL umPrint(umMessage, src=RoutineName)
      WRITE(umMessage,'(A,E16.8)') ' Layer thickness: ',height_domain*delta_eta
      CALL umPrint(umMessage, src=RoutineName)
    END IF
  END IF


CASE(vert_quadratic_theta)
! Quadratic grid

  delta_eta = 1.0/model_levels
  eta_theta_levels(0) = 0.0
  DO k = 1, model_levels
    eta_theta_levels(k) = (k * delta_eta) * (k * delta_eta)
    eta_rho_levels(k)   = 0.5 * ( eta_theta_levels(k-1) + eta_theta_levels(k) )
  END DO

  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
    CALL umPrint('Quadratic vertical grid',src=RoutineName)
    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE(umMessage,'(2(A,E16.8))')                                      &
     ' eta_theta_levels(model_levels): ', eta_theta_levels(model_levels),  &
     ', height_domain (metres): ', height_domain
      CALL umPrint(umMessage, src=RoutineName)
    END IF
  END IF

CASE(vert_bi_quadratic)
! Quadratic in theta only (no modification of rho levels)

  delta_eta = 1.0/model_levels
  eta_theta_levels(0) = 0.0
  DO k = 1, model_levels
    eta_theta_levels(k) = (k * delta_eta) * (k * delta_eta)
  END DO

  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
    CALL umPrint('Quadratic vertical grid for u and theta levels', &
                 src=RoutineName)
    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE(umMessage,'(2(A,E16.8))')                                      &
     ' eta_theta_levels(model_levels): ', eta_theta_levels(model_levels),  &
     ', height_domain (metres): ', height_domain
      CALL umPrint(umMessage, src=RoutineName)
    END IF
  END IF


CASE(vert_quadratic_uvtheta)
! Quadratic in both theta and rho

  delta_eta = 1.0/(2.0*model_levels)
  eta_theta_levels(0) = 0.0
  DO k = 1, model_levels
    eta_theta_levels(k) = (2.0 * k * delta_eta) * (2.0 * k * delta_eta)
    eta_rho_levels(k)   = (((2*k)-1.0) * delta_eta) * (((2*k)-1.0) * delta_eta)
  END DO

  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
    CALL umPrint('Common quadratic vertical grid for u and theta levels', &
                 src=RoutineName)
    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE(umMessage,'(2(A,E16.8))')                                      &
     ' eta_theta_levels(model_levels): ', eta_theta_levels(model_levels),  &
     ', height_domain (metres): ', height_domain
      CALL umPrint(umMessage, src=RoutineName)
    END IF
  END IF


CASE(vert_schar)
! Schar grid - eta_theta_levels are approx. those used by Schar.

  eta_theta_levels(0) = 0.0
  DO k = 1, model_levels
    eta_theta_levels(k) = k**1.2 / model_levels**1.2
    eta_rho_levels(k)   = 0.5 * (eta_theta_levels(k-1) + eta_theta_levels(k))
  END DO

  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
    CALL umPrint('Schar vertical grid',src=RoutineName)
    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE(umMessage,'(2(A,E16.8))')                                      &
     ' eta_theta_levels(model_levels): ', eta_theta_levels(model_levels),  &
     ', height_domain (metres): ', height_domain
      CALL umPrint(umMessage, src=RoutineName)
    END IF
  END IF


CASE(vert_dwd)
! DWD stretched grid - eta_theta_levels are similar to DWD.

  ! Set first levels to actual height values:
  eta_theta_levels(0) = 0.0
  eta_theta_levels(1) = 40.0
  eta_theta_levels(2) = 100.0
  delta_eta = 100.0
  ! Remaining levels equate to a stretched grid
  DO k = 3, model_levels
    eta_theta_levels(k) = eta_theta_levels(k-1) + delta_eta
    delta_eta = delta_eta + 40.0
  END DO
  ! Namelist input must be overridden:
  height_domain = eta_theta_levels(model_levels)

  ! Normalise theta levels and set rho levels:
  DO k = 1, model_levels
    eta_theta_levels(k) = eta_theta_levels(k) / eta_theta_levels(model_levels)
    eta_rho_levels(k)   = 0.5 * (eta_theta_levels(k-1) + eta_theta_levels(k))
  END DO

  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
    CALL umPrint('DWD stretched vertical grid',src=RoutineName)
    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE(umMessage,'(2(A,E16.8))')                                      &
     ' eta_theta_levels(model_levels): ', eta_theta_levels(model_levels),  &
     ', height_domain (metres): ', height_domain
      CALL umPrint(umMessage, src=RoutineName)
    END IF
  END IF


CASE(vert_stretch_plus_regular)
! Regular grid followed by stretched grid (of big_layers levels)

  ! Number of regular layers:
  regular_layers = model_levels - transit_layers - big_layers
  IF (regular_layers < 0) THEN
    icode = ABS(regular_layers)
    WRITE(umMessage, '(3(A,I0),A)') 'Size of stretched grid (big_layers (',   &
      big_layers, ') + transit_layers (', transit_layers, '))' // newline //  &
      'is greater than the total number of model levels (',model_levels, ')'
    CALL ereport(RoutineName, icode, umMessage)
  END IF

  ! magnification factor is applied over transit_layers levels:
  mag_factor = big_factor**1.0/MAX(1,transit_layers)

  ! eta levels interval is chosen to scale to 1 over the domain
  ! NB  mag_factor * (1.0 - big_factor) / (1.0 - mag_factor) is the
  ! first transit layer magnified
  delta_eta = 1.0 / (regular_layers +                                  &
                    (1.0 - big_factor) / (1.0 - mag_factor) +          &
                    big_factor * big_layers)

  ! Set height domain to be the full height domain
  ! (the input value height_domain is for regular levels)
  height_domain = height_domain / (regular_layers * delta_eta)

  eta_theta_levels(0) = 0.0
  ! Regular portion of grid
  DO k = 1, regular_layers
    eta_theta_levels(k) = k * delta_eta
  END DO
  ! Transition zone between regular and fully stretched grid
  DO k = regular_layers + 1, regular_layers + transit_layers
    eta_theta_levels(k) = eta_theta_levels(k-1) + delta_eta
    delta_eta = delta_eta * mag_factor
  END DO
  ! Fully stretched grid
  DO k = regular_layers + transit_layers + 1, model_levels
    eta_theta_levels(k) = eta_theta_levels(k-1) + delta_eta
  END DO
  ! Rho levels
  DO k = 1, model_levels
    eta_rho_levels(k) = 0.5 * (eta_theta_levels(k-1) + eta_theta_levels(k))
  END DO

  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
    CALL umPrint('Regular then stretched vertical grid',src=RoutineName)
    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE(umMessage,'(A,E16.8,A,I0,A)') ' regular grid: ',               &
       height_domain*delta_eta, ' metres, ', regular_layers, ' layers'
      CALL umPrint(umMessage, src=RoutineName)
      WRITE(umMessage,'(A,I0,A,E16.8)') ' stretching over ',               &
        transit_layers, ' levels, magnification factor: ',mag_factor
      CALL umPrint(umMessage, src=RoutineName)
      WRITE(umMessage,'(A,E12.5,A,I0,A)') 'stretched grid: ', big_factor,  &
        ' * regular depth over ', big_layers, ' layers'
      CALL umPrint(umMessage, src=RoutineName)
      WRITE(umMessage,'(2(A,E16.8))')                                      &
     ' eta_theta_levels(model_levels): ', eta_theta_levels(model_levels),  &
     ', height_domain (metres): ', height_domain
      CALL umPrint(umMessage, src=RoutineName)
    END IF
  END IF


CASE(vert_quad_stretch_thin)
! Quadratic stretched grid, thin near surface
! Hardwired assumption of 38 levels in the parameters used.

  delta_eta = 1.0/(model_levels - 1)
  eta_theta_levels(0) = 0.0
  ! Quadratic grid:
  DO k = 1, model_levels -1
    eta_theta_levels(k) = (k * delta_eta) * (k * delta_eta)
  END DO

  ! Calculate ratio of successive eta differences
  ! Store the first level at which the ratio is less than some value
  crit_level = model_levels
  DO k = 1, model_levels-2
    eta_ratio = ( eta_theta_levels(k+1) - eta_theta_levels(k) )        &
              / ( eta_theta_levels(k)   - eta_theta_levels(k-1) )
    IF (eta_ratio < mag_value1) THEN
      crit_level = k
      EXIT
    END IF
  END DO

  IF (crit_level /= model_levels) THEN
    delta_eta = eta_theta_levels(crit_level) - eta_theta_levels(crit_level-1)
    mag1 = mag_value1
    mag3 = mag_factor3
    DO k = crit_level+1, model_levels-1
      delta_eta = mag1 * delta_eta
      eta_theta_levels(k) = eta_theta_levels(k-1) + delta_eta
      mag1 = mag1 + mag_factor1
      IF (mag1 > mag_value3) THEN
        mag1 = mag1 - mag_factor1 + mag3
        mag3 = 2.0 * mag3
        IF (mag1 > 2.0) THEN
          mag1 = 1.3
        END IF
      ELSE IF (mag1 > mag_value2) THEN
        mag1 = mag1 - mag_factor1 + mag_factor2
      END IF
    END DO
  END IF

  ! Re-normalise eta_theta_levels, starting with top of quadratic grid
  eta_theta_levels(model_levels) = eta_theta_levels(model_levels-1)
  DO k = model_levels-1, 2, -1
    eta_theta_levels(k) = eta_theta_levels(k-1) / eta_theta_levels(model_levels)
  END DO
  eta_theta_levels(model_levels) = 1.0

  ! Reset height domain to give r_theta = first_theta_height at level 1
  height_domain = first_theta_height / eta_theta_levels(2)

  ! Add thin bottom layer
  eta_theta_levels(1) = thin_theta_height / height_domain

  ! Set eta_rho_levels
  DO k = 1, model_levels
    eta_rho_levels(k) = 0.5 * (eta_theta_levels(k-1) + eta_theta_levels(k))
  END DO

  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
    CALL umPrint('Quadratic stretched vertical grid, thin near surface',   &
                 src=RoutineName)
    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE(umMessage,'(A,E16.8)') ' height_domain (metres): ', height_domain
      CALL umPrint(umMessage, src=RoutineName)
    END IF
  END IF


CASE(vert_geometric_theta)
! Geometrically spaced grid

  delta_eta = 1.0 / model_levels
  eta_theta_levels(0) = 0.0

  DO k = 1, model_levels
    eta_theta_levels(k) = eta_theta_levels(k-1) + delta_eta
    eta_rho_levels(k)   = 0.5 * ( eta_theta_levels(k-1) + eta_theta_levels(k) )
    delta_eta = delta_eta * vert_grid_ratio
  END DO
  ! Renormalise. The order of the two calculations matters!
  DO k = 1, model_levels
    eta_rho_levels(k) = eta_rho_levels(k) / eta_theta_levels(model_levels)
    eta_theta_levels(k) = eta_theta_levels(k) / eta_theta_levels(model_levels)
  END DO

  IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
    CALL umPrint('Semi-geometric vertical grid', src=RoutineName)
    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE(umMessage,'(A,E16.8)') ' height_domain (metres): ', height_domain
      CALL umPrint(umMessage, src=RoutineName)
      IF (PrintStatus >= PrStatus_Diag) THEN
        WRITE(umMessage,'(A5,A17,A17)') 'Level', '  eta_theta_level',      &
                                        '    eta_rho_level'
        DO k = 1, model_levels
          WRITE(umMessage,'(I5,2(A1,E16.8))') k,' ', eta_theta_levels(k),  &
                                                ' ', eta_rho_levels(k)
        END DO
      END IF  ! Diag
    END IF  ! Oper
  END IF  ! Norm & PE0


CASE(vert_idl_um_grid)
! Height gen smoothed grid
! No special treatment of vertical grids.
! Do nothing.

CASE(vert_dump)
! Vertical grid taken from input dump or ancillary
! Should already have been handled in rcf_readnl_vertical.
! Do nothing.

CASE DEFAULT

  icode = 100
  WRITE(umMessage,'(A,I0)') 'Unsupported value for grid_number: ', grid_number
  CALL ereport(RoutineName, icode, umMessage)

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE  rcf_ideal_set_eta_levels
END MODULE rcf_ideal_set_eta_levels_mod
