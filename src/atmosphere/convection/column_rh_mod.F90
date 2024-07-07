! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! To calculate the column integrated relative humidity
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
MODULE column_rh_mod


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

! Model level heights from centre of Earth
USE level_heights_mod, ONLY: &
  r_theta_levels             &  ! Radii on theta levels (m)
 ,r_rho_levels                  ! Radii on rho levels (m)

! Number of model levels
USE nlsizes_namelist_mod, ONLY: model_levels

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_new         => qsat,                           &
                    qsat_mix_new     => qsat_mix,                       &
                    l_new_qsat_conv !Currently defaults to FALSE

USE gen_phys_inputs_mod, ONLY: l_mr_physics

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='COLUMN_RH_MOD'

CONTAINS

SUBROUTINE calc_column_rh(                 &
   npnts, nlcl                             &
   , row_length,rows                       &
   , index_i, index_j                      &
   , p, t, q, rho_only                     &
   , z_half, zmax                          &
   , qsat, column_rh, column_rh_bl         &
   , column_q )

IMPLICIT NONE

INTEGER, INTENT(IN) ::     &
  npnts                    & ! Number of points
, nlcl(npnts)                ! Level of LCL

INTEGER, INTENT(IN) ::     &
  row_length               & ! Local number of points on a row
, rows                       ! Local number of rows in a theta field

INTEGER, INTENT(IN) ::     &
  index_i(npnts)           & ! Column number of unstable points
, index_j(npnts)             ! Row number of unstable points

REAL, INTENT(IN) ::          &
   p(npnts, model_levels)    & ! pressure (Pa)
   , t(npnts, model_levels)  & ! temperature (K)
   , q(npnts, model_levels)    ! water vapour (kg/kg)

REAL, INTENT(IN) ::        &
   rho_only(row_length,rows,model_levels)
                             ! Density (kg/m3)

REAL, INTENT(IN) ::        &
   z_half(npnts, model_levels) ! Height on half levels

REAL, INTENT(IN) ::        &
   zmax                      ! Top of column

REAL, INTENT(OUT) ::       &
   qsat(npnts, model_levels) ! q saturation value

REAL, INTENT(OUT) ::       &
   column_rh(npnts)        & ! Column integrated RH (fraction)
   , column_rh_bl(npnts)     ! Column integrated RH up to LCL

REAL, INTENT(OUT) ::       &
   column_q(npnts)           ! Column integrated q

!------------------------------
! Local variables
!------------------------------

REAL :: column_qsat          ! column_qsat
REAL :: temp_mass            ! column_qsat

INTEGER :: k,ii,i,j          ! loop counters

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_COLUMN_RH'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO  k = 1,model_levels
  IF ( l_new_qsat_conv ) THEN
    IF ( l_mr_physics ) THEN
      CALL qsat_mix_new(qsat(:,k),t(:,k),p(:,k),npnts)
    ELSE
      CALL qsat_new(qsat(:,k),t(:,k),p(:,k),npnts)
    END IF
  ELSE
    ! DEPENDS ON: qsat_mix
    CALL qsat_mix(qsat(:,k),t(:,k),p(:,k),npnts,l_mr_physics)
  END IF
END DO

DO ii=1,npnts
  i = index_i(ii)
  j = index_j(ii)

  column_q(ii) = 0.0
  column_qsat  = 0.0

  DO k = 2,model_levels-1

    IF (z_half(ii,k) <= zmax) THEN

      temp_mass = rho_only(i,j,k)                                         &
                *  r_theta_levels(i,j,k) * r_theta_levels(i,j,k)          &
                * (r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k))

      column_q(ii) = column_q(ii) + temp_mass * q(ii,k)
      column_qsat  = column_qsat  + temp_mass * qsat(ii,k)

      IF ( k == nlcl(ii) ) THEN
        column_rh_bl(ii) = column_q(ii)/column_qsat
      END IF
    END IF

  END DO

  column_rh(ii)=column_q(ii)/column_qsat

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN


END SUBROUTINE calc_column_rh

END MODULE column_rh_mod
