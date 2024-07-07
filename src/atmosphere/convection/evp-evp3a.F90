! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the evaporation of precipitation
!
MODULE evp_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'EVP_MOD'
CONTAINS

SUBROUTINE evp(npnts, iphase, bevap, area_fac                              &
               , precip, tevp, cca, rho, delq, delpkm1, pkm1               &
               , evap, full_evap)

USE planet_constants_mod, ONLY: g

USE cv_param_mod, ONLY:                                                    &
   p_lq1, p_lq2, p_ice1, p_ice2, rho_lqp1, rho_lqp2, rho_lqa, rho_lqb,     &
   rho_icp1, rho_icp2, rho_icea, rho_iceb,                                 &
   lq_a, lq_b, lq_c, ice_a, ice_b, ice_c, ice_d

USE science_fixes_mod, ONLY: l_fix_conv_precip_evap

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
!
! Description: Calculates the evaporation of precipitation
!
! Method: UM documentation paper 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.2.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------
INTEGER, INTENT(IN) ::  &
  npnts                 & ! Vector length
 ,iphase                  ! Indication for rain (1), or snow (2)

LOGICAL, INTENT(IN) ::  &
  bevap(npnts)            ! Mask for points where evaporation takes place

REAL, INTENT(IN) ::     &
  area_fac                ! Fraction of convective cloud amount to give
                          ! local cloud area
REAL, INTENT(IN) ::     &
  precip(npnts)         & ! Amount of precipitation (kg/m**2/s)

 ,tevp(npnts)           & ! Temperature of layer k (K)

 ,cca(npnts)            & ! Convective cloud amount (fraction)

 ,rho(npnts)            & ! Density of air

 ,delq(npnts)           & ! Change in humidity mixing ratio across layer k
                          ! (kg/kg)
 ,delpkm1(npnts)        & ! Change in pressure across layer k-1 (Pa)

 ,pkm1(npnts)             ! Pressure at level k-1 (Pa)

REAL, INTENT(OUT) ::    &
  evap(npnts)             ! Evaporation

LOGICAL, INTENT(OUT) :: &
  full_evap(npnts)        ! True if all precip evaporated

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER ::       &
  i               ! Loop counter


REAL ::        &
  econ         & ! Quadratic term
 ,c1           & ! Constant
 ,c2           & ! Constant
 ,lrate        & ! Local rate of precipitation
 ,area         & ! Fractional area occupied by precipitation / downdraught
 , tl1,ti1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EVP'

!-----------------------------------------------------------------------
! Start of routine
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

full_evap(:) = .FALSE.

tl1=0.5*p_lq1
ti1=0.5*p_ice1

IF (iphase == 1) THEN        ! Rain

  DO i=1,npnts
    IF (bevap(i)) THEN
      IF (precip(i)  >   0.0) THEN

        ! The following code implements the formulae in section 2.8.4.1b
        ! of UM Documentation Paper 027: Convection Schemes:
        econ = ((lq_a*tevp(i)+lq_b)*tevp(i)+lq_c)* (100000.0/pkm1(i))
        area = area_fac*cca(i)
        lrate = precip(i)/area

        IF ( l_fix_conv_precip_evap ) THEN

          ! Correct formula doesn't include scaling by the area
          ! fraction here, as this is done afterwards in the routines that
          ! call evp (pevp_bcb and devap)
          c1 = rho_lqa * (lrate * lrate * rho(i))**tl1
          c2 = rho_lqb * (lrate**p_lq2) * (rho(i)**rho_lqp2)
          evap(i) = econ * (c1 + c2) * delq(i) * delpkm1(i) / g

          ! Limit evaporation to not exceed the total precipitation available
          IF ( evap(i) >= lrate ) THEN
            evap(i) = lrate
            full_evap(i) = .TRUE.
          END IF

        ELSE

          ! Original formula with erroneous duplicated factor of area
          c1 = rho_lqa * area * (lrate * lrate * rho(i))**tl1
          c2 = rho_lqb * area * (lrate**p_lq2) * (rho(i)**rho_lqp2)
          evap(i) = MIN( econ*(c1+c2)*delq(i)*delpkm1(i)/g, lrate )
          ! Possibly numerically dodgy way of assessing whether full evap
          IF (evap(i) == lrate) full_evap(i) = .TRUE.

        END IF

      ELSE
        evap(i) = 0.0
      END IF
    END IF
  END DO

ELSE IF (iphase == 2) THEN        ! Snow

  DO i=1,npnts
    IF (bevap(i)) THEN
      IF (precip(i)  >   0.0) THEN

        ! The following code implements the formulae in section 2.8.4.1b
        ! of UM Documentation Paper 027: Convection Schemes:
        IF (tevp(i) <= 243.58) THEN
          econ = ice_d*(100000.0/pkm1(i))
        ELSE
          econ = ((ice_a*tevp(i)+ice_b)*tevp(i)+ice_c)*(100000.0/pkm1(i))
        END IF
        area = area_fac*cca(i)
        lrate = precip(i)/area

        IF ( l_fix_conv_precip_evap ) THEN

          ! Correct formula doesn't include scaling by the area
          ! fraction here, as this is done afterwards in the routines that
          ! call evp (pevp_bcb and devap)
          c1 = rho_icea * (lrate * lrate * rho(i))**ti1
          c2 = rho_iceb * (lrate**p_ice2) * (rho(i)**rho_icp2)
          evap(i) = econ * (c1 + c2) * delq(i) * delpkm1(i) / g

          ! Limit evaporation to not exceed the total precipitation available
          IF ( evap(i) >= lrate ) THEN
            evap(i) = lrate
            full_evap(i) = .TRUE.
          ! Also don't allow negative evaporation (vapour subliming onto falling
          ! convective snow when supersaturated w.r.t. ice).
          ELSE IF ( evap(i) < 0.0 ) THEN
            evap(i) = 0.0
          END IF

        ELSE

          ! Original formula with erroneous duplicated factor of area
          c1 = rho_icea * area * (lrate * lrate * rho(i))**ti1
          c2 = rho_iceb * area * (lrate**p_ice2) * (rho(i)**rho_icp2)
          evap(i)=MAX( 0.0, MIN( econ*(c1+c2)*delq(i)*delpkm1(i)/g, lrate ) )
          ! Possibly numerically dodgy way of assessing whether full evap
          IF (evap(i) == lrate) full_evap(i) = .TRUE.

        END IF

      ELSE
        evap(i) = 0.0
      END IF
    END IF
  END DO

END IF     ! test on iphase

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE evp

END MODULE evp_mod
