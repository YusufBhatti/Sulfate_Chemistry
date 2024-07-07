! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Boundary Layer

! PURPOSE: To calculate buoyancy parameters on p,T,q-levels

! CODE DESCRIPTION:
!   LANGUAGE: FORTRAN 95
!   THIS CODE IS WRITTEN TO UMDP 3 PROGRAMMING STANDARDS.

MODULE buoy_tq_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'BUOY_TQ_MOD'
CONTAINS

SUBROUTINE buoy_tq (                                                    &
! IN dimensions/logicals
 bl_levels,                                                             &
! IN fields
 p,t,q,qcf,qcl,cf_bulk,                                                 &
! OUT fields
 bt,bq,bt_cld,bq_cld,bt_gb,bq_gb,a_qs,a_dqsdt,dqsdt                     &
 )

USE atm_fields_bounds_mod, ONLY: tdims, tdims_l
USE gen_phys_inputs_mod, ONLY: l_mr_physics
USE planet_constants_mod, ONLY: r, repsilon, c_virtual, etar, lcrcp, ls, lsrcp
USE water_constants_mod, ONLY: lc, lf, tm
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_new         => qsat,                           &
                    qsat_mix_new     => qsat_mix,                       &
                    l_new_qsat_bl !Currently defaults to FALSE

IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.
INTEGER, INTENT(IN) ::                                                  &
 bl_levels              ! IN No. of atmospheric levels for which

REAL, INTENT(IN) ::                                                     &
 p(tdims%i_start:tdims%i_end,                                           &
   tdims%j_start:tdims%j_end,0:bl_levels+1),                            &
                                    ! IN Pressure at pressure points.
 t(tdims%i_start:tdims%i_end,                                           &
   tdims%j_start:tdims%j_end,bl_levels),                                &
                                    ! IN Temperature (K). At P points
 q(tdims_l%i_start:tdims_l%i_end,                                       &
   tdims_l%j_start:tdims_l%j_end,tdims_l%k_start:bl_levels),            &
                              ! IN Sp humidity (kg water per kg air).
 qcl(tdims_l%i_start:tdims_l%i_end,                                     &
     tdims_l%j_start:tdims_l%j_end,tdims_l%k_start:bl_levels),          &
                              ! IN Cloud liq water (kg per kg air).
 qcf(tdims_l%i_start:tdims_l%i_end,                                     &
     tdims_l%j_start:tdims_l%j_end,tdims_l%k_start:bl_levels),          &
                              ! IN Cloud liq water (kg per kg air).
 cf_bulk(tdims%i_start:tdims%i_end,                                     &
         tdims%j_start:tdims%j_end, bl_levels)! IN Cloud fraction (decimal).

! ARGUMENTS WITH INTENT OUT. IE: OUTPUT VARIABLES.
REAL, INTENT(OUT) ::                                                    &
 bq(tdims%i_start:tdims%i_end,                                          &
    tdims%j_start:tdims%j_end,bl_levels),                               &
                              ! OUT A buoyancy parameter for clear air
 bt(tdims%i_start:tdims%i_end,                                          &
    tdims%j_start:tdims%j_end,bl_levels),                               &
                              ! OUT A buoyancy parameter for clear air
 bq_cld(tdims%i_start:tdims%i_end,                                      &
        tdims%j_start:tdims%j_end,bl_levels),                           &
                            ! OUT A buoyancy parameter for cloudy air
 bt_cld(tdims%i_start:tdims%i_end,                                      &
        tdims%j_start:tdims%j_end,bl_levels),                           &
                            ! OUT A buoyancy parameter for cloudy air
 bq_gb(tdims%i_start:tdims%i_end,                                       &
       tdims%j_start:tdims%j_end,bl_levels),                            &
                              ! OUT A grid-box mean buoyancy parameter
 bt_gb(tdims%i_start:tdims%i_end,                                       &
       tdims%j_start:tdims%j_end,bl_levels),                            &
                              ! OUT A grid-box mean buoyancy parameter
 a_qs(tdims%i_start:tdims%i_end,                                        &
       tdims%j_start:tdims%j_end,bl_levels),                            &
                            ! OUT Saturated lapse rate factor
 a_dqsdt(tdims%i_start:tdims%i_end,                                     &
         tdims%j_start:tdims%j_end,bl_levels),                          &
                            ! OUT Saturated lapse rate factor
 dqsdt(tdims%i_start:tdims%i_end,                                       &
       tdims%j_start:tdims%j_end,bl_levels)
                            ! OUT Derivative of q_SAT w.r.t. T

! LOCAL VARIABLES.
REAL ::                                                                 &
 qs(tdims%i_start:tdims%i_end,                                          &
    tdims%j_start:tdims%j_end), & ! WORK Saturated mixing ratio.

 tmp1(tdims%i_start:tdims%i_end), & ! TEMP array to contain lc or ls

 tmp2(tdims%i_start:tdims%i_end) ! TEMP array to contain lcrcp or lsrcp

INTEGER ::                                                              &
  i,j,                                                                  &
  k

REAL ::                                                                 &
  bc

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BUOY_TQ'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-----------------------------------------------------------------------
! 1.  Loop round levels.
!-----------------------------------------------------------------------

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP PRIVATE(i, j, k, bc, qs, tmp1, tmp2)                                    &
!$OMP SHARED(bl_levels, p, t, q, qcf, qcl, cf_bulk, bt, bq, bt_cld, bq_cld,   &
!$OMP        bt_gb, bq_gb, a_qs, a_dqsdt, dqsdt, tdims, l_mr_physics, r,     &
!$OMP        repsilon, c_virtual, etar, lcrcp, ls, lsrcp, l_new_qsat_bl)

DO k = 1, bl_levels

  !-----------------------------------------------------------------------
  ! 1.1 Calculate saturated specific humidity at pressure and
  !     temperature of current level.
  !-----------------------------------------------------------------------

!No halos in these variables
  IF ( l_new_qsat_bl ) THEN
    IF ( l_mr_physics ) THEN
      CALL qsat_mix_new(qs,t(:,:,k),p(:,:,k),tdims%i_end,tdims%j_end)
    ELSE
      CALL qsat_new(qs,t(:,:,k),p(:,:,k),tdims%i_end,tdims%j_end)
    END IF
  ELSE
! DEPENDS ON: qsat_mix
    CALL qsat_mix(qs,t(1,1,k),p(1,1,k),tdims%i_end*tdims%j_end,l_mr_physics)
  END IF

  !   Using the temp arrays and splitting the i index helps vectorisation,
  !   Before the index split, no vectorisation was taking place as each
  !   conditional was computationally heavy with only minor differences
  !   between each side (differences being the lc, ls, lcrcp & lsrcp variables)

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      IF (t(i,j,k)  >   tm) THEN
        tmp1(i) = lc
        tmp2(i) = lcrcp
        !            ...  (Clausius-Clapeyron) for T above freezing
      ELSE
        tmp1(i) = ls
        tmp2(i) = lsrcp
        !            ...  (Clausius-Clapeyron) for T below freezing
      END IF
    END DO ! p_points,i

      !---------------------------------------------------------------
      ! 1.2 Calculate buoyancy parameters BT and BQ, required for the
      !     calculation of stability.
      !---------------------------------------------------------------

    DO i = tdims%i_start, tdims%i_end
      bt(i,j,k) = 1.0/t(i,j,k)
      bq(i,j,k) =                                                       &
        c_virtual/(1.0+c_virtual*q(i,j,k)-qcl(i,j,k)-qcf(i,j,k))


      dqsdt(i,j,k) = (repsilon * tmp1(i) * qs(i,j))                     &
                   / ( r * t(i,j,k) * t(i,j,k) )


      a_qs(i,j,k) = 1.0 / (1.0 + tmp2(i)*dqsdt(i,j,k))

      a_dqsdt(i,j,k) = a_qs(i,j,k) * dqsdt(i,j,k)

      bc = tmp2(i)*bt(i,j,k) - etar*bq(i,j,k)

      !--------------------------------------------------------------
      ! 1.3 Calculate in-cloud buoyancy parameters.
      !--------------------------------------------------------------

      bt_cld(i,j,k) = bt(i,j,k) - a_dqsdt(i,j,k) * bc
      bq_cld(i,j,k) = bq(i,j,k) + a_qs(i,j,k) * bc

      !--------------------------------------------------------------
      ! 1.4 Calculate grid-box mean buoyancy parameters.
      !--------------------------------------------------------------

      bt_gb(i,j,k) = bt(i,j,k) +                                        &
                     cf_bulk(i,j,k)*( bt_cld(i,j,k) - bt(i,j,k) )
      bq_gb(i,j,k) = bq(i,j,k) +                                        &
                     cf_bulk(i,j,k)*( bq_cld(i,j,k) - bq(i,j,k) )

    END DO ! p_points,i
  END DO ! p_points,j
END DO ! bl_levels

!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE buoy_tq
END MODULE buoy_tq_mod
