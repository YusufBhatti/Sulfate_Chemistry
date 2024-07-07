! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE qt_bal_cld_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'QT_BAL_CLD_MOD'
CONTAINS

SUBROUTINE qt_bal_cld                                                   &
    (p_star,p_theta_levels,p,                                           &
     theta,exner_theta_levels,                                          &
     q,qcl,qcf,qcf2,                                                    &
     rhcpt, rhc_row_length, rhc_rows, bl_levels,                        &
     fv_cos_theta_latitude,                                             &
     l_mcr_qcf2, l_mixing_ratio, ntml, cumulus,                         &
     area_cloud_fraction,  bulk_cloud_fraction,                         &
     cloud_fraction_liquid,  cloud_fraction_frozen,                     &
     mype)

! Purpose:
!        reset q, t and the cloud fields to be consistent at the
!        end of the timestep
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Cloud
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

USE planet_constants_mod, ONLY: cp
USE atm_fields_bounds_mod, ONLY: pdims, pdims_s,          &
                                 tdims, tdims_s, tdims_l
USE water_constants_mod, ONLY: lc
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE ls_arcld_mod, ONLY: ls_arcld
IMPLICIT NONE

INTEGER ::                                                            &
  mype,               &  ! My processor number
  rhc_row_length,     &  ! Array size for RHcrit array
  rhc_rows,           &  ! Array size for RHcrit array
  bl_levels


LOGICAL ::                                                            &
  l_mcr_qcf2,                                                         &
                         ! true if second cloud ice variable in use
  l_mixing_ratio         ! true if using mixing ratio formulation

REAL ::                                                               &
  p(                     pdims_s%i_start:pdims_s%i_end,               &
                         pdims_s%j_start:pdims_s%j_end,               &
                         pdims_s%k_start:pdims_s%k_end+1),            &
  p_theta_levels(        pdims_s%i_start:pdims_s%i_end,               &
                         pdims_s%j_start:pdims_s%j_end,               &
                         pdims_s%k_start:pdims_s%k_end),              &
  p_star(                  pdims%i_start:pdims%i_end,                 &
                           pdims%j_start:pdims%j_end),                &
  theta(                 tdims_s%i_start:tdims_s%i_end,               &
                         tdims_s%j_start:tdims_s%j_end,               &
                         tdims_s%k_start:tdims_s%k_end),              &
  exner_theta_levels(    tdims_s%i_start:tdims_s%i_end,               &
                         tdims_s%j_start:tdims_s%j_end,               &
                         tdims_s%k_start:tdims_s%k_end),              &
  q(                     tdims_l%i_start:tdims_l%i_end,               &
                         tdims_l%j_start:tdims_l%j_end,               &
                         tdims_l%k_start:tdims_l%k_end),              &
  qcl(                   tdims_l%i_start:tdims_l%i_end,               &
                         tdims_l%j_start:tdims_l%j_end,               &
                         tdims_l%k_start:tdims_l%k_end),              &
  qcf(                   tdims_l%i_start:tdims_l%i_end,               &
                         tdims_l%j_start:tdims_l%j_end,               &
                         tdims_l%k_start:tdims_l%k_end),              &
  qcf2(                  tdims_l%i_start:tdims_l%i_end,               &
                         tdims_l%j_start:tdims_l%j_end,               &
                         tdims_l%k_start:tdims_l%k_end),              &
  rhcpt(rhc_row_length, rhc_rows, tdims%k_end),                       &
! trig arrays
    fv_cos_theta_latitude (tdims_s%i_start:tdims_s%i_end,               &
                           tdims_s%j_start:tdims_s%j_end)

! Diagnostic variables

REAL ::                                                               &
  area_cloud_fraction(     tdims%i_start:tdims%i_end,                 &
                           tdims%j_start:tdims%j_end,                 &
                                       1:tdims%k_end),                &
  bulk_cloud_fraction(   tdims_l%i_start:tdims_l%i_end,               &
                         tdims_l%j_start:tdims_l%j_end,               &
                         tdims_l%k_start:tdims_l%k_end),              &
  cloud_fraction_liquid( tdims_l%i_start:tdims_l%i_end,               &
                         tdims_l%j_start:tdims_l%j_end,               &
                         tdims_l%k_start:tdims_l%k_end),              &
  cloud_fraction_frozen( tdims_l%i_start:tdims_l%i_end,               &
                         tdims_l%j_start:tdims_l%j_end,               &
                         tdims_l%k_start:tdims_l%k_end)

INTEGER ::                                                            &
  ntml (                   pdims%i_start:pdims%i_end,                 &
                           pdims%j_start:pdims%j_end)

LOGICAL ::                                                            &
  cumulus (                tdims%i_start:tdims%i_end,                 &
                           tdims%j_start:tdims%j_end)
                         ! bl convection flag

! Local variables
REAL ::                                                               &
 p_layer_centres(          pdims%i_start:pdims%i_end,                 &
                           pdims%j_start:pdims%j_end,                 &
                                       0:pdims%k_end),                &
            ! pressure at layer centres. Same as p_theta_levels
            ! except bottom level = p_star, and at top = 0.
 p_layer_boundaries(       pdims%i_start:pdims%i_end,                 &
                           pdims%j_start:pdims%j_end,                 &
                                       0:pdims%k_end),                &
 tl(                       tdims%i_start:tdims%i_end,                 &
                           tdims%j_start:tdims%j_end,                 &
                                       1:tdims%k_end),                &
 qt(                       tdims%i_start:tdims%i_end,                 &
                           tdims%j_start:tdims%j_end,                 &
                                       1:tdims%k_end),                &
 qcf_in(                   tdims%i_start:tdims%i_end,                 &
                           tdims%j_start:tdims%j_end,                 &
                                       1:tdims%k_end),                &
 qcl_out(                  tdims%i_start:tdims%i_end,                 &
                           tdims%j_start:tdims%j_end,                 &
                                       1:tdims%k_end),                &
 cf_inout(                 tdims%i_start:tdims%i_end,                 &
                           tdims%j_start:tdims%j_end,                 &
                                       1:tdims%k_end),                &
 cfl_inout(                tdims%i_start:tdims%i_end,                 &
                           tdims%j_start:tdims%j_end,                 &
                                       1:tdims%k_end),                &
 cff_inout(                tdims%i_start:tdims%i_end,                 &
                           tdims%j_start:tdims%j_end,                 &
                                       1:tdims%k_end)

INTEGER ::                                                            &
  large_levels,                                                       &
  levels_per_level

INTEGER ::                                                            &
  i,j,k,errorstatus     ! loop variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='QT_BAL_CLD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, j, k)                          &
!$OMP SHARED(pdims, p_layer_centres,  p_layer_boundaries, p_star,      &
!$OMP    p_theta_levels, p, tl, theta, exner_theta_levels, qcl, qt,    &
!$OMP    qcf_in, q, qcf, cf_inout, bulk_cloud_fraction, cfl_inout,     &
!$OMP    cloud_fraction_liquid, cff_inout, cloud_fraction_frozen,      &
!$OMP    tdims, cp, qcf2, l_mcr_qcf2)

!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    p_layer_centres(i,j,0) = p_star(i,j)
    p_layer_boundaries(i,j,0) = p_star(i,j)
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k = pdims%k_start, pdims%k_end-1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      p_layer_centres(i,j,k) = p_theta_levels(i,j,k)
      p_layer_boundaries(i,j,k) = p(i,j,k+1)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

k=pdims%k_end
!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    p_layer_centres(i,j,k) = p_theta_levels(i,j,k)
    p_layer_boundaries(i,j,k) = 0.0
  END DO
END DO
!$OMP END DO NOWAIT

! ----------------------------------------------------------------------
! Section  Convert qT and Tl for input to cloud scheme.
! ----------------------------------------------------------------------

! Create Tl and qT
!$OMP DO SCHEDULE(STATIC)
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      tl(i,j,k) =                                                     &
             (theta(i,j,k) *exner_theta_levels(i,j,k))                &
                     - (lc * qcl(i,j,k)) / cp
      qt(i,j,k) = q(i,j,k) + qcl(i,j,k)
      qcf_in(i,j,k)=qcf(i,j,k)
      cf_inout(i,j,k)= bulk_cloud_fraction(i,j,k)
      cfl_inout(i,j,k)= cloud_fraction_liquid(i,j,k)
      cff_inout(i,j,k)= cloud_fraction_frozen(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

    ! If second cloud ice variable in use then add to qcf
    ! for the cloud scheme call
IF (l_mcr_qcf2) THEN
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        qcf_in(i,j,k) = qcf_in(i,j,k) + qcf2(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

!$OMP END PARALLEL

! ----------------------------------------------------------------------
! Section BL.4b Call cloud scheme to convert Tl and qT to T, q and qcl
!              calculate bulk_cloud fields from qT and qcf
!               and calculate area_cloud fields.
! ----------------------------------------------------------------------


! Determine number of sublevels for vertical gradient area cloud
! Want an odd number of sublevels per level: 3 is hardwired in do loops
levels_per_level = 3
large_levels = ((tdims%k_end - 2)*levels_per_level) + 2

CALL ls_arcld( p_layer_centres, rhcpt, p_layer_boundaries,            &
                 rhc_row_length, rhc_rows, bl_levels,                 &
                 levels_per_level, large_levels,                      &
                 fv_cos_theta_latitude,                               &
                 ntml, cumulus, l_mixing_ratio, qcf_in,               &
                 tl, qt, qcl_out,                                     &
                 area_cloud_fraction,  cf_inout,                      &
                 cfl_inout, cff_inout ,                               &
                 errorstatus)

! qt holds q (no halos), tl holds t(no halos),
! qcl_out holds qcl(no halos)

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i, j, k)     &
!$OMP SHARED(tdims, theta, tl, exner_theta_levels, q, qt, qcl,        &
!$OMP    qcl_out, bulk_cloud_fraction, cloud_fraction_liquid,         &
!$OMP    cf_inout, cfl_inout, cloud_fraction_frozen, cff_inout)
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      theta(i,j,k) = tl(i,j,k)/exner_theta_levels(i,j,k)
      q(i,j,k) = qt(i,j,k)
      qcl(i,j,k) =qcl_out(i,j,k)
      bulk_cloud_fraction(i,j,k) = cf_inout(i,j,k)
      cloud_fraction_liquid(i,j,k) = cfl_inout(i,j,k)
      cloud_fraction_frozen(i,j,k) = cff_inout(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE qt_bal_cld
END MODULE qt_bal_cld_mod
