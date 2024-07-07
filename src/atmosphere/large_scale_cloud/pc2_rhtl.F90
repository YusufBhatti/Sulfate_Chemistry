! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  PC2 Cloud Scheme: Calculation of RH(TL) for later use by initiation

MODULE pc2_rhtl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'PC2_RHTL_MOD'
CONTAINS

SUBROUTINE pc2_rhtl(                                                    &
!      Parallel variables
  halo_i, halo_j, off_x, off_y,                                         &
!      Array dimensions
 row_length, rows,                                                      &
!      Output value of RH(TL)
 rhts,                                                                  &
!      Logical control
 l_mixing_ratio                                                         &
 )

USE planet_constants_mod,  ONLY: lcrcp
USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim
USE atm_fields_bounds_mod, ONLY: pdims, tdims, tdims_s, tdims_l, pdims_s
USE atm_fields_mod, ONLY: q, qcl, theta, exner_theta_levels, p_theta_levels
USE nlsizes_namelist_mod,  ONLY: model_levels

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_wat_new     => qsat_wat,                       &
                    qsat_wat_mix_new => qsat_wat_mix,                   &
                    l_new_qsat_lsc !Currently defaults to FALSE

IMPLICIT NONE

! Purpose:
!   Calculate total relative humidity wrt T_L
!
! Method:
!   Straight calculation using RH_T(TL) = (q+qcl)/qsatwat(TL)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: PC2 Cloud Scheme Documentation

!  Global Variables:----------------------------------------------------
!
!  Subroutine Arguments:------------------------------------------------
INTEGER           , INTENT(IN) ::                                       &
!      Model dimensions
   row_length, rows,                                                    &
!      Parallel setup variables
    halo_i,                                                             &
!      Size of halo in i direction.
    halo_j,                                                             &
!      Size of halo in j direction.
    off_x,                                                              &
!      Size of small halo in i.
    off_y
!      Size of small halo in j.

LOGICAL           , INTENT(IN) ::                                     &
 l_mixing_ratio     ! Use mixing ratio formulation

REAL              , INTENT(OUT) ::                                    &
!      RH_T(T_L)
   rhts(row_length,rows,model_levels)

!  External functions:

!  Local arrays---------------------------------------------------------

REAL ::                                                               &
 tl(                 tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end),                      &
!      Liquid water temperature (K)
   qsl_tl(             tdims%i_start:tdims%i_end,                       &
                       tdims%j_start:tdims%j_end),                      &
!      Qsat wrt liquid water at temp TL
   p_no_halos(         pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end)
!      Pressure without halo values (Pa)

!  Local scalars

INTEGER :: k,i,j ! Loop counters:   K - vertical level index
!                                   I,J - horizontal position index

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PC2_RHTL'

! ==Main Block==--------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,tl,p_no_halos,qsl_tl)    &
!$OMP& SHARED(tdims,theta,exner_theta_levels,lcrcp,qcl,p_theta_levels, &
!$OMP& l_mixing_ratio,rhts,q,l_new_qsat_lsc)
DO k = 1, tdims%k_end

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      !       Calculate liquid temperature TL
      tl(i,j) = theta(i,j,k) * exner_theta_levels(i,j,k)              &
              -lcrcp*qcl(i,j,k)
      p_no_halos(i,j) = p_theta_levels(i,j,k)
    END DO
  END DO

  ! Calculate qsat(TL) with respect to liquid water
  IF ( l_new_qsat_lsc ) THEN
    IF ( l_mixing_ratio ) THEN
      CALL qsat_wat_mix_new(qsl_tl,tl,p_no_halos,tdims%i_len,tdims%j_len)
    ELSE
      CALL qsat_wat_new(qsl_tl,tl,p_no_halos,tdims%i_len,tdims%j_len)
    END IF
  ELSE
    ! DEPENDS ON: qsat_wat_mix
    CALL qsat_wat_mix(qsl_tl,tl,p_no_halos, tdims%i_len*tdims%j_len,          &
                      l_mixing_ratio)
  END IF


  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      !       Calculate RH_T(TL)
      rhts(i,j,k) = ( q(i,j,k) + qcl(i,j,k) ) / qsl_tl(i,j)
    END DO
  END DO

END DO  ! k
!$OMP END PARALLEL DO
! End of the subroutine

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pc2_rhtl
END MODULE pc2_rhtl_mod
