! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Service

! Overarching module for calculating saturation vapour pressure.
! Calculates:
! Saturation Specific Humidity (Qsat): Vapour to Liquid/Ice.
! Saturation Specific Humidity (Qsat_Wat): Vapour to Liquid.

! All return a saturation mixing ratio given a temperature and pressure
! using saturation vapour pressures calculated using the Goff-Gratch
! formulae, adopted by the WMO as taken from Landolt-Bornstein, 1987
! Numerical Data and Functional Relationships in Science and Technolgy.
! Group V/vol 4B meteorology. Phyiscal and Chemical properties or air, P35

! Note regarding values stored in the lookup tables (es):
! _wat versions     : over water above and below 0 deg c.
! non-_wat versions : over water above 0 deg c
!                     over ice below 0 deg c

! Method:
! Uses lookup tables to find eSAT, calculates qSAT directly from that.

! Documentation: UMDP No.29

MODULE qsat_mod

! Module-wide USE statements
USE um_types, ONLY: real64, real32

IMPLICIT NONE

! Set everything as private, and only expose the interfaces we want to
PRIVATE
PUBLIC :: qsat, qsat_wat, qsat_mix, qsat_wat_mix, l_new_qsat_bl,              &
          l_new_qsat_chem_aero, l_new_qsat_atm_serv, l_new_qsat_mphys,        &
          l_new_qsat_cntl, l_new_qsat_misc, l_new_qsat_lsc,                   &
          l_new_qsat_acassim, l_new_qsat_lsprec, l_new_qsat_conv,             &
          l_new_qsat_ideal, l_new_qsat_jules, init_qsat_switches

! This will migrate to science_fixes_mod. These are changes in 
! init_qsat_switches below, according to setting in science_fixes_mod.
! If you want  to experiment with fine control, set what you want where
! indicated in init_qsat_switches
LOGICAL :: l_new_qsat_bl         = .FALSE. !Boundary Layer
LOGICAL :: l_new_qsat_chem_aero  = .FALSE. !Chemistry, aerosols
LOGICAL :: l_new_qsat_atm_serv   = .FALSE. !Atmosphere service
LOGICAL :: l_new_qsat_mphys      = .FALSE. !Micro physics, not lsp
LOGICAL :: l_new_qsat_lsprec     = .FALSE. !Large scale precip
LOGICAL :: l_new_qsat_cntl       = .FALSE. !Control level
LOGICAL :: l_new_qsat_misc       = .FALSE. !Misc others
LOGICAL :: l_new_qsat_lsc        = .FALSE. !Large scale cloud
LOGICAL :: l_new_qsat_acassim    = .FALSE. !AC assimilation
LOGICAL :: l_new_qsat_conv       = .FALSE. !Convection
LOGICAL :: l_new_qsat_ideal      = .FALSE. !Idealised
LOGICAL :: l_new_qsat_jules      = .FALSE. !JULES

INTERFACE qsat
  MODULE PROCEDURE qsat_real64_1D, qsat_real32_1D,                            &
                   qsat_real64_2D, qsat_real32_2D,                            &
                   qsat_real64_3D, qsat_real32_3D,                            &
                   qsat_real64_scalar, qsat_real32_scalar
END INTERFACE

INTERFACE qsat_wat
  MODULE PROCEDURE qsat_wat_real64_1D, qsat_wat_real32_1D,                    &
                   qsat_wat_real64_2D, qsat_wat_real32_2D,                    &
                   qsat_wat_real64_3D, qsat_wat_real32_3D,                    &
                   qsat_wat_real64_scalar, qsat_wat_real32_scalar
END INTERFACE

INTERFACE qsat_mix
  MODULE PROCEDURE qsat_mix_real64_1D, qsat_mix_real32_1D,                    &
                   qsat_mix_real64_2D, qsat_mix_real32_2D,                    &
                   qsat_mix_real64_3D, qsat_mix_real32_3D,                    &
                   qsat_mix_real64_scalar, qsat_mix_real32_scalar
END INTERFACE

INTERFACE qsat_wat_mix
  MODULE PROCEDURE qsat_wat_mix_real64_1D, qsat_wat_mix_real32_1D,            &
                   qsat_wat_mix_real64_2D, qsat_wat_mix_real32_2D,            &
                   qsat_wat_mix_real64_3D, qsat_wat_mix_real32_3D,            &
                   qsat_wat_mix_real64_scalar, qsat_wat_mix_real32_scalar
END INTERFACE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='QSAT_NEW_MOD'

CONTAINS

! Set up the qsat switches according to what was read in by gen_phys_inputs_mod

SUBROUTINE init_qsat_switches()

USE gen_phys_inputs_mod, ONLY: l_new_qsat

IMPLICIT NONE

!Physics
l_new_qsat_bl         = l_new_qsat
l_new_qsat_mphys      = l_new_qsat
l_new_qsat_lsprec     = l_new_qsat
l_new_qsat_lsc        = l_new_qsat
l_new_qsat_conv       = l_new_qsat
l_new_qsat_jules      = l_new_qsat

!Chemistry
l_new_qsat_chem_aero  = l_new_qsat

!Others
l_new_qsat_atm_serv   = l_new_qsat
l_new_qsat_cntl       = l_new_qsat
l_new_qsat_misc       = l_new_qsat
l_new_qsat_acassim    = l_new_qsat
l_new_qsat_ideal      = l_new_qsat

!If you wish to manually override these for your own experiments, then here 
!is a good place to do it

!When running the 32-bit lsp scheme, we need to use the new precision-aware
!qsat family. Not doing so causes the model to fail at bcgstab.
!Obviously, this will disappear when the old qsat is retired.
#if defined(LSPREC_32B)
l_new_qsat_lsprec     = .TRUE.
#endif

RETURN
END SUBROUTINE init_qsat_switches

! Do all the 64-bit routines first, then all the 32-bit ones

! ------------------------------------------------------------------------------
! 1D subroutines.
! These are called by all the other dimensionalities.
! There are 4 routines that provide every permutation of with/without
! _wat and _mix. The _mix versions apply a slightly different algorithm,
! and _wat uses different data tables.

SUBROUTINE qsat_real64_1D(qs, t, p, npnts)
! Remember to use the correct flavour of es
USE qsat_data_real64_mod, ONLY:  t_low, t_high, delta_t, es
! Use in required info. Make sure it's the right KIND
USE planet_constants_mod, ONLY: repsilon, one_minus_epsilon
USE conversions_mod,      ONLY: zerodegc
IMPLICIT NONE
INTEGER, PARAMETER :: prec = real64
#include "qsat_mod_qsat.h"
END SUBROUTINE qsat_real64_1D

SUBROUTINE qsat_wat_real64_1D(qs, t, p, npnts)
! Remember to use the correct flavour of es
USE qsat_data_real64_mod, ONLY:  t_low, t_high, delta_t, es => es_wat
! Use in required info. Make sure it's the right KIND
USE planet_constants_mod, ONLY: repsilon, one_minus_epsilon
USE conversions_mod,      ONLY: zerodegc
IMPLICIT NONE
INTEGER, PARAMETER :: prec = real64
#include "qsat_mod_qsat.h"
END SUBROUTINE qsat_wat_real64_1D

SUBROUTINE qsat_mix_real64_1D(qs, t, p, npnts)
! Remember to use the correct flavour of es
USE qsat_data_real64_mod, ONLY:  t_low, t_high, delta_t, es
! Use in required info. Make sure it's the right KIND
USE planet_constants_mod, ONLY: repsilon, one_minus_epsilon
USE conversions_mod,      ONLY: zerodegc
IMPLICIT NONE
INTEGER, PARAMETER :: prec = real64
#include "qsat_mod_qsat_mix.h"
END SUBROUTINE qsat_mix_real64_1D

SUBROUTINE qsat_wat_mix_real64_1D(qs, t, p, npnts)
! Remember to use the correct flavour of es
USE qsat_data_real64_mod, ONLY:  t_low, t_high, delta_t, es => es_wat
! Use in required info. Make sure it's the right KIND
USE planet_constants_mod, ONLY: repsilon, one_minus_epsilon
USE conversions_mod,      ONLY: zerodegc
IMPLICIT NONE
INTEGER, PARAMETER :: prec = real64
#include "qsat_mod_qsat_mix.h"
RETURN
END SUBROUTINE qsat_wat_mix_real64_1D

! ------------------------------------------------------------------------------
! 2D subroutines

SUBROUTINE qsat_real64_2D(qs, t, p, npntsi, npntsj)
IMPLICIT NONE
INTEGER,               INTENT(IN)  :: npntsi, npntsj
REAL (KIND=real64),    INTENT(IN)  :: t(npntsi,npntsj), p(npntsi,npntsj)
REAL (KIND=real64),    INTENT(OUT) :: qs(npntsi,npntsj)
INTEGER                            :: j
DO j=1, npntsj
  CALL qsat_real64_1D(qs(1,j),t(1,j),p(1,j),npntsi)
END DO
RETURN
END SUBROUTINE qsat_real64_2D

SUBROUTINE qsat_wat_real64_2D(qs, t, p, npntsi, npntsj)
IMPLICIT NONE
INTEGER,               INTENT(IN)  :: npntsi, npntsj
REAL (KIND=real64),    INTENT(IN)  :: t(npntsi,npntsj), p(npntsi,npntsj)
REAL (KIND=real64),    INTENT(OUT) :: qs(npntsi,npntsj)
INTEGER                            :: j
DO j=1, npntsj
  CALL qsat_wat_real64_1D(qs(1,j),t(1,j),p(1,j),npntsi)
END DO
RETURN
END SUBROUTINE qsat_wat_real64_2D

SUBROUTINE qsat_mix_real64_2D(qs, t, p, npntsi, npntsj)
IMPLICIT NONE
INTEGER,               INTENT(IN)  :: npntsi, npntsj
REAL (KIND=real64),    INTENT(IN)  :: t(npntsi,npntsj), p(npntsi,npntsj)
REAL (KIND=real64),    INTENT(OUT) :: qs(npntsi,npntsj)
INTEGER                            :: j
DO j=1, npntsj
  CALL qsat_mix_real64_1D(qs(1,j),t(1,j),p(1,j),npntsi)
END DO
RETURN
END SUBROUTINE qsat_mix_real64_2D

SUBROUTINE qsat_wat_mix_real64_2D(qs, t, p, npntsi, npntsj)
IMPLICIT NONE
INTEGER,               INTENT(IN)  :: npntsi, npntsj
REAL (KIND=real64),    INTENT(IN)  :: t(npntsi,npntsj), p(npntsi,npntsj)
REAL (KIND=real64),    INTENT(OUT) :: qs(npntsi,npntsj)
INTEGER                            :: j
DO j=1, npntsj
  CALL qsat_wat_mix_real64_1D(qs(1,j),t(1,j),p(1,j),npntsi)
END DO
RETURN
END SUBROUTINE qsat_wat_mix_real64_2D

! ------------------------------------------------------------------------------
! 3D subroutines

SUBROUTINE qsat_real64_3D(qs, t, p, npntsi, npntsj, npntsk)
IMPLICIT NONE
INTEGER,               INTENT(IN)  :: npntsi, npntsj, npntsk
REAL (KIND=real64),    INTENT(IN)  :: t(npntsi,npntsj,npntsk),                &
                                      p(npntsi,npntsj,npntsk)
REAL (KIND=real64),    INTENT(OUT) :: qs(npntsi,npntsj,npntsk)
INTEGER                            :: j, k
DO k=1, npntsk
  DO j=1, npntsj
    CALL qsat_real64_1D(qs(1,j,k),t(1,j,k),p(1,j,k),npntsi)
  END DO
END DO
RETURN
END SUBROUTINE qsat_real64_3D

SUBROUTINE qsat_wat_real64_3D(qs, t, p, npntsi, npntsj, npntsk)
IMPLICIT NONE
INTEGER,               INTENT(IN)  :: npntsi, npntsj, npntsk
REAL (KIND=real64),    INTENT(IN)  :: t(npntsi,npntsj,npntsk),                &
                                      p(npntsi,npntsj,npntsk)
REAL (KIND=real64),    INTENT(OUT) :: qs(npntsi,npntsj,npntsk)
INTEGER                            :: j, k
DO k=1, npntsk
  DO j=1, npntsj
    CALL qsat_wat_real64_1D(qs(1,j,k),t(1,j,k),p(1,j,k),npntsi)
  END DO
END DO
RETURN
END SUBROUTINE qsat_wat_real64_3D

SUBROUTINE qsat_mix_real64_3D(qs, t, p, npntsi, npntsj, npntsk)
IMPLICIT NONE
INTEGER,               INTENT(IN)  :: npntsi, npntsj, npntsk
REAL (KIND=real64),    INTENT(IN)  :: t(npntsi,npntsj,npntsk),                &
                                      p(npntsi,npntsj,npntsk)
REAL (KIND=real64),    INTENT(OUT) :: qs(npntsi,npntsj,npntsk)
INTEGER                            :: j, k
DO k=1, npntsk
  DO j=1, npntsj
    CALL qsat_mix_real64_1D(qs(1,j,k),t(1,j,k),p(1,j,k),npntsi)
  END DO
END DO
RETURN
END SUBROUTINE qsat_mix_real64_3D

SUBROUTINE qsat_wat_mix_real64_3D(qs, t, p, npntsi, npntsj, npntsk)
IMPLICIT NONE
INTEGER,               INTENT(IN)  :: npntsi, npntsj, npntsk
REAL (KIND=real64),    INTENT(IN)  :: t(npntsi,npntsj,npntsk),                &
                                      p(npntsi,npntsj,npntsk)
REAL (KIND=real64),    INTENT(OUT) :: qs(npntsi,npntsj,npntsk)
INTEGER                            :: j, k
DO k=1, npntsk
  DO j=1, npntsj
    CALL qsat_wat_mix_real64_1D(qs(1,j,k),t(1,j,k),p(1,j,k),npntsi)
  END DO
END DO
RETURN
END SUBROUTINE qsat_wat_mix_real64_3D

! ------------------------------------------------------------------------------
! scalar subroutines

SUBROUTINE qsat_real64_scalar(qs, t, p)
IMPLICIT NONE
REAL (KIND=real64),    INTENT(IN)  :: t, p
REAL (KIND=real64),    INTENT(OUT) :: qs
REAL (KIND=real64)                 :: t_arr(1), p_arr(1), qs_arr(1)
t_arr(1) = t
p_arr(1) = p
CALL qsat_real64_1D(qs_arr,t_arr,p_arr,1)
qs = qs_arr(1)
RETURN
END SUBROUTINE qsat_real64_scalar

SUBROUTINE qsat_wat_real64_scalar(qs, t, p)
IMPLICIT NONE
REAL (KIND=real64),    INTENT(IN)  :: t, p
REAL (KIND=real64),    INTENT(OUT) :: qs
REAL (KIND=real64)                 :: t_arr(1), p_arr(1), qs_arr(1)
t_arr(1) = t
p_arr(1) = p
CALL qsat_wat_real64_1D(qs_arr,t_arr,p_arr,1)
qs = qs_arr(1)
RETURN
END SUBROUTINE qsat_wat_real64_scalar

SUBROUTINE qsat_mix_real64_scalar(qs, t, p)
IMPLICIT NONE
REAL (KIND=real64),    INTENT(IN)  :: t, p
REAL (KIND=real64),    INTENT(OUT) :: qs
REAL (KIND=real64)                 :: t_arr(1), p_arr(1), qs_arr(1)
t_arr(1) = t
p_arr(1) = p
CALL qsat_mix_real64_1D(qs_arr,t_arr,p_arr,1)
qs = qs_arr(1)
RETURN
END SUBROUTINE qsat_mix_real64_scalar

SUBROUTINE qsat_wat_mix_real64_scalar(qs, t, p)
IMPLICIT NONE
REAL (KIND=real64),    INTENT(IN)  :: t, p
REAL (KIND=real64),    INTENT(OUT) :: qs
REAL (KIND=real64)                 :: t_arr(1), p_arr(1), qs_arr(1)
t_arr(1) = t
p_arr(1) = p
CALL qsat_wat_mix_real64_1D(qs_arr,t_arr,p_arr,1)
qs = qs_arr(1)
RETURN
END SUBROUTINE qsat_wat_mix_real64_scalar

!Now for all the 32-bit routines

! ------------------------------------------------------------------------------
! 1D subroutines.
! These are called by all the other dimensionalities.
! There are 4 routines that provide every permutation of with/without
! _wat and _mix. The _mix versions apply a slightly different algorithm,
! and _wat uses different data tables.

SUBROUTINE qsat_real32_1D(qs, t, p, npnts)
! Remember to use the correct flavour of es
USE qsat_data_real32_mod, ONLY:  t_low, t_high, delta_t, es
! Use in required info. Make sure it's the right KIND
USE planet_constants_mod, ONLY: repsilon          => repsilon_32b,            &
                                one_minus_epsilon => one_minus_epsilon_32b
USE conversions_mod,      ONLY: zerodegc   => zerodegc_32
IMPLICIT NONE
INTEGER, PARAMETER :: prec = real32
#include "qsat_mod_qsat.h"
END SUBROUTINE qsat_real32_1D

SUBROUTINE qsat_wat_real32_1D(qs, t, p, npnts)
! Remember to use the correct flavour of es
USE qsat_data_real32_mod, ONLY:  t_low, t_high, delta_t, es => es_wat
! Use in required info. Make sure it's the right KIND
USE planet_constants_mod, ONLY: repsilon          => repsilon_32b,            &
                                one_minus_epsilon => one_minus_epsilon_32b
USE conversions_mod,      ONLY: zerodegc   => zerodegc_32
IMPLICIT NONE
INTEGER, PARAMETER :: prec = real32
#include "qsat_mod_qsat.h"
END SUBROUTINE qsat_wat_real32_1D

SUBROUTINE qsat_mix_real32_1D(qs, t, p, npnts)
! Remember to use the correct flavour of es
USE qsat_data_real32_mod, ONLY:  t_low, t_high, delta_t, es
! Use in required info. Make sure it's the right KIND
USE planet_constants_mod, ONLY: repsilon          => repsilon_32b,            &
                                one_minus_epsilon => one_minus_epsilon_32b
USE conversions_mod,      ONLY: zerodegc   => zerodegc_32
IMPLICIT NONE
INTEGER, PARAMETER :: prec = real32
#include "qsat_mod_qsat_mix.h"
END SUBROUTINE qsat_mix_real32_1D

SUBROUTINE qsat_wat_mix_real32_1D(qs, t, p, npnts)
! Remember to use the correct flavour of es
USE qsat_data_real32_mod, ONLY:  t_low, t_high, delta_t, es => es_wat
! Use in required info. Make sure it's the right KIND
USE planet_constants_mod, ONLY: repsilon          => repsilon_32b,            &
                                one_minus_epsilon => one_minus_epsilon_32b
USE conversions_mod,      ONLY: zerodegc   => zerodegc_32
IMPLICIT NONE
INTEGER, PARAMETER :: prec = real32
#include "qsat_mod_qsat_mix.h"
RETURN
END SUBROUTINE qsat_wat_mix_real32_1D

! ------------------------------------------------------------------------------
! 2D subroutines

SUBROUTINE qsat_real32_2D(qs, t, p, npntsi, npntsj)
IMPLICIT NONE
INTEGER,               INTENT(IN)  :: npntsi, npntsj
REAL (KIND=real32),    INTENT(IN)  :: t(npntsi,npntsj), p(npntsi,npntsj)
REAL (KIND=real32),    INTENT(OUT) :: qs(npntsi,npntsj)
INTEGER                            :: j
DO j=1, npntsj
  CALL qsat_real32_1D(qs(1,j),t(1,j),p(1,j),npntsi)
END DO
RETURN
END SUBROUTINE qsat_real32_2D

SUBROUTINE qsat_wat_real32_2D(qs, t, p, npntsi, npntsj)
IMPLICIT NONE
INTEGER,               INTENT(IN)  :: npntsi, npntsj
REAL (KIND=real32),    INTENT(IN)  :: t(npntsi,npntsj), p(npntsi,npntsj)
REAL (KIND=real32),    INTENT(OUT) :: qs(npntsi,npntsj)
INTEGER                            :: j
DO j=1, npntsj
  CALL qsat_wat_real32_1D(qs(1,j),t(1,j),p(1,j),npntsi)
END DO
RETURN
END SUBROUTINE qsat_wat_real32_2D

SUBROUTINE qsat_mix_real32_2D(qs, t, p, npntsi, npntsj)
IMPLICIT NONE
INTEGER,               INTENT(IN)  :: npntsi, npntsj
REAL (KIND=real32),    INTENT(IN)  :: t(npntsi,npntsj), p(npntsi,npntsj)
REAL (KIND=real32),    INTENT(OUT) :: qs(npntsi,npntsj)
INTEGER                            :: j
DO j=1, npntsj
  CALL qsat_mix_real32_1D(qs(1,j),t(1,j),p(1,j),npntsi)
END DO
RETURN
END SUBROUTINE qsat_mix_real32_2D

SUBROUTINE qsat_wat_mix_real32_2D(qs, t, p, npntsi, npntsj)
IMPLICIT NONE
INTEGER,               INTENT(IN)  :: npntsi, npntsj
REAL (KIND=real32),    INTENT(IN)  :: t(npntsi,npntsj), p(npntsi,npntsj)
REAL (KIND=real32),    INTENT(OUT) :: qs(npntsi,npntsj)
INTEGER                            :: j
DO j=1, npntsj
  CALL qsat_wat_mix_real32_1D(qs(1,j),t(1,j),p(1,j),npntsi)
END DO
RETURN
END SUBROUTINE qsat_wat_mix_real32_2D

! ------------------------------------------------------------------------------
! 3D subroutines

SUBROUTINE qsat_real32_3D(qs, t, p, npntsi, npntsj, npntsk)
IMPLICIT NONE
INTEGER,               INTENT(IN)  :: npntsi, npntsj, npntsk
REAL (KIND=real32),    INTENT(IN)  :: t(npntsi,npntsj,npntsk),                &
                                      p(npntsi,npntsj,npntsk)
REAL (KIND=real32),    INTENT(OUT) :: qs(npntsi,npntsj,npntsk)
INTEGER                            :: j, k
DO k=1, npntsk
  DO j=1, npntsj
    CALL qsat_real32_1D(qs(1,j,k),t(1,j,k),p(1,j,k),npntsi)
  END DO
END DO
RETURN
END SUBROUTINE qsat_real32_3D

SUBROUTINE qsat_wat_real32_3D(qs, t, p, npntsi, npntsj, npntsk)
IMPLICIT NONE
INTEGER,               INTENT(IN)  :: npntsi, npntsj, npntsk
REAL (KIND=real32),    INTENT(IN)  :: t(npntsi,npntsj,npntsk),                &
                                      p(npntsi,npntsj,npntsk)
REAL (KIND=real32),    INTENT(OUT) :: qs(npntsi,npntsj,npntsk)
INTEGER                            :: j, k
DO k=1, npntsk
  DO j=1, npntsj
    CALL qsat_wat_real32_1D(qs(1,j,k),t(1,j,k),p(1,j,k),npntsi)
  END DO
END DO
RETURN
END SUBROUTINE qsat_wat_real32_3D

SUBROUTINE qsat_mix_real32_3D(qs, t, p, npntsi, npntsj, npntsk)
IMPLICIT NONE
INTEGER,               INTENT(IN)  :: npntsi, npntsj, npntsk
REAL (KIND=real32),    INTENT(IN)  :: t(npntsi,npntsj,npntsk),                &
                                      p(npntsi,npntsj,npntsk)
REAL (KIND=real32),    INTENT(OUT) :: qs(npntsi,npntsj,npntsk)
INTEGER                            :: j, k
DO k=1, npntsk
  DO j=1, npntsj
    CALL qsat_mix_real32_1D(qs(1,j,k),t(1,j,k),p(1,j,k),npntsi)
  END DO
END DO
RETURN
END SUBROUTINE qsat_mix_real32_3D

SUBROUTINE qsat_wat_mix_real32_3D(qs, t, p, npntsi, npntsj, npntsk)
IMPLICIT NONE
INTEGER,               INTENT(IN)  :: npntsi, npntsj, npntsk
REAL (KIND=real32),    INTENT(IN)  :: t(npntsi,npntsj,npntsk),                &
                                      p(npntsi,npntsj,npntsk)
REAL (KIND=real32),    INTENT(OUT) :: qs(npntsi,npntsj,npntsk)
INTEGER                            :: j, k
DO k=1, npntsk
  DO j=1, npntsj
    CALL qsat_wat_mix_real32_1D(qs(1,j,k),t(1,j,k),p(1,j,k),npntsi)
  END DO
END DO
RETURN
END SUBROUTINE qsat_wat_mix_real32_3D

! ------------------------------------------------------------------------------
! scalar subroutines

SUBROUTINE qsat_real32_scalar(qs, t, p)
IMPLICIT NONE
REAL (KIND=real32),    INTENT(IN)  :: t, p
REAL (KIND=real32),    INTENT(OUT) :: qs
REAL (KIND=real32)                 :: t_arr(1), p_arr(1), qs_arr(1)
t_arr(1) = t
p_arr(1) = p
CALL qsat_real32_1D(qs_arr,t_arr,p_arr,1)
qs = qs_arr(1)
RETURN
END SUBROUTINE qsat_real32_scalar

SUBROUTINE qsat_wat_real32_scalar(qs, t, p)
IMPLICIT NONE
REAL (KIND=real32),    INTENT(IN)  :: t, p
REAL (KIND=real32),    INTENT(OUT) :: qs
REAL (KIND=real32)                 :: t_arr(1), p_arr(1), qs_arr(1)
t_arr(1) = t
p_arr(1) = p
CALL qsat_wat_real32_1D(qs_arr,t_arr,p_arr,1)
qs = qs_arr(1)
RETURN
END SUBROUTINE qsat_wat_real32_scalar

SUBROUTINE qsat_mix_real32_scalar(qs, t, p)
IMPLICIT NONE
REAL (KIND=real32),    INTENT(IN)  :: t, p
REAL (KIND=real32),    INTENT(OUT) :: qs
REAL (KIND=real32)                 :: t_arr(1), p_arr(1), qs_arr(1)
t_arr(1) = t
p_arr(1) = p
CALL qsat_mix_real32_1D(qs_arr,t_arr,p_arr,1)
qs = qs_arr(1)
RETURN
END SUBROUTINE qsat_mix_real32_scalar

SUBROUTINE qsat_wat_mix_real32_scalar(qs, t, p)
IMPLICIT NONE
REAL (KIND=real32),    INTENT(IN)  :: t, p
REAL (KIND=real32),    INTENT(OUT) :: qs
REAL (KIND=real32)                 :: t_arr(1), p_arr(1), qs_arr(1)
t_arr(1) = t
p_arr(1) = p
CALL qsat_wat_mix_real32_1D(qs_arr,t_arr,p_arr,1)
qs = qs_arr(1)
RETURN
END SUBROUTINE qsat_wat_mix_real32_scalar

END MODULE qsat_mod
