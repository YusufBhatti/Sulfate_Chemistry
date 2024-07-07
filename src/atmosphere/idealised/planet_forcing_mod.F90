! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Control subroutines for the suite of idealised dry planet simulations
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE planet_forcing_mod
  
  USE tforce_mod
  USE trelax_mod
  
  IMPLICIT NONE

CONTAINS

  FUNCTION l_planet_forcing(tforce_number)
    
    IMPLICIT NONE
!
! Description:
!     Returns .TRUE. if this module should be used for the forcing
!     profile, .FALSE. otherwise.

!   Input
    INTEGER :: tforce_number

!   Output
    LOGICAL :: l_planet_forcing
  
    l_planet_forcing =                                                  &
        tforce_number == tf_TLE .OR.                                    &
        tforce_number == tf_EL .OR.                                     &
        tforce_number == tf_SHJ .OR.                                    &
        tforce_number == tf_HD209458b_Heng .OR.                         &
        tforce_number == tf_HD209458b_Heng_smooth .OR.                  &
        tforce_number == tf_HD209458b_iro .OR.                          &
        tforce_number == tf_Y_Dwarf .OR.                                &
        tforce_number == tf_GJ1214b .OR.                                &
        tforce_number == tf_GJ1214b_dT800
  
  END FUNCTION l_planet_forcing
  
  FUNCTION l_planet_relax(trelax_number)
  
    IMPLICIT NONE
!
! Description:
!     Returns .TRUE. if this module should be used for the forcing
!     relaxation timescale, .FALSE. otherwise.

!   Input
    INTEGER :: trelax_number

!   Output
    LOGICAL :: l_planet_relax
  
    l_planet_relax =                                                    &
        trelax_number == tr_EL .OR.                                     &
        trelax_number == tr_SHJ .OR.                                    &
        trelax_number == tr_HD209458b_Iro .OR.                          &
        trelax_number == tr_Y_Dwarf .OR.                                &
        trelax_number == tr_GJ1214b
  
  END FUNCTION l_planet_relax

  FUNCTION planet_forcing_theta(tforce_number, i, j, k,                 &
            exner_theta_levels) RESULT(theta_out)
    
    USE atm_fields_bounds_mod
!   Functions to calculate the forcing profile
    USE tidally_locked_earth_forcing_mod,                               &
          ONLY: tidally_locked_earth_theta
    USE earth_like_forcing_mod, ONLY: earth_like_theta
    USE shallow_hot_jupiter_forcing_mod, ONLY: shallow_hot_jupiter_theta
    USE hd209458b_forcing_mod, ONLY: hd209458b_theta
    USE gj1214b_forcing_mod, ONLY: gj1214b_theta
    USE y_dwarf_forcing_mod, ONLY: y_dwarf_theta
  
    IMPLICIT NONE
!
! Description:
!       Wrapper for call to calculation of forcing theta profile

    ! Input
    INTEGER :: tforce_number
    INTEGER :: i,j,k
    REAL    :: exner_theta_levels(tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end) 
    
!   Output
    REAL :: theta_out
    
    IF (tforce_number == tf_TLE) THEN
       theta_out=tidally_locked_earth_theta(i, j, k, exner_theta_levels)
    ELSE IF (tforce_number == tf_EL) THEN
       theta_out=earth_like_theta(i, j, k, exner_theta_levels)
    ELSE IF (tforce_number == tf_SHJ) THEN
       theta_out=shallow_hot_jupiter_theta(i, j, k, exner_theta_levels)
    ELSE IF (tforce_number == tf_HD209458b_Heng .OR.     &
         tforce_number == tf_HD209458b_Heng_smooth .OR.  &
         tforce_number == tf_HD209458b_iro) THEN
       theta_out=hd209458b_theta(i, j, k, exner_theta_levels, &
            tforce_number)
    ELSE IF (tforce_number == tf_GJ1214b .OR.                                &
             tforce_number == tf_GJ1214b_dT800) THEN
       theta_out=gj1214b_theta(i, j, k, exner_theta_levels,tforce_number)
    ELSE IF (tforce_number == tf_Y_Dwarf) THEN
       theta_out=y_dwarf_theta(i, j, k, exner_theta_levels)
    END IF
  
  END FUNCTION planet_forcing_theta
  
  FUNCTION planet_recip_newt(trelax_number, i, j, k,                    &
      exner_theta_levels) RESULT(recip_tscale_out)
    
    USE atm_fields_bounds_mod
!   Functions to calculate the relaxation timescale
    USE earth_like_forcing_mod, ONLY: earth_like_recip_newt
    USE shallow_hot_jupiter_forcing_mod,                                &
          ONLY: shallow_hot_jupiter_recip_newt
    USE hd209458b_forcing_mod, ONLY: hd209458b_recip_newt
    USE gj1214b_forcing_mod, ONLY: gj1214b_recip_newt
    USE y_dwarf_forcing_mod, ONLY: y_dwarf_recip_newt
  
    IMPLICIT NONE
!
! Description:
!        Wrapper for call to calculation of relaxation time scale

! Input
    INTEGER :: trelax_number
    INTEGER :: i,j,k
    REAL    :: exner_theta_levels(tdims_s%i_start:tdims_s%i_end,        &
         tdims_s%j_start:tdims_s%j_end, tdims_s%k_start:tdims_s%k_end) 
    
! Output
    REAL :: recip_tscale_out
  
    IF (trelax_number == tr_EL) THEN
       recip_tscale_out=earth_like_recip_newt()
    ELSE IF (trelax_number == tr_SHJ) THEN
       recip_tscale_out=shallow_hot_jupiter_recip_newt()
    ELSE IF (trelax_number == tr_HD209458b_Iro) THEN
       recip_tscale_out=hd209458b_recip_newt(i, j, k,                    &
            exner_theta_levels, trelax_number) 
    ELSE IF (trelax_number == tr_GJ1214b) THEN 
       recip_tscale_out=gj1214b_recip_newt(i, j, k,                      &
            exner_theta_levels)
    ELSE IF (trelax_number == tr_Y_Dwarf) THEN
       recip_tscale_out=y_dwarf_recip_newt(i, j, k, exner_theta_levels)
    END IF
  
  END FUNCTION planet_recip_newt
  
END MODULE planet_forcing_mod
