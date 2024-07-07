! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Global standard conversions

MODULE conversions_mod

! Description:
!   Module to pull in, and rename, the invariant physical constants used for
!   conversions from one scale to another from the external shared library.
!
! This file should not have any constants permanently defined in it.
! However, it is acceptable to add new constants here and notify the UM
! System team you have a new constant to be moved into Shumlib before
! requesting a review
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Constants

USE f_shum_conversions_mod, ONLY:                                             &
    ! Time conversions.
    rsec_per_day         => shum_rsec_per_day_const,                          &
    isec_per_day         => shum_isec_per_day_const,                          &
    rsec_per_hour        => shum_rsec_per_hour_const,                         &
    isec_per_hour        => shum_isec_per_hour_const,                         &
    isec_per_min         => shum_isec_per_min_const,                          &
    rhour_per_day        => shum_rhour_per_day_const,                         &
    ihour_per_day        => shum_ihour_per_day_const,                         &
    rhour_per_sec        => shum_rhour_per_sec_const,                         &
    rday_per_hour        => shum_rday_per_hour_const,                         &
    ! Angle conversions.
    pi                   => shum_pi_const,                                    &
    pi_over_180          => shum_pi_over_180_const,                           &
    recip_pi_over_180    => shum_180_over_pi_const,                           &
    ! Temperature conversions.
    zerodegc             => shum_zerodegc_const,                              &
    ! Speed conversions.
    kt2ms                => shum_kt2ms_const,                                 &
    ! Length conversions.
    ft2m                 => shum_ft2m_const,                                  &
    ! Time conversions in 32Bit form.
    rsec_per_day_32      => shum_rsec_per_day_const_32,                       &
    isec_per_day_32      => shum_isec_per_day_const_32,                       &
    rsec_per_hour_32     => shum_rsec_per_hour_const_32,                      &
    isec_per_hour_32     => shum_isec_per_hour_const_32,                      &
    isec_per_min_32      => shum_isec_per_min_const_32,                       &
    rhour_per_day_32     => shum_rhour_per_day_const_32,                      &
    ihour_per_day_32     => shum_ihour_per_day_const_32,                      &
    rhour_per_sec_32     => shum_rhour_per_sec_const_32,                      &
    rday_per_hour_32     => shum_rday_per_hour_const_32,                      &
    ! Angle conversions in 32Bit form.
    pi_32                => shum_pi_const_32,                                 &
    pi_over_180_32       => shum_pi_over_180_const_32,                        &
    recip_pi_over_180_32 => shum_180_over_pi_const_32,                        &
    ! Temperature conversions in 32Bit form.
    zerodegc_32          => shum_zerodegc_const_32,                           &
    ! Speed conversions in 32Bit form.
    kt2ms_32             => shum_kt2ms_const_32,                              &
    ! Length conversions in 32Bit form.
    ft2m_32              => shum_ft2m_const_32

IMPLICIT NONE

! Define new conversions constants here :

END MODULE conversions_mod
