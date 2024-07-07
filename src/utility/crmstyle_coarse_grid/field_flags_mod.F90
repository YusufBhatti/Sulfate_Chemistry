! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Flags to indicate whether field found for next time

MODULE field_flags_mod

IMPLICIT NONE
SAVE

! Description:
!   Module containing flags for fields required
!
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CRMstyle_coarse_grid
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3  programming standards.
!
! Declarations:

! Prognotsics required
LOGICAL ::         & ! Stashcode - field
  l_u              & !  2      U wind
 ,l_v              & !  3      v wind
 ,l_w              & !  150    w wind
 ,l_theta          & !  4      potential temperature
 ,l_q              & !  10     water vapour
 ,l_qcl            & !  254    cloud water
 ,l_qcf            & !  12     cloud ice
 ,l_qrain          & !  272    rain
 ,l_qgraup         & !  273    graupel
 ,l_ptheta           !  408    pressure on theta levels

! Tendencies
LOGICAL ::         & ! Stashcode - field
  l_dt1            & ! 1181    dT SW radiation
 ,l_dt2            & ! 2181    dT LW radiation
 ,l_dt4            & ! 4181    dT microphysics
 ,l_dt9            & ! 9181    dT cloud & BL
 ,l_dt12           & ! 12181   dT advection
 ,l_dt30           & ! 30181   dT total - not in all runs ?
 ,l_dq4            & ! 4181    dq microphysics
 ,l_dq9            & ! 9181    dq cloud & BL
 ,l_dq12           & ! 12181   dq advection
 ,l_dq30           & ! 30181   dq total - not in all runs ?
 ,l_dqcl4          & ! 4181    dqcl microphysics
 ,l_dqcl9          & ! 9181    dqcl cloud & BL
 ,l_dqcl12         & ! 12181   dqcl advection
 ,l_dqcl30         & ! 30181   dqcl total - not in all runs ?
 ,l_dqcf4          & ! 4181    dqcf microphysics
 ,l_dqcf3          & ! 3181    dqcf BL
 ,l_dqcf12         & ! 12181   dqcf advection
 ,l_dqcf30           ! 30181   dqcf total - not in all runs ?

! single fields
LOGICAL ::         & ! Stashcode - field
  l_rain           & ! 4203   surface rain
 ,l_snow           & ! 4204   surface snow
 ,l_precip         & ! 5216   - indicate both of above
 ,l_sh             & ! 3217   sensible heat
 ,l_lh             & ! 3234   latent heat
 ,l_zh             & ! 25     Boundary layer depth
 ,l_tstar          & ! 24     surface temperature
 ,l_pstar          & ! 409    surface pressure
 ,l_orog           & ! 33  orography
 ,l_landsea          ! 505 fractional land


LOGICAL ::         & !
  l_got_fields       ! indicate whether found all required fields

END MODULE field_flags_mod
