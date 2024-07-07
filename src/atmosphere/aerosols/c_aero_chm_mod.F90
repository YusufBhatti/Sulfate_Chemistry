! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Aerosols
MODULE c_aero_chm_mod

IMPLICIT NONE
! C_AERO_CHM start
! Contains constants required for aerosol conversion and nucleation
! scavenging by cloud droplets.
! same values apply for soot, smoke, ammonium nitrate and ocff
!
!     air parcel lifetime in cloud
REAL,PARAMETER:: cloudtau = 1.08e4            ! secs (=3 hours)
!
!     timescale for suspended aerosol to evaporate
REAL,PARAMETER:: evaptau = 300.0              ! secs  (=5 mins)
!
!     timescale for accumulation mode particles
REAL,PARAMETER:: nuctau = 30.0                ! secs
!
!     Cloud liquid water threshold for nucleation scavenging to occur.
REAL,PARAMETER:: thold = 1.0e-8               ! kg/kg
!
! C_AERO_CHM end

END MODULE c_aero_chm_mod
