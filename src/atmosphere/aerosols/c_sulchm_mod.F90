! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! A module containing constants/parameters for Sulphur Cycle Chemistry
!
MODULE c_sulchm_mod

IMPLICIT NONE

! Description:
!   This module contains constants for the
!   Sulphur Cycle Chemistry
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Aerosols
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 standards.
!----------------------------------------------------------------------

! Timescale for dissolved SO4 to evaporate (5 mins)
REAL, PARAMETER :: evaptau = 300.0 ! (s)

! Timescale for accumulation mode particles to nucleate once they enter a cloud 
REAL, PARAMETER :: nuctau = 30.0 ! (s)

! Sulphur/2*Nitrogen
REAL, PARAMETER :: relm_s_2n = 3.206 / 2.80

! Reaction rate for DMS+OH  (cc/mcl/s)
REAL, PARAMETER :: k_dms_oh = 9.1e-12

! Rate coeff for CH3SO2+O3 -> CH3SO3+O2 (cc/mcl/s)
REAL, PARAMETER :: k4_ch3so2_o3 = 1.0e-14 

! Rate coeff for CH3SO3+HO2 -> MSA+O2
REAL, PARAMETER :: k5_ch3so3_ho2 = 4.0e-11 

! High pressure reaction rate limit (cc/mcl/s) from STOCHEM model
REAL, PARAMETER :: k_so2oh_hi = 2.0e-12

! Branching ratio
REAL, PARAMETER :: brat_so2 = 0.9            ! SO2 in DMS oxidn
REAL, PARAMETER :: brat_msa = 1.0 - brat_so2 ! MSA in DMS oxidn

! Relative atomic mass
REAL, PARAMETER :: relm_s_h2o2 = 3.206 / 3.40 ! Sulphur/RMM_H2O2


! Power of temp dependence of K_SO2OH_LO
REAL, PARAMETER :: parh = 3.3

! Parameters for calcn of low pressure reaction rate of SO2 with OH 
REAL, PARAMETER :: k1 = 4.0e-31   ! (cc/mcl)2/s
REAL, PARAMETER :: t1 = 300.0     ! K

! Parameters for calcn of HO2 + HO2 reaction rate
REAL, PARAMETER :: k2 = 2.2e-13
REAL, PARAMETER :: t2 = 600.0     ! K
REAL, PARAMETER :: k3 = 1.9e-33
REAL, PARAMETER :: t3 = 890.0     ! K
REAL, PARAMETER :: k4 = 1.4e-21
REAL, PARAMETER :: t4 = 2200.0    ! K


! Parameters for interpolation between LO and HI reaction rate limits
REAL, PARAMETER :: fc = 0.45
REAL, PARAMETER :: fac1 = 1.1904 ! ( 0.75-1.27*LOG10(FC)

! Air parcel lifetime in cloud (3 hours)
REAL, PARAMETER :: cloudtau = 1.08e4 ! s

! Chem lifetime in cloud before oxidn (15 mins)
REAL, PARAMETER :: chemtau = 9.0e2 ! s

! Min mmr of O3 required for oxidn (kg/kg, equiv. 10ppbv)
REAL, PARAMETER :: o3_min = 1.6e-8

! Threshold for cloud liquid water (kg/kg)
REAL, PARAMETER :: thold = 1.0e-8 

! Threshold concn of accu mode particles below which PSI=1 (m-3)
REAL, PARAMETER :: num_star = 1.0e6

! Geometric standard devn of the Aitken mode distn
REAL, PARAMETER :: sigma_ait = 1.30

! Mean radius of particles (m)
REAL, PARAMETER :: rad_ait = 6.5e-9  ! Aitken mode
REAL, PARAMETER :: rad_acc = 95.0e-9 ! Acccumulation mode

! Mole fraction of S in particle
REAL, PARAMETER :: chi = 32.0 / 132.0

! Standard devn of particle size distn for accum mode
REAL, PARAMETER :: sigma = 1.4

END MODULE c_sulchm_mod

