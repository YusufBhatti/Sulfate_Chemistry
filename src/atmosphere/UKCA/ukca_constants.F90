! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module to contain constants used in UKCA
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_constants

USE planet_constants_mod, ONLY: vkarmn => vkman, gg => g, ra => r

USE water_constants_mod,  ONLY: rho_water

USE missing_data_mod,     ONLY: rmdi, imdi

USE chemistry_constants_mod, ONLY: avogadro, boltzmann, rho_so4

USE c_rmol, ONLY: rmol


IMPLICIT NONE


! UKCA_ MODE constants

! MODE constants already defined in the UM
! 1) Original values
!      REAL, PARAMETER :: AVC=6.022E23
!      REAL, PARAMETER :: ZBOLTZ=1.3807E-23
!      REAL, PARAMETER :: VKARMN=0.4
!      REAL, PARAMETER :: RA=287.05
!      REAL, PARAMETER :: GG=9.80665
!      REAL, PARAMETER :: RAD_E=6371229.0
!      REAL, PARAMETER :: MM_DA=AVC*ZBOLTZ/RA
!      REAL, PARAMETER :: RHOSUL=1800.0E0      ! UM is slightly different
!      REAL, PARAMETER :: MMSUL=0.098E0
!      REAL, PARAMETER :: RR=8.314
!      PARAMETER(BCONST=3.44E13)        ! Volume_mode
!      PARAMETER(NU_H2SO4=3.0)          ! Volume_mode
!      PARAMETER(CONVERT=1.0e-18)       ! Volume_mode
!      PARAMETER(RHOW=1000.0,MMW=0.018) ! Volume_mode
!      PARAMETER(RHOSUL=1800.0)`        ! Volume_mode

! 2) UM definitions of MODE constants
REAL, PARAMETER :: avc=Avogadro
REAL, PARAMETER :: recip_avogad=1.0/Avogadro
REAL, PARAMETER :: zboltz=Boltzmann
REAL, PARAMETER :: rr=rmol
REAL, PARAMETER :: rhosul=rho_so4       ! 1769 or 1800 ?
REAL, PARAMETER :: rhow=rho_water       ! Water density (kg/m^3)

! 3) MODE/aerosol chemistry constants not already in UM
REAL, PARAMETER :: nmol=1.0e2          !
REAL, PARAMETER :: tdays=0.0           !
REAL, PARAMETER :: ems_eps=1.0e-8      !
REAL, PARAMETER :: conc_eps=1.0e-8     !
REAL, PARAMETER :: dn_eps=1.0e-8       !
REAL, PARAMETER :: bconst=3.44e13      ! Volume_mode
REAL, PARAMETER :: nu_h2so4=3.0        ! Volume_mode
REAL, PARAMETER :: convert=1.0e-18     ! Volume_mode
REAL, PARAMETER :: H_plus=1.0e-5       ! cloud/rain H+ concentration


! Dobson Unit (molecules per m2)
REAL, PARAMETER :: dobson = 2.685e20

! molecular masses in kg/mol
LOGICAL, PARAMETER  ::  l_ukca_diurnal_isopems   = .TRUE.
REAL, PARAMETER :: m_air        = 0.02897       ! Air
REAL, PARAMETER :: mmsul        = 0.09808       ! H2SO4
REAL, PARAMETER :: mmw          = 0.0180154     ! H2O


!  -------------------------------------------------------------------
!  Conversion factor from vmr to mmr for each species.
!             vmr*c_species = mmr
!             c_species = m_species/m_air  (m_air = 28.97)
!
!  Followed by molecular masses for each species (e.g. m_ch4) in g/mol
!  -------------------------------------------------------------------

REAL, PARAMETER :: c_o3p        = 0.5523
REAL, PARAMETER :: c_o1d        = 0.5523
REAL, PARAMETER :: c_o3         = 1.657
REAL, PARAMETER :: c_no         = 1.036
REAL, PARAMETER :: c_no3        = 2.140
REAL, PARAMETER :: c_no2        = 1.588
REAL, PARAMETER :: c_n2o5       = 3.728
REAL, PARAMETER :: c_ho2no2     = 2.727
REAL, PARAMETER :: c_hono2      = 2.175
REAL, PARAMETER :: c_hno3       = c_hono2   ! used in nitrate.F90
REAL, PARAMETER :: c_oh         = 0.5868
REAL, PARAMETER :: c_ho2        = 1.139
REAL, PARAMETER :: c_h2         = 0.06904
REAL, PARAMETER :: c_h2o2       = 1.174
REAL, PARAMETER :: c_ch4        = 0.5523
REAL, PARAMETER :: c_c          = 0.4142
REAL, PARAMETER :: c_co         = 0.9665
REAL, PARAMETER :: c_co2        = 1.5188
REAL, PARAMETER :: c_hcho       = 1.036
REAL, PARAMETER :: C_MeOO       = 1.622
REAL, PARAMETER :: c_h2o        = 0.6213
REAL, PARAMETER :: c_h2os       = 0.6213
REAL, PARAMETER :: C_MeOOH      = 1.657
REAL, PARAMETER :: c_hono       = 1.622
REAL, PARAMETER :: c_o2         = 1.105
REAL, PARAMETER :: c_n2         = 0.9665
REAL, PARAMETER :: c_c2h6       = 1.036
REAL, PARAMETER :: C_EtOO       = 2.106
REAL, PARAMETER :: C_EtOOH      = 2.140
REAL, PARAMETER :: C_MeCHO      = 1.519
REAL, PARAMETER :: c_toth       = 1.000
REAL, PARAMETER :: C_MeCO3      = 2.589
REAL, PARAMETER :: c_pan        = 4.177
REAL, PARAMETER :: c_c3h8       = 1.519
REAL, PARAMETER :: C_PrOO       = 2.589
REAL, PARAMETER :: C_PrOOH      = 2.623
REAL, PARAMETER :: C_EtCHO      = 2.002
REAL, PARAMETER :: C_EtCO3      = 3.072
REAL, PARAMETER :: C_Me2CO      = 2.002
REAL, PARAMETER :: C_MeCOCH2OO  = 3.072
REAL, PARAMETER :: C_MeCOCH2OOH = 3.107
REAL, PARAMETER :: c_ppan       = 4.660
REAL, PARAMETER :: C_MeONO2     = 2.658
REAL, PARAMETER :: c_n          = 0.48325
REAL, PARAMETER :: c_h          = 0.03452
REAL, PARAMETER :: c_n2o        = 1.5188
REAL, PARAMETER :: C_CFCl3      = 4.7480
REAL, PARAMETER :: C_CF2Cl2     = 4.1783
REAL, PARAMETER :: C_ClO        = 1.7784
REAL, PARAMETER :: C_HCl        = 1.2604
REAL, PARAMETER :: C_ClONO2     = 3.3668
REAL, PARAMETER :: C_HOCl       = 1.8129
REAL, PARAMETER :: C_OClO       = 2.3309
REAL, PARAMETER :: C_BrO        = 3.315
REAL, PARAMETER :: C_BrONO2     = 4.9034
REAL, PARAMETER :: C_HBr        = 2.7970
REAL, PARAMETER :: C_HOBr       = 3.3495
REAL, PARAMETER :: C_BrCl       = 3.9884
REAL, PARAMETER :: C_MeBr       = 3.2805
REAL, PARAMETER :: c_so2        = 2.2112
REAL, PARAMETER :: c_so3        = 2.7615
REAL, PARAMETER :: C_Me2S       = 2.145
REAL, PARAMETER :: c_dms        = 2.145
REAL, PARAMETER :: c_dmso       = 2.6965
REAL, PARAMETER :: c_ocs        = 2.0711
REAL, PARAMETER :: c_cos        = 2.0711
REAL, PARAMETER :: c_h2s        = 1.1766
REAL, PARAMETER :: c_cs2        = 2.6282
REAL, PARAMETER :: c_msia       = 2.7650 !LER, UC, Feb 2019
REAL, PARAMETER :: c_sad        = 4.1255
REAL, PARAMETER :: c_msa        = 3.317
REAL, PARAMETER :: c_s          = 1.1046
REAL, PARAMETER :: c_h2so4      = 3.385
REAL, PARAMETER :: C_CF2ClCFCl2 = 6.4722
REAL, PARAMETER :: C_CHF2Cl     = 2.9858
REAL, PARAMETER :: C_MeCCl3     = 4.6082
REAL, PARAMETER :: C_CCl4       = 5.3158
REAL, PARAMETER :: C_MeCl       = 1.7432
REAL, PARAMETER :: C_CF2ClBr    = 5.7128
REAL, PARAMETER :: C_CF3Br      = 5.1432
REAL, PARAMETER :: C_Cl         = 1.2261
REAL, PARAMETER :: C_Cl2O2      = 3.5568
REAL, PARAMETER :: C_Br         = 2.7627
REAL, PARAMETER :: C_CH2Br2     = 6.0013
REAL, PARAMETER :: c_mecf2cl    = 3.4673
REAL, PARAMETER :: c_cf2br2     = 7.2489
REAL, PARAMETER :: c_cf2brcf2br = 8.9748
REAL, PARAMETER :: c_cf2clcf3   = 5.3314
REAL, PARAMETER :: c_cf2clcf2cl = 5.8992
REAL, PARAMETER :: c_mecfcl2    = 4.0352
REAL, PARAMETER :: c_c5h8       = 2.3473
REAL, PARAMETER :: c_iso2       = 4.0387
REAL, PARAMETER :: c_isooh      = 4.0732
REAL, PARAMETER :: c_ison       = 5.3504  ! Might be revised to 5.0052
!                       for RAQ chem where ISON = (NO3)C4H6CHO: C5H7NO4: 145
REAL, PARAMETER :: c_macr       = 2.4163
REAL, PARAMETER :: c_macro2     = 4.1077
REAL, PARAMETER :: c_macrooh    = 4.1422
REAL, PARAMETER :: c_mpan       = 5.0742
REAL, PARAMETER :: c_hacet      = 2.5544
REAL, PARAMETER :: c_mgly       = 2.4853
REAL, PARAMETER :: c_nald       = 3.6244
REAL, PARAMETER :: c_hcooh      = 1.5878
REAL, PARAMETER :: c_meco3h     = 2.6234
REAL, PARAMETER :: c_meco2h     = 2.0711
REAL, PARAMETER :: c_nh3        = 0.5879
REAL, PARAMETER :: c_monoterp   = 4.7034
REAL, PARAMETER :: c_sec_org    = 5.1782    ! Molecular weight=150.

!     species required for ukca_scenario_rcp
REAL, PARAMETER :: C_cf3chf2    = 4.1429 ! HFC125
REAL, PARAMETER :: C_ch2fcf3    = 3.5219 ! HFC134a
REAL, PARAMETER :: C_chbr3      = 8.7238


!     Extra species for RAQ chemistry
REAL, PARAMETER :: c_c4h10      = 2.0021
REAL, PARAMETER :: c_c2h4       = 0.9665
REAL, PARAMETER :: c_c3h6       = 1.4498
REAL, PARAMETER :: c_rnc2h4     = 3.6244
REAL, PARAMETER :: c_rnc3h6     = 4.1077
!     rnc2h4 & rnc3h6 are CH2(NO3)CHO & CH3CH(NO3)CHO
REAL, PARAMETER :: c_ch3oh      = 1.1046
REAL, PARAMETER :: c_meoh       = 1.1046
REAL, PARAMETER :: c_toluene    = 3.1757
REAL, PARAMETER :: c_oxylene    = 3.6590
REAL, PARAMETER :: c_memald     = 3.3828
!     MEMALDIAL is CH3-CO-CH=CH-CHO, ring fragmentation
!     product from degradation of toluene and o-xylene
REAL, PARAMETER :: c_buooh      = 3.1067
REAL, PARAMETER :: c_mek        = 2.4853
REAL, PARAMETER :: c_mvk        = 2.4163
REAL, PARAMETER :: c_mvkooh     = 4.1422
REAL, PARAMETER :: c_gly        = 2.0021
REAL, PARAMETER :: c_rnc5h8     = 5.0052
REAL, PARAMETER :: c_orgnit     = 5.5230
REAL, PARAMETER :: c_buoo       = 3.0721
REAL, PARAMETER :: c_meko2      = 3.5554
REAL, PARAMETER :: c_hoc2h4o2   = 2.6579
REAL, PARAMETER :: c_hoc3h6o2   = 3.1412
REAL, PARAMETER :: c_hoipo2     = 4.0387
REAL, PARAMETER :: c_tolp1      = 3.7280
REAL, PARAMETER :: c_oxyl1      = 4.2113
REAL, PARAMETER :: c_memald1    = 3.9351
REAL, PARAMETER :: c_homvko2    = 4.1077

!     Extra species for EXTTC chemistry
REAL, PARAMETER :: c_isoo       = 4.0387
REAL, PARAMETER :: c_macroo     = 4.1077
REAL, PARAMETER :: c_apin       = 4.9645
REAL, PARAMETER :: C_PropeOO    = 2.5198
REAL, PARAMETER :: C_PropeOOH   = 2.5544
REAL, PARAMETER :: C_OnitU      = 3.5554
REAL, PARAMETER :: c_mekoo      = 3.5554
REAL, PARAMETER :: c_mekooh     = 3.5889
REAL, PARAMETER :: C_EteOO      = 2.0366
REAL, PARAMETER :: c_alka       = 2.0021
REAL, PARAMETER :: c_alkaoo     = 3.0721
REAL, PARAMETER :: c_alkaooh    = 3.1067
REAL, PARAMETER :: c_arom       = 3.4173
REAL, PARAMETER :: c_aromoo     = 4.4874
REAL, PARAMETER :: c_aromooh    = 4.5219
REAL, PARAMETER :: c_bsvoc1     = 4.9645   ! as APIN
REAL, PARAMETER :: c_bsvoc2     = 4.9645   ! as APIN
REAL, PARAMETER :: C_B2ndry     = 4.9645   ! as APIN
REAL, PARAMETER :: c_asvoc1     = 3.4173   ! as AROM
REAL, PARAMETER :: c_asvoc2     = 3.4173   ! as AROM
REAL, PARAMETER :: C_A2ndry     = 3.4173   ! as AROM
REAL, PARAMETER :: c_bsoa       = 5.1778   ! 150.0
REAL, PARAMETER :: c_asoa       = 5.1778   ! 150.0
REAL, PARAMETER :: c_isosvoc1   = 2.3473   ! as C5H8
REAL, PARAMETER :: c_isosvoc2   = 2.3473   ! as C5H8
REAL, PARAMETER :: c_isosoa     = 4.4874   ! 130.0

!     molecular masses in g/mol of emitted species,
!     for budget calculations
REAL, PARAMETER :: m_ho2     =  33.007
REAL, PARAMETER :: m_ch4     =  16.0
REAL, PARAMETER :: m_co      =  28.0
REAL, PARAMETER :: m_hcho    =  30.0
REAL, PARAMETER :: m_c2h6    =  30.0
REAL, PARAMETER :: m_c3h8    =  44.0
REAL, PARAMETER :: m_mecho   =  44.0
REAL, PARAMETER :: m_no2     =  46.0
REAL, PARAMETER :: m_n2o5    = 108.01
REAL, PARAMETER :: m_me2co   =  58.0
REAL, PARAMETER :: m_isop    =  68.0
REAL, PARAMETER :: m_no      =  30.0
REAL, PARAMETER :: m_n       =  14.0
REAL, PARAMETER :: m_c       =  12.0
REAL, PARAMETER :: m_monoterp =  136.24

!     molecular masses of stratospheric species, for which surface
!     mmrs are prescribed
REAL, PARAMETER :: m_hcl        =  36.5
REAL, PARAMETER :: m_n2o        =  44.0
REAL, PARAMETER :: m_clo        =  51.5
REAL, PARAMETER :: m_hocl       =  52.5
REAL, PARAMETER :: m_oclo       =  67.5
REAL, PARAMETER :: m_clono2     =  97.5
REAL, PARAMETER :: m_cf2cl2     = 121.0
REAL, PARAMETER :: m_cfcl3      = 137.5
REAL, PARAMETER :: m_hbr        =  81.0
REAL, PARAMETER :: m_mebr       =  95.0
REAL, PARAMETER :: m_bro        =  96.0
REAL, PARAMETER :: m_hobr       =  97.0
REAL, PARAMETER :: m_brcl       = 115.5
REAL, PARAMETER :: m_brono2     = 142.0
REAL, PARAMETER :: m_cf2clcfcl2 = 187.5
REAL, PARAMETER :: m_chf2cl     =  86.5
REAL, PARAMETER :: m_meccl3     = 133.5
REAL, PARAMETER :: m_ccl4       = 154.0
REAL, PARAMETER :: m_mecl       =  50.5
REAL, PARAMETER :: m_cf2clbr    = 165.5
REAL, PARAMETER :: m_cf3br      = 149.0
REAL, PARAMETER :: m_ch2br2     = 173.835

! sulphur containing, etc.
REAL, PARAMETER :: m_ocs        =  60.0
REAL, PARAMETER :: m_cos        =  60.0
REAL, PARAMETER :: m_h2s        =  34.086
REAL, PARAMETER :: m_cs2        =  76.14
REAL, PARAMETER :: m_dms        =  62.1
REAL, PARAMETER :: m_dmso       =  78.13
REAL, PARAMETER :: m_me2s       =  62.1
REAL, PARAMETER :: m_msa        =  96.1
REAL, PARAMETER :: m_sec_org    =  150.0
REAL, PARAMETER :: m_s          =  32.07
REAL, PARAMETER :: m_so2        =  64.06
REAL, PARAMETER :: m_so3        =  80.06
REAL, PARAMETER :: m_so4        =  96.06
REAL, PARAMETER :: m_h2so4      =  98.07
REAL, PARAMETER :: m_nh3        =  17.03
REAL, PARAMETER :: m_nh42so4    =  132.16

REAL, PARAMETER :: m_cl         = 35.5
REAL, PARAMETER :: m_cl2o2      = 103.0
REAL, PARAMETER :: m_br         = 80.0
REAL, PARAMETER :: m_h2         = 2.016
REAL, PARAMETER :: m_h2o        = 18.0
REAL, PARAMETER :: m_mecoch2ooh = 90.0
REAL, PARAMETER :: m_isooh      = 118.0
REAL, PARAMETER :: m_mpan       = 147.0
REAL, PARAMETER :: m_ppan       = 135.0

!     Extra masses for RAQ or other chemistries
REAL, PARAMETER :: m_c5h8    =  68.0
REAL, PARAMETER :: m_c4h10   =  58.0
REAL, PARAMETER :: m_c2h4    =  28.0
REAL, PARAMETER :: m_c3h6    =  42.0
REAL, PARAMETER :: m_toluene =  92.0
REAL, PARAMETER :: m_oxylene = 106.0
REAL, PARAMETER :: m_ch3oh   =  32.0
REAL, PARAMETER :: m_meoh    =  32.0
REAL, PARAMETER :: m_buooh   =  90.0
REAL, PARAMETER :: m_mvkooh  = 120.0
REAL, PARAMETER :: m_orgnit  = 160.0
REAL, PARAMETER :: m_macrooh = 120.0
REAL, PARAMETER :: m_gly     =  58.0

REAL, PARAMETER :: m_hono = 47.0
!     Extra masses for Wesely scheme
REAL, PARAMETER :: m_macr    = 70.0
REAL, PARAMETER :: m_etcho   = 58.0
REAL, PARAMETER :: m_nald    = 105.0
REAL, PARAMETER :: m_mgly    = 72.0
REAL, PARAMETER :: m_hacet  = 74.0


!     Extra masses for EXTTC chemistry
REAL, PARAMETER :: m_apin     =  136.0
REAL, PARAMETER :: m_mvk      =  70.0
REAL, PARAMETER :: m_mek      =  72.0
REAL, PARAMETER :: m_alka     =  58.0        ! as butane
REAL, PARAMETER :: m_arom     =  99.0        ! (toluene + xylene)/2
REAL, PARAMETER :: m_bsvoc1   = 144.0
REAL, PARAMETER :: m_bsvoc2   = 144.0
REAL, PARAMETER :: m_asvoc1   = 99.0
REAL, PARAMETER :: m_asvoc2   = 99.0
REAL, PARAMETER :: m_isosvoc1 = 68.0
REAL, PARAMETER :: m_isosvoc2 = 68.0
REAL, PARAMETER :: m_onitu    = 102.0
REAL, PARAMETER :: m_bsoa     = 150.0
REAL, PARAMETER :: m_asoa     = 150.0
REAL, PARAMETER :: m_isosoa   = 130.0
REAL, PARAMETER :: m_alkaooh  = 90.0
REAL, PARAMETER :: m_aromooh  = 130.0
REAL, PARAMETER :: m_mekooh   = 104.0

!     The mass of organic nitrate is an approximation,
!     calculated as the average of ORGNIT formed by two
!     reacs. in UKCA_CHEMCO_RAQ:
!      NO2 + TOLP1 --> ORGNIT (A)
!      NO2 + OXYL1 --> ORGNIT (B)
!      * TOL  = methylbenzene       = C6H5(CH3)
!        OXYL = 1,2-dimethylbenzene = C6H4(CH3)2
!      * TOL  + OH --> TOLP1: C6H4(OH)(CH3)  = methyl phenol
!        OXYL + OH --> OXYL1: C6H3(OH)(CH3)2 = dimethyl phenol
!      * ORGNIT A: TOLP1 + NO2 ~ C6H3(CH3)(OH)NO2  ~
!                  C7H7NO3: methyl nitrophenol   -> 153
!        ORGNIT B: OXYL1 + NO2 ~ C6H2(CH3)2(OH)NO2 ~
!                  C8H9NO3: dimethyl nitrophenol -> 167
!  -------------------------------------------------------------------

END MODULE ukca_constants
