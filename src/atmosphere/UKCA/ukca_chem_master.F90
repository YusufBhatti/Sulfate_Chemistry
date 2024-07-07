! *****************************COPYRIGHT*******************************
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module providing UKCA with inputs describing chemistry
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
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
MODULE ukca_chem_master_mod

IMPLICIT NONE

PRIVATE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_CHEM_MASTER_MOD'

PUBLIC  :: ukca_chem_master

CONTAINS

SUBROUTINE ukca_chem_master

USE ukca_chem_defs_mod, ONLY: chch_t1, ratb_t1, rath_t1, ratj_t1, ratt_t1,&
                              depvel_t, wetdep,                           & 
                              chch_defs_new, ratb_defs_new, rath_defs_new,&
                              ratj_defs_new, ratt_defs_new, depvel_defs_new,&
                              henry_defs_new
USE ukca_option_mod,    ONLY: l_ukca_strattrop, l_ukca_strat, l_ukca_trop,&
                              l_ukca_achem, l_ukca_raq, l_ukca_offline, &
                              l_ukca_tropisop,l_ukca_trophet, L_ukca_het_psc, &
                              l_ukca_stratcfc, &
                              jpspec, jpbk, jptk, jphk, jppj, jpdd, jpdw,jpnr,&
                              jpctr
USE asad_mod,           ONLY: jddepc, jddept
USE carbon_options_mod, ONLY: l_co2_interactive
USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE umprintmgr

IMPLICIT NONE


! In the below items, item_no identifies different versions of the same entry
! (e.g., updates of the same reaction). If e.g. in different schemes, the 
! branching ratios or products are different, these should be identified 
! through different settings of chemistry. If multiple versions are defined, 
! the highest version <= the selected version is used. Versions of items need 
! to be listed together.


! Define chemistry scheme types
! Entries will only be selected if the active chemistry scheme flag matches
! with one of the values in the "scheme" column.
! Only one chemistry scheme flag can be active at a time.
INTEGER, PARAMETER :: ST  = 1   ! stratosphere-troposphere scheme
INTEGER, PARAMETER :: T   = 2   ! troposphere scheme
INTEGER, PARAMETER :: S   = 4   ! stratosphere scheme
INTEGER, PARAMETER :: R   = 8   ! RAQ scheme
INTEGER, PARAMETER :: OL  = 16  ! offline scheme
INTEGER, PARAMETER :: TI  = 32  ! troposphere-isoprene scheme

! Define qualifiers
! Entries with the corresponding qualifier flag in their 'qualifier' column
! will be selected, whilst those with the corresponding flag in their
! 'disqualifier' column will be deselected. Any entry with no value in
! the qualifier and disqualifier columns will be automatically selected.
! Multiple qualifier flags can be used at the same time. However, a match
! in the disqualifier column will take precedent over a match in the
! qualifier column.
INTEGER, PARAMETER :: A  = 1          ! aerosol chemistry
INTEGER, PARAMETER :: TH = 2          ! tropospheric heterogeneous reactions
INTEGER, PARAMETER :: HP = 4          ! heterogeneous PSC chemistry
INTEGER, PARAMETER :: ES = 8          ! extended stratospheric reactions
INTEGER, PARAMETER :: CI = 16         ! 3-D CO2


! define size of master chemistry
INTEGER, PARAMETER :: n_chch_master = 149 ! number of known species
INTEGER, PARAMETER :: n_het_master  =  15 ! number of heterogeneous reactions #Updated from 10 to 15; LER, UC, Apr2019
INTEGER, PARAMETER :: n_dry_master  =  57 ! number of dry deposition reactions
INTEGER, PARAMETER :: n_wet_master  =  52 ! number of wet deposition reactions #Updated from 49 to 52; LER, UC, Apr2019
INTEGER, PARAMETER :: n_bimol_master= 398 ! number of bimolecular reactions #Updated from 400 to 398; LER, UC, Feb2019
INTEGER, PARAMETER :: n_ratj_master =  76 ! number of photolysis reactions
INTEGER, PARAMETER :: n_ratt_master =  49 ! number of termolecular reactions

! Define below the chemistry schemes

! ST scheme follows ATA NLA CheST Chemistry v1.2
! These reactions generally have comments identifying where they come from 
! (e.g., "T001" is the first termolecular reaction. Reactions not labelled are 
! taken from the previous "trop", "tropisop", "strat", or offline schemes, 
! where they were not labelled (hence I'm not sure about the origins of the 
! data). I suggest that eventually these reactions are superseded by the ST 
! scheme.
 
!The columns take the following meanings:
! item number, species name, (unused), species classification, family (unused),
! dry deposition flag, wet deposition flag, emission flag (unused), scheme,
! selection qualifier, deselection / disqualifier, (UM) version when entered 
! into the master mechanism.  
TYPE(chch_t1), PARAMETER :: chch_defs_master(1:n_chch_master)=(/   &
!   1
chch_t1( 1,'O(3P)     ',1,'TR        ','Ox        ',0,0,0,S+ST,0,0,107),&
!   1
chch_t1( 1,'O(3P)     ',1,'SS        ','Ox        ',0,0,0,TI+T+R,0,0,107),&
!   2
chch_t1( 2,'O(1D)     ',1,'SS        ','Ox        ',0,0,0,TI+S+T+ST+R,0,0,107),&
!   3 DD: 1,
chch_t1( 3,'O3        ',1,'TR        ','Ox        ',1,0,0,TI+S+T+ST+R,0,A,107),&
!   3  Wet deposition of ozone was introduced for aerosol chemistry
!      The GLOMAP-mode aerosol scheme needs to know about dissolved O3 in the
!      aqueous phase to be able to perform SO2 in-cloud oxidatuion.
!      Wet deposition is not a loss process for O3.
chch_t1( 3,'O3        ',1,'TR        ','Ox        ',1,1,0,TI+ST+S+T+R,A,0,107),&
!   3
chch_t1( 3,'O3        ',1,'CF        ','Ox        ',0,0,0,OL,0,0,107),&
!   4
chch_t1( 4,'N         ',1,'TR        ','NOx       ',0,0,0,S+ST,0,0,107),&
! NOx is emitted as NO2 for all schemes except RAQ which causes ozone 
! depletion (NO + O3 -> NO2 + O2). 
! This could be made consistent. Also RAQ scheme does not do dry deposition of 
! NO.
!   5 DD: 2,
chch_t1( 5,'NO        ',1,'TR        ','NOx       ',1,0,0,S+ST,0,0,107),&
!   5 DD: 2,
chch_t1( 5,'NO        ',1,'TR        ','NOx       ',1,0,1,TI+T,0,0,107),&
!   5
chch_t1( 5,'NO        ',1,'TR        ','NOx       ',0,0,1,R,0,0,107),&
!   6 DD: 3,WD: 1,
! No dry deposition for NO3 in RAQ/S/T schemes.
chch_t1( 6,'NO3       ',1,'TR        ','NOx       ',0,1,0,R+S+T,0,0,107),&
chch_t1( 6,'NO3       ',1,'TR        ','NOx       ',1,1,0,TI+ST,0,0,107),&
chch_t1( 6,'NO3       ',1,'CF        ','NOx       ',0,0,0,OL,0,0,107),&
!   7 DD: 4,      EM: 1
chch_t1( 7,'NO2       ',1,'TR        ','NOx       ',1,0,1,TI+S+ST,0,0,107),&
chch_t1( 7,'NO2       ',1,'TR        ','NOx       ',1,0,0,T+R,0,0,107),&
! No dry deposition for N2O5 in R/S/T schemes
!   8 DD: 5,WD: 2,
chch_t1( 8,'N2O5      ',2,'TR        ','          ',0,1,0,R+S+T,0,0,107),&
chch_t1( 8,'N2O5      ',2,'TR        ','          ',1,1,0,TI+ST,0,0,107),&
!   9 DD: 6,WD: 3,
! No dry deposition for HO2NO2 in R and T schemes
chch_t1( 9,'HO2NO2    ',1,'TR        ','          ',0,1,0,R+T,0,0,107),&
chch_t1( 9,'HO2NO2    ',1,'TR        ','          ',1,1,0,TI+S+ST,A,0,107),&!A added, LER Apr2019
!  10 DD: 7,WD: 4,
chch_t1(10,'HONO2     ',1,'TR        ','          ',1,1,0,TI+S+T+ST+R,0,0,107),&
!  11 DD: 8,WD: 5,
chch_t1(11,'H2O2      ',1,'TR        ','          ',1,1,0,TI+S+T+ST+OL+R,A,0,&!A added, LER Apr2019
   107),&
!  12             EM: 2
! Dry deposition of CH4 not enabled in S and ST schemes.
chch_t1(12,'CH4       ',1,'TR        ','          ',0,0,2,S+ST,0,0,107),&
!  13 DD: 9,      EM:  3
chch_t1(12,'CH4       ',1,'TR        ','          ',1,0,1,TI+T+R,0,0,107),&
!  13 DD:10,      EM:  3
chch_t1(13,'CO        ',1,'TR        ','          ',1,0,1,TI+S+T+ST+R,0,0,107),&
!  14 DD:11,WD: 6,EM:  4
! No dry deposition for HCHO in S/T/R schemes
chch_t1(14,'HCHO      ',1,'TR        ','          ',0,1,1,S+T+R,0,0,107),&
chch_t1(14,'HCHO      ',1,'TR        ','          ',1,1,4,TI+ST,0,0,107),&
!  15       WD: 7,
chch_t1(15,'MeOO      ',1,'TR        ','          ',0,1,0,TI+S+ST,0,0,107),&
!  15       WD: 7,
chch_t1(15,'MeOO      ',1,'SS        ','          ',0,1,0,T+R,0,0,107),&
!  16 DD:12,WD: 8,
chch_t1(16,'MeOOH     ',1,'TR        ','          ',1,1,0,TI+S+T+ST+R,0,0,107),&
!  17
chch_t1(17,'H         ',1,'TR        ','HOx       ',0,0,0,S+ST,0,0,107),&
!  18
! No water in RAQ scheme? I presume H2O is implicit.
chch_t1(18,'H2O       ',1,'TR        ','          ',0,0,0,S+ST,0,0,107),&
!  18
chch_t1(18,'H2O       ',1,'CF        ','          ',0,0,0,TI+T+R,0,0,107),&
!  19
chch_t1(19,'OH        ',1,'TR        ','HOx       ',0,0,0,TI+S+ST,0,0,107),&
!  20       WD: 9,
chch_t1(19,'OH        ',1,'SS        ','HOx       ',0,0,0,T+R,0,0,107),&
!  20       WD: 9,
chch_t1(19,'OH        ',1,'CF        ','HOx       ',0,0,0,OL,0,0,107),&
!  20       WD: 9,
chch_t1(20,'HO2       ',1,'TR        ','HOx       ',0,1,0,TI+S+ST,0,0,107),&
!  20       WD: 9,
chch_t1(20,'HO2       ',1,'SS        ','HOx       ',0,1,0,T+R,0,0,107),&
!  20       WD: 9,
chch_t1(20,'HO2       ',1,'CF        ','HOx       ',0,0,0,OL,0,0,107),&
!  21
chch_t1(21,'Cl        ',1,'TR        ','Clx       ',0,0,0,S+ST,0,0,107),&
!  22
chch_t1(22,'Cl2O2     ',1,'TR        ','Clx       ',0,0,0,S+ST,0,0,107),&
!  23
chch_t1(23,'ClO       ',1,'TR        ','Clx       ',0,0,0,S+ST,0,0,107),&
!  24
chch_t1(24,'OClO      ',1,'TR        ','          ',0,0,0,S+ST,0,0,107),&
!  25
chch_t1(25,'Br        ',1,'TR        ','Brx       ',0,0,0,S+ST,0,0,107),&
!  26
chch_t1(26,'BrO       ',1,'TR        ','Brx       ',0,0,0,S+ST,0,0,107),&
!  27
chch_t1(27,'BrCl      ',1,'TR        ','          ',0,0,0,S+ST,0,0,107),&
!  28       WD:10,
chch_t1(28,'BrONO2    ',1,'TR        ','          ',0,1,0,S+ST,0,0,107),&
!  29
chch_t1(29,'N2O       ',1,'TR        ','          ',0,0,5,S+ST,0,0,107),&
!  30 DD:13,WD:11,
chch_t1(30,'HCl       ',1,'TR        ','          ',1,1,0,S+ST,0,0,107),&
!  31 DD:14,WD:12,
chch_t1(31,'HOCl      ',1,'TR        ','          ',1,1,0,S+ST,0,0,107),&
!  32 DD:15,WD:13,
chch_t1(32,'HBr       ',1,'TR        ','          ',1,1,0,S+ST,0,0,107),&
!  33 DD:15,WD:14,
chch_t1(33,'HOBr      ',1,'TR        ','          ',1,1,0,S+ST,A,0,107),&!A added, LER Apr2019
!  34       WD:15,
chch_t1(34,'ClONO2    ',1,'TR        ','          ',0,1,0,S+ST,0,0,107),&
!  35
chch_t1(35,'CFCl3     ',1,'TR        ','          ',0,0,6,S+ST,0,0,107),&
!  36
chch_t1(36,'CF2Cl2    ',1,'TR        ','          ',0,0,7,S+ST,0,0,107),&
!  37
chch_t1(37,'MeBr      ',1,'TR        ','          ',0,0,8,S+ST,0,0,107),&
!  38 DD:17,WD:16,
chch_t1(38,'HONO      ',1,'TR        ','          ',1,1,0,TI+T+ST,0,0,107),&
!  39             EM: 5
chch_t1(39,'C2H6      ',1,'TR        ','          ',0,0,8,TI+T+ST+R,0,0,107),&
!  40
chch_t1(40,'EtOO      ',1,'TR        ','          ',0,0,0,TI+ST,0,0,107),&
!  40
chch_t1(40,'EtOO      ',1,'SS        ','          ',0,0,0, T+R,0,0,107),&
!  41 DD:18,WD: 17,
chch_t1(41,'EtOOH     ',1,'TR        ','          ',1,1,0,TI+T+ST+R,0,0,107),&
!  42 DD:19,       EM: 6
! No dry deposition for MeCHO in R and T schemes 
chch_t1(42,'MeCHO     ',1,'TR        ','          ',0,0,1,R+T,0,0,107),&
chch_t1(42,'MeCHO     ',1,'TR        ','          ',1,0,1,TI+ST,0,0,107),&
!  43
chch_t1(43,'MeCO3     ',1,'TR        ','          ',0,0,0,TI+ST,0,0,107),&
!  43
chch_t1(43,'MeCO3     ',1,'SS        ','          ',0,0,0,T+R,0,0,107),&
!  44 DD:20,
chch_t1(44,'PAN       ',1,'TR        ','          ',1,0,0,TI+T+ST+R,0,0,107),&
!  45              EM: 7
chch_t1(45,'C3H8      ',1,'TR        ','          ',0,0,1,TI+T+ST+R,0,0,107),&
!  46
chch_t1(46,'n-PrOO    ',1,'TR        ','          ',0,0,0,TI+ST,0,0,107),&
!  46
chch_t1(46,'n-PrOO    ',1,'SS        ','          ',0,0,0,T,0,0,107),&
!  47
chch_t1(47,'i-PrOO    ',1,'TR        ','          ',0,0,0,TI+ST,0,0,107),&
!  47
chch_t1(47,'i-PrOO    ',1,'SS        ','          ',0,0,0,T+R,0,0,107),&
!  48 DD:21,WD:18,
chch_t1(48,'n-PrOOH   ',1,'TR        ','          ',1,1,0,TI+T+ST,0,0,107),&
!  49 DD:22,WD:19,
chch_t1(49,'i-PrOOH   ',1,'TR        ','          ',1,1,0,TI+T+ST+R,0,0,107),&
!  50 DD:23,
chch_t1(50,'EtCHO     ',1,'TR        ','          ',1,0,0,TI+T+ST,0,0,107),&
!  51
chch_t1(51,'EtCO3     ',1,'TR        ','          ',0,0,0,TI+ST,0,0,107),&
!  51
chch_t1(51,'EtCO3     ',1,'SS        ','          ',0,0,0,T,0,0,107),&
!  52             EM: 8
chch_t1(52,'Me2CO     ',1,'TR        ','          ',0,0,1,TI+T+ST+R,0,0,107),&
!  53
chch_t1(53,'MeCOCH2OO ',1,'TR        ','          ',0,0,0,TI+ST,0,0,107),&
!  53
chch_t1(53,'MeCOCH2OO ',1,'SS        ','          ',0,0,0,T+R,0,0,107),&
!  54 DD:24,WD:20,
chch_t1(54,'MeCOCH2OOH',1,'TR        ','          ',1,1,0,TI+T+ST,0,0,107),&
!  55 DD:25,
chch_t1(55,'PPAN      ',1,'TR        ','          ',1,0,0,TI+T+ST,0,0,107),&
!  56
chch_t1(56,'MeONO2    ',1,'TR        ','          ',0,0,0,TI+T+ST,0,0,107),&
!  57             EM: 9
chch_t1(57,'C5H8      ',1,'TR        ','          ',0,0,1,TI+ST+R,0,0,107),&
!  58
chch_t1(58,'ISO2      ',1,'TR        ','          ',0,0,0,TI+ST,0,0,107),&
!  59 DD:26,WD:21,
chch_t1(59,'ISOOH     ',1,'TR        ','          ',1,1,0,TI+ST+R,0,0,107),&
!  60 DD:27,WD:22,
! No DD for ISON in R scheme
chch_t1(60,'ISON      ',1,'TR        ','          ',0,1,0,R,0,0,107),&
chch_t1(60,'ISON      ',1,'TR        ','          ',1,1,0,TI+ST,0,0,107),&
!  61 DD:28,
chch_t1(61,'MACR      ',1,'TR        ','          ',1,0,0,TI+ST,0,0,107),&
!  62
chch_t1(62,'MACRO2    ',1,'TR        ','          ',0,0,0,TI+ST,0,0,107),&
!  63 DD:29,WD:23,
chch_t1(63,'MACROOH   ',1,'TR        ','          ',1,1,0,TI+ST,0,0,107),&
!  64 DD:30,
chch_t1(64,'MPAN      ',1,'TR        ','          ',1,0,0,TI+ST,0,0,107),&
!  65 DD:31,WD:24,
chch_t1(65,'HACET     ',1,'TR        ','          ',1,1,0,TI+ST,0,0,107),&
!  66 DD:32,WD:25,
! No dry deposition for methyl glyoxal in R scheme
chch_t1(66,'MGLY      ',1,'TR        ','          ',0,1,0,R,0,0,107),&
chch_t1(66,'MGLY      ',1,'TR        ','          ',1,1,0,TI+ST,0,0,107),&
!  67 DD:33,
chch_t1(67,'NALD      ',1,'TR        ','          ',1,0,0,TI+ST,0,0,107),&
!  68 DD:34,WD:26,
chch_t1(68,'HCOOH     ',1,'TR        ','          ',1,1,0,TI+ST,0,0,107),&
!  69 DD:35,WD:27,
chch_t1(69,'MeCO3H    ',1,'TR        ','          ',1,1,0,TI+ST,0,0,107),&
!  70 DD:36,WD:28,
chch_t1(70,'MeCO2H    ',1,'TR        ','          ',1,1,0,TI+ST,0,0,107),&
!  71 DD:37,
! Dry deposition for H2 in R scheme but not in other schemes. 
! H2 is prescribed at the surface. Could it become interactive? 
chch_t1(71,'H2        ',1,'TR        ','          ',0,0,0,S+ST,0,0,107),&
chch_t1(71,'H2        ',1,'TR        ','          ',1,0,1,R,0,0,107),  &
!  71
chch_t1(71,'H2        ',1,'CT        ','          ',0,0,0,TI+T,0,0,107),&
!  72 DD:38,WD:29,
! No dry deposition for MeOH in R scheme. No emissions of MeOH in other schemes.
chch_t1(72,'MeOH      ',1,'TR        ','          ',1,1,0,TI+ST,0,0,107),&
chch_t1(72,'MeOH      ',1,'TR        ','          ',0,1,1,R,0,0,107),  &
!  73
chch_t1(73,'CO2       ',1,'CT        ','          ',0,0,0,TI+T+S+ST,0,CI,107),&
!  73
! Stratospheric chemical schemes are enabled to take 3D CO2 as input
chch_t1(73,'CO2       ',1,'CF        ','          ',0,0,0,S+ST,CI,0,111),&
!  74
chch_t1(74,'O2        ',1,'CT        ','          ',0,0,0,TI+S+T+ST,0,0,107),&
!  75
chch_t1(75,'N2        ',1,'CT        ','          ',0,0,0,TI+S+T+ST,0,0,107),&
chch_t1(76,'DMS       ',1,'TR        ','          ',0,1,1,TI+S+ST+OL+R,A,0,&!Made soluble, LER, Mar2019
   107), &
!  77 DD:39,WD:31,EM:10
chch_t1(77,'SO2       ',1,'TR        ','          ',1,1,1,TI+S+ST+OL+R,A,0,&
                                                                     107), &
!  78 DD:40,
chch_t1(78,'H2SO4     ',1,'TR        ','          ',1,0,0,TI+S+ST+OL+R,A,0,&
                                                                     107), &
!  79 DD:41,
chch_t1(79,'MSA       ',1,'TR        ','          ',0,0,0,S+ST+R,A,0,107),&
chch_t1(79,'MSA       ',1,'TR        ','          ',1,0,0,TI,A,0,107),&
!  80 DD:42,WD:32
! No dry deposition for DMSO in R scheme.
chch_t1(80,'DMSO      ',1,'TR        ','          ',0,1,0,R,A,0,107),&
chch_t1(80,'DMSO      ',1,'TR        ','          ',1,1,0,TI+ST+OL,A,0,107),&
!  81 DD:43,WD:33,EM:11
chch_t1(81,'NH3       ',1,'TR        ','          ',1,1,1,TI+ST+R,A,0,107),&
!  82
chch_t1(82,'CS2       ',1,'TR        ','          ',0,0,1,TI+S+ST,A,0,107),&
!  83
chch_t1(83,'COS       ',1,'TR        ','          ',0,0,1,TI+S+ST,A,0,107),&
!  84
chch_t1(84,'H2S       ',1,'TR        ','          ',0,0,1,TI+S+ST,A,0,107),&
!  85
chch_t1(85,'SO3       ',1,'TR        ','          ',0,0,0,S+ST,A,0,107),&
!  86 DD:44,      EM: 12
chch_t1(86,'Monoterp  ',1,'TR        ','          ',1,0,1,TI+ST+OL+R,A,0,107),&
!  87 DD:45,WD:34
chch_t1(87,'Sec_Org   ',1,'TR        ','          ',1,1,0,TI+ST+OL+R,A,0,107),&
chch_t1(88,'O(3P)S    ',1,'SS        ','          ',0,0,0,T,0,0,107),&
chch_t1(89,'O(1D)S    ',1,'SS        ','          ',0,0,0,T,0,0,107),&
! O3S does not have wet deposition but O3 does for A chemistry!
!  88 DD:46,
chch_t1(90,'O3S       ',1,'TR        ','          ',1,0,0,T+R,0,0,107),&
chch_t1(91,'OHS       ',1,'SS        ','          ',0,0,0,T,0,0,107),&
!  89       WD:35
chch_t1(92,'HO2S      ',1,'SS        ','          ',0,1,0,T,0,0,107),&
chch_t1(93,'s-BuOO    ',1,'SS        ','          ',0,0,0,R,0,0,107),&
chch_t1(94,'MVK       ',1,'TR        ','          ',0,0,0,R,0,0,107),&
!  90 DD:47,WD:36
chch_t1(95,'MVKOOH    ',1,'TR        ','          ',1,1,0,R,0,0,107),&
chch_t1(96,'MEKO2     ',1,'SS        ','          ',0,0,0,R,0,0,107),&
chch_t1(97,'HOC2H4O2  ',1,'SS        ','          ',0,0,0,R,0,0,107),&
!  91 DD:48,WD:37
chch_t1(98,'ORGNIT    ',1,'TR        ','          ',1,1,0,R,0,0,107),&
chch_t1(99,'HOC3H6O2  ',1,'SS        ','          ',0,0,0,R,0,0,107),&
chch_t1(100,'OXYL1     ',1,'SS        ','          ',0,0,0,R,0,0,107),&
chch_t1(101,'MEMALD1   ',1,'SS        ','          ',0,0,0,R,0,0,107),&
chch_t1(102,'RNC2H4    ',1,'TR        ','          ',0,0,0,R,0,0,107),&
chch_t1(103,'HOIPO2    ',1,'SS        ','          ',0,0,0,R,0,0,107),&
chch_t1(104,'RNC3H6    ',1,'TR        ','          ',0,0,0,R,0,0,107),&
chch_t1(105,'HOMVKO2   ',1,'SS        ','          ',0,0,0,R,0,0,107),&
!  92             EM: 13
chch_t1(106,'C2H4      ',1,'TR        ','          ',0,0,1,R,0,0,107),&
!  93             EM: 14
chch_t1(107,'C3H6      ',1,'TR        ','          ',0,0,1,R,0,0,107),&
!  94             EM: 15
chch_t1(108,'C4H10     ',1,'TR        ','          ',0,0,1,R,0,0,107),&
!  95 DD:49,WD:38
chch_t1(109,'s-BuOOH   ',1,'TR        ','          ',1,1,0,R,0,0,107),&
chch_t1(110,'MEK       ',1,'TR        ','          ',0,0,0,R,0,0,107),&
!  96             EM: 16
chch_t1(111,'TOLUENE   ',1,'TR        ','          ',0,0,1,R,0,0,107),&
chch_t1(112,'TOLP1     ',1,'SS        ','          ',0,0,0,R,0,0,107),&
chch_t1(113,'MEMALD    ',1,'TR        ','          ',0,0,0,R,0,0,107),&
!  97       WD:39
chch_t1(114,'GLY       ',1,'TR        ','          ',0,1,0,R,0,0,107),&
chch_t1(115,'oXYLENE   ',1,'TR        ','          ',0,0,1,R,0,0,107),&
chch_t1(116,'MSIA      ',1,'TR        ','          ',0,1,0,ST,A,0,107)/) !Made soluble, LER, Mar2019

! Heterogeneous chemistry
! Columns take the following meanings:
! Item number, reactant1, reactant2, product1, product2, product3, product4,
! product yield x 4, chemistry scheme, qualifier, disqualifier, version
TYPE(rath_t1), PARAMETER :: rath_defs_master(1:n_het_master)=(/         &
rath_t1(1,'ClONO2    ','H2O       ','HOCl      ','HONO2     ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000, S+ST,HP,0,107), &
rath_t1(2,'ClONO2    ','HCl       ','Cl        ','Cl        ','HONO2     ', &
'          ', 0.000, 0.000, 0.000, 0.000, S+ST,HP,0,107), &
rath_t1(3,'HOCl      ','HCl       ','Cl        ','Cl        ','H2O       ', &
'          ', 0.000, 0.000, 0.000, 0.000, S+ST,HP,0,107), &
rath_t1(4,'N2O5      ','H2O       ','HONO2     ','HONO2     ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000, S+ST,HP,0,107), &
rath_t1(5,'N2O5      ','HCl       ','Cl        ','NO2       ','HONO2     ', &
'          ', 0.000, 0.000, 0.000, 0.000, S+ST,HP,0,107), &

! Aerosol chemistry: there are no gas phase products, the 'NULLx' products
!  identify the reactions in asad_hetero

! NULL0-NULL2 are existing reactions that also form part of the Chen et al. scheme. LER 28.03.19
!HSO3+H2O2(aq)
rath_t1(6,'SO2       ','H2O2      ','NULL0     ','          ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000, TI+S+ST+OL+R,A,0,107),             &
!HSO3+O3(aq)
rath_t1(7,'SO2       ','O3        ','NULL1     ','          ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000, TI+S+ST+OL+R,A,0,107),             &
!HSO3+O3(aq)
rath_t1(8,'SO2       ','O3        ','NULL2     ','          ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000, TI+S+ST+OL+R,A,0,107),             &

! Chen et al (2018) aqueous phase chemistry, LER 28.03.19
!DMS(aq)+O3(aq)
rath_t1(9,'DMS       ','O3        ','NULL3     ','          ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000, ST,A,0,107), &
!MSIA(aq)+O3(aq)
rath_t1(10,'MSIA      ','O3       ','NULL4     ','          ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000, ST,A,0,107), &
!MSI-+O3(aq)
rath_t1(11,'MSIA      ','O3       ','NULL5     ','          ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000, ST,A,0,107), &
!HSO3-+HOBr(aq)
rath_t1(12,'SO2       ','HOBr     ','NULL6     ','          ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000, ST,A,0,107), &
!SO3--+HOBr(aq)
rath_t1(13,'SO2       ','HOBr     ','NULL7     ','          ','          ', &
'          ', 0.000, 0.000, 0.000, 0.000, ST,A,0,107), &

! Tropospheric heterogenous reactions
rath_t1(14,'N2O5      ','          ','HONO2     ','          ','          ', &
'          ', 2.000, 0.000, 0.000, 0.000, TI+R,TH,0,107),               &
! Heterogenous
rath_t1(15,'HO2       ','          ','H2O2      ','          ','          ', &
'          ', 0.500, 0.000, 0.000, 0.000, TI,TH,0,107)  /)


! photolysis reactions
! Item number, reactant, 'PHOTON', product x4,
! fractional production coefficient x4, quantum yield (%),
! photolysis rate label, chemistry scheme, qualifier, disqualifier, version 
TYPE(ratj_t1), PARAMETER :: ratj_defs_master(1:n_ratj_master)=(/&
! 1
ratj_t1(1,'EtOOH     ','PHOTON    ','MeCHO     ','HO2       ','OH        ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',TI+T+ST+R,0,0,107) ,&
! 2
ratj_t1(2,'H2O2      ','PHOTON    ','OH        ','OH        ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jh2o2     ',TI+S+T+ST+R,0,0,107) ,&
! 3
! This should produce H+ CHO -> H + HO2 + CO in ST scheme.
ratj_t1(3,'HCHO      ','PHOTON    ','HO2       ','HO2       ','CO        ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhchoa    ',TI+T+ST+R,0,0,107) ,&
ratj_t1(3,'HCHO      ','PHOTON    ','H         ','CO        ','HO2       ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhchoa    ',S,0,0,107) ,&
! 4
ratj_t1(4,'HCHO      ','PHOTON    ','H2        ','CO        ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhchob    ',TI+S+T+ST+R,0,0,107) ,&
! 5
! In S and ST, branching is assumed following JPL 2011. Update T, TI, and R?
ratj_t1(5,'HO2NO2    ','PHOTON    ','HO2       ','NO2       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jpna67    ',ST,0,0,107) ,&
ratj_t1(5,'HO2NO2    ','PHOTON    ','NO2       ','HO2       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jpna67    ',S,0,0,107) ,&
ratj_t1(5,'HO2NO2    ','PHOTON    ','HO2       ','NO2       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jpna      ',TI+T+R,0,0,107) ,&
! 6
ratj_t1(6,'HONO2     ','PHOTON    ','OH        ','NO2       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhono2    ',TI+T+ST+R,0,0,107) ,&
ratj_t1(6,'HONO2     ','PHOTON    ','NO2       ','OH        ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhono2    ',S,0,0,107) ,&
! 7
ratj_t1(7,'MeCHO     ','PHOTON    ','MeOO      ','HO2       ','CO        ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jaceta    ',TI+T+ST+R,0,0,107) ,&
! 8
! Missing in R.
ratj_t1(8,'MeCHO     ','PHOTON    ','CH4       ','CO        ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jacetb    ',TI+T+ST,0,0,107) ,&
! 9
! THIS REACTION IS MISSING IN THE TROPOSPHERIC CHEMISTRY SCHEME.
ratj_t1(9,'MeOOH     ','PHOTON    ','HO2       ','HCHO      ','OH        ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',TI+ST+R,0,0,107) ,&
ratj_t1(9,'MeOOH     ','PHOTON    ','HCHO      ','HO2       ','OH        ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',S,0,0,107) ,&
! 10
ratj_t1(10,'N2O5      ','PHOTON    ','NO3       ','NO2       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jn2o5     ',TI+T+ST+R,0,0,107),&
ratj_t1(10,'N2O5      ','PHOTON    ','NO2       ','NO3       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jn2o5     ',S,0,0,107) ,&
! 11
ratj_t1(11,'NO2       ','PHOTON    ','NO        ','O(3P)     ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jno2      ',TI+S+T+ST+R,0,0,107) ,&
! 12
ratj_t1(12,'NO3       ','PHOTON    ','NO        ','O2        ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jno3a     ',TI+S+T+ST+R,0,0,107) ,&
! 13
ratj_t1(13,'NO3       ','PHOTON    ','NO2       ','O(3P)     ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jno3b     ',TI+S+T+ST+R,0,0,107) ,&
! 14
! Missing in R
ratj_t1(14,'O2        ','PHOTON    ','O(3P)     ','O(3P)     ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jo2       ',TI+S+T+ST,0,0,107) ,&
! 15
ratj_t1(15,'O3        ','PHOTON    ','O2        ','O(1D)     ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jo3a      ',TI+S+T+ST+R,0,0,107) ,&
! 16
ratj_t1(16,'O3        ','PHOTON    ','O2        ','O(3P)     ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jo3b      ',TI+S+T+ST+R,0,0,107) ,&
! 17
ratj_t1(17,'PAN       ','PHOTON    ','MeCO3     ','NO2       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jpan      ',TI+T+ST+R,0,0,107) ,&
! 18
! HONO not modelled in R
ratj_t1(18,'HONO      ','PHOTON    ','OH        ','NO        ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhono     ',TI+T+ST,0,0,107) ,&
! 19
ratj_t1(19,'EtCHO     ','PHOTON    ','EtOO      ','HO2       ','CO        ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jetcho    ',TI+T+ST,0,0,107) ,&
! 20
ratj_t1(20,'Me2CO     ','PHOTON    ','MeCO3     ','MeOO      ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jaceto    ',TI+ST+R,0,0,107) ,&
! The naming of 'jacetone' is probably inconsistent with FAST-JX. Error?
ratj_t1(20,'Me2CO     ','PHOTON    ','MeCO3     ','MeOO      ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jacetone  ',T,0,0,107) ,&
! 21 missing in A
ratj_t1(21,'n-PrOOH   ','PHOTON    ','EtCHO     ','HO2       ','OH        ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',T+ST,0,0,107) ,&
ratj_t1(21,'n-PrOOH   ','PHOTON    ','HO2       ','EtCHO     ','OH        ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',TI,0,0,107) ,&
! 22 missing in A
ratj_t1(22,'i-PrOOH   ','PHOTON    ','Me2CO     ','HO2       ','OH        ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',T+ST+R,0,0,107) ,&
ratj_t1(22,'i-PrOOH   ','PHOTON    ','HO2       ','Me2CO     ','OH        ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',TI,0,0,107) ,&
! 23
! reaction missing in tropospheric chemistry scheme
ratj_t1(23,'MeCOCH2OOH','PHOTON    ','MeCO3     ','HCHO      ','OH        ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',ST,0,0,107) ,&
ratj_t1(23,'MeCOCH2OOH','PHOTON    ','HCHO      ','MeCO3      ','OH        ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',TI,0,0,107) ,&
! 24
ratj_t1(24,'PPAN      ','PHOTON    ','EtCO3     ','NO2       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jpan      ',TI+T+ST,0,0,107) ,&
! 25
ratj_t1(25,'MeONO2    ','PHOTON    ','HO2       ','HCHO      ','NO2       ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmena     ',TI+T+ST,0,0,107) ,&
! 26 Different products in TI and ST versus R.
ratj_t1(26,'ISOOH     ','PHOTON    ','OH        ','MACR      ','HCHO      ',&
     'HO2       ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',TI+ST,0,0,107) ,&
ratj_t1(26,'ISOOH     ','PHOTON    ','OH        ','MVK      ','HCHO      ',&
     'HO2       ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',R,0,0,107) ,&
! 27
ratj_t1(27,'ISON      ','PHOTON    ','NO2       ','MACR      ','HCHO      ',&
     'HO2       ', 0.0,0.0,0.0,0.0, 100.000,'jiprn     ',TI+ST,0,0,107) ,&
! 28
ratj_t1(28,'MACR      ','PHOTON    ','MeCO3     ','HCHO      ','CO        ',&
     'HO2       ', 0.0,0.0,0.0,0.0, 100.000,'jmacr     ',TI+ST,0,0,107) ,&
! 29
ratj_t1(29,'MPAN      ','PHOTON    ','MACRO2    ','NO2       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jpan      ',TI+ST,0,0,107) ,&
! 30
ratj_t1(30,'MACROOH   ','PHOTON    ','OH        ','HO2       ','OH        ',&
     'HO2       ', 0.0,0.0,0.0,0.0, 100.000,'jmacro    ',TI+ST,0,0,107) ,&
! 31
ratj_t1(31,'MACROOH   ','PHOTON    ','HACET     ','CO        ','MGLY      ',&
     'HCHO      ', 0.0,0.0,0.0,0.0, 100.000,'jmacro    ',TI+ST,0,0,107) ,&
! 32
ratj_t1(32,'HACET     ','PHOTON    ','MeCO3     ','HCHO      ','HO2       ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhacet    ',TI+ST,0,0,107) ,&
! 33
ratj_t1(33,'MGLY      ','PHOTON    ','MeCO3     ','CO        ','HO2       ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmkal     ',TI+ST+R,0,0,107) ,&
! 34
ratj_t1(34,'NALD      ','PHOTON    ','HCHO      ','CO        ','NO2       ',&
     'HO2       ', 0.0,0.0,0.0,0.0, 100.000,'jaceta    ',TI+ST,0,0,107) ,&
! 35
ratj_t1(35,'MeCO3H    ','PHOTON    ','MeOO      ','OH        ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmeco3h   ',TI+ST,0,0,107) ,&
! 36
ratj_t1(36,'BrCl      ','PHOTON    ','Br        ','Cl        ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jbrcl     ',S+ST,0,0,107) ,&
! 37
ratj_t1(37,'BrO       ','PHOTON    ','Br        ','O(3P)     ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jbro      ',S+ST,0,0,107) ,&
! 38
ratj_t1(38,'BrONO2    ','PHOTON    ','Br        ','NO3       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jbrna     ',S+ST,0,0,107) ,&
! 39
ratj_t1(39,'BrONO2    ','PHOTON    ','BrO       ','NO2       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jbrnb     ',S+ST,0,0,107) ,&
! 40
ratj_t1(40,'O2        ','PHOTON    ','O(3P)     ','O(1D)     ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jo2b      ',S+ST,0,0,107) ,&
! 41
ratj_t1(41,'OClO      ','PHOTON    ','O(3P)     ','ClO       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'joclo     ',S+ST,0,0,107) ,&
! 42
ratj_t1(42,'NO        ','PHOTON    ','N         ','O(3P)     ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jno       ',S+ST,0,0,107) ,&
! 43
ratj_t1(43,'HOBr      ','PHOTON    ','OH        ','Br        ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhobr     ',S+ST,0,0,107) ,&
! 44
ratj_t1(44,'N2O       ','PHOTON    ','N2        ','O(1D)     ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jn2o      ',S+ST,0,0,107) ,&
! 45
ratj_t1(45,'H2O       ','PHOTON    ','OH        ','H         ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jh2o      ',S+ST,0,0,107) ,&
! 46
ratj_t1(46,'ClONO2    ','PHOTON    ','Cl        ','NO3       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jclna     ',S+ST,0,0,107) ,&
! 47
ratj_t1(47,'ClONO2    ','PHOTON    ','ClO       ','NO2       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jclnb     ',S+ST,0,0,107) ,&
! 48
ratj_t1(48,'HCl       ','PHOTON    ','H         ','Cl        ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhcl      ',S+ST,0,0,107) ,&
! 49
ratj_t1(49,'HOCl      ','PHOTON    ','OH        ','Cl        ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jhocl     ',S+ST,0,0,107) ,&
! 50
ratj_t1(50,'Cl2O2     ','PHOTON    ','Cl        ','Cl        ','O2        ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jcl2o2    ',S+ST,0,0,107) ,&
! 51
ratj_t1(51,'CFCl3     ','PHOTON    ','Cl        ','Cl        ','Cl        ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jcfcl3    ',S+ST,0,0,107) ,&
! 52
ratj_t1(52,'CF2Cl2    ','PHOTON    ','Cl        ','Cl        ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jcfc2     ',S+ST,0,0,107) ,&
! 53
ratj_t1(53,'MeBr      ','PHOTON    ','Br        ','H         ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmebr     ',S+ST,0,0,107) ,&
! 54
ratj_t1(54,'CH4       ','PHOTON    ','MeOO      ','H         ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jch4      ',S+ST,0,0,107) ,&
! 55
ratj_t1(55,'CO2       ','PHOTON    ','CO        ','O(3P)     ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jco2      ',S+ST,0,0,107),&
! 56
! Reaction missing in the R scheme.
ratj_t1(56,'HO2NO2    ','PHOTON    ','OH        ','NO3       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jpna33    ',ST,0,0,107),&
ratj_t1(56,'HO2NO2    ','PHOTON    ','NO3       ','OH       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jpna33    ',S,0,0,107),&
! 57
ratj_t1(57,'CS2       ','PHOTON    ','COS        ','SO2       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jcs2      ',S+ST,A,0,107), &
! 58
ratj_t1(58,'COS       ','PHOTON    ','CO         ','SO2       ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jcos      ',S+ST,A,0,107), &
! 59
ratj_t1(59,'H2SO4     ','PHOTON    ','SO3        ','OH        ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jh2so4    ',S+ST,A,0,107), &
! 60
ratj_t1(60,'SO3       ','PHOTON    ','SO2        ','O(3P)     ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jso3      ',S+ST,A,0,107),&
ratj_t1(61,'MEK       ','PHOTON    ','MeCO3     ','EtOO      ','          ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jaceto    ',R,0,0,107),&
ratj_t1(62,'GLY       ','PHOTON    ','CO        ','CO        ','HO2       ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmkal     ',R,0,0,107),&
ratj_t1(63,'s-BuOOH   ','PHOTON    ','MEK       ','HO2       ','OH        ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',R,0,0,107),&
ratj_t1(64,'MVKOOH    ','PHOTON    ','MGLY      ','HCHO      ','OH        ',&
     '          ', 0.0,0.0,0.0,0.0, 100.000,'jmhp      ',R,0,0,107) /)

! Termolecular reactions
! Rates taken from  file "UKCA_Reaction_Rates - Termolecular Reactions.csv"
!on 11:41 2012/05/01

! Item number, reactant x2, product x2, broadening coefficient parameter F,
! low-pressure limit coefficients (k, alpha and beta), infinite-pressure limit
! coefficients (k, alpha and beta), product yields x 2, chemistry scheme,
! qualifier, disqualifier, version 
TYPE(ratt_t1), PARAMETER :: ratt_defs_master(1:n_ratt_master)=(/   &
ratt_t1(1,'O(3P)     ','O2        ','O3        ','m         ',     0.0,&
  5.70e-34, -2.60, 0.00,  0.00e+00,  0.0, 0.0, 0.0, 0.0, TI+T,0,0,107),&
! T001 JPL 2011
ratt_t1(1,'O(3P)     ','O2        ','O3        ','m         ',     0.0,&
  6.00e-34, -2.50, 0.00,  0.00e+00,  0.0, 0.0, 0.0, 0.0, ST+R,0,0,107),&
! C
ratt_t1(1,'O(3P)     ','O2        ','O3        ','m         ',     0.0,&
  5.60e-34, -2.60, 0.00,  0.00e+00,  0.0, 0.0, 0.0, 0.0, S,0,0,107) ,& 
! T002 JPL 2011
! Missing in T and TI schemes. Note the discrepancy in low-temp rate
! coefficient
ratt_t1(2,'O(3P)     ','NO        ','NO2       ','m         ',     0.6,&
  9.00e-32, -1.50, 0.00,  3.00e-11,  0.0, 0.0, 0.0, 0.0, ST+R,0,0,107),&
ratt_t1(2,'O(3P)     ','NO        ','NO2       ','m         ',     0.6,&
  9.00e-31, -1.50, 0.00,  3.00e-11,  0.0, 0.0, 0.0, 0.0, S,0,0,107),&
! T003 JPL 2011, Missing in T and TI scheme
ratt_t1(3,'O(3P)     ','NO2       ','NO3       ','m         ',     0.6,&
  2.50e-31, -1.80, 0.00,  2.20e-11, -0.7, 0.0, 0.0, 0.0, ST,0,0,107),&
ratt_t1(3,'O(3P)     ','NO2       ','NO3       ','m         ',     0.6,&
  1.30e-31, -1.50, 0.00,  2.30e-11, 0.24, 0.0, 0.0, 0.0, S,0,0,107) ,& 
! T004 JPL 2011
ratt_t1(4,'O(1D)     ','N2        ','N2O       ','m         ',     0.0,&
  2.80e-36, -0.90, 0.00,  0.00e+00,  0.0, 0.0, 0.0, 0.0, ST,0,0,107),&
! JPL 2003
ratt_t1(4,'O(1D)     ','N2        ','N2O       ','m         ',     0.0,&
  3.50e-37, -0.60, 0.00,  0.00e+00,  0.0, 0.0, 0.0, 0.0, S,0,0,107) ,& 
! T005 JPL 2011
ratt_t1(5,'BrO       ','NO2       ','BrONO2    ','m         ',     0.6,&
  5.20e-31, -3.20, 0.00,  6.90e-12,  0.0, 0.0, 0.0, 0.0, ST,0,0,107),&
ratt_t1(5,'BrO       ','NO2       ','BrONO2    ','m         ',   327.0,&
  4.70e-31, -3.10, 0.00,  1.40e-11, -1.2, 0.0, 0.0, 0.0, S,0,0,107) ,& 
! T006 JPL 2011
ratt_t1(6,'ClO       ','ClO       ','Cl2O2     ','m         ',     0.6,&
  1.60e-32, -4.50, 0.00,  3.00e-12, -2.0, 0.0, 0.0, 0.0, ST,0,0,107),&
ratt_t1(6,'ClO       ','ClO       ','Cl2O2     ','m         ',     0.6,&
  1.70e-32, -4.00, 0.00,  5.40e-12,  0.0, 0.0, 0.0, 0.0, S,0,0,107) ,& 
! T007 IUPAC 2005
ratt_t1(7,'Cl2O2     ','m         ','ClO       ','ClO       ',     0.5,&
  3.70e-07,  0.00, 7690., 1.80e+14,  0.0,7690., 0.0, 0.0, ST,0,0,107),&
ratt_t1(7,'Cl2O2     ','m         ','ClO       ','ClO       ',     0.6,&
  1.00e-06,  1.00, 8000., 4.80e+15,  0.00,8820., 0.0, 0.0, S,0,0,107) ,& 
! T008 JPL 2011
ratt_t1(8,'ClO       ','NO2       ','ClONO2    ','m         ',     0.6,&
  1.80e-31, -3.40, 0.00,  1.50e-11,  0.0, 0.0, 0.0, 0.0, ST,0,0,107),&
ratt_t1(8,'ClO       ','NO2       ','ClONO2    ','m         ',   430.0,&
  1.60e-31, -3.40, 0.00,  1.50e-11,  0.0, 0.0, 0.0, 0.0, S,0,0,107) ,& 
! T009 JPL 2011
ratt_t1(9,'H         ','O2        ','HO2       ','m         ',     0.6,&
  4.40e-32, -1.30, 0.00,  7.50e-11,  0.0, 0.0, 0.0, 0.0, ST,0,0,107),&
! JPL2003
ratt_t1(9,'H         ','O2        ','HO2       ','m         ',     0.6,&
  5.70e-32, -1.60, 0.00,  7.50e-11,  0.0, 0.0, 0.0, 0.0, S,0,0,107),& 
! T010 JPL 2011 see also asad_trimol.F90
ratt_t1(10,'HO2       ','HO2       ','H2O2      ','O2        ',     0.0,&
  2.10e-33,  0.00,-920.,  0.00e+00,  0.0, 0.0, 0.0, 0.0, ST+R,0,0,107),&
ratt_t1(10,'HO2       ','HO2       ','H2O2      ','O2        ',     0.0,&
  1.90e-33,  0.00,-980.,  0.00e+00,  0.0, 0.0, 0.0, 0.0, TI+S+T,0,0,107),&
! T011 JPL 2011
ratt_t1(11,'HO2       ','NO2       ','HO2NO2    ','m         ',     0.6,&
  2.00e-31, -3.40, 0.00,  2.90e-12,  0.0, 0.0, 0.0, 0.0, ST+R,0,0,107),&
ratt_t1(11,'HO2       ','NO2       ','HO2NO2    ','m         ',     0.6,&
  1.80e-31, -3.20, 0.00,  4.70e-12,  0.0, 0.0, 0.0, 0.0, TI+S+T,0,0,107),&
! T012 IUPAC 2001
ratt_t1(12,'HO2NO2    ','m         ','HO2       ','NO2       ',     0.5,&
  4.10e-05,  0.00,10650., 4.80e+15,  0.0,11170.,0.0, 0.0, T+ST+R,0,0,107),&
ratt_t1(12,'HO2NO2    ','m         ','HO2       ','NO2       ',     0.6,&
  4.10e-05,  0.00,10650., 4.80e+15,  0.0,11170.,0.0, 0.0, TI,0,0,107),&
! Note the different rate
ratt_t1(12,'HO2NO2    ','m         ','HO2       ','NO2       ',     0.6,&
  4.10e-06,  0.00,10650., 4.80e+15,  0.0,11170.,0.0, 0.0, S,0,0,107) ,& 
ratt_t1(13,'OH        ','NO        ','HONO      ','m         ',  1420.0,&
  7.40e-31, -2.40, 0.00,  3.30e-11, -0.3, 0.0, 0.0, 0.0, TI+T,0,0,107),&
! T013 JPL 2011
ratt_t1(13,'OH        ','NO        ','HONO      ','m         ',     0.6,&
  7.00e-31, -2.60, 0.00,  3.60e-11, -0.10, 0.0, 0.0, 0.0, ST,0,0,107),&
! T014 JPL 2011
ratt_t1(14,'OH        ','NO2       ','HONO2     ','m         ',     0.6,&
  1.80e-30, -3.00, 0.00,  2.80e-11,  0.0, 0.0, 0.0, 0.0, ST+R,0,0,107),&
ratt_t1(14,'OH        ','NO2       ','HONO2     ','m         ',     0.4,&
  3.30e-30, -3.00, 0.00,  4.10e-11,  0.0, 0.0, 0.0, 0.0, TI+S+T,0,0,107),&
! T015 JPL 2011
ratt_t1(15,'OH        ','OH        ','H2O2      ','m         ',     0.6,&
  6.90e-31, -1.00, 0.00,  2.60e-11,  0.0, 0.0, 0.0, 0.0, ST,0,0,107),&
ratt_t1(15,'OH        ','OH        ','H2O2      ','m         ',     0.5,&
  6.90e-31, -0.80, 0.00,  2.60e-11,  0.0, 0.0, 0.0, 0.0, TI+S+T,0,0,107),&
! T016 MCMv3.2
ratt_t1(16,'MeCO3     ','NO2       ','PAN       ','m         ',     0.3,&
  2.70e-28, -7.10, 0.00,  1.20e-11, -0.9, 0.0, 0.0, 0.0, TI+T+ST+R,0,0,107),&
! T017 MCMv3.2
ratt_t1(17,'PAN       ','m         ','MeCO3     ','NO2       ',     0.3,&
  4.90e-03,  0.00,12100., 5.40e+16,  0.0,13830.,0.0, 0.0, TI+T+ST+R,0,0,107),&
! T018 MCMv3.2
ratt_t1(18,'EtCO3     ','NO2       ','PPAN      ','m         ',     0.3,&
  2.70e-28, -7.10, 0.00,  1.20e-11, -0.9, 0.0, 0.0, 0.0, TI+T+ST,0,0,107),&
ratt_t1(19,'PPAN      ','m         ','EtCO3     ','NO2       ',     0.4,&
  1.70e-03,  0.00,11280., 8.30e+16,  0.0,13940.,0.0, 0.0, TI+T,0,0,107),&
! T019 MCMv3.2
ratt_t1(19,'PPAN      ','m         ','EtCO3     ','NO2       ',     0.3,&
  4.90e-03,  0.00,12100., 5.40e+16,  0.0,13830.,0.0, 0.0, ST,0,0,107),&
! T020 MVMv3.2
ratt_t1(20,'MACRO2    ','NO2       ','MPAN      ','m         ',     0.3,&
  2.70e-28, -7.10, 0.00,  1.20e-11, -0.9, 0.0, 0.0, 0.0, TI+ST,0,0,107),&
! T021 MCMv3.2
ratt_t1(21,'MPAN      ','m         ','MACRO2    ','NO2       ',     0.3,&
  4.90e-03,  0.00,12100., 5.40e+16,  0.0,13830.,0.0, 0.0, TI+ST,0,0,107),&
! T022 IUPAC 2002
! Missing in S schemes
ratt_t1(22,'NO2       ','NO3       ','N2O5      ','m         ',     0.3,&
  3.60e-30, -4.10, 0.00,  1.90e-12,  0.2, 0.0, 0.0, 0.0, TI+T+S+ST+R,0,0,107),&
! T023 IUPAC 2002
ratt_t1(23,'N2O5      ','m         ','NO2       ','NO3       ',     0.3,&
  1.30e-03, -3.50,11000., 9.70e+14,  0.1,11080.,0.0, 0.0, TI+S+T+ST+R,0,0,107),&
! T024 IUPAC 2001
ratt_t1(24,'NO        ','NO        ','NO2       ','NO2       ',     0.0,&
  3.30e-39,  0.00,-530.,  0.00e+00,  0.0, 0.0, 0.0, 0.0, T+ST,0,0,107),&
! B; not in TI/TI scheme
ratt_t1(24,'NO        ','NO        ','NO2       ','NO2       ',     0.0,&
  6.93e-40,  0.00,-530.,  0.00e+00,  0.0, 0.0, 0.0, 0.0, S,0,0,107),&
! ! Original reaction:
!ratt_t1(25,'SO2       ','OH        ','SO3       ','HO2       ',     0.6,&
!  3.00e-31, -3.30, 0.00,  1.50e-12,  0.0, 0.0, 0.0, 0.0, ST+S,A,0,107),&
!Replace with SO2+OH reaction from Chen et al (2018), LER, UC, Feb 2019:
ratt_t1(25,'SO2       ','OH        ','H2SO4     ','HO2       ',     0.6,&
  3.30e-31, -4.30, 0.00,  1.60e-12,  0.0, 0.0, 0.0, 0.0, ST+S,A,0,107),&
ratt_t1(25,'SO2       ','OH        ','HO2       ','H2SO4     ',     0.6,&
  3.00e-31, -3.30, 0.00,  1.50e-12,  0.0, 0.0, 0.0, 0.0, TI,A,0,107),&
ratt_t1(25,'SO2       ','OH        ','H2SO4     ','          ',     0.6,&
  3.00e-31, -3.30, 0.00,  1.50e-12,  0.0, 0.0, 0.0, 0.0, OL,A,0,107),&
ratt_t1(26,'O(3P)S    ','O2        ','O3S       ','m         ',     0.0,&
  5.70e-34, -2.60, 0.00,  0.00e+00,  0.0, 0.0, 0.0, 0.0, T,0,0,107),&
ratt_t1(27,'OH        ','C2H4      ','HOC2H4O2  ','          ',     0.0,&
  8.60e-29, -3.10, 0.00,  9.00e-12, -0.85, 0.0, 0.0, 0.0, R,0,0,107),&
ratt_t1(28,'OH        ','C3H6      ','HOC3H6O2  ','          ',     0.0,&
  8.00e-27, -3.50, 0.00,  3.00e-11, -1.00, 0.0, 0.0, 0.0, R,0,0,107) /)

!----------------------------------------------------------------------
! NOTES: CheST Termolecular Reactions
!----------------------------------------------------------------------
! T001 O(3P)+O2 -> O3 m JPL 2011
! T001 IUPAC 2002 recommend k = 5.623E-34*(T/300)^-2.6 (based on
! T001 weighted mean of kO2+kN2)
!----------------------------------------------------------------------
! T002 O(3P)+NO -> NO2 m JPL 2011
! T002 IUPAC 2002 recommend k0 = 1.0E-3*(T/300)-1.6*[N2]
!----------------------------------------------------------------------
! T003 O(3P)+NO2 -> NO3 m JPL 2011
! T003 IUPAC 2002 recommend k0 = 1.3E-31*(T/300)^-1.5 kinf =
! T003 2.3E-11*(T/300)^0.24 Fc = 0.6
!----------------------------------------------------------------------
! T004 O(1D)+N2 -> N2O m JPL 2011
! T004 IUPAC 2002 k0=2.8E-36*[N2]
!----------------------------------------------------------------------
! T005 BrO+NO2 -> BrONO2 m JPL 2011
! T005 IUPAC Fc = 0.55 k0 = 4.2E-31*exp(T/300)^-2.4*[N2] kinf=2.7E-11
!----------------------------------------------------------------------
! T006 ClO+ClO -> Cl2O2 m JPL 2011
! T006 IUPAC recommend k0 = 2.0E-32*(T/300)^-4*[N2] kinf = 1.0E-11.
! T006 Fc=0.45 In general this is a problematic rate
!----------------------------------------------------------------------
! T007 Cl2O2+m -> ClO ClO IUPAC 2005
! T007 No JPL data
!----------------------------------------------------------------------
! T008 ClO+NO2 -> ClONO2 m JPL 2011
! T008 IUPAC recommend k0 = 1.6E-31*(T/300)^-3.4*[N2] kinf = 7.0E-11
! T008 Fc=0.4
!----------------------------------------------------------------------
! T009 H+O2 -> HO2 m JPL 2011
! T009 IUPAC 2009 recommend k0 = 4.3E-32*(T/300)^-1.2*[N2] kinf =
! T009 9.6E-11 and FC(ent) = 0.5
!----------------------------------------------------------------------
! T010 HO2+HO2 -> H2O2 O2 JPL 2011 see also asad_trimol.F90
! T010 IUPAC (2001) k = 1.9E-33*exp(980/T)*[N2]. Note that this reaction
! T010 is special and in the presence of H2O the rate constant needs to
! T010 be adjusted by: {1 + 1.4E-21*[H2O]*exp(2200/T)} (same H2O
! T010 expression for JPL). JPL also include the formation of HO2-H2O as
! T010 a separate species.
!----------------------------------------------------------------------
! T011 HO2+NO2 -> HO2NO2 m JPL 2011
! T011 IUPAC 2009 k0 = 1.4E-31*(T/300)^-3.1*[N2] kinf = 4.0E-12 Fc = 0.4
!----------------------------------------------------------------------
! T012 HO2NO2+m -> HO2 NO2 IUPAC 2001
! T012 No JPL data.. NOTE should multiply the expression by [N2] NOT [M]
!----------------------------------------------------------------------
! T013 OH+NO -> HONO m JPL 2011
! T013 IUPAC 2002 recommend k0 = 7.40E-31*(T/300)^-2.4 kinf =
! T013 3.3E-11*(T/300)^-0.3 Fc = 0.81
!----------------------------------------------------------------------
! T014 OH+NO2 -> HONO2 m JPL 2011
! T014 IUPCA 2009 recommend k0 = 3.3E-30*(T/300)^-3.0 kinf = 6.0E-11 Fc
! T014 = 0.4
!----------------------------------------------------------------------
! T015 OH+OH -> H2O2 m JPL 2011
! T015 IUPAC 2009 recommend k0 = 6.9E-31*(T/300)^-0.8 kinf =
! T015 3.9E-11*(T/300)^-0.47
!----------------------------------------------------------------------
! T016 MeCO3+NO2 -> PAN m MCMv3.2
! T016 Based on IUPAC 2003. JPL recommend k0 = 9.7E-29*(T/300)^-5.6 kinf
! T016 = 9.3E-12*(T/300)^-1.5 Fc = 0.6
!----------------------------------------------------------------------
! T017 PAN+m -> MeCO3 NO2 MCMv3.2
! T017 IUPAC 2003
!----------------------------------------------------------------------
! T018 EtCO3+NO2 -> PPAN m MCMv3.2
! T018 JPL recommend k0 = 9.0E-28*(T/300)^-8.9 kinf =
! T018 7.7E-12*(T/300)^-0.2 Fc = 0.6. IUPAC 2003 k=
! T018 1.70E-3*exp(-11280/T) kinf = 8.3E16*exp(-13940/T) Fc = 0.36
!----------------------------------------------------------------------
! T019 PPAN+m -> EtCO3 NO2 MCMv3.2
! T019 IUPAC 2003 k0 = 1.70E-03*exp(-11280/T) kinf =
! T019 8.30E+16*exp(-13940/T) Fc = 0.36
!----------------------------------------------------------------------
! T020 MACRO2+NO2 -> MPAN m MVMv3.2
! T020 -
!----------------------------------------------------------------------
! T021 MPAN+m -> MACRO2 NO2 MCMv3.2
! T021 IUPAC recommend k = 1.6E16*exp(-13500/T). No JPL data
!----------------------------------------------------------------------
! T022 NO2+NO3 -> N2O5 m IUPAC 2002
! T022 JPL recommend k0 = 2.0E-30*(T/300)^-4.4 kinf =
! T022 1.4E-12*(T/300)^-0.7 Fc = 0.6
!----------------------------------------------------------------------
! T023 N2O5+m -> NO2 NO3 IUPAC 2002
! T023 No JPL data. NOTE there is code in asad_trimol.F90 which looks
! T023 for this reaction and modifies the rate by an extra factor. I can
! T023 NOT find why it does this in the literature.
!----------------------------------------------------------------------
! T024 NO+NO -> NO2 NO2 IUPAC 2001
! T024 No JPL data..
!----------------------------------------------------------------------

! Dry deposition velocities. Defined in executable block because of length
! limit.
TYPE(depvel_t) :: depvel_defs_master(1:n_dry_master)

! The following formula is used to calculate the effective Henry's Law 
! coefficient, which takes the effects of dissociation and complex formation 
! on a species' solubility into account: (see Giannakopoulos, 1998)
!
!       H(eff) = K(298)exp{[-deltaH/R]x[(1/T)-(1/298)]}
!
! The data in columns 1 & 2 above give the data for this gas-aqueous transfer,
!       Column 1 = K(298) [M/atm]
!       Column 2 = -deltaH/R [K-1]
!
! If the species dissociates in the aqueous phase, the above term is multiplied
! by 1+{K(aq)/[H+]}, where
!       K(aq) = K(298)exp{[-deltaH/R]x[(1/T)-(1/298)]}
! The data in columns 3 & 4 give the data for this aqueous-phase dissociation,
!       Column 3 = K(298) [M]
!       Column 4 = -deltaH/R [K-1]
! The data in columns 5 and 6 give the data for a second dissociation,
! e.g for SO2, HSO3^{-}, and SO3^{2-}
!       Column 5 = K(298) [M]
!       Column 6 = -deltaH/R [K-1]

! Item number, species name, Henry's coefficients x6, chemistry scheme,
! qualifier, disqualifier, version
TYPE(wetdep), PARAMETER :: henry_defs_master(1:n_wet_master)=(/  &
!    1  
wetdep(1,'NO3       ',&
(/0.20e+01,0.20e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+T+ST+R,0,0,107),&
wetdep(1,'NO3       ',&
(/0.60e+00,0.00e+00,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),S,0,0,107),&
!    2 
wetdep(2,'N2O5      ',&
(/0.21e+06,0.87e+04,0.20e+02,0.00e+00,0.00e+00,0.00e+00/),TI+T+ST+R,0,0,107),&
wetdep(2,'N2O5      ',&
(/0.21e+06,0.87e+04,0.157e+02,0.00e+00,0.00e+00,0.00e+00/),S,0,0,107),&
!    3  
wetdep(3,'HO2NO2    ',&
(/0.13e+05,0.69e+04,0.10e-04,0.00e+00,0.00e+00,0.00e+00/),TI+T+ST+R,0,0,107),&
wetdep(3,'HO2NO2    ',&
(/0.20e+05,0.00e+00,0.10e-04,0.00e+00,0.00e+00,0.00e+00/),S,0,0,107),&  
!    4  
wetdep(4,'HONO2     ',&
(/0.21e+06,0.87e+04,0.20e+02,0.00e+00,0.00e+00,0.00e+00/),TI+T+ST+R,A,0,107),&!added A flag, LER, Apr2019
wetdep(4,'HONO2     ',&
(/0.21e+06,0.87e+04,0.157e+02,0.00e+00,0.00e+00,0.00e+00/),S,0,0,107),&
!    5
wetdep(5,'H2O2      ',&
(/0.83e+05,0.74e+04,0.24e-11,-0.373e+04,0.e+00,0.e+00/),TI+T+ST+OL+R,0,0,107),&
wetdep(5,'H2O2      ',&
(/0.83e+05,0.74e+04,0.22e-11,-0.373e+04,0.00e+00,0.00e+00/),S,0,0,107),&
!    6
wetdep(6,'HCHO      ',&
(/0.33e+04,0.65e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+T+ST+R,0,0,107),&
wetdep(6,'HCHO      ',&
(/0.30e+04,0.72e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),S,0,0,107),&
!    7  
wetdep(7,'MeOO      ',&
(/0.20e+04,0.66e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+T+ST+R,0,0,107),&
wetdep(7,'MeOO      ',&
(/0.20e+04,0.564e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),S,0,0,107),&
!    8
wetdep(8,'MeOOH     ',&
(/0.31e+03,0.50e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+T+ST+R,0,0,107),&
wetdep(8,'MeOOH     ',&
(/0.31e+03,0.52e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),S,0,0,107),&
!    9  
wetdep(9,'HO2       ',&
(/0.40e+04,0.59e+04,0.20e-04,0.00e+00,0.00e+00,0.00e+00/),TI+T+ST+R,0,0,107),&
wetdep(9,'HO2       ',&
(/0.40e+04,0.59e+04,0.16e-04,0.00e+00,0.00e+00,0.00e+00/),S,0,0,107),&
!   10  
wetdep(10,'BrONO2    ',&
(/0.21e+06,0.87e+04,0.157e+03,0.00e+00,0.00e+00,0.00e+00/),ST,0,0,107),&
wetdep(10,'BrONO2    ',&
(/0.00e+05,0.00e+00,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),S,0,0,107),&
!   11  
wetdep(11,'HCl       ',&
(/0.19e+02,0.60e+03,0.10e+05,0.00e+00,0.00e+00,0.00e+00/),S+ST,0,0,107),&
!   12
wetdep(12,'HOCl      ',&
(/0.93e+03,0.59e+04,0.32e-07,0.00e+00,0.00e+00,0.00e+00/),S+ST,0,0,107),&
!   13
wetdep(13,'HBr       ',&
(/0.13e+01,0.102e+05,0.10e+10,0.00e+00,0.00e+00,0.00e+00/),S+ST,0,0,107),&
!   14
wetdep(14,'HOBr      ',&
(/0.61e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),S+ST,A,0,107),&!added A flag, LER, Apr2019
!   15  
wetdep(15,'ClONO2    ',&
(/0.21e+06,0.87e+04,0.157e+02,0.00e+00,0.00e+00,0.00e+00/),S+ST,0,0,107),&
!   16  
wetdep(16,'HONO      ',&
(/0.50e+02,0.49e+04,0.56e-03,-0.126e+04,0.00e+00,0.00e+00/),TI+T+ST,0,0,107),&
!   17  
wetdep(17,'EtOOH     ',&
(/0.34e+03,0.57e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+T+ST+R,0,0,107),&
!   18
wetdep(18,'n-PrOOH   ',&
(/0.34e+03,0.57e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+T+ST,0,0,107),&
!   19  
wetdep(19,'i-PrOOH   ',&
(/0.34e+03,0.57e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+T+ST+R,0,0,107),&
!   20  
wetdep(20,'MeCOCH2OOH',&
(/0.34e+03,0.57e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+T+ST,0,0,107),&
!   21
wetdep(21,'ISOOH     ',&
(/0.17e+07,0.97e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+ST+R,0,0,107),&
!   22  
wetdep(22,'ISON      ',&
(/0.30e+04,0.74e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+ST+R,0,0,107),&
!   23  
wetdep(23,'MACROOH   ',&
(/0.17e+07,0.97e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+ST,0,0,107),&
!   24  
wetdep(24,'HACET     ',&
(/0.14e+03,0.72e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+ST,0,0,107),&
!   25  
wetdep(25,'MGLY      ',&
(/0.35e+04,0.72e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+ST+R,0,0,107),&
!   26
wetdep(26,'HCOOH     ',&
(/0.69e+04,0.56e+04,0.18e-03,-0.151e+04,0.00e+00,0.00e+00/),TI+ST,0,0,107),&
!   27  
wetdep(27,'MeCO3H    ',&
(/0.75e+03,0.53e+04,0.63e-08,0.00e+00,0.00e+00,0.00e+00/),TI+ST,0,0,107),&
!   28  
wetdep(28,'MeCO2H    ',&
(/0.47e+04,0.60e+04,0.18e-04,0.00e+00,0.00e+00,0.00e+00/),TI+ST,0,0,107),&
!   29  
wetdep(29,'MeOH      ',&
(/0.23e+03,0.49e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+ST+R,0,0,107),&
!   30  
wetdep(30,'O3        ',&
(/0.113e-01,0.23e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+S+ST+OL,A,0,107),&!Added 'A' flag, LER, Mar2019
!   31  
wetdep(31,'SO2       ',&
(/0.123e+01,0.302e+04,0.123e-01,0.201e+04,0.60e-07,0.112e+04/),TI+S+ST+OL,A,0,&
    107),&
!   32  
!wetdep(32,'DMSO      ',&
!(/0.50e+05,0.6425e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+ST+OL,A,0,107),&
!Updating DMSO Henry's Law coefficients as per Chen et al. (2018), LER, Mar2019
wetdep(32,'DMSO      ',&
(/1.00e+07,0.2580e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+ST+OL,A,0,107),&
!   33  
wetdep(33,'NH3       ',&
(/0.10e+07,0.00e+00,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),TI+ST,A,0,107),&
!   34  
wetdep(34,'Sec_Org   ',&
(/0.10e+06,0.12e+02,0.00e-00,0.00e+00,0.00e+00,0.00e+00/),TI+ST+OL,A,0,107), &
!   35
wetdep(35,'HO2S      ',&
(/0.40e+04,0.59e+04,0.20e-04,0.00e+00,0.00e+00,0.00e+00/),T,0,0,107),&
wetdep(36,'MVKOOH    ',&
(/0.17e+07,0.97e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),R,0,0,107),& 
wetdep(37,'ORGNIT    ',&
(/0.13e+03,0.00e+00,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),R,0,0,107), &
wetdep(38,'s-BuOOH   ',&
(/0.34e+03,0.57e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),R,0,0,107), & 
wetdep(39,'GLY       ',&
(/0.36e+06,0.00e+00,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),R,0,0,107), &
!Adding DMS, MSA and MSIA as per Chen et al. (2018), LER, UC, Mar2019
wetdep(40,'DMS       ',&
(/0.56,0.4480e+04,0.00e+00,0.00e+00,0.00e+00,0.00e+00/),ST,A,0,107), &
wetdep(41,'MSIA      ',&
(/1.00e+08,0.1760e+04,-1.00,-1.00,0.00e+00,0.00e+00/),ST,A,0,107), &
wetdep(42,'MSA       ',&
(/1.00e+09,0.1760e+04,-2.00,-2.00,0.00e+00,0.00e+00/),ST,A,0,107) /)
!For MSIA and MSA, -1 and -2 are used in the second dissociation position to flag
!that a pH correction needs to be applied to calculate the effective Henry's Law
!constants. See ukca_fracdiss.F90 for more details. LER, UC, Mar2019

! Bimolecular reactions are too many to define here in one statement.
TYPE(ratb_t1) :: ratb_defs_master(1:n_bimol_master) 


! To initialise chemistry files from specified species and rate arrays
! using the defined logicals. Composes chemistry definition arrays from
! components using logical definitions from RUN_UKCA namelist held in
! UKCA module ukca_option_mod


INTEGER            :: ierr
CHARACTER (LEN=errormessagelength) :: cmessage = ''

INTEGER :: chem_scheme
INTEGER :: qualifier
INTEGER :: i
INTEGER :: j
INTEGER :: k
LOGICAL :: found
INTEGER, PARAMETER :: version = 111 ! This is a version identifier. Set here
                    ! to the UM version when rates were added to the scheme.
                    ! Any entries can be considered whose version identifier is
                    ! <= this global version number. If there are multiple of 
                    ! those, take the one with the largest version number.
LOGICAL :: take_this
TYPE(CHCH_T1) :: last_chch  
TYPE(RATB_T1) :: last_bimol
TYPE(RATT_T1) :: last_termol
TYPE(RATJ_T1) :: last_photol
TYPE(RATH_T1) :: last_het
TYPE(DEPVEL_T):: last_depvel
TYPE(WETDEP)  :: last_wetdep

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CHEM_MASTER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Rates taken from  file "UKCA_Reaction_Rates - Bimolecular Reactions.csv"
! on 12:08 2012/05/01

! Item number, reactants x2, products x4, Arrhenius coefficients (k0, alpha
! and beta), fractional production rates x4, chemistry scheme, qualifier,
! disqualifier, version
ratb_defs_master(1:50) = (/ &
! B001 JPL2011
ratb_t1(1,'Br        ','Cl2O2     ','BrCl      ','Cl        ','O2        ',&
'          ',5.90e-12,  0.00,  170.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(1,'Br        ','Cl2O2     ','BrCl      ','Cl        ','O2        ',&
'          ',3.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! 3 B002 JPL2011
ratb_t1(2,'Br        ','HCHO      ','HBr       ','CO        ','HO2       ',&
'          ',1.70e-11,  0.00,  800.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B003 JPL2011
ratb_t1(3,'Br        ','HO2       ','HBr       ','O2        ','          ',&
'          ',4.80e-12,  0.00,  310.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(3,'Br        ','HO2       ','HBr       ','O2        ','          ',&
'          ',1.40e-11,  0.00,  590.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B004 JPL2011
ratb_t1(4,'Br        ','O3        ','BrO       ','O2        ','          ',&
'          ',1.60e-11,  0.00, 780.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(4,'Br        ','O3        ','BrO       ','O2        ','          ',&
'          ',1.70e-11,  0.00,  800.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B005 JPL2011
ratb_t1(5,'Br        ','OClO      ','BrO       ','ClO       ','          ',&
'          ',2.60e-11,  0.00, 1300.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(5,'Br        ','OClO      ','BrO       ','ClO       ','          ',&
'          ',2.70e-11,  0.00, 1300.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B006 JPL2011
ratb_t1(6,'BrO       ','BrO       ','Br        ','Br        ','O2        ',&
'          ',2.40e-12,  0.00,  -40.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B007 JPL2011
ratb_t1(7,'BrO       ','ClO       ','Br        ','Cl        ','O2        ',&
'          ',2.30e-12,  0.00, -260.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(7,'BrO       ','ClO       ','Br        ','Cl        ','O2        ',&
'          ',2.90e-12,  0.00, -220.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B008 JPL2011
ratb_t1(8,'BrO       ','ClO       ','Br        ','OClO      ','          ',&
'          ',9.50e-13,  0.00, -550.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(8,'BrO       ','ClO       ','Br        ','OClO      ','          ',&
'          ',1.60e-12,  0.00, -430.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B009 JPL2011
ratb_t1(9,'BrO       ','ClO       ','BrCl      ','O2        ','          ',&
'          ',4.10e-13,  0.00, -290.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(9,'BrO       ','ClO       ','BrCl      ','O2        ','          ',&
'          ',5.80e-13,  0.00, -170.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B010 JPL2011
ratb_t1(10,'BrO       ','HO2       ','HBr       ','O3        ','          ',&
'          ',0.00e+00,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(10,'BrO       ','HO2       ','HBr       ','O3        ','          ',&
'          ',7.00e-14,  0.00, -540.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B011 JPL2011
ratb_t1(11,'BrO       ','HO2       ','HOBr      ','O2        ','          ',&
'          ',4.50e-12,  0.00, -460.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(11,'BrO       ','HO2       ','HOBr      ','O2        ','          ',&
'          ',3.30e-12,  0.00, -540.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B012 JPL2011
ratb_t1(12,'BrO       ','NO        ','Br        ','NO2       ','          ',&
'          ',8.80e-12,  0.00, -260.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(12,'BrO       ','NO        ','Br        ','NO2       ','          ',&
'          ',8.70e-12,  0.00, -260.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B013 JPL2011
ratb_t1(13,'BrO       ','OH        ','Br        ','HO2       ','          ',&
'          ',1.70e-11,  0.00, -250.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(13,'BrO       ','OH        ','Br        ','HO2       ','          ',&
'          ',7.50e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B014 JPL2011
ratb_t1(14,'CF2Cl2    ','O(1D)     ','Cl        ','ClO       ','          ',&
'          ',1.40e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B015 JPL2011
ratb_t1(15,'CFCl3     ','O(1D)     ','Cl        ','Cl        ','ClO       ',&
'          ',2.30e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B016 JPL2011
ratb_t1(16,'Cl        ','CH4       ','HCl       ','MeOO      ','          ',&
'          ',7.30e-12,  0.00, 1280.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(16,'Cl        ','CH4       ','HCl       ','MeOO      ','          ',&
'          ',6.60e-12,  0.00, 1240.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B017 IUPAC2006
ratb_t1(17,'Cl        ','Cl2O2     ','Cl        ','Cl        ','Cl        ',&
'          ',7.60e-11,  0.00,  -65.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(17,'Cl        ','Cl2O2     ','Cl        ','Cl        ','Cl        ',&
'          ',1.00e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B018 JPL2011
ratb_t1(18,'Cl        ','ClONO2    ','Cl        ','Cl        ','NO3       ',&
'          ',6.50e-12,  0.00, -135.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B019 JPL2011
ratb_t1(19,'Cl        ','H2        ','HCl       ','H         ','          ',&
'          ',3.05e-11,  0.00, 2270.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(19,'Cl        ','H2        ','HCl       ','H         ','          ',&
'          ',3.90e-11,  0.00, 2310.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B020 JPL2011
ratb_t1(20,'Cl        ','H2O2      ','HCl       ','HO2       ','          ',&
'          ',1.10e-11,  0.00,  980.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B021 JPL2011
ratb_t1(21,'Cl        ','HCHO      ','HCl       ','CO        ','HO2       ',&
'          ',8.10e-11,  0.00,   30.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(21,'Cl        ','HCHO      ','HCl       ','CO        ','HO2       ',&
'          ',8.20e-11,  0.00,   34.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B022 JPL2011
ratb_t1(22,'Cl        ','HO2       ','ClO       ','OH        ','          ',&
'          ',3.65e-11,  0.00,  375.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(22,'Cl        ','HO2       ','ClO       ','OH        ','          ',&
'          ',4.10e-11,  0.00,  450.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B023 JPL2011
ratb_t1(23,'Cl        ','HO2       ','HCl       ','O2        ','          ',&
'          ',1.40e-11,  0.00, -270.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(23,'Cl        ','HO2       ','HCl       ','O2        ','          ',&
'          ',1.80e-11,  0.00, -170.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B024 JPL2011
ratb_t1(24,'Cl        ','HOCl      ','Cl        ','Cl        ','OH        ',&
'          ',3.40e-12,  0.00,  130.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(24,'Cl        ','HOCl      ','Cl        ','Cl        ','OH        ',&
'          ',2.50e-12,  0.00,  130.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B025 JPL2011
ratb_t1(25,'Cl        ','MeOOH     ','HCl       ','MeOO      ','          ',&
'          ',5.70e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B026 JPL2011
ratb_t1(26,'Cl        ','NO3       ','ClO       ','NO2       ','          ',&
'          ',2.40e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B027 JPL2011
ratb_t1(27,'Cl        ','O3        ','ClO       ','O2        ','          ',&
'          ',2.30e-11,  0.00,  200.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(27,'Cl        ','O3        ','ClO       ','O2        ','          ',&
'          ',2.90e-11,  0.00,  260.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B028 JPL2011
ratb_t1(28,'Cl        ','OClO      ','ClO       ','ClO       ','          ',&
'          ',3.40e-11,  0.00, -160.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(28,'Cl        ','OClO      ','ClO       ','ClO       ','          ',&
'          ',3.20e-11,  0.00, -170.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B029 JPL2011
ratb_t1(29,'ClO       ','ClO       ','Cl        ','Cl        ','O2        ',&
'          ',1.00e-12,  0.00, 1590.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B030 JPL2011
ratb_t1(30,'ClO       ','ClO       ','Cl        ','Cl        ','O2        ',&
'          ',3.00e-11,  0.00, 2450.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107) /)

ratb_defs_master(51:100) = (/&
! B031 JPL2011
ratb_t1(31,'ClO       ','ClO       ','Cl        ','OClO      ','          ',&
'          ',3.50e-13,  0.00, 1370.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B032 JPL2011
ratb_t1(32,'ClO       ','HO2       ','HCl       ','O3        ','          ',&
'          ',0.00e+00,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
ratb_t1(33,'ClO       ','HO2       ','HOCl      ','O2        ','          ',&
'          ',4.60e-13,  0.00, -710.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B033 JPL2011
ratb_t1(33,'ClO       ','HO2       ','HOCl      ','O2        ','          ',&
'          ',2.60e-12,  0.00, -290.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
! B034 JPL2011
ratb_t1(34,'ClO       ','MeOO      ','Cl        ','HCHO      ','HO2       ',&
'          ',3.30e-12,  0.00,  115.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B035 JPL2011
ratb_t1(35,'ClO       ','NO        ','Cl        ','NO2       ','          ',&
'          ',6.40e-12,  0.00, -290.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(35,'ClO       ','NO        ','Cl        ','NO2       ','          ',&
'          ',6.20e-12,  0.00, -295.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B036 JPL2011
ratb_t1(36,'ClO       ','NO3       ','Cl        ','O2        ','NO2       ',&
'          ',4.60e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B037 JPL2011
ratb_t1(37,'ClO       ','NO3       ','OClO      ','NO2       ','          ',&
'          ',0.00e+00,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B038 IUPAC2002
ratb_t1(38,'EtCO3     ','NO        ','EtOO      ','CO2       ','NO2       ',&
'          ',6.70e-12,  0.00, -340.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B039 MCMv3.2
ratb_t1(39,'EtCO3     ','NO3       ','EtOO      ','CO2       ','NO2       ',&
'          ',4.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B040 IUPAC2002
ratb_t1(40,'EtOO      ','MeCO3     ','MeCHO     ','HO2       ','MeOO      ',&
'          ',4.40e-13,  0.00,-1070.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B041 IUPAC2005
ratb_t1(41,'EtOO      ','NO        ','MeCHO     ','HO2       ','NO2       ',&
'          ',2.55e-12,  0.00, -380.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(41,'EtOO      ','NO        ','MeCHO     ','HO2       ','NO2       ',&
'          ',2.60e-12,  0.00, -380.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B042 IUPAC2008
ratb_t1(42,'EtOO      ','NO3       ','MeCHO     ','HO2       ','NO2       ',&
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B043 JPL2011
ratb_t1(43,'H         ','HO2       ','H2        ','O2        ','          ',&
'          ',6.90e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(43,'H         ','HO2       ','H2        ','O2        ','          ',&
'          ',2.35e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B044 JPL2011
ratb_t1(44,'H         ','HO2       ','O(3P)     ','H2O       ','          ',&
'          ',1.62e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B045 JPL2011
ratb_t1(45,'H         ','HO2       ','OH        ','OH        ','          ',&
'          ',7.20e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(45,'H         ','HO2       ','OH        ','OH        ','          ',&
'          ',5.59e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B046 JPL2011
ratb_t1(46,'H         ','NO2       ','OH        ','NO        ','          ',&
'          ',4.00e-10,  0.00,  340.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B047 JPL2011
ratb_t1(47,'H         ','O3        ','OH        ','O2        ','          ',&
'          ',1.40e-10,  0.00,  470.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
ratb_t1(48,'HO2       ','EtCO3     ','O2        ','EtCO3H    ','          ',&
'          ',3.05e-13,  0.00,  -1040.00, 0.00, 0.00, 0.00, 0.00,TI,0,0,107),&
! B048 MCMv3.2*
ratb_t1(48,'HO2       ','EtCO3     ','O2        ','EtCO3H    ','          ',&
'          ',4.40e-13,  0.00, -980.00, 0.00, 0.00, 0.00, 0.00,ST+T,0,0,107),&
ratb_t1(49,'HO2       ','EtCO3     ','O3        ','EtCO2H    ','          ',&
'          ',1.25e-13,  0.00,  -1040.00, 0.00, 0.00, 0.00, 0.00,TI,0,0,107),&
! B049 MCMv3.2
ratb_t1(49,'HO2       ','EtCO3     ','O3        ','EtCO2H    ','          ',&
'          ',7.80e-14,  0.00, -980.00, 0.00, 0.00, 0.00, 0.00,T+ST,0,0,107),&
! B050 IUPAC2011
ratb_t1(50,'HO2       ','EtOO      ','EtOOH     ','          ','          ',&
'          ',6.40e-13,  0.00, -710.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(50,'HO2       ','EtOO      ','EtOOH     ','          ','          ',&
'          ',3.80e-13,  0.00, -900.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B051 JPL2011 see also asad_bimol
ratb_t1(51,'HO2       ','HO2       ','H2O2      ','          ','          ',&
'          ',3.00e-13,  0.00, -460.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(51,'HO2       ','HO2       ','H2O2      ','          ','          ',&
'          ',2.20e-13,  0.00, -600.00, 0.00, 0.00, 0.00, 0.00,TI+T+OL,0,0,107),&
ratb_t1(51,'HO2       ','HO2       ','H2O2      ','O2        ','          ',&
'          ',2.20e-13,  0.00, -600.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B052 Poschl00
ratb_t1(52,'HO2       ','ISO2      ','ISOOH     ','          ','          ',&
'          ',2.05e-13,  0.00,  -1300.00, 0.00, 0.00, 0.00, 0.00,TI+ST,0,0,107),&
! B053 Poschl00
ratb_t1(53,'HO2       ','MACRO2    ','MACROOH   ','          ','          ',&
'          ',1.82e-13,  0.00,  -1300.00, 0.00, 0.00, 0.00, 0.00,TI+ST,0,0,107),&
! B054 IUPAC2009
ratb_t1(54,'HO2       ','MeCO3     ','MeCO2H    ','O3        ','          ',&
'          ',7.80e-14,  0.00, -980.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(54,'HO2       ','MeCO3     ','MeCO2H    ','O3        ','          ',&
'          ',1.04e-13,  0.00, -980.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B055 IUPAC2009
ratb_t1(55,'HO2       ','MeCO3     ','MeCO3H    ','          ','          ',&
'          ',2.13e-13,  0.00, -980.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(55,'HO2       ','MeCO3     ','MeCO3H    ','          ','          ',&
'          ',2.08e-13,  0.00, -980.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B056 IUPAC2009
ratb_t1(56,'HO2       ','MeCO3     ','OH        ','MeOO      ','          ',&
'          ',2.29e-13,  0.00, -980.00, 0.00, 0.00, 0.00, 0.00,T+ST,0,0,107),&
ratb_t1(56,'HO2       ','MeCO3     ','OH        ','MeOO      ','          ',&
'          ',2.08e-13,  0.00, -980.00, 0.00, 0.00, 0.00, 0.00,TI,0,0,107),&

ratb_t1(56,'HO2       ','MeCO3     ','O3        ','MeOO      ','MeCO3     ',&
'          ',5.20e-13,  0.00, -980.00, 0.300, 0.800, 0.200, 0.00,R,0,0,107),&
! B057 IUPAC2009
ratb_t1(57,'HO2       ','MeCOCH2OO ','MeCOCH2OOH','          ','          ',&
'          ',9.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(57,'HO2       ','MeCOCH2OO ','MeCOCH2OOH','          ','          ',&
'          ',1.36e-13,  0.00,  -1250.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B058 IUPAC2009 see also asad_bimol
ratb_t1(58,'HO2       ','MeOO      ','HCHO      ','          ','          ',&
'          ',3.80e-13,  0.00, -780.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B059 IUPAC2009 see also asad_bimol
! Is branching not handled consistently? No branching in R.
ratb_t1(59,'HO2       ','MeOO      ','MeOOH     ','          ','          ',&
'          ',3.80e-13,  0.00, -780.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
ratb_t1(59,'HO2       ','MeOO      ','O2        ','MeOOH     ','          ',&
'          ',3.80e-13,  0.00, -780.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B060 JPL2011
ratb_t1(60,'HO2       ','NO        ','OH        ','NO2       ','          ',&
'          ',3.30e-12,  0.00, -270.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(60,'HO2       ','NO        ','OH        ','NO2       ','          ',&
'          ',3.50e-12,  0.00, -250.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
ratb_t1(60,'HO2       ','NO        ','OH        ','NO2       ','          ',&
'          ',3.60e-12,  0.00, -270.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B061 JPL2011
ratb_t1(61,'HO2       ','NO3       ','OH        ','NO2       ','          ',&
'          ',3.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(61,'HO2       ','NO3       ','OH        ','NO2       ','          ',&
'          ',4.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107) /)

ratb_defs_master(101:150) = (/&
ratb_t1(61,'HO2       ','NO3       ','OH        ','NO2       ','O2        ',&
'          ',2.15e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
ratb_t1(500,'HO2       ','NO3       ','O2        ','HONO2     ','          ',&
'          ',2.15e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B062 IUPAC2001
ratb_t1(62,'HO2       ','O3        ','OH        ','O2        ','          ',&
'          ',2.03e-16,  4.57, -693.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
ratb_t1(62,'HO2       ','O3        ','OH        ','O2        ','O2        ',&
'          ',2.03e-16,  4.57, -693.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B063 MCMv3.2
ratb_t1(63,'HO2       ','i-PrOO    ','i-PrOOH   ','          ','          ',&
'          ',1.51e-13,  0.00,-1300.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
! B064 MCMv3.2
ratb_t1(64,'HO2       ','n-PrOO    ','n-PrOOH   ','          ','          ',&
'          ',1.51e-13,  0.00,-1300.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B065 Poschl00
ratb_t1(65,'ISO2      ','ISO2      ','MACR      ','MACR      ','HCHO      ',&
'HO2       ',2.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+ST,0,0,107),&
! B066 Poschl00
ratb_t1(66,'MACRO2    ','MACRO2    ','HACET     ','MGLY      ','HCHO      ',&
'CO        ',1.00e-12,  0.00,    0.00, 2.00, 2.00, 1.00, 1.00,TI+ST,0,0,107),&
! B067 Poschl00
ratb_t1(67,'MACRO2    ','MACRO2    ','HO2       ','          ','          ',&
'          ',1.00e-12,  0.00,    0.00, 2.00, 0.00, 0.00, 0.00,TI+ST,0,0,107),&
! B068 JPL2011
ratb_t1(68,'MeBr      ','Cl        ','Br        ','HCl       ','          ',&
'          ',1.40e-11,  0.00, 1030.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(68,'MeBr      ','Cl        ','Br        ','HCl       ','          ',&
'          ',1.70e-11,  0.00, 1080.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B069 JPL2011
ratb_t1(69,'MeBr      ','O(1D)     ','Br        ','OH        ','          ',&
'          ',1.80e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B070 JPL2011
ratb_t1(70,'MeBr      ','OH        ','Br        ','H2O       ','          ',&
'          ',2.35e-12,  0.00, 1300.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B071 IUPAC2002
ratb_t1(71,'MeCO3     ','NO        ','MeOO      ','CO2       ','NO2       ',&
'          ',7.50e-12,  0.00, -290.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
! B072 IUPAC2008
! Missing in R
ratb_t1(72,'MeCO3     ','NO3       ','MeOO      ','CO2       ','NO2       ',&
'          ',4.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B073 MCMv3.2
ratb_t1(73,'MeCOCH2OO ','NO        ','MeCO3     ','HCHO      ','NO2       ',&
'          ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(73,'MeCOCH2OO ','NO        ','MeCO3     ','HCHO      ','NO2       ',&
'          ',2.80e-12,  0.00, -300.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B074 MCMv3.2
ratb_t1(74,'MeCOCH2OO ','NO3       ','MeCO3     ','HCHO      ','NO2       ',&
'          ',2.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(74,'MeCOCH2OO ','NO3       ','MeCO3     ','HCHO      ','NO2       ',&
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B075 IUPAC2002
ratb_t1(75,'MeOO      ','MeCO3     ','HO2       ','HCHO      ','MeOO      ',&
'          ',1.80e-12,  0.00, -500.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
! B076 IUPAC2002
ratb_t1(76,'MeOO      ','MeCO3     ','MeCO2H    ','HCHO      ','          ',&
'          ',2.00e-13,  0.00, -500.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
! B077 IUPAC2002 see also asad_bimol
ratb_t1(77,'MeOO      ','MeOO      ','HO2       ','HO2       ','HCHO      ',&
'HCHO      ',1.03e-13,  0.00, -365.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
ratb_t1(77,'MeOO      ','MeOO      ','HCHO      ','HO2       ','O2        ',&
'          ',9.50e-14,  0.00, -390.00, 2.000, 2.000, 1.000, 0.00,S,0,0,107),&
! B078 IUPAC2002 see also asad_bimol
ratb_t1(78,'MeOO      ','MeOO      ','MeOH      ','HCHO      ','          ',&
'          ',1.03e-13,  0.00, -365.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
! B079 IUPAC2005
ratb_t1(79,'MeOO      ','NO        ','HO2       ','HCHO      ','NO2       ',&
'          ',2.30e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(79,'MeOO      ','NO        ','HO2       ','HCHO      ','NO2       ',&
'          ',2.95e-12,  0.00, -285.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
ratb_t1(79,'MeOO      ','NO        ','HCHO      ','HO2       ','NO2       ',&
'          ',2.95e-12,  0.00, -285.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B080 IUPAC2005
ratb_t1(80,'MeOO      ','NO        ','MeONO2    ','          ','          ',&
'          ',2.30e-15,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(80,'MeOO      ','NO        ','MeONO2    ','          ','          ',&
'          ',2.95e-15,  0.00, -285.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B081 MCMv3.2
ratb_t1(81,'MeOO      ','NO3       ','HO2       ','HCHO      ','NO2       ',&
'          ',  1.20e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(81,'MeOO      ','NO3       ','HO2       ','HCHO      ','NO2       ',&
'          ',1.30e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B082 JPL2011
ratb_t1(82,'N         ','NO        ','N2        ','O(3P)     ','          ',&
'          ',2.10e-11,  0.00, -100.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B083 JPL2011
ratb_t1(83,'N         ','NO2       ','N2O       ','O(3P)     ','          ',&
'          ',5.80e-12,  0.00, -220.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B084 JPL2011
ratb_t1(84,'N         ','O2        ','NO        ','O(3P)     ','          ',&
'          ',1.50e-11,  0.00, 3600.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B085 ????
! Not in strat scheme
ratb_t1(85,'N2O5      ','H2O       ','HONO2     ','HONO2     ','          ',&
'          ',2.50e-22,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B086 Poschl00
ratb_t1(86,'NO        ','ISO2      ','ISON      ','          ','          ',&
'          ',1.12e-13,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,TI+ST,0,0,107),&
! B087 Poschl00
ratb_t1(87,'NO        ','ISO2      ','NO2       ','MACR      ','HCHO      ',&
'HO2       ',2.43e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,TI+ST,0,0,107),&
! B088 Poschl00
ratb_t1(88,'NO        ','MACRO2    ','MGLY      ','HCHO      ','HO2       ',&
'          ',1.27e-12,  0.00, -360.00, 1.00, 1.50, 1.50, 0.00,TI+ST,0,0,107),&
! B089 Poschl00
ratb_t1(89,'NO        ','MACRO2    ','NO2       ','MeCO3     ','HACET     ',&
'CO        ',1.27e-12,  0.00, -360.00, 2.00, 0.50, 0.50, 0.50,TI+ST,0,0,107),&
! B090 JPL2011
ratb_t1(90,'NO        ','NO3       ','NO2       ','NO2       ','          ',&
'          ',1.50e-11,  0.00, -170.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(90,'NO        ','NO3       ','NO2       ','NO2       ','          ',&
'          ',1.80e-11,  0.00, -110.00, 0.00, 0.00, 0.00, 0.00,TI+T+S,0,0,107),&
! B091 JPL2011
ratb_t1(91,'NO        ','O3        ','NO2       ','          ','          ',&
'          ',3.00e-12,  0.00, 1500.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(91,'NO        ','O3        ','NO2       ','O2        ','          ',&
'          ',1.40e-12,  0.00, 1310.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
ratb_t1(91,'NO        ','O3        ','NO2       ','          ','          ',&
'          ',1.40e-12,  0.00, 1310.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B092 JPL2011
ratb_t1(92,'NO2       ','NO3       ','NO        ','NO2       ','O2        ',&
'          ',4.50e-14,  0.00, 1260.00, 0.0, 0.0, 0.0, 0.0,S+T+ST+R,0,0,107),&
! B093 JPL2011
ratb_t1(93,'NO2       ','O3        ','NO3       ','          ','          ',&
'          ',1.20e-13,  0.00, 2450.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(93,'NO2       ','O3        ','NO3       ','O2        ','          ',&
'          ',1.40e-13,  0.00, 2470.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
ratb_t1(93,'NO2       ','O3        ','NO3       ','          ','          ',&
'          ',1.40e-13,  0.00, 2470.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B094 JPL2011
ratb_t1(94,'NO3       ','Br        ','BrO       ','NO2       ','          ',&
'          ',1.60e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B095 IUPAC2007
ratb_t1(95,'NO3       ','C5H8      ','ISON      ','          ','          ',&
'          ',3.15e-12,  0.00,  450.00, 0.00, 0.00, 0.00, 0.00,TI+ST,0,0,107) /)

ratb_defs_master(151:200) = (/&
! I don't know how to map (NO3)C4H6CHO. The string is too long for ASAD.
ratb_t1(95,'NO3       ','C5H8      ','(NO3)C4H6CHO','HO2       ','          ',&
'          ',3.03e-12,  0.00,  446.00, 0.00, 0.00, 0.00, 0.00,R,0,0,107),&
! B096 IUPAC2007
ratb_t1(96,'NO3       ','EtCHO     ','HONO2     ','EtCO3     ','          ',&
'          ',6.30e-15,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(96,'NO3       ','EtCHO     ','HONO2     ','EtCO3     ','          ',&
'          ',3.46e-12,  0.00, 1862.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B097 IUPAC2007
ratb_t1(97,'NO3       ','HCHO      ','HONO2     ','HO2       ','CO        ',&
'          ',2.00e-12,  0.00, 2440.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
ratb_t1(97,'NO3       ','HCHO      ','HONO2     ','CO        ','HO2       ',&
'          ',5.60e-16,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B098 MCMv3.2
ratb_t1(98,'NO3       ','MGLY      ','MeCO3     ','CO        ','HONO2     ',&
'          ',3.36e-12,  0.00, 1860.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(98,'NO3       ','MGLY      ','MeCO3     ','CO        ','HONO2     ',&
'          ',3.46e-12,  0.00, 1860.00, 0.00, 0.00, 0.00, 0.00,TI,0,0,107),&
! B099 IUPAC2007
ratb_t1(99,'NO3       ','Me2CO     ','HONO2     ','MeCOCH2OO ','          ',&
'          ',3.00e-17,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B100 IUPAC2007
ratb_t1(100,'NO3       ','MeCHO     ','HONO2     ','MeCO3     ','          ',&
'          ',1.40e-12,  0.00, 1860.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
! B101 JPL2011
ratb_t1(101,'O(1D)     ','CH4       ','HCHO      ','H2        ','          ',&
'          ',9.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(101,'O(1D)     ','CH4       ','HCHO      ','H2        ','          ',&
'          ',7.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
ratb_t1(101,'O(1D)     ','CH4       ','HCHO      ','H2        ','          ',&
'          ',1.50e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B102 JPL2011
ratb_t1(102,'O(1D)     ','CH4       ','HCHO      ','HO2       ','HO2       ',&
'          ',3.45e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B103 JPL2011
ratb_t1(103,'O(1D)     ','CH4       ','OH        ','MeOO      ','          ',&
'          ',1.31e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(103,'O(1D)     ','CH4       ','OH        ','MeOO      ','          ',&
'          ',1.35e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
ratb_t1(103,'O(1D)     ','CH4       ','OH        ','MeOO      ','          ',&
'          ',1.05e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B104 JPL2011
! Not in TI scheme
ratb_t1(104,'O(1D)     ','CO2       ','O(3P)     ','CO2       ','          ',&
'          ',7.50e-11,  0.00, -115.00, 0.00, 0.00, 0.00, 0.00,T+ST,0,0,107),&
ratb_t1(104,'O(1D)     ','CO2       ','O(3P)     ','CO2       ','          ',&
'          ',7.40e-11,  0.00, -120.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B105 IUPAC2008
! Not in TI scheme
ratb_t1(105,'O(1D)     ','H2        ','OH        ','H         ','          ',&
'          ',1.20e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,T+ST,0,0,107),&
ratb_t1(105,'O(1D)     ','H2        ','OH        ','H         ','          ',&
'          ',1.10e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B106 JPL2011
ratb_t1(106,'O(1D)     ','H2O       ','OH        ','OH        ','          ',&
'          ',1.63e-10,  0.00,  -60.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(106,'O(1D)     ','H2O       ','OH        ','OH        ','          ',&
'          ',2.20e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+S+T,0,0,107),&
! B107 JPL2011
ratb_t1(107,'O(1D)     ','HBr       ','HBr       ','O(3P)     ','          ',&
'          ',3.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B108 JPL2011
ratb_t1(108,'O(1D)     ','HBr       ','OH        ','Br        ','          ',&
'          ',1.20e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B109 JPL2011
ratb_t1(109,'O(1D)     ','HCl       ','H         ','ClO       ','          ',&
'          ',3.60e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B110 JPL2011
ratb_t1(110,'O(1D)     ','HCl       ','O(3P)     ','HCl       ','          ',&
'          ',1.35e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B111 JPL2011
ratb_t1(111,'O(1D)     ','HCl       ','OH        ','Cl        ','          ',&
'          ',1.01e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B112 JPL2011
ratb_t1(112,'O(1D)     ','N2        ','O(3P)     ','N2        ','          ',&
'          ',2.15e-11,  0.00, -110.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(112,'O(1D)     ','N2        ','O(3P)     ','N2        ','          ',&
'          ',2.10e-11,  0.00, -115.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
ratb_t1(112,'O(1D)     ','N2        ','O(3P)     ','N2        ','          ',&
'          ',1.80e-11,  0.00, -110.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B113 JPL2011
ratb_t1(113,'O(1D)     ','N2O       ','N2        ','O2        ','          ',&
'          ',4.60e-11,  0.00,  -20.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(113,'O(1D)     ','N2O       ','N2        ','O2        ','          ',&
'          ',4.90e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B114 JPL2011
ratb_t1(114,'O(1D)     ','N2O       ','NO        ','NO        ','          ',&
'          ',7.30e-11,  0.00,  -20.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(114,'O(1D)     ','N2O       ','NO        ','NO        ','          ',&
'          ',6.70e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B115 JPL2011
ratb_t1(115,'O(1D)     ','O2        ','O(3P)     ','O2        ','          ',&
'          ',3.30e-11,  0.00,  -55.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(115,'O(1D)     ','O2        ','O(3P)     ','O2        ','          ',&
'          ',3.20e-11,  0.00,  -67.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
ratb_t1(115,'O(1D)     ','O2        ','O(3P)     ','O2        ','          ',&
'          ',3.20e-11,  0.00,  -70.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B116 JPL2011
! Not in TI scheme
ratb_t1(116,'O(1D)     ','O3        ','O2        ','O(3P)     ','O(3P)     ',&
'          ',1.20e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+T+ST,0,0,107),&
! B117 JPL2011
! Not in TI scheme
ratb_t1(117,'O(1D)     ','O3        ','O2        ','O2        ','          ',&
'          ',1.20e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+T+ST,0,0,107),&
! B118 JPL2011
ratb_t1(118,'O(3P)     ','BrO       ','O2        ','Br        ','          ',&
'          ',1.90e-11,  0.00, -230.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B119 JPL2011
ratb_t1(119,'O(3P)     ','ClO       ','Cl        ','O2        ','          ',&
'          ',2.80e-11,  0.00,  -85.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(119,'O(3P)     ','ClO       ','Cl        ','O2        ','          ',&
'          ',3.80e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B120 JPL2011
ratb_t1(120,'O(3P)     ','ClONO2    ','ClO       ','NO3       ','          ',&
'          ',3.60e-12,  0.00,  840.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(120,'O(3P)     ','ClONO2    ','ClO       ','NO3       ','          ',&
'          ',2.90e-12,  0.00,  800.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B121 ????
! Not in TI scheme
ratb_t1(121,'O(3P)     ','H2        ','OH        ','HO2       ','          ',&
'          ',  9.00e-18,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,T,0,0,107),&
ratb_t1(121,'O(3P)     ','H2        ','OH        ','H         ','          ',&
'          ',9.00e-18,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B122 JPL2011
! Not in T/TI/R schemes
ratb_t1(122,'O(3P)     ','H2O2      ','OH        ','HO2       ','          ',&
'          ',1.40e-12,  0.00, 2000.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B123 JPL2011
ratb_t1(123,'O(3P)     ','HBr       ','OH        ','Br        ','          ',&
'          ',5.80e-12,  0.00, 1500.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B124 JPL2011
! Not in TI/R schemes
ratb_t1(124,'O(3P)     ','HCHO      ','OH        ','CO        ','HO2       ',&
'          ',3.40e-11,  0.00, 1600.00, 0.00, 0.00, 0.00, 0.00,S+T+ST,0,0,107),&
! B125 JPL2011
ratb_t1(125,'O(3P)     ','HCl       ','OH        ','Cl        ','          ',&
'          ',1.00e-11,  0.00, 3300.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107)/)

ratb_defs_master(201:250) = (/&
! B126 IUPAC2001
! Not in TI/R schemes
ratb_t1(126,'O(3P)     ','HO2       ','OH        ','O2        ','          ',&
'          ',2.70e-11,  0.00, -224.00, 0.00, 0.00, 0.00, 0.00,S+T+ST,0,0,107),&
! B127 JPL2011
ratb_t1(127,'O(3P)     ','HOCl      ','OH        ','ClO       ','          ',&
'          ',1.70e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B128 JPL2011
ratb_t1(128,'O(3P)     ','NO2       ','NO        ','O2        ','          ',&
'          ',5.10e-12,  0.00, -210.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(128,'O(3P)     ','NO2       ','O2        ','NO        ','          ',&
'          ',5.60e-12,  0.00, -180.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
ratb_t1(128,'O(3P)     ','NO2       ','NO        ','O2        ','          ',&
'          ',5.50e-12,  0.00, -188.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B129 IUPAC2009
! Not in TI/R schemes
ratb_t1(129,'O(3P)     ','NO3       ','O2        ','NO2       ','          ',&
'          ',1.70e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+T+ST,0,0,107),&
! B130 JPL2011
ratb_t1(130,'O(3P)     ','O3        ','O2        ','O2        ','          ',&
'          ',8.00e-12,  0.00, 2060.00, 0.0, 0.0, 0.0, 0.0,TI+S+T+ST,0,0,107),&
! B131 JPL2011
ratb_t1(131,'O(3P)     ','OClO      ','O2        ','ClO       ','          ',&
'          ',2.40e-12,  0.00,  960.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B132 JPL2011
! Not in TI/R scheme
ratb_t1(132,'O(3P)     ','OH        ','O2        ','HO2         ','          ',&
'          ',1.80e-11,  0.00, -180.00, 0.00, 0.00, 0.00, 0.00,T,0,0,107),&
ratb_t1(132,'O(3P)     ','OH        ','O2        ','H         ','          ',&
'          ',1.80e-11,  0.00, -180.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(132,'O(3P)     ','OH        ','O2        ','H         ','          ',&
'          ',2.40e-11,  0.00, -110.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B133 IUPAC2007*
ratb_t1(133,'O3        ','C5H8      ','HO2       ','OH        ','          ',&
'          ',3.33e-15,  0.00, 1995.00, 0.75, 0.75, 0.00, 0.00,TI+ST,0,0,107),&
! B134 IUPAC2007*
ratb_t1(134,'O3        ','C5H8      ','MACR      ','HCHO      ','MACRO2    ',&
'MeCO3     ',3.33e-15,  0.00, 1995.00, 1.95, 1.74, 0.30, 0.30,TI+ST,0,0,107),&
! B135 IUPAC2007*
ratb_t1(135,'O3        ','C5H8      ','MeOO      ','HCOOH     ','CO        ',&
'H2O2      ',3.33e-15,  0.00, 1995.00, 0.24, 0.84, 0.42, 0.27,TI+ST,0,0,107),&
! Don't know how to map CH2O2.
ratb_t1(135,'O3        ','C5H8      ','MVK       ','CO        ','CH2O2     ',&
'HO2       ',3.93e-15,  0.00, 1913.00, 2.00, 1.56, 0.44, 0.54,R,0,0,107),&
ratb_t1(400,'O3        ','C5H8      ','OH        ','          ','          ',&
'          ',3.93e-15,  0.00, 1913.00, 0.54, 0.00, 0.00, 0.00,R,0,0,107),&
! B136 IUPAC2007*
ratb_t1(136,'O3        ','MACR      ','MGLY      ','HCOOH     ','HO2       ',&
'CO        ',2.13e-16,  0.00, 1520.00, 1.80, 0.90, 0.64, 0.44,TI+ST,0,0,107),&
! B137 IUPAC2007*
ratb_t1(137,'O3        ','MACR      ','MGLY      ','HCOOH     ','HO2       ',&
'CO        ',3.50e-16,  0.00, 2100.00, 1.80, 0.90, 0.64, 0.44,TI+ST,0,0,107),&
! B138 IUPAC2007*
ratb_t1(138,'O3        ','MACR      ','OH        ','MeCO3     ','          ',&
'          ',2.13e-16,  0.00, 1520.00, 0.38, 0.20, 0.00, 0.00,TI+ST,0,0,107),&
! B139 IUPAC2007*
ratb_t1(139,'O3        ','MACR      ','OH        ','MeCO3     ','          ',&
'          ',3.50e-16,  0.00, 2100.00, 0.38, 0.20, 0.00, 0.00,TI+ST,0,0,107),&
! B140 JPL2011
ratb_t1(140,'OClO      ','NO        ','NO2       ','ClO       ','          ',&
'          ',2.50e-12,  0.00,  600.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(140,'OClO      ','NO        ','NO2       ','ClO       ','          ',&
'          ',3.40e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B141 IUPAC2007
ratb_t1(141,'OH        ','C2H6      ','H2O       ','EtOO      ','          ',&
'          ',6.90e-12,  0.00, 1000.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
! B142 IUPAC2007 see also asad_bimol
ratb_t1(142,'OH        ','C3H8      ','i-PrOO    ','H2O       ','          ',&
'          ',7.60e-12,  0.00,  585.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
! RAQ does not do n-PrOO production; this would get only half the loss rate 
! for C3H8.
! B143 IUPAC2007 see also asad_bimol
ratb_t1(143,'OH        ','C3H8      ','n-PrOO    ','H2O       ','          ',&
'          ',7.60e-12,  0.00,  585.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B144 IUPAC2009
ratb_t1(144,'OH        ','C5H8      ','ISO2      ','          ','          ',&
'          ',2.70e-11,  0.00, -390.00, 0.00, 0.00, 0.00, 0.00,TI+ST,0,0,107),&
ratb_t1(144,'OH        ','C5H8      ','HOIPO2    ','H2O       ','          ',&
'          ',2.54e-11,  0.00, -410.00, 0.00, 0.00, 0.00, 0.00,R,0,0,107),&
! B145 JPL2011
ratb_t1(145,'OH        ','CH4       ','H2O       ','MeOO      ','          ',&
'          ',2.45e-12,  0.00, 1775.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(145,'OH        ','CH4       ','H2O       ','MeOO      ','          ',&
'          ',1.85e-12,  0.00, 1690.00, 0.00, 0.00, 0.00, 0.00,TI+S+T,0,0,107),&
! B146 IUPAC2005 see also asad_bimol
ratb_t1(146,'OH        ','CO        ','HO2       ','          ','          ',&
'          ',1.44e-13,  0.00,    0.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
ratb_t1(146,'OH        ','CO        ','H         ','CO2       ','          ',&
'          ',1.30e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B147 JPL2011
ratb_t1(147,'OH        ','ClO       ','HCl       ','O2        ','          ',&
'          ',6.00e-13,  0.00, -230.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B148 JPL2011
ratb_t1(148,'OH        ','ClO       ','HO2       ','Cl        ','          ',&
'          ',7.40e-12,  0.00, -270.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B149 JPL2011
ratb_t1(149,'OH        ','ClONO2    ','HOCl      ','NO3       ','          ',&
'          ',1.20e-12,  0.00,  330.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B150 IUPAC2007
ratb_t1(150,'OH        ','EtCHO     ','H2O       ','EtCO3     ','          ',&
'          ',4.90e-12,  0.00, -405.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(150,'OH        ','EtCHO     ','H2O       ','EtCO3     ','          ',&
'          ',5.10e-12,  0.00, -405.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B151 MCMv3.2
ratb_t1(151,'OH        ','EtOOH     ','H2O       ','EtOO      ','          ',&
'          ',1.90e-12,  0.00, -190.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B152 MCMv3.2
ratb_t1(152,'OH        ','EtOOH     ','H2O       ','MeCHO     ','OH        ',&
'          ',8.01e-12,  0.00,    0.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
! B153 JPL2011
! Products should be H2O + H not H2O + HO2 in ST chemistry.
ratb_t1(153,'OH        ','H2        ','H2O       ','HO2       ','          ',&
'          ',2.80e-12,  0.00, 1800.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(153,'OH        ','H2        ','H2O       ','HO2       ','          ',&
'          ',7.70e-12,  0.00, 2100.00, 0.00, 0.00, 0.00, 0.00,TI+T+R,0,0,107),&
ratb_t1(153,'OH        ','H2        ','H2O       ','H         ','          ',&
'          ',7.70e-12,  0.00, 2100.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B154 IUPAC2001
ratb_t1(154,'OH        ','H2O2      ','H2O       ','HO2       ','          ',&
'          ',2.90e-12,  0.00,  160.00, 0.0, 0.0, 0.0, 0.0,TI+S+T+ST+R,0,0,107),&
ratb_t1(154,'OH        ','H2O2      ','H2O       ','          ','          ',&
'          ',2.90e-12,  0.00,  160.00, 0.00, 0.00, 0.00, 0.00,OL,0,0,107),&
! B155 IUPAC2007
ratb_t1(155,'OH        ','HACET     ','MGLY      ','HO2       ','          ',&
'          ',1.60e-12,  0.00, -305.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(155,'OH        ','HACET     ','MGLY      ','HO2       ','          ',&
'          ',3.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI,0,0,107),&
! B156 JPL2011
ratb_t1(156,'OH        ','HBr       ','H2O       ','Br        ','          ',&
'          ',5.50e-12,  0.00, -200.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(156,'OH        ','HBr       ','H2O       ','Br        ','          ',&
'          ',1.10e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B157 IUPAC2007
ratb_t1(157,'OH        ','HCHO      ','H2O       ','HO2       ','CO        ',&
'          ',5.40e-12,  0.00, -135.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
ratb_t1(157,'OH        ','HCHO      ','H2O       ','CO        ','HO2       ',&
'          ',8.20e-12,  0.00,  -40.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B158 IUPAC2007
! more products needed here: OH + HCOOH -> H2O + CO2 + H
ratb_t1(158,'OH        ','HCOOH     ','HO2       ','          ','          ',&
'          ',4.50e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+ST,0,0,107)/)

ratb_defs_master(251:300) = (/&
! B159 JPL2011
ratb_t1(159,'OH        ','HCl       ','H2O       ','Cl        ','          ',&
'          ',1.80e-12,  0.00,  250.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(159,'OH        ','HCl       ','H2O       ','Cl        ','          ',&
'          ',1.80e-12,  0.00,  240.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B160 JPL2011
ratb_t1(160,'OH        ','HO2       ','H2O       ','          ','          ',&
'          ',4.80e-11,  0.00, -250.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
ratb_t1(160,'OH        ','HO2       ','H2O       ','O2        ','          ',&
'          ',4.80e-11,  0.00, -250.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B161 IUPAC2007
ratb_t1(161,'OH        ','HO2NO2    ','H2O       ','NO2       ','O2        ',&
'          ',3.20e-13,  0.00, -690.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(161,'OH        ','HO2NO2    ','H2O       ','NO2       ','O2        ',&
'          ',1.30e-12,  0.00, -380.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
ratb_t1(161,'OH        ','HO2NO2    ','H2O       ','NO2       ','          ',&
'          ',1.90e-12,  0.00, -270.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B162 JPL2011
ratb_t1(162,'OH        ','HOCl      ','ClO       ','H2O       ','          ',&
'          ',3.00e-12,  0.00,  500.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
! B163 IUPAC2004
ratb_t1(163,'OH        ','HONO      ','H2O       ','NO2       ','          ',&
'          ',2.50e-12,  0.00, -260.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B164 IUPAC2004 see also asad_bimol
ratb_t1(164,'OH        ','HONO2     ','H2O       ','NO3       ','          ',&
'          ',2.40e-14,  0.00, -460.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(164,'OH        ','HONO2     ','H2O       ','NO3       ','          ',&
'          ',1.50e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T+S,0,0,107),&
! B165 Poschl00
ratb_t1(165,'OH        ','ISON      ','HACET     ','NALD      ','          ',&
'          ',1.30e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+ST,0,0,107),&
! B166 Poschl00
ratb_t1(166,'OH        ','ISOOH     ','MACR      ','OH        ','          ',&
'          ',1.00e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+ST,0,0,107),&
! Very different rates & products
ratb_t1(166,'OH        ','ISOOH     ','MVK       ','HCHO      ','OH        ',&
'          ',4.20e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,R,0,0,107),&
! B167 IUPAC2007
ratb_t1(167,'OH        ','MACR      ','MACRO2    ','          ','          ',&
'          ',1.30e-12,  0.00, -610.00, 0.00, 0.00, 0.00, 0.00,TI+ST,0,0,107),&
! B168 IUPAC2007
ratb_t1(168,'OH        ','MACR      ','MACRO2    ','          ','          ',&
'          ',4.00e-12,  0.00, -380.00, 0.00, 0.00, 0.00, 0.00,TI+ST,0,0,107),&
! B169 MCMv3.2
ratb_t1(169,'OH        ','MACROOH   ','MACRO2    ','          ','          ',&
'          ',3.77e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(169,'OH        ','MACROOH   ','MACRO2    ','          ','          ',&
'          ',3.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI,0,0,107),&
! B170 IUPAC2008
ratb_t1(170,'OH        ','MGLY      ','MeCO3     ','CO        ','          ',&
'          ',1.90e-12,  0.00, -575.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(170,'OH        ','MGLY      ','MeCO3     ','CO        ','          ',&
'          ',1.50e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI,0,0,107),&
! Note the different rates
ratb_t1(170,'OH        ','MGLY      ','MeCO3     ','CO        ','          ',&
'          ',1.72e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,R,0,0,107), &
! B171 IUPAC2006
ratb_t1(171,'OH        ','MPAN      ','HACET     ','NO2       ','          ',&
'          ',2.90e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+ST,0,0,107),&
! B172 IUPAC2007
ratb_t1(172,'OH        ','Me2CO     ','H2O       ','MeCOCH2OO ','          ',&
'          ',1.70e-14,  0.00, -423.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(172,'OH        ','Me2CO     ','H2O       ','MeCOCH2OO ','          ',&
'          ',1.70e-14,  0.00, -420.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B173 IUPAC2007
ratb_t1(173,'OH        ','Me2CO     ','H2O       ','MeCOCH2OO ','          ',&
'          ',8.80e-12,  0.00, 1320.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
! B174 IUPAC2009
ratb_t1(174,'OH        ','MeCHO     ','H2O       ','MeCO3     ','          ',&
'          ',4.70e-12,  0.00, -345.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
ratb_t1(174,'OH        ','MeCHO     ','H2O       ','MeCO3     ','          ',&
'          ',4.40e-12,  0.00, -365.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B175 MCMv3.2
ratb_t1(175,'OH        ','MeCO2H    ','MeOO      ','          ','          ',&
'          ',8.00e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(175,'OH        ','MeCO2H    ','MeOO      ','          ','          ',&
'          ',4.00e-13,  0.00, -200.00, 0.00, 0.00, 0.00, 0.00,TI,0,0,107),&
! B176 MCMv3.2
ratb_t1(176,'OH        ','MeCO3H    ','MeCO3     ','          ','          ',&
'          ',3.70e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+ST,0,0,107),&
! B177 MCMv3.2
ratb_t1(177,'OH        ','MeCOCH2OOH','H2O       ','MeCOCH2OO ','          ',&
'          ',1.90e-12,  0.00, -190.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B178 MCMv3.2
ratb_t1(178,'OH        ','MeCOCH2OOH','OH        ','MGLY      ','          ',&
'          ',8.39e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B179 IUPAC2006
ratb_t1(179,'OH        ','MeOH      ','HO2       ','HCHO      ','          ',&
'          ',2.85e-12,  0.00,  345.00, 0.00, 0.00, 0.00, 0.00,TI+ST+R,0,0,107),&
! B180 IUPAC2006
ratb_t1(180,'OH        ','MeONO2    ','HCHO      ','NO2       ','H2O       ',&
'          ',4.00e-13,  0.00,  845.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B060a added Alex
! Not in TI, S, ST and R schemes
ratb_t1( 60,'HO2       ','NO        ','HONO2     ','          ','          ',&
'          ',3.60e-12,  0.00, -270.00, 0.00, 0.00, 0.00, 0.00,T,0,0,107),&
! B181 IUPAC2007
! branching not in S scheme
ratb_t1(181,'OH        ','MeOOH     ','H2O       ','HCHO      ','OH        ',&
'          ',2.12e-12,  0.00, -190.00, 0.00, 0.00, 0.00, 0.00,ST+R,0,0,107),&
! Note the different rates
ratb_t1(181,'OH        ','MeOOH     ','H2O       ','HCHO      ','OH        ',&
'          ',1.02e-12,  0.00, -190.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B182 IUPAC2007
ratb_t1(182,'OH        ','MeOOH     ','H2O       ','MeOO      ','          ',&
'          ',1.89e-12,  0.00, -190.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
ratb_t1(182,'OH        ','MeOOH     ','H2O       ','MeOO      ','          ',&
'          ',2.70e-12,  0.00, -200.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B183 IUPAC2009
ratb_t1(183,'OH        ','NALD      ','HCHO      ','CO        ','NO2       ',&
'          ',4.70e-12,  0.00, -345.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(183,'OH        ','NALD      ','HCHO      ','CO        ','NO2       ',&
'          ',4.40e-12,  0.00, -365.00, 0.00, 0.00, 0.00, 0.00,TI,0,0,107),&
! B184 JPL2011
ratb_t1(184,'OH        ','NO3       ','HO2       ','NO2       ','          ',&
'          ',2.20e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,S+ST,0,0,107),&
ratb_t1(184,'OH        ','NO3       ','HO2       ','NO2       ','          ',&
'          ',2.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B185 JPL2011
ratb_t1(185,'OH        ','O3        ','HO2       ','O2        ','          ',&
'          ',1.70e-12,  0.00,  940.00, 0.0, 0.0, 0.0, 0.0,TI+T+S+ST+R,0,0,107),&
! B186 JPL2011
ratb_t1(186,'OH        ','OClO      ','HOCl      ','O2        ','          ',&
'          ',1.40e-12,  0.00, -600.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(186,'OH        ','OClO      ','HOCl      ','O2        ','          ',&
'          ',4.50e-13,  0.00, -800.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B187 IUPAC2001
! Note the different signs of temperature dependence
ratb_t1(187,'OH        ','OH        ','H2O       ','O(3P)     ','          ',&
'          ',6.31e-14,  2.60, -945.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
ratb_t1(187,'OH        ','OH        ','H2O       ','O(3P)     ','          ',&
'          ',4.20e-12,  0.00,  240.00, 0.00, 0.00, 0.00, 0.00,S,0,0,107),&
! B188 MCMv3.2
ratb_t1(188,'OH        ','PAN       ','HCHO      ','NO2       ','H2O       ',&
'          ',3.00e-14,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
ratb_t1(188,'OH        ','PAN       ','HCHO      ','NO3       ','          ',&
'          ',9.50e-13,  0.00,  650.00, 0.00, 0.00, 0.00, 0.00,R,0,0,107) /)

ratb_defs_master(301:350) = (/&
! B189 MCMv3.2
ratb_t1(189,'OH        ','PPAN      ','MeCHO     ','NO2       ','H2O       ',&
'          ',1.27e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B190 MCMv3.2
ratb_t1(190,'OH        ','i-PrOOH   ','Me2CO     ','OH        ','          ',&
'          ',1.66e-11,  0.00,    0.00, 0.00, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
! B191 MCMv3.2
ratb_t1(191,'OH        ','i-PrOOH   ','i-PrOO    ','H2O       ','          ',&
'          ',1.90e-12,  0.00, -190.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B192 MCMv3.2
ratb_t1(192,'OH        ','n-PrOOH   ','EtCHO     ','H2O       ','OH        ',&
'          ',1.10e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B193 MCMv3.2
ratb_t1(193,'OH        ','n-PrOOH   ','n-PrOO    ','H2O       ','          ',&
'          ',1.90e-12,  0.00, -190.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B194 IUPAC2005
ratb_t1(194,'i-PrOO    ','NO        ','Me2CO     ','HO2       ','NO2       ',&
'          ',2.70e-12,  0.00, -360.00, 0.0, 0.0, 0.0, 0.0,TI+T+ST+R,0,0,107),&
! B195 MCMv3.2
ratb_t1(195,'i-PrOO    ','NO3       ','Me2CO     ','HO2       ','NO2       ',&
'          ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(195,'i-PrOO    ','NO3       ','Me2CO     ','HO2       ','NO2       ',&
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,TI+T,0,0,107),&
! B196 IUPAC2005
ratb_t1(196,'n-PrOO    ','NO        ','EtCHO     ','HO2       ','NO2       ',&
'          ',2.90e-12,  0.00, -350.00, 0.00, 0.00, 0.00, 0.00,TI+T+ST,0,0,107),&
! B197 MCM3.2
ratb_t1(197,'n-PrOO    ','NO3       ','EtCHO     ','HO2       ','NO2       ',&
'          ',2.70e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00,ST,0,0,107),&
ratb_t1(197,'n-PrOO    ','NO3       ','EtCHO     ','HO2       ','NO2       ',&
'          ',2.50e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00,T+TI,0,0,107),&  
! Missing in TI scheme    
ratb_t1(198,'CS2       ','O(3P)     ','COS       ','SO2       ','CO        ',&
'          ',3.20e-11,  0.00,  650.00, 0.00, 0.00, 0.00, 0.00, S+ST,A,0,107),&
ratb_t1(199,'CS2       ','OH        ','COS       ','SO2       ','          ',&
'          ',1.25e-16,  0.00,-4550.00, 0.00, 0.00, 0.00, 0.00, S+ST+R,A,0,107),&
ratb_t1(199,'CS2       ','OH        ','SO2       ','COS       ','          ',&
'          ',8.80e-16,  0.00,  -2300.00, 0.00, 0.00, 0.00, 0.00, TI,A,0,107),&
!! Commenting out existing DMS scheme to test the gas-phase DMS scheme of Chen et al. (2018)
!! LER, UC, Feb 2019
!ratb_t1(200,'DMS       ','OH        ','SO2       ','          ','          ',&
!'          ',1.20e-11,  0.00,  260.00, 0.00, 0.00, 0.00, 0.00, S+ST+R,A,0,107),&
!! Fractional production of MeOO can't be 0
!ratb_t1(200,'DMS       ','OH        ','SO2       ','DMSO      ','MeOO      ',&
!'          ',3.04e-12,  0.00, -350.00, 0.60, 0.40, 0.00, 0.00, TI,A,0,107),&
!ratb_t1(200,'DMS       ','OH        ','SO2       ','          ','          ',&
!'          ',9.60e-12,  0.00,  240.00, 0.00, 0.00, 0.00, 0.00, OL,A,0,107),&
!!ratb_t1(2001,'DMS       ','OH        ','SO2       ','MeOO      ','HCHO      ',&Commenting to keep code blocks same
!!'          ',9.60e-12,  0.00,  240.00, 0.00, 0.00, 0.00, 0.00, TI,A,0,107),&LER UC 19.02.19
!! Difference in fractional products. Inconsistent?
!ratb_t1(201,'DMS       ','OH        ','MSA       ','SO2       ','          ',&
!'          ',3.04e-12,  0.00, -350.00, 0.00, 0.00, 0.00, 0.00, S+ST+R,A,0,107),&
!!ratb_t1(201,'DMS       ','OH        ','SO2       ','DMSO      ','          ',&Commenting to keep code blocks same
!!'          ',3.04e-12,  0.00, -350.00, 0.60, 0.40, 0.00, 0.00, OL,A,0,107),&LER UC 19.02.19
!! Incomplete products
!ratb_t1(202,'DMS       ','NO3       ','SO2       ','          ','          ',&
!'          ',1.90e-13,  0.00, -500.00, 0.0, 0.0, 0.0, 0.0, S+ST+OL+R,A,0,107),&
!ratb_t1(202,'DMS       ','NO3       ','SO2       ','HONO2     ','MeOO      ',&
!'HCHO      ',1.90e-13,  0.00, -500.00, 0.00, 0.00, 0.00, 0.00, TI,A,0,107),&
!! Not in TIA scheme
!ratb_t1(203,'DMS       ','O(3P)     ','SO2       ','          ','          ',&
!'          ',1.30e-11,  0.00, -410.00, 0.00, 0.00, 0.00, 0.00, S+ST,A,0,107),&
!ratb_t1(204,'DMSO      ','OH        ','SO2       ','          ','          ',&
!'          ',5.80e-11,  0.00,    0.00, 0.60, 0.00, 0.00, 0.00, OL+R,A,0,107),&
!ratb_t1(204,'DMSO      ','OH        ','SO2       ','MSA       ','          ',&
!'          ',5.80e-11,  0.00,    0.00, 0.60, 0.40, 0.00, 0.00, TI,A,0,107),&
!Scheme of Chen et al. (2018), ACP, 13617-13637, follows: (LER, UC, Feb 2019)**
!ratb_t1(200,'DMS       ','OH        ','SO2       ','DMSO      ','MeOO      ',&
!'          ',3.04e-12,  0.00, -350.00, 0.60, 0.40, 0.00, 0.00, TI+ST,A,0,107),&
!ratb_t1(201,'DMS       ','OH        ','SO2       ','MeOO      ','HCHO      ',&
!'          ',4.69e-12,  0.00, -280.00, 0.00, 0.00, 0.00, 0.00, ST,A,0,107),&
!ratb_t1(202,'DMS       ','NO3       ','SO2       ','HONO2     ','MeOO      ',&
!'HCHO      ',1.13e-12,  0.00,  530.00, 0.00, 0.00, 0.00, 0.00, ST,A,0,107),&
!ratb_t1(203,'DMS       ','BrO       ','DMSO      ','Br        ','          ',&
!'          ',3.39e-13,  0.00,  950.00, 0.00, 0.00, 0.00, 0.00, ST,A,0,107),&
!ratb_t1(204,'DMS       ','O3        ','SO2       ','          ','          ',&
!'          ',1.00e-19,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, ST,A,0,107),&
!ratb_t1(205,'DMS       ','Cl        ','SO2       ','DMSO      ','HCl       ',&
!'ClO       ',3.40e-10,  0.00,    0.00, 0.50, 0.50, 0.50, 0.50, ST,A,0,107),&
!ratb_t1(206,'DMSO      ','OH        ','MSIA      ','SO2       ','          ',&
!'          ',8.94e-11,  0.00, 800.00, 0.95, 0.05, 0.00, 0.00, ST,A,0,107),&
!ratb_t1(207,'MSIA      ','OH        ','SO2       ','MSA       ','          ',&
!'          ',9.00e-11,  0.00,    0.00, 0.90, 0.10, 0.00, 0.00, ST,A,0,107),&
!ratb_t1(208,'MSIA      ','O3        ','MSA       ','          ','          ',&
!'          ',2.00e-18,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, ST,A,0,107),& 
!Scheme of Chen et al. (2018), with correct E/R. LER, UC, June 2019
ratb_t1(200,'DMS       ','OH        ','SO2       ','DMSO      ','MeOO      ',&! # need to confirm (but assume no change)
'          ',3.04e-12,  0.00, -350.00, 0.60, 0.40, 0.00, 0.00, TI+ST,A,0,107),&
ratb_t1(201,'DMS       ','OH        ','SO2       ','MeOO      ','HCHO      ',&
'          ',1.2e-11,  0.00,  280.00, 0.00, 0.00, 0.00, 0.00, ST,A,0,107),&! # need to confirm (4.69e-12 to 1.2e-11)
ratb_t1(202,'DMS       ','NO3       ','SO2       ','HONO2     ','MeOO      ',&
'HCHO      ',1.9e-13,  0.00, -530.00, 0.00, 0.00, 0.00, 0.00, ST,A,0,107),& !# YAB Changed from 1.13e-12 to 1.9e-12
ratb_t1(203,'DMS       ','BrO       ','DMSO      ','Br        ','          ',&
'          ',1.4e-14,  0.00, -950.00, 0.00, 0.00, 0.00, 0.00, ST,A,0,107),& !# YAB changed from 3.39e-13 to 1.4e-14
ratb_t1(204,'DMS       ','O3        ','SO2       ','          ','          ',&
'          ',1.00e-19,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, ST,A,0,107),&
ratb_t1(205,'DMS       ','Cl        ','SO2       ','DMSO      ','HCl       ',&
'ClO       ',3.40e-10,  0.00,    0.00, 0.50, 0.50, 0.50, 0.50, ST,A,0,107),&
ratb_t1(206,'DMSO      ','OH        ','MSIA      ','SO2       ','          ',&
'          ',6.1e-12,  0.00, -800.00, 0.95, 0.05, 0.00, 0.00, ST,A,0,107),& !# YAB changed from 8.94e-11 to  6.1e-12 
ratb_t1(207,'MSIA      ','OH        ','SO2       ','MSA       ','          ',&
'          ',9.00e-11,  0.00,    0.00, 0.90, 0.10, 0.00, 0.00, ST,A,0,107),&
ratb_t1(208,'MSIA      ','O3        ','MSA       ','          ','          ',&
'          ',2.00e-18,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, ST,A,0,107),&
!!^End of Chen et al. (2018) bimolecular reactions. LER, UC, Feb 2019
! Not in TIA scheme. Incomplete products
ratb_t1(209,'H2S       ','O(3P)     ','OH        ','SO2       ','          ',&
'          ',9.20e-12,  0.00, 1800.00, 0.00, 0.00, 0.00, 0.00, S+ST,A,0,107),&
ratb_t1(210,'H2S       ','OH        ','SO2       ','H2O       ','          ',&
'          ',6.00e-12,  0.00,   75.00, 0.00, 0.00, 0.00, 0.00, S+ST,A,0,107),&
ratb_t1(210,'H2S       ','OH        ','SO2       ','          ','          ',&
'          ',6.00e-12,  0.00,   75.00, 0.00, 0.00, 0.00, 0.00, TI,A,0,107),&
! Not in TIA scheme
ratb_t1(211,'COS       ','O(3P)     ','CO        ','SO2       ','          ',&
'          ',2.10e-11,  0.00, 2200.00, 0.00, 0.00, 0.00, 0.00, S+ST,A,0,107),&
ratb_t1(212,'COS       ','OH        ','CO2       ','SO2       ','          ',&
'          ',1.10e-13,  0.00, 1200.00, 0.00, 0.00, 0.00, 0.00, S+ST,A,0,107),&
ratb_t1(212,'COS       ','OH        ','SO2       ','          ','          ',&
'          ',1.10e-13,  0.00, 1200.00, 0.00, 0.00, 0.00, 0.00, TI,A,0,107),&
! Not in TI scheme
ratb_t1(213,'SO2       ','O3        ','SO3       ','          ','          ',&
'          ',3.00e-12,  0.00, 7000.00, 0.00, 0.00, 0.00, 0.00, S+ST,A,0,107),&
! Not in TI scheme
ratb_t1(214,'SO3       ','H2O       ','H2SO4     ','H2O       ','          ',&
'          ',8.50e-41,  0.00,-6540.00, 0.00, 0.00, 0.00, 0.00, S+ST,A,0,107),&
ratb_t1(215,'Monoterp  ','OH        ','Sec_Org   ','          ','          ',&
'          ',1.20e-11,  0.00, -444.0, 0.13, 0.0, 0.0, 0.0, TI+ST+OL+R,A,0,107),&
ratb_t1(216,'Monoterp  ','O3        ','Sec_Org   ','          ','          ',&
'          ',1.01e-15,  0.00,  732.00, 0.13, 0.0, 0.0, 0.0, TI+ST+OL,A,0,107),&
ratb_t1(217,'Monoterp  ','NO3       ','Sec_Org   ','          ','          ',&
'          ',1.19e-12,  0.00, -925.00, 0.13, 0.0, 0.0, 0.0, TI+ST+OL,A,0,107),&
! Make sure the below five reactions are in sync with the corresponding 
! reactions listed above. 
ratb_t1(218,'HO2S      ','O3S       ','HO2       ','O2        ','          ',&
'          ',2.03e-16,  4.57, -693.00, 0.00, 0.00, 0.00, 0.00, T,0,0,107),&
ratb_t1(219,'OHS       ','O3S       ','OH        ','O2        ','          ',&
'          ',1.70e-12,  0.00,  940.00, 0.00, 0.00, 0.00, 0.00, T,0,0,107),&
ratb_t1(220,'O(1D)S    ','H2O       ','H2O       ','          ','          ',&
'          ',2.20e-10,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, T,0,0,107),&
ratb_t1(221,'O(1D)S    ','N2        ','O(3P)S    ','N2        ','          ',&
'          ',2.10e-11,  0.00, -115.00, 0.00, 0.00, 0.00, 0.00, T,0,0,107),&
ratb_t1(222,'O(1D)S    ','O2        ','O(3P)S    ','O2        ','          ',&
'          ',3.20e-11,  0.00,  -67.00, 0.00, 0.00, 0.00, 0.00, T,0,0,107),&
! These are reactions that are in R but not in the other mechs. Should they be?
ratb_t1(223,'NO3       ','NO3       ','NO2       ','NO2       ','O2        ',&
'          ',8.50e-13,  0.00, 2450.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(224,'OH        ','NH3       ','NH2       ','H2O       ','          ',&
'          ',3.50e-12,  0.00,  925.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(225,'NO3       ','HO2       ','HONO2     ','O2        ','          ',&
'          ',4.20e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(226,'EtOO      ','MeOO      ','MeCHO     ','HO2       ','HO2       ',&
'HCHO      ',2.00e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),& 
ratb_t1(227,'C4H10     ','OH        ','s-BuOO    ','H2O       ','          ',&
'          ',7.90e-13,  2.00, -300.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(228,'s-BuOO    ','NO        ','MEK       ','HO2       ','NO2       ',&
'          ',2.54e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(229,'s-BuOO    ','MeOO      ','MEK       ','HO2       ','HO2       ',&
'HCHO      ',2.50e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(230,'OH        ','MEK       ','MEKO2     ','H2O       ','          ',&
'          ',1.30e-12,  0.00,   25.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
! don't know how to map the products onto the RAC species
ratb_t1(231,'EtOO      ','EtOO      ','EtO       ','EtO       ','O2        ',&
'          ',6.40e-14,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(232,'MeCO3     ','MeCO3     ','MeOO      ','MeOO      ','CO2       ',&
'CO2       ',2.90e-12,  0.00, -500.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(233,'MeOO      ','MeCOCH2OO ','MeCO3     ','HCHO      ','HCHO      ',&
'HO2       ',3.80e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107) /)

ratb_defs_master(351:n_bimol_master) = (/&
ratb_t1(234,'i-PrOO    ','MeOO      ','Me2CO     ','HCHO      ','HO2       ',&
'HO2       ',4.00e-14,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(235,'HO2       ','s-BuOO    ','s-BuOOH   ','O2        ','          ',&
'          ',1.82e-13,  0.00,  -1300.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(236,'OH        ','s-BuOOH   ','MEK       ','OH        ','          ',&
'          ',2.15e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
! how does MeCOCH2OOMe map onto the species?
ratb_t1(237,'NO        ','MeCOCHO2Me','MeCO3     ','MeCHO     ','NO2       ',&
'          ',2.54e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(238,'MeOO      ','MeCOCHO2Me','HCHO      ','HO2       ','MeCO3     ',&
'MeCHO     ',8.80e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(239,'HOC2H4O2  ','NO        ','HCHO      ','HCHO      ','HO2       ',&
'NO2       ',9.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(240,'MeOO      ','HOC2H4O2  ','HCHO      ','HO2       ','O2        ',&
'          ',2.00e-12,  0.00,    0.00, 3.00, 2.00, 1.00, 0.00, R,0,0,107),&
! products: HCHO, 0.47 CH2O2, 0.31 CO, 0.22 CO2, 0.31 H2O, 0.13 H2 and 0.2 HO2.
! How does CH2O2 map onto the RAQ species? If it is dropped, and CO2 and H2O 
! are dropped, this can be combined into one.
ratb_t1(241,'O3        ','C2H4      ','HCHO      ','CH2O2     ','CO        ',&
'CO2       ',6.00e-15,  0.00, 2630.00, 2.00, 0.94, 0.62, 0.44, R,0,0,107),&
ratb_t1(242,'O3        ','C2H4      ','H2O       ','H2        ','HO2       ',&
'CO2       ',6.00e-15,  0.00, 2630.00, 0.62, 0.26, 0.40, 0.00, R,0,0,107),&
ratb_t1(243,'O3        ','C3H6      ','HCHO      ','CH4       ','CO        ',&
'CO2       ',1.375e-15, 0.00, 1878.00, 2.00, 0.60, 0.80, 1.20, R,0,0,107),&
ratb_t1(244,'O3        ','C3H6      ','OH        ','MeOH      ','HO2       ',&
'MeOO      ',1.375e-15, 0.00, 1878.00, 0.56, 0.24, 0.60, 1.16, R,0,0,107),&
! products 0.420 CO2 and 0.580 H2O are dropped here
ratb_t1(245,'O3        ','C3H6      ','MeCHO     ','H2        ','CO        ',&
'HO2       ',2.75e-15,  0.00, 1878.00, 1.00, 0.24, 1.16, 0.18, R,0,0,107),&
ratb_t1(246,'HOC3H6O2  ','NO        ','HCHO      ','HO2       ','MeCHO     ',&  
'NO2       ',2.54e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(247,'MeOO      ','HOC3H6O2  ','HCHO      ','HO2       ','MeCHO     ',&
'O2        ',6.00e-13,  0.00,    0.00, 2.00, 2.00, 1.00, 1.00, R,0,0,107),&
ratb_t1(248,'O3        ','MVK       ','MGLY      ','CO        ','CH2O2      ',&
'HO2       ',3.78e-15,  0.00,  -1521.00, 2.00, 1.52, 0.48, 0.72, R,0,0,107),&
ratb_t1(249,'O3        ','MVK       ','OH        ','          ','          ',&
'          ',3.78e-15,  0.00,  -1521.00, 0.72, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(250,'HOIPO2    ','MeOO      ','HO2       ','HO2       ','HCHO      ',&
'MVK       ',5.00e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(251,'HOMVKO2   ','MeOO      ','HO2       ','HO2       ','HCHO      ',&
'MGLYOX    ',2.00e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(252,'HOIPO2    ','HO2       ','ISOOH     ','O2        ','          ',&
'          ',2.45e-13,  0.00,  -1250.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(253,'HOMVKO2   ','HO2       ','MVKOOH    ','O2        ','          ',&
'          ',2.23e-13,  0.00,  -1250.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(254,'MVKOOH    ','OH        ','MGLY      ','HCHO      ','OH        ',&
'          ',5.77e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(255,'GLY       ','OH        ','HO2       ','CO        ','CO        ',&
'          ',1.14e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(256,'HOC5H8O2  ','NO        ','MVK       ','HO2       ','HCHO      ',&
'NO2       ',2.08e-12,  0.00, -180.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(257,'OH        ','MVK       ','HOMVKO2   ','H2O       ','          ',&
'          ',4.13e-12,  0.00, -452.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),& 
! don't know how to map MeCOCHO
ratb_t1(258,'HOMVKO2   ','NO        ','MGLY      ','HCHO      ','HO2       ',&
'NO2       ',2.50e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),& 
ratb_t1(259,'NO3       ','C2H6      ','EtOO      ','HONO2     ','          ',&
'          ',5.70e-12,  0.00, 4426.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),& 
ratb_t1(260,'NO3       ','C4H10     ','s-BuOO    ','HONO2     ','          ',&
'          ',2.80e-12,  0.00, 3280.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
! don't know how to map CH2(NO3)CHO. Too long for ASAD.
ratb_t1(261,'NO3       ','C2H4      ','CH2(NO3)CHO','HO2       ','          ',&
'          ',3.30e-12,  2.00, 2880.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
ratb_t1(262,'CH2(NO3)CHO','OH        ','HCHO       ','NO2       ','CO2       ',&
'          ',4.95e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
! don't know how to map CH3CH(NO3)CHO. Too long for ASAD.
ratb_t1(263,'NO3       ','C3H6      ','CH3CH(NO3)CHO','HO2       ',&
            '          ',&
'          ',4.59e-13,  0.00, 1156.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),& 
ratb_t1(264,'CH3CH(NO3)CHO','OH        ','MeCHO     ','NO2       ',&
            'CO2       ',&
'          ',5.25e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
! '(NO3)C4H6CHO' is too long for ASAD. 
ratb_t1(265,'(NO3)C4H6CHO','OH        ','MVK       ','NO2       ','          ',&
'          ',4.16e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),& 
ratb_t1(266,'OH        ','oXYLENE   ','HO2       ','MEMALD    ','MGLY      ',&
'          ',1.36e-11,  0.00,    0.00, 1.00, 0.80, 0.80, 0.00, R,0,0,107),&
ratb_t1(267,'OXYL1     ','NO2       ','ORGNIT    ','          ','          ',&
'          ',1.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),& 
ratb_t1(268,'OH        ','MEMALD    ','MEMALD1   ','          ','          ',&
'          ',5.60e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&  
! MGLY = MeCOCHO; GLY =  CHOCHO                    
ratb_t1(269,'MEMALD1   ','NO        ','HO2       ','MGLY      ','GLY       ',&
'NO2       ',2.54e-12,  0.00, -360.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&  
ratb_t1(270,'OH        ','TOLUENE   ','MEMALD    ','GLY       ','HO2       ',&
'          ',1.18e-12,  0.00, -338.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&  
ratb_t1(271,'OH        ','TOLUENE   ','TOLP1     ','          ','          ',&
'          ',3.60e-13,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&  
ratb_t1(272,'OH        ','oXYLENE   ','OXYL1     ','          ','          ',&
'          ',1.36e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&  
ratb_t1(273,'OH        ','ORGNIT    ','MEMALD    ','GLY       ','NO2       ',&
'          ',2.70e-12,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&  
ratb_t1(274,'NO3       ','ORGNIT    ','MEMALD    ','GLY       ','NO2       ',&
'NO2       ',7.00e-14,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&  
ratb_t1(275,'OXYL1     ','HO2       ','MGLY      ','MEMALD    ','          ',&
'          ',2.50e-13,  0.00,  -1300.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&  
ratb_t1(276,'MEMALD1   ','MeOO      ','HO2       ','HCHO      ','MGLY      ',&
'GLY       ',1.00e-13,  0.00,    0.00, 2.00, 1.00, 1.00, 1.00, R,0,0,107),&
ratb_t1(277,'HO2       ','TOLP1     ','MEMALD    ','GLY       ','OH        ',&
'          ',1.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&  
ratb_t1(278,'NO2       ','TOLP1     ','ORGNIT    ','          ','          ',&
'          ',1.00e-11,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,0,0,107),&
! not sure about the rate or products
ratb_t1(279,'SO2       ','OH        ','          ','          ','          ',&
'          ',0.00e+00,  0.00,    0.00, 0.00, 0.00, 0.00, 0.00, R,A,0,107),&
! The below reaction is missing in ST.
ratb_t1(280,'MACRO2    ','MeOO      ','MGLY      ','HACET     ','MeCO3     ',& 
'HCHO      ',1.00e-12,  0.00,    0.00, 1.00, 0.75, 0.25, 2.75, TI,0,0,107),&
ratb_t1(281,'MACRO2    ','MeOO      ','HO2       ','CO        ','          ',& 
'          ',1.00e-12,  0.00,    0.00, 1.17, 0.25, 0.00, 0.00, TI,0,0,107) /)

!----------------------------------------------------------------------
! NOTES: CheST Bimolecular Reactions
!----------------------------------------------------------------------
! B001 Br+Cl2O2 -> BrCl Cl O2 JPL2011
! B001 k(298)=3.3E-12. Both JPL2011 and IUPAC2006 agree. Reaction forms
! B001 ClOO rather than Cl and O2.
!----------------------------------------------------------------------
! B002 Br+HCHO -> HBr CO HO2 JPL2011
! B002 k(298)=1.1E-12. Unchanged since JPL2009. No reference by IUPAC.
! B002 JPL cite two studies of which this value is derived from by a
! B002 weighted mean.
!----------------------------------------------------------------------
! B003 Br+HO2 -> HBr O2  JPL2011
! B003 k(298)=1.7E-12. Unchanged since JPL2006. IUPAC2003 recommend
! B003 k=7.7E-12*exp(-450/temp).
!----------------------------------------------------------------------
! B004 Br+O3 -> BrO O2  JPL2011
! B004 k(298)=1.2E-12. Unchanged since JPL2009. IUPAC2003 recommend k=1
! B004 .7E-11exp(-800/temp)
!----------------------------------------------------------------------
! B005 Br+OClO -> BrO ClO  JPL2011
! B005 k(298)=3.4E-13. Originally from JPL1990
!----------------------------------------------------------------------
! B006 BrO+BrO -> Br Br O2 JPL2011
! B006 Originally from JPL1997
!----------------------------------------------------------------------
! B007 BrO+ClO -> Br Cl O2 JPL2011
! B007 k(298)=5.5E-12. Unchanged since JPL2009. IUPAC2004 recommend
! B007 k=2.9E-12*exp(220/temp) with k(298)=6.1E-12. Both JPL and IUPAC
! B007 recommended products are Br + ClOO.
!----------------------------------------------------------------------
! B008 BrO+ClO -> Br OClO  JPL2011
! B008 k(298)=6.0E-12. Unchanged since JPL2009. IUPAC2004 recommend
! B008 k=1.60E-13*exp(430/temp) with k(298)=6.8E-12 (at 220 K k =
! B008 1.13E-12 c.f 1.16E-12). This channel accounts for around 60% of
! B008 the total for this reaction.
!----------------------------------------------------------------------
! B009 BrO+ClO -> BrCl O2  JPL2011
! B009 k(298)=1.1E-12. Unchanged since JPL2009. IUPAC2004 recommend
! B009 k=5.8E-13*exp(170/temp) with k(298)=1.0E-12. Branching ratios put
! B009 this channel at ~8%.
!----------------------------------------------------------------------
! B010 BrO+HO2 -> HBr O3  JPL2011
! B010 Remove this reaction since both JPL2006 and IUPAC2003 suggest
! B010 that the branching ratio for this channel is 0.0
!----------------------------------------------------------------------
! B011 BrO+HO2 -> HOBr O2  JPL2011
! B011 k(298)=2.1E-11. Unchanged since JPL2006. IUPAC2003 recommend k=4
! B011 .5E-12exp(500/temp) with k(298)=2.4E-11. Both JPL and IUPAC
! B011 suggest that essentially this is the only channel for this
! B011 reaction.
!----------------------------------------------------------------------
! B012 BrO+NO -> Br NO2  JPL2011
! B012 k(298)=2.1E-11. Originally from JPL1982
!----------------------------------------------------------------------
! B013 BrO+OH -> Br HO2  JPL2011
! B013 k(298)=3.9E-11. Unchanged since JPL2006. IUPAC2003 recommend
! B013 k=1.8E-11*exp(250/temp) with k(298)=4.1E-11. The temperature
! B013 dependence comes from the study by Bedjanian et al.
!----------------------------------------------------------------------
! B014 CF2Cl2+O(1D) -> Cl ClO  JPL2011
! B014 k(298)=1.4E-10. Unchanged since JPL1992. Both JPL1992 and
! B014 IUPAC2011 agree. The ClO channel is ~87%.
!----------------------------------------------------------------------
! B015 CFCl3+O(1D) -> Cl Cl ClO JPL2011
! B015 k(298)=2.3E-10. Unchanged since JPL1992. Both JPL1992 and
! B015 IUPAC2011 agree. The ClO channel is ~80%.
!----------------------------------------------------------------------
! B016 Cl+CH4 -> HCl MeOO  JPL2011
! B016 k(298)=1.0E-13. IUPAC recommend k=6.6E-12*exp(-1240/temp).
!----------------------------------------------------------------------
! B017 Cl+Cl2O2 -> Cl Cl Cl IUPAC2006
! B017 k(298)=1.0E-10. Probable error in JPL2011 gives beta=9.4E-11!
! B017 IUPAC2006 recommend k = 7.6E-11*exp(65/temp) based on study by
! B017 Ingham et al.
!----------------------------------------------------------------------
! B018 Cl+ClONO2 -> Cl Cl NO3 JPL2011
! B018 IUPAC similar (6.2 and -145). Actual reaction forms Cl2 rather
! B018 than 2Cl.
!----------------------------------------------------------------------
! B019 Cl+H2 -> HCl H  JPL2011
! B019 IUPAC similar (3.9 and 2310)
!----------------------------------------------------------------------
! B020 Cl+H2O2 -> HCl HO2  JPL2011
! B020 IUPAC identical.
!----------------------------------------------------------------------
! B021 Cl+HCHO -> HCl CO HO2 JPL2011
! B021 IUPAC similar (8.1 and 34).
!----------------------------------------------------------------------
! B022 Cl+HO2 -> ClO OH  JPL2011
! B022 IUPAC gives total rate + k2/k. Hard to compare
!----------------------------------------------------------------------
! B023 Cl+HO2 -> HCl O2  JPL2011
! B023 IUPAC gives total rate + k2/k. Hard to compare
!----------------------------------------------------------------------
! B024 Cl+HOCl -> Cl Cl OH JPL2011
! B024 No IUPAC rate
!----------------------------------------------------------------------
! B025 Cl+MeOOH -> HCl MeOO  JPL2011
! B025 IUPAC rate similar (5.9E-11)
!----------------------------------------------------------------------
! B026 Cl+NO3 -> ClO NO2  JPL2011
! B026 IUPAC rates identical
!----------------------------------------------------------------------
! B027 Cl+O3 -> ClO O2  JPL2011
! B027 IUPAC rates similar (2.8 and 250)
!----------------------------------------------------------------------
! B028 Cl+OClO -> ClO ClO  JPL2011
! B028 IUPAC rates similar (3.2 and -170)
!----------------------------------------------------------------------
! B029 ClO+ClO -> Cl Cl O2 JPL2011
! B029 IUPAC rates identical
!----------------------------------------------------------------------
! B030 ClO+ClO -> Cl Cl O2 JPL2011
! B030 IUPAC rates identical
!----------------------------------------------------------------------
! B031 ClO+ClO -> Cl OClO  JPL2011
! B031 IUPAC rates identical
!----------------------------------------------------------------------
! B032 ClO+HO2 -> HCl O3  JPL2011
! B032 No rate given in IUPAC or JPL. Reaction should probably be
! B032 removed.
!----------------------------------------------------------------------
! B033 ClO+HO2 -> HOCl O2  JPL2011
! B033 IUPAC rates similar (2.2 and -340)
!----------------------------------------------------------------------
! B034 ClO+MeOO -> Cl HCHO HO2 JPL2011
! B034 IUPAC provides a more detailed breakdown of products. Would
! B034 require a modification of asad_bimol to treat properly. Plus
! B034 retain policy of using JPL for ClOx reactions.
!----------------------------------------------------------------------
! B035 ClO+NO -> Cl NO2  JPL2011
! B035 IUPAC rates very similar (6.2 and -295)
!----------------------------------------------------------------------
! B036 ClO+NO3 -> Cl O2 NO2 JPL2011
! B036 Assume all ClO + NO3 forms ClOO. Note IUPAC still has 1.2/4.6
! B036 going to other channel
!----------------------------------------------------------------------
! B037 ClO+NO3 -> OClO NO2  JPL2011
! B037 Assume all ClO + NO3 forms ClOO. Note IUPAC still has 1.2/4.6
! B037 going to other channel. If use JPL approach reaction could be
! B037 removed.
!----------------------------------------------------------------------
! B038 EtCO3+NO -> EtOO CO2 NO2 IUPAC2002
! B038 No JPL rate
!----------------------------------------------------------------------
! B039 EtCO3+NO3 -> EtOO CO2 NO2 MCMv3.2
! B039 No IUPAC/JPL rate
!----------------------------------------------------------------------
! B040 EtOO+MeCO3 -> MeCHO HO2 MeOO IUPAC2002
! B040 No JPL rate
!----------------------------------------------------------------------
! B041 EtOO+NO -> MeCHO HO2 NO2 IUPAC2005
! B041 No JPL rate
!----------------------------------------------------------------------
! B042 EtOO+NO3 -> MeCHO HO2 NO2 IUPAC2008
! B042 No JPL rate
!----------------------------------------------------------------------
! B043 H+HO2 -> H2 O2  JPL2011
! B043 No IUPAC rate
!----------------------------------------------------------------------
! B044 H+HO2 -> O(3P) H2O  JPL2011
! B044 No IUPAC rate
!----------------------------------------------------------------------
! B045 H+HO2 -> OH OH  JPL2011
! B045 No IUPAC rate
!----------------------------------------------------------------------
! B046 H+NO2 -> OH NO  JPL2011
! B046 No IUPAC rate
!----------------------------------------------------------------------
! B047 H+O3 -> OH O2  JPL2011
! B047 No IUPAC rate
!----------------------------------------------------------------------
! B048 HO2+EtCO3 -> O2 EtCO3H  MCMv3.2*
! B048 No IUPAC/JPL rate. No exact correspondence in MCM. Take total
! B048 rate and subtract O3 channel.
!----------------------------------------------------------------------
! B049 HO2+EtCO3 -> O3 EtCO2H  MVMv3.2
! B049 No IUPAC/JPL rate
!----------------------------------------------------------------------
! B050 HO2+EtOO -> EtOOH   IUPAC2011
! B050 -
!----------------------------------------------------------------------
! B051 HO2+HO2 -> H2O2   JPL2011 see also asad_bimol
! B051 IUPAC slightly different (2.2 & -600) . Rate modified in the
! B051 presence of H2O in asad bimol.
!----------------------------------------------------------------------
! B052 HO2+ISO2 -> ISOOH   Poschl00
! B052 -
!----------------------------------------------------------------------
! B053 HO2+MACRO2 -> MACROOH   Poschl00
! B053 -
!----------------------------------------------------------------------
! B054 HO2+MeCO3 -> MeCO2H O3  IUPAC2009
! B054 -
!----------------------------------------------------------------------
! B055 HO2+MeCO3 -> MeCO3H   IUPAC2009
! B055 -
!----------------------------------------------------------------------
! B056 HO2+MeCO3 -> OH MeOO  IUPAC2009
! B056 -
!----------------------------------------------------------------------
! B057 HO2+MeCOCH2OO -> MeCOCH2OOH   IUPAC2009
! B057 Use measured value rather than the expression from the MCM (which
! B057 does include T dependence)
!----------------------------------------------------------------------
! B058 HO2+MeOO -> HCHO H2O  IUPAC2009 see also asad_bimol
! B058 Total rate given. Split between channels controlled by asad_bimol
!----------------------------------------------------------------------
! B059 HO2+MeOO -> MeOOH O2  IUPAC2009 see also asad_bimol
! B059 Total rate given. Split between channels controlled by asad_bimol
!----------------------------------------------------------------------
! B060 HO2+NO -> OH NO2  JPL2011
! B060 IUPAC values similar (3.45 and 270)
!----------------------------------------------------------------------
! B061 HO2+NO3 -> OH NO2 O2 JPL2011
! B061 IUPAC value slightly larger (4 cf 3.5)
!----------------------------------------------------------------------
! B062 HO2+O3 -> OH O2  IUPAC2001
! B062 Use more sophisticated expression of IUPAC. JPL have simple term
! B062 1*10^-14 and 490.
!----------------------------------------------------------------------
! B063 HO2+i-PrOO -> i-PrOOH   MCMv3.2
! B063 No rate in IUPAC or JPL
!----------------------------------------------------------------------
! B064 HO2+n-PrOO -> n-PrOOH   MCMv3.2
! B064 No rate in IUPAC or JPL
!----------------------------------------------------------------------
! B065 ISO2+ISO2 -> MACR MACR HCHOHO2 Poschl00
! B065 -
!----------------------------------------------------------------------
! B066 MACRO2+MACRO2 -> HACET MGLY HCHOCO Poschl00
! B066 -
!----------------------------------------------------------------------
! B067 MACRO2+MACRO2 -> HO2   Poschl00
! B067 -
!----------------------------------------------------------------------
! B068 MeBr+Cl -> Br HCl  JPL2011
! B068 -
!----------------------------------------------------------------------
! B069 MeBr+O(1D) -> Br OH  JPL2011
! B069 -
!----------------------------------------------------------------------
! B070 MeBr+OH -> Br H2O  JPL2011
! B070 -
!----------------------------------------------------------------------
! B071 MeCO3+NO -> MeOO CO2 NO2 IUPAC2002
! B071 -
!----------------------------------------------------------------------
! B072 MeCO3+NO3 -> MeOO CO2 NO2 IUPAC2008
! B072 -
!----------------------------------------------------------------------
! B073 MeCOCH2OO+NO -> MeCO3 HCHO NO2 MCMv3.2
! B073 No rate in IUPAC or JPL
!----------------------------------------------------------------------
! B074 MeCOCH2OO+NO3 -> MeCO3 HCHO NO2 MCMv3.2
! B074 No rate in IUPAC or JPL
!----------------------------------------------------------------------
! B075 MeOO+MeCO3 -> HO2 HCHO MeOO IUPAC2002
! B075 JPL only gives total rates
!----------------------------------------------------------------------
! B076 MeOO+MeCO3 -> MeCO2H HCHO  IUPAC2002
! B076 JPL only gives total rates
!----------------------------------------------------------------------
! B077 MeOO+MeOO -> HO2 HO2 HCHOHCHO IUPAC2002 see also asad_bimol
! B077 IUPAC gives total k and k2 as f(T) resulting in complicated
! B077 expression evaluated in asad_bimol. JPL only gives a total rate
! B077 (9.5 and -390).
!----------------------------------------------------------------------
! B078 MeOO+MeOO -> MeOH HCHO  IUPAC2002 see also asad_bimol
! B078 IUPAC gives total k and k2 as f(T) resulting in complicated
! B078 expression evaluated in asad_bimol. JPL only gives a total rate
! B078 (9.5 and -390).
!----------------------------------------------------------------------
! B079 MeOO+NO -> HO2 HCHO NO2 IUPAC2005
! B079 IUPAC gives rates to CH3O. Assume 99.9% of the CH3O goes to HCHO
! B079 and HO2 and 0.1% goes to MeONO2. Where does this partition come
! B079 from?
!----------------------------------------------------------------------
! B080 MeOO+NO -> MeONO2   IUPAC2005
! B080 IUPAC gives rates to CH3O. Assume 99.9% of the CH3O goes to HCHO
! B080 and HO2 and 0.1% goes to MeONO2. Where does this partition come
! B080 from?
!----------------------------------------------------------------------
! B081 MeOO+NO3 -> HO2 HCHO NO2 MCMv3.2
! B081 No JPL value.
!----------------------------------------------------------------------
! B082 N+NO -> N2 O(3P)  JPL2011
! B082 No IUPAC value.
!----------------------------------------------------------------------
! B083 N+NO2 -> N2O O(3P)  JPL2011
! B083 No IUPAC value.
!----------------------------------------------------------------------
! B084 N+O2 -> NO O(3P)  JPL2011
! B084 IUPAC  values slightly different (5.1 and -198). Different
! B084 choices of results and fitting.
!----------------------------------------------------------------------
! B085 N2O5+H2O -> HONO2 HONO2  ????
! B085 Where is this from? Both JPL and IUPAC only give upper limits
! B085 (which are lower than this value)
!----------------------------------------------------------------------
! B086 NO+ISO2 -> ISON   Poschl00
! B086 -
!----------------------------------------------------------------------
! B087 NO+ISO2 -> NO2 MACR HCHOHO2 Poschl00
! B087 -
!----------------------------------------------------------------------
! B088 NO+MACRO2 -> MGLY HCHO HO2 Poschl00
! B088 Note to accommodate all the products split into two reactions and
! B088 half the rates.
!----------------------------------------------------------------------
! B089 NO+MACRO2 -> NO2 MeCO3 HACETCO Poschl00
! B089 Note to accommodate all the products split into two reactions and
! B089 half the rates.
!----------------------------------------------------------------------
! B090 NO+NO3 -> NO2 NO2  JPL2011
! B090 Discrepancy with IUPAC (1.8E-11 & beta =170) related to
! B090 differences in averaging.
!----------------------------------------------------------------------
! B091 NO+O3 -> NO2   JPL2011
! B091 Considerable discrepancy with IUPAC (1.4E-12 & beta=1310)
! B091 resulting in differences ~10% between the two. Use JPL as include
! B091 the extra (room temperature only though) of Stedman & Niki and
! B091 Bemand et al. Believe difference arises from fitting & averaging
! B091 vs averaging & fitting.
!----------------------------------------------------------------------
! B092 NO2+NO3 -> NO NO2 O2 JPL2011
! B092 This is a termolecular reaction. Why is it here? IUPAC don't
! B092 recognise the reaction and JPL note its existence is not
! B092 definitively proven. However they do use studies that infer the
! B092 rates from thermal decomposition of N2O5.
!----------------------------------------------------------------------
! B093 NO2+O3 -> NO3   JPL2011
! B093 Slight disagreement with IUPAC rates (1.4E-11 & -2470). Use JPL
! B093 by default.
!----------------------------------------------------------------------
! B094 NO3+Br -> BrO NO2  JPL2011
! B094 Both IUPAC & JPL recommend approach of Mellouki et al (1989)
!----------------------------------------------------------------------
! B095 NO3+C5H8 -> ISON   IUPAC2007
! B095 -
!----------------------------------------------------------------------
! B096 NO3+EtCHO -> HONO2 EtCO3  IUPAC2007
! B096 No JPL rates. No T dependence given buy IUPAC.
!----------------------------------------------------------------------
! B097 NO3+HCHO -> HONO2 HO2 CO IUPAC2007
! B097 No direct measurements of T dependence. Infer from T dependence
! B097 of MeCHO + NO3.  JPL values at 298K same as IUPAC.
!----------------------------------------------------------------------
! B098 NO3+MGLY -> MeCO3 CO HONO2 MCMv3.2
! B098 No rate in JPL/IUPAC. 2.4*IUPAC MeCHO + NO3.
!----------------------------------------------------------------------
! B099 NO3+Me2CO -> HONO2 MeCOCH2OO  IUPAC2007
! B099 Number is actually an upper limit. Use it as no other number
! B099 available (JPL don't provide one either)
!----------------------------------------------------------------------
! B100 NO3+MeCHO -> HONO2 MeCO3  IUPAC2007
! B100 -
!----------------------------------------------------------------------
! B101 O(1D)+CH4 -> HCHO H2  JPL2011
! B101 IUPAC values slightly lower (7.50E-12)
!----------------------------------------------------------------------
! B102 O(1D)+CH4 -> HCHO HO2 HO2 JPL2011
! B102 IUPAC values same
!----------------------------------------------------------------------
! B103 O(1D)+CH4 -> OH MeOO  JPL2011
! B103 IUPAC values slightly lower (1.05E-10)
!----------------------------------------------------------------------
! B104 O(1D)+CO2 -> O(3P) CO2  JPL2011
! B104 -
!----------------------------------------------------------------------
! B105 O(1D)+H2 -> OH H  IUPAC2008
! B105 JPL values identical
!----------------------------------------------------------------------
! B106 O(1D)+H2O -> OH OH  JPL2011
! B106 JPL value dominated by the work of Dunlea & Ravishankra. Slightly
! B106 lower than the IUPAC values (5% at room temperature)
!----------------------------------------------------------------------
! B107 O(1D)+HBr -> HBr O(3P)  JPL2011
! B107 Use Wine al for BR to O3P
!----------------------------------------------------------------------
! B108 O(1D)+HBr -> OH Br  JPL2011
! B108 Use Wine et al O3P BR and assume rest of HBR gos to OH + Br
!----------------------------------------------------------------------
! B109 O(1D)+HCl -> H ClO  JPL2011
! B109 Use Wine et al for branching ratios.
!----------------------------------------------------------------------
! B110 O(1D)+HCl -> O(3P) HCl  JPL2011
! B110 Use Wine et al for branching ratios.
!----------------------------------------------------------------------
! B111 O(1D)+HCl -> OH Cl  JPL2011
! B111 Use Wine et al for branching ratios.
!----------------------------------------------------------------------
! B112 O(1D)+N2 -> O(3P) N2  JPL2011
! B112 IUPAC identical
!----------------------------------------------------------------------
! B113 O(1D)+N2O -> N2 O2  JPL2011
! B113 Use JPL values (slightly higher than IUPAC 4.3e-11) as slightly
! B113 newer and include T dependence
!----------------------------------------------------------------------
! B114 O(1D)+N2O -> NO NO  JPL2011
! B114 Use JPL values (slightly lower than IUPAC 7.6e-11) as slightly
! B114 newer and include T dependence
!----------------------------------------------------------------------
! B115 O(1D)+O2 -> O(3P) O2  JPL2011
! B115 IUPAC values slightly different (3.2 and -67)
!----------------------------------------------------------------------
! B116 O(1D)+O3 -> O2 O(3P) O(3P) JPL2011
! B116 IUPAC values identical
!----------------------------------------------------------------------
! B117 O(1D)+O3 -> O2 O2  JPL2011
! B117 IUPAC values identical
!----------------------------------------------------------------------
! B118 O(3P)+BrO -> O2 Br  JPL2011
! B118 -
!----------------------------------------------------------------------
! B119 O(3P)+ClO -> Cl O2  JPL2011
! B119 Updated JPL 10-6
!----------------------------------------------------------------------
! B120 O(3P)+ClONO2 -> ClO NO3  JPL2011
! B120 Updated JPL 10-6
!----------------------------------------------------------------------
! B121 O(3P)+H2 -> OH H  ????
! B121 Not defined in IUPAC or JPL. Where does this come from?
!----------------------------------------------------------------------
! B122 O(3P)+H2O2 -> OH HO2  JPL2011
! B122 IUPAC values identical
!----------------------------------------------------------------------
! B123 O(3P)+HBr -> OH Br  JPL2011
! B123 -
!----------------------------------------------------------------------
! B124 O(3P)+HCHO -> OH CO HO2 JPL2011
! B124 Value not given in IUPAC (are we sure?)
!----------------------------------------------------------------------
! B125 O(3P)+HCl -> OH Cl  JPL2011
! B125 -
!----------------------------------------------------------------------
! B126 O(3P)+HO2 -> OH O2  IUPAC2001
! B126 JPL are the same to one significant figure.
!----------------------------------------------------------------------
! B127 O(3P)+HOCl -> OH ClO  JPL2011
! B127 -
!----------------------------------------------------------------------
! B128 O(3P)+NO2 -> NO O2  JPL2011
! B128 Note T dependence in IUPAC is slightly different (198 cf 210)
!----------------------------------------------------------------------
! B129 O(3P)+NO3 -> O2 NO2  IUPAC2009
! B129 Rate smaller in JPL2011 (1.1e-11) who prefer the results of
! B129 Graham & Johnston (1978) to Canosa-Mas et al. (1989). Though
! B129 IUPAC note that the results agree within uncertainties.
!----------------------------------------------------------------------
! B130 O(3P)+O3 -> O2 O2  JPL2011
! B130 IUPAC identical
!----------------------------------------------------------------------
! B131 O(3P)+OClO -> O2 ClO  JPL2011
! B131 -
!----------------------------------------------------------------------
! B132 O(3P)+OH -> O2 H  JPL2011
! B132 IUPAC slightly different (2.4 & -110). JPL include slightly more
! B132 results and newer results. However main difference probably
! B132 derives from fitting.
!----------------------------------------------------------------------
! B133 O3+C5H8 -> HO2 OH  IUPAC2007*
! B133 Total rates from IUPAC. Arbitrarily split into three reactions.
! B133 Fractional products given by Poschl et al (2000).
!----------------------------------------------------------------------
! B134 O3+C5H8 -> MACR HCHO MACRO2MeCO3 IUPAC2007*
! B134 Total rates from IUPAC. Arbitrarily split into three reactions.
! B134 Fractional products given by Poschl et al (2000).
!----------------------------------------------------------------------
! B135 O3+C5H8 -> MeOO HCOOH COH2O2 IUPAC2007*
! B135 Total rates from IUPAC. Arbitrarily split into three reactions.
! B135 Fractional products given by Poschl et al (2000).
!----------------------------------------------------------------------
! B136 O3+MACR -> MGLY HCOOH HO2CO IUPAC2007*
! B136 Complicated expression. Rate is average of IUPAC MACR + O3 and
! B136 MVK + O3. Split into further two reactions to allow all products
! B136 to be included. Hence rates are IUPAC/4. Fractional Products from
! B136 Poschl et al (2000).
!----------------------------------------------------------------------
! B137 O3+MACR -> MGLY HCOOH HO2CO IUPAC2007*
! B137 Complicated expression. Rate is average of IUPAC MACR + O3 and
! B137 MVK + O3. Split into further two reactions to allow all products
! B137 to be included. Hence rates are IUPAC/4. Fractional Products from
! B137 Poschl et al (2000).
!----------------------------------------------------------------------
! B138 O3+MACR -> OH MeCO3  IUPAC2007*
! B138 Complicated expression. Rate is average of IUPAC MACR + O3 and
! B138 MVK + O3. Split into further two reactions to allow all products
! B138 to be included. Hence rates are IUPAC/4. Fractional Products from
! B138 Poschl et al (2000).
!----------------------------------------------------------------------
! B139 O3+MACR -> OH MeCO3  IUPAC2007*
! B139 Complicated expression. Rate is average of IUPAC MACR + O3 and
! B139 MVK + O3. Split into further two reactions to allow all products
! B139 to be included. Hence rates are IUPAC/4. Fractional Products from
! B139 Poschl et al (2000).
!----------------------------------------------------------------------
! B140 OClO+NO -> NO2 ClO  JPL2011
! B140 Introduce T dependence of Bemand et al. New measurements by Li et
! B140 al in good agreement with these
!----------------------------------------------------------------------
! B141 OH+C2H6 -> H2O EtOO  IUPAC2007
! B141 JPL slightly different (7.22 & 1020). Difference arises from
! B141 fitting. Use IUPAC by default for organics.
!----------------------------------------------------------------------
! B142 OH+C3H8 -> i-PrOO H2O  IUPAC2007 see also asad_bimol
! B142 Total k rate (given) from IUPAC2007. Split between the two
! B142 channels is temperature dependent using measurements of Droege
! B142 and Tully (1986) and covered by asad_bimol. Note JPL values
! B142 slightly higher
!----------------------------------------------------------------------
! B143 OH+C3H8 -> n-PrOO H2O  IUPAC2007 see also asad_bimol
! B143 Total k rate (given) from IUPAC2007. Split between the two
! B143 channels is temperature dependent using measurements of Droege
! B143 and Tully (1986) and covered by asad_bimol. Note JPL values
! B143 slightly higher
!----------------------------------------------------------------------
! B144 OH+C5H8 -> ISO2   IUPAC2009
! B144 No JPL rate.
!----------------------------------------------------------------------
! B145 OH+CH4 -> H2O MeOO  JPL2011
! B145 Small difference with IUPAC results (1.85 and 1690)
!----------------------------------------------------------------------
! B146 OH+CO -> HO2   IUPAC2005 see also asad_bimol
! B146 Base rate modified by density dependence term in asad_bimol. JPL
! B146 treat as a termolecular reaction.
!----------------------------------------------------------------------
! B147 OH+ClO -> HCl O2  JPL2011
! B147 -
!----------------------------------------------------------------------
! B148 OH+ClO -> HO2 Cl  JPL2011
! B148 -
!----------------------------------------------------------------------
! B149 OH+ClONO2 -> HOCl NO3  JPL2011
! B149 -
!----------------------------------------------------------------------
! B150 OH+EtCHO -> H2O EtCO3  IUPAC2007
! B150 Latest value includes measurements of Le Crane et al (2005)>
!----------------------------------------------------------------------
! B151 OH+EtOOH -> H2O EtOO  MCMv3.2
! B151 IUPAC provides a limit
!----------------------------------------------------------------------
! B152 OH+EtOOH -> H2O MeCHO OH MCMv3.2
! B152 IUPAC provides a limit
!----------------------------------------------------------------------
! B153 OH+H2 -> H2O HO2  JPL2011
! B153 IUPAC different parameterisation (7.7 and 2100). Same at 298
!----------------------------------------------------------------------
! B154 OH+H2O2 -> H2O HO2  IUPAC2001
! B154 -
!----------------------------------------------------------------------
! B155 OH+HACET -> MGLY HO2  IUPAC2007
! B155 Includes temperature dependence of Dillon et al (2006)
!----------------------------------------------------------------------
! B156 OH+HBr -> H2O Br  JPL2011
! B156 -
!----------------------------------------------------------------------
! B157 OH+HCHO -> H2O HO2 CO IUPAC2007
! B157 -
!----------------------------------------------------------------------
! B158 OH+HCOOH -> HO2   IUPAC2007
! B158 -
!----------------------------------------------------------------------
! B159 OH+HCl -> H2O Cl  JPL2011
! B159 -
!----------------------------------------------------------------------
! B160 OH+HO2 -> H2O   JPL2011
! B160 IUPAC identical
!----------------------------------------------------------------------
! B161 OH+HO2NO2 -> H2O NO2  IUPAC2007
! B161 The preferred values are based on the recent and extensive
! B161 absolute rate study of JimÃÂ©nez et al. (2004) i
!----------------------------------------------------------------------
! B162 OH+HOCl -> ClO H2O  JPL2011
! B162 -
!----------------------------------------------------------------------
! B163 OH+HONO -> H2O NO2  IUPAC2004
! B163 IUPAC includes the more modern results of Burkholder. JPL values
! B163 slightly different (1.8 and 390)
!----------------------------------------------------------------------
! B164 OH+HONO2 -> H2O NO3  IUPAC2004 see also asad_bimol
! B164 Include rate with no density dependence. Density dependence
! B164 calculated using asad_bimol
!----------------------------------------------------------------------
! B165 OH+ISON -> HACET NALD  Poschl00
! B165 Use original MIM rate of Poschl et al (2000)
!----------------------------------------------------------------------
! B166 OH+ISOOH -> MACR OH  Poschl00
! B166 Use original MIM rate of Poschl et al (2000)
!----------------------------------------------------------------------
! B167 OH+MACR -> MACRO2   IUPAC2007
! B167 MACR rate is average of MACR+OH and MVK+OH. This is the IUPAC MVK
! B167 rate.
!----------------------------------------------------------------------
! B168 OH+MACR -> MACRO2   IUPAC2007
! B168 MACR rate is average of MACR+OH and MVK+OH. This is the IUPAC
! B168 MACR rate.
!----------------------------------------------------------------------
! B169 OH+MACROOH -> MACRO2   MCMv3.2
! B169  No data in IUPAC or JPL. Rate k is based on MCM v3.2
!----------------------------------------------------------------------
! B170 OH+MGLY -> MeCO3 CO  IUPAC2008
! B170 -
!----------------------------------------------------------------------
! B171 OH+MPAN -> HACET NO2  IUPAC2006
! B171 -
!----------------------------------------------------------------------
! B172 OH+Me2CO -> H2O MeCOCH2OO  IUPAC2007
! B172 Two part reaction due to IUPAC recommendation of the form k1+k2
!----------------------------------------------------------------------
! B173 OH+Me2CO -> H2O MeCOCH2OO  IUPAC2007
! B173  Two part reaction due to IUPAC recommendation of the form k1+k2
!----------------------------------------------------------------------
! B174 OH+MeCHO -> H2O MeCO3  IUPAC2009
! B174 Use total cross section.
!----------------------------------------------------------------------
! B175 OH+MeCO2H -> MeOO   MCMv3.2
! B175 No data in IUPAC or JPL. Rate k is based on MCM v3.2. Can't find
! B175 the data used for previous versions of the code.
!----------------------------------------------------------------------
! B176 OH+MeCO3H -> MeCO3   MCMv3.2
! B176 No data in IUPAC or JPL. Rate k is based on MCM v3.2
!----------------------------------------------------------------------
! B177 OH+MeCOCH2OOH -> H2O MeCOCH2OO  MCMv3.2
! B177  No data in IUPAC or JPL. Rate k is based on MCM v3.2
!----------------------------------------------------------------------
! B178 OH+MeCOCH2OOH -> OH MGLY  MCMv3.2
! B178  No data in IUPAC or JPL. Rate k is based on MCM v3.2
!----------------------------------------------------------------------
! B179 OH+MeOH -> HO2 HCHO  IUPAC2006
! B179 JPL rates same.
!----------------------------------------------------------------------
! B180 OH+MeONO2 -> HCHO NO2 H2O IUPAC2006
! B180 JPL rates slightly different 8.0E-13 and 1000.  Large
! B180 uncertainties on measurement.
!----------------------------------------------------------------------
! B181 OH+MeOOH -> H2O HCHO OH IUPAC2007
! B181 Use total rate from IUPAC * 0.4 (IUPAC BR to HCHO). Note JPL rate
! B181 is considerably lower. There are discrepancies between measured
! B181 values of this cross section which JPL and IUPAC reconcile
! B181 differently. We employ the IUPAC values.
!----------------------------------------------------------------------
! B182 OH+MeOOH -> H2O MeOO  IUPAC2007
! B182 Use total rate from IUPAC * 0.6 (IUPAC BR to MeOO). Note JPL rate
! B182 is considerably lower.
!----------------------------------------------------------------------
! B183 OH+NALD -> HCHO CO NO2 IUPAC2009
! B183 As in Poschl et al use MeCHO + OH rates (but updated).
!----------------------------------------------------------------------
! B184 OH+NO3 -> HO2 NO2  JPL2011
! B184 Note small discrepancy with IUPAC values (JPL rates 10% higher).
! B184 Same measurements used
!----------------------------------------------------------------------
! B185 OH+O3 -> HO2 O2  JPL2011
! B185 Same values in IUPAC.
!----------------------------------------------------------------------
! B186 OH+OClO -> HOCl O2  JPL2011
! B186 Updated based on the recommended value reported by Gierczak et
! B186 al.
!----------------------------------------------------------------------
! B187 OH+OH -> H2O O(3P)  IUPAC2001
! B187 No T dependence given by JPL
!----------------------------------------------------------------------
! B188 OH+PAN -> HCHO NO2 H2O MCMv3.2
! B188 No data in IUPAC or JPL. Rate k is based on MCM v3.2
!----------------------------------------------------------------------
! B189 OH+PPAN -> MeCHO NO2 H2O MCMv3.2
! B189 No data in IUPAC or JPL. Rate k is based on MCM v3.2
!----------------------------------------------------------------------
! B190 OH+i-PrOOH -> Me2CO OH  MCMv3.2
! B190 No data in IUPAC or JPL. Rate k is based on MCM v3.2
!----------------------------------------------------------------------
! B191 OH+i-PrOOH -> i-PrOO H2O  MCMv3.2
! B191 No data in IUPAC or JPL. Rate k is based on MCM v3.2
!----------------------------------------------------------------------
! B192 OH+n-PrOOH -> EtCHO H2O OH MCMv3.2
! B192 Rate k is based on MCM v3.2. IUPAC only provides total cross
! B192 sections.
!----------------------------------------------------------------------
! B193 OH+n-PrOOH -> n-PrOO H2O  MCMv3.2
! B193 Rate k is based on MCM v3.2. IUPAC only provides total cross
! B193 sections.
!----------------------------------------------------------------------
! B194 i-PrOO+NO -> Me2CO HO2 NO2 IUPAC2005
! B194 -
!----------------------------------------------------------------------
! B195 i-PrOO+NO3 -> Me2CO HO2 NO2 MCMv3.2
! B195 No data in IUPAC or JPL. Take rates from MCM v3.2.
!----------------------------------------------------------------------
! B196 n-PrOO+NO -> EtCHO HO2 NO2 IUPAC2005
! B196 The recommendation accepts the Arrhenius expression of Eberhard
! B196 and Howard (1996).  Neglect n-propyl nitrate formation.
!----------------------------------------------------------------------
! B197 n-PrOO+NO3 -> EtCHO HO2 NO2 MCM3.2
! B197 No data in IUPAC or JPL. Take rates from MCM v3.2.
!----------------------------------------------------------------------


! Dry deposition coefficients. Only used for off-line dry-deposition scheme.
! Item number, species name, dry dep velocities x30, chemistry scheme, 
! qualifier, disqualifier, version

! Split to avoid warning "Too many continuation lines"
depvel_defs_master(1:29) = (/   &
! 1
depvel_t(1,'SO2       ',&
(/0.80,  0.80,  0.80,  0.80,  0.80,  0.80,& !     39.1
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,& !     39.2
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,& !     39.3
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,& !     39.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10/),& !   39.5
  TI+S+ST+OL+R,A,0,107),&
! 2
! no DD of DMSO in R scheme
depvel_t(2,'DMSO      ',&
(/1.00,  1.00,  1.00,  1.00,  1.00,  1.00,& !     42.1
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !     42.2
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !     42.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !     42.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00/),& !   42.5
  TI+ST+OL,A,0,107),&
! 3
depvel_t(3,'NH3       ',&
(/0.80,  0.80,  0.80,  0.80,  0.80,  0.80,& !     43.1
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,& !     43.2
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,& !     43.3
  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,& !     43.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10/),& !   43.5
  TI+ST+R,A,0,107),&
! 4 
depvel_t(4,'H2SO4     ',&
(/0.50,  0.50,  0.50,  0.50,  0.50,  0.50,& !     40.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,& !     40.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,& !     40.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,& !     40.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50/),& !   40.1
  TI+OL+S+ST+R,A,0,107),&
! 5
! No DD of Monoterp in R scheme
depvel_t(5,'Monoterp  ',&
(/0.20,  0.20,  0.20,  0.20,  0.20,  0.20,& !     44.1
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,& !     44.2
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,& !     44.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,& !     44.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10/),& !   44.5
  TI+ST+OL,A,0,107),&
! 6
! No DD of Sec_Org in R scheme
depvel_t(6,'Sec_Org   ',&
(/0.20,  0.20,  0.20,  0.20,  0.20,  0.20,& !     45.1
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,& !     45.2
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,& !     45.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,& !     45.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10/),& !   45.5
  TI+ST+OL,A,0,107),&
! 7
! R and T are at older revision than S and ST. Make consistent
depvel_t(7,'O3        ',& ! (Ganzeveld & Lelieveld (1995) note 1)
                          ! (modified to be the same as Guang version)
(/0.05,  0.05,  0.05,  0.05,  0.05,  0.05,& !      1.1
  0.85,  0.30,  0.65,  0.65,  0.25,  0.45,& !      1.2
  0.65,  0.25,  0.45,  0.65,  0.25,  0.45,& !      1.3
  0.18,  0.18,  0.18,  0.18,  0.18,  0.18,& !      1.4
  0.05,  0.05,  0.05,  0.05,  0.05,  0.05/),& !    1.5
  TI+ST,0,0,107),&
! 8*
depvel_t(7,'O3        ',&
(/0.05,  0.05,  0.05,  0.05,  0.05,  0.05,&
  1.00,  0.11,  0.56,  0.26,  0.11,  0.19,&
  1.00,  0.37,  0.69,  0.59,  0.46,  0.53,&
  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,&
  0.05,  0.05,  0.05,  0.05,  0.05,  0.05/),&
  T+R,0,0,107),&
! 9**
! O3 (Ganzeveld& Lelieveld (1995) - note 1)
depvel_t(7,'O3        ',&
(/0.07,  0.07,  0.07,  0.07,  0.07,  0.07,&
  1.00,  0.11,  0.56,  0.26,  0.11,  0.19,&
  1.00,  0.37,  0.69,  0.59,  0.46,  0.53,&
  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,&
  0.07,  0.07,  0.07,  0.07,  0.07,  0.07/),&
  S,0,0,107),&
! 10
! No DD of NO in R scheme
depvel_t(10,'NO        ',& ! (inferred from NO2 - see Giannakopoulos (1998))
(/0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !      2.1
  0.14,  0.01,  0.07,  0.01,  0.01,  0.01,& !      2.2
  0.10,  0.01,  0.06,  0.01,  0.01,  0.01,& !      2.3
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,& !      2.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00/),& !    2.5
  TI+T+S+ST,0,0,107),&
! 11
! No DD of NO3 in R scheme 
depvel_t(11,'NO3     ',& ! (as NO2)
(/0.02,  0.02,  0.02,  0.02,  0.02,  0.02,& !      3.1
  0.83,  0.04,  0.44,  0.06,  0.04,  0.05,& !      3.2
  0.63,  0.06,  0.35,  0.07,  0.06,  0.07,& !      3.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !      3.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01/),& !    3.5
  TI+T+S+ST,0,0,107),&
! 12
depvel_t(12,'NO2       ',&
(/0.02,  0.02,  0.02,  0.02,  0.02,  0.02,& !      4.1
  0.83,  0.04,  0.44,  0.06,  0.04,  0.05,& !      4.2
  0.63,  0.06,  0.35,  0.07,  0.06,  0.07,& !      4.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !      4.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01/),& !    4.5
  TI+T+S+ST+R,0,0,107),&
! 13
! No DD of N2O5 in R scheme
depvel_t(13,'N2O5      ',& ! (as HNO3)
(/1.00,  1.00,  1.00,  1.00,  1.00,  1.00,& !      5.1
  4.00,  3.00,  3.50,  2.00,  1.00,  1.50,& !      5.2
  2.50,  1.50,  2.00,  1.00,  1.00,  1.00,& !      5.3
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,& !      5.4
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50/),& !    5.5
  TI+T+S+ST,0,0,107),&
! 14
! No DD of HO2NO2 in R scheme
depvel_t(14,'HO2NO2    ',& ! (as HNO3)
(/1.00,  1.00,  1.00,  1.00,  1.00,  1.00,& !      6.1
  3.20,  1.40,  2.30,  2.00,  1.00,  1.50,& !      6.2
  1.80,  0.80,  1.30,  1.00,  1.00,  1.00,& !      6.3
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,& !      6.4
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50/),& !    6.5
  TI+ST,0,0,107),&
! 15*
depvel_t(14,'HO2NO2    ',& ! (as HNO3)
(/1.00,  1.00,  1.00,  1.00,  1.00,  1.00,&
  4.00,  3.00,  3.50,  2.00,  1.00,  1.50,&
  2.50,  1.50,  2.00,  1.00,  1.00,  1.00,&
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,&
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50/),&
  T+S,0,0,107),&  
! 16
! T and R are outdated versus S, ST
depvel_t(16,'HONO2     ',& ! (Zhang et al., 2003)
(/1.00,  1.00,  1.00,  1.00,  1.00,  1.00,& !      7.1
  3.20,  1.40,  2.30,  2.00,  1.00,  1.50,& !      7.2
  1.80,  0.80,  1.30,  1.00,  1.00,  1.00,& !      7.3
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,& !      7.4
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50/),& !    7.5
  TI+ST,0,0,107),&
! 17*
depvel_t(16,'HONO2     ',&
(/1.00,  1.00,  1.00,  1.00,  1.00,  1.00,&
  4.00,  3.00,  3.50,  2.00,  1.00,  1.50,&
  2.50,  1.50,  2.00,  1.00,  1.00,  1.00,&
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,&
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50/),&
  S+T+R,0,0,107),&
! 18
depvel_t(18,'H2O2      ',&
(/1.00,  1.00,  1.00,  1.00,  1.00,  1.00,& !      8.1
  1.25,  0.16,  0.71,  0.28,  0.12,  0.20,& !      8.2
  1.25,  0.53,  0.89,  0.83,  0.78,  0.81,& !      8.3
  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,& !      8.4
  0.32,  0.32,  0.32,  0.32,  0.32,  0.32/),& !    8.5
  TI+T+S+ST+OL+R,0,0,107),&
! 19
depvel_t(19,'CH4       ',&
(/0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !      9.1
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !      9.2
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !      9.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !      9.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00/),& !    9.5
  TI+T+R,0,0,107),&
! 20
depvel_t(20,'CO        ',& ! (see Giannokopoulos (1998))
(/0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !     10.1
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     10.2
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     10.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     10.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00/),& !   10.5
  TI+T+S+ST+R,0,0,107),&
! 21
! No DD for HCHO in R scheme
depvel_t(21,'HCHO      ',& ! (Zhang et al., 2003)
(/1.00,  0.60,  0.80,  1.00,  0.60,  0.80,& !     11.1
  1.60,  0.40,  1.00,  0.01,  0.01,  0.01,& !     11.2
  0.71,  0.20,  0.45,  0.03,  0.03,  0.03,& !     11.3
  0.40,  0.40,  0.40,  0.00,  0.00,  0.00,& !     11.4
  0.30,  0.30,  0.30,  0.30,  0.30,  0.30/),& !   11.5
  TI+ST,0,0,107),&
! 22*
depvel_t(21,'HCHO      ',&
(/1.00,  1.00,  1.00,  1.00,  1.00,  1.00,&
  1.00,  0.01,  0.50,  0.01,  0.01,  0.01,&
  0.71,  0.03,  0.37,  0.03,  0.03,  0.03,&
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,&
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04/),&
  S+T,0,0,107),&
! 23
depvel_t(23,'MeOOH     ',&
(/0.25,  0.25,  0.25,  0.10,  0.10,  0.10,& !     12.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,& !     12.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,& !     12.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     12.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01/),& !   12.5
  TI+T+S+ST+R,0,0,107),&
! 24
depvel_t(24,'HCl       ',& ! (same as HBr)
(/0.50,  0.50,  0.50,  0.50,  0.50,  0.50,& !     13.1
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,& !     13.2
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,& !     13.3
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,& !     13.4
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20/),& !   13.5
  S+ST,0,0,107),&
! 25
depvel_t(25,'HOCl      ',& ! (same as HOBr)
(/0.20,  0.20,  0.20,  0.20,  0.20,  0.20,& !     14.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,& !     14.2
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,& !     14.3
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,& !     14.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10/),& !   14.5
  S+ST,0,0,107),&
! 26
depvel_t(26,'HBr       ',& ! (note 3)
(/0.50,  0.50,  0.50,  0.50,  0.50,  0.50,& !     15.1
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,& !     15.2
  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,& !     15.3
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,& !     15.4
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20/),& !   15.5
  S+ST,0,0,107),&
! 27
depvel_t(27,'HOBr      ',& ! (note 3)
(/0.20,  0.20,  0.20,  0.20,  0.20,  0.20,& !     16.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,& !     16.2
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,& !     16.3
  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,& !     16.4
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10/),& !   16.5
  S+ST,0,0,107),&
! 28
depvel_t(28,'HONO      ',&
(/0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !     17.1
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !     17.2
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !     17.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !     17.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00/),& !   17.5
  TI+T+ST,0,0,107),&
! 29
depvel_t(29,'EtOOH     ',& !(as MeOOH)
(/0.25,  0.25,  0.25,  0.10,  0.10,  0.10,& !     18.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,& !     18.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,& !     18.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     18.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01/),& !   18.5
  TI+T+ST+R,0,0,107) /)

! Split to avoid warning "Too many continuation lines"
depvel_defs_master(30:n_dry_master) = (/&
! 30
! No DD in R scheme for MeCHO
depvel_t(30,'MeCHO     ',&
(/0.02,  0.02,  0.02,  0.02,  0.02,  0.02,& !     19.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,& !     19.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,& !     19.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !     19.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00/),& !   19.5
  TI+T+ST,0,0,107),&
! 31
! Make consistent 
depvel_t(31,'PAN       ',& ! (Zhang et al., 2003)
(/0.01,  0.01,  0.01,  0.01,  0.01,  0.01,& !     20.1
  0.63,  0.14,  0.39,  0.06,  0.05,  0.06,& !     20.2
  0.42,  0.14,  0.24,  0.07,  0.06,  0.07,& !     20.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,& !     20.4
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00/),& !   20.5
  TI+ST,0,0,107),&
! 32
depvel_t(32,'PAN       ',&
(/0.01,  0.01,  0.01,  0.01,  0.01,  0.01,&
  0.53,  0.04,  0.29,  0.06,  0.05,  0.06,&
  0.42,  0.06,  0.24,  0.07,  0.06,  0.07,&
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00/),&
  T+R,0,0,107),&
! 33
depvel_t(33,'n-PrOOH   ',& ! (as MeOOH)
(/0.25,  0.25,  0.25,  0.10,  0.10,  0.10,& !     21.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,& !     21.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,& !     21.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     21.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01/),& !   21.5
  TI+T+ST,0,0,107),&
! 34
depvel_t(34,'i-PrOOH   ',& ! (as MeOOH)
(/0.25,  0.25,  0.25,  0.10,  0.10,  0.10,& !     22.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,& !     22.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,& !     22.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     22.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01/),& !   22.5
  TI+T+ST+R,0,0,107),&
! 35
depvel_t(35,'EtCHO     ',& ! (as MeCHO)
(/0.02,  0.02,  0.02,  0.02,  0.02,  0.02,& !     23.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,& !     23.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,& !     23.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !     23.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00/),& !   23.5
  TI+T+ST,0,0,107),& 
! 36
depvel_t(36,'MeCOCH2OOH',& !(as MeCO3H - where do the desert data come from!?)
(/0.36,  0.36,  0.36,  0.19,  0.19,  0.19,& !     24.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,& !     24.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,& !     24.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     24.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02/),& !   24.5
  TI+T+ST,0,0,107),&
! 37
depvel_t(37,'PPAN      ',& ! (as PAN)
(/0.01,  0.01,  0.01,  0.01,  0.01,  0.01,& !     25.1
  0.63,  0.14,  0.39,  0.06,  0.05,  0.06,& !     25.2
  0.42,  0.14,  0.24,  0.07,  0.06,  0.07,& !     25.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,& !     25.4
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00/),& !   25.5
  TI+ST,0,0,107),&
! 38*
depvel_t(37,'PPAN      ',&
(/0.01,  0.01,  0.01,  0.01,  0.01,  0.01,&
  0.53,  0.04,  0.29,  0.06,  0.05,  0.06,&
  0.42,  0.06,  0.24,  0.07,  0.06,  0.07,&
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00/),&
  T,0,0,107),&
! 39
! Inconsistent DD for ISOOH
depvel_t(39,'ISOOH     ',& ! (as MeCO3H) (MIM)
(/0.36,  0.36,  0.36,  0.19,  0.19,  0.19,& !     26.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,& !     26.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,& !     26.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     26.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02/),& !   26.5
  TI+ST,0,0,107),&
! 40*
depvel_t(39,'ISOOH     ',&
(/0.25,  0.25,  0.25,  0.10,  0.10,  0.10,& !ISOOH
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,&
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,&
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,&
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01/),&
  R,0,0,107),&
! 41
depvel_t(41,'ISON      ',& ! (as MeCO2H) (MIM) - note 2
(/0.50,  0.50,  0.50,  0.25,  0.25,  0.25,& !     27.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,& !     27.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,& !     27.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,& !     27.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03/),& !   27.5
  TI+ST,0,0,107),&
! 42
depvel_t(42,'MACR      ',& ! (as MeCHO) (MIM)
(/0.02,  0.02,  0.02,  0.02,  0.02,  0.02,& !     28.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,& !     28.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,& !     28.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !     28.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00/),& !   28.5
  TI+ST,0,0,107),&
! 43
depvel_t(43,'MACROOH   ',& ! (as MeCO3H) (MIM)
(/0.36,  0.36,  0.36,  0.19,  0.19,  0.19,& !     29.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,& !     29.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,& !     29.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     29.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02/),& !   29.5
  TI+ST,0,0,107),&
! 44
depvel_t(44,'MPAN      ',& ! (as PAN) (MIM)
(/0.01,  0.01,  0.01,  0.01,  0.01,  0.01,& !     30.1
  0.63,  0.14,  0.39,  0.06,  0.05,  0.06,& !     30.2
  0.42,  0.14,  0.24,  0.07,  0.06,  0.07,& !     30.3
  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,& !     30.4
  0.01,  0.00,  0.01,  0.01,  0.00,  0.00/),& !   30.5
  TI+ST,0,0,107),&
! 45
depvel_t(45,'HACET     ',& ! (as MeCO3H) (MIM)
(/0.36,  0.36,  0.36,  0.19,  0.19,  0.19,& !     31.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,& !     31.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,& !     31.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     31.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02/),& !   31.5
  TI+ST,0,0,107),&
! 46
depvel_t(46,'MGLY      ',& ! (as MeCHO) (MIM)
(/0.02,  0.02,  0.02,  0.02,  0.02,  0.02,& !     32.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,& !     32.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,& !     32.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !     32.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00/),& !   32.5
  TI+ST,0,0,107),&
! 47
depvel_t(47,'NALD      ',& ! (as MeCHO) (MIM)
(/0.02,  0.02,  0.02,  0.02,  0.02,  0.02,& !     33.1
  0.31,  0.00,  0.16,  0.00,  0.00,  0.00,& !     33.2
  0.26,  0.00,  0.13,  0.00,  0.00,  0.00,& !     33.3
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !     33.4
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00/),& !   33.5
  TI+ST,0,0,107),&
! 48
depvel_t(48,'HCOOH     ',& ! (from Sander & Crutzen (1996) - note 3) (MIM)
(/0.50,  0.50,  0.50,  0.25,  0.25,  0.25,& !     34.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,& !     34.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,& !     34.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,& !     34.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03/),& !   34.5
  TI+ST,0,0,107),&
! 49
depvel_t(49,'MeCO3H    ',& ! (from Giannakopoulos (1998)) (MIM)
(/0.36,  0.36,  0.36,  0.19,  0.19,  0.19,& !     35.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,& !     35.2
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,& !     35.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     35.4
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02/),& !   35.5
  TI+ST,0,0,107),&
! 50
depvel_t(50,'MeCO2H    ',& ! (as HCOOH) (MIM)
(/0.50,  0.50,  0.50,  0.25,  0.25,  0.25,& !     36.1
  1.00,  0.06,  0.53,  0.08,  0.07,  0.08,& !     36.2
  0.74,  0.10,  0.42,  0.11,  0.10,  0.11,& !     36.3
  0.04,  0.04,  0.04,  0.04,  0.04,  0.04,& !     36.4
  0.03,  0.02,  0.03,  0.03,  0.02,  0.03/),& !   36.5
  TI+ST,0,0,107),&
! 51
depvel_t(51,'MeOH      ',& ! (as MeOOH)
(/0.25,  0.25,  0.25,  0.10,  0.10,  0.10,& !     38.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,& !     38.2
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,& !     38.3
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     38.4
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01/),& !   38.5
  TI+ST,0,0,107),&
! 52
depvel_t(52,'O3S       ',&
(/0.07,  0.07,  0.07,  0.07,  0.07,  0.07,& !     46.1
  1.00,  0.11,  0.56,  0.26,  0.11,  0.19,& !     46.1
  1.00,  0.37,  0.69,  0.59,  0.46,  0.53,& !     46.1
  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,& !     46.1
  0.07,  0.07,  0.07,  0.07,  0.07,  0.07/),& !   46.1
  T+R,0,0,107),&
! 53
depvel_t(53,'MVKOOH    ',&
(/0.36,  0.36,  0.36,  0.19,  0.19,  0.19,& !     47.1
  0.71,  0.04,  0.38,  0.06,  0.05,  0.06,& !     47.1
  0.53,  0.07,  0.30,  0.08,  0.07,  0.08,& !     47.1
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     47.1
  0.02,  0.01,  0.02,  0.02,  0.01,  0.02/),& !   47.1
  R,0,0,107),&
! 54
depvel_t(54,'ORGNIT    ',&
(/0.25,  0.25,  0.25,  0.10,  0.10,  0.10,& !     48.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,& !     48.1
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,& !     48.1
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     48.1
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01/),& !   48.1
  R,0,0,107),&
! 55
depvel_t(55,'H2        ',&
(/0.00,  0.00,  0.00,  0.00,  0.00,  0.00,& !     37.1
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     37.1
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     37.1
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     37.1
  0.00,  0.00,  0.00,  0.00,  0.00,  0.00/),& !   37.1
  R,0,0,107),&
! 56
depvel_t(56,'s-BuOOH   ',&
(/0.25,  0.25,  0.25,  0.10,  0.10,  0.10,& !     49.1
  0.83,  0.04,  0.44,  0.06,  0.05,  0.06,& !     49.1
  0.63,  0.06,  0.35,  0.08,  0.06,  0.07,& !     49.1
  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,& !     49.1
  0.01,  0.01,  0.01,  0.01,  0.01,  0.01/),& !   49.1
  R,0,0,107),&
! 57
depvel_t(57,'MSA       ',&
(/0.50,  0.50,  0.50,  0.50,  0.50,  0.50,& !     41.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,& !     41.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,& !     41.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,& !     41.1
  0.50,  0.50,  0.50,  0.50,  0.50,  0.50/),& !   41.1
  TI,A,0,107) /)

! determine which chemistry is to be used. Test here that only one scheme is 
! selected.
ierr = 0
chem_scheme = 0
i = 0
IF (l_ukca_strat) THEN
  chem_scheme = S
  i = i + 1
END IF
IF (l_ukca_strattrop) THEN
  chem_scheme = ST
  i = i + 1
END IF
IF (l_ukca_trop) THEN
  chem_scheme = T
  i = i + 1
END IF
IF (l_ukca_tropisop) THEN
  chem_scheme = TI
  i = i + 1
END IF
IF (l_ukca_raq         ) THEN
  chem_scheme = R
  i = i + 1
END IF
IF (l_ukca_offline     ) THEN
  chem_scheme = OL
  i = i + 1
END IF

! Determine which qualifiers to use, based on namelist options.
qualifier = 0
IF (l_ukca_achem .OR. l_ukca_offline) qualifier = A
IF (l_ukca_het_psc) qualifier = qualifier + HP
IF (l_ukca_trophet) qualifier = qualifier + TH
IF (l_ukca_stratcfc)qualifier = qualifier + ES
IF (l_co2_interactive) qualifier = qualifier + CI

! No chemistry scheme selected / chemistry scheme not implemented here
IF (chem_scheme == 0) THEN
  ierr = 11
  cmessage = 'No chemistry scheme selected'
END IF
IF (i > 1) THEN
  ierr = 12
  cmessage = "Can't cope with more than one chemistry scheme."
END IF

! calculate sizes of arrays to be extracted
jpspec = 0
jpbk = 0
jppj = 0
jptk = 0
jphk = 0
jpdd = 0
jpdw = 0

! select species and work out number of dry- and wet-deposited species. 
IF (ANY(IAND(chch_defs_master%qual,chch_defs_master%disqual) > 0)) THEN
  ierr=30
  cmessage = 'Do not select and deselect the same species.'
END IF

last_chch = chch_t1( 0,'',  0,'','',  0,  0,0, 0,0,0,0)
DO i=1,n_chch_master
! only use the entry if the version is less than or equal to the selected 
! version
  IF (chch_defs_master(i)%version <= version) THEN
! only use if it is part of the chosen chemistry scheme 
    IF (IAND(chch_defs_master(i)%scheme,chem_scheme) > 0) THEN
      IF (match(chch_defs_master(i)%qual,chch_defs_master(i)%disqual,&
          qualifier)) THEN
! only use if it is of a higher version than the last chosen item, if that's 
! for the same item
        take_this = .FALSE.
        IF (chch_defs_master(i)%item_no /= last_chch%item_no) THEN
          take_this = .TRUE.
        ELSE IF (chch_defs_master(i)%version == last_chch%version) THEN
! Can't have two entries for the same item number and version
          ierr = 1
          WRITE(cmessage,'(A68,I4)') &
          'Two entries in chch_defs_master with same item and version numbers.'&
          ,chch_defs_master(i)%item_no
        ELSE IF (chch_defs_master(i)%version > last_chch%version) THEN
! This is a higher version than one previously selected.
          take_this = .TRUE.
        END IF
        IF (take_this) THEN
! Fresh item, no previous version had been selected
          IF (chch_defs_master(i)%item_no /= last_chch%item_no) THEN
! the case of this being the first entry of this item number. Increase 
! counts jpspec, jpdd, and jpdw by 1 as needed.
            jpspec = jpspec + 1
            jpdd = jpdd + chch_defs_master(i)%switch1
            jpdw = jpdw + chch_defs_master(i)%switch2
          ELSE 
! This is the case of this being a higher version than one previously taken. 
! Do not increase jpspec (done that for the older version). Increase or 
! decrease jpdd and jpdw if there are differences between the two versions.
            jpdd = jpdd + chch_defs_master(i)%switch1 - last_chch%switch1
            jpdw = jpdw + chch_defs_master(i)%switch2 - last_chch%switch2
          END IF
! memorize entry for next round
          last_chch = chch_defs_master(i)          
        END IF
      END IF
    END IF
  END IF
END DO  


! work out number of bimolecular reactions
last_bimol = ratb_t1(0,'','','','','','',0.,0.,0., 0., 0., 0., 0.,0,0,0,0)

IF (ANY(IAND(ratb_defs_master%qual,ratb_defs_master%disqual) > 0)) THEN
  ierr=40
  cmessage = 'Do not select and deselect the same bimol reactions.'
END IF
DO i=1,n_bimol_master
! only count if part of selected chemistry scheme
  IF (IAND(ratb_defs_master(i)%scheme,chem_scheme) > 0) THEN
    IF (match(ratb_defs_master(i)%qual,ratb_defs_master(i)%disqual,qualifier))&
       THEN
! only count if version is less than or equal to chosen global version
      IF (ratb_defs_master(i)%version <= version) THEN
! if previously chosen item has same item number, ignore this
        take_this = .FALSE.
        IF (ratb_defs_master(i)%item_no /= last_bimol%item_no) THEN
          take_this = .TRUE.
        ELSE IF (ratb_defs_master(i)%version == last_bimol%version) THEN
! we can't have two entries with the same item number and version.
          ierr = 2
          WRITE (cmessage,'(A68,I4)') &
   'Two entries in ratb_defs_master with same item and version numbers.',&
            ratb_defs_master(i)%item_no
        ELSE IF (ratb_defs_master(i)%version > last_bimol%version) THEN
! higher version
          take_this = .TRUE.
        END IF
        IF (take_this) THEN
          IF (ratb_defs_master(i)%item_no /= last_bimol%item_no) THEN
            jpbk = jpbk + 1
          END IF
          last_bimol = ratb_defs_master(i)
        END IF
      END IF
    END IF
  END IF
END DO 

! find number of termolecular reactions
last_termol = ratt_t1(0,'','','','',0.,0.,0.,0.,0.,0.,0.,0., 0.,0,0,0,0)
IF (ANY(IAND(ratt_defs_master%qual,ratt_defs_master%disqual) > 0)) THEN
  ierr=50
  cmessage = 'Do not select and deselect the same termol reactions.'
END IF
DO i=1,n_ratt_master
  IF (IAND(ratt_defs_master(i)%scheme,chem_scheme) > 0) THEN
    IF (match(ratt_defs_master(i)%qual,ratt_defs_master(i)%disqual,qualifier))&
    THEN
      IF (ratt_defs_master(i)%version <= version) THEN
        take_this = .FALSE.
        IF (ratt_defs_master(i)%item_no /= last_termol%item_no) THEN
          take_this = .TRUE.
        ELSE IF (ratt_defs_master(i)%version == last_termol%version) THEN
          ierr = 3
          cmessage = &
          'Two items in ratt_defs_master with same item and version numbers.'
        ELSE IF (ratt_defs_master(i)%version > last_termol%version) THEN
          take_this = .TRUE.
        END IF
        IF (take_this) THEN
          IF (ratt_defs_master(i)%item_no /= last_termol%item_no) THEN
            jptk = jptk + 1
          END IF
          last_termol = ratt_defs_master(i)
        END IF
      END IF
    END IF
  END IF
END DO  

! find number of photolysis reactions
last_photol = ratj_t1(0,'','','','','','',0.,0.,0.,0.,0.,'',0,0,0,0)
IF (ANY(IAND(ratj_defs_master%qual,ratj_defs_master%disqual) > 0)) THEN
  ierr=60
  cmessage = 'Do not select and deselect the same photol reactions.'
END IF
DO i=1,n_ratj_master
  IF (IAND(ratj_defs_master(i)%scheme,chem_scheme) > 0) THEN
    IF (match(ratj_defs_master(i)%qual,ratj_defs_master(i)%disqual,qualifier))&
    THEN
      IF (ratj_defs_master(i)%version <= version) THEN
        take_this = .FALSE.
        IF (ratj_defs_master(i)%item_no /= last_photol%item_no) THEN
          take_this = .TRUE.
        ELSE IF (ratj_defs_master(i)%version == last_photol%version) THEN
          ierr = 4
          WRITE(cmessage,'(A58,I4)') &
            'Two items in ratj_defs_master with same item and version ',&
            ratj_defs_master(i)%item_no
        ELSE IF (ratj_defs_master(i)%version > last_photol%version) THEN
          take_this = .TRUE.
        END IF
        IF (take_this) THEN
          IF (ratj_defs_master(i)%item_no /= last_photol%item_no) &
            jppj = jppj + 1
          last_photol = ratj_defs_master(i)
        END IF 
      END IF
    END IF
  END IF
END DO

! find number of heterogeneous reactions
last_het = rath_t1(0,'','','','','','',0.,0.,0.,0.,0,0,0,0)
IF (ANY(IAND(rath_defs_master%qual,rath_defs_master%disqual) > 0)) THEN
  ierr=70
  cmessage = 'Do not select and deselect the same hetero reactions.'
END IF
DO i=1,n_het_master
  IF (IAND(rath_defs_master(i)%scheme,chem_scheme) > 0) THEN
    IF (match(rath_defs_master(i)%qual,rath_defs_master(i)%disqual,qualifier)) &
    THEN
      IF (rath_defs_master(i)%version <= version) THEN
        take_this = .FALSE.
        IF (rath_defs_master(i)%item_no /= last_het%item_no) THEN
          take_this = .TRUE.
        ELSE IF (rath_defs_master(i)%version == last_het%version) THEN
          ierr = 5
          cmessage = &
          'Two items in rath_defs_master with same item and version numbers.'
        ELSE IF (rath_defs_master(i)%version > last_het%version) THEN
          take_this = .TRUE.
        END IF
        IF (take_this) THEN
          IF (rath_defs_master(i)%item_no /= last_het%item_no) THEN
            jphk = jphk + 1
          END IF
          last_het = rath_defs_master(i)
        END IF
      END IF
    END IF
  END IF
END DO  

! calculate total number of reactions
jpnr = jpbk + jptk + jppj + jphk

! allocate final arrays for species and reactions
ALLOCATE(chch_defs_new(jpspec))
ALLOCATE(ratb_defs_new(jpbk))
ALLOCATE(ratj_defs_new(jppj))
ALLOCATE(ratt_defs_new(jptk))
ALLOCATE(rath_defs_new(jphk))
ALLOCATE(depvel_defs_new(jddept,jddepc,jpdd))
ALLOCATE(henry_defs_new(6,jpdw))

! pick out species needed
j = 0
last_chch = chch_t1( 0,'',  0,'','',  0,  0,0, 0,0,0,0)
DO i=1,n_chch_master
! only use the entry if the version is less than or equal to the selected 
! version
  IF ((chch_defs_master(i)%version <= version) .AND. &
   (IAND(chch_defs_master(i)%scheme,chem_scheme) > 0) .AND. &
   (match(chch_defs_master(i)%qual,chch_defs_master(i)%disqual,qualifier))) &
   THEN
! increase count only if fresh item not just different version of previous 
! item. If not, overwrite previously chosen item.
    IF (chch_defs_master(i)%item_no /= last_chch%item_no) j = j + 1
! If new item or higher version of previously selected item, copy into array
    IF ((chch_defs_master(i)%item_no /= last_chch%item_no) .OR. &
        (chch_defs_master(i)%version > last_chch%version)) THEN
! memorize entry for next round
      chch_defs_new(j) = chch_defs_master(i)     
      last_chch = chch_defs_master(i)    
    END IF
  END IF
END DO

WRITE (umMessage,'(A9,I4)') 'jpspec = ',jpspec
CALL umPrint(umMessage,src=RoutineName)
WRITE (umMessage,'(A4)') 'CHCH'
CALL umPrint(umMessage,src=RoutineName)
DO i=1,jpspec
  WRITE (umMessage,'(I2,1X,A10,1X,A2,3(1X,I1),1X,I3)') &
         chch_defs_new(i)%item_no,&
         chch_defs_new(i)%speci,chch_defs_new(i)%ctype,&
         chch_defs_new(i)%switch1,chch_defs_new(i)%switch2,&
         chch_defs_new(i)%switch3,chch_defs_new(i)%version
  CALL umPrint(umMessage,src=RoutineName)
END DO

! Pick out bimolecular reactions needed
j=0
last_bimol = ratb_t1(0,'','','','','','',0.,0.,0., 0., 0., 0., 0.,0,0,0,0)
DO i=1,n_bimol_master
  IF ((IAND(ratb_defs_master(i)%scheme,chem_scheme) > 0) .AND. &
    (ratb_defs_master(i)%version <= version) .AND. &
    (match(ratb_defs_master(i)%qual,ratb_defs_master(i)%disqual,qualifier))) &
  THEN
    IF (ratb_defs_master(i)%item_no /= last_bimol%item_no) j = j + 1
    IF ((ratb_defs_master(i)%item_no /= last_bimol%item_no) .OR. &
        (ratb_defs_master(i)%version > last_bimol%version)) THEN
      ratb_defs_new(j) = ratb_defs_master(i)
      last_bimol = ratb_defs_master(i)
    END IF
  END IF
END DO 
 
! Pick out termolecular reactions needed
j=0
last_termol = ratt_t1(0,'','','','',0.,0.,0.,0.,0.,0.,0.,0., 0.,0,0,0,0)
DO i=1,n_ratt_master
  IF ((IAND(ratt_defs_master(i)%scheme,chem_scheme) > 0) .AND. &
    (ratt_defs_master(i)%version <= version) .AND. &
    (match(ratt_defs_master(i)%qual,ratt_defs_master(i)%disqual,qualifier))) &
  THEN
    IF (ratt_defs_master(i)%item_no /= last_termol%item_no) j = j + 1
    IF ((ratt_defs_master(i)%item_no /= last_termol%item_no) .OR. &
      (ratt_defs_master(i)%version > last_termol%version)) THEN
      ratt_defs_new(j) = ratt_defs_master(i)
      last_termol = ratt_defs_master(i)
    END IF
  END IF
END DO 

! Pick out photolysis reactions needed
j=0
last_photol = ratj_t1(0,'','','','','','',0.,0.,0.,0.,0.,'',0,0,0,0)
DO i=1,n_ratj_master
  IF ((IAND(ratj_defs_master(i)%scheme,chem_scheme) > 0) .AND. &
    (ratj_defs_master(i)%version <= version) .AND. &
    (match(ratj_defs_master(i)%qual,ratj_defs_master(i)%disqual,qualifier))) &
  THEN
    IF (ratj_defs_master(i)%item_no /= last_photol%item_no) j = j + 1
    IF ((ratj_defs_master(i)%item_no /= last_photol%item_no) .OR. &
        (ratj_defs_master(i)%version > last_photol%version)) THEN
      ratj_defs_new(j) = ratj_defs_master(i)
      last_photol = ratj_defs_master(i)
    END IF
  END IF
END DO 

! Pick out heterogeneous reactions needed
j=0
last_het = rath_t1(0,'','','','','','',0.,0.,0.,0.,0,0,0,0)
DO i=1,n_het_master
  IF ((IAND(rath_defs_master(i)%scheme,chem_scheme) > 0) .AND. &
    (rath_defs_master(i)%version <= version) .AND. &
    (match(rath_defs_master(i)%qual,rath_defs_master(i)%disqual,qualifier))) &
  THEN
    IF (rath_defs_master(i)%item_no /= last_het%item_no) j = j + 1
    IF ((rath_defs_master(i)%item_no /= last_het%item_no) .OR. &
        (rath_defs_master(i)%version > last_het%version)) THEN
      rath_defs_new(j) = rath_defs_master(i)
      last_het = rath_defs_master(i)
    END IF
  END IF
END DO 

 
! pick out species with dry deposition
j = 0
last_depvel = depvel_t(0,'',(/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,&
                              0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,&
                              0.,0.,0.,0.,0.,0.,0.,0.,0.,0./),0,0,0,0)
DO i=1,jpspec
  IF (chch_defs_new(i)%switch1 == 1) THEN
    found = .FALSE.
    j = j + 1
    DO k=1,n_dry_master
! check for matching species, chemistry scheme, and version
      IF (depvel_defs_master(k)%speci == chch_defs_new(i)%speci) THEN
        IF (IAND(depvel_defs_master(k)%scheme,chem_scheme) > 0) THEN
          IF (depvel_defs_master(k)%version <= version) THEN
            IF (match(depvel_defs_master(k)%qual,depvel_defs_master(k)%disqual,&
                qualifier)) THEN
              take_this = .FALSE.
              IF (depvel_defs_master(k)%item_no == last_depvel%item_no) THEN
! can't have two entries at the same version
                IF (depvel_defs_master(k)%version == last_depvel%version) THEN
                  ierr = 7
                  WRITE (cmessage,'(A70,I4)') &
    'Two entries in depvel_defs_master with same item and version numbers.',&
                    depvel_defs_master(k)%speci
                END IF
                IF (depvel_defs_master(k)%version > last_depvel%version) &
                  take_this = .TRUE.
              ELSE
                take_this = .TRUE.
              END IF
              IF (take_this) THEN
                depvel_defs_new(:,:,j) = RESHAPE(&
                  depvel_defs_master(k)%velocity,(/jddept,jddepc/))
                last_depvel =  depvel_defs_master(k)
                found = .TRUE.
              END IF
            END IF
          END IF
        END IF
      END IF
    END DO
! check for whether dry deposition parameters have been located.
    IF (.NOT.(found)) THEN
      ierr = 8
      cmessage = 'Dry deposition selected for '//chch_defs_new(i)%speci//&
        ' but dry deposition parameters missing.'
    END IF
  END IF
END DO 

! pick out species with wet deposition
j = 0
last_wetdep = wetdep(0,'',(/0.,0.,0.,0.,0.,0./),0,0,0,0)
DO i=1,jpspec
  IF (chch_defs_new(i)%switch2 == 1) THEN
    found = .FALSE.
    j = j + 1
    DO k=1,n_wet_master
! check for matching species, chemistry scheme, and version
      IF (henry_defs_master(k)%speci == chch_defs_new(i)%speci) THEN
        IF (IAND(henry_defs_master(k)%scheme,chem_scheme) > 0) THEN
          IF (henry_defs_master(k)%version <= version) THEN
            IF (match(henry_defs_master(k)%qual,henry_defs_master(k)%disqual,&
                qualifier)) THEN
              take_this = .FALSE.
              IF (henry_defs_master(k)%item_no == last_wetdep%item_no) THEN
! can't have two entries at the same version
                IF (henry_defs_master(k)%version == last_wetdep%version) THEN
                  ierr = 9
                  cmessage = &
    'Two entries in henry_defs_master with same item and version numbers.'
                END IF
                IF (henry_defs_master(k)%version > last_wetdep%version) &
                  take_this = .TRUE.
              ELSE
                take_this = .TRUE.
              END IF
              IF (take_this) THEN
                henry_defs_new(:,j) = henry_defs_master(k)%params
                last_wetdep =  henry_defs_master(k)
                found = .TRUE.
              END IF
            END IF
          END IF
        END IF
      END IF
    END DO
! check for whether dry deposition parameters have been located
    IF (.NOT.(found)) THEN
      ierr = 10
      cmessage = 'Wet deposition selected for '//chch_defs_new(i)%speci//&
        ' but wet deposition parameters missing.'      
    END IF
  END IF
END DO

! calculate total number of tracers
jpctr = 0
DO i=1,jpspec
  IF (TRIM(chch_defs_new(i)%ctype) == 'TR') jpctr = jpctr + 1
END DO

WRITE(umMessage,'(A8,I4)') 'jpctr  = ',jpctr
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A8,I4)') 'jpbk   = ',jpbk
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A4)') 'RATB'
CALL umPrint(umMessage,src=RoutineName)
DO i=1,jpbk
  WRITE (umMessage,'(I3,1X,6(A10,1X))') &
                       ratb_defs_new(i)%item_no,ratb_defs_new(i)%react1,&
                       ratb_defs_new(i)%react2,ratb_defs_new(i)%prod1,&
                       ratb_defs_new(i)%prod2,ratb_defs_new(i)%prod3,&
                       ratb_defs_new(i)%prod4
  CALL umPrint(umMessage,src=RoutineName)
  WRITE (umMessage,'(I3,3(1X,E9.2))') &
                       ratb_defs_new(i)%version,&
                       ratb_defs_new(i)%k0   ,ratb_defs_new(i)%alpha,&
                       ratb_defs_new(i)%beta
  CALL umPrint(umMessage,src=RoutineName)
END DO
WRITE (umMessage,'(A8,I4)') 'jppj   = ',jppj
CALL umPrint(umMessage,src=RoutineName)
WRITE (umMessage,'(A4)') 'RATJ'
CALL umPrint(umMessage,src=RoutineName)
DO i=1,jppj
  WRITE (umMessage,'(I2,1X,7(A10,1X),I3)') ratj_defs_new(i)%item_no,&
    ratj_defs_new(i)%react1,ratj_defs_new(i)%react2,ratj_defs_new(i)%prod1,&
    ratj_defs_new(i)%prod2,ratj_defs_new(i)%prod3,ratj_defs_new(i)%prod4,&
    ratj_defs_new(i)%fname,ratj_defs_new(i)%version
  CALL umPrint(umMessage,src=RoutineName)
END DO
WRITE (umMessage,'(A8,I4)') 'jptk   = ',jptk
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A4)') 'RATT'
CALL umPrint(umMessage,src=RoutineName)
DO i=1,jptk
  WRITE(umMessage,'(I2,1X,4(A10,1X))') ratt_defs_new(i)%item_no,&
    ratt_defs_new(i)%react1,ratt_defs_new(i)%react2,ratt_defs_new(i)%prod1,&
    ratt_defs_new(i)%prod2
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(I3,6(1X,E9.2))') &
    ratt_defs_new(i)%version,ratt_defs_new(i)%k1,&
    ratt_defs_new(i)%alpha1,ratt_defs_new(i)%beta1  ,ratt_defs_new(i)%k2,&
    ratt_defs_new(i)%alpha2,ratt_defs_new(i)%beta2
  CALL umPrint(umMessage,src=RoutineName)
END DO
WRITE(umMessage,'(A8,I4)') 'jphk   = ',jphk
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A4)') 'RATH'
CALL umPrint(umMessage,src=RoutineName)
DO i=1,jphk
  WRITE(umMessage,'(I2,1X,6(A10,1X),I3)') rath_defs_new(i)%item_no,&
    rath_defs_new(i)%react1,rath_defs_new(i)%react2,rath_defs_new(i)%prod1,&
    rath_defs_new(i)%prod2,rath_defs_new(i)%prod3,rath_defs_new(i)%prod4,&
    rath_defs_new(i)%version
  CALL umPrint(umMessage,src=RoutineName)
END DO

WRITE(umMessage,'(A8,I4)') 'jpdd   = ',jpdd
CALL umPrint(umMessage,src=RoutineName)
j=0
DO i=1,jpspec 
  IF (chch_defs_new(i)%switch1 == 1) THEN
    j=j+1
    WRITE(umMessage,'(A10)') chch_defs_new(i)%speci
    CALL umPrint(umMessage,src=RoutineName)
    DO k=1,jddepc
      WRITE(umMessage,'(6(F5.2))') depvel_defs_new(:,k,j)
      CALL umPrint(umMessage,src=RoutineName)
    END DO
  END IF
END DO

j=0
WRITE(umMessage,'(A8,I4)') 'jpdw   = ',jpdw
CALL umPrint(umMessage,src=RoutineName)
DO i=1,jpspec
  IF (chch_defs_new(i)%switch2 == 1) THEN
    j=j+1
    WRITE(umMessage,'(A10)') chch_defs_new(i)%speci
    CALL umPrint(umMessage,src=RoutineName)
    WRITE(umMessage,'(6(E10.3))') henry_defs_new(:,j)
    CALL umPrint(umMessage,src=RoutineName)
  END IF
END DO

IF (ierr > 1) THEN 
  CALL ereport(RoutineName,ierr,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_chem_master


! ######################################################################

LOGICAL FUNCTION match(qual, disqual, qualifier)

! If qual == 0 and disqual == 0, select unconditionally.
! If qual > 0, only select if the corresponding flag match is True.
! If disqual > 0, select if the corresponding flag is False, 
!                 deselect if the corresponding flag is True.
IMPLICIT NONE

INTEGER :: qual, disqual, qualifier

match = ((qual == 0) .AND. (disqual == 0.))

IF (qual > 0) THEN
  IF (IAND(qual,  qualifier) > 0) match = .TRUE.
END IF

IF (disqual > 0) THEN
  IF (IAND(disqual,qualifier) > 0) THEN
    match = .FALSE.
  ELSE
    match = .TRUE.
  END IF
END IF

RETURN
END FUNCTION match


END MODULE ukca_chem_master_mod
