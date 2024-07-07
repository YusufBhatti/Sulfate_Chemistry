! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module containing the reactions which define the chemical reaction
!  fluxes for STASH output.
!
!  Part of the UKCA model. UKCA is a community model supported
!  by The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!    This code is written to UMDP3 v8 programming standards.
!
! ######################################################################
!
! TP       = RXN for reaction
!            DEP for deposition
!            EMS for emission
!            NET for a tendency of a particular species
!            STE for the Strat-Trop exchange across the 2PVU+380K Tropopause
!                for a particular speices
!            MAS for the mass of air per gridcell
!            PSC for location of type 1/2 PSCs
! STASH#   = First 2 digits are section number, last 3 are item number,
!            5 digits in total, e.g. 34412 etc.
! WHICH    = IF TP==RXN: B = bimolecular         3D
!                        H = heterogeneous       3D
!                        T = termolecular        3D
!                        J = photolysis          3D
!            IF TP==DEP: D = dry deposition      2D
!                        W = wet deposition      3D
!            IF TP==EMS: S = surface emissions   2D
!                        A = aircraft emissions  3D
!                        L = lightning emissions 3D
!                        V = volcanic emissions  3D (only available if have SO2 tracer)
!                        S = 3D SO2   emissions  3D (only available if have SO2 tracer)
!            IF TP==NET  X = not used
!            IF TP==STE  X = not used
!            IF TP==MAS  X = not used
!            IF TP==PSC  1/2 3D
! RXN#     = If there are multiple reactions with the same reactants AND products
!            then this number specifies which of those (in order) you wish to budget.
!            Set to 0 otherwise.
! #SPECIES = number of species. IF TP==RXN this is total of reactants + products
!                               IF TP!=RXN this should be 1 ONLY
! R1       = IF TP!=RXN this is the species being deposited/emitted
!
!
! NOTE: If the same stash number is specified for multiple diagnostics then these
!       will be summed up before outputting into stash for processing. Specifing
!       the same reaction multiple times is allowed, but only works for diagnostics
!       of the same dimensionallity, i.e. 2D *or* 3D
!
! UNITS OF DIAGNOSTICS ARE: moles/gridcell/s
!
! ######################################################################

MODULE asad_flux_dat

USE missing_data_mod,      ONLY: imdi

IMPLICIT NONE

PRIVATE

! Describes the bimolecular reaction rates
TYPE asad_flux_defn
  CHARACTER(LEN=3)  :: diag_type         ! which type of flux:RXN,DEP,EMS
  INTEGER           :: stash_number      ! stash number, e.g. 50001 etc.
  CHARACTER(LEN=1)  :: rxn_type          ! which rxn type: B,H,T,J,D,W,S,A,L,V
  LOGICAL           :: tropospheric_mask ! T or F
  INTEGER           :: rxn_location      ! for multiple reacions with the same
                                         ! reactants and products, default=0
  INTEGER           :: num_species       ! number of reactants+products
  CHARACTER(LEN=10) :: reactants(2)      ! list of reactants
  CHARACTER(LEN=10) :: products(4)       ! list of products
END TYPE asad_flux_defn
PUBLIC asad_flux_defn

PUBLIC :: asad_chemical_fluxes

CHARACTER(LEN=10) :: blank0 = '          '   ! Defines null product

! Number of chemical fluxes defined below
INTEGER, PARAMETER :: n_chemical_fluxes = 310

TYPE(asad_flux_defn), ALLOCATABLE, SAVE :: asad_chemical_fluxes(:)

TYPE(asad_flux_defn), PARAMETER :: asad_trop_ox_budget_prod(20) = &
! Production of Ox
       (/ asad_flux_defn('RXN',50001,'B',.TRUE.,0,4,                    &
       (/'HO2       ','NO        '/),                                   &
       (/'OH        ','NO2       ','          ','          '/)),        &
       asad_flux_defn('RXN',50002,'B',.TRUE.,0,5,                       &
       (/'MeOO      ','NO        '/),                                   &
       (/'HO2       ','HCHO      ','NO2       ','          '/)),        &

! NO + RO2 reactions: sum into STASH section 50 item 3
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                       &
       (/'EtOO      ','NO        '/),                                   &
       (/'MeCHO     ','HO2       ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                       &
       (/'MeCO3     ','NO        '/),                                   &
       (/'MeOO      ','CO2       ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                       &
       (/'n-PrOO    ','NO        '/),                                   &
       (/'EtCHO     ','HO2       ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                       &
       (/'i-PrOO    ','NO        '/),                                   &
       (/'Me2CO     ','HO2       ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                       &
       (/'EtCO3     ','NO        '/),                                   &
       (/'EtOO      ','CO2       ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,5,                       &
       (/'MeCOCH2OO ','NO        '/),                                   &
       (/'MeCO3     ','HCHO      ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                       &
       (/'NO        ','ISO2      '/),                                   &
       (/'NO2       ','MACR      ','HCHO      ','HO2       '/)),        &
       asad_flux_defn('RXN',50003,'B',.TRUE.,0,6,                       &
       (/'NO        ','MACRO2    '/),                                   &
       (/'NO2       ','MeCO3     ','HACET     ','CO        '/)),        &

! OH + inorganic acid reactions: sum into STASH section 50 item 4
       asad_flux_defn('RXN',50004,'B',.TRUE.,0,4,                       &
       (/'OH        ','HONO2     '/),                                   &
       (/'H2O       ','NO3       ','          ','          '/)),        &
       asad_flux_defn('RXN',50004,'B',.TRUE.,0,4,                       &
       (/'OH        ','HONO      '/),                                   &
       (/'H2O       ','NO2       ','          ','          '/)),        &

! OH + organic nitrate reactions: sum into STASH section 50 item 5
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                       &
       (/'OH        ','MeONO2    '/),                                   &
       (/'HCHO      ','NO2       ','H2O       ','          '/)),        &
       asad_flux_defn('RXN',50005,'B',.TRUE.,0,5,                       &
       (/'OH        ','NALD      '/),                                   &
       (/'HCHO      ','CO        ','NO2       ','          '/)),        &

! Organic nitrate photolysis: sum into STASH section 50 item 6
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,5,                       &
       (/'MeONO2    ','PHOTON    '/),                                   &
       (/'HO2       ','HCHO      ','NO2       ','          '/)),        &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,6,                       &
       (/'NALD      ','PHOTON    '/),                                   &
       (/'HCHO      ','CO        ','NO2       ','HO2       '/)),        &
       asad_flux_defn('RXN',50006,'J',.TRUE.,0,6,                       &
       (/'ISON      ','PHOTON    '/),                                   &
       (/'NO2       ','MACR      ','HCHO      ','HO2       '/)),        &

! OH + PAN-type reactions: sum into STASH section 50 item 7
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,5,                       &
       (/'OH        ','PAN       '/),                                   &
       (/'HCHO      ','NO2       ','H2O       ','          '/)),        &
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,5,                       &
       (/'OH        ','PPAN      '/),                                   &
       (/'MeCHO     ','NO2       ','H2O       ','          '/)),        &
       asad_flux_defn('RXN',50007,'B',.TRUE.,0,4,                       &
       (/'OH        ','MPAN      '/),                                   &
       (/'HACET     ','NO2       ','          ','          '/))         &
       /)


TYPE(asad_flux_defn), PARAMETER :: asad_trop_ox_budget_loss01(23) &
! Loss of Ox
       = (/ asad_flux_defn('RXN',50011,'B',.TRUE.,0,4,                  &
       (/'O(1D)     ','H2O       '/),                                   &
       (/'OH        ','OH        ','          ','          '/)),        &

! Minor reactions, should have negligble impact:
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                       &
       (/'O(1D)     ','CH4       '/),                                   &
       (/'OH        ','MeOO      ','          ','          '/)),        &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                       &
       (/'O(1D)     ','CH4       '/),                                   &
       (/'HCHO      ','H2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,5,                       &
       (/'O(1D)     ','CH4       '/),                                   &
       (/'HCHO      ','HO2       ','HO2       ','          '/)),        &

! Include with above. Each is loss of 2xOx, so include twice to sum
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                       &
       (/'O(3P)     ','O3        '/),                                   &
       (/'O2        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                       &
       (/'O(3P)     ','O3        '/),                                   &
       (/'O2        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                       &
       (/'O(3P)     ','NO2       '/),                                   &
       (/'NO        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50012,'B',.TRUE.,0,4,                       &
       (/'O(3P)     ','NO2       '/),                                   &
       (/'NO        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50013,'B',.TRUE.,0,4,                       &
       (/'HO2       ','O3        '/),                                   &
       (/'OH        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50014,'B',.TRUE.,0,4,                       &
       (/'OH        ','O3        '/),                                   &
       (/'HO2       ','O2        ','          ','          '/)),        &

! O3 + alkene:
! O3 + isoprene reactions
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                       &
       (/'O3        ','C5H8      '/),                                   &
       (/'MACR      ','HCHO      ','MACRO2    ','MeCO3     '/)),        &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,6,                       &
       (/'O3        ','C5H8      '/),                                   &
       (/'MeOO      ','HCOOH     ','CO        ','H2O2      '/)),        &
       asad_flux_defn('RXN',50015,'B',.TRUE.,0,4,                       &
       (/'O3        ','C5H8      '/),                                   &
       (/'HO2       ','OH        ','          ','          '/)),        &

! O3 + MACR reactions. ratb_defs specifies 2 different rates for each of these,
!  so need to select each one in turn
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,6,                       &
       (/'O3        ','MACR      '/),                                   &
       (/'MGLY      ','HCOOH     ','HO2       ','CO        '/)),        &
       asad_flux_defn('RXN',50015,'B',.TRUE.,1,4,                       &
       (/'O3        ','MACR      '/),                                   &
       (/'OH        ','MeCO3     ','          ','          '/)),        &
       asad_flux_defn('RXN',50015,'B',.TRUE.,2,6,                       &
       (/'O3        ','MACR      '/),                                   &
       (/'MGLY      ','HCOOH     ','HO2       ','CO        '/)),        &
       asad_flux_defn('RXN',50015,'B',.TRUE.,2,4,                       &
       (/'O3        ','MACR      '/),                                   &
       (/'OH        ','MeCO3     ','          ','          '/)),        &

! N2O5 + H20 reaction
       asad_flux_defn('RXN',50016,'B',.TRUE.,0,4,                       &
       (/'N2O5      ','H2O       '/),                                   &
       (/'HONO2     ','HONO2     ','          ','          '/)),        &
! Heterogenious reactions are not included yet
!RXN    50016     H       4       N2O5    H2O    HONO2    HONO2

! NO3 chemical loss
! sink of 2xOx:
       asad_flux_defn('RXN',50017,'J',.TRUE.,0,4,                       &
       (/'NO3       ','PHOTON    '/),                                   &
       (/'NO        ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50017,'J',.TRUE.,0,4,                       &
       (/'NO3       ','PHOTON    '/),                                   &
       (/'NO        ','O2        ','          ','          '/)),        &
! these are sinks of 1xOx
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                       &
       (/'HO2       ','NO3       '/),                                   &
       (/'OH        ','NO2       ','          ','          '/)),        &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                       &
       (/'OH        ','NO3       '/),                                   &
       (/'HO2       ','NO2       ','          ','          '/)),        &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                       &
       (/'MeOO      ','NO3       '/),                                   &
       (/'HO2       ','HCHO      ','NO2       ','          '/))         &
       /)

TYPE(asad_flux_defn), PARAMETER :: asad_trop_ox_budget_loss02(12) &
 = (/ asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                  &
 (/'EtOO      ','NO3       '/),                                   &
 (/'MeCHO     ','HO2       ','NO2       ','          '/)),        &
 asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                       &
 (/'MeCO3     ','NO3       '/),                                   &
 (/'MeOO      ','CO2       ','NO2       ','          '/)),        &
 asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                       &
 (/'n-PrOO    ','NO3       '/),                                   &
 (/'EtCHO     ','HO2       ','NO2       ','          '/)),        &
 asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                       &
 (/'i-PrOO    ','NO3       '/),                                   &
 (/'Me2CO     ','HO2       ','NO2       ','          '/)),        &
 asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                       &
 (/'EtCO3     ','NO3       '/),                                   &
 (/'EtOO      ','CO2       ','NO2       ','          '/)),        &
 asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                       &
 (/'MeCOCH2OO ','NO3       '/),                                   &
 (/'MeCO3     ','HCHO      ','NO2       ','          '/)),        &
 asad_flux_defn('RXN',50017,'B',.TRUE.,0,5,                       &
 (/'NO3       ','HCHO      '/),                                   &
 (/'HONO2     ','HO2       ','CO        ','          '/)),        &
 asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                       &
 (/'NO3       ','MeCHO     '/),                                   &
 (/'HONO2     ','MeCO3     ','          ','          '/)),        &
 asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                       &
 (/'NO3       ','EtCHO     '/),                                   &
 (/'HONO2     ','EtCO3     ','          ','          '/)),        &
 asad_flux_defn('RXN',50017,'B',.TRUE.,0,4,                       &
 (/'NO3       ','Me2CO     '/),                                   &
 (/'HONO2     ','MeCOCH2OO ','          ','          '/)),        &
! sink of 2xOx
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                       &
       (/'NO3       ','C5H8      '/),                                   &
       (/'ISON      ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50017,'B',.TRUE.,0,3,                       &
       (/'NO3       ','C5H8      '/),                                   &
       (/'ISON      ','          ','          ','          '/))         &
       /)

TYPE(asad_flux_defn), PARAMETER :: asad_trop_ox_budget_drydep(11) &
! Dry deposition:
! O3 dry deposition
       = (/asad_flux_defn('DEP',50021,'D',.TRUE.,0,1,                   &
       (/'O3        ','          '/),                                   &
       (/'          ','          ','          ','          '/)),&
! NOy dry deposition
! sink of 2xOx
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'NO3       ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'NO3       ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
! sink on 3xOx
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'HO2NO2    ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'HONO2     ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'PAN       ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'PPAN      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50022,'D',.TRUE.,0,1,                       &
       (/'MPAN      ','          '/),                                   &
       (/'          ','          ','          ','          '/))         &
       /)

TYPE(asad_flux_defn), PARAMETER :: asad_trop_ox_budget_wetdep(7)= &
! Wet deposition:
! NOy wet deposition
! sink of 2xOx
       (/ asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                    &
       (/'NO3       ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                       &
       (/'NO3       ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
! sink on 3xOx
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                       &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                       &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                       &
       (/'N2O5      ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                       &
       (/'HO2NO2    ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('DEP',50031,'W',.TRUE.,0,1,                       &
       (/'HONO2     ','          '/),                                   &
       (/'          ','          ','          ','          '/))         &
       /)

TYPE(asad_flux_defn), PARAMETER :: asad_trop_other_fluxes(16) =   &
! Extra fluxes of interest
       (/ asad_flux_defn('RXN',50042,'B',.TRUE.,0,3,                    &
       (/'NO3       ','C5H8      '/),                                   &
       (/'ISON      ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50043,'B',.TRUE.,0,3,                       &
       (/'NO        ','ISO2      '/),                                   &
       (/'ISON      ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50044,'T',.TRUE.,0,4,                       &
       (/'HO2       ','HO2       '/),                                   &
       (/'H2O2      ','O2        ','          ','          '/)),        &
       asad_flux_defn('RXN',50044,'B',.TRUE.,0,3,                       &
       (/'HO2       ','HO2       '/),                                   &
       (/'H2O2      ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                       &
       (/'HO2       ','MeOO      '/),                                   &
       (/'MeOOH     ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                       &
       (/'HO2       ','EtOO      '/),                                   &
       (/'EtOOH     ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                       &
       (/'HO2       ','MeCO3     '/),                                   &
       (/'MeCO3H    ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,4,                       &
       (/'HO2       ','MeCO3     '/),                                   &
       (/'MeCO2H    ','O3        ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                       &
       (/'HO2       ','n-PrOO    '/),                                   &
       (/'n-PrOOH   ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                       &
       (/'HO2       ','i-PrOO    '/),                                   &
       (/'i-PrOOH   ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,4,                       &
       (/'HO2       ','EtCO3     '/),                                   &
       (/'O2        ','EtCO3H    ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,4,                       &
       (/'HO2       ','EtCO3     '/),                                   &
       (/'EtCO2H    ','O3        ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                       &
       (/'HO2       ','MeCOCH2OO '/),                                   &
       (/'MeCOCH2OOH','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                       &
       (/'HO2       ','ISO2      '/),                                   &
       (/'ISOOH     ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50045,'B',.TRUE.,0,3,                       &
       (/'HO2       ','MACRO2    '/),                                   &
       (/'MACROOH   ','          ','          ','          '/)),        &
       asad_flux_defn('RXN',50046,'T',.TRUE.,0,4,                       &
       (/'OH        ','NO2       '/),                                   &
       (/'HONO2     ','m         ','          ','          '/))         &
       /)

TYPE(asad_flux_defn), PARAMETER :: asad_general_interest(8) = (/  &
! For CH4 lifetime
       asad_flux_defn('RXN',50041,'B',.TRUE.,0,4,                       &
       (/'OH        ','CH4       '/),                                   &
       (/'H2O       ','MeOO      ','          ','          '/)),        &
! Ozone STE
       asad_flux_defn('STE',50051,'X',.TRUE.,0,1,                       &
       (/'O3        ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
! Ozone tendency in the troposphere
       asad_flux_defn('NET',50052,'X',.TRUE.,0,1,                       &
       (/'O3        ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
! Ozone in the troposphere
       asad_flux_defn('OUT',50053,'X',.TRUE.,0,1,                       &
       (/'O3        ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
! Ozone tendency - whole atmos
       asad_flux_defn('NET',50054,'X',.FALSE.,0,1,                      &
       (/'O3        ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
! Tropospheric mass
       asad_flux_defn('MAS',50061,'X',.TRUE.,0,1,                       &
       (/'X         ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
! Tropospheric mask
       asad_flux_defn('TPM',50062,'X',.TRUE.,0,1,                       &
       (/'X         ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
! Total atmospheric mass
       asad_flux_defn('MAS',50063,'X',.FALSE.,0,1,                      &
       (/'X         ','          '/),                                   &
       (/'          ','          ','          ','          '/))         &
       /)


TYPE(asad_flux_defn), PARAMETER :: asad_trop_co_budget(21) = (/   &
! CO budget
! CO loss
       asad_flux_defn('RXN',50071,'B',.TRUE.,0,3,                       &
       (/'OH        ','CO        '/),                                   &
       (/'HO2       ','          ','          ','          '/)),        &
! CO prod - bimol rxns
! HCHO + OH/NO3
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                       &
       (/'NO3       ','HCHO      '/),                                   &
       (/'HONO2     ','HO2       ','CO        ','          '/)),        &
       asad_flux_defn('RXN',50072,'B',.TRUE.,0,5,                       &
       (/'OH        ','HCHO      '/),                                   &
       (/'H2O       ','HO2       ','CO        ','          '/)),        &
! MGLY + OH/NO3
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,4,                       &
       (/'OH        ','MGLY      '/),                                   &
       (/'MeCO3     ','CO        ','          ','          '/)),        &
       asad_flux_defn('RXN',50073,'B',.TRUE.,0,5,                       &
       (/'NO3       ','MGLY      '/),                                   &
       (/'MeCO3     ','CO        ','HONO2     ','          '/)),        &
! O3 + MACR/ISOP & OTHER FLUXES
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                       &
       (/'O3        ','C5H8      '/),                                   &
       (/'MeOO      ','HCOOH     ','CO        ','H2O2      '/)),        &
       asad_flux_defn('RXN',50074,'B',.TRUE.,1,6,                       &
       (/'O3        ','MACR      '/),                                   &
       (/'MGLY      ','HCOOH     ','HO2       ','CO        '/)),        &
       asad_flux_defn('RXN',50074,'B',.TRUE.,2,6,                       &
       (/'O3        ','MACR      '/),                                   &
       (/'MGLY      ','HCOOH     ','HO2       ','CO        '/)),        &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                       &
       (/'NO        ','MACRO2    '/),                                   &
       (/'NO2       ','MeCO3     ','HACET     ','CO        '/)),        &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,6,                       &
       (/'MACRO2    ','MACRO2    '/),                                   &
       (/'HACET     ','MGLY      ','HCHO      ','CO        '/)),        &
       asad_flux_defn('RXN',50074,'B',.TRUE.,0,5,                       &
       (/'OH        ','NALD      '/),                                   &
       (/'HCHO      ','CO        ','NO2       ','          '/)),        &
! CO prod - photol rxns
! HCHO photolysis: RADICAL
       asad_flux_defn('RXN',50075,'J',.TRUE.,0,5,                       &
       (/'HCHO      ','PHOTON    '/),                                   &
       (/'HO2       ','HO2       ','CO        ','          '/)),        &
! HCHO photolysis: MOLECULAR
       asad_flux_defn('RXN',50076,'J',.TRUE.,0,4,                       &
       (/'HCHO      ','PHOTON    '/),                                   &
       (/'H2        ','CO        ','          ','          '/)),        &
! MGLY photolysis
       asad_flux_defn('RXN',50077,'J',.TRUE.,0,5,                       &
       (/'MGLY      ','PHOTON    '/),                                   &
       (/'MeCO3     ','CO        ','HO2       ','          '/)),        &
! OTHER CO PROD PHOTOLYSIS RXNS
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                       &
       (/'MeCHO     ','PHOTON    '/),                                   &
       (/'MeOO      ','HO2       ','CO        ','          '/)),        &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,4,                       &
       (/'MeCHO     ','PHOTON    '/),                                   &
       (/'CH4       ','CO        ','          ','          '/)),        &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,5,                       &
       (/'EtCHO     ','PHOTON    '/),                                   &
       (/'EtOO      ','HO2       ','CO        ','          '/)),        &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                       &
       (/'MACR      ','PHOTON    '/),                                   &
       (/'MeCO3     ','HCHO      ','CO        ','HO2       '/)),        &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                       &
       (/'MACROOH   ','PHOTON    '/),                                   &
       (/'HACET     ','CO        ','MGLY      ','HCHO      '/)),        &
       asad_flux_defn('RXN',50078,'J',.TRUE.,0,6,                       &
       (/'NALD      ','PHOTON    '/),                                   &
       (/'HCHO      ','CO        ','NO2       ','HO2       '/)),        &
! CO drydep
       asad_flux_defn('DEP',50079,'D',.TRUE.,0,1,                       &
       (/'CO        ','          '/),                                   &
       (/'          ','          ','          ','          '/))         &
       /)

TYPE(asad_flux_defn), PARAMETER :: asad_lightning_diags(6) = (/   &
! WARNING: Lightning emissions are calculated as NO2, but emitted as NO
       asad_flux_defn('EMS',50081,'L',.FALSE.,0,1,                      &
       (/'NO        ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('LGT',50082,'T',.FALSE.,0,1,                      &
       (/'X         ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('LGT',50083,'G',.FALSE.,0,1,                      &
       (/'X         ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('LGT',50084,'C',.FALSE.,0,1,                      &
       (/'X         ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('LGT',50085,'N',.FALSE.,0,1,                      &
       (/'X         ','          '/),                                   &
       (/'          ','          ','          ','          '/)),        &
       asad_flux_defn('LIN',50086,'X',.FALSE.,0,1,                      &
       (/'X         ','          '/),                                   &
       (/'          ','          ','          ','          '/))         &
       /)

TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                       &
                                     asad_strat_oh_prod(31) = (/ &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,5,                      &
(/'Cl        ','HOCl      '/),                                   &
(/'Cl        ','Cl        ','OH        ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'O(3P)     ','HBr       '/),                                  &
(/'OH        ','Br        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'O(3P)     ','HCl       '/),                                  &
(/'OH        ','Cl        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'O(3P)     ','HOCl      '/),                                  &
(/'OH        ','ClO       ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'O(1D)     ','HBr       '/),                                  &
(/'OH        ','Br        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'O(1D)     ','HCl       '/),                                  &
(/'OH        ','Cl        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'O(1D)     ','MeBr      '/),                                  &
(/'OH        ','Br        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'Cl        ','HO2       '/),                                  &
(/'ClO       ','OH        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'Cl        ','HO2       '/),                                  &
(/'ClO       ','OH        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'NO3       ','HO2       '/),                                  &
(/'NO2       ','OH        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'O3        ','HO2       '/),                                  &
(/'O2        ','OH        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'O(3P)     ','HO2       '/),                                  &
(/'OH        ','O2        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'O(1D)     ','CH4       '/),                                  &
(/'OH        ','MeOO      ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'O(1D)     ','H2        '/),                                  &
(/'OH        ','H         ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'O(1D)     ','H2O       '/),                                  &
(/'OH        ','OH        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'O(1D)     ','H2O       '/),                                  &
(/'OH        ','OH        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'H         ','HO2       '/),                                  &
(/'OH        ','OH        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'H         ','HO2       '/),                                  &
(/'OH        ','OH        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'H         ','O3        '/),                                  &
(/'OH        ','O2        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'H         ','NO2       '/),                                  &
(/'OH        ','NO        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'O(3P)     ','H2        '/),                                  &
(/'OH        ','H         ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'O(3P)     ','H2O2      '/),                                  &
(/'OH        ','HO2       ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,5,                      &
 (/'O(3P)     ','HCHO      '/),                                  &
(/'OH        ','CO        ','HO2       ','          '/)),        &
asad_flux_defn('RXN',50091,'B',.FALSE.,0,4,                      &
 (/'NO        ','HO2       '/),                                  &
(/'NO2       ','OH        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                      &
  (/'HOBr      ','PHOTON    '/),                                 &
(/'OH        ','Br        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                      &
  (/'HOCl      ','PHOTON    '/),                                 &
(/'OH        ','Cl        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                      &
  (/'HONO2     ','PHOTON    '/),                                 &
(/'NO2       ','OH        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                      &
  (/'HO2NO2    ','PHOTON    '/),                                 &
(/'NO3       ','OH        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                      &
  (/'H2O       ','PHOTON    '/),                                 &
(/'OH        ','H         ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,4,                      &
  (/'H2O2      ','PHOTON    '/),                                 &
(/'OH        ','OH        ','          ','          '/)),        &
asad_flux_defn('RXN',50091,'J',.FALSE.,0,5,                      &
  (/'MeOOH     ','PHOTON    '/),                                 &
(/'HCHO      ','HO2       ','OH        ','          '/))         &
/)

! OH loss reactions
TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                         &
                                       asad_strat_oh_loss(26) = (/ &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
  (/'BrO       ','OH        '/),                                   &
  (/'Br        ','HO2       ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
  (/'ClO       ','OH        '/),                                   &
  (/'Cl        ','HO2       ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
  (/'OH        ','ClO       '/),                                   &
  (/'HCl       ','O2        ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
  (/'OH        ','HBr       '/),                                   &
  (/'H2O       ','Br        ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
  (/'OH        ','HCl       '/),                                   &
  (/'H2O       ','Cl        ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
  (/'OH        ','HOCl      '/),                                   &
  (/'ClO       ','H2O       ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
  (/'OH        ','OClO      '/),                                   &
  (/'HOCl      ','O2        ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
  (/'OH        ','ClONO2    '/),                                   &
  (/'HOCl      ','NO3       ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
  (/'OH        ','MeBr      '/),                                   &
  (/'Br        ','H2O       ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
   (/'OH        ','O3        '/),                                  &
  (/'HO2       ','O2        ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
  (/'O(3P)     ','OH        '/),                                   &
  (/'O2        ','H         ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
  (/'OH        ','CH4       '/),                                   &
  (/'H2O       ','MeOO      ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,3,                      &
  (/'OH        ','CO        '/),                                   &
  (/'HO2       ','          ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
   (/'OH        ','NO3       '/),                                  &
  (/'HO2       ','NO2       ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'T',.FALSE.,0,4,                      &
  (/'OH        ','NO2       '/),                                   &
  (/'HONO2     ','m         ','          ','          '/)),        &
  Asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
  (/'OH        ','HONO2     '/),                                   &
  (/'NO3       ','H2O       ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,5,                      &
  (/'OH        ','HO2NO2    '/),                                   &
  (/'H2O       ','NO2       ','O2        ','          '/)),        &
  asad_flux_defn('RXN',50092,'T',.FALSE.,0,4,                      &
   (/'OH        ','OH        '/),                                  &
  (/'H2O2      ','m         ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'T',.FALSE.,0,4,                      &
  (/'OH        ','OH        '/),                                   &
  (/'H2O2      ','m         ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
  (/'OH        ','OH        '/),                                   &
  (/'H2O       ','O(3P)     ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
  (/'OH        ','OH        '/),                                   &
  (/'H2O       ','O(3P)     ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,3,                      &
  (/'OH        ','HO2       '/),                                   &
  (/'H2O       ','          ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
   (/'OH        ','H2O2      '/),                                  &
  (/'H2O       ','HO2       ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,5,                      &
   (/'OH        ','HCHO      '/),                                  &
  (/'H2O       ','CO        ','HO2       ','          '/)),        &
  Asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
  (/'OH        ','H2        '/),                                   &
  (/'H2O       ','H02       ','          ','          '/)),        &
  asad_flux_defn('RXN',50092,'B',.FALSE.,0,4,                      &
  (/'OH        ','MeOOH     '/),                                   &
  (/'MeOO      ','H2O       ','          ','          '/))         &
  /)

! Simple strat ozone budget
TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                       &
                                   asad_strat_o3_budget(20) = (/ &
! production
asad_flux_defn('RXN',50101,'J',.FALSE.,0,4,                      &
(/'O2        ','PHOTON    '/),                                   &
(/'O(3P)     ','O(3P)     ','          ','          '/)),        &
asad_flux_defn('RXN',50101,'J',.FALSE.,0,4,                      &
(/'O2        ','PHOTON    '/),                                   &
(/'O(3P)     ','O(1D)     ','          ','          '/)),        &
asad_flux_defn('RXN',50102,'B',.FALSE.,0,4,                      &
(/'HO2       ','NO        '/),                                   &
(/'OH        ','NO2       ','          ','          '/)),        &
asad_flux_defn('RXN',50103,'B',.FALSE.,0,5,                      &
(/'MeOO      ','NO        '/),                                   &
(/'HO2       ','HCHO      ','NO2       ','          '/)),        &
asad_flux_defn('RXN',50104,'B',.FALSE.,0,4,                      &
(/'OH        ','HONO2     '/),                                   &
(/'H2O       ','NO3       ','          ','          '/)),        &
! loss
asad_flux_defn('RXN',50111,'J',.FALSE.,0,5,                      &
(/'Cl2O2     ','PHOTON    '/),                                   &
(/'Cl        ','Cl        ','O2        ','          '/)),        &
asad_flux_defn('RXN',50112,'B',.FALSE.,0,4,                      &
(/'BrO       ','ClO       '/),                                   &
(/'Br        ','OClO      ','          ','          '/)),        &
asad_flux_defn('RXN',50113,'B',.FALSE.,0,4,                      &
(/'HO2       ','O3        '/),                                   &
(/'OH        ','O2        ','          ','          '/)),        &
asad_flux_defn('RXN',50114,'B',.FALSE.,0,4,                      &
(/'ClO       ','HO2       '/),                                   &
(/'HOCl      ','O2        ','          ','          '/)),        &
asad_flux_defn('RXN',50115,'B',.FALSE.,0,4,                      &
(/'BrO       ','HO2       '/),                                   &
(/'HOBr      ','O2        ','          ','          '/)),        &
asad_flux_defn('RXN',50116,'B',.FALSE.,0,4,                      &
(/'O(3P)     ','ClO       '/),                                   &
(/'Cl        ','O2        ','          ','          '/)),        &
asad_flux_defn('RXN',50117,'B',.FALSE.,0,4,                      &
(/'O(3P)     ','NO2       '/),                                   &
(/'NO        ','O2        ','          ','          '/)),        &
asad_flux_defn('RXN',50118,'B',.FALSE.,0,4,                      &
(/'O(3P)     ','HO2       '/),                                   &
(/'OH        ','O2        ','          ','          '/)),        &
asad_flux_defn('RXN',50119,'B',.FALSE.,0,4,                      &
(/'O3        ','H         '/),                                   &
(/'OH        ','O2        ','          ','          '/)),        &
asad_flux_defn('RXN',50120,'B',.FALSE.,0,4,                      &
(/'O(3P)     ','O3        '/),                                   &
(/'O2        ','O2        ','          ','          '/)),        &
asad_flux_defn('RXN',50121,'J',.FALSE.,0,4,                      &
(/'NO3       ','PHOTON    '/),                                   &
(/'NO        ','O2        ','          ','          '/)),        &
asad_flux_defn('RXN',50122,'B',.FALSE.,0,4,                      &
(/'O(1D)     ','H2O       '/),                                   &
(/'OH        ','OH        ','          ','          '/)),        &
asad_flux_defn('RXN',50123,'B',.FALSE.,0,5,                      &
(/'HO2       ','NO3       '/),                                   &
(/'OH        ','NO2       ','O2        ','          '/)),        &
asad_flux_defn('RXN',50124,'B',.FALSE.,0,4,                      &
(/'OH        ','NO3       '/),                                   &
(/'HO2       ','NO2       ','          ','          '/)),        &
asad_flux_defn('RXN',50125,'B',.FALSE.,0,5,                      &
(/'NO3       ','HCHO      '/),                                   &
(/'HONO2     ','HO2       ','CO        ','          '/))         &
/)

TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                        &
                                      asad_strat_o3_misc(15) = (/ &
 asad_flux_defn('DEP',50131,'D',.FALSE.,0,1,                      &
 (/'O3        ','          '/),                                   &
 (/'          ','          ','          ','          '/)),        &
 asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                      &
 (/'NO3       ','          '/),                                   &
 (/'          ','          ','          ','          '/)),        &
 asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                      &
 (/'NO3       ','          '/),                                   &
 (/'          ','          ','          ','          '/)),        &
 asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                      &
 (/'N2O5      ','          '/),                                   &
 (/'          ','          ','          ','          '/)),        &
 asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                      &
 (/'N2O5      ','          '/),                                   &
 (/'          ','          ','          ','          '/)),        &
 asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                      &
 (/'N2O5      ','          '/),                                   &
 (/'          ','          ','          ','          '/)),        &
 asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                      &
 (/'HO2NO2    ','          '/),                                   &
 (/'          ','          ','          ','          '/)),        &
 asad_flux_defn('DEP',50132,'D',.FALSE.,0,1,                      &
 (/'HONO2     ','          '/),                                   &
 (/'          ','          ','          ','          '/)),        &
 asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                      &
 (/'NO3       ','          '/),                                   &
 (/'          ','          ','          ','          '/)),        &
 asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                      &
 (/'NO3       ','          '/),                                   &
 (/'          ','          ','          ','          '/)),        &
 asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                      &
 (/'N2O5      ','          '/),                                   &
 (/'          ','          ','          ','          '/)),        &
 asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                      &
 (/'N2O5      ','          '/),                                   &
 (/'          ','          ','          ','          '/)),        &
 asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                      &
 (/'N2O5      ','          '/),                                   &
 (/'          ','          ','          ','          '/)),        &
 asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                      &
 (/'HO2NO2    ','          '/),                                   &
 (/'          ','          ','          ','          '/)),        &
 asad_flux_defn('DEP',50133,'W',.FALSE.,0,1,                      &
 (/'HONO2     ','          '/),                                   &
 (/'          ','          ','          ','          '/))         &
 /)

! RC(O)O2 + NO2 --> *PAN production reactions, output to 50-248 
TYPE(asad_flux_defn), PARAMETER :: asad_rco2no2_pan_prod(3) = (/ &
asad_flux_defn('RXN',50248,'T',.FALSE.,0,4,                      & ! T016
(/'MeCO3     ','NO2       '/),                                   &
(/'PAN       ','m         ','          ','          '/) ),       &
asad_flux_defn('RXN',50248,'T',.FALSE.,0,4,                      & ! T018
(/'EtCO3     ','NO2       '/),                                   &
(/'PPAN      ','m         ','          ','          '/) ),       &
asad_flux_defn('RXN',50248,'T',.FALSE.,0,4,                      & ! T020
(/'MACRO2    ','NO2       '/),                                   &
(/'MPAN      ','m         ','          ','          '/) )        &
  /)

! RO2 + HO2 type reactions, output to 50-249 
TYPE(asad_flux_defn), PARAMETER :: asad_ro2ho2_reacn(13) = (/    &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,4,                      & ! B048
(/'HO2       ','EtCO3     '/),                                   &
(/'O2        ','EtCO3H    ','          ','          '/) ),       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,4,                      & ! B049
(/'HO2       ','EtCO3     '/),                                   &
(/'O3        ','EtCO2H    ','          ','          '/) ),       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B050
(/'HO2       ','EtOO      '/),                                   &
(/'EtOOH     ','          ','          ','          '/) ),       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B052
(/'HO2       ','ISO2      '/),                                   &
(/'ISOOH     ','          ','          ','          '/) ),       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B053
(/'HO2       ','MACRO2    '/),                                   &
(/'MACROOH   ','          ','          ','          '/) ),       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,4,                      & ! B054
(/'HO2       ','MeCO3     '/),                                   &
(/'MeCO2H    ','O3        ','          ','          '/) ),       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B055
(/'HO2       ','MeCO3     '/),                                   &
(/'MeCO3H    ','          ','          ','          '/) ),       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,4,                      & ! B056
(/'HO2       ','MeCO3     '/),                                   &
(/'OH        ','MeOO      ','          ','          '/) ),       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B057
(/'HO2       ','MeCOCH2OO '/),                                   &
(/'MeCOCH2OOH','          ','          ','          '/) ),       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,4,                      & ! B058
(/'HO2       ','MeOO      '/),                                   &
(/'HCHO      ','          ','          ','          '/) ),       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,4,                      & ! B059
(/'HO2       ','MeOO      '/),                                   &
(/'MeOOH     ','          ','          ','          '/) ),       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B063
(/'HO2       ','i-PrOO    '/),                                   &
(/'i-PrOOH   ','          ','          ','          '/) ),       &
asad_flux_defn('RXN',50249,'B',.FALSE.,0,3,                      & ! B064
(/'HO2       ','n-PrOO    '/),                                   &
(/'n-PrOOH   ','          ','          ','          '/) )        &
  /)

! RO2 + NO3 type reactions, output to 50-250 
TYPE(asad_flux_defn), PARAMETER :: asad_ro2no3_reacn(7) = (/     &
asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B039
(/'EtCO3     ','NO3       '/),                                   &
(/'EtOO      ','CO2       ','NO2       ','          '/) ),       &
asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B042
(/'EtOO      ','NO3       '/),                                   &
(/'MeCHO     ','HO2       ','NO2       ','          '/) ),       &
asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B072
(/'MeCO3     ','NO3       '/),                                   &
(/'MeOO      ','CO2       ','NO2       ','          '/) ),       &
asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B074
(/'MeCOCH2OO ','NO3       '/),                                   &
(/'MeCO3     ','HCHO      ','NO2       ','          '/) ),       &
asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B081
(/'MeOO      ','NO3       '/),                                   &
(/'HO2       ','HCHO      ','NO2       ','          '/) ),       &
asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B195
(/'i-PrOO    ','NO3       '/),                                   &
(/'Me2CO     ','HO2       ','NO2       ','          '/) ),       &
asad_flux_defn('RXN',50250,'B',.FALSE.,0,5,                      & ! B197
(/'n-PrOO    ','NO3       '/),                                   &
(/'EtCHO     ','HO2       ','NO2       ','          '/) )        &
  /)

! RO2 + RO2 type reactions, output to 50-251 
TYPE(asad_flux_defn), PARAMETER :: asad_ro2ro2_reacn(8) = (/      &
asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B040
(/'EtOO      ','MeCO3     '/),                                   &
(/'MeCHO     ','HO2       ','MeOO      ','          '/) ),       &
asad_flux_defn('RXN',50251,'B',.FALSE.,0,6,                      & ! B065
(/'ISO2      ','ISO2      '/),                                   &
(/'MACR      ','MACR      ','HCHO      ','HO2       '/) ),       &
asad_flux_defn('RXN',50251,'B',.FALSE.,0,6,                      & ! B066
(/'MACRO2    ','MACRO2    '/),                                   &
(/'HACET     ','MGLY      ','HCHO      ','CO        '/) ),       &
asad_flux_defn('RXN',50251,'B',.FALSE.,0,3,                      & ! B067
(/'MACRO2    ','MACRO2    '/),                                   &
(/'HO2       ','          ','          ','          '/) ),       &
asad_flux_defn('RXN',50251,'B',.FALSE.,0,5,                      & ! B075
(/'MeOO      ','MeCO3     '/),                                   &
(/'HO2       ','HCHO      ','MeOO      ','          '/) ),       &
asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B076
(/'MeOO      ','MeCO3     '/),                                   &
(/'MeCO2H    ','HCHO      ','          ','          '/) ),       &
asad_flux_defn('RXN',50251,'B',.FALSE.,0,6,                      & ! B077
(/'MeOO      ','MeOO      '/),                                   &
(/'HO2       ','HO2       ','HCHO      ','HCHO      '/) ),       &
asad_flux_defn('RXN',50251,'B',.FALSE.,0,4,                      & ! B078
(/'MeOO      ','MeOO      '/),                                   &
(/'MeOH      ','HCHO      ','          ','          '/) )        &
  /)

! Methane Oxidation reactions, output to 50-247
TYPE(asad_flux_defn), PARAMETER :: asad_ch4_oxidn(6) = (/        &
asad_flux_defn('RXN',50247,'J',.FALSE.,0,4,                      & ! RATJ_T
(/'CH4       ','PHOTON    '/),                                   &
(/'MeOO      ','H         ','          ','          '/)),        &
asad_flux_defn('RXN',50247,'B',.FALSE.,0,4,                      & ! B016
(/'Cl        ','CH4       '/),                                   &
(/'HCl       ','MeOO      ','          ','          '/) ),       &
asad_flux_defn('RXN',50247,'B',.FALSE.,0,4,                      & ! B101
(/'O(1D)     ','CH4       '/),                                   &
(/'HCHO      ','H2        ','          ','          '/) ),       &
asad_flux_defn('RXN',50247,'B',.FALSE.,0,5,                      & ! B102
(/'O(1D)     ','CH4       '/),                                   &
(/'HCHO      ','HO2       ','HO2       ','          '/) ),       &
asad_flux_defn('RXN',50247,'B',.FALSE.,0,4,                      & ! B103
(/'O(1D)     ','CH4       '/),                                   &
(/'OH        ','MeOO      ','          ','          '/) ),       &
asad_flux_defn('RXN',50247,'B',.FALSE.,0,4,                      & ! B145
(/'OH        ','CH4       '/),                                   &
(/'H2O       ','MeOO      ','          ','          '/) )        &
  /)

! Chemical Production of O(1D), output to 50-254
TYPE(asad_flux_defn), PARAMETER :: asad_o1d_prod(3) = (/         &
asad_flux_defn('RXN',50254,'J',.FALSE.,0,4,                      &
(/'O3        ','PHOTON    '/),                                   &
(/'O2        ','O(1D)     ','          ','          '/)),        &
asad_flux_defn('RXN',50254,'J',.FALSE.,0,4,                      &
(/'O2        ','PHOTON    '/),                                   &
(/'O(3P)     ','O(1D)     ','          ','          '/)),        &
asad_flux_defn('RXN',50254,'J',.FALSE.,0,4,                      &
(/'N2O       ','PHOTON    '/),                                   &
(/'N2        ','O(1D)     ','          ','          '/))         &
  /)

! Tropospheric sulphur chemistry for online oxidants
TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                       &
                            asad_aerosol_chem_online(30) = (/    &
asad_flux_defn('RXN',50135,'B',.FALSE.,0,4,                      &
(/'MSIA      ','OH        '/),                                   &
(/'SO2       ','MSA       ','          ','          '/)),        &
asad_flux_defn('RXN',50136,'B',.FALSE.,0,4,                      &
(/'DMSO      ','OH        '/),                                   &
(/'MSIA      ','SO2       ','          ','          '/)),        &
asad_flux_defn('RXN',50137,'B',.FALSE.,0,6,                      &
(/'DMS       ','Cl        '/),                                   &
(/'SO2       ','DMSO      ','HCl       ','ClO       '/)),        &
asad_flux_defn('RXN',50138,'B',.FALSE.,0,4,                      &
(/'DMS       ','BrO       '/),                                   &
(/'DMSO      ','Br        ','          ','          '/)),        &
asad_flux_defn('RXN',50139,'B',.FALSE.,0,3,                      &
(/'MSIA      ','O3        '/),                                   &
(/'MSA       ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50140,'B',.FALSE.,0,5,                      &
(/'DMS       ','OH        '/),                                   &
(/'SO2       ','MeOO      ','HCHO      ','          '/)),        &
asad_flux_defn('RXN',50141,'B',.FALSE.,0,5,                      &
(/'DMS       ','OH        '/),                                   &
(/'SO2       ','DMSO      ','MeOO      ','          '/)),        &
asad_flux_defn('RXN',50142,'B',.FALSE.,0,6,                      &
(/'DMS       ','NO3       '/),                                   &
(/'SO2       ','HONO2     ','MeOO      ','HCHO      '/)),        &
asad_flux_defn('RXN',50143,'B',.FALSE.,0,4,                      &
(/'DMSO      ','OH        '/),                                   &
(/'SO2       ','MSA       ','          ','          '/)),        &
asad_flux_defn('RXN',50144,'B',.FALSE.,0,4,                      &
(/'CS2       ','OH        '/),                                   &
(/'SO2       ','COS       ','          ','          '/)),        &
asad_flux_defn('RXN',50145,'B',.FALSE.,0,3,                      &
(/'H2S       ','OH        '/),                                   &
(/'SO2       ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50146,'B',.FALSE.,0,3,                      &
(/'COS       ','OH        '/),                                   &
(/'SO2       ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50147,'B',.FALSE.,0,3,                      &
(/'Monoterp  ','OH        '/),                                   &
(/'Sec_Org   ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50148,'B',.FALSE.,0,3,                      &
(/'Monoterp  ','O3        '/),                                   &
(/'Sec_Org   ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50149,'B',.FALSE.,0,3,                      &
(/'Monoterp  ','NO3       '/),                                   &
(/'Sec_Org   ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50150,'T',.FALSE.,0,4,                      &
(/'SO2       ','OH        '/),                                   &
(/'HO2       ','H2SO4     ','          ','          '/)),        &
asad_flux_defn('RXN',50151,'H',.FALSE.,0,3,                      &
(/'SO2       ','H2O2      '/),                                   &
(/'NULL0     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50152,'H',.FALSE.,0,3,                      &
(/'SO2       ','O3        '/),                                   &
(/'NULL1     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50153,'H',.FALSE.,0,3,                      &
(/'SO2       ','O3        '/),                                   &
(/'NULL2     ','          ','          ','          '/)),        &
asad_flux_defn('DEP',50154,'D',.TRUE.,0,1,                       &
(/'SO2       ','          '/),                                   &
(/'          ','          ','          ','          '/)),        &
asad_flux_defn('DEP',50155,'W',.TRUE.,0,1,                       &
(/'SO2       ','          '/),                                   &
(/'          ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50331,'H',.FALSE.,0,3,                      &
(/'DMS       ','O3        '/),                                   &
(/'NULL3     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50332,'H',.FALSE.,0,3,                      &
(/'MSIA      ','O3        '/),                                   &
(/'NULL4     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50333,'H',.FALSE.,0,3,                      &
(/'MSIA      ','O3        '/),                                   &
(/'NULL5     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50334,'H',.FALSE.,0,3,                      &
(/'SO2       ','HOBr      '/),                                   &
(/'NULL6     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50335,'H',.FALSE.,0,3,                      &
(/'SO2       ','HOBr      '/),                                   &
(/'NULL7     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50336,'B',.FALSE.,0,5,                      &
(/'DMS       ','OH        '/),                                   &
(/'SO2       ','MeOO      ','HCHO      ','          '/)),        &
asad_flux_defn('RXN',50337,'B',.FALSE.,0,5,                      &
(/'DMS       ','OH        '/),                                   &
(/'SO2       ','DMSO      ','MeOO      ','          '/)),        &
asad_flux_defn('RXN',50338,'B',.FALSE.,0,6,                      &
(/'DMS       ','NO3       '/),                                   &
(/'SO2       ','HONO2     ','MeOO      ','HCHO      '/)),        &
asad_flux_defn('RXN',50339,'B',.FALSE.,0,3,                      &
(/'DMS       ','O3        '/),                                   &
(/'SO2       ','          ','          ','          '/))         &

/)

! Tropospheric sulphur chemistry for Offline oxidants chemistry
TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                       &
                                   asad_aerosol_chem_offline(30) &
                                                         = (/    &
asad_flux_defn('RXN',50135,'B',.FALSE.,0,4,                      &
(/'MSIA      ','OH        '/),                                   &
(/'SO2       ','MSA       ','          ','          '/)),        &
asad_flux_defn('RXN',50136,'B',.FALSE.,0,4,                      &
(/'DMSO      ','OH        '/),                                   &
(/'MSIA      ','SO2       ','          ','          '/)),        &
asad_flux_defn('RXN',50137,'B',.FALSE.,0,6,                      &
(/'DMS       ','Cl        '/),                                   &
(/'SO2       ','DMSO      ','HCl       ','ClO       '/)),        &
asad_flux_defn('RXN',50138,'B',.FALSE.,0,4,                      &
(/'DMS       ','BrO       '/),                                   &
(/'DMSO      ','Br        ','          ','          '/)),        &
asad_flux_defn('RXN',50139,'B',.FALSE.,0,3,                      &
(/'MSIA      ','O3        '/),                                   &
(/'MSA       ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50140,'B',.FALSE.,0,3,                      &
(/'DMS       ','OH        '/),                                   &
(/'SO2       ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50141,'B',.FALSE.,0,4,                      &
(/'DMS       ','OH        '/),                                   &
(/'SO2       ','DMSO      ','          ','          '/)),        &
asad_flux_defn('RXN',50142,'B',.FALSE.,0,3,                      &
(/'DMS       ','NO3       '/),                                   &
(/'SO2       ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50143,'B',.FALSE.,0,3,                      &
(/'DMSO      ','OH        '/),                                   &
(/'SO2       ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50144,'B',.FALSE.,0,4,                      &
(/'CS2       ','OH        '/),                                   &
(/'SO2       ','COS       ','          ','          '/)),        &
asad_flux_defn('RXN',50145,'B',.FALSE.,0,3,                      &
(/'H2S       ','OH        '/),                                   &
(/'SO2       ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50146,'B',.FALSE.,0,3,                      &
(/'COS       ','OH        '/),                                   &
(/'SO2       ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50147,'B',.FALSE.,0,3,                      &
(/'Monoterp  ','OH        '/),                                   &
(/'Sec_Org   ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50148,'B',.FALSE.,0,3,                      &
(/'Monoterp  ','O3        '/),                                   &
(/'Sec_Org   ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50149,'B',.FALSE.,0,3,                      &
(/'Monoterp  ','NO3       '/),                                   &
(/'Sec_Org   ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50150,'T',.FALSE.,0,3,                      &
(/'SO2       ','OH        '/),                                   &
(/'H2SO4     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50151,'H',.FALSE.,0,3,                      &
(/'SO2       ','H2O2      '/),                                   &
(/'NULL0     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50152,'H',.FALSE.,0,3,                      &
(/'SO2       ','O3        '/),                                   &
(/'NULL1     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50153,'H',.FALSE.,0,3,                      &
(/'SO2       ','O3        '/),                                   &
(/'NULL2     ','          ','          ','          '/)),        &
asad_flux_defn('DEP',50154,'D',.TRUE.,0,1,                       &
(/'SO2       ','          '/),                                   &
(/'          ','          ','          ','          '/)),        &
asad_flux_defn('DEP',50155,'W',.TRUE.,0,1,                       &
(/'SO2       ','          '/),                                   &
(/'          ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50331,'H',.FALSE.,0,3,                      &
(/'DMS       ','O3        '/),                                   &
(/'NULL3     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50332,'H',.FALSE.,0,3,                      &
(/'MSIA      ','O3        '/),                                   &
(/'NULL4     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50333,'H',.FALSE.,0,3,                      &
(/'MSIA      ','O3        '/),                                   &
(/'NULL5     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50334,'H',.FALSE.,0,3,                      &
(/'SO2       ','HOBr      '/),                                   &
(/'NULL6     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50335,'H',.FALSE.,0,3,                      &
(/'SO2       ','HOBr      '/),                                   &
(/'NULL7     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50336,'B',.FALSE.,0,5,                      &
(/'DMS       ','OH        '/),                                   &
(/'SO2       ','MeOO      ','HCHO      ','          '/)),        &
asad_flux_defn('RXN',50337,'B',.FALSE.,0,5,                      &
(/'DMS       ','OH        '/),                                   &
(/'SO2       ','DMSO      ','MeOO      ','          '/)),        &
asad_flux_defn('RXN',50338,'B',.FALSE.,0,6,                      &
(/'DMS       ','NO3       '/),                                   &
(/'SO2       ','HONO2     ','MeOO      ','HCHO      '/)),        &
asad_flux_defn('RXN',50339,'B',.FALSE.,0,3,                      &
(/'DMS       ','O3        '/),                                   &
(/'SO2       ','          ','          ','          '/))         &

/)

! Water production: sum into section 50, item 238
! Reaction identifiers from ukca_chem_strattrop are shown
TYPE(asad_flux_defn), PARAMETER :: asad_h2o_budget(38) = (/      &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'H         ','HO2       '/),                           & ! B044
(/'O(3P)     ','H2O       ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'MeBr      ','OH        '/),                           & ! B070
(/'Br        ','H2O       ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','C2H6      '/),                           & ! B141 
(/'H2O       ','EtOO      ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','C3H8      '/),                           & ! B142
(/'i-PrOO    ','H2O       ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','C3H8      '/),                           & ! B143
(/'n-PrOO    ','H2O       ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','CH4       '/),                           & ! B145
(/'H2O       ','MeOO      ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','EtCHO     '/),                           & ! B150
(/'H2O       ','EtCO3     ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','EtOOH     '/),                           & ! B151
(/'H2O       ','EtOO      ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                      &
(/'OH        ','EtOOH     '/),                           & ! B152
(/'H2O       ','MeCHO     ','OH        ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','H2        '/),                           & ! B153
(/'H2O       ','HO2       ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','H2O2      '/),                           & ! B154
(/'H2O       ','HO2       ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','HBr       '/),                           & ! B156
(/'H2O       ','Br        ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                      &
(/'OH        ','HCHO      '/),                           & ! B157
(/'H2O       ','HO2       ','CO        ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','HCl       '/),                           & ! B159
(/'H2O       ','Cl        ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,3,                      &
(/'OH        ','HO2       '/),                           & ! B160
(/'H2O       ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','HONO2     '/),                           & ! B164
(/'H2O       ','NO3       ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                      &
(/'OH        ','HO2NO2    '/),                           & ! B161
(/'H2O       ','NO2       ','O2        ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','HOCl      '/),                           & ! B162
(/'ClO       ','H2O       ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','HONO      '/),                           & ! B163
(/'H2O       ','NO2       ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,1,4,                      &
(/'OH        ','Me2CO     '/),                           & ! B172
(/'H2O       ','MeCOCH2OO ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,2,4,                      &
(/'OH        ','Me2CO     '/),                           & ! B173
(/'H2O       ','MeCOCH2OO ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','MeCHO     '/),                           & ! B174
(/'H2O       ','MeCO3     ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','MeCOCH2OOH'/),                           & ! B177
(/'H2O       ','MeCOCH2OO ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                      &
(/'OH        ','MeONO2    '/),                           & ! B180
(/'HCHO      ','NO2       ','H2O       ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                      &
(/'OH        ','MeOOH     '/),                           & ! B181
(/'H2O       ','HCHO      ','OH        ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','MeOOH     '/),                           & ! B182
(/'H2O       ','MeOO      ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','OH        '/),                           & ! B187
(/'H2O       ','O(3P)     ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                      &
(/'OH        ','PAN       '/),                           & ! B188
(/'HCHO      ','NO2       ','H2O       ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                      &
(/'OH        ','PPAN      '/),                           & ! B189
(/'MeCHO     ','NO2       ','H2O       ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','i-PrOOH   '/),                           & ! B191
(/'i-PrOO    ','H2O       ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,4,                      &
(/'OH        ','n-PrOOH   '/),                           & ! B192
(/'n-PrOO    ','H2O       ','          ','          '/)),        &
asad_flux_defn('RXN',50238,'B',.FALSE.,0,5,                      &
(/'OH        ','n-PrOOH   '/),                           & ! B193
(/'EtCHO     ','H2O       ','OH        ','          '/)),        &
asad_flux_defn('RXN',50238,'H',.FALSE.,0,5,                      &
(/'HOCl      ','HCl       '/),                           & ! Het PSC
(/'Cl        ','Cl        ','H2O       ','          '/)),        &

! Water loss: sum into section 50, item 239
asad_flux_defn('RXN',50239,'B',.FALSE.,0,4,                      &
(/'O(1D)     ','H2O       '/),                            & ! B106
(/'OH        ','OH        ','          ','          '/)),        &
asad_flux_defn('RXN',50239,'B',.FALSE.,0,4,                      &
(/'N2O5      ','H2O       '/),                            & ! B085
(/'HONO2     ','HONO2     ','          ','          '/)),        &
asad_flux_defn('RXN',50239,'J',.FALSE.,0,4,                      &
(/'H2O       ','PHOTON    '/),                            & ! Photol
(/'OH        ','H         ','          ','          '/)),        &
asad_flux_defn('RXN',50239,'H',.FALSE.,0,4,                      &
(/'ClONO2    ','H2O       '/),                            & ! Het PSC
(/'HOCl      ','HONO2     ','          ','          '/)),        &
asad_flux_defn('RXN',50239,'H',.FALSE.,0,4,                      &
(/'N2O5      ','H2O       '/),                            & ! Het PSC
(/'HONO2     ','HONO2     ','          ','          '/))         &
/)

! Strat-Trop sulphur chemistry (contains explicit SO3)
TYPE(asad_flux_defn), PARAMETER, PUBLIC ::                       &
                         asad_aerosol_chem_strattrop(30) = (/    &
asad_flux_defn('RXN',50135,'B',.FALSE.,0,4,                      &
(/'MSIA      ','OH        '/),                                   &
(/'SO2       ','MSA       ','          ','          '/)),        &
asad_flux_defn('RXN',50136,'B',.FALSE.,0,4,                      &
(/'DMSO      ','OH        '/),                                   &
(/'MSIA      ','SO2       ','          ','          '/)),        &
asad_flux_defn('RXN',50137,'B',.FALSE.,0,6,                      &
(/'DMS       ','Cl        '/),                                   &
(/'SO2       ','DMSO      ','HCl       ','ClO       '/)),        &
asad_flux_defn('RXN',50138,'B',.FALSE.,0,4,                      &
(/'DMS       ','BrO       '/),                                   &
(/'DMSO      ','Br        ','          ','          '/)),        &
asad_flux_defn('RXN',50139,'B',.FALSE.,0,3,                      &
(/'MSIA      ','O3        '/),                                   &
(/'MSA       ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50140,'B',.FALSE.,0,5,                      &
(/'DMS       ','OH        '/),                                   &
(/'SO2       ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50141,'B',.FALSE.,0,5,                      &
(/'DMS       ','OH        '/),                                   &
(/'SO2       ','MSA       ','          ','          '/)),        &
asad_flux_defn('RXN',50142,'B',.FALSE.,0,6,                      &
(/'DMS       ','NO3       '/),                                   &
(/'SO2       ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50143,'B',.FALSE.,0,4,                      &
(/'DMSO      ','OH        '/),                                   &
(/'SO2       ','MSA       ','          ','          '/)),        &
asad_flux_defn('RXN',50144,'B',.FALSE.,0,4,                      &
(/'CS2       ','OH        '/),                                   &
(/'SO2       ','COS       ','          ','          '/)),        &
asad_flux_defn('RXN',50145,'B',.FALSE.,0,4,                      &
(/'H2S       ','OH        '/),                                   &
(/'SO2       ','H2O       ','          ','          '/)),        &
asad_flux_defn('RXN',50146,'B',.FALSE.,0,3,                      &
(/'COS       ','OH        '/),                                   &
(/'SO2       ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50147,'B',.FALSE.,0,3,                      &
(/'Monoterp  ','OH        '/),                                   &
(/'Sec_Org   ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50148,'B',.FALSE.,0,3,                      &
(/'Monoterp  ','O3        '/),                                   &
(/'Sec_Org   ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50149,'B',.FALSE.,0,3,                      &
(/'Monoterp  ','NO3       '/),                                   &
(/'Sec_Org   ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50150,'T',.FALSE.,0,4,                      &
(/'SO2       ','OH        '/),                                   &
(/'H2SO4     ','HO2       ','          ','          '/)),        &
asad_flux_defn('RXN',50151,'H',.FALSE.,0,3,                      &
(/'SO2       ','H2O2      '/),                                   &
(/'NULL0     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50152,'H',.FALSE.,0,3,                      &
(/'SO2       ','O3        '/),                                   &
(/'NULL1     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50153,'H',.FALSE.,0,3,                      &
(/'SO2       ','O3        '/),                                   &
(/'NULL2     ','          ','          ','          '/)),        &
asad_flux_defn('DEP',50154,'D',.TRUE.,0,1,                       &
(/'SO2       ','          '/),                                   &
(/'          ','          ','          ','          '/)),        &
asad_flux_defn('DEP',50155,'W',.TRUE.,0,1,                       &
(/'SO2       ','          '/),                                   &
(/'          ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50331,'H',.FALSE.,0,3,                      &
(/'DMS       ','O3        '/),                                   &
(/'NULL3     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50332,'H',.FALSE.,0,3,                      &
(/'MSIA      ','O3        '/),                                   &
(/'NULL4     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50333,'H',.FALSE.,0,3,                      &
(/'MSIA      ','O3        '/),                                   &
(/'NULL5     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50334,'H',.FALSE.,0,3,                      &
(/'SO2       ','HOBr      '/),                                   &
(/'NULL6     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50335,'H',.FALSE.,0,3,                      &
(/'SO2       ','HOBr      '/),                                   &
(/'NULL7     ','          ','          ','          '/)),        &
asad_flux_defn('RXN',50336,'B',.FALSE.,0,5,                      &
(/'DMS       ','OH        '/),                                   &
(/'SO2       ','MeOO      ','HCHO      ','          '/)),        &
asad_flux_defn('RXN',50337,'B',.FALSE.,0,5,                      &
(/'DMS       ','OH        '/),                                   &
(/'SO2       ','DMSO      ','MeOO      ','          '/)),        &
asad_flux_defn('RXN',50338,'B',.FALSE.,0,6,                      &
(/'DMS       ','NO3       '/),                                   &
(/'SO2       ','HONO2     ','MeOO      ','HCHO      '/)),        &
asad_flux_defn('RXN',50339,'B',.FALSE.,0,3,                      &
(/'DMS       ','O3        '/),                                   &
(/'SO2       ','          ','          ','          '/))         &
/)

TYPE(asad_flux_defn), PUBLIC :: asad_aerosol_chem(30)

PUBLIC :: asad_load_default_fluxes

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ASAD_FLUX_DAT'

CONTAINS

! #####################################################################

SUBROUTINE asad_load_default_fluxes
! To load flux information into asad_chemical_fluxes array for use in
!  asad flux diagnostics scheme. Note that both the total number of
!  fluxes and the number in each section (currently: reaction fluxes;
!  photolytic fluxes; dry deposition; and wet deposition) are limited.
!  Reactions etc are either taken from the definition arrays above
!  or are specified in a user-STASHmaster file, and read in using
!  the chemical definition arrays.
!  NOT currently handling the Strat-trop exchange as this should be
!  in section 38.
!  Called from UKCA_MAIN1.

USE ukca_option_mod,     ONLY: l_ukca_offline, l_ukca_offline_be, &
                               l_ukca_strattrop
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim
USE umPrintMgr

USE umPrintMgr

USE version_mod, ONLY: nitemp

IMPLICIT NONE

INTEGER :: n_max_diags     ! max No of diagnostics
INTEGER :: idiag           ! item number for diagnostic
INTEGER :: i, j, k         ! counters
INTEGER :: ierr            ! error code

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_LOAD_DEFAULT_FLUXES'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

n_max_diags = nitemp

! Pick up correct version of the aerosol chemistry
IF (l_ukca_offline .OR. l_ukca_offline_be) THEN
  asad_aerosol_chem = asad_aerosol_chem_offline
ELSE IF (l_ukca_strattrop) THEN
  asad_aerosol_chem = asad_aerosol_chem_strattrop
ELSE
  asad_aerosol_chem = asad_aerosol_chem_online
END IF

! Use pre-determined reaction list

ALLOCATE(asad_chemical_fluxes(n_chemical_fluxes))
asad_chemical_fluxes=(/        &
   asad_trop_ox_budget_prod,   & ! 20
   asad_trop_ox_budget_loss01, & ! 23 43
   asad_trop_ox_budget_loss02, & ! 12 55
   asad_trop_ox_budget_drydep, & ! 11 66
   asad_trop_ox_budget_wetdep, & ! 7  73
   asad_trop_other_fluxes,     & ! 16 89
   asad_general_interest,      & ! 8  97
   asad_trop_co_budget,        & ! 21 118
   asad_lightning_diags,       & ! 6  124
   asad_strat_oh_prod,         & ! 31 155
   asad_strat_oh_loss,         & ! 26 181
   asad_strat_o3_budget,       & ! 20 201
   asad_strat_o3_misc,         & ! 15 216
   asad_aerosol_chem,          & ! 16 232
   asad_rco2no2_pan_prod,      & ! 3  235
   asad_ro2ho2_reacn,          & ! 13 248
   asad_ro2no3_reacn,          & ! 7  255
   asad_ro2ro2_reacn,          & ! 8  263
   asad_ch4_oxidn,             & ! 6  269
   asad_o1d_prod,              & ! 3  272
   asad_h2o_budget             & ! 38 307
/)

IF (printstatus > PrStatus_Normal) THEN  
  WRITE(umMessage,'(A)') 'LIST OF DEFINED ASAD DIAGNOSTICS:'
  CALL umPrint(umMessage,src='asad_flux_dat')
  WRITE(umMessage,'(A)') '========================================='
  CALL umPrint(umMessage,src='asad_flux_dat')
  WRITE(umMessage,'(A)') ' '
  CALL umPrint(umMessage,src='asad_flux_dat')
  WRITE(umMessage,'(A)') 'N=STASH=TYPE=R1=======R2========P1========P2'//     &
             '========P3========P4========No.='
  CALL umPrint(umMessage,src='asad_flux_dat')
  j = 0
  DO i = 1, SIZE(asad_chemical_fluxes)
    IF (asad_chemical_fluxes(i)%stash_number /= imdi) THEN
      WRITE(umMessage,'(I3,1x,I5,1x,A3,1x,6A10,I3)') i,           &
                  asad_chemical_fluxes(i)%stash_number,           &
                  asad_chemical_fluxes(i)%diag_type,              &
                  asad_chemical_fluxes(i)%reactants(1),           &
                  asad_chemical_fluxes(i)%reactants(2),           &
                  asad_chemical_fluxes(i)%products(1),            &
                  asad_chemical_fluxes(i)%products(2),            &
                  asad_chemical_fluxes(i)%products(3),            &
                  asad_chemical_fluxes(i)%products(4),            &
                  asad_chemical_fluxes(i)%num_species
      CALL umPrint(umMessage,src='asad_flux_dat')
      j = j + 1
    END IF
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_load_default_fluxes

! ######################################################################
FUNCTION prcount(i)
! To sum total number of reactants (set to 2) and products

USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN) :: i                      ! Reaction index

INTEGER, PARAMETER  :: n_reactants = 2        ! Default no of reactants
INTEGER             :: prcount                ! No of reactants + products
INTEGER             :: n                      ! counter

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRCOUNT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

n = n_reactants
IF (asad_chemical_fluxes(i)%products(1) /= blank0) n = n+1
IF (asad_chemical_fluxes(i)%products(2) /= blank0) n = n+1
IF (asad_chemical_fluxes(i)%products(3) /= blank0) n = n+1
IF (asad_chemical_fluxes(i)%products(4) /= blank0) n = n+1

prcount = n

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION prcount


END MODULE asad_flux_dat
