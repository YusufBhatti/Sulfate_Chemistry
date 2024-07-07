! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  A routine to set up values such as the number of reactions used
!  by each chemistry scheme
!
! Method:
!  i_ukca_chem, l_ukca_chem_aero and l_ukca_trophet are read
!  in from namelists
!  This routine uses these values to set the number or reactants, and the
!  number of each type of reaction, plus the number of dry and wet
!  deposited species
!  Also includes error checking for schemes which have not yet
!  been set up
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
!   Language:  FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------
!
MODULE ukca_setup_chem_mod

IMPLICIT NONE

! Item numbers for HO2 and OH when non-transported prognostics
! Ensure consistency by using same parameter here and in ukca_ntp_mod
INTEGER, PARAMETER :: i_ho2_be = 993
INTEGER, PARAMETER :: i_oh_be  = 995

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_SETUP_CHEM_MOD'

CONTAINS

SUBROUTINE ukca_setup_chem

USE ukca_option_mod, ONLY:  i_ukca_chem, l_ukca_chem_aero,              &
   l_ukca_achem, l_ukca_aerchem, l_ukca_raq, l_ukca_trop,               &
   l_ukca_tropisop, l_ukca_strattrop, l_ukca_stratcfc,                  &
   l_ukca_strat, l_ukca_chem, l_ukca_het_psc,                           &
   l_ukca_mode, l_ukca_advh2o, l_ukca_trophet, ukca_int_method,         &
   jpctr, jpspec, jpbk, jptk, jphk, jppj, jpdd, jpdw, jpnr,             &
   i_ukca_dms_flux, chem_timestep,                                      &
   l_ukca_offline, l_ukca_offline_be, l_ukca_nr_aqchem, l_ukca_rado3,   &
   l_ukca_radch4, l_ukca_radn2o, l_ukca_radf11, l_ukca_radf12,          &
   l_ukca_radf113, l_ukca_radf22, l_ukca_classic_hetchem,               &
   l_ukca_raqaero

USE ukca_chem_schemes_mod, ONLY: int_method_be_explicit, int_method_nr, &
   i_ukca_chem_off, i_ukca_chem_trop, i_ukca_chem_raq,                  &
   i_ukca_chem_tropisop, i_ukca_chem_strattrop, i_ukca_chem_strat,      &
   i_ukca_chem_offline, i_ukca_chem_offline_be

USE ukca_d1_defs,    ONLY: ukca_item_sulpc
USE run_aerosol_mod, ONLY: l_sulpc_dms
USE dms_flux_mod_4a, ONLY: i_liss_merlivat,i_wanninkhof,i_nightingale,i_blomquist !08/11/22 Added
USE missing_data_mod, ONLY: imdi
USE ereport_mod,     ONLY: ereport
USE umPrintMgr
USE parkind1,        ONLY: jprb, jpim
USE yomhook,         ONLY: lhook, dr_hook

USE ukca_chem_offline, ONLY: ndry_offline, nwet_offline
USE ukca_chem_aer,     ONLY: ndry_aer_aer, nwet_aer_aer
USE ukca_chem1_dat,    ONLY: ndry_trop, nwet_trop
USE ukca_chem_raq,     ONLY: ndry_raq, nwet_raq
USE ukca_chem_raqaero_mod, ONLY: ndry_raqaero, nwet_raqaero

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

CHARACTER (LEN=errormessagelength) :: cmessage        ! Error message
INTEGER            :: errcode         ! Variable passed to ereport

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SETUP_CHEM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


SELECT CASE (i_ukca_chem)

CASE (i_ukca_chem_off)
  ! Chemistry off completely
  l_ukca_chem     = .FALSE.
  ukca_int_method = 0
  jpctr           = 0
  jpspec          = 0
  jpbk            = 0
  jptk            = 0
  jppj            = 0
  jphk            = 0
  jpdd            = 0
  jpdw            = 0

CASE (i_ukca_chem_trop)

  l_ukca_chem     = .TRUE.
  l_ukca_trop     = .TRUE.
  ukca_int_method = int_method_BE_explicit
  IF (l_ukca_chem_aero) THEN
    ! Tropospheric chemistry with Aerosols (BE)
    l_ukca_aerchem  = .TRUE.
    jpctr           = 33
    jpspec          = 53
    jpbk            = 88
    jptk            = 14
    jppj            = 20
    jphk            = 0
    jpdd            = ndry_aer_aer
    jpdw            = nwet_aer_aer

  ELSE
    ! Tropopspheric chemistry (BE) without aerosols
    jpctr            = 26
    jpspec           = 46
    jpbk             = 88
    jptk             = 14
    jppj             = 20
    jphk             = 0
    jpdd             = ndry_trop
    jpdw             = nwet_trop

  END IF

CASE (i_ukca_chem_raq)
  ! Regional Air Quality (BE)
  l_ukca_chem     = .TRUE.
  l_ukca_raq      = .TRUE.
  l_ukca_raqaero  = .FALSE.
  ukca_int_method = int_method_BE_explicit
  jpctr           = 40
  jpspec          = 58
  jpbk            = 113
  jptk            = 12
  jppj            = 23
  jphk            = 0   ! By default zero heterogeneous reactions
  jpdd            = ndry_raq
  jpdw            = nwet_raq

  ! Increase number of heterogeneous reactions if needed. For the moment
  ! only the heterog. hydrolysis of N2O5 is used, but more reactions
  ! can be added in the future.
  IF ( l_ukca_classic_hetchem ) THEN
    jphk = jphk + 2
  END IF

  IF (l_ukca_chem_aero) THEN
    IF (l_ukca_classic_hetchem) THEN
      cmessage='RAQ-AERO is not compatibile with the CLASSIC ' //            &
               'heterogeneous chemistry'
      errcode=113
      CALL ereport('ukca_setup_chem',errcode,cmessage)      
    END IF
    l_ukca_raq     = .FALSE.
    l_ukca_raqaero = .TRUE.
    jpctr          = jpctr + 8
    jpspec         = jpspec + 8
    jpbk           = jpbk + 10
    jptk           = jptk + 2
    jppj           = jppj + 0
    jphk           = jphk + 2
    jpdd           = ndry_raqaero
    jpdw           = nwet_raqaero

  END IF

CASE (i_ukca_chem_offline_be)
! Offline oxidants aerosol precursor chemistry with backward-Euler solver
  l_ukca_chem     = .TRUE.
  l_ukca_offline_be = .TRUE.
  ukca_int_method = int_method_BE_explicit
  jpctr           = 7
  jpspec          = 11
  jpbk            = 9
  jptk            = 1
  jppj            = 0
  jphk            = 3
  jpdd            = ndry_offline
  jpdw            = nwet_offline

  IF (l_ukca_rado3 .OR. l_ukca_radch4 .OR. l_ukca_radn2o .OR.         &
      l_ukca_radf11 .OR. l_ukca_radf12 .OR. l_ukca_radf113 .OR.       &
      l_ukca_radf22) THEN
    CALL umPrint( 'Offline chemistry does not support any radiative ' &
      //'feedbacks from UKCA trace gases', src='ukca_setup_chem_mod')
    cmessage='Unsupported option choice'
    errcode=ABS(i_ukca_chem)
    CALL ereport('ukca_setup_chem',errcode,cmessage)
  END IF

CASE (i_ukca_chem_tropisop)
  l_ukca_chem     = .TRUE.
  l_ukca_tropisop = .TRUE.
  ukca_int_method = int_method_nr

  ! If aerosol chemistry on, add additional reactions
  IF (l_ukca_chem_aero) THEN
    l_ukca_achem     = .TRUE.
    l_ukca_nr_aqchem = .TRUE.

  ELSE
     ! Can't have trophet reactions without aerosol chemistry
    IF (l_ukca_trophet) THEN
      CALL umPrint( 'trop het chem requires aerosol chemistry on', &
          src='ukca_setup_chem_mod')
      cmessage='Unsupported option choice'
      errcode=ABS(i_ukca_chem)
      CALL ereport('ukca_setup_chem',errcode,cmessage)
    END IF
  END IF

CASE (i_ukca_chem_strattrop)
  l_ukca_chem      = .TRUE.
  l_ukca_strattrop = .TRUE.
  l_ukca_advh2o    = .TRUE.
  ukca_int_method = int_method_nr

  ! If aerosol chemistry on, add additional reactions
  IF (l_ukca_chem_aero) THEN
    l_ukca_achem     =.TRUE.
    l_ukca_nr_aqchem = .TRUE.

  ELSE
     ! Can't have trophet reactions without aerosol chemistry
    IF (l_ukca_trophet) THEN
      CALL umPrint( 'trop het chem requires aerosol chemistry on', &
          src='ukca_setup_chem_mod')
      cmessage='Unsupported option choice'
      errcode=ABS(i_ukca_chem)
      CALL ereport('ukca_setup_chem',errcode,cmessage)
    END IF
  END IF

CASE (i_ukca_chem_strat)
  l_ukca_chem     = .TRUE.
  l_ukca_strat    = .TRUE.
  l_ukca_advh2o    = .TRUE.
  ukca_int_method = int_method_nr
!  ! If aerosol chemistry on, add additional reactions
  IF (l_ukca_chem_aero) THEN
    l_ukca_achem     =.TRUE.
    l_ukca_nr_aqchem = .TRUE.

  END IF

  IF (l_ukca_trophet) THEN
    CALL umPrint( 'Strat chem does not support trop het chem', &
        src='ukca_setup_chem_mod')
    cmessage='Unsupported option choice'
    errcode=ABS(i_ukca_chem)
    CALL ereport('ukca_setup_chem',errcode,cmessage)
  END IF

CASE (i_ukca_chem_offline)
  l_ukca_chem     = .TRUE.
  l_ukca_offline  = .TRUE.
  l_ukca_nr_aqchem = .TRUE.
  l_ukca_advh2o    = .FALSE.
  ukca_int_method = int_method_nr

  IF (l_ukca_rado3 .OR. l_ukca_radch4 .OR. l_ukca_radn2o .OR.         &
      l_ukca_radf11 .OR. l_ukca_radf12 .OR. l_ukca_radf113 .OR.       &
      l_ukca_radf22) THEN
    CALL umPrint( 'Offline chemistry does not support any radiative ' &
      //'feedbacks from UKCA trace gases', src='ukca_setup_chem_mod')
    cmessage='Unsupported option choice'
    errcode=ABS(i_ukca_chem)
    CALL ereport('ukca_setup_chem',errcode,cmessage)
  END IF

CASE DEFAULT
  ! If we get here, we don't know how to deal with this chemistry scheme,
  ! stop
  WRITE(umMessage,'(A,1X,I6)') 'Unknown chemistry scheme: ',i_ukca_chem
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  cmessage='Unknown chemistry scheme'
  errcode=ABS(i_ukca_chem)
  CALL ereport('ukca_setup_chem',errcode,cmessage)

END SELECT

! calculate jpnr from the sum of the reactions
jpnr = jpbk + jptk + jppj + jphk

! Set item numbers for UKCA oxidants used in the CLASSIC sulphur cycle:
! Settings for B-E solvers with OH and HO2 as non-transported species
!       001 : O3 MASS MIXING RATIO AFTER TSTEP
!       007 : HONO2 MASS MIXING RATIO AFTER TSTEP
!       008 : H2O2 MASS MIXING RATIO AFTER TSTEP
!       995 : OH MASS MIXING RATIO  AFTER TSTEP
!       993 : HO2 MASS MIXING RATIO AFTER TSTEP

! Settings for N-R solvers with OH and HO2 as transported species
!        81 : OH MASS MIXING RATIO  AFTER TSTEP
!        82 : HO2 MASS MIXING RATIO AFTER TSTEP

IF (ukca_int_method == int_method_nr) THEN
  ukca_item_sulpc(:) = (/1,7,8,81,82/)   ! N-R solver
ELSE
  ukca_item_sulpc(:) = (/1,7,8,i_oh_be,i_ho2_be/)
END IF


! Error checking
IF (ukca_int_method == int_method_nr .OR.                     &
    i_ukca_chem ==  i_ukca_chem_offline_be) THEN
  ! The chemical timestep must be set for N-R and offline oxidants (BE) schemes
  IF (chem_timestep == imdi) THEN
    cmessage = ' The chemical timestep has not been set'
    errcode = 1
    CALL ereport('ukca_setup_chem',errcode,cmessage)
  END IF
END IF

IF (ukca_int_method  == int_method_BE_explicit .AND.           &
      i_ukca_chem /= i_ukca_chem_offline_be .AND.              &
      i_ukca_chem /= i_ukca_chem_raq) THEN 
    ! If aerosol chemistry or tropospheric het chem on abort - not supported
  IF (l_ukca_chem_aero .OR. l_ukca_trophet) THEN
    CALL umPrint( 'This BE scheme does not support add-on aerosol '  &
       // 'chemistry or tropospheric heterogeneous chemistry',       &
       src='ukca_setup_chem_mod')
    cmessage='Unsupported option combination'
    errcode=ABS(i_ukca_chem)
    CALL ereport('ukca_setup_chem',errcode,cmessage)
  END IF
END IF

! Tropospheric heterogeneous chemistry requires GLOMAP-mode
IF (l_ukca_trophet .AND. .NOT. l_ukca_mode) THEN
  CALL umPrint( 'Need Glomap MODE on for ' &
     // 'tropospheric heterogeneous chemistry',src='ukca_setup_chem_mod')
  cmessage='Unsupported option combination'
  errcode=1
  CALL ereport('ukca_setup_chem',errcode,cmessage)
END IF

! Check chemistry supports GLOMAP-mode
IF (l_ukca_mode) THEN

   ! Newton Raphson
   ! Require that l_ukca_achem true or offline chemistry
  IF (ukca_int_method  == int_method_nr) THEN
    IF (.NOT. l_ukca_nr_aqchem) THEN
      CALL umPrint( 'Need aerosol chemistry on for ' &
     // 'Glomap MODE',src='ukca_setup_chem_mod')
      cmessage='Unsupported option combination'
      errcode=2
      CALL ereport('ukca_setup_chem',errcode,cmessage)
    END IF

    ! Backward Euler - only Trop + Aerosols, offline and RAQ-Aero support
    ! GLOMAP-mode
  ELSE
    IF (.NOT. (l_ukca_aerchem .OR. l_ukca_offline_be .OR. l_ukca_raqaero)) THEN 
      CALL umPrint( 'Need aerosol chemistry on for ' &
     // 'Glomap MODE',src='ukca_setup_chem_mod')
      cmessage='Unsupported option combination'
      errcode=3
      CALL ereport('ukca_setup_chem',errcode,cmessage)
    END IF
  END IF

END IF

! Check specified DMS scheme logicals, only if CLASSIC DMS Off
! DMS flux routine not called from UKCA if CLASSIC is ON
IF ( l_ukca_chem_aero .AND. .NOT. l_sulpc_dms ) THEN  
  SELECT CASE (i_ukca_dms_flux )
  CASE (i_liss_merlivat)
    !    Do nothing
  CASE (i_wanninkhof)
    !    Do nothing
  CASE (i_nightingale)
    !    Do nothing
  CASE (i_blomquist)
    !    Do nothing
  CASE DEFAULT
    ! If not set or set to unknown value, throw up an error
    CALL umPrint( 'CLASSIC DMS is OFF but UKCA DMS scheme not selected ' &
      // 'i_ukca_dms_flux should be 1,2,3,4') ! YAB 08/11 added 4
    errcode = 4
    cmessage = 'RUN_UKCA: DMS flux scheme not specified'
    CALL ereport('UKCA_SETUP_CHEM_MOD',errcode,cmessage)
  END SELECT
END IF              ! chem_aero and no sulpc_dms

! Print the values we have set
IF (PrintStatus >= PrStatus_Diag) THEN
  CALL umPrint(       'Logicals and parameters set by ukca_setup_chem:', &
      src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_chem      = ', l_ukca_chem
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_trop      = ', l_ukca_trop
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_achem   = ',   l_ukca_achem
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_nr_aqchem = ', l_ukca_nr_aqchem
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_aerchem   = ', l_ukca_aerchem
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_tropisop  = ', l_ukca_tropisop
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_strattrop = ', l_ukca_strattrop
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_strat     = ', l_ukca_strat
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,L7)') 'l_ukca_offline   = ', l_ukca_offline
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'ukca_int_method  = ', ukca_int_method
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jpctr            = ', jpctr
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jpspec           = ', jpspec
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jpbk             = ', jpbk
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jptk             = ', jptk
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jppj             = ', jppj
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jphk             = ', jphk
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jpnr             = ', jpnr
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jpdd             = ', jpdd
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I5)') 'jpdw             = ', jpdw
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
  WRITE(umMessage,'(A,1X,I2)') 'i_ukca_dms_flux  = ', i_ukca_dms_flux
  CALL umPrint(umMessage,src='ukca_setup_chem_mod')
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ukca_setup_chem

END MODULE ukca_setup_chem_mod
