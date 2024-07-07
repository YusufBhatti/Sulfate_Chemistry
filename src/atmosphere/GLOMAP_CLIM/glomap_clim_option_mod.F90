! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Module to hold GLOMAP-mode climatology variables in 
!          RUN_GLOMAP_AEROCLIM namelist
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
MODULE glomap_clim_option_mod

USE errormessagelength_mod, ONLY: &
    errormessagelength

USE filenamelength_mod,     ONLY: &
    filenamelength

USE missing_data_mod,       ONLY: &
    imdi

USE parkind1,               ONLY: &
    jpim,                         &
    jprb

USE yomhook,                ONLY: &
    dr_hook,                      &
    lhook

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

! Declarations for GLOMAP_CLIM sub-model
! -----------------------------------------------------------------------------

! Namelist items:

! Namelist logicals:
LOGICAL :: l_glomap_mode_clim    = .FALSE. ! True when GLOMAP climatology is on

LOGICAL :: l_glomap_clim_arg_act = .FALSE. ! True when using Abdul-Razzak and
                                           !  Ghan Activation Method

LOGICAL :: l_glomap_clim_aie1    = .FALSE. ! True when first aerosol indirect
                                           !  effect required (on radiation)

LOGICAL :: l_glomap_clim_aie2    = .FALSE. ! True when second aerosol indirect
                                           !  effect required (on precip.)

LOGICAL :: l_glomap_clim_radaer  = .FALSE. ! True when GLOMAP climatology is 
                                           ! used in radaer

LOGICAL :: l_glomap_clim_radaer_sustrat=.FALSE. ! Use H2SO4 for stratospheric 
                                                ! sulphate

! Namelist integers:
INTEGER :: i_glomap_clim_setup = imdi   ! Controls aerosol scheme
                                        ! e.g. SUSSOCBC_5MODE when ==2
                                        ! SUSS_4mode (1) not currently
                                        ! available
                                        ! SUSSBCOCDU_7mode (8) not currently
                                        ! available

! Namelist characters:

! RADAER lookup tables and optical properties namelists.
CHARACTER(LEN=filenamelength) :: gclmaclw = 'gclmaclw is unset' 
                               !  Aitken + Insoluble acc mode (LW)
CHARACTER(LEN=filenamelength) :: gclmacsw = 'gclmacsw is unset' 
                               !  Aitken + Insoluble acc mode (SW)
CHARACTER(LEN=filenamelength) :: gclmanlw = 'gclmanlw is unset' 
                               !  Soluble accum mode (LW)
CHARACTER(LEN=filenamelength) :: gclmansw = 'gclmansw is unset' 
                               !  Soluble accum mode (SW)
CHARACTER(LEN=filenamelength) :: gclmcrlw = 'gclmcrlw is unset' 
                               !  Coarse mode (LW)
CHARACTER(LEN=filenamelength) :: gclmcrsw = 'gclmcrsw is unset' 
                               !  Coarse mode (SW)
CHARACTER(LEN=filenamelength) :: gclmprec = 'gclmprec is unset' 
                               !  Precomputed values

!Define the RUN_GLOMAP_AEROCLIM namelist
NAMELIST/RUN_GLOMAP_AEROCLIM/ l_glomap_mode_clim,                             &
                              l_glomap_clim_arg_act,                          &
                              l_glomap_clim_aie1,                             &
                              l_glomap_clim_aie2,                             &
                              l_glomap_clim_radaer,                           &
                              l_glomap_clim_radaer_sustrat,                   &
                              i_glomap_clim_setup,                            &
                              gclmaclw, gclmacsw, gclmanlw, gclmansw,         &
                              gclmcrlw, gclmcrsw, gclmprec

!===========================================================================
! start of items not set in namelist
!===========================================================================
! Glomap-mode integers:
INTEGER, PARAMETER :: i_gc_sussocbc_5mode = 2 ! CASE(i_glomap_clim_setup == 2)
!===========================================================================
! end of items not set in namelist
!===========================================================================

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GLOMAP_CLIM_OPTION_MOD'

CONTAINS

!---------------------------------------------------------------------------
SUBROUTINE print_nlist_run_glomap_clim ()

USE umPrintMgr, ONLY: &
    umPrint,          &
    umMessage

IMPLICIT NONE

! Local variables
REAL(KIND=jprb)             :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_GLOMAP_CLIM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

WRITE(umMessage,'(A)')    'Contents of namelist run_glomap_aeroclim'
CALL umPrint(umMessage, src=ModuleName)

WRITE(umMessage,'(A,L1)') ' l_glomap_mode_clim = ', l_glomap_mode_clim
CALL umPrint(umMessage, src=ModuleName)

WRITE(umMessage,'(A,L1)') ' l_glomap_clim_arg_act = ', l_glomap_clim_arg_act
CALL umPrint(umMessage, src=ModuleName)

WRITE(umMessage,'(A,L1)') ' l_glomap_clim_aie1 = ', l_glomap_clim_aie1
CALL umPrint(umMessage, src=ModuleName)

WRITE(umMessage,'(A,L1)') ' l_glomap_clim_aie2 = ', l_glomap_clim_aie2
CALL umPrint(umMessage, src=ModuleName)

WRITE(umMessage,'(A,L1)') ' l_glomap_clim_radaer = ', l_glomap_clim_radaer
CALL umPrint(umMessage, src=ModuleName)

WRITE(umMessage,'(A,L1)') ' l_glomap_clim_radaer_sustrat = ',                 &
                            l_glomap_clim_radaer_sustrat
CALL umPrint(umMessage, src=ModuleName)

WRITE(umMessage,'(A,I6)') ' i_glomap_clim_setup = ', i_glomap_clim_setup
CALL umPrint(umMessage, src=ModuleName)

WRITE(umMessage,"(A,A)")  ' gclmaclw = ', TRIM(gclmaclw)
CALL umprint(umMessage,src=ModuleName)

WRITE(umMessage,"(A,A)")  ' gclmacsw = ', TRIM(gclmacsw)
CALL umprint(umMessage,src=ModuleName)

WRITE(umMessage,"(A,A)")  ' gclmanlw = ', TRIM(gclmanlw)
CALL umprint(umMessage,src=ModuleName)

WRITE(umMessage,"(A,A)")  ' gclmansw = ', TRIM(gclmansw)
CALL umprint(umMessage,src=ModuleName)

WRITE(umMessage,"(A,A)")  ' gclmcrlw = ', TRIM(gclmcrlw)
CALL umprint(umMessage,src=ModuleName)

WRITE(umMessage,"(A,A)")  ' gclmcrsw = ', TRIM(gclmcrsw)
CALL umprint(umMessage,src=ModuleName)

WRITE(umMessage,"(A,A)")  ' gclmprec = ', TRIM(gclmprec)
CALL umprint(umMessage,src=ModuleName)

WRITE(umMessage,'(A)')    '- - - - - - end of namelist - - - - - -'
CALL umPrint(umMessage, src = ModuleName)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE print_nlist_run_glomap_clim

!---------------------------------------------------------------------------
! Description:
!   Subroutine to apply logic checks based on the
!   options selected in the run_glomap_aeroclim namelist.

SUBROUTINE check_glomap_clim_options()

USE chk_opts_mod, ONLY: &
    chk_var,            &
    def_src

IMPLICIT NONE

! Local variables
CHARACTER(LEN=*), PARAMETER       :: RoutineName='CHECK_GLOMAP_CLIM_OPTIONS'
CHARACTER(LEN=errormessagelength) :: cmessage
REAL(KIND=jprb)                   :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

def_src = RoutineName

IF (l_glomap_mode_clim) THEN
  ! check if valid value of i_glomap_clim_setup
  ! currently i_glomap_clim_setup must equal 2
  CALL chk_var ( i_glomap_clim_setup, 'i_glomap_clim_setup',                  &
                                      [i_gc_sussocbc_5mode],                  &
                 cmessage='Invalid value of i_glomap_clim_setup, ' //         &
                          'this must == 2 )' )
END IF

def_src = ''

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE check_glomap_clim_options

!---------------------------------------------------------------------------
! Description:
!   Subroutine to apply logic checks based on the
!   options selected in the run_glomap_aeroclim namelist.

SUBROUTINE check_run_glomap_clim()

USE ereport_mod,      ONLY: &
    ereport

USE mphys_inputs_mod, ONLY: &
    l_mcr_arcl,             &
    l_autoconv_murk

USE ukca_option_mod,  ONLY: &
    l_ukca,                 &
    l_ukca_radaer,          &
    l_ukca_arg_act,         &
    l_ukca_aie1,            &
    l_ukca_aie2

IMPLICIT NONE

! Local variables
CHARACTER(LEN=*), PARAMETER       :: RoutineName='CHECK_RUN_GLOMAP_CLIM'
CHARACTER(LEN=errormessagelength) :: cmessage
INTEGER                           :: errcode
REAL(KIND=jprb)                   :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Use of ukca and glomap_mode climatology
! This is not currently tested - it may be desirable to use both in future
IF (l_ukca .AND. l_glomap_mode_clim) THEN
  cmessage = 'Cannot set both l_ukca & l_glomap_mode_clim to .true.'
  errcode  = 54
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Direct aerosol effects
IF (l_ukca_radaer .AND. l_glomap_clim_radaer) THEN
  cmessage = 'Cannot supply RADAER with aerosol fields from both online ' //  &
             'GLOMAP-mode (UKCA) and GLOMAP-mode climatology aerosols'
  errcode  = 54
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Indirect aerosol effects
IF ( (l_glomap_clim_aie1 .OR. l_glomap_clim_aie2) .AND. .NOT.                 &
     l_glomap_mode_clim ) THEN
  cmessage = 'Cannot use AIE without GLOMAP-mode aerosols'
  errcode  = 54
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Activate
IF (l_ukca_arg_act .AND. l_glomap_clim_arg_act) THEN
  cmessage = 'Cannot feed both ukca GLOMAP-mode and climatology aerosols'
  errcode  = 54
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Aerosol indirect effect
IF ( ( l_glomap_clim_aie1 .OR. l_glomap_clim_aie2 )                           &
       .AND. .NOT. l_glomap_clim_arg_act ) THEN
  cmessage = 'Currently only Abdul-Razzak & Ghan method calculates cdnc. ' // &
             'Cannot set l_glomap_clim_aie1 or l_glomap_clim_aie2 '        // &
             ' without l_glomap_clim_arg_act.'
  errcode=54
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Aerosol indirect effect
IF (l_ukca_aie1 .AND. l_glomap_clim_aie1) THEN
  cmessage='Cannot set both l_ukca_aie1 and l_glomap_clim_aie1 to .true.'
  errcode=54
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Aerosol indirect effect
IF (l_ukca_aie2 .AND. l_glomap_clim_aie2) THEN
  cmessage='Cannot set both l_ukca_aie2 and l_glomap_clim_aie2 to .true.'
  errcode=54
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Aerosol indirect effect
IF (l_autoconv_murk .AND. l_glomap_clim_aie2) THEN
  cmessage='Cannot set both l_autoconv_murk and l_glomap_clim_aie2 to .true.'
  errcode=54
  CALL ereport(RoutineName,errcode,cmessage)
END IF

! Aerosol indirect effect
IF (l_mcr_arcl .AND. l_glomap_clim_aie2) THEN
  cmessage='Cannot set both l_mcr_arcl and l_glomap_clim_aie2 to .true.'
  errcode=54
  CALL ereport(RoutineName,errcode,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE check_run_glomap_clim

!---------------------------------------------------------------------------

SUBROUTINE read_nml_run_glomap_clim(unit_in)

USE check_iostat_mod, ONLY: &
    check_iostat

USE setup_namelist,   ONLY: &
    setup_nml_type

USE um_parcore,       ONLY: &
    mype

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in

! Local variables
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER(LEN=*), PARAMETER       :: RoutineName='READ_NML_RUN_GLOMAP_CLIM'
REAL(KIND=jprb)                   :: zhook_handle

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int       = 1
INTEGER, PARAMETER :: n_log       = 6
INTEGER, PARAMETER :: n_chars     = (7 * filenamelength)

TYPE my_namelist
  SEQUENCE
  INTEGER :: i_glomap_clim_setup
  LOGICAL :: l_glomap_mode_clim
  LOGICAL :: l_glomap_clim_arg_act
  LOGICAL :: l_glomap_clim_aie1
  LOGICAL :: l_glomap_clim_aie2
  LOGICAL :: l_glomap_clim_radaer
  LOGICAL :: l_glomap_clim_radaer_sustrat
  CHARACTER (LEN=filenamelength) :: gclmaclw 
  CHARACTER (LEN=filenamelength) :: gclmacsw 
  CHARACTER (LEN=filenamelength) :: gclmanlw
  CHARACTER (LEN=filenamelength) :: gclmansw
  CHARACTER (LEN=filenamelength) :: gclmcrlw
  CHARACTER (LEN=filenamelength) :: gclmcrsw
  CHARACTER (LEN=filenamelength) :: gclmprec
END TYPE my_namelist
TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,                &
                    n_log_in=n_log, n_chars_in=n_chars)

IF (mype == 0) THEN
  READ(UNIT=unit_in, NML=run_glomap_aeroclim, IOSTAT=ErrorStatus,             &
       IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_GLOMAP_AEROCLIM", iomessage)
  
  my_nml % i_glomap_clim_setup          = i_glomap_clim_setup
  my_nml % l_glomap_mode_clim           = l_glomap_mode_clim
  my_nml % l_glomap_clim_arg_act        = l_glomap_clim_arg_act
  my_nml % l_glomap_clim_aie1           = l_glomap_clim_aie1
  my_nml % l_glomap_clim_aie2           = l_glomap_clim_aie2
  my_nml % l_glomap_clim_radaer         = l_glomap_clim_radaer
  my_nml % l_glomap_clim_radaer_sustrat = l_glomap_clim_radaer_sustrat
  my_nml % gclmaclw                     = gclmaclw
  my_nml % gclmacsw                     = gclmacsw
  my_nml % gclmanlw                     = gclmanlw
  my_nml % gclmansw                     = gclmansw
  my_nml % gclmcrlw                     = gclmcrlw
  my_nml % gclmcrsw                     = gclmcrsw
  my_nml % gclmprec                     = gclmprec
END IF

CALL mpl_bcast(my_nml, 1, mpl_nml_type, 0, my_comm, icode)

IF (mype /= 0) THEN
  i_glomap_clim_setup          = my_nml % i_glomap_clim_setup
  l_glomap_mode_clim           = my_nml % l_glomap_mode_clim
  l_glomap_clim_arg_act        = my_nml % l_glomap_clim_arg_act
  l_glomap_clim_aie1           = my_nml % l_glomap_clim_aie1
  l_glomap_clim_aie2           = my_nml % l_glomap_clim_aie2
  l_glomap_clim_radaer         = my_nml % l_glomap_clim_radaer
  l_glomap_clim_radaer_sustrat = my_nml % l_glomap_clim_radaer_sustrat
  gclmaclw                     = my_nml % gclmaclw
  gclmacsw                     = my_nml % gclmacsw
  gclmanlw                     = my_nml % gclmanlw
  gclmansw                     = my_nml % gclmansw
  gclmcrlw                     = my_nml % gclmcrlw
  gclmcrsw                     = my_nml % gclmcrsw
  gclmprec                     = my_nml % gclmprec
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE read_nml_run_glomap_clim

END MODULE glomap_clim_option_mod
