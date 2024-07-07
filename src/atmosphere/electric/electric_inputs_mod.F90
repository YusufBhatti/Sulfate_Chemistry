! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

! Purpose:
!   Module containing runtime options/data used by the electrification scheme

! Method:
!   Switches and associated data values used by the electrification scheme
!   are defined here and assigned default values. These may be overridden
!   by namelist input.

!   A description of what each switch or number refers to is provided
!   with the namelist

!   Any routine wishing to use these options may do so with the 'USE'
!   statement.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Electric

MODULE electric_inputs_mod

USE missing_data_mod,       ONLY: imdi, rmdi
USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
!===========================================================================
! Logical options for the electric namelist (default is .false.)
!===========================================================================
LOGICAL :: l_use_electric = .FALSE.

!===========================================================================
! Integer options for the electric namelist (default is imdi)
!===========================================================================
INTEGER :: electric_method = imdi

! Parameters for electric_method, used to define which scheme to use
! in flash_rate_mod but not used in namelist
INTEGER, PARAMETER :: em_gwp    = 1 
INTEGER, PARAMETER :: em_mccaul = 2

!===========================================================================
! Real options for the electric namelist (default is rmdi)
!===========================================================================

! Settings for Mccaul et al (2009) scheme
REAL :: k1 = rmdi

REAL :: k2 = rmdi

! Settings for Graupel water path based scheme

REAL :: g1 = rmdi

REAL :: g2 = rmdi


! Define the RUN_ELECTRIC namelist

NAMELIST/run_electric/                                                        &
l_use_electric, electric_method, k1, k2, g1, g2

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ELECTRIC_INPUTS_MOD'

CONTAINS

SUBROUTINE print_nlist_run_electric()

USE umPrintMgr, ONLY: umPrint

IMPLICIT NONE

! Local variables

CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_ELECTRIC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Start of subroutine printing

CALL umPrint('Contents of namelist run_electric',                        &
              src='electric_inputs_mod')

WRITE(lineBuffer,*)'l_use_electric = ',l_use_electric
CALL umPrint(lineBuffer,src='electric_inputs_mod')

WRITE(lineBuffer,*)'electric_method =', electric_method
CALL umPrint(lineBuffer,src='electric_inputs_mod')

WRITE(lineBuffer,*)'k1 =', k1
CALL umPrint(lineBuffer,src='electric_inputs_mod')

WRITE(lineBuffer,*)'k2 =', k2
CALL umPrint(lineBuffer,src='electric_inputs_mod')

WRITE(lineBuffer,*)'g1 =', g1
CALL umPrint(lineBuffer,src='electric_inputs_mod')

WRITE(lineBuffer,*)'g2 =', g2
CALL umPrint(lineBuffer,src='electric_inputs_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -',                  &
              src='electric_inputs_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_electric

#if !defined(LFRIC)
SUBROUTINE read_nml_run_electric(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_ELECTRIC'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 1
INTEGER, PARAMETER :: n_real = 4
INTEGER, PARAMETER :: n_log = 1

TYPE my_namelist
  SEQUENCE
  INTEGER :: electric_method
  REAL :: k1
  REAL :: k2
  REAL :: g1
  REAL :: g2
  LOGICAL :: l_use_electric
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                    n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=RUN_Electric, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_Electric", iomessage)

  my_nml % electric_method = electric_method
  my_nml % k1 = k1
  my_nml % k2 = k2
  my_nml % g1 = g1
  my_nml % g2 = g2
  my_nml % l_use_electric = l_use_electric

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  electric_method = my_nml % electric_method
  k1 = my_nml % k1
  k2 = my_nml % k2
  g1 = my_nml % g1
  g2 = my_nml % g2
  l_use_electric = my_nml % l_use_electric

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_electric
#endif

SUBROUTINE check_run_electric

USE chk_opts_mod,           ONLY: chk_var, def_src

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_RUN_ELECTRIC'

CHARACTER(LEN=errormessagelength) :: comments
CHARACTER(LEN=100) :: ChkStr

REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_use_electric) THEN
  ! Only bother checking the namelist if the electric scheme is in use.
  ! This IF test should perhaps be in readlsta, but we'll leave it here
  ! for now in early testing.

  ! 1. Check electric_method is sensible. This must fail if it is not
  ! a sensible number as everything else in the scheme is dependent
  ! on it

  comments = 'Electric method should be selected for scheme to run'

  WRITE(ChkStr,'(2(A,I1),A)') '[',em_gwp,',',em_mccaul,']'

  CALL chk_var( electric_method, 'electric_method', ChkStr, cmessage=comments)
                
  ! If we get this far, electric_method is OK, so can check other variables
  ! can now be checked based on the value of electric_method

  SELECT CASE (electric_method)

  CASE (em_gwp)

    ! Need to check g1 and g2
    ! These should mirror exactly what is in the Rose metadata
    CALL chk_var( g1, 'g1', '[0.0:100.0]' )
    CALL chk_var( g2, 'g2', '[-1000.0:1000.0]' )

  CASE (em_mccaul)

    ! Need to check k1 and k2
    ! These should mirror exactly what is in the Rose metadata
    CALL chk_var( k1, 'k1', '[0.0:50.0]' )
    CALL chk_var( k2, 'k2', '[0.0:50.0]' )

  END SELECT

END IF ! l_use_electric

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)


END SUBROUTINE check_run_electric

END MODULE electric_inputs_mod
