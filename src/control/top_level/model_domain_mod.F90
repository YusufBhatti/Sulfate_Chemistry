! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Purpose: Read in and hold information related to model domain. This
!          relates to the domain properties should not be confused
!          with grid properties, i.e. number of rows/columns etc.
!

MODULE model_domain_mod

USE missing_data_mod,       ONLY: rmdi, imdi
USE errormessagelength_mod, ONLY: errormessagelength
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!-----------------------------------------------------
! Parameters
!-----------------------------------------------------
! Possible values for namelist variable - model_type
! 0 => Small execs don't use this, so we have a placeholder value
! 1 => Global Model (i.e. cyclic EW, and NS overpole - a sphere)
! 2 => Limited Area Model (Classic style)
! 3 => Limited Area Model (Cyclic boundary conditions - EW only)
! 4 => Limited Area Model (Cyclic boundary conditions - EW and NS, (i.e. torus)
! 5 => Single Column Model (i.e. one horizontal point)

INTEGER, PARAMETER :: mt_smexe         = 0
INTEGER, PARAMETER :: mt_global        = 1
INTEGER, PARAMETER :: mt_lam           = 2
INTEGER, PARAMETER :: mt_cyclic_lam    = 3
INTEGER, PARAMETER :: mt_bi_cyclic_lam = 4
INTEGER, PARAMETER :: mt_single_column = 5

! DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1


! Item not in the model domain namelist
! It is read in and set up by the reconfiguration within the horizont namelist.
INTEGER :: output_grid_stagger = 6

INTEGER,PARAMETER :: FH_GridStagger_C         = 3
INTEGER,PARAMETER :: FH_GridStagger_Endgame   = 6

!-----------------------------------------------------
! Items set in namelist
!-----------------------------------------------------

INTEGER :: model_type = imdi
LOGICAL :: l_regular  = .FALSE. ! TRUE if NOT variable resolution
LOGICAL :: l_cartesian = .FALSE.   ! True = Cartesian co-ordinates

NAMELIST /model_domain/ model_type, l_regular, l_cartesian


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MODEL_DOMAIN_MOD'

CONTAINS

SUBROUTINE check_nml_model_domain()

! Description:
!   Subroutine to apply checks based on options read in model_domain namelist.

USE chk_opts_mod, ONLY: chk_var, def_src
USE ereport_mod,  ONLY: ereport

IMPLICIT NONE

CHARACTER (LEN=*), PARAMETER       :: RoutineName = 'CHECK_NML_MODEL_DOMAIN'
REAL(KIND=jprb)                    :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

! Check available model domain type is chosen
CALL chk_var( model_type, 'model_type'                             &
            , [mt_global, mt_lam, mt_cyclic_lam, mt_bi_cyclic_lam, &
               mt_single_column] )

def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_nml_model_domain

#if !defined(LFRIC)
SUBROUTINE read_nml_model_domain(unit_in)

USE um_parcore,       ONLY: mype
USE check_iostat_mod, ONLY: check_iostat
USE setup_namelist,   ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'READ_NML_MODEL_DOMAIN'
REAL(KIND=jprb)   :: zhook_handle

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 1
INTEGER, PARAMETER :: n_log = 2

TYPE my_namelist
  SEQUENCE
  INTEGER :: model_type
  LOGICAL :: l_regular
  LOGICAL :: l_cartesian
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=model_domain, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist MODEL_DOMAIN", iomessage)

  my_nml % model_type = model_type
  my_nml % l_regular  = l_regular
  my_nml % l_cartesian  = l_cartesian

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  model_type = my_nml % model_type
  l_regular  = my_nml % l_regular
  l_cartesian  = my_nml % l_cartesian

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE  read_nml_model_domain
#endif

SUBROUTINE print_nlist_model_domain()

USE umPrintMgr, ONLY: umPrint

IMPLICIT NONE

CHARACTER(LEN=50000) :: lineBuffer
CHARACTER(LEN=*), PARAMETER    :: RoutineName = 'PRINT_NLIST_MODEL_DOMAIN'
REAL(KIND=jprb)      :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist model_domain', src='model_domain_mod')

WRITE(lineBuffer,'(A,I2)') ' model_type = ', model_type
CALL umPrint(lineBuffer, src='model_domain_mod')
WRITE(lineBuffer,'(A,L1)') ' l_regular  = ', l_regular
CALL umPrint(lineBuffer, src='model_domain_mod')
WRITE(lineBuffer,'(A,L1)') ' l_cartesian  = ', l_cartesian
CALL umPrint(lineBuffer, src='model_domain_mod')
CALL umPrint('- - - - - - end of namelist - - - - - -', src='model_domain_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_model_domain

END MODULE model_domain_mod
