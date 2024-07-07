! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Vertical grid information

MODULE vertnamelist_mod

USE Atmos_Max_Sizes, ONLY:   model_levels_max

USE missing_data_mod, ONLY:  rmdi, imdi

USE errormessagelength_mod, ONLY:                                        &
                             errormessagelength

USE Ereport_Mod, ONLY:       Ereport

USE nlsizes_namelist_mod, ONLY:                                          &
                             model_levels

IMPLICIT NONE

PRIVATE

PUBLIC ::                                                                &
                             vertlevs,                                   &
                             first_constant_r_rho_level,                 &
                             z_top_of_model,                             &
                             eta_theta,                                  &
                             eta_rho,                                    &
                             read_nml_vertlevs,                          &
                             check_nml_vertlevs,                         &
                             print_nlist_vertlevs

! Description:
!   Contains the vertical grid namelist for variable resolution, shared
!   between LBC generation, CreateBC, the SCM and the reconfiguration
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Grids
!
! Code Description:
!   Language: Fortran 95
!   This code is written to UMDP3 v10.3 programming standards.
!

INTEGER :: first_constant_r_rho_level  ! Lowest level where rho level has
                                       ! constant radius
                                       ! (i.e. not terrain-following)

REAL :: z_top_of_model                 ! Top of model (metres)
REAL :: eta_theta (model_levels_max+1) ! Theta levels (as fraction of whole)
REAL :: eta_rho (model_levels_max)     ! Rho levels (as fraction of whole)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='VERTNAMELIST_MOD'

NAMELIST /vertlevs/                                                            &
  first_constant_r_rho_level, z_top_of_model, eta_theta, eta_rho

CONTAINS

SUBROUTINE read_nml_vertlevs(unit_in)

USE check_iostat_mod, ONLY:           check_iostat
USE setup_namelist, ONLY:             setup_nml_type
USE um_parcore, ONLY:                 mype

! DrHook
USE yomhook,   ONLY:                  lhook,           &
                                      dr_hook

USE parkind1,  ONLY:                  jprb,            &
                                      jpim

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
    
CHARACTER(LEN=errormessagelength) :: iomessage

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 1
INTEGER, PARAMETER :: n_real = 2 + 2*model_levels_max

TYPE my_namelist
  SEQUENCE
  INTEGER :: first_constant_r_rho_level
  REAL :: z_top_of_model
  REAL :: eta_theta (model_levels_max+1)
  REAL :: eta_rho (model_levels_max)
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

CHARACTER (LEN=errormessagelength)           :: Cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'READ_NML_VERTLEVS'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set VERTLEVS defaults
eta_theta(:)                  = rmdi
eta_rho(:)                    = rmdi
z_top_of_model                = rmdi
first_constant_r_rho_level    = imdi

! Quick error check to make sure parameter model_levels_max isn't too small
IF ( model_levels_max < model_levels ) THEN
  ErrorStatus = 10
  Cmessage = 'Internal parameter model_levels_max is too small!'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int, &
                    n_real_in=n_real)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=vertlevs, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist vertlevs", iomessage)

  my_nml % first_constant_r_rho_level = first_constant_r_rho_level
  my_nml % z_top_of_model             = z_top_of_model
  my_nml % eta_theta                  = eta_theta
  my_nml % eta_rho                    = eta_rho

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  first_constant_r_rho_level = my_nml % first_constant_r_rho_level
  z_top_of_model             = my_nml % z_top_of_model
  eta_theta                  = my_nml % eta_theta
  eta_rho                    = my_nml % eta_rho

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE read_nml_vertlevs

SUBROUTINE print_nlist_vertlevs()

USE umPrintMgr, ONLY:        umPrint

! DrHook
USE yomhook,   ONLY:         lhook,                                      &
                             dr_hook

USE parkind1,  ONLY:         jprb,                                       &
                             jpim

IMPLICIT NONE
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'PRINT_NLIST_VERTLEVS'

CHARACTER(LEN=50000) :: lineBuffer
INTEGER              :: i

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist vertlevs',         &
                        src='vertnamelist_mod')

WRITE(linebuffer,'(A,F16.4)' ) 'z_top_of_model ',z_top_of_model
CALL umPrint(linebuffer,src='vertnamelist_mod')
WRITE(linebuffer,'(A,I0)' ) 'first_constant_r_rho_level ',    &
                             first_constant_r_rho_level
CALL umPrint(linebuffer,src='vertnamelist_mod')
WRITE(linebuffer,'(A)' ) 'Eta_Theta'
CALL umPrint(linebuffer,src='vertnamelist_mod')

DO i=1,model_levels + 1
  WRITE(linebuffer,'(F10.7)' ) eta_theta(i)
  CALL umPrint(linebuffer,src='vertnamelist_mod')
END DO
WRITE(linebuffer,'(A)' ) 'Eta_Rho'
CALL umPrint(linebuffer,src='vertnamelist_mod')
DO i=1,model_levels
  WRITE(linebuffer,'(F10.7)' ) eta_rho(i)
  CALL umPrint(linebuffer,src='vertnamelist_mod')
END DO

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='vertnamelist_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE print_nlist_vertlevs

SUBROUTINE check_nml_vertlevs()
! Description:
!   Subroutine to check variables based in recon_vertical namelist.

USE chk_opts_mod, ONLY: chk_var, def_src

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'CHECK_NML_VERTLEVS'

INTEGER                      :: levels_theta  ! number of levels of
INTEGER                      :: levels_rho    ! values read in.
INTEGER                      :: i             ! looper
INTEGER                      :: ErrorStatus
REAL                         :: last_value    ! used in checking code
CHARACTER (LEN=errormessagelength)           :: Cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

! Check that the number of eta values correspond to the number
! of model levels
levels_theta = 0
levels_rho   = 0
DO i = 1, model_levels_max
  IF (eta_theta(i) /= rmdi) THEN
    levels_theta = levels_theta + 1
  END IF

  IF (eta_rho(i) /= rmdi) THEN
    levels_rho = levels_rho + 1
  END IF
END DO

CALL chk_var( levels_theta, 'no. of theta levels', [model_levels + 1] )

CALL chk_var( levels_rho, 'no. of rho levels', [model_levels] )

! Check that eta_rho levels and eta_theta levels are monotone
! ascending.
last_value = rmdi
DO i = 1, model_levels + 1
  IF ( eta_theta(i) <= last_value ) THEN
    Cmessage = 'Eta_Theta values are not monotone ascending'
    ErrorStatus = 70
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  last_value = eta_theta(i)
END DO

last_value = rmdi
DO i = 1, model_levels
  IF ( eta_rho(i) <= last_value ) THEN
    Cmessage = 'Eta_Rho values are not monotone ascending'
    ErrorStatus = 80
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  last_value = eta_rho(i)
END DO

! z_top_of_model should be > twice the max height of orography - which we've
! chosen not to code up at this time - as we're not sure we've read orography
! in by this point in execution....
CALL chk_var( z_top_of_model, 'z_top_of_model', '[>0.0]' )

! like z_top_of_model, first_constant_r_rho_level should be 
!   > twice the max height of orography
CALL chk_var( first_constant_r_rho_level, 'first_constant_r_rho_level','[>0]')

def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE check_nml_vertlevs

END MODULE vertnamelist_mod
