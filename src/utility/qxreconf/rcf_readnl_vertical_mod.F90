! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  reads the vertical levels namelists

MODULE rcf_readnl_vertical_Mod

USE Rcf_V_Int_Ctl_Mod, ONLY: v_int_order

USE umPrintMgr, ONLY:        umPrint

IMPLICIT NONE

PRIVATE

PUBLIC ::                                                                &
        Rcf_Readnl_Vertical,    &
        print_nlist_vertical,   &
        check_nml_vertical

NAMELIST /recon_vertical/ v_int_order

!  Subroutine Rcf_Readnl_vertical - reads the vertical levels namelists
!
! Description:
! Module to read in the VERTICAL levels namelist
! Note that it *must* be called after readnl_recona as
! some sizing is required!
!
! Method:
!  Contains subroutines to read, check and print the VERTICAL namelist.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_READNL_VERTICAL_MOD'

CONTAINS

SUBROUTINE Rcf_Readnl_Vertical( unit_in )

USE check_iostat_mod, ONLY:   check_iostat

USE setup_namelist, ONLY:     setup_nml_type

USE um_parcore, ONLY:         mype

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 1
INTEGER, PARAMETER :: n_int = 1

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_READNL_VERTICAL'

TYPE my_namelist
  SEQUENCE
  INTEGER :: v_int_order                 
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int )

IF (mype == 0) THEN
  READ(UNIT=unit_in, NML=recon_vertical, IOSTAT=ErrorStatus)
  my_nml % v_int_order = v_int_order
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN
  v_int_order = my_nml % v_int_order
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE Rcf_Readnl_Vertical


SUBROUTINE print_nlist_vertical()

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

CHARACTER(LEN=50000) :: lineBuffer
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'PRINT_NLIST_VERTICAL'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist recon_vertical',         &
                        src='rcf_readnl_vertical_mod')

WRITE(lineBuffer,'(A,I0)')' v_int_order = ',v_int_order
CALL umPrint(lineBuffer,src='rcf_readnl_vertical_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='rcf_readnl_vertical_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE print_nlist_vertical


SUBROUTINE check_nml_vertical()
! Description:
!   Subroutine to check variables based in recon_vertical namelist.

USE chk_opts_mod, ONLY: chk_var, def_src

USE interpor_mod, ONLY:       &
  interp_order_linear,        &
  interp_order_linear_noex,   &
  interp_order_cubic,         &
  interp_order_quintic

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'CHECK_NML_VERTICAL'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

CALL chk_var( v_int_order, 'v_int_order', [              &
              interp_order_linear,                       &
              interp_order_linear_noex,                  &
              interp_order_cubic,                        &
              interp_order_quintic] )

def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_nml_vertical

END MODULE Rcf_Readnl_Vertical_Mod
