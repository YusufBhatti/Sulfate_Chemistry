! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Read the HORIZONT namelist

MODULE rcf_readnl_horizont_mod

!  Subroutine Rcf_Readnl_Horizont - Read the HORIZONT namelist
!
! Description:
!   Reads the HORIZONT namelist controlling horizontal interpolation
!   and domains.
!
! Method:
!   Data  read and Output_Grid set accordingly.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE missing_data_mod, ONLY: imdi, rmdi

USE model_domain_mod, ONLY: output_grid_stagger

USE Rcf_Interp_Weights_Mod, ONLY: &
    h_int_method,              &
    bilinear,                  &
    area_weighted,             &
    nearest_neighbour,         &
    smcp_int_nearest_neighbour,&
    l_limit_rotations

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

PRIVATE :: lhook, dr_hook, jpim, jprb

! Data for namelist reads - mostly for LAM
INTEGER, PARAMETER     :: iproj = imdi    ! 'Used' to populate lookups and FH
INTEGER, SAVE          :: orog_blend_width = imdi     !} For orography blending
REAL, ALLOCATABLE      :: blend_weights(:)            !} zone

! storage for namelist information
INTEGER, PARAMETER               :: orog_blend_max = 40
REAL                             :: orog_blend_weights(orog_blend_max)

NAMELIST /horizont/ h_int_method,                                      &
                    orog_blend_weights,                                &
                    smcp_int_nearest_neighbour, l_limit_rotations,     &
                    output_grid_stagger

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_READNL_HORIZONT_MOD'

CONTAINS

SUBROUTINE rcf_readnl_horizont( nft )

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)              :: nft             ! File Unit

! Local vars/params
CHARACTER (LEN=*), PARAMETER     :: RoutineName = 'RCF_READNL_HORIZONT'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set defaults
orog_blend_weights(:)  = rmdi

! Read horizont Namelist
CALL read_nml_horizont(nft)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_readnl_horizont

SUBROUTINE print_nlist_horizont()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_HORIZONT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist horizont', &
    src='rcf_readnl_horizont_mod')

WRITE(lineBuffer,*)'h_int_method = ',h_int_method
CALL umPrint(lineBuffer,src='rcf_readnl_horizont_mod')
WRITE(lineBuffer,*)'orog_blend_weights = ',orog_blend_weights
CALL umPrint(lineBuffer,src='rcf_readnl_horizont_mod')
WRITE(lineBuffer,*)'smcp_int_nearest_neighbour = ', &
    smcp_int_nearest_neighbour
CALL umPrint(lineBuffer,src='rcf_readnl_horizont_mod')
WRITE(lineBuffer,*)'l_limit_rotations = ',l_limit_rotations
CALL umPrint(lineBuffer,src='rcf_readnl_horizont_mod')
WRITE(lineBuffer,'(A,I10)')'output_grid_stagger == ',output_grid_stagger
CALL umPrint(lineBuffer,src='rcf_readnl_horizont_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='rcf_readnl_horizont_mod')
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE print_nlist_horizont

SUBROUTINE read_nml_horizont(unit_in)

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

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 2
INTEGER, PARAMETER :: n_real = orog_blend_max
INTEGER, PARAMETER :: n_log = 2

TYPE my_namelist
  SEQUENCE
  INTEGER :: h_int_method
  INTEGER :: output_grid_stagger
  REAL    :: orog_blend_weights(orog_blend_max)
  LOGICAL :: smcp_int_nearest_neighbour
  LOGICAL :: l_limit_rotations
END TYPE my_namelist

TYPE (my_namelist) :: my_nml
CHARACTER (LEN=*), PARAMETER  :: RoutineName='READ_NML_HORIZONT'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                    n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=horizont, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist horizont", iomessage)

  my_nml % h_int_method               = h_int_method
  my_nml % output_grid_stagger        = output_grid_stagger
  my_nml % orog_blend_weights         = orog_blend_weights
  my_nml % smcp_int_nearest_neighbour = smcp_int_nearest_neighbour
  my_nml % l_limit_rotations          = l_limit_rotations

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  h_int_method                = my_nml % h_int_method
  output_grid_stagger         = my_nml % output_grid_stagger
  orog_blend_weights          = my_nml % orog_blend_weights
  smcp_int_nearest_neighbour  = my_nml % smcp_int_nearest_neighbour
  l_limit_rotations           = my_nml % l_limit_rotations

END IF

CALL mpl_type_free(mpl_nml_type,icode)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_nml_horizont

SUBROUTINE check_nml_horizont()

USE chk_opts_mod, ONLY:                chk_var, def_src

USE Rcf_HeadAddress_Mod, ONLY:         FH_GridStagger_C,           &
                                       FH_GridStagger_Endgame

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE

CHARACTER (LEN=*), PARAMETER       :: RoutineName = 'CHECK_NML_HORIZONT'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER      :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER      :: zhook_out = 1
REAL(KIND=jprb)                    :: zhook_handle

INTEGER                            :: icount
CHARACTER (LEN=errormessagelength) :: Cmessage

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)
def_src = RoutineName

CALL chk_var( h_int_method, 'h_int_method', [ bilinear, area_weighted,     &
                                              nearest_neighbour ]  )
! find the apparent size of the blending width set in the inputs.
orog_blend_width=0
DO icount = 1, orog_blend_max
  IF (orog_blend_weights(icount) /= rmdi ) THEN
    orog_blend_width=orog_blend_width+1
  END IF
END DO
IF (orog_blend_width > 0 .AND. model_type /= mt_global ) THEN
  DO icount = 1,orog_blend_width
    WRITE (Cmessage, '(A,I0,A)') "Element ",icount," of orog_blend_weights" // &
                                 " appears to be out of range in HORIZONT"
    CALL chk_var( orog_blend_weights(icount), 'orog_blend_weights',       &
                 '[0.0:1.0]', cmessage=cmessage )
  END DO
END IF

CALL chk_var( output_grid_stagger, 'output_grid_stagger',                 &
             [FH_GridStagger_C, FH_GridStagger_Endgame] )


def_src = ''
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE check_nml_horizont

END MODULE rcf_readnl_horizont_mod
