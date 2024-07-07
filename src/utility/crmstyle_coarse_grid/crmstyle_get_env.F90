! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  Get environmental filenames for program

MODULE crmstyle_get_env_mod

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_GET_ENV_MOD'

CONTAINS

! ------------------------------------------------------------------------------
! Description:
! Gets various environmental filenames
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
! ------------------------------------------------------------------------------

SUBROUTINE crmstyle_get_env(icall)


USE crmstyle_cntl_mod, ONLY:                                            &
  num_ff, l_bcu_mask, l_class_col, l_no_orog, l_all_sea

USE crmstyle_grid_info_mod, ONLY:                                       &
  nprocs, nproc_x, nproc_y, thread_level_set

USE crmstyle_filenames_mod, ONLY:                                        &
  ProgName, input_file, orogfile, landseafile, lev_nl_file, pp_file      &
 ,all_file, acc_file, acu_file, bcu_file, wg1_file, ppd_file             &
 ,nbd_file, nid_file, adu_file, acw_file, bcw_file, bcu_mask_file        &
 ,ff_env_name

USE ereport_mod, ONLY: ereport

USE get_env_var_mod, ONLY: get_env_var

USE umPrintMgr, ONLY: umPrint, umMessage    ! for writing output

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) ::  &
  icall                    ! call number - detemines what is done


!---------------------------------------------------------------------------
! local variables
!---------------------------------------------------------------------------
INTEGER ::         &
   length          & ! length of returned string
 , i

CHARACTER(LEN=8) :: c_nproc            ! to get nproc_x and nproc_y from
CHARACTER(LEN=errormessagelength) :: cmessage            ! error message
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'CRMSTYLE_GET_ENV'


! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
! Processor/node info plus control namelist input filename

IF (icall == 1) THEN

  CALL get_env_var('NPROCX',c_nproc, allow_missing=.TRUE., length=length)
  IF (length < 0) THEN
    cmessage = 'Warning: Environment variable NPROCX has '              &
          //'not been set. Setting nproc_x to 1.'
    CALL ereport(RoutineName, length, cmessage)
    nproc_x=1
  ELSE
    ! convert to integer
    READ(c_nproc,'(I4)') nproc_x
  END IF

  CALL get_env_var('NPROCY',c_nproc, allow_missing=.TRUE., length=length)
  IF (length < 0) THEN
    cmessage = 'Warning: Environment variable NPROCY has '              &
          //'not been set. Setting nproc_y to 1.'
    CALL ereport(RoutineName, length, cmessage)
    nproc_y=1
  ELSE
    ! convert to integer
    READ(c_nproc,'(I4)') nproc_y
  END IF

  ! total number of processors/nodes
  CALL get_env_var('NPROC',c_nproc, allow_missing=.TRUE., length=length)
  IF (length < 0) THEN
    cmessage = 'Warning: Environment variable NPROC has '              &
        //'not been set. Setting nprocs to 1.'
    CALL ereport(RoutineName, length, cmessage)
    nprocs=1
  ELSE
    ! convert to integer
    READ(c_nproc,'(I4)') nprocs
  END IF

  CALL get_env_var( "CRMSTYLE_INPUT", input_file)

END IF

!------------------------------------------------------------------------------
! Second call now know num_ff and l_bcu_mask
!------------------------------------------------------------------------------

IF (icall == 2) THEN

  IF (.NOT. l_no_orog) THEN
    CALL get_env_var( "OROG_FILE", orogfile )
  END IF

  IF (.NOT. l_all_sea) THEN
    CALL get_env_var( "MASK_FILE", landseafile )
  END IF

  CALL get_env_var( "VERT_NAMELIST", lev_nl_file )

  ! PP input filenames
  ! First file is compulsory:
  CALL get_env_var(ff_env_name(1), pp_file(1))
  WRITE(umMessage,'(A9,I2,A,A)') ' pp_file ',i,ff_env_name(1),pp_file(1)
  CALL umPrint(umMessage,src=RoutineName)

  ! Files 2-24 are optional:
  DO i = 2,num_ff

    CALL get_env_var(ff_env_name(i), pp_file(i), allow_missing=.TRUE., &
                     length=length)

    IF (length > 0) THEN
      WRITE(umMessage,'(A9,I2,A,A)') ' pp_file ',i,ff_env_name(i),pp_file(i)
      CALL umPrint(umMessage,src=RoutineName)
    END IF

  END DO

  ! Output filenames

  CALL get_env_var( "ALL_FILE", all_file )

  CALL get_env_var( "ACC_FILE", acc_file )

  CALL get_env_var( "ACU_FILE", acu_file )

  CALL get_env_var( "BCU_FILE", bcu_file )

  CALL get_env_var( "WG1_FILE", wg1_file )

  CALL get_env_var( "PPD_FILE", ppd_file )

  CALL get_env_var( "NBD_FILE", nbd_file )

  CALL get_env_var( "NID_FILE", nid_file )

  CALL get_env_var( "ADU_FILE", adu_file )

  CALL get_env_var( "ACW_FILE", acw_file )

  CALL get_env_var( "BCW_FILE", bcw_file )

  ! Only required if l_bcu_mask but not read namelist at this stage
  ! Also required for output of column classification full fields
  IF (l_bcu_mask .OR. l_class_col) THEN
    CALL get_env_var( "BCU_MASK", bcu_mask_file )
  END IF


END IF
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_get_env

END MODULE crmstyle_get_env_mod
