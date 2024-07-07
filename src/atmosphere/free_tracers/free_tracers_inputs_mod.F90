! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************

!  Data module for free tracers related switches and settings
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Free Tracers
MODULE free_tracers_inputs_mod

USE missing_data_mod, ONLY: rmdi, imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

!-----------------------------------------------------
! Items set in module
!-----------------------------------------------------

INTEGER,PARAMETER:: a_tracer_first = 1
! First atmospheric tracer (STASH No)

INTEGER,PARAMETER:: a_tracer_last = 150
! Last atmospheric tracer  (STASH No)

INTEGER, PARAMETER :: a_max_trvars = (a_tracer_last - a_tracer_first) + 1
! Max number of tracers allowed

LOGICAL :: tracer_a (0:a_max_trvars)

INTEGER :: a_tr_index(a_max_trvars)
! Index to relative position.
! a_tr_index(n) gives position in jtracer for tracer number N.
! Set in set_atm_pointers.
! a_tr_index(n) is the position, in the list of tracers
! actually present in D1, that tracer number N (in the list
! of all tracers selectable from the user interface) occupies,
! if it is present.
! If tracer number N is absent then A_TR_INDEX(N) is -1.

INTEGER :: a_tr_stashitem(a_max_trvars)
! a_tr_stashitem is set up in set_atm_pointers

INTEGER :: a_tr_lbc_stashitem(a_max_trvars)
! a_tr_lbc_stashitem is set up in inbounda and is only
! referenced if LBC code is active.

INTEGER :: a_tr_active_lbc_index(a_max_trvars)

!-----------------------------------------------------
! Items set in namelist
!-----------------------------------------------------

LOGICAL  ::  l_free_tracer = .FALSE.
! Include tracers in the atmosphere

LOGICAL  ::  l_free_tracer_lbc = .FALSE.
! Turn on free tracer LBCs

INTEGER  ::  i_free_tracer(a_max_trvars) = 0
! Specifies which STASH tracer items to be included,
! these will have STASHmaster entries from section 33
!  0 - Do not include
!  1 - Include from dump

INTEGER  ::  i_free_tracer_lbc(a_max_trvars) = 0
! Specifies which tracers have lateral boundary condition
! data in the LBC input file.
! 0 - No LBC data
! 1 - LBC data present

LOGICAL :: l_bl_tracer_mix = .FALSE.
! Boundary layer tracer mixing

! ---- PV tracer variables
LOGICAL :: l_pv_tracer = .FALSE.

LOGICAL :: l_pv_dyn = .FALSE.

LOGICAL :: l_pv_split_rad = .FALSE.

LOGICAL :: l_pv_adv_only = .FALSE.

LOGICAL :: l_calc_pv_full = .FALSE.

INTEGER :: num_dPV = imdi

! ---- Diabatic tracer variables

LOGICAL :: l_diab_tracer= .FALSE.

LOGICAL :: l_diab_tr_rad= .FALSE.

LOGICAL :: l_diab_tr_bl= .FALSE.

INTEGER :: num_dtheta = imdi

NAMELIST /run_free_tracers/ l_free_tracer, l_free_tracer_lbc,     &
     i_free_tracer, i_free_tracer_lbc, l_bl_tracer_mix,           &
     l_pv_tracer, l_pv_dyn, l_pv_split_rad,l_pv_adv_only,         &
     l_calc_pv_full,                                              &
     l_diab_tracer,l_diab_tr_rad,l_diab_tr_bl


!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FREE_TRACERS_INPUTS_MOD'

CONTAINS

SUBROUTINE print_nlist_run_free_tracers()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_FREE_TRACERS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_free_tracers', &
     src='free_tracers_inputs_mod')

WRITE(lineBuffer,*) ' l_free_tracer = ',l_free_tracer
CALL umPrint(lineBuffer,src='free_tracers_inputs_mod')
WRITE(lineBuffer,*) ' l_free_tracer_lbc  = ',l_free_tracer_lbc
CALL umPrint(lineBuffer,src='free_tracers_inputs_mod')
WRITE(lineBuffer,*) ' i_free_tracer = ',i_free_tracer
CALL umPrint(lineBuffer,src='free_tracers_inputs_mod')
WRITE(lineBuffer,*) ' i_free_tracer_lbc = ',i_free_tracer_lbc
CALL umPrint(lineBuffer,src='free_tracers_inputs_mod')
WRITE(lineBuffer,*) ' l_bl_tracer_mix = ',l_bl_tracer_mix
CALL umPrint(lineBuffer,src='free_tracers_inputs_mod')
WRITE(lineBuffer,'(A,L1)') ' l_pv_tracer = ',l_pv_tracer
CALL umPrint(lineBuffer,src='free_tracers_inputs_mod')
WRITE(lineBuffer,'(A,L1)') ' l_pv_dyn = ',l_pv_dyn
CALL umPrint(lineBuffer,src='free_tracers_inputs_mod')
WRITE(lineBuffer,'(A,L1)') ' l_pv_split_rad = ',l_pv_split_rad
CALL umPrint(lineBuffer,src='free_tracers_inputs_mod')
WRITE(lineBuffer,'(A,L1)') ' l_pv_adv_only = ',l_pv_adv_only
CALL umPrint(lineBuffer,src='free_tracers_inputs_mod')
WRITE(lineBuffer,'(A,L1)') ' l_calc_pv_full = ',l_calc_pv_full
CALL umPrint(lineBuffer,src='free_tracers_inputs_mod')
WRITE(lineBuffer,'(A,L1)') ' l_diab_tracer = ',l_diab_tracer
CALL umPrint(lineBuffer,src='free_tracers_inputs_mod')
WRITE(lineBuffer,'(A,L1)') ' l_diab_tr_rad = ',l_diab_tr_rad
CALL umPrint(lineBuffer,src='free_tracers_inputs_mod')
WRITE(lineBuffer,'(A,L1)') ' l_diab_tr_bl = ',l_diab_tr_bl
CALL umPrint(lineBuffer,src='free_tracers_inputs_mod')

WRITE(lineBuffer,'(A)') '- - - - - - end of namelist - - - - - -'
CALL umPrint(lineBuffer,src='free_tracers_inputs_mod')


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_free_tracers

#if !defined(LFRIC)
SUBROUTINE read_nml_run_free_tracers(unit_in)

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_FREE_TRACERS'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 2 * a_max_trvars
INTEGER, PARAMETER :: n_log = 11


TYPE my_namelist
  SEQUENCE
  INTEGER :: i_free_tracer(a_max_trvars)
  INTEGER :: i_free_tracer_lbc(a_max_trvars)
  LOGICAL :: l_free_tracer
  LOGICAL :: l_free_tracer_lbc
  LOGICAL :: l_bl_tracer_mix
  LOGICAL :: l_pv_tracer
  LOGICAL :: l_pv_dyn  
  LOGICAL :: l_pv_split_rad
  LOGICAL :: l_pv_adv_only
  LOGICAL :: l_calc_pv_full
  LOGICAL :: l_diab_tracer
  LOGICAL :: l_diab_tr_rad
  LOGICAL :: l_diab_tr_bl
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                    n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=run_free_tracers, IOSTAT=ErrorStatus,     &
       IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist run_free_tracers", iomessage)

  my_nml % i_free_tracer     = i_free_tracer
  my_nml % i_free_tracer_lbc = i_free_tracer_lbc
  my_nml % l_free_tracer     =  l_free_tracer
  my_nml % l_free_tracer_lbc =  l_free_tracer_lbc
  my_nml % l_bl_tracer_mix   = l_bl_tracer_mix
  my_nml % l_pv_tracer       = l_pv_tracer
  my_nml % l_pv_dyn          = l_pv_dyn
  my_nml % l_pv_split_rad    = l_pv_split_rad
  my_nml % l_pv_adv_only     = l_pv_adv_only
  my_nml % l_calc_pv_full    = l_calc_pv_full
  my_nml % l_diab_tracer     = l_diab_tracer
  my_nml % l_diab_tr_rad     = l_diab_tr_rad
  my_nml % l_diab_tr_bl      = l_diab_tr_bl
  
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  i_free_tracer     = my_nml % i_free_tracer
  i_free_tracer_lbc = my_nml % i_free_tracer_lbc
  l_free_tracer     = my_nml % l_free_tracer
  l_free_tracer_lbc = my_nml % l_free_tracer_lbc
  l_bl_tracer_mix   = my_nml % l_bl_tracer_mix
  l_pv_tracer       = my_nml % l_pv_tracer
  l_pv_dyn          = my_nml % l_pv_dyn
  l_pv_split_rad    = my_nml % l_pv_split_rad
  l_pv_adv_only     = my_nml % l_pv_adv_only
  l_calc_pv_full    = my_nml % l_calc_pv_full
  l_diab_tracer     = my_nml % l_diab_tracer
  l_diab_tr_rad     = my_nml % l_diab_tr_rad
  l_diab_tr_bl      = my_nml % l_diab_tr_bl

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_free_tracers
#endif

END MODULE free_tracers_inputs_mod
