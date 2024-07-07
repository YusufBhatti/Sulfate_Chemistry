MODULE oasis_operations_mod
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='OASIS_OPERATIONS_MOD'

CONTAINS

LOGICAL FUNCTION oasis_get_ok(recv_code)
!
! Description:
! Checks the incoming recv_code against a list of valid OASIS return codes 
! to determine if a successful receive (get) operation has been completed. 
! Returns TRUE if the code indicates a successful get, FALSE otherwise.  
!  
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Coupling
!=======================================================================
USE um_types, ONLY: integer32

! OASIS3-MCT interface return codes
USE mod_prism, ONLY:       prism_recvout,     &
                           prism_fromrestout, &
                           prism_recvd,       &
                           prism_fromrest

IMPLICIT NONE

INTEGER (KIND=integer32), INTENT(IN) :: recv_code

OASIS_get_ok = ((recv_code == prism_recvout)    .OR. &
                (recv_code == prism_fromrestout).OR. &
                (recv_code == prism_recvd)      .OR. &
                (recv_code == prism_fromrest))

END FUNCTION oasis_get_ok

!=======================================================================

SUBROUTINE read_nml_transcount(unit_in, errorstatus)

USE um_types, ONLY: integer64
USE um_parcore, ONLY: mype
USE oasis_atm_data_mod, ONLY: transcount, &
                              transient_a2o_count, transient_o2a_count, &
                              transient_a2c_count, transient_c2a_count

USE setup_namelist, ONLY: setup_nml_type

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE mpl, ONLY: mpl_integer

IMPLICIT NONE

INTEGER,INTENT(IN) :: unit_in
INTEGER,INTENT(OUT) :: ErrorStatus
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_TRANSCOUNT'

! set number of types and number of each type of variable in my_namelist type

INTEGER, PARAMETER :: no_of_types = 1 ! We only deal with one type (integer).
INTEGER, PARAMETER :: n_int = 4       ! The number of integers to broadcast.

! setup namelist type - variables of the same type in the namelist
! should be together
! variables should be ordered - INTEGER, REAL, LOGICAL, CHARACTER
! each variable in the namelist should be declared in the TYPE statement

TYPE my_namelist
  SEQUENCE
  INTEGER (KIND=integer64) :: transient_a2o_count
  INTEGER (KIND=integer64) :: transient_o2a_count
  INTEGER (KIND=integer64) :: transient_a2c_count
  INTEGER (KIND=integer64) :: transient_c2a_count
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL gc_get_communicator(my_comm, icode)

!remove optional arguments that correspond to n_type=0
CALL setup_nml_type(no_of_types, mpl_nml_type,                      &
                    n_int_in=n_int)

! read in namelist on processor 0, populate namelist type variable (my_nml)
! for convenience do in same order as defined in TYPE above

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=transcount, IOSTAT=ErrorStatus)

  my_nml %  transient_a2o_count = transient_a2o_count
  my_nml %  transient_o2a_count = transient_o2a_count
  my_nml %  transient_a2c_count = transient_a2c_count
  my_nml %  transient_c2a_count = transient_c2a_count

END IF

! broadcast namelist type variable (my_nml) to all other processors

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)
CALL mpl_bcast(Errorstatus,1,mpl_integer,0,my_comm,icode)

! all processors except 0 populate namelist using namelist type
! variable (my_nml)

IF (mype /= 0) THEN

  transient_a2o_count = my_nml %  transient_a2o_count
  transient_o2a_count = my_nml %  transient_o2a_count
  transient_a2c_count = my_nml %  transient_a2c_count
  transient_c2a_count = my_nml %  transient_c2a_count

END IF

! free up memory by freeing the mpl type

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_transcount

!=======================================================================

SUBROUTINE read_nml_transfld(unit_in,ErrorStatus)

USE um_parcore, ONLY: mype
USE mpl, ONLY: mpl_integer
USE oasis_atm_data_mod, ONLY: transfld, & 
                              i_number, &
                              f_id,     &
                              map_type, &
                              max_lev,  &
                              c_out,    &
                              c_in,     &
                              n_out,    &
                              n_in

USE setup_namelist, ONLY: setup_nml_type

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! if unit number is hard-wired can just hard-wire later in routine
INTEGER,INTENT(IN) :: unit_in
INTEGER,INTENT(OUT) :: ErrorStatus
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: icode

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_TRANSFLD'

! set number of types and number of each type of variable in my_namelist type

INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 4
INTEGER, PARAMETER :: n_chars = 4 + 4 + 8 + 8

! setup namelist type - variables of the same type in the namelist
! should be together
! variables should be ordered - INTEGER, REAL, LOGICAL, CHARACTER
! each variable in the namelist should be declared in the TYPE statement

TYPE my_namelist
  SEQUENCE
  INTEGER :: i_number
  INTEGER :: f_id
  INTEGER :: map_type 
  INTEGER :: max_lev
  CHARACTER(LEN=4) :: c_out
  CHARACTER(LEN=4) :: c_in
  CHARACTER(LEN=8) :: n_out
  CHARACTER(LEN=8) :: n_in
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL gc_get_communicator(my_comm, icode)

!remove optional arguments that correspond to n_type=0
CALL setup_nml_type(no_of_types, mpl_nml_type,                      &
                    n_int_in=n_int,n_chars_in=n_chars)

! read in namelist on processor 0, populate namelist type variable (my_nml)
! for convenience do in same order as defined in TYPE above

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=transfld, IOSTAT=ErrorStatus)

  my_nml % i_number = i_number
  my_nml % f_id     = f_id
  my_nml % map_type = map_type
  my_nml % max_lev  = max_lev
  my_nml % c_out    = c_out
  my_nml % c_in     = c_in
  my_nml % n_out    = n_out
  my_nml % n_in     = n_in

END IF

! broadcast namelist type variable (my_nml) to all other processors

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)
CALL mpl_bcast(Errorstatus,1,mpl_integer,0,my_comm,icode)

! all processors except 0 populate namelist using namelist type
! variable (my_nml)

IF (mype /= 0) THEN

  i_number = my_nml % i_number         
  f_id     = my_nml % f_id             
  map_type = my_nml % map_type
  max_lev  = my_nml % max_lev        
  c_out    = my_nml % c_out
  c_in     = my_nml % c_in 
  n_out    = my_nml % n_out
  n_in     = my_nml % n_in 

END IF

! free up memory by freeing the mpl type

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_transfld

END MODULE oasis_operations_mod
