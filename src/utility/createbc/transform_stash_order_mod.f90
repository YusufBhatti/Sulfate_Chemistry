! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE transform_stash_order_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.


! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! Variable to contain the acceptable list of stash codes in order
INTEGER, ALLOCATABLE :: valid_stash_code_order(:)

! Number of potential stash codes
INTEGER :: len_stash_code

! Flag to confirm valid_stash_code_order has been initialised
LOGICAL :: initialised = .FALSE.

! Non-tracer stash codes, in order
INTEGER :: master_stash(43) = (/     &
                                33,  & ! Orogaphy
                                2,   & ! u-wind 
                                3,   & ! v-wind
                                150, & ! w-wind 
                                253, & ! rho
                                4,   & ! theta
                                10,  & ! q
                                254, & ! qcl
                                12,  & ! qcf
                                255, & ! exner
                                256, & ! adv u-wind
                                257, & ! adv v-wind
                                258, & ! adv w-wind
                                271, & ! qcf2
                                272, & ! qrain
                                273, & ! qgraup
                                266, & ! bulk_cf
                                267, & ! liquid_cf
                                268, & ! frozen_cf
                                90,  & ! total_aero
                                431, & ! dust mmr bin 1
                                432, & ! dust mmr bin 2
                                433, & ! dust mmr bin 3
                                434, & ! dust mmr bin 4
                                435, & ! dust mmr bin 5
                                436, & ! dust mmr bin 6
                                101, & ! so2
                                102, & ! dms
                                103, & ! mmr so4_aitken
                                104, & ! mmr so4_accum
                                105, & ! mmr so4_diss
                                107, & ! mmr nh3
                                108, & ! mmr bc_fr
                                109, & ! mmr bc_ag
                                110, & ! mmr bc_cl
                                111, & ! mmr smoke fr
                                112, & ! mmr smoke ag
                                113, & ! mmr smoke cl
                                114, & ! mmr ocff_fr
                                115, & ! mmr ocff_ag
                                116, & ! mmr ocff_cl
                                117, & ! mmr nitr acc
                                118  & ! mmr nitr diss
                                /)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'TRANSFORM_STASH_ORDER_MOD'

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE initialise_master_stash()
! Fill the valid_stash_code_order variable, based on the master_stash list
! above and the number of potential tracers
USE um_stashcode_mod,           ONLY: stashcode_tracer_sec, stashcode_ukca_sec
USE lbc_stashcode_mapping_mod, ONLY: num_ukca_tracers, num_free_tracers
IMPLICIT NONE

! Loop counter
INTEGER :: i
CHARACTER(LEN=*), PARAMETER :: routinename = 'INITIALISE_MASTER_STASH'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set number of potential stash codes  
len_stash_code = 43 + num_ukca_tracers + num_free_tracers
  
ALLOCATE(valid_stash_code_order(len_stash_code))
  
! first 43 are the master_stash items
valid_stash_code_order(1:43) = master_stash(1:43)
  
! Fill remainder with free tracers, then UKCA tracers
DO i = 1, num_free_tracers
  valid_stash_code_order(43+i) = i + stashcode_tracer_sec * 1000
END DO
  
DO i = 1, num_ukca_tracers
  valid_stash_code_order(43+num_free_tracers+i) = i + stashcode_ukca_sec * 1000
END DO

initialised = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE initialise_master_stash

!-------------------------------------------------------------------------------

SUBROUTINE transform_stash_order(num, stash_codes)
! Transform a list of stash_codes into the order the UM requires to read in
! successfully.
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: num
INTEGER, INTENT(INOUT) :: stash_codes(num)

! Counters
INTEGER :: i, j

! Temporary stash codes array
INTEGER :: tmp_stash_codes(len_stash_code)

! Error handling variables
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: routinename = 'TRANSFORM_STASH_ORDER'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise valid_stash_code_order if not already done
IF (.NOT. initialised) CALL initialise_master_stash()

! Create a temporary copy of the array
tmp_stash_codes = valid_stash_code_order

! For every item in the list of potential stash codes, if it's not present
! in the list sent as a subroutine argument, set it to be -1
DO i = 1, len_stash_code
  IF (.NOT. is_in(valid_stash_code_order(i), stash_codes)) THEN
    tmp_stash_codes(i) = -1
  END IF
END DO

! For each stash code in the tmp_stash_codes which isn't -1, set the stash_codes
! argument to be equal to that code.
j = 0
DO i = 1, len_stash_code
  IF (tmp_stash_codes(i) == -1) THEN
    CYCLE
  END IF
  j = j + 1
  stash_codes(j) = tmp_stash_codes(i)
END DO

! Check that the number of stash codes ordered equals the number of stash codes
! sent down as an argument.
IF (j /= num) THEN
  cmessage = 'Error ordering stash codes'
  icode = 110
  CALL ereport(routinename, icode, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE transform_stash_order

!-------------------------------------------------------------------------------

! Return true is value is in list
LOGICAL FUNCTION is_in(val, list)
IMPLICIT NONE
INTEGER :: val
INTEGER :: list(:)
CHARACTER(LEN=*), PARAMETER :: routinename = 'IS_IN'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

is_in = .FALSE.
IF (ANY(list == val)) THEN
  is_in = .TRUE.
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION is_in

END MODULE transform_stash_order_mod
