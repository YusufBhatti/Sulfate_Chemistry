! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Load the STASHmaster_A file 

MODULE stashmaster_utils_mod

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

USE ppxlook_mod,        ONLY: read_atmos_stashmaster, exppxi
USE filenamelength_mod, ONLY: filenamelength

! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'STASHMASTER_UTILS_MOD'

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE load_stashmaster(stashmaster_path)
! Given the path to the STASHmaster directory, loads the STASHmaster file
IMPLICIT NONE
CHARACTER(LEN=filenamelength), INTENT(IN) :: stashmaster_path
CHARACTER(LEN=*), PARAMETER :: routinename = 'LOAD_STASHMASTER'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CALL read_atmos_stashmaster(stashmaster_path)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE load_stashmaster

!-------------------------------------------------------------------------------

INTEGER FUNCTION query_stashmaster(stashcode, element)
! Given a stashcode and element (e.g. bottom level code), returns the value of
! that element for that stashcode.
IMPLICIT NONE
INTEGER, INTENT(IN) :: stashcode, element
INTEGER :: section, item
INTEGER, PARAMETER ::  model = 1   ! Atmosphere model only
CHARACTER(LEN=*), PARAMETER :: routinename = 'QUERY_STASHMASTER'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
section = stashcode / 1000
item = MOD(stashcode, 1000)
query_stashmaster = exppxi(model,section,item,element)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION query_stashmaster

!-------------------------------------------------------------------------------

INTEGER FUNCTION number_of_dust_fields(list)
! Count the number of dust bin stash codes in a list of stash codes
USE um_stashcode_mod, ONLY: stashcode_dust1_mmr, stashcode_dust2_mmr, stashcode_dust3_mmr, &
                            stashcode_dust4_mmr, stashcode_dust5_mmr, stashcode_dust6_mmr
IMPLICIT NONE
INTEGER :: list(:)
INTEGER :: i
LOGICAL :: found_bin1, found_bin2, found_bin3, found_bin4, found_bin5, found_bin6
CHARACTER(LEN=*), PARAMETER :: routinename = 'NUMBER_OF_DUST_FIELDS'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
found_bin1 = .FALSE.  
found_bin2 = .FALSE.  
found_bin3 = .FALSE.  
found_bin4 = .FALSE.  
found_bin5 = .FALSE.  
found_bin6 = .FALSE.  
number_of_dust_fields = 0

DO i = 1, SIZE(list)
  IF (list(i) == stashcode_dust1_mmr .AND. .NOT. found_bin1) THEN
    found_bin1 = .TRUE.
    number_of_dust_fields = number_of_dust_fields + 1
  ELSE IF (list(i) == stashcode_dust2_mmr .AND. .NOT. found_bin2) THEN
    found_bin2 = .TRUE.
    number_of_dust_fields = number_of_dust_fields + 1
  ELSE IF (list(i) == stashcode_dust3_mmr .AND. .NOT. found_bin3) THEN
    found_bin3 = .TRUE.
    number_of_dust_fields = number_of_dust_fields + 1
  ELSE IF (list(i) == stashcode_dust4_mmr .AND. .NOT. found_bin4) THEN
    found_bin4 = .TRUE.
    number_of_dust_fields = number_of_dust_fields + 1
  ELSE IF (list(i) == stashcode_dust5_mmr .AND. .NOT. found_bin5) THEN
    found_bin5 = .TRUE.
    number_of_dust_fields = number_of_dust_fields + 1
  ELSE IF (list(i) == stashcode_dust6_mmr .AND. .NOT. found_bin6) THEN
    found_bin6 = .TRUE.
    number_of_dust_fields = number_of_dust_fields + 1
  END IF
END DO      

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION number_of_dust_fields


END MODULE stashmaster_utils_mod

