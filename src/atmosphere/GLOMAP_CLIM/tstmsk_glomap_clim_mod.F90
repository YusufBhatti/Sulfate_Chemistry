! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Checks version mask and option code for section 54 
!  (GLOMAP-mode climatology aerosols)

MODULE tstmsk_glomap_clim_mod

IMPLICIT NONE

! Description:
!   This subroutine checks if an GLOMAP-mode climatology aerosol variable 
!   (section 54) is turned on or not. For variables which are on, it 
!   sets the logical lmask to true.
!
! Method:
!   lmask is set to false at the start of the routine.
!
!   If GLOMAP-mode climatology is off, the code exits.
!   If GLOMAP-mode climatology is on, then the option codes which are
!   inputs to this subroutine are tested.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: Fortran 95.
!   This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TSTMSK_GLOMAP_CLIM_MOD'

CONTAINS

SUBROUTINE tstmsk_glomap_clim(  n1, n2, n3, n4, n5, n6, n7, n8, n9,n10,        &
                               n11,n12,n13,n14,n15,n16,n17,n18,n19,n20,        &
                               n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,lmask )

USE ereport_mod,            ONLY: &
    ereport

USE errormessagelength_mod, ONLY: &
    errormessagelength

USE glomap_clim_option_mod, ONLY: &
    l_glomap_mode_clim,           &
    i_glomap_clim_setup,          &
    i_gc_sussocbc_5mode

USE parkind1,               ONLY: &
    jpim,                         &
    jprb

USE yomhook,                ONLY: &
    lhook,                        &
    dr_hook

! Note - before adding to the UM we need to add DrHook

IMPLICIT NONE

! Control is by means of 30 option codes n1 to n30 for each
! STASHmaster entry

! 30 option codes are inputs to the routine
INTEGER, INTENT(IN) ::  n1, n2, n3, n4, n5, n6, n7, n8, n9,n10
INTEGER, INTENT(IN) :: n11,n12,n13,n14,n15,n16,n17,n18,n19,n20
INTEGER, INTENT(IN) :: n21,n22,n23,n24,n25,n26,n27,n28,n29,n30

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Local variables
INTEGER                       :: errcode   ! error code for ereport
INTEGER                       :: ntotal    ! sum of option codes

CHARACTER(LEN=*), PARAMETER       :: RoutineName='TSTMSK_GLOMAP_CLIM'
CHARACTER(LEN=errormessagelength) :: cmessage

LOGICAL, INTENT(OUT)          :: lmask

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! By default GLOMAP-mode climatology variables are off
lmask=.FALSE.

! If all option codes are zero then we set lmask to TRUE
! for backward compatibility with other tstmsk sections

ntotal=  n1 +  n2 +  n3 +  n4 +  n5 +  n6 +  n7 +  n8 +  n9 + n10 +            &
        n11 + n12 + n13 + n14 + n15 + n16 + n17 + n18 + n19 + n20 +            &
        n21 + n22 + n23 + n24 + n25 + n26 + n27 + n28 + n29 + n30

IF (ntotal == 0) THEN
  lmask=.TRUE.
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! If GLOMAP-mode climatology is off, all related items are turned off so return
IF (.NOT. l_glomap_mode_clim) THEN
   IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
   RETURN
END IF

SELECT CASE (n30)
CASE (1)
  ! If n30=1 it is controlled by the value of i_glomap_clim_setup
  
  ! we now call check_glomap with the correct option code.
  SELECT CASE (i_glomap_clim_setup)
  CASE (i_gc_sussocbc_5mode) ! SUSSOCBC_5MODE
    CALL check_glomapclim(n2, lmask)
    ! Check option code 2
  CASE DEFAULT
    ! This is an error as there is no code for this glomap-mode scheme
    errcode  = 54000 + i_glomap_clim_setup
    cmessage = 'Unexpected value of i_glomap_clim_setup. Must == 2. In ' //   &
                routinename
    CALL ereport(routinename,errcode,cmessage)
  END SELECT
CASE DEFAULT
   ! This is an error, call ereport
   errcode  = 54999
   cmessage = 'Unexpected value of n30 in tstmsk_glomap_clim'
   CALL ereport(routinename,errcode,cmessage)
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE tstmsk_glomap_clim

! -----------------------------------------------------------------------------
! The same logic is used for all values of i_glomap_clim_setup
! but check a different option code 
! -----------------------------------------------------------------------------
SUBROUTINE check_glomapclim(n, lmask)

USE ereport_mod,            ONLY: &
    ereport

USE errormessagelength_mod, ONLY: &
    errormessagelength

USE glomap_clim_option_mod, ONLY: &
    l_glomap_clim_radaer,         &
    l_glomap_clim_arg_act,        &
    l_glomap_clim_aie1,           &
    l_glomap_clim_aie2

USE parkind1,               ONLY: &
    jpim,                         &
    jprb

USE yomhook,                ONLY: &
    lhook,                        &
    dr_hook

IMPLICIT NONE

! option code being tested
INTEGER,INTENT(IN)    :: n
! logical value returned by this subroutine
LOGICAL,INTENT(INOUT) :: lmask

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Local variables
INTEGER                           :: errcode   ! error code for ereport
CHARACTER(LEN=errormessagelength) :: cmessage  ! used for ereport
CHARACTER(LEN=*), PARAMETER       :: RoutineName='CHECK_GLOMAPCLIM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Check option code n
SELECT CASE (n)
CASE (0)
  ! If n is zero, this is always off
  lmask = .FALSE.
CASE (1)
  ! If n is 1 then it is always on
  lmask = .TRUE.
CASE (3)
  ! If n is 3 then it is on if ACTIVATE is on
  lmask = l_glomap_clim_arg_act
CASE (4)
  ! If n is 4 then it is on if ACTIVATE is on or aie1 or aie2 is on
  lmask = l_glomap_clim_arg_act .OR. l_glomap_clim_aie1 .OR. l_glomap_clim_aie2
CASE (5)
  ! If n is 5 then it is on if RADAER is on
  lmask = l_glomap_clim_radaer
CASE DEFAULT
    ! This is an error, call ereport
  errcode = 54777
  cmessage='Unknown option code in '//routinename
  CALL ereport(routinename,errcode,cmessage)
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE check_glomapclim

END MODULE tstmsk_glomap_clim_mod
