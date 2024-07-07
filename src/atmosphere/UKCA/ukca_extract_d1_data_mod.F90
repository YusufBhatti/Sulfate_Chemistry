MODULE ukca_extract_d1_data_mod

IMPLICIT NONE

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Subroutines to extract data from D1 array. Converted from
!   REDIST_STOCHEM.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!  Use of an interface block means that a generic call can select
!  the appropriate one of these routines automatically from the
!  input variables. The use of the PRIVATE statement means that
!  the only way to access these routines is via the generic interface.
!
!  Called from UKCA_MAIN1.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------

! Default is private, only generic interface public
PRIVATE
PUBLIC ukca_extract_d1_data

INTERFACE ukca_extract_d1_data
MODULE PROCEDURE                  &
  ukca_extract_d1_data1d,         &
  ukca_extract_d1_data2d,         &
  ukca_extract_d1_integer_data2d, &
  ukca_extract_d1_logical_data2d, &
  ukca_extract_d1_data3d
END INTERFACE ukca_extract_d1_data

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_EXTRACT_D1_DATA_MOD'

CONTAINS

! ---------------------------------------------------------------------
! Subroutine to extract 1D array of REAL data from D1
! ---------------------------------------------------------------------
SUBROUTINE ukca_extract_d1_data1d(                            &
first,n,x)

USE ukca_d1_defs
USE d1_array_mod, ONLY: d1
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY: &
    len_tot, n_obj_d1_max

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN)   :: n        ! id of array
LOGICAL, INTENT(IN)   :: first
REAL, INTENT(INOUT)   :: x(:)     ! extracted array

!       Local variables

CHARACTER (LEN=errormessagelength)    :: cmessage
INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EXTRACT_D1_DATA1D'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (first .AND. SIZE(x) /= UkcaD1Codes(n)%length) THEN
  errcode =  (1000*UkcaD1Codes(n)%section) + UkcaD1Codes(n)%item
  cmessage='Array sizes in local variable and D1 do not agree'
  WRITE(umMessage,'(A50,A13,I8,A5,I8)') cmessage, ' Error code: ', &
    errcode,' PE: ',mype
  CALL umPrint(umMessage,src='ukca_extract_d1_data1d')
  WRITE(umMessage,'(A15,I8,A15,I8)') 'Expected size: ',SIZE(x),    &
    ' Length in D1: ',UkcaD1Codes(n)%length
  CALL umPrint(umMessage,src='ukca_extract_d1_data1d')

  CALL ereport('extract_d1_data2d',errcode,cmessage)
ELSE
  x = d1(UkcaD1Codes(n)%address:UkcaD1Codes(n)%address         &
          +UkcaD1Codes(n)%length-1)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_extract_d1_data1d

! ---------------------------------------------------------------------
! Subroutine to extract 2D array of REAL data from D1
! ---------------------------------------------------------------------
SUBROUTINE ukca_extract_d1_data2d(                             &
first,n,x)

USE ukca_d1_defs
USE d1_array_mod, ONLY: d1
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY: &
    len_tot, n_obj_d1_max
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN)  :: n       ! id of array
LOGICAL, INTENT(IN)  :: first
REAL, INTENT(INOUT)  :: x(:,:)  ! extracted array

REAL, ALLOCATABLE  :: data1(:)
CHARACTER (LEN=errormessagelength) :: cmessage
INTEGER            :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EXTRACT_D1_DATA2D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (first .AND. SIZE(x) /= UkcaD1Codes(n)%length) THEN
  errcode =  (1000*UkcaD1Codes(n)%section)+                    &
       UkcaD1Codes(n)%item
  cmessage='Array sizes in local variable and D1 do not agree'
  WRITE(umMessage,'(A50,A13,I8,A5,I8)') cmessage, ' Error code: ', &
    errcode,' PE: ',mype
  CALL umPrint(umMessage,src='ukca_extract_d1_data2d')
  WRITE(umMessage,'(A15,I8,A15,I8)') 'Expected size: ',SIZE(x),    &
    ' Length in D1: ',UkcaD1Codes(n)%length
  CALL umPrint(umMessage,src='ukca_extract_d1_data2d')
  WRITE(umMessage,'(A19,I8,A2,I8)') 'Error: UKCAD1CODES(',n,')=',  &
    UkcaD1Codes(n)
  CALL umPrint(umMessage,src='ukca_extract_d1_data2d')
  CALL ereport('extract_d1_data2d',errcode,cmessage)
ELSE
  ALLOCATE(data1(UkcaD1Codes(n)%length))
  data1 = d1(UkcaD1Codes(n)%address:UkcaD1Codes(n)%address     &
          +UkcaD1Codes(n)%length-1)
  x = RESHAPE(data1,(/SIZE(x,DIM=1),SIZE(x,DIM=2)/))
  DEALLOCATE(data1)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_extract_d1_data2d

! ---------------------------------------------------------------------
! Subroutine to extract 2D array of INTEGER data from D1
! ---------------------------------------------------------------------
SUBROUTINE ukca_extract_d1_integer_data2d(                     &
 first,n,x)

USE ukca_d1_defs
USE d1_array_mod, ONLY: d1
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY: &
    len_tot, n_obj_d1_max
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN)    :: n       ! id of array
LOGICAL, INTENT(IN)    :: first
INTEGER, INTENT(INOUT) :: x(:,:)  ! extracted array

INTEGER, ALLOCATABLE :: data1(:)
CHARACTER(LEN=errormessagelength)    :: cmessage
INTEGER              :: errcode
INTEGER              :: int_val   ! used for transfers

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EXTRACT_D1_INTEGER_DATA2D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (first .AND. SIZE(x) /= UkcaD1Codes(n)%length) THEN
  errcode = n
  cmessage='Array sizes in local variable and D1 do not agree'
  WRITE(umMessage,'(A50,A13,I8,A5,I8)') cmessage, ' Error code: ', &
    errcode,' PE: ',mype
  CALL umPrint(umMessage,src='ukca_extract_d1_integer_data2d')
  WRITE(umMessage,'(A15,I8,A15,I8)') 'Expected size: ',SIZE(x),    &
    ' Length in D1: ',UkcaD1Codes(n)%length
  CALL umPrint(umMessage,src='ukca_extract_d1_integer_data2d')

  CALL ereport('extract_d1_data2d',errcode,cmessage)
ELSE
  ALLOCATE(data1(UkcaD1Codes(n)%length))
  data1 = TRANSFER( d1(UkcaD1Codes(n)%address:UkcaD1Codes(n)%address    &
                       +UkcaD1Codes(n)%length-1), int_val,              &
                    SIZE=UkcaD1Codes(n)%length )
  x = RESHAPE(data1,(/SIZE(x,DIM=1),SIZE(x,DIM=2)/))
  DEALLOCATE(data1)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_extract_d1_integer_data2d

! ---------------------------------------------------------------------
! Subroutine to extract 2D array of LOGICAL data from D1
! ---------------------------------------------------------------------
SUBROUTINE ukca_extract_d1_logical_data2d(                     &
 first,n,x)

USE ukca_d1_defs
USE d1_array_mod, ONLY: d1
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY: &
    len_tot, n_obj_d1_max

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN)    :: n       ! id of array
LOGICAL, INTENT(IN)    :: first
LOGICAL, INTENT(INOUT) :: x(:,:)  ! extracted array

LOGICAL, ALLOCATABLE   :: data1(:)
CHARACTER(LEN=errormessagelength)      :: cmessage
INTEGER                :: errcode
LOGICAL                :: log_val ! used for transfers

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EXTRACT_D1_LOGICAL_DATA2D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (first .AND. SIZE(x) /= UkcaD1Codes(n)%length) THEN
  errcode = n
  cmessage='Array sizes in local variable and D1 do not agree'
  WRITE(umMessage,'(A50,A13,I8,A5,I8)') cmessage, ' Error code: ', &
    errcode,' PE: ',mype
  CALL umPrint(umMessage,src='ukca_extract_d1_logical_data2d')
  WRITE(umMessage,'(A15,I8,A15,I8)') 'Expected size: ',SIZE(x),    &
    ' Length in D1: ',UkcaD1Codes(n)%length
  CALL umPrint(umMessage,src='ukca_extract_d1_logical_data2d')

  WRITE(cmessage,'(A)')'Check if domain profiles DTILE and DPFT'   &
      //'match with NTILE and NPFT'
  CALL ereport('extract_d1_data2d',errcode,cmessage)
ELSE
  ALLOCATE(data1(UkcaD1Codes(n)%length))
  data1 = TRANSFER( d1(UkcaD1Codes(n)%address:UkcaD1Codes(n)%address    &
                       +UkcaD1Codes(n)%length-1), log_val,              &
                    SIZE=UkcaD1Codes(n)%length )
  x = RESHAPE(data1,(/SIZE(x,DIM=1),SIZE(x,DIM=2)/))
  DEALLOCATE(data1)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_extract_d1_logical_data2d

! ---------------------------------------------------------------------
! Subroutine to extract 3D array of REAL data from D1
! ---------------------------------------------------------------------
SUBROUTINE ukca_extract_d1_data3d(                             &
 first,n,x)

USE ukca_d1_defs
USE d1_array_mod, ONLY: d1
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE nlsizes_namelist_mod, ONLY: &
    len_tot, n_obj_d1_max
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER, INTENT(IN) :: n         ! id of array
LOGICAL, INTENT(IN) :: first
REAL, INTENT(INOUT) :: x(:,:,:)  ! extracted array

REAL, ALLOCATABLE   :: data1(:)
CHARACTER(LEN=errormessagelength)   :: cmessage
INTEGER             :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_EXTRACT_D1_DATA3D'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (first .AND. SIZE(x) /= UkcaD1Codes(n)%length) THEN
  errcode =  (1000*UkcaD1Codes(n)%section)+&
       UkcaD1Codes(n)%item
  cmessage='Array sizes in local variable and D1 do not agree'
  WRITE(umMessage,'(A50,A13,I8,A5,I8)') cmessage, ' Error code: ', &
    errcode,' PE: ',mype
  CALL umPrint(umMessage,src='ukca_extract_d1_data3d')
  WRITE(umMessage,'(A15,I8,A15,I8)') 'Expected size: ',SIZE(x),    &
    ' Length in D1: ',UkcaD1Codes(n)%length
  CALL umPrint(umMessage,src='ukca_extract_d1_data3d')

  WRITE(cmessage,'(A)')'Check if domain profiles DTILE and DPFT'   &
      //' match with NTILE and NPFT'
  CALL ereport('extract_d1_data3d',errcode,cmessage)
ELSE
  ALLOCATE(data1(UkcaD1Codes(n)%length))
  data1 = d1(UkcaD1Codes(n)%address:UkcaD1Codes(n)%address     &
          +UkcaD1Codes(n)%length-1)
  x=RESHAPE(data1,(/SIZE(x,DIM=1),SIZE(x,DIM=2),SIZE(x,DIM=3)/))
  DEALLOCATE(data1)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_extract_d1_data3d
END MODULE ukca_extract_d1_data_mod
