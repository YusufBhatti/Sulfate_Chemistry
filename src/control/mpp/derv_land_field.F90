! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine DERV_LAND_FIELD : Computes no of land points in MPP jobs

! Subroutine Interface :

SUBROUTINE derv_land_field (dump_unit,icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE io
USE Model_File
USE UM_ParParams
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE nlsizes_namelist_mod, ONLY: &
    global_land_field, land_field, local_land_field

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Description : Calculates the no of land points on each PE.

! Method : Call CALC_LAND_FIELD to read in Land-Sea Mask from
!          Atmosphere Dump and then calculate no of land points.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP

! Code Description :
! Language : FORTRAN90

! Declarations :

!     Arguments

! unit attached to input dump file. This may be astart in the first time 
! step, or checkpoint_dump_im for later time steps at the start of
! continuation runs (CRUNs).
INTEGER, INTENT(IN) :: dump_unit 

INTEGER, INTENT(OUT) :: icode          ! OUT Error Code
CHARACTER(LEN=errormessagelength) :: cmessage ! OUT Error message

!     Local variables
INTEGER :: ilen1_lookup   ! First dimesion of look-up table
INTEGER :: ilen2_lookup   ! Second dimension of look-up table
INTEGER :: fixhd(256)     ! Fixed header
INTEGER :: unit_no        ! File unit for dump file

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DERV_LAND_FIELD'



IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!     land_field is the global no of land-points.

!     Initialise global_land_field
global_land_field = land_field

WRITE(umMessage,*) ' global_land_field set to ',land_field
CALL umPrint(umMessage,src='derv_land_field')

!     Link to atmos input dump

unit_no = dump_unit  ! unit attached to input dump file

CALL setpos(unit_no,0,icode)

IF (icode >  0) THEN
  WRITE(umMessage,*) 'Error in SETPOS called from DERV_LAND_FIELD.'
  CALL umPrint(umMessage,src='derv_land_field')
  WRITE(umMessage,*) 'Trying to reset pointer to start of file.'
  CALL umPrint(umMessage,src='derv_land_field')
  GO TO 9999   !  Return
END IF

!     Read fixed header
! DEPENDS ON: read_flh
CALL read_flh(unit_no,fixhd,256,icode,cmessage)

!     Check error code from read_flh
IF (icode >  0) THEN
  WRITE(umMessage,*) 'Error in READ_FLH called from DERV_LAND_FIELD.'
  CALL umPrint(umMessage,src='derv_land_field')
  WRITE(umMessage,*) 'Trying to read fixed header from atmos dump.'
  CALL umPrint(umMessage,src='derv_land_field')
  GO TO 9999   !  Return
END IF

!     Get dimensions of look-up table
ilen1_lookup=fixhd(151)
ilen2_lookup=fixhd(152)

!     Proceed to calculate no of land points on each PE.
! DEPENDS ON: calc_land_field
CALL calc_land_field (unit_no,fixhd,ilen1_lookup,ilen2_lookup,    &
                      icode,cmessage)

!     land_field now contains the no of land_points for this PE.

!     Initialise local_land_field
local_land_field = land_field

WRITE(umMessage,*) ' local_land_field set to ',land_field
CALL umPrint(umMessage,src='derv_land_field')

9999 CONTINUE

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE derv_land_field
