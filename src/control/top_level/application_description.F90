! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine: Application_Description -----------------------------------
!
!  Purpose: Provides runtime registry and enquiry functions allowing
!           UM code to determine the executable target (program main),
!           build type, sizes, and parallel nature. This facilitates
!           reduction in preprocessor based compilation in favour of
!           software switching
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Top Level

MODULE Application_Description

USE application_description_runtypes

IMPLICIT NONE

CONTAINS

SUBROUTINE setApplicationDesc(appType) ! Call once from main.
USE UM_Types
USE ereport_mod

IMPLICIT NONE
INTEGER, INTENT(IN) :: appType
INTEGER             :: rtype=4
INTEGER             :: itype=4
INTEGER             :: ptype=4
INTEGER             :: errorCode
CHARACTER(LEN=132)  :: message

errorCode=1

IF (appType>0 .AND. appType <= num_exe_types ) THEN
  exe_type=appType
ELSE
  WRITE(message,'(A,I3)') &
      'Tried to set application to unknown value:',appType
  CALL ereport('applcation_description:setApplicationDescription', &
      errorCode,message)
END IF

rtype= umFortranRealSize()
itype= umFortranIntegerSize()
ptype= umFortranPointerSize()

IF (rtype==8 .AND. itype==8) THEN
  exe_data_64=.TRUE.
ELSE IF (rtype==4 .AND. itype==4) THEN
  exe_data_64=.FALSE.
ELSE
  WRITE(message,'(A)') &
      'Application has different lengths for reals and integers'
  CALL ereport('applcation_description:setApplicationDescription', &
      errorCode,message)
END IF

IF (ptype==8) THEN
  exe_addr_64=.TRUE.
ELSE IF (ptype==4) THEN
  exe_data_64=.FALSE.
ELSE
  WRITE(message,'(A)') &
      'Application has an unknown pointer size'
  CALL ereport('applcation_description:setApplicationDescription', &
      errorCode,message)
END IF

END SUBROUTINE setApplicationDesc

! Provide the application type
FUNCTION getExeType() RESULT(r)
IMPLICIT NONE
INTEGER :: r
r=exe_type
END FUNCTION getExeType

! Is this an MPI program (ie I can call MPL, rather than nproc>1)
FUNCTION isParallel() RESULT(r)
USE ereport_mod
IMPLICIT NONE
LOGICAL :: r
INTEGER :: istat
INTEGER :: val
INTEGER :: errorCode
CHARACTER(LEN=132)  :: message
INTEGER, PARAMETER :: gc_parallel=1
INTEGER, PARAMETER :: gc_serial=2
INTEGER, PARAMETER :: gc_is_parallel=3
CALL GC_GetOpt(gc_is_parallel,val,istat)
IF (val==gc_serial) THEN
  r=.FALSE.
ELSE IF (val==gc_parallel) THEN
  r=.TRUE.
ELSE
  errorCode=2
  WRITE(message,'(A,I3)') &
      'GCOM returned an unknown value for the GC_IS_PARALLEL check:', &
      val
  CALL ereport('applcation_description:isParallel', &
      errorCode,message)
END IF
END FUNCTION isParallel

! Convenience routines, to taste.
FUNCTION isSmallExec() RESULT(r)
IMPLICIT NONE
LOGICAL :: r
r=.FALSE.
IF (exe_type==1) THEN
  r=.TRUE.
ELSE IF (exe_type>2) THEN
  r=.TRUE.
END IF

END FUNCTION isSmallExec

FUNCTION UM_NameStr(i) RESULT(r)
IMPLICIT NONE
INTEGER, INTENT(IN)           :: i
CHARACTER(LEN=UM_Name_Len)    :: r
WRITE(r,'(A)')UM_name((exe_type-1)*UM_Name_Len+1:i*UM_Name_len)
END FUNCTION UM_NameStr


END MODULE Application_Description
