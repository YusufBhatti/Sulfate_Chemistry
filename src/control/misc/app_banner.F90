! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine: Application_Description -----------------------------------
!
!  Purpose:
!    Report the contents of the application description module.
!    The routine can't be inside that module because of printing dependencies.
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Misc

MODULE App_Banner

USE Application_Description
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
IMPLICIT NONE

CONTAINS

SUBROUTINE reportApplication()
IMPLICIT NONE
INTEGER           :: i
CHARACTER(LEN=2)  :: ds
CHARACTER(LEN=2)  :: ps
CALL umPrint(' ******************************************',src='app_banner')
WRITE(umMessage,'(A,I2,A,A)')    &
    ' App ID:',exe_type, &
    ', Name:',TRIM(UM_NameStr(exe_type))
CALL umPrint(umMessage,src='app_banner')
CALL umPrint(' ------------------------------------------',src='app_banner')

WRITE(ds,'(A)')'32'
WRITE(ps,'(A)')'32'
IF (exe_data_64) WRITE(ds,'(A)')'64'
IF (exe_addr_64) WRITE(ps,'(A)')'64'
WRITE(umMessage,'(A,A,A,A,A)')'  - Data size is ',ds, &
    ' bit. Program is ',ps,' bit.'
CALL umPrint(umMessage,src='app_banner')
IF (isParallel()) THEN
  WRITE(umMessage,'(A,I2)')'  - Program is parallel.'
  CALL umPrint(umMessage,src='app_banner')
ELSE
  WRITE(umMessage,'(A,I2)')'  - Program is serial.'
  CALL umPrint(umMessage,src='app_banner')
END IF
CALL umPrint(' ******************************************',src='app_banner')

END SUBROUTINE reportApplication

END MODULE App_Banner
