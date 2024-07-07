! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Description:  Dummy routines to break dependencies in SSTPert library
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: SSTpert library

SUBROUTINE gc_ibcast(arg1, arg2, arg3, arg4, arg5, arg6)
IMPLICIT NONE 
INTEGER, INTENT(IN) :: arg1, arg2, arg3, arg4, arg5
INTEGER, INTENT(IN) :: arg6(:)
END SUBROUTINE gc_ibcast

MODULE um_parvars
IMPLICIT NONE 
INTEGER, PARAMETER :: mype = 0
INTEGER, PARAMETER :: nproc = 1
END MODULE um_parvars

MODULE dump_headers_mod
IMPLICIT NONE 
INTEGER, PARAMETER :: ih_sp_seed = 32 
INTEGER, ALLOCATABLE :: a_inthd(:)
INTEGER, PARAMETER :: ih_stph_n1 = 30
INTEGER, PARAMETER :: ih_stph_n2 = 31
INTEGER, PARAMETER :: ih_stochastic_flag = 29
END MODULE dump_headers_mod

MODULE umPrintMgr
USE ISO_FORTRAN_ENV, ONLY: OUTPUT_UNIT
IMPLICIT NONE
INTEGER, PARAMETER :: maxLineLen = 1024
CHARACTER (LEN=maxLineLen)    :: umMessage
INTEGER, PARAMETER :: prstatus_normal = 1
INTEGER :: printstatus = 0
LOGICAL :: printeractive = .FALSE.
!$OMP THREADPRIVATE (umMessage)
CHARACTER(LEN=1)              :: newline
CONTAINS
RECURSIVE SUBROUTINE umPrint(line,level,pe,src,model,                          &
                             UsrPrefix,HangIndent,stdErrorToo)
IMPLICIT NONE 
CHARACTER(LEN=*)            :: line
INTEGER, OPTIONAL           :: level
INTEGER, OPTIONAL           :: pe
CHARACTER(LEN=*), OPTIONAL  :: src
CHARACTER(LEN=*), OPTIONAL  :: model
CHARACTER(LEN=*), OPTIONAL  :: UsrPrefix  
CHARACTER(LEN=*), OPTIONAL  :: HangIndent 
LOGICAL, OPTIONAL           :: stdErrorToo
WRITE(OUTPUT_UNIT, "(A)") TRIM(line)
END SUBROUTINE umPrint
SUBROUTINE umprintFlush()
IMPLICIT NONE 
END SUBROUTINE umprintFlush
END MODULE umPrintMgr

MODULE stochastic_physics_run_mod
USE missing_data_mod, ONLY: imdi
IMPLICIT NONE 
LOGICAL :: l_x_eq_sin_x= .TRUE.
INTEGER, PARAMETER :: firsttimestep_true  = 1
INTEGER, PARAMETER :: firsttimestep_crun  = -1
INTEGER :: stph_n1 = imdi
INTEGER :: stph_n2 = imdi
REAL, ALLOCATABLE :: Ymn(:,:,:)
INTEGER :: stphseed = 1
INTEGER, PARAMETER :: stph_seed_present          = 1
INTEGER, PARAMETER :: stph_spt_data_present      = 2
INTEGER, PARAMETER :: stph_skeb2_data_present    = 4
INTEGER, PARAMETER :: stph_rp2_data_present      = 8
INTEGER :: stph_skeb2_data_check = imdi
INTEGER :: stph_header_flag = stph_seed_present  
END MODULE stochastic_physics_run_mod
