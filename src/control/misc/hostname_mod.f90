! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Purpose:
!    Return the current hostname using the C function gethostname()
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Misc


MODULE hostname_mod

IMPLICIT NONE


CONTAINS

FUNCTION get_hostname()

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_SIZE_T
  USE f_shum_string_conv_mod, ONLY: f_shum_c2f_string
  IMPLICIT NONE

  INTERFACE 
    SUBROUTINE gethostname(name, length) BIND(C, NAME='gethostname')
        USE ISO_C_BINDING, ONLY: C_CHAR, C_SIZE_T
        INTEGER (C_SIZE_T), VALUE, INTENT(IN) :: length
        CHARACTER(KIND=C_CHAR, LEN=1), INTENT(OUT) :: name(*)
    END SUBROUTINE gethostname
  END INTERFACE

  INTEGER, PARAMETER :: hostname_length = 256
  CHARACTER(KIND=C_CHAR, LEN=1) :: hostname(hostname_length+1)
  CHARACTER(LEN=hostname_length) :: get_hostname
  
  CALL gethostname(hostname, INT(hostname_length+1,KIND=C_SIZE_T))

  get_hostname = TRIM(f_shum_c2f_string(hostname))
  
END FUNCTION get_hostname

END MODULE hostname_mod
