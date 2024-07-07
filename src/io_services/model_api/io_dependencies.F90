! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! A module to contain downstream data dependencies
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: IO Services

MODULE io_dependencies

! Older c layer
#if defined(C95_2A)
! DEPENDS ON: portio2a
! DEPENDS ON: pio_io_timer
#endif

! Newer C99 based layer
#if defined(C95_2B)
USE umPrintMgr
USE thread_utils
! DEPENDS ON: portio2b
! DEPENDS ON: c_io
! DEPENDS ON: c_io_timing
! DEPENDS ON: c_io_byteswap
! DEPENDS ON: c_io_rbuffering
! DEPENDS ON: c_io_wbuffering
! DEPENDS ON: c_io_unix
! DEPENDS ON: c_io_libc
! DEPENDS ON: c_io_trace
#endif

! DEPENDS ON: pio_data_conv
! DEPENDS ON: pio_byteswap
! DEPENDS ON: portutils

#if defined (__linux__)
! DEPENDS ON: c_affinity
#endif

IMPLICIT NONE

END MODULE io_dependencies
