! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   -----------------------------------------------------------
!
!   Purpose: Generate extra data for the timeseries.
!   This extra data provides information about what processing was done
!   to produce the timeseries. This information will hopefully be of 
!   use to users doing further processing of the timeseries data.
!   This subroutine
!   EXTRA_MAKE_VECTOR:  computes the long/lat ht domain info
!                     : and puts that into the correct place in the
!                     : extra data
!
!   To some extent this routine has much in common with the
!   multi_spatial routine but as it has a different function
!   viz generate info on timeseries rather than generating a single time
!   for the timeseries it is coded separately.
!   When modifying multi_spatial be sure also to modify this routine and
!   vice versa
!   Programming Standard: UM DOC Paper3, Verion 4 (05/02/92)
!
!   System Task:C4
!
!   Interface and arguments ------------------------------------
!
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: STASH


!   Interface and arguments: -----------------------------
SUBROUTINE extra_make_vector(control,control_len,record_cnt,      &
  no_records,extra_data,extra_data_len,                           &
   bzx,bzy,bdx,bdy)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE stparam_mod, ONLY: st_south_code, st_west_code, st_north_code,&
                       st_east_code, st_input_bottom, st_input_top

IMPLICIT NONE
INTEGER :: control_len ! IN size of control record
INTEGER :: control(control_len) ! IN stash control record
INTEGER :: record_cnt ! IN record that is being processed
INTEGER :: no_records ! IN total number of records
INTEGER :: extra_data_len ! IN size of extra data
REAL :: extra_data(extra_data_len) !IN/OUT extra data
REAL :: bdx,bdy,bzx,bzy ! IN grid descriptors
!  ------------------------------------------------------
!
!   Local variables
INTEGER :: addr ! what address in extra data are we at
INTEGER :: record_len ! how many words in a block ?

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EXTRA_MAKE_VECTOR'
!
IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
record_len=no_records+1
addr=1+record_cnt
!  put in the first latitude
extra_data(addr)=control(st_south_code)*bdy+bzy
addr=addr+record_len
!  put in the first long
extra_data(addr)=control(st_west_code)*bdx+bzx
addr=addr+record_len
!  put in the second lat
extra_data(addr)=control(st_north_code)*bdy+bzy
addr=addr+record_len
!  put in the second long
extra_data(addr)=control(st_east_code)*bdx+bzx
addr=addr+record_len
!  put in the lowest level
extra_data(addr)=control(st_input_bottom)
addr=addr+record_len
!  and now the highest level
extra_data(addr)=control(st_input_top)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE extra_make_vector

