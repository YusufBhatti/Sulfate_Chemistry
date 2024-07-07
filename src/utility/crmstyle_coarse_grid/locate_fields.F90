! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Routine to locate required fields in header

MODULE locate_fields_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LOCATE_FIELDS_MOD'

CONTAINS

SUBROUTINE locate_fields(num_ff,num_want, stash_list,lev_list, date_typ,     &
                          date1, date2, ff_hdr,                              &
                          fstart_num,fend_num,file_num)


USE IO_Mod, ONLY:         &
  PP_Header_type, UM_Header_type
USE crmstyle_cntl_mod, ONLY:     &
   l_sect30, num_want_no30


USE umPrintMgr                              ! for writing output
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Description:
!   This subroutine reads a given number of fields, matching given
!   criteria from an open UM fieldsfile.  The fields are then store in
!   given locations in the 'Fields' array.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utilty - crmstyle_coarse_grid
!
! Code Description:
!   Language:           Fortran 90
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments:

INTEGER, INTENT(IN) ::     &
  num_ff                   & ! Number of fieldsfiles
 ,num_want                 & ! required number of different stash fields
 ,date1(6)                 & ! date time for prognostics (time T)
 ,date2(6)                 & ! date time for prognostics (time T+1 step)
 ,stash_list(num_want)     & ! required stash_lists
 ,lev_list(num_want)       & ! Number of levels required
 ,date_typ(num_want)         ! use date 1 and date 2 to check

TYPE(UM_Header_type), INTENT(IN) :: ff_hdr(num_ff)  ! UM Headers: fieldsfile


INTEGER, INTENT(OUT) :: &
  fstart_num(num_want)  & ! start position of required field
 ,fend_num(num_want)    & ! end position of required fields
 ,file_num(num_want)      ! file number

!-----------------------------------------------------------------------
! Local Variables:
INTEGER :: i, j, k                      ! local counter
INTEGER :: count_lev
INTEGER :: num_get

CHARACTER(LEN=*), PARAMETER :: RoutineName = "LOCATE_FIELDS"

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------

DO i = 1,num_want

  fstart_num(i) = 0
  fend_num(i)   = 0
  file_num(i)   = 0

END DO

! Format ok provided num_want remains below 40

WRITE(umMessage,'(A)') ' Routine locate fields - stashlist'
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(40I6)') (stash_list(i),i=1,num_want)
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A)') ' Routine locate fields - lev_list'
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(40I6)') (lev_list(i),i=1,num_want)
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(A)') ' Routine locate fields - date_typ'
CALL umPrint(umMessage,src=RoutineName)
WRITE(umMessage,'(40I6)') (date_typ(i),i=1,num_want)
CALL umPrint(umMessage,src=RoutineName)

IF (l_sect30) THEN
  num_get = num_want
ELSE
  num_get = num_want_no30
END IF

!-------------------------------
! Loop through required fields
!-------------------------------
DO i = 1,num_get
  count_lev = 0
  IF (date_typ(i) == 1) THEN

    DO k=1, num_ff
      DO j=1, ff_hdr(k)%numflds

        IF (stash_list(i) == ff_hdr(k)%lookup(j)%STCode .AND.                &
           fstart_num(i) == 0 .AND.                                          &
           date1(6) == ff_hdr(k)%lookup(j)%validsec .AND.                    &
           date1(5) == ff_hdr(k)%lookup(j)%validmin .AND.                    &
           date1(4) == ff_hdr(k)%lookup(j)%validhour  .AND.                  &
           date1(3) == ff_hdr(k)%lookup(j)%validdate .AND.                   &
           date1(2) == ff_hdr(k)%lookup(j)%validmonth .AND.                  &
           date1(1) == ff_hdr(k)%lookup(j)%validyear ) THEN
          fstart_num(i) = j
          file_num(i)  = k
          count_lev = 1
          IF (lev_list(i) == 1) THEN
            fend_num(i) = fstart_num(i)

          END IF
        ELSE IF (stash_list(i) == ff_hdr(k)%lookup(j)%STCode                 &
                        .AND. count_lev >= 1  .AND.                          &
           date1(6) == ff_hdr(k)%lookup(j)%validsec .AND.                    &
           date1(5) == ff_hdr(k)%lookup(j)%validmin .AND.                    &
           date1(4) == ff_hdr(k)%lookup(j)%validhour  .AND.                  &
           date1(3) == ff_hdr(k)%lookup(j)%validdate .AND.                   &
           date1(2) == ff_hdr(k)%lookup(j)%validmonth .AND.                  &
           date1(1) == ff_hdr(k)%lookup(j)%validyear ) THEN
          count_lev = count_lev+1
          IF (lev_list(i) == count_lev) THEN
            fend_num(i) = j
          END IF
        END IF
      END DO
    END DO

  ELSE    ! date typ 2

    DO k=1, num_ff
      DO j=1, ff_hdr(k)%numflds
        IF (stash_list(i) == ff_hdr(k)%lookup(j)%STCode .AND.                &
           fstart_num(i) == 0 .AND.                                          &
           date2(6) == ff_hdr(k)%lookup(j)%validsec .AND.                    &
           date2(5) == ff_hdr(k)%lookup(j)%validmin .AND.                    &
           date2(4) == ff_hdr(k)%lookup(j)%validhour  .AND.                  &
           date2(3) == ff_hdr(k)%lookup(j)%validdate .AND.                   &
           date2(2) == ff_hdr(k)%lookup(j)%validmonth .AND.                  &
           date2(1) == ff_hdr(k)%lookup(j)%validyear ) THEN
          fstart_num(i) = j
          file_num(i)  = k
          count_lev = 1
          IF (lev_list(i) == 1) THEN
            fend_num(i) = fstart_num(i)
          END IF
        ELSE IF (stash_list(i) == ff_hdr(k)%lookup(j)%STCode                 &
                        .AND. count_lev >= 1 .AND.                           &
           date2(6) == ff_hdr(k)%lookup(j)%validsec .AND.                    &
           date2(5) == ff_hdr(k)%lookup(j)%validmin .AND.                    &
           date2(4) == ff_hdr(k)%lookup(j)%validhour  .AND.                  &
           date2(3) == ff_hdr(k)%lookup(j)%validdate .AND.                   &
           date2(2) == ff_hdr(k)%lookup(j)%validmonth .AND.                  &
           date2(1) == ff_hdr(k)%lookup(j)%validyear ) THEN
          count_lev = count_lev+1
          IF (lev_list(i) == count_lev) THEN
            fend_num(i) = j
          END IF
        END IF
      END DO
    END DO

  END IF   ! end date_typ
END DO

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE locate_fields

END MODULE locate_fields_mod
