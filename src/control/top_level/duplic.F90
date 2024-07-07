! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Deletes duplicate diags & times; checks for overlap levs & times.
!
! Subroutine Interface:

SUBROUTINE duplic(nrecs,ntimes,nlevels)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE version_mod, ONLY:                                            &
    nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,           &
    nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
    nlevp, npslevp, npslistp
USE submodel_mod, ONLY: n_internal_model_max
USE stparam_mod, ONLY: st_freq_code, st_output_bottom, st_macrotag,&
    st_model_code, st_proc_no_code, st_start_time_code,            &
    st_output_code, st_output_top, st_end_time_code, st_input_code,&
    st_period_code, st_gridpoint_code, st_weight_code,             &
    st_north_code, st_south_code, st_west_code, st_east_code,      &
    st_series_ptr, st_pseudo_out, st_output_type
USE stextend_mod, ONLY: llistty, indx_s, list_s, itim_s,  &
                        rlevlst_s, levlst_s

IMPLICIT NONE
!
! Description:
!   Deletes duplicate diagnostic entries from STASH list; deletes
!   duplicate STASH times;
!   checks for overlap of levels and times, to some extent.
!   Called by STPROC.
!   Input : NRECS   No. of STASH list records
!           NTIMES  No. of STASH times
!           NLEVELS No. of STASH levels
!           LIST_S  STASH list array with prelim. pointer system
!           ITIM_S  STASH times array
!   Output: NRECS   Reduced no. of STASH list records
!           NTIMES  Reduced no. of STASH times
!           NLEVEL  Reduced no. of STASH levels
!           ITIM_S  Reduced STASH times array
!           LIST_S  Reduced STASH list with prelim. pointers,
!                           consistent with STASH times.
!
! Method:
!
!   (a) STASH times tables in ITIM_S.
!       The times at which STASH processing occurs for a diagnostic
!   IREC may be specified by the entries (start_time,end_time,period)
!   in LIST_S.
!       Alternatively, if LIST_S(st_freq_code,IREC) has value '-n',
!   then STASH processing times for this diagnostic are given by a
!   'times table' in ITIM_S.
!   (In such a case, the above 3 entries in LIST_S are ignored).
!   The times table is given by column 'n' of ITIM_S, i.e.,
!   ITIM_S(time,n).
!   In this routine, the logical array entry LTUSE(n) is set to
!   .TRUE. if col 'n' of ITIM_S contains a times table. Any column
!   of ITIM_S which does not contain a times table is filled with
!   -1's. The cols which contain times tables are then shuffled along,
!   so that they occupy the first NTIMES cols of ITIM_S. The pointers
!   in LIST_S(st_freq_code,IREC) are altered accordingly.
!
!   (b) STASH levels lists in (R)LEVLST_S.
!       The levels on which STASH processing occurs for a diagnostic
!   IREC is specified by the entries (output_bottom, output_top) in
!   LIST_S.
!     If LIST_S(bot)=m, then output is on a range of model levels,
!   with level m as the bottom level, and LIST_S(top) points to the
!   top output model level.
!     If LIST_S(bot)=-n, then there is a levels list in col 'n' of
!   LEVLST_S, and LIST_S(top) contains a code value indicating the
!   type of levels (model, pressures, heights or theta). Each levels
!   list also has a corresponding entry in LLISTTY, indicating whether
!   the list is real or integer.
!     In this routine, the cols of LEVLST_S which contain levels lists
!   are shuffled along so that they occupy the first NLEVELS cols of
!   LEVLST_S. The pointers in LIST_S(output_bottom,IREC), and the
!   entries in LLISTTY, are altered accordingly.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
! Subroutine Arguments:
!
!   Scalar arguments with intent(InOut):

INTEGER :: nrecs      ! No. of STASH list records
INTEGER :: ntimes     ! No. of STASH times
INTEGER :: nlevels    ! No. of STASH levels

! Local scalars:

LOGICAL :: ltrpt
LOGICAL :: testflag
INTEGER :: i
INTEGER :: i1
INTEGER :: i2
INTEGER :: iend
INTEGER :: iitm
INTEGER :: il
INTEGER :: irec
INTEGER :: isec
INTEGER :: istr
INTEGER :: it
INTEGER :: itags1
INTEGER :: itags2
INTEGER :: itagu1
INTEGER :: itagu2
INTEGER :: itend1
INTEGER :: itend2
INTEGER :: modl
INTEGER :: nlevsw
INTEGER :: nrecsw
INTEGER :: ntimesw

! Local arrays:

LOGICAL :: ltuse(2*nproftp+2)  ! LTUSE(n) set to .T. if column n
                            ! in ITIM_S contains a STASH times
                            ! table.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DUPLIC'

!- End of Header -----------------------------------------------------


! Initialise LTUSE array

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
DO i=1,ntimes
  ltuse(i)=.FALSE.
END DO


! Blank out unused STASH times

DO irec=1,nrecs
  IF (list_s(st_freq_code,irec) <  0) THEN   ! STASH times table
    ltuse(-list_s(st_freq_code,irec))=.TRUE. !  exists for IREC
  END IF
END DO

DO i=1,ntimes
  IF (.NOT. ltuse(i)) itim_s(1,i)=-1  ! Fill unused columns in
END DO                              ! ITIM_S with -1 in each row.


! Delete blank STASH times

ntimesw=1

DO it=1,ntimes

  ! If col 'IT' contains a times table, find
  ! corresponding record IREC in LIST_S, and replace entry '-IT'
  ! by '-NTIMESW'. In each case, NTIMESW <= IT.

  IF (itim_s(1,it) /= -1) THEN
    DO irec=1,nrecs

      IF (list_s(st_freq_code,irec) == -it) THEN
        list_s(st_freq_code,irec)=-ntimesw
      END IF

    END DO

    IF (it /= ntimesw) THEN
      ! Move times table in col 'IT' to col 'NTIMESW'. Hence array
      ! ITIM_S is compressed.
      DO i=1,ntimep
        itim_s(i,ntimesw)=itim_s(i,it)
      END DO
    END IF

    ntimesw=ntimesw+1
  END IF
END DO

ntimes=ntimesw-1  ! No. of STASH-times tables remaining, so far


! Delete blank STASH levels

nlevsw=1

DO il=1,nlevels

  ! If col 'IL' of LEVLST_S contains a levs list, then find corresponding
  !  record IREC in LIST_S and replace entry '-IL' by '-NLEVSW'. In each
  ! case, NLEVSW <= IL.

  IF (levlst_s(1,il) /= 0) THEN
    DO irec=1,nrecs
      IF (list_s(st_output_bottom,irec) == -il) THEN
        list_s(st_output_bottom,irec)=-nlevsw
      END IF
    END DO
    IF (il /= nlevsw) THEN
      ! Move levels list in col 'IL' to col 'NLEVSW'. Hence array
      ! LEVLST_S is compressed.
      DO i=1,nlevp_s
        levlst_s(i,nlevsw)=levlst_s(i,il)
        rlevlst_s(i,nlevsw)=rlevlst_s(i,il)
      END DO
      llistty(nlevsw)=llistty(il) ! Move corresponding entry in
    END IF                        ! LLISTTY

    nlevsw=nlevsw+1
  END IF
END DO

nlevels=nlevsw-1


! Check for duplication/overlap of STASH levels

nrecsw=nrecs

DO modl  = 1,n_internal_model_max
  DO isec  = 0,nsectp
    DO iitm  = 1,nitemp

      IF (indx_s(2,modl,isec,iitm) >= 2) THEN  !More than one STASH rec
                                              !  for (model,sec,item)
        istr=     indx_s(1,modl,isec,iitm)    !1st record with m,s,i
        iend=istr+indx_s(2,modl,isec,iitm)-1  !Last record with m,s,i

        DO i1=istr,iend-1

          itags1=list_s(st_macrotag,i1)/1000              ! System tag
          itagu1=list_s(st_macrotag,i1)-1000*itags1       ! User tag

          IF (list_s(st_model_code,i1) <= n_internal_model_max) THEN
            !              Not flagged redundant
            DO i2=i1+1,iend

              itags2=list_s(st_macrotag,i2)/1000           ! System tag
              itagu2=list_s(st_macrotag,i2)-1000*itags2    ! User tag

              IF ((list_s(st_proc_no_code,i1) ==                         &
                  list_s(st_proc_no_code,i2)) .AND.                      &
                 (list_s(st_freq_code,i1) ==                            &
                  list_s(st_freq_code,i2)) .AND.                         &
                 (list_s(st_period_code,i1) ==                          &
                  list_s(st_period_code,i2)) .AND.                       &
                 (list_s(st_gridpoint_code,i1) ==                       &
                  list_s(st_gridpoint_code,i2)) .AND.                    &
                 (list_s(st_weight_code,i1) ==                          &
                  list_s(st_weight_code,i2)) .AND.                       &
                 (list_s(st_north_code,i1) ==                           &
                  list_s(st_north_code,i2)) .AND.                        &
                 (list_s(st_south_code,i1) ==                           &
                  list_s(st_south_code,i2)) .AND.                        &
                 (list_s(st_west_code,i1) ==                            &
                  list_s(st_west_code,i2)) .AND.                         &
                 (list_s(st_east_code,i1) ==                            &
                  list_s(st_east_code,i2)) .AND.                         &
                 (list_s(st_input_code,i1) ==                           &
                  list_s(st_input_code,i2)) .AND.                        &
                 (list_s(st_output_code,i1) ==                          &
                  list_s(st_output_code,i2)) .AND.                       &
                 (list_s(st_series_ptr,i1) ==                           &
                  list_s(st_series_ptr,i2)) .AND.                        &
                 (list_s(st_pseudo_out,i1) ==                           &
                  list_s(st_pseudo_out,i2)) .AND.                        &
                 (list_s(st_output_type,i1) ==                          &
                  list_s(st_output_type,i2) .OR.                         &
                  list_s(st_output_type,i1) == 0 .OR.                   &
                  list_s(st_output_type,i2) == 0) .AND.                  &
              ((itags1 == itags2) .OR. (itags1 == 0) .OR.                  &
              (itags2 == 0)) .AND.                                       &
              ((itagu1 == itagu2) .OR. (itagu1 == 0) .OR.                  &
              (itagu2 == 0)) .AND.                                       &
            (list_s(st_model_code,i2) <= n_internal_model_max)) THEN
                !            Not flagged redundant

                ! If they are the same in all but time and level

                itend1=list_s(st_end_time_code,i1)
                itend2=list_s(st_end_time_code,i2)

                IF (itend1 == -1) itend1=                                &
                   list_s(st_start_time_code,i2)+1     ! Force overlap

                IF (itend2 == -1) itend2=                                &
                   list_s(st_start_time_code,i1)+1     ! Force overlap

                ! Where period_code is zero we have to prevent second part from
                ! being evaluated so break this OR statement into two parts:
                testflag=.FALSE.
                IF ((list_s(st_period_code,i1) == 0) .OR.                 &
                  (list_s(st_period_code,i1) == -1)) THEN
                  testflag=.TRUE.
                ELSE IF ((MOD(list_s(st_start_time_code,i2)-              &
                     list_s(st_start_time_code,i1),                     &
                     list_s(st_period_code,i1)) == 0)) THEN
                  testflag=.TRUE.
                END IF

                IF ((.NOT. ((list_s(st_start_time_code,i1)                &
                                                  >  itend2) .OR.        &
                  (itend1 <  list_s(st_start_time_code,i2))) .OR.        &
                    (list_s(st_output_code,i1) >  0)) .AND.              &
                (MOD(list_s(st_start_time_code,i2)-                     &
                     list_s(st_start_time_code,i1),                     &
                     list_s(st_freq_code,i1)) == 0) .AND.                &
                     testflag .AND.                                      &
                    (list_s(st_output_bottom,i2) ==                     &
                     list_s(st_output_bottom,i1)) .AND.                  &
                    (list_s(st_output_top,i2) ==                        &
                     list_s(st_output_top,i1))) THEN

                  ! (Times overlap or in dump) and overlay in freq & period
                  !                            and levels the same

                  IF (itagu1 == 0) THEN
                    list_s(st_macrotag,i1)=itagu2
                    itagu1=itagu2
                  ELSE
                    list_s(st_macrotag,i1)=itagu1
                  END IF

                  IF (itags1 == 0) THEN
                    list_s(st_macrotag,i1)=                             &
                    itags2*1000+list_s(st_macrotag,i1)
                    itags1=itags2
                  ELSE
                    list_s(st_macrotag,i1)=                             &
                    itags1*1000+list_s(st_macrotag,i1)
                  END IF

                  list_s(st_start_time_code,i1)=                        &
                  MIN(list_s(st_start_time_code,i1),                    &
                      list_s(st_start_time_code,i2))

                  IF ((list_s(st_end_time_code,i1) == -1) .OR.            &
                     (list_s(st_end_time_code,i2) == -1)) THEN
                    list_s(st_end_time_code,i1)=-1
                  ELSE
                    list_s(st_end_time_code,i1)=                        &
                    MAX(list_s(st_end_time_code,i1),                    &
                        list_s(st_end_time_code,i2))
                  END IF

                  list_s(st_model_code,i2)=n_internal_model_max+1
                  ! Sets model id to be greater than no of models,
                  ! so that this diag is put at the end of any sorted list.

                  nrecsw=nrecsw-1

                  DO i=istr,iend                      !Change pointers
                    IF (list_s(st_input_code,i )   ==                    &
                      -list_s(nelemp+1     ,i2)) THEN
                      list_s(st_input_code,i ) =                       &
                     -list_s(nelemp+1     ,i1)
                    END IF
                  END DO
                END IF

              END IF   ! I1,I2 comparison
            END DO     ! I2
          END IF      ! I1 Not flagged redundant
        END DO       ! I1

      END IF   ! More than one STASH record for m,s,i
    END DO     ! Items
  END DO     ! Sections
END DO     ! Models

! Remove unwanted records (i.e., those flagged redundant)

! DEPENDS ON: order
CALL order(nrecs)
nrecs=nrecsw
! DEPENDS ON: sindx
CALL sindx(nrecs)
!
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE duplic
