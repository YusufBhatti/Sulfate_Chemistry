! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine COPYDIAG_03D -------------------------------------------
!
!   Purpose : To copy a diagnostic field from secondary space to the
!             main data array for stash processing, and to extend the
!             data to a full horizontal field. Input data of multilevel
!             fields is assumed to be on all model levels. Output data
!             is on the levels required.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH

SUBROUTINE copydiag_03d(                                          &
     diagout,diagin,                                              &
     row_length,rows,levels,                                      &
     offx_out, offy_out,                                          &
     offx_in, offy_in,                                            &
     at_extremity,                                                &
     stlist,len_stlist,stash_levels,                              &
     len_stashlevels,                                             &
     im,IS,ie,                                                    &
     icode,cmessage)

USE um_parparams, ONLY: pnorth, peast, psouth, pwest
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE

INTEGER ::                                                        &
       row_length,                                                &
                    ! Number of points in a row
       rows,                                                      &
                    ! Number of rows
       offx_out, offy_out,                                        &
                           ! offset dimensions of diagout
       offx_in, offy_in,                                          &
                           ! offset dimensions of diagin
       levels,                                                    &
                    ! Number of levels in input data
       len_stlist,                                                &
                    !
       stlist(len_stlist),                                        &
                           ! Stash list
       len_stashlevels,                                           &
                        !
       stash_levels(len_stashlevels,*),                           &
                                        ! Stash levels list.
       im,IS,ie,                                                  &
                    ! Model, section, item
       icode        ! Return code =0 Normal exit
!                                       >1 Error message

LOGICAL :: at_extremity(4)
!! logicals indicating if a processor is at the edge of the LPG

CHARACTER(LEN=errormessagelength) :: cmessage

REAL ::                                                           &
   diagin(1-offx_in:row_length+offx_in, 1-offy_in:rows+offy_in    &
          ,0:levels)                                              &
                          ! Output data
,  diagout(1-offx_out:row_length+offx_out,1-offy_out:rows+offy_out&
          ,*)             ! Input data

!     Local variables

INTEGER ::                                                        &
   i,j,                                                           &
                         !
   k,                                                             &
                         !
   kout, nout,                                                    &
                         !
   levels_list(levels+1)         ! list of levels to be copied

LOGICAL ::                                                        &
   list(0:levels) !

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPYDIAG_03D'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
icode = 0
cmessage = ""

! DEPENDS ON: set_zero_levels_list
CALL set_zero_levels_list(levels, len_stlist, stlist, list, &
  stash_levels, len_stashlevels, icode, cmessage)
IF (icode > 0) GO TO 9999

!  Move data from DIAGIN to DIAGOUT at levels requested

nout = 0
DO k = 0, levels
  IF (list(k)) THEN
    nout = nout + 1
    levels_list(nout) = k
  END IF
END DO

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( nout, levels_list, rows, row_length, diagout, diagin,    &
!$OMP         model_type, offx_out, offy_out, at_extremity )           &
!$OMP PRIVATE ( i, j, k, kout )
!DIR$ IVDEP
DO kout = 1, nout

  k = levels_list(kout)
  !           Copy fields
  DO j = 1, rows
    DO i = 1, row_length
      diagout(i,j,kout) = diagin(i,j,k)
    END DO
  END DO


  IF (model_type /= mt_global) THEN

    ! Marker only for possible extensions of code. Since diagnostics are
    ! so far defined as having no halos, extra code to populate halos is
    ! not used at present, so bypassed here on in normal case, when
    ! offx_out=offy_out=0 .

    IF (offx_out /= 0 .OR. offy_out /= 0) THEN ! check no halos

      IF (at_extremity(PNorth)) THEN
        DO i = 1-offx_out, row_length+offx_out
          DO j = rows, rows+offy_out
            diagout(i,j,kout) = diagout(i,rows,kout)
          END DO
        END DO
      END IF

      IF (at_extremity(PSouth)) THEN
        DO i = 1-offx_out, row_length+offx_out
          DO j = 1-offy_out, 1
            diagout(i,j,kout) = diagout(i,1,kout)
          END DO
        END DO
      END IF

      ! copy diagnostic information to e and w boundaries

      IF (at_extremity(PEast)) THEN
        DO i = 1-offx_out, 1
          DO j = 1-offy_out, rows+offy_out
            diagout(i,j,kout) = diagout(1,j,kout)
          END DO
        END DO
      END IF


      IF (at_extremity(PWest)) THEN
        DO i = row_length, row_length+offx_out
          DO j = 1-offy_out, rows+offy_out
            diagout(i,j,kout) = diagout(row_length,j,kout)
          END DO
        END DO
      END IF

    END IF      ! check no halos
  END IF     ! .not. GLOBAL

END DO
!$OMP END PARALLEL DO

9999 CONTINUE
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE copydiag_03d
