! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Change pointer system for "child" diags to final STASH list version
!
! Subroutine Interface:

SUBROUTINE pointr(nrecs)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE version_mod, ONLY:                                            &
    nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,           &
    nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
    nlevp, npslevp, npslistp
USE submodel_mod, ONLY: n_internal_model_max
USE stextend_mod, ONLY: indx_s, list_s
USE stparam_mod, ONLY: st_input_code

IMPLICIT NONE

! Description:
!   The stash list with preliminary pointer system, and the stash
!   index, are input to this routine. The stash list with final
!   pointer system is output.
!   Called by STPROC.
!
!   Fuller explanation:
!   Any diag in the stash list which has a processing code in the range
!   2-7 (i.e., accumulate, time mean, append time series, max, min,
!   trajectory) has one or more "child records". A child is another
!   diag, with the same m,s,i. The output from the parent diag is used
!   as input to the child diag, which is then processed to produce
!   further output. Each child record has an entry which points to its
!   parent record - the st_input_code entry. In routine PRELIM, a
!   preliminary pointer system is set up, involving the use
!   of the "extra entry", NELEMP+1. In each record, entry NELEMP+1 is
!   set to the current value of NRECS, i.e., the position of that
!   record in the prelim stash list. The value of the st_input_code
!   entry for a child record is set to the negative of the NRECS value
!   for its parent. Note that, in the prelim stash list, the children
!   of a particular parent appear immediately after the parent.
!   So, after PRELIM, each record in the stash list identifies
!   itself by its NELEMP+1 entry, and each child record identifies
!   its parent by its st_input_code entry. The final position of
!   each record in the stash list is given by the INDX_S array.
!   This subroutine therefore changes the  st_input_code entry of
!   each child record so that it agrees with INDX_S.
!   The NELEMP+1 entry is then no longer relevant.
!
! Method:
!   Uses INDX_S array to identify parent records (i.e., diagnostics
!   which have more than one entry in the stash list).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
! Subroutine arguments:

!   Scalar arguments with intent(in):

INTEGER :: nrecs    ! No. of records in stash list

! Local scalars

INTEGER :: modl  ! Loop counter for internal models
INTEGER :: isec  ! Do. sections
INTEGER :: iitm  ! Do. items
INTEGER :: istr  ! Position of parent record in stash list
INTEGER :: iend  ! Position of final child record in stash list
INTEGER :: i1
INTEGER :: i2
INTEGER :: i3

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='POINTR'

!- End of Header ----------------------------------------------------

! Loop over models, section, items

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
DO modl=1,n_internal_model_max
  DO isec=0,nsectp
    DO iitm=1,nitemp

      ! Examine INDX_S entry to find out whether there are child record(s)

      IF (indx_s(2,modl,isec,iitm) >= 2) THEN

        istr=     indx_s(1,modl,isec,iitm)
        iend=istr+indx_s(2,modl,isec,iitm)-1

        DO i1=istr,iend-1
          DO i2=i1+1,iend
            IF (list_s(st_input_code,i2) ==                            &
              -list_s(nelemp+1     ,i1)) THEN    
              list_s(st_input_code,i2)=-i1-nrecs
            END IF
          END DO
        END DO

        DO i3=istr,iend
          IF (list_s(st_input_code,i3) <  0) THEN
            list_s(st_input_code,i3)=                                &
            list_s(st_input_code,i3)+nrecs
          END IF
        END DO

      END IF

    END DO  ! Items
  END DO  ! Sections
END DO  ! Internal models

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pointr
