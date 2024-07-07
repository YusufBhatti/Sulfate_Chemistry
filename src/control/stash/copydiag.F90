! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Purpose : To copy a single diagnostic field from secondary space to
!             the main data array for stash processing
!             Option to extend data to a full horizontal field.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH

SUBROUTINE copydiag(                                              &
     diagout,diagin,                                              &
     row_length,rows,                                             &
     offx_out, offy_out,                                          &
     offx_in, offy_in,                                            &
     at_extremity,                                                &
     im,IS,ie,                                                    &
     icode,cmessage)

USE um_parparams, ONLY: pnorth, peast, psouth, pwest
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_global

IMPLICIT NONE

INTEGER ::                                                        &
     row_length                                                   &
                          ! Number of points in a row
,    rows                                                         &
                          ! Number of rows
,    offx_out, offy_out                                           &
                          ! offset dimensions of diagout
,    offx_in, offy_in                                             &
                          ! offset dimensions of diagin
,    im,IS,ie                                                     &
                          ! Model, section, item
,    icode                ! Return code  =0 Normal exit  >1 Error

LOGICAL ::                                                        &
     at_extremity(4)      ! Indicates if this processor is at
                          !  north, south east or west of the
                          !  processor grid
CHARACTER(LEN=errormessagelength) :: cmessage

! ARGPPX arguments:


REAL ::                                                           &
 diagin(1-offx_in:row_length+offx_in, 1-offy_in:rows+offy_in)     &
                          ! Output field
,diagout(1-offx_out:row_length+offx_out,1-offy_out:rows+offy_out)
                          ! Input field

!     Local variables

INTEGER ::                                                        &
   i,j  ! loop bound

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='COPYDIAG'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
icode = 0
cmessage = ""

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(i,j)                                                            &
!$OMP SHARED(rows, row_length, diagout, diagin, offx_out, offy_out,           &
!$OMP        model_type,at_extremity)

!     Copy fields
!$OMP DO SCHEDULE(STATIC)
DO j = 1, rows
  DO i = 1, row_length
    diagout(i,j) = diagin(i,j)
  END DO
END DO
!$OMP END DO

IF (model_type /= mt_global) THEN

  ! Marker only for possible extensions of code. Since diagnostics are
  ! so far defined as having no halos, extra code to populate halos is
  ! not used at present, so bypassed here on in normal case, when
  ! offx_out=offy_out=0 .

  IF (offx_out /= 0 .OR. offy_out /= 0) THEN ! check no halos

    IF (at_extremity(PNorth)) THEN
!$OMP DO SCHEDULE(STATIC)
      DO i = 1-offx_out, row_length+offx_out
        DO j = rows, rows+offy_out
          diagout(i,j) = diagout(i,rows)
        END DO
      END DO
!$OMP END DO
    END IF

    IF (at_extremity(PSouth)) THEN
!$OMP DO SCHEDULE(STATIC)
      DO i = 1-offx_out, row_length+offx_out
        DO j = 1-offy_out, 1
          diagout(i,j) = diagout(i,1)
        END DO
      END DO
!$OMP END DO
    END IF

    ! copy diagnostic information to e and w boundaries

    IF (at_extremity(PEast)) THEN
!$OMP DO SCHEDULE(STATIC)
      DO i = 1-offx_out, 1
        DO j = 1-offy_out, rows+offy_out
          diagout(i,j) = diagout(1,j)
        END DO
      END DO
!$OMP END DO
    END IF


    IF (at_extremity(PWest)) THEN
!$OMP DO SCHEDULE(STATIC)
      DO i = row_length, row_length+offx_out
        DO j = 1-offy_out, rows+offy_out
          diagout(i,j) = diagout(row_length,j)
        END DO
      END DO
!$OMP END DO
    END IF

  END IF      ! check no halos
END IF     ! .not. GLOBAL
!$OMP END PARALLEL

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

END SUBROUTINE copydiag
