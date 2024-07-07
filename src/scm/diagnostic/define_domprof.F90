! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Set up a domain profile for SCM diagnostics

SUBROUTINE define_domprof                                                     &
  ( tindex, tname, rowa1, rowa2, rowb1, rowb2, lev1, lev2, SCMop )


USE UM_types
USE s_scmop_mod, ONLY: SCMop_type, maxndomprof
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY: umPrint, umMessage, newline

IMPLICIT NONE


! Description:
! Set up a domain profile, i.e. a spatial region over which
! diagnostics may be defined. This will be stored in SCMop
! and in future can be referred to using tindex.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model

! Code Description:
! Language: Fortran 90

INTEGER ::           &
  tindex              ! In An index for the profile

CHARACTER (LEN=*) :: &
  tname               ! In A name for the profile

INTEGER ::           &
  rowa1              &! In The horizontal area over which
, rowa2              &! diagnostics of this type are defined
, rowb1              &
, rowb2

INTEGER ::           &
  lev1               &! In The vertical range over which " " "
, lev2

TYPE(SCMop_type) :: SCMop ! InOut The derived-type structure
                          !       containing all the diagnostic
                          !       information

! Dr Hook
!==============================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DEFINE_DOMPROF'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Check that tindex is not greater than maxndomprof (i.e.
! avoid array out-of-bounds errors)
IF (tindex >  maxndomprof) THEN
  WRITE(umMessage,'(A)')                                                       &
    '========================================================='//     newline//&
    '| Define_Domprof ERROR:'                                  //     newline//&
    '| Trying to define a domain'                              //     newline//&
    '| profile with an index above the maximum. Increase'      //     newline//&
    '| maxndomprof or use a lower index. Domain profile'       //     newline//&
    '| not defined: ' // TRIM(ADJUSTL(tname))                  //     newline//&
    '========================================================='
  CALL umPrint(umMessage,src='define_domprof')

  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
  RETURN

END IF

! For safety, enforce a rule that no diagnostic may be defined
! unless *every* domain has been defined. This is to ensure
! SCMop%nelements(x) can never become incorrect by a domain being
! re-defined.
IF (SCMop%nentries >  0) THEN
  WRITE(umMessage,'(A,I0)')                                                    &
    'Define_Domprof ERROR: nentries is non-zero, setting'//           newline//&
    'to zero now. Some diagnostics may be discarded ',                         &
    SCMop%nentries
  CALL umPrint(umMessage,src='define_domprof')
  SCMop%nentries = 0
END IF

! Fill the relevant section of the SCMop structure
SCMop%d_name(tindex)  = tname
SCMop%d_rowa1(tindex) = rowa1
SCMop%d_rowa2(tindex) = rowa2
SCMop%d_rowb1(tindex) = rowb1
SCMop%d_rowb2(tindex) = rowb2
SCMop%d_lev1(tindex)  = lev1
SCMop%d_lev2(tindex)  = lev2

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE define_domprof

