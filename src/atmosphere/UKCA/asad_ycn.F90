! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Purpose: This routine is intended for use with stiff integrators which
!     require the evaluation of df/dt at any point in time.
!
!     This routine is intended to be passed as an argument to the stiff
!     integrators such as the NAG and SVODE drivers.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_DIFFUN
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!     Interface
!     ---------
!        tz     - on entry, specifies the time (unchanged).
!        f1     - on entry, contains species at time TZ (unchanged).
!        dfdt   - on exit, must contain time derivative of the gridpt.
!
!     As all the chemistry arrays are in modules, we have to copy the
!     input to the asad_mod arrays and copy the tendency to the output
!     argument at the end.
!
!     Method
!     ------
!     IMPORTANT!!! This subroutine assumes that only a single gridpt
!     is being worked on by the integrator calling this subroutine.
!
!     The first step is to copy the passed values of the
!     species back into ASAD_MOD so that the production
!     and loss can be computed. Since the integrator is only integrating
!     the variables, f, passed between the model and ASAD we only copy t
!     the species that will change during the timestep. Note that we onl
!     copy to the first element in the species array, y.
!
!     When no families are in use, we can use the index array nlf since
!     this stores a list of all the species of type TR and nf = jpctr.
!
!     The next step is to compute the rates of change. The routine diffu
!     is used but we tell it only to compute the first gridpt. On exit
!     from diffun, the fdot array will have been assigned.
!
!     Externals
!     ---------
!       diffun    - to compute rates of change.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
MODULE asad_ycn_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'ASAD_YCN_MOD'

CONTAINS

SUBROUTINE asad_ycn(tz,f1,dfdt)

USE asad_mod,           ONLY: nf, nlf, y, fdot
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE errormessagelength_mod, ONLY: errormessagelength
USE asad_diffun_mod, ONLY: asad_diffun
IMPLICIT NONE



REAL, INTENT(IN)  :: tz
REAL, INTENT(IN)  :: f1(jpctr)

REAL, INTENT(OUT) :: dfdt(jpctr)

!       Local variables

INTEGER :: j    ! Loop variable
INTEGER :: js

LOGICAL, SAVE :: gfirst = .TRUE.

CHARACTER(LEN=errormessagelength) :: cmessage  ! Error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_YCN'


!       1. Initialise species array; ONLY COPY TO FIRST ELEMENT!
!          ---------- ------- ------ ---- ---- -- ----- --------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF ( gfirst ) THEN
  IF ( method  >=  10 .AND. jpctr  /=  nf ) THEN
    WRITE(umMessage,*) '** INTERNAL ASAD ERROR: jpctr  /=  nf in ycn',&
    ' There should not be any families in use with the stiff', &
    ' integrators.'
    CALL umPrint(umMessage,src='asad_ycn')
    cmessage = 'jpctr  /=  nf in ycn'

    CALL ereport('ASAD_YCN',jpctr,cmessage)
  END IF
  gfirst = .FALSE.
END IF

!       n.b. nf = jpctr when no families are in use

DO j = 1, nf
  js  = nlf(j)
  y(1,js) = f1(j)
END DO

!       2.   Compute rates of change.
!            ------- ----- -- -------

CALL asad_diffun( 1 )

DO j = 1, jpctr
  dfdt(j) = fdot(1,j)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_ycn
!--------------------------------------------------------------------
END MODULE asad_ycn_mod
