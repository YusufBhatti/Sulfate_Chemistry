! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINES VANMOPS_MIXED_PHASE and CC_TO_RHTOT-------------------
!
!    Purpose : Performs vertical analysis of MOPS cloud data
!              for (2A) mixed phase cloud microphysics scheme
!       When IPASS=1,
!       model field is interpolated to ob locations and increments
!       calculated. Preliminary weight normalisation factors are also
!       calculated.
!
!       When IPASS=2,
!       data density is interpolated to ob locations and
!       final weight normalisation factors are calculated.
!
!       Documentation in a Working Paper by B Macpherson & D Wilson is
!       on Metweb at URL
!       http://fr2010/~frff/papers/MOPS_for_new_micro.ps.gz
!
!    For use on Cray
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!
!
!
!    Arguments:---------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation


MODULE cc_to_rhtot_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CC_TO_RHTOT_MOD'

CONTAINS

SUBROUTINE cc_to_rhtot  (cc,npts,rhc,rhtot,k,bl_levels)

USE cloud_inputs_mod, ONLY: i_eacf, all_clouds, not_mixph
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER ::                                                        &
 npts                                                             &
                     ! IN No. of points on level.
,k                                                                &
                     ! IN: Level no.
,bl_levels           ! IN: Number of BL levels

REAL ::                                                           &
 rhc                                                              &
                     ! IN Critical relative humidity (fraction).
,cc(npts)                                                         &
                     ! IN Cloud cover (fraction).
,rhtot(npts)                                                      &
                     ! OUT Total relative humidity (fraction).
,tempnum

!-------------------------------------------------------------------
!     External subroutine calls NONE
!-------------------------------------------------------------------

! Local variables---------------------------------------------------

INTEGER :: i     ! Do loop index

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CC_TO_RHTOT'


! Code in calling routine restricts CC to range 0-1
! but check to be safe

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO i=1,npts
  IF (cc(i) >  1.0) cc(i) = 1.0
  IF (cc(i) <  0.0) cc(i) = 0.0

  IF (i_eacf == all_clouds .OR. i_eacf == not_mixph) THEN
    ! Use empirically adjusted cloud fraction

    IF (k <= bl_levels) THEN
      IF (cc(i) <= 0.5) THEN
        tempnum = (-1.0+SQRT(2.0*cc(i)))*0.816 - 0.184
        rhtot(i) = tempnum*(1.0-rhc) + 1.0
      ELSE
        tempnum = (1.0-SQRT(2.0*(1.0-cc(i))))*0.816 - 0.184
        rhtot(i) = tempnum*(1.0-rhc) + 1.0
      END IF
    ELSE
      IF (cc(i) <= 0.5) THEN
        tempnum = (-1.0+SQRT(2.0*cc(i)))*0.9045 - 0.0955
        rhtot(i) = tempnum*(1.0-rhc) + 1.0
      ELSE
        tempnum = (1.0-SQRT(2.0*(1.0-cc(i))))*0.9045 - 0.0955
        rhtot(i) = tempnum*(1.0-rhc) + 1.0
      END IF
    END IF

  ELSE  ! i_eacf
    ! Use normal method to calculate RHtot

    IF (cc(i) >  0.5) THEN
      !  working paper eqn 9b
      rhtot(i) = 1.0 + (1.0-rhc) * (1.0 - SQRT(2.0*(1.0-cc(i))  )  )
    ELSE
      !  working paper eqn 9a
      rhtot(i) = 1.0 + (1.0-rhc) * ( SQRT(2.0*cc(i))-1.0 )
    END IF

  END IF  ! i_eacf
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE cc_to_rhtot

END MODULE cc_to_rhtot_mod
