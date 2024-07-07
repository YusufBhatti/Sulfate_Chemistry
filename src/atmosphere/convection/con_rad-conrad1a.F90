! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Convective Cloud Top, Base and Amount Calculation Scheme.
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection

MODULE con_rad_4a5a_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CON_RAD_4A5A_MOD'
CONTAINS

!-----------------------------------------------------------------
!  SUBROUTINE con_rad_4a5a
!
!  PURPOSE : Calculates convective cloud top, base and amount
!
!  Suitable for single column model use
!
!  PROGRAMMING STANDARDS : Unified Model Documentation Paper NO. 4
!  VERSION NO. 1
!
!  LOGICAL COMPONENT NUMBER: P27
!
!  SYSTEM TASK : P27
!
!  DOCUMENTATION : Unified Model Documentation Paper P27
!                  Section (9)
!
!-----------------------------------------------------------------


SUBROUTINE con_rad_4a5a (                                        &
          k,xpk,xpkp1,flxkp1,bterm,cca,iccb,icct,start_lev,      &
          tcw,ccw,cclwp,delpkp1,lcca,lcbase,lctop,npnts,         &
          l_q_interact)

USE science_fixes_mod,      ONLY: l_fix_ccb_cct
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr,             ONLY: PrintStatus, PrStatus_Diag

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE planet_constants_mod, ONLY: g

IMPLICIT NONE



!
!----------------------------------------------------------------------
! VECTOR LENGTH AND LOOP VARIABLES
!----------------------------------------------------------------------
!
INTEGER :: npnts        ! IN VECTOR LENGTH
!
INTEGER :: k            ! IN PRESENT MODEL LAYER
!
INTEGER :: i            ! LOOP COUNTER
!
! ----------------------------------------------------------------------
! Arguments with Intent IN.  ie:  Input variables.
! ----------------------------------------------------------------------
!
REAL :: xpk(npnts)      ! IN PARCEL CLOUD WATER IN LAYER K (KG/KG)
!
REAL :: xpkp1(npnts)    ! IN PARCEL CLOUD WATER IN LAYER K+1 (KG/KG)
!
LOGICAL :: bterm(npnts) ! IN MASK FOR POINTS WHERE CONVECTION
                     !    IS ENDING
!
REAL :: flxkp1(npnts)   ! IN PARCEL MASSFLUX IN LAYER K+1 (PA/S)
!
REAL :: delpkp1(npnts)  ! IN PRESSURE DIFFERENCE ACROSS LAYER K+1
!
REAL :: ccw(npnts)      ! IN PARCEL CLOUD WATER AS CALCULATED BEFORE
                     !    PRECIPITATION. LAYER K+1 (KG/KG)
!
!
INTEGER :: start_lev(npnts)  ! IN LEVEL AT WHICH CONVECTION INITIATED
!
LOGICAL :: l_q_interact ! IN Switch (PC2) allows overwriting of
                     ! parcel variables
!
! ----------------------------------------------------------------------
! Arguments with Intent IN/OUT. ie: input variables changed on output.
! ----------------------------------------------------------------------
!
REAL :: tcw(npnts)      ! INOUT
                     ! IN  TOTAL CONDENSED WATER SUMMED TO
                     !     LAYER K (KG/M**2/S)
                     ! OUT TOTAL CONDENSED WATER SUMMED TO
                     !     LAYER K+1 OR IF CONVECTION HAS
                     !     TERMINATED ZEROED (KG/M**2/S)
!
REAL :: cclwp(npnts)    ! INOUT
                     ! IN  TOTAL CLOUD LIQUID WATER PATH
                     !     SUMMED TO LAYER K  (KG/M**2)
                     ! OUT TOTAL CLOUD LIQUID WATER PATH
                     !     SUMMED TO LAYER K+1 (KG/M**2)
REAL :: lcca(npnts)      ! INOUT LOWEST CONV.CLOUD AMOUNT (%)
!
INTEGER :: lcbase(npnts) ! INOUT LOWEST CONV.CLOUD BASE LEVEL
!
INTEGER :: lctop(npnts)  ! INOUT LOWEST CONV.CLOUD TOP LEVEL
!
! ----------------------------------------------------------------------
! Arguments with Intent OUT. ie: Output variables.
! ----------------------------------------------------------------------
!
REAL :: cca(npnts)      ! OUT CONVECTIVE CLOUD AMOUNT (%)
!
INTEGER :: iccb(npnts)   ! OUT CONVECTIVE CLOUD BASE LEVEL
!
INTEGER :: icct(npnts)   ! OUT CONVECTIVE CLOUD TOP LEVEL
!

INTEGER :: error  ! Error code for ereport
CHARACTER (LEN=errormessagelength) :: cmessage  ! String to store error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CON_RAD_4A5A'

!
!  External subroutine calls: ------------------------------------------
!
!     EXTERNAL None.
!
!- End of Header
!
! ==Main Block==--------------------------------------------------------
!
! ----------------------------------------------------------------------
!  CALCULATE Cloud Base and Cloud top, and Lowest Base and Top
!
!  WHEN CLOUD BASE SET ZERO TOTAL CONDENSED WATER
! ----------------------------------------------------------------------
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( l_fix_ccb_cct ) THEN
  ! Corrected version of the ccb / cct calculation; consistently set the
  ! base and top to be the first and last level in the ascent where the
  ! parcel contains cloud.

  DO i = 1, npnts
    ! Set cloud-base if cloud just appeared at k+1 (but no cloud at k),
    ! or if we're at the first level of ascent and already cloudy.
    IF ( ( (.NOT.xpk(i) > 0.0) .AND. xpkp1(i) > 0.0 ) .OR.           &
         ( xpk(i) > 0.0 .AND. k  ==  start_lev(i) )       ) THEN
      IF ( xpk(i) > 0.0 ) THEN
        ! Cloudy at start-level: ccb corresponds to the rho-level below k
        iccb(i) = k
      ELSE
        ! Not cloudy at k: ccb corresponds to the rho-level between k and k+1
        iccb(i) = k + 1
      END IF
      cclwp(i)  = 0.0
      tcw(i)   = 0.0
      IF ( lcbase(i) == 0 )  lcbase(i) = iccb(i)
    END IF
    ! Set cloud-top if cloud just vanished at k+1 (but cloud present at k),
    ! or if we're at the last level of ascent and still cloudy.
    IF ( ( xpk(i) > 0.0 .AND. (.NOT.xpkp1(i) > 0.0) ) .OR.           &
         ( xpkp1(i) > 0.0 .AND. bterm(i) )                ) THEN
      ! cct corresponds to the rho-level between k and k+1
      icct(i) = k + 1
      IF ( lctop(i) == 0 )  lctop(i) = icct(i)
    END IF
  END DO

ELSE
  ! Original, confusing code; sets the base to be the first cloudy level
  ! in the ascent, but then the top is only set if still cloudy at
  ! the termination level.  If the parcel forms cloud (so base gets set)
  ! but then the condensate evaporates before reaching the termination
  ! level, the cloud-top level will be left unset, leading to
  ! inconsistent CCrad fields.

  DO i = 1, npnts
  
    IF ( xpk(i)  <=  0.0 .AND. ccw(i)  >   0.0 ) THEN
      !       Assuming initial parcel condensate is zero
      iccb(i)=k+1
      cclwp(i)=0.0
      !
      IF ( lcbase(i)  ==  0 ) THEN
        !         Lowest cloud base
        lcbase(i) = k+1
      END IF  ! If_lcbase_1

      ! Note that, if not l_q_interact, xpk must always be zero
      ! in the 1st level of ascent, so if any condensation occurs,
      ! the above check will set the cloud-base.  However, if using PC2
      ! (l_q_interact), the initial parcel can contain condensate (xpk>0),
      ! so need an additional check...

    ELSE IF ( l_q_interact .AND. xpk(i)  >   0.0 .AND.              &
             k  ==  start_lev(i) ) THEN

      ! ...If the initial parcel contains condensate, then cloud-base level
      ! is the parcel initial level, not k+1.

      ! Non-zero initial parcel condensate (initialized to environment)
      iccb(i)  = k
      cclwp(i) = 0.0
      !
      IF ( lcbase(i)  ==  0 ) THEN
        !         Lowest cloud base
        lcbase(i) = k
      END IF  ! If_lcbase_2

    END IF

    IF (bterm(i) .AND.                                              &
        ((ccw(i) >  0.0) .OR. (xpk(i) >  0.0)) ) icct(i) = k+1

    IF (bterm(i) .AND.  lctop(i) == 0 .AND.                         &
        ((ccw(i) >  0.0) .OR. (xpk(i) >  0.0)) ) THEN
      lctop(i) = k+1
    END IF
  
  END DO  ! i = 1, npnts

END IF  ! l_fix_ccb_cct

! If full diagnostic output requested, or not using the fix...
IF ( PrintStatus==PrStatus_Diag .OR. (.NOT.l_fix_ccb_cct) ) THEN
  ! ...print a warning if cloud-base and cloud-top are not set consistently
  error = -111
  DO i = 1, npnts
    IF ( bterm(i) ) THEN
      IF (      ( iccb(i)>0   .AND. icct(i)<iccb(i)    )                  &
           .OR. ( icct(i)>0   .AND. iccb(i)==0         )                  &
           .OR. ( lcbase(i)>0 .AND. lctop(i)<lcbase(i) )                  &
           .OR. ( lctop(i)>0  .AND. lcbase(i)==0       )                  &
           .OR. ( iccb(i)>0   .AND. lcbase(i)==0       )                  &
           .OR. ( iccb(i)==0  .AND. lcbase(i)>0        )  ) THEN
        ! - If the base has been set, the top must also be set (and vice versa)
        ! - The top mustn't be below the base.
        ! - If the current base/top has been set, the lowest base/top must
        !   also be set (and vice versa).
        WRITE(cmessage,"(A,4(A,I4))")                                     &
          "Convective cloud-base and cloud-top levels are inconsistent!", &
          " iccb = ",   iccb(i),    " icct = ",  icct(i),                 &
          " lcbase = ", lcbase(i),  " lctop = ", lctop(i)
        CALL ereport( RoutineName, error, cmessage )
      END IF
    END IF
  END DO
END IF


DO i = 1, npnts
  !
  IF ( flxkp1(i)  >   0.0) THEN
    !
    !---------------------------------------------------------------------
    ! SUM TOTAL CONDENSED WATER PER SECOND - ASSUMES THAT THE INITIAL
    ! CONVECTIVE LAYER IS UNSATURATED
    !---------------------------------------------------------------------
    !
    tcw(i) = tcw(i) + flxkp1(i) * ccw(i) / g
    !
    !---------------------------------------------------------------------
    ! SUM CONV CONDENSED WATER PATH - ASSUMES THAT THE INITIAL
    ! CONVECTIVE LAYER IS UNSATURATED
    !---------------------------------------------------------------------
    !
    cclwp(i) = cclwp(i) + xpkp1(i) * delpkp1(i) / g

  END IF
  !
  !---------------------------------------------------------------------
  ! CALCULATE CONVECTIVE CLOUD AMOUNT IF CONVECTION TERMINATES IN
  ! LAYER K AND TOTAL CONDENSED WATER PATH OVER A TIME STEP
  !
  ! UM DOCUMENTATION PAPER P27
  ! SECTION (9), EQUATION (37)
  !---------------------------------------------------------------------
  !
  IF ( bterm(i) .AND. tcw(i) >  0.0 ) THEN
    !
    IF ( tcw(i)  <   2.002e-6 ) tcw(i) = 2.002e-6
    !
    cca(i) = 0.7873 + 0.06 * LOG(tcw(i))
    IF (cca(i)  >   1.0) cca(i) = 1.0
    !
    IF (lcca(i) <= 0.0) THEN
      lcca(i) = 0.7873 + 0.06 * LOG(tcw(i))
      IF (lcca(i)  >   1.0) lcca(i) = 1.0
    END IF
    !
    tcw(i) = 0.0
    !
  END IF
END DO ! I loop over NPNTS
!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE con_rad_4a5a
END MODULE con_rad_4a5a_mod
