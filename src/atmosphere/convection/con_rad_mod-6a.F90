! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Convective cloud microphysics routine
!
MODULE con_rad_6a_mod

IMPLICIT NONE

! Description:
!   Calculates convective cloud top, base and amount
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CON_RAD_6A_MOD'

CONTAINS

! Subroutine Interface:

SUBROUTINE con_rad_6a (k, npnts, start_lev,                        &
                       ccwk, ccwkp1, tmp_ccwkp1, flxkp1, delpkp1,  &
                       l_q_interact, bterm,                        &
                       tcw, cclwp, lcca, lcbase, lctop,            &
                       iccb, icct, cca)

USE science_fixes_mod,      ONLY: l_fix_ccb_cct
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr,             ONLY: PrintStatus, PrStatus_Diag

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE planet_constants_mod, ONLY: g

IMPLICIT NONE



!----------------------------------------------------------------------
! Arguments with INTENT(IN)
!----------------------------------------------------------------------
INTEGER,INTENT(IN) :: k               ! Present model layer number
INTEGER,INTENT(IN) :: npnts           ! Vector length
INTEGER,INTENT(IN) :: start_lev(npnts)! Level at which convection initiated

REAL,INTENT(IN) :: ccwk(npnts)      ! Total condensate in level k (kg/kg)
REAL,INTENT(IN) :: ccwkp1(npnts)    ! Total condensate in level k+1 (kg/kg)
REAL,INTENT(IN) :: tmp_ccwkp1(npnts)! Total condensate in level k+1 (kg/kg)
                                    ! before precipitation
REAL,INTENT(IN) :: flxkp1(npnts)    ! Parcel mass flux in layer k+1 (Pa/s)
REAL,INTENT(IN) :: delpkp1(npnts)   ! pressure difference across layer k+1 (Pa)

LOGICAL,INTENT(IN) :: l_q_interact  ! True if PC2 is switched on
LOGICAL,INTENT(IN) :: bterm(npnts)  ! Mask for parcels which terminate
                                    ! in layer k+1

! ----------------------------------------------------------------------
! Arguments with INTENT(IN/OUT)
! ----------------------------------------------------------------------
REAL,INTENT(INOUT) :: tcw(npnts)    ! Total condensed water (kg/m**2/s)
                                    ! IN summed to layer k
                                    ! OUT summed to layer k+1
REAL,INTENT(INOUT) :: cclwp(npnts)  ! Condensed water path (kg/m**2)
                                    ! IN summed to layer k
                                    ! OUT summed to layer k+1
REAL,INTENT(INOUT) :: lcca(npnts)   ! Lowest conv. cloud amount (%)

INTEGER,INTENT(INOUT) :: lcbase(npnts)! Lowest conv. cloud base level
INTEGER,INTENT(INOUT) :: lctop(npnts) ! Lowest conv. cloud top level

INTEGER,INTENT(INOUT) :: iccb(npnts)  ! convective cloud base_level
INTEGER,INTENT(INOUT) :: icct(npnts)  ! convective cloud top level
! Need to be inout as need to be remembered between subsequent calls to convec2

REAL,INTENT(OUT) :: cca(npnts)      ! convective cloud amount (%)

!-------------------------------------------------------------------------------
! Local variables
!-------------------------------------------------------------------------------
INTEGER :: i          ! Loop counter

INTEGER :: error  ! Error code for ereport
CHARACTER (LEN=errormessagelength) :: cmessage  ! String to store error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CON_RAD_6A'

! ----------------------------------------------------------------------
!  Calculate cloud base and lowest cloud base
!  when cloud base set zero total condensed water
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
!  Calculate cloud base and cloud top (and lowest base and top)
! ----------------------------------------------------------------------

IF ( l_fix_ccb_cct ) THEN
  ! Corrected version of the ccb / cct calculation; consistently
  ! set the convective base and top to just be the start-level and
  ! termination level of the ascent, ignoring cloud.
  ! This is consistent with how iccb and icct are used in the 6A
  ! convection scheme, to determine the depth of the ascent, not just
  ! the cloudy part.

  DO i=1, npnts
    IF ( k == start_lev(i) ) THEN
      iccb(i) = k
      cclwp(i)  = 0.0
      tcw(i)   = 0.0
      IF ( lcbase(i) == 0 )  lcbase(i) = iccb(i)
    END IF
    IF (bterm(i)) THEN
      icct(i)  = k+1
      IF ( lctop(i) == 0 )  lctop(i) = icct(i)
    END IF
  END DO

ELSE
  ! Original, confusing code; sets the base to be the first
  ! cloudy level in the ascent, except if using PC2 (l_q_interact),
  ! in which case the base is forced to be at the ascent start level
  ! even if the parcel contains no cloud at k or k+1.
  ! Then the conv cloud top is set to just be the top of the ascent,
  ! regardless of whether the parcel contained cloud, or if using PC2.

  DO i = 1, npnts
  
    IF ( ccwk(i)  <=  0.0 .AND. tmp_ccwkp1(i)  >   0.0 ) THEN
      !       Assuming initial parcel condensate is zero
      iccb(i)   = k+1
      cclwp(i)  = 0.0

      IF ( lcbase(i)  ==  0 ) THEN
        !         Lowest cloud base
        lcbase(i) = k+1
      END IF  ! If_lcbase_1

      ! Note that, if not l_q_interact, ccwk must always be zero
      ! in the 1st level of ascent, so if any condensation occurs,
      ! the above check will set the cloud-base.  However, if using PC2
      ! (l_q_interact), the initial parcel can contain condensate (ccwk>0),
      ! so need an additional check...

    ELSE IF ( l_q_interact .AND. k  ==  start_lev(i) ) THEN

      ! ...If the initial parcel contains condensate, then cloud-base level
      ! is the parcel initial level, not k+1.

      ! Non-zero initial parcel condensate (initialized to environment)
      iccb(i)  = k
      cclwp(i) = 0.0
      tcw(i)   = 0.0

      IF ( lcbase(i)  ==  0 ) THEN
        !         Lowest cloud base
        lcbase(i) = k
      END IF  ! If_lcbase_2

    END IF

    IF (bterm(i)) THEN
      icct(i)  = k+1
    END IF

    IF (bterm(i) .AND.  lctop(i) == 0 ) THEN
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

  IF ( flxkp1(i)  >   0.0) THEN

    !---------------------------------------------------------------------
    ! Sum total condensed water per second - assumes that the initial
    ! convective layer is unsaturated.
    ! Uses the CCW before precipitation
    !---------------------------------------------------------------------
    tcw(i)   = tcw(i)   + flxkp1(i) * tmp_ccwkp1(i) / g

    !---------------------------------------------------------------------
    ! Sum conv condensed water path - assumes that the initial
    ! convective layer is unsaturated
    ! Uses the CCW after precipitation
    !---------------------------------------------------------------------
    cclwp(i) = cclwp(i) + ccwkp1(i) * delpkp1(i) / g

  END IF

  !---------------------------------------------------------------------
  ! calculate convective cloud amount if convection terminates in
  ! layer k and total condensed water path over a time step
  !
  ! UM documentation paper p27
  ! section (9), equation (37)
  !---------------------------------------------------------------------
  IF ( bterm(i) .AND. tcw(i) >  0.0 ) THEN

    IF ( tcw(i)  <   2.002e-6 ) tcw(i) = 2.002e-6

    cca(i) = 0.7873 + 0.06 * LOG(tcw(i))
    IF (cca(i)  >   1.0) cca(i) = 1.0

    IF (lcca(i) <= 0.0) THEN
      lcca(i) = 0.7873 + 0.06 * LOG(tcw(i))
      IF (lcca(i)  >   1.0) lcca(i) = 1.0
    END IF

    tcw(i) = 0.0

  END IF
END DO ! I loop over NPNTS

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE con_rad_6a
END MODULE con_rad_6a_mod
