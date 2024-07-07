! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This subroutine Gathers the winds from 30201 and 30202 stash
! and saves them in separate arrays to be used by the FV-TRACK code.
! It loads the logical to activate different outputs (850mbar, TC
! levels, new TC method), user requested levels for TC and arrays
! to output the vorticity.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Climate Diagnostics

MODULE winds_for_track_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='WINDS_FOR_TRACK_MOD'

CONTAINS

SUBROUTINE winds_for_track(u_p_levs, v_p_levs, TC_track_p_levs,         &
                           tc_track_press, u_press, v_press,            &
                           u_p, v_p,                                    &
                           q_track, q_track_TC, q_track_TC_new,         &
                           u_850, v_850, u_tc_track, v_tc_track)

USE atm_fields_bounds_mod, ONLY: udims, vdims
USE track_mod,             ONLY: nlevs_avg, u_tc_new, v_tc_new
USE ereport_mod,           ONLY: ereport
USE errormessagelength_mod,ONLY: errormessagelength
USE umPrintMgr,            ONLY: umPrint, umMessage

USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim

IMPLICIT NONE

! Argument variables
INTEGER , INTENT(IN)::                                                  &
    u_p_levs, v_p_levs                                                  &
   ! Number of levels for U and V at pressure levels (30201,30202)
,   TC_track_p_levs
   ! Number of levels for Tropical cylone analysis

REAL, INTENT(IN) ::                                                     &
    tc_track_press(TC_track_p_levs)                                     &
   ! Pressure levels for TC tracking
,   u_press(u_p_levs)                                                   &
   ! U pressure levels as chosen for 30201
,   v_press(v_p_levs)
   ! V pressure levels as chosen for 30202

REAL, INTENT(IN)::                                                      &
    u_p(udims%i_start:udims%i_end, vdims%j_start:vdims%j_end, u_p_levs) &
   ! u at selected pressures levels chosen for 30201
,   v_p(udims%i_start:udims%i_end, vdims%j_start:vdims%j_end, v_p_levs)
   ! v at selected pressure levels chosen for 30202

LOGICAL, INTENT(IN) ::                                                  &
   q_track                                                              &
   ! q_track Switches the track output each 6 hours
,  q_track_TC                                                           &
   ! q_track Switches the track output each 6 hours
   ! for Tropical cyclones levels for warm core diagnosis.
,  q_track_TC_new
   ! q_track Switches the track output each 6 hours
   ! for the new method to track TC.

REAL,INTENT(INOUT) ::                                                   &
    u_850(udims%i_start:udims%i_end, vdims%j_start:vdims%j_end, 1)      &
   ! u at 850hPa
,   v_850(udims%i_start:udims%i_end, vdims%j_start:vdims%j_end,1)       &
   ! v at 850hPa
,   u_tc_track(udims%i_start:udims%i_end,                               &
             vdims%j_start:vdims%j_end,TC_track_p_levs)                 &
   ! u at selected pressures for TC warm-core identification

,   v_tc_track(udims%i_start:udims%i_end,                               &
             vdims%j_start:vdims%j_end,TC_track_p_levs)
   ! v at selected pressures for TC warm-core identification

! Local variables
INTEGER :: i,j,k,k_track ! indexes for do loops

REAL ::                                                                 &
 TC_new_Press(nlevs_avg)
 ! Pressure Levels to average for vorticity in new TC tracking

LOGICAL ::                                                              &
    u_tc_lev_flag(tc_track_P_levs)                                      &
,   v_tc_lev_flag(tc_track_P_levs)                                      &
   ! logical flags to activate if TC levels are not present
,   u_850_flag, v_850_flag                                              &
   ! Logical flag to activate if 850 level is not present
,   u_tc_new_flag(nlevs_avg),v_tc_new_flag(nlevs_avg)
   ! Logical flags to activate 'TC-new' levels if not present

!Error message handlers
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER  :: RoutineName='WINDS_FOR_TRACK'


! DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!GET 850 level
IF (q_track) THEN

  ! Declare 850 variables as TRUE
  u_850_flag=.TRUE.
  v_850_flag=.TRUE.

  ! Check there is the 850 level for U
  DO  k=1,u_p_levs

    ! If 850hPa level found
    IF (u_press(k) == 850.0) THEN
      ! Copy it
      DO j = vdims%j_start, vdims%j_end
        DO i = udims%i_start, udims%i_end
          u_850(i,j,1)=u_p(i,j,k)
        END DO
      END DO
      ! Switch off the flag
      u_850_flag=.FALSE.
    END IF
  END DO

  ! Check there is 850 level for V
  DO  k=1,v_p_levs
    ! If 850hPa level found
    IF (v_press(k) == 850.0) THEN
      DO j = vdims%j_start, vdims%j_end
        DO i = udims%i_start, udims%i_end
          v_850(i,j,1)=v_p(i,j,k)
        END DO
      END DO
      v_850_flag=.FALSE.
    END IF
  END DO

  ! If levels are not present stop the model:
  IF (u_850_flag) THEN
    WRITE(umMessage,'(A)')'**ERROR in TRACK**: U at 850hPa (30201) not chosen'
    CALL umPrint(umMessage,src='winds_for_track')
    icode   = 1
    WRITE (cmessage,'(A)') 'Please request U at 850hPa sect 30201'
    CALL ereport(routinename, icode, cmessage)
  END IF

  IF (v_850_flag) THEN
    WRITE(umMessage,'(A)')'**ERROR in TRACK**: V at 850hPa (30202) not chosen'
    CALL umPrint(umMessage,src='winds_for_track')
    icode   = 1
    WRITE (cmessage,'(A)') 'Please request V at 850hPa sect 30202'
    CALL ereport(routinename, icode, cmessage)
  END IF

  ! End if over 850hPa level
END IF

!Get TC levels
IF (q_track_TC) THEN

  ! Initialize tc_lev_flag
  DO k=1,TC_track_p_levs
    u_tc_lev_flag(k)=.TRUE.
    v_tc_lev_flag(k)=.TRUE.
  END DO

  ! Loop for levels where u_tc_track is requested
  DO  k=1,u_p_levs
    DO  k_track=1,TC_track_p_levs
      IF (u_press(k) == TC_track_PRESS(k_track)) THEN
        DO j = vdims%j_start, vdims%j_end
          DO i = udims%i_start, udims%i_end
            u_tc_track(i,j,k_track)=u_p(i,j,k)
          END DO
        END DO

        u_tc_lev_flag(k_track)=.FALSE.
      END IF
    END DO
  END DO

  DO  k=1,v_p_levs
    DO  k_track=1,TC_track_P_LEVS
      IF (v_press(k) == TC_track_PRESS(k_track)) THEN
        DO j = vdims%j_start, vdims%j_end
          DO i = udims%i_start, udims%i_end
            v_tc_track(i,j,k_track)=v_p(i,j,k)
          END DO
        END DO
        v_tc_lev_flag(k_track)=.FALSE.
      END IF
    END DO
  END DO

  ! If any level has not been found stop the model!
  DO k=1,TC_track_p_levs
    IF (u_tc_lev_flag(k)) THEN
      WRITE(umMessage,'(A)')                                             &
      '**ERROR TRACK **: U for TC levels have been not selected at 30201'
      CALL umPrint(umMessage,src='winds_for_track')
      WRITE(umMessage,'(A,F3.0)')'U Pressure level needed at 30201 ',    &
      tc_track_press(k)
      CALL umPrint(umMessage,src='winds_for_track')
      icode   = 1
      WRITE (cmessage,'(A)')                                             &
      'TRACK ERROR (TC). There are levels not selected at stash:30201'
      CALL ereport(routinename, icode, cmessage)
    END IF

    IF (v_tc_lev_flag(k)) THEN
      WRITE(umMessage,'(A)')                                             &
      '**ERROR TRACK **: V for TC levels have been not selected at 30202'
      CALL umPrint(umMessage,src='winds_for_track')
      WRITE(umMessage,'(A,F3.0)')'V Pressure level needed at 30202: ',   &
      tc_track_press(k)
      CALL umPrint(umMessage,src='winds_for_track')
      icode   = 1
      WRITE (cmessage,'(A)')                                             &
      'TRACK ERROR (TC). There are levels not selected at stash:30202'
      CALL ereport(routinename, icode, cmessage)
    END IF

  END DO

  ! End if over TC levels
END IF

! Get levels for new TC tracking
! It gathers 850,700 and 600 to be meaned and tracked for
! nbot_tc-ntop_tc filtered fields
IF (q_track_TC_new) THEN

  ! Allocate u_tc_new and v_tc_new
  IF (.NOT. ALLOCATED(u_tc_new)) THEN
    ALLOCATE(u_tc_new(udims%i_start:udims%i_end,                         &
             vdims%j_start:vdims%j_end,nlevs_avg) )
  END IF

  IF (.NOT. ALLOCATED(v_tc_new)) THEN
    ALLOCATE(v_tc_new(udims%i_start:udims%i_end,                         &
             vdims%j_start:vdims%j_end,nlevs_avg) )
  END IF

  ! Set up levels to average
  TC_new_Press(1) = 850.0
  TC_new_Press(2) = 700.0
  TC_new_Press(3) = 600.0

  ! Initialize tc_new_flag
  DO k=1,nlevs_avg
    u_tc_new_flag(k)=.TRUE.
    v_tc_new_flag(k)=.TRUE.
  END DO

  DO  k=1,u_p_levs
    DO  k_track=1,nlevs_avg
      IF (u_press(k) == TC_new_PRESS(k_track)) THEN
        DO j = vdims%j_start, vdims%j_end
          DO i = udims%i_start, udims%i_end
            u_tc_new(i,j,k_track)=u_p(i,j,k)
          END DO
        END DO
        u_tc_new_flag(k_track)=.FALSE.
      END IF
    END DO
  END DO

  DO  k=1,v_p_levs
    DO  k_track=1,nlevs_avg
      IF (u_press(k) == TC_new_PRESS(k_track)) THEN
        DO j = vdims%j_start, vdims%j_end
          DO i = udims%i_start, udims%i_end
            v_tc_new(i,j,k_track)=v_p(i,j,k)
          END DO
        END DO
        v_tc_new_flag(k_track)=.FALSE.
      END IF
    END DO
  END DO

  ! If any level has not been found stop the model!
  DO k=1,nlevs_avg
    IF (u_tc_new_flag(k)) THEN
      WRITE(umMessage,'(A)')                                             &
      '**ERROR TRACK **: U for new-TC levels not selected at 30201'
      CALL umPrint(umMessage,src='winds_for_track')
      WRITE(umMessage,'(A,F3.0)')'U Pressure level needed at 30201 ',    &
      TC_new_PRESS(k)
      CALL umPrint(umMessage,src='winds_for_track')
      icode   = 1
      WRITE (cmessage,'(A)')                                             &
      'TRACK ERROR (TC). There are levels not selected at stash:30201'
      CALL ereport(routinename, icode, cmessage)
    END IF

    IF (v_tc_new_flag(k)) THEN
      WRITE(umMessage,'(A)')                                             &
      '**ERROR TRACK **: V for new-TC levels not selected at 30202'
      CALL umPrint(umMessage,src='winds_for_track')
      WRITE(umMessage,'(A,F3.0)')'V Pressure level needed at 30202 ',    &
      TC_new_PRESS(k)
      CALL umPrint(umMessage,src='winds_for_track')
      icode   = 1
      WRITE (cmessage,'(A)')                                             &
      'TRACK ERROR (TC). There are levels not selected at stash:30202'
      CALL ereport(routinename, icode, cmessage)
    END IF

  END DO

  ! End if over q_track_TC_new
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE winds_for_track

END MODULE winds_for_track_mod
