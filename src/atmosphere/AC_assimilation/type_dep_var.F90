! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE TYPE_DEP_VAR -------------------------------------------
!
!    Purpose : Process &ACP Namelist arrays which are
!              Observation Type Dependent.
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!
!
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: AC Assimilation
MODULE type_dep_var_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TYPE_DEP_VAR_MOD'

CONTAINS

SUBROUTINE type_dep_var(timeb,timea,tgetobb,tgetoba,radinf,obthin,&
                         cscale_start,cscale_obtime,cscale_end,   &
                         icode,cmessage)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE UM_ParCore, ONLY: mype
USE comobs_mod, ONLY: nobtypmx
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE ac_control_mod

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE



INTEGER :: timea   (nobtypmx)         !IN time window after ob time
INTEGER :: timeb   (nobtypmx)         !IN time window before ob time
INTEGER :: tgetoba (nobtypmx)         !IN window for getting obs (a)
INTEGER :: tgetobb (nobtypmx)         !IN window for getting obs (b)
INTEGER :: radinf  (nobtypmx)         !IN horiz infl radius
INTEGER :: obthin  (nobtypmx)         !IN obs thinning factors
INTEGER :: cscale_start  (nobtypmx)   !IN horiz cor scale at start
INTEGER :: cscale_obtime (nobtypmx)   !IN   "    "    "   "  ob time
INTEGER :: cscale_end    (nobtypmx)   !IN   "    "    "   "   end
INTEGER :: icode                      !OUT error code and message
CHARACTER(LEN=errormessagelength) :: cmessage

!     Local arrays/variables.

INTEGER :: iscfact(nobtypmx)
INTEGER :: jobt,jobt2,ac_type,j
INTEGER :: first_type,last_type,n_types
LOGICAL :: lfound

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TYPE_DEP_VAR'


30  FORMAT (a,f5.2,a,i4)

!     Process TIMEB
!     -------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
DO jobt=1,nobtypmx
  IF (timeb(jobt) >  0) THEN
    ac_type = timeb(jobt)/1000
    timeb(jobt) = MOD(timeb(jobt),1000)
    IF (timeb(jobt) <  1 .OR. timeb(jobt) >  999) THEN
      icode=1
      cmessage = 'TYPE_DEP_VAR : Invalid value given for TIMEB'
      GO TO 999
    END IF
    IF (ac_type == 999) THEN
      !           Use new value for ALL Obs types
      DO jobt2=1,nobtypmx
        def_timeb(jobt2)=timeb(jobt)
      END DO
      WRITE(umMessage,*) 'TIMEB : Defaults changed to ',                 &
          timeb(jobt),' mins for all obs types.'
      CALL umPrint(umMessage,src='type_dep_var',pe=0)
    ELSE
      lfound = .FALSE.
      DO jobt2=1,nobtypmx
        IF (ac_type == master_ac_types(jobt2)) THEN
          lfound = .TRUE.
          def_timeb(jobt2)=timeb(jobt)
          WRITE(umMessage,*) 'TIMEB : Default changed to ',               &
              timeb(jobt),' mins for obs type ',ac_type
          CALL umPrint(umMessage,src='type_dep_var',pe=0)
        END IF
      END DO
      IF (.NOT. lfound) THEN
        icode=1
        cmessage =                                                &
        'TYPE_DEP_VAR : Invalid AC Obs type used in TIMEB'
        GO TO 999
      END IF
    END IF
  END IF
END DO

!     Process TIMEA
!     -------------
DO jobt=1,nobtypmx
  IF (timea(jobt) >  0) THEN
    ac_type = timea(jobt)/1000
    timea(jobt) = MOD(timea(jobt),1000)
    IF (timea(jobt) <  1 .OR. timea(jobt) >  999) THEN
      icode=1
      cmessage = 'TYPE_DEP_VAR : Invalid value given for TIMEA'
      GO TO 999
    END IF
    IF (ac_type == 999) THEN
      !           Use new value for ALL Obs types
      DO jobt2=1,nobtypmx
        def_timea(jobt2)=timea(jobt)
      END DO
      WRITE(umMessage,*) 'TIMEA : Defaults changed to ',                 &
          timea(jobt),' mins for all obs types.'
      CALL umPrint(umMessage,src='type_dep_var',pe=0)
    ELSE
      lfound = .FALSE.
      DO jobt2=1,nobtypmx
        IF (ac_type == master_ac_types(jobt2)) THEN
          lfound = .TRUE.
          def_timea(jobt2)=timea(jobt)
          WRITE(umMessage,*) 'TIMEA : Default changed to ',              &
              timea(jobt),' mins for obs type ',ac_type
          CALL umPrint(umMessage,src='type_dep_var',pe=0)
        END IF
      END DO
      IF (.NOT. lfound) THEN
        icode=1
        cmessage =                                                &
        'TYPE_DEP_VAR : Invalid AC Obs type used in TIMEA'
        GO TO 999
      END IF
    END IF
  END IF
END DO

!     Process TGETOBB
!     ---------------
DO jobt=1,nobtypmx
  IF (tgetobb(jobt) >  0) THEN
    ac_type = tgetobb(jobt)/1000
    tgetobb(jobt) = MOD(tgetobb(jobt),1000)
    IF (tgetobb(jobt) <  1 .OR. tgetobb(jobt) >  999) THEN
      icode=1
      cmessage = 'TYPE_DEP_VAR : Invalid value given in TGETOBB'
      GO TO 999
    END IF
    IF (ac_type == 999) THEN
      !           Use new value for ALL Obs types
      DO jobt2=1,nobtypmx
        def_tgetobb(jobt2)=tgetobb(jobt)
      END DO
      WRITE(umMessage,*) 'TGETOBB : Default changed to ',               &
          tgetobb(jobt),' mins for all obs types'
      CALL umPrint(umMessage,src='type_dep_var',pe=0)
    ELSE
      lfound = .FALSE.
      DO jobt2=1,nobtypmx
        IF (ac_type == master_ac_types(jobt2)) THEN
          lfound = .TRUE.
          def_tgetobb(jobt2)=tgetobb(jobt)
          WRITE(umMessage,*) 'TGETOBB : Default changed to ',            &
              tgetobb(jobt),' mins for obs type ',ac_type
          CALL umPrint(umMessage,src='type_dep_var',pe=0)
        END IF
      END DO
      IF (.NOT. lfound) THEN
        icode=1
        cmessage =                                                &
        'TYPE_DEP_VAR : Invalid AC Obs type used in TGETOBB'
        GO TO 999
      END IF
    END IF
  END IF
END DO

!     Process TGETOBA
!     ---------------
DO jobt=1,nobtypmx
  IF (tgetoba(jobt) >  0) THEN
    ac_type = tgetoba(jobt)/1000
    tgetoba(jobt) = MOD(tgetoba(jobt),1000)
    IF (tgetoba(jobt) <  1 .OR. tgetoba(jobt) >  999) THEN
      icode=1
      cmessage = 'TYPE_DEP_VAR : Invalid value given in TGETOBA'
      GO TO 999
    END IF
    IF (ac_type == 999) THEN
      !           Use new value for ALL Obs types
      DO jobt2=1,nobtypmx
        def_tgetoba(jobt2)=tgetoba(jobt)
      END DO
      WRITE(umMessage,*) 'TGETOBA : Default changed to ',                &
          tgetoba(jobt),' mins for all obs types'
      CALL umPrint(umMessage,src='type_dep_var',pe=0)
    ELSE
      lfound = .FALSE.
      DO jobt2=1,nobtypmx
        IF (ac_type == master_ac_types(jobt2)) THEN
          lfound = .TRUE.
          def_tgetoba(jobt2)=tgetoba(jobt)
          WRITE(umMessage,*) 'TGETOBA : Default changed to ',            &
              tgetoba(jobt),' mins for obs type',ac_type
          CALL umPrint(umMessage,src='type_dep_var',pe=0)
        END IF
      END DO
      IF (.NOT. lfound) THEN
        icode=1
        cmessage =                                                &
        'TYPE_DEP_VAR : Invalid AC Obs type used in TGETOBA'
        GO TO 999
      END IF
    END IF
  END IF
END DO

!     Process RADINF
!     --------------
DO jobt=1,nobtypmx
  IF (radinf(jobt) >  0) THEN
    ac_type = radinf(jobt)/1000
    radinf(jobt) = MOD(radinf(jobt),1000)
    IF (radinf(jobt) <  1 .OR. radinf(jobt) >  999) THEN
      icode=1
      cmessage = 'TYPE_DEP_VAR : Invalid value given in RADINF'
      GO TO 999
    END IF
    IF (ac_type == 999) THEN
      !           Use new value for ALL Obs types
      DO jobt2=1,nobtypmx
        def_radinf(jobt2)=radinf(jobt)/100.0
      END DO
      IF (mype == 0)                                               &
      PRINT 30, ' RADINF  : Defaults changed to ',                &
                  def_radinf(1),' for all obs types'
    ELSE
      lfound = .FALSE.
      DO jobt2=1,nobtypmx
        IF (ac_type == master_ac_types(jobt2)) THEN
          lfound = .TRUE.
          def_radinf(jobt2)=radinf(jobt)/100.0
          IF (mype == 0)                                           &
          PRINT 30, ' RADINF  : Default changed to ',             &
                      def_radinf(jobt2),' for obs type ',ac_type
        END IF
      END DO
      IF (.NOT. lfound) THEN
        icode=1
        cmessage =                                                &
        'TYPE_DEP_VAR : Invalid AC Obs type used in RADINF'
        GO TO 999
      END IF
    END IF
  END IF
END DO


!     Process OBTHIN
!     --------------
DO jobt=1,nobtypmx
  IF (obthin(jobt) >  0) THEN
    ac_type = obthin(jobt)/1000
    obthin(jobt) = MOD(obthin(jobt),1000)
    IF (obthin(jobt) <  1 .OR. obthin(jobt) >  999) THEN
      icode=1
      cmessage = 'TYPE_DEP_VAR : Invalid value given in OBTHIN'
      GO TO 999
    END IF
    IF (ac_type == 999) THEN
      !           Use new value for ALL Obs types
      DO jobt2=1,nobtypmx
        def_obthin(jobt2)=obthin(jobt)
      END DO
      WRITE(umMessage,*) ' OBTHIN  : Defaults changed to ',              &
          def_obthin(1),' for all obs types'
      CALL umPrint(umMessage,src='type_dep_var',pe=0)
    ELSE
      lfound = .FALSE.
      DO jobt2=1,nobtypmx
        IF (ac_type == master_ac_types(jobt2)) THEN
          lfound = .TRUE.
          def_obthin(jobt2)=obthin(jobt)
          WRITE(umMessage,*) ' OBTHIN  : Default changed to ',           &
              def_obthin(jobt2),' for obs type ',ac_type
          CALL umPrint(umMessage,src='type_dep_var',pe=0)
        END IF
      END DO
      IF (.NOT. lfound) THEN
        icode=1
        cmessage =                                                &
        'TYPE_DEP_VAR : Invalid AC Obs type used in OBTHIN'
        GO TO 999
      END IF
    END IF
  END IF
END DO

!     Process CSCALE_START
!     --------------------
DO jobt=1,nobtypmx
  IF (cscale_start(jobt) >  0) THEN
    ac_type = cscale_start(jobt)/1000
    cscale_start(jobt) = MOD(cscale_start(jobt),1000)
    IF (cscale_start(jobt) <  1 .OR.                              &
        cscale_start(jobt) >  999) THEN
      icode=1
      cmessage =                                                  &
      'TYPE_DEP_VAR : Invalid value given in CSCALE_START'
      GO TO 999
    END IF
    IF (ac_type == 999) THEN
      !           Use new value for ALL Obs types
      DO jobt2=1,nobtypmx
        def_cscale_start(jobt2)=cscale_start(jobt)
      END DO
      WRITE(umMessage,*) 'CSCALE_START  : Defaults changed to ',   &
          cscale_start(jobt),' km for all obs types.'
      CALL umPrint(umMessage,src='type_dep_var',pe=0)
    ELSE
      lfound = .FALSE.
      DO jobt2=1,nobtypmx
        IF (ac_type == master_ac_types(jobt2)) THEN
          lfound = .TRUE.
          def_cscale_start(jobt2)=cscale_start(jobt)
          WRITE(umMessage,*) 'CSCALE_START  : Default changed to ',      &
              cscale_start(jobt),' km for obs type ',ac_type
          CALL umPrint(umMessage,src='type_dep_var',pe=0)
        END IF
      END DO
      IF (.NOT. lfound) THEN
        icode=1
        cmessage =                                                &
        'TYPE_DEP_VAR : Invalid AC Obs type used in CSCALE_START'
        GO TO 999
      END IF
    END IF
  END IF
END DO

!     Process CSCALE_OBTIME
!     ---------------------
DO jobt=1,nobtypmx
  IF (cscale_obtime(jobt) >  0) THEN
    ac_type = cscale_obtime(jobt)/1000
    cscale_obtime(jobt) = MOD(cscale_obtime(jobt),1000)
    IF (cscale_obtime(jobt) <  1 .OR.                             &
        cscale_obtime(jobt) >  999) THEN
      icode=1
      cmessage =                                                  &
      'TYPE_DEP_VAR : Invalid value given in CSCALE_OBTIME'
      GO TO 999
    END IF
    IF (ac_type == 999) THEN
      !           Use new value for ALL Obs types
      DO jobt2=1,nobtypmx
        def_cscale_obtime(jobt2)=cscale_obtime(jobt)
      END DO
      WRITE(umMessage,*) 'CSCALE_OBTIME : Defaults changed to ',         &
          cscale_obtime(jobt),' km for all obs types.'
      CALL umPrint(umMessage,src='type_dep_var',pe=0)
    ELSE
      lfound = .FALSE.
      DO jobt2=1,nobtypmx
        IF (ac_type == master_ac_types(jobt2)) THEN
          lfound = .TRUE.
          def_cscale_obtime(jobt2)=cscale_obtime(jobt)
          WRITE(umMessage,*) 'CSCALE_OBTIME : Default changed to ',      &
              cscale_obtime(jobt),' km for obs type ',ac_type
          CALL umPrint(umMessage,src='type_dep_var',pe=0)
        END IF
      END DO
      IF (.NOT. lfound) THEN
        icode=1
        cmessage =                                                &
        'TYPE_DEP_VAR : Invalid Obs type used in CSCALE_OBTIME'
        GO TO 999
      END IF
    END IF
  END IF
END DO

!     Process CSCALE_END
!     ------------------
DO jobt=1,nobtypmx
  IF (cscale_end(jobt) >  0) THEN
    ac_type = cscale_end(jobt)/1000
    cscale_end(jobt) = MOD(cscale_end(jobt),1000)
    IF (cscale_end(jobt) <  1 .OR.                                &
        cscale_end(jobt) >  999) THEN
      icode=1
      cmessage =                                                  &
      'TYPE_DEP_VAR : Invalid value given in CSCALE_END'
      GO TO 999
    END IF
    IF (ac_type == 999) THEN
      !           Use new value for ALL Obs types
      DO jobt2=1,nobtypmx
        def_cscale_end(jobt2)=cscale_end(jobt)
      END DO
      WRITE(umMessage,*) 'CSCALE_END    : Defaults changed to ',         &
          cscale_end(jobt),' km for all obs types.'
      CALL umPrint(umMessage,src='type_dep_var',pe=0)
    ELSE
      lfound = .FALSE.
      DO jobt2=1,nobtypmx
        IF (ac_type == master_ac_types(jobt2)) THEN
          lfound = .TRUE.
          def_cscale_end(jobt2)=cscale_end(jobt)
          WRITE(umMessage,*) 'CSCALE_END    : Default changed to ',      &
              cscale_end(jobt),' km for obs type ',ac_type
          CALL umPrint(umMessage,src='type_dep_var',pe=0)
        END IF
      END DO
      IF (.NOT. lfound) THEN
        icode=1
        cmessage =                                                &
        'TYPE_DEP_VAR : Invalid AC Obs type used in CSCALE_END'
        GO TO 999
      END IF
    END IF
  END IF
END DO

!     Process NO_SCFACT
!     -----------------
DO jobt=1,nobtypmx
  iscfact(jobt)   = no_scfact(jobt)
  no_scfact(jobt) = 1   !  Apply scale factor to all types.
END DO
DO jobt=1,nobtypmx
  IF (iscfact(jobt) >  0) THEN
    lfound = .FALSE.
    DO jobt2=1,nobtypmx
      IF (iscfact(jobt) == master_ac_types(jobt2)) THEN
        lfound = .TRUE.
        no_scfact(jobt2) = 0   !  No scale factor for this type
      END IF
    END DO
    IF (.NOT. lfound) THEN
      icode=1
      cmessage =                                                  &
      'TYPE_DEP_VAR : Invalid Obs type used in NO_SCFACT'
      GO TO 999
    END IF
  END IF
END DO

n_types = 0
DO jobt=1,nobtypmx
  IF (master_ac_types(jobt) >  0) THEN
    n_types = n_types+1
  END IF
END DO

DO j=1,(n_types+6)/7
  first_type=(j-1)*7+1
  last_type =MIN(j*7,n_types)
  IF (mype == 0) THEN
    WRITE(umMessage,*) ' '
    CALL umPrint(umMessage,src='type_dep_var')
    WRITE(umMessage,'(A,7I8)')   ' AC Obs Types      ',                  &
        (master_ac_types(jobt),jobt=first_type,last_type)
    CALL umPrint(umMessage,src='type_dep_var')
    WRITE(umMessage,'(A,7F8.1)') ' TIMEB (mins)      ',                  &
        (def_timeb(jobt),jobt=first_type,last_type)
    CALL umPrint(umMessage,src='type_dep_var')
    WRITE(umMessage,'(A,7F8.1)') ' TIMEA (mins)      ',                  &
        (def_timea(jobt),jobt=first_type,last_type)
    CALL umPrint(umMessage,src='type_dep_var')
    WRITE(umMessage,'(A,7F8.1)') ' TGETOBB (mins)    ',                  &
        (def_tgetobb(jobt),jobt=first_type,last_type)
    CALL umPrint(umMessage,src='type_dep_var')
    WRITE(umMessage,'(A,7F8.1)') ' TGETOBA (mins)    ',                  &
        (def_tgetoba(jobt),jobt=first_type,last_type)
    CALL umPrint(umMessage,src='type_dep_var')
    WRITE(umMessage,'(A,7F8.2)') ' RADINF            ',                  &
        (def_radinf(jobt),jobt=first_type,last_type)
    CALL umPrint(umMessage,src='type_dep_var')
    WRITE(umMessage,'(A,7I8)') ' OBTHIN            ',                    &
        (def_obthin(jobt),jobt=first_type,last_type)
    CALL umPrint(umMessage,src='type_dep_var')
    WRITE(umMessage,'(A,7F8.1)') ' CSCALE_START (km) ',                  &
        (def_cscale_start(jobt),jobt=first_type,last_type)
    CALL umPrint(umMessage,src='type_dep_var')
    WRITE(umMessage,'(A,7F8.1)') ' CSCALE_OBTIME (km)',                  &
        (def_cscale_obtime(jobt),jobt=first_type,last_type)
    CALL umPrint(umMessage,src='type_dep_var')
    WRITE(umMessage,'(A,7F8.1)') ' CSCALE_END (km)   ',                  &
        (def_cscale_end(jobt),jobt=first_type,last_type)
    CALL umPrint(umMessage,src='type_dep_var')
    WRITE(umMessage,'(A,7I8)')   ' NO_SCFACT         ',                  &
        (no_scfact(jobt),jobt=first_type,last_type)
    CALL umPrint(umMessage,src='type_dep_var')
  END IF
END DO

!     Convert Horizontal Correlation Scales from Km to metres.
DO jobt=1,nobtypmx
  def_cscale_start(jobt)  = def_cscale_start(jobt)  * 1000.0
  def_cscale_obtime(jobt) = def_cscale_obtime(jobt) * 1000.0
  def_cscale_end(jobt)    = def_cscale_end(jobt)    * 1000.0
END DO

999  CONTINUE
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE type_dep_var
END MODULE type_dep_var_mod
