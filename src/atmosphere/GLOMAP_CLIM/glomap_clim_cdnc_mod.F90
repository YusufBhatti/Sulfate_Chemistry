! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!   To obtain GLOMAP-MODE input to the structure UKCA_CDNC
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
MODULE glomap_clim_cndc_mod

IMPLICIT NONE

TYPE stashcodetype
  INTEGER :: section        ! section code
  INTEGER :: item           ! item code
  LOGICAL :: put_stash      ! item has to be written to stash
ENDTYPE stashcodetype

! List of drydiam fields to write to stash (when not l_glomap_clim_radaer)
TYPE(stashcodetype), ALLOCATABLE, PRIVATE :: drydiamfields(:)

! List of cdnc fields to write to stash
TYPE(stashcodetype), ALLOCATABLE, PRIVATE :: cdncfields(:)

! Length of cdncfields array
INTEGER, PARAMETER, PRIVATE :: n_fields_cdnc = 4

! Length of drydiamfields array
INTEGER, PARAMETER, PRIVATE :: n_fields_drydiam = 6

! Field to write to stash
REAL, ALLOCATABLE, PRIVATE :: outfield(:,:,:)

! Temperature
REAL, ALLOCATABLE, PUBLIC  :: t_theta_levels(:,:,:)

! Pressure on theta levels without halos
REAL, ALLOCATABLE, PUBLIC  :: pr_theta_levels(:,:,:)

! Sat vap pressure with respect to liquid water irrespective of temperature
REAL, ALLOCATABLE, PRIVATE :: qsvp(:,:,:)

! sat mixing ratio with respect to liquid water irrespective of temperature
REAL, ALLOCATABLE, PRIVATE :: qsmr(:,:,:)

! weighted cdnc = total cdnc * cldflg [m-3]
REAL, ALLOCATABLE, PRIVATE :: cdnc(:,:,:)

! weighted <cdnc^-1/3> = <cdnc^-1/3> * cldflg [m-3]
REAL, ALLOCATABLE, PRIVATE :: cdnc3(:,:,:)

! weighted cdnc = total cdnc * liq_cloud_frac [m-3]
REAL, ALLOCATABLE, PRIVATE :: cdncwt(:,:,:)

! to pass as argument to ukca_activate
REAL, ALLOCATABLE, PRIVATE :: bl_tke(:,:,:)

! to pass as argument to ukca_activate
REAL, ALLOCATABLE, PRIVATE :: mode_tracers(:,:,:,:)

! to pass as argument to ukca_activate
REAL, ALLOCATABLE, PRIVATE :: mode_diags(:,:,:,:)

! Local buffer read from D1 and its size.
REAL, ALLOCATABLE, PRIVATE :: buffer(:,:,:)
INTEGER, PRIVATE           :: buffer_size

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName ='GLOMAP_CLIM_CNDC_MOD'

PUBLIC :: glomap_clim_arg_act_get_cdnc

PRIVATE :: glomap_clim_cdnc_identify_fields,    &
           glomap_clim_drydiam_identify_fields, &
           set_field_cdnc,                      &
           set_field_drydiam,                   &
           glomap_clim_cdnc_bl_tke

CONTAINS

SUBROUTINE glomap_clim_arg_act_get_cdnc ( ukca_cdnc, stashwork54 )

USE atm_fields_real_mod,          ONLY: &
    p_theta_levels,                     &
    exner_theta_levels,                 &
    theta,                              &
    q,                                  &
    w,                                  &
    cf_liquid,                          &
    qcl

USE atm_step_local,               ONLY: &
    first_atmstep_call

USE ereport_mod,                  ONLY: &
    ereport

USE errormessagelength_mod,       ONLY: &
    errormessagelength

USE glomap_clim_calc_drydiam_mod, ONLY: &
    glomap_clim_calc_drydiam

USE glomap_clim_option_mod,       ONLY: &
    l_glomap_clim_radaer

USE nlsizes_namelist_mod,         ONLY: &
    row_length,                         &
    rows,                               &
    model_levels,                       &
    bl_levels

USE parkind1,                     ONLY: &
    jpim,                               &
    jprb

USE planet_constants_mod,         ONLY: &
    repsilon

USE qsat_mod,                     ONLY: &
    qsat_wat_mix

USE stash_array_mod,              ONLY: &
    stash_maxlen,                       &
    sf,                                 &
    len_stlist,                         &
    num_stash_levels,                   &
    stlist,                             &
    si,                                 &
    stindex,                            &
    stash_levels

USE submodel_mod,                 ONLY: &
    atmos_im

USE ukca_activate_mod,            ONLY: &
    ukca_activate

USE ukca_cdnc_mod,                ONLY: &
    ukca_cdnc_struct

USE um_parvars,                   ONLY: &
    at_extremity

USE um_stashcode_mod,             ONLY: &
    stashcode_glomap_clim_sec,          &
    stashcode_gc_cdnc3,                 &
    stashcode_gc_cdnc

USE umPrintMgr,                   ONLY: &
    umPrint,                            &
    umMessage

USE yomhook,                      ONLY: &
    lhook,                              &
    dr_hook

IMPLICIT NONE

! Arguments

! Structure for UKCA/cdnc interaction
TYPE (ukca_cdnc_struct), INTENT(INOUT) :: ukca_cdnc

! stashwork54 array
REAL, INTENT(INOUT) :: stashwork54 ( stash_maxlen ( stashcode_glomap_clim_sec, &
                                                    atmos_im ) )

! Local variables

INTEGER, PARAMETER :: zero     = 0
INTEGER, PARAMETER :: one      = 1
INTEGER, PARAMETER :: im_index = 1  ! internal model index

INTEGER :: k, l     ! loop variables
INTEGER :: section  ! stash section
INTEGER :: item     ! stash item
INTEGER :: ireturn  ! error code from set_field
INTEGER :: ierrcode ! error code
INTEGER :: icode    ! return code

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=errormessagelength) :: cmessage

CHARACTER(LEN=*), PARAMETER :: RoutineName='GLOMAP_CLIM_ARG_ACT_GET_CDNC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

! Populate cdncfields structure
CALL glomap_clim_cdnc_identify_fields()

IF (.NOT. ALLOCATED(t_theta_levels))                                           &
                 ALLOCATE(t_theta_levels( row_length, rows, model_levels ) )
IF (.NOT. ALLOCATED(pr_theta_levels))                                          &
                ALLOCATE(pr_theta_levels( row_length, rows, model_levels ) )
IF (.NOT. ALLOCATED(qsmr)) ALLOCATE(qsmr( row_length, rows, model_levels ) )
IF (.NOT. ALLOCATED(qsvp)) ALLOCATE(qsvp( row_length, rows, model_levels ) )

t_theta_levels(:,:,:)= exner_theta_levels(1:row_length,1:rows,1:model_levels)* &
                                    theta(1:row_length,1:rows,1:model_levels)

pr_theta_levels(:,:,:)  =  p_theta_levels(1:row_length,1:rows,1:model_levels)

DO k = 1, model_levels
  CALL qsat_wat_mix( qsmr(:,:,k), t_theta_levels(:,:,k),                       &
                     pr_theta_levels(:,:,k), row_length, rows )
  
  ! Derive saturated vapour pressure from saturated mixing ratio
  qsvp(:,:,k) = (           qsmr(1:row_length,1:rows,k)  *                     &
                 pr_theta_levels(1:row_length,1:rows,k)) /                     &
             ( repsilon + ( qsmr(1:row_length,1:rows,k)  * (1.0 - repsilon) ) )
END DO

IF ( .NOT. l_glomap_clim_radaer ) THEN
  ! Populate drydiamfields structure
  CALL glomap_clim_drydiam_identify_fields()
  
  ! dry diameter will not have been set during radaer
  CALL glomap_clim_calc_drydiam( t_theta_levels, pr_theta_levels )
  
  ! the following allows drydiam to be written to stash
  IF (.NOT. ALLOCATED(outfield))                                               &
             ALLOCATE(outfield(row_length,rows,model_levels))
  
  stash_54_dd:DO l=1,n_fields_drydiam
    IF ( .NOT. drydiamfields(l)%put_stash ) CYCLE stash_54_dd
    
    item    = drydiamfields(l)%item
    section = drydiamfields(l)%section
    
    IF ( section /= stashcode_glomap_clim_sec ) CYCLE stash_54_dd
    
    IF ( .NOT. sf(item,section) ) CYCLE stash_54_dd
    
    CALL set_field_drydiam ( item, section, outfield, ireturn )
    
    IF (ireturn /= 0) THEN
      WRITE(umMessage,'(A,I4,A,I4)') 'item: ', item, ' section: ', section
      CALL umPrint(umMessage, src=RoutineName)
      
      ierrcode = 1
      cmessage = 'outfield contains negative values'
      CALL ereport(Modulename//':'//RoutineName,ierrcode,cmessage)
    END IF
    
    icode = 0
    
    ! Copy requested fields into STASHwork array
    
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d ( stashwork54(si(item,section,im_index)),                 &
                       outfield(:,:,:),                                        &
                       row_length, rows, model_levels,                         &
                       zero, zero, zero, zero, at_extremity,                   &
                       stlist( one, stindex(one, item, section, im_index) ),   &
                       len_stlist,                                             &
                       stash_levels,num_stash_levels+1,                        &
                       atmos_im,section,item,icode,cmessage)
    
    IF (icode >  0) THEN
      ierrcode = ABS( section*1000 + item )
      
      WRITE(umMessage,'(A,I6)')'Error from copydiag_3d, stashcode: ', ierrcode
      CALL umPrint( umMessage, src=RoutineName )
      
      CALL ereport( Modulename//':'//RoutineName, ierrcode, cmessage )
    END IF
  END DO stash_54_dd
  
  IF (ALLOCATED(outfield))        DEALLOCATE(outfield)
END IF

IF (.NOT. ALLOCATED(bl_tke)) ALLOCATE( bl_tke(row_length,rows,bl_levels) )

IF ( first_atmstep_call ) THEN
  ! bl_tke is uninitialised at time step one
  ! a minimum updraght velocity field is used instead of calculated from bl_tke
  bl_tke = 0.0
ELSE
  CALL glomap_clim_cdnc_bl_tke ( bl_tke )
END IF

IF (.NOT. ALLOCATED(cdnc))   ALLOCATE(  cdnc( row_length, rows, model_levels ))
IF (.NOT. ALLOCATED(cdnc3))  ALLOCATE( cdnc3( row_length, rows, model_levels ))
IF (.NOT. ALLOCATED(cdncwt)) ALLOCATE(cdncwt( row_length, rows, model_levels ))

IF (.NOT. ALLOCATED(mode_tracers))                                             &
           ALLOCATE(mode_tracers( row_length, rows, model_levels, one ) )

IF (.NOT. ALLOCATED(mode_diags))                                               &
           ALLOCATE(mode_diags(   row_length, rows, model_levels, one ) )

! note that vertical velocity (w) is required from the surface level
CALL ukca_activate ( row_length,                                               &
                     rows,                                                     &
                     bl_levels,                                                &
                     (row_length*rows),                                        &
                     one,                                                      &
                     one,                                                      &
                     pr_theta_levels(1:row_length,1:rows,1:model_levels),      &
                     t_theta_levels(1:row_length,1:rows,1:model_levels),       &
                     q(1:row_length,1:rows,1:model_levels),                    &
                     qsvp(1:row_length,1:rows,1:model_levels),                 &
                     bl_tke(1:row_length,1:rows,1:bl_levels),                  &
                     w(1:row_length,1:rows,0:model_levels-1),                  &
                     cf_liquid(1:row_length,1:rows,1:model_levels),            &
                     qcl(1:row_length,1:rows,1:model_levels),                  &
                     cdnc,                                                     &
                     cdnc3,                                                    &
                     cdncwt,                                                   &
                     mode_tracers,                                             &
                     mode_diags)

IF ( ALLOCATED(mode_diags) )   DEALLOCATE(mode_diags)
IF ( ALLOCATED(mode_tracers) ) DEALLOCATE(mode_tracers)

! Cloud droplet number concentration
ukca_cdnc%cdnc(:,:,:)  = cdnc(:,:,:)

! <Cloud droplet number concentration^-1/3>
ukca_cdnc%cdnc3(:,:,:) = cdnc3(:,:,:)

IF (.NOT. ALLOCATED(outfield)) ALLOCATE(outfield(row_length,rows,model_levels))

stash_54:DO l=1,n_fields_cdnc
  IF ( .NOT. cdncfields(l)%put_stash ) CYCLE stash_54
  
  item    = cdncfields(l)%item
  section = cdncfields(l)%section
  
  IF ( section /= stashcode_glomap_clim_sec ) CYCLE stash_54
  
  IF ( .NOT. sf(item,section) ) CYCLE stash_54
  
  CALL set_field_cdnc ( item, section, cdnc, cdnc3, cdncwt,                    &
                        cf_liquid(1:row_length,1:rows,1:model_levels),         &
                        outfield, ireturn )
  
  IF (ireturn /= 0) THEN
    WRITE(umMessage,'(A,I4,A,I4)') 'item: ', item, ' section: ', section
    CALL umPrint(umMessage, src=RoutineName)
    
    ierrcode = 1
    cmessage = 'outfield contains negative values'
    CALL ereport(Modulename//':'//RoutineName,ierrcode,cmessage)
  END IF
  
  icode = 0
  
  ! Copy requested fields into STASHwork array
  
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d ( stashwork54(si(item,section,im_index)),                   &
                     outfield(:,:,:),                                          &
                     row_length, rows, model_levels,                           &
                     zero, zero, zero, zero, at_extremity,                     &
                     stlist( one, stindex(one, item, section, im_index) ),     &
                     len_stlist,                                               &
                     stash_levels,num_stash_levels+1,                          &
                     atmos_im,section,item,icode,cmessage)
  
  IF (icode >  0) THEN
    ierrcode = ABS( section*1000 + item )
    
    WRITE(umMessage,'(A,I6)')'Error from copydiag_3d, stashcode: ', ierrcode
    CALL umPrint( umMessage, src=RoutineName )
    
    CALL ereport( Modulename//':'//RoutineName, ierrcode, cmessage )
  END IF
END DO stash_54

IF (ALLOCATED(outfield))        DEALLOCATE(outfield)
IF (ALLOCATED(cdncfields))      DEALLOCATE(cdncfields)
IF (ALLOCATED(cdncwt))          DEALLOCATE(cdncwt)
IF (ALLOCATED(cdnc3))           DEALLOCATE(cdnc3)
IF (ALLOCATED(cdnc))            DEALLOCATE(cdnc)
IF (ALLOCATED(bl_tke))          DEALLOCATE(bl_tke)
IF (ALLOCATED(qsvp))            DEALLOCATE(qsvp)
IF (ALLOCATED(qsmr))            DEALLOCATE(qsmr)
IF (ALLOCATED(pr_theta_levels)) DEALLOCATE(pr_theta_levels)
IF (ALLOCATED(t_theta_levels))  DEALLOCATE(t_theta_levels)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
END SUBROUTINE glomap_clim_arg_act_get_cdnc

SUBROUTINE glomap_clim_cdnc_identify_fields

USE glomap_clim_option_mod,  ONLY: &
    l_glomap_clim_arg_act,         &
    l_glomap_clim_aie1,            &
    l_glomap_clim_aie2

USE parkind1,                ONLY: &
    jprb,                          &
    jpim

USE um_stashcode_mod,        ONLY: &
    stashcode_glomap_clim_sec,     &
    stashcode_gc_cf_liquid,        &
    stashcode_gc_cdncwt,           &
    stashcode_gc_cdnc3,            &
    stashcode_gc_cdnc

USE yomhook,                 ONLY: &
    lhook,                         &
    dr_hook

IMPLICIT NONE

INTEGER, PARAMETER :: msect = 1000*stashcode_glomap_clim_sec  ! 1000 * section

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'GLOMAP_CLIM_CDNC_IDENTIFY_FIELDS'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! List of Section 54 items to be copied to stash
IF (.NOT. ALLOCATED(cdncfields)) ALLOCATE(cdncfields(n_fields_cdnc))

cdncfields(1)%section   = stashcode_glomap_clim_sec
cdncfields(1)%item      = stashcode_gc_cf_liquid - msect
! Only write cf_liquid to stash when l_glomap_clim_arg_act
IF (l_glomap_clim_arg_act) THEN
  cdncfields(1)%put_stash = .TRUE.
ELSE
  cdncfields(1)%put_stash = .FALSE.
END IF

cdncfields(2)%section   = stashcode_glomap_clim_sec
cdncfields(2)%item      = stashcode_gc_cdncwt - msect
! Only write cdncwt to stash when when l_glomap_clim_arg_act
IF (l_glomap_clim_arg_act) THEN
  cdncfields(2)%put_stash = .TRUE.
ELSE
  cdncfields(2)%put_stash = .FALSE.
END IF

cdncfields(3)%section   = stashcode_glomap_clim_sec
cdncfields(3)%item      = stashcode_gc_cdnc3 - msect
! Only write cdnc3 to stash when l_glomap_clim_arg_act
IF (l_glomap_clim_arg_act) THEN
  cdncfields(3)%put_stash = .TRUE.
ELSE
  cdncfields(3)%put_stash = .FALSE.
END IF

cdncfields(4)%section   = stashcode_glomap_clim_sec
cdncfields(4)%item      = stashcode_gc_cdnc - msect
! Only write cdnc to stash when a required logical is set
IF (l_glomap_clim_arg_act .OR. l_glomap_clim_aie1 .OR. l_glomap_clim_aie2) THEN
  cdncfields(4)%put_stash = .TRUE.
ELSE
  cdncfields(4)%put_stash = .FALSE.
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE glomap_clim_cdnc_identify_fields

SUBROUTINE glomap_clim_drydiam_identify_fields

USE glomap_clim_option_mod,  ONLY: &
    l_glomap_clim_arg_act,         &
    l_glomap_clim_radaer

USE parkind1,                ONLY: &
    jprb,                          &
    jpim

USE um_stashcode_mod,        ONLY: &
    stashcode_glomap_clim_sec,     &
    stashcode_gc_dryd_ait_sol,     &
    stashcode_gc_dryd_acc_sol,     &
    stashcode_gc_dryd_cor_sol,     &
    stashcode_gc_dryd_ait_ins,     &
    stashcode_gc_dryd_acc_ins,     &
    stashcode_gc_dryd_cor_ins

USE yomhook,                 ONLY: &
    lhook,                         &
    dr_hook

IMPLICIT NONE

INTEGER, PARAMETER :: msect = 1000*stashcode_glomap_clim_sec  ! 1000 * section

CHARACTER(LEN=*), PARAMETER :: RoutineName='GLOMAP_CLIM_DRYDIAM_IDENTIFY_FIELDS'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! List of Section 54 items to be copied to stash
IF (.NOT. ALLOCATED(drydiamfields)) ALLOCATE(drydiamfields(n_fields_drydiam))

drydiamfields(1)%section   = stashcode_glomap_clim_sec
drydiamfields(1)%item      = stashcode_gc_dryd_ait_sol - msect
!Write dryd_ait_sol to stash when l_glomap_clim_arg_act not l_glomap_clim_radaer
IF (l_glomap_clim_arg_act .AND. (.NOT. l_glomap_clim_radaer)) THEN
  drydiamfields(1)%put_stash = .TRUE.
ELSE
  drydiamfields(1)%put_stash = .FALSE.
END IF

drydiamfields(2)%section   = stashcode_glomap_clim_sec
drydiamfields(2)%item      = stashcode_gc_dryd_acc_sol - msect
!Write dryd_acc_sol to stash when l_glomap_clim_arg_act not l_glomap_clim_radaer
IF (l_glomap_clim_arg_act .AND. (.NOT. l_glomap_clim_radaer)) THEN
  drydiamfields(2)%put_stash = .TRUE.
ELSE
  drydiamfields(2)%put_stash = .FALSE.
END IF

drydiamfields(3)%section   = stashcode_glomap_clim_sec
drydiamfields(3)%item      = stashcode_gc_dryd_cor_sol - msect
!Write dryd_cor_sol to stash when l_glomap_clim_arg_act not l_glomap_clim_radaer
IF (l_glomap_clim_arg_act .AND. (.NOT. l_glomap_clim_radaer)) THEN
  drydiamfields(3)%put_stash = .TRUE.
ELSE
  drydiamfields(3)%put_stash = .FALSE.
END IF

drydiamfields(4)%section   = stashcode_glomap_clim_sec
drydiamfields(4)%item      = stashcode_gc_dryd_ait_ins - msect
!Write dryd_ait_ins to stash when l_glomap_clim_arg_act not l_glomap_clim_radaer
IF (l_glomap_clim_arg_act .AND. (.NOT. l_glomap_clim_radaer)) THEN
  drydiamfields(4)%put_stash = .TRUE.
ELSE
  drydiamfields(4)%put_stash = .FALSE.
END IF

drydiamfields(5)%section   = stashcode_glomap_clim_sec
drydiamfields(5)%item      = stashcode_gc_dryd_acc_ins - msect
!Write dryd_acc_ins to stash when l_glomap_clim_arg_act not l_glomap_clim_radaer
IF (l_glomap_clim_arg_act .AND. (.NOT. l_glomap_clim_radaer)) THEN
  drydiamfields(5)%put_stash = .TRUE.
ELSE
  drydiamfields(5)%put_stash = .FALSE.
END IF

drydiamfields(6)%section   = stashcode_glomap_clim_sec
drydiamfields(6)%item      = stashcode_gc_dryd_cor_ins - msect
!Write dryd_cor_ins to stash when l_glomap_clim_arg_act not l_glomap_clim_radaer
IF (l_glomap_clim_arg_act .AND. (.NOT. l_glomap_clim_radaer)) THEN
  drydiamfields(6)%put_stash = .TRUE.
ELSE
  drydiamfields(6)%put_stash = .FALSE.
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE glomap_clim_drydiam_identify_fields

SUBROUTINE set_field_cdnc ( item, section, cdnc, cdnc3, cdncwt, gc_cf_liquid,  &
                            outfield, ireturn)

! To fill outfield array with correct 3D field using section, item

USE ereport_mod,             ONLY: &
    ereport

USE errormessagelength_mod,  ONLY: &
    errormessagelength

USE nlsizes_namelist_mod,    ONLY: &
    row_length,                    &
    rows,                          &
    model_levels

USE parkind1,                ONLY: &
    jprb,                          &
    jpim

USE um_stashcode_mod,        ONLY: &
    stashcode_glomap_clim_sec,     &
    stashcode_gc_cf_liquid,        &
    stashcode_gc_cdncwt,           &
    stashcode_gc_cdnc3,            &
    stashcode_gc_cdnc

USE yomhook,                 ONLY: &
    lhook,                         &
    dr_hook

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN)  :: item                                   ! Stash item
INTEGER, INTENT(IN)  :: section                                ! Stash section
REAL,    INTENT(IN)  :: cdnc(row_length,rows,model_levels)     ! cdnc
REAL,    INTENT(IN)  :: cdnc3(row_length,rows,model_levels)    ! cdnc3
REAL,    INTENT(IN)  :: cdncwt(row_length,rows,model_levels)   ! cdncwt
REAL,    INTENT(IN)  :: gc_cf_liquid(row_length,rows,model_levels) ! cf_liquid
REAL,    INTENT(OUT) :: outfield(row_length,rows,model_levels) ! output field
INTEGER, INTENT(OUT) :: ireturn                                ! success/fail

! Local variables

INTEGER, PARAMETER                :: msect = stashcode_glomap_clim_sec * 1000

INTEGER                           :: ierrcode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER       :: RoutineName='SET_FIELD_CDNC'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

outfield = -0.1
ireturn  = 0

SELECT CASE (item)
! cf_liquid
CASE ( stashcode_gc_cf_liquid - msect )
  outfield = gc_cf_liquid
! cdncwt
CASE ( stashcode_gc_cdncwt - msect )
  outfield = cdncwt
! cdnc3
CASE ( stashcode_gc_cdnc3 - msect )
  outfield = cdnc3
! cdnc
CASE ( stashcode_gc_cdnc  - msect )
  outfield = cdnc
! Default
CASE DEFAULT
  cmessage = 'Item not found in case statement'
  ierrcode = ABS(item)
  IF ( ierrcode == 0 ) THEN
    ierrcode = 54000
  END IF
  CALL ereport(ModuleName//':'//RoutineName,ierrcode,cmessage)
END SELECT

! Set flag if field not filled
IF (ANY(outfield < 0.0)) ireturn = 1

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_field_cdnc

SUBROUTINE set_field_drydiam ( item, section, outfield, ireturn)

! To fill outfield array with correct 3D field using section, item

USE ereport_mod,             ONLY: &
    ereport

USE errormessagelength_mod,  ONLY: &
    errormessagelength

USE nlsizes_namelist_mod,    ONLY: &
    row_length,                    &
    rows,                          &
    model_levels

USE parkind1,                ONLY: &
    jprb,                          &
    jpim

USE um_stashcode_mod,        ONLY: &
    stashcode_glomap_clim_sec,     &
    stashcode_gc_dryd_ait_sol,     &
    stashcode_gc_dryd_acc_sol,     &
    stashcode_gc_dryd_cor_sol,     &
    stashcode_gc_dryd_ait_ins,     &
    stashcode_gc_dryd_acc_ins,     &
    stashcode_gc_dryd_cor_ins

USE ukca_aero_ctl_mod,       ONLY: &
    drydiam

USE ukca_mode_setup,         ONLY: &
    mode_ait_sol,                  &
    mode_acc_sol,                  &
    mode_cor_sol,                  &
    mode_ait_insol,                &
    mode_acc_insol,                &
    mode_cor_insol

USE yomhook,                 ONLY: &
    lhook,                         &
    dr_hook

IMPLICIT NONE

! Arguments

INTEGER, INTENT(IN)  :: item                                   ! Stash item
INTEGER, INTENT(IN)  :: section                                ! Stash section
REAL,    INTENT(OUT) :: outfield(row_length,rows,model_levels) ! output field
INTEGER, INTENT(OUT) :: ireturn                                ! success/fail

! Local variables

INTEGER, PARAMETER                :: msect = stashcode_glomap_clim_sec * 1000

INTEGER                           :: ierrcode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER       :: RoutineName='SET_FIELD_DRYDIAM'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

outfield = -0.1
ireturn  = 0

SELECT CASE (item)
! dryd_ait_sol
CASE ( stashcode_gc_dryd_ait_sol - msect )
  outfield(:,:,:) = drydiam(:,:,:,mode_ait_sol)
! dryd_acc_sol
CASE ( stashcode_gc_dryd_acc_sol - msect )
  outfield(:,:,:) = drydiam(:,:,:,mode_acc_sol)
! dryd_cor_sol
CASE ( stashcode_gc_dryd_cor_sol - msect )
  outfield(:,:,:) = drydiam(:,:,:,mode_cor_sol)
! dryd_ait_ins
CASE ( stashcode_gc_dryd_ait_ins - msect )
  outfield(:,:,:) = drydiam(:,:,:,mode_ait_insol)
! dryd_acc_ins
CASE ( stashcode_gc_dryd_acc_ins - msect )
  outfield(:,:,:) = drydiam(:,:,:,mode_acc_insol)
! dryd_cor_ins
CASE ( stashcode_gc_dryd_cor_ins - msect )
  outfield(:,:,:) = drydiam(:,:,:,mode_cor_insol)
! Default
CASE DEFAULT
  cmessage = 'Item not found in case statement'
  ierrcode = ABS(item)
  IF ( ierrcode == 0 ) THEN
    ierrcode = 54000
  END IF
  CALL ereport(ModuleName//':'//RoutineName,ierrcode,cmessage)
END SELECT

! Set flag if field not filled
IF (ANY(outfield < 0.0)) ireturn = 1

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_field_drydiam

SUBROUTINE glomap_clim_cdnc_bl_tke ( bl_tke )

USE d1_array_mod,           ONLY: &
    d1,                           &
    d1_addr,                      &
    no_obj_d1,                    &
    d1_item,                      &
    d1_section,                   &
    d1_address,                   &
    d1_length,                    &
    d1_no_levels

USE ereport_mod,            ONLY: &
    ereport

USE errormessagelength_mod, ONLY: &
    errormessagelength

USE nlsizes_namelist_mod,   ONLY: &
    row_length,                   &
    rows,                         &
    bl_levels

USE parkind1,               ONLY: &
    jpim,                         &
    jprb

USE submodel_mod,           ONLY: &
    submodel_for_sm,              &
    atmos_im

USE umPrintMgr,             ONLY: &
    umPrint,                      &
    umMessage

USE um_stashcode_mod,       ONLY: &
    stashcode_bl_sec,             &
    stashcode_bl_tke

USE yomhook,                ONLY: &
    lhook,                        &
    dr_hook

IMPLICIT NONE

! Arguments

REAL, INTENT(OUT) :: bl_tke ( row_length, rows, bl_levels )

! Local variables

INTEGER :: i_obj          ! Counter
INTEGER :: m_atm_modl     ! Atmosphere submodel index

INTEGER :: section        ! section code
INTEGER :: item           ! item code
INTEGER :: n_levels       ! number of levels
INTEGER :: address        ! address in D1
INTEGER :: length         ! length of field

LOGICAL :: l_is_bl_tke

INTEGER                           :: ierr
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER       :: RoutineName='GLOMAP_CLIM_CDNC_BL_TKE'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

m_atm_modl = submodel_for_sm(atmos_im)

section   = stashcode_bl_sec 
item      = stashcode_bl_tke

l_is_bl_tke = .FALSE.

! Go through D1 and retain some information
location_in_d1:DO i_obj=1,no_obj_d1(m_atm_modl)
  IF (      (  d1_addr(    d1_item,   i_obj, m_atm_modl ) == item )     .AND.  &
            (  d1_addr( d1_section,   i_obj, m_atm_modl ) == section ) ) THEN
    
    address  = d1_addr( d1_address,   i_obj, m_atm_modl )
    length   = d1_addr( d1_length,    i_obj, m_atm_modl )
    n_levels = d1_addr( d1_no_levels, i_obj, m_atm_modl )
    
    l_is_bl_tke = .TRUE.
    WRITE(umMessage,'(A)') 'Found bl_tke'
    CALL umPrint(umMessage,src=RoutineName)
    EXIT location_in_d1
  END IF
END DO location_in_d1

IF ( .NOT. l_is_bl_tke ) THEN
  ierr = 3473
  WRITE(umMessage,'(A)') 'bl_tke not found in d1. Include diagnostic in suite.'
  CALL umPrint(umMessage,src=RoutineName)
  cmessage = RoutineName // ' : bl_tke (03473) not found in d1. '
  CALL ereport( Modulename//':'//RoutineName, ierr, cmessage )
END IF

IF ( n_levels < bl_levels ) THEN
  ierr = 1
  WRITE(umMessage,'(A,I10)') 'n_levels = ', n_levels
  CALL umPrint(umMessage,src=RoutineName)
  WRITE(umMessage,'(A,I10)') 'bl_levels = ', bl_levels
  CALL umPrint(umMessage,src=RoutineName)
  
  cmessage = RoutineName // ' : n_levels < bl_levels'
  CALL ereport( Modulename//':'//RoutineName, ierr, cmessage )
END IF

buffer_size = row_length * rows * n_levels

IF ( length /= buffer_size ) THEN
  ierr = 1
  WRITE(umMessage,'(2(A,I20))') 'buffer_size ', buffer_size, ' length ', length
  CALL umPrint(umMessage,src=RoutineName)
  
  cmessage = RoutineName // ' : Unexpected total size of D1 diag.'
  CALL ereport( Modulename//':'//RoutineName, ierr, cmessage )
END IF

IF (.NOT. ALLOCATED(buffer)) ALLOCATE( buffer ( row_length, rows, n_levels ) )

buffer = RESHAPE( d1(address:address+length-1), (/row_length,rows,n_levels/) )

! levels above bl_levels initialised to zero
bl_tke = 0.0

bl_tke(:,:,1:bl_levels) =  buffer(:,:,1:bl_levels)

IF (ALLOCATED(buffer)) DEALLOCATE( buffer )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE glomap_clim_cdnc_bl_tke

END MODULE glomap_clim_cndc_mod
