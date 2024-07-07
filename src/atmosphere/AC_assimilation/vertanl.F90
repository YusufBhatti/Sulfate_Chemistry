! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE VERTANL ------------------------------------------------
!
!    Purpose : High level vertical analysis routine.
!              Calls lower level VAN--- routines which are specific
!              to particular AC Observation types.
!
!     Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!     Project Task : P3
!
!
!
!     Called by: AC2
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: AC Assimilation
MODULE vertanl_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='VERTANL_MOD'

CONTAINS

SUBROUTINE vertanl (pstgs,                                        &
                    kact,ipass,lenobt,ndv,obs_no,                 &
                    model_obs_type,obs_lat,obs_long,obs,inobs,    &
                    exner,pstar,theta,rh,qcl,qcf,                 &
                    conv_cld,layer_cloud,pressure,                &
                    ls_rain,ls_snow,conv_rain,conv_snow,          &
                    rhcrit,                                       &
                    p_field,wtsfld,lenwts,nlevwt,                 &
                    cf1pt,cf2pt,cf3pt,cf4pt,                      &
                    np1pt,np2pt,np3pt,np4pt,                      &
                    obs_incr,normf,lmissd,                        &
                    bl_levels,                                    &
                    row_length,p_rows,                            &
                    lenob,nptobt,no_anal_levs,no_anal_var,        &
                    icode,cmessage)
!
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Atmos_Max_Sizes
USE UM_ParParams
USE comobs_mod, ONLY: nobtypmx, obs_info
USE getob3_mod, ONLY: getob3
USE vanmops_mixed_phase_mod, ONLY: vanmops_mixed_phase
USE vanrain_mod, ONLY: vanrain
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
USE ac_control_mod
USE errormessagelength_mod, ONLY: errormessagelength
USE nlsizes_namelist_mod, ONLY: model_levels

IMPLICIT NONE

INTEGER :: bl_levels,row_length,p_rows
INTEGER :: p_field,lenwts,nlevwt
INTEGER :: lenobt,ndv,inobs
INTEGER :: lenob,no_anal_levs,no_anal_var

INTEGER :: np1pt(lenobt),np2pt(lenobt),np3pt(lenobt)
INTEGER :: np4pt(lenobt)
INTEGER :: obs_no(lenobt)
INTEGER :: model_obs_type(lenobt)

REAL ::                                                           &
               exner    (p_field,model_levels),                   &
               pressure (p_field,model_levels),                   &
               layer_cloud(p_field,model_levels),                 &
               pstgs    (p_field),                                &
               pstar    (p_field),                                &
               theta    (p_field,model_levels),                   &
               rh       (p_field,model_levels),                   &
               qcl      (p_field,model_levels),                   &
               qcf      (p_field,model_levels),                   &
               conv_cld (p_field,model_levels),                   &
               ls_rain(p_field),                                  &
               ls_snow(p_field),                                  &
               conv_rain(p_field),                                &
               conv_snow(p_field),                                &
               wtsfld   (lenwts,nlevwt),                          &
               rhcrit(model_levels),                              &
               cf1pt(lenobt),cf2pt(lenobt),cf3pt(lenobt),         &
               cf4pt(lenobt),                                     &
               obs_incr(lenob+1,no_anal_levs,no_anal_var),        &
               normf(lenob+1,no_anal_levs),                       &
               obs_lat(lenobt),obs_long(lenobt),                  &
               obs(inobs,*)

LOGICAL :: lmissd (lenob+1,no_anal_levs)
!
INTEGER :: kact, ipass, ktype, nptobt, nlevob
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
!
!     DYNAMIC ALLOCATION
!
REAL ::        obdata(lenobt,ndv)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='VERTANL'

!
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (ltimer_ac) CALL timer ('VERTANL ',3)

!  ** 1       Get obs and error values using GETOB3
ktype  = lact(kact)
nlevob = obs_info % noblev(kact)

CALL getob3 (kact,obs,obdata,                                     &
             obs_lat,obs_long,                                    &
             lenobt,inobs,ndv,obs_no,icode,cmessage)
IF (icode >  0) GO TO 999
!
!  ** 2      Call appropriate vertical analysis routine VAN---
!     NB. ICODE is not checked here, because only one VAN routine is
!     called followed by a jump to the end of the routine.
!
IF (ktype == 406) THEN
  ! for MOPS cloud data

  IF (ipass == 1) THEN

    CALL vanmops_mixed_phase (                                 &
                 kact,ipass,theta,exner,conv_cld,              &
                 layer_cloud,pressure,                         &
                 rh,qcl,qcf,p_field,no_anal_levs,rhcrit,       &
                 obdata,cf1pt,cf2pt,cf3pt,cf4pt,               &
                 np1pt,np2pt,np3pt,np4pt,                      &
                 obs_incr,normf,obs_lat,obs_long,              &
                 obs_no,lmissd,nptobt,lenobt,ndv,lenob,        &
                 no_anal_levs,no_anal_var,                     &
                 bl_levels,icode,cmessage)
  ELSE IF (ipass == 2) THEN

    CALL vanmops_mixed_phase (                                 &
                 kact,ipass,theta,exner,conv_cld,              &
                 layer_cloud,pressure,                         &
                 wtsfld,qcl,qcf,lenwts,nlevwt,rhcrit,          &
                 obdata,cf1pt,cf2pt,cf3pt,cf4pt,               &
                 np1pt,np2pt,np3pt,np4pt,                      &
                 obs_incr,normf,obs_lat,obs_long,              &
                 obs_no,lmissd,nptobt,lenobt,ndv,lenob,        &
                 no_anal_levs,no_anal_var,                     &
                 bl_levels,icode,cmessage)
  END IF


ELSE IF (ktype == 506) THEN
  !
  IF (ipass == 1) THEN

    CALL vanrain (kact,ipass,ls_rain,                             &
      ls_snow,conv_rain,conv_snow,p_field,                        &
      obdata,cf1pt,cf2pt,cf3pt,cf4pt,                             &
      np1pt,np2pt,np3pt,np4pt,                                    &
      obs_incr,normf,obs_lat,obs_long,                            &
      obs_no,lmissd,nptobt,                                       &
      lenobt,ndv,lenob,                                           &
      no_anal_levs,no_anal_var,                                   &
      icode,cmessage)

  ELSE IF (ipass == 2) THEN

    CALL vanrain (kact,ipass,wtsfld,                              &
      wtsfld,wtsfld,wtsfld,lenwts,                                &
      obdata,cf1pt,cf2pt,cf3pt,cf4pt,                             &
      np1pt,np2pt,np3pt,np4pt,                                    &
      obs_incr,normf,obs_lat,obs_long,                            &
      obs_no,lmissd,nptobt,                                       &
      lenobt,ndv,lenob,                                           &
      no_anal_levs,no_anal_var,                                   &
      icode,cmessage)

  END IF


ELSE
  icode=1
  cmessage = 'VERTANL : Obs Type not known'
  WRITE(umMessage,*) 'VERTANL KTYPE not processed ',ktype
  CALL umPrint(umMessage,src='vertanl')
END IF
!
999  CONTINUE

IF (icode >  0) THEN
  WRITE(umMessage,*) 'VERTANL Error code',icode
  CALL umPrint(umMessage,src='vertanl')
  WRITE(umMessage,*) cmessage
  CALL umPrint(umMessage,src='vertanl')
END IF

IF (ltimer_ac) CALL timer ('VERTANL ',4)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE vertanl
END MODULE vertanl_mod
