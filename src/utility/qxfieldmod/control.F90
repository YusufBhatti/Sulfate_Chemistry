! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Small execs
SUBROUTINE control(ppunit1,ppunit2,len1_lookup,len2_lookup,       &
                   lookup,pp_inthd,len_inthd,                     &
                   pp_fixhd,len_fixhd,icode,cmessage,nent)

USE check_iostat_mod
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
INTEGER ::                                                        &
     len_fixhd                                                    &
    ,len_inthd                                                    &
    ,len_lookup                                                   &
    ,len1_lookup                                                  &
    ,len2_lookup                                                  &
    ,lookup(len1_lookup,len2_lookup)                              &
    ,looknew(len1_lookup,len2_lookup)                             &
    ,pp_fixhd(len_fixhd)                                          &
    ,pp_inthd(len_inthd)                                          &
    ,len_io                                                       &
    ,icode                                                        &
    ,ppunit1                                                      &
    ,ppunit2                                                      &
    ,nent
INTEGER ::                                                        &
     row_length                                                   &
    ,p_rows                                                       &
    ,p_field                                                      &
    ,lenbuf                                                       &
    ,i                                                            &
    ,j
REAL ::                                                           &
     a_io

INTEGER ::                                                        &
    stime_mod                                                     &
,   etime_mod                                                     &
,   nfields_mod                                                   &
,   mtype_mod(500)                                                &
,   mlevs_mod(500)                                                &
,   stime_sel                                                     &
,   etime_sel                                                     &
,   nfields_sel                                                   &
,   mtype_sel(500)                                                &
,   mlevs_sel(500)                                                &
,   stime_rej                                                     &
,   etime_rej                                                     &
,   nfields_rej                                                   &
,   mtype_rej(500)                                                &
,   mlevs_rej(500)                                                &
,   ppunit_orog                                                   &
,   stime_thi                                                     &
,   etime_thi                                                     &
,   nfields_thi                                                   &
,   mtype_thi(500)                                                &
,   mlevs_thi(500)                                                &
,   ixxstep_thi(500)                                              &
,   iyystep_thi(500)

INTEGER ::                                                        &
 ErrorStatus      ! Return code : 0 Normal Exit : >0 Error

CHARACTER ::                                                      &
    output_pack_type*6

REAL ::                                                           &
    amult(500)                                                    &
,   wind_10m_orog                                                 &
                           !  LEVEL ABOVE WHICH 10M WIND FIXED
,   wind_10m_scale         !  SCALE APPLIED TO LEVEL 1 WINDS

LOGICAL ::                                                        &
    modify                                                        &
,   reject                                                        &
,   SELECT                                                        &
,   wind_10m                                                      &
,   thin

NAMELIST /mods/                                                   &
  modify,stime_mod,etime_mod,nfields_mod,                         &
                                      mtype_mod,mlevs_mod,amult,  &
  SELECT,stime_sel,etime_sel,nfields_sel,mtype_sel,mlevs_sel,     &
  reject,stime_rej,etime_rej,nfields_rej,mtype_rej,mlevs_rej,     &
  wind_10m,wind_10m_scale,wind_10m_orog,ppunit_orog,              &
  thin,stime_thi,etime_thi,nfields_thi,mtype_thi,mlevs_thi,       &
                                        ixxstep_thi,iyystep_thi,  &
  output_pack_type

!-----------------------------------------------------------------------
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=errormessagelength) :: iomessage
!
! ---------------------------------------------------------------
!      init namelist
! ---------------------------------------------------------------
modify   = .FALSE.
reject   = .FALSE.
SELECT   = .FALSE.
wind_10m = .FALSE.
thin=.FALSE.
stime_mod = -99
etime_mod = -99
nfields_mod=0
stime_sel = -99
etime_sel = -99
nfields_sel=0
stime_rej = -99
etime_rej = -99
nfields_rej=0
stime_thi = -99
etime_thi = -99
nfields_thi=0
DO i=1,500
  mtype_mod(i)=0
  mlevs_mod(i)=0
  amult(i)=1.0
  mtype_sel(i)=0
  mlevs_sel(i)=0
  mtype_rej(i)=0
  mlevs_rej(i)=0
  mtype_thi(i)=0
  mlevs_thi(i)=0
  ixxstep_thi(i)=2
  iyystep_thi(i)=2
END DO
wind_10m_orog  = -9999.0
wind_10m_scale = .7
ppunit_orog    = 12
output_pack_type='WGDOS '
!
! ---------------------------------------------------------------
!      read namelist
! ---------------------------------------------------------------
READ (UNIT=5, NML=mods, IOSTAT=ErrorStatus, IOMSG=iomessage)
CALL check_iostat(errorstatus, "namelist MODS", iomessage)
WRITE(7,mods)

! ---------------------------------------------------------------
!      Set up constants
! ---------------------------------------------------------------
row_length=pp_inthd(6)
p_rows=pp_inthd(7)
p_field=row_length*p_rows
lenbuf=p_field + 512

! DEPENDS ON: fields
CALL fields(pp_fixhd,len_fixhd,lenbuf,p_field,                    &
             lookup,lookup,len1_lookup,len2_lookup,nent,          &
             stime_mod,etime_mod,nfields_mod,                     &
                                       mtype_mod,mlevs_mod,amult, &
             stime_sel,etime_sel,nfields_sel,mtype_sel,mlevs_sel, &
             stime_rej,etime_rej,nfields_rej,mtype_rej,mlevs_rej, &
             stime_thi,etime_thi,nfields_thi,mtype_thi,mlevs_thi, &
                                         ixxstep_thi,iyystep_thi, &
             modify,SELECT,reject,thin,output_pack_type,          &
             wind_10m,wind_10m_orog,wind_10m_scale,ppunit_orog,   &
             ppunit1,ppunit2,icode,cmessage)
9999 CONTINUE
RETURN
END SUBROUTINE control

