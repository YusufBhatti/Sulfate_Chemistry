! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  write out CRM style sampling to ff

MODULE crmstyle_write_crm_ff_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_WRITE_CRM_FF_MOD'

CONTAINS

SUBROUTINE crmstyle_write_crm_ff(icall_type,maxfldsout, example_hdr,        &
  ref_modlev_field,                                                         &
  crm_w, crm_u, crm_v, crm_th, crm_thv, crm_rho, crm_rh, crm_a,             &
  crm_dpx, crm_dpy, crm_q,   crm_qcl, crm_qcf, crm_qrain, crm_qgraup,       &
  crm_dt1, crm_dt2, crm_dt4, crm_dt9, crm_dt12, crm_dt30,                   &
  crm_dq4, crm_dq9, crm_dq12, crm_dq30,                                     &
  crm_dqcl4, crm_dqcl9, crm_dqcl12, crm_dqcl30,                             &
  crm_dqcf4, crm_dqcf3, crm_dqcf12, crm_dqcf30,                             &
  crm_dqrain30, crm_dqgr30, crm_drho,                                       &
  crm_thw, crm_thvw, crm_qw, crm_qclw, crm_qcfw, crm_qrainw, crm_qgraupw,   &
  crm_uw, crm_vw, crm_ww, crm_w3, crm_vv, crm_uu, crm_uv,                   &
  crm_uth, crm_vth, crm_uthv, crm_vthv, crm_uq, crm_vq, crm_wp,             &
  crm_hdr, error_code                )


USE crmstyle_grid_info_mod, ONLY:                                       &
  local_new_x,local_new_y

USE crmstyle_cntl_mod, ONLY:                                              &
  num_x,num_y, mlevs, l_sect30, l_qgraup

USE missing_data_mod, ONLY: rmdi, imdi
USE word_sizes_mod, ONLY: iwp,wp

USE io, ONLY: file_close

USE io_constants, ONLY: ioNoDelete

USE filenamelength_mod, ONLY: filenamelength

USE file_manager, ONLY: release_file_unit

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type,          &
  UM_Header_type,         &
  lenfixhd

USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal,            &
  EndofFile

USE gather_cp_to_pp_struct_mod, ONLY: gather_cp_to_pp_struct
USE ereport_mod, ONLY: ereport, ereport_finalise

USE umPrintMgr                              ! for writing output
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!  Writeout all the fields averaged for the gievn partition
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utility - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! ------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::    &
  icall_type              & ! call type
 ,MaxFldsout                ! max output fields


TYPE(UM_Header_type), INTENT(IN)     :: &
  example_hdr                             ! Header to copy
TYPE(PP_Field_type), INTENT(INOUT) ::      &
  ref_modlev_field(mlevs)                 ! pp fields on model levels


REAL(wp), INTENT(IN) ::                    &
  crm_a(local_new_x,local_new_y,mlevs)     &
 ,crm_w(local_new_x,local_new_y,mlevs)     &
 ,crm_u(local_new_x,local_new_y,mlevs)     &
 ,crm_v(local_new_x,local_new_y,mlevs)     &
 ,crm_th(local_new_x,local_new_y,mlevs)    &
 ,crm_thv(local_new_x,local_new_y,mlevs)   &
 ,crm_rh(local_new_x,local_new_y,mlevs)    &
 ,crm_rho(local_new_x,local_new_y,mlevs)   &
 ,crm_dpx(local_new_x,local_new_y,mlevs)   &
 ,crm_dpy(local_new_x,local_new_y,mlevs)   &
 ,crm_q(local_new_x,local_new_y,mlevs)      &
 ,crm_qcl(local_new_x,local_new_y,mlevs)    &
 ,crm_qcf(local_new_x,local_new_y,mlevs)    &
 ,crm_qrain(local_new_x,local_new_y,mlevs)  &
 ,crm_qgraup(local_new_x,local_new_y,mlevs) &
 ,crm_dt1(local_new_x,local_new_y,mlevs)    &
 ,crm_dt2(local_new_x,local_new_y,mlevs)    &
 ,crm_dt4(local_new_x,local_new_y,mlevs)    &
 ,crm_dt9(local_new_x,local_new_y,mlevs)    &
 ,crm_dt12(local_new_x,local_new_y,mlevs)   &
 ,crm_dt30(local_new_x,local_new_y,mlevs)   &
 ,crm_dq4(local_new_x,local_new_y,mlevs)    &
 ,crm_dq9(local_new_x,local_new_y,mlevs)    &
 ,crm_dq12(local_new_x,local_new_y,mlevs)   &
 ,crm_dq30(local_new_x,local_new_y,mlevs)   &
 ,crm_dqcl4(local_new_x,local_new_y,mlevs)    &
 ,crm_dqcl9(local_new_x,local_new_y,mlevs)    &
 ,crm_dqcl12(local_new_x,local_new_y,mlevs)   &
 ,crm_dqcl30(local_new_x,local_new_y,mlevs)   &
 ,crm_dqcf4(local_new_x,local_new_y,mlevs)    &
 ,crm_dqcf3(local_new_x,local_new_y,mlevs)    &
 ,crm_dqcf12(local_new_x,local_new_y,mlevs)   &
 ,crm_dqcf30(local_new_x,local_new_y,mlevs)   &
 ,crm_dqrain30(local_new_x,local_new_y,mlevs) &
 ,crm_dqgr30(local_new_x,local_new_y,mlevs)   &
 ,crm_drho(local_new_x,local_new_y,mlevs)     &
 ,crm_thw(local_new_x,local_new_y,mlevs)      &
 ,crm_thvw(local_new_x,local_new_y,mlevs)     &
 ,crm_qw(local_new_x,local_new_y,mlevs)       &
 ,crm_qclw(local_new_x,local_new_y,mlevs)     &
 ,crm_qcfw(local_new_x,local_new_y,mlevs)     &
 ,crm_qrainw(local_new_x,local_new_y,mlevs)   &
 ,crm_qgraupw(local_new_x,local_new_y,mlevs)  &
 ,crm_uw(local_new_x,local_new_y,mlevs)   &
 ,crm_vw(local_new_x,local_new_y,mlevs)   &
 ,crm_ww(local_new_x,local_new_y,mlevs)   &
 ,crm_w3(local_new_x,local_new_y,mlevs)   &
 ,crm_uu(local_new_x,local_new_y,mlevs)   &
 ,crm_vv(local_new_x,local_new_y,mlevs)   &
 ,crm_uv(local_new_x,local_new_y,mlevs)   &
 ,crm_uth(local_new_x,local_new_y,mlevs)  &
 ,crm_vth(local_new_x,local_new_y,mlevs)  &
 ,crm_uthv(local_new_x,local_new_y,mlevs) &
 ,crm_vthv(local_new_x,local_new_y,mlevs) &
 ,crm_uq(local_new_x,local_new_y,mlevs)   &
 ,crm_vq(local_new_x,local_new_y,mlevs)   &
 ,crm_wp(local_new_x,local_new_y,mlevs)

TYPE(UM_Header_type), INTENT(INOUT)     :: &
  crm_hdr                                   ! Header for output file
INTEGER, INTENT(OUT) ::    &
  error_code                 ! return code
  ! error_code is intent out, but is never checked in parent routine.
  ! Should it be deleted ?

!-------------------------------------------------------------------------
! Local variables
INTEGER ::               &
  i, j, k, ij,ii,jj                  ! loop counters

INTEGER ::               &
  ErrorStatus             ! Error code from operations on file


CHARACTER(LEN=*), PARAMETER :: RoutineName = "CRMSTYLE_WRITE_CRM_FF"

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
! Opens and sets up header
ErrorStatus = 0

error_code = 0 ! setting error_code prevents compiler warnings, but it isn't
               ! used. Should it be deleted ?

IF (icall_type == 1) THEN
  ! Setups new header and allocates lookup table
  ! DEPENDS ON: new_umhdr
  CALL new_umhdr( example_hdr, MaxFldsOut, crm_hdr, ErrorStatus )

  IF ( ErrorStatus /= StatusOK ) THEN
    ErrorStatus = StatusFatal
    CALL EReport( "crmstyle_write_crm_ff", ErrorStatus,    &
               "Error Setting up Output File header." )
  END IF

END IF

! Copy field to structure and set stashcode and processing code before
! writing to output fields file.

WRITE(umMessage,'(A,I4,2A)') ' Unit ',crm_hdr %UnitNum,                    &
                             ' File ',TRIM(crm_hdr % FileName)
CALL umPrint(umMessage,src=RoutineName)


CALL gather_cp_to_pp_struct( 4, 0, local_new_x,local_new_y,                &
                             num_x,num_y,mlevs,maxfldsout,                 &
                       crm_th,ref_modlev_field, crm_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(211, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                  &
                       crm_thv,ref_modlev_field, crm_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct( 2, 0, local_new_x,local_new_y,                &
                             num_x,num_y,mlevs,maxfldsout,                 &
                            crm_u,ref_modlev_field, crm_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct( 3, 0, local_new_x,local_new_y,                &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            crm_v,ref_modlev_field, crm_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(150, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            crm_w,ref_modlev_field, crm_hdr, ErrorStatus )


CALL gather_cp_to_pp_struct(10, 0, local_new_x,local_new_y,                &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            crm_q,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(254, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            crm_qcl,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(12, 0, local_new_x,local_new_y,                &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            crm_qcf,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(272, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            crm_qrain,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(273, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            crm_qgraup,ref_modlev_field, crm_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(30113, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            crm_rh,ref_modlev_field, crm_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(266, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            crm_a,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(253, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            crm_rho,ref_modlev_field, crm_hdr, ErrorStatus )

! Tendencies
CALL gather_cp_to_pp_struct(1181, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dt1,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(2181, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dt2,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(4181, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dt4,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(9181, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dt9,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(12181, 0, local_new_x,local_new_y,           &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dt12,ref_modlev_field, crm_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(4182, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dq4,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(9182, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dq9,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(12182, 0, local_new_x,local_new_y,           &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dq12,ref_modlev_field, crm_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(4183, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dqcl4,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(9183, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dqcl9,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(12183, 0, local_new_x,local_new_y,           &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dqcl12,ref_modlev_field, crm_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(4184, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dqcf4,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(3184, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dqcf3,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(12184, 0, local_new_x,local_new_y,           &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dqcf12,ref_modlev_field, crm_hdr, ErrorStatus )

IF (l_sect30) THEN
  CALL gather_cp_to_pp_struct(30181, 0, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dt30,ref_modlev_field, crm_hdr, ErrorStatus )
  CALL gather_cp_to_pp_struct(30182, 0, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dq30,ref_modlev_field, crm_hdr, ErrorStatus )
  CALL gather_cp_to_pp_struct(30183, 0, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dqcl30,ref_modlev_field, crm_hdr, ErrorStatus )
  CALL gather_cp_to_pp_struct(30184, 0, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dqcf30,ref_modlev_field, crm_hdr, ErrorStatus )
  CALL gather_cp_to_pp_struct(30188, 0, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_drho,ref_modlev_field, crm_hdr, ErrorStatus )
  CALL gather_cp_to_pp_struct(30189, 0, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dqrain30,ref_modlev_field, crm_hdr, ErrorStatus )

END IF

IF (l_sect30 .AND. l_qgraup) THEN
  CALL gather_cp_to_pp_struct(30190, 0, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                &
                       crm_dqgr30,ref_modlev_field, crm_hdr, ErrorStatus )
END IF

! Fluxes

CALL gather_cp_to_pp_struct(5292, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                            crm_thw,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5293, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                            crm_thvw,ref_modlev_field, crm_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(5290, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                            crm_qw,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5294, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                            crm_qclw,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5295, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                            crm_qcfw,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5296, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                            crm_qrainw,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5297, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                           crm_qgraupw,ref_modlev_field, crm_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(30233, 0, local_new_x,local_new_y,           &
                            num_x,num_y,mlevs,maxfldsout,                &
                            crm_ww,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30255, 0, local_new_x,local_new_y,           &
                            num_x,num_y,mlevs,maxfldsout,                &
                            crm_w3,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30213, 0, local_new_x,local_new_y,           &
                            num_x,num_y,mlevs,maxfldsout,                &
                            crm_uw,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30223, 0, local_new_x,local_new_y,           &
                            num_x,num_y,mlevs,maxfldsout,                &
                            crm_vw,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30211, 0, local_new_x,local_new_y,           &
                            num_x,num_y,mlevs,maxfldsout,                &
                            crm_uu,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30222, 0, local_new_x,local_new_y,           &
                            num_x,num_y,mlevs,maxfldsout,                &
                            crm_vv,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30212, 0, local_new_x,local_new_y,           &
                            num_x,num_y,mlevs,maxfldsout,                &
                            crm_uv,ref_modlev_field, crm_hdr, ErrorStatus )
! pressure gradients
CALL gather_cp_to_pp_struct(5001, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                            crm_dpx,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5002, 0, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                &
                            crm_dpy,ref_modlev_field, crm_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(30214, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            crm_uth,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30224, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            crm_vth,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30215, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            crm_uq,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30225, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            crm_vq,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5026, 0, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            crm_wp,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5027, 0, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            crm_uthv,ref_modlev_field, crm_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5028, 0, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            crm_vthv,ref_modlev_field, crm_hdr, ErrorStatus )

!Close output file

IF (icall_type == 3) THEN
  ! Write out updated lookup table
  ! DEPENDS ON: writelookup
  CALL WriteLookup( crm_hdr, ErrorStatus )
  WRITE(umMessage,'(I7,A)') crm_hdr % NumFlds, " fields written."
  CALL umPrint(umMessage,src=RoutineName)

  IF ( ASSOCIATED( crm_hdr % Lookup ) ) THEN
    CALL File_Close ( crm_hdr % UnitNum,               &
                      crm_hdr % FileName,              &
                      filenamelength,                  &
                      delete=ioNoDelete, error=ErrorStatus )

    CALL release_file_unit(crm_hdr % UnitNum, handler="portio")
  END IF

END IF
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_write_crm_ff

END MODULE crmstyle_write_crm_ff_mod

