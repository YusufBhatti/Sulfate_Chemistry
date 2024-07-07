! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  write out coarse grid means and SD to ff

MODULE crmstyle_write_all_ff_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_WRITE_ALL_FF_MOD'

CONTAINS

SUBROUTINE crmstyle_write_all_ff(icall_type, maxfldsout, example_hdr,       &
                                 ref_modlev_field, all_hdr,  error_code  )

USE crmstyle_sample_arrays_mod, ONLY:                                       &
  all_w, all_u,all_v, all_t, all_th, all_thv, all_rho,                      &
  all_q, all_qcl, all_qcf, all_qrain, all_qgraup,                           &
  all_ptheta, all_rh, all_a, all_dpx, all_dpy,                              &
  all_dt1, all_dt2, all_dt4, all_dt9, all_dt12,all_dt30,                    &
  all_dq4, all_dq9, all_dq12, all_dq30,                                     &
  all_dqcl4, all_dqcl9, all_dqcl12, all_dqcl30,                             &
  all_dqcf4, all_dqcf3, all_dqcf12, all_dqcf30,                             &
  all_dqrain30, all_dqgr30, all_drho,                                       &
  all_thw, all_thvw,                                                        &
  all_qw, all_qclw, all_qcfw, all_qrainw, all_qgraupw,                      &
  all_uw, all_vw, all_ww, all_w3, all_vv,                                   &
  all_uu, all_uv,                                                           &
  all_uth, all_vth, all_uthv, all_vthv, all_uq, all_vq, all_wp,all_wp_hydro,&
  all_sd_w, all_sd_u, all_sd_v, all_sd_th, all_sd_thv, all_sd_t, all_sd_rh, &
  all_sd_rho, all_sd_ptheta, all_sd_q,   all_sd_qcl,                        &
  all_sd_qcf, all_sd_qrain, all_sd_qgraup,                                  &
  all_sd_dt1, all_sd_dt2, all_sd_dt4, all_sd_dt9, all_sd_dt12, all_sd_dt30, &
  all_sd_dq4, all_sd_dq9, all_sd_dq12, all_sd_dq30,                         &
  all_sd_dqcl4, all_sd_dqcl9, all_sd_dqcl12, all_sd_dqcl30,                 &
  all_sd_dqcf4, all_sd_dqcf3, all_sd_dqcf12, all_sd_dqcf30,                 &
  all_sd_dqrain30, all_sd_dqgr30, all_sd_drho

USE crmstyle_grid_info_mod, ONLY:                                          &
  local_new_x,local_new_y

USE crmstyle_cntl_mod, ONLY:                                                &
  num_x,num_y, mlevs, l_sect30, l_qgraup

USE word_sizes_mod, ONLY: iwp,wp

USE missing_data_mod, ONLY: rmdi, imdi

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
INTEGER, INTENT(IN) ::              &
  icall_type                        & ! call type
 ,MaxFldsout                          ! max output fields

TYPE(UM_Header_type), INTENT(IN) :: &
  example_hdr                         ! Header to copy

TYPE(PP_Field_type), INTENT(INOUT)  :: &
  ref_modlev_field(mlevs)             ! pp fields on model levels


TYPE(UM_Header_type), INTENT(INOUT)     :: &
  all_hdr                                    ! Header to for file

INTEGER, INTENT(OUT) ::    &
  error_code                 ! return code
  ! error_code is intent out, but is never checked in parent routine.
  ! Should it be deleted ?
!-------------------------------------------------------------------------
! Local variables
INTEGER ::               &
  i                        ! loop counters

INTEGER ::               &
  ErrorStatus             ! Error code from operations on file

CHARACTER(LEN=*), PARAMETER :: RoutineName = "CRMSTYLE_WRITE_ALL_FF"

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------


!=============================================================================
! Note currently outputting ~71 model level fields
!=============================================================================
error_code = 0 ! setting error_code prevents compiler warnings, but it isn't
               ! used. Should it be deleted ?

ErrorStatus = 0

IF (icall_type == 1) THEN
  ! Opens and sets up header
  ! Setups new header and allocates lookup table
  ! DEPENDS ON: new_umhdr
  CALL new_umhdr( example_hdr, MaxFldsOut, all_hdr, ErrorStatus )

  WRITE(umMessage,'(A,3I6)') ' setup lookup header ',all_hdr % NumFlds,     &
                 all_hdr % Len2Lookup, ErrorStatus
  CALL umPrint(umMessage,src=RoutineName)
  IF ( ErrorStatus /= StatusOK ) THEN
    ErrorStatus = StatusFatal
    CALL EReport( "crmstyle_write_all_ff", ErrorStatus,    &
               "Error Setting up Output File header." )
  END IF

END IF


! Copy field to structure and set stashcode and processing code before
! writing to output fields file.

WRITE(umMessage,'(A,I6,2A)') ' Unit ',all_hdr %UnitNum, ' File ',    &
                        TRIM(all_hdr % FileName)
CALL umPrint(umMessage,src=RoutineName)

CALL gather_cp_to_pp_struct( 4, 0, local_new_x,local_new_y,                 &
                             num_x,num_y,mlevs,maxfldsout,                  &
                             all_th,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(211, 0, local_new_x,local_new_y,                &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_thv,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(16004, 0, local_new_x,local_new_y,              &
                             num_x,num_y,mlevs,maxfldsout,                  &
                             all_t,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct( 2, 0, local_new_x,local_new_y,                 &
                             num_x,num_y,mlevs,maxfldsout,                  &
                             all_u,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct( 3, 0, local_new_x,local_new_y,                 &
                             num_x,num_y,mlevs,maxfldsout,                  &
                             all_v,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(150, 0, local_new_x,local_new_y,                &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_w,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(10, 0, local_new_x,local_new_y,                 &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_q,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(254, 0, local_new_x,local_new_y,                &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_qcl,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(12, 0, local_new_x,local_new_y,                 &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_qcf,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(272, 0, local_new_x,local_new_y,                &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_qrain,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(273, 0, local_new_x,local_new_y,                &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_qgraup,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(30113, 0, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_rh,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(408, 0, local_new_x,local_new_y,                &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_ptheta,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(266, 0, local_new_x,local_new_y,                &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_a, ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(253, 0, local_new_x,local_new_y,                &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_rho, ref_modlev_field, all_hdr, ErrorStatus )

! S.d. of fields
CALL gather_cp_to_pp_struct( 4, 512, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_th,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(211, 512, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_thv,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(16004, 512, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_t,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct( 2, 512, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_u,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct( 3, 512, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_v,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(150, 512, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_w,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(10, 512, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_q,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(254, 512, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_qcl,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(12, 512, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_qcf,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(272, 512, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_qrain,ref_modlev_field, all_hdr, ErrorStatus)
CALL gather_cp_to_pp_struct(273, 512, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                   &
                          all_sd_qgraup,ref_modlev_field, all_hdr, ErrorStatus)

CALL gather_cp_to_pp_struct(30113, 512, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_rh,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(408, 512, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                   &
                          all_sd_ptheta,ref_modlev_field, all_hdr, ErrorStatus)

CALL gather_cp_to_pp_struct(253, 512, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_rho, ref_modlev_field, all_hdr, ErrorStatus )

! Tendencies
CALL gather_cp_to_pp_struct(1181, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_dt1,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(2181, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_dt2,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(4181, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_dt4,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(9181, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_dt9,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(12181, 0, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_dt12,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(4182, 0, local_new_x,local_new_y,num_x,         &
                            num_y,mlevs,maxfldsout,                         &
                            all_dq4,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(9182, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                   &
                       all_dq9,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(12182, 0, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_dq12,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(4183, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_dqcl4,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(9183, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_dqcl9,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(12183, 0, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_dqcl12,ref_modlev_field, all_hdr, ErrorStatus)

CALL gather_cp_to_pp_struct(4184, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_dqcf4,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(3184, 0, local_new_x,local_new_y,               &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_dqcf3,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(12184, 0, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_dqcf12,ref_modlev_field, all_hdr, ErrorStatus )

! Tendencies s.d
CALL gather_cp_to_pp_struct(1181, 512, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_dt1,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(2181, 512, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_dt2,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(4181, 512, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_dt4,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(9181, 512, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_dt9,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(12181, 512, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_dt12,ref_modlev_field, all_hdr, ErrorStatus)

CALL gather_cp_to_pp_struct(4182, 512, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_dq4,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(9182, 512, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_dq9,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(12182, 512, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_dq12,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(4183, 512, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_dqcl4,ref_modlev_field, all_hdr, ErrorStatus)
CALL gather_cp_to_pp_struct(9183, 512, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_dqcl9,ref_modlev_field, all_hdr, ErrorStatus)
CALL gather_cp_to_pp_struct(12183, 512, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_dqcl12,ref_modlev_field,all_hdr,ErrorStatus)

CALL gather_cp_to_pp_struct(4184, 512, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_dqcf4,ref_modlev_field, all_hdr, ErrorStatus)
CALL gather_cp_to_pp_struct(3184, 512, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_dqcf3,ref_modlev_field, all_hdr, ErrorStatus)
CALL gather_cp_to_pp_struct(12184, 512, local_new_x,local_new_y,            &
                            num_x,num_y,mlevs,maxfldsout,                   &
                            all_sd_dqcf12,ref_modlev_field,all_hdr,ErrorStatus)


IF (l_sect30) THEN
  CALL gather_cp_to_pp_struct(30181, 0, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                &
                       all_dt30,ref_modlev_field, all_hdr, ErrorStatus )
  CALL gather_cp_to_pp_struct(30182, 0, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                &
                       all_dq30,ref_modlev_field, all_hdr, ErrorStatus )
  CALL gather_cp_to_pp_struct(30183, 0, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                &
                       all_dqcl30,ref_modlev_field, all_hdr, ErrorStatus )
  CALL gather_cp_to_pp_struct(30184, 0, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                &
                       all_dqcf30,ref_modlev_field, all_hdr, ErrorStatus )
  CALL gather_cp_to_pp_struct(30188, 0, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                &
                       all_drho,ref_modlev_field, all_hdr, ErrorStatus )
  CALL gather_cp_to_pp_struct(30189, 0, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                &
                       all_dqrain30,ref_modlev_field, all_hdr, ErrorStatus )

  CALL gather_cp_to_pp_struct(30181, 512, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                  &
                       all_sd_dt30,ref_modlev_field, all_hdr, ErrorStatus )
  CALL gather_cp_to_pp_struct(30182, 512, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                  &
                       all_sd_dq30,ref_modlev_field, all_hdr, ErrorStatus )
  CALL gather_cp_to_pp_struct(30183, 512, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                  &
                       all_sd_dqcl30,ref_modlev_field, all_hdr, ErrorStatus )
  CALL gather_cp_to_pp_struct(30184, 512, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                  &
                       all_sd_dqcf30,ref_modlev_field, all_hdr, ErrorStatus )
  CALL gather_cp_to_pp_struct(30188, 512, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                  &
                       all_sd_drho,ref_modlev_field, all_hdr, ErrorStatus )
  CALL gather_cp_to_pp_struct(30189, 512, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                  &
                     all_sd_dqrain30,ref_modlev_field, all_hdr, ErrorStatus )
END IF

IF (l_sect30 .AND. l_qgraup) THEN
  CALL gather_cp_to_pp_struct(30190, 0, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                &
                       all_dqgr30,ref_modlev_field, all_hdr, ErrorStatus )
  CALL gather_cp_to_pp_struct(30190, 512, local_new_x,local_new_y,         &
                            num_x,num_y,mlevs,maxfldsout,                  &
                       all_sd_dqgr30,ref_modlev_field, all_hdr, ErrorStatus )
END IF
! Fluxes

CALL gather_cp_to_pp_struct(5292, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                 &
                            all_thw,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5293, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                 &
                            all_thvw,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(5290, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                 &
                            all_qw,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5294, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                 &
                            all_qclw,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5295, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                 &
                            all_qcfw,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5296, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                 &
                            all_qrainw,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5297, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                 &
                            all_qgraupw,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(30233, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            all_ww,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30255, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            all_w3,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30213, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            all_uw,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30223, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            all_vw,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30211, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            all_uu,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30222, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            all_vv,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30212, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            all_uv,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(5001, 0, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            all_dpx,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5002, 0, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            all_dpy,ref_modlev_field, all_hdr, ErrorStatus )

CALL gather_cp_to_pp_struct(30214, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            all_uth,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30224, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            all_vth,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30215, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            all_uq,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(30225, 0, local_new_x,local_new_y,             &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            all_vq,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5026, 0, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                  &
                          all_wp_hydro,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5027, 0, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            all_uthv,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5028, 0, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            all_vthv,ref_modlev_field, all_hdr, ErrorStatus )
CALL gather_cp_to_pp_struct(5029, 0, local_new_x,local_new_y,              &
                            num_x,num_y,mlevs,maxfldsout,                  &
                            all_wp,ref_modlev_field, all_hdr, ErrorStatus )

! Close output file

IF (icall_type == 3) THEN
  ! Write out updated lookup table
  ! DEPENDS ON: writelookup
  CALL WriteLookup( all_hdr, ErrorStatus )
  WRITE(umMessage,'(I7,A)') all_hdr % NumFlds, " fields written."
  CALL umPrint(umMessage,src=RoutineName)

  IF ( ASSOCIATED( all_hdr % Lookup ) ) THEN
    CALL File_Close ( all_hdr % UnitNum,               &
                      all_hdr % FileName,              &
                      filenamelength,                  &
                      delete=ioNoDelete, error=ErrorStatus )

    CALL release_file_unit(all_hdr % UnitNum, handler="portio")
  END IF
END IF
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_write_all_ff

END MODULE crmstyle_write_all_ff_mod
