! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  write out coarse grid means and SD to ff

MODULE crmstyle_write_single_ff_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='CRMSTYLE_WRITE_SINGLE_FF_MOD'

CONTAINS

SUBROUTINE crmstyle_write_single_ff(icall_type, maxfldsout, example_hdr,  &
  ref_modlev_field,                                                       &
  g1_mean_zh, g1_mean_sh, g1_mean_lh,  g1_mean_pstar, g1_mean_tstar,      &
  g1_mean_rain,                                                           &
  g1_mean_snow, g1_mean_precip, g1_mean_orog,  g1_mean_land,              &
  g1_sd_zh, g1_sd_sh, g1_sd_lh,  g1_sd_pstar, g1_sd_tstar, g1_sd_rain,    &
  g1_sd_snow, g1_sd_precip, g1_sd_orog,  g1_sd_land,                      &
  fract_conv, fract_strat, prec_conv, prec_strat, cape, cin, zlcl, zfree, &
  zneutral, bcu_pcape, bcu_dpcapedt, bcu_dpcapedt_bl,                     &
  bcw_pcape, bcw_dpcapedt, bcw_dpcapedt_bl, ppd_dpcapedt_bl,              &
  bcu_dilcape, bcu_dcapedt, bcw_dilcape, bcw_dcapedt,                     &
  single_hdr, error_code                )


USE crmstyle_grid_info_mod, ONLY:                                        &
  local_new_x, local_new_y

USE crmstyle_cntl_mod, ONLY:                                              &
  num_x,num_y, mlevs, l_class_col, l_cape, l_pcape, l_bcu, l_bcw, l_ppd

USE word_sizes_mod, ONLY: iwp,wp

USE UM_ParCore, ONLY: mype, nproc_max
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
! This file belongs in section: Utilty - crmstyle_coarse_grid
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
  ref_modlev_field(1)                      ! pp fields on model levels


REAL(wp), INTENT(IN) ::       &
  g1_mean_zh(local_new_x,local_new_y)     &
 ,g1_mean_sh(local_new_x,local_new_y)     &
 ,g1_mean_lh(local_new_x,local_new_y)     &
 ,g1_mean_pstar(local_new_x,local_new_y)  &
 ,g1_mean_tstar(local_new_x,local_new_y)  &
 ,g1_mean_rain(local_new_x,local_new_y)   &
 ,g1_mean_snow(local_new_x,local_new_y)   &
 ,g1_mean_precip(local_new_x,local_new_y) &
 ,g1_mean_orog(local_new_x,local_new_y)   &
 ,g1_mean_land(local_new_x,local_new_y)   &
 ,g1_sd_zh(local_new_x,local_new_y)       &
 ,g1_sd_sh(local_new_x,local_new_y)       &
 ,g1_sd_lh(local_new_x,local_new_y)       &
 ,g1_sd_pstar(local_new_x,local_new_y)    &
 ,g1_sd_tstar(local_new_x,local_new_y)    &
 ,g1_sd_rain(local_new_x,local_new_y)     &
 ,g1_sd_snow(local_new_x,local_new_y)     &
 ,g1_sd_precip(local_new_x,local_new_y)   &
 ,g1_sd_orog(local_new_x,local_new_y)     &
 ,g1_sd_land(local_new_x,local_new_y)     &
 ,fract_conv(local_new_x,local_new_y)     &
 ,fract_strat(local_new_x,local_new_y)    &
 ,prec_conv(local_new_x,local_new_y)      &
 ,prec_strat(local_new_x,local_new_y)     &
 ,cape(local_new_x,local_new_y)           &
 ,cin(local_new_x,local_new_y)            &
 ,zlcl(local_new_x,local_new_y)           &
 ,zfree(local_new_x,local_new_y)          &
 ,zneutral(local_new_x,local_new_y)       &
 ,bcu_pcape(local_new_x,local_new_y)      &
 ,bcu_dpcapedt(local_new_x,local_new_y)   &
 ,bcu_dpcapedt_bl(local_new_x,local_new_y)&
 ,bcw_pcape(local_new_x,local_new_y)      &
 ,bcw_dpcapedt(local_new_x,local_new_y)   &
 ,bcw_dpcapedt_bl(local_new_x,local_new_y)&
 ,ppd_dpcapedt_bl(local_new_x,local_new_y)&
 ,bcu_dilcape(local_new_x,local_new_y)    &
 ,bcu_dcapedt(local_new_x,local_new_y)    &
 ,bcw_dilcape(local_new_x,local_new_y)    &
 ,bcw_dcapedt(local_new_x,local_new_y)

TYPE(UM_Header_type), INTENT(INOUT)     :: &
  single_hdr                                  ! Header to for file

INTEGER, INTENT(OUT) ::    &
  error_code                 ! return code
  ! error_code is intent out, but is never checked in parent routine.
  ! Should it be deleted ?
!-------------------------------------------------------------------------
! Local variables
INTEGER ::               &
  i,j                      ! loop counters

INTEGER ::               &
  ErrorStatus              ! Error code from operations on file
INTEGER ::               &
  itest                    ! =1 setup and output test field

REAL(wp) ::                &
  test_field(local_new_x,local_new_y)    ! Use to check gathering working

CHARACTER(LEN=*), PARAMETER :: RoutineName = "CRMSTYLE_WRITE_SINGLE_FF"

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------


!=============================================================================
error_code = 0 ! setting error_code prevents compiler warnings, but it isn't
               ! used. Should it be deleted ?

ErrorStatus = 0
itest = 0      ! set to 1 if you wish to test gathering of field is working

IF (icall_type == 1) THEN
  ! Opens and sets up header
  ! Setups new header and allocates lookup table
  ! DEPENDS ON: new_umhdr
  CALL new_umhdr( example_hdr, MaxFldsOut, single_hdr,ErrorStatus)

  IF ( ErrorStatus /= StatusOK ) THEN
    ErrorStatus = StatusFatal
    CALL EReport( "crmstyle_write_single_ff", ErrorStatus,    &
               "Error Setting up Output File header." )
  END IF

END IF


! Setup test field to check gather 
IF (itest == 1) THEN

  DO j=1,local_new_y
    DO i=1,local_new_x
      test_field(i,j) = REAL(mype)*10.0+REAL((j-1)*local_new_x+i)
    END DO
  END DO

  CALL gather_cp_to_pp_struct(  1, 0, local_new_x, local_new_y,               &
                             num_x,num_y, 1, maxfldsout,                      &
                          test_field,ref_modlev_field, single_hdr,ErrorStatus)

END IF
! Copy field to structure and set stashcode and processing code before
! writing to output fields file.

WRITE(umMessage,'(A,I6,2A)') ' Unit ',single_hdr %UnitNum,                   &
                    ' File ',TRIM(single_hdr % FileName)
CALL umPrint(umMessage,src=RoutineName)


CALL gather_cp_to_pp_struct( 25, 0, local_new_x, local_new_y,                 &
                             num_x,num_y, 1, maxfldsout,                      &
                       g1_mean_zh,ref_modlev_field, single_hdr,ErrorStatus)
CALL gather_cp_to_pp_struct( 3217, 0, local_new_x, local_new_y,               &
                             num_x,num_y, 1, maxfldsout,                      &
                       g1_mean_sh,ref_modlev_field, single_hdr,ErrorStatus)
CALL gather_cp_to_pp_struct( 3234, 0, local_new_x, local_new_y,               &
                             num_x,num_y, 1, maxfldsout,                      &
                       g1_mean_lh,ref_modlev_field, single_hdr,ErrorStatus)
CALL gather_cp_to_pp_struct( 409, 0, local_new_x, local_new_y,                &
                             num_x,num_y, 1, maxfldsout,                      &
                       g1_mean_pstar,ref_modlev_field, single_hdr,ErrorStatus)
CALL gather_cp_to_pp_struct( 24, 0, local_new_x, local_new_y,                 &
                             num_x,num_y, 1, maxfldsout,                      &
                       g1_mean_tstar,ref_modlev_field, single_hdr,ErrorStatus)
CALL gather_cp_to_pp_struct( 4203, 0, local_new_x, local_new_y,               &
                             num_x,num_y, 1, maxfldsout,                      &
                       g1_mean_rain,ref_modlev_field, single_hdr,ErrorStatus)
CALL gather_cp_to_pp_struct( 4204, 0, local_new_x, local_new_y,               &
                             num_x,num_y, 1, maxfldsout,                      &
                       g1_mean_snow,ref_modlev_field, single_hdr,ErrorStatus)
CALL gather_cp_to_pp_struct( 5216, 0, local_new_x, local_new_y,               &
                             num_x,num_y, 1, maxfldsout,                      &
                       g1_mean_precip,ref_modlev_field, single_hdr, ErrorStatus)
CALL gather_cp_to_pp_struct( 33, 0, local_new_x, local_new_y,                 &
                             num_x,num_y, 1, maxfldsout,                      &
                       g1_mean_orog,ref_modlev_field, single_hdr,ErrorStatus)
CALL gather_cp_to_pp_struct( 505, 0, local_new_x, local_new_y,                &
                             num_x,num_y, 1, maxfldsout,                      &
                       g1_mean_land,ref_modlev_field, single_hdr,ErrorStatus)

! Tendencies s.d

CALL gather_cp_to_pp_struct( 25, 512, local_new_x, local_new_y,               &
                             num_x,num_y, 1, maxfldsout,                      &
                       g1_sd_zh,ref_modlev_field, single_hdr,ErrorStatus)
CALL gather_cp_to_pp_struct( 3217, 512, local_new_x, local_new_y,             &
                              num_x,num_y, 1, maxfldsout,                     &
                       g1_sd_sh,ref_modlev_field, single_hdr,ErrorStatus)
CALL gather_cp_to_pp_struct( 3234, 512, local_new_x, local_new_y,             &
                              num_x,num_y, 1, maxfldsout,                     &
                       g1_sd_lh,ref_modlev_field, single_hdr,ErrorStatus)
CALL gather_cp_to_pp_struct( 409, 512, local_new_x, local_new_y,              &
                              num_x,num_y, 1, maxfldsout,                     &
                       g1_sd_pstar,ref_modlev_field, single_hdr,ErrorStatus)
CALL gather_cp_to_pp_struct( 24, 512, local_new_x, local_new_y,               &
                              num_x,num_y, 1, maxfldsout,                     &
                       g1_sd_tstar,ref_modlev_field, single_hdr,ErrorStatus)
CALL gather_cp_to_pp_struct( 4203, 512, local_new_x, local_new_y,             &
                              num_x,num_y, 1, maxfldsout,                     &
                       g1_sd_rain,ref_modlev_field, single_hdr,ErrorStatus)
CALL gather_cp_to_pp_struct( 4204, 512, local_new_x, local_new_y,             &
                              num_x,num_y, 1, maxfldsout,                     &
                       g1_sd_snow,ref_modlev_field, single_hdr,ErrorStatus)
CALL gather_cp_to_pp_struct( 5216, 512, local_new_x, local_new_y,             &
                              num_x,num_y, 1, maxfldsout,                     &
                       g1_sd_precip,ref_modlev_field, single_hdr,ErrorStatus)
CALL gather_cp_to_pp_struct( 33, 512, local_new_x, local_new_y,               &
                             num_x,num_y, 1, maxfldsout,                      &
                       g1_sd_orog,ref_modlev_field, single_hdr,ErrorStatus)
CALL gather_cp_to_pp_struct( 505, 512, local_new_x, local_new_y,              &
                              num_x,num_y, 1, maxfldsout,                     &
                       g1_sd_land,ref_modlev_field, single_hdr,ErrorStatus)

IF (l_class_col) THEN
  ! Using code for deep convection to indicate convective
  CALL gather_cp_to_pp_struct(5269,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                       fract_conv,ref_modlev_field, single_hdr,ErrorStatus)
  ! Using code for mid convection to indicate stratiform
  CALL gather_cp_to_pp_struct(5272,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                       fract_strat,ref_modlev_field, single_hdr,ErrorStatus)

  ! Using code for deep precip to indicate convective
  CALL gather_cp_to_pp_struct(5277,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                       prec_conv,ref_modlev_field, single_hdr,ErrorStatus)
  ! Using code for mid convection to indicate stratiform
  CALL gather_cp_to_pp_struct(5279,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                       prec_strat,ref_modlev_field, single_hdr,ErrorStatus)

END IF

IF (l_cape .OR. l_pcape) THEN

  ! LCL  (like cloud base)
  CALL gather_cp_to_pp_struct(5210,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                       zlcl,ref_modlev_field, single_hdr,ErrorStatus)

END IF
! CAPE & CIN info

IF (l_cape) THEN
  ! CAPE - undilute
  CALL gather_cp_to_pp_struct(5233,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                       cape,ref_modlev_field, single_hdr,ErrorStatus)
  ! CIN undilute
  CALL gather_cp_to_pp_struct(5234,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                       cin,ref_modlev_field, single_hdr,ErrorStatus)

  ! Level of free convection use 5009
  CALL gather_cp_to_pp_struct(5009,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                       zfree,ref_modlev_field, single_hdr,ErrorStatus)
  ! zneutral (like cloud top)
  CALL gather_cp_to_pp_struct(5211,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                       zneutral,ref_modlev_field, single_hdr,ErrorStatus)

END IF

! BCu PCAPE and dPCAPEdt  - made up stash codes

IF (l_pcape .AND. l_bcu) THEN  
  ! PCAPE - dilute
  CALL gather_cp_to_pp_struct(5010,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                       bcu_pcape,ref_modlev_field, single_hdr,ErrorStatus)
  ! Rate of change of PCAPE over cloud depth
  CALL gather_cp_to_pp_struct(5011,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                       bcu_dpcapedt,ref_modlev_field, single_hdr,ErrorStatus)

  ! Rate of change of PCAPE below cloud 
  CALL gather_cp_to_pp_struct(5012,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                      bcu_dpcapedt_bl,ref_modlev_field, single_hdr,ErrorStatus)

  IF (l_ppd) THEN
    ! Rate of change of Tv in DD below cloud
    CALL gather_cp_to_pp_struct(5013,  0, local_new_x, local_new_y,           &
                             num_x,num_y, 1, maxfldsout,                      &
                      ppd_dpcapedt_bl,ref_modlev_field, single_hdr,ErrorStatus)
  END IF

  CALL gather_cp_to_pp_struct(5014,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                       bcu_dilcape,ref_modlev_field, single_hdr,ErrorStatus)
  ! Rate of change of CAPE over cloud depth
  CALL gather_cp_to_pp_struct(5015,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                       bcu_dcapedt,ref_modlev_field, single_hdr,ErrorStatus)
END IF

! BCw PCAPE and dPCAPEdt  - made up stash codes

IF (l_pcape .AND. l_bcw) THEN  
  ! PCAPE - dilute
  CALL gather_cp_to_pp_struct(5020,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                       bcw_pcape,ref_modlev_field, single_hdr,ErrorStatus)
  ! Rate of change of PCAPE over cloud depth
  CALL gather_cp_to_pp_struct(5021,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                       bcw_dpcapedt,ref_modlev_field, single_hdr,ErrorStatus)

  ! Rate of change of PCAPE below cloud 
  CALL gather_cp_to_pp_struct(5022,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                      bcw_dpcapedt_bl,ref_modlev_field, single_hdr,ErrorStatus)

  CALL gather_cp_to_pp_struct(5024,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                       bcw_dilcape,ref_modlev_field, single_hdr,ErrorStatus)
  ! Rate of change of CAPE over cloud depth
  CALL gather_cp_to_pp_struct(5025,  0, local_new_x, local_new_y,             &
                             num_x,num_y, 1, maxfldsout,                      &
                       bcw_dcapedt,ref_modlev_field, single_hdr,ErrorStatus)
END IF



! Close output file if last write

IF (icall_type == 3) THEN
  ! Write out updated lookup table
  ! DEPENDS ON: writelookup
  CALL WriteLookup( single_hdr,ErrorStatus)
  WRITE(umMessage,'(I7,A)') single_hdr % NumFlds, " fields written."
  CALL umPrint(umMessage,src=RoutineName)

  IF ( ASSOCIATED( single_hdr % Lookup ) ) THEN
    CALL File_Close ( single_hdr % UnitNum,               &
                      single_hdr % FileName,              &
                      filenamelength,                     &
                      delete=ioNoDelete, error=ErrorStatus)

    CALL release_file_unit(single_hdr % UnitNum, handler="portio")
  END IF
END IF
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_write_single_ff

END MODULE crmstyle_write_single_ff_mod
