! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  read in landsea mask and grid info and calculate cos(latitude)
MODULE crmstyle_read_landsea_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CRMSTYLE_READ_LANDSEA_MOD'

CONTAINS

SUBROUTINE crmstyle_read_landsea(in_cols, in_rows, mype, all_proc_group, &
                                 local_row_len,local_rows,               &
                                 bzy,bzx,bdy,bdx,pseudo_lat, pseudo_lon, &
                                 mask_hdr )

USE hires_data_mod , ONLY:                                              &
   landsea

USE crmwork_arrays_mod,  ONLY:         &
  xcoslat_full

! UM constants

USE conversions_mod,     ONLY: pi_over_180

! set_unit_bcast_flag etc
USE io

USE errormessagelength_mod, ONLY: errormessagelength

USE file_manager, ONLY: assign_file_unit

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type,          &
  UM_Header_type

USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal,            &
  EndofFile

USE UM_ParParams, ONLY: halo_type_no_halo, halo_type_single, fld_type_p

USE ereport_mod, ONLY: ereport, ereport_finalise

USE get_env_var_mod, ONLY: get_env_var

! Subroutine
USE get_anc_flds_mod, ONLY: get_anc_flds

USE umPrintMgr                              ! for writing output
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   Read Land sea mask ancillary file and grid info and calculate cos(latitude)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! ------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::   &
  in_cols                & ! Number of columns
 ,in_rows                & ! Number of rows
 ,mype                   & ! Processor number
 ,all_proc_group         & ! Group id for all processors
 ,local_row_len          & ! local row length
 ,local_rows               ! local rows

REAL, INTENT(OUT)  ::  &
  bzy                  & ! Full Grid info
 ,bzx                  &
 ,bdy                  &
 ,bdx                  &
 ,pseudo_lat           &
 ,pseudo_lon

TYPE(UM_Header_type),INTENT(OUT) :: mask_hdr       ! UM Headers:  land sea


!-------------------------------------------------------------------------
! Local variables
INTEGER ::               &
  i, j, k                  ! loop counters

INTEGER ::               &
  unit_mask              & ! unit number for land sea mask
 ,ErrorStatus            & ! Error code from operations on file
 ,numflds                & ! number of fields
 ,MaxFlds                & ! Maximum number of fields
 ,itest                    ! 0 - normal , 1 use test array to check MPP
                           ! decomposition works

INTEGER ::               &
  Store(1)               & !
 ,STCode(1)                ! stashcode to search for

LOGICAL :: PPHdrMod = .FALSE.
LOGICAL :: l_append = .FALSE.

TYPE(PP_Field_type)  :: maskfield(1)    ! land sea mask

REAL ::                                  &
  landsea_temp(local_row_len,local_rows) & ! fractional landsea
 ,landsea_full(in_cols,in_rows)          & ! full landsea fraction
 ,latitude(in_rows)                        ! full latitude

REAL, ALLOCATABLE :: &
  test_full(:,:)     & ! full test grid
 ,test_local(:,:)      ! local test grid data

CHARACTER(LEN=*), PARAMETER :: RoutineName = "CRMSTYLE_READ_LANDSEA"
CHARACTER(LEN=errormessagelength) :: cmessage  ! error message

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
itest = 0   ! no longer testing decompooistion works

ErrorStatus = StatusOK

CALL get_env_var("MASK_FILE", mask_hdr % FileName)

CALL assign_file_unit(mask_hdr % FileName, mask_hdr % UnitNum, handler="portio")

! The read_umdr reads in the information and sends to all PEs
! DEPENDS ON: read_umhdr
CALL read_umhdr( mask_hdr, l_append, ErrorStatus )

IF ( ErrorStatus == StatusOK ) THEN

  ! want grid size info from header

  ! Note all PE need to know this information
  bdx = mask_hdr%realc(1)
  bdy = mask_hdr%realc(2)
  bzx = mask_hdr%realc(3)
  bzy = mask_hdr%realc(4)
  pseudo_lat = mask_hdr%realc(5)
  pseudo_lon = mask_hdr%realc(6)
  WRITE(umMessage,'(A,6f10.5)') ' Full grid ',bzx, bdx, bzy, bdy, pseudo_lat, &
                        pseudo_lon
  CALL umPrint(umMessage,src=RoutineName)

ELSE ! problems reading file

  WRITE(umMessage,'(A,I10)') 'Problems reading header ',ErrorStatus
  CALL umPrint(umMessage,src=RoutineName)
  ErrorStatus = StatusFatal
  CALL EReport( RoutineName, ErrorStatus,"Error Reading land fraction Header" )

END IF

maxflds = 1
NumFlds       =       1
Store     (1) =       1
STCode    (1) =      505 ! fractional land

CALL set_unit_bcast_flag(mask_hdr % UnitNum)
IF (mype == 0) THEN
  CALL get_anc_flds( NumFlds, MaxFlds, STCode,Store, PPHdrMod, mask_hdr,     &
                   bzy,bzx,bdy,bdx, maskfield, ErrorStatus )
  IF ( ErrorStatus /= StatusOK ) THEN
    WRITE(umMessage,'(A,I10)') ' failed to read landsea ',ErrorStatus
    CALL umPrint(umMessage,src=RoutineName)
    ErrorStatus = StatusFatal
    CALL EReport( RoutineName, ErrorStatus,              &
                 "Error Reading land sea fraction field from ancillary" )
  ELSE
    WRITE(umMessage,'(A)') ' Read in Land sea mask'
    CALL umPrint(umMessage,src=RoutineName)
  END IF

!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                             &
!$OMP& SHARED(in_rows, in_cols, latitude, bzy, bdy, xcoslat_full,        &
!$OMP& landsea_full, maskField )
  DO j=1,in_rows

    ! latitude of grid points - full grid
    latitude(j) = bzy + REAL(j)*bdy

    DO i= 1,in_cols
      landsea_full(i,j) = maskField(1) % RData(i,j)
      xcoslat_full(i,j) = COS(latitude(j)*pi_over_180)
    END DO
  END DO
!$OMP END PARALLEL DO

    ! Release space
  DEALLOCATE( maskfield(1) % RData )
  NULLIFY( maskField(1) % RData )

END IF   ! mype == 0
CALL clear_unit_bcast_flag(mask_hdr % UnitNum)

! Full land sea mask on PE 0.  Want to scatter back to all PEs

! Note scatter field only fills the fields values not any halo values so cannot
! be used to send a field with halo values.

! DEPENDS ON: scatter_field
CALL scatter_field(landsea_temp,landsea_full,                     &
                  local_row_len,local_rows,                       &
                  in_cols,in_rows,                                &
                  fld_type_p,halo_type_no_halo,                   &
                  0,all_proc_group)

! 64 bit to 32 bit

!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                             &
!$OMP& SHARED(local_rows, local_row_len, landsea, landsea_temp )
DO j= 1,local_rows
  DO i= 1,local_row_len
    landsea(i,j) =  landsea_temp(i,j)
  END DO
END DO
!$OMP END PARALLEL DO

!------------------------------------------------------------------------------
! Test decompostion of fields and scattering
!------------------------------------------------------------------------------

IF (itest == 1) THEN
  ALLOCATE(test_full(in_cols,in_rows) )
  ALLOCATE(test_local(local_row_len,local_rows))
  ! test field
  DO j= 1,in_rows
    DO i= 1,in_cols
      test_full(i,j) = REAL(j)*10000.0+REAL(i)
    END DO
  END DO

  ! DEPENDS ON: scatter_field
  CALL scatter_field(test_local,test_full,                        &
                  local_row_len,local_rows,                       &
                  in_cols,in_rows,                                &
                  fld_type_p,halo_type_no_halo,                   &
                  0,all_proc_group)

  WRITE(umMessage,'(A,I10)') ' test no halo',mype
  CALL umPrint(umMessage,src=RoutineName)
  DO j = 1,local_rows
    WRITE(umMessage,'(1200f10.0)') (test_local(i,j), i=1,local_row_len)
    CALL umPrint(umMessage,src=RoutineName)
  END DO

  DEALLOCATE(test_local)
  DEALLOCATE(test_full)
END IF   ! test decomposition

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE crmstyle_read_landsea

END MODULE crmstyle_read_landsea_mod
