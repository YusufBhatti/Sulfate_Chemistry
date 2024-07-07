! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
!  gather and copy a 2d or 3d field into a set of pp structure for writeout

MODULE gather_cp_to_pp_struct_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GATHER_CP_TO_PP_STRUCT_MOD'

CONTAINS

SUBROUTINE gather_cp_to_pp_struct(stashcode,lbproc,cols,rows,cols_out,rows_out,&
                             mlevs,maxfldsout,                                 &
                             infield,pp_struct, um_hdr, ErrorStatus )


USE word_sizes_mod, ONLY: iwp,wp

USE UM_ParVars, ONLY: gc_all_proc_group
USE UM_ParParams, ONLY: halo_type_no_halo, halo_type_single, fld_type_p
USE io

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type,          &
  UM_Header_type

USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal,            &
  EndofFile

USE ereport_mod, ONLY: ereport, ereport_finalise
USE errormessagelength_mod, ONLY: errormessagelength

USE umPrintMgr                              ! for writing output
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
! Gather coarse resolution output field from all nodes/processors, copy to
! a pp structure and output field.
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
INTEGER, INTENT(IN)  ::  &
  stashcode              &  ! stashcode for field
 ,lbproc                 &  ! processing code
 ,cols                   &  ! columns - local
 ,rows                   &  ! rows  -local
 ,cols_out               &  ! output columns
 ,rows_out               &  ! output rows
 ,mlevs                  &  ! number of model levels
 ,MaxFldsOut                ! Maximum size of lookup table for Fieldsfile

REAL(wp), INTENT(IN)  ::  &
  infield(cols,rows,mlevs)        ! input field for this processor

TYPE(PP_Field_type), INTENT(INOUT) :: &
  pp_struct(mlevs)                     ! pp fields on model levels

TYPE(UM_Header_type), INTENT(INOUT) :: um_hdr   ! Fields file header

INTEGER, INTENT(INOUT) :: ErrorStatus

!-------------------------------------------------------------------------
! Local variables
INTEGER ::               &
  i, j, k, ij              ! loop counters
INTEGER ::               &
  start_field_num
INTEGER ::               &
  loc_to_copy(mlevs)       ! locations of fields to writeout (all)

REAL  ::                                &
  packacc(mlevs)                        & ! not set as not using packing
 ,copy_infield(cols,rows)               & ! input field
 ,full_field(cols_out,rows_out,mlevs)     ! full output field

CHARACTER(LEN=*), PARAMETER :: RoutineName = "GATHER_CP_TO_PP_STRUCT"

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
! Gather field from all PEs for output

DO k=1,mlevs

!$OMP PARALLEL DO PRIVATE(i,j) DEFAULT(NONE)                 &
!$OMP& SHARED(k, rows, cols, copy_infield, infield)
  DO j=1,rows
    DO i=1,cols
      copy_infield(i,j) = infield(i,j,k)
    END DO
  END DO
!$OMP END PARALLEL DO

  ! DEPENDS ON: gather_field
  CALL gather_field(copy_infield,full_field(1,1,k),                  &
                  cols,rows,                                         &
                  cols_out,rows_out,                                 &
                  fld_type_p,halo_type_no_halo,                      &
                  0,gc_all_proc_group )


END DO

!------------------------------------------------------------------------------
! Current used fields in lookup table
start_field_num = um_hdr%numflds

!$OMP PARALLEL DO PRIVATE(i,j,ij,k) DEFAULT(NONE)                 &
!$OMP& SHARED(mlevs, rows_out, cols_out, pp_struct, full_field,   &
!$OMP&   start_field_num, loc_to_copy, stashcode, lbproc)

DO k=1,mlevs

  ! Set stashcode
  pp_struct(k) % Hdr % stcode   = stashcode

  ! Set processing code
  pp_struct(k) % Hdr % lbproc   = lbproc

  loc_to_copy(k) = k  ! locations in input pp structure
  pp_struct(k) %LookupPos = start_field_num +k

  ! copy field to rdata

  DO j=1,rows_out
    DO i=1,cols_out
      ij=(j-1)*cols_out+i
      pp_struct(k) % rdata(ij,1) = full_field(i,j,k)
    END DO
  END DO

END DO  ! level loop

!$OMP END PARALLEL DO

WRITE(umMessage,'(A9,I6,A14,I6)') ' Stash : ',stashcode,                  &
                                  ' total fields ',um_hdr%numflds
CALL umPrint(umMessage,src=RoutineName)

errorstatus = 0   ! no error before call

! writes 2d/3d fields to fields file
! DEPENDS ON: writeflds
CALL writeflds(mlevs,mlevs,maxfldsout,"NONE",loc_to_copy,packacc, &
                 pp_struct, um_hdr, ErrorStatus )

IF ( ErrorStatus /= StatusOK ) THEN
  ErrorStatus = StatusFatal
  CALL EReport( "gather_cp_to_pp_struct", ErrorStatus,    &
               "Error copying fields to output file" )
END IF
!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE gather_cp_to_pp_struct

END MODULE gather_cp_to_pp_struct_mod
