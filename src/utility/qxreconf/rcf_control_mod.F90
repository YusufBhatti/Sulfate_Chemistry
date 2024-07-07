! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  The main control routine for the reconfiguration

MODULE Rcf_Control_Mod

!  Subroutine Rcf_Control - main control routine
!
! Description:
!    Top level control of the reconfiguration
!
! Method:
!    Headers, Fields land-sea masks and interpolation weights are
!    set up. Then ancillaries are read and main dump creation done.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!    This file belongs in section: Reconfiguration
!
! Code Description:
!    Language: FORTRAN 90
!    This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_CONTROL_MOD'

CONTAINS

SUBROUTINE Rcf_Control( hdr_in, hdr_out )

USE Rcf_Setup_Header_Mod, ONLY: &
    Rcf_Setup_Header

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_Setup_Field_mod, ONLY: &
    Rcf_Setup_Field

USE Rcf_init_h_interp_mod, ONLY: &
    Rcf_init_h_interp

USE Rcf_ReadLSMIn_mod, ONLY: &
    Rcf_ReadLSMIn

USE Rcf_Setup_LSM_out_mod, ONLY: &
    Rcf_setup_lsm_out

USE Rcf_Grid_Type_Mod, ONLY: &
    Input_Grid,           &
    Output_Grid

USE Rcf_Ancil_mod, ONLY: &
    Rcf_Ancil

USE rcf_netcdf_atmos_mod, ONLY: &
    rcf_netcdf_atmos

USE rcf_netcdf_init_mod, ONLY: &
    rcf_netcdf_init

USE Rcf_Create_Dump_Mod, ONLY: &
    Rcf_Create_Dump

USE Rcf_WriteUMhdr_Mod, ONLY: &
    Rcf_WriteUMhdr

USE rcf_set_data_source_mod, ONLY: &
    rcf_initialise_data_source

USE rcf_data_source_mod, ONLY: &
    data_source_type

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (um_header_type), INTENT(INOUT)   :: hdr_in
TYPE (um_header_type), INTENT(INOUT)   :: hdr_out

! Local variables
TYPE (field_type), POINTER :: fields_in(:)
TYPE (field_type), POINTER :: fields_out(:)
INTEGER                    :: field_count_in
INTEGER                    :: field_count_out
CHARACTER (LEN=20)         :: title
TYPE (data_source_type), POINTER :: data_source(:)

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_CONTROL'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!---------------------------------------------------------------
! Initialise fields data-type(s) and field source array
!---------------------------------------------------------------
NULLIFY( fields_in  )
NULLIFY( fields_out )
NULLIFY( data_source )
!----------------------------------------------------------------
! Setup the output header information (including lookups)
!----------------------------------------------------------------
CALL Rcf_Setup_Header( hdr_in, hdr_out )

!---------------------------------------------------------------
! Write out the header
!---------------------------------------------------------------
CALL Rcf_WriteUMhdr( Hdr_Out )

!----------------------------------------------------------------
! Setup the field data-types for both the input and output grids
!----------------------------------------------------------------
title = 'Input grid'
CALL Rcf_Setup_Field( fields_in, Hdr_In, Input_Grid,                &
                      field_count_in, title)

title = 'Output grid'
CALL Rcf_Setup_Field( fields_out, Hdr_Out, Output_Grid,             &
                      field_count_out, title )

!------------------------------------------------------------------
! Setup horizontal interpolation weights
!------------------------------------------------------------------
CALL Rcf_init_h_interp( Input_Grid, Output_Grid, Hdr_In, Hdr_Out )

!-------------------------------------------------------------------
! Setup input and output land-sea masks and coastal adjustment
! indexes etc.
!-------------------------------------------------------------------
CALL Rcf_ReadLSMIn( Hdr_In, fields_in, field_count_in )
CALL Rcf_Setup_LSM_Out( hdr_out, fields_in, field_count_in, fields_out, &
                        field_count_out, Output_Grid )

!-------------------------------------------------------------------
! Setup data source
!-------------------------------------------------------------------
! Initialise the source of the data for output dump, further data source
! processing will be performed from within rcf_create_dump

CALL rcf_initialise_data_source(data_source, fields_in, fields_out,   &
                                field_count_in, field_count_out,      &
                                hdr_in, hdr_out)
!-------------------------------------------------------------------
! Do the ancillary processing
!-------------------------------------------------------------------
CALL Rcf_Ancil ( Hdr_In, Hdr_Out,              &
                 Fields_In, Field_Count_In,    &
                 Fields_Out, Field_Count_Out,  &
                 data_source )

!-------------------------------------------------------------------
! Do the NetCDF climatology processing
! Unlike ancil files there is no need to get unit number for NetCDF files.
!-------------------------------------------------------------------
CALL rcf_netcdf_init ( )

CALL rcf_netcdf_atmos ( hdr_out, fields_out, field_count_out )

!-------------------------------------------------------------------
! Create the output dump
!-------------------------------------------------------------------
CALL Rcf_Create_Dump( hdr_in, hdr_out, fields_in, fields_out,       &
                      field_count_in, field_count_out, data_source )


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Control
END MODULE Rcf_Control_Mod
